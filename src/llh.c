#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "edible.h"
#include "new_models.h"


/*  Main routine to calculate the likelihood and
 * partial derivatives matrix of the tree*/
void calcllhs(struct treenode *node_p,struct treenode *tree2,unsigned int a, double **matrix){
  int b,c;

  extern int branches;
  extern int leaves;
  extern int mode;
  extern int nodecount;
  extern int individual;
  extern int interesting_branches[];
  extern FILE *prob_file_p;
  extern double kappa;
  extern int is_kappa;
  extern double **rootedexpect;
  
  double *partial;
  double *rootedpartial;
  double temp;
  double **cache_hit;
  
  /*  Makeleaf array point to tree starting at node_p*/
  leaves=0;
  doleaf(node_p,0);

  rootedpartial=NULL;
  
  if(ISMODE(PROBS)){
    if(ISMODE(ROOTED) && ISMODE(INDIVIDUAL) && individual==1){
      b=nodecount+((ISMODE(NODEASROOT))?1:2);
      rootedpartial=calloc(b,sizeof(double));
      if(rootedpartial==NULL)
	nomemory();
    }
  }
  
  /*  Code to work out nucleotide codes from the number.
   *  If we are using a boot strap, then we need to generate
   * the nucleotides using the tree evolve them. This is done
   * in the llhcalc routine - not here*/
  if(ISMODE(BOOTSTRAP) || ISMODE(PERCENTILE)){
    /*  "Evolve" the sequence we're going to calculate and
     * check to see if the result has been cached*/
    cache_hit=pcalc(node_p);
    /*  If we have a cache, copy the result, output information
     * to probability file if necessary and the return*/
    if(cache_hit!=NULL){
      if(ISMODE(HKY) && NOTMODE(NOKAPPA))
	is_kappa=1;
      for(b=0;b<branches+is_kappa+1;b++)
	for(c=0;c<branches+is_kappa+1;c++)
	  matrix[b][c]=cache_hit[b][c];
      is_kappa=0;
      
      /*  Dump necessary information to files*/
      if(ISMODE(PROBS)){
	fprintf(prob_file_p,"\t\t");
	fprintf(prob_file_p,"%e\t",matrix[branches][branches]);
	if(ISMODE(INDIVIDUAL) && individual==1){
	  if(ISMODE(ROOTED)){
	    planttree(matrix,rootedexpect);
	    rooted_derivative(matrix[branches],rootedpartial);
	    fprintf(prob_file_p,"%e\t%e\t",rootedpartial[interesting_branches[0]]*rootedpartial[interesting_branches[0]]
		    /(matrix[branches][branches]*matrix[branches][branches])
		    -rootedexpect[interesting_branches[0]][interesting_branches[0]]/matrix[branches][branches]
		    ,rootedpartial[interesting_branches[0]]);
	    fprintf(prob_file_p,"%e",rootedexpect[interesting_branches[0]][interesting_branches[0]]);
	  }
	  else{
	    fprintf(prob_file_p,"%e\t%e\t",matrix[branches][interesting_branches[0]]*matrix[branches][interesting_branches[0]]
		    /(matrix[branches][branches]*matrix[branches][branches])
		    -matrix[interesting_branches[0]][interesting_branches[0]]/matrix[branches][branches]
		    ,matrix[branches][interesting_branches[0]]);
	    fprintf(prob_file_p,"%e",matrix[interesting_branches[0]][interesting_branches[0]]);
	  }
	}
	fprintf(prob_file_p,"\n");
      }
      if(rootedpartial!=NULL)
	free(rootedpartial);
      
      return;
    }
  }
  /*  If not sampling, add sequence to the tree*/
  else
    tree_sequence(a);
  
  /*  Calculate likelihood and partial derivatives.
   * also generates a sequence in boot strap case*/
  matrix[branches][branches]=llhcalc(node_p);

  /*   Make a copy of the tree and calculate the partial 
   * derivatives for each branch and store in newly acquired array*/
  partial=calloc(branches+branches+1,sizeof(double));
  if(partial==NULL)
    nomemory();
  for(b=0;b<branches;b++){   
    filltree(node_p,tree2,0);
    partialcalc(tree2,partial,b);
  }

  /*  If we wish to calculate the information WRT kappa, then
   * it is more efficient to use FD methods than exact calculation*/
  if(ISMODE(HKY) && NOTMODE(NOKAPPA)){
    is_kappa=1;

    /*  Do the first evaluation llh(x+dx)*/
    kappa+=DELTA;
    do_rate();
    matrix[branches][branches+1]=llhcalc(node_p);
    matrix[branches+1][branches+1]=matrix[branches][branches+1];

    /*  The second evaluation llh(x-dx)*/
    kappa-=2*DELTA;
    do_rate();
    temp=llhcalc(node_p);

    /*  dllh/dx~[f(x+dx)-f(x-dx)]/(2*DELTA)*/
    matrix[branches][branches+1]-=temp;
    matrix[branches][branches+1]/=2*DELTA;
    
    /*  d2llh/dx2~[f(x+dx)+f(x-dx)-2*f(x)]/(DELTA*DELTA)*/
    matrix[branches+1][branches+1]+=temp;
    matrix[branches+1][branches+1]-=2*matrix[branches][branches];
    matrix[branches+1][branches+1]/=DELTA*DELTA;

    /*  Reset the rates and so on*/
    kappa+=DELTA;
    do_rate();
    is_kappa=0;
  }
  
  /*  Calculate the second derivatives. Unless we are sampling or wish
   * to calculate the variance, the off diagonal second partial 
   * derivatives are returned as 0*/
  for(b=0;b<branches;b++){
    for(c=0;c<(b+1);c++){
    matrix[b][c]=second_derivative(tree2,partial,b,c);
    /*  We have symmetry*/
    matrix[c][b]=matrix[b][c];
    }
    /*  Fill in the first order partial derivatives*/
    matrix[branches][b]=partial[b];
  }

  /*  Free the additional memory used*/
  free(partial);

  /*  Dump necessary information to files*/
  if(ISMODE(PROBS)){
    fprintf(prob_file_p,"\t\t");
    fprintf(prob_file_p,"%e\t",matrix[branches][branches]);
    if(ISMODE(INDIVIDUAL) && individual==1){
      if(ISMODE(ROOTED)){
        planttree(matrix,rootedexpect);
        rooted_derivative(matrix[branches],rootedpartial);
        fprintf(prob_file_p,"%e\t%e\t",rootedpartial[interesting_branches[0]]*rootedpartial[interesting_branches[0]]
                /(matrix[branches][branches]*matrix[branches][branches])
                -rootedexpect[interesting_branches[0]][interesting_branches[0]]/matrix[branches][branches]
                ,rootedpartial[interesting_branches[0]]);
        fprintf(prob_file_p,"%e",rootedexpect[interesting_branches[0]][interesting_branches[0]]);
      }
      else{
        fprintf(prob_file_p,"%e\t%e\t",matrix[branches][interesting_branches[0]]*matrix[branches][interesting_branches[0]]
	        /(matrix[branches][branches]*matrix[branches][branches])
	        -matrix[interesting_branches[0]][interesting_branches[0]]/matrix[branches][branches]
	        ,matrix[branches][interesting_branches[0]]);
        fprintf(prob_file_p,"%e",matrix[interesting_branches[0]][interesting_branches[0]]);
      }
    }
    fprintf(prob_file_p,"\n");
  }
  if(rootedpartial!=NULL)
    free(rootedpartial);
  
  /*  See if it's worth adding the result generated to the cache*/
  if((ISMODE(BOOTSTRAP) || ISMODE(PERCENTILE)) && ISMODE(CACHE))
    update_cache(matrix);
}


/*  Starts off calculation of likelihood of given
 * tree.*/
double llhcalc(struct treenode *node_p){

  int a,b,n,r;
  double l,llh, llhr;
  double tot;
  
  extern int mode;
  extern double p[4];
  extern int is_kappa;
  extern int nrates;
  extern double * rate_prob;
  
  assert(NULL!=node_p);
  assert(NULL!=rate_prob);
  assert(nrates>0);
  #ifndef NDEBUG
    tot = 0;
    for ( r=0 ; r<nrates ; r++){
      tot += rate_prob[r];
    }
    assert(fabs(1.-tot)<DBL_EPSILON*nrates);
  #endif
  
  
  llh = 0.;
  for ( r=0 ; r<nrates ; r++){
    /*  Function to recurse down the tree*/
    n=0;
    while(n<DOODAH && CHILD(node)!=NULL){
      llhcalc_branch(CHILD(node),r);
      n++;
    }

    a=n; llhr = 0;
    for(b=0;b<4;b++){/*  Loop to calculate the partial likelihood of this tree */
      l=1;           /* given the partial likelihoods already calculated.      */
      for(n=0;n<a;n++)
        l*=CHILD(node)->llh[b];
      llhr += p[b]*l;
    }

    llh += rate_prob[r] * llhr;    
  }
    
  assert(llh>=0. && llh<=1.);
  return llh;
}



/*  Given the nucleotides of the leaves of the tree,
 * this routine should calculate the likelihood of
 * a bit of the tree evolving the data*/

void llhcalc_branch(struct treenode *node_p, int r){
  extern int kudge;
  extern int branch1;
  extern int branch2;
  extern struct treenode *branch[];
  extern FILE *prob_file_p;
  extern int mode;
  extern int is_kappa;
  
  int a,n;

  
  /*  If we are on a leaf, do initial calculation of prob.*/
  n=1;
  if(CHILD(node)==NULL){
    for(a=0;a<4;a++){
      if(kudge==1 && (branch[branch1]==node_p || branch[branch2]==node_p)){
	node_p->llh[a]=prob_calcd_new(a,node_p->nucleotide,node_p->bnum,r);
      } else{
	node_p->llh[a]=prob_calc_new(a,node_p->nucleotide,node_p->bnum,r);
      }
    }
    return;
  }

/*  Now - can assume that we are on an internal node
 * Recurse down all the (non-parent) branches        */
  while(n<DOODAH && CHILD(node)!=NULL){
    llhcalc_branch(CHILD(node),r);
    n++;
  }
/*  We can assume the llh's of all the child nodes
 * have been calculated - and thus can calculate the
 * llh vector for the current node*/
  a=n;
/*  'a' is now either DOODAH, or the first 'branch'
 * with nothing coming off it                        */
  if(kudge==1 && (branch[branch1]==node_p || branch[branch2]==node_p))
    llhvector(node_p,a,1,r);
  else
    llhvector(node_p,a,0,r);
}

/*  Routine for just the calculation of the llh
 * vectors.
 *  Assumes internal node, with all the vectors
 * from child nodes already calculated.
 * a is the either DOODAH or the maximum number of
 * non-null branches from the node.*/
void llhvector(struct treenode *node_p,int a,int derivative, int r){
  int b,c,n;
  double llhmatrix[4][4];
  double l,llh;
  double (*prob_function)(int,int,int,int)=NULL;


  switch (derivative) {
    case 0:
      prob_function=prob_calc_new;
      break;
    case 1:
      prob_function=prob_calcd_new;
      break;
    case 2:
      prob_function=prob_calcd2_new;
      break;
    default:
      abort();
  }
/*  Matrix of the probabilities [x][y] for
 * nucleotide x to change to y
 * ** Nucleotide codes **   0 == A     1 == C
 *                          2 == G     3 == T*/

  for(b=0;b<4;b++)
    for(c=0;c<4;c++){
      llhmatrix[c][b]=prob_function(c,b,node_p->bnum,r);
    }
/*  Main code to calculate the partial likelihoods
 * and save in the llh vector at the node*/
  for(b=0;b<4;b++){
    llh=0;
    for(c=0;c<4;c++){
      l=llhmatrix[b][c];
      for(n=1;n<a;n++)
	l*=CHILD(node)->llh[c];
      llh+=l;
    }
    node_p->llh[b]=llh;
  }
}



/*  Routine to calculate probability of change along
 * length. Currently using JC model*/
double prob_calc(int from, int to, double length){
  extern int mode;
  extern double p[4];
  extern double rate;
  extern double kappa;

  double pi_j;
  double a,d;
    
  /*  Decide whether we are using HKY85 model or JC model of
   * nucleotide substitution*/
  if(NOTMODE(HKY))
    return (double)0.25+(DELTAF(from,to)-0.25)*exp(-FACTOR*length);
  else{
    if(to==0 || to==2)
      pi_j=p[0]+p[2];
    else
      pi_j=p[1]+p[3];
    a=1+pi_j*(kappa-1);
    
    /*  We have case from and to are the same nucleotide*/
    if(from==to){
      d=(1-p[to]/pi_j)*exp(-rate*length*a);
      d+=p[to]*(1-pi_j)*exp(-rate*length)/pi_j;
      d+=p[to];
    }
    /*  Case where a transition has occurred A to G, G to A, C to T to T to C*/
    else if(((from&1)==0 && (to&1)==0) || ((from&1)!=0 && (to&1)!=0)){
      d=p[to]*(1-pi_j)*exp(-rate*length)/pi_j;
      d-=p[to]*exp(-rate*length*a)/pi_j;
      d+=p[to];
    }
    /*  A transversion has occurred A to C or T, G to C or T or
     * there reversions in time*/
    else{
      d=p[to]*(1-exp(-rate*length));
    }
    return d;
  }
}  


/*  Routine to calculate probability of change along
 * length. Currently using JC model*/
double prob_calcd(int from, int to, double length){
  extern int mode;
  extern double p[4];
  extern double kappa;
  extern double rate;

  double pi_j,a,d;

  
  /*  Case, using the JC model*/
  if(NOTMODE(HKY))
    return (double)(0.25-DELTAF(from,to))*FACTOR*exp(-FACTOR*length);
  /*  Else give the derivative by time in the HKY85 model*/
  else{
    if(to==0 || to==2)
      pi_j=p[0]+p[2];
    else
      pi_j=p[1]+p[3];
    a=1+pi_j*(kappa-1);
  
    
    /*  We have case from and to are the same nucleotide*/
    if(from==to){
      d=-rate*a*(1-p[to]/pi_j)*exp(-rate*a*length);
      d-=rate*p[to]*(1-pi_j)*exp(-rate*length)/pi_j;
    }
    /*  Case where a transition has occurred*/
    else if(((from&1)==0 && (to&1)==0) || ((from&1)!=0 && (to&1)!=0)){
      d=rate*a*p[to]*exp(-rate*a*length)/pi_j;
      d-=rate*p[to]*(1-pi_j)*exp(-rate*length)/pi_j;
    }
    /*  Otherwise a transversion has occurred*/
    else{
      d=p[to]*rate*exp(-rate*length);
    }
    return d;
  }
}


/*  Routine to calculate probability of change along
 * length. Currently using JC model*/
double prob_calcd2(int from, int to, double length){
  extern int mode;
  extern double p[4];
  extern double kappa;
  extern double rate;

  double pi_j,a,d;

  /*  Are we using the JC model?*/
  if(NOTMODE(HKY))
    return (double)(DELTAF(from,to)-0.25)*FACTOR*FACTOR*exp(-FACTOR*length);
  else{
    if(to==0 || to==2)
      pi_j=p[0]+p[2];
    else
      pi_j=p[1]+p[3];
    a=1+pi_j*(kappa-1);
    /*  We have case from and to are the same nucleotide*/
    if(from==to){
      d=rate*rate*p[to]*(1-pi_j)*exp(-rate*length)/pi_j;
      d+=rate*rate*a*a*(1-p[to]/pi_j)*exp(-rate*a*length);
    }
    /*  Case where a transition has occurred A to G, G to A, C to T to T to C*/
    else if(((from&1)==0 && (to&1)==0) || ((from&1)!=0 && (to&1)!=0)){
      d=rate*rate*p[to]*(1-pi_j)*exp(-rate*length)/pi_j;
      d-=rate*rate*a*a*p[to]*exp(-rate*a*length)/pi_j;
    }
    /*  Otherwise a transversion has occurred*/
    else{
      d=-p[to]*rate*rate*exp(-rate*length);
    }
    return d;
  }
}

/*  Routine to calculate the second derivatives of
 * the likelihood*/
double second_derivative(struct treenode *node_p,double *matrix,int b1,int b2){
  extern int branches;
  extern int individual;
  extern int kudge;
  extern int branch1;
  extern int branch2;
  extern int mode;
  double deriv;

  
  if(b1==b2)
    return matrix[b1+branches];
  /*  If we are not sampling, the expectation of the second derivatives is zero
   * and we can just ignore them. On the other hand if we are sampling of doing variance
   * calculations, all the second derivatives are needed*/
  if(NOTMODE(PERCENTILE) && NOTMODE(VARIANCE) && NOTMODE(BOOTSTRAP))
    if(!(ISMODE(PROBS) && ISMODE(INDIVIDUAL) && ISMODE(ROOTED) && individual==1))
      return 0;
  
  branch1=b1;
  branch2=b2;
  kudge=1;
 
  deriv=llhcalc(node_p);

  kudge=0;
  return deriv;
}

/*  Routine to evolve a sequence down the tree and then check if
 * the evolved sequence is in the cache*/
double **pcalc(struct treenode *node_p){
  extern double p[4];
  int a,n;
  double l;
  double **temp;
  
  /*  Generate the original ancestor from the base frequencies*/
  l=randomd();
  a=-1;
  while(l>0){
    l-=p[++a];
  }
  node_p->nucleotide=a;
  
  /*  Recurse down the tree, evolving nucleotides*/
  n=0;
  while(n<DOODAH && CHILD(node)!=NULL){
    pcalc_branch(CHILD(node));
    n++;
  }
  
  /*  Check the result cache*/
  temp=check_cache();
  
  return temp;
}

/*  Routine to evolve a nucleotide down a branch*/
void pcalc_branch(struct treenode *node_p){
  extern int mode;
  extern FILE *prob_file_p;
  int n;
  
  /*  Generate new nucleotide from parent*/
  node_p->nucleotide=evolve_nucleotide(node_p->length[0],(node_p->node[0])->nucleotide);
  if(ISMODE(PROBS) && node_p->node[1]==NULL)
    print_nucleotide(node_p->nucleotide,prob_file_p);
  
  /*  Recurse down tree until leaf*/
  n=1;
  while(n<DOODAH && CHILD(node)!=NULL){
    pcalc_branch(CHILD(node));
    n++;
  }
}
