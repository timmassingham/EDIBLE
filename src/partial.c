#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include "edible.h"
#include "new_models.h"
#include <float.h>


/*  Starts off calculation of likelihood of given
 * tree.*/
void partialcalc(struct treenode *node_p, double (*matrix)[], int twig){
  extern int branches;
  int a,b,n,r;
  double l,llh, llhr;
  double tot;
  
  extern double p[4];
  extern int nrates;
  extern double * rate_cat, * rate_prob;
  
  assert(NULL!=node_p);
  assert(NULL!=rate_cat);
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
      partialcalc_branch(CHILD(node),twig, 1, r);
      n++;
    }

    a=n; llhr = 0;
    for(b=0;b<4;b++){/*  Loop to calculate the partial likelihood of this tree */
      l=1;           /* given the partial likelihoods already calculated.      */
      for(n=0;n<a;n++)
        l*=CHILD(node)->llh[b];
      llhr += p[b]*l;
    }

    llh += rate_prob[r] * rate_cat[r] * llhr;    
  }
  (*matrix)[twig] = llh;
  
  llh = 0.;
  for ( r=0 ; r<nrates ; r++){
    /*  Function to recurse down the tree*/
    n=0;
    while(n<DOODAH && CHILD(node)!=NULL){
      partialcalc_branch(CHILD(node),twig, 2, r);
      n++;
    }

    a=n; llhr = 0;
    for(b=0;b<4;b++){/*  Loop to calculate the partial likelihood of this tree */
      l=1;           /* given the partial likelihoods already calculated.      */
      for(n=0;n<a;n++)
        l*=CHILD(node)->llh[b];
      llhr += p[b]*l;
    }

    llh += rate_prob[r] * rate_cat[r] * rate_cat[r] * llhr;    
  }
  (*matrix)[twig+branches] = llh;  
}

void partialcalc_branch(struct treenode *node_p, int twig, int deriv, int r){
  extern struct treenode *branch[];
  
  int a,n;
  
  assert(deriv==1 || deriv==2);
  
  /*  If we are on a leaf, do initial calculation of prob.*/
  n=1;
  if(CHILD(node)==NULL){
    for(a=0;a<4;a++){
      if(branch[twig]==node_p){
        if ( deriv==1){
	  node_p->llh[a]=prob_calcd_new(a,node_p->nucleotide,node_p->bnum,r);
        } else if (deriv==2)
          node_p->llh[a]=prob_calcd2_new(a,node_p->nucleotide,node_p->bnum,r);
      } else {
	node_p->llh[a]=prob_calc_new(a,node_p->nucleotide,node_p->bnum,r);
      }
    }
    return;
  }

/*  Now - can assume that we are on an internal node
 * Recurse down all the (non-parent) branches        */
  while(n<DOODAH && CHILD(node)!=NULL){
    partialcalc_branch(CHILD(node),twig,deriv,r);
    n++;
  }
/*  We can assume the llh's of all the child nodes
 * have been calculated - and thus can calculate the
 * llh vector for the current node*/
  a=n;
/*  'a' is now either DOODAH, or the first 'branch'
 * with nothing coming off it                        */
  if(branch[twig]==node_p){
    if ( deriv==1)
      llhvector(node_p,a,1,r);
    else if (deriv==2)
      llhvector(node_p,a,2,r);
  } else {
    llhvector(node_p,a,0,r);
  }
}


/*  Routine to calculate the conversion matrix between
 * the rooted and unrooted forms of the tree.*/
int (*(*findspade(struct treenode *node_p))[])[]{
  extern int branches;
  extern int nodecount;
  extern int root;
  extern int mode;
  extern struct treenode *branch[];
  int (*(*conv_matrix)[])[];
  int a,c,end,is_k=0;
  struct treenode *node_c;

  c=(ISMODE(NODEASROOT))?1:2;
  
  /*  Have extra row and column is we are using HKY85 model*/
  if(ISMODE(HKY) && NOTMODE(NOKAPPA)){
    c++;
    is_k=1;
  }

  /*  Get memory from system for the conversion matrix*/
  conv_matrix=calloc((nodecount+c),sizeof(int *));
  if(conv_matrix==NULL)
    nomemory();
  for(a=0;a<(nodecount+c);a++){
    (*conv_matrix)[a]=calloc(branches+is_k,sizeof(int));
    if((*conv_matrix)[a]==NULL)
      nomemory();
  }

/*  Using the made-up node names to differentiate the
 * the parameters in the rooted model.
 * Doing basic set up for the matrix*/
  for(a=0;a<branches;a++){
    sscanf((branch[a]->node[0])->name,"Node-%d",&end);
    (*(*conv_matrix)[end])[a]=1;
    end=-1;
    sscanf(branch[a]->name,"Node-%d",&end);
    if(end!=-1)
      (*(*conv_matrix)[end])[a]=-1;
  }
  node_c=branch[root];

  /*  Since the (code) root and the (genetic) root may not
   * necessarily coincide, we must recurse up the tree 
   * from the genetic root correcting the matrix as we go*/
  while(node_c->node[0]!=node_p){
    sscanf((node_c->node[0])->name,"Node-%d",&end);
    a=0;
    while(branch[a]!=node_c){a++;}
    (*(*conv_matrix)[end])[a]=-1;
    a=0;
    while(branch[a]!=node_c->node[0]){a++;}
    (*(*conv_matrix)[end])[a]=1;
    node_c=node_c->node[0];
  }


  /* Deal with end case if we have "floating" genetic root*/
  if(NOTMODE(NODEASROOT)){
  /*  Root node has different behaviour from the others*/
    (*(*conv_matrix)[nodecount+1])[root]=2;
    a=0;
    while(branch[a]!=node_c){a++;}
    (*(*conv_matrix)[0])[a]=-1;
  }
  else{
    sscanf((branch[root]->node[0])->name,"Node-%d",&end);   
    (*(*conv_matrix)[end])[root]=1;
  }
  
  /*  If we are using the HKY85 model, the information relating to
   * kappa doesn't change in the rooted case*/
  if(ISMODE(HKY) && NOTMODE(NOKAPPA))
    (*(*conv_matrix)[nodecount+c-1])[branches]=1;
  
  return conv_matrix;
}
