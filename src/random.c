#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "edible.h"
#include <float.h>

/*  Routine to sample from data set and return an estimate
 * of the amount of information that can be expected in
 * (100-x)% of experiments*/
double sample_percentile(struct treenode *node_p, struct treenode *tree2,unsigned int e,int factor_flag,double factor){
  int a;
  extern FILE *sample_file_p;
  extern int branches;
  extern int leaves;
  extern int sample_size;
  extern int sequence_length;
  extern int mode;
  int percent_point;
  double information;
  extern double percentile;
  extern int is_kappa;
  double (*percent)[];
  double (*(*partial)[])[];

  /*  Make copy of our tree (so we can reuse some of the already
   * calculate values*/
  filltree(node_p,tree2,0);
  leaves=0;
  doleaf(tree2,0);/*  Change the leaf pointers so they point to the
                  * new tree*/

  /*  If we are using the HKY85 model, then we need a bit
   * more memory to calculate the derivative WRT kappa*/
  if(ISMODE(HKY) && NOTMODE(NOKAPPA))
    is_kappa=1;
  
  /*  Get memory for array of first and second derivatives*/
  partial=calloc(branches+1+is_kappa,sizeof(double *));
  if(partial==NULL)
    nomemory();
  for(a=0;a<(branches+1+is_kappa);a++){
    (*partial)[a]=calloc(branches+1+is_kappa,sizeof(double));
    if((*partial)[a]==NULL)
      nomemory();
  }
    
  /*  Want the median to be between the last two entries of
   * the array.*/
  percent_point=(int)(percentile*sample_size+1.5);
  percent=calloc(sample_size,sizeof(double));   
  if(percent==NULL)
    nomemory();
  for(a=0;a<percent_point;a++)
    (*percent)[a]=DBL_MAX;

  /*  Sample from distribution, keeping only the bottom so many 
   * results (depending on what percentile point we are trying 
   * to calculate*/ 
  for(a=0;a<sample_size;a++){
    information=randomresult(node_p,tree2,partial,e,sequence_length,factor_flag,factor);
    reorder(percent,information,percent_point);
  }

/*  Calculate the weighted median of the data set -
 * interpolates between the last two values to get a
 * "more accurate" estimate.*/
  information=(*percent)[percent_point-2]
             +(percentile*sample_size+1.5-percent_point)
             *((*percent)[percent_point-1]-(*percent)[percent_point-2]);

  for(a=0;a<(branches+1+is_kappa);a++)
    free((*partial)[a]);
  free(partial);
  if(sample_file_p!=NULL){
    for(a=0;a<sample_size;a++)
      fprintf(sample_file_p,"%e\n",(*percent)[a]);
  }
  free(percent);

  return information;
}


/*  Routine to return a result from an experiment generated
 * randomly.*/
double randomresult(struct treenode *node_p, struct treenode *tree,double (*(*partial)[])[],unsigned int e,int sequence,int factor_flag,double factor){
  int a,c,d;
  double (*data)[];
  double b;
  extern int branches;
  extern double (*(*expect)[])[];
  extern int mode;
  
  for(a=0;a<branches;a++)
    for(c=0;c<(a+1);c++) 
      (*(*expect)[a])[c]=0;
 
  for(d=0;d<sequence;d++){
    /*   Get a new sample from the distribution*/
    calcllhs(node_p,tree,e,partial);
    
    for(a=0;a<branches;a++)  /* Add new entry to expectation matrix*/
      for(c=0;c<(a+1);c++){
        (*(*expect)[a])[c]+=evaluate_information(partial,a,c);
        (*(*expect)[c])[a]=(*(*expect)[a])[c];
      }
    if(ISMODE(HKY) && NOTMODE(NOKAPPA))
      (*(*expect)[branches])[branches]+=evaluate_information(partial,branches+1,branches+1);
  }

  /*  The estimate of the expected information is the sum of the
   * observed informations/number of results*/
  for(a=0;a<branches;a++)
    for(c=0;c<branches;c++)
      (*(*expect)[a])[c]/=sequence;
  if(ISMODE(HKY) && NOTMODE(NOKAPPA))
    (*(*expect)[branches])[branches]/=sequence;

  /*  Not having to worry about variance calculation (ATM;-)
   * scaling of information about kappa is more complicated and so
   * happens separately */
  scale_tree(factor_flag,factor,expect,branches);
  
  data=determinant();

  if(ISMODE(HKY) && NOTMODE(NOKAPPA))
     b=(*data)[1];
  else
    b=(*data)[0];
  
  free(data);
  return b;
}


/*  Routine to return random double between 0.0 and 1.0
 * fine enough for data*/
#define IM 2147483647

double randomd(void){
  long a,b;
  unsigned long c;
  double d;   
 
  a=randoms();
  b=randoms();
  c=IM;
  c++;
  d=(double)a/IM;
  d+=b;
  d/=c;

  return d;   
}

/*  Routine to return random number based on minimum
 * standard generator.
 * Seed _cannot_ be zero.*/
#define IA 16807  
#define IQ 127773
#define IR 2836

long randoms(void){
  long k;
  extern int seed;
  
  k=seed/IQ; 
  seed=IA*(seed-k*IQ)-IR*k;
  if(seed<0)
    seed+=IM;
  return seed; 
}


/*  Routine to place element in its correct place in an array
 * knocking off elements off the end if necessary*/
void reorder(double (*percent)[],double information,int percent_point){
  extern int sample_size;
  int a,b;
  double d,e;

  a=0;
  while(a<sample_size && (*percent)[a]<information){a++;}

  /*  Else move all the elements of the array along one
   * and dumps the one which falls off the end*/
  d=(*percent)[a];
  (*percent)[a]=information;
  for(b=++a;b<sample_size;b++){
    e=(*percent)[b];
    (*percent)[b]=d;
    d=e;
  }
}

/*  Routine to evolve a new nucleotide from old, using the model
 * given.*/
int evolve_nucleotide(double length, int start_nucleotide){
  double p;
  double l;
  int nucleo=0;

  l=randomd();

  do{
    p=prob_calc(start_nucleotide, nucleo, length);
    if(l<p)
      return nucleo;

    l-=p;
  }while(++nucleo<4);

  printf("Nucleotide evolves into something closely resembling an error\n"
	 "Measure not one - probability theory failure. Terminating\n");
  exit(0);
}

/*  Check to see if the current tree is in the cache and return 
 * a pointer to it matrix if it is. NULL if not*/
double (*(*check_cache(void))[])[]{
  extern struct treenode *leaf[];
  extern int leaves;
  extern struct crecord (*cache)[];
  extern int cache_size;
  extern int mode;
  
  int a,b,match=0;
  double (*(*temp)[])[];
  
  temp=NULL;
  if(NOTMODE(CACHE))
    return temp;
  
  for(a=0;a<cache_size && match==0 && (*(*cache)[a].leaf_nucleotide)[0]!=-1;a++){
    match=1;
    for(b=0;b<leaves && match==1;b++)
      if((*(*cache)[a].leaf_nucleotide)[b]!=leaf[b]->nucleotide)
	match=0;
    if(match==1)
      temp=(*cache)[a].matrix;
  }
  
  return temp;
}

/*  See if we can fit the result generated into cache*/
void update_cache(double (*(*matrix)[])[]){
  extern struct crecord (*cache)[];
  extern int cache_size;
  extern struct treenode *leaf[];
  extern int mode;
  extern int is_kappa;
  extern int branches;
  extern int leaves;
  
  int a,b,i,j;
  
  if(ISMODE(HKY) && NOTMODE(NOKAPPA))
    is_kappa=1;
  
  b=0;
  for(a=0;a<cache_size && (*(*cache)[a].leaf_nucleotide)[0]!=-1;a++)
    if((*(*(*cache)[a].matrix)[branches])[branches]
       <(*(*(*cache)[b].matrix)[branches])[branches])
      b=a;

  if(a!=cache_size)
    b=a;
  else
    if((*(*(*cache)[b].matrix)[branches])[branches]>(*(*matrix)[branches])[branches])
      return;
  
  for(i=0;i<leaves;i++)
    (*(*cache)[b].leaf_nucleotide)[i]=leaf[i]->nucleotide;
  for(i=0;i<branches+is_kappa+1;i++)
    for(j=0;j<branches+is_kappa+1;j++)
      (*(*(*cache)[b].matrix)[i])[j]=(*(*matrix)[i])[j];
}

/*  Clean the cache for another round of results to be generated*/
void wipe_cache(void){
  extern struct crecord (*cache)[];
  extern int cache_size;
  int a;
  
  for(a=0;a<cache_size;a++)
    (*(*cache)[a].leaf_nucleotide)[0]=-1;
}
