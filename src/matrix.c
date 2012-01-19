#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "edible.h"

/*  Routine to calculate the matrix of expected information
 * from the first derivatives. Returns pointer to matrix*/
void expectation(unsigned int e,double factor,int factor_flag, struct treenode *tree, struct treenode *tree2){

  extern double **expect;
  extern double **var;
  extern double **var2;
  int a,b,d=0;
  unsigned int c;
  extern int branches;
  extern int nodecount;
  extern int mode;
  extern int individual;
  extern int is_kappa;
  extern int interesting_branches[];
  extern int nodecount;
  
  double **matrix;
  double **matrix2;
  double temp;

  /*  If we are using HKY85, we need a bit more memory for the
   * derivative WRT kappa*/
  if(ISMODE(HKY) && NOTMODE(NOKAPPA))
    is_kappa=1;
  
  /*  Zero the relevant arrays.*/
  for(a=0;a<branches+is_kappa;a++)
    for(b=0;b<branches+is_kappa;b++){
      expect[a][b]=0;
    }

  if(ISMODE(VARIANCE)){
    d=branches+is_kappa;
    if(ISMODE(ROOTED))
      d=nodecount+((ISMODE(NODEASROOT))?1:2)+is_kappa;
    
    for(a=0;a<d;a++)
      for(b=0;b<d;b++)
        var[a][b]=0;
  }

  
  /*  Get array to do all our intermediate calculations in*/
  matrix=calloc(branches+1+is_kappa,sizeof(double *));
  if(matrix==NULL)
    nomemory();
  for(a=0;a<(branches+1+is_kappa);a++){
    matrix[a]=calloc(branches+1+is_kappa,sizeof(double));
    if(matrix[a]==NULL)
      nomemory();
  }
  matrix2=calloc(branches+1+is_kappa,sizeof(double *));
  if(matrix2==NULL)
    nomemory();
  for(a=0;a<(branches+1+is_kappa);a++){
    matrix2[a]=calloc(branches+1+is_kappa,sizeof(double));
    if(matrix2[a]==NULL)
      nomemory();
  }

  is_kappa=0;
  
  /*  If we want information about one parameter only and
   * are not in a rooted tree, only calculate that
   * information*/
  if(ISMODE(INDIVIDUAL) && individual==1 && NOTMODE(ROOTED)){
    for(c=0;c<e;c++){
      calcllhs(tree,tree2,c,matrix);

      /*  If we are sampling then we have to take into account the
       * partial second derivatives.
       *  Also must divide by the probability since we are sampling rather than
       * calculating the expectation*/
      temp=evaluate_information(matrix,interesting_branches[0],interesting_branches[0]);
      expect[interesting_branches[0]][interesting_branches[0]]+=temp;
      if(ISMODE(VARIANCE)){
        temp*=temp/matrix[branches][branches];
        if(ISMODE(BOOTSTRAP))
          temp*=matrix[branches][branches];
        var[interesting_branches[0]][interesting_branches[0]]+=temp;
      }

      /*  If we are using HKY85 model and wish information about kappa*/
      if(ISMODE(HKY) && NOTMODE(NOKAPPA)){
        temp=evaluate_information(matrix,branches+1,branches+1);
        expect[branches][branches]+=temp;
        if(ISMODE(VARIANCE)){
          temp*=temp/matrix[branches][branches];
          if(ISMODE(BOOTSTRAP))
            temp*=matrix[branches][branches];
          var[branches][branches]+=temp;
	}
      }
    }

      
    /*  If we are sampling for the expected information, then we
     * must divide by the number of samples*/
    if(ISMODE(BOOTSTRAP)){
      expect[interesting_branches[0]][interesting_branches[0]]/=e;
      if(ISMODE(VARIANCE))
	var[interesting_branches[0]][interesting_branches[0]]/=e;
      /*  Do same if we are sampling for Kappa*/
      if(ISMODE(HKY) && NOTMODE(NOKAPPA)){
	expect[branches][branches]/=e;
	if(ISMODE(VARIANCE))
	  var[branches][branches]/=e;
      }
    }
  }
  else{
  /*  Need to calculate entire expectation matrix*/
    for(c=0;c<e;c++){
      calcllhs(tree,tree2,c,matrix);

      for(a=0;a<branches;a++)
        for(b=0;b<(a+1);b++){  /*  Expectation matrix is symmetric  */
          /*  If we are sampling then we have to take into account the 
	   * second order partial derivatives.
	   *  Also must divide by the probability since we are sampling
           * rather than calculating the expectation*/
	  matrix2[a][b]=evaluate_information(matrix,a,b);
          expect[a][b]+=matrix2[a][b];
          matrix2[b][a]=matrix2[a][b];
        }
      if(ISMODE(HKY) && NOTMODE(NOKAPPA)){
	matrix2[branches][branches]=evaluate_information(matrix,branches+1,branches+1);
	expect[branches][branches]+=matrix2[branches][branches];
	is_kappa=1;
      }

      /*  If we want the variance of the expected information to be
       * calculated*/
      if(ISMODE(VARIANCE)){
        /*  Rooted case, we must convert the information just calculated
         * the rooted form.*/
        if(ISMODE(ROOTED)){
          planttree(matrix2,var2);
          for(a=0;a<nodecount+((ISMODE(NODEASROOT))?1:2)+is_kappa;a++)
            for(b=0;b<nodecount+((ISMODE(NODEASROOT))?1:2)+is_kappa;b++){
              temp=var2[a][b]*var2[a][b]/matrix[branches][branches];
              if(ISMODE(BOOTSTRAP))
                temp*=matrix[branches][branches];
              var[a][b]+=temp;
            }
        }
	/*  Or the far easier unrooted form - Look, no matrix
	 * multiplication!*/
        else{
          for(a=0;a<branches+is_kappa;a++)
            for(b=0;b<branches+is_kappa;b++){
              temp=matrix2[a][b]*matrix2[a][b]/matrix[branches][branches];
              if(ISMODE(BOOTSTRAP))
                temp*=matrix[branches][branches];
              var[a][b]+=temp;
            }
        }
      }
      is_kappa=0;
    }
    
    if(ISMODE(HKY) && NOTMODE(NOKAPPA))
      is_kappa=1;
    /*  If we are sampling for the expected information, then we must
     * divide by the number of samples*/
    if(ISMODE(BOOTSTRAP)){
      for(a=0;a<branches+is_kappa;a++)
	for(b=0;b<(a+1);b++)
	  expect[a][b]/=e;
      if(ISMODE(VARIANCE))
	for(a=0;a<d;a++)
	  for(b=0;b<d;b++)
	    var[a][b]/=e;
    }
  }

  /*  Exploit symmetry*/
  for(a=0;a<branches+is_kappa;a++)
    for(b=0;b<(a+1);b++){
      expect[b][a]=expect[a][b];
    }

  /*  Scale the tree.*/
  scale_tree(factor_flag,factor,expect,branches);
  /*  We also need to scaling the variance matrix,
   * the scaling has twice the effect on the variance matrix
   * since it is, in effect, an expectation squared*/
  if(ISMODE(VARIANCE)){
    a=branches;
    if(ISMODE(ROOTED))
      a=nodecount+((ISMODE(NODEASROOT))?1:2);
    scale_tree(factor_flag,factor,var,a);
    scale_tree(factor_flag,factor,var,a);
  }
  
  for(a=0;a<(branches+1+is_kappa);a++){
    free(matrix[a]);
    free(matrix2[a]);
  }
  free(matrix);
  free(matrix2);
  
  is_kappa=0;
}



/*  Routine to calculate the matrix of expectation of the
 * rooted tree from the unrooted*/
void planttree(double **original_matrix,double **new_matrix){
  double **matrix;
  extern int **conv_matrix;
  int a,b,c,d,is_k=0;
  extern int branches;
  extern int nodecount;
  extern int mode;

  d=(ISMODE(NODEASROOT))?1:2;
  
  /*  Have an extra row/column if we want information about kappa*/
  if(ISMODE(HKY) && NOTMODE(NOKAPPA)){
    d++;
    is_k=1;
  }

/*  Get memory for necessary matrix*/
  matrix=calloc(nodecount+d,sizeof(double *));
  if(matrix==NULL)
    nomemory();
  for(a=0;a<nodecount+d;a++){
    matrix[a]=calloc(branches+is_k,sizeof(double));
    if(matrix[a]==NULL)
      nomemory();
  }

/*  Zero rootedexpect matrix*/
  for(a=0;a<nodecount+d;a++)
    for(b=0;b<nodecount+d;b++)
      new_matrix[a][b]=0;

/*  Do first matrix multiplication, conv_matrix*expectation*/
  for(a=0;a<nodecount+d;a++)
    for(b=0;b<branches+is_k;b++)
      for(c=0;c<branches+is_k;c++)
        matrix[a][b]+=(double)conv_matrix[a][c]*original_matrix[c][b];

/*  Do second matrix multiplication, result*conv_matrix^{t}*/
  for(a=0;a<nodecount+d;a++)
    for(b=0;b<nodecount+d;b++)
      for(c=0;c<branches+is_k;c++)
       new_matrix[a][b]+=matrix[a][c]*conv_matrix[b][c];

/*  Free memory used, since we don't need the intermediate
 * matrix any more*/
  for(a=0;a<nodecount+d;a++)
    free(matrix[a]);
  free(matrix);
}

void rooted_derivative(double *original_vector,double *new_vector){
  extern int **conv_matrix;
  int a,b,d;
  extern int branches;
  extern int nodecount;
  extern int mode;

  d=(ISMODE(NODEASROOT))?1:2;

  /*  Not worried about kappa in this case*/

  /*  Zero the new vector*/
  for(a=0;a<nodecount+d;a++)
    new_vector[a]=0;

  /*  Multiply the 1st partial derivative vector by the conversion matrix
   * to get the new  (rooted) form*/
  for(a=0;a<nodecount+d;a++)
    for(b=0;b<branches;b++)
      new_vector[a]+=conv_matrix[a][b]*original_vector[b];
}


