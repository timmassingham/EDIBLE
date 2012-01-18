#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "edible.h"
#include <string.h>
#include <time.h>
#include "Meschach/matrix.h"
#include "Meschach/matrix2.h"


/* Calculate x^y, returning an integer.
 * No overflow checking*/
int ipow(int x, int y){
  int a,b;
  b=1;
  for(a=0;a<y;a++)
    b*=x;
  return b;
}

/* If memory can't be allocated - crash*/
void nomemory(void){
  printf("Out of memory\n*Plunk*\n");
  exit(2);
}



/*  Convert +ve number to string
 * Kernighan and Ritchie can take the blame for this if
 * it doesn't work properly*/
char *itotext(int n,char *s){
  int c,i,j;

  i=0;

  do{
    s[i++]=n%10+'0';
  } while((n/=10)>0);
  s[i]='\0';

  for(i=0,j=strlen(s)-1;i<j;i++,j--){
    c=s[i];
    s[i]=s[j];
    s[j]=c;
  }
  return s;
}


/*  Routine to take the matrix given and calculate the log-
 * determinant, by calling LU decomposition routine and
 * then multiplying down diagonals. Returns the
 * log-determinant calculated.*/
double * determinant(void){
  int a,max,c;
  extern int branches;
  double *det;
  extern int mode;
  extern int nodecount;
  extern double (*(*expect)[])[];
  extern double (*(*rootedexpect)[])[];
  extern int individual;
  extern int interesting_branches[];
  extern int is_kappa;
  double (*(*matrix)[])[];
  MAT * matrix2;

  is_kappa=0;
  if(ISMODE(HKY) && NOTMODE(NOKAPPA))
    is_kappa=1;
  matrix=expect;
  max=branches;
  if(ISMODE(ROOTED)){ /*  If want rooted tree then create new*/
    planttree(expect,rootedexpect);   /* matrix*/
    matrix=rootedexpect;
    max=nodecount+2;
    if(ISMODE(NODEASROOT))
      max=nodecount+1;
  }

  if(ISMODE(MATRICES)){  /* If want intermediate matrices dumped*/
    dump(matrix,max+is_kappa,"Full matrix");
  }

  if(ISMODE(INDIVIDUAL)){ /*  We want information about some, but
                           * not all of the elements*/
    if(NOTMODE(DETINDIV)){
      det=calloc(individual+is_kappa,sizeof(double));
      for(a=0;a<individual;a++)
        det[a]=(*(*matrix)[interesting_branches[a]])[interesting_branches[a]];
      if(is_kappa==1)
	det[individual]=(*(*matrix)[max])[max];
      is_kappa=0;
      return det;
    }

    /*  Case - we want the determinate of the sub-matrix formed 
     * by several parameters*/
    /*  Get memory for new matrix*/
    matrix2 = m_get(individual+is_kappa,individual+is_kappa);
    if(NULL==matrix2){
	    nomemory();
    }
    m_zero(matrix2);


    /*  Creates the sub-matrix from the original expected information
     * matrix*/
    for(a=0;a<individual;a++)
      for(c=0;c<individual;c++)
	matrix2->me[a][c]=(*(*matrix)[interesting_branches[a]])[interesting_branches[c]];
    if(is_kappa==1){
      matrix2->me[individual][individual]=(*(*matrix)[max])[max];
    }
    
    max=individual;
    if(ISMODE(MATRICES))
      dump(matrix2->me,max,"Sub-matrix to be calculated");
  } else {
      matrix2 = m_get(max,max);
      if(NULL==matrix2){
          nomemory();
      }
      m_zero(matrix2);
      for ( a=0 ; a<max ; a++){
          for ( c=0 ; c<max ; c++){
              matrix2->me[a][c] = (*(*matrix)[a])[c];
          }
      }
  }
 
  /*  Perform LU decomposition on whichever matrix we've been handed*/
  det=calloc(1+is_kappa,sizeof(double));
  matrix2=CHfactor(matrix2);

  /*  The determinant of the matrix is the product of
   * the diagonal elements of the decomposed form*/
  for(a=0;a<max;a++){
    det[0] += 2.0 * log(matrix2->me[a][a]);
  }
  if(is_kappa==1){
    det[1] = 2.0 * log(matrix2->me[max][max]);
  }

  M_FREE(matrix2);

  return det;
}


/*  Routine to take given matrix and perform LU
 * factorisation with row pivoting on it.
 *  Returns the number of pivots used in calculation.
 * Will fail for extremely degenerate matrices (all
 * possible choices of pivot down rows are zero.) In which
 * case, expect divide by zero error.
 *  Assumes matrix is symmetric - is the case when we will 
 * be using the routine*/
int ludecomp(double (*(*matrix)[])[],int max){
  int level;
  int pivots;
  int a,b;
  double (*row)[];
  double elt;

  level=0;
  pivots=0;

  /*  Go through all the rows from the top to bottom*/
  while(level<(max-1)){
    b=level;
    /*  In general we are best using the entry in the
     * row with the greatest absolute value as a pivot.
     * Doing this switch multiplies the determinant by -1*/
    for(a=level;a<max;a++)
      b=(fabs((*(*matrix)[a])[level])>fabs((*(*matrix)[b])[level]))?a:b;
    if(b!=level){
      row=(*matrix)[b];
      (*matrix)[b]=(*matrix)[level];
      (*matrix)[level]=row;
      pivots++;
    }

    /*  Pivot element*/
    elt=(*(*matrix)[level])[level];
    /*  If pivot is very small, then all the elements must be
     * almost zero and the method is unreliable at best - die
     * with determinant 0*/
    if(fabs(elt)<1e-50)
      return -1;

    /*  Scale all the row by the pivot element*/
    for(a=(level+1);a<max;a++)
      (*(*matrix)[a])[level]/=elt;
    /*  Subtract new row from all remaining rows*/
    for(a=(level+1);a<max;a++)
      for(b=(level+1);b<=a;b++){  /*  Use the symmetry of the matrix*/
        (*(*matrix)[a])[b]-=(*(*matrix)[a])[level]*(*(*matrix)[level])[b];
        (*(*matrix)[b])[a]=(*(*matrix)[a])[b];
      }
    level++;
  }
  return pivots;
}

/*  Modification to getc to ignore white space*/
char getnextc(FILE *fp){
  char a;
  int b;

  b=0;
  while(b==0){
    a=getc(fp);
    if(a!=' ' && a!='\n' && a!='\t')
      b++;
  }
  return a;
}

/*  Fairly simple function to print output nucleotides*/
void print_nucleotide(int nucleo, FILE *fp){

  switch(nucleo){
  case 0:fprintf(fp,"A");
         break;
  case 1:fprintf(fp,"C");
         break;
  case 2:fprintf(fp,"G");
         break;
  case 3:fprintf(fp,"T");
         break;
  }
}

/*  Initialise the random number generator*/
void initialise_rg(void){
  extern int seed;
  time_t *me;

  me=NULL;                                                             
  /*  Use time to generate seed                                                 
   * (time = seconds since 00:00.00 1 Jan '70) On UNIX systems...*/             
           seed=(time(me)); /*  Initialise random number generator.*/           
           seed=(seed==0)?1:seed; /* 0 can never be a seed or the*/
}


/*  Put sequence from array onto a tree and print it out if necessary*/
void tree_sequence(unsigned int a){
  extern int leaves;
  extern struct treenode *leaf[];
  extern int mode;
  extern FILE *prob_file_p;

  int b,d;

  for(b=leaves-1;b>-1;b--){
    d=a&(3<<(2*b));
    d>>=(2*b);
    leaf[leaves-1-b]->nucleotide=d;
    if(ISMODE(PROBS))
      print_nucleotide(leaf[leaves-1-b]->nucleotide,prob_file_p);
  }
}

/*  Do the information calculation, so the code isn't cluttered with
 * the formula in 8 or so places*/
double evaluate_information(double (*(*matrix)[])[],const int a,const int b){
  double v;
  extern int mode;
  extern int branches;

  if ( (*(*matrix)[branches])[branches]==0.)
    return 0.;
  v=(*(*matrix)[branches])[a]*(*(*matrix)[branches])[b]/(*(*matrix)[branches])[branches];
  if(ISMODE(PERCENTILE) || ISMODE(VARIANCE) || ISMODE(BOOTSTRAP))
    v=(v-(*(*matrix)[a])[b]);
  if(ISMODE(BOOTSTRAP) || ISMODE(PERCENTILE))
    v/=(*(*matrix)[branches])[branches];

  return v;
}

/*  Calculate the new rate if Kappa changes*/
void do_rate(void){
  extern double rate;
  extern double kappa;
  extern double p[4];

  rate=2*kappa*(p[0]*p[2]+p[1]*p[3])
                   +2*(p[0]*p[1]+p[0]*p[3]+p[1]*p[2]+p[2]*p[3]);
  rate=1/rate;
}

/*  Lovely O(n^2) algorithm for sorting an array, just like
 * young computer scientists get slapped round the wrist for creating.
 * Seriously, it isn't going to be sorting any great amount of data,
 * so having a stupid algorithm doesn't matter*/
void thick_sort(int array[],int length){
  int a,b,sh=1;
  int max;
  
  while(sh==1){
    sh=0;
    b=0;
    max=array[0];
    
    for(a=1;a<length;a++)
      if(array[a]>max){
	max=array[a];
	sh=1;
	b=a;
      }
    
    length--;
    if(array[length]<max){
      array[b]=array[length];
      array[length]=max;
    }
  }
}

