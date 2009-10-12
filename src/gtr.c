#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "edible.h"
#include "Meschach/matrix.h"
#include "Meschach/matrix2.h"


static double U[16],V[16];
static double E[4];
static double Q[16];
static double GTRscale;
/*  Q = U*E*V, note that V is stored as its transpose  */

int Factorize ( double * A, double * val, int n);

void CreateGTR (void){
  extern double * gtr;
  extern double p[];
  extern int mode;
  double tot;
  int row,col,ent;
 
  assert(NULL!=gtr);
  assert(NULL!=p);
  assert(ISMODE(GTR));

  
  Q[1] = Q[4] = gtr[3];
  Q[2] = Q[8] = 1.;
  Q[3] = Q[12] = gtr[1];
  Q[6] = Q[9] = gtr[4];
  Q[7] = Q[13] = gtr[0];
  Q[11] = Q[14] = gtr[2];

  for ( row=0 ; row<4 ; row++){
    for ( col=0 ; col<4 ; col++){
      Q[row*4+col] *= p[col];
    }
  }

  for ( row=0 ; row<4 ; row++){
    Q[row*4+row] = tot = 0.;
    for ( col=0 ; col<4 ; col++){
       tot -= Q[row*4+col];
    }
    Q[row*4+row] = tot;
  }


  for ( ent=0 ; ent<16 ; ent++){
    U[ent] = Q[ent];
  }

  /*  Make symmetric*/
  for ( row=0 ; row<4 ; row++){
    for ( col=0 ; col<4 ; col++){
      U[row*4+col] *= sqrt(p[row])/sqrt(p[col]);
    }
  }
  Factorize (U,E,4);
  for ( row=0 ; row<4 ; row++){
    if ( E[row]>0.)
      E[row] = 0.;
  }
  /*  U is now eigen vectors of symetricified matrix, unsymmetrify*/
  for ( row=0 ; row<4 ; row++){
    for ( col=0 ; col<4 ; col++){
      V[col*4+row] = U[row*4+col];
    }
  }
  for ( ent=0 ; ent<16 ; ent++){
    U[ent] = V[ent];
  }
  for ( row=0 ; row<4 ; row++){
    for ( col=0 ; col<4 ; col++){
      U[row*4+col] /= sqrt(p[row]);
      V[row*4+col] *= sqrt(p[row]);
    }
  }

  /*  Now have Q = U*E*V^T  */
  GTRscale = 0.;
  for ( row=0 ; row<4 ; row++){
    GTRscale -= Q[row*4+row] * p[row];
  }
  GTRscale = 1. / GTRscale;
}


void GetGTRmat ( double * mat, double length){
  int row,col,k;
  double expev[4];
  double tmp[16];
  double tot;

  assert(NULL!=mat);
  assert(length>=0);

  
  for ( row=0 ; row<4 ; row++){
    expev[row] = exp(length*GTRscale*E[row]);
  }
  for ( row=0 ; row<4 ; row++){
    for ( col=0 ; col<4 ; col++){
      tmp[row*4+col] = U[row*4+col] * expev[col];
    }
  }
  
  for ( row=0 ; row<4 ; row++){
    for ( col=0 ; col<4 ; col++){
      mat[row*4+col] = 0.;
      for ( k=0 ; k<4 ; k++){
	mat[row*4+col] += tmp[row*4+k] * V[col*4+k];
      }
    }
  }

  #ifndef NDEBUG
  for ( row=0 ; row<4 ; row++){
    tot = 0.;
    for ( col=0 ; col<4 ; col++){
      tot+= mat[row*4+col];
      assert(mat[row*4+col]>=0.);
    }
    assert(fabs(tot-1.)<4*DBL_EPSILON);
  }
  #endif
}

void GetQP ( double * mat, double * pmat){
  int row,col,k;
  
  assert(NULL!=mat);
  assert(NULL!=pmat);

  for ( row=0 ; row<4 ; row++){
    for ( col=0 ; col<4 ; col++){
      mat[row*4+col] = 0.;
      for ( k=0 ; k<4 ; k++){
	mat[row*4+col] += Q[row*4+k] * pmat[k*4+col];
      }
      mat[row*4+col] *= GTRscale;
    }
  }
}
	  

int Factorize ( double * A, double * val, int n){
  MAT * B, * Q;
  VEC * out;
  int i,j;

  assert(NULL!=A);
  assert(NULL!=val);
  assert(n>0);

  out = malloc(sizeof(VEC));
  out->max_dim = out->dim = n;
  out->ve = val;

  B = malloc(sizeof(MAT));
  B->m = B->n = B->max_m = B->max_n = n;
  B->max_size = n*n;
  B->base = A;
  B->me = malloc(n*sizeof(double *));
  for ( i=0 ; i<n ; i++){
    B->me[i] = &A[i*n];
  }

  Q = m_get (n,n);

  symmeig (B,Q,out);
  for ( i=0 ; i<n ; i++){
    for ( j=0 ; j<n ; j++){
      A[j*n+i] = Q->me[i][j];
    }
  }

  m_free(Q);
  free(B->me);
  free(B);
  free(out);

  return 1;
}

