#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "edible.h"
#include "gtr.h"


#define OOM(A) { if (NULL== (A) ) {puts("Out of memory!"); exit (EXIT_FAILURE);} }

static double *** ProbMatArray=NULL,***DProbMatArray=NULL, *** DDProbMatArray=NULL;



double prob_calc_new ( int from, int to, int b, int r){
  extern int branches;
  extern int nrates;
  double result;
  int bran,rc;

  assert(from>=0 && from<4);
  assert(to>=0 && to<4);
  assert(b>=0 && b<branches);
  assert(r>=0 && r<nrates);
  assert(NULL!=ProbMatArray);
  #ifndef NDEBUG
    for ( bran=0 ; bran<branches ; bran++){
      assert(NULL!=ProbMatArray[bran]);
      for ( rc=0 ; rc<nrates ; rc++){
	assert(NULL!=ProbMatArray[bran][rc]);
      }
    }
  #endif

  result = ProbMatArray[b][r][from*4+to];


  assert(result>=0 && result<=1);
  return result;
}


double prob_calcd_new ( int from, int to, int b, int r){
  extern int branches;
  extern int nrates;
  double result;
  int bran,rc;

  assert(from>=0 && from<4);
  assert(to>=0 && to<4);
  assert(b>=0 && b<branches);
  assert(r>=0 && r<nrates);
  assert(NULL!=DProbMatArray);
  #ifndef NDEBUG
    for ( bran=0 ; bran<branches ; bran++){
      assert(NULL!=DProbMatArray[bran]);
      for ( rc=0 ; rc<nrates ; rc++){
        assert(NULL!=DProbMatArray[bran][rc]);
      }
    }
  #endif


  result = DProbMatArray[b][r][from*4+to];


  return result;
}

double prob_calcd2_new ( int from, int to, int b, int r){
  extern int branches;
  extern int nrates;
  double result;
  int bran,rc;

  assert(from>=0 && from<4);
  assert(to>=0 && to<4);
  assert(b>=0 && b<branches);
  assert(r>=0 && r<nrates);
  assert(NULL!=DDProbMatArray);
  #ifndef NDEBUG
    for ( bran=0 ; bran<branches ; bran++){
      assert(NULL!=DDProbMatArray[bran]);
      for ( rc=0 ; rc<nrates ; rc++){
        assert(NULL!=DDProbMatArray[bran][rc]);
      }
    }
  #endif


  result = DDProbMatArray[b][r][from*4+to];


  return result;
}


void CreatePMats ( void ){
  extern int branches;
  extern struct treenode *branch[];
  extern int nrates;
  extern double *  rate_cat;
  extern int mode;
  int bran,r,from,to;
  double length;

  assert(branches>0);
  assert(nrates>0);
  assert(NULL!=rate_cat);
  assert( (NULL==ProbMatArray && NULL==DProbMatArray && NULL==DDProbMatArray)
        ||(NULL!=ProbMatArray && NULL!=DProbMatArray && NULL!=DDProbMatArray) );

  if ( NULL==ProbMatArray ){
    ProbMatArray = calloc (branches,sizeof(double **));		OOM(ProbMatArray);
    DProbMatArray = calloc (branches,sizeof(double **));        OOM(DProbMatArray);
    DDProbMatArray = calloc (branches,sizeof(double **));       OOM(DDProbMatArray);
    for ( bran=0 ; bran<branches ; bran++){
      ProbMatArray[bran] = calloc(nrates,sizeof(double *));	OOM(ProbMatArray[bran]);
      DProbMatArray[bran] = calloc(nrates,sizeof(double *));    OOM(DProbMatArray[bran]
);
      DDProbMatArray[bran] = calloc(nrates,sizeof(double *));   OOM(DDProbMatArray[bran]);
      for ( r=0 ; r<nrates ; r++){
	ProbMatArray[bran][r] = calloc(16,sizeof(double));     OOM(ProbMatArray[bran][r]);
	DProbMatArray[bran][r] = calloc(16,sizeof(double));    OOM(DProbMatArray[bran][
r]);
	DDProbMatArray[bran][r] = calloc(16,sizeof(double));   OOM(DDProbMatArray[bran][
r]);
      }
    }
  }

  if ( NOTMODE(GTR) ){
    for ( bran=0 ; bran<branches ; bran++){
      length = branch[bran]->length[0];
      for ( r=0 ; r<nrates ; r++){
	for ( from=0 ; from<4 ; from++){
	  for ( to=0 ; to<4 ; to++){
	    ProbMatArray[bran][r][from*4+to] = prob_calc(from,to,rate_cat[r] * length);
	    DProbMatArray[bran][r][from*4+to] = prob_calcd(from,to,rate_cat[r] * length);
	    DDProbMatArray[bran][r][from*4+to] = prob_calcd2(from,to,rate_cat[r] * length);
	  }
	}
      }
    }
  } else {
    for ( bran=0 ; bran<branches ; bran++){
      length = branch[bran]->length[0];
      for ( r=0 ; r<nrates ; r++){
	GetGTRmat(ProbMatArray[bran][r],rate_cat[r]*length);
	GetQP(DProbMatArray[bran][r],ProbMatArray[bran][r]);
	GetQP(DDProbMatArray[bran][r],DProbMatArray[bran][r]);
      }
    }
  }
    


}
