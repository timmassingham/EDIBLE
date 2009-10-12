/*  Program to calculate the likelihood of a given tree
 * evolving the nucleotides entered.
 *  WARNING: If DOODAH is exceeded - the program may still
 * run, giving an erroneous answer*/


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "edible.h"
#include "gtr.h"
#include <string.h>
#include <time.h>
#include <float.h>

#define FOOBAR main /*  Calling an important function FOOBAR
                     * is apparently a shooting offence*/

int nodecount=0;
int leaves=0;
int branches=0;
int is_kappa=0;
struct treenode *leaf[LEAFMAX];
struct treenode *branch[BRANCHMAX];
char *out_file;
FILE *matrix_file_p;
FILE *prob_file_p;
FILE *variance_file_p;
FILE *sample_file_p;
int sequence_length, sample_size;
int boot_strap_size;
int cache_size=0;
struct crecord (*cache)[];

double percentile;
int mode=0;
/*  Value automatically altered by program
 *  Modes: +1 dumps intermediates to file, +2 turns on sampling
 *         +4 allows rooted tree,          +8 want individual value
 *         +16 dump probabilities to file, +32 rooted behaviour=node
 *         +64 calculate variances         +128 print all trees
 *         +256 use bootstrap method       +512 use HKY85 model
 *         +1024 don't calculate kappa     +2048 determinate of param.
 *         +4096 use cache when sampling*/
int root=0;
int individual=0;
int (*(*conv_matrix)[])[];
double (*(*expect)[])[];
double (*(*var)[])[];
double (*(*var2)[])[];
double (*(*rootedexpect)[])[];
char *outstring;
int seed;
int kudge=0,branch1,branch2;
int interesting_branches[BRANCHMAX];

double p[4];
double kappa;
double rate;

int nrates=0;
double * rate_prob=NULL, *rate_cat=NULL;
double *** ProbMatArray=NULL, *** DProbMatArray=NULL, *** DDProbMatArray=NULL;
double * gtr=NULL;

int FOOBAR(argc, argv)
int argc;
char *argv[];
{
struct treenode snode;
unsigned int e;
int a,b;
extern int boot_strap_size;
extern int sample_size;
extern int sequence_length;
extern int cache_size;
extern struct crecord (*cache)[];
extern int leaves;
extern double percentile;
extern char *out_file;
extern FILE *matrix_file_p;
extern FILE *prob_file_p;
extern FILE *variance_file_p;
extern FILE *sample_file_p;
extern double (*(*expect)[])[];
extern double (*(*var)[])[];
extern double (*(*var2)[])[];
extern double (*(*rootedexpect)[])[];
extern char *outstring;
extern struct treenode *leaf[];
extern int is_kappa;
extern int interesting_branches[];

extern double p[4];
extern double kappa;

char *matrix_file;
char *prob_file;
char *variance_file;
char *string;
char c;
FILE *file_p;
int r;
double tot;

a=1;

/*  Set up parameters from the JC model*/
p[0]=0.25;
p[1]=0.25;
p[2]=0.25;
p[3]=0.25;
kappa=1;
do_rate();

/*  Try to read all command line options and fail if there
 * are some we don't understand*/
while(a<argc && *argv[a]=='-'){
  argv[a]++;
  
  switch(*argv[a]){

  /*  Calculate variances.*/
  case 'v':a++;
           mode+=VARIANCE;
           variance_file=argv[a++];
           variance_file_p=fopen(variance_file,"w");
           if(variance_file_p==NULL){
             printf("Can't open file to dump variances to\n");
             exit(2);
           }
           break;

  /*  Dump intermediate matrices to a file*/
  case 'm':a++;
           matrix_file=argv[a++];
           matrix_file_p=fopen(matrix_file,"w");
           if(matrix_file_p==NULL){
             printf("Can't open file to dump intermediate matrices to\n");
             exit(2);
           }   
           mode+=MATRICES; 
           break;

  /*  Dump the probabilities for each nucleotide sequence
   * to a file*/
  case 'p':a++;  
           prob_file=argv[a++];
           prob_file_p=fopen(prob_file,"w");
           if(prob_file_p==NULL){
             printf("Can't open file to dump probabilities to\n");
             exit(2);
           }
           mode+=PROBS;
           fprintf(prob_file_p,"Probabilities for tree\n");
           fprintf(prob_file_p,"Nucleotides\tProbabillity\tInformation\t1st Partial\t2nd Partial\n");
           break;

  /*  Turn on sampling of the distribution to estimate the nth
   * percentile.*/
  case 's':a++;
           c=argv[a][0];
           sample_file_p=NULL;
           /* sample_file may not start with a digit or '-'*/
           if(isdigit(c)==0 && c!='-'){
	     sample_file_p=fopen(argv[a++],"w");
	     if(sample_file_p==NULL){
	       printf("Can't open file to dump sampled distribution to\n");
	       exit(2);
	     }
	   }
           sscanf(argv[a++],"%d",&sample_size);
           sscanf(argv[a++],"%d",&sequence_length);
           sscanf(argv[a++],"%le",&percentile);
           if(sample_size<1 || sequence_length<1 || percentile<0 || percentile>100){
             printf("Percentile parameters are invalid!\n");
             exit(0);
           }
           percentile/=100;
           mode+=PERCENTILE;
           initialise_rg();
          
           if(sample_file_p!=NULL){
	     fprintf(sample_file_p,"#Sampled distribution, size %d based on a sequence of length %d\n",sample_size,sequence_length);
	   }
           break;

  /* Print out all the trees we use (useful when performing several
   * calculations while manipulating the tree*/
  case 't':a++;
           mode+=TREES;
           break;
  /*  We want to use boot strap techniques - equivalent to straight
   * forward estimation in all cases apart from trying to estimate 
   * the percentile*/
  case 'b':a++;
           sscanf(argv[a++],"%d",&boot_strap_size);
	   mode+=BOOTSTRAP;
           initialise_rg();
	   break;
  /*  Enable HKY85 model of nucleotide substitution, along with
   * information calculations about kappa. The 'k' option 
   * disables calculations about kappa*/
  case 'h':a++;
	   mode+=HKY;
	   is_kappa=1;
	   sscanf(argv[a++],"%le",&p[0]);
	   sscanf(argv[a++],"%le",&p[1]);
	   sscanf(argv[a++],"%le",&p[2]);
	   sscanf(argv[a++],"%le",&p[3]);
	   sscanf(argv[a++],"%le",&kappa);
	   
	   p[3]=1-(p[0]+p[1]+p[2]);
	   if(p[0]<0 || p[1]<0 || p[2]<0 || p[3]<0){
	     printf("Probabilities are negative - terminating\n");
	     exit(2);
	   }
           do_rate();
	   break;
  /*  Disable calculations about Kappa in the HKY model.
   * Does nothing in the JC case*/
  case 'k':a++;
           mode+=NOKAPPA;
           break;
  /*  Work with the determinate of all the values specified,
   * rather than outputting individual informations*/
  case 'd':a++;
           mode+=DETINDIV;
           break;
  /*  Add a result cache for sampling purposes*/
  case 'c':a++;
           mode+=CACHE;
           sscanf(argv[a++],"%d",&cache_size);
           if(cache_size<0){
	     printf("Negative cache! Exiting\n");
	     exit(0);
	   }
           cache=calloc(cache_size,sizeof(struct crecord));
           if(cache==NULL)
             nomemory();
           break;
  case 'r': a++;
	    mode+=HAVE_RATES;
            sscanf(argv[a++],"%d",&nrates);
            if (nrates<1){
              printf ("Less than one rate category!");
              exit(EXIT_FAILURE);
            }  
            rate_prob = calloc(nrates,sizeof(double));
            rate_cat = calloc(nrates,sizeof(double));
            for ( r=0 ; r<nrates ; r++){
              sscanf(argv[a++],"%le",&rate_cat[r]);
              if ( rate_cat[r]<=0.){
                printf ("Rate is less than zero (%e).\n",rate_cat[r]);
                exit(EXIT_FAILURE);
              }
            }
            tot = 0.;
            for ( r=0 ; r<nrates ; r++){
              sscanf(argv[a++],"%le",&rate_prob[r]);
              tot += rate_prob[r];
              if ( rate_prob[r]<=0. || rate_prob[r]>=1.){
                printf ("Error in prior probability of rate %e. Probability is %e\n",rate_cat[r],rate_prob[r]);
                exit(EXIT_FAILURE);
              }
            }
            if ( fabs(1.-tot)>DBL_EPSILON*nrates){
              printf ("Rate prior probabilities do not sum to one.\n");
              exit(EXIT_FAILURE);
            }
            
            break;
  case 'g': a++;
	    mode += GTR;
	    gtr = calloc(5,sizeof(double));
	    for ( r=0 ; r<4 ; r++){
	      sscanf(argv[a++],"%le",&p[r]);
	    }
	    for ( r=0 ; r<5 ; r++){
	      sscanf(argv[a++],"%le",&gtr[r]);
	    }
	    CreateGTR();
            break;
  }
}

/* We have arguments left over, something nasty has happened!*/
if(a!=argc-2){
  printf("\tWrong number of arguments!\n"
         "(Non-optimized version)\n\n"
         "  Usage:  edible <-s [sample_file] sample_size sequence_length percentile>"   
         "\n\t<-m matrix_file> <-p prob_file> <-t> <-b boot_strap_size> <-d>\n"
	 "\t <-h pi_a pi_c pi_g pi_t kappa> <-g pi_a pi_c pi_g pi_t GTR_parameters>\n"
    "\t<-r ncat rates probs> <-k> <-c cache_size> tree_file output_file \n\n"
         "  matrix_file is file to dump expectation matrices to.\n"
         "  prob_file is file to dump probabilities to.\n"
         "  sample_file is file to dump sampled distribution to. Optional.\n"
         "  sample_size for use when sampling for x%% level.\n"
         "  sequence_length: length of sequences used in sampling.\n"
         "  percentile: percentile of distribution to calculate, x%%.\n"
         "  option t causes intermediate trees to be printed\n"
	 "  boot_strap_size is the size of the sample to use in estimations\n"
	 "  -d want determinate of information about several parameters specified\n"
	 "  -h uses the HKY85 model rather than the JC model of nucleotide substitution\n"
	 "pi_? is the base frequency of the nucleotide ?. Kappa is the transition bias\n"
    "  -g use GTR model of nucleotide substitution, see documentation for details\n"
    "  -r use ncat categories of rate varation, rates and probabilities given\n"
    "  -k don't calculate the information WRT kappa in the HKY85 model\n"
	 "  cache_size is the amount of sequences to store results for when sampling\n\n");
  exit(2);
}

/*  If rate not given, use one rate category*/
if ( NULL==rate_cat && NULL==rate_prob){
  rate_cat = malloc(sizeof(double));
  rate_cat[0] = 1.;
  rate_prob = malloc(sizeof(double));
  rate_prob[0] = 1.;
  nrates = 1;
}
/*  Many options have been broken by adding rate code*/
if (ISMODE(BOOTSTRAP) || ISMODE(PERCENTILE) || (ISMODE(HKY) && NOTMODE(NOKAPPA))){
  puts ("Rates are not compatible with bootstraps, percentiles, or HKY with kappa");
  exit(EXIT_FAILURE);
}

  

if(ISMODE(PERCENTILE) && ISMODE(BOOTSTRAP))
  printf("Don't want to sample twice - doing percentile calculations\n");
  
out_file=argv[a+1];
readtree(argv[a],&snode);

/*  If we want to cache results then get memory for the cache*/
if(ISMODE(CACHE)){
  for(a=0;a<cache_size;a++){
    (*cache)[a].leaf_nucleotide=calloc(leaves,sizeof(int));
    for(b=0;b<leaves;b++)
      (*(*cache)[a].leaf_nucleotide)[b]=-1;
    (*cache)[a].matrix=calloc(branches+1+is_kappa,sizeof(double *));
    if((*cache)[a].matrix==NULL)
      nomemory();
    for(b=0;b<branches+1+is_kappa;b++){
      (*(*cache)[a].matrix)[b]=calloc(branches+1+is_kappa,sizeof(double));
      if((*(*cache)[a].matrix)[b]==NULL)
	nomemory();
    }
  }
}
  
/*  If we are sampling, we can only order by one parameter
 * complain if we have more*/
if(ISMODE(PERCENTILE) && ISMODE(INDIVIDUAL) && individual>1 && NOTMODE(DETINDIV)){
  printf("  Argh - total breakdown of order. A partial order just isn't good\n"
         "enough to keep the peace during percentile calculations. Please specify one\n"
         "parameter (or determinate) so we can tame likelihood land\n"
         "  This town is just too dangerous, I'm leaving.\n");
  exit(0);
}

/* Offer choice of tree manipulations*/
printf("\nOptions:  1. Calculate the expected information of the tree\n"
       "          2. Scale tree by a given range of factors\n"
       "          3. Scale a single branch by a given range of factors\n"
       "          4. Slide a given branch along to others for given lengths\n"
       "            Please choose:");

  /*  Consider how much memory we need - if we are calculating all
   * the probabilities then we need lots. If we are using boot strap
   * techniques, we are only going to store the boot strap sample*/
  e=ipow(4,leaves);
  if(ISMODE(BOOTSTRAP))
      e=boot_strap_size;

  /*  Get memory for the expectation matrix*/
  expect=calloc(branches+is_kappa,sizeof(double *));
  if(expect==NULL)
    nomemory();
  for(a=0;a<branches+is_kappa;a++){
    (*expect)[a]=calloc(branches+is_kappa,sizeof(double));
    if((*expect)[a]==NULL)
      nomemory();
  }
  /*  Get memory for the variance matrix, if we are doing variance
   *  calculations*/
  if(ISMODE(VARIANCE) && NOTMODE(ROOTED)){
    var=calloc(branches+is_kappa,sizeof(double *));
    if(var==NULL)
      nomemory();
    for(a=0;a<branches+is_kappa;a++){
      (*var)[a]=calloc(branches+is_kappa,sizeof(double));
      if((*var)[a]==NULL)
	nomemory();
    }
  }

  /*  Open the output file*/
  file_p=fopen(out_file,"w");

  /*  Depending on options, get any additional memory needed
   * and write correct comments to files*/
  if(ISMODE(ROOTED)){
    conv_matrix=findspade(&snode); /*  memory needed for conversion matrix
                                    * if tree is rooted*/
    b=((mode&32)!=0)?1:2;

    /*  Get memory for the new expectation matrix after the conversion
     * matrix has been applied*/
    rootedexpect=calloc(nodecount+b+is_kappa,sizeof(double *));
    if(rootedexpect==NULL)
      nomemory();
    for(a=0;a<nodecount+b+is_kappa;a++){
      (*rootedexpect)[a]=calloc(nodecount+b+is_kappa,sizeof(double));
      if((*rootedexpect)[a]==NULL)
        nomemory();
    }

    /*  If we are doing a variance calculation, then get the memory for
     * the rooted case (we have to convert each entry as we go along
     * - slow!*/
    if(ISMODE(VARIANCE)){
      var=calloc(nodecount+b+is_kappa,sizeof(double* ));
      if(var==NULL)
        nomemory();
      for(a=0;a<nodecount+b+is_kappa;a++){
        (*var)[a]=calloc(nodecount+b+is_kappa,sizeof(double));
        if((*var)[a]==NULL)
          nomemory();
      }
      var2=calloc(nodecount+b+is_kappa,sizeof(double* ));
      if(var2==NULL)
        nomemory();
      for(a=0;a<nodecount+b+is_kappa;a++){
        (*var2)[a]=calloc(nodecount+b+is_kappa,sizeof(double));
        if((*var2)[a]==NULL)
          nomemory();
      }
    }

    /*  Reset the is_kappa variable for use in the main loops*/
    is_kappa=0;

    /*  Output some useful information to the open files*/
    for(a=0;a<individual;a++)
      sscanf((branch[interesting_branches[a]]->node[0])->name,"node-%d",&interesting_branches[a]);
    printnode(&snode,0,file_p);
    if(ISMODE(MATRICES)){
      if(ISMODE(INDIVIDUAL))
        fprintf(matrix_file_p,"#Zeros are uncalculated!\n");
      printnode(&snode,0,matrix_file_p);
      for(a=0;a<nodecount+1;a++)
        fprintf(matrix_file_p,"#Parameter %d measures evolutionary time from "
                              "node-%d to present day\n",a,a);
      if(NOTMODE(NODEASROOT)){
        fprintf(matrix_file_p,"#Parameter %d measures evolutionary time from"
                              " the root node to present day\n",nodecount+1);
      }
    }
      for(a=0;a<nodecount+1;a++)  
        fprintf(file_p,"#Parameter %d measures evolutionary time from "
                       "node-%d to present day\n",a,a);
      if(NOTMODE(NODEASROOT)){
        fprintf(file_p,"#Parameter %d measures evolutionary time from "
                       "the root node to present day\n",nodecount+1);
      }
  }
  else {
    printnode(&snode,0,file_p);
    if(ISMODE(MATRICES)){
      if(ISMODE(INDIVIDUAL))
        fprintf(matrix_file_p,"#Zeros are uncalculated!\n");
      for(a=0;a<branches;a++)
        fprintf(matrix_file_p,"#Branch %d is between %s and %s\n",a
                      ,branch[a]->name,(branch[a]->node[0])->name);
    }
    for(a=0;a<branches;a++)
      fprintf(file_p,"#Branch %d is between %s and %s\n",a,branch[a]->name
                    ,(branch[a]->node[0])->name);
  }

  /*  Add a bit more information about the options chosen to the output file*/
  if(ISMODE(INDIVIDUAL))
    for(a=0;a<individual;a++)
      fprintf(file_p,"#Doing parameter/branch %d\n",interesting_branches[a]);
  if(ISMODE(PERCENTILE))
    fprintf(file_p,"#Sampling for the %e percentile of the distribution\n",percentile*100);
  if(ISMODE(PROBS)){
    for(a=0;a<leaves;a++)
      fprintf(prob_file_p,"#Leaf %d is %s\n",a,leaf[a]->name);
    if(ISMODE(INDIVIDUAL) && individual==1){
      if(ISMODE(ROOTED))
        fprintf(prob_file_p,"Derivatives WRT parameter %d\n",interesting_branches[0]);
      else
        fprintf(prob_file_p,"Derivatives WRT branch %d\n",interesting_branches[0]);
    }
  }

  /*  Depending on the options chosen, we are looking at different parameters*/
  switch((mode&(DETINDIV+INDIVIDUAL+ROOTED))){
  case 0:
  case ROOTED:
  case DETINDIV:outstring=calloc(34,sizeof(char));
	    string=strcpy(outstring,"Determinant of expectation matrix");
	    break;
  case INDIVIDUAL:   outstring=calloc(25,sizeof(char));
	    string=strcpy(outstring,"Information about branch");
	    break;
  case ROOTED+INDIVIDUAL:  outstring=calloc(28,sizeof(char));
	    string=strcpy(outstring,"Information about parameter");
	    break;
  case DETINDIV+INDIVIDUAL:outstring=calloc(32,sizeof(char));
	    string=strcpy(outstring,"Determinate of several branches");
	    break;
  case DETINDIV+INDIVIDUAL+ROOTED:outstring=calloc(34,sizeof(char));
	    string=strcpy(outstring,"Determinate of several parameters");
	    break;
  }

/*  Choose option required.*/
a=0;
scanf("%d",&a);  
while(a<1 || a>4){
  printf("Invalid. Please choose again: ");
  scanf("%d",&a);
}
 
/*  Go to the relevant procedure, depending on option*/
switch(a){
  case 1: standard(&snode,file_p,e); 
          break;
  case 2: growtree(&snode,file_p,e);
          break;
  case 3: growbranch(&snode,file_p,e);
          break;
  case 4: greasebranch(&snode,file_p,e);
          break;
}

/*  Close all the files we are using (unnecessary, since the compiler should
 * do it for us, but...*/
fclose(file_p);
if(ISMODE(MATRICES))
  fclose(matrix_file_p);
if(ISMODE(PROBS))
  fclose(prob_file_p);
if(ISMODE(VARIANCE))
  fclose(variance_file_p);

return 0;
}
