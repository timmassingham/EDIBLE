#define BUGS FEATURES
#define FEATURES on
#define CHILD(A) (node_p->A[n])
#define PARENT(B) (CHILD(node)->B[0])
#define ISMODE(A) (mode&A)!=0
#define NOTMODE(A) (mode&A)==0

/*  Fix complilation problems on MS platform*/
#define OFFSPRING 1

/*  Compiler directives - CHILD(node) is address of n th
 * child from a particular node.
 *                        CHILD(length) is length of branch
 * to child node n.
 *                        PARENT(node) is the "parent" node
 * for a particular node (except for the one originally
 * started at.
 *                        PARENT(length) is the length of
 * the branch to the particular node's "parent"*/


/*  Define the modes we are going to use (see edible.c)*/
#define MATRICES         1
#define PERCENTILE       2
#define ROOTED           4
#define INDIVIDUAL       8
#define PROBS           16
#define NODEASROOT      32
#define VARIANCE        64
#define TREES          128
#define BOOTSTRAP      256
#define HKY            512
#define NOKAPPA       1024
#define DETINDIV      2048
#define CACHE         4096
#define HAVE_RATES    8192
#define GTR	     16384

#define DOODAH 6 /* Maximum number of branches from any node
                  * eg: 3 for binary, 4 for trinary trees
                  * and so on till you run out of memory*/
#define LEAFMAX 30
#define BRANCHMAX 100
#define DELTAF(A,B) (double)((A==B)?1:0)
#define FACTOR (4.0/3)
#define DELTA 0.0001  /*  The value kappa is perturbed by in the FD
                       * approximation used to calculate the first and
                       * second derivatives of P by kappa*/

struct treenode{
  char name[16];   /*Node names must be less than 16 char.*/
  struct treenode *node[DOODAH];
  double length[DOODAH];
  double llh[4];
  int nucleotide;
  int bnum;
};

struct crecord{
  double **matrix;
  int *leaf_nucleotide;
};

/* Functions from read.c*/
void readtree(char *file, struct treenode *snodes);
void makenode(FILE *fp,struct treenode *node_p,int flag);
void printnode(struct treenode *node_p,int flag,FILE *fp);
void dump(double **matrix,int maxi,char *s);

/* Functions from llh.c*/
void calcllhs(struct treenode *node_p,struct treenode *tree2,unsigned int a, double **matrix);
double llhcalc(struct treenode *node_p);
void llhcalc_branch(struct treenode *node_p, int r);
void llhvector(struct treenode *node_p,int a,int derivative, int r);
double prob_calc(int from, int to, double length);
double prob_calcd(int from, int to, double length);
double prob_calcd2(int from, int to, double length);
double second_derivative(struct treenode *node_p,double *matrix,int b1,int b2);
double **pcalc(struct treenode *node_p);
void pcalc_branch(struct treenode *node_p);

/* Functions from utility.c*/
int ipow(int x, int y);
void nomemory(void);
char *itotext(int n,char *s);
double * determinant(void);
int ludecomp(double **matrix,int max);
char getnextc(FILE *fp);
void print_nucleotide(int nucleo, FILE *fp);
void initialise_rg(void);
void tree_sequence(unsigned int a);
double evaluate_information(double **matrix,const int a,const int b);
void do_rate(void);
void thick_sort(int array[],int length);

/* Functions from partial.c*/
void partialcalc(struct treenode *node_p, double *matrix, int twig);
void partialcalc_branch(struct treenode *node_p, int twig, int deriv, int r);
int **findspade(struct treenode *node_p);

/* Functions from matrix.c*/
void expectation(unsigned int e,double factor,int factor_flag,struct treenode *tree, struct treenode *tree2);
double *variance(unsigned int e);
void planttree(double **original_matrix,double **new_matrix);
void rooted_derivative(double *original_vector,double *new_vector);

/* Functions from options.c*/
void standard(struct treenode *node_p,FILE *file_p,unsigned int e);
void growtree(struct treenode *node_p,FILE *file_p,unsigned int e);
void growbranch(struct treenode *node_p,FILE *file_p,unsigned int e);
void greasebranch(struct treenode *node_p,FILE *file_p,unsigned int e);
double *find_information(struct treenode *tree, struct treenode *tree2, unsigned int e, int factor_flag, double factor);

/* Functions from tree.c*/
struct treenode *treecopy(struct treenode *node2_p,int flag);
void filltree(struct treenode *from, struct treenode *node_p,int flag);
void doleaf(struct treenode *node_p,int flag);
int findnode(struct treenode *node1_p,struct treenode *node2_p);
void scale_tree(int factor_flag,double factor, double **expect,int rows);
void print_tree(struct treenode *node_p, FILE *fp, int flag);

/* Functions from random.c*/
double sample_percentile(struct treenode *node_p,struct treenode *tree2, unsigned int e,int factor_flag,double factor);
double randomresult(struct treenode *node_p,struct treenode *tree,double **partial,unsigned int e,int sequence, int factor_flag, double factor);
double randomd(void);
long randoms(void);
void reorder(double *percent,double information,int percent_point);
int evolve_nucleotide(double length, int start_nucleotide);
double **check_cache(void);
void update_cache(double **matrix);
void wipe_cache(void);
