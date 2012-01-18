  struct treenode *tree2;
/* double (*(*P)[])[] = P is pointer to an array of pointers to doubles
 *        (*P)[a] = is a pointer to an array of doubles (a is int)
 *        (*(*P)[a])[b] = is a double, (a,b) in the array*/

  double *det;
  double factor;
  int factor_flag;
  
/*  0 don't include factor in calculations.
 *  1 scaling entire tree, include factor in all calculations.
 * about that branch.*/

  extern int leaves;
  extern int mode;
  extern int branches;
  extern struct treenode *branch[];
  extern FILE *matrix_file_p;
  extern FILE *variance_file_p;
  extern FILE *sample_file_p;
  extern FILE *prob_file_p;
  extern char *outstring;
  extern int individual;
  extern int interesting_branches[];
