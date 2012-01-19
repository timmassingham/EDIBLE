#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include "edible.h"

/*  Opens file for input and calls makenode(*) to start
 * creating the tree*/
void readtree(char *file, struct treenode *snodes)
{
  FILE *input_file;
  char c;
  char *string;
  extern int individual;
  extern int interesting_branches[];
  extern int mode;
  extern int nodecount;
  extern struct treenode *branch[];
  extern int root;
  int a,b;
  
  /*  Try to open file*/
  input_file=fopen(file,"r");
  if(input_file==NULL){
    printf("Can't open file!!\n*Plink*\n");
    exit(2);
  }
  c=getnextc(input_file);

  /*  Call the (code) root node-0 and then recurse through the
   * file constructing the tree in memory*/
  string=strcpy(snodes->name,"Node-0");
  makenode(input_file,snodes,0);

  /*  If we are asked for information about an individual node,
   * set the correct flag and variable.
   * First do case in (rooted) tree when we want information about
   * the root node */
  if(ISMODE(ROOTED) && (c=getnextc(input_file))=='*'){
    if(ISMODE(NODEASROOT))
      interesting_branches[individual++]=root;
    else
      interesting_branches[individual++]=nodecount+1;
    if(NOTMODE(INDIVIDUAL))
      mode+=INDIVIDUAL;
  }
  if((mode&12)==12){ /* 12=individual information and rooted*/
    for(a=0;a<individual-1;a++)
      sscanf((branch[interesting_branches[a]]->node[0])->name,"Node-%d",&interesting_branches[a]);
  /*  Allow for information about the root node (different
   * number to the above scheme*/
    if(c!='*' || ISMODE(NODEASROOT))
      sscanf((branch[interesting_branches[individual-1]]->node[0])->name,"Node-%d",&interesting_branches[individual-1]);
  }

  /*  If we have a rooted tree, we may have multiple entries for each
   * parameter. Sort the array and remove duplicates.*/
  thick_sort(interesting_branches,individual);

  b=1;
  a=1;
  while(a<individual){
    while(a<individual && interesting_branches[a-1]==interesting_branches[a]){a++;}
    if(a<individual){
      interesting_branches[b]=interesting_branches[a++];
      b++;
    }
  }
  individual=b;

  fclose(input_file);
}



/*  Code to read data from the open file and fill out all
 * the information for one node. Calls itself recursively
 * when it meets the start of another node - creating
 * the memory space first.*/
void makenode(FILE *fp,struct treenode *node_p,int flag){
char l[16],m[16];
char c;
char *string;
int a, n;
struct treenode *cnode;
extern int nodecount;
extern int leaves;
extern struct treenode *leaf[];
extern struct treenode *branch[];
extern int branches;
extern int root;
extern int individual;
extern int mode;
extern int interesting_branches[];

n=flag;
c=' '; /*  Initialise c, so it is not ')'*/

/*  Read branch from file and recurse down it until we reach the
 * end, creating the tree as we go*/
while(c!=')'){   /*  Until this branch ends {file format has ')'}*/
  c=getnextc(fp);

  /*  Get memory*/
  cnode=malloc(sizeof(struct treenode));
  if(cnode==NULL)
    nomemory();

  if(c=='('){              /*  Check if new node is attached to current one*/
    makenode(fp,cnode,1);  /* and call recursively.*/
    nodecount++;
    c='0';
  }
  CHILD(node)=cnode;

  if(c=='0'){    /*  Think of suitable name for internal node*/
    string=strcpy(m,"Node-");
    string=itotext(nodecount,l);
    string=strcat(m,l);
    string=strcpy(CHILD(node)->name,m);
    c=getnextc(fp);
  }
  else{          /*  Or read name from the file*/
    l[0]=c;
    a=1;
    while((c=getnextc(fp))!=':' && c!='(' && c!=')' && c!='*')
      if(a<15)
        l[a++]=c;
    l[a]='\0';
    sscanf(l,"%s",CHILD(node)->name);
    leaf[leaves]=CHILD(node);
    leaves++;
    for(a=1;a<DOODAH;a++)
      CHILD(node)->node[a]=NULL;
  }

  a=0;
  if(c!=':' && c!='*'){    /*  Be pedantic about file format*/
    printf("Error in file format %c - can't find length\n*Plink*\n",c);
    exit(2);
  }

  /*  If we want information about this branch only, set the flag*/
  if(c=='*'){
    if(NOTMODE(INDIVIDUAL))
      mode+=INDIVIDUAL;
    interesting_branches[individual++]=branches;
  }

  /*  Read in length of branch to new node*/
  while((c=getnextc(fp))!=',' && c!='R' && c!='r' && c!='(' && c!=')')
    l[a++]=c;
  l[a]='\0';
  sscanf(l,"%le",&CHILD(length));

/*  Add parent to node just created
 * and also save length, since memory has been set aside
 * and the code becomes less obscure*/
  PARENT(length)=CHILD(length);
  PARENT(node)=node_p;
  if(c=='R' || c=='r'){
    root=branches;
    mode+=ROOTED;
    if(c=='r') /*  If the root is an actual node, rather than a*/
      mode+=NODEASROOT;/* branch then flag*/
    c=getnextc(fp);
  }
  CHILD(node)->bnum = branches;
  branch[branches++]=CHILD(node);
  n++;
}
/*  Blank out the remaining pointers, so they aren't
 * checked at a later date*/
a=n;
for(n=a;n<DOODAH;n++){
  CHILD(length)=-1;
  CHILD(node)=NULL;
}
}



/*     For testing(?) - prints out tree     */
void printnode(struct treenode *node_p,int flag,FILE *fp){
  int n;

  n=flag;
  while(n<DOODAH && CHILD(node)!=NULL){
    fprintf(fp,"# %s\tConnected to %15s by Length %e\n",node_p->name
                ,CHILD(node)->name,CHILD(length));

/*  When calling parent from child - be careful of
 * infinite recursion!
 *  The first branch always the parent, except for
 * the tree root, in which case, the first branch
 * is also a child.*/ 
    printnode(CHILD(node),1,fp);
    n++;
  }
}

/*  Routine to dump a matrix to matrix output file   */
void dump(double **matrix,int max,char *s){
  int a,b;
  extern FILE *matrix_file_p;

  fprintf(matrix_file_p,"%s\n",s);
  for(a=0;a<max;a++){
    for(b=0;b<max;b++)
      fprintf(matrix_file_p,"%1.11E,",matrix[a][b]);
    fprintf(matrix_file_p,"\n");
  }
}
