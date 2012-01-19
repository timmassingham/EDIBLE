#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "edible.h"
#include <string.h>

/*  Make copy of tree and return pointer to start node
 *  Again just recurses down the tree, copying
 * as it goes. Flag=0 at (code) root, flag=1 everywhere else*/
struct treenode *treecopy(struct treenode *node2_p,int flag){
  struct treenode *node_p;
  int a,n;
  char *string;
  
  /*  Get memory for new node of tree*/
  node_p=malloc(sizeof(struct treenode));
  if(node_p==NULL)
    nomemory();

  /* Copy over name*/
  string=strcpy(node_p->name,node2_p->name);
  node_p->bnum = node2_p->bnum;

  /* Cycle through all the branches leading to children,
   * recursing down if the branch actually exists*/
  for(n=flag;n<DOODAH && node2_p->node[n]!=NULL;n++){
    if((node2_p->node[n])!=NULL){
      CHILD(node)=treecopy(node2_p->node[n],1);
      PARENT(node)=node_p;
    }
  } 

  /* Blank out all remaining branches. ie: the ones which 
   * don't lead anywhere*/
  a=n;
  for(n=a;n<DOODAH;n++){
      CHILD(node)=NULL; 
      CHILD(length)=-1; 
  }
  return node_p;
}
 
/*  Fills in copy of tree which already exists.
 *  Work by recursing down tree. Flags as in treecopy
 * procedure*/
void filltree(struct treenode *from, struct treenode *node_p,int flag){
  int n;
  extern struct treenode *branch[];
  extern int branches;
 
  /*  If at root, then reset the branch counter*/ 
  if(flag==0)
    branches=0;
  /*  Copy nucleotide*/
  node_p->nucleotide=from->nucleotide;
  /*  Copy the "likelihood" vector across for reuse*/
  for(n=0;n<4;n++)
    node_p->llh[n]=from->llh[n];
  node_p->bnum = from->bnum;

  /*  If the node has children, then recurse down them
   * filling as we go*/
  for(n=flag;n<DOODAH && CHILD(node)!=NULL;n++){
    CHILD(length)=from->length[n];
    filltree(from->node[n],CHILD(node),1);
    PARENT(length)=CHILD(length);
    /*  Update the branch array to point at the new tree*/
    branch[branches]=CHILD(node);
    branches++;
  }
}  
   
/*  Routine to find all the leaves on a tree and fill
 * in the leaf[] array. Recurses down tree*/
void doleaf(struct treenode *node_p,int flag){
  extern int leaves;
  extern struct treenode *leaf[];
  int n;

  n=1;
  if(CHILD(node)==NULL)
    leaf[leaves++]=node_p;
  else
    for(n=flag;n<DOODAH && CHILD(node)!=NULL;n++)
      doleaf(CHILD(node),1);
}
 
 
/*  Routine to return n such that
 * (node_p->node[0])->node[n]=node_p*/
int findnode(struct treenode *node1_p,struct treenode *node2_p){
  int n;

  for(n=0;n<DOODAH && node1_p->node[n]!=node2_p;n++){}
  if(node1_p==node2_p)
    n=0;
  return n;
}


/*  Does the correct scaling to the information matrix, given 
which rates on the tree are changing*/
void scale_tree(int factor_flag, double factor, double **expect,int rows){
  int a,b;

  /*  Multiply necessary informations by correct factor */
  if(factor_flag==-1)/* Multiply all informations*/
    for(a=0;a<rows;a++)
      for(b=0;b<rows;b++)
        expect[a][b]*=factor*factor;
  else
    /*  We only require factor to be applied to one
     * information. If we are taking the determinant
     * factor applied to all in same row and column 
     * as information*/
    for(a=0;a<rows;a++){
      expect[factor_flag][a]*=factor;
      expect[a][factor_flag]*=factor;
    }
}


void print_tree(struct treenode *node_p, FILE *fp, int flag){
  int a,n;
  
  n=flag;
  while(n<DOODAH && CHILD(node)!=NULL){n++;}
  a=n;

  if(flag==0)
    fprintf(fp,"#");
  fprintf(fp,"(");
  for(n=flag;n<a;n++){
    if(CHILD(node)->node[1]!=NULL)
      print_tree(CHILD(node),fp,1);
    else
      fprintf(fp,"%s:%1.5f",CHILD(node)->name,CHILD(node)->length[0]);
    if(n<a-1)
      fprintf(fp,",");
  }

  fprintf(fp,")");
  if(flag!=0)
    fprintf(fp,":%1.5f",node_p->length[0]);
  else
    fprintf(fp,"\n");
}
  
