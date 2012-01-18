#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "edible.h"
#include "new_models.h"
#include <time.h>

/*  Standard option - calculate the probability, the matrix of
 * expected information and the diagonal values of the variance
 * Dump these to disk and calculate the determinant of the
 * expectation matrix.*/
void standard(struct treenode *node_p,FILE *file_p,unsigned int e){
#include "variables.h"
  int a,b;

  /*  Make a copy of the tree, fill and calculate the likelihoods*/
  tree2=treecopy(node_p,0);
  filltree(node_p,tree2,0);

  if(sample_file_p!=NULL)
    fprintf(sample_file_p,"#Points 1\n\n");
  
  /*  We aren't altering the rate, so factor is 1*/
  factor=1.0;
  factor_flag=0;

  det=find_information(node_p,tree2,e,factor_flag,factor);

  /*  Print the required output - may be more than one value*/
  b=1;
  if(NOTMODE(DETINDIV) && ISMODE(INDIVIDUAL))
    b=individual;
  
  for(a=0;a<b;a++)
    fprintf(file_p,"%s = %1.11E\n",outstring,det[a]);
  if(ISMODE(HKY) && NOTMODE(NOKAPPA) && NOTMODE(PERCENTILE))
    fprintf(file_p,"Information about Kappa = %1.11E\n",det[b]);

  free(det);
}




/*  Routine to vary the length of all branches by the same factor
 * and dump the expected information to disk. Since we are altering 
 * the rate of change, we need to include a factor in all the information
 * calculations*/
void growtree(struct treenode *node_p,FILE *file_p,unsigned int e){
#include "variables.h"
  int steps,a,b,c;
  
  double f_min,f_max,f;

  struct treenode *tree;


/*  Code to get the factors from the users*/
  printf("Altering length of all branches by same factor\n\n");
  printf("Please enter min. factor, max factor and number of points:\n");
  do{
    printf("\tmin. factor:");
    scanf("%lf",&f_min);
  }while(f_min<=0);
  do{
    printf("\tmax. factor:");
    scanf("%lf",&f_max);
  }while(f_min>=f_max);
  do{
    printf("\tnumber of points:");
    scanf("%d",&steps);
  }while(steps<=0);

  /*  Dump text to file explaining what we are about to do*/
  fprintf(file_p,"#Results from Varying length of all branches by factor\n");
  fprintf(file_p,"#Factor\t\t\t%s\n",outstring);
  if(ISMODE(MATRICES))
    fprintf(matrix_file_p,"#Matrices from Varying length of all branches by factor\n");
  if(sample_file_p!=NULL){
    fprintf(sample_file_p,"#Points %d\n",steps);
    fprintf(sample_file_p,"#Samples from varying length of all branches by factor\n");
  }
  if(ISMODE(PROBS))
    fprintf(prob_file_p,"#Probabilities from varying length of all branches by factor\n");
  
/*  Start main loop*/

  /*  Make two copies of the tree*/
  tree2=treecopy(node_p,0);
  tree=treecopy(node_p,0);

  /*  Calculate what factor we need to include in the calculations*/
  for(c=0;c<(steps);c++){
    if(steps==1)  /* Prevent divide by zero if only one step*/
      f=f_min;
    else
      f=f_min+c*(f_max-f_min)/(steps-1);
    factor=f;
    factor_flag=-1; /*  Flag = -1 means that we want to include the
                     * factor on all branches*/

    /*  Add all the information to one copy from the original tree*/
    filltree(node_p,tree,0);

/*  Change all the length in the tree to the scaled version
 * The array branch[] currently points to those of 'tree'*/
    for(b=0;b<branches;b++){
      branch[b]->length[0]=branch[b]->length[0]*f;
      a=findnode(branch[b]->node[0],branch[b]);
      (branch[b]->node[0])->length[a]=branch[b]->length[0];
    }
    /*  Make copy of the tree with scaled branches*/
    filltree(tree,tree2,0);

    /*  Branch[] now points to those in tree2
     * Make leaves point to those on tree*/
    leaves=0;
    doleaf(tree,0);

    printf("Doing factor %E\n",f);

    /*  If we've been asked to print out the intermediate trees, then
     * dump them to all the open files*/
    if(ISMODE(TREES))
      print_tree(tree,file_p,0);
    fprintf(file_p,"%E\n",f);
    if(ISMODE(MATRICES)){
      fprintf(matrix_file_p,"\n#Factor %E\n",f);
      if(ISMODE(TREES))
        print_tree(tree,matrix_file_p,0);
    }
    if(ISMODE(VARIANCE)){
      fprintf(variance_file_p,"\n#Factor %E\n",f);
      if(ISMODE(TREES))
        print_tree(tree,variance_file_p,0);
    }
    if(sample_file_p!=NULL){
      fprintf(sample_file_p,"\n#Factor %E\n",f);
      if(ISMODE(TREES))
	print_tree(tree,sample_file_p,0);
    }
    if(ISMODE(PROBS)){
      fprintf(prob_file_p,"\n#Factor %E\n",f);
      if(ISMODE(TREES))
	print_tree(tree,prob_file_p,0);
    }
    det=find_information(tree,tree2,e,factor_flag,factor);

    /*  Possible may want to print out more than one result*/
    b=1;
    if(ISMODE(INDIVIDUAL))
      b=individual;
    
    if(ISMODE(INDIVIDUAL)){
      if(ISMODE(DETINDIV)){
	fprintf(file_p,"\t\t\t%E\tD(",det[0]);
	for(a=0;a<b-1;a++)
	  fprintf(file_p,"%d,",interesting_branches[a]);
	fprintf(file_p,"%d)\n",interesting_branches[a]);
      }else
	for(a=0;a<b;a++)
	  fprintf(file_p,"\t\t\t%E\t%d\n",det[a],interesting_branches[a]);
    }else
      fprintf(file_p,"\t\t\t%E\tD\n",det[0]);
    
    if(ISMODE(HKY) && NOTMODE(NOKAPPA) && NOTMODE(PERCENTILE))
      fprintf(file_p,"\t\t\t%E\tKappa\n",det[b]);
    free(det);
  }
}



/*  Routine to vary the length of one branch though several
 * multipliers. This has little meaning when considering a clock-like
 * (rooted) tree*/
void growbranch(struct treenode *node_p,FILE *file_p,unsigned int e){
#include "variables.h"
  int steps,a,b,c;
  int elastic;

  double f_min,f_max,f;

  struct treenode *tree;

/*  Test to see whether tree is rooted - changing branch 
 * length in unrooted trees doesn't have much experimental
 * meaning                                                  */
  if(((mode&4)|(mode&32))!=0)
    printf("\a\n**Changing branch length in a rooted tree doesn't"
           " have meaning\n**in an experimental design problem!\n\n\a");

  /*  Make two copies of the tree*/
  tree2=treecopy(node_p,0);
  tree=treecopy(node_p,0);
  /*  Fill in the first*/
  filltree(node_p,tree,0);

  /*  Choose branch to stretch - better be elastic*/
  printf("\n");
  for(a=0;a<branches;a++)
    printf("Branch %d goes from %s to %s\n",a
           ,(branch[a]->node[0])->name,branch[a]->name);
  do{
    printf("\nWhich branch would you like to alter?");
    scanf("%d",&elastic);
  }while(elastic<0 || elastic>branches);

  /*  Get length scaling factors from user*/
  printf("Altering length of a branch leading to"
                    " %s by a varying factor\n\n",branch[elastic]->name);
  printf("Please enter min. factor, max factor and number of points:\n");
  do{
    printf("\tmin. factor:");
    scanf("%lf",&f_min);
  }while(f_min<=0);
  do{
    printf("\tmax. factor:");
    scanf("%lf",&f_max);
  }while(f_min>=f_max);
  do{
    printf("\tnumber of points:");
    scanf("%d",&steps);
  }while(steps<=0);

  /*  Dump the necessary information to file*/
  fprintf(file_p,"#Results from varying length of branch"
	  " between %s and %s by factor\n"
	  ,(branch[elastic]->node[0])->name,branch[elastic]->name);
  fprintf(file_p,"#Factor\t\t\t%s\n",outstring);
  if(ISMODE(MATRICES))
    fprintf(matrix_file_p,"#Matrices from varying length of branch"
	    " between %s and %s by factor\n"
	    ,(branch[elastic]->node[0])->name,branch[elastic]->name);
  if(sample_file_p!=NULL){
    fprintf(sample_file_p,"#Points %d\n",steps);
    fprintf(sample_file_p,"#Samples from varying length of branch"
	    " between %s and %s by factor\n"
	    ,(branch[elastic]->node[0])->name,branch[elastic]->name);
  }
  if(ISMODE(PROBS))
    fprintf(prob_file_p,"#Probabilities from varying length of branch"
	    " between %s and %s by factor\n"
	    ,(branch[elastic]->node[0])->name,branch[elastic]->name);
  
  /*  Main loop to stretch the branch*/
  for(c=0;c<steps;c++){
    /*  Calculate the factor wanted and set the flag - we only
     * want to apply the factor to one of the informations this
     * time*/
    if(steps==1)
      f=f_min;
    else
      f=f_min+c*(f_max-f_min)/(steps-1);

    factor=f;
    factor_flag=elastic;

    /*  Fill in first copy from the original*/
    filltree(node_p,tree,0);

   /*  Branch[] points to those of 'tree'
    * Change length of the branch by factor - need to find the node
    * that points from the other side of the branch.*/
    branch[elastic]->length[0]=branch[elastic]->length[0]*f;
    a=findnode(branch[elastic]->node[0],branch[elastic]);
    (branch[elastic]->node[0])->length[a]=branch[elastic]->length[0];

    /*  Fill in the second copy from the first (with the new branch 
     * length and set all the leaves to point to "tree"*/
    filltree(tree,tree2,0);
    leaves=0;
    doleaf(tree,0);

    printf("Doing Factor %E\n",f);
    /*  If we are printing out the intermediate trees, then dump them
     * to every file*/
    if(ISMODE(TREES))
      print_tree(tree,file_p,0);
    fprintf(file_p,"%E\n",f);
    if(ISMODE(MATRICES)){
      fprintf(matrix_file_p,"\n#Factor %E\n",f);
      if(ISMODE(TREES))
        print_tree(tree,matrix_file_p,0);
    }
    if(ISMODE(VARIANCE)){
      fprintf(variance_file_p,"\n#Factor %E\n",f);
      if(ISMODE(TREES))
        print_tree(tree,variance_file_p,0);
    }
    if(sample_file_p!=NULL){
      fprintf(sample_file_p,"\n#Factor %E\n",f);
      if(ISMODE(TREES))
	print_tree(tree,sample_file_p,0);
    }
    if(ISMODE(PROBS)){
      fprintf(prob_file_p,"\n#Factor %E\n",f);
      if(ISMODE(TREES))
	print_tree(tree,prob_file_p,0);
    }
    
    det=find_information(tree,tree2,e,factor_flag,factor);

    /*  Possible may want to print out more than one result*/
    b=1;
    if(ISMODE(INDIVIDUAL))
      b=individual;
    
    if(ISMODE(INDIVIDUAL)){
      if(ISMODE(DETINDIV)){
	fprintf(file_p,"\t\t\t%E\tD(",det[0]);
	for(a=0;a<b-1;a++)
	  fprintf(file_p,"%d,",interesting_branches[a]);
	fprintf(file_p,"%d)\n",interesting_branches[a]);
      }else
	for(a=0;a<b;a++)
	  fprintf(file_p,"\t\t\t%E\t%d\n",det[a],interesting_branches[a]);
    }else
      fprintf(file_p,"\t\t\t%E\tD\n",det[0]);
    
    if(ISMODE(HKY) && NOTMODE(NOKAPPA) && NOTMODE(PERCENTILE))
      fprintf(file_p,"\t\t\t%E\tKappa\n",det[b]);
    free(det);
  }
}


/*  Routine to slid one branch of a node along the a branch
 *  The moving node is assumed to be simple!*/
void greasebranch(struct treenode *node_p,FILE *file_p,unsigned int e){
#include "variables.h"
  int steps,a,b,c;
  int slippe,to,from,to_root,to_false_root;
  int sp_case;

  double f_min,f_max,f,l;

  extern int root;
  struct treenode *movingnode,*node_c;
  struct treenode *tree;

  /*  Make two copies of the tree*/
  tree2=treecopy(node_p,0);
  tree=treecopy(node_p,0);
  filltree(node_p,tree,0);

  /*  Print out the tree branches and respective numbers*/
  printf("\n");
  for(a=0;a<branches;a++)
    printf("Branch %d goes from %s to %s\n",a
           ,(branch[a]->node[0])->name,branch[a]->name);

  /*  Get the necessary parameters from the user*/
  do{
    printf("\nWhich branch would you like to move?");
    scanf("%d",&slippe);
  }while(slippe<0 || slippe>branches);

  do{
    printf("\nWhich branch would you like to move branch %d along(towards)?"
                                                                    ,slippe);
    scanf("%d",&to);
    printf("\nWhich branch would you like to move branch %d along(away)?"
                                                                 ,slippe);
    scanf("%d",&from);
  }while(to<0 || to>branches || from<0 || from>branches);

/*  Out of the three parameters given, two of them must point
 * the common node.*/
  movingnode=(branch[to]->node[0]==branch[from]->node[0])
             ?branch[to]->node[0]:branch[slippe]->node[0];

/*  Find which branch in the common node points to each of the
 * parameter branches.*/
  from=findnode(movingnode,branch[from]);
  to=findnode(movingnode,branch[to]);
  slippe=findnode(movingnode,branch[slippe]);

/*  Get the amount of variation from the user*/
  printf("Sliding branch leading to %s\n\n",branch[slippe]->name);
  printf("Please enter min. & max. distance to move and number of points:\n");
  do{
    printf("\tmin. distance:");
    scanf("%lf",&f_min);
  }while(f_min<=-movingnode->length[from] || f_min>=movingnode->length[to]);
  do{
    printf("\tmax. distance:");
    scanf("%lf",&f_max);
  }while(f_min>=f_max || f_max>=movingnode->length[to]);
  do{
    printf("\tnumber of points:");
    scanf("%d",&steps);
  }while(steps<=0);

/*  Are we moving the branch towards the (code) root?*/
  to_false_root=(to==0)?0:1;

/*  Work out where root node is relative to movingnode.
 * Needed to work out whether to increase or decrease branch
 * lengths when moving node*/
  to_root=to_false_root; /*  Correct unless we are between the
                          * two roots in the tree */
  /*  If we have a genetic root...*/
  if(ISMODE(ROOTED)){
    /*  Decide depending on whether the root is floating or not,
     * which node it refers to*/
    node_c=(ISMODE(NODEASROOT))?branch[root]->node[0]:branch[root];

    /*  Iterate from the genetic root to the code root and set
     * "to_root" as appropriate. to_root only differs from to_false_root
     * if we are moving a branch between the two roots*/
    while(node_c!=tree){

      /*  If we reach movingnode without passing either the to or from
       * node first, then we must have come down on of the other branches
       * or our root is the movingnode itself. In either case, we can't 
       * keep evolutionary times fixed - problem*/
      if(node_c==movingnode){
        printf("Moving branch containing root - can't keep evolutionary times fixed\n");
        exit(0);
      }

      /*  We encounter the "to" node first and so must be moving 
       * towards the genetic root*/
      if(node_c==movingnode->node[to]){
        to_root=0;
        break;
      }

      /*  We encounter the from root first and so must be moving
       * away from the genetic node*/
      if(node_c==movingnode->node[from]){
        to_root=1;
        break;
      }

      /*  Iterate up one step*/
      node_c=node_c->node[0];
    }

  /*  If the root is on the branch we are moving, bad things
   * happen (TM) Moving the branch can't be done if evolutionary
   * times are fixed.
   *  We have deal with the case when root is on the moving bit of
   * tree but not the false root.*/
    if(slippe==0){
      /*  Decide where the root is depending on whether it is fixed
       * or floating*/
      node_c=(ISMODE(NODEASROOT))?branch[root]->node[0]:branch[root];
      a=0;

      /*  We may assume that the code root is on the moving branch of
       * the tree. If we have to pass through movingnode to get from
       * the genetic root to the code root then we are fine*/
      while(node_c!=tree){
	if(node_c==movingnode){
          a=1;
          break;
        }
        node_c=node_c->node[0];
      }

      /*  One special case if the movingnode is the fixed root*/
      if(ISMODE(NODEASROOT) && branch[root]->node[0]==movingnode)
        a=0;

      /*  If both the code and genetic roots are on the moving section
       * of tree, complain*/
      if(a==0){
        printf("Moving branch containing root - can't keep evolutionary"
               " times fixed\n");
        exit(0);
      }
    }

/*   If the root of the tree is on the same branch as we are moving across
 * then we'll have trouble when we cross it - give warning*/
    if(NOTMODE(NODEASROOT)) /*  genetic root is floating*/
      if(to_root==0) /*  we are heading towards it*/

        /*  Depending on whether we are going away or towards the
         * code root, see if the genetic root is on the branch we are
         * moving along*/
        if((to_false_root==0 && branch[root]==movingnode)
         ||(to_false_root==1 && branch[root]->node[0]==movingnode))
          fprintf(file_p,"# *** Root in branch moving across. Beware, here be dragons ***\n");
  }

  /*  Dump some information to file informing user what we are about
   * to do*/
  fprintf(file_p,"#Results from sliding branch from %s to"
	  " %s along branch from %s to %s\n"
	  ,movingnode->name,(movingnode->node[slippe])->name
	  ,(movingnode->node[from])->name,(movingnode->node[to])->name);
  fprintf(file_p,"#Distance moved towards %s\t%s\n",(movingnode->node[to])->name,outstring);
  if(ISMODE(MATRICES))
    fprintf(matrix_file_p,"#Matrices from sliding branch from %s to"
	    " %s along branch from %s to %s\n"
	    ,movingnode->name,(movingnode->node[slippe])->name
	    ,(movingnode->node[from])->name,(movingnode->node[to])->name);
  if(sample_file_p!=NULL){
    fprintf(sample_file_p,"#Points %d\n",steps);
    fprintf(sample_file_p,"#Samples from sliding branch from %s to"
	    " %s along branch from %s to %s\n"
	    ,movingnode->name,(movingnode->node[slippe])->name
	    ,(movingnode->node[from])->name,(movingnode->node[to])->name);
  }
  if(ISMODE(PROBS))
    fprintf(prob_file_p,"#Probabilities from sliding branch from %s to"
	    " %s along branch from %s to %s\n"
	    ,movingnode->name,(movingnode->node[slippe])->name
	    ,(movingnode->node[from])->name,(movingnode->node[to])->name);
  
  /*  l is the total length of the "branch" we are moving across*/
  l=movingnode->length[to]+movingnode->length[from];
  
  if(steps==1)
    f=f_min;
  else
    f=(f_max-f_min)/(steps-1);

  /* We aren't changing any rates, so we don't need to apply the
   * factor to any information (so factor=1)*/
  factor=1.0;
  factor_flag=0;
  
  /* Initialise the tree...*/
  movingnode->length[to]-=f_min;
  if(to==0){
    b=findnode(movingnode->node[0],movingnode);
    (movingnode->node[to])->length[b]=movingnode->length[to];
  }
  else
    (movingnode->node[to])->length[0]=movingnode->length[to];

  movingnode->length[from]+=f_min;
  if(from==0){
    b=findnode(movingnode->node[0],movingnode);
    (movingnode->node[from])->length[b]=movingnode->length[from];
  }
  else
    (movingnode->node[from])->length[0]=movingnode->length[from];
  /*  If the tree is rooted, then we must preserve evolutionary 
   * times. This requires adding bits onto any other branches
   * connected to the moving node, bar the to and from branches.*/
  if(ISMODE(ROOTED)){
    /*  Case moving towards the evolutionary root - all branches
     * connected to the moving node, bar from and to need to have 
     * the distance added onto them*/
    sp_case=(movingnode==tree)?1:0;
    if(to_root==0){
      if(sp_case==1 && to!=0 && from!=0){
	b=findnode(movingnode->node[0],movingnode);
	movingnode->length[0]+=f_min;
	(movingnode->node[0])->length[b]=movingnode->length[0];
      }
      for(a=sp_case;a<DOODAH && movingnode->node[a]!=NULL;a++)
	if(a!=to && a!=from){
	  (movingnode->node[a])->length[0]+=f_min;
	  movingnode->length[a]=(movingnode->node[a])->length[0];
        }
    }
    /*  Else we are moving away from it and all branches connected to
     * the moving node need the distance subtracted from them. Same
     * special case as above*/
    else{
      if(sp_case==1 && to!=0 && from!=0){
	b=findnode(movingnode->node[0],movingnode);
	movingnode->length[0]-=f_min;
	(movingnode->node[0])->length[b]=movingnode->length[0];
      }
      for(a=sp_case;a<DOODAH && movingnode->node[a]!=NULL;a++)
	if(a!=to && a!=from){
	(movingnode->node[a])->length[0]-=f_min;
	movingnode->length[a]=(movingnode->node[a])->length[0];
      }
    }
  }  

  
  for(c=0;c<steps;c++){

/*  Adding new lengths - if any of the markers are zero then
 * we have to search for the correct pointer in the parent node*/

    movingnode->length[to]-=(c!=0)*f;
    /*  Also need to update length of branch that leads to this
     * one. If to=0 then we must look up the tree*/
    if(to==0){
      b=findnode(movingnode->node[0],movingnode);
      (movingnode->node[to])->length[b]=movingnode->length[to];
    }
    else
      (movingnode->node[to])->length[0]=movingnode->length[to];

    movingnode->length[from]+=(c!=0)*f;
    /*  Also need to update length of branch that leads to this                 
     * one. If from=0 then we must look up the tree*/
    if(from==0){
      b=findnode(movingnode->node[0],movingnode);
      (movingnode->node[from])->length[b]=movingnode->length[from];
    }
    else
      (movingnode->node[from])->length[0]=movingnode->length[from];

    /*  If the tree is rooted, then we must preserve evolutionary 
     * times. This requires adding bits onto any other branches
     * connected to the moving node, bar the to and from branches.*/
    if(ISMODE(ROOTED)){
      /*  Case moving towards the evolutionary root - all branches
       * connected to the moving node, bar from and to need to have 
       * the distance added onto them*/
      sp_case=(movingnode==tree)?1:0; /*  Is the movingnode the code root?*/

      if(to_root==0){
        if(sp_case==1 && to!=0 && from!=0){
          b=findnode(movingnode->node[0],movingnode);
          movingnode->length[0]+=(c!=0)*f;
          (movingnode->node[0])->length[b]=movingnode->length[0];
        }
	for(a=sp_case;a<DOODAH && movingnode->node[a]!=NULL;a++)
	  if(a!=to && a!=from){
 	    (movingnode->node[a])->length[0]+=(c!=0)*f;
	    movingnode->length[a]=(movingnode->node[a])->length[0];
	  }
      }
      /*  Else we are moving away from it and all branches connected to
       * the moving node need the distance subtracted from them. Same
       * special case as above*/
      else{
        if(sp_case==1 && to!=0 && from!=0){
          b=findnode(movingnode->node[0],movingnode);
          movingnode->length[0]-=(c!=0)*f;
          (movingnode->node[0])->length[b]=movingnode->length[0];
        }
	for(a=sp_case;a<DOODAH && movingnode->node[a]!=NULL;a++)
	  if(a!=to && a!=from){
	    (movingnode->node[a])->length[0]-=(c!=0)*f;
	    movingnode->length[a]=(movingnode->node[a])->length[0];
	  }
      }
    }

    /*  Fill in second copy of tree from first and update the leaf
     * array to point to tree*/
    filltree(tree,tree2,0);
    leaves=0;
    doleaf(tree,0);

    /*  Print a bit of information to all open files*/
    printf("Doing distance %E\n",c*f+f_min);
    if(ISMODE(TREES))
      print_tree(tree,file_p,0);
    fprintf(file_p,"%E\n",c*f+f_min);
    if(ISMODE(MATRICES)){
      fprintf(matrix_file_p,"\n#Distance %E\n",c*f+f_min);
      if(ISMODE(TREES))
        print_tree(tree,matrix_file_p,0);
    }
    if(ISMODE(VARIANCE)){
      fprintf(variance_file_p,"\n#Distance %E\n",c*f+f_min);
      if(ISMODE(TREES))
        print_tree(tree,variance_file_p,0);
    }
    if(sample_file_p!=NULL){
      fprintf(sample_file_p,"\n#Distance %E\n",c*f+f_min);
      if(ISMODE(TREES))
	print_tree(tree,sample_file_p,0);
    }
    if(ISMODE(PROBS)){
      fprintf(prob_file_p,"\n#Distance %E\n",c*f+f_min);
      if(ISMODE(TREES))
	print_tree(tree,prob_file_p,0);
    }
    
    det=find_information(tree,tree2,e,factor_flag,factor);

    /*  Possible may want to print out more than one result*/
    b=1;
    if(ISMODE(INDIVIDUAL))
      b=individual;
    
    if(ISMODE(INDIVIDUAL)){
      if(ISMODE(DETINDIV)){
	fprintf(file_p,"\t\t\t%E\tD(",det[0]);
	for(a=0;a<b-1;a++)
	  fprintf(file_p,"%d,",interesting_branches[a]);
	fprintf(file_p,"%d)\n",interesting_branches[a]);
      }else
	for(a=0;a<b;a++)
	  fprintf(file_p,"\t\t\t%E\t%d\n",det[a],interesting_branches[a]);
    }else
      fprintf(file_p,"\t\t\t%E\tD\n",det[0]);
    
    if(ISMODE(HKY) && NOTMODE(NOKAPPA) && NOTMODE(PERCENTILE))
      fprintf(file_p,"\t\t\t%E\tKappa\n",det[b]);
    free(det);
  }
}


double *find_information(struct treenode *tree,struct treenode *tree2,unsigned int e, int factor_flag, double factor){

  extern double (*(*var)[])[];
  extern double (*(*expect)[])[];
  extern double (*(*rootedexpect)[])[];
  extern int mode;
  extern int nodecount;
  extern int branches;
  extern FILE *variance_file_p;
  extern int is_kappa;
 
  int a,b,d;
  double *det;

  CreatePMats ();
  
  /*  Work out the information the user requires
   * If sampling then sample to estimate the nth percentile*/
  if(ISMODE(PERCENTILE)){
    det=calloc(1,sizeof(double));
    det[0]=sample_percentile(tree,tree2,e,factor_flag,factor);
    if(ISMODE(CACHE))
      wipe_cache();
  }
  /*  Else calculate the expected information matrix*/
  else{ 
    expectation(e,factor,factor_flag,tree,tree2);
    if(ISMODE(CACHE))
      wipe_cache();
    
    /*  If the variances are wanted, then calculate and dump
       * the results to a file*/
    if(ISMODE(VARIANCE)){
      if(ISMODE(HKY) && NOTMODE(NOKAPPA))
        is_kappa=1;
      d=branches+is_kappa;
      /*  If we are working with rooted trees, the variance matrix
       * will already be in rooted form and so of a different size.
       * We must use the rooted expectation tree to create the variance
       * from E(X^2)*/
      if(ISMODE(ROOTED)){
        d=nodecount+is_kappa+((ISMODE(NODEASROOT))?1:2);
        planttree(expect,rootedexpect);
        for(a=0;a<d;a++)
          for(b=0;b<d;b++)
            (*(*var)[a])[b]-=(*(*rootedexpect)[a])[b]*(*(*rootedexpect)[a])[b];
      }
      /*  Same, but not rooted.*/
      else{
        for(a=0;a<d;a++)
          for(b=0;b<d;b++)
            (*(*var)[a])[b]-=(*(*expect)[a])[b]*(*(*expect)[a])[b];
      }

      /*  Dump the calculated variances to a file*/
      for(a=0;a<d;a++){
	for(b=0;b<d;b++)
	  fprintf(variance_file_p,"%e,",(*(*var)[a])[b]);
	fprintf(variance_file_p,"\n");
      }
      fprintf(variance_file_p,"\n");

      is_kappa=0;
    }

    /*  Work out the information required (either returns determinant
     * of expected information / estimated expected information or
     * an individual element, depending on what was required -
     * Note: recalculated rooted-expect again if doing variance calc.
     * Slightly wasteful*/
    det=determinant();
  }
  return det;
}
