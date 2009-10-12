/*  Program to calculate the nth percentile and mean from the output
 * of EDIBLE.*/
#include <stdio.h>
#include <stdlib.h>

struct percentile{
  double point;
  struct percentile *next;
};

int grt(double *one, double *two);

void main(int argc, char *argv[]){
  FILE *input_file_p;
  FILE *output_file_p;
  int sample_size,sequence_length,pt;
  double (*points)[];
  double (*results)[];
  struct percentile *percent, *temp;
  int max_percentiles,j,k,min_percentiles,m;
  double a;
  int (*wibble)();
  char ch;
  
  /*  Check number of arguments, exit if wrong*/
  if(argc!=3){
    printf("Wrong number of arguments!\n"
	   "Usage: pcalc sample_file output_file\n");
    exit(0);
  }
  
  /*  Open the necessary files for use. Die if they can't be.*/
  input_file_p=fopen(argv[1],"r");
  output_file_p=fopen(argv[2],"w");
  if(input_file_p==NULL){
    printf("Can't open input file.\n");
    exit(0);
  }
  if(output_file_p==NULL){
    printf("Can't open output file.\n");
    exit(0);
  }
 
  /*  Read in some values from file*/
  fscanf(input_file_p,"%*s %*s %*s %d ",&sample_size);
  fscanf(input_file_p,"%*s %*s %*s %*s %*s %*s %d\n",&sequence_length);
  fscanf(input_file_p,"%*s %d",&pt);
  if(sample_size<=0 || sequence_length<=0 || pt<=0){
    printf("Problems with file format.\n");
    exit(0);
  }
  
  /*  Read in the percentiles from the user*/
  printf("Which percentiles would you like calculated?\n");
  printf("Enter value between 0 and 100 followed by return\n");
  printf("Enter a value outside this range to start the calculation\n");
  max_percentiles=0;
  percent=malloc(sizeof(struct percentile));
  temp=percent;
  do{
    printf("?");
    scanf("%le",&a);
    if(a>0 && a<100){
      max_percentiles++;
      temp->next=malloc(sizeof(struct percentile));
      temp=temp->next;
      temp->point=a/100;
      temp->next=NULL;
    }
  }while(a>0 && a<100);

  printf("Going through data\n");
  temp=percent;
  percent=percent->next;
  free(temp);
  
  /*  Percentile now points to the first result wanted, 'i' is the 
   * number of percentile points wanted*/
  points=calloc(max_percentiles,sizeof(double));
  results=calloc(sample_size,sizeof(double));
  if(points==NULL || results==NULL){
    printf("Out of memory!\n");
    exit(2);
  }

  /*  Tranfer the results from the structure to an array
   * (easier to sort and use)*/
							  
  for(j=0;j<max_percentiles;j++){
    (*points)[j]=percent->point;
    temp=percent;
    percent=percent->next;
    free(temp);
  }

  /*  Sort the percentile points to be calculated*/
  wibble=&grt;
  qsort(points,max_percentiles,sizeof(double),wibble);

  /*  Check to see if any of the points can't be calculated, must be 
   * inbetween two of the sampled points to be used*/
  while((int)(sample_size*(*points)[max_percentiles-1]+0.5)>sample_size-1){
    printf("%e too large for interpolation.\n",(*points)[max_percentiles-1]*100);
    max_percentiles--;
  }
  min_percentiles=0;
  while((int)(sample_size*(*points)[min_percentiles]+0.5)<1){
    printf("%e too small for interpolation.\n",(*points)[min_percentiles]*100);
    min_percentiles++;
  }
  
  /*  Go to end of line*/
  while((ch=getc(input_file_p))!='\n'){}
  /*  Copy next line into the output file*/
  while((ch=getc(input_file_p))!='\n'){
    fprintf(output_file_p,"%c",ch);
  }
  fprintf(output_file_p,"\n");

  /*  Starting the main loop*/
  for(j=0;j<pt;j++){
    /*  Find the start of this data and scan in the factor*/
    ch=getc(input_file_p);
    ch=getc(input_file_p);
    if(ch=='#'){
      fprintf(output_file_p,"#");
      while((ch=getc(input_file_p))!=' ')
	fprintf(output_file_p,"%c",ch);
      fscanf(input_file_p,"%le",&a);
      fprintf(output_file_p," %e\n",a);
      ch=getc(input_file_p);
      ch=getc(input_file_p);
    }

    /*  Have we printed the tree in the file? If so, then 
     * copy it into the output*/
    if(ch=='#'){
      fprintf(output_file_p,"#");
      while((ch=getc(input_file_p))!='\n')
	fprintf(output_file_p,"%c",ch);
      fprintf(output_file_p,"\n");
    }
    else
      ch=ungetc(ch,input_file_p);

    /*  Scan in the sampled values into an array*/
    a=0;
    for(k=0;k<sample_size;k++){
      fscanf(input_file_p,"%le",&(*results)[k]);
      a+=(*results)[k];
    }
    /*  Print the mean and then interpolate all the percentile points*/
    fprintf(output_file_p,"#Sample mean=%e\n",(double)(a/sample_size));
    for(k=min_percentiles;k<max_percentiles;k++){
      m=(int)((*points)[k]*sample_size+1.5);
      a=(*results)[m-2]+((*points)[k]*sample_size+1.5-m)
	*((*results)[k-1]-(*results)[k-2]);
      fprintf(output_file_p,"%e\t%e\n",(*points)[k],a);
    } 
    fprintf(output_file_p,"\n");
    ch=getc(input_file_p);
  }
}

/*  Routine to compare two doubles and return +1 if the first is 
 * greater than the second, else -1*/
int grt(double *one, double *two){
  if((*one-*two)<0)
    return -1;
  else
    return 1;
}
