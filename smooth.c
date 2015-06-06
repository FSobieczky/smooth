/*
smooth.c  :  Perform a smoothing of a PGM - file in ASCII (non-raw!) format
                (Note: use pnmnoraw file1.PGM file2.PGM to transform into
                       ASCII format.)
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

char *strcpy(char*, char*);
char *strcat(char*, char*);

output_f(float *w, int B, int H, char *filename)
{
  int i, j, temp;
  FILE *fpo;

  if((fpo=fopen(filename,"w"))==NULL)
    {
      printf("Couldn't open %s for writing !\n", filename);
      exit(-1);
    }

  fprintf(fpo,"P2\n%d %d\n%d\n", B, H, 255);
  for(j = 0; j<H; j++)
    {
      for(i = 0; i<B; i++)
	{
	  temp = (int)w[j*B+i];
	  fprintf(fpo, "%d ", temp);
	}
      fprintf(fpo, "\n");
    }
  fclose(fpo);
}


main(int argc, char **argv)
{
  int B, H, N, NN, temp1, temp2, n, t, s, m, i, j, k;
  int ngreyvalues;
  FILE *fp;
  char line[100], outfile[120];
  int *f, *range, *newrange, *deg;
  float *u, *v, *w, W;
  int *nr, *nu, *nl, *nd;

  if(argc!=4 && argc!=5)
    {
      printf("Usage: smooth <file> n s [<outfile>]\n\n`file' is plain PGM (Magic number P2)\
\n\nWrites smoothed file in PGM ASCII format to `file.out' or `outfile', if specified.\
  \n\n n is an integer greater or equal to 1 = number of iterates \n\
     in the smooothing (e.g. 3)\n\n\
 s is the difference in greyvalue between two neighboring greyvalues \n\
   (range: 0 - 255) to be sufficiently different so that the edge connecting \n\
   the two vertices is removed from the lattice and no smoothing across it \n\
   takes place..\n ");
      exit(-1);
    }

  t = atoi(argv[2]);   /* Number of smoothing steps */
  s = atoi(argv[3]);   /* threshold in greyvalues, when to cut an edge */
  if(argc==5)          /* then last argument is output filename.  */
    {          
      strcpy(outfile, argv[4]);
    }
  else  /* Otherwise it is `file.out' where file is input-file name.  */
    {
      strcpy(outfile, argv[1]);  
      strcat(outfile,".out");
    }
  if((fp=fopen(argv[1],"r"))==NULL)
    {
      printf("Couldn't open %s for reading !\n", argv[1]);
      exit(-1);
    }

  fscanf(fp,"%s", line);
  if(strcmp(line, "P2")!=0)
    {
      printf("'%s' is not a proper pbm-file.\n\n", argv[1]);
      exit(-1);
      }  

  fscanf(fp,"%s", line);
  if(line[0]=='#')
    {
      printf("Remove comment line (2nd line) from '%s'.\n\n", argv[1]);
      exit(-1);
      }
  fclose(fp);    /* Initial checking of inputfile finished. */

  if((fp=fopen(argv[1],"r"))==NULL)
    {
      printf("Couldn't open %s for reading !\n", argv[1]);
      exit(-1);
    }  /* file is reopened, so that fp* jumps to beginning of it. */

  fscanf(fp,"%s", line);
  temp1 = fscanf(fp,"%d %d", &B, &H);
  temp2 = fscanf(fp,"%d", &ngreyvalues);
  if(ngreyvalues!=255)
    { 
      printf("File doesn't have 0-255 greyvalues! (B:%d, H:%d)\n\n", B, H);
      exit(-1);
    }
  if(temp1!=2 || temp2!=1)  /* in this case fscanf didn't read in 3 values successfully, s.a. */
    {
      printf("'%s' is not a proper file.\n\n", argv[1]);
      exit(-1);
    }         /* Second checking of inputfile finished. If no `exit(-1)' so far, file ok. */
  

  N = H*B;  /* Number of Pixels */

  f = (int *)malloc((sizeof(int)*N+1));  /* The input image - stays fixed */
  nr = (int *)malloc((sizeof(int)*N+1));  /* edge to right present? */
  nu = (int *)malloc((sizeof(int)*N+1));  /* edge up present? */
  nl = (int *)malloc((sizeof(int)*N+1));  /* edge to left present? */
  nd = (int *)malloc((sizeof(int)*N+1));  /* edge downward present? */
  u = (float *)malloc((sizeof(float)*N+1));  /* The image which is changed from step to step */
  v = (float *)malloc((sizeof(float)*N+1));  /* The image one time step back  -  the last image */
  w = (float *)malloc((sizeof(float)*N+1)); /* Eventually the value of the normalized trace */
  deg = (int *)malloc((sizeof(int)*N+1));  /* number of neighbors of each of the vertices */

  for(j = 0; j<H; j++)
    {
      for(i = 0; i<B; i++)
	{
	  fscanf(fp, "%d", &(f[j*B+i]));
	}
    }
  fclose(fp);   
  /* Input finished. f[0..N] should have values between 0 and 255, now */

  
  for(i=0; i<N; i++)  /* Initializing deg[], and nr[],nu[],nl[],nd[] */
    {                 /* This is the abbreviation for neighborright, neighborup, 			 neighborleft, neighbordown */
      deg[i] = 0;
      nr[i] = 0;
      nu[i] = 0;
      nl[i] = 0;
      nd[i] = 0;
      if((i+1)%B!=0)  /* not on right edge */
	if(f[i]>=f[i+1]-s && f[i]<=f[i+1]+s)                      
	  {         
	         /* a large value of s here will make it */
	    /* more likely that smoothing occurs between f[i] and f[i+1] etc. */
	    deg[i]++;
	    nr[i] = 1;
	  }
      if(i%B !=0)    /* not on left edge */
	if(f[i]>=f[i-1]-s && f[i]<=f[i-1]+s)
	  {
	    deg[i]++;
	    nl[i]=1;
	  }
      if(i<N-B)     /* not on upper edge */
	if(f[i]>=f[i+B]-s && f[i]<=f[i+B]+s)
	  {
	    deg[i]++;
	    nu[i]=1;
	  }
      if(i>B-1)     /* not on lower edge */
	if(f[i]>=f[i-B]-s && f[i]<=f[i-B]+s)
	  {
	    deg[i]++;
	    nd[i]=1;
	  }
    } 
      /* Degree-vector deg[] set up  -  Note: maximum is 4; */
      /* Also nr[], nu[], nl[], nd[] also set up, now*/

      /*  vertex around which smoothing occurs: */
      /*     m = j*B + i  - goes through whole picture: */

  for(i=0; i<N; i++)
    {
      v[i]=f[i];  /* Picture in the beginning - no smoothing yet */
    }

  /* Main Loop: */
  
  for(n = 0; n<t; n++)
    {	  
      for(i=0; i<N; i++)
	{
	  u[i] = 0.0;
	  if(nr[i]==1)
	    u[i] += v[i+1];
	  if(nl[i]==1)
	    u[i] += v[i-1];
	  if(nu[i]==1)
	    u[i] += v[i+B];
	  if(nd[i]==1)
	    u[i] += v[i-B];
	  
	  u[i] += v[i]*(4-deg[i]);
	  
	  u[i] = 0.25*u[i];
	}
      
      for(i=0; i<N; i++)  /* Put new picture into old */
	v[i] = u[i];      /*  after every step */
    }

  output_f(u, B, H, outfile);   /* Output */
}
