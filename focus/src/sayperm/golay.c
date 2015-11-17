/*#define _GOLAYMAIN_*/

void GenerateGolay(int *golay)
 {
 
 /* Generates the so called Golay code[24,12] from a generator matrix, 
    z.B bricogne, Efficient Sampling Methods for Combination of Signs,
    Phases, Hyperphases, and molecular Orientations 
    or http://web.usna.navy.mil/~wdj/codes0.htm
 */

  int count;
  int i[12], j, k;  
  
  
const int gm[] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                  0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 
                  0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 
                  0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 
                  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 
                  0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 
                  0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 
                  0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 
                  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1};

/*


const int gm[] = {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 
                  1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 
                  1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 
                  1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 
                  1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 
                  1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 
                  1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 
                  1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 
                  1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 
                  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 
                  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
*/

/* These two generator matrices should be equivalent. They are from different sources */

  count = 0;

  for(i[0]=0; i[0]<=1; i[0]++)
  for(i[1]=0; i[1]<=1; i[1]++)
  for(i[2]=0; i[2]<=1; i[2]++)
  for(i[3]=0; i[3]<=1; i[3]++)
  for(i[4]=0; i[4]<=1; i[4]++)
  for(i[5]=0; i[5]<=1; i[5]++)
  for(i[6]=0; i[6]<=1; i[6]++)
  for(i[7]=0; i[7]<=1; i[7]++)
  for(i[8]=0; i[8]<=1; i[8]++)
  for(i[9]=0; i[9]<=1; i[9]++)
  for(i[10]=0; i[10]<=1; i[10]++)
  for(i[11]=0; i[11]<=1; i[11]++)
     {
        for(j=0; j<24; j++)
           {
             golay[24*count+j] = 0;
   
           for(k=0; k<12; k++)
              {
                 golay[24*count+j] += gm[24*k+j] * i[k];
              }
            
	    
            if(golay[24*count+j]&1) golay[24*count+j] = 0;
            else  golay[24*count+j] = 1;
             
           }
        count++;
    }
  
}

void GolayToCentrPhases(int *golay, double *golayphase)
  {
   int i;
  
   for(i=0; i<98304; i++)
     {
      if(golay[i]==0) golayphase[i] = 0.0;
      if(golay[i]==1) golayphase[i] = 0.5;
     } 
  }

void GolayToAcentrPhases(int *golay, double *golayphase)
  {
   int i,j;

     for(i=0; i<4096; i++)
     for(j=0; j<12; j++)
        {
	 if ((golay[24*i+j]==0)&&(golay[24*i+j+11]==0)) {golayphase[12*i+j] = 0.125; continue;}
	 if ((golay[24*i+j]==0)&&(golay[24*i+j+11]==1)) {golayphase[12*i+j] = 0.625; continue;}
	 if ((golay[24*i+j]==1)&&(golay[24*i+j+11]==0)) {golayphase[12*i+j] = 0.375; continue;}
	 if ((golay[24*i+j]==1)&&(golay[24*i+j+11]==1)) {golayphase[12*i+j] = 0.875; continue;}
        }
    return;
   }
  

  
#ifdef _GOLAYMAIN_  

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#define PI 3.1415927

void main(void)
{

 int *golay, i, j;
 double *golayphase;
 golay = (int *)malloc(98304*sizeof(int));
 golayphase = (double *)malloc(98304*sizeof(double));
 GenerateGolay(golay);
 GolayToCentrPhases(golay, golayphase);
 
 for(i=0; i<4096; i++) {
   fprintf(stdout,"%5d", i);
   for(j=0; j<24; j++) 
       fprintf(stdout,"%5.1f",golayphase[12*i+j]); 
       fprintf(stdout,"\n");}

    
 return;
}

#endif /*GOLAYMAIN*/
