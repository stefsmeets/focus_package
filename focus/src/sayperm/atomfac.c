#include "function.h"
#include "atominfodecl.h"
#include <math.h>

void DetScatFacVarWidth(double *SfacSquare, double *Sfac, double fac);
int ReadPseudoAtomShape(double *x, double *y);
double CalcScatteringFac(double *x, double *y, double s, int number);
void CalcGaussAtomShape(double *x, double *y, double *z, double fac);
int CalcGaussAtomShapeCount(double fac);

void DetScatFacVarWidth(double *SfacSquare, double *Sfac, double fac)
 {
  double *x, *y, *z;
  int points, i;
  double H;

  H = 0.0;
  i=0;
  points = CalcGaussAtomShapeCount(fac);
  x = malloc((points+1) * sizeof(double));
  y = malloc((points+1) * sizeof(double));
  z = malloc((points+1) * sizeof(double));
  
  CalcGaussAtomShape(x, y, z, fac);
  /* while(i<=points){fprintf(stdout,"%9.3f %9.3f %9.3f\n",x[i],y[i], z[i]); i++;} */
  i = 0;

  for (i=1; i<= NumIndepRefl; i++)
    {
     H = CalculateH(Eq_hkl[i].h[0], 
              Eq_hkl[i].k[0], Eq_hkl[i].l[0]);
  
     Sfac[i] =  ehkl[i] * CalcScatteringFac(x, y, H, points);
     SfacSquare[i] = ehkl[i] * CalcScatteringFac(x, z, H, points);
     fprintf(stdout,"%4d %4d %4d %9.5f %9.3f %9.3f %9.3f %9.3f\n", Eq_hkl[i].h[0],  Eq_hkl[i].k[0], Eq_hkl[i].l[0], ehkl[i], H, Sfac[i], SfacSquare[i],
CalcScatteringFac(x, z, H, points)); 
    }
  free(x);
  free(y);
  free(z);

  }
/*
int ReadPseudoAtomShape(double *x, double *y)
  {

     FILE *po;
     int i;
     char buff[500];
      

     if ((po = fopen ( "shape.inp","r")) == NULL) 
          {fprintf(stdout, "Sorry, can not open shape.inp\n"); exit (1);}

     i = 0;
     while( EOF != fgetline(po, buff))
     if ( buff[0] != '#')
       {
         sscanf(buff,"%lf %lf",&x[i], &y[i]);
         i++;
       }
   return(i);
   fclose(po);
   }

int CountAtomShapeInpLine(void)
{
  FILE  *po;
  int valid_lines;
  char buff[500];
 
  if ((po = fopen ( "shape.inp","r")) == NULL) 
      {fprintf(stdout, "Sorry, can not open shape.inp\n"); exit (1);}
      valid_lines = 0;
      while( EOF != fgetline(po, buff))
                        if ( buff[0] != '#') valid_lines++;
                     
  fclose(po);
  return (valid_lines-1);
 }

*/

double CalcScatteringFac(double *x, double *y, double s, int number)
  {
    double A;
    int i;
    A = 0.000;
    
       
    for(i=1; i<=number; i++)
      {
          if (x[i-1] > x[i]) {
           fprintf(stdout," %4d %9.3f %9.3f \n",i, x[i-1],x[i]);   
            fprintf(stdout,"Check the pseudoatom shape table\n");
            exit(1); }
          
         /* Summation over the product of delta x,cos(2Pi
            xi/d), and the yi, this Sum is cubed and multiplied 
            by 4/3 Pi to expand this fourier from one dimension
            to a three dim sphere, A is the Scattering factor
            for a particular d value (1/s)*/
         
           
          A += (x[i]-x[i-1]) * cos(2*PI*x[i]*s)*y[i];
         
        } 
     A =  4.1887902*A*A*A;
     
     return(A);

   }   


int CalcGaussAtomShapeCount(double fac)
  {
   double x[2];
   double y;
   int i;
    
   x[0] = 0.0000;
   y = 1.0000;
   i = 0;
   
   while(y > 0.0001)
     {
       i++;
       x[1] = x[0] + 0.02;
       y = exp(-fac * PI * x[1] * x[1]);
       x[0] = x[1];
     }
   return(i);
   }    

void CalcGaussAtomShape(double *x, double *y, double *z, double fac)
  {
   int i;
   x[0] = 0.0000;
   y[0] = 1.0000;
   z[0] = 1.0000;
   
   i = 0;
   while(y[i] > 0.0001)
     {
       i++;
       x[i] = x[i-1] + 0.02;
       y[i] = exp(-fac * PI * x[i] * x[i]);
       z[i] = y[i] * y[i];
        /*fprintf(stdout," %9.3f %9.3f %9.3f\n" , x[i], y[i], fac);*/    
     }
   return;
   }
