#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "sayre.h"
#include "golay.h"
#include "function.h"

double  MakeSayreCentr(double escale, double volume, double *SfacSquare,
                      double *Sfac,
		      int Reflexstat1, int Reflexstat2, int Reflexstat3,
		      int Reflexstat4)
{


  
 int count, count0, value, count_trip;
 int i, j, n;
 double E_Sum, E_SumA,SumDeltaE;
 double SumDeltaObs;
 double *imm_sumA;
 double *Ahkl, *PhaseSayre; 
 double *PhaseIntern;
 int *stat_intern;

  
 double SumDeltaEBefore;
 int convergence;
 int ClustMem;
 double *ObsE;
 double *SayreE;
 double E_Sum_buff; 
 
 ObsE = malloc(2 * sizeof(double)); 
 SayreE = malloc(2 * sizeof(double)); 
 PhaseIntern = malloc((NumIndepRefl+1) * sizeof(double));
 imm_sumA = malloc((NumDepRefl+1) * (NumDepRefl+1) * sizeof(double));
 Ahkl = malloc((NumDepRefl+1) * sizeof(double));
 PhaseSayre = malloc((NumDepRefl+1) * sizeof(double));
 stat_intern = malloc((NumDepRefl+1) * sizeof(int)); 
 
 count = 0;
 count0 = 0;
 count_trip=0;
 SumDeltaE = 0;
 SumDeltaObs = 1.0;
 E_SumA  = 0.0;
 convergence = 0;
 ClustMem = 0;
  
  
  
  for (i=1; i<= NumIndepRefl; i++) PhaseIntern[i] = phase0[i];
    
 
  for (i=1; i<= NumIndepRefl; i++)
      for (j=0; j<Eq_hkl[i].N; j++)
        {
         count++;
         Ahkl[count] = Sfac[i] * cos((PhaseIntern[i]+Eq_hkl[i].TH[j]/STBF)*TWO_PI);
	 stat_intern[count] = stat[i];
        }
    count = 0;     
   





  for (i=1; i<=NumIndepRefl; i++)
    if ((stat[i] != Reflexstat1) && (stat[i] != Reflexstat2) && 
            (stat[i] != Reflexstat3) && (stat[i] != Reflexstat4))
    {   
    
    for (j=1; j<=NumDepRefl; j++)
     if ((stat_intern[j] != Reflexstat1) && (stat_intern[j] != Reflexstat2) && 
               (stat_intern[j] != Reflexstat3)&& (stat[i] != Reflexstat4))
      {	       
    for (n=1; n<=NumDepRefl; n++)
       if ((stat_intern[n] != Reflexstat1) && (stat_intern[n] != Reflexstat2) &&
                      (stat_intern[n] != Reflexstat3)&& (stat[i] != Reflexstat4))
       
          {
             
            if   ( IsTriplATripl( Eq_hkl[i].h[0], Eq_hkl[i].k[0], Eq_hkl[i].l[0],
                                  h[j], k[j], l[j],
                                  h[n], k[n], l[n]) == 1)

	   
	       {
	       
                count0++;
                imm_sumA[count0] = Ahkl[j]*Ahkl[n];

                E_SumA = E_SumA + imm_sumA[count0];
                count_trip++;
                 
                value = 0;  
                        
                }
                
             }   
 
         }
     
     
     
     if (count_trip != 0) 
       {
           
        
         E_Sum  = E_SumA;
	
         if ((Eq_hkl[i].h[0]==0 )&& (Eq_hkl[i].k[0]==0) && (Eq_hkl[i].l[0]==0))
           {  
	     E_Sum_buff = (E_Sum*2 - Ahkl[1]*Ahkl[1])/(escale * volume);
             escale = escale * (E_Sum_buff/SfacSquare[i]);
	     E_Sum = (E_Sum*2 - Ahkl[1]*Ahkl[1])/(escale * volume);
             
	   }  
	     
	 
	 else
              E_Sum  = E_Sum / (escale * volume); 
	 
        }

     
       /* force to be centrosymmetric*/
     
       if (E_SumA < 0) PhaseSayre[i] = 0.5;
       else  PhaseSayre[i] = 0.0; 
       
     /*fprintf(stdout," E(%2d %4d %4d) %4d %9.3f %9.3f %9.3f %9.3f %9.3f %4d %4d\n",  
               Eq_hkl[i].h[0],  Eq_hkl[i].k[0], Eq_hkl[i].l[0], count_trip,
                Sfac[i], SfacSquare[i], E_Sum, phase0[i], PhaseSayre[i], stat[i],
                i);*/      
     

     SumDeltaE = SumDeltaE + fabs(fabs(SfacSquare[i]) - fabs(E_Sum));
     SumDeltaObs = SumDeltaObs + fabs(SfacSquare[i]);


     E_Sum  = 0.0;
     E_SumA  = 0.0;
     count0 = 0; 
     count_trip = 0;   
    
    }
    
    
     SumDeltaEBefore = SumDeltaE;
     SumDeltaE = SumDeltaE/SumDeltaObs;         

  

free (stat_intern);
free (imm_sumA);
free (Ahkl);
free (ObsE);
free (SayreE);
free (PhaseIntern);
free (PhaseSayre);
return SumDeltaE;     

}
void SayrePhaseExtCentr(double escale, double volume, double *SfacSquare, double *Sfac,
                        int Reflexstat1, int Reflexstat2, int Reflexstat3,
                        int Reflexstat4)
                 
{


  
  int count, count0, value, count_trip, *stat_intern;
  int i, j, n;
  double E_Sum, E_SumA;
  double *imm_sumA;
  double *Ahkl, *PhaseSayre; 


 double *SayreE;
 SayreE = malloc(100 * sizeof(double));
 imm_sumA = malloc((NumDepRefl+1) * (NumDepRefl+1) * sizeof(double));
 Ahkl = malloc((NumDepRefl+1) * sizeof(double));
 PhaseSayre = malloc((NumDepRefl+1) * sizeof(double));
 stat_intern = malloc((NumDepRefl+1) * sizeof(int));

  count = 0;
  count0 = 0;
  count_trip=0;
  E_SumA  = 0.0;

  
 
 for (i=1; i<= NumIndepRefl; i++)
      for (j=0; j<Eq_hkl[i].N; j++)
        {
         count++;
         if ((stat_intern[j] == Reflexstat1) || (stat_intern[j] == Reflexstat2) ||
            (stat_intern[j] == Reflexstat3) || (stat_intern[j] == Reflexstat4))
            {
              Ahkl[count] = Sfac[i] * cos((phase0[i]+Eq_hkl[i].TH[j]/STBF)*TWO_PI);
	      stat_intern[count] = stat[i];
            }
        }
     count = 0;     
  


/*
   Search for triplets of independent not marked reflections
   only in use of the marked reflections
*/



  for (i=1; i<=NumIndepRefl; i++)
    if ((stat[i] != Reflexstat1) && (stat[i] != Reflexstat2) &&
        (stat[i] != Reflexstat3) && (stat[i] != Reflexstat4))
    {
       
    for (j=1; j<=NumDepRefl; j++)
      if ((stat_intern[j] == Reflexstat1) || (stat_intern[j] == Reflexstat2) ||
          (stat_intern[j] == Reflexstat3) || (stat_intern[j] == Reflexstat4))
      {
       
      for (n=1; n<=NumDepRefl; n++)
        if ((stat_intern[n] == Reflexstat1) || (stat_intern[n] == Reflexstat2) ||
	 (stat_intern[n] == Reflexstat3) || (stat_intern[n] == Reflexstat4))
         {
         if   (( IsTriplATripl( Eq_hkl[i].h[0], Eq_hkl[i].k[0], Eq_hkl[i].l[0],
                                h[j], k[j], l[j],
                                h[n], k[n], l[n]) == 1))

	  {
          count0++;
          imm_sumA[count0] = Ahkl[j]*Ahkl[n];

          E_SumA = E_SumA + imm_sumA[count0];
          count_trip++;
                 
          value = 0;  
                        
           }
                
         }
             
       }      
 
    
    
     if (count_trip != 0) 
       {
           
        
         E_Sum  = E_SumA;
	
         if ((Eq_hkl[i].h[0]==0 )&& (Eq_hkl[i].k[0]==0) && (Eq_hkl[i].l[0]==0))
             E_Sum = (E_Sum*2 - Ahkl[1]*Ahkl[1])/(escale * volume);
	 
	 else
              E_Sum  = E_Sum / (escale * volume); 
	 
        }
  

        if (E_SumA < 0) {PhaseSayre[i] = 0.5; E_Sum = -E_Sum;}
	else PhaseSayre[i] = 0.0;
        
     
     /*  fprintf(stdout," E(%2d %4d %4d) %4d %9.3f %9.3f %9.3f %9.3f %9.3f %4d %4d\n" ,  
               Eq_hkl[i].h[0],  Eq_hkl[i].k[0], Eq_hkl[i].l[0], count_trip,
                Sfac[i], SfacSquare[i],E_Sum, phase0[i], PhaseSayre[i], stat[i],
                i); */                  
     
     E_Sum  = 0.0;
     E_SumA  = 0.0;
     count0 = 0; 
     count_trip = 0;   
    
    }
     
     
    for(i=1;i<=NumIndepRefl;i++)               
       if ((stat[i] != Reflexstat1) && (stat[i] != Reflexstat2) &&
        (stat[i] != Reflexstat3) && (stat[i] != Reflexstat4 ))
            
             phase0[i] = PhaseSayre[i];  

free (PhaseSayre);
free (imm_sumA);
free (Ahkl);
free (SayreE);
free (stat_intern);
return;     

}


double  MakeSayreCycleCentr(double escale, double volume, double *SfacSquare,
                           double *Sfac)
		      
{


  
  int count, count0, value, count_trip, count_cycl;
  int i, j, n;
  double E_Sum, E_SumA, SumDeltaE;
  double SumDeltaObs;
  double *imm_sumA;
  double *Ahkl, *PhaseSayre; 
  double *PhaseIntern;
  int *stat_intern;

  
 double SumDeltaEBefore;
 int convergence;
 int ClustMem;
 double *ObsE;
 double *SayreE;
 double E_Sum_buff;
 double *SfacSquareBuff; 
 
 ObsE = malloc((NumIndepRefl+1) * sizeof(double)); 
 SayreE = malloc((NumIndepRefl+1) * sizeof(double)); 
 PhaseIntern = malloc((NumIndepRefl+1) * sizeof(double));
 imm_sumA = malloc((NumDepRefl+1) * (NumDepRefl+1) * sizeof(double));
 Ahkl = malloc((NumDepRefl+1) * sizeof(double));
 PhaseSayre = malloc((NumDepRefl+1) * sizeof(double));
 stat_intern = malloc((NumDepRefl+1) * sizeof(int));
 SfacSquareBuff = malloc((NumIndepRefl+1) * sizeof(double)); 
 
  
count = 0;
  count0 = 0;
  count_trip=0;
  SumDeltaE = 0;
  SumDeltaObs = 1.0;
  E_SumA  = 0.0;
  convergence = 0;
  ClustMem = 0;
  count_cycl=0;
  
  
  for (i=1; i<= NumIndepRefl; i++)
      {
        PhaseIntern[i] = phase0[i];
        SfacSquareBuff[i] = SfacSquare[i];
        
      }
 
  for (i=1; i<= NumIndepRefl; i++)
      for (j=0; j<Eq_hkl[i].N; j++)
        {
         count++;
         Ahkl[count] = Sfac[i] * cos((PhaseIntern[i]+Eq_hkl[i].TH[j]/STBF)*TWO_PI);
	 stat_intern[count] = stat[i];
        }
         
   

do 
   {

  SumDeltaE = 0.00;
  SumDeltaObs = 0.00;
  count = 0;
  count0 = 0;
  
   
  for (i=1; i<=NumIndepRefl; i++)
    
    {   
    
    for (j=1; j<=NumDepRefl; j++)
     
      {	       
    for (n=1; n<=NumDepRefl; n++)
       
       
          {
             
            if   ( IsTriplATripl( Eq_hkl[i].h[0], Eq_hkl[i].k[0], Eq_hkl[i].l[0],
                                  h[j], k[j], l[j],
                                  h[n], k[n], l[n]) == 1)

	   
	       {
	       
                count0++;
                imm_sumA[count0] = Ahkl[j]*Ahkl[n];

                E_SumA = E_SumA + imm_sumA[count0];
                count_trip++;
                 
                value = 0;  
                        
                }
                
             }   
 
         }
     
     
     
     if (count_trip != 0) 
       {
           
        
         E_Sum  = E_SumA;
	
         if ((Eq_hkl[i].h[0]==0 )&& (Eq_hkl[i].k[0]==0) && (Eq_hkl[i].l[0]==0))
           {  
	     E_Sum_buff = (E_Sum*2 - Ahkl[1]*Ahkl[1])/(escale * volume);
	     escale = escale*(E_Sum_buff/SfacSquareBuff[i]);
	     E_Sum = (E_Sum*2 - Ahkl[1]*Ahkl[1])/(escale * volume);
	   }  
	     
	 
	 else
              E_Sum  = E_Sum / (escale * volume); 
	 
        }

     if (E_SumA<0.00) PhaseSayre[i] = 0.5;
     else PhaseSayre[i] = 0.0;

    
    /* fprintf(stdout," E(%2d %4d %4d) %4d %9.3f %9.3f %9.3f %9.3f %9.3f %4d %4d\n",  
               Eq_hkl[i].h[0],  Eq_hkl[i].k[0], Eq_hkl[i].l[0], count_trip,
                Sfac[i], SfacSquareBuff[i],E_Sum, phase0[i], PhaseSayre[i], stat[i],
                i); */   
     
     
     if (stat[i] > 4) 
      { 
         
         ClustMem++; 
         ObsE[ClustMem] = SfacSquare[i];
         PhaseIntern[i] = PhaseSayre[i];
         SayreE[ClustMem] = E_Sum;
   
         
        if  (stat[i+1] != stat[i] ) 
          {
            MakeNormalization(SayreE, ObsE, ClustMem);

            for ( n = 1; n <= ClustMem; n++)
                  {
                   
                   SfacSquareBuff[i-ClustMem+n] = ObsE[n];
                  }
                   
            ClustMem = 0; 
          }
       }                  

     SumDeltaE = SumDeltaE + fabs(fabs(SfacSquare[i]) - fabs(E_Sum));
     SumDeltaObs = SumDeltaObs + fabs(SfacSquare[i]);

     E_Sum  = 0.0;
     E_SumA  = 0.0;
     count0 = 0; 
     count_trip = 0;   
    
    }
    
    
    for(i=1;i<=NumIndepRefl;i++)
        {                  
          if (stat[i] != 0)    phase0[i] = PhaseSayre[i];
               
        }
        
     
     
     
     
     
     SumDeltaE = SumDeltaE/SumDeltaObs;
     fprintf(stdout, "%4.2f SumDeltaE = %4.2f  %4d\n",SumDeltaEBefore, SumDeltaE,
     convergence);
     if ( fabs(SumDeltaEBefore - SumDeltaE) < 0.005) convergence++;
     SumDeltaEBefore = SumDeltaE;

     /* AtomsToPoints(SfacSquareBuff, Sfac); */
     
     for (i=1; i<= NumIndepRefl; i++)
       for (j=0; j<Eq_hkl[i].N; j++)
        {
         count++; 
         Ahkl[count] = Sfac[i] * cos((phase0[i]+Eq_hkl[i].TH[j]/STBF)*TWO_PI);
        } 
     count_cycl++; if (count_cycl>10) {convergence=3; count_cycl=0;}
  }
  
  while (convergence < 3 );     
              

  
free (SfacSquareBuff);
free (stat_intern);
free (imm_sumA);
free (Ahkl);
free (ObsE);
free (SayreE);
free (PhaseIntern);
free (PhaseSayre);
return SumDeltaE;     

}

double  MakeSayreAcentr(double escale, double volume, double *SfacSquare,
                  double *Sfac,
		  int Reflexstat1, int Reflexstat2, int Reflexstat3,
		  int Reflexstat4)
{


  
  int count, count0, value, count_trip;
  int i, j, n;
  double E_Sum, E_SumA, E_SumB, SumDeltaE;
  double SumDeltaObs;
  double *imm_sumA, *imm_sumB ;
  double *Ahkl, *Bhkl, *PhaseSayre; 
  double *PhaseIntern;
  int *stat_intern;

  
 double SumDeltaEBefore;
 int convergence;
 int ClustMem;
 double *ObsE;
 double *SayreE;
 double E_Sum_buff; 
 
 ObsE = malloc(2 * sizeof(double)); 
 SayreE = malloc(2 * sizeof(double)); 
 PhaseIntern = malloc((NumIndepRefl+1) * sizeof(double));
 imm_sumA = malloc((NumDepRefl+1) * (NumDepRefl+1) * sizeof(double));
 imm_sumB = malloc((NumDepRefl+1) * (NumDepRefl+1) * sizeof(double));
 Ahkl = malloc((NumDepRefl+1) * sizeof(double));
 Bhkl = malloc((NumDepRefl+1) * sizeof(double));
 PhaseSayre = malloc((NumDepRefl+1) * sizeof(double));
 stat_intern = malloc((NumDepRefl+1) * sizeof(int)); 
 
  count = 0;
  count0 = 0;
  count_trip=0;
  SumDeltaE = 0;
  SumDeltaObs = 1.0;
  E_SumA  = 0.0;
  E_SumB  = 0.0;
  convergence = 0;
  ClustMem = 0;
  
  
  
  for (i=1; i<= NumIndepRefl; i++) PhaseIntern[i] = phase0[i];
    
 
  for (i=1; i<= NumIndepRefl; i++)
      for (j=0; j<Eq_hkl[i].N; j++)
        {
         count++;
         Ahkl[count] = Sfac[i] * cos((PhaseIntern[i]+Eq_hkl[i].TH[j]/STBF)*TWO_PI);
         Bhkl[count] = Sfac[i] * sin((PhaseIntern[i]+Eq_hkl[i].TH[j]/STBF)*TWO_PI);         
	 stat_intern[count] = stat[i];
        }
    count = 0;     
   





  for (i=1; i<=NumIndepRefl; i++)
    if ((stat[i] != Reflexstat1) && (stat[i] != Reflexstat2) && 
            (stat[i] != Reflexstat3) && (stat[i] != Reflexstat4))
    {   
    
    for (j=1; j<=NumDepRefl; j++)
     if ((stat_intern[j] != Reflexstat1) && (stat_intern[j] != Reflexstat2) && 
               (stat_intern[j] != Reflexstat3)&& (stat[i] != Reflexstat4))
      {	       
    for (n=1; n<=NumDepRefl; n++)
       if ((stat_intern[n] != Reflexstat1) && (stat_intern[n] != Reflexstat2) &&
                      (stat_intern[n] != Reflexstat3)&& (stat[i] != Reflexstat4))
       
          {
             
            if   ( IsTriplATripl( Eq_hkl[i].h[0], Eq_hkl[i].k[0], Eq_hkl[i].l[0],
                                  h[j], k[j], l[j],
                                  h[n], k[n], l[n]) == 1)

	   
	       {
	       
                count0++;
                imm_sumA[count0] = Ahkl[j]*Ahkl[n] - Bhkl[j]*Bhkl[n];
                imm_sumB[count0] = Ahkl[j]*Bhkl[n] + Bhkl[j]*Ahkl[n];

                E_SumA = E_SumA + imm_sumA[count0];
                E_SumB = E_SumB + imm_sumB[count0];
                count_trip++;
                 
                value = 0;  
                        
                }
                
             }   
 
         }
     
     
     
     if (count_trip != 0) 
       {
           
        
         E_Sum  = sqrt(E_SumA*E_SumA + E_SumB*E_SumB);
	
         if ((Eq_hkl[i].h[0]==0 )&& (Eq_hkl[i].k[0]==0) && (Eq_hkl[i].l[0]==0))
           {  
	     E_Sum_buff = (E_Sum*2 - Ahkl[1]*Ahkl[1])/(escale * volume);
             escale = escale * (E_Sum_buff/SfacSquare[i]);
	     E_Sum = (E_Sum*2 - Ahkl[1]*Ahkl[1])/(escale * volume);
             
	   }  
	     
	 
	 else
              E_Sum  = E_Sum / (escale * volume); 
	 
        }

      if ((E_SumA > 0.000001)||(E_SumA < -0.000001))            
              {
               PhaseSayre[i] = atan2(E_SumB,E_SumA) / TWO_PI;
               if (E_SumB < 0) PhaseSayre[i] = PhaseSayre[i] + 1; 
              }  
      else PhaseSayre[i] = 0.25;
       
      /* fprintf(stdout," E(%2d %4d %4d) %4d %9.3f %9.3f %9.5f %9.5f %9.3f %4d %4d\n",  
               Eq_hkl[i].h[0],  Eq_hkl[i].k[0], Eq_hkl[i].l[0], count_trip,
                Sfac[i], SfacSquare[i]*1000, E_Sum*1000, phase0[i], PhaseSayre[i], stat[i],
                i); */     
     
     if(i>1)
         {
           SumDeltaE = SumDeltaE + fabs(fabs(SfacSquare[i]) - fabs(E_Sum));
           SumDeltaObs = SumDeltaObs + fabs(SfacSquare[i]);
         }
 
     E_Sum  = 0.0;
     E_SumA  = 0.0;
     E_SumB  = 0.0; 
     count0 = 0; 
     count_trip = 0;   
    
    }
    
    
     SumDeltaEBefore = SumDeltaE;
     SumDeltaE = SumDeltaE/SumDeltaObs;         

  

free (stat_intern);
free (imm_sumA);
free (imm_sumB);
free (Ahkl);
free (Bhkl);
free (ObsE);
free (SayreE);
free (PhaseIntern);
free (PhaseSayre);
return SumDeltaE;     

}

void SayrePhaseExtAcentr(double escale, double volume, double *SfacSquare, double *Sfac,
                    int Reflexstat1, int Reflexstat2, int Reflexstat3,
                    int Reflexstat4)
                 
{


  
  int count, count0, value, count_trip, *stat_intern;
  int i, j, n;
  double E_Sum, E_SumA, E_SumB;
  double *imm_sumA, *imm_sumB ;
  double *Ahkl, *Bhkl, *PhaseSayre; 


 double *SayreE;
 SayreE = malloc(100 * sizeof(double));
 imm_sumA = malloc((NumDepRefl+1) * (NumDepRefl+1) * sizeof(double));
 imm_sumB = malloc((NumDepRefl+1) * (NumDepRefl+1) * sizeof(double));
 Ahkl = malloc((NumDepRefl+1) * sizeof(double));
 Bhkl = malloc((NumDepRefl+1) * sizeof(double));
 PhaseSayre = malloc((NumDepRefl+1) * sizeof(double));
 stat_intern = malloc((NumDepRefl+1) * sizeof(int));

  count = 0;
  count0 = 0;
  count_trip=0;
  E_SumA  = 0.0;
  E_SumB  = 0.0;

  
 
 for (i=1; i<= NumIndepRefl; i++)
      for (j=0; j<Eq_hkl[i].N; j++)
        {
         count++;
         if ((stat_intern[j] == Reflexstat1) || (stat_intern[j] == Reflexstat2) ||
            (stat_intern[j] == Reflexstat3) || (stat_intern[j] == Reflexstat4))
            {
              Ahkl[count] = Sfac[i] * cos((phase0[i]+Eq_hkl[i].TH[j]/STBF)*TWO_PI);
              Bhkl[count] = Sfac[i] * sin((phase0[i]+Eq_hkl[i].TH[j]/STBF)*TWO_PI);
	      stat_intern[count] = stat[i];
            }
        }
     count = 0;     
  


/*
   Search for triplets of independent not marked reflections
   only in use of the marked reflections
*/



  for (i=1; i<=NumIndepRefl; i++)
    if ((stat[i] != Reflexstat1) && (stat[i] != Reflexstat2) &&
        (stat[i] != Reflexstat3) && (stat[i] != Reflexstat4))
    {
       
    for (j=1; j<=NumDepRefl; j++)
      if ((stat_intern[j] == Reflexstat1) || (stat_intern[j] == Reflexstat2) ||
          (stat_intern[j] == Reflexstat3) || (stat_intern[j] == Reflexstat4))
      {
      
      for (n=1; n<=NumDepRefl; n++)
        if ((stat_intern[n] == Reflexstat1) || (stat_intern[n] == Reflexstat2) ||
	 (stat_intern[n] == Reflexstat3) || (stat_intern[n] == Reflexstat4))
         {
           
         if   (( IsTriplATripl( Eq_hkl[i].h[0], Eq_hkl[i].k[0], Eq_hkl[i].l[0],
                                h[j], k[j], l[j],
                                h[n], k[n], l[n]) == 1))

	  {
	       
          count0++;
          imm_sumA[count0] = Ahkl[j]*Ahkl[n] - Bhkl[j]*Bhkl[n];
          imm_sumB[count0] = Ahkl[j]*Bhkl[n] + Bhkl[j]*Ahkl[n];

          E_SumA = E_SumA + imm_sumA[count0];
          E_SumB = E_SumB + imm_sumB[count0];
          count_trip++;
                 
          value = 0;  
                        
           }
                
         }
             
       }      
 
    
    
     if (count_trip != 0) 
       {
           
        
         E_Sum  = sqrt(E_SumA*E_SumA + E_SumB*E_SumB);
	
         if ((Eq_hkl[i].h[0]==0 )&& (Eq_hkl[i].k[0]==0) && (Eq_hkl[i].l[0]==0))
             E_Sum = (E_Sum*2 - Ahkl[1]*Ahkl[1])/(escale * volume);
	 
	 else
              E_Sum  = E_Sum / (escale * volume); 
	 
        }
  
        if ((E_SumA > 0.000001)||(E_SumA < -0.000001))            
              {
               PhaseSayre[i] = atan2(E_SumB,E_SumA) / TWO_PI;
               if (E_SumB < 0) PhaseSayre[i] = PhaseSayre[i] + 1; 
              }  
         else PhaseSayre[i] = 0.25;
        
     
     /* fprintf(stdout," E(%2d %4d %4d) %4d %9.3f %9.3f %9.3f %9.3f %9.3f %4d %4d\n" ,  
               Eq_hkl[i].h[0],  Eq_hkl[i].k[0], Eq_hkl[i].l[0], count_trip,
                Sfac[i], SfacSquare[i],E_Sum, phase0[i], PhaseSayre[i], stat[i],
                i); */                 
     
     E_Sum  = 0.0;
     E_SumA  = 0.0;
     E_SumB  = 0.0; 
     count0 = 0; 
     count_trip = 0;   
    
    }
     
     
    for(i=1;i<=NumIndepRefl;i++)               
       if ((stat[i] != Reflexstat1) && (stat[i] != Reflexstat2) &&
        (stat[i] != Reflexstat3) && (stat[i] != Reflexstat4 ))
            
             phase0[i] = PhaseSayre[i];  

free (PhaseSayre);
free (imm_sumA);
free (imm_sumB);
free (Ahkl);
free (Bhkl);
free (SayreE);
free (stat_intern);
return;     

}



double  MakeSayreCycleAcentr(double escale, double volume, double *SfacSquare,
                      double *Sfac)
		      
{


  
  int count, count0, value, count_trip, count_cycl;
  int i, j, n;
  double E_Sum, E_SumA, E_SumB, SumDeltaE;
  double SumDeltaObs;
  double *imm_sumA, *imm_sumB ;
  double *Ahkl, *Bhkl, *PhaseSayre; 
  double *PhaseIntern;
  int *stat_intern;

  
 double SumDeltaEBefore;
 int convergence;
 int ClustMem;
 double *ObsE;
 double *SayreE;
 double E_Sum_buff;
 double *SfacSquareBuff; 
 
 ObsE = malloc((NumIndepRefl+1) * sizeof(double)); 
 SayreE = malloc((NumIndepRefl+1) * sizeof(double)); 
 PhaseIntern = malloc((NumIndepRefl+1) * sizeof(double));
 imm_sumA = malloc((NumDepRefl+1) * (NumDepRefl+1) * sizeof(double));
 imm_sumB = malloc((NumDepRefl+1) * (NumDepRefl+1) * sizeof(double));
 Ahkl = malloc((NumDepRefl+1) * sizeof(double));
 Bhkl = malloc((NumDepRefl+1) * sizeof(double));
 PhaseSayre = malloc((NumDepRefl+1) * sizeof(double));
 stat_intern = malloc((NumDepRefl+1) * sizeof(int));
 SfacSquareBuff = malloc((NumIndepRefl+1) * sizeof(double)); 
 
  
count = 0;
  count0 = 0;
  count_trip=0;
  SumDeltaE = 0;
  SumDeltaObs = 1.0;
  E_SumA  = 0.0;
  E_SumB  = 0.0;
  convergence = 0;
  ClustMem = 0;
  count_cycl=0;
  
  
  for (i=1; i<= NumIndepRefl; i++)
      {
        PhaseIntern[i] = phase0[i];
        SfacSquareBuff[i] = SfacSquare[i];
        
      }
 
  for (i=1; i<= NumIndepRefl; i++)
      for (j=0; j<Eq_hkl[i].N; j++)
        {
         count++;
         Ahkl[count] = Sfac[i] * cos((PhaseIntern[i]+Eq_hkl[i].TH[j]/STBF)*TWO_PI);
         Bhkl[count] = Sfac[i] * sin((PhaseIntern[i]+Eq_hkl[i].TH[j]/STBF)*TWO_PI);         
	 stat_intern[count] = stat[i];
        }
         
   

do 
   {

  SumDeltaE = 0.00;
  SumDeltaObs = 0.00;
  count = 0;
  count0 = 0;
  
   
  for (i=1; i<=NumIndepRefl; i++)
    
    {   
    
    for (j=1; j<=NumDepRefl; j++)
     
      {	       
    for (n=1; n<=NumDepRefl; n++)
       
       
          {
             
            if   ( IsTriplATripl( Eq_hkl[i].h[0], Eq_hkl[i].k[0], Eq_hkl[i].l[0],
                                  h[j], k[j], l[j],
                                  h[n], k[n], l[n]) == 1)

	   
	       {
	       
                count0++;
                imm_sumA[count0] = Ahkl[j]*Ahkl[n] - Bhkl[j]*Bhkl[n];
                imm_sumB[count0] = Ahkl[j]*Bhkl[n] + Bhkl[j]*Ahkl[n];

                E_SumA = E_SumA + imm_sumA[count0];
                E_SumB = E_SumB + imm_sumB[count0];
                count_trip++;
                 
                value = 0;  
                        
                }
                
             }   
 
         }
     
     
     
     if (count_trip != 0) 
       {
           
        
         E_Sum  = sqrt(E_SumA*E_SumA + E_SumB*E_SumB);
	
         if ((Eq_hkl[i].h[0]==0 )&& (Eq_hkl[i].k[0]==0) && (Eq_hkl[i].l[0]==0))
           {  
	     E_Sum_buff = (E_Sum*2 - Ahkl[1]*Ahkl[1])/(escale * volume);
	     escale = escale*(E_Sum_buff/SfacSquareBuff[i]);
	     E_Sum = (E_Sum*2 - Ahkl[1]*Ahkl[1])/(escale * volume);
	   }  
	     
	 
	 else
              E_Sum  = E_Sum / (escale * volume); 
	 
        }

   
   
   
     if ((E_SumA > 0.000001)||(E_SumA < -0.000001))            
              {
               PhaseSayre[i] = atan2(E_SumB,E_SumA) / TWO_PI;
               if (E_SumB < 0) PhaseSayre[i] = PhaseSayre[i] + 1; 
              }  
       else PhaseSayre[i] = 0.25;

    
      if ((PhaseSayre[i] <= 0.0625) || (PhaseSayre[i] > 0.9375)) PhaseSayre[i] = 0.000; 
      if ((PhaseSayre[i] >= 0.0625) && (PhaseSayre[i] < 0.1875)) PhaseSayre[i] = 0.125;
      if ((PhaseSayre[i] >= 0.1875) && (PhaseSayre[i] < 0.3125)) PhaseSayre[i] = 0.250;
      if ((PhaseSayre[i] >= 0.3125) && (PhaseSayre[i] < 0.4375)) PhaseSayre[i] = 0.375;
      if ((PhaseSayre[i] >= 0.4375) && (PhaseSayre[i] < 0.5625)) PhaseSayre[i] = 0.500;
      if ((PhaseSayre[i] >= 0.5625) && (PhaseSayre[i] < 0.6875)) PhaseSayre[i] = 0.625;
      if ((PhaseSayre[i] >= 0.6875) && (PhaseSayre[i] < 0.8125)) PhaseSayre[i] = 0.750;
      if ((PhaseSayre[i] >= 0.8125) && (PhaseSayre[i] < 0.9375)) PhaseSayre[i] = 0.875; 


    
    /* fprintf(stdout," E(%2d %4d %4d) %4d %9.3f %9.3f %9.3f %9.3f %9.3f %4d %4d\n",  
               Eq_hkl[i].h[0],  Eq_hkl[i].k[0], Eq_hkl[i].l[0], count_trip,
                Sfac[i], SfacSquareBuff[i],E_Sum, phase0[i], PhaseSayre[i], stat[i],
                i); */   
     
     
     if (stat[i] > 4) 
      { 
         
         ClustMem++; 
         ObsE[ClustMem] = SfacSquare[i];
         PhaseIntern[i] = PhaseSayre[i];
         SayreE[ClustMem] = E_Sum;
   
         
        if  (stat[i+1] != stat[i] ) 
          {
            MakeNormalization(SayreE, ObsE, ClustMem);

            for ( n = 1; n <= ClustMem; n++)
                  {
                   
                   SfacSquareBuff[i-ClustMem+n] = ObsE[n];
                  }
                   
            ClustMem = 0; 
          }
       }                  

     SumDeltaE = SumDeltaE + fabs(fabs(SfacSquare[i]) - fabs(E_Sum));
     SumDeltaObs = SumDeltaObs + fabs(SfacSquare[i]);

     E_Sum  = 0.0;
     E_SumA  = 0.0;
     E_SumB  = 0.0; 
     count0 = 0; 
     count_trip = 0;   
    
    }
    
    
    for(i=1;i<=NumIndepRefl;i++)
        {                  
          if (stat[i] != 0)    phase0[i] = PhaseSayre[i];
               
        }
        
     
     
     
     
     
     SumDeltaE = SumDeltaE/SumDeltaObs;
     fprintf(stdout, "%4.2f SumDeltaE = %4.2f  %4d\n",SumDeltaEBefore, SumDeltaE,
     convergence);
     if ( fabs(SumDeltaEBefore - SumDeltaE) < 0.005) convergence++;
     SumDeltaEBefore = SumDeltaE;

     
     for (i=1; i<= NumIndepRefl; i++)
       for (j=0; j<Eq_hkl[i].N; j++)
        {
         count++; 
         Ahkl[count] = Sfac[i] * cos((phase0[i]+Eq_hkl[i].TH[j]/STBF)*TWO_PI);
         Bhkl[count] = Sfac[i] * sin((phase0[i]+Eq_hkl[i].TH[j]/STBF)*TWO_PI);
        } 
     count_cycl++; if (count_cycl>10) {convergence=3; count_cycl=0;}
  }
  
  while (convergence < 3 );     
              

  
free (SfacSquareBuff);
free (stat_intern);
free (imm_sumA);
free (imm_sumB);
free (Ahkl);
free (Bhkl);
free (ObsE);
free (SayreE);
free (PhaseIntern);
free (PhaseSayre);
return SumDeltaE;     


}



void PhaseRatioAcentr(int *h_m_phase)
 {
   int i;
    
   for(i=0; i<8; i++) h_m_phase[i] = 0;
   for(i=1; i<=NumIndepRefl; i++)
    { 
      if ((phase0[i] <= 0.0625) || (phase0[i] > 0.9375)) h_m_phase[0]++; 
      if ((phase0[i] >= 0.0625) && (phase0[i] < 0.1875)) h_m_phase[1]++;
      if ((phase0[i] >= 0.1875) && (phase0[i] < 0.3125)) h_m_phase[2]++;
      if ((phase0[i] >= 0.3125) && (phase0[i] < 0.4375)) h_m_phase[3]++;
      if ((phase0[i] >= 0.4375) && (phase0[i] < 0.5625)) h_m_phase[4]++;
      if ((phase0[i] >= 0.5625) && (phase0[i] < 0.6875)) h_m_phase[5]++;
      if ((phase0[i] >= 0.6875) && (phase0[i] < 0.8125)) h_m_phase[6]++;
      if ((phase0[i] >= 0.8125) && (phase0[i] < 0.9375)) h_m_phase[7]++; 
    }
 return;
}


void PhaseExtWoolfCodeCentr(double *Sfac, double *SfacSquare)
  {
    int n1, n2, i, m, count;
    double R1, R2, R3, R4; 
    double *phase_comb, *perm_phase;
    int phase_number;
    


  phase_comb = malloc( 136 * sizeof(double));
  phase_number = CountPhaseInpLine();        
  perm_phase = malloc((phase_number+5)  * sizeof(double));
  phase_number = ReadPhaseInp(perm_phase);
  ReadPhaseCombInp(phase_comb);

  count = 0;

  for(n1=0; n1<16; n1++)
  for(n2=0; n2<16; n2++)
    {
      m=0;
      for(i=1; i<=NumIndepRefl; i++)
       if (stat[i] == 1) 
            {
              if (m<7) 
                 {phase0[i] = phase_comb[16*m+n1];fprintf(stdout,"%4d %4d %4.1f\n", m, n1, phase_comb[16*m+n1]); m++;}
              else
                 {phase0[i] = phase_comb[16*(m-7)+n2];fprintf(stdout,"%4d %4d %4.1f\n", m, n2, phase_comb[16*(m-7)+n2]); m++;}
            }            
       
    /* Transform of one dim array into k; l array 
       by num_columm*line_number+columm_number */   
                                                 
             count++;
       
            
             SayrePhaseExtCentr(escale, volume, SfacSquare, Sfac, 0, 1, 1, 1);
             R1 = MakeSayreCentr(escale, volume, SfacSquare, Sfac, 3, 4, 5, 6);
             if (R1 < MaxRvalue[1])
             
               {
	         SayrePhaseExtCentr(escale, volume, SfacSquare, Sfac, 0, 1, 2, 2 );
                 R2 = MakeSayreCentr(escale, volume, SfacSquare, Sfac, 4, 5, 6 ,7);
		 if (R2 < MaxRvalue[2])
		   {
		     
		     SayrePhaseExtCentr(escale, volume, SfacSquare, Sfac, 0, 1, 2, 3);
                     R3 = MakeSayreCentr(escale, volume, SfacSquare, Sfac, 22, 22, 22, 22);
		     if (R3 < MaxRvalue[3])
		       {
		       R4 = MakeSayreCycleCentr(escale, volume, SfacSquare, Sfac);
		       WriteNewHklInp(CellParam, Sfac);
		       fprintf(stdout,"#%-9d R1 %6.3f;1 R2 %6.3f;  R3 %6.3f; R4 %6.3f %6d %4d\n" 
		               , count, R1, R2, R3, R4, num_plus, num_minus);

		       fprintf(stdout,"\n\n\n");
		       fflush(stdout);
                       	       
		       }
		       
		   }  
		 
               }  
   
     
   } /* loop over phase combinations */
return;

}




void PhaseExtGolayCodeCentr(double *Sfac, double *SfacSquare)
  {
    int n1,i, m, count;
    double R1, R2, R3, R4; 
 
 
    int *golay;
    double *golayphase;
    fprintf(stdout,"in PhaseExtGolayCodeCentr1 \n" );
    golay = (int *)malloc(98304*sizeof(int));
    golayphase = (double *)malloc(98304*sizeof(double));
    fprintf(stdout,"in PhaseExtGolayCodeCentr2 \n" );
    GenerateGolay(golay);
    fprintf(stdout,"in PhaseExtGolayCodeCentr3 \n" );
    GolayToCentrPhases(golay, golayphase);
    fprintf(stdout,"in PhaseExtGolayCodeCentr4 \n" );

  count = 0;

  for(n1=0; n1<4096; n1++)
    {
      m=0;
      fprintf(stdout,"%5d ", n1);
      for(i=1; i<=NumIndepRefl; i++)
       if (stat[i] == 1) 
            {
                 {phase0[i] = golayphase[24*n1+m];fprintf(stdout," %4.1f", golayphase[24*n1+m]); m++;}
            }            
       fprintf(stdout,"\n"); 
    /* Transform of one dim array into k; l array 
       by num_columm*line_number+columm_number */   
                                                 
             count++;
       
            
             SayrePhaseExtCentr(escale, volume, SfacSquare, Sfac, 0, 1, 1, 1);
             R1 = MakeSayreCentr(escale, volume, SfacSquare, Sfac, 3, 4, 5, 6);
             if (R1 < MaxRvalue[1])
             
               {
	         SayrePhaseExtCentr(escale, volume, SfacSquare, Sfac, 0, 1, 2, 2 );
                 R2 = MakeSayreCentr(escale, volume, SfacSquare, Sfac, 4, 5, 6 ,7);
		 if (R2 < MaxRvalue[2])
		   {
		     
		     SayrePhaseExtCentr(escale, volume, SfacSquare, Sfac, 0, 1, 2, 3);
                     R3 = MakeSayreCentr(escale, volume, SfacSquare, Sfac, 22, 22, 22, 22);
		     if (R3 < MaxRvalue[3])
		       {
		       R4 = MakeSayreCycleCentr(escale, volume, SfacSquare, Sfac);
		       WriteNewHklInp(CellParam, Sfac);
		       fprintf(stdout,"#%-9d R1 %6.3f;1 R2 %6.3f;  R3 %6.3f; R4 %6.3f %6d %4d\n" 
		               , count, R1, R2, R3, R4, num_plus, num_minus);

		       fprintf(stdout,"\n\n\n");
		       fflush(stdout);
                       	       
		       }
		       
		   }  
		 
               }  
   
     
   } /* loop over phase combinations */
return;

}


void PhaseExtGolayCodeAcentr(double *Sfac, double *SfacSquare)
  {
    int n1,i, m, count;
    double R1, R2, R3, R4; 
 
 
    int *golay;
    double *golayphase;
    
    golay = (int *)malloc(98304*sizeof(int));
    golayphase = (double *)malloc(98304*sizeof(int));
    
    GenerateGolay(golay);
    GolayToAcentrPhases(golay, golayphase);


  count = 0;

  for(n1=0; n1<4096; n1++)
    {
      m=0;
      fprintf(stdout,"%5d ", n1);
      for(i=1; i<=NumIndepRefl; i++)
       if (stat[i] == 1) 
            {
                 {phase0[i] = golayphase[12*n1+m];fprintf(stdout," %5.2f", golayphase[12*n1+m]); m++;}
            }            
       fprintf(stdout,"\n"); 
    /* Transform of one dim array into k; l array 
       by num_columm*line_number+columm_number */   
                                                 
             count++;
       
            
             SayrePhaseExtAcentr(escale, volume, SfacSquare, Sfac, 0, 1, 1, 1);
             R1 = MakeSayreAcentr(escale, volume, SfacSquare, Sfac, 3, 4, 5, 6);
             if (R1 < MaxRvalue[1])
             
               {
	         SayrePhaseExtAcentr(escale, volume, SfacSquare, Sfac, 0, 1, 2, 2 );
                 R2 = MakeSayreAcentr(escale, volume, SfacSquare, Sfac, 4, 5, 6 ,7);
		 if (R2 < MaxRvalue[2])
		   {
		     
		     SayrePhaseExtAcentr(escale, volume, SfacSquare, Sfac, 0, 1, 2, 3);
                     R3 = MakeSayreAcentr(escale, volume, SfacSquare, Sfac, 22, 22, 22, 22);
		     if (R3 < MaxRvalue[3])
		       {
		       WriteNewHklInp(CellParam, Sfac);
		       R4 = MakeSayreCycleAcentr(escale, volume, SfacSquare, Sfac);
		       fprintf(stdout,"#%-9d R1 %6.3f;1 R2 %6.3f;  R3 %6.3f; R4 %6.3f %6d %4d\n" 
		               , count, R1, R2, R3, R4, num_plus, num_minus);

		       fprintf(stdout,"\n\n\n");
		       fflush(stdout);
                       	       
		       }
		       
		   }  
		 
               }  
   
     
   } /* loop over phase combinations */
return;

}


void GeneratePermSynth(int *grid_dim)
{

 int i,m,n, n1;
 double *perm_phase, *phase_comb;
 int phase_number;
 double *roxyz;


fprintf(stdout," %4d %4d %4d %4d %4d %4d\n ", 
                  grid_dim[1], 
                  grid_dim[2], grid_dim[3],grid_dim[4],
                  grid_dim[5], grid_dim[6]);


 NumIndepRefl = CountHklInpLine();

 phase_comb = malloc( 136 * sizeof(double));
 roxyz = malloc(grid_dim[1] * grid_dim[2] * grid_dim[3] *
                 grid_dim[4] * grid_dim[5] * grid_dim[6] * sizeof(double));
 
       
  phase_number = CountPhaseInpLine();        
  perm_phase = malloc((phase_number+5)  * sizeof(double));
  
  phase_number = ReadPhaseInp(perm_phase);
  ReadPhaseCombInp(phase_comb);
  fprintf(stdout,"%4d \n",phase_number); 

for(n=0; n<16; n++)
    {
     fprintf(stdout,"\n");   
     for(m=0; m<7; m++)
       fprintf(stdout,"%4.1f", phase_comb[16*m+n]);
    }
fprintf(stdout,"\n");

for(n1=0; n1<16; n1++)
    {
      m=0;
      for(i=1; i<=NumIndepRefl; i++)
       if (stat[i] == 1) 
            {
              if (m<7) 
                 {phase0[i] = phase_comb[16*m+n1];fprintf(stdout,"%4d %4d %4.1f\n", m, n1, phase_comb[16*m+n1]); m++;}
             
            }            
       
    		     WriteNewHklInp(CellParam, ehkl); 
                     fourier(grid_dim, roxyz, 0, 0);
                     OutputFourier(roxyz, CellParam, grid_dim, n1);
  
    } /* loop over phase combinations */
return;

}



void GeneratePermSynthExt(int *grid_dim, double *Sfac, double *SfacSquare)
{

 int i,m,n, n1;
 double *perm_phase, *phase_comb;
 int phase_number;
 double *roxyz;


fprintf(stdout," %4d %4d %4d %4d %4d %4d\n ", 
                  grid_dim[1], 
                  grid_dim[2], grid_dim[3],grid_dim[4],
                  grid_dim[5], grid_dim[6]);


 NumIndepRefl = CountHklInpLine();

 phase_comb = malloc( 136 * sizeof(double));
 roxyz = malloc(grid_dim[1] * grid_dim[2] * grid_dim[3] *
                 grid_dim[4] * grid_dim[5] * grid_dim[6] * sizeof(double));
 
       
  phase_number = CountPhaseInpLine();        
  perm_phase = malloc((phase_number+5)  * sizeof(double));
  
  phase_number = ReadPhaseInp(perm_phase);
  ReadPhaseCombInp(phase_comb);
  fprintf(stdout,"%4d \n",phase_number); 

for(n=0; n<16; n++)
    {
     fprintf(stdout,"\n");   
     for(m=0; m<7; m++)
       fprintf(stdout,"%4.1f", phase_comb[16*m+n]);
    }
fprintf(stdout,"\n");

for(n1=0; n1<16; n1++)
    {
      m=0;
      for(i=1; i<=NumIndepRefl; i++)
       if (stat[i] == 1) 
            {
              if (m<7) 
                 {phase0[i] = phase_comb[16*m+n1];fprintf(stdout,"%4d %4d %4.1f\n", m, n1, phase_comb[16*m+n1]); m++;}
             
            }            
                     SayrePhaseExtCentr(escale, volume, SfacSquare, Sfac, 0, 1, 1, 1); 
    		     WriteNewHklInp(CellParam, ehkl); 
                     fourier(grid_dim, roxyz, 0, 0);
                     OutputFourier(roxyz, CellParam, grid_dim, n1);
  
    } /* loop over phase combinations */
return;

}

void fou_focus(int *grid_dim)
  {
   double *roxyz;
   roxyz = malloc(grid_dim[1] * grid_dim[2] * grid_dim[3] *
                 grid_dim[4] * grid_dim[5] * grid_dim[6] * sizeof(double));
    fprintf(stdout, "output at sayperm 1446 \n");  
   fourier(grid_dim, roxyz, 0, 0);
    fprintf(stdout, "output at sayperm 1448 \n");  
   OutputFourier(roxyz, CellParam, grid_dim, 1);
    fprintf(stdout, "output at fourier 1450 \n");  
   return;
  } 





double  MakeSayreCentrTable(double escale, double volume, double *SfacSquare,
                      double *Sfac,
		      int Reflexstat1, int Reflexstat2, int Reflexstat3,
		      int Reflexstat4)
{


  
 int count, count0, value, count_trip;
 int i, j, n;
 double E_Sum, E_SumA,SumDeltaE;
 double SumDeltaObs;
 double *imm_sumA;
 double *Ahkl, *PhaseSayre; 
 double *PhaseIntern;
 int *stat_intern;

  
 double SumDeltaEBefore;
 int convergence;
 int ClustMem;
 double *ObsE;
 double *SayreE;
 double E_Sum_buff; 
 
 ObsE = malloc(2 * sizeof(double)); 
 SayreE = malloc(2 * sizeof(double)); 
 PhaseIntern = malloc((NumIndepRefl+1) * sizeof(double));
 imm_sumA = malloc((NumDepRefl+1) * (NumDepRefl+1) * sizeof(double));
 Ahkl = malloc((NumDepRefl+1) * sizeof(double));
 PhaseSayre = malloc((NumDepRefl+1) * sizeof(double));
 stat_intern = malloc((NumDepRefl+1) * sizeof(int)); 
 
 count = 0;
 count0 = 0;
 count_trip=0;
 SumDeltaE = 0;
 SumDeltaObs = 1.0;
 E_SumA  = 0.0;
 convergence = 0;
 ClustMem = 0;
  
  
  
  for (i=1; i<= NumIndepRefl; i++) PhaseIntern[i] = phase0[i];
    
 
  for (i=1; i<= NumIndepRefl; i++)
      for (j=0; j<Eq_hkl[i].N; j++)
        {
         count++;
         Ahkl[count] = Sfac[i] * cos((PhaseIntern[i]+Eq_hkl[i].TH[j]/STBF)*TWO_PI);
	 stat_intern[count] = stat[i];
        }
    count = 0;     
   





  for (i=1; i<=NumIndepRefl; i++)
    if ((stat[i] != Reflexstat1) && (stat[i] != Reflexstat2) && 
            (stat[i] != Reflexstat3) && (stat[i] != Reflexstat4))
    {   
    
    for (j=1; j<=NumDepRefl; j++)
     if ((stat_intern[j] != Reflexstat1) && (stat_intern[j] != Reflexstat2) && 
               (stat_intern[j] != Reflexstat3)&& (stat[i] != Reflexstat4))
      {	       
    for (n=1; n<=NumDepRefl; n++)
       if ((stat_intern[n] != Reflexstat1) && (stat_intern[n] != Reflexstat2) &&
                      (stat_intern[n] != Reflexstat3)&& (stat[i] != Reflexstat4))
       
          {
             
            if   ( IsTriplATripl( Eq_hkl[i].h[0], Eq_hkl[i].k[0], Eq_hkl[i].l[0],
                                  h[j], k[j], l[j],
                                  h[n], k[n], l[n]) == 1)

	   
	       {
	       
                count0++;
                imm_sumA[count0] = Ahkl[j]*Ahkl[n];

                E_SumA = E_SumA + imm_sumA[count0];
                count_trip++;
                 
                value = 0;  
                        
                }
                
             }   
 
         }
     
     
     
     if (count_trip != 0) 
       {
           
        
         E_Sum  = E_SumA;
	
         if ((Eq_hkl[i].h[0]==0 )&& (Eq_hkl[i].k[0]==0) && (Eq_hkl[i].l[0]==0))
           {  
	     E_Sum_buff = (E_Sum*2 - Ahkl[1]*Ahkl[1])/(escale * volume);
             escale = escale * (E_Sum_buff/SfacSquare[i]);
	     E_Sum = (E_Sum*2 - Ahkl[1]*Ahkl[1])/(escale * volume);
             
	   }  
	     
	 
	 else
              E_Sum  = E_Sum / (escale * volume); 
	 
        }

     
       /* force to be centrosymmetric*/
     
       if (E_SumA < 0) PhaseSayre[i] = 0.5;
       else  PhaseSayre[i] = 0.0; 
       
     fprintf(stdout," E(%2d %4d %4d) %4d %9.3f %9.3f %9.3f %9.3f %9.3f %4d %4d\n",  
               Eq_hkl[i].h[0],  Eq_hkl[i].k[0], Eq_hkl[i].l[0], count_trip,
                Sfac[i], SfacSquare[i], E_Sum, phase0[i], PhaseSayre[i], stat[i],
                i);      
     

     SumDeltaE = SumDeltaE + fabs(fabs(SfacSquare[i]) - fabs(E_Sum));
     SumDeltaObs = SumDeltaObs + fabs(SfacSquare[i]);


     E_Sum  = 0.0;
     E_SumA  = 0.0;
     count0 = 0; 
     count_trip = 0;   
    
    }
    
    
     SumDeltaEBefore = SumDeltaE;
     SumDeltaE = SumDeltaE/SumDeltaObs;         

  

free (stat_intern);
free (imm_sumA);
free (Ahkl);
free (ObsE);
free (SayreE);
free (PhaseIntern);
free (PhaseSayre);
return SumDeltaE;     

}
