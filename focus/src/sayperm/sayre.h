#include "function.h"

double  MakeSayreAcentr(double escale, double volume, double *SfacSquare, double *Sfac,
		  int Reflexstat1, int Reflexstat2, int Reflexstat3,
		  int Reflexstat4);
void SayrePhaseExtAcentr(double escale, double volume, double *SfacSquare, double *Sfac,
                    int Reflexstat1, int Reflexstat2, int Reflexstat3,
                    int Reflexstat4);
double  MakeSayreCycleAcentr(double escale, double volume, double *SfacSquare,
                      double *Sfac);
double  MakeSayreCentr(double escale, double volume, double *SfacSquare,double *Sfac, 
		      int Reflexstat1, int Reflexstat2, int Reflexstat3,
		      int Reflexstat4);
void SayrePhaseExtCentr(double escale, double volume, double *SfacSquare,
                        double *Sfac,
                        int Reflexstat1, int Reflexstat2, int Reflexstat3,
                        int Reflexstat4);
double  MakeSayreCycleCentr(double escale, double volume, double *SfacSquare,
                           double *Sfac);
void PhaseExtWoolfCodeCentr(double *Sfac, double *SfacSquare);		  
void PhaseExtGolayCodeCentr(double *Sfac, double *SfacSquare);
void PhaseExtGolayCodeAcentr(double *Sfac, double *SfacSquare);
void PhaseRatioAcentr(int *h_m_phase);
void GeneratePermSynth(int *grid_dim);
void GeneratePermSynthExt(int *grid_dim, double *Sfac, double *SfacSquare);
void fou_focus(int *grid_dim);
double  MakeSayreCentrTable(double escale, double volume, double *SfacSquare, double *Sfac,
		      int Reflexstat1, int Reflexstat2, int Reflexstat3,
		      int Reflexstat4);
