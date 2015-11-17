#ifndef _FUNCTION
#define _FUNCTION

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "sginfo.h"
#include "atominfodecl.h"



#define MAX_LINE_SIZE  255
#define TWO_PI         6.2831853
#define PI 3.1415927
#define SQRT_2 1.4142136
#define  RAD  0.01745329251994


typedef struct
  {
    double a;
    double b;
    double c;
    double alpha;
    double beta;
    double gamma;
  }
T_CellParam;



extern T_Eq_hkl *Eq_hkl;
extern T_SgInfo SgInfo;
extern T_CellParam *CellParam;
extern T_ChemicalComp *Comp;


extern double *ehkl, *fhkl;
extern double *norm_ehkl;
extern int *stat;
extern double *phase0;
extern int NumDepRefl, NumIndepRefl;
extern double *MaxRvalue;
extern int *h, *k, *l;
extern int countperm;
extern int num_plus, num_minus;
extern double escale, volume;
extern double B, Scale;

int ReadFouInp (char *space_group, T_CellParam * CellParam, int *grid_dim,
		int dev, int *limit,
		double *iso_value, double *escale,
		double *Scale, double *B,
		int *perm_limit, int *atom_number, double *f0,
		double *weight, double *PAWidth, char *code_name);
int ReadHklInp (void);
int ReadPhaseInp (double *perm_phase);
void ReadPhaseCombInp (double *phase_comb);
int CountPhaseInpLine (void);
void fourier (int *grid_dim, double *roxyz, int begin, int end);
void OutputFourier (double *roxyz, T_CellParam * CellParam, int *grid_dim, int map_num);
int CountSymmDepRefl (void);
void MakeNormalization (double *SayreE, double *ObsE, int ClustMem);
double CalculateH (int h, int k, int l);
void PointToAtomsOneDim (double *SfacSquare, double *Sfac);
void PointToAtoms (double *SfacSquare, double *Sfac, double PAWidth);
void RemoveTempScattPart (double B, double Scale, int ENum, int AtKind);
void FillHklAlloc (void);
void PositivityProbab (int atom_number);
void WriteNewHklInp (T_CellParam * CellParam, double *Sfac);
double Peakiness (int *grid_dim, double *roxyz);
void PointToAtomsWoolfson (double *SfacP, double *SfacPSquare, double *SfacPCub, double weightP,
      double *SfacQ, double *SfacQSquare, double *SfacQCub, double weightQ);
void MakeWoolfson (double escale, double volume,
		   double *fp, double *fp2, double *fp3,
		   double *fq, double *fq2, double *fq3,
		   int Reflexstat1, int Reflexstat2,
		   int Reflexstat3, int Reflexstat4);
int IsQuartAQuart (int h1, int k1, int l1,
		   int h2, int k2, int l2,
		   int h3, int k3, int l3,
		   int h4, int k4, int l4);
const T_PSE *PSEInfo (char *AtomLabel);
double CalcVolume (T_CellParam * CellParam);
int CountHklInpLine (void);int
IsTriplATripl (int h1, int k1, int l1,
               int h2, int k2, int l2,
               int h3, int k3, int l3);
int BuildSgInfo (T_SgInfo * SgInfo, const char *SgName);


#endif
