/* Sfacsquare should be reseted in every main loop !! */



#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "sginfo.h"
#include "function.h"
#include "sayre.h"


T_Eq_hkl *Eq_hkl;
T_SgInfo SgInfo;
T_CellParam *CellParam;
T_ChemicalComp *Comp;


double *ehkl, *fhkl;
double *norm_ehkl;
int *stat;
double *phase0;
int NumDepRefl, NumIndepRefl;
double *MaxRvalue;
int *h, *k, *l;
int countperm;
int num_plus, num_minus;
double escale, volume;
double B, Scale;

int 
main (int argc, char *argv[])

{

  int i;

  int count;
  int CountAtomKind, ElectronNumber;
  char *space_group, *code_name;
  int set;
  double *Sfac, *SfacSquare, *SfacBuff;
  int NumReflPerm, atom_number;
  int count_exit, perm_reflections;
  int countprint;
  double f0;
  int *grid_dim;

  double R3, PAWidth;

  countperm = 0;
  count_exit = 0;
  set = 0;
  countprint = 0;
  count = 0;
  perm_reflections = 0;
  CountAtomKind = 0;
  ElectronNumber = 0;
  space_group = malloc (40 * sizeof (char));
  CellParam = malloc (sizeof (T_CellParam));

  Comp = malloc (sizeof (T_ChemicalComp));
  for (i = 0; i < 10; i++)
    Comp->atom[i] = malloc (20 * sizeof (char));

  MaxRvalue = malloc (7 * sizeof (double));
  grid_dim = malloc (7 * sizeof (int));
  code_name = malloc (40 * sizeof (char));

  strcpy (space_group, "P1");
  strcpy (code_name, "woolfson");
  escale = 1.0;
  Scale = 1.0;
  B = 0.0;
  volume = 1.0;
  PAWidth = 0.000;

  for (i = 0; i < argc; i++)
    fprintf (stdout, "argv[%3d]     %s\n", i, argv[i]);

  NumIndepRefl = CountHklInpLine ();
  Eq_hkl = malloc ((NumIndepRefl + 1) * sizeof (T_Eq_hkl));
  stat = malloc ((NumIndepRefl + 1) * sizeof (int));
  ehkl = malloc ((NumIndepRefl + 1) * sizeof (double));
  phase0 = malloc ((NumIndepRefl + 1) * sizeof (double));
  norm_ehkl = malloc ((NumIndepRefl + 1) * sizeof (double));
  Sfac = malloc ((NumIndepRefl + 1) * sizeof (double));
  SfacBuff = malloc ((NumIndepRefl + 1) * sizeof (double));
  SfacSquare = malloc ((NumIndepRefl + 1) * sizeof (double));
  fhkl = malloc ((NumIndepRefl + 1) * sizeof (double));




/* if (BuildSgInfo(&SgInfo, space_group) != 0)
   fprintf(stderr, "%s\n", SgError);
   else
   {
   ListSgInfo(&SgInfo, 1, 0, stdout);
   if (SgError)
   fprintf(stderr, "%s\n", SgError);
   }
   fprintf(stdout," acentric? %3d \n", SgInfo.Centric); */



  if (ReadFouInp (space_group, CellParam, grid_dim, 0, NULL, NULL, &escale, &Scale,
	       &B, NULL, &atom_number, &f0, NULL, &PAWidth, code_name) != 1)
    fprintf (stdout, "Error during reading of sp.inp");
  volume = CalcVolume (CellParam);

  fprintf (stdout, "\n\n\n");
  fprintf (stdout, "UnitCell         %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", CellParam->a,
	   CellParam->b, CellParam->c, CellParam->alpha, CellParam->beta,
	   CellParam->gamma);
  fprintf (stdout, "SpaceGroup       %s\n", space_group);
  fprintf (stdout, "MaxRvalue        %5.3f %5.3f %5.3f\n", MaxRvalue[1], MaxRvalue[2], MaxRvalue[3]);
  fprintf (stdout, "FScale           %5.3f\n", Scale);
  fprintf (stdout, "TempFac          %5.3f\n", B);
  fprintf (stdout, "PseudoAtom       %5.3f \n", PAWidth);
  fprintf (stdout, "PermCode         %s\n", code_name);

  fprintf (stdout, "volume (calc):    %5.3f\n", volume);
  fprintf (stdout, "\n\n\n");

  for (i = 1; i <= NumIndepRefl; i++)
    fhkl[i] = ehkl[i];

  NumDepRefl = CountSymmDepRefl ();
  fprintf (stdout, "NumIndepRefl :%4d\n", NumIndepRefl);
  fprintf (stdout, "NumDepRefl :%4d\n", NumDepRefl);
  fprintf (stdout, "%6.2f  %6.2f  %6.2f %6.2f %6.2f %6.2f %4d\n", CellParam->a,
	   CellParam->b, CellParam->c, CellParam->alpha, CellParam->beta,
	   CellParam->gamma, NumDepRefl);


  h = (int *) malloc ((NumDepRefl + 1) * sizeof (int));
  k = (int *) malloc ((NumDepRefl + 1) * sizeof (int));
  l = (int *) malloc ((NumDepRefl + 1) * sizeof (int));


/*  FillHklAlloc();
   RemoveTempScattPart(B, Scale, atom_number, f0);  
   fac = 0.55;
   DetScatFacVarWidth(SfacSquare, Sfac, fac);
 */

  FillHklAlloc ();
  RemoveTempScattPart (B, Scale, ElectronNumber, CountAtomKind);

  PointToAtoms (SfacSquare, Sfac, PAWidth);
  for (i = 1; i <= NumIndepRefl; i++)
    SfacBuff[i] = Sfac[i];
  NumReflPerm = 0;


  for (i = 1; i <= NumIndepRefl; i++)
    {
      SfacBuff[i] = Sfac[i];
      if (stat[i] == 1)
	perm_reflections++;
    }


  NumReflPerm = 0;

  if ((strcmp (code_name, "hamming") == 0) && (SgInfo.Centric == 0))
    fprintf (stdout, "Combination of Hamming Code and Acentric SG impossible, sorry!\n ");
  else if ((strcmp (code_name, "hamming") == 0) && (SgInfo.Centric == 1))
    {
      if (perm_reflections != 14)
	{
	  fprintf (stdout, "Can permute forteen reflections, only\n");
	  exit (0);
	}
      PhaseExtWoolfCodeCentr (Sfac, SfacSquare);
    }
  else if ((strcmp (code_name, "golay") == 0) && (SgInfo.Centric == 0))
    {
      if (perm_reflections != 12)
	{
	  fprintf (stdout, "Can permute twelve reflections, only\n");
	  exit (0);
	}
      PhaseExtGolayCodeAcentr (Sfac, SfacSquare);
    }
  else if ((strcmp (code_name, "golay") == 0) && (SgInfo.Centric == 1))
    {
      if (perm_reflections != 24)
	{
	  fprintf (stdout, "Can permute 24 reflections, only\n");
	  exit (0);
	}
      PhaseExtGolayCodeCentr (Sfac, SfacSquare);
    }
  else if (strcmp (code_name, "permsynth") == 0)
    {
      if (perm_reflections != 7)
	{
	  fprintf (stdout, "Can permute 7 reflections, only\n");
	  exit (0);
	}
      GeneratePermSynth (grid_dim);
    }
  else if (strcmp (code_name, "permsynthext") == 0)
    {
      if (perm_reflections != 7)
	{
	  fprintf (stdout, "Can permute 7 reflections, only\n");
	  exit (0);
	}
      GeneratePermSynthExt (grid_dim, Sfac, SfacSquare);
    }
  else if (strcmp (code_name, "fourier") == 0)
    {
      fou_focus (grid_dim);
    }
  else if (strcmp (code_name, "test") == 0)
    {
      R3 = MakeSayreCentrTable (escale, volume, SfacSquare, Sfac, 22, 22, 22, 22);
      fprintf (stdout, "R-value : %4.3f\n ", R3);
    }

  return (0);

}
