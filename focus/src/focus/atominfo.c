#include <stdio.h>
#include <ctype.h>
#include <math.h>


#define ATOMINFO_C__


#include "atominfo.h"

/*
http://www.ccp14.ac.uk/ccp/web-mirrors/objcryst/ObjCryst/atominfo.html

AtomInfo - Scattering factors etc. for ANSI C
Ralf W. Grosse-Kunstleve 

Compilation: 
cc -O -o atominfo -DSTAND_ALONE atominfo.c -lm
Usage: 
usage: atominfo [-IT92|WK95] [-ExactLabel] [-Anomalous] [-AnomalousHenke] AtomLabel
                [Start_sinTovL [End_sinTovL [StepSize]]]
References: 
IT92
 International Tables for Crystallography, Volume C 1992, Table 6.1.1.4 (pp. 500-502)
WK95
 Waasmaier & Kirfel; Acta Cryst. (1995), A51, 416-431.

ftp://wrzx02.rz.uni-wuerzburg.de/pub/local/Crystallography/sfac.dat 
 Element names and atomic weights
 CRC Handbook of Chemistry & Physics, 63rd edition, 1982-1983
 CRC Handbook of Chemistry & Physics, 70th edition, 1989-1990
 Atom radii
 Inorganic Crystal Structure Database (ICSD) User's Manual 1986
 Bond Scattering Lengths and cross-sections (Added V. Favre-Nicolin 14Sep2000) 
            Neutron News, Vol. 3, No. 3, 1992, pp. 29-37. 
           http://www.ncnr.nist.gov/resources/n-lengths/list.html 

Anomalous Scattering Factors: (taken from the DABAX database) 
 Henke Tables (Added V. Favre-Nicolin 14Sep2000):

B. L. Henke, E. M. Gullikson, and J. C. Davis, Atomic Data and Nuclear Data Tables Vol. 54 No. 2 (July 1993). 
ftp://grace.lbl.gov/pub/sf/ 
  
 Sasaki Tables (Added V. Favre-Nicolin 13Feb2001):

S.Sasaki (1989),KEK Report, 88-14, 1-136 
S.Sasaki (1990) , KEK Report, 90-16, 1-143 
ftp://pfweis.kek.jp/pub/Sasaki-table/ 
(these include long-range and detailed data around the absorption edges, all in only one data thanks to Bruce Ravel's reformatting).


For more detailed references open atominfo.h. 

Concise documentation: 
Missing. Sorry.
Example 1: 
% atominfo Si

14 "Si" "silicon" 28.086
Scattering Factor Label (IT92): "Si"
  IT92-CAA a1 b1 a2 b2 a3 b3 a4 b4 c
    6.2915 2.4386 3.0353 32.3337 1.9891 0.6785 1.541 81.6937 1.1407
Scattering Factor Label (WK95): "Si"
  WK95-CAA a1 b1 a2 b2 a3 b3 a4 b4 a5 b5 c
    5.27533 2.63134 3.19104 33.7307 1.51151 0.081119 1.35685 86.2886 2.51911 1.17009 0.145073
Atom Radius 1.34 "Si"
Example 2: 
% atominfo C 0 1 0.05

6 "C" "carbon" 12.011
Scattering Factor Label (IT92): "C"
  IT92-CAA a1 b1 a2 b2 a3 b3 a4 b4 c
    2.31 20.8439 1.02 10.2075 1.5886 0.5687 0.865 51.6512 0.2156
Scattering Factor Label (WK95): "C"
  WK95-CAA a1 b1 a2 b2 a3 b3 a4 b4 a5 b5 c
    2.65751 14.7808 1.07808 0.776775 1.49091 42.0868 -4.24107 -0.000294 0.713791 0.239535 4.29798
Atom Radius 0.86 "C"
   stol     IT92     WK95
  0.000    5.999    5.997
  0.050    5.749    5.749
  0.100    5.108    5.110
  0.150    4.310    4.310
  0.200    3.560    3.557
  0.250    2.950    2.949
  0.300    2.494    2.497
  0.350    2.171    2.173
  0.400    1.948    1.947
  0.450    1.794    1.791
  0.500    1.686    1.683
  0.550    1.604    1.603
  0.600    1.537    1.539
  0.650    1.479    1.483
  0.700    1.425    1.430
  0.750    1.373    1.377
  0.800    1.321    1.324
  0.850    1.270    1.271
  0.900    1.218    1.218
  0.950    1.167    1.166
  1.000    1.115    1.113
Source code: 
atominfo.c 
atominfo.h
*/

static void StripLabel(const char *Label, int Exact, char *buf, int mbuf)
{
  int         i;
  const char  *digit;


  mbuf--;

  while (*Label && isspace(*Label)) Label++;

  digit = NULL;

  for (i = 0; i < mbuf && *Label; Label++)
  {
    if (isspace(*Label))
      break;

    if (isdigit(*Label))
    {
      if (digit) break;
      digit = Label;
    }
    else if (*Label == '+' || *Label == '-')
    {
      if (i + 1 < mbuf)
      {
        if (digit == NULL)
            digit = "1";

        buf[i++] = *digit;
      }

      buf[i++] = *Label++;

      break;
    }
    else if (digit)
      break;
    else
      buf[i++] = toupper(*Label);
  }

  if (Exact && *Label && ! isspace(*Label))
    i = 0;

  buf[i] = '\0';
}


const T_PSE *FindInPSE(const char *Label, int Exact)
{
  int          i, m;
  char         buf[3];
  const T_PSE  *ePSE, *mPSE;


  StripLabel(Label, Exact, buf, sizeof buf / sizeof (*buf));

  if (buf[0] == '\0')
    return NULL;

  m    = 0;
  mPSE = NULL;

  for (ePSE = PSE; ePSE->Z; ePSE++)
  {
    for (i = 0; buf[i] && ePSE->Symbol[i]; i++)
      if (buf[i] != toupper(ePSE->Symbol[i])) break;

    if (buf[i] == '\0' && ePSE->Symbol[i] == '\0')
      return ePSE;

    if (i == 1 && isalpha(ePSE->Symbol[1]))
      i = 0;

    if (i > m)
    {
      m    = i;
      mPSE = ePSE;
    }
  }

  if (Exact)
    return NULL;

  return mPSE;
}

//////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// International tables, section C 2011, Table 4.3.2.2 //
//////////////////////////////////////////////////////////////////////////////////////////

const T_SF_ALL_CAA *FindSF_IT4322(const char *Label, int Exact)
// Takes string *label (label of scatterer) => returns SFC struct (SF_IT4322 @ atominfo.h)
// Only way to retrieve data from atominfo.h
{
  int                  i, m;
  char                 buf[6];
  const T_SF_ALL_CAA  *SFC, *mSFC;


  StripLabel(Label, Exact, buf, sizeof buf / sizeof (*buf));

  if (buf[0] == '\0')
    return 0;

  m    = 0;
  mSFC = NULL;

  for (SFC = SF_ALL_CAA; SFC->Label; SFC++)
  {
    for (i = 0; buf[i] && SFC->Label[i]; i++)
      if (buf[i] != toupper(SFC->Label[i])) break;

    if (buf[i] == '\0' && SFC->Label[i] == '\0')
    {
      //if (SFC->c) return SFC;

      for (i = 0; i < (sizeof SFC->e / sizeof (*SFC->e)); i++)
        if (SFC->e[i]) return SFC;

      for (i = 0; i < (sizeof SFC->d / sizeof (*SFC->d)); i++)
        if (SFC->d[i]) return SFC;

      return NULL;
    }

    if (i == 1 && isalpha(SFC->Label[1]))
      i = 0;

    if (i > m && ! isdigit(SFC->Label[i - 1]))
    {
      m    = i;
      mSFC = SFC;
    }
  }

  if (Exact)
    return NULL;

  return mSFC;
}

double CalcSF_IT4322(const T_SF_ALL_CAA *CAA, double stol2)
{
  int     i;
  double  sf;


  // sf = CAA->c;
  sf = 0.0;


  for (i = 0; i < (sizeof CAA->d / sizeof (*CAA->d)); i++)
    sf += CAA->d[i] * exp(-CAA->e[i] * stol2);

  return sf;
}

//////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// International tables, section C 2011, Table 4.3.2.3 //
//////////////////////////////////////////////////////////////////////////////////////////

const T_SF_ALL_CAA *FindSF_IT4323(const char *Label, int Exact)
{
  int                  i, m;
  char                 buf[6];
  const T_SF_ALL_CAA  *SFC, *mSFC;

  StripLabel(Label, Exact, buf, sizeof buf / sizeof (*buf));

  if (buf[0] == '\0')
    return 0;

  m    = 0;
  mSFC = NULL;

  for (SFC = SF_ALL_CAA; SFC->Label; SFC++)
  {
    for (i = 0; buf[i] && SFC->Label[i]; i++)
      if (buf[i] != toupper(SFC->Label[i])) break;

    if (buf[i] == '\0' && SFC->Label[i] == '\0')
    {
      //if (SFC->c) return SFC;

      for (i = 0; i < (sizeof SFC->g / sizeof (*SFC->g)); i++)
        if (SFC->g[i]) return SFC;

      for (i = 0; i < (sizeof SFC->f / sizeof (*SFC->f)); i++)
        if (SFC->f[i]) return SFC;

      return NULL;
    }

    if (i == 1 && isalpha(SFC->Label[1]))
      i = 0;

    if (i > m && ! isdigit(SFC->Label[i - 1]))
    {
      m    = i;
      mSFC = SFC;
    }
  }

  if (Exact)
    return NULL;

  return mSFC;
}

double CalcSF_IT4323(const T_SF_ALL_CAA *CAA, double stol2)
{
  int     i;
  double  sf;

  // sf = CAA->c;
  sf = 0.0;

  for (i = 0; i < (sizeof CAA->f / sizeof (*CAA->f)); i++)
    sf += CAA->f[i] * exp(-CAA->g[i] * stol2);

  return sf;
}

//////////////////////////////////////////////////////////////////////

const T_SF_IT92_CAA *FindSF_IT92_CAA(const char *Label, int Exact)
{
  int                  i, m;
  char                 buf[6];
  const T_SF_IT92_CAA  *SFC, *mSFC;


  StripLabel(Label, Exact, buf, sizeof buf / sizeof (*buf));

  if (buf[0] == '\0')
    return 0;

  m    = 0;
  mSFC = NULL;

  for (SFC = SF_IT92_CAA; SFC->Label; SFC++)
  {
    for (i = 0; buf[i] && SFC->Label[i]; i++)
      if (buf[i] != toupper(SFC->Label[i])) break;

    if (buf[i] == '\0' && SFC->Label[i] == '\0')
    {
      if (SFC->c) return SFC;

      for (i = 0; i < (sizeof SFC->b / sizeof (*SFC->b)); i++)
        if (SFC->b[i]) return SFC;

      for (i = 0; i < (sizeof SFC->a / sizeof (*SFC->a)); i++)
        if (SFC->a[i]) return SFC;

      return NULL;
    }

    if (i == 1 && isalpha(SFC->Label[1]))
      i = 0;

    if (i > m && ! isdigit(SFC->Label[i - 1]))
    {
      m    = i;
      mSFC = SFC;
    }
  }

  if (Exact)
    return NULL;

  return mSFC;
}


double CalcSF_IT92_CAA(const T_SF_IT92_CAA *CAA, double stol2)
{
  int     i;
  double  sf;


  sf = CAA->c;

  for (i = 0; i < (sizeof CAA->a / sizeof (*CAA->a)); i++)
    sf += CAA->a[i] * exp(-CAA->b[i] * stol2);

  return sf;
}


const T_SF_ALL_CAA *FindSF_WK95_CAA(const char *Label, int Exact)
{
  int                  i, m;
  char                 buf[6];
  const T_SF_ALL_CAA  *SFC, *mSFC;


  StripLabel(Label, Exact, buf, sizeof buf / sizeof (*buf));

  if (buf[0] == '\0')
    return 0;

  m    = 0;
  mSFC = NULL;

  for (SFC = SF_ALL_CAA; SFC->Label; SFC++)
  {
    for (i = 0; buf[i] && SFC->Label[i]; i++)
      if (buf[i] != toupper(SFC->Label[i])) break;

    if (buf[i] == '\0' && SFC->Label[i] == '\0')
    {
      if (SFC->c) return SFC;

      for (i = 0; i < (sizeof SFC->b / sizeof (*SFC->b)); i++)
        if (SFC->b[i]) return SFC;

      for (i = 0; i < (sizeof SFC->a / sizeof (*SFC->a)); i++)
        if (SFC->a[i]) return SFC;

      return NULL;
    }

    if (i == 1 && isalpha(SFC->Label[1]))
      i = 0;

    if (i > m && ! isdigit(SFC->Label[i - 1]))
    {
      m    = i;
      mSFC = SFC;
    }
  }

  if (Exact)
    return NULL;

  return mSFC;
}


double CalcSF_WK95_CAA(const T_SF_ALL_CAA *CAA, double stol2)
{
  int     i;
  double  sf;
  //printf("RAWR\n");

  sf = CAA->c;

  for (i = 0; i < (sizeof CAA->a / sizeof (*CAA->a)); i++)
    sf += CAA->a[i] * exp(-CAA->b[i] * stol2);

  return sf;
}


const T_AtomRadius *FindAtomRadius(const char *Label, int Exact)
{
  int                 i, m;
  char                buf[5];
  const T_AtomRadius  *LAR, *mLAR;


  StripLabel(Label, Exact, buf, sizeof buf / sizeof (*buf));

  if (buf[0] == '\0')
    return NULL;

  m    = 0;
  mLAR = NULL;

  for (LAR = ListAtomRadius; LAR->Label; LAR++)
  {
    for (i = 0; buf[i] && LAR->Label[i]; i++)
      if (buf[i] != toupper(LAR->Label[i])) break;

    if (buf[i] == '\0' && LAR->Label[i] == '\0')
      return LAR;

    if (i == 1 && isalpha(LAR->Label[1]))
      i = 0;

    if (i > m && ! isdigit(LAR->Label[i - 1]))
    {
      m    = i;
      mLAR = LAR;
    }
  }

  if (Exact)
    return NULL;

  return mLAR;
}


const T_ChXrayWaveLength *ChXrayWaveLengthOf(const char *Label)
{
  const char                *s, *t;
  const T_ChXrayWaveLength  *cxw;

  // only called if <<Lambda $label>> is given
  // returns wavelength struct from atominfo.h: const T_ChXrayWaveLength ListChXrayWaveLengths[]

  for (cxw = ListChXrayWaveLengths; cxw->Label; cxw++)
  {
    for (s = cxw->Label, t = Label; *s || *t; s++, t++)
      if (toupper(*s) != toupper(*t)) break;

    if (! (*s || *t)) 
      {
      return cxw;
      }
  }

  return NULL;
}







#if STAND_ALONE  /* main function for debugging purposes */


#include <stdlib.h>
#include <string.h>


static char *progn = "atominfo";


static int str_icmp(const char *s, const char *t)
{
  char  cs, ct;


  while (*s || *t)
  {
    cs = toupper(*s++);
    ct = toupper(*t++);
    if (cs < ct) return -1;
    if (cs > ct) return  1;
  }

  return 0;
}


static void ShowSF_IT92_CAA(const T_SF_IT92_CAA *CAA)
{
  int     i, n;


  n = sizeof CAA->a / sizeof (*CAA->a);

  fputs("  IT92-CAA", stdout);

  for (i = 0; i < n; i++)
    fprintf(stdout, " a%1$d b%1$d", i + 1);

  fputs(" c\n   ", stdout);

  for (i = 0; i < (sizeof CAA->a / sizeof (*CAA->a)); i++)
    fprintf(stdout, " %.6g %.6g", CAA->a[i], CAA->b[i]);

  fprintf(stdout, " %.6g\n", CAA->c);
}


static void ShowSF_WK95_CAA(const T_SF_ALL_CAA *CAA)
{
  int     i, n;


  n = sizeof CAA->a / sizeof (*CAA->a);

  fputs("  WK95-CAA", stdout);

  for (i = 0; i < n; i++)
    fprintf(stdout, " a%1$d b%1$d", i + 1);

  fputs(" c\n   ", stdout);

  for (i = 0; i < (sizeof CAA->a / sizeof (*CAA->a)); i++)
    fprintf(stdout, " %.6g %.6g", CAA->a[i], CAA->b[i]);

  fprintf(stdout, " %.6g\n", CAA->c);
}


static void ShowSF_IT4322(const T_SF_ALL_CAA *CAA)
{
  int     i, n;


  n = sizeof CAA->d / sizeof (*CAA->d);

  fputs("  IT4322", stdout);

  for (i = 0; i < n; i++)
    fprintf(stdout, " a%1$d b%1$d", i + 1);

  fputs("\n", stdout);

  for (i = 0; i < (sizeof CAA->d / sizeof (*CAA->d)); i++)
    fprintf(stdout, " %.6g %.6g", CAA->d[i], CAA->e[i]);
  fputs("\n", stdout);

  //fprintf(stdout, " %.6g\n", CAA->c);
}

static void ShowSF_IT4323(const T_SF_ALL_CAA *CAA)
{
  int     i, n;


  n = sizeof CAA->f / sizeof (*CAA->f);

  fputs("  IT4323", stdout);

  for (i = 0; i < n; i++)
    fprintf(stdout, " a%1$d b%1$d", i + 1);

  fputs("\n", stdout);

  for (i = 0; i < (sizeof CAA->f / sizeof (*CAA->f)); i++)
    fprintf(stdout, " %.6g %.6g", CAA->f[i], CAA->g[i]);
  fputs("\n", stdout);

  //fprintf(stdout, " %.6g\n", CAA->c);
}



static void usage(void)
{
  fprintf(stderr,
    "usage: %s [-IT92|WK95|IT4322|IT4323] [-ExactLabel] AtomLabel\n"
    "                [Start_sinTovL [End_sinTovL [StepSize]]]\n",
    progn);

  exit(1);
}


int main(int argc, char *argv[])
{
  int     iarg, i, n;
  int     iSES;
  int     F_IT92, F_WK95, F_ExactLabel, F_IT4322, F_IT4323;
  char    *AtomLabel, xtrac;
  double  Start_stol, End_stol, StepSize;
  double  stol, sfIT92, sfWK95, val, sfIT4322, sfIT4323;

  const T_PSE          *ePSE;
  const T_SF_IT92_CAA  *SF_IT92;
  const T_SF_ALL_CAA   *SF_ALL_CAA;
  const T_AtomRadius   *AR;

  const T_SF_ALL_CAA   *SF_WK95;
  const T_SF_ALL_CAA   *SF_IT4322;
  const T_SF_ALL_CAA   *SF_IT4323;


  F_IT92 = 0;
  F_WK95 = 0;
  F_IT4322 = 0; 
  F_IT4323 = 0; 
  F_ExactLabel = 0;
  AtomLabel = NULL;

  Start_stol = 0.;
  End_stol = 1.;
  StepSize = .05;

  iSES = -1;

  for (iarg = 1; iarg < argc; iarg++)
  {
    if (argv[iarg][0] == '-' && isdigit(argv[iarg][1]) == 0)
    {
      if (str_icmp(argv[iarg], "-ExactLabel") == 0)
        F_ExactLabel = 1;
      else if (str_icmp(argv[iarg], "-IT92") == 0)
        F_IT92 = 1;
      else if (str_icmp(argv[iarg], "-WK95") == 0)
        F_WK95 = 1;
      else if (str_icmp(argv[iarg], "-IT4322") == 0)
        F_IT4322 = 1;
      else if (str_icmp(argv[iarg], "-IT4323") == 0)
        F_IT4323 = 1;
      else
        usage();
    }
    else if (AtomLabel == NULL)
    {
      AtomLabel = argv[iarg];
      iSES = 0;
    }
    else if (iSES >= 0 && iSES < 3)
    {
          n = sscanf(argv[iarg], "%lf", &val, &xtrac);
      if (n != 1 && ! (n == 2 && isspace(xtrac)))
        usage();

      if      (iSES == 0)
      {
            Start_stol = val;
        if (Start_stol < 0.) {
          fprintf(stderr, "%s: illegal parameter: Start_sinTovL\n", progn);
          usage();
        }
      }
      else if (iSES == 1)
      {
            End_stol = val;
        if (End_stol < 0. || End_stol >= 100.) {
          fprintf(stderr, "%s: illegal parameter: End_sinTovL\n", progn);
          usage();
        }
      }
      else
      {
            StepSize = val;
        if (StepSize <= 0.) {
          fprintf(stderr, "%s: illegal parameter: StepSize\n", progn);
          usage();
        }
      }

      iSES++;
    }
    else
      usage();
  }

  if (AtomLabel == NULL)
    usage();

  putc('#', stdout);
  for (i = 0; i < argc;)
    fprintf(stdout, " %s", argv[i++]);
  putc('\n', stdout);

  ePSE    = FindInPSE(AtomLabel, F_ExactLabel);
  SF_IT92 = FindSF_IT92_CAA(AtomLabel, F_ExactLabel);
  SF_WK95 = FindSF_WK95_CAA(AtomLabel, F_ExactLabel);
  SF_IT4322 = FindSF_IT4322(AtomLabel, F_ExactLabel);
  SF_IT4323 = FindSF_IT4323(AtomLabel, F_ExactLabel);
  AR      = FindAtomRadius(AtomLabel, F_ExactLabel);

  if (ePSE)
    fprintf(stdout, "%d \"%s\" \"%s\" %.3f\n",
      ePSE->Z, ePSE->Symbol, ePSE->Name, ePSE->Mass);
  else
    fprintf(stdout, "not found in PSE\n");

  fprintf(stdout,
    "\nScattering Factor Label (IT92):");

  if (SF_IT92) {
    fprintf(stdout, " \"%s\"\n", SF_IT92->Label);
    ShowSF_IT92_CAA(SF_IT92);
  }
  else
    fprintf(stdout, " not found\n");

  fprintf(stdout,
    "\nScattering Factor Label (WK95):");

  if (SF_WK95) {
    fprintf(stdout, " \"%s\"\n", SF_WK95->Label);
    ShowSF_WK95_CAA(SF_WK95);
  }
  else
    fprintf(stdout, " not found\n");
  
  fprintf(stdout,
    "\nScattering Factor Label (IT4322):");

  if (SF_IT4322) {
    fprintf(stdout, " \"%s\"\n", SF_IT4322->Label);
    ShowSF_IT4322(SF_IT4322);
  }
  else
    fprintf(stdout, " not found\n");
  
  fprintf(stdout,
    "\nScattering Factor Label (IT4323):");

  if (SF_IT4323) {
    fprintf(stdout, " \"%s\"\n", SF_IT4323->Label);
    ShowSF_IT4323(SF_IT4323);
  }
  else
    fprintf(stdout, " not found\n");

  fprintf(stdout, "\nAtom Radius");

  if (AR)
    fprintf(stdout, " %.2f \"%s\"\n", AR->Radius, AR->Label);
  else
    fprintf(stdout, " not found\n");

  if (F_IT92 == 0 && F_WK95 == 0 && F_IT4322 == 0 && F_IT4323 == 0) {
    F_IT92   = 1;
    F_WK95   = 1;
    F_IT4322= 1;
    F_IT4323= 1;
  }

  if ((SF_IT92 || SF_WK95 || SF_IT4322 || SF_IT4323) && iSES > 0)
  {
      fprintf(stdout, "   stol");

    if (SF_IT92 && F_IT92 && SF_IT4322 && SF_IT4322)
      fprintf(stdout, "     IT92");

    if (SF_IT92 && F_IT92 && SF_IT4322 && SF_IT4322)
      fprintf(stdout, "     WK95");

    if (SF_IT92 && F_IT92 && SF_IT4322 && SF_IT4322)
      fprintf(stdout, "   IT4322");
    
    if (SF_IT92 && F_IT92 && SF_IT4322 && SF_IT4322)
      fprintf(stdout, "   IT4323");


    putc('\n', stdout);

    sfIT92 = 0.;
    sfWK95 = 0.;
    sfIT4322 = 0.;
    sfIT4323 = 0.;

    for (i = 0;                                         /* ARBITRARY */
         (stol = Start_stol + i * StepSize) <= (End_stol + 1.e-4);
         i++
        )
    {
      fprintf(stdout, "  %5.3f", stol);

      if (SF_IT92 && F_IT92 && SF_IT4322 && SF_IT4323) {
        sfIT92 = CalcSF_IT92_CAA(SF_IT92, stol * stol);
        fprintf(stdout, "  %7.3f", sfIT92);
      }

      if (SF_IT92 && F_IT92 && SF_IT4322 && SF_IT4323) {
        sfWK95 = CalcSF_WK95_CAA(SF_WK95, stol * stol);
        fprintf(stdout, "  %7.3f", sfWK95);
      }

      if (SF_IT92 && F_IT92 && SF_IT4322 && SF_IT4323) {
        sfIT4322 = CalcSF_IT4322(SF_IT4322, stol * stol);
        fprintf(stdout, "  %7.3f", sfIT4322);
      }
      
      if (SF_IT92 && F_IT92 && SF_IT4322 && SF_IT4323) {
        sfIT4323 = CalcSF_IT4323(SF_IT4323, stol * stol);
        fprintf(stdout, "  %7.3f", sfIT4323);
      }


      putc('\n', stdout);
    }
  }

  return 0;
}


#endif  /* STAND_ALONE */
