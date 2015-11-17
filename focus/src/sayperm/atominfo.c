#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "atominfodecl.h"
#include "atominfo.h"


#define ATOMINFO_C__i




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
  int i = 0;
  double sf = 0.0;
  double sf1 = 0.0;
  double b = 0.0;
  double stol = 0.0;
 double prod = 0.0;


  sf = CAA->c;
  fprintf(stdout," in CalcSF_IT92 %5.3lf\n", sf);


  for (i = 0; i < (sizeof CAA->a / sizeof (*CAA->a)); i++){
    b = CAA->b[i];
    stol = stol2;
    /*sf += CAA->a[i] * exp(-CAA->b[i] * stol2);*/
    prod = -b * stol;
    sf1 = exp(prod);
    sf += CAA->a[i] * sf1;
  fprintf(stdout," in CalcSF_IT92 %3d %3d %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf %5.3lf\n", sizeof CAA->a, sizeof (*CAA->a), CAA->a[i], CAA->b[i], stol2, sf, sf1, prod);
 b = 0.0;
 stol = 0.0;
}

  return sf;
}


const T_SF_WK95_CAA *FindSF_WK95_CAA(const char *Label, int Exact)
{
  int                  i, m;
  char                 buf[6];
  const T_SF_WK95_CAA  *SFC, *mSFC;


  StripLabel(Label, Exact, buf, sizeof buf / sizeof (*buf));

  if (buf[0] == '\0')
    return 0;

  m    = 0;
  mSFC = NULL;

  for (SFC = SF_WK95_CAA; SFC->Label; SFC++)
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


double CalcSF_WK95_CAA(const T_SF_WK95_CAA *CAA, double stol2)
{
  int     i;
  double  sf;


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


  for (cxw = ListChXrayWaveLengths; cxw->Label; cxw++)
  {
    for (s = cxw->Label, t = Label; *s || *t; s++, t++)
      if (toupper(*s) != toupper(*t)) break;

    if (! (*s || *t)) return cxw;
  }

  return NULL;
}




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


static void ShowSF_WK95_CAA(const T_SF_WK95_CAA *CAA)
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


static void usage(void)
{
  fprintf(stderr,
    "usage: %s [-IT92|WK95] [-ExactLabel] AtomLabel\n"
    "                [Start_sinTovL [End_sinTovL [StepSize]]]\n",
    progn);

  exit(1);
}

#if STAND_ALONE  /* main function for debugging purposes */


int main(int argc, char *argv[])
{
  int     iarg, i, n;
  int     iSES;
  int     F_IT92, F_WK95, F_ExactLabel;
  char    *AtomLabel, xtrac;
  double  Start_stol, End_stol, StepSize;
  double  stol, sfIT92, sfWK95, val;

  const T_PSE          *ePSE;
  const T_SF_IT92_CAA  *SF_IT92;
  const T_SF_WK95_CAA  *SF_WK95;
  const T_AtomRadius   *AR;


  F_IT92 = 0;
  F_WK95 = 0;
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
  AR      = FindAtomRadius(AtomLabel, F_ExactLabel);

  if (ePSE)
    fprintf(stdout, "%d \"%s\" \"%s\" %.3f\n",
      ePSE->Z, ePSE->Symbol, ePSE->Name, ePSE->Mass);
  else
    fprintf(stdout, "not found in PSE\n");

  fprintf(stdout,
    "Scattering Factor Label (IT92):");

  if (SF_IT92) {
    fprintf(stdout, " \"%s\"\n", SF_IT92->Label);
    ShowSF_IT92_CAA(SF_IT92);
  }
  else
    fprintf(stdout, " not found\n");

  fprintf(stdout,
    "Scattering Factor Label (WK95):");

  if (SF_WK95) {
    fprintf(stdout, " \"%s\"\n", SF_WK95->Label);
    ShowSF_WK95_CAA(SF_WK95);
  }
  else
    fprintf(stdout, " not found\n");

  fprintf(stdout, "Atom Radius");

  if (AR)
    fprintf(stdout, " %.2f \"%s\"\n", AR->Radius, AR->Label);
  else
    fprintf(stdout, " not found\n");

  if (F_IT92 == 0 && F_WK95 == 0) {
    F_IT92 = 1;
    F_WK95 = 1;
  }

  if ((SF_IT92 || SF_WK95) && iSES > 0)
  {
      fprintf(stdout, "   stol");

    if (SF_IT92 && F_IT92)
      fprintf(stdout, "     IT92");

    if (SF_WK95 && F_WK95)
      fprintf(stdout, "     WK95");

    putc('\n', stdout);

    sfIT92 = 0.;
    sfWK95 = 0.;

    for (i = 0;                                         /* ARBITRARY */
         (stol = Start_stol + i * StepSize) <= (End_stol + 1.e-4);
         i++
        )
    {
      fprintf(stdout, "  %5.3f", stol);

      if (SF_IT92 && F_IT92) {
        sfIT92 = CalcSF_IT92_CAA(SF_IT92, stol * stol);
        fprintf(stdout, "  %7.3f", sfIT92);
      }

      if (SF_WK95 && F_WK95) {
        sfWK95 = CalcSF_WK95_CAA(SF_WK95, stol * stol);
        fprintf(stdout, "  %7.3f", sfWK95);
      }

      putc('\n', stdout);
    }
  }

  return 0;
}


#endif  /* STAND_ALONE */
