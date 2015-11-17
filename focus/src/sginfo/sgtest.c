/*
  Space Group Info (c) 1994-97 Ralf W. Grosse-Kunstleve
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>


#define SGCOREDEF__
#include "sginfo.h"


/* sghklvfy.c */

int Verify_sghkl(const T_SgInfo *SgInfo, int FriedelSym,
                 int Maxh, int Maxk, int Maxl);


#define Fprintf (void) fprintf
#define Fflush  (void) fflush


static int str_icmp(char *s, char *t)
{
  char     cs, ct;

  while (*s || *t)
  { cs = toupper(*s++);
    ct = toupper(*t++);
    if (cs < ct) return -1;
    if (cs > ct) return  1;
  }
  return 0;
}


static char *SgInfoCpy(T_SgInfo *Target, T_SgInfo *Source)
{
  int          ErrCount;
  T_RTMx       *ListSeitzMx;
  T_RotMxInfo  *ListRotMxInfo;
  T_hklCond    *ReflCond;
  T_hklCond    *RestCond;
  int          *SysEnhanced;


  ErrCount = 0;

  ListSeitzMx   = Target->ListSeitzMx;
  ListRotMxInfo = Target->ListRotMxInfo;
  ReflCond      = Target->ReflCond;
  RestCond      = Target->RestCond;
  SysEnhanced   = Target->SysEnhanced;

  (void) memcpy(Target, Source, sizeof (*Target));

  Target->ListSeitzMx   = ListSeitzMx;
  Target->ListRotMxInfo = ListRotMxInfo;
  Target->ReflCond      = ReflCond;
  Target->RestCond      = RestCond;
  Target->SysEnhanced   = SysEnhanced;

  (void) memcpy(Target->ListSeitzMx, Source->ListSeitzMx,
                Source->nList * sizeof (*(Target->ListSeitzMx)));

  if (Target->ListRotMxInfo != NULL)
  {
    if (Source->ListRotMxInfo == NULL) {
      (void) memset(Target->ListRotMxInfo, 0,
                    Source->MaxList * sizeof (*(Target->ListRotMxInfo)));
      ErrCount++;
    }
    else
      (void) memcpy(Target->ListRotMxInfo, Source->ListRotMxInfo,
                    Source->MaxList * sizeof (*(Target->ListRotMxInfo)));
  }

  if (Target->ReflCond != NULL)
  {
    if (Source->ReflCond == NULL) {
      (void) memset(Target->ReflCond, 0,
                    Source->MaxList * sizeof (*(Target->ReflCond)));
      ErrCount++;
    }
    else
      (void) memcpy(Target->ReflCond, Source->ReflCond,
                    Source->MaxList * sizeof (*(Target->ReflCond)));
  }

  if (Target->RestCond != NULL)
  {
    if (Source->RestCond == NULL) {
      (void) memset(Target->RestCond, 0,
                    Source->MaxList * sizeof (*(Target->RestCond)));
      ErrCount++;
    }
    else
      (void) memcpy(Target->RestCond, Source->RestCond,
                    Source->MaxList * sizeof (*(Target->RestCond)));
  }

  if (Target->SysEnhanced != NULL)
  {
    if (Source->SysEnhanced == NULL) {
      (void) memset(Target->SysEnhanced, 0,
                    Source->MaxList * sizeof (*(Target->SysEnhanced)));
      ErrCount++;
    }
    else
      (void) memcpy(Target->SysEnhanced, Source->SysEnhanced,
                    Source->MaxList * sizeof (*(Target->SysEnhanced)));
  }

  if (ErrCount)
    return "Internal Error: SgInfoCpy()";

  return NULL;
}


static const char *progn = "sgtest";


static void progerror(const char *message)
{
  Fflush(stdout);
  Fprintf(stderr, "%s: %s\n", progn, message);
  exit(1);
}


static void NotEnoughCore(void)
{
  progerror("Not enough core");
}


static void AllocSgInfo(T_SgInfo *SgInfo)
{
  SgInfo->MaxList = 192;

  SgInfo->ListSeitzMx
    = malloc(SgInfo->MaxList * sizeof (*SgInfo->ListSeitzMx));
  if (SgInfo->ListSeitzMx == NULL) NotEnoughCore();

  SgInfo->ListRotMxInfo
    = malloc(SgInfo->MaxList * sizeof (*SgInfo->ListRotMxInfo));
  if (SgInfo->ListRotMxInfo == NULL) NotEnoughCore();

  SgInfo->ReflCond
    = malloc(SgInfo->MaxList * sizeof (*SgInfo->ReflCond));
  if (SgInfo->ReflCond == NULL) NotEnoughCore();

  SgInfo->RestCond
    = malloc(SgInfo->MaxList * sizeof (*SgInfo->RestCond));
  if (SgInfo->RestCond == NULL) NotEnoughCore();

  SgInfo->SysEnhanced
    = malloc(SgInfo->MaxList * sizeof (*SgInfo->SysEnhanced));
  if (SgInfo->SysEnhanced == NULL) NotEnoughCore();
}


static void usage(void)
{
  Fprintf(stderr,
    "usage: %s [-v] [-Verify] [-light] [#Start [#End]]\n",
    progn);
  exit (1);
}


int main(int argc, char *argv[])
{
  int                F_Verbose, F_Verify, F_light;
  int                SgNumStart, SgNumEnd;
  int                i, n, pos_hsym, OrSh[3];
  int                nNextBasis, iNextBasis;
  char               xtrac;
  T_SgInfo           SpgrInfo0, SpgrInfoN, SpgrInfo, BC_SgInfo;
  const T_TabSgName  *TSgN, *ConvTSgN;
  T_RTMx             CBMx, InvCBMx;


  F_Verbose = 0;
  F_Verify  = 0;
  F_light   = 0;

  SgNumStart = -1;
  SgNumEnd   = -1;

  for (i = 1; i < argc; i++)
  {
    if      (str_icmp(argv[i], "-v") == 0)
      F_Verbose = 1;
    else if (str_icmp(argv[i], "-Verify") == 0)
      F_Verify = 1;
    else if (str_icmp(argv[i], "-light") == 0)
      F_light = 1;
    else if (SgNumStart == -1)
    {
          n = sscanf(argv[i], "%d %c", &SgNumStart, &xtrac);
      if (n != 1 || SgNumStart < 1) usage();
    }
    else if (SgNumEnd   == -1)
    {
          n = sscanf(argv[i], "%d %c", &SgNumEnd, &xtrac);
      if (n != 1 || SgNumEnd < 1) usage();
    }
    else
      usage();
  }

  if (SgNumEnd < 0) {
    if (SgNumStart < 0) {
      SgNumStart =   1;
      SgNumEnd   = 230;
    }
    else
      SgNumEnd   = SgNumStart;
  }

  AllocSgInfo(&SpgrInfo0);
  AllocSgInfo(&SpgrInfoN);
  AllocSgInfo(&SpgrInfo);
  AllocSgInfo(&BC_SgInfo);

  for (TSgN = TabSgName; TSgN->HallSymbol; TSgN++)
  {
    if (TSgN->SgNumber < SgNumStart)
      continue;

    if (TSgN->SgNumber > SgNumEnd)
      break;

    ResetSgInfo(&SpgrInfo0);

    pos_hsym = ParseHallSymbol(TSgN->HallSymbol, &SpgrInfo0);

    if (SgError != NULL)
    {
      Fprintf(stdout, "    %s\n", TSgN->HallSymbol);
      for (i = 0; i < pos_hsym; i++) putc('-', stdout);
      Fprintf(stdout, "---^\n");
      Fprintf(stdout, "%s\n", SgError);
      SgError = NULL;
      continue;
    }
                     i = VolAPointGroups[TSgN->SgNumber];
    if (   PG_Number(i) >= PG_Number(PG_4)
        && PG_Number(i) <= PG_Number(PG_6_mmm)
        && TSgN->Extension[0] != 'R')
      nNextBasis = 3;
    else
      nNextBasis = 1;

    for (iNextBasis = 0; iNextBasis < nNextBasis; iNextBasis++)
    {
      Fprintf(stdout, "TSgN  ");
      PrintTabSgNameEntry(TSgN, 0, 0, 0, stdout);

      i = PG_Index(VolAPointGroups[TSgN->SgNumber]);

      Fprintf(stdout, " : %s", PG_Names[i]);
      Fprintf(stdout, " : %s", PG_Names[PG_Index(LG_Code_of_PG_Index[i])]);

      if      (iNextBasis == 0)
      {
        if (nNextBasis != 1)
          Fprintf(stdout, " : z->z");

        SgError = SgInfoCpy(&SpgrInfoN, &SpgrInfo0);
      }
      else if (iNextBasis == 1)
      {
        Fprintf(stdout, " : z->x");

        for (i = 0; i < 9;  i++) {
             CBMx.s.R[i] = CRBF * RMx_3_111[i];
          InvCBMx.s.R[i] = CRBF * RMx_3i111[i];
        }

        for (i = 0; i < 3; i++) {
             CBMx.s.T[i] = 0;
          InvCBMx.s.T[i] = 0;
        }

        ResetSgInfo(&SpgrInfoN);
        (void) TransformSgInfo(&SpgrInfo0, &CBMx, &InvCBMx, &SpgrInfoN);
      }
      else if (iNextBasis == 2)
      {
        Fprintf(stdout, " : z->y");

        for (i = 0; i < 9;  i++) {
             CBMx.s.R[i] = CRBF * RMx_3i111[i];
          InvCBMx.s.R[i] = CRBF * RMx_3_111[i];
        }

        for (i = 0; i < 3; i++) {
             CBMx.s.T[i] = 0;
          InvCBMx.s.T[i] = 0;
        }

        ResetSgInfo(&SpgrInfoN);
        (void) TransformSgInfo(&SpgrInfo0, &CBMx, &InvCBMx, &SpgrInfoN);
      }

      putc('\n', stdout);
      Fflush(stdout);

      if (SgError)
      {
        Fprintf(stdout, "@ (%d) %s => %s\n",
          TSgN->SgNumber, TSgN->HallSymbol, SgError);

        SgError = NULL;
        goto Next_TSgN;
      }

      for (OrSh[0] = 0; OrSh[0] < 12; OrSh[0]++)
      for (OrSh[1] = 0; OrSh[1] < 12; OrSh[1]++)
      for (OrSh[2] = 0; OrSh[2] < 12; OrSh[2]++)
      {
            SgError = SgInfoCpy(&SpgrInfo, &SpgrInfoN);
        if (SgError)
        {
          Fprintf(stdout, "@ (%d) %s => %s\n",
            TSgN->SgNumber, TSgN->HallSymbol, SgError);

          exit(1);
        }

        if (F_light == 0)
          for (i = 0; i < 3; i++)
            SpgrInfo.OriginShift[i] = OrSh[i];

        if (CompleteSgInfo(&SpgrInfo) != 0)
        {
          if (SgError != NULL)
          {
            Fprintf(stdout, "@ (%d) %s (%d %d %d) => %s\n",
              TSgN->SgNumber, TSgN->HallSymbol,
              OrSh[0], OrSh[1], OrSh[2],
              SgError);

            SgError = NULL;
            goto Next_TSgN;
          }
        }

        if (F_Verbose)
          Fprintf(stdout, "Hall Symbol  %s\n", SpgrInfo.HallSymbol);

        if (F_Verify && SpgrInfo.ReflCond)
        {
          if (SetReflCond(&SpgrInfo) >= 0)
          {
            if (   SpgrInfo.RestCond
                && SetRestCond(&SpgrInfo) < 0)
            {
              Fprintf(stdout, "@ (%d) %s (%d %d %d) => %s\n",
                TSgN->SgNumber, SpgrInfo.HallSymbol,
                OrSh[0], OrSh[1], OrSh[2],
                SgError);

              SgError = NULL;
            }

            if (SpgrInfo.SysEnhanced)
              (void) SetSysEnhanced(&SpgrInfo);
          }

          if (SgError != NULL)
          {
            Fprintf(stdout, "@ (%d) %s (%d %d %d) => %s\n",
              TSgN->SgNumber, SpgrInfo.HallSymbol,
              OrSh[0], OrSh[1], OrSh[2],
              SgError);

            SgError = NULL;
          }

                   i = Verify_sghkl(&SpgrInfo, 4, abs(4), abs(4), abs(4));
          if      (i < 0)
          {
            Fprintf(stdout, "@ (%d) %s (%d %d %d) => %s\n",
              TSgN->SgNumber, SpgrInfo.HallSymbol,
              OrSh[0], OrSh[1], OrSh[2],
              SgError);

            SgError = NULL;
          }
          else if (i == 0 && F_Verbose)
            Fprintf(stdout, "Verify o.k.: sghkl\n");
        }

        ConvTSgN = FindReferenceSpaceGroup(&SpgrInfo, &CBMx, &InvCBMx);

        if (SgError)
        {
          Fprintf(stdout, "@ (%d) %s (%d %d %d) => %s => %s\n",
            TSgN->SgNumber, TSgN->HallSymbol,
            OrSh[0], OrSh[1], OrSh[2],
            SpgrInfo.HallSymbol, SgError);

          SgError = NULL;
          goto Next_TSgN;
        }
        else
        {
          if (F_Verbose) {
            PrintMapleRTMx(   &CBMx, CRBF, CTBF,    "CBMx", stdout);
            PrintMapleRTMx(&InvCBMx, CRBF, CTBF, "InvCBMx", stdout);
          }

          ResetSgInfo(&BC_SgInfo);

              i = TransformSgInfo(&SpgrInfo, &CBMx, &InvCBMx, &BC_SgInfo);
          if (i == 0)
              i = CompleteSgInfo(&BC_SgInfo);
          if (i != 0) {
            Fprintf(stdout, "@ (%d) %s (%d %d %d) => %s => (%d) %s => %s\n",
              TSgN->SgNumber, TSgN->HallSymbol,
              OrSh[0], OrSh[1], OrSh[2],
              SpgrInfo.HallSymbol,
              ConvTSgN->SgNumber, ConvTSgN->HallSymbol,
              SgError);

            SgError = NULL;
            goto Next_TSgN;
          }
          else if (BC_SgInfo.TabSgName == NULL) {
            Fprintf(stdout, "@ (%d) %s (%d %d %d) => %s => (%d) %s => %s\n",
              TSgN->SgNumber, TSgN->HallSymbol,
              OrSh[0], OrSh[1], OrSh[2],
              SpgrInfo.HallSymbol,
              ConvTSgN->SgNumber, ConvTSgN->HallSymbol,
              BC_SgInfo.HallSymbol);

            goto Next_TSgN;
          }
          else
          {
            if (   BC_SgInfo.TabSgName->SgNumber != ConvTSgN->SgNumber
                || BC_SgInfo.TabSgName->SgNumber != TSgN->SgNumber)
            {
              Fprintf(stdout, "@ (%d) %s (%d %d %d) => %s => (%d) %s => ",
                TSgN->SgNumber, TSgN->HallSymbol,
                OrSh[0], OrSh[1], OrSh[2],
                SpgrInfo.HallSymbol,
                ConvTSgN->SgNumber, ConvTSgN->HallSymbol);

              PrintTabSgNameEntry(BC_SgInfo.TabSgName, 0, 0, 0, stdout);
              putc('\n', stdout);

              goto Next_TSgN;
            }

            if (F_Verbose)
            {
              Fprintf(stdout, "Change of Basis => ");
              PrintTabSgNameEntry(BC_SgInfo.TabSgName, 0, 0, 0, stdout);
              putc('\n', stdout);
            }
          }
        }

        if (F_light)
          goto Next_iNextBase;
      }

      Next_iNextBase:;
    }

    Next_TSgN:;
  }

  return 0;
}
