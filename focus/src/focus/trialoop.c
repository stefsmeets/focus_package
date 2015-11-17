#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "main.h"
#include "xtal.h"
#include "trialoop.h"
#include "io.h"
#include "ranmar.h"
#include "nodsurf.h"


#define MemCpy(t, s, n)  memcpy((t), (s), (n) * sizeof (*(t)))


#ifdef JUNK
static void PrintPhaseCode(T_PhaseCode *code)
{
  int      iCT, iTab;

  static int TabSymbol[] =
   {
     '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
     '-', 'i', 'h', 'g', 'f', 'e', 'd', 'c', 'b', 'a',
     '0'
   };


  Fprintf(stdout, "PhaseCode(");

  for (iCT = 0; iCT < nActivePhase; iCT++)
  {
    if (iCT % 72 == 0) putc('\n', stdout);

        iTab = (code[iCT] + 9) / 18;
    if (iTab >= 0 && iTab < (sizeof TabSymbol / sizeof (*TabSymbol)))
      putc(TabSymbol[iTab], stdout);
    else
      putc('~', stdout);
  }

  Fprintf(stdout, "\n)");
}
#endif


static int NextPhaseCode(T_PhaseCode *Code, T_PhaseCode *FixedCode)
{
  int           nCR, iCT;
  T_CodeTransl  *CT;
  static int    iRandomBlindCalls = 0;


  nCR = nCallRanmar;

  if      (F_LoadDensityFileName)
  {
    for (iCT = 0; iCT < nActivePhase; iCT++) Code[iCT] = 0;
  }
  else if (F_AllPhaseCodes0)
  {
    if (FixedCode) {
      (void) MemCpy(Code, FixedCode, nActivePhase);
      for (iCT = 0; iCT < nActivePhase; iCT++)
        if (Code[iCT] < 0) Code[iCT] = 0;
    }
    else
      for (iCT = 0; iCT < nActivePhase; iCT++) Code[iCT] = 0;
  }
  else
  {
    if (iRandomBlindCalls < nRandomBlindCalls)
    {
      if (nCallRanmar > RandomBlindCalls[iRandomBlindCalls])
        DoRandomInitialization();

      while (nCallRanmar < RandomBlindCalls[iRandomBlindCalls])
        (void) ranmar();

      nCR = nCallRanmar;

      iRandomBlindCalls++;
    }

    if (FixedCode)
      (void) MemCpy(Code, FixedCode, nActivePhase);
    else
      for (iCT = 0; iCT < nActivePhase; iCT++) Code[iCT] = -1;

    CT = CodeTransl;

    for (iCT = 0; iCT < nActivePhase; iCT++, CT++)
    {
      if (Code[iCT] >= 0) continue;

      if (CT->FobsRaw->PhaseRestriction < 0)
        Code[iCT] = (int) (ranmar() * 360.) % 360;
      else
      {
        if (ranmar() < .5) Code[iCT] =   0;
        else               Code[iCT] = 180;
      }
    }
  }

  return nCR;
}


void TriaLoop(int MaximumTrials, int nTrialsPerPhaseSet)
{
  T_PhaseCode          *PhaseCode, *FixedCode;
  T_FourierParameters  FP[1];
  T_PeakFlags          *PeakFlags;
  Fprec                *Surface;
  int                  StatLoadNextPhaseSet, iTrialsPerPhaseSet;
  long                 PosLoadNextPhaseSet;

  T_Ticks  Ticks;
  long     CheckTime = -SignalFileCheckInterval;


  nTrials = 0;
  nConvergedPhases = 0;
  nNotConvergedPhases = 0;

  StatLoadNextPhaseSet = 0;
  PosLoadNextPhaseSet = 1;
  iTrialsPerPhaseSet = 0;

  CheckMalloc(PhaseCode, nActivePhase);

  FixedCode = NULL;

  if (F_PhaseSetsFileName)
    CheckMalloc(FixedCode, nActivePhase);

  InitFourierParameters(FP, CodeTransl, nActivePhase);

  CheckMalloc(PeakFlags, Nx * Ny * Nz);
  DS_MarkEquiv(&SpgrInfo, PeakFlags, Nx, Ny, Nz);

  Surface = NULL;

  if (F_SurfaceFileName) {
    CheckMalloc(Surface, Nx * Ny * Nz);
    LoadSurface(Surface, Nx, Ny, Nz, PeakFlags);
  }

  /* initialize random number generator */
  DoRandomInitialization();

     (void) GetTicks(&Ticks);
  PrintTicks(stdout, &Ticks, "# Time ", "  Start of TriaLoop\n");
  Fflush(stdout);

  for (;;)
  {
    if (F_SignalFileName)
    {
      if (Ticks.sec_since_some_day - CheckTime >= SignalFileCheckInterval)
      {
        CheckSignalFile();
        CheckTime = Ticks.sec_since_some_day;
      }
    }

    if (QuitProgram)
      break;

    if (   FixedCode
        && (   StatLoadNextPhaseSet == 0
            || iTrialsPerPhaseSet >= nTrialsPerPhaseSet))
    {
      PosLoadNextPhaseSet = LoadNextPhaseSet(NULL) + 1;

      Fprintf(stdout,
        "# Looking for next phase set, starting at line #%ld in file %s\n",
        PosLoadNextPhaseSet, F_PhaseSetsFileName);

          StatLoadNextPhaseSet = LoadNextPhaseSet(FixedCode);
      if (StatLoadNextPhaseSet == -1)
        break;

      iTrialsPerPhaseSet = 0;
    }

    Fprintf(stdout, "WoT %d", nTrials);
    if (FixedCode)
      Fprintf(stdout, " PosPhaseSet %ld Fixed %ld Random %d Trial %d",
        PosLoadNextPhaseSet,
        (long) nActivePhase - StatLoadNextPhaseSet,
        StatLoadNextPhaseSet,
        iTrialsPerPhaseSet);
    putc('\n', stdout);
    Fflush(stdout);

    CurrPhaseCode_nCallRanmar = NextPhaseCode(PhaseCode, FixedCode);
    RecyclePhases(FP, PeakFlags, Surface, PhaseCode);

    nTrials++;
    if (FixedCode) iTrialsPerPhaseSet++;

       (void) GetTicks(&Ticks);
    PrintTicks(stdout, &Ticks, "# Time ", "\n");
    Fflush(stdout);

    if (MaximumTrials && MaximumTrials <= nTrials)
      break;
  }

  if (Surface) AppFree(Surface, Nx * Ny * Nz);
  AppFree(PeakFlags, Nx * Ny * Nz);
  FreeFourierParameters(FP);

  PrintTicks(stdout, NULL, "# Time ", "  End of TriaLoop\n");
}
