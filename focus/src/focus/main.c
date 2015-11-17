#define VersionID "01.92.1"

#ifndef A_AppGlobal__
#define A_AppGlobal__
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#ifdef __TURBOC__

#include <dos.h>
unsigned _stklen = 16 * 1024;
#include <alloc.h>
#define malloc farmalloc
#define realloc farrealloc
#define free farfree

#elif ! (defined(__ALPHA) && defined(__VMS)) \
   && ! defined(__APPLE__)
#include <malloc.h>

#endif

#include "unixstd.h"

#include "main.h"
#include "lib.h"
#include "xtal.h"
#include "trialoop.h"
#include "matrix.h"
#include "io.h"
#include "ranmar.h"
#include "nodsurf.h"


static char *AppHostname(void)
{
  static int  First = 1;
  static char hostname[64];


  if (First)
  {
#if ! defined(__TURBOC__)
    if (gethostname(hostname, sizeof hostname) != 0)
#endif
      (void) strcpy(hostname, "unknown");

    First = 0;
  }

  return hostname;
}


static long AppPID(void)
{
  long PID = 0;

#if ! (defined(__TURBOC__) || (defined(__ALPHA) && defined(__VMS)))
    PID = getpid();
#endif

  return PID;
}


static void InitializeParameters(void)
{
  progn = "focus";
  Debug0 = 0; /* use "variables" to avoid compiler message */
  Debug1 = 1; /* indicating unreachable code               */
  MemoryUsed = 0;
  MaxMemoryUsed = 0;
  QuitProgram = 0;
  ReportOnExit = 0;
  CountWarnings = 0;

  RandomInitialization = 0;
  FlagTimeRandomInitialization = 1;
  nRandomBlindCalls = 0;
   RandomBlindCalls = NULL;
  LargestFwFragment.nAsyN    = 0;
  LargestFwFragment.nSymN    = 0;
  LargestFwFragment.mBpN     = 0.;
  LargestFwFragment.Low_iEPL = 0;
  LargestFwFragment.Top_iEPL = 0;
  LargestFwFragment.LnB  = NULL;
  nFeedBackCycles = 0;
   FeedBackCycles = NULL;
  FeedBackBreakIf.PhaseDiff.Value = .05;
  FeedBackBreakIf.PhaseDiff.Percent = 1;
  FeedBackBreakIf.PhaseDiff.Weighted = 0;
  FeedBackBreakIf.DeltaR.Value = 0.01;
  FeedBackBreakIf.DeltaR.Percent = 1;
  FeedBackBreakIf.DeltaR.Weighted = 0;
  FeedBackBreakIf.Branch = BrkIf_Branch_PhaseDiff | BrkIf_Branch_DeltaR;
  nTrials = 0;
  nConvergedPhases = 0;
  nNotConvergedPhases = 0;

  TopTitle = NULL;

  Ilinec = 0;
  Imode = 0;
  Itabsize = 8;
  *Iline = '\0';

  FobsMaxQ = 1.;
  nFobsRaw = 0;
  FobsRaw = NULL;
  ResetSgInfo(&SpgrInfo);
  fTrVector = NULL;
  List_fSeitzMx = NULL;
  nInputLatConD = 6;
  LatConD.a = LatConD.b = LatConD.c = LatConD.v = 1.;
  LatConD.alpha = LatConD.beta = LatConD.gamma = 90. * PIover180;
  LatConR.a = LatConR.b = LatConR.c = LatConR.v = 1.;
  LatConR.alpha = LatConR.beta = LatConR.gamma = 90. * PIover180;
  MinLatticeTr2 = 1.;
  MaxLatticeTr2 = sqrt(3.);
  LambdaName   = "";
  LambdaLength = 1.;
   WyckoffList = NULL;
  nWyckoffList = 0;
  MaxPotentialAtoms = 0;
  MaxRecycledAtoms = 0;
  MaxPeaksFwSearch = 0;
  MaxPeaksFwFragmentSearch = 0;
  MaxRawPeaks = 0;
  NextPeakMx = NULL;
  MinDistance2Mx = NULL;
  PeakSearchLevel = 3;
  eDensityCutOff.Value = 0.;
  eDensityCutOff.Percent = 1;
  MinPfI = 17;
  eD_PeaksSortElement = PSE_Grid_eD;
  nAtomType = 0;
  uAtomType = 0;
   AtomType = NULL;
  MinConsecutiveNodes = -1;
  IndicateFw = 0;
  FwSearchMethod = FwSM_FwTracking;
  MinLoopSize =  3;
  MaxLoopSize = 24;
  EvenLoopSizesOnly = 0;
  nNodeTypes = 0;
   NodeTypes = NULL;
  NCNmin = 0;
  NCNmax = 0;
  MinSymNodes = 0;
  MaxSymNodes = -1;
  Check3DimConnectivity = 1;
  MinNodeDist2            = Square(2.8);
  MaxNodeDist2            = Square(3.4);
  IdealT_NodeDist2        = Square(3.1);
  IdealTetrahedralEdge2   = -1.;
  IdealTetrahedralVol2    = -1.;
  IdealTetrahedralVolFrac = .2;
  ModeCheckTetrahedra = 1;
  nSite = 0;
   Site = NULL;
  nSF_Tables = 0;
   SF_Tables = NULL;
  Nx = 24;
  Ny = 24;
  Nz = 24;
  CatchDistance2 = Square(.5);
  FobsScale = 1.;
  SigmaCutOff = 0.;
  AbsentRedivisionLimit.Value = 0.;
  AbsentRedivisionLimit.Type  = ARLT_FWHM;
  AbsentRedivisionLimit.Frac_I_ignored = 0.;
  AbsentRedivisionLimit.Frac_I_w_moved = 0.;
  OverlapFactor = 0.;
  OverlapAction = OvlAct_NoAction;
  GenerateFWHM.Valid = 0;
  GenerateFWHM.U = 0.;
  GenerateFWHM.V = 0.;
  GenerateFWHM.W = 0.;
  ReflectionUsage.Value = 100.;
  ReflectionUsage.Percent = 1;
  SumFmrg = 0.;
  MaxFmrg = 0.;
  BogFmrg = 0.;
  eDminRho = 0.;
  eDmaxRho = 0.;
  IRho3dVscale = 0.;
  IRho3dV = 0.;
  IRius = 0.;
  tRius = 0.;
  eD_PeakList = NULL;
  NeD_PeakList = 0;
  nList__RawSymEquiv = 0;
  List_xRawSymEquiv = NULL;
  List_cRawSymEquiv = NULL;
  nCodeTransl = 0;
  CodeTransl = NULL;

  ProfileStart = 0.;
  ProfileEnd   = 0.;
  ProfileStep  = 0.01;
  ProfileGenEnd = 0.;
  ProfilePOLRA = 1.;
  ProfileFWHM.Valid = 0;
  ProfileFWHM.U = 0.;
  ProfileFWHM.V = 0.;
  ProfileFWHM.W = 0.;
  ProfileAsym.Valid = 0;
  ProfileAsym.a1 = 0.;
  ProfileAsym.a2 = 0.;
  ProfileAsym.a3 = 0.;
  ProfilePeakShape = PPS_PseudoVoigt;
  PseudoVoigtPeakRange = 10.;
  PseudoVoigtFracLorentz = .5;
  ProfileBackground = 0.;
  ProfileReferenceRefl.Mode = PRRM_Free;
  ProfileReferenceRefl.h = 0;
  ProfileReferenceRefl.k = 0;
  ProfileReferenceRefl.l = 0;
  ProfileReferenceRefl.TTheta = 0.;
  ProfileReferenceRefl.Index = -1;
  ProfileReferenceMax = 0.;

  FT_TimeStart.second = -1;
  FT_SumTicks = 0;
  nFourierTransform = 0;
  nCatchRawPeak = 0;
  CountInterSectionII = 0;
  CountInterSectionIP = 0;
  CountInterSectionPI = 0;
  CountInterSectionPP = 0;
  IllPeakInterpolation = 0;
  SucPeakInterpolation = 0;
  FwSearchTimeStart.second = -1;
  FwSearchSumTicks = 0;
  iFramework = 0;
  nBadTetrahedraFW = 0;
  nNo3DimConFW = 0;
  nSmallLoopsFW = 0;
  nRejectedOddLoopsFW = 0;

  nNormalizeSurface = 0;
  vNormalizeSurface[0] = 0.;
  vNormalizeSurface[1] = 0.;
  vNormalizeSurface[2] = 0.;

  CurrCorrelationCoefficient = 0.;
  CurrPhaseCode_nCallRanmar = 0;

  ModeScatteringFactorTable = 0; // 0 = WK95  2 = IT4322  3 = IT4323 
}


void p_rogerror(const char *message, const char *file, const int line)
{
  Fflush(stdout);
  Fprintf(stderr, "\n%s(%s:%d)", progn, file, line);
  if (message) Fprintf(stderr, ": %s", message);
  putc('\n', stderr);
  AppExit(1);
}


void I_nternalError(const char *message, const char *file, const int line)
{
  Fflush(stdout);
  Fprintf(stderr, "\n%s(%s:%d): Internal Error", progn, file, line);
  if (message) Fprintf(stderr, ": %s", message);
  putc('\n', stderr);
  AppExit(3);
}


void N_otEnoughCore(const char *file, const int line)
{
  p_rogerror("Not enough core", file, line);
}


void I_llegalLine(const char *fnin, long lcount,
                  const char *file, const int line)
{
  Fflush(stdout);
  Fprintf(stderr, "%s(%s:%d): Illegal line #%ld in %s\n",
    progn, file, line, lcount, fnin);
  exit(1);
}


void *App_Malloc(int n, size_m size, char *var_name)
{
  size_m  nsize;
  void    *ptr;


      nsize = n * size;
  if (nsize != 0)
  {
    ptr = malloc(nsize);
    if (ptr != NULL) MemoryUsed += nsize;
    if (MaxMemoryUsed < MemoryUsed) MaxMemoryUsed = MemoryUsed;
  }

  if (Debug0)
    Fprintf(stdout, "# App_Malloc: %s(%d * %ld): +%ld => %ld %ld\n",
      var_name, n, (long) size,
      (long) nsize, (long) MemoryUsed, (long) MaxMemoryUsed);

  if (nsize == 0)
#ifdef NO_MALLOC_0_BYTES
    progerror("App_Malloc: attempt to allocate 0 bytes");
#else
    ptr = NULL;
#endif

  return ptr;
}


void *App_Realloc(void *ptr, int nn, int n, size_m size,
                  char *nptr_name, char *ptr_name)
{
  size_m  nnsize, nsize;
  void    *nptr;


   nsize =  n * size;
  nnsize = nn * size;

  if (nnsize != 0)
  {
    nptr = realloc(ptr, nnsize);
    if (nptr != NULL) MemoryUsed += nnsize - nsize;
    if (MaxMemoryUsed < MemoryUsed) MaxMemoryUsed = MemoryUsed;
  }

  if (Debug0)
  {
    Fprintf(stdout,
      "# App_Realloc: %s(%d * %ld) - %s(%d * %ld): %ld - %ld = %ld",
      nptr_name, nn, (long) size,
       ptr_name,  n, (long) size,
      (long) nnsize, (long) nsize, (long)(nnsize - nsize));

    Fprintf(stdout, " => %ld %ld\n",
      (long) MemoryUsed, (long) MaxMemoryUsed);

    if (nptr != ptr)
      Fprintf(stdout, "# App_Realloc: %s != %s\n", nptr_name, ptr_name);
  }

  if (nnsize == 0)
#ifdef NO_MALLOC_0_BYTES
    progerror("App_Realloc: attempt to allocate 0 bytes");
#else
    nptr = NULL;
#endif

  return nptr;
}


void App_Free(void *ptr, size_m n, size_m size, char *var_name)
{
  size_m  nsize;

  free(ptr);
                nsize = n * size;
  MemoryUsed -= nsize;

  if (Debug0)
    Fprintf(stdout, "# App_Free:   %s(%ld * %ld): -%ld => %ld %ld\n",
      var_name, (long) n, (long) size,
      (long) nsize, (long) MemoryUsed, (long) MaxMemoryUsed);
}


char *App_Strdup(char *s, char *var_name)
{
  size_m  size;
  char    *sdup;


  size = (strlen(s) + 1);

  sdup = (char *) malloc(size);
  if (sdup != NULL)
  {
    MemoryUsed += size;
    if (MaxMemoryUsed < MemoryUsed) MaxMemoryUsed = MemoryUsed;
    (void) strcpy(sdup, s);
  }

  if (Debug0)
    Fprintf(stdout, "# App_Strdup: %s(%d * %ld): +%ld => %ld %ld\n",
      var_name, 1, (long) size,
      (long) size, (long) MemoryUsed, (long) MaxMemoryUsed);

  return sdup;
}


void AppExit(int status)
{
  if (status == 2) Fprintf(stdout, "\n\n");

  if      (QuitProgram == 1)
    Fprintf(stdout, "# Program stopped due to signal\n\n");
  else if (QuitProgram == 2)
    Fprintf(stdout, "# Program stopped due to signal file\n\n");

  if (ReportOnExit == 1)
  {
    Fprintf(stdout, "#Trials               = %d\n",  nTrials);
    Fprintf(stdout, "#ConvergedPhases      = %d\n",  nConvergedPhases);
    Fprintf(stdout, "#NotConvergedPhases   = %d\n",  nNotConvergedPhases);

    Fprintf(stdout, "#FourierTransforms    = %d\n",  nFourierTransform);
    Fprintf(stdout, "#IllPeakInterpolation = %ld\n", IllPeakInterpolation);
    Fprintf(stdout, "#SucPeakInterpolation = %ld\n", SucPeakInterpolation);
  }

  if (nCatchRawPeak ||
      CountInterSectionII ||
      CountInterSectionIP ||
      CountInterSectionPI ||
      CountInterSectionPP ||
      FwSearchSumTicks || FwSearchTimeStart.second != -1)
  {
    Fprintf(stdout, "#CatchRawPeak         = %ld\n", nCatchRawPeak);
    Fprintf(stdout, "#InterSectionII       = %ld\n", CountInterSectionII);
    Fprintf(stdout, "#InterSectionIP       = %ld\n", CountInterSectionIP);
    Fprintf(stdout, "#InterSectionPI       = %ld\n", CountInterSectionPI);
    Fprintf(stdout, "#InterSectionPP       = %ld\n", CountInterSectionPP);
  }

  if (FwSearchSumTicks || FwSearchTimeStart.second != -1)
  {
    Fprintf(stdout, "#Frameworks           = %d\n",  iFramework);
    Fprintf(stdout, "#BadTetrahedraFW      = %d\n",  nBadTetrahedraFW);
    Fprintf(stdout, "#No3DimConFW          = %d\n",  nNo3DimConFW);
    Fprintf(stdout, "#SmallLoopsFW         = %d\n",  nSmallLoopsFW);
    Fprintf(stdout, "#RejectedOddLoopsFW   = %d\n",  nRejectedOddLoopsFW);
  }

  Fprintf(stdout, "#Warnings = %ld\n", (long) CountWarnings);
  Fprintf(stdout, "# MemoryUsed = %ld %ld\n",
                  (long) MemoryUsed, (long) MaxMemoryUsed);
  PrintTicks(stdout, NULL, "# Time ", "  Exit\n");

  if (FT_SumTicks || FT_TimeStart.second != -1)
  {
    (void) GetTicks(&FT_TimeEnd);

    if (FT_TimeStart.second != -1)
    {
      FT_SumTicks +=   FT_TimeEnd.cpu_user
                     - FT_TimeStart.cpu_user;
      FT_SumTicks +=   FT_TimeEnd.cpu_system
                     - FT_TimeStart.cpu_system;
    }

    Fprintf(stdout, "# Time Fourier Transform = %.2f (%.6g s / transform)\n",
      (double) FT_SumTicks / FT_TimeEnd.ticks_per_second,
      (double) FT_SumTicks / FT_TimeEnd.ticks_per_second
                           / (double) nFourierTransform);
  }

  if (FwSearchSumTicks || FwSearchTimeStart.second != -1)
  {
    (void) GetTicks(&FwSearchTimeEnd);

    if (FwSearchTimeStart.second != -1)
    {
      FwSearchSumTicks +=   FwSearchTimeEnd.cpu_user
                          - FwSearchTimeStart.cpu_user;
      FwSearchSumTicks +=   FwSearchTimeEnd.cpu_system
                          - FwSearchTimeStart.cpu_system;
    }

    Fprintf(stdout, "# Time Framework Search = %.2f\n",
      (double) FwSearchSumTicks / FwSearchTimeEnd.ticks_per_second);
  }

#ifdef PROVOKE_SEG_FAULT
  if (status == 3) /* provoke a segmentation fault */
  {
    int  *ip = NULL;
    Fflush(stdout);
    Fflush(stderr);
    status = *ip;
  }
#endif

  exit(status);
}


static int NT_SortFunction(const T_NodeType *a,
                           const T_NodeType *b)
{
  if (a->CN < b->CN) return -1;
  if (a->CN > b->CN) return  1;
  return 0;
}


static void SetMaxRawPeaks(void)
{
  int  iAT, nApUC;


  nApUC = 0;

  for (iAT = 0; iAT < nAtomType; iAT++)
    if (AtomType[iAT].On)
      nApUC += AtomType[iAT].nPerUnitCell;

  if (MaxPotentialAtoms < MaxRecycledAtoms)
      MaxPotentialAtoms = MaxRecycledAtoms;

  if (MaxPotentialAtoms == 0)
      MaxPotentialAtoms = nApUC * 6 / 5;

  if (MaxRecycledAtoms == 0)
      MaxRecycledAtoms = nApUC;

  if (MaxPeaksFwSearch == 0)
      MaxPeaksFwSearch = SpgrInfo.OrderL * 50;

  if (MaxPeaksFwFragmentSearch == 0)
      MaxPeaksFwFragmentSearch = SpgrInfo.OrderL * 20;

      MaxRawPeaks = MaxPotentialAtoms;
  if (MaxRawPeaks < MaxRecycledAtoms)
      MaxRawPeaks = MaxRecycledAtoms;
  if (MaxRawPeaks < MaxPeaksFwSearch)
      MaxRawPeaks = MaxPeaksFwSearch;
  if (MaxRawPeaks < MaxPeaksFwFragmentSearch)
      MaxRawPeaks = MaxPeaksFwFragmentSearch;
}


static void ReworkNodeTypes(void)
{
  int  i, unlimited;


  if (nNodeTypes == 0)
  {                        nNodeTypes = 1;
    CheckMalloc(NodeTypes, nNodeTypes);

    NodeTypes[0].CN      = 4;
    NodeTypes[0].MaxNpAU = 0;

    for (i = -6; i <= 6; i++)
      NodeTypes[0].NoSPO[i + 6] = TetrahedraNoSPO[i + 6];
  }


  NCNmin = NCNmax = NodeTypes[0].CN;
  unlimited = 0;

  for (i = 0; i < nNodeTypes; i++)
  {
    if (NCNmin > NodeTypes[i].CN)
        NCNmin = NodeTypes[i].CN;

    if (NCNmax < NodeTypes[i].CN)
        NCNmax = NodeTypes[i].CN;

    if (NodeTypes[i].MaxNpAU == 0)
      unlimited = 1;
  }

  if (unlimited == 0)
  {
    Fprintf(stdout, "\n# WARNING: No NodeType with an unlimited"
                    " number per asymmetric unit\n\n");
    CountWarnings++;
  }

  if (NCNmax > MaxNCNmax)
    InternalError("NCNmax out of range");

  if (nNodeTypes > 1)
    qsort((void *) NodeTypes, nNodeTypes, sizeof (*NodeTypes),
          (SortFunction) NT_SortFunction);
}


static void OverlapStatistics(void)
{
  int        iFR, nTotal, nOverl;
  T_FobsRaw  *FR;
  Fprec      d, MFtotal, MFoverl, RatioMF;


  if (OverlapFactor == 0.) return;

  Fprintf(stdout, ">Begin OverlapStatistics\n");

  MFtotal = 0.;
  MFoverl = 0.;
  nTotal  = 0;
  nOverl  = 0;

  for (iFR = 0, FR = FobsRaw; iFR < nFobsRaw; iFR++, FR++)
  {
    if (FR->Status != FRS_F000)
    {
      MFtotal += FR->M * FR->Fmrg;
       nTotal++;

      if (FR->Overlap) {
        MFoverl += FR->M * FR->Fmrg;
         nOverl++;
      }

      if (FR->Q > 0.) d = 1. / AppSqrt(FR->Q);
      else            d = 0.;

      if (MFtotal) RatioMF = 100. * MFoverl / MFtotal;
      else         RatioMF = 0.;

      Fprintf(stdout, "%.6g %.6g %.6g\n",
        d, 100. * nOverl / nTotal, RatioMF);
    }
  }

  Fprintf(stdout, ">End OverlapStatistics\n\n");
}


static void CompleteAtomTypeSF_Info(void)
{
  int         iAT;
  T_AtomType  *AT;


  AT = AtomType;

  for (iAT = 0; iAT < nAtomType; iAT++, AT++)
  {
    if (CompleteSF_Info(&AT->SF_Info, 1, 1) != 0)
    {
      char  buf[128];

      Sprintf(buf,
        "AtomType  %s  %-4.60s: Unknown Scattering Factor Label",
        EchoAtomTypeClass(NULL, AT->Class), AT->SF_Info.Lbl);

      progerror(buf);
    }

    AppFree((char *) AT->SF_Info.Lbl, strlen(AT->SF_Info.Lbl) + 1);

    if (AT->SF_Info.SFT)
      AT->SF_Info.Lbl = AT->SF_Info.SFT->Label;
    else
      AT->SF_Info.Lbl = AT->SF_Info.CAA->Label;
  }
}


static void CompleteInputChemistrySF_Info(void)
{
  int               i, n;
  T_SF_Info         *SFI;
  int               iIC;
  T_InputChemistry  *IC;


  IC = InputChemistry;

  for (iIC = 0; iIC < nInputChemistry; iIC++, IC++)
  {
    n = 2; if (IC->Type == ICT_AddBondAtom) n++;

    for (i = 0; i < n; i++)
    {
      SFI = &IC->AtomID[i].SF_Info;

      if (CompleteSF_Info(SFI, 1, 0) == 0)
      {
        AppFree((char *) SFI->Lbl, strlen(SFI->Lbl) + 1);

        if (SFI->SFT)
          SFI->Lbl = SFI->SFT->Label;
        else
          SFI->Lbl = SFI->CAA->Label;
      }
    }
  }
}


static int AT_SortFunction(const T_AtomType *a, const T_AtomType *b)
{
  Fprec       a_occ_charge, b_occ_charge;


  if (a->On != 0 && b->On == 0) return -1;
  if (a->On == 0 && b->On != 0) return  1;

  a_occ_charge = a->SF_Info.f_stol_0 * a->OccDefault;
  b_occ_charge = b->SF_Info.f_stol_0 * b->OccDefault;

  if (a_occ_charge > b_occ_charge) return -1;
  if (a_occ_charge < b_occ_charge) return  1;

  if (a->Class == ATC_Node       && b->Class != ATC_Node      ) return -1;
  if (a->Class != ATC_Node       && b->Class == ATC_Node      ) return  1;

  if (a->Class == ATC_NodeBridge && b->Class != ATC_NodeBridge) return -1;
  if (a->Class != ATC_NodeBridge && b->Class == ATC_NodeBridge) return  1;

  return 0;
}


static int HarmonizeSgLatCon(T_SgInfo *SgInfo, T_LatticeConstants *lc, int np)
{
  switch(SgInfo->XtalSystem)
  {
    case XS_Triclinic:
      if (np != 6) goto IllUnitCell;
      break;
    case XS_Monoclinic:
      if (np != 4 && np != 6) goto IllUnitCell;
      switch (SgInfo->UniqueRefAxis)
      {
        case 'x': lc->beta  = lc->gamma = 90. * PIover180; break;
        case 'y': if (np != 6) lc->beta  = lc->alpha;
                  lc->alpha = lc->gamma = 90. * PIover180; break;
        case 'z': if (np != 6) lc->gamma = lc->alpha;
                  lc->alpha = lc->beta  = 90. * PIover180; break;
        default:
          goto IntErr;
      }
      break;
    case XS_Orthorhombic:
      if (np != 3 && np != 6) goto IllUnitCell;
      lc->alpha = lc->beta = lc->gamma = 90. * PIover180;
      break;
    case XS_Tetragonal:
      if (np != 2 && np != 6) goto IllUnitCell;
      switch (SgInfo->UniqueRefAxis)
      {
        case 'x': lc->c = lc->b; break;
        case 'y': lc->c = lc->a; break;
        case 'z': if (np != 6) lc->c = lc->b;
                  lc->b = lc->a; break;
        default:
          goto IntErr;
      }
      lc->alpha = lc->beta = lc->gamma = 90. * PIover180;
      break;
    case XS_Trigonal:
      if (np != 2 && np != 6) goto IllUnitCell;
      if (SgInfo->UniqueDirCode == '*')
      {
        if (np != 6) lc->alpha = lc->b * PIover180;
        lc->c = lc->b = lc->a;
        lc->gamma = lc->beta = lc->alpha;
        break;
      }
    case XS_Hexagonal:
      if (np != 2 && np != 6) goto IllUnitCell;
      switch (SgInfo->UniqueRefAxis)
      {
        case 'x': lc->c = lc->b;
                  lc->alpha = 120. * PIover180;
                  lc->beta  = lc->gamma = 90. * PIover180; break;
        case 'y': lc->c = lc->a;
                  lc->beta  = 120. * PIover180;
                  lc->alpha = lc->gamma = 90. * PIover180; break;
        case 'z': if (np != 6) lc->c = lc->b;
                  lc->b = lc->a;
                  lc->gamma = 120. * PIover180;
                  lc->alpha = lc->beta  = 90. * PIover180; break;
        default:
          goto IntErr;
      }
      break;
    case XS_Cubic:
      if (np != 1 && np != 6) goto IllUnitCell;
      lc->c = lc->b = lc->a;
      lc->alpha = lc->beta = lc->gamma = 90. * PIover180;
      break;
    default:
      goto IntErr;
  }

  return  0;

  IntErr: SetSgError("Internal Error: HarmonizeSgLatCon()");
  return -1;

  IllUnitCell: SetSgError("Error: Illegal UnitCell or SpaceGroup");
  return -1;
}


static void DoGenerateFWHM(T_FobsRaw *FR, int nFR)
{
  Fprec  tanT, FWHM2;


  if (! GenerateFWHM.Valid)
    return;

  while (nFR--)
  {
    if (FR->Status != FRS_F000)
    {
      tanT = AppTan(TwoThetaDeg(FR->Q) * .5 * PIover180);

      FWHM2 = GenerateFWHM.U + tanT * (GenerateFWHM.V + GenerateFWHM.W * tanT);

      if (FWHM2 < 0.)
        progerror("Illegal GenerateFWHM parameters");

      FR->FWHM = AppSqrt(FWHM2);
    }

    FR++;
  }
}


static void AllAsyPoints(void)
{
  int          ix, iy, iz, iD, iSite;
  char         buf[256], *cp;
  T_PeakFlags  *PeakFlags, *PF;
  Fprec        *Surface;


  for (iSite = 0; iSite < nSite; iSite++)
  {
        cp = Site[iSite].Label;
    if (cp) AppFree(cp, strlen(cp) + 1);
        cp = (char *) Site[iSite].SF_Info.Lbl;
    if (cp) AppFree(cp, strlen(cp) + 1);
  }

  AppFree(Site, nSite);

  Site = NULL;
  nSite = 0;

  Surface = NULL;

  if (F_SurfaceFileName) {
    CheckMalloc(Surface, Nx * Ny * Nz);
    LoadSurface(Surface, Nx, Ny, Nz, NULL);
  }

  CheckMalloc(PeakFlags, Nx * Ny * Nz);
  DS_MarkEquiv(&SpgrInfo, PeakFlags, Nx, Ny, Nz);

  PF = PeakFlags;

  for (iD = 0; iD < Nx * Ny * Nz; iD++, PF++)
    if ((*PF) & PF_HighBit && (Surface == NULL || Surface[iD] > 0.))
      nSite++;

  if (nSite != 0)
  {
    CheckMalloc(Site, nSite);

    iSite = 0;

    iD = 0;
    PF = PeakFlags;

    for (ix = 0; ix < Nx; ix++)
    for (iy = 0; iy < Ny; iy++)
    for (iz = 0; iz < Nz; iz++, iD++, PF++)
    {
      if ((*PF) & PF_HighBit && (Surface == NULL || Surface[iD] > 0.))
      {
          Sprintf(buf, "N%5.5d", iSite);
        AppStrdup(buf, Site[iSite].Label);
        if (NULL    == Site[iSite].Label) NotEnoughCore();
        Site[iSite].SF_Info.Lbl      = NULL;
        Site[iSite].SF_Info.SFT      = NULL;
        Site[iSite].SF_Info.CAA      = NULL;
        Site[iSite].SF_Info.f_stol_0 = -1.;
        Site[iSite].x                = (Fprec) ix / Nx;
        Site[iSite].y                = (Fprec) iy / Ny;
        Site[iSite].z                = (Fprec) iz / Nz;
        Site[iSite].Occ              = 0.;
        Site[iSite].F_Occ            = 0;
        Site[iSite].Uiso             = 0.;
        Site[iSite].F_Uiso           = 0;
             iSite++;
      }
    }

    if (iSite != nSite)
      InternalError("Corrupt PeakFlags");
  }

  AppFree(PeakFlags, Nx * Ny * Nz);
  if (Surface) AppFree(Surface, Nx * Ny * Nz);
}


void DoRandomInitialization(void)
{
  long  ij, kl;


  ij = labs(RandomInitialization);

  kl = (long)(1.234 * (double) ij);
  if (ij < 0) ij = -ij;
  if (kl < 0) kl = -kl;
  ij %= 31328;
  kl %= 30081;

  if (rmarin((int) ij, (int) kl) != 0)
    InternalError("call rmarin()");
}


static void PrintMin_dMax_hkl(const char *Lbl, T_MaxQ_MaxH *MQH)
{

  Fprintf(stdout, ">Begin %s\n", Lbl);

  Fprintf(stdout, "Min d       = ");
    if (MQH->Q == 0.)
      Fprintf(stdout, "Infinity\n");
    else
      Fprintf(stdout, "%.6g\n", 1. / AppSqrt(MQH->Q));
  Fprintf(stdout, "Max sinTovL = %.6g\n",
    AppSqrt(MQH->Q) * .5);
  if (LambdaLength > 0.)
    Fprintf(stdout, "Max 2-Theta = %.6g\n",
      TwoThetaDeg(MQH->Q));

  Fprintf(stdout, "Max |h| = %3d\n", MQH->h);
  Fprintf(stdout, "Max |k| = %3d\n", MQH->k);
  Fprintf(stdout, "Max |l| = %3d\n", MQH->l);

  Fprintf(stdout, ">End %s\n\n", Lbl);
}


static void SetMaxQ_MaxH(void)
{
  int        iFR, i;
  T_FobsRaw  *FR;
  T_Eq_hkl   Eq_hkl;
  int        Eq_h, Eq_k, Eq_l;


  Ac_MaxQ_MaxH.Q = 0.;
  Ac_MaxQ_MaxH.h = Ac_MaxQ_MaxH.k = Ac_MaxQ_MaxH.l = 0;
  AS_MaxQ_MaxH.Q = 0.;
  AS_MaxQ_MaxH.h = AS_MaxQ_MaxH.k = AS_MaxQ_MaxH.l = 0;

  FR = FobsRaw;

  for (iFR = 0; iFR < nFobsRaw; iFR++, FR++)
  {
    if (FR->Status != FRS_Active && FR->Status != FRS_Sleeping) continue;

    if (AS_MaxQ_MaxH.Q < FR->Q)
        AS_MaxQ_MaxH.Q = FR->Q;

    if (   FR->Status == FRS_Active
        && Ac_MaxQ_MaxH.Q < FR->Q)
           Ac_MaxQ_MaxH.Q = FR->Q;

    (void) BuildEq_hkl(&SpgrInfo, 1, &Eq_hkl, FR->h, FR->k, FR->l);

    for (i = 0; i < Eq_hkl.N; i++)
    {
      Eq_h = abs(Eq_hkl.h[i]);
      Eq_k = abs(Eq_hkl.k[i]);
      Eq_l = abs(Eq_hkl.l[i]);

      if (AS_MaxQ_MaxH.h < Eq_h) AS_MaxQ_MaxH.h = Eq_h;
      if (AS_MaxQ_MaxH.k < Eq_k) AS_MaxQ_MaxH.k = Eq_k;
      if (AS_MaxQ_MaxH.l < Eq_l) AS_MaxQ_MaxH.l = Eq_l;

      if (FR->Status != FRS_Active) continue;

      if (Ac_MaxQ_MaxH.h < Eq_h) Ac_MaxQ_MaxH.h = Eq_h;
      if (Ac_MaxQ_MaxH.k < Eq_k) Ac_MaxQ_MaxH.k = Eq_k;
      if (Ac_MaxQ_MaxH.l < Eq_l) Ac_MaxQ_MaxH.l = Eq_l;
    }
  }
}


#define usage() u_sage(__LINE__)

static void u_sage(int source_code_line)
{
  Fprintf(stderr,
    "usage(%d): %s\n"
    "  [-DebugAll]\n"
    "  [-SignalFile=FileName]\n"
    "  [-SiteFcal]\n"
    "  [-List_hkl]\n"
    "  [-Collection_hkl]\n"
    "  [-PowderStepScan]\n"
    "  [-StdPeak=FileName]\n"
    "  [-PrintFmrg]\n"
    "  [-SiteFrame]\n"
    "  [-SiteLabel]\n"
    "  [-AllSitesOn]\n"
    "  [-SitePhases]\n"
    "  [-AllPhaseCodes0]\n"
    "  [-Coseq[=Value]]\n"
    "  [-CoseqSelect=#]\n"
    "  [-CoseqUnitCell]\n"
    "  [-CoseqProtocol=FileName]\n"
    "  [-Surface=FileName]\n"
    "  [-CoseqSave=FileName]\n"
    "  [-CoseqSaveTime=Minutes]\n"
    "  [-CoseqKeep[=#]]\n"
    "  [-CoseqReload=FileName]\n"
    "  [-CoseqSplit=#x,#y,#z]\n"
    "  [-SetAllFmrg[=Value]]\n"
    "  [-Put_strudat[=FileName]]\n"
    "  [-PutAllPeaks]\n"
    "  [-ShowLargestFwFragment]\n"
    "  [-eD_Histogram[=FileName]]\n"
    "  [-Put_eDmap[=FileName]]\n"
    "  [-AllAsyPoints[=#]]\n"
    "  [-RandomSites[=#]]\n"
    "  [-TestMoreSpecial[=#]]\n"
    "  [-LoadDensity=FileName]\n"
    "  [-PhaseSets=FileName]\n"
    "  [-Surface=FileName]\n"
    "  focus.inp [Maximum#Trials]\n",
    source_code_line, progn);
  exit(1);
}


int main(int argc, char *argv[])
{
  int   i, n;
  char  *fnin;
  FILE  *fpin;
  char  *cp, xtrac;

  size_m  m;

  int  MaximumTrials, nTrialsPerPhaseSet;

  int     OverlapGroups;
  double  dbuf;

  int    F_PrintFmrg;
  int    F_SiteFcal, F_List_hkl, F_Collection_hkl, F_PowderStepScan;
  int    F_RandomSites, F_TestMoreSpecial;
  char  *F_StdPeakFileName;


  InitTicks(); /* initialize cputime routines */


  putc('#', stdout);
  for (i = 0; i < argc;)
    Fprintf(stdout, " %s", argv[i++]);
  putc('\n', stdout);

  Fprintf(stdout, "# Version %s  %s  %s\n",
    VersionID, __DATE__, __TIME__);
  Fprintf(stdout, "#    - Modified for electron diffraction\n");
  Fprintf(stdout, "# %s\n", SgInfoVersion());
  Fprintf(stdout, "# Host: %s\n", AppHostname());
  Fprintf(stdout, "# PID: %ld\n", AppPID());
  putc('\n', stdout);

  InitializeParameters();

  fnin = NULL;
  fpin = NULL;
  MaximumTrials = -1;

  F_SiteFcal = 0;
  F_List_hkl = 0;
  F_Collection_hkl = 0;
  F_PowderStepScan = 0;
  F_StdPeakFileName = NULL;
  F_PrintFmrg = 0;

  F_SiteFrame = 0;
  F_SiteLabel = 0;
  F_AllSitesOn = 0;
  F_SitePhases = 0;
  F_SignalFileName = NULL;
  F_AllPhaseCodes0 = 0;
  F_SetAllFmrg = 0;
  F_SetAllFmrgValue = 0.;
  F_CoseqValue = 10;
  F_CoseqSelect = -1;
  F_CoseqUnitCell = 0;
  F_CoseqProtocolFileName = NULL;
  F_CoseqSaveFileName = NULL;
  F_CoseqSaveTime = -1;
  F_CoseqKeep = -1;
  F_CoseqReloadFileName = NULL;
  F_CoseqSplit = NULL;
  F_nCoseqSplit = 0;
  F_CoseqSplit_x = -1;
  F_CoseqSplit_y = -1;
  F_CoseqSplit_z = -1;
  F_Put_strudat = 0;
  F_Put_strudatFileName = NULL;
  F_PutAllPeaks = 0;
  F_ShowLargestFwFragment = 0;
  F_eD_Histogram = 0;
  F_eD_HistogramFileName = NULL;
  F_Put_eDmap = 0;
  F_Put_eDmapFileName = NULL;
  F_nAllAsyPoints = -1;
   F_AllAsyPoints = NULL;
  F_RandomSites = -1;
  F_TestMoreSpecial = -1;
  F_LoadDensityFileName = NULL;
  F_PhaseSetsFileName = NULL;
  F_SurfaceFileName = NULL;

  for (i = 1; i < argc; i++)
  {
    if      (str_icmp(argv[i], "-DebugAll") == 0)
      Debug0 = 1;
    else if (str_icmp(argv[i], "-SiteFcal") == 0)
      F_SiteFcal = 1;
    else if (str_icmp(argv[i], "-List_hkl") == 0)
      F_List_hkl = 1;
    else if (str_icmp(argv[i], "-Collection_hkl") == 0)
      F_Collection_hkl = 1;
    else if (str_icmp(argv[i], "-PowderStepScan") == 0)
      F_PowderStepScan = 1;
    else if (str_ibegin(argv[i], "-StdPeak=") == 0)
    {
      if (F_StdPeakFileName) usage();
                 cp = argv[i] + 9;
      if (strlen(cp) == 0) usage();
      AppStrdup(cp, F_StdPeakFileName);
                if (F_StdPeakFileName == NULL) NotEnoughCore();
    }
    else if (str_icmp(argv[i], "-PrintFmrg") == 0)
      F_PrintFmrg = 1;
    else if (str_icmp(argv[i], "-SitePhases") == 0)
      F_SitePhases = 1;
    else if (str_ibegin(argv[i], "-SignalFile") == 0)
    {
      if (F_SignalFileName) usage();
      cp = argv[i] + 11;
      if (*cp++ != '=') usage();
      if (strlen(cp) == 0) usage();
      AppStrdup(cp, F_SignalFileName);
      if (F_SignalFileName == NULL) NotEnoughCore();
          fpin = fopen(F_SignalFileName, "r");
      if (fpin != NULL)
      {
        Fprintf(stderr,
          "%s: signal file %s already existing: remove first\n",
          progn, F_SignalFileName);
        exit(1);
      }
    }
    else if (str_icmp(argv[i], "-AllPhaseCodes0") == 0)
      F_AllPhaseCodes0 = 1;
    else if (str_ibegin(argv[i], "-SetAllFmrg") == 0)
    {
           cp = argv[i] + 11;
      if (*cp)
      {
        if (*cp++ != '=') usage();
        n = sscanf(cp, "%lf", &dbuf);
        if (n != 1) usage();
        F_SetAllFmrgValue = (Fprec) dbuf;
      }
      else
        F_SetAllFmrgValue = 1.;
      F_SetAllFmrg = 1;
    }
    else if (str_icmp(argv[i], "-SiteFrame") == 0)
      F_SiteFrame = 1;
    else if (str_icmp(argv[i], "-SiteLabel") == 0)
      F_SiteLabel = 1;
    else if (str_icmp(argv[i], "-AllSitesOn") == 0)
      F_AllSitesOn = 1;
    else if (str_ibegin(argv[i], "-CoseqProtocol=") == 0)
    {
      if (F_CoseqProtocolFileName) usage();
                 cp = argv[i] + 15;
      if (strlen(cp) == 0) usage();
      AppStrdup(cp, F_CoseqProtocolFileName);
      if (F_CoseqProtocolFileName == NULL) NotEnoughCore();
          fpin = fopen(F_CoseqProtocolFileName, "r");
      if (fpin != NULL)
      {
        Fprintf(stderr,
          "%s: coseq protocol file %s already existing: remove first\n",
          progn, F_CoseqProtocolFileName);
        exit(1);
      }
          fpin = fopen(F_CoseqProtocolFileName, "w");
      if (fpin == NULL)
      {
        Fprintf(stderr,
          "%s: Cannot write to file %s\n",
          progn, F_CoseqProtocolFileName);
        exit(1);
      }
      Fclose(fpin), fpin = NULL;
    }
    else if (str_ibegin(argv[i], "-CoseqSelect=") == 0)
    {
      if (F_CoseqSelect >= 0) usage();
                     cp = argv[i] + 13;
          n = sscanf(cp, "%d%c", &F_CoseqSelect, &xtrac);
      if (n != 1 || F_CoseqSelect < 0) usage();
    }
    else if (str_icmp(argv[i], "-CoseqUnitCell") == 0)
    {
      if (F_CoseqUnitCell != 0) usage();
          F_CoseqUnitCell = 1;
    }
    else if (str_ibegin(argv[i], "-CoseqSaveTime=") == 0)
    {
      if (F_CoseqSaveTime >= 0) usage();
                     cp = argv[i] + 15;
          n = sscanf(cp, "%d%c", &F_CoseqSaveTime, &xtrac);
      if (n != 1 || F_CoseqSaveTime < 0) usage();
    }
    else if (str_ibegin(argv[i], "-CoseqSave=") == 0)
    {
      if (F_CoseqSaveFileName) usage();
                 cp = argv[i] + 11;
      if (strlen(cp) == 0) usage();
      AppStrdup(cp, F_CoseqSaveFileName);
                if (F_CoseqSaveFileName == NULL) NotEnoughCore();
    }
    else if (str_ibegin(argv[i], "-CoseqKeep") == 0)
    {
      if (F_CoseqKeep >= 0) usage();
                cp = argv[i] + 10;
      if      (*cp == '\0')
        F_CoseqKeep = 0;
      else if (*cp == '=')
      {                cp++;
            n = sscanf(cp, "%d%c", &F_CoseqKeep, &xtrac);
        if (n != 1 || F_CoseqKeep < 0) usage();
      }
      else
        usage();
    }
    else if (str_ibegin(argv[i], "-CoseqReload=") == 0)
    {
      if (F_CoseqReloadFileName) usage();
                 cp = argv[i] + 13;
      if (strlen(cp) == 0) usage();
      AppStrdup(cp, F_CoseqReloadFileName);
                if (F_CoseqReloadFileName == NULL) NotEnoughCore();
    }
    else if (str_ibegin(argv[i], "-CoseqSplit=") == 0)
    {
      if (F_CoseqSplit)
        usage();

      F_nCoseqSplit = getnum(argv[i] + 12, &F_CoseqSplit, 3);

      for (n = 0; n < F_nCoseqSplit; n++)
        if (F_CoseqSplit[n] < 1)
          usage();

      if      (F_nCoseqSplit == 3)
      {
        F_CoseqSplit_x = F_CoseqSplit[0];
        F_CoseqSplit_y = F_CoseqSplit[1];
        F_CoseqSplit_z = F_CoseqSplit[2];
      }
      else if (F_nCoseqSplit != 1)
        usage();
    }
    else if (str_ibegin(argv[i], "-Coseq") == 0)
    {
           cp = argv[i] + 6;
      if (*cp)
      {
        if (*cp++ != '=') usage();
        n = sscanf(cp, "%d", &F_CoseqValue);
        if (n != 1 || F_CoseqValue < 0) usage();
      }
    }
    else if (str_ibegin(argv[i], "-Put_strudat") == 0)
    {
      cp = argv[i] + 12;
      if (*cp)
      {
        if (*cp++ != '=') usage();
        if (strlen(cp) == 0) usage();
        AppStrdup(cp, F_Put_strudatFileName);
        if (F_Put_strudatFileName == NULL) NotEnoughCore();
      }
      F_Put_strudat = 1;
    }
    else if (str_icmp(argv[i], "-PutAllPeaks") == 0)
      F_PutAllPeaks = 1;
    else if (str_icmp(argv[i], "-ShowLargestFwFragment") == 0)
      F_ShowLargestFwFragment = 1;
    else if (str_ibegin(argv[i], "-eD_Histogram") == 0)
    {
           cp = argv[i] + 13;
      if (*cp)
      {
        if (*cp++ != '=') usage();
        if (strlen(cp) == 0) usage();
        AppStrdup(cp, F_eD_HistogramFileName);
        if (F_eD_HistogramFileName == NULL) NotEnoughCore();
      }
      F_eD_Histogram = 1;
    }
    else if (str_ibegin(argv[i], "-Put_eDmap") == 0)
    {
      cp = argv[i] + 10;
      if (*cp)
      {
        if (*cp++ != '=') usage();
        if (strlen(cp) == 0) usage();
        AppStrdup(cp, F_Put_eDmapFileName);
        if (F_Put_eDmapFileName == NULL) NotEnoughCore();
      }
      F_Put_eDmap = 1;
    }
    else if (str_ibegin(argv[i], "-AllAsyPoints") == 0)
    {
      if (F_nAllAsyPoints >= 0) usage();

                cp = argv[i] + 13;
      if      (*cp == '\0')
        F_nAllAsyPoints = 0;
      else if (*cp == '=')
      {                          cp++;
        F_nAllAsyPoints = getnum(cp, &F_AllAsyPoints, 2);

        for (n = 0; n < F_nAllAsyPoints; n++)
          if (F_AllAsyPoints[n] < 0)
            usage();

        if (F_nAllAsyPoints < 1) usage();
      }
      else
        usage();
    }
    else if (str_ibegin(argv[i], "-RandomSites") == 0)
    {
      if (F_RandomSites >= 0) usage();
                cp = argv[i] + 12;
      if      (*cp == '\0')
        F_RandomSites = 0;
      else if (*cp == '=')
      {                cp++;
            n = sscanf(cp, "%d %c", &F_RandomSites, &xtrac);
        if (n != 1 || F_RandomSites < 0) usage();
      }
      else
        usage();
    }
    else if (str_ibegin(argv[i], "-TestMoreSpecial") == 0)
    {
      if (F_TestMoreSpecial >= 0) usage();
                cp = argv[i] + 16;
      if      (*cp == '\0')
        F_TestMoreSpecial = 0;
      else if (*cp == '=')
      {                cp++;
            n = sscanf(cp, "%d %c", &F_TestMoreSpecial, &xtrac);
        if (n != 1 || F_TestMoreSpecial < 0) usage();
      }
      else
        usage();
    }
    else if (str_ibegin(argv[i], "-LoadDensity=") == 0)
    {
      if (F_LoadDensityFileName) usage();
                 cp = argv[i] + 13;
      if (strlen(cp) == 0) usage();
      AppStrdup(cp, F_LoadDensityFileName);
      if (F_LoadDensityFileName == NULL) NotEnoughCore();
    }
    else if (str_ibegin(argv[i], "-PhaseSets=") == 0)
    {
      if (F_PhaseSetsFileName) usage();
                 cp = argv[i] + 11;
      if (strlen(cp) == 0) usage();
      AppStrdup(cp, F_PhaseSetsFileName);
      if (F_PhaseSetsFileName == NULL) NotEnoughCore();
    }
    else if (str_ibegin(argv[i], "-Surface=") == 0)
    {
      if (F_SurfaceFileName) usage();
                 cp = argv[i] + 9;
      if (strlen(cp) == 0) usage();
      AppStrdup(cp, F_SurfaceFileName);
      if (F_SurfaceFileName == NULL) NotEnoughCore();
    }
    else if (*argv[i] != '-')
    {
      if      (fnin == NULL)
        fnin = argv[i];
      else if (MaximumTrials < 0)
      {
            n = sscanf(argv[i], "%d%c", &MaximumTrials, &xtrac);
        if (n != 1 || MaximumTrials < 0) usage();
      }
      else
        usage();
    }
    else
      usage();
  }

  if (F_CoseqKeep < 0)
      F_CoseqKeep = 1;
  if (F_CoseqSaveTime < 0)
      F_CoseqSaveTime = 30;
  if (F_CoseqSplit == NULL)
  {
      F_nCoseqSplit = 3;
      F_CoseqSplit_x = 1;
      F_CoseqSplit_y = 1;
      F_CoseqSplit_z = 1;
  }

  if (fnin != NULL)
  {
        fpin = fopen(fnin, "r");
    if (fpin == NULL)
    {
      Fprintf(stderr, "%s: Can't open %s\n", progn, fnin);
      exit(1);
    }

    if (ReadCmdFile(fpin) != RCF_OK)
      IllegalLine(fnin, Ilinec);
  }

  if (FlagTimeRandomInitialization)
    RandomInitialization = GetTicks(NULL);

  if (MinLoopSize > MaxLoopSize)
    progerror("MinLoopSize > MaxLoopSize");

  if (nFeedBackCycles == 0 && argc < 2)
  {                             nFeedBackCycles = 10;
    CheckMalloc(FeedBackCycles, nFeedBackCycles);

    for (i = 0; i < 10; i++)
      FeedBackCycles[i] = 1;
  }

  if (TopTitle == NULL) TopTitle = fnin;

  if (SpgrInfo.nList == 0)
    if (BuildSpgrInfo(&SpgrInfo, "P1") != 0)
      InternalError("Can't set default space group");

  if (HarmonizeSgLatCon(&SpgrInfo, &LatConD, nInputLatConD) != 0)
    progerror(SgError);

  cp = "Illegal lattice constants";

  if (Lc2RLc(&LatConD, &LatConR) != 0) progerror(cp);
  if (SetupTrCrystCarte(&LatConD, &LatConR, 0, TMx_xc, 0) != 0) progerror(cp);
  if (SetupTrCrystCarte(&LatConD, &LatConR, 0, TMx_cx, 1) != 0) progerror(cp);
  CalcLatticeTr2(&LatConD, SpgrInfo.LatticeInfo,
                 &MinLatticeTr2, &MaxLatticeTr2);

  SetMaxRawPeaks();
  ReworkNodeTypes();

  SortSF_Tables();
  CompleteAtomTypeSF_Info();
  CompleteInputChemistrySF_Info();

  if (   F_SiteFcal
      || F_List_hkl
      || F_Collection_hkl
      || F_PowderStepScan
      || F_SitePhases)
    CompleteSite(Site, nSite, AtomType, nAtomType);

  CheckMalloc(fTrVector, SpgrInfo.LatticeInfo->nTrVector);
  CheckMalloc(List_fSeitzMx, SpgrInfo.nList);
  fSymOps(&SpgrInfo, List_fSeitzMx, fTrVector);

  IdealTetrahedron();

  m = MemoryUsed;

  if (F_SiteFcal)
    DoSiteFcal(0);

  if (F_List_hkl)
    DoSiteFcal(1);

  if (F_Collection_hkl || F_PowderStepScan)
    CollectionLists(F_Collection_hkl, F_PowderStepScan, F_StdPeakFileName);

  if (F_SiteFcal || F_List_hkl || F_Collection_hkl || F_PowderStepScan) {
    if (m != MemoryUsed)
      InternalError("Corrupt memory management");
    AppExit(0);
  }

                       SortFobsRaw(FobsRaw, nFobsRaw);
                 SetStatus1FobsRaw(FobsRaw, nFobsRaw);
                    DoGenerateFWHM(FobsRaw, nFobsRaw);
                          CalcFmrg(FobsRaw, nFobsRaw);
                    RedivideAbsent(FobsRaw, nFobsRaw);
  OverlapGroups = SetOverlapGroups(FobsRaw, nFobsRaw);
                      TreatOverlap(FobsRaw, nFobsRaw);
                 SetStatus2FobsRaw(FobsRaw, nFobsRaw);

  SetF000();

  if (nCodeTransl > 0)
  {
    CheckMalloc(CodeTransl, nCodeTransl);
    MaxFmrg = InitCodeTransl(CodeTransl, nCodeTransl, FobsRaw, nFobsRaw);
    BogFmrg = MaxFmrg * 1.e-5; /* ARBITRARY */
  }

  SetSleeping(CodeTransl, nCodeTransl);

  nSysAbsent = 0;
  nOmit = 0;
  nOffGrid = 0;
  nSymEquivFobs = 0;
  nActivePhase = 0;
  nSleeping = 0;

  for (i = 0; i < nFobsRaw; i++)
  {
    switch (FobsRaw[i].Status)
    {
      case FRS_F000:      nF000++;         break;
      case FRS_SysAbsent: nSysAbsent++;    break;
      case FRS_SymEquiv:  nSymEquivFobs++; break;
      case FRS_Omit:      nOmit++;         break;
      case FRS_OffGrid:   nOffGrid++;      break;
      case FRS_Active:    nActivePhase++;  break;
      case FRS_Sleeping:  nSleeping++;     break;
      case FRS_Undefined:
        InternalError("FobsRaw->Status Undefined");
        break;
      default:
        InternalError("Corrupt FobsRaw->Status");
        break;
    }
  }

  if (nActivePhase + nSleeping != nCodeTransl)
    InternalError("nActivePhase + nSleeping != nCodeTransl");

  if (argc < 2) {
    EchoInput(stdout, 0);
    exit(0);
  }

  EchoInput(stdout, 1);
  putc('\n', stdout);

  if (F_PrintFmrg)
    PrintFmrg(stdout);

  SetMaxQ_MaxH();

  if (Debug1) /* Debug: Print #Undefined... */
  {
    Fprintf(stdout, "#F000  %d\n", nF000);
    Fprintf(stdout, "#SysAbsent  %d\n", nSysAbsent);
    Fprintf(stdout, "#Omit  %d\n", nOmit);
    Fprintf(stdout, "#OffGrid  %d\n", nOffGrid);
    Fprintf(stdout, "#SymEquiv  %d\n", nSymEquivFobs);
    Fprintf(stdout, "#Active  %d\n", nActivePhase);
    Fprintf(stdout, "#Sleeping  %d\n", nSleeping);
    putc('\n', stdout);

    Fprintf(stdout, "# AbsentRedivision Frac I ignored = %.6g\n",
      AbsentRedivisionLimit.Frac_I_ignored);
    Fprintf(stdout, "# AbsentRedivision Frac I w moved = %.6g\n",
      AbsentRedivisionLimit.Frac_I_w_moved);
    putc('\n', stdout);

    Fprintf(stdout, "#OverlapGroups  %d\n", OverlapGroups);
    putc('\n', stdout);

    PrintMin_dMax_hkl("ActiveAndSleepingReflections", &AS_MaxQ_MaxH);
    PrintMin_dMax_hkl("ActiveReflections",            &Ac_MaxQ_MaxH);
  }

  if (Debug1) /* Debug: */ OverlapStatistics();

  if (Debug1) /* Debug: Print TMx_xc TMx_cx */
  {
    MxDump(TMx_xc, 3, 3, "TMx_xc");
    MxDump(TMx_cx, 3, 3, "TMx_cx");
    putc('\n', stdout);
  }

  if (SpgrInfo.HallSymbol[0])
  {
    Fprintf(stdout, "HallSymbol  %s\n", SpgrInfo.HallSymbol);
    ListSgInfo(&SpgrInfo, 1, 1, stdout);
    putc('\n', stdout);
  }

  BuildWyckoffList(&SpgrInfo);

  if (F_RandomSites >= 0 || F_TestMoreSpecial >= 0)
  {
    PrintWyckoffList();

#ifndef DisableAppSignal
    AppSignal();
#endif

    if (F_RandomSites >= 0)
      RandomSites(F_RandomSites);

    if (F_TestMoreSpecial >= 0)
      TestMoreSpecial(F_TestMoreSpecial);

    AppExit(0);
  }

  SetCanBeCN(&SpgrInfo);

  if (nAtomType)
  {
    qsort((void *) AtomType, nAtomType, sizeof (*AtomType),
          (SortFunction) AT_SortFunction);

    for (i = 0; i < nAtomType; i++)
      if (AtomType[i].On == 0)
        break;

    uAtomType = i;
  }
  else
    uAtomType = 0;

  if (MaxSymNodes < 0)
  {
    MaxSymNodes = 0;
    for (i = 0; i < uAtomType; i++)
      if (AtomType[i].Class == ATC_Node)
        MaxSymNodes += AtomType[i].nPerUnitCell;

    if (MaxSymNodes != 0)
      Fprintf(stdout, "# MaxSymNodes set to %d\n\n", MaxSymNodes);
  }

  if (F_SiteFrame == 0 || F_nAllAsyPoints >= 0 || nSite > 0)
    SetupXtalInternals();

  if (F_nAllAsyPoints >= 0)
         AllAsyPoints();

#ifndef DisableAppSignal
    AppSignal();
#endif

  if      (F_SiteFrame)
             SiteFrame();

  else if (nFobsRaw > 0)
  {
    if (F_SitePhases)
      DoSitePhases();

    if (F_SetAllFmrg != 0)
    {
      SumFmrg = 0.;

      for (i = 0; i < nActivePhase; i++)
      {
        CodeTransl[i].FobsRaw->Fmrg = F_SetAllFmrgValue;
        CodeTransl[i].M_Fmrg = CodeTransl[i].FobsRaw->M * F_SetAllFmrgValue;
        SumFmrg += CodeTransl[i].M_Fmrg;
      }

      Fprintf(stdout, "# SetAllFmrg = %.2f  =>  SumFmrg = %.2f\n\n",
        F_SetAllFmrgValue, SumFmrg);
    }

    if (nActivePhase > 0 && MaximumTrials != 0)
    {
      nTrialsPerPhaseSet = 0;

      if (F_PhaseSetsFileName)
      {
        nTrialsPerPhaseSet = 1;

        if (MaximumTrials > 0) {
          nTrialsPerPhaseSet = MaximumTrials;
          MaximumTrials = -1;
        }
      }

      if (MaximumTrials < 0)
      {
        if (     RandomBlindCalls != NULL
            &&  nRandomBlindCalls != 0
            && (nRandomBlindCalls > 1 || RandomBlindCalls[0] != 0))
          MaximumTrials = nRandomBlindCalls;
        else
          MaximumTrials = 0;
      }

      ReportOnExit = 1;
      TriaLoop(MaximumTrials, nTrialsPerPhaseSet);
    }
  }

  AppExit(0);

  return 0;
}
