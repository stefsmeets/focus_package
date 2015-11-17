#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "main.h"
#include "xtal.h"
#include "io.h"
#include "trialoop.h"
#include "ranmar.h"
#include "matrix.h"
#include "nodsurf.h"


static int    stride_x0;
static int    stride_y0;
static int    stride_z0;
static int    stride_x1;
static int    stride_y1;
static int    stride_z1;
static Fprec  interpfx;
static Fprec  interpfy;
static Fprec  interpfz;


void InitCalcIRius(void)
{
  Fprec  tovd; /* tRius over delta */


  tovd = tRius / LatConD.a * Nx;
  stride_x0 = (int) tovd;
  interpfx = tovd - (Fprec) stride_x0;
  stride_x0 %= Nx;
  stride_x1 = (stride_x0 + 1) % Nx;

  tovd = tRius / LatConD.b * Ny;
  stride_y0 = (int) tovd;
  interpfy = tovd - (Fprec) stride_y0;
  stride_y0 %= Ny;
  stride_y1 = (stride_y0 + 1) % Ny;

  tovd = tRius / LatConD.c * Nz;
  stride_z0 = (int) tovd;
  interpfz = tovd - (Fprec) stride_z0;
  stride_z0 %= Nz;
  stride_z1 = (stride_z0 + 1) % Nz;
}


static void CalcIRius(const Fprec *eD)
{
  int    ix, iy, iz;
  int    sx0, sy0, sz0;
  int    sx1, sy1, sz1;
  int    sy0Nz, sy1Nz;
  Fprec  add, sub;

  const Fprec  *r0, *rix, *rixiy;
  const Fprec  *rtx0, *rty0, *rtz0;
  const Fprec  *rtx1, *rty1, *rtz1;


  sx0 = stride_x0;
  sx1 = stride_x1;
  sy0 = stride_y0;
  sy1 = stride_y1;
  sz0 = stride_z0;
  sz1 = stride_z1;

  sy0Nz = sy0 * Nz;
  sy1Nz = sy1 * Nz;

  r0 = eD;
  rix = eD;
  rixiy = eD;

  IRius = 0.;

  rtx0 = eD + (sx0 * Ny * Nz);
  rtx1 = eD + (sx1 * Ny * Nz);

  for (ix = 0; ix < Nx; ix++)
  {
    rty0 = rix + sy0Nz;
    rty1 = rix + sy1Nz;

    for (iy = 0; iy < Ny; iy++)
    {
      rtz0 = rixiy + sz0;
      rtz1 = rixiy + sz1;

      for (iz = 0; iz < Nz; iz++)
      {
                               sub = (*rtx0) + interpfx * ((*rtx1) - (*rtx0));
        add =  (*r0) - AppFabs(sub);
                               sub = (*rty0) + interpfy * ((*rty1) - (*rty0));
        add *= (*r0) - AppFabs(sub);
                               sub = (*rtz0) + interpfz * ((*rtz1) - (*rtz0));
        add *= (*r0) - AppFabs(sub);

        IRius += add;

        rtz0 =  rtz1++;
        rty0++; rty1++;
        rtx0++; rtx1++;
        r0++;

        sz0 = sz1++;
        if (sz1 == Nz) sz1 = 0, rtz1 = rixiy;
      }
      sy0 = sy1++;
      if      (sy0 == 0)           rty0 = rix;
      else if (sy1 == Ny) sy1 = 0, rty1 = rix;

      rixiy = r0;
    }
    sx0 = sx1++;
    if      (sx0 == 0)           rtx0 = eD;
    else if (sx1 == Nx) sx1 = 0, rtx1 = eD;

    rix = r0;
  }

  IRius *= IRho3dVscale;
}


static void eD_Statistics(Fprec *eD)
{
  int    ieD;
  Fprec  rho;

  const Fprec *eD0;


  eD0 = eD;

            rho = *eD;
  IRho3dV = rho * rho * rho;

  eDminRho = eDmaxRho = rho;

  eD++;

  for (ieD = 1; ieD < Nx * Ny * Nz; ieD++, eD++)
  {
               rho = *eD;
    IRho3dV += rho * rho * rho;

    if      (eDminRho > rho) eDminRho = rho;
    else if (eDmaxRho < rho) eDmaxRho = rho;
  }

  IRho3dV *= IRho3dVscale;

  if (tRius != 0.)
    CalcIRius(eD0);
  else
        IRius = 0.;
}


static void CalcCT_Fcal(T_eD_PeakList *EPList, int nEPList,
                        T_CodeTransl *CT0, int nCT)
{
  int            iCT, iEPL, iAT, iPos;
  T_CodeTransl   *CT;
  Fprec          h, k, l;
  T_eD_PeakList  *EPL;
  T_fVector      *SE;
  Fprec          *FF, fj, phase;
  Fprec          Fcr, Fci;


  for (iCT = 0, CT = CT0; iCT < nCT; iCT++, CT++)
  {
    h = CT->FobsRaw->h;
    k = CT->FobsRaw->k;
    l = CT->FobsRaw->l;
    FF = CT->FF;
    Fcr = 0.;
    Fci = 0.;

    for (iEPL = 0, EPL = EPList; iEPL < nEPList; iEPL++, EPL++)
    {
          iAT = EPL->iAtomType;
      if (iAT >= 0)
      {
        fj = FF[iAT];
        SE = EPL->xSymEquiv;
               iPos = EPL->WL_Entry->nPositions;
        while (iPos--)
        {
          phase = (Fprec) TwoPI * (h * SE->x + k * SE->y + l * SE->z); SE++;
          Fcr += fj * AppCos(phase);
          Fci += fj * AppSin(phase);
        }
      }
    }

    CT->Fcal.r = Fcr;
    CT->Fcal.i = Fci;
    CT->FcalAbs = AppSqrt(Fcr * Fcr + Fci * Fci);
    if (CT->FcalAbs > BogFmrg)
      CT->FcalPhi = AppAtan2(Fci, Fcr);
    else
      CT->FcalPhi = 0.;
  }


  if (Debug0) /* Debug: Print CT_Fc */
  {
    Fprintf(stdout, ">Begin CT_Fc\n");

    for (iCT = 0, CT = CT0; iCT < nCT; iCT++, CT++)
    {
      Fprintf(stdout, "CT_Fc  %3d %3d %3d  %7.2f %4.0f  %7.2f %7.2f\n",
        CT->FobsRaw->h, CT->FobsRaw->k, CT->FobsRaw->l,
        CT->FcalAbs,
        (CT->FcalPhi < 0. ? CT->FcalPhi + TwoPI : CT->FcalPhi) / PIover180,
        CT->Fcal.r, CT->Fcal.i);
    }
    Fprintf(stdout, ">End CT_Fc\n");
  }
}


static int CountConsecutiveNodes(void)
{
  int            iEPL, nConsecutiveNodes;
  T_eD_PeakList  *EPL;


  nConsecutiveNodes = 0;

        EPL = eD_PeakList;
  for (iEPL = 0; iEPL < NeD_PeakList; iEPL++, EPL++)
  {
    if (EPL->iAtomType >= 0)
    {
      if (AtomType[EPL->iAtomType].Class != ATC_Node)
        break;

      nConsecutiveNodes += EPL->WL_Entry->nPositions;
    }
  }

  return nConsecutiveNodes;
}


void RecyclePhases(T_FourierParameters *FP,
                   T_PeakFlags *PeakFlags, Fprec *Surface,
                   T_PhaseCode *PhaseCode)
{
  int           PutHeader;
  char          buf[256], *cp;
  int           iCT, iD;
  int           FeedBack, FlagFirstR, FlagConverged;
  int           iCycle, iFBC, jFBC, FwFragWasUsed;
  int           Dphi360;
  Fprec         wSumDiff, R, LastR, DiffR;
  T_CodeTransl  *CT;
  T_PhaseCode   *PC;


  if (eD_PeakList) InternalError(NULL);

  PutHeader = 1;

  CurrCorrelationCoefficient = 0.;
  FlagFirstR = 1;
  LastR = 0.;
  FlagConverged = 0;
  FeedBack = 1;
  Re_i_Cycle = 0;
      iCycle = 0;
        iFBC = 0;

  for (;;)
  {
    DoCodeTranslation(PhaseCode,
                      FP->Density, FP->Mx, FP->My, FP->Mz, FP->Friedel);

    if (Re_i_Cycle || F_LoadDensityFileName == NULL)
      FourierTransform(FP);
    else
      LoadDensity(FP->Density);

    if (Surface)
      for (iD = 0; iD < Nx * Ny * Nz; iD++) FP->Density[iD] *= Surface[iD];

    if (F_Put_eDmap) Put_eDmap(FP->Density, NULL);

    eD_Statistics(FP->Density);

    jFBC = iFBC;

    while (iFBC < nFeedBackCycles)
    {
      if (FeedBackCycles[iFBC] == iCycle) {
        iCycle = 0;
        iFBC++;
      }
      else
        break;
    }

    if (iFBC < nFeedBackCycles)
      jFBC = iFBC;

    FwFragWasUsed = Build_eD_PeakList(FP->Density, PeakFlags, jFBC % 2);

    if (nAsyPeaks <= 0)
    {
      R = 9.99;
      wSumDiff =
      FeedBack = 0;
    }
    else
    {
      CalcCT_Fcal(eD_PeakList, NeD_PeakList, CodeTransl, nActivePhase);
      CalcCCandR(CodeTransl, nActivePhase, &CurrCorrelationCoefficient, &R);

      wSumDiff = 0.;

      CT = CodeTransl;
      PC = PhaseCode;

      for (iCT = 0; iCT < nActivePhase; iCT++, CT++, PC++)
      {
        if (CT->FcalAbs > BogFmrg)
        {
          Phi2Phi360(CT->FcalPhi, CT->Phi360Fcal);

          Dphi360 = CT->Phi360Fmrg - CT->Phi360Fcal;
          if (Dphi360 <   0) Dphi360 += 360;
          if (Dphi360 > 180) Dphi360 = 360 - Dphi360;

          wSumDiff += Dphi360 * CT->M_Fmrg / (Fprec) 180.;
        }
        else
        {
          CT->Phi360Fcal = CT->Phi360Fmrg;
          Dphi360 = 0;
        }

        if (Debug0) /* Debug: Print DPhi */
          Fprintf(stdout, "DPhi  %3d  %3d  %3d\n",
            CT->Phi360Fmrg, CT->Phi360Fcal, Dphi360);
      }

      if (wSumDiff > CF0)
      {
        switch (FeedBackBreakIf.Branch)
        {
          case 0:
            break;
          case 1:
#define CondDeltaR \
            { \
              if (FlagFirstR == 0) \
              { \
                DiffR = LastR - R; \
                FlagConverged \
                  = AppFabs(DiffR) <= FeedBackBreakIf.DeltaR.Value; \
              } \
              LastR = R; \
              FlagFirstR = 0; \
            }
            CondDeltaR;
            break;
          case 4:
#define CondPhaseDiff \
            { \
              if (FeedBackBreakIf.PhaseDiff.Percent) \
                FlagConverged \
                  = wSumDiff / SumFmrg <= FeedBackBreakIf.PhaseDiff.Value; \
              else \
                FlagConverged \
                  = wSumDiff <= FeedBackBreakIf.PhaseDiff.Value; \
            }
            CondPhaseDiff;
            break;
          case 5:
            CondDeltaR;
            if (FlagConverged == 1)
              CondPhaseDiff;
            break;
          case 7:
            CondDeltaR;
            if (FlagConverged == 0)
              CondPhaseDiff;
            break;
          default:
            InternalError("switch(FeedBackBreakIf.Branch)");
            break;
        }

        if (FlagConverged)
          FeedBack = 0;
      }
      else
      {
        FlagConverged = 1;
        FeedBack = 0;
      }
    }

    if (iFBC == nFeedBackCycles)
      FeedBack = 0;

    if (FeedBack)
      Sprintf(buf, "@ %3d ",  Re_i_Cycle);
    else
      Sprintf(buf, "@ %+3d ", Re_i_Cycle);

    cp = buf; while (*cp) cp++;

    if (FwFragWasUsed) *cp++ = 'F';
    else               *cp++ = 'A';

    Sprintf(cp,
      " %5.3f %4.2f %4.1f %5.1f %5.1f %7.0f %7.0f %3d %3d %3d %ld %d",
      CurrCorrelationCoefficient, R, wSumDiff  / SumFmrg * 100.,
      eDminRho, eDmaxRho, IRho3dV, IRius,
      nPeaks, nAsyPeaks, NeD_PeakList,
      CurrPhaseCode_nCallRanmar, nTrials);

    if (F_Put_strudat)
      DoPut_strudat(eD_PeakList, NeD_PeakList, Re_i_Cycle,
                    buf, F_PutAllPeaks);

    if (   MinConsecutiveNodes >= 0
        && MinConsecutiveNodes <= CountConsecutiveNodes())
    {
      if (PutHeader)
      {
        Fprintf(stdout, "   #C M  CC    R     %%D eDmin eDmax");
        Fprintf(stdout, " Irho3dV   IRius  #U  #A  #T #RC #T\n");
        PutHeader = 0;
      }

      Fprintf(stdout, "%s\n", buf);
      Fflush(stdout);
    }

    if (FeedBack == 0)
      break;

    CT = CodeTransl;
    PC = PhaseCode;

    for (iCT = 0; iCT < nActivePhase; iCT++, CT++, PC++)
    {
      *PC = CT->Phi360Fcal - CT->CTData->TH;
      if (*PC < 0) *PC += 360;

      if (CT->FobsRaw->PhaseRestriction >= 0)
      {
        if (*PC != 0 && *PC != 180)
          CheckRestrictedPhase(CT, PC);
      }
    }

    Re_i_Cycle++;
        iCycle++;
  }

  CurrCorrelationCoefficient = 0.;

  if (FlagConverged)    nConvergedPhases++;
  else               nNotConvergedPhases++;

  Free_eD_PeakList();
}
