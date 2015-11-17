#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "main.h"
#include "xtal.h"
#include "io.h"


typedef struct
  {
    T_SymNodes  *SN;
    T_PUC       PUC;
    int         Count;
  }
  T_LoopConfPool;


static int LoopCountSF(const T_LoopCount *a, const T_LoopCount *b)
{
  if (a->Loop != 0 && b->Loop == 0) return -1;
  if (a->Loop == 0 && b->Loop != 0) return  1;
  if (a->Loop < b->Loop) return -1;
  if (a->Loop > b->Loop) return  1;
  if (a->Count < b->Count) return -1;
  if (a->Count > b->Count) return  1;
  return 0;
}


#define IxLoopCount(i, j, n)\
  (((n) * ((n) - 1) - ((n) - (i)) * ((n) - (i) - 1)) / 2 + (j) - (i) - 1)


static void BuildLoopConf(T_ListSymNodes *LSymN)
{
  int  n, iEPL, iBond, jBond, SymN_Offset;

  T_SymNodeBonds  *NB;
  T_LoopConfPool  *Pool, NewN, *PMi, *PMj;
  T_iVector       PMi_UC, NewN_UC;
  int             sPool, AllocIncrement;

  int  iOffs, OffsCurr, OffsNew, OffsEnd, NewOffsEnd;
  int  iSphere, Loop;

  T_SymNodes  *PivotSN;
  int         PiBond, WtdHits;

  T_LoopConf   *LCf;
  int          iAng, LS, iAsyN;


  CheckMalloc(LSymN->LCf, LSymN->nAsyN);

  LCf         = LSymN->LCf;
  SymN_Offset = 0;

  for (iEPL = LSymN->Low_iEPL; iEPL <= LSymN->Top_iEPL; iEPL++)
  {
    if (LSymN->LnB[iEPL] == 0)
      continue;

    if ((int)(LCf - LSymN->LCf) == LSymN->nAsyN)
      InternalError(NULL);

    PivotSN = LSymN->SymN + SymN_Offset;

    SymN_Offset += eD_PeakList[iEPL].WL_Entry->nPositions;

    LCf->nAngles = PivotSN->nBonds * (PivotSN->nBonds - 1) / 2;

    if (LCf->nAngles == 0)
      LCf->LC = NULL;
    else
    {
      CheckMalloc(LCf->LC, LCf->nAngles);

      for (iAng = 0; iAng < LCf->nAngles; iAng++) {
        LCf->LC[iAng].Loop  = 0;
        LCf->LC[iAng].Count = 0;
      }
    }

    LCf++;
  }

  AllocIncrement = 10000; /* ARBITRARY */

                    sPool = 1000; /* ARBITRARY */
  CheckMalloc(Pool, sPool);

  LCf         = LSymN->LCf;
  SymN_Offset = 0;

  for (iEPL = LSymN->Low_iEPL; iEPL <= LSymN->Top_iEPL; iEPL++)
  {
    if (LSymN->LnB[iEPL] == 0)
      continue;

    if ((int)(LCf - LSymN->LCf) == LSymN->nAsyN)
      InternalError(NULL);

    PivotSN = LSymN->SymN + SymN_Offset;

    SymN_Offset += eD_PeakList[iEPL].WL_Entry->nPositions;

    for (PiBond = 0; PiBond < PivotSN->nBonds - 1; PiBond++)
    {
      WtdHits = PivotSN->nBonds - (PiBond + 1);

      Pool[0].SN = PivotSN->NB[PiBond].NextSymN;

      PackUC(Pool[0].PUC, &(PivotSN->NB[PiBond].UC_Offset));

      Pool[0].Count = 1;

      OffsCurr = 0;
      OffsNew  =
      OffsEnd  = 1;

      for (iSphere = 0;; iSphere++)
      {
        Loop = iSphere + 3;

        if (Loop > LSymN->MaxLoopSize)
          break;

        PMi = Pool + OffsCurr;

        for (iOffs = OffsCurr; iOffs < OffsNew; iOffs++, PMi++)
        {
          UnpackUC(PMi->PUC, &PMi_UC);

          NB = PMi->SN->NB;

          for (iBond = 0; iBond < PMi->SN->nBonds; iBond++)
          {
            NewN.SN   = NB->NextSymN;
            NewN_UC.x = PMi_UC.x + NB->UC_Offset.x;
            NewN_UC.y = PMi_UC.y + NB->UC_Offset.y;
            NewN_UC.z = PMi_UC.z + NB->UC_Offset.z;

            if (   NewN.SN   == PivotSN
                && NewN_UC.x == 0
                && NewN_UC.y == 0
                && NewN_UC.z == 0)
              goto NextBond;

            if (CheckUC(&NewN_UC))
              progerror("|NewN_UC| overflow");

            PackUC(NewN.PUC, &NewN_UC);

            PMj = Pool;

            for (n = OffsEnd; n--; PMj++)
            {
              if (   PMj->SN  == NewN.SN
                  && PMj->PUC == NewN.PUC)
              {
                if (PMj - Pool >= OffsNew)
                  PMj->Count += PMi->Count;

                goto CheckHit;
              }
            }

            if (OffsEnd == sPool)
            {
                                       n = sPool + AllocIncrement;
              CheckRealloc(Pool, Pool, n,  sPool);

              sPool = n;

              PMi = Pool + iOffs;

              if (OffsCurr != 0)
              {
                int  NewAllocIncr;

                NewAllocIncr = OffsNew - OffsCurr;
                NewAllocIncr *= NewAllocIncr;
                NewAllocIncr = (int)((Fprec) NewAllocIncr / OffsCurr);

                if (AllocIncrement < NewAllocIncr)
                    AllocIncrement = NewAllocIncr;
              }
            }

            PMj = Pool + OffsEnd;

            PMj->SN    = NewN.SN;
            PMj->PUC   = NewN.PUC;
            PMj->Count = PMi->Count;

            OffsEnd++;

            CheckHit:

            for (jBond = PiBond + 1; jBond < PivotSN->nBonds; jBond++)
            {
              if (   PivotSN->NB[jBond].NextSymN     == NewN.SN
                  && PivotSN->NB[jBond].UC_Offset.x  == NewN_UC.x
                  && PivotSN->NB[jBond].UC_Offset.y  == NewN_UC.y
                  && PivotSN->NB[jBond].UC_Offset.z  == NewN_UC.z)
              {
                iAng = IxLoopCount(PiBond, jBond, PivotSN->nBonds);

                if      (LCf->LC[iAng].Loop == Loop)
                  LCf->LC[iAng].Count += PMi->Count;

                else if (LCf->LC[iAng].Loop ==    0)
                {
                  WtdHits--;

                  LCf->LC[iAng].Loop  = Loop;
                  LCf->LC[iAng].Count = PMi->Count;
                }

                if (Debug0 && LCf->LC[iAng].Loop == Loop)
                  Fprintf(stdout, "Loop  %d - %3.3d - %d  %d:%d\n",
                    PiBond, iEPL, jBond,
                    LCf->LC[iAng].Loop,
                    LCf->LC[iAng].Count);
              }
            }

            NextBond:

            NB++;
          }
        }

        if (WtdHits == 0)
          break;

        NewOffsEnd = OffsEnd - OffsCurr;

        if (OffsCurr && NewOffsEnd)
          (void) memmove(&Pool[0], &Pool[OffsCurr],
                         sizeof (*Pool) * NewOffsEnd);

        OffsCurr = OffsNew - OffsCurr;
        OffsNew  =
        OffsEnd  = NewOffsEnd;
      }
    }

    if (LCf->nAngles)
      qsort((void *) LCf->LC, LCf->nAngles, sizeof (*LCf->LC),
            (SortFunction) LoopCountSF);

    LCf++;
  }

  AppFree(Pool, sPool);

  CheckMalloc(LSymN->LS_Flags, LSymN->MaxLoopSize + 1);

  for (LS = 0; LS <= LSymN->MaxLoopSize; LS++)
    LSymN->LS_Flags[LS] = 0;

  LCf = LSymN->LCf;

  for (iAsyN = 0; iAsyN < LSymN->nAsyN; iAsyN++, LCf++)
    for (iAng = 0; iAng < LCf->nAngles; iAng++)
      LSymN->LS_Flags[LCf->LC[iAng].Loop]++;
}


static void FreeLoopConf(T_ListSymNodes *LSymN)
{
  int         iAsyN;
  T_LoopConf  *LCf;


  LCf = LSymN->LCf;

  for (iAsyN = 0; iAsyN < LSymN->nAsyN; iAsyN++, LCf++)
    if (LCf->LC)
      AppFree(LCf->LC, LCf->nAngles);

  AppFree(LSymN->LCf, LSymN->nAsyN);
          LSymN->LCf = NULL;

  AppFree(LSymN->LS_Flags, LSymN->MaxLoopSize + 1);
          LSymN->LS_Flags = NULL;
}


static void PrintNodeLabel(const T_ListSymNodes *LSymN, int iEPL)
{
  if (F_SiteFrame && F_SiteLabel) {
    Fprintf(stdout, "%-4s", Site[eD_PeakList[iEPL].Index].Label);
  }
  else
  {
    if (LSymN->Colored == 0 || NodeColor(iEPL) == 0)
      putc('T', stdout);
    else
      putc('U', stdout);

    Fprintf(stdout, "%3.3d", eD_PeakList[iEPL].Index);
  }
}


static void PrintLoopConf(const T_ListSymNodes *LSymN)
{
  int         LS, iAng, iEPL;
  T_LoopConf  *LCf;


  Fprintf(stdout, ">Begin LoopConf");
  if (IndicateFw) Fprintf(stdout, ":%4.4d", iFramework);
  putc('\n', stdout);

  LCf = LSymN->LCf;

  for (iEPL = LSymN->Low_iEPL; iEPL <= LSymN->Top_iEPL; iEPL++)
  {
    if (LSymN->LnB[iEPL] == 0)
      continue;

    if ((int)(LCf - LSymN->LCf) == LSymN->nAsyN)
      InternalError(NULL);

    PrintNodeLabel(LSymN, iEPL);

    for (iAng = 0; iAng < LCf->nAngles; iAng++) {
      if (LCf->LC[iAng].Loop)
        Fprintf(stdout, " %d:%d",
          LCf->LC[iAng].Loop,
          LCf->LC[iAng].Count);
      else
        Fprintf(stdout, " >%d", LSymN->MaxLoopSize);
    }

    putc('\n', stdout);

    LCf++;
  }

  Fprintf(stdout, "Total");

  for (LS = 1; LS <= LSymN->MaxLoopSize; LS++)
    if (LSymN->LS_Flags[LS])
      Fprintf(stdout, " %d", LS);

  if (LSymN->LS_Flags[0])
    Fprintf(stdout, " >%d", LSymN->MaxLoopSize);

  putc('\n', stdout);

  Fprintf(stdout, ">End LoopConf\n\n");
}


static int IntegerSortFunction(const int *a, const int *b)
{
  if (*a < *b) return -1;
  if (*a > *b) return  1;

  return 0;
}


static void TopologyInfo(int *SphereStorage, int nAsyN, int nSphereStorage)
{
  int  i;
  int  *SSiAsyN, iAsyN;
  int  *SSjAsyN, jAsyN;
  int  *Flags, nPairs, nDistinct;

  int  iEPL, TDn;
  int  M, SumM, SumNk, SumMSumNk;


  CheckMalloc(Flags, nAsyN);

  nPairs = 0;

  if (FwSearchMethod == FwSM_AltFwTracking)
  {
    for (iAsyN = 0; iAsyN < nAsyN; iAsyN++)
      Flags[iAsyN] = 0;

    SSiAsyN = SphereStorage + nAsyN;

    for (iAsyN = 0; iAsyN < nAsyN - 1; iAsyN++, SSiAsyN += nSphereStorage)
    {
      if (Flags[iAsyN])
        continue;

      SSjAsyN = SphereStorage + nAsyN + (iAsyN + 1) * nSphereStorage;

      for (jAsyN = iAsyN + 1; jAsyN < nAsyN; jAsyN++, SSjAsyN += nSphereStorage)
      {
        if (   Flags[jAsyN]
            ||    NodeColor(SphereStorage[iAsyN])
               == NodeColor(SphereStorage[jAsyN]))
          continue;

        for (i = 0; i < nSphereStorage; i++)
          if (SSiAsyN[i] != SSjAsyN[i])
            break;

        if (i == nSphereStorage)
        {
          Flags[iAsyN] = 1;
          Flags[jAsyN] = 1;
          nPairs++;
          break;
        }
      }
    }
  }

  for (iAsyN = 0; iAsyN < nAsyN; iAsyN++)
    Flags[iAsyN] = 0;

  nDistinct = 0;

  SSiAsyN = SphereStorage + nAsyN;

  for (iAsyN = 0; iAsyN < nAsyN; iAsyN++, SSiAsyN += nSphereStorage)
  {
    if (Flags[iAsyN])
      continue;

    Flags[iAsyN]++;
    nDistinct++;

    SSjAsyN = SphereStorage + nAsyN + (iAsyN + 1) * nSphereStorage;

    for (jAsyN = iAsyN + 1; jAsyN < nAsyN; jAsyN++, SSjAsyN += nSphereStorage)
    {
      if (Flags[jAsyN])
        continue;

      for (i = 0; i < nSphereStorage; i++)
        if (SSiAsyN[i] != SSjAsyN[i])
          break;

      if (i == nSphereStorage)
      {
        Flags[iAsyN]++;
        Flags[jAsyN] = -1;
        break;
      }
    }
  }

  if (nAsyN > 1)
    qsort((void *) Flags, nAsyN, sizeof (*Flags),
          (SortFunction) IntegerSortFunction);

  Fprintf(stdout, "# TopologyInfo");
  if (IndicateFw) Fprintf(stdout, ":%4.4d", iFramework);

  Fprintf(stdout, " %d : %d :", nAsyN, nDistinct);

  for (iAsyN = nAsyN - 1; iAsyN >= 0; iAsyN--)
  {
    if (Flags[iAsyN] < 0)
      break;

    Fprintf(stdout, " %d", Flags[iAsyN]);
  }

  if (FwSearchMethod == FwSM_AltFwTracking)
    Fprintf(stdout, " : %d", nPairs);

  putc('\n', stdout);

  AppFree(Flags, nAsyN);

      TDn = 10;
  if (TDn > nSphereStorage)
      TDn = nSphereStorage;

  SumM      = 0;
  SumMSumNk = 0;

  SSiAsyN = SphereStorage + nAsyN;

  for (iAsyN = 0; iAsyN < nAsyN; iAsyN++, SSiAsyN += nSphereStorage)
  {
    SumNk = 1;

    for (i = 0; i < TDn; i++)
      SumNk += SSiAsyN[i];
                                 iEPL = SphereStorage[iAsyN];
                 M = eD_PeakList[iEPL].WL_Entry->nPositions;
    SumM      += M;
    SumMSumNk += M * SumNk;
  }

  Fprintf(stdout, "# TD%d", TDn);
  if (IndicateFw) Fprintf(stdout, ":%4.4d", iFramework);

  if (SumM)
    Fprintf(stdout, " %.6g\n", (double) SumMSumNk / SumM);
  else
    Fprintf(stdout, " 0\n");
}


static void CalcCoseqSplit(double a, double b, double c,
                           int nCellC, int nUnits,
                           int *F_x, int *F_y, int *F_z)
{
  double  fx, fy, fz;


  fx = pow(b * c * nUnits / (a * a * nCellC), 1./3.);
  fy = a / b * fx;
  fz = a / c * fx;

  *F_x = (int)(fx + 2./3.); if (*F_x < 1) *F_x = 1;
  *F_y = (int)(fy + 2./3.); if (*F_y < 1) *F_y = 1;
  *F_z = (int)(fz + 2./3.); if (*F_z < 1) *F_z = 1;

}


static void AllocPUC(T_PUC_Buf *RootPB, int nSNFxFyFz)
{
  int        iPB, sPUC;
  T_PUC_Buf  *PBi;


  if (RootPB->PUC)
    InternalError("AllocPUC(): RootPB->PUC != NULL");

              /* ARBITRARY */
      sPUC = (int)(10000. / nSNFxFyFz + .5);
  if (sPUC < 10) sPUC = 10; /* ARBITRARY */

  PBi = RootPB;

  for (iPB = 0; iPB < nSNFxFyFz; iPB++, PBi++)
  {
                          PBi->sPUC = sPUC;
    CheckMalloc(PBi->PUC, PBi->sPUC);

#ifdef COSEQ_TRACE                             /* ARBITRARY */
                           PBi->sNgen = (int)(sPUC / 2.5 + .5);
    CheckMalloc(PBi->Ngen, PBi->sNgen);
#endif
  }
}


static void CoseqRun(const T_ListSymNodes *LSymN)
{
  int  i, n, New_oEnd, iSN, iEPL, iBond, SymN_Offset;
  int  nNew, LastnNew;

  T_SymNodes      *SN;
  T_SymNodeBonds  *NB;
  T_iVector       NewUC, PBiUC;
  T_PUC           NewPUC, *PUCi, *PUCj;

  int  iOffs;
  int  LastSphere, iSphere;

  int  *SphereStorage, *SSiAsyN, iAsyN;
  int  nSphereStorage;

  T_Ticks  LastSaveTime, CurrentTicks;
  double   TI;
  int      CoseqSaveDone;

  int        FxFyFz, nSNFxFyFz;
  int        iPB, jPB;
  T_PUC_Buf  *RootPB, *PBi, *PBj;

#ifdef COSEQ_TRACE
  int  *Slot;
  int  nSlots;

  nSlots = 0;

  SN = LSymN->SymN;

  for (iSN = 0; iSN < LSymN->nSymN; iSN++, SN++)
    if (nSlots < SN->nBonds)
        nSlots = SN->nBonds;

  if (nSlots == 0)
    InternalError("CoseqRun(): Corrupt nSlots");

                    nSlots++;
  CheckMalloc(Slot, nSlots);
#endif


  if (F_CoseqSaveFileName) (void) GetTicks(&LastSaveTime);

  LastSphere = F_CoseqValue;

  if (F_CoseqUnitCell)
    nSphereStorage = 0;
  else
  {
        nSphereStorage = 50; /* ARBITRARY */
    if (nSphereStorage > LastSphere)
        nSphereStorage = LastSphere;
  }

  CheckMalloc(SphereStorage, LSymN->nAsyN + LSymN->nAsyN * nSphereStorage);

  iAsyN = 0;

  for (iEPL = LSymN->Low_iEPL; iEPL <= LSymN->Top_iEPL; iEPL++)
  {
    if (LSymN->LnB[iEPL])
    {
      if (iAsyN >= LSymN->nAsyN)
        InternalError("CoseqRun(): Corrupt LSymN->nAsyN");

      SphereStorage[iAsyN++] = iEPL;
    }
  }

  for (i = 0; i < LSymN->nAsyN * nSphereStorage; i++)
    SphereStorage[LSymN->nAsyN + i] = 0;

  if      (F_nCoseqSplit == 1)
  {
    CalcCoseqSplit(LatConD.a, LatConD.b, LatConD.c,
                   LSymN->nSymN, F_CoseqSplit[0],
                   &F_CoseqSplit_x,
                   &F_CoseqSplit_y,
                   &F_CoseqSplit_z);
  }
  else if (F_nCoseqSplit != 3)
    InternalError("CoseqRun(): Corrupt F_nCoseqSplit");

     FxFyFz = F_CoseqSplit_x * F_CoseqSplit_y * F_CoseqSplit_z;
  nSNFxFyFz = LSymN->nSymN * FxFyFz;

  if (F_nCoseqSplit == 1)
    Fprintf(stdout, "# CoseqSplit %d * %d * %d * %d = %d (* %d = %d)\n\n",
      F_CoseqSplit_x, F_CoseqSplit_y, F_CoseqSplit_z, LSymN->nSymN,
      nSNFxFyFz, (int) sizeof (*RootPB), (int)(nSNFxFyFz * sizeof (*RootPB)));

  CheckMalloc(RootPB, nSNFxFyFz);

  SN  = LSymN->SymN;
  PBi = RootPB;

  for (iSN = 0; iSN < LSymN->nSymN; iSN++, SN++)
  {
    SN->PB = PBi;

    for (iPB = 0; iPB < FxFyFz; iPB++, PBi++)
    {
      PBi->PUC  = NULL;
      PBi->sPUC = 0;
#ifdef COSEQ_TRACE
      PBi->Ngen  = NULL;
      PBi->sNgen = 0;
#endif
    }
  }

  Fprintf(stdout, ">Begin coseq");
  if (IndicateFw) Fprintf(stdout, ":%4.4d", iFramework);
  putc('\n', stdout);

      SSiAsyN = SphereStorage + LSymN->nAsyN;
        iAsyN = 0;
  SymN_Offset = 0;

  for (iEPL = LSymN->Low_iEPL; iEPL <= LSymN->Top_iEPL; iEPL++)
  {
    if (LSymN->LnB[iEPL] == 0)
      continue;

    PBi = RootPB;

    for (iPB = 0; iPB < nSNFxFyFz; iPB++, PBi++)
    {
      PBi->oCurr = 0;
      PBi->oNew  = 0;
      PBi->oEnd  = 0;
    }

    if (F_CoseqReloadFileName)
    {
      nNew =
        CoseqReload(LSymN, FxFyFz,
                    &iEPL, &SymN_Offset, &iSphere,
                    SphereStorage + LSymN->nAsyN, nSphereStorage, &iAsyN);

      PrintTicks(stdout, NULL, "# Time Coseq Reload done  ", "\n");
      Fprintf(stdout, "# MaxMemoryUsed = %ld\n", (long) MaxMemoryUsed);
      Fflush(stdout);

      if (   iEPL < LSymN->Low_iEPL
          || iEPL > LSymN->Top_iEPL
          || LSymN->LnB[iEPL] == 0)
        progerror("CoseqReload mismatch");

      if (   F_CoseqSelect >= 0
          && F_CoseqSelect != eD_PeakList[iEPL].Index)
        break;

      SSiAsyN = &SphereStorage[LSymN->nAsyN + iAsyN * nSphereStorage];

      PrintNodeLabel(LSymN, iEPL);

      Fprintf(stdout, " (%d missing, last was %d)", iSphere, nNew);
      Fflush(stdout);

      AppFree(F_CoseqReloadFileName, strlen(F_CoseqReloadFileName) + 1);
              F_CoseqReloadFileName = NULL;
    }
    else if (F_CoseqUnitCell == 0)
    {
      SN = LSymN->SymN + SymN_Offset;

      SymN_Offset += eD_PeakList[iEPL].WL_Entry->nPositions;

      if (   F_CoseqSelect >= 0
          && F_CoseqSelect != eD_PeakList[iEPL].Index)
      {
        SSiAsyN += nSphereStorage;
          iAsyN++;

        continue;
      }
                             NewUC.x = NewUC.y = NewUC.z = 0;
             PackUC(NewPUC, &NewUC);
      Calc_iPUC_Buf(iPB,    &NewUC);

      if (RootPB->PUC == NULL)
        AllocPUC(RootPB, nSNFxFyFz);

      PBi = &(SN->PB[iPB]);
      PBi->PUC[0] = NewPUC;
      PBi->oNew =
      PBi->oEnd = 1;

      nNew = 0;
      iSphere = 0;

      PrintNodeLabel(LSymN, iEPL);
    }
    else
    {
                             NewUC.x = NewUC.y = NewUC.z = 0;
             PackUC(NewPUC, &NewUC);
      Calc_iPUC_Buf(iPB,    &NewUC);

      AllocPUC(RootPB, nSNFxFyFz);

      SN = LSymN->SymN;

      for (iSN = 0; iSN < LSymN->nSymN; iSN++, SN++)
      {
        PBi = &(SN->PB[iPB]);
        PBi->PUC[0] = NewPUC;
        PBi->oNew =
        PBi->oEnd = 1;
      }

      nNew = 0;
      iSphere = 0;

      Fprintf(stdout, "CSUC %3d", LSymN->nSymN);
    }

    if (F_CoseqProtocolFileName)
      CoseqProtocol(-1, 0, 0); /* set protocol file "time 0" */

    for (; iSphere < LastSphere; iSphere++)
    {
      LastnNew = nNew;
          nNew = 0;

      SN  = LSymN->SymN;
      PBi = RootPB;

      for (iSN = 0; iSN < LSymN->nSymN; iSN++, SN++)
      {
        for (iPB = 0; iPB < FxFyFz; iPB++, PBi++)
        {
                            iOffs = PBi->oCurr;
          PUCi = PBi->PUC + iOffs;

          for (; iOffs < PBi->oNew; iOffs++, PUCi++)
          {
            UnpackUC(*PUCi, &PBiUC);

            NB = SN->NB;

            for (iBond = 0; iBond < SN->nBonds; iBond++)
            {
              NewUC.x = PBiUC.x + NB->UC_Offset.x;
              NewUC.y = PBiUC.y + NB->UC_Offset.y;
              NewUC.z = PBiUC.z + NB->UC_Offset.z;

              if (CheckUC(&NewUC))
                progerror("CoseqRun(): |NewUC| overflow");

                  PackUC(NewPUC, &NewUC);
              Calc_iPUC_Buf(jPB, &NewUC);

                     PBj = &NB->NextSymN->PB[jPB];
              PUCj = PBj->PUC;
                 n = PBj->oEnd;

#ifndef COSEQ_TRACE
              while (n--)
                if (*PUCj++ == NewPUC)
                  goto NextBond;
#else
              for (i = 0; i < n; i++)
              {
                if (*PUCj++ == NewPUC)
                {
                      i -= PBj->oNew;
                  if (i >= 0)
                    PBj->Ngen[i]++;

                  goto NextBond;
                }
              }
#endif

              if (PBj->oEnd == PBj->sPUC)
              {                                           /* ARBITRARY */
                    n = (int)((double) LastnNew / nSNFxFyFz * 3.1 + .5);
                if (n <= PBj->sPUC)
                    n = (int)(PBj->sPUC * 1.1 + .5); /* ARBITRARY */
                if (n - PBj->sPUC < 10) /* ARBITRARY */
                    n = PBj->sPUC + 10;

                CheckRealloc(PBj->PUC, PBj->PUC, n, PBj->sPUC);

                PBj->sPUC = n;

                PUCi = PBi->PUC + iOffs;
                PUCj = PBj->PUC + PBj->oEnd;
              }

#ifdef COSEQ_TRACE
                  i = PBj->oEnd - PBj->oNew;
              if (i == PBj->sNgen)
              {
                    n = PBj->sPUC / 2.5; /* ARBITRARY */
                if (n <= PBj->sNgen)
                    n = (int)(PBj->sNgen * 1.1 + .5); /* ARBITRARY */
                if (n - PBj->sNgen < 10) /* ARBITRARY */
                    n = PBj->sNgen + 10;

                CheckRealloc(PBj->Ngen, PBj->Ngen, n, PBj->sNgen);

                PBj->sNgen = n;
              }

              PBj->Ngen[i] = 1;
#endif

              *PUCj = NewPUC;
              PBj->oEnd++;
              nNew++;

              NextBond:

              NB++;
            }
          }
        }
      }

#ifdef COSEQ_TRACE
      for (i = 1; i < nSlots; i++)
        Slot[i] = 0;

      PBi = RootPB;

      for (iPB = 0; iPB < nSNFxFyFz; iPB++, PBi++)
      {
        int *Ngen = PBi->Ngen;

               n = PBi->oEnd - PBi->oNew;
        while (n--)
        {
          if (0 >= *Ngen || *Ngen >= nSlots)
            InternalError("CoseqRun(): Corrupt *Ngen");

          Slot[*Ngen++]++;
        }
      }

      Fprintf(stdout, "\n[");

      n = 0;

      for (i = 1; i < nSlots; i++)
      {
        Fprintf(stdout, " %d", Slot[i]);
        n += i * Slot[i];
      }

      Fprintf(stdout, " : %d : %d ]", n, n - nNew);
#endif

      if (iSphere < 10)
        Fprintf(stdout, " %4d", nNew);
      else
      {
#ifndef COSEQ_TRACE
        if (iSphere % 10 == 0) Fprintf(stdout, "\n    ");
#endif
        Fprintf(stdout, " %d", nNew);
#if ! (defined(__ALPHA) && defined(__VMS))
        Fflush(stdout);
#endif
      }

      if (F_CoseqProtocolFileName)
        CoseqProtocol(iEPL, iSphere, nNew);

      if (iSphere < nSphereStorage)
      {
        SSiAsyN[iSphere] = nNew;

        if (F_CoseqSelect < 0 && iSphere + 1 == nSphereStorage)
        {
          int  *SSj, j, k;

          SSj = SphereStorage + LSymN->nAsyN;

          for (j = 0; j < iAsyN; j++, SSj += nSphereStorage)
          {
            for (k = 0; k < nSphereStorage; k++)
              if (SSiAsyN[k] != SSj[k]) break;

            if (k == nSphereStorage)
              break;
          }

          if (j < iAsyN)
          {
            if (nSphereStorage != LastSphere)
              Fprintf(stdout, " ...");

            Fprintf(stdout, "\n    N1 to N%d equal those of ",
              nSphereStorage);

            PrintNodeLabel(LSymN, SphereStorage[j]);

            break;
          }
        }
      }

      SN  = LSymN->SymN;
      PBi = RootPB;

      for (iSN = 0; iSN < LSymN->nSymN; iSN++, SN++)
      {
        for (iPB = 0; iPB < FxFyFz; iPB++, PBi++)
        {
          New_oEnd = PBi->oEnd - PBi->oCurr;

          if (New_oEnd && PBi->oCurr)
            (void) memmove(&(PBi->PUC[0]), &(PBi->PUC[PBi->oCurr]),
                           sizeof (*(PBi->PUC)) * New_oEnd);

          PBi->oCurr = PBi->oNew - PBi->oCurr;
          PBi->oNew  =
          PBi->oEnd  = New_oEnd;
        }
      }

      CoseqSaveDone = 0;

      if (F_CoseqSaveFileName)
      {
        (void) GetTicks(&CurrentTicks);

        TI =  (double)(CurrentTicks.cpu_user   - LastSaveTime.cpu_user);
        TI += (double)(CurrentTicks.cpu_system - LastSaveTime.cpu_system);
        TI /= CurrentTicks.ticks_per_second * 60.;
        TI += .5;

        if ((int) TI >= F_CoseqSaveTime || iSphere + 1 == LastSphere)
        {
          CoseqSave(LSymN, FxFyFz,
                    iEPL, SymN_Offset, iSphere + 1,
                    SphereStorage + LSymN->nAsyN, nSphereStorage, iAsyN);

          CoseqSaveDone = 1;

          (void) GetTicks(&LastSaveTime);
        }
      }

      if (iSphere >= 20 - 1) /* ARBITRARY */
      {
        if (F_SignalFileName) CheckSignalFile();

        if (QuitProgram)
        {
          putc('\n', stdout);
          Fprintf(stdout, ">End coseq\n\n");

          if (F_CoseqSaveFileName && CoseqSaveDone == 0)
            CoseqSave(LSymN, FxFyFz,
                      iEPL, SymN_Offset, iSphere + 1,
                      SphereStorage + LSymN->nAsyN, nSphereStorage, iAsyN);
          AppExit(1);
        }
      }
    }

    putc('\n', stdout);
    Fflush(stdout);

    if (F_CoseqUnitCell)
      break;

    SSiAsyN += nSphereStorage;
      iAsyN++;
  }

  Fprintf(stdout, ">End coseq\n\n");

  if (F_CoseqSelect < 0 && F_CoseqUnitCell == 0)
  {
    TopologyInfo(SphereStorage, LSymN->nAsyN, nSphereStorage);
    putc('\n', stdout);
  }

  Fflush(stdout);

  PBi = RootPB;

  for (iPB = 0; iPB < nSNFxFyFz; iPB++, PBi++)
  {
    AppFree(PBi->PUC, PBi->sPUC);
#ifdef COSEQ_TRACE
    AppFree(PBi->Ngen, PBi->sNgen);
#endif
  }

  AppFree(RootPB, nSNFxFyFz);
  AppFree(SphereStorage, LSymN->nAsyN + LSymN->nAsyN * nSphereStorage);

#ifdef COSEQ_TRACE
  AppFree(Slot, nSlots);
#endif
}


static int ThreeDimConCheck(const T_ListSymNodes *LSymN)
{
  int             *IsCon; /* IsCon[UCx + 1][UCy + 1][UCz + 1][iSymN] */
  int             nIsCon, iIC, nInCell, Is3D;
  int                 UCx,     UCy,     UCz,     iSN, nSN;
  long            NextUCx, NextUCy, NextUCz, iNextSN, iNextIC;
  int             nNB, iNB;
  T_SymNodeBonds  *NB;

#define IxIC(x, y, z, i) (((((x + 1) * 3) + (y + 1)) * 3 + (z + 1)) * nSN + i)


      nSN = LSymN->nSymN;
  if (nSN == 0)
    return 0;

  nIsCon = nSN * 27;

  CheckMalloc(IsCon, nIsCon);

  for (iIC = 0; iIC < nIsCon; iIC++)
    IsCon[iIC] = 0;

  IsCon[IxIC(0, 0, 0, 0)] = -1;

  for (;;)
  {
    iIC = 0;

    for (UCx = -1; UCx <= 1; UCx++)
    for (UCy = -1; UCy <= 1; UCy++)
    for (UCz = -1; UCz <= 1; UCz++)
    for (iSN = 0; iSN < nSN; iSN++, iIC++)
      if (IsCon[iIC] == -1) goto FollowBonds;

    break;

    FollowBonds:

    nNB = LSymN->SymN[iSN].nBonds;
     NB = LSymN->SymN[iSN].NB;

    for (iNB = 0; iNB < nNB; iNB++, NB++)
    {
#define                                    NotValidUC(x) (x < -1 || x > 1)
      NextUCx = UCx + NB->UC_Offset.x; if (NotValidUC(NextUCx)) continue;
      NextUCy = UCy + NB->UC_Offset.y; if (NotValidUC(NextUCy)) continue;
      NextUCz = UCz + NB->UC_Offset.z; if (NotValidUC(NextUCz)) continue;
      iNextSN = NB->NextSymN - LSymN->SymN;
#undef                                     NotValidUC

                iNextIC = IxIC(NextUCx, NextUCy, NextUCz, iNextSN);
      if (IsCon[iNextIC] == 0)
          IsCon[iNextIC] = -1;
    }

    IsCon[iIC] = 1;
  }

  Is3D = 0;

  iIC = 0;

  for (UCx = -1; UCx <= 1; UCx++)
  for (UCy = -1; UCy <= 1; UCy++)
  for (UCz = -1; UCz <= 1; UCz++)
  {
    nInCell = 0;

    for (iSN = 0; iSN < nSN; iSN++, iIC++)
      if (IsCon[iIC]) nInCell++;

    if (UCx || UCy || UCz) {
      if (nInCell == 0)   goto Clean;
    }
    else {
      if (nInCell != nSN) goto Clean;
    }
  }

  Is3D = 1;

  Clean:

  AppFree(IsCon, nIsCon);

  return Is3D;

#undef IxIC
}


static void BuildListSymNodes(const T_ListSymNodes *LSymN)
{
  int     iEPL, jEPL;
  int     nPos, iPos, jPos, kPos;
  T_RTMx  smxa, smxb, smxab;

  int                nPNNB, iPNNB, Sum_nPNNB, iLPNNB;
  T_PotentialNNBond  *PNNB;

  int             jAsyN, SymN_Offset;
  int             *AsyN_Offset;
  int             iSN;
  T_SymNodes      *SN;
  int             iNB;
  T_SymNodeBonds  *NB;
  int             j_m0pSh;
  T_fVector       fUC_Shift_ii;
  T_fVector       xjjb, xjkb, xjkb_xjk;

#include "mindist2.h"


  CheckMalloc(AsyN_Offset, LSymN->nAsyN);

  jAsyN = iSN = 0;

  for (iEPL = LSymN->Low_iEPL; iEPL <= LSymN->Top_iEPL; iEPL++) {
    if (LSymN->LnB[iEPL] != 0) {
      AsyN_Offset[jAsyN++] = iSN;
                             iSN += eD_PeakList[iEPL].WL_Entry->nPositions;
    }
  }

   SN = LSymN->SymN;
  iSN = 0;
   NB = LSymN->SymNB;
  iNB = 0;

  for (iEPL = LSymN->Low_iEPL; iEPL <= LSymN->Top_iEPL; iEPL++)
  {
    if (LSymN->LnB[iEPL] == 0) continue;

    nPos = eD_PeakList[iEPL].WL_Entry->nPositions;

    for (iPos = 0; iPos < nPos; iPos++)
    {
      Reconstruct_fUC_Shift(iEPL, iPos, &smxa, &fUC_Shift_ii);

      SN->nBonds = LSymN->LnB[iEPL];
      SN->NB     = NB;

      Sum_nPNNB = jAsyN = 0;

      for (jEPL = LSymN->Low_iEPL; jEPL <= LSymN->Top_iEPL; jEPL++)
      {
        if (LSymN->LnB[jEPL] == 0)
          continue;

         PNNB = NULL;
        nPNNB = 0;

        if (LSymN->Colored == 0 || NodeColor(iEPL) != NodeColor(jEPL))
        {
          for (iLPNNB = 0; iLPNNB < eD_PeakList[iEPL].lLPNNB; iLPNNB++)
          {
            if (jEPL == eD_PeakList[iEPL].LPNNB[iLPNNB].iEPL)
            {
               PNNB = eD_PeakList[iEPL].LPNNB[iLPNNB].PNNB;
              nPNNB = eD_PeakList[iEPL].LPNNB[iLPNNB].nPNNB;

              break;
            }
          }
        }

        for (iPNNB = 0; iPNNB < nPNNB; iPNNB++, PNNB++)
        {
          jPos = abs(PNNB->Info) - 1;
          j_m0pSh = jPos % 27;
          jPos /= 27;

          SeitzMx_of_iPos(jPos, eD_PeakList[jEPL].WL_Entry->WL_Flag,
                                                        &smxb, NULL);
          SeitzMxMultiply(&smxab, &smxa, &smxb);

          kPos = iPos_of_SeitzMx(&smxab, eD_PeakList[jEPL].WL_Entry->WL_Flag,
                                        &eD_PeakList[jEPL].Position);

          xjjb.x =   eD_PeakList[jEPL].xSymEquiv[jPos].x
                   + List_m0pShift[j_m0pSh].x / fSTBF;
          xjjb.y =   eD_PeakList[jEPL].xSymEquiv[jPos].y
                   + List_m0pShift[j_m0pSh].y / fSTBF;
          xjjb.z =   eD_PeakList[jEPL].xSymEquiv[jPos].z
                   + List_m0pShift[j_m0pSh].z / fSTBF;

          iRTMx_t_fVector(&xjkb, &smxa, &xjjb);

          xjkb.x += fUC_Shift_ii.x;
          xjkb.y += fUC_Shift_ii.y;
          xjkb.z += fUC_Shift_ii.z;

          xjkb_xjk.x = xjkb.x - eD_PeakList[jEPL].xSymEquiv[kPos].x;
          xjkb_xjk.y = xjkb.y - eD_PeakList[jEPL].xSymEquiv[kPos].y;
          xjkb_xjk.z = xjkb.z - eD_PeakList[jEPL].xSymEquiv[kPos].z;

          if (   IsUnitTr(xjkb_xjk.x) == 0
              || IsUnitTr(xjkb_xjk.y) == 0
              || IsUnitTr(xjkb_xjk.z) == 0)
            InternalError("Coseq(): Corrupt xjkb_xjk");

              SymN_Offset = AsyN_Offset[jAsyN] + kPos;
          if (SymN_Offset >= LSymN->nSymN)
            InternalError("Coseq(): SymN_Offset >= LSymN->nSymN");

          if (Sum_nPNNB >= SN->nBonds)
            InternalError("Coseq(): Sum_nPNNB >= SN->nBonds");

          if (iNB >= LSymN->nSymNB)
            InternalError("Coseq(): iNB >= LSymN->nSymNB");

          NB->NextSymN     = LSymN->SymN + SymN_Offset;
          NB->UC_Offset.x  = INT_ROUNDED(xjkb_xjk.x);
          NB->UC_Offset.y  = INT_ROUNDED(xjkb_xjk.y);
          NB->UC_Offset.z  = INT_ROUNDED(xjkb_xjk.z);

           NB++;
          iNB++;
          Sum_nPNNB++;
        }

        jAsyN++;
      }

      if (Sum_nPNNB != SN->nBonds)
        InternalError("Coseq(): Sum_nPNNB != SN->nBonds");

       SN++;
      iSN++;
    }
  }

  if (iSN != LSymN->nSymN) InternalError("Coseq(): iSN != LSymN->nSymN");

  if (F_CoseqSaveFileName)
  {
    Fprintf(stdout, ">Begin ListSymNodes\n");
    Fprintf(stdout, "#SymNodes = %d\n",     LSymN->nSymN);
    Fprintf(stdout, "#SymNodeBonds = %d\n", LSymN->nSymNB);

    for (iSN = 0; iSN < LSymN->nSymN; iSN++)
    {
      NB = LSymN->SymN[iSN].NB;

      for (iNB = 0; iNB < LSymN->SymN[iSN].nBonds; iNB++, NB++)
      {
        if (iNB == 0)
          Fprintf(stdout, " %5d", iSN);
        else
          Fprintf(stdout, "      ");

        Fprintf(stdout, "  %5ld  %2d %2d %2d\n",
        (long)(NB->NextSymN - LSymN->SymN),
        NB->UC_Offset.x,
        NB->UC_Offset.y,
        NB->UC_Offset.z);
      }
    }
    Fprintf(stdout, ">End ListSymNodes\n\n");
  }

  AppFree(AsyN_Offset, LSymN->nAsyN);
}


int CheckFragmentLoops(int Colored, const int *LnB, int Low_iEPL, int Top_iEPL)
{
  int             iEPL, LS, RejectedLoop;
  T_ListSymNodes  LSymN;


  LSymN.Colored  = Colored;
  LSymN.LnB      = LnB;
  LSymN.Low_iEPL = Low_iEPL;
  LSymN.Top_iEPL = Top_iEPL;

  LSymN.nAsyN  = 0;
  LSymN.nSymN  = 0;
  LSymN.nSymNB = 0;

  for (iEPL = Low_iEPL; iEPL <= Top_iEPL; iEPL++)
  {
    if (LnB[iEPL] != 0)
    {
      LSymN.nAsyN++;
      LSymN.nSymN  += eD_PeakList[iEPL].WL_Entry->nPositions;
      LSymN.nSymNB += eD_PeakList[iEPL].WL_Entry->nPositions * LnB[iEPL];
    }
  }

  CheckMalloc(LSymN.SymN,  LSymN.nSymN);
  CheckMalloc(LSymN.SymNB, LSymN.nSymNB);

  if (EvenLoopSizesOnly)
    LSymN.MaxLoopSize = MaxLoopSize;
  else
    LSymN.MaxLoopSize = MinLoopSize - 1;

  LSymN.LCf         = NULL;
  LSymN.LS_Flags    = NULL;

  BuildListSymNodes(&LSymN);

  BuildLoopConf(&LSymN);

  RejectedLoop = 0;

  for (LS = 1; LS < MinLoopSize; LS++) {
    if (LSymN.LS_Flags[LS]) {
      RejectedLoop = -1;
      goto Clean;
    }
  }

  if (EvenLoopSizesOnly)
  {
    for (LS = 1; LS <= LSymN.MaxLoopSize; LS += 2) {
      if (LSymN.LS_Flags[LS]) {
        RejectedLoop = -1;
        goto Clean;
      }
    }
  }

  Clean:

  FreeLoopConf(&LSymN);
  AppFree(LSymN.SymNB, LSymN.nSymNB);
  AppFree(LSymN.SymN,  LSymN.nSymN);

  return RejectedLoop;
}


void EvalFw(int Colored, const int *LnB, int Low_iEPL, int Top_iEPL)
{
  int             iEPL, Is3DimCon, LS;
  T_ListSymNodes  LSymN;


  if (   ModeCheckTetrahedra
      && CheckTetrahedra(Colored, LnB, Low_iEPL, Top_iEPL, NULL) != 0)
  {
    nBadTetrahedraFW++;

    if (F_nAllAsyPoints >= 0 || Debug0) { /* Debug: Print BadTetrahedraFW */
      Fprintf(stdout, "# BadTetrahedraFW");
      for (iEPL = Low_iEPL; iEPL <= Top_iEPL; iEPL++)
        if (LnB[iEPL]) Fprintf(stdout, " %d", iEPL);
      putc('\n', stdout);
      Fflush(stdout);
    }

    return;
  }

  LSymN.Colored  = Colored;
  LSymN.LnB      = LnB;
  LSymN.Low_iEPL = Low_iEPL;
  LSymN.Top_iEPL = Top_iEPL;

  LSymN.nAsyN  = 0;
  LSymN.nSymN  = 0;
  LSymN.nSymNB = 0;

  for (iEPL = Low_iEPL; iEPL <= Top_iEPL; iEPL++)
  {
    if (LnB[iEPL] != 0)
    {
      LSymN.nAsyN++;
      LSymN.nSymN  += eD_PeakList[iEPL].WL_Entry->nPositions;
      LSymN.nSymNB += eD_PeakList[iEPL].WL_Entry->nPositions * LnB[iEPL];
    }
  }

  CheckMalloc(LSymN.SymN,  LSymN.nSymN);
  CheckMalloc(LSymN.SymNB, LSymN.nSymNB);

  LSymN.MaxLoopSize = MaxLoopSize;
  LSymN.LCf         = NULL;
  LSymN.LS_Flags    = NULL;

  BuildListSymNodes(&LSymN);

  Is3DimCon = 0;

  if (Check3DimConnectivity)
  {
        Is3DimCon = ThreeDimConCheck(&LSymN);
    if (Is3DimCon == 0) {
      nNo3DimConFW++;
      goto Clean;
    }
  }

  BuildLoopConf(&LSymN);

  for (LS = 1; LS < MinLoopSize; LS++) {
    if (LSymN.LS_Flags[LS]) {
      FreeLoopConf(&LSymN);
      nSmallLoopsFW++;
      goto Clean;
    }
  }

  if (EvenLoopSizesOnly)
  {
    for (LS = 1; LS <= LSymN.MaxLoopSize; LS += 2) {
      if (LSymN.LS_Flags[LS]) {
        FreeLoopConf(&LSymN);
        nRejectedOddLoopsFW++;
        goto Clean;
      }
    }
  }

  PrintFramework(Colored, LnB, Low_iEPL, Top_iEPL, NULL);
  PrintLoopConf(&LSymN);
   FreeLoopConf(&LSymN);

  if (F_CoseqValue > 0)
  {
    CoseqRun(&LSymN);

    if (F_CoseqProtocolFileName)
      CoseqProtocol(-1, 0, 0); /* close protocol file (if open) */

    if (F_CoseqSaveFileName)
      AdminCoseqSave(NULL, 0); /* free file name buffer (if necessary) */
  }

  iFramework++;

  Clean:

  AppFree(LSymN.SymNB, LSymN.nSymNB);
  AppFree(LSymN.SymN,  LSymN.nSymN);
}
