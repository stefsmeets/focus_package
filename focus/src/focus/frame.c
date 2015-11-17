#include <stdio.h>
#include <stdlib.h>

#include "main.h"
#include "xtal.h"
#include "io.h"


typedef struct T_LnBC /* List n Bonds Chain */
  {
    int            *LnB;
    int            *WasPivot;
    struct T_LnBC  *Next;
  }
  T_LnBC;


static int Red_nPL = 0;


static int riEPL_SortFunction(const int *a, const int *b)
{
  if (*a < *b) return -1;
  if (*a > *b) return  1;

  return 0;
}


static void Collect_riEPL(int Li, int Ti, int *LnB, int nAsyN, int *riEPL)
{
  int  iEPL, iAsyN;


  iAsyN = 0;

  for (iEPL = Li; iEPL <= Ti; iEPL++)
  {
    if (LnB[iEPL] == 0) continue;
    if (iAsyN == nAsyN) InternalError("Corrupt LnB");

    riEPL[iAsyN] = eD_PeakList[iEPL].Index;
          iAsyN++;
  }

  if (iAsyN != nAsyN)
    InternalError("Corrupt LnB");

  qsort((void *) riEPL, nAsyN, sizeof (*riEPL),
        (SortFunction) riEPL_SortFunction);
}


static int CmpFwFiEPL(int La, int Ta, int *LnBa,
                      int Lb, int Tb, int *LnBb,
                      int nAsyN)
{
  int  Update, *riEPLa, *riEPLb, iAsyN;


  CheckMalloc(riEPLa, nAsyN);
  CheckMalloc(riEPLb, nAsyN);

  Collect_riEPL(La, Ta, LnBa, nAsyN, riEPLa);
  Collect_riEPL(Lb, Tb, LnBb, nAsyN, riEPLb);

  Update = 0;

  for (iAsyN = 0; iAsyN < nAsyN; iAsyN++)
  {
    if (riEPLa[iAsyN] < riEPLb[iAsyN])               break;
    if (riEPLa[iAsyN] > riEPLb[iAsyN]) { Update = 1; break; }
  }

  AppFree(riEPLa, nAsyN);
  AppFree(riEPLb, nAsyN);

  return Update;
}


static int CouldBecomeFw(int Li, int Ti, const int *LnB, int *NCpAU,
                         const int *CanNotBeTetrahedron)
{
  int  CN;


  for (; Li <= Ti; Li++)
    if (LnB[Li]) break;

  if (Li > Ti) return 1;

  for (CN = LnB[Li]; CN <= NCNmax; CN++)
  {
    if (CN == 4 && CanNotBeTetrahedron && CanNotBeTetrahedron[Li])
      continue;

    if (NCpAU[CN] && eD_PeakList[Li].WL_Entry->CanBeCN[CN])
    {
      NCpAU[CN]--;
      if (CouldBecomeFw(Li + 1, Ti, LnB, NCpAU, CanNotBeTetrahedron)) return 1;
      NCpAU[CN]++;
    }
  }

  return 0;
}


static void InitAllowedNCpAU(int *ANCpAU)
{
  int  CN, iNT;

  T_NodeType  *NT;


  for (CN = 0; CN <= MaxNCNmax; CN++) ANCpAU[CN] = 0;

  NT = NodeTypes;

  for (iNT = 0; iNT < nNodeTypes; iNT++, NT++)
  {
    if (NT->MaxNpAU == 0)
      ANCpAU[NT->CN] = Red_nPL;
    else
      ANCpAU[NT->CN] = NT->MaxNpAU;
  }
}


static int UpdateLargestFwFragment(int Colored,
                                   int *LnB, int Low_iEPL, int Top_iEPL,
                                   int nSymN)
{
  int    iEPL, *LnBi, Update;
  int    SumLnB, nAsyN;
  float  mBpN;
  int    *CanNotBeTetrahedron;
  int    NCpAU[MaxNCNmax + 1];


  nAsyN  = 0;
  SumLnB = 0;

  LnBi = &LnB[Low_iEPL];

  for (iEPL = Low_iEPL; iEPL <= Top_iEPL; iEPL++, LnBi++) {
    if (*LnBi) {
      nAsyN++;
      SumLnB += *LnBi;
    }
  }

  mBpN = (float) SumLnB / nAsyN;

  if (Debug0) /* Debug: Print FwFragment */
  {
    T_LargestFwFragment  FwF;
    FwF.nAsyN    = nAsyN;
    FwF.nSymN    = nSymN;
    FwF.mBpN     = mBpN;
    FwF.Low_iEPL = Low_iEPL;
    FwF.Top_iEPL = Top_iEPL;
    FwF.LnB      = LnB;
    PrintFramework(Colored, NULL, 0, 0, &FwF);
  }

  Update = 0;

  if      (LargestFwFragment.mBpN <  mBpN)
    Update = 1;
  else if (LargestFwFragment.mBpN == mBpN)
  {
    if      (LargestFwFragment.nAsyN <  nAsyN)
      Update = 1;
    else if (LargestFwFragment.nAsyN == nAsyN)
    {
      if      (LargestFwFragment.nSymN <  nSymN)
        Update = 1;
      else if (LargestFwFragment.nSymN == nSymN)
        Update = CmpFwFiEPL(LargestFwFragment.Low_iEPL,
                            LargestFwFragment.Top_iEPL,
                            LargestFwFragment.LnB,
                            Low_iEPL,
                            Top_iEPL,
                            LnB,
                            nAsyN);
    }
  }

  if (Update && ModeCheckTetrahedra != 0)
  {
    CheckMalloc(CanNotBeTetrahedron, Top_iEPL + 1);

    if (CheckTetrahedra(Colored, LnB, Low_iEPL, Top_iEPL,
                        CanNotBeTetrahedron) != 0)
    {
      InitAllowedNCpAU(NCpAU);

      if (CouldBecomeFw(Low_iEPL, Top_iEPL, LnB, NCpAU,
                        CanNotBeTetrahedron) == 0)
        Update = -1;
    }

    AppFree(CanNotBeTetrahedron, Top_iEPL + 1);

    if (Update == -1) return -1;
  }

  if (   (MinLoopSize > 3 || EvenLoopSizesOnly)
      && CheckFragmentLoops(Colored, LnB, Low_iEPL, Top_iEPL) != 0)
    return -1;

  if (Debug0) /* Debug: Indicate Update LargestFwFragment */
  {
    if (Update)
      Fprintf(stdout, "# Update LargestFwFragment\n");
    else
      Fprintf(stdout, "# Smaller than LargestFwFragment\n");
  }

  if (Update)
  {
    LargestFwFragment.nAsyN    = nAsyN;
    LargestFwFragment.nSymN    = nSymN;
    LargestFwFragment.mBpN     = mBpN;
    LargestFwFragment.Low_iEPL = Low_iEPL;
    LargestFwFragment.Top_iEPL = Top_iEPL;

    LnBi = &LnB[Low_iEPL];

    if (Colored == 0)
    {
      for (iEPL = Low_iEPL; iEPL <= Top_iEPL; iEPL++, LnBi++)
        LargestFwFragment.LnB[iEPL] = *LnBi;
    }
    else
    {
      for (iEPL = Low_iEPL; iEPL <= Top_iEPL; iEPL++, LnBi++)
        if (NodeColor(iEPL) == 0)
          LargestFwFragment.LnB[iEPL] =   *LnBi;
        else
          LargestFwFragment.LnB[iEPL] = -(*LnBi);
    }
  }

  return 0;
}


static void CompleteNodeConnectivity(int *AllowedNCpAU,
                                     int CNmin, int CNmax,
                                     T_LnBC *pLnBC,
                                     int Low_iEPL, int Top_iEPL,
                                     int nSymN,
                                     int pEPL)
{
  int    *LnB, *pLnB, *LnBj, pLnB_pEPL;
  int    *WasPivot, *pWasPivot;
  int    iLPNNB[MaxNCNmax], iLPNNB_iL, iLmin, iLnext, iL, jL;
  int    SumBonds[MaxNCNmax], pCN, aCN, NCpAU[MaxNCNmax + 1];
  int    iEPL, jEPL;
  int    NewTop_iEPL, New_nSymN, NewCNmin, NewCNmax;
  int    ValidFw, Colored;

  int             lLPNNBp, jLPNNB;
  T_ListPNNBonds  *LPNNBp, *LPNNBj;
  int             CanBeCNp[MaxNCNmax + 1];

  T_eD_PeakList  *EPLj;


  if (pLnBC->Next == NULL)
  {
    CheckMalloc(pLnBC->Next, 1);
    CheckMalloc(pLnBC->Next->LnB,      Red_nPL);
    CheckMalloc(pLnBC->Next->WasPivot, Red_nPL);
                pLnBC->Next->Next = NULL;
  }

  pLnB      = pLnBC->LnB;
   LnB      = pLnBC->Next->LnB;

  if (FwSearchMethod != FwSM_AltFwTracking)
    Colored = 0;
  else
    Colored = 1;

  pWasPivot = pLnBC->WasPivot;
   WasPivot = pLnBC->Next->WasPivot;

  for (iEPL = 0; iEPL < Low_iEPL; iEPL++)
    LnB[iEPL] = 0;

  for (iEPL = 0; iEPL < Low_iEPL; iEPL++)
    WasPivot[iEPL] = 0;

   LPNNBp = eD_PeakList[pEPL].LPNNB;
  lLPNNBp = eD_PeakList[pEPL].lLPNNB;

  for (pCN = NCNmin; pCN <= NCNmax; pCN++) /* copy to stack for efficiency */
    CanBeCNp[pCN] = eD_PeakList[pEPL].WL_Entry->CanBeCN[pCN];

  pLnB_pEPL = pLnB[pEPL];

  iLmin = 0;

  if (nSymN == 0)
  {
    if (LPNNBp[0].iEPL == pEPL && Colored == 0)
    {
      iLmin++;
      iLPNNB[1] = 0;
    }
    else
    {
      nSymN = eD_PeakList[pEPL].WL_Entry->nPositions;
      pLnB[pEPL] = -1;
    }
  }

  iLnext = iLmin;
  iL = 0;
  iLPNNB_iL = iLPNNB[0] = 0;

  for (;;)
  {
    iEPL = LPNNBp[iLPNNB_iL].iEPL;

    if (iEPL < Low_iEPL) goto NextCombination;
    if (pLnB[iEPL])      goto NextCombination;

    SumBonds[iL] = LPNNBp[iLPNNB_iL].nPNNB;

    if (iL == 0)
      SumBonds[iL] += pLnB_pEPL;
    else
      SumBonds[iL] += SumBonds[iL - 1];

    pCN = SumBonds[iL];

    if (pCN > CNmax) goto NextCombination;

    if (LargestFwFragment.LnB == NULL)
    {
      if (AllowedNCpAU[pCN] == 0) goto MoreBonds;
      if (CanBeCNp[pCN] == 0)     goto MoreBonds;
    }
    else if (pCN == CNmax)
    {
      if (CanBeCNp[pCN] == 0)     goto MoreBonds;
    }

    for (iEPL = Low_iEPL; iEPL < Red_nPL; iEPL++)
      LnB[iEPL] = pLnB[iEPL];

    for (jL = 0; jL <= iL; jL++)
      LnB[LPNNBp[iLPNNB[jL]].iEPL] = -1;

    if (Colored)
      for (jL = 0; jL <= iL; jL++)
        eD_PeakList[LPNNBp[iLPNNB[jL]].iEPL].Color = ! eD_PeakList[pEPL].Color;

    NewTop_iEPL = Top_iEPL;
    New_nSymN = nSymN;

    for (jL = 0; jL <= iL; jL++)
    {
      jEPL = LPNNBp[iLPNNB[jL]].iEPL;

      if (NewTop_iEPL < jEPL)
          NewTop_iEPL = jEPL;

      EPLj = &eD_PeakList[jEPL];

          New_nSymN += EPLj->WL_Entry->nPositions;
      if (New_nSymN > MaxSymNodes)
        goto NextCombination;

      LPNNBj = EPLj->LPNNB;

      for (jLPNNB = 0; jLPNNB < EPLj->nLPNNB; jLPNNB++, LPNNBj++)
      {
                iEPL = LPNNBj->iEPL;
        if (LnB[iEPL] == 0) continue;

        if (LPNNBj->nPNNB < 0) goto NextCombination;

        if (Colored && EPLj->Color == eD_PeakList[iEPL].Color) continue;

        if (pWasPivot[iEPL])   goto NextCombination;

        if (LPNNBj->ReturnLPNNB && pLnB[iEPL])
        {
          if (LnB[iEPL] == -1)
              LnB[iEPL] = 0;

              LnB[iEPL] += LPNNBj->ReturnLPNNB->nPNNB;
          if (LnB[iEPL] > CNmax)
            goto NextCombination;
        }

#ifndef NO_EXTRA_CHECKS
        if (LPNNBj->ReturnLPNNB == NULL && iEPL != jEPL)
          InternalError("Corrupt ReturnLPNNB");
#endif
        if (LnB[jEPL] == -1)
            LnB[jEPL] = 0;

            LnB[jEPL] += LPNNBj->nPNNB;
        if (LnB[jEPL] > CNmax)
          goto NextCombination;
      }
    }

#ifndef NO_EXTRA_CHECKS
    if (LnB[pEPL] != pCN) InternalError("Corrupt LnB");
#endif

    if (Debug0) /* Debug: Trace LnB */
    {
      int nN = 0; /* number of nodes (at beginning of site list) */
      static int CountOK = 0;
      int OK = 1;
      Fprintf(stdout, "LnB[%d]", pEPL);
      for (iEPL = 0; iEPL <= NewTop_iEPL; iEPL++) {
        Fprintf(stdout, " %d", LnB[iEPL]);
        if (eD_PeakList[iEPL].Index < nN) putc('*', stdout);
        if (LnB[iEPL] && eD_PeakList[iEPL].Index >= nN) OK = 0;
      }
      if (OK) Fprintf(stdout, " OK%d", ++CountOK);
      putc('\n', stdout);
      Fflush(stdout);
    }

    ValidFw = 0;

    if (New_nSymN >= MinSymNodes)
    {
      InitAllowedNCpAU(NCpAU);

      LnBj = &LnB[Low_iEPL];

      for (jEPL = Low_iEPL; jEPL <= NewTop_iEPL; jEPL++, LnBj++)
      {
        if (*LnBj == 0) continue;
        if (NCpAU[*LnBj]-- == 0) goto AfterEvalFw;
        if (eD_PeakList[jEPL].WL_Entry->CanBeCN[*LnBj] == 0) goto AfterEvalFw;
      }

      ValidFw = 1;

      EvalFw(Colored, LnB, Low_iEPL, NewTop_iEPL);

      if (F_SignalFileName) CheckSignalFile();
      if (QuitProgram) AppExit(1);

      AfterEvalFw:;
    }

    if (LargestFwFragment.LnB)
    {
      if (ValidFw == 0)
      {
        if (New_nSymN == MaxSymNodes) goto NextCombination;

        if (NCNmin != NCNmax)
        {
          InitAllowedNCpAU(NCpAU);

          if (CouldBecomeFw(Low_iEPL, NewTop_iEPL, LnB, NCpAU, NULL) == 0)
            goto NextCombination;
        }
      }

      if (UpdateLargestFwFragment(Colored,
                                  LnB, Low_iEPL, NewTop_iEPL, New_nSymN) != 0)
        goto NextCombination;
    }

    if (ModeCheckTetrahedra != 0)
    {
      int  *CanNotBeTetrahedron;

      CheckMalloc(CanNotBeTetrahedron, NewTop_iEPL + 1);

      if (CheckTetrahedra(Colored, LnB, Low_iEPL, NewTop_iEPL,
                          CanNotBeTetrahedron) != 0)
      {
        InitAllowedNCpAU(NCpAU);

        if (CouldBecomeFw(Low_iEPL, NewTop_iEPL, LnB, NCpAU,
                          CanNotBeTetrahedron) == 0) {
          AppFree(CanNotBeTetrahedron, NewTop_iEPL + 1);
          goto NextCombination;
        }
      }

      AppFree(CanNotBeTetrahedron, NewTop_iEPL + 1);
    }

    for (jEPL = Low_iEPL; jEPL < Red_nPL; jEPL++)
      WasPivot[jEPL] = pWasPivot[jEPL];

    WasPivot[pEPL] = 1;

    for (jEPL = Low_iEPL + 1; jEPL <= NewTop_iEPL; jEPL++)
    {
      if (WasPivot[jEPL] == 0 && LnB[jEPL] && LnB[jEPL] < CNmax)
      {
        aCN = pCN;

        if (LargestFwFragment.LnB)
          while (AllowedNCpAU[aCN] == 0) aCN++;

#ifndef NO_EXTRA_CHECKS
        if (aCN > NCNmax) InternalError("Corrupt aCN");
#endif
        AllowedNCpAU[aCN]--;

        NewCNmin = CNmin;
        NewCNmax = CNmax;
        while (NewCNmin <= CNmax && ! AllowedNCpAU[NewCNmin]) NewCNmin++;
        while (NewCNmax >= CNmin && ! AllowedNCpAU[NewCNmin]) NewCNmax--;

        if (NewCNmin <= NewCNmax)
          CompleteNodeConnectivity(AllowedNCpAU,
                                   NewCNmin, NewCNmax,
                                   pLnBC->Next,
                                   Low_iEPL, NewTop_iEPL,
                                   New_nSymN,
                                   jEPL);
        AllowedNCpAU[aCN]++;
        if (LnB[jEPL] < CNmin && LargestFwFragment.LnB == NULL) break;
        WasPivot[jEPL] = 1;
      }
    }

    MoreBonds:

    if (pCN == CNmax) goto NextCombination;

    iLnext = ++iL;
    iLPNNB[iLnext] = iLPNNB_iL;
    goto Ne_tCombination;

    NextCombination:

    iL = iLnext;

    Ne_tCombination:

        iLPNNB_iL = ++iLPNNB[iL];
    if (iLPNNB_iL == lLPNNBp)
    {
      if (iL == iLmin) break;
      iLnext = --iL;
      goto Ne_tCombination;
    }
  }
}


static void RemoveLnBC(T_LnBC *LnBC)
{
  T_LnBC  *LnBC_Next;


  while (LnBC)
  {
    LnBC_Next = LnBC->Next;
        AppFree(LnBC->LnB,      Red_nPL);
        AppFree(LnBC->WasPivot, Red_nPL);
        AppFree(LnBC, 1);

    LnBC = LnBC_Next;
  }
}


static void FwTracking(void)
{
  int     AllowedNCpAU[MaxNCNmax + 1];
  int     iEPL, fEPL, lEPL, sEPL, nSeedNodes;
  T_LnBC  *LnBC;


  if (F_nAllAsyPoints >= 0 || Debug0) /* Debug: Print Time FwTracking Start */
    PrintTicks(stdout, NULL, "# Time FwTracking Start  ", "\n");

  nSeedNodes = 0;

  for (iEPL = 0; iEPL < nCeD_PeakList; iEPL++)
  {
    if (eD_PeakList[iEPL].lLPNNB == 0)
      break;

    if (eD_PeakList[iEPL].IsAOS < 0)
      nSeedNodes++;
  }

  Red_nPL = iEPL;

  if (F_nAllAsyPoints >= 0 || Debug0) { /* Debug: Print Red_nPL */
    Fprintf(stdout, "# nCeD_PeakList = %d\n", nCeD_PeakList);
    Fprintf(stdout, "# Red_nPL       = %d\n", Red_nPL);
    Fprintf(stdout, "# nSeedNodes    = %d\n", nSeedNodes);
    putc('\n', stdout);
    Fflush(stdout);
  }

  fEPL = 0;
  lEPL = Red_nPL - 1;

  if (F_nAllAsyPoints > 0) fEPL = F_AllAsyPoints[0];
  if (F_nAllAsyPoints > 1) lEPL = F_AllAsyPoints[1];

  if (Red_nPL && (fEPL < Red_nPL || lEPL < Red_nPL))
  {
    InitAllowedNCpAU(AllowedNCpAU);

    CheckMalloc(LnBC, 1);
    CheckMalloc(LnBC->LnB,      Red_nPL);
    CheckMalloc(LnBC->WasPivot, Red_nPL);
                LnBC->Next = NULL;

    for (iEPL = 0; iEPL < Red_nPL; iEPL++) LnBC->WasPivot[iEPL] = 0;
    for (iEPL = 0; iEPL < Red_nPL; iEPL++) LnBC->LnB[iEPL]      = 0;

    if (fEPL >= Red_nPL) fEPL = Red_nPL - 1;
    if (lEPL >= Red_nPL) lEPL = Red_nPL - 1;

    if (fEPL <= lEPL) sEPL = 1; else sEPL = -1;

    for (iEPL = fEPL; iEPL != lEPL + sEPL; iEPL += sEPL)
    {
      if (eD_PeakList[iEPL].IsAOS >= 0) continue;

      if (FwSearchMethod == FwSM_AltFwTracking)
        eD_PeakList[iEPL].Color = 0;

      if (F_nAllAsyPoints >= 0) {
        Fprintf(stdout, "# CompleteNodeConnectivity[%d] ", iEPL);
        PrintTicks(stdout, NULL, NULL, "\n");
        Fflush(stdout);
      }

      CompleteNodeConnectivity(AllowedNCpAU,
                               NCNmin, NCNmax,
                               LnBC,
                               iEPL, iEPL,
                               0,
                               iEPL);
      LnBC->LnB[iEPL] = 0;
    }

    RemoveLnBC(LnBC);
  }

  Red_nPL = 0;

  if (F_nAllAsyPoints >= 0 || Debug0) /* Debug: Print Time FwTracking End */
    PrintTicks(stdout, NULL, "# Time FwTracking End  ", "\n");
}


static void AllSitesOn(void)
{
  int  iEPL, iLPNNB, nPNNB;
  int  *List_nBonds;


  CheckMalloc(List_nBonds, NeD_PeakList);

  for (iEPL = 0; iEPL < NeD_PeakList; iEPL++)
  {
    List_nBonds[iEPL] = 0;

    for (iLPNNB = 0; iLPNNB < eD_PeakList[iEPL].lLPNNB; iLPNNB++)
    {
          nPNNB = eD_PeakList[iEPL].LPNNB[iLPNNB].nPNNB;
      if (nPNNB > 0)
        List_nBonds[iEPL] += nPNNB;
    }
  }

  EvalFw(0, List_nBonds, 0, NeD_PeakList - 1);

  AppFree(List_nBonds, NeD_PeakList);
}


static void ReIndexLargestFwFragment(void)
{
  int            iEPL, *LnBi, *RLnB;
  T_eD_PeakList  *EPLi;


  for (iEPL = 0; iEPL < LargestFwFragment.Low_iEPL; iEPL++)
    LargestFwFragment.LnB[iEPL] = 0;

  for (iEPL = LargestFwFragment.Top_iEPL + 1; iEPL < nCeD_PeakList; iEPL++)
    LargestFwFragment.LnB[iEPL] = 0;

  CheckMalloc(RLnB, nCeD_PeakList);

  EPLi = eD_PeakList;
  LnBi = LargestFwFragment.LnB;

  for (iEPL = 0; iEPL < nCeD_PeakList; iEPL++, EPLi++, LnBi++)
    RLnB[EPLi->Index] = *LnBi;

  LnBi = LargestFwFragment.LnB;

  for (iEPL = 0; iEPL < nCeD_PeakList; iEPL++, LnBi++)
    *LnBi = RLnB[iEPL];

  AppFree(RLnB, nCeD_PeakList);

  for (iEPL = 0; iEPL < nCeD_PeakList; iEPL++)
    if (LargestFwFragment.LnB[iEPL]) {
      LargestFwFragment.Low_iEPL = iEPL;
      break;
    }

  for (iEPL = nCeD_PeakList - 1; iEPL >= 0; iEPL--)
    if (LargestFwFragment.LnB[iEPL]) {
      LargestFwFragment.Top_iEPL = iEPL;
      break;
    }
}


void FwSearch(int FindLargestFwFragment)
{
  size_m  PrevMemoryUsed, MU;


  PrevMemoryUsed = MemoryUsed;

  (void) GetTicks(&FwSearchTimeStart);

  IndicateFw = (F_nAllAsyPoints >= 0 || F_SiteFrame == 0);

  if (F_nAllAsyPoints >= 0) {
    Fprintf(stdout, "# MemoryUsed before BuildLPNNB(): %ld\n",
      (long) MemoryUsed);
    Fflush(stdout);
    MU = MemoryUsed;
  }

  BuildLPNNB();

  if (F_nAllAsyPoints >= 0) {
    Fprintf(stdout, "# MemoryUsed after  BuildLPNNB(): %ld = %ld + %ld\n\n",
      (long) MemoryUsed, (long) MU, (long)(MemoryUsed - MU));
    Fflush(stdout);
  }

  if (F_nAllAsyPoints >= 0 || F_SiteFrame || Debug0) /* Debug: PrintLPNNB */
    PrintLPNNB();

  if (F_SignalFileName) CheckSignalFile();
  if (QuitProgram) AppExit(1);

  if (F_AllSitesOn)
    AllSitesOn();
  else
    FwTracking();

  if (LargestFwFragment.nAsyN) ReIndexLargestFwFragment();

  FreeLPNNB();

  if (F_ShowLargestFwFragment)
  {
    if (LargestFwFragment.nAsyN) {
      Fprintf(stdout, "# LargestFwFragment\n");
      PrintFramework(0, NULL, 0, 0, &LargestFwFragment);
    }
    else if (FindLargestFwFragment)
      Fprintf(stdout, "# No LargestFwFragment\n");
  }

  (void) GetTicks(&FwSearchTimeEnd);

  FwSearchSumTicks +=   FwSearchTimeEnd.cpu_user
                      - FwSearchTimeStart.cpu_user;
  FwSearchSumTicks +=   FwSearchTimeEnd.cpu_system
                      - FwSearchTimeStart.cpu_system;

  FwSearchTimeStart.second = -1;

  if (MemoryUsed != PrevMemoryUsed)
  {
    Fprintf(stdout, "# PrevMemoryUsed = %ld MemoryUsed = %ld\n",
                    (long) PrevMemoryUsed, (long) MemoryUsed);

    InternalError("Corrupt memory management");
  }
}
