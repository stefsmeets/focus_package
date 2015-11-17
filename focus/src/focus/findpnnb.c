#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "main.h"
#include "xtal.h"


static Fprec EpsDeltaDist2 = 0.;


static void Identity_fSeitzMx_of_iSymOp(int iSymOp, T_fVector *Pos,
                                        T_fRTMx *fSMx)
{
  int            iList, iLoopInv, i;
  T_fVector      DiffPos;
  const T_fRTMx  *lfsmx;


  if (SpgrInfo.Centric == -1)
  {
    iList    = iSymOp / 2;
    iLoopInv = iSymOp % 2;
  }
  else
  {
    iList    = iSymOp;
    iLoopInv = 0;
  }

  lfsmx = &List_fSeitzMx[iList];

  if (iLoopInv == 0)
    for (i = 0; i < 12; i++) fSMx->a[i] =  lfsmx->a[i];
  else
    for (i = 0; i < 12; i++) fSMx->a[i] = -lfsmx->a[i];

  fRTMx_t_fVector(&DiffPos, fSMx, Pos);

  DiffPos.x -= Pos->x;
  DiffPos.y -= Pos->y;
  DiffPos.z -= Pos->z;

  fSMx->s.T[0] -= DiffPos.x;
  fSMx->s.T[1] -= DiffPos.y;
  fSMx->s.T[2] -= DiffPos.z;

#ifndef NO_EXTRA_CHECKS
  {
    int        nTrV, iTrV, TrialTrV[3];
    const int  *TrV;

#define Discretize(fTrV_)\
    ((int)((fTrV_) * fSTBF + (Fprec)((fTrV_) < CF0 ? -.5 : .5)))

    TrialTrV[0] = iModPositive(Discretize(DiffPos.x), STBF);
    TrialTrV[1] = iModPositive(Discretize(DiffPos.y), STBF);
    TrialTrV[2] = iModPositive(Discretize(DiffPos.z), STBF);

#undef Discretize

     TrV = SpgrInfo.LatticeInfo-> TrVector;
    nTrV = SpgrInfo.LatticeInfo->nTrVector;

    for (iTrV = 0; iTrV < nTrV; iTrV++, TrV += 3)
      if ( TrialTrV[0] == TrV[0]
        && TrialTrV[1] == TrV[1]
        && TrialTrV[2] == TrV[2]) break;

    if (iTrV == nTrV)
      InternalError("Cannot find TrV");

    if (   IsUnitTr(DiffPos.x - (Fprec) TrV[0] / fSTBF) == 0
        || IsUnitTr(DiffPos.y - (Fprec) TrV[1] / fSTBF) == 0
        || IsUnitTr(DiffPos.z - (Fprec) TrV[2] / fSTBF) == 0)
      InternalError("Cannot find TrV");
  }
#endif /* NO_EXTRA_CHECKS */
}


static void BuildListDist(int iEPL, int jEPL, T_PotentialNNBond *LD)
{
  int                j_nPos;
  Fprec              dx0, dy0, dz0, dx, dy, dz;
  T_fVector          *i_cSymEquiv, *j_cSymEquiv;
  T_PotentialNNBond  *LD0;

  Fprec  MinNodeDist2_minus;
  Fprec  MaxNodeDist2_plus;


  MinNodeDist2_minus = MinNodeDist2 * (Fprec)(1. - .02); /* ARBITRARY */
  MaxNodeDist2_plus  = MaxNodeDist2 * (Fprec)(1. + .02); /* ARBITRARY */

  LD0 = LD;

  i_cSymEquiv = eD_PeakList[iEPL].cSymEquiv;
  j_cSymEquiv = eD_PeakList[jEPL].cSymEquiv;
  j_nPos      = eD_PeakList[jEPL].WL_Entry->nPositions;

  while (j_nPos--)
  {
#define EvalDist2(dx_, dy_, dz_)\
  {\
    LD->Vect.x = dx_;\
    LD->Vect.y = dy_;\
    LD->Vect.z = dz_;\
    LD->Abs2 = dx_ * dx_ + dy_ * dy_ + dz_ * dz_;\
      LD->Info = 0;\
    if (   MinNodeDist2_minus > LD->Abs2\
        ||                      LD->Abs2 > MaxNodeDist2_plus)\
      LD->Info = -1;\
    LD++;\
  }

    dx0 = j_cSymEquiv->x - i_cSymEquiv->x;
    dy0 = j_cSymEquiv->y - i_cSymEquiv->y;
    dz0 = j_cSymEquiv->z - i_cSymEquiv->z;

    EvalDist2(dx0, dy0, dz0);

#include "mindist2.c" /* hand-made inline code */
#undef EvalDist2

    j_cSymEquiv++;
  }

  if (iEPL == jEPL) LD0->Info = -1;
}


static int VerifyPNNB(int iEPL, int jEPL, int j_iPos, int j_m0pSh,
                      T_PotentialNNBond **ListDist_ji)
{
  T_PotentialNNBond  *ListDist;
  int                i_nPos, i_iPos, i_m0pSh, iLD;
  T_RTMx             SMx_ii;
  T_fVector          fUC_Shift_ii;
  T_fVector          *xSymEquiv, ji, Pot_j0;
  Fprec              dx, dy, dz, Delta[3], DeltaDist2;

#include "mindist2.h"


  i_nPos = eD_PeakList[iEPL].WL_Entry->nPositions;

  if (*ListDist_ji == NULL)
  {
    CheckMalloc(*ListDist_ji, i_nPos * 27);

    BuildListDist(jEPL, iEPL, *ListDist_ji);
  }

  ListDist = *ListDist_ji;

  iLD = 0;

  for (i_iPos  = 0; i_iPos  < i_nPos; i_iPos++)
  for (i_m0pSh = 0; i_m0pSh <     27; i_m0pSh++, iLD++)
  {
    if (   ListDist[iLD].Info == 0
        && MinNodeDist2 <= ListDist[iLD].Abs2
        &&                 ListDist[iLD].Abs2 <= MaxNodeDist2)
    {
      Reconstruct_fUC_Shift(iEPL, i_iPos, &SMx_ii, &fUC_Shift_ii);

             xSymEquiv = &eD_PeakList[jEPL].xSymEquiv[j_iPos];
      ji.x = xSymEquiv->x + List_m0pShift[j_m0pSh].x / fSTBF;
      ji.y = xSymEquiv->y + List_m0pShift[j_m0pSh].y / fSTBF;
      ji.z = xSymEquiv->z + List_m0pShift[j_m0pSh].z / fSTBF;

      iRTMx_t_fVector(&Pot_j0, &SMx_ii, &ji);

      Pot_j0.x += fUC_Shift_ii.x + List_m0pShift[i_m0pSh].x / fSTBF;
      Pot_j0.y += fUC_Shift_ii.y + List_m0pShift[i_m0pSh].y / fSTBF;
      Pot_j0.z += fUC_Shift_ii.z + List_m0pShift[i_m0pSh].z / fSTBF;

      dx = eD_PeakList[jEPL].Position.x - Pot_j0.x;
      dy = eD_PeakList[jEPL].Position.y - Pot_j0.y;
      dz = eD_PeakList[jEPL].Position.z - Pot_j0.z;

      Tform_xc(dx, dy, dz, Delta[0], Delta[1], Delta[2]);

          DeltaDist2 = DotC(Delta, Delta);
      if (DeltaDist2 < EpsDeltaDist2)
        return 1;
    }
  }

  return 0;
}


static int Set_PNNB(int iEPL, int jEPL, T_ListPNNBonds *eLPNNB)
{
  int        j_nPos, j_iPos;
  int        j_m0pSh, jLD, k_m0pSh, kLD;
  Fprec      dx, dy, dz, Delta[3], DeltaDist2;
  int        iSymOp, k_iPos;
  T_fRTMx    fSMx;
  T_fVector  AsyBonded, Pot_ki;

  int                nPNNB;
  T_PotentialNNBond  *PNNB;
  T_PotentialNNBond   PNNB_Buf[MaxNCNmax];
  T_PotentialNNBond  *ListDist, *ListDist_ji;

#include "mindist2.h"


  j_nPos = eD_PeakList[jEPL].WL_Entry->nPositions;

  CheckMalloc(ListDist, j_nPos * 27);

  BuildListDist(iEPL, jEPL, ListDist);

  ListDist_ji = NULL;

   PNNB = PNNB_Buf;
  nPNNB = 0;

  jLD = 0;

  for (j_iPos  = 0; j_iPos  < j_nPos; j_iPos++)
  for (j_m0pSh = 0; j_m0pSh <     27; j_m0pSh++, jLD++)
  {
    if (   ListDist[jLD].Info == 0
        && MinNodeDist2 <= ListDist[jLD].Abs2
        &&                 ListDist[jLD].Abs2 <= MaxNodeDist2
        && (   iEPL == jEPL
            || VerifyPNNB(iEPL, jEPL, j_iPos, j_m0pSh, &ListDist_ji)))
    {
#define AddPNNB(iLD, SignInfo)\
      if (nPNNB == NCNmax) {\
        nPNNB = PNNB_TooMany;\
        goto Clean;\
      }\
      \
      PNNB->Info   = SignInfo(iLD + 1);\
      PNNB->Vect.x = ListDist[iLD].Vect.x;\
      PNNB->Vect.y = ListDist[iLD].Vect.y;\
      PNNB->Vect.z = ListDist[iLD].Vect.z;\
      PNNB->Abs2   = ListDist[iLD].Abs2;\
      PNNB++;\
      nPNNB++;\
      \
      ListDist[iLD].Info = 1

      AddPNNB(jLD, +);

      AsyBonded.x = eD_PeakList[jEPL].xSymEquiv[j_iPos].x
                               + List_m0pShift[j_m0pSh].x / fSTBF;
      AsyBonded.y = eD_PeakList[jEPL].xSymEquiv[j_iPos].y
                               + List_m0pShift[j_m0pSh].y / fSTBF;
      AsyBonded.z = eD_PeakList[jEPL].xSymEquiv[j_iPos].z
                               + List_m0pShift[j_m0pSh].z / fSTBF;

      /* skip first SymOp which is always the identity (and WL_F_Active)
       */
      for (iSymOp = 1; iSymOp < SpgrInfo.OrderP; iSymOp++)
      {
        if (eD_PeakList[iEPL].WL_Entry->WL_Flag[iSymOp] == WL_F_Identity)
        {
          Identity_fSeitzMx_of_iSymOp(
            iSymOp, &eD_PeakList[iEPL].xSymEquiv[0], &fSMx);

          fRTMx_t_fVector(&Pot_ki, &fSMx, &AsyBonded);

          kLD = 0;

          for (k_iPos  = 0; k_iPos  < j_nPos; k_iPos++)
          for (k_m0pSh = 0; k_m0pSh <     27; k_m0pSh++, kLD++)
          {
            if (ListDist[kLD].Info == 0)
            {
              dx = eD_PeakList[jEPL].xSymEquiv[k_iPos].x
                              + List_m0pShift[k_m0pSh].x / fSTBF;
              dy = eD_PeakList[jEPL].xSymEquiv[k_iPos].y
                              + List_m0pShift[k_m0pSh].y / fSTBF;
              dz = eD_PeakList[jEPL].xSymEquiv[k_iPos].z
                              + List_m0pShift[k_m0pSh].z / fSTBF;

              dx -= Pot_ki.x;
              dy -= Pot_ki.y;
              dz -= Pot_ki.z;

              Tform_xc(dx, dy, dz, Delta[0], Delta[1], Delta[2]);

                  DeltaDist2 = DotC(Delta, Delta);
              if (DeltaDist2 < EpsDeltaDist2) {
                AddPNNB(kLD, -);
              }
            }
          }
        }
      }
#undef AddPNNB
    }
  }

  Clean:

  if (ListDist_ji)
    AppFree(ListDist_ji, eD_PeakList[iEPL].WL_Entry->nPositions * 27);

  AppFree(ListDist, j_nPos * 27);

  if (nPNNB > 0)
  {
    CheckMalloc(PNNB, nPNNB);
    (void) memmove(PNNB, PNNB_Buf, sizeof (*PNNB) * nPNNB);
  }
  else
    PNNB = NULL;

  eLPNNB->iEPL  = jEPL;
  eLPNNB->nPNNB = nPNNB;
  eLPNNB->PNNB  =  PNNB;

  return nPNNB;
}


static void RemovePNNB(T_eD_PeakList *eEPL, int iLPNNB, int New_nPNNB)
{
  int  jLPNNB;


  if (eEPL->LPNNB[iLPNNB].nPNNB > 0)
    AppFree(eEPL->LPNNB[iLPNNB].PNNB, eEPL->LPNNB[iLPNNB].nPNNB);

  if (New_nPNNB == 0)
  {
    for (jLPNNB = iLPNNB + 1; jLPNNB < eEPL->nLPNNB; iLPNNB++, jLPNNB++)
    {
      eEPL->LPNNB[iLPNNB].iEPL  = eEPL->LPNNB[jLPNNB].iEPL;
      eEPL->LPNNB[iLPNNB].nPNNB = eEPL->LPNNB[jLPNNB].nPNNB;
      eEPL->LPNNB[iLPNNB].PNNB  = eEPL->LPNNB[jLPNNB].PNNB;
    }

    eEPL->nLPNNB--;
  }

  eEPL->LPNNB[iLPNNB].PNNB  = NULL;
  eEPL->LPNNB[iLPNNB].nPNNB = New_nPNNB;
}


static int TidyLPNNB(void)
{
  int            iEPL,  jEPL;
  T_eD_PeakList  *EPLi, *EPLj;
  int            iLPNNB, jLPNNB;
  int            Have_ij, Have_ji;
  int            NeedAnotherTidy, SumPNNB;


  NeedAnotherTidy = 0;

  EPLi = eD_PeakList;

  for (iEPL = 0; iEPL < nCeD_PeakList; iEPL++, EPLi++)
  {
    if (EPLi->nLPNNB > 0)
    {
      EPLj = eD_PeakList;

      for (jEPL = 0; jEPL < nCeD_PeakList; jEPL++, EPLj++)
      {
        if (iEPL == jEPL) continue;

        Have_ij = 0;

        for (iLPNNB = 0; iLPNNB < EPLi->nLPNNB; iLPNNB++)
          if (EPLi->LPNNB[iLPNNB].iEPL == jEPL) {
            Have_ij = 1;
            break;
        }

        Have_ji = 0;

        if (EPLj->nLPNNB > 0)
        {
          for (jLPNNB = 0; jLPNNB < EPLj->nLPNNB; jLPNNB++)
            if (EPLj->LPNNB[jLPNNB].iEPL == iEPL) {
              Have_ji = 1;
              break;
          }
        }

        if (Have_ij != Have_ji)
        {
          if (Have_ij) RemovePNNB(EPLi, iLPNNB, 0);
          else         RemovePNNB(EPLj, jLPNNB, 0);

          NeedAnotherTidy = 1;
        }
        else if (Have_ij)
        {
          Have_ij = (EPLi->LPNNB[iLPNNB].nPNNB > 0);
          Have_ji = (EPLj->LPNNB[jLPNNB].nPNNB > 0);

          if (Have_ij != Have_ji)
          {
            if (Have_ij) RemovePNNB(EPLi, iLPNNB, EPLj->LPNNB[jLPNNB].nPNNB);
            else         RemovePNNB(EPLj, jLPNNB, EPLi->LPNNB[iLPNNB].nPNNB);

            NeedAnotherTidy = 1;
          }
        }
      }

      SumPNNB = 0;

      for (iLPNNB = 0; iLPNNB < EPLi->nLPNNB; iLPNNB++)
        if (EPLi->LPNNB[iLPNNB].nPNNB > 0)
          SumPNNB += EPLi->LPNNB[iLPNNB].nPNNB;

      if (SumPNNB == 0 || (SumPNNB < NCNmin && F_AllSitesOn == 0))
      {
        for (iLPNNB = 0; iLPNNB < EPLi->nLPNNB; iLPNNB++)
          if (EPLi->LPNNB[iLPNNB].nPNNB > 0)
            AppFree(EPLi->LPNNB[iLPNNB].PNNB, EPLi->LPNNB[iLPNNB].nPNNB);

        AppFree(EPLi->LPNNB, EPLi->mLPNNB);

        EPLi->mLPNNB = 0;
        EPLi->nLPNNB = 0;
        EPLi->LPNNB  = NULL;

        NeedAnotherTidy = 1;
      }
    }
  }

  return NeedAnotherTidy;
}


static int LPNNBiSortFunction(const T_ListPNNBonds *a,
                              const T_ListPNNBonds *b)
{
  if (a->nPNNB > 0 && b->nPNNB <= 0) return -1;
  if (a->nPNNB < 0 && b->nPNNB >= 0) return  1;
  if (a->iEPL < b->iEPL) return -1;
  if (a->iEPL > b->iEPL) return  1;
  return 0;
}


static void SortLPNNB(void)
{
  int             iEPL, Offs;
  T_eD_PeakList   *EPLi;
  int             iLPNNB;
  T_ListPNNBonds  *LPNNBi;


  EPLi = eD_PeakList;

  for (iEPL = 0; iEPL < nCeD_PeakList; iEPL++, EPLi++)
  {
    if (EPLi->nLPNNB == 0) continue;

    LPNNBi = EPLi->LPNNB;

      Offs = 0;
    if (LPNNBi->iEPL == iEPL)
      Offs = 1;

    qsort((void *)(LPNNBi + Offs), EPLi->nLPNNB - Offs, sizeof (*LPNNBi),
          (SortFunction) LPNNBiSortFunction);

    for (iLPNNB = 0; iLPNNB < EPLi->nLPNNB; iLPNNB++, LPNNBi++)
      if (LPNNBi->nPNNB < 1) break;

    if (EPLi->lLPNNB != iLPNNB) InternalError("Corrupt EPLi->lLPNNB");
  }
}


static void SetReturnLPNNB(void)
{
  int             iEPL,  jEPL;
  T_eD_PeakList   *EPLi, *EPLj;
  int             iLPNNB,  jLPNNB;
  T_ListPNNBonds  *LPNNBi, *LPNNBj;


  EPLi = eD_PeakList;

  for (iEPL = 0; iEPL < nCeD_PeakList; iEPL++, EPLi++)
  {
    LPNNBi = EPLi->LPNNB;

    for (iLPNNB = 0; iLPNNB < EPLi->nLPNNB; iLPNNB++, LPNNBi++)
    {
          jEPL = LPNNBi->iEPL;
      if (jEPL < iEPL) continue;

      if (iEPL == jEPL)
      {
        LPNNBi->ReturnLPNNB = NULL;
        continue;
      }

      EPLj = &eD_PeakList[jEPL];

      LPNNBj = EPLj->LPNNB;

      for (jLPNNB = 0; jLPNNB < EPLj->nLPNNB; jLPNNB++, LPNNBj++)
      {
        if (LPNNBj->iEPL == iEPL)
        {
          LPNNBi->ReturnLPNNB = LPNNBj;
          LPNNBj->ReturnLPNNB = LPNNBi;
          break;
        }
      }

      if (jLPNNB == EPLj->nLPNNB)
        InternalError("jLPNNB == EPLj->nLPNNB");
    }
  }
}


void BuildLPNNB(void)
{
  int             iEPL, jEPL, ijEPL, i;
  T_eD_PeakList   *EPLi;
  Fprec           MinDist2;
  int             nLPNNB, iLPNNB, SumPNNB;
  T_ListPNNBonds  *LPNNB;
  int             *NoOiEPL;


      EpsDeltaDist2 = MinNodeDist2 * 1.e-6; /* ARBITRARY */
  if (EpsDeltaDist2 <= 0.)
      EpsDeltaDist2 = 1.e-6;

  EPLi = eD_PeakList;

  for (iEPL = 0; iEPL < nCeD_PeakList; iEPL++, EPLi++)
  {
    EPLi->mLPNNB = 0;
    EPLi->nLPNNB = 0;
    EPLi->LPNNB  = NULL;

    if (   NextPeakMx[iEPL * (NeD_PeakList + 1)] < MinNodeDist2
        || (F_AllSitesOn == 0 && EPLi->WL_Entry->CanBeNode == 0))
      EPLi->nLPNNB = -1;
  }

  CheckMalloc(LPNNB, nCeD_PeakList);

  EPLi = eD_PeakList;

  for (iEPL = 0; iEPL < nCeD_PeakList; iEPL++, EPLi++)
  {
    if (EPLi->nLPNNB >= 0)
    {
      SumPNNB = 0;
      nLPNNB = 0;

      for (ijEPL = 0; ijEPL < nCeD_PeakList; ijEPL++)
      {
        if      (ijEPL == 0)    jEPL = iEPL;
        else if (ijEPL <= iEPL) jEPL = ijEPL - 1;
        else                    jEPL = ijEPL;

        if (eD_PeakList[jEPL].nLPNNB >= 0)
        {
              MinDist2 = NextPeakMx[iEPL * NeD_PeakList + jEPL];
          if (MinDist2 <= MaxNodeDist2)
          {
            if (MinDist2 <= MinNodeDist2)
            {
              if (ijEPL == 0) {
                EPLi->nLPNNB = -1; /* self distance too small */
                goto Next_iEPL;
              }

              LPNNB[nLPNNB].nPNNB = PNNB_TooClose;
              LPNNB[nLPNNB].PNNB  = NULL;
              LPNNB[nLPNNB].iEPL  = jEPL;
                    nLPNNB++;
            }
            else if (Set_PNNB(iEPL, jEPL, &LPNNB[nLPNNB]))
            {
              if (ijEPL == 0 && LPNNB[nLPNNB].nPNNB < 0) {
                EPLi->nLPNNB = -1; /* too many self-bonds */
                goto Next_iEPL;
              }

              if (LPNNB[nLPNNB].nPNNB > 0)
                SumPNNB += LPNNB[nLPNNB].nPNNB;

              nLPNNB++;
            }
          }
        }
      }

      if (SumPNNB >= NCNmin || (F_AllSitesOn && nLPNNB > 0))
      {
        EPLi->mLPNNB = nLPNNB;
        EPLi->nLPNNB = nLPNNB;
        CheckMalloc(EPLi->LPNNB, nLPNNB);
        (void) memmove(EPLi->LPNNB, LPNNB, sizeof (*LPNNB) * nLPNNB);
      }
      else if (nLPNNB)
      {
        for (i = 0; i < nLPNNB; i++)
          if (LPNNB[i].nPNNB > 0)
            AppFree(LPNNB[i].PNNB, LPNNB[i].nPNNB);
      }

      Next_iEPL:;
    }
  }

  AppFree(LPNNB, nCeD_PeakList);

  while (TidyLPNNB());

  EPLi = eD_PeakList;

  for (iEPL = 0; iEPL < nCeD_PeakList; iEPL++, EPLi++)
  {
    if (EPLi->nLPNNB < 0)
        EPLi->nLPNNB = 0;

    EPLi->lLPNNB = 0;

    if (EPLi->nLPNNB == 0) continue;

    LPNNB = EPLi->LPNNB;

    for (iLPNNB = 0; iLPNNB < EPLi->nLPNNB; iLPNNB++, LPNNB++)
      if (LPNNB->nPNNB > 0) EPLi->lLPNNB++;

    if (EPLi->lLPNNB == 0) InternalError("EPLi->lLPNNB == 0");
  }

  Sort_eD_PeakList(PSE_lLPNNB, nCeD_PeakList);

  CheckMalloc(NoOiEPL, nCeD_PeakList);

  EPLi = eD_PeakList;

  for (iEPL = 0; iEPL < nCeD_PeakList; iEPL++, EPLi++)
    NoOiEPL[EPLi->Index] = iEPL;

  EPLi = eD_PeakList;

  for (iEPL = 0; iEPL < nCeD_PeakList; iEPL++, EPLi++)
  {
    LPNNB = EPLi->LPNNB;

    for (iLPNNB = 0; iLPNNB < EPLi->nLPNNB; iLPNNB++, LPNNB++)
      LPNNB->iEPL = NoOiEPL[LPNNB->iEPL];
  }

  AppFree(NoOiEPL, nCeD_PeakList);

  SortLPNNB();

  SetReturnLPNNB();

  SetAOS();
}


void FreeLPNNB(void)
{
  int            iEPL, i;
  T_eD_PeakList  *EPLi;


  EPLi = eD_PeakList;

  for (iEPL = 0; iEPL < nCeD_PeakList; iEPL++, EPLi++)
  {
    if (EPLi->mLPNNB > 0)
    {
      for (i = 0; i < EPLi->nLPNNB; i++)
        if (EPLi->LPNNB[i].nPNNB > 0)
          AppFree(EPLi->LPNNB[i].PNNB, EPLi->LPNNB[i].nPNNB);

      AppFree(EPLi->LPNNB, EPLi->mLPNNB);
    }

    EPLi->mLPNNB = 0;
    EPLi->nLPNNB = 0;
    EPLi->LPNNB  = NULL;
  }

  Sort_eD_PeakList(PSE_Index, nCeD_PeakList);
}


void PrintLPNNB(void)
{
  int               iEPL, iLPNNB, iPNNB;
  T_eD_PeakList     *EPLi;
  T_ListPNNBonds    *LPNNB;
  T_PotentialNNBond  *PNNB;


  Fprintf(stdout, ">Begin LPNNB\n");

  EPLi = eD_PeakList;

  for (iEPL = 0; iEPL < nCeD_PeakList; iEPL++, EPLi++)
  {
    Fprintf(stdout, "%d %d %d (%d)",
      iEPL, EPLi->lLPNNB, EPLi->nLPNNB, EPLi->Index);

    if (EPLi->IsAOS >= 0)
      Fprintf(stdout, " IsAOS %d", EPLi->IsAOS);

    putc('\n', stdout);

    if (EPLi->nLPNNB > 0)
    {
      LPNNB = EPLi->LPNNB;

      for (iLPNNB = 0; iLPNNB < EPLi->nLPNNB; iLPNNB++, LPNNB++)
      {
        Fprintf(stdout, "  %d %d", LPNNB->iEPL, LPNNB->nPNNB);

        PNNB = LPNNB->PNNB;

        for (iPNNB = 0; iPNNB < LPNNB->nPNNB; iPNNB++, PNNB++)
        {
          Fprintf(stdout, " %6.4f", AppSqrt(PNNB->Abs2));

          if (PNNB->Info < 0)
            putc('*', stdout);
          else if (iPNNB + 1 < LPNNB->nPNNB)
            putc(' ', stdout);
        }

        if (LPNNB->ReturnLPNNB)
          Fprintf(stdout, " -> %d\n", LPNNB->ReturnLPNNB->iEPL);
        else if (LPNNB->nPNNB > 0)
          Fprintf(stdout, " -><-\n");
      }
    }
  }

  Fprintf(stdout, ">End LPNNB\n\n");
  Fflush(stdout);
}
