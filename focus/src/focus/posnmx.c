#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "main.h"
#include "xtal.h"


void SeitzMx_of_iPos(int iPos, int *WL_Flag, T_RTMx *SMx, T_fRTMx *fSMx)
{
  int            nSymOp, nTrV;
  int            nActive, iSymOp, iList, iLoopInv, iTrV, i;
  const int      *TrV;
  const T_RTMx   *lsmx;
  const T_fRTMx  *lfsmx;


  nSymOp = SpgrInfo.OrderP;
  nTrV = SpgrInfo.LatticeInfo->nTrVector;

  nActive = iPos / nTrV;
  iTrV    = iPos % nTrV;

  for (iSymOp = 0; iSymOp < nSymOp; iSymOp++)
    if (*WL_Flag++ == WL_F_Active)
      if (nActive-- == 0)
        break;

  if (iSymOp == nSymOp)
    InternalError("iSymOp == nSymOp");

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

  if (SMx != NULL)
  {
    lsmx = &SpgrInfo.ListSeitzMx[iList];

    if (iLoopInv == 0)
      for (i = 0; i < 12; i++) SMx->a[i] =  lsmx->a[i];
    else
      for (i = 0; i < 12; i++) SMx->a[i] = -lsmx->a[i];

    TrV = &SpgrInfo.LatticeInfo->TrVector[iTrV * 3];

    for (i = 0; i < 3; i++)
      SMx->s.T[i] += TrV[i];
  }

  if (fSMx != NULL)
  {
    lfsmx = &List_fSeitzMx[iList];

    if (iLoopInv == 0)
      for (i = 0; i < 12; i++) fSMx->a[i] =  lfsmx->a[i];
    else
      for (i = 0; i < 12; i++) fSMx->a[i] = -lfsmx->a[i];

    fSMx->s.T[0] += fTrVector[iTrV].x;
    fSMx->s.T[1] += fTrVector[iTrV].y;
    fSMx->s.T[2] += fTrVector[iTrV].z;
  }
}


int IsUnitTr(Fprec Tr)
{
  Tr = AppFabs(Tr);
  Tr = AppFmod(Tr, CF1);

#define Eps 1.e-4 /* ARBITRARY */
  if (Tr < (Fprec)      Eps ) return 1;
  if (Tr > (Fprec)(1. - Eps)) return 1;
#undef Eps

  return 0;
}


int iPos_of_SeitzMx(T_RTMx *SMx, int *WL_Flag, T_fVector *Pos)
{
  int           iList, nLoopInv, iLoopInv, nTrV, iTrV;
  const T_RTMx  *lsmx;
  const int     *TrV;
  int           iSymOp, iPos, sign, i;
  T_RTMx        DfSMx;
  T_fVector     RT;


  nLoopInv = Sg_nLoopInv(&SpgrInfo);
  nTrV = SpgrInfo.LatticeInfo->nTrVector;

  lsmx = SpgrInfo.ListSeitzMx;

  iPos =
  iSymOp = 0;

  for (iList = 0; iList < SpgrInfo.nList; iList++, lsmx++)
  {
    for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
    {
      if (WL_Flag[iSymOp++] == WL_F_Active)
      {
        if (iLoopInv == 0) sign =  1;
        else               sign = -1;

        for (i = 0; i < 9; i++)
          if (SMx->s.R[i] != sign * lsmx->s.R[i]) break;

        if (i == 9) /* rotation parts match */
        {
          TrV = SpgrInfo.LatticeInfo->TrVector;

          for (iTrV = 0; iTrV < nTrV; iTrV++, TrV += 3)
          {
            for (i = 0; i < 3; i++)
              if (   SMx->s.T[i]
                  != iModPositive(sign * lsmx->s.T[i] + TrV[i], STBF))
                break;

            if (i == 3)
              return iPos + iTrV;
          }
        }

        /* maybe SMx is inactive: find appropriate active entry
         */

        for (i = 0; i < 12; i++)
          DfSMx.a[i] = SMx->a[i] - sign * lsmx->a[i];

        iRTMx_t_fVector(&RT, &DfSMx, Pos);

        for (iTrV = 0; iTrV < nTrV; iTrV++)
        {
          if (   IsUnitTr(RT.x - fTrVector[iTrV].x)
              && IsUnitTr(RT.y - fTrVector[iTrV].y)
              && IsUnitTr(RT.z - fTrVector[iTrV].z))
            return iPos + iTrV;
        }

        iPos += nTrV;
      }
    }
  }

  InternalError("iPos not found");
  return -1; /* never reached */
}


void Reconstruct_fUC_Shift(int iEPL, int iPos,
                           T_RTMx *SMx, T_fVector *fUC_Shift_ii)
{
  T_fVector  e, *xSymEquiv;


  SeitzMx_of_iPos(iPos, eD_PeakList[iEPL].WL_Entry->WL_Flag, SMx, NULL);

  if (iPos == 0)
  {
    fUC_Shift_ii->x = 0.;
    fUC_Shift_ii->y = 0.;
    fUC_Shift_ii->z = 0.;
  }
  else
  {
    iRTMx_t_fVector(&e, SMx, &eD_PeakList[iEPL].Position);

                      xSymEquiv = &eD_PeakList[iEPL].xSymEquiv[iPos];
    fUC_Shift_ii->x = xSymEquiv->x - e.x;
    fUC_Shift_ii->y = xSymEquiv->y - e.y;
    fUC_Shift_ii->z = xSymEquiv->z - e.z;

#if ! defined(NO_EXTRA_CHECKS)
    if (   IsUnitTr(fUC_Shift_ii->x) == 0
        || IsUnitTr(fUC_Shift_ii->y) == 0
        || IsUnitTr(fUC_Shift_ii->z) == 0)
      InternalError("Corrupt fUC_Shift_ii");
#endif
                        xSymEquiv = &eD_PeakList[iEPL].xSymEquiv[0];
    fUC_Shift_ii->x -= (xSymEquiv->x - eD_PeakList[iEPL].Position.x);
    fUC_Shift_ii->y -= (xSymEquiv->y - eD_PeakList[iEPL].Position.y);
    fUC_Shift_ii->z -= (xSymEquiv->z - eD_PeakList[iEPL].Position.z);
  }

#if ! defined(NO_EXTRA_CHECKS)
  {
    /* 08/31/99: this test is not really consistent with the
       code above.
       Running on Red Hat Linux 6.0 resulted in the
       internal error below to be issued.
       "Fix": NormOf() introduced after peak interpolation.
     */

    T_fVector  Pot_ii;
    Fprec      dx, dy, dz, Delta[3], DeltaDist2;

    iRTMx_t_fVector(&Pot_ii, SMx, &eD_PeakList[iEPL].Position);

    Pot_ii.x += fUC_Shift_ii->x;
    Pot_ii.y += fUC_Shift_ii->y;
    Pot_ii.z += fUC_Shift_ii->z;

    dx = eD_PeakList[iEPL].xSymEquiv[iPos].x - Pot_ii.x;
    dy = eD_PeakList[iEPL].xSymEquiv[iPos].y - Pot_ii.y;
    dz = eD_PeakList[iEPL].xSymEquiv[iPos].z - Pot_ii.z;

    Tform_xc(dx, dy, dz, Delta[0], Delta[1], Delta[2]);

        DeltaDist2 = DotC(Delta, Delta);
    if (DeltaDist2 >= 1.e-4) /* ARBITRARY */
    {
      Fprintf(stdout, "# Want %8.6g %8.6g %8.6g\n",
        eD_PeakList[iEPL].xSymEquiv[iPos].x,
        eD_PeakList[iEPL].xSymEquiv[iPos].y,
        eD_PeakList[iEPL].xSymEquiv[iPos].z);

      Fprintf(stdout, "# Have %8.6g %8.6g %8.6g\n",
        Pot_ii.x,
        Pot_ii.y,
        Pot_ii.z);

      InternalError("DeltaDist2 too large");
    }
  }
#endif
}
