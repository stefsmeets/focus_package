#include <stdio.h>
#include <stdlib.h>

#include "main.h"
#include "xtal.h"


#define iRotMx_t_fVector(RotMxV_, RotMx_, V_)\
  {\
    (RotMxV_)->x =   (RotMx_)[0] * (V_)->x\
                   + (RotMx_)[1] * (V_)->y\
                   + (RotMx_)[2] * (V_)->z;\
    (RotMxV_)->y =   (RotMx_)[3] * (V_)->x\
                   + (RotMx_)[4] * (V_)->y\
                   + (RotMx_)[5] * (V_)->z;\
    (RotMxV_)->z =   (RotMx_)[6] * (V_)->x\
                   + (RotMx_)[7] * (V_)->y\
                   + (RotMx_)[8] * (V_)->z;\
  }


static int IsAOS(T_fVector *ShV)
{
  int              iList, iLoopInv, nLoopInv, nTrV, iTrV, i;
  const T_RTMx     *lsmx;
  const T_fVector  *fTrV;
  int              R_I[9];
  T_fVector        R_I_ShV, VT;


  nLoopInv = Sg_nLoopInv(&SpgrInfo);
  nTrV = SpgrInfo.LatticeInfo->nTrVector;

  for (iList = 0; iList < SpgrInfo.nList; iList++)
  {
    lsmx = &SpgrInfo.ListSeitzMx[iList];

    for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
    {
      if (iLoopInv == 0)
        for (i = 0; i < 9; i++)
        {
          if (i % 4) R_I[i] =  lsmx->s.R[i];
          else       R_I[i] =  lsmx->s.R[i] - 1;
        }
      else
        for (i = 0; i < 9; i++)
        {
          if (i % 4) R_I[i] = -lsmx->s.R[i];
          else       R_I[i] = -lsmx->s.R[i] - 1;
        }

      iRotMx_t_fVector(&R_I_ShV, R_I, ShV);

      fTrV = fTrVector;

      for (iTrV = 0; iTrV < nTrV; iTrV++, fTrV++)
      {
        VT.x = R_I_ShV.x + fTrV->x;
        VT.y = R_I_ShV.y + fTrV->y;
        VT.z = R_I_ShV.z + fTrV->z;

        if (   IsUnitTr(VT.x)
            && IsUnitTr(VT.y)
            && IsUnitTr(VT.z))
          break;
      }

      if (iTrV == nTrV) return 0;
    }
  }

  return 1;
}


void SetAOS(void)
{
  int            iEPL,  jEPL;
  T_eD_PeakList  *EPLi, *EPLj;
  int            j_nPos, j_iPos;
  T_fVector      *xSymEquiv, ShV;


  EPLi = eD_PeakList;

  for (iEPL = 0; iEPL < nCeD_PeakList; iEPL++, EPLi++)
  {
    if (EPLi->IsAOS >= 0) continue;

    EPLj = EPLi + 1;

    for (jEPL = iEPL + 1; jEPL < nCeD_PeakList; jEPL++, EPLj++)
    {
      xSymEquiv = EPLj->xSymEquiv;
         j_nPos = EPLj->WL_Entry->nPositions;

      for (j_iPos = 0; j_iPos < j_nPos; j_iPos++, xSymEquiv++)
      {
        ShV.x = xSymEquiv->x - EPLi->Position.x;
        ShV.y = xSymEquiv->y - EPLi->Position.y;
        ShV.z = xSymEquiv->z - EPLi->Position.z;

        if (IsAOS(&ShV))
        {
          EPLj->IsAOS = iEPL;
          goto Next_iEPL;
        }
      }
    }

    Next_iEPL:;
  }
}
