/*
  Space Group Info (c) 1994-97 Ralf W. Grosse-Kunstleve
 */

#include <stdio.h>
#include <stdlib.h>


#include "sginfo.h"


static void ReverseCols(int nRows, int nCols, int *Mx)
{
  int  iRow, iCol, jCol, itmp;


  for (iRow = 0; iRow < nRows; iRow++, Mx += nCols)
    for (iCol = 0, jCol = nCols - 1; iCol < jCol; iCol++, jCol--) {
      itmp = Mx[iCol];
             Mx[iCol] = Mx[jCol];
                        Mx[jCol] = itmp;
    }
}


static void SetupMNp(int *MNp, int nCols, int *V, const T_RTMx *M, int iM)
{
  /*
    u0  u1  u2  u3  u4  u5  u6  u7  u8  v0  v1  v2       i  j  k        v
    ---------------------------------------------------------------------
    r0,  0,  0, r1,  0,  0, r2,  0,  0,  0,  0,  0, ..., 0, 0, 0, ... = 0
    r3,  0,  0, r4,  0,  0, r5,  0,  0,  0,  0,  0, ..., 0, 0, 0, ... = 0
    r6,  0,  0, r7,  0,  0, r8,  0,  0,  0,  0,  0, ..., 0, 0, 0, ... = 0
     0, r0,  0,  0, r1,  0,  0, r2,  0,  0,  0,  0, ..., 0, 0, 0, ... = 0
     0, r3,  0,  0, r4,  0,  0, r5,  0,  0,  0,  0, ..., 0, 0, 0, ... = 0
     0, r6,  0,  0, r7,  0,  0, r8,  0,  0,  0,  0, ..., 0, 0, 0, ... = 0
     0,  0, r0,  0,  0, r1,  0,  0, r2,  0,  0,  0, ..., 0, 0, 0, ... = 0
     0,  0, r3,  0,  0, r4,  0,  0, r5,  0,  0,  0, ..., 0, 0, 0, ... = 0
     0,  0, r6,  0,  0, r7,  0,  0, r8,  0,  0,  0, ..., 0, 0, 0, ... = 0
     0,  0,  0,  0,  0,  0,  0,  0,  0, r0, r1, r2, ..., 1, 0, 0, ... = -t0
     0,  0,  0,  0,  0,  0,  0,  0,  0, r3, r4, r5, ..., 0, 1, 0, ... = -t1
     0,  0,  0,  0,  0,  0,  0,  0,  0, r6, r7, r8, ..., 0, 0, 1, ... = -t2
   */

  int  ix, ib, ir, ic;


  for (ix = 0; ix < 12 * nCols; ix++) MNp[ix] = 0;
  for (ix = 0; ix < 12;         ix++)   V[ix] = 0;

  for (ib = 0; ib < 3; ib++)
    for (ir = 0; ir < 3; ir++)
      for (ic = 0; ic < 3; ic++)
        MNp[ib * (3 * nCols + 1) + ir * nCols + ic * 3] = M->s.R[ir * 3 + ic];

  for (ir = 0; ir < 3; ir++)
    for (ic = 0; ic < 3; ic++)
      MNp[(ir + 9) * nCols + 9 + ic] = M->s.R[ir * 3 + ic];

  for (ir = 0; ir < 3; ir++)
    MNp[(ir + 9) * nCols + 12 + iM * 3 + ir] = 1;

  for (ir = 0; ir < 3; ir++)
    V[ir + 9] = -M->s.T[ir];

  ReverseCols(12, nCols, MNp);

  MathMx(MNp, NULL, 0, 12, nCols, 2, 1);
  putc('\n', stdout);
  MathMx(V, NULL, 0, 1, 12, 2, 0);
  putc('\n', stdout);
}


static void SetupNpM(int *NpM, int nCols, int *V, const T_RTMx *M, int iM)
{
  /*
    u0  u1  u2  u3  u4  u5  u6  u7  u8  v0  v1  v2       i  j  k        v
    ---------------------------------------------------------------------
    r0, r3, r6,  0,  0,  0,  0,  0,  0,  0,  0,  0, ..., 0, 0, 0, ... = 0
     0,  0,  0, r0, r3, r6,  0,  0,  0,  0,  0,  0, ..., 0, 0, 0, ... = 0
     0,  0,  0,  0,  0,  0, r0, r3, r6,  0,  0,  0, ..., 0, 0, 0, ... = 0
    r1, r4, r7,  0,  0,  0,  0,  0,  0,  0,  0,  0, ..., 0, 0, 0, ... = 0
     0,  0,  0, r1, r4, r7,  0,  0,  0,  0,  0,  0, ..., 0, 0, 0, ... = 0
     0,  0,  0,  0,  0,  0, r1, r4, r7,  0,  0,  0, ..., 0, 0, 0, ... = 0
    r2, r5, r8,  0,  0,  0,  0,  0,  0,  0,  0,  0, ..., 0, 0, 0, ... = 0
     0,  0,  0, r2, r5, r8,  0,  0,  0,  0,  0,  0, ..., 0, 0, 0, ... = 0
     0,  0,  0,  0,  0,  0, r2, r5, r8,  0,  0,  0, ..., 0, 0, 0, ... = 0
    t0, t1, t2,  0,  0,  0,  0,  0,  0,  1,  0,  0, ..., 1, 0, 0, ... = 0
     0,  0,  0, t0, t1, t2,  0,  0,  0,  0,  1,  0, ..., 0, 1, 0, ... = 0
     0,  0,  0,  0,  0,  0, t0, t1, t2,  0,  0,  1, ..., 0, 0, 1, ... = 0
   */

  int  ix, ib, ir, ic;


  for (ix = 0; ix < 12 * nCols; ix++) NpM[ix] = 0;
  for (ix = 0; ix < 12;         ix++)   V[ix] = 0;

  for (ib = 0; ib < 3; ib++)
    for (ir = 0; ir < 4; ir++)
      for (ic = 0; ic < 3; ic++)
        NpM[ib * (nCols + 3) + ir * 3 * nCols + ic] = M->a[ir * 3 + ic];

  for (ir = 0; ir < 3; ir++)
    NpM[(ir + 9) * nCols + 9 + ir] = 1;

  for (ir = 0; ir < 3; ir++)
    NpM[(ir + 9) * nCols + 12 + iM * 3 + ir] = 1;

  ReverseCols(12, nCols, NpM);

  MathMx(NpM, NULL, 0, 12, nCols, 2, 1);
  putc('\n', stdout);
  MathMx(V, NULL, 0, 1, 12, 2, 0);
  putc('\n', stdout);
}


static void SetupAXp(int *AXp, int nCols, const int *A, int s)
{
  /*
    x0  x1  x2      is js ks
    ---------------------------------
    a0, a1, a2, ..., 1, 0, 0, ... = 0
    a3, a4, a5, ..., 0, 1, 0, ... = 0
    a6, a7, a8, ..., 0, 0, 1, ... = 0
   */

  int  ix, ir, ic;


  for (ix = 0; ix < 3 * nCols; ix++) AXp[ix] = 0;

  for (ir = 0; ir < 3; ir++)
    for (ic = 0; ic < 3; ic++)
      AXp[ir * nCols + ic] = A[ir * 3 + ic];

  for (ir = 0; ir < 3; ir++)
    for (ic = 0; ic < 3; ic++)
      if (A[ir * 3 + ic] != 0) {
        AXp[ir * nCols + 3 + 3 * s + ir] = 1;
        break;
      }

  ReverseCols(3, nCols, AXp);

  MathMx(AXp, NULL, 0, 3, nCols, 2, 1);
  putc('\n', stdout);
}


int sgjunk(const T_SgInfo *SgInfo)
{
  int           nS, iS, nRows, nCols, nRmIXp, *RmIXp, i;
  int           Stat, Rank;
  int           iList, iLoopInv, nLoopInv;
  int           RmI[9];
  const T_RTMx  *lsmx;
  T_RTMx        CBMx[1], InvCBMx[1], PrimSMx[1];


  RmIXp  = NULL;
  Stat = -1;

  for (i = 0; i < 9; i++) CBMx->s.R[i] = CRBF * SgInfo->CCMx_LP[i];
  for (i = 0; i < 3; i++) CBMx->s.T[i] = 0;

  if (InverseRTMx(CBMx, InvCBMx, CRBF) == 0) goto CleanExit;

  nS = SgInfo->OrderP;

  nRows = 3     * nS;
  nCols = 3 + 3 * nS;
  nRmIXp = nRows * nCols;

      RmIXp = malloc(nRmIXp * sizeof (*RmIXp));
  if (RmIXp == NULL) {
    SetSgError("Not enough core");
    goto CleanExit;
  }

  nLoopInv = Sg_nLoopInv(SgInfo);
  iS = 0;
  lsmx = SgInfo->ListSeitzMx;

  for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
  {
    if (CB_SMx(PrimSMx, CBMx, lsmx, InvCBMx) != 0) goto CleanExit;

    for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++, iS++)
    {
      SetRminusI(PrimSMx->s.R, RmI, iLoopInv);
      SetupAXp(&RmIXp[iS * 3 * nCols], nCols, RmI, iS);
    }
  }

  Rank = iReducedRowEchelon(RmIXp, nRows, nCols, NULL);
#ifdef JUNK
  ReverseCols(nRows, nCols, RmIXp);
#endif

  fprintf(stdout, "Rank = %d\n", Rank);
  MathMx(RmIXp, NULL, 0, nRows, nCols, 2, 1);
  putc('\n', stdout);

  Stat = 0;

  CleanExit:

  if (RmIXp) free(RmIXp);

  if (Stat != 0) SetSgError("Internal Error: sgjunk()");
  return Stat;
}
