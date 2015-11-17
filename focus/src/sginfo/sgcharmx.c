/*
  Space Group Info (c) 1994-97 Ralf W. Grosse-Kunstleve
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define SGCOREDEF__
#include "sginfo.h"


#define MemCmp(t, s, n)  memcmp((t), (s), (n) * sizeof (*(t)))


int iReducedRowEchelon(int *Mx, int nr, int nc, int *V, int RemoveGCD)
{
  /* M.B. Boisen, Jr. & G.V. Gibbs
     Mathematical Crystallography, Revised Edition 1990
     pp. 317-320
   */

  int  pr, pc, ir, ic, piv, mir, lcm, gcd, fp, fi;


#define Mx(i, j) Mx[i * nc + j]

  pc = 0;

  for (pr = 0; pr < nr; pr++, pc++)
  {
    for (; pc < nc; pc++)
      for (ir = pr; ir < nr; ir++)
        if ((piv = Mx(ir, pc)) != 0) goto DoRowOp;

    break;

    DoRowOp:

    if (ir != pr) {
      for (ic = 0; ic < nc; ic++) {
        fi = Mx(ir, ic);
             Mx(ir, ic) = Mx(pr, ic);
                          Mx(pr, ic) = fi;
      }
      if (V) { fi = V[ir]; V[ir] = V[pr]; V[pr] = fi; }
    }

    for (ir = 0; ir < nr; ir++)
    {
      if (ir == pr) continue;
      if ((mir = Mx(ir, pc)) == 0) continue;

      lcm = iLCM(piv, mir);
      fp = lcm / piv;
      fi = lcm / mir;

      for (ic = 0; ic < nc; ic++) {
        Mx(ir, ic) *= fi;
        Mx(ir, ic) -= fp * Mx(pr, ic);
      }

      if (V) {
        V[ir] *= fi;
        V[ir] -= fp * V[pr];
      }
    }
  }

  for (ir = 0; ir < pr; ir++)
  {
    if (RemoveGCD) {
      gcd = (V ? V[ir] : 0);
      for (ic = ir; ic < nc; ic++) gcd = iGCD(gcd, Mx(ir, ic));
    }
    else
      gcd = 1;

    for (ic = ir; ic < nc; ic++)
      if ((fi = Mx(ir, ic)) != 0) {
        if (fi < 0) gcd *= -1;
        break;
      }

    for (; ic < nc; ic++)
      Mx(ir, ic) /= gcd;

    if (V) V[ir] /= gcd;
  }

#undef Mx

  return pr;
}


int RRE2EigenVector(const int *Mx, int *EV)
{
  /* Mx must be in iReducedRowEchelon() form with Rank 2.
   */

  int  i;


  if      (Mx[0] == 0)               /*      [ 0 a 0 ]    x = 1        */
  {                                  /* Mx = [ 0 0 b ] => y = 0        */
    EV[0] = 1;                       /*      [ 0 0 0 ]    z = 0        */
    EV[1] = 0;
    EV[2] = 0;
  }
  else if (Mx[1] == 0 && Mx[4] != 0) /*      [ a 0 p ]    x = -p/a*z   */
  {                                  /* Mx = [ 0 b q ] => y = -q/b*z   */
    EV[2] = iLCM(Mx[0], Mx[4]);      /*      [ 0 0 0 ]    z = lcm(a,b) */
    EV[0] = (-Mx[2] * EV[2]) / Mx[0];
    EV[1] = (-Mx[5] * EV[2]) / Mx[4];
  }
  else if (Mx[5] != 0)               /*      [ a p 0 ]    x = -p/a*y   */
  {                                  /* Mx = [ 0 q b ] => y = lcm(a,b) */
    EV[1] = iLCM(Mx[0], Mx[5]);      /*      [ 0 0 0 ]    z = -q/b*y   */
    EV[2] = (-Mx[4] * EV[1]) / Mx[5];
    EV[0] = (-Mx[1] * EV[1]) / Mx[0];
  }
  else
  {
    EV[0] = 0;
    EV[1] = 0;
    EV[2] = 0;
    SetSgError("Internal Error: Corrupt iReducedRowEchelon Matrix");
    return -1;
  }

  if (SignHemisphere(EV[0], EV[1], EV[2]) < 0)
    for (i = 0; i < 3; i++) EV[i] *= -1;

  return 0;
}


static int RRE2Tr(const int *Mx, const int *wl, int *Tr)
{
  /* Mx must be in iReducedRowEchelon() form with Rank 2.
   */

  int  ix[3], iwl, i;


  if      (Mx[0] == 0)               /*        [ 0 a 0 | r ]    x = 0   */
  {                                  /* Mx|V = [ 0 0 b | s ] => y = r/a */
    ix[0] = -1;                      /*        [ 0 0 0 | 0 ]    z = s/b */
    ix[1] =  1;
    ix[2] =  5;
  }
  else if (Mx[1] == 0 && Mx[4] != 0) /*        [ a 0 p | r ]    x = r/a */
  {                                  /* Mx|V = [ 0 b q | s ] => y = s/b */
    ix[0] =  0;                      /*        [ 0 0 0 | 0 ]    z = 0   */
    ix[1] =  4;
    ix[2] = -1;
  }
  else if (Mx[5] != 0)               /*        [ a p 0 | r ]    x = r/a */
  {                                  /* Mx|V = [ 0 q b | s ] => y = 0   */
    ix[0] =  0;                      /*        [ 0 0 0 | 0 ]    z = s/b */
    ix[1] = -1;
    ix[2] =  5;
  }
  else
    goto ReturnInternalError;

  iwl = 0;

  for (i = 0; i < 3; i++)
  {
    if (ix[i] < 0)
      Tr[i] = 0;
    else
    {
      if (Mx[ix[i]] == 0) goto ReturnInternalError;

          Tr[i] = wl[iwl++];
      if (Tr[i] %  Mx[ix[i]]) return -1;
          Tr[i] /= Mx[ix[i]];
    }
  }

  return 0;

  ReturnInternalError:
  SetSgError("Internal Error: Corrupt iReducedRowEchelon Matrix");
  return -1;
}


static void RRE1BVAB(const int *Mx, const int *BVC, int *BVA, int *BVB,
                     const int *Mbv)
{
  int       iPerm, mPerm, Det, MinDet, f, i;
  const int rPerm[3][3] = {{ 0, 1, 2}, { 2, 0, 1}, { 1, 2, 0}};


#define pMx(ir)   Mx[rPerm[iPerm][ir]]
#define pBVA(ir) BVA[rPerm[iPerm][ir]]
#define pBVB(ir) BVB[rPerm[iPerm][ir]]
#define pBVC(ir) BVC[rPerm[iPerm][ir]]

  for (iPerm = 0;;)
  {
    if      (pMx(0) == 0)
    {
      if      (pMx(1) == 0)
      {
        f = 1; if (pBVC(2) < 0) f = -1;
        pBVA(0) = f;  pBVB(0) = 0;
        pBVA(1) = 0;  pBVB(1) = 1;
        pBVA(2) = 0;  pBVB(2) = 0;
        return;
      }
      else if (pMx(2) != 0)
      {
        f = 1; if (pMx(2) * pBVC(2) + pMx(1) * pBVC(1) < 0) f = -1;
        pBVA(0) = f;  pBVB(0) =       0;
        pBVA(1) = 0;  pBVB(1) =  pMx(2);
        pBVA(2) = 0;  pBVB(2) = -pMx(1);
        return;
      }

      iPerm = 1;
    }
    else if (pMx(1) == 0)
      iPerm = 2;
    else if (pMx(2) == 0)
      iPerm = 1;
    else
      break;
  }

  MinDet = 0;
  mPerm = 0;

  for (iPerm = 0; iPerm < 3; iPerm++)
  {
    if (Mbv == NULL)
    {
      Det =   pBVC(0) * pMx(0) * pMx(1)
            + pBVC(1) * pMx(1) * pMx(1)
            + pBVC(2) * pMx(1) * pMx(2);
    }
    else
    {
      pBVA(0) =  pMx(1);
      pBVA(1) = -pMx(0);
      pBVA(2) = 0;

      RotMx_t_Vector(BVB, Mbv, BVA, 0);

      Det =   BVC[0] * (BVA[1] * BVB[2] - BVA[2] * BVB[1])
            - BVC[1] * (BVA[0] * BVB[2] - BVA[2] * BVB[0])
            + BVC[2] * (BVA[0] * BVB[1] - BVA[1] * BVB[0]);
    }

    if ((abs(MinDet) > abs(Det) || MinDet == 0) && Det != 0) {
      MinDet = Det;
      mPerm = iPerm;
    }
  }

  iPerm = mPerm;

  pBVA(0) =  pMx(1);
  pBVA(1) = -pMx(0);
  pBVA(2) = 0;

  if (Mbv == NULL)
  {
    pBVB(0) = 0;
    pBVB(1) =  pMx(2);
    pBVB(2) = -pMx(1);
  }
  else
    RotMx_t_Vector(BVB, Mbv, BVA, 0);

  if (MinDet < 0)
    for (i = 0; i < 3; i++) BVA[i] *= -1;

#undef pMx
#undef pBVA
#undef pBVB
#undef pBVC

  return;
}


static int VerifyRotMxOrder(const int *ProperR, int AbsOrder, int *CumMx)
{
  int  MxA[9], MxB[9];
  int  *RR, *RRR, *Swp, iO, i;

  const int IdentityMx[] = { 1, 0, 0,
                             0, 1, 0,
                             0, 0, 1 };


  if (CumMx) for (i = 0; i < 9; i++) CumMx[i] = ProperR[i];

  RR = (int *) ProperR;
  RRR = MxA;

  for (iO = 1; iO < AbsOrder; iO++)
  {
    for (i = 0; i < 9; i++)
      if (-1 > RR[i] || RR[i] > 1)
        return -1;

    if (MemCmp(IdentityMx, RR, 9) == 0)
      return -1;

    RotMxMultiply(RRR, ProperR, RR);
    if (RR == ProperR) RR = MxB;
    Swp = RR; RR = RRR; RRR = Swp;

    if (CumMx) for (i = 0; i < 9; i++) CumMx[i] += RR[i];
  }

  if (MemCmp(IdentityMx, RR, 9) != 0)
    return -1;

  return 0;
}


static int SetBasis(const int ProperR[9], int AbsOrder, int CumMx[9],
                    int Basis[3][3])
{
  int        i, j;
  int        RmI[9], MbvBuf[9];
  const int  *Mbv;


  if (AbsOrder == 1)
  {
    for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
        Basis[i][j] = (i == j ? 1 : 0);

    return 0;
  }

  SetRminusI(ProperR, RmI, 0);
  if (iReducedRowEchelon(RmI, 3, 3, NULL, 1) != 2) return -1;
  if (RRE2EigenVector(RmI, Basis[2]) != 0) return -1;

  if (iReducedRowEchelon(CumMx, 3, 3, NULL, 1) != 1) return -1;

  Mbv = NULL;

  if (AbsOrder > 2)
  {
    if (AbsOrder > 4) {
      RotMxMultiply(MbvBuf, ProperR, ProperR);
      Mbv = MbvBuf;
    }
    else
      Mbv = ProperR;
  }

  RRE1BVAB(CumMx, Basis[2], Basis[0], Basis[1], Mbv);

  if (deterRotMx((int *) Basis) < 1) {
    SetSgError("Internal Error: SetBasis(): Det < 1");
    return -1;
  }

  return 0;
}


int SenseOfRotation(const int *R, int Order, const int *EV)
{
  /* M.B. Boisen, Jr. & G.V. Gibbs
     Mathematical Crystallography, Revised Edition 1990
     pp. 348-349, 354-356
   */

  int  f, trace;


  f = 1; if (Order < 0) f = -1;
  trace = f * (R[0] + R[4] + R[8]);

  if (trace == 3 || trace == -1) return 0; /* 1-fold or 2-fold */

  if (EV[1] == 0 && EV[2] == 0) {
    if (EV[0] * f * R[7] > 0)
      return 1;
  }
  else
    if (f * (R[3] * EV[2] - R[6] * EV[1]) > 0)
      return 1;

  return -1;
}


int GetRotMxOrder(const int *RotMx)
{
  int deter = deterRotMx(RotMx);

  if (deter == -1 || deter == 1)
  {
    switch (traceRotMx(RotMx))
    {
      case -3:                  return -1;
      case -2:                  return -6;
      case -1: if (deter == -1) return -4;
               else             return  2;
      case  0: if (deter == -1) return -3;
               else             return  3;
      case  1: if (deter == -1) return -2;
               else             return  4;
      case  2:                  return  6;
      case  3:                  return  1;
    }
  }

  return 0;
}


void ResetSeitzMxInfo(T_SMxI *SI, int KeepRMxI)
{
  int  i;


  if (KeepRMxI == 0)
  {
    SI->Order = 0;
    for (i = 0; i < 9; i++) SI->Basis[0][i] = 0;
    SI->SenseOfRotation = 0;
  }

  for (i = 0; i < 3; i++) SI->wg[i] = 0;
  for (i = 0; i < 3; i++) SI->Tr[i] = 0;
}


int SetRotMxInfo(const int *R, T_SMxI *SI)
{
  int        Order, AbsOrder, i;
  int        M_ProperR[9];
  const int   *ProperR;
  int        M_CumMx[9], *CumMx;


  if (SI) ResetSeitzMxInfo(SI, 0);

      Order = GetRotMxOrder(R);
  if (Order == 0)
    return 0;

  ProperR = R;

      AbsOrder = Order;
  if (AbsOrder < 0) {
      AbsOrder *= -1;
    for (i = 0; i < 9; i++) M_ProperR[i] = -R[i];
    ProperR = M_ProperR;
  }

  CumMx = NULL; if (SI != NULL) CumMx = M_CumMx;

  if (VerifyRotMxOrder(ProperR, AbsOrder, CumMx) != 0)
    return 0;

  if (SI)
  {
    if (SetBasis(ProperR, AbsOrder, CumMx, SI->Basis) != 0) {
      ResetSeitzMxInfo(SI, 0);
      return 0;
    }

    SI->SenseOfRotation = SenseOfRotation(R, Order, SI->Basis[2]);
    SI->Order = Order;
  }

  return Order;
}


int MakeCumRMx(const int *R, int Order, int *CumRMx)
{
  int  MxA[9], MxB[9];
  int  *RR, *RRR, *Swp, iO, i;


  if (Order < 0) {
    Order *= -1;
    if (Order % 2) Order *= 2;
  }

  InitRotMx(CumRMx, 1);

  if (Order > 1)
  {
    RR = (int *) R;
    RRR = MxA;

    for (iO = 1;;)
    {
      for (i = 0; i < 9; i++) CumRMx[i] += RR[i];

      if (++iO == Order)
        break;

      RotMxMultiply(RRR, R, RR);
      if (RR == R) RR = MxB;
      Swp = RR; RR = RRR; RRR = Swp;
    }
  }

  return Order;
}


int SetSeitzMxInfo(const T_RTMx *SMx, T_SMxI *SI)
{
  int  Mul, Mx[9], wl[3], Rank, i;


  if (SI->Order == 0) {
    if (SetRotMxInfo(SMx->s.R, SI) == 0)
      return 0;
  }
  else
    ResetSeitzMxInfo(SI, 1);

  Mul = MakeCumRMx(SMx->s.R, SI->Order, Mx);

  RotMx_t_Vector(SI->wg, Mx, SMx->s.T, 0);

  for (i = 0; i < 3; i++) {
    if (SI->wg[i] %  Mul) goto ReturnError;
        SI->wg[i] /= Mul;
  }

  for (i = 0; i < 3; i++)
    wl[i] = -(SMx->s.T[i] - SI->wg[i]) * (CTBF / STBF);

  SetRminusI(SMx->s.R, Mx, 0);
  Rank = iReducedRowEchelon(Mx, 3, 3, wl, 1);

  for (i = 2; i >= Rank; i--)
    if (iModPositive(wl[i], STBF) != 0) goto ReturnError;

  if      (Rank == 1)
  {
    for (i = 0; i < 3; i++) {
      if (Mx[i] != 0) break;
      SI->Tr[i] = 0;
    }

    if (i < 3) {
          SI->Tr[i] = wl[0];
      if (SI->Tr[i] %  Mx[i]) goto ReturnError;
          SI->Tr[i] /= Mx[i];
                 i++;
    }

    for (; i < 3; i++) SI->Tr[i] = 0;
  }
  else if (Rank == 2)
  {
    if (RRE2Tr(Mx, wl, SI->Tr) != 0)
      goto ReturnError;
  }
  else if (Rank == 3)
  {
    for (i = 0; i < 3; i++) {
          SI->Tr[i] = wl[i];
      if (SI->Tr[i] %  Mx[4 * i]) goto ReturnError;
          SI->Tr[i] /= Mx[4 * i];
    }
  }

  return SI->Order;

  ReturnError:
  ResetSeitzMxInfo(SI, 1);
  return 0;
}
