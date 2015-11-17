/*
  Space Group Info (c) 1994-97 Ralf W. Grosse-Kunstleve
 */

#include <stdio.h>
#include <stdlib.h>


#define SGCOREDEF__
#include "sginfo.h"


#define Fprintf (void) fprintf


int iModShort(int ix, int iy)
{
      ix = iModPositive(ix, iy);
  if (ix > iy / 2)
      ix -= iy;

  return ix;
}


int SignHemisphere(int h, int k, int l)
{
  if (l >  0) return  1;
  if (l == 0) {
    if (k >  0) return  1;
    if (k == 0) {
      if (h >  0) return  1;
      if (h == 0)
        return 0;
    }
  }

  return -1;
}


int iGCD(const int a, const int b)
{
  int  ri, rj, rk;


      ri = a;
  if (ri < 0) ri = -ri;

  if ((rj = b) != 0)
  {
    for (;;)
    {
      rk = ri % rj; if (rk == 0) { ri = rj; break; }
      ri = rj % rk; if (ri == 0) { ri = rk; break; }
      rj = rk % ri; if (rj == 0) {          break; }
    }

    if (ri < 0) ri = -ri;
  }

  return ri;
}


int FindGCD(const int *S, int nS)
{
  int  ri, rj, rk;


  if (nS-- == 0) return 0;

      ri = *S++;
  if (ri < 0) ri = -ri;

  while (nS--)
  {
    if ((rj = *S++) != 0)
    {
      for (;;)
      {
        rk = ri % rj; if (rk == 0) { ri = rj; break; }
        ri = rj % rk; if (ri == 0) { ri = rk; break; }
        rj = rk % ri; if (rj == 0) {          break; }
      }

      if (ri < 0) ri = -ri;

      if (ri == 1) break;
    }
  }

  return ri;
}


int iLCM(const int a, const int b)
{
  int  ri, rj, rk;


      ri = a;
  if (ri == 0) ri = 1;

  if ((rj = b) != 0)
  {
    for (;;)
    {
      rk = ri % rj; if (rk == 0) { ri = rj; break; }
      ri = rj % rk; if (ri == 0) { ri = rk; break; }
      rj = rk % ri; if (rj == 0) {          break; }
    }

    ri = a / ri * b;
  }

  if (ri < 0) return -ri;
              return  ri;
}


int FindLCM(const int *S, int nS)
{
  int  a, b, ri, rj, rk;


  if (nS-- == 0) return 1;

      ri = *S++;
  if (ri == 0) ri = 1;

  while (nS--)
  {
    if ((rj = *S++) != 0)
    {
      a = ri; b = rj;

      for (;;)
      {
        rk = ri % rj; if (rk == 0) { ri = rj; break; }
        ri = rj % rk; if (ri == 0) { ri = rk; break; }
        rj = rk % ri; if (rj == 0) {          break; }
      }

      ri = a / ri * b;
    }
  }

  if (ri < 0) return -ri;
              return  ri;
}


void SimplifyFraction(int nume, int deno, int *o_nume, int *o_deno)
{
  int gcd = iGCD(nume, deno);
  if (gcd)
  {
    *o_nume = nume / gcd;
    *o_deno = deno / gcd;

    if (*o_deno < 0) {
      *o_nume *= -1;
      *o_deno *= -1;
    }
  }
}


const int *PseudoMetricalMatrix(int XtalSystem,
                                int UniqueDirCode, int UniqueRefAxis,
                                int Reciprocal)
{
  static const int PseudoGzD[] = {  2, -1,  0,  /* a = b = c = sqrt(2) */
                                   -1,  2,  0,  /* alpha = beta = 90   */
                                    0,  0,  2   /* gamma = 120         */
                                 };
  static const int PseudoGzR[] = {  2,  1,  0,  /* a = b = c = sqrt(2) */
                                    1,  2,  0,  /* alpha = beta = 90   */
                                    0,  0,  2   /* gamma =  60         */
                                 };
  static const int PseudoGyD[] = {  2,  0, -1,
                                    0,  2,  0,
                                   -1,  0,  2
                                 };
  static const int PseudoGyR[] = {  2,  0,  1,
                                    0,  2,  0,
                                    1,  0,  2
                                 };
  static const int PseudoGxD[] = {  2,  0,  0,
                                    0,  2, -1,
                                    0, -1,  2
                                 };
  static const int PseudoGxR[] = {  2,  0,  0,
                                    0,  2,  1,
                                    0,  1,  2
                                 };
  const int  *PseudoG = NULL;

  if (   (XtalSystem == XS_Trigonal && UniqueDirCode != '*')
      ||  XtalSystem == XS_Hexagonal)
  {
    if      (UniqueRefAxis == 'z') {
      if (Reciprocal == 0) PseudoG = PseudoGzD;
      else                 PseudoG = PseudoGzR;
    }
    else if (UniqueRefAxis == 'y') {
      if (Reciprocal == 0) PseudoG = PseudoGyD;
      else                 PseudoG = PseudoGyR;
    }
    else if (UniqueRefAxis == 'x') {
      if (Reciprocal == 0) PseudoG = PseudoGxD;
      else                 PseudoG = PseudoGxR;
    }
  }

  return PseudoG;
}


#define GmulV(GV, G, V) \
  (GV)[0] = (G)[0] * (V)[0] + (G)[1] * (V)[1] + (G)[2] * (V)[2]; \
  (GV)[1] = (G)[3] * (V)[0] + (G)[4] * (V)[1] + (G)[5] * (V)[2]; \
  (GV)[2] = (G)[6] * (V)[0] + (G)[7] * (V)[1] + (G)[8] * (V)[2]


int iScalProd(const int *u, const int *v, const int *PseudoG)
{
  int  Prod, Gv[3];


  if (PseudoG)
  {
    GmulV(Gv, PseudoG, v);
    v = Gv;
  }

  Prod =   u[0] * v[0]
         + u[1] * v[1]
         + u[2] * v[2];

  return Prod;
}


void iCrossProd(int *rxs, const int *r, const int *s, const int *PseudoG)
{
  int  Gr[3], Gs[3];


  if (PseudoG)
  {
    GmulV(Gr, PseudoG, r);
    GmulV(Gs, PseudoG, s);
    r = Gr;
    s = Gs;
  }

  rxs[0] = r[1] * s[2] - r[2] * s[1];
  rxs[1] = r[2] * s[0] - r[0] * s[2];
  rxs[2] = r[0] * s[1] - r[1] * s[0];
}


#undef GmulV


int ExpandListSeitzMx(T_SgInfo *SgInfo)
{
  int           iList, f, i;
  int           nTrV, iTrV, nLoopInv, iLoopInv;
  const int     *TrV;
  const T_RTMx  *lsmx;
  T_RTMx        *tsmx;
  T_RotMxInfo   *trmxi;


  if (SgInfo->MaxList < SgInfo->OrderL) {
    SetSgError(
      "Internal Error: ExpandListSeitzMx(): SgInfo->MaxList < SgInfo->OrderL");
    return -1;
  }

  tsmx = &SgInfo->ListSeitzMx[SgInfo->nList];

  nLoopInv = Sg_nLoopInv(SgInfo);

  nTrV = SgInfo->LatticeInfo->nTrVector;
   TrV = SgInfo->LatticeInfo->TrVector;

  for (iTrV = 0; iTrV < nTrV; iTrV++, TrV += 3)
  {
    for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
    {
      if (iTrV == 0 && iLoopInv == 0)
        continue;

      if (iLoopInv == 0) f =  1;
      else               f = -1;

      lsmx = SgInfo->ListSeitzMx;

      for (iList = 0; iList < SgInfo->nList; iList++, lsmx++, tsmx++)
      {
        for (i = 0; i < 9; i++)
          tsmx->s.R[i] =              f * lsmx->s.R[i];

        for (i = 0; i < 3; i++)
          tsmx->s.T[i] = iModPositive(f * lsmx->s.T[i] + TrV[i], STBF);
      }
    }
  }

  if (SgInfo->ListRotMxInfo)
  {
    iList = SgInfo->nList;

    tsmx  = &SgInfo->ListSeitzMx[iList];
    trmxi = &SgInfo->ListRotMxInfo[iList];

    for (; iList < SgInfo->OrderL; iList++, tsmx++, trmxi++) {
      if (GetRotMxInfo(tsmx->s.R, trmxi, NULL) == 0) {
        SetSgError("Internal Error: ExpandListSeitzMx(): Corrupt SeitzMx");
        return -1;
      }
    }
  }

  return 0;
}


void SetRminusI(const int *R, int *RmI, int Inv)
{
  int  i;


  if (Inv == 0)
    for (i = 0; i < 9; i++)
      RmI[i] =  R[i];
  else
    for (i = 0; i < 9; i++)
      RmI[i] = -R[i];

  for (i = 0; i < 9; i += 4)
    RmI[i] -= 1;
}


int BuildCumRMx(const T_SgInfo *SgInfo, int iList, int iLoopInv, int *CumRMx)
{
  int           nO, iO, i;
  const T_RTMx  *lsmx;
  int           *rmxab, rmxa[9], *rmxb, BufMx1[9], BufMx2[9];


  lsmx = &SgInfo->ListSeitzMx[iList];

  if (SgInfo->ListRotMxInfo)
    nO = SgInfo->ListRotMxInfo[iList].Order;
  else
    nO = GetRotMxOrder(lsmx->s.R);

  if (nO < 0)
  {
    nO = -nO;
    if (nO % 2 && iLoopInv == 0) nO *= 2;
  }
  else
  {
    if (nO % 2 && iLoopInv != 0) nO *= 2;
  }

  if (iLoopInv == 0)
    for (i = 0; i < 9; i++) rmxa[i] =  lsmx->s.R[i];
  else
    for (i = 0; i < 9; i++) rmxa[i] = -lsmx->s.R[i];

  InitRotMx(CumRMx, 1);

  if (nO > 1)
  {
    rmxab = rmxa;
    rmxb = NULL;

    for (iO = 1;;)
    {
      for (i = 0; i < 9; i++)
        CumRMx[i] += rmxab[i];

      if (++iO == nO)
        break;

      rmxb = rmxab;
      if (iO % 2) rmxab = BufMx1;
      else        rmxab = BufMx2;

      RotMxMultiply(rmxab, rmxa, rmxb);
    }
  }

  return nO;
}


int InverseRTMx(const T_RTMx *RTMx, T_RTMx *InvRTMx, int RBF)
{
  int  DetF, i;


  /*
     InvRTMx->s.R =  Inverse[RTMx->s.R]
     InvRTMx->s.T = -Inverse[RTMx->s.R] * RTMx->s.T
   */

  if (RBF == 0) RBF = 1;

      DetF = deterRotMx(RTMx->s.R);
  if (DetF %  (RBF * RBF) || DetF == 0) return 0;
      DetF /= (RBF * RBF);

  iCoFactorMxTp(RTMx->s.R, InvRTMx->s.R);

  for (i = 0; i < 9; i++) {
    if (InvRTMx->s.R[i] %  DetF) return 0;
        InvRTMx->s.R[i] /= DetF;
  }

  RotMx_t_Vector(InvRTMx->s.T, InvRTMx->s.R, RTMx->s.T, 0);

  for (i = 0; i < 3; i++) {
    if (InvRTMx->s.T[i] %  RBF) return 0;
        InvRTMx->s.T[i] /= RBF;
        InvRTMx->s.T[i] *= -1;
  }

  return DetF;
}


int SetPgInfo(const T_SgInfo *SgInfo, T_SgInfo *PgInfo)
{
  int                iList, i;
  T_RTMx             SMx[1];
  const T_RTMx       *lsmx;
  T_RotMxInfo        RotMxInfo;
  const T_RotMxInfo  *lrmxi;


  ResetSgInfo(PgInfo);

  PgInfo->GenOption = 1; /* faster group generation */

  SMx->s.T[0] =
  SMx->s.T[1] =
  SMx->s.T[2] = 0;

  lsmx = SgInfo->ListSeitzMx;

  for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
  {
        lrmxi = ListOrBufRotMxInfo(SgInfo, iList, &RotMxInfo);
    if (lrmxi == NULL)
      return -1;

    /* Since GenOption == 1 and Patterson groups are always
       centro-symmetric, convert improper to proper symmetry
       operations
     */

    if (lrmxi->Order > 0)
      for (i = 0; i < 9; i++)
        SMx->s.R[i] =  lsmx->s.R[i];
    else
      for (i = 0; i < 9; i++)
        SMx->s.R[i] = -lsmx->s.R[i];

    if (Add2ListSeitzMx(PgInfo, SMx) != 0)
      return -1;
  }

  if (AddInversion2ListSeitzMx(PgInfo))
    return -1;

  if (AddLatticeTr2ListSeitzMx(PgInfo, SgInfo->LatticeInfo))
    return -1;

  if (CompleteSgInfo(PgInfo) != 0)
    return -1;

  return 0;
}


void MathMx(const int *MxN, const int *MxD, int D,
            int nr, int nc, int Width, int NewLine)
{
  int  ir, ic;


  if (Width <= 0) Width = 4;

  if (nr > 1) putc('{', stdout);
  for (ir = 0; ir < nr; ir++)
  {
    putc('{', stdout);
    for (ic = 0; ic < nc; ic++)
    {
      if (ic) putc(',', stdout);

      if (MxD)
        Fprintf(stdout, "%*s", Width,
          FormatFraction(*MxN++, *MxD++, 0, NULL, 0));
      else if (D)
        Fprintf(stdout, "%*s", Width,
          FormatFraction(*MxN++,    D,   0, NULL, 0));
      else
        Fprintf(stdout, "%*d", Width,
          *MxN++);
    }
    putc('}', stdout);
    if (ir + 1 < nr) {
      putc(',', stdout);
      if (NewLine) {
        putc('\n', stdout);
        putc(' ', stdout);
      }
    }
  }
  if (nr > 1) putc('}', stdout);
}
