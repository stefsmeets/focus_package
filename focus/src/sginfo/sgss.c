/*
  Space Group Info (c) 1994-97 Ralf W. Grosse-Kunstleve
 */

#include <stdio.h>
#include <stdlib.h>


#include "sginfo.h"


int Is_ss(const T_SgInfo *SgInfo, int h, int k, int l)
{
  int           i_ssVM, u;
  const T_ssVM  *ssVM;


  ssVM = SgInfo->ssVM;

  for (i_ssVM = 0; i_ssVM < SgInfo->n_ssVM; i_ssVM++, ssVM++)
  {
    u =  ssVM->V[0] * h;
    u += ssVM->V[1] * k;
    u += ssVM->V[2] * l;

    if (ssVM->M) {
      if (u % ssVM->M) return 0; }
    else {
      if (u)           return 0; }
  }

  return 1;
}


void Set_uvw(const T_SgInfo *SgInfo, int h, int k, int l, int *uvw)
{
  int           i_ssVM, u;
  const T_ssVM  *ssVM;


  ssVM = SgInfo->ssVM;

  for (i_ssVM = 0; i_ssVM < SgInfo->n_ssVM; i_ssVM++, ssVM++)
  {
    u =  ssVM->V[0] * h;
    u += ssVM->V[1] * k;
    u += ssVM->V[2] * l;

    if (ssVM->M) u %= ssVM->M;

    uvw[i_ssVM] = u;
  }
}


static int Test_ssVM(const T_SgInfo *SgInfo, int *Shift, int Mod)
{
  int           V[3], m[3], iv;
  int           iList, iLoopInv, nLoopInv;
  int           RmI[9];
  const T_RTMx  *lsmx;
  int           nTrV, iTrV;
  const int     *TrV;


  nLoopInv = Sg_nLoopInv(SgInfo);
  nTrV = SgInfo->LatticeInfo->nTrVector;

  for (iList = 0; iList < SgInfo->nList; iList++)
  {
    lsmx = &SgInfo->ListSeitzMx[iList];

    for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
    {
      SetRminusI(lsmx->s.R, RmI, iLoopInv);

      RotMx_t_Vector(V, RmI, Shift, 0);

      if (Mod == 0)
      {
        for (iv = 0; iv < 3; iv++)
          if (V[iv] != 0)
            return 0;
      }
      else
      {
        for (iv = 0; iv < 3; iv++)
        {
              V[iv] *= STBF;
          if (V[iv] % Mod != 0)
            return 0;

          V[iv] /= Mod;
        }

        TrV = SgInfo->LatticeInfo->TrVector;

        for (iTrV = 0; iTrV < nTrV; iTrV++, TrV += 3)
        {
          for (iv = 0; iv < 3; iv++)
          {
                m[iv] = (V[iv] + TrV[iv]) % STBF;
            if (m[iv] != 0)
              break;
          }

          if (iv == 3)
            break;
        }

        if (iTrV == nTrV)
          return 0;
      }
    }
  }

  return 1;
}


static const int Trial_ssVM[] =
  {
    3, 2,   1,  0,  0,   0, 2,
    3, 2,   1,  1,  0,   4, 2,
    3, 1,   1, -1,  0,   3,
    3, 1,   2,  4,  3,   6,
    3, 1,   2,  0,  3,   6,
    3, 1,   2, -1,  0,   4,
    3, 1,   1,  1, -1,   4,
    1, 3,   1,  1,  1,   0, 4, 2,
    0
  };


static int Collect_ssVM(const T_SgInfo *SgInfo, T_ssVM *Vs)
{
  const int  *TVM;
  int        iv, nVs, nPerm, iPerm, nMods, iMods;
  int        V[3], M[3];


  nVs = 0;

          TVM = Trial_ssVM;
  while (*TVM)
  {
    nPerm = *TVM++;
    nMods = *TVM++;

    for (iv = 0; iv < 3; iv++)
      V[iv] = *TVM++;

    for (iMods = 0; iMods < nMods; iMods++)
      M[iMods] = *TVM++;

    for (iPerm = 0; iPerm < nPerm; iPerm++)
    {
      if (iPerm)
      {
          iv = V[2];
        V[2] = V[1];
        V[1] = V[0];
        V[0] = iv;
      }

      for (iMods = 0; iMods < nMods; iMods++)
      {
        if (Test_ssVM(SgInfo, V, M[iMods]) == 1)
        {
          for (iv = 0; iv < 3; iv++)
            Vs[nVs].V[iv] = V[iv];

          Vs[nVs].M  = M[iMods];
             nVs++;

          break;
        }
      }
    }
  }

  return nVs;
}


static int SetNextPt(const int *V, int M, int f,
                     const int *PivotPt, int *NextPt)
{
  int  iv;


  for (iv = 0; iv < 3; iv++)
  {
        NextPt[iv] = V[iv] * f * 12;
    if (NextPt[iv] % M) {
      SetSgError("Internal Error: SetNextPt()");
      return -1;
    }

    NextPt[iv] = iModPositive(PivotPt[iv] + NextPt[iv] / M, 12);
  }

  return (NextPt[0] * 12 + NextPt[1]) * 12 + NextPt[2];
}


static int MarkCoveredOriginShifts(const T_ssVM *BS, int nBS, int *Field)
{
  int        iBS, iF, iv;
  int        PivotPt[7][3], f[6];


  if (nBS > sizeof f / sizeof (*f)) {
    SetSgError("Internal Error: MarkCoveredOriginShifts()");
    return -1;
  }

  if (nBS == 0) {
    Field[0] = 1;
    return 0;
  }

  for (iv = 0; iv < 3; iv++) PivotPt[0][iv] = 0;

  f[0] = 0;

  iBS = 0;

  for (;;)
  {
        iF = SetNextPt(BS[iBS].V, BS[iBS].M, f[iBS],
                       PivotPt[iBS], PivotPt[iBS + 1]);
    if (iF < 0)
      return -1;

    Field[iF] = 1;

    if (iBS + 1 < nBS)
    {
      iBS++;
      f[iBS] = 0;
    }
    else for (;;)
    {
      if (f[iBS] + 1 < BS[iBS].M) {
        f[iBS]++;
        break;
      }

      if (iBS == 0)
        goto Done;

      iBS--;
    }
  }

  Done: return 0;
}


static int IsCoveredVM(const int *Field, const int *V, int M)
{
  int  f, iF, jF;
  int  PivotPt[3], NextPt[3];


  if (M == 0)
      M = 12;

  iF = 0;

  for (PivotPt[0] = 0; PivotPt[0] < 12; PivotPt[0]++)
  for (PivotPt[1] = 0; PivotPt[1] < 12; PivotPt[1]++)
  for (PivotPt[2] = 0; PivotPt[2] < 12; PivotPt[2]++, iF++)
  {
    if (Field[iF] == 0)
      continue;

    for (f = 0; f < M; f++)
    {
          jF = SetNextPt(V, M, f, PivotPt, NextPt);
      if (jF < 0)
        return -1;

      if (Field[jF] == 0)
        return 0;
    }
  }

  return 1;
}


static int IsBasicSet(const T_SgInfo *SgInfo,
                      const T_ssVM *Vs, int nVs,
                      const int *Ix, int nIx)
{
  int        nBS, iIx, iVs, iv, iF, Stat;
  T_ssVM     BS[6];
  int        nTrV, iTrV, gcd;
  const int  *TrV;
  int        *Field;


      Field = malloc((12 * 12 * 12) * sizeof (*Field));
  if (Field == NULL) {
    SetSgError("Not enough core");
    return -1;
  }

  nBS = 0;

  for (iIx = 0; iIx < nIx; iIx++)
  {
    for (iv = 0; iv < 3; iv++)
      BS[nBS].V[iv] = Vs[Ix[iIx]].V[iv];

        BS[nBS].M = Vs[Ix[iIx]].M;
    if (BS[nBS].M == 0)
        BS[nBS].M = 12;

    nBS++;
  }

  nTrV = SgInfo->LatticeInfo->nTrVector;
   TrV = SgInfo->LatticeInfo->TrVector + 3;

  for (iTrV = 1; iTrV < nTrV; iTrV++, TrV += 3)
  {
    gcd = FindGCD(TrV, 3);

    for (iv = 0; iv < 3; iv++)
      BS[nBS].V[iv] = TrV[iv] / gcd;

    BS[nBS].M = STBF / gcd;
       nBS++;
  }

  for (iF = 0; iF < (12 * 12 * 12); iF++)
    Field[iF] = 0;

  if (MarkCoveredOriginShifts(BS, nBS, Field) != 0) {
    free(Field);
    return -1;
  }

  for (iVs = 0; iVs < nVs; iVs++)
  {
    for (iIx = 0; iIx < nIx; iIx++)
      if (iVs == Ix[iIx]) goto Next_iVs;

        Stat = IsCoveredVM(Field, Vs[iVs].V, Vs[iVs].M);
    if (Stat != 1) {
      free(Field);
      return Stat;
    }

    Next_iVs:;
  }

  free(Field);

  return 1;
}


static int NextOf_n_from_m(int m, int n, int *ix)
{
  int  p, l;


  p = l = n - 1;

  for (; p >= 0;)
  {
        ix[p]++;
    if (ix[p] == m - l + p)
      p--;
    else if (p < l)
    {
      ix[p + 1] = ix[p];
         p++;
    }
    else
      return 1;
  }

  return 0;
}


static int Select_ssVM(T_SgInfo *SgInfo, T_ssVM *Vs, int nVs)
{
  int  Ix[3], nIx, iIx;
  int  Stat, iv;


  for (nIx = 0; nIx <= nVs && nIx <= 3; nIx++)
  {
    for (iIx = 0; iIx < nIx; iIx++)
      Ix[iIx] = iIx;

    do
    {
          Stat = IsBasicSet(SgInfo, Vs, nVs, Ix, nIx);
      if (Stat < 0)
        return -1;
      if (Stat)
      {
        SgInfo->n_ssVM = nIx;

        for (iIx = 0; iIx < nIx; iIx++)
        {
          for (iv = 0; iv < 3; iv++)
            SgInfo->ssVM[iIx].V[iv] = Vs[Ix[iIx]].V[iv];

          SgInfo->ssVM[iIx].M = Vs[Ix[iIx]].M;
        }

        return 0;
      }
    }
    while (NextOf_n_from_m(nVs, nIx, Ix) != 0);
  }

  SetSgError("Internal Error: Select_ssVM()");
  return -1;
}


int Set_ss(T_SgInfo *SgInfo)
{
  int     nVs;
  T_ssVM  Vs[22];


  SgInfo->n_ssVM = -1;

      nVs = Collect_ssVM(SgInfo, Vs);
  if (nVs > 22)
  {
    SetSgError("Internal Error: nVs > 22");
    return -1;
  }

  if (nVs == 22) /* shortcut for space group 1 */
    nVs = 3;

  if (Select_ssVM(SgInfo, Vs, nVs) != 0)
    return -1;

  return SgInfo->n_ssVM;
}
