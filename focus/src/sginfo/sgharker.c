/*
  Space Group Info (c) 1994-97 Ralf W. Grosse-Kunstleve
 */

#include <stdio.h>
#include <stdlib.h>


#include "sginfo.h"


static int CoreSetNormVect(const int *Ox, int *N, const int *PseudoG)
{
  int  i, j;
  int  Prime[3][3], Pt[3][3], CP[3][3];
  int  Dim, NonNullCP[3], nNonNullCP;


  Prime[0][0] = 11; Prime[0][1] = 13; Prime[0][2] = 17;
  Prime[1][0] = 19; Prime[1][1] = 23; Prime[1][2] = 29;
  Prime[2][0] = 31; Prime[2][1] = 37; Prime[2][2] = 41;

  for (i = 0; i < 3; i++)
    RotMx_t_Vector(Pt[i], Ox, Prime[i], 0);

  nNonNullCP = 0;

  for (i = 0; i < 3; i++)
  {
    iCrossProd(CP[i], Pt[i], Pt[(i + 1) % 2], PseudoG);

    NonNullCP[i] = 0;

    for (j = 0; j < 3; j++)
      if (CP[i][j]) {
         NonNullCP[i] = 1;
        nNonNullCP++;
        break;
      }
  }

  if (nNonNullCP == 0)
  {
    Dim = 0;

    for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++) {
            N[j] = Pt[i][j];
        if (N[j])
          Dim = 1;
      }

      if (Dim)
        break;
    }

    return Dim;
  }

  if (nNonNullCP == 1) {
    N[0] = N[1] = N[2] = 0;
    SetSgError("Internal Error: CoreSetNormVect()");
    return -1;
  }

  for (i = 0; i < 3; i++)
    if (NonNullCP[i]) break;

  for (j = 0; j < 3; j++)
    N[j] = CP[i][j];

  if (nNonNullCP == 2)
    return 2;

  if (iScalProd(N, Pt[2], PseudoG) == 0)
    return 2;

  N[0] = N[1] = N[2] = 0;
  return 3;
}


static int SetNormVect(const int *Ox, int *N, const int *PseudoG)
{
  int  Dim, GCD, i;
  int  Cp, Cm;


      Dim = CoreSetNormVect(Ox, N, PseudoG);
  if (Dim < 0)
    return -1;

  GCD = FindGCD(N, 3);

  if (GCD)
    for (i = 0; i < 3; i++)
      N[i] /= GCD;

  Cp = Cm = 0;

  for (i = 0; i < 3; i++) {
    if ( N[i] < 0) Cp++;
    if (-N[i] < 0) Cm++;
  }

  if (Cp > Cm || (Cp == Cm && N[0] < 0))
    for (i = 0; i < 3; i++)
      N[i] = -N[i];

  return Dim;
}


static int FindHarkerInfo(T_HarkerInfo *HarkerInfo, int nHI, T_RTMx *Ox)
{
  int  iHI, i;


  for (iHI = 0; iHI < nHI; iHI++)
  {
    for (i = 0; i < 12; i++)
      if (Ox->a[i] != HarkerInfo[iHI].Ox.a[i])
        goto NextOx;

    return iHI;

    NextOx:;
  }

  return -1;
}


static int HarkerInfoSortFunction(const T_HarkerInfo *a, const T_HarkerInfo *b)
{
  int                 i;


  if (a->m < b->m) return -1;
  if (a->m > b->m) return  1;

  if (a->Dim > b->Dim) return -1;
  if (a->Dim < b->Dim) return  1;

  if (a->Dim == 1 || a->Dim == 2)
  {
        i = CmpEigenVectors(a->N, b->N);
    if (i != 0)
      return i;

    for (i = 0; i < 3; i++)
    {
      if (a->Ox.s.T[i] < b->Ox.s.T[i]) return -1;
      if (a->Ox.s.T[i] > b->Ox.s.T[i]) return  1;
    }
  }

  for (i = 0; i < 12; i++)
  {
    if (a->Ox.a[i] < b->Ox.a[i]) return -1;
    if (a->Ox.a[i] > b->Ox.a[i]) return  1;
  }

  return 0;
}


int SetHarkerInfo(T_SgInfo *SgInfo,
                  T_SgInfo *PgInfo,
                  T_HarkerInfo *HarkerInfo, int mHI,
                  int Grid[3])
{
  int        iSymTr, jSymTr, pSymTr;
  T_RTMx     SmS[1], PSmS[1];
  int        nHI, iHI, i;
  const int  *PseudoG = NULL;


  if (ExpandListSeitzMx(SgInfo) != 0)
    return -1;
  if (ExpandListSeitzMx(PgInfo) != 0)
    return -1;

  if (Grid) for (i = 0; i < 3; i++) Grid[i] = 1;

  if (HarkerInfo)
    PseudoG = PseudoMetricalMatrix(SgInfo->XtalSystem,
                                   SgInfo->UniqueDirCode,
                                   SgInfo->UniqueRefAxis,
                                   0);
  nHI = 0;

  for (iSymTr = 0; iSymTr < SgInfo->OrderP; iSymTr++)
  {
    for (jSymTr = 0; jSymTr < SgInfo->OrderL; jSymTr++)
    {
      for (i = 0; i < 12; i++)
        SmS->a[i] =   SgInfo->ListSeitzMx[iSymTr].a[i]
                    - SgInfo->ListSeitzMx[jSymTr].a[i];

      for (i = 0; i < 3; i++)
        SmS->s.T[i] = iModPositive(SmS->s.T[i], STBF);

      if (Grid)
        for (i = 0; i < 3; i++)
          Grid[i] = iLCM(Grid[i], STBF / iGCD(SmS->s.T[i], STBF));

      if (HarkerInfo)
      {
            iHI = FindHarkerInfo(HarkerInfo, nHI, SmS);
        if (iHI >= 0)
          HarkerInfo[iHI].m++;
        else
        {
          for (pSymTr = 1; pSymTr < PgInfo->OrderL; pSymTr++)
          {
            SeitzMxMultiply(PSmS, &PgInfo->ListSeitzMx[pSymTr], SmS);

                iHI = FindHarkerInfo(HarkerInfo, nHI, PSmS);
            if (iHI >= 0)
              break;
          }

          if (pSymTr == PgInfo->OrderL)
          {
            if (nHI == mHI) {
              SetSgError(
                "Internal Error: Allocated space for HarkerInfo too small");
              return -1;
            }

            for (i = 0; i < 9; i++)
              HarkerInfo[nHI].Ox.s.R[i] = SmS->s.R[i];

            for (i = 0; i < 3; i++)
              HarkerInfo[nHI].Ox.s.T[i] = iModPositive(SmS->s.T[i], STBF);

            HarkerInfo[nHI].m = 1;

                HarkerInfo[nHI].Dim = SetNormVect(HarkerInfo[nHI].Ox.s.R,
                                                  HarkerInfo[nHI].N,
                                                  PseudoG);
            if (HarkerInfo[nHI].Dim < 0)
              return -1;

            nHI++;
          }
        }
      }
    }
  }

  if (nHI > 1)
    qsort((void *) HarkerInfo, nHI, sizeof (*HarkerInfo),
          (int (*)(const void *, const void *)) HarkerInfoSortFunction);

  return nHI;
}
