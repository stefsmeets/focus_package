/*
  Space Group Info (c) 1994-97 Ralf W. Grosse-Kunstleve
 */

#include <stdio.h>
#include <stdlib.h>


#include "sginfo.h"


#define Fprintf (void) fprintf


static int iPlaneMeetsLine(const int a[3], const int p[3],
                           const int b[3], const int q[3],
                           int m[3], const int *PseudoG)
{
  /* Plane: p.(x-a)=0
     Line:  b+k*q=x   => p.(b+k*q-a)=0
                      => k*q.p+p.(b-a)=0
                      => k=-(p.(b-a)/(q.p))
   */

  int  qp, b_a[3], b_ap, i;


      qp = iScalProd(q, p, PseudoG);
  if (qp == 0) {
    for (i = 0; i < 3; i++) m[i] = 0;
    return 1;
  }

  for (i = 0; i < 3; i++) b_a[i] = b[i] - a[i];

  b_ap = iScalProd(b_a, p, PseudoG);

  for (i = 0; i < 3; i++) {
        m[i] = -b_ap * q[i];
    if (m[i] %  qp) return -1;
        m[i] /= qp;
        m[i] %= CTBF;
  }

  for (i = 0; i < 3; i++) m[i] += b[i];

  return 0;
}


static void RefineGrid(int Grid[3], const int m[3])
{
  int  i, gcd, lcm;


  for (i = 0; i < 3; i++)
  {
    gcd = iGCD(m[i], CTBF);
    lcm = iLCM(Grid[i], CTBF / gcd);
    Grid[i] = lcm;
  }
}


static void ShowSpecialPoint(int iSymTr, int i, int m[3])
{
  Fprintf(stdout, "(%d,%d)", iSymTr + 1, i);

  for (i = 0; i < 3; i++)
    Fprintf(stdout, " %6s", FormatFraction(m[i], CTBF, 0, NULL, 0));

  putc('\n', stdout);
}


static void ShowSI(const T_SMxI *SI)
{
  int  i;


  Fprintf(stdout, "%2d", SI->Order);
  Fprintf(stdout, " [ %2d %2d %2d ] ",
    SI->Basis[2][0], SI->Basis[2][1], SI->Basis[2][2]);
  if (SI->SenseOfRotation < 0)
    putc('c', stdout);
  else
    putc(' ', stdout);
  putc(' ', stdout);
  Fprintf(stdout, "wg =");
  for (i = 0; i < 3; i++)
    Fprintf(stdout, " %6s", FormatFraction(SI->wg[i], STBF, 0, NULL, 0));
  Fprintf(stdout, ", Tr =");
  for (i = 0; i < 3; i++)
    Fprintf(stdout, " %6s", FormatFraction(SI->Tr[i], CTBF, 0, NULL, 0));
  putc('\n', stdout);
}


int SgGrid(T_SgInfo *SgInfo, int Grid[3][3], int Show)
{
  int        iSymTr, i, j, stat, Active_wg;
  T_SMxI     SI[1];
  T_RTMx     SMx[1];
  int        NAxes[3][3], NFaces[3][3];
  int        a[3], m[3], LTr[3], Tr[3];
  const int  *PseudoG;


  if (ExpandListSeitzMx(SgInfo) != 0)
    return -1;

  PseudoG = PseudoMetricalMatrix(SgInfo->XtalSystem,
                                 SgInfo->UniqueDirCode, SgInfo->UniqueRefAxis,
                                 0);
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      NAxes[i][j] = (j == i ? 1 : 0);

  for (i = 0; i < 3; i++)
    iCrossProd(NFaces[i], NAxes[i], NAxes[(i + 1) % 3], PseudoG);

  for (i = 0; i < 3; i++) a[i] = 0;

  for (i = 0; i < 3; i++)
  for (j = 0; j < 3; j++)
    Grid[i][j] = 1;

  for (iSymTr = 0; iSymTr < SgInfo->OrderL; iSymTr++)
  {
    ResetSeitzMxInfo(SI, 0);

    for (i = 0; i < 9; i++)
      SMx->s.R[i] = SgInfo->ListSeitzMx[iSymTr].s.R[i];

#define Loop(i) for (LTr[i] = -1; LTr[i] <= 1; LTr[i]++)
    Loop(0) Loop(1) Loop(2)
#undef Loop
    {
      for (i = 0; i < 3; i++)
        SMx->s.T[i] = SgInfo->ListSeitzMx[iSymTr].s.T[i] + LTr[i] * STBF;

      if (SetSeitzMxInfo(SMx, SI) == 0)
      {
        if (SgError != NULL) return -1;
        SetSgError("Corrupt ListSeitzMx");
        return -1;
      }

      if (Show) ShowSI(SI);

      Active_wg = 0;

      for (i = 0; i < 3; i++)
        if (iModPositive(SI->wg[i], STBF) != 0) {
          Active_wg = 1;
          break;
        }

      for (i = 0; i < 3; i++) Tr[i] = SI->wg[i] * (CTBF / STBF);
      RefineGrid(Grid[2], Tr);

      if      (SI->Order == 1)
      {
        if (Show) ShowSpecialPoint(iSymTr, 0, Tr);
        RefineGrid(Grid[1], Tr);
      }
      else if (SI->Order == -2)
      {
        for (i = 0; i < 3; i++)
        {
          stat = iPlaneMeetsLine(SI->Tr, SI->Basis[2],
                                 a, NAxes[i], m, PseudoG);
          if (stat < 0) return -1;
          if (stat == 0) {
            if (Show) ShowSpecialPoint(iSymTr, -(i + 1), m);
            RefineGrid(Grid[1], m);
            if (Active_wg == 0)
              RefineGrid(Grid[0], m);
          }
        }
      }
      else
      {
        if (SI->Order != -1)
        {
          for (i = 0; i < 3; i++)
          {
            stat = iPlaneMeetsLine(a, NFaces[i],
                                   SI->Tr, SI->Basis[2], m, PseudoG);
            if (stat < 0) return -1;
            if (stat == 0) {
              if (Show) ShowSpecialPoint(iSymTr, (i + 1), m);
              RefineGrid(Grid[1], m);
              if (Active_wg == 0)
                RefineGrid(Grid[0], m);
            }
          }
        }

        if (SI->Order < 0) {
          if (Show) ShowSpecialPoint(iSymTr, 0, SI->Tr);
          if (Active_wg != 0) { SetSgError("Corrupt SMxI->wg"); return -1; }
          RefineGrid(Grid[0], SI->Tr);
          RefineGrid(Grid[1], SI->Tr);
        }
      }
    }
  }

  return 0;
}
