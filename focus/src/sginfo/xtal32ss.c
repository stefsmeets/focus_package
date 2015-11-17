/*
  Space Group Info (c) 1994-97 Ralf W. Grosse-Kunstleve
 */

#include <stdio.h>
#include <stdlib.h>


#include "sginfo.h"


#define Fprintf (void) fprintf


static int PrimitiveRotMx(const int *CCMx_LP, int *RotMx, const int *CCMx_PL,
                          int deterCCMx_LP)
{
  int       i;
  int       BufMx[9];


  /* Mp = Tlp . Mz . Tpl */

  RotMxMultiply(BufMx, RotMx, CCMx_PL);
  RotMxMultiply(RotMx, CCMx_LP, BufMx);

  for (i = 0; i < 9; i++)
  {
    if (RotMx[i] % deterCCMx_LP) {
      SetSgError("Internal Error: PrimitiveRotMx()");
      return -1;
    }
  }

  for (i = 0; i < 9; i++)
    RotMx[i] /= deterCCMx_LP;

  return 0;
}


int xtal32ss(T_SgInfo *SgInfo)
{
  static const int Tab_ss_Vector[] =
    {
       1,  0,  0,   0, /*  h      */
       0,  1,  0,   1, /*  k      */
       0,  0,  1,   2, /*  l      */
       1,  1,  0,   0, /*  h+k    */
       1, -1,  0,   0, /*  h-k    */
       0,  1,  1,   1, /*  k+l    */
       0,  1, -1,   1, /*  k-l    */
       1,  0,  1,   1, /*  h+l    */
       1,  0, -1,   1, /*  h-l    */
       1,  1,  1,   0, /*  h+k+l  */
       1,  1, -1,   0, /*  h+k-l  */
       1, -1,  1,   0, /*  h-k+l  */
      -1,  1,  1,   0, /* -h+k+l  */
       2,  1, -1,   0, /*  2h+k-l */
       2, -1,  1,   0, /*  2h-k+l */
      -1,  2,  1,   0, /* -h+2k+l */
       1,  2, -1,   0, /*  h+2k-l */
      -1,  1,  2,   0, /* -h+k+2l */
       1, -1,  2,   0  /*  h-k+2l */
    };

  static int nTab_ss_Vector
     = sizeof Tab_ss_Vector / sizeof (*Tab_ss_Vector) / 4;

  int        deterCCMx_LP, CCMx_PL[9];
  int        i, itabsiv;
  int        nLoopInv, iLoopInv, n_ssVM, i_ssVM;
  int        n, m, l;
  int        IsFine;
  int        item[3];
  int        RmI[9], ss_Buf[9];
  int        iList;
  T_RTMx     *lsmx;
  const int  *tabsiv;


  if (SgInfo->LatticeInfo->Code != 'P')
  {
    deterCCMx_LP = deterRotMx(SgInfo->CCMx_LP);
                iCoFactorMxTp(SgInfo->CCMx_LP, CCMx_PL);

    if (deterCCMx_LP < 1)
      goto ReturnError;
  }

  nLoopInv = Sg_nLoopInv(SgInfo);

  SgInfo->n_ssVM = n_ssVM = 0;

  for (i_ssVM = 0; i_ssVM < 3; i_ssVM++)
  {
    for (i = 0; i < 3; i++)
      SgInfo->ssVM[i_ssVM].V[i] = 0;

    SgInfo->ssVM[i_ssVM].M = 1;
    item[i_ssVM] = 1;
  }

  tabsiv = Tab_ss_Vector;

  for (itabsiv = 0; itabsiv < nTab_ss_Vector; itabsiv++, tabsiv += 4)
  {
    IsFine = 1;
    m = -1;

    for (iList = 0; IsFine && iList < SgInfo->nList; iList++)
    {
      lsmx = &SgInfo->ListSeitzMx[iList];

      for (iLoopInv = 0; IsFine && iLoopInv < nLoopInv; iLoopInv++)
      {
        if (iLoopInv == 0)
          for (i = 0; i < 9; i++)
          {
            if (i % 4) RmI[i] =  lsmx->s.R[i];
            else       RmI[i] =  lsmx->s.R[i] - 1;
          }
        else
          for (i = 0; i < 9; i++)
          {
            if (i % 4) RmI[i] = -lsmx->s.R[i];
            else       RmI[i] = -lsmx->s.R[i] - 1;
          }

        if (SgInfo->LatticeInfo->Code != 'P')
        {
          if (PrimitiveRotMx(SgInfo->CCMx_LP, RmI, CCMx_PL,
                                deterCCMx_LP) < 0)
            return -1;
        }

        for (i = 0; IsFine && i < 3; i++)
        {
          n =  tabsiv[0] * RmI[i * 3 + 0];
          n += tabsiv[1] * RmI[i * 3 + 1];
          n += tabsiv[2] * RmI[i * 3 + 2];
          n = abs(n);

          if (n == 1)
            IsFine = 0;
          else if (m < 2)
            m = n;
          else if (n > 0 && n != m)
            IsFine = 0;
        }
      }
    }

    if (IsFine)
    {
#if DEBUG_xtal32ss
      Fprintf(stdout, "H-Kt %2d %2d %2d   %d\n",
        tabsiv[0], tabsiv[1], tabsiv[2], m);
#endif

      l = tabsiv[3];

      while (item[l] > 1) /* just "if", see break's */
      {
        if (m == item[l]) break;

        if (m == 3 && (   SgInfo->XtalSystem != XS_Trigonal
                       || SgInfo->UniqueDirCode != '=')) break;

        if (m == 4 && (   SgInfo->XtalSystem == XS_Triclinic
                       || SgInfo->XtalSystem == XS_Monoclinic)) break;

        if (m == 2) break;

        /* if (m > 1 || m != 4) break; */

        n_ssVM--;
        item[l] = 1;
        break;
      }

      if (item[l] == 1)
      {
        if (itabsiv > 12)
          n_ssVM = 0;

        item[l] = m;
        SgInfo->ssVM[n_ssVM].M = m;

        for (i = 0; i < 3; i++)
          SgInfo->ssVM[n_ssVM].V[i] = tabsiv[i];

        n_ssVM++;
      }
    }
  }

#if DEBUG_xtal32ss
  Fprintf(stdout, "H-Kt\n");
#endif

  if (SgInfo->LatticeInfo->Code != 'P')
  {
#if DEBUG_xtal32ss
    for (i_ssVM = 0; i_ssVM < n_ssVM; i_ssVM++)
      Fprintf(stdout, "H-Kp %2d %2d %2d   %d\n",
        SgInfo->ssVM[i_ssVM].V[0],
        SgInfo->ssVM[i_ssVM].V[1],
        SgInfo->ssVM[i_ssVM].V[2],
        SgInfo->ssVM[i_ssVM].M);
    Fprintf(stdout, "H-Kp\n");
#endif

    for (i_ssVM = 0; i_ssVM < n_ssVM; i_ssVM++)
    {
      n = i_ssVM * 3;

      for (i = 0; i < 3; i++)
      {
        ss_Buf[n + i]
          =   SgInfo->ssVM[i_ssVM].V[0] * CCMx_PL[i * 3 + 0]
            + SgInfo->ssVM[i_ssVM].V[1] * CCMx_PL[i * 3 + 1]
            + SgInfo->ssVM[i_ssVM].V[2] * CCMx_PL[i * 3 + 2];
      }
    }

    for (i_ssVM = 0; i_ssVM < n_ssVM; i_ssVM++)
    {
      n = i_ssVM * 3;

      for (i = 0; i < 3; i++)
      {
        if (ss_Buf[n + i] % deterCCMx_LP)
        {
          Fprintf(stdout, " %3d %3d %3d\n",
            ss_Buf[n + 0], ss_Buf[n + 1], ss_Buf[n + 2]);
          goto ReturnError;
        }

        SgInfo->ssVM[i_ssVM].V[i] = ss_Buf[n + i] / deterCCMx_LP;
      }
    }
  }

  SgInfo->n_ssVM = n_ssVM;
  return n_ssVM;

  ReturnError:

  SetSgError("Internal Error: xtal32ss()");
  return -1;
}
