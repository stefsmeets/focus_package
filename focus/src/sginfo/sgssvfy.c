/*
  Space Group Info (c) 1994-97 Ralf W. Grosse-Kunstleve
 */

#include <stdio.h>
#include <stdlib.h>


#include "sginfo.h"


#define Fprintf (void) fprintf


#define nxs_malloc(ptr, n) (ptr) = malloc((n) * sizeof (*(ptr)))


static void MarkLegalOrigins(const T_SgInfo *SgInfo, int *TestField)
{
  int           O[3], V[3], mx, my, mz, i;
  int           IsFine, iList, iLoopInv, nLoopInv;
  int           BufMx[9];
  const T_RTMx  *lsmx;
  int           nTrV, iTrV;
  const int     *TrV;


  nLoopInv = Sg_nLoopInv(SgInfo);

  nTrV = SgInfo->LatticeInfo->nTrVector;

  for (O[0] = 0; O[0] < 12; O[0]++)
  for (O[1] = 0; O[1] < 12; O[1]++)
  for (O[2] = 0; O[2] < 12; O[2]++)
  {
    IsFine = 1;

    for (iList = 0; IsFine && iList < SgInfo->nList; iList++)
    {
      lsmx = &SgInfo->ListSeitzMx[iList];

      for (iLoopInv = 0; IsFine && iLoopInv < nLoopInv; iLoopInv++)
      {
        if (iLoopInv == 0)
          for (i = 0; i < 9; i++)
          {
            if (i % 4) BufMx[i] =  lsmx->s.R[i];
            else       BufMx[i] =  lsmx->s.R[i] - 1;
          }
        else
          for (i = 0; i < 9; i++)
          {
            if (i % 4) BufMx[i] = -lsmx->s.R[i];
            else       BufMx[i] = -lsmx->s.R[i] - 1;
          }

        RotMx_t_Vector(V, BufMx, O, 12);

        TrV = SgInfo->LatticeInfo->TrVector;

        for (iTrV = 0; iTrV < nTrV; iTrV++)
        {
          mx = (V[0] * (STBF / 12) + *TrV++) % STBF;
          my = (V[1] * (STBF / 12) + *TrV++) % STBF;
          mz = (V[2] * (STBF / 12) + *TrV++) % STBF;

          if (mx == 0 && my == 0 && mz == 0)
            break;
        }

        if (iTrV == nTrV) IsFine = 0;
      }
    }

    *TestField++ = IsFine;

#if DEBUG_MarkLegalOrigins
    if (IsFine != 0)
      Fprintf(stdout, " %2d %2d %2d\n", O[0], O[1], O[2]);
#endif
  }
}


#define IsArbitraryShift(iShift) \
  (    (iShift) == 1 || (iShift) ==  5 \
    || (iShift) == 7 || (iShift) == 11)


static int Verify_ss(int h, int k, int l, const int *TestField)
{
  int    O[3], TH;


  for (O[0] = 0; O[0] < 12; O[0]++)
  for (O[1] = 0; O[1] < 12; O[1]++)
  for (O[2] = 0; O[2] < 12; O[2]++)
  {
    if (*TestField++)
    {
          TH = h * O[0] + k * O[1] + l * O[2];
          TH %= 12;
      if (TH) return 0;

      if (IsArbitraryShift(O[0])) TH += h;
      if (IsArbitraryShift(O[1])) TH += k;
      if (IsArbitraryShift(O[2])) TH += l;
      if (TH) return 0;
    }
  }

  return 1;
}


int Check_ssVM(T_SgInfo *SgInfo)
{
  int       h, k, l, iList;
  int       Maxh, Maxk, Maxl;
  int       Minh, Mink, Minl;
  int       nTestField, *TestField;
  int       nProperty, *Property, *pp;
  int       RetVal, would_be, is;


  TestField = NULL;
  Property  = NULL;
  RetVal    = -1;
                        nTestField = 12 * 12 * 12;
  nxs_malloc(TestField, nTestField);
  if (TestField == NULL) {
    SetSgError("Not enough core");
    goto CleanAndReturn;
  }

  MarkLegalOrigins(SgInfo, TestField);

  Maxh = Maxk = Maxl =  7;
  Minh = Mink = Minl = -7;

  nProperty =   (Maxh - Minh + 1)
              * (Maxk - Mink + 1)
              * (Maxl - Minl + 1);

  nxs_malloc(Property, nProperty);
  if (Property == NULL) {
    SetSgError("Not enough core");
    goto CleanAndReturn;
  }

  pp = Property;
  for (h = Minh; h <= Maxh; h++)
  for (k = Mink; k <= Maxk; k++)
  for (l = Minl; l <= Maxl; l++)
  {
    iList = IsSysAbsent_hkl(SgInfo, h, k, l, NULL);
    if (SgError != NULL)
      goto CleanAndReturn;

    if (iList == 0)
      *pp++ = Verify_ss(h, k, l, TestField);
    else
      *pp++ = -1;
  }

  pp = Property;
  for (h = Minh; h <= Maxh; h++)
  for (k = Mink; k <= Maxk; k++)
  for (l = Minl; l <= Maxl; l++)
  {
    is = *pp++;

    if (is >= 0) {
                would_be = Is_ss(SgInfo, h, k, l);
      if (is != would_be) {
        RetVal = 0;
        goto CleanAndReturn;
      }
    }
  }

  RetVal = 1;

  CleanAndReturn:

  if (Property)  free(Property);
  if (TestField) free(TestField);

  return RetVal;
}
