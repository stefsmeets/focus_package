/*
  Space Group Info (c) 1994-97 Ralf W. Grosse-Kunstleve
 */

#include <stdio.h>
#include <stdlib.h>


#include "sginfo.h"


#define Fprintf (void) fprintf


int Verify_sghkl(const T_SgInfo *SgInfo, int FriedelSym,
                 int Maxh, int Maxk, int Maxl)
{
  int  h, k, l;
  int  Minh, Mink, Minl;
  int  iList, iEq, HKLabsent, Restr, EpsCond, EpsLoop, nNotHidden, ErrCount;
  int  nLoopInv, iLoopInv;

  int                              nRefl, iRefl, jRefl;
  struct { int h, k, l, M, R, f; } *Refl, *Ri, *Rj;

  T_Eq_hkl  Eq_hkl[1];


  nRefl =   (2 * Maxh + 1)
          * (2 * Maxk + 1)
          * (2 * Maxl + 1);

      Refl = malloc(nRefl * sizeof (*Refl));
  if (Refl == NULL) {
    SetSgError("Internal Error: Not enough core");
    return -1;
  }

  if (SetListMin_hkl(SgInfo, FriedelSym,
                      Maxh,  Maxk,  Maxl,
                     &Minh, &Mink, &Minl) != 0)
    return -1;

  Ri = Refl;

  for (h = -Maxh; h <= Maxh; h++)
  for (k = -Maxk; k <= Maxk; k++)
  for (l = -Maxl; l <= Maxl; l++, Ri++)
  {
    Ri->h = h;
    Ri->k = k;
    Ri->l = l;

    Ri->M = Mult_hkl(SgInfo, FriedelSym, h, k, l);
    if (SgError != NULL)
      return -1;

    Ri->f = 0;

    iList = LoopIsSysAbsent_hkl(SgInfo, h, k, l, &Ri->R);
    if (SgError != NULL)
      return -1;

    if (iList != 0)
      Ri->f = -1;
    else
    {
      iList = IsHidden_hkl(SgInfo, FriedelSym,
                           Minh, Mink, Minl,
                           Maxh, Maxk, Maxl,
                              h,    k,    l);
      if (iList == 0)
        Ri->f = 1;
    }
  }

  ErrCount = 0;

  Ri = Refl;

  for (iRefl = 0; iRefl < nRefl; iRefl++, Ri++)
  {
    if (SgInfo->nReflCond >= 0)
    {
      HKLabsent = CondIsSysAbsent_hkl(SgInfo->ReflCond, SgInfo->nReflCond,
                                      Ri->h, Ri->k, Ri->l);

      if ((Ri->f == -1) != (HKLabsent != 0)) {
        Fprintf(stdout,
          "Error: %3d %3d %3d %s(%d) vs. %s(%d) mismatch\n",
          Ri->h, Ri->k, Ri->l,
          "LoopIsSysAbsent_hkl", Ri->f,
          "CondIsSysAbsent_hkl", HKLabsent);
        ErrCount++;
      }

      if (Ri->f != -1 && SgInfo->nRestCond >= 0)
      {
        Restr = Get_hklRestriction(SgInfo->RestCond, SgInfo->nRestCond,
                                   Ri->h, Ri->k, Ri->l);
        if (Ri->R != Restr) {
          Fprintf(stdout,
            "Error: %3d %3d %3d %s(%d) vs. %s(%d) mismatch\n",
            Ri->h, Ri->k, Ri->l,
            "LoopIsSysAbsent_hkl", Ri->R,
            "Get_hklRestriction", Restr);
          ErrCount++;
        }
      }

      if (Ri->f != -1)
      {
        EpsLoop = Epsilon_hkl(SgInfo, Ri->h, Ri->k, Ri->l);

        EpsCond = Get_hklEpsilon(SgInfo->ReflCond, SgInfo->nReflCond,
                                 SgInfo->SysEnhanced,
                                 Ri->h, Ri->k, Ri->l);

        if (EpsLoop != EpsCond) {
          Fprintf(stdout,
            "Error: %3d %3d %3d %s(%d) vs. %s(%d) mismatch\n",
            Ri->h, Ri->k, Ri->l,
            "Epsilon_hkl", EpsLoop,
            "Get_hklEpsilon", EpsCond);
          ErrCount++;
        }
      }
    }

    if (Ri->f == 1)
    {
      if (Ri->h < Minh || Ri->k < Mink || Ri->l < Minl) {
        Fprintf(stdout, "Error: %3d %3d %3d not hidden\n",
          Ri->h, Ri->k, Ri->l);
        ErrCount++;
      }
    }

    (void) BuildEq_hkl(SgInfo, FriedelSym, Eq_hkl, Ri->h, Ri->k, Ri->l);
    if (SgError != NULL)
      return -1;

    if (Ri->M != Eq_hkl->M) {
      Fprintf(stdout,
        "Error: %3d %3d %3d Mult_hkl(%d) vs. BuildEq_hkl(%d) mismatch\n",
        Ri->h, Ri->k, Ri->l, Ri->M, Eq_hkl->M);
      ErrCount++;
    }

    nLoopInv = (Eq_hkl->Centric == 0 ? 1 : 2);

    nNotHidden = 0;

    for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
    {
      for (iEq = 0; iEq < Eq_hkl->N; iEq++)
      {
        h = Eq_hkl->h[iEq]; if (abs(h) > Maxh) continue;
        k = Eq_hkl->k[iEq]; if (abs(k) > Maxk) continue;
        l = Eq_hkl->l[iEq]; if (abs(l) > Maxl) continue;

        if (iLoopInv) {
          h *= -1;
          k *= -1;
          l *= -1;
        }

        Rj = Refl;

        for (jRefl = 0; jRefl < nRefl; jRefl++, Rj++)
          if (   h == Rj->h
              && k == Rj->k
              && l == Rj->l)
            break;

        if (jRefl == nRefl) {
          Fprintf(stdout,
            "Error: %3d %3d %3d sym. equiv. %3d %3d %3d not found\n",
            Ri->h, Ri->k, Ri->l, h, k, l);
          ErrCount++;
        }
        else
        {
          if (AreSymEquivalent_hkl(SgInfo, FriedelSym,
              Ri->h, Ri->k, Ri->l,
              Rj->h, Rj->k, Rj->l) == 0)
          {
            if (SgError != NULL)
              return -1;

            Fprintf(stdout,
              "Error: %3d %3d %3d vs. %3d %3d %3d %s failed\n",
              Ri->h, Ri->k, Ri->l,
              Rj->h, Rj->k, Rj->l, "AreSymEquivalent_hkl()");
            ErrCount++;
          }

          if (Ri->M != Rj->M) {
            Fprintf(stdout,
              "Error: %3d %3d %3d vs. %3d %3d %3d multiplicity mismatch\n",
              Ri->h, Ri->k, Ri->l,
              Rj->h, Rj->k, Rj->l);
            ErrCount++;
          }

          if ((Ri->f == -1) != (Rj->f == -1)) {
            Fprintf(stdout,
              "Error: %3d %3d %3d vs. %3d %3d %3d sys. abs. mismatch\n",
              Ri->h, Ri->k, Ri->l,
              Rj->h, Rj->k, Rj->l);
            ErrCount++;
          }

          if (iRefl != jRefl && Ri->f == 1 && Rj->f == 1) {
            Fprintf(stdout,
              "Error: %3d %3d %3d vs. %3d %3d %3d sym. equiv. not hidden\n",
              Ri->h, Ri->k, Ri->l,
              Rj->h, Rj->k, Rj->l);
            ErrCount++;
          }

          if (Rj->f == 1)
            nNotHidden++;
        }
      }
    }

    if (Ri->f != -1 && nNotHidden == 0) {
      Fprintf(stdout, "Error: %3d %3d %3d all sym. equiv. hidden\n",
        Ri->h, Ri->k, Ri->l);
      ErrCount++;
    }
  }

  free(Refl);

  return ErrCount;
}
