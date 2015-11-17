/*
  Space Group Info (c) 1994-97 Ralf W. Grosse-Kunstleve
 */

#include <stdio.h>
#include <stdlib.h>


#include "sginfo.h"


#define HmulR(hm, km, lm, h, k, l, R) \
  (hm) = (R)[0] * (h) + (R)[3] * (k) + (R)[6] * (l); \
  (km) = (R)[1] * (h) + (R)[4] * (k) + (R)[7] * (l); \
  (lm) = (R)[2] * (h) + (R)[5] * (k) + (R)[8] * (l)


static const char *IErr_Inc_SymMx =
  "Internal Error: Inconsistent symmetry matrices";


int LoopIsSysAbsent_hkl(const T_SgInfo *SgInfo,
                        int h, int k, int l, int *TH_Restriction)
{
  int           iTrV, nTrV;
  const int     *TrV;
  int           iList, hm, km, lm;
  int           TH, THr, FlagMismatch;
  const T_RTMx  *lsmx;


  /* check list of symmetry operations
     take care of lattice type and "centric" flag */

  THr = -1;
  if (TH_Restriction != NULL) *TH_Restriction = THr;
  FlagMismatch = 0;

  nTrV = SgInfo->LatticeInfo->nTrVector;
  lsmx = SgInfo->ListSeitzMx;

  for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
  {
    HmulR(hm, km, lm, h, k, l, lsmx->s.R);

    TrV = SgInfo->LatticeInfo->TrVector;

    for (iTrV = 0; iTrV < nTrV; iTrV++)
    {
      TH =  (lsmx->s.T[0] + *TrV++) * h;
      TH += (lsmx->s.T[1] + *TrV++) * k;
      TH += (lsmx->s.T[2] + *TrV++) * l;
      TH %= STBF; if (TH < 0) TH += STBF;

      if      (-h == hm && -k == km && -l == lm)
      {
        if (TH != 0 && SgInfo->Centric == -1)
          return -(iList + 1 + iTrV * SgInfo->nList);

        if (THr < 0)
          THr = TH;
        else if (THr != TH)
          FlagMismatch = 1; /* must be systematic absent */
                            /* will check later ...      */
      }
      else if ( h == hm &&  k == km &&  l == lm)
      {
        if (TH != 0)
          return  (iList + 1 + iTrV * SgInfo->nList);
      }
      else
        break;
    }
  }

  if (THr >= 0 && FlagMismatch) /* ... consistency check */
    SetSgError(IErr_Inc_SymMx);

  if (TH_Restriction != NULL)
  {
    if (SgInfo->Centric == -1) *TH_Restriction = 0;
    else                       *TH_Restriction = THr;
  }

  return 0;
}


int IsSysAbsent_hkl(const T_SgInfo *SgInfo,
                    int h, int k, int l, int *TH_Restriction)
{
  if (    SgInfo->nReflCond < 0
      || (SgInfo->nRestCond < 0 && TH_Restriction != NULL))
    return LoopIsSysAbsent_hkl(SgInfo, h, k, l, TH_Restriction);

  if (TH_Restriction != NULL)
  {
     *TH_Restriction = Get_hklRestriction(SgInfo->RestCond, SgInfo->nRestCond,
                                          h, k, l);
    if (SgError != NULL)
      return -1;
  }

  return CondIsSysAbsent_hkl(SgInfo->ReflCond, SgInfo->nReflCond, h, k, l);
}


int BuildEq_hkl(const T_SgInfo *SgInfo, int FriedelSym,
                T_Eq_hkl *Eq_hkl, int h, int k, int l)
{
  int       iList, hm, km, lm, i;
  T_RTMx    *lsmx;
  T_Eq_hkl  BufEq_hkl;


  if (Eq_hkl == NULL)
    Eq_hkl = &BufEq_hkl;

  Eq_hkl->Centric =    SgInfo->Centric == -1
                    || SgInfo->InversionOffOrigin != 0
                    || FriedelSym != 0;
  Eq_hkl->M = 0;
  Eq_hkl->N = 0;

  lsmx = SgInfo->ListSeitzMx;

  for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
  {
    HmulR(hm, km, lm, h, k, l, lsmx->s.R);

    for (i = 0; i < Eq_hkl->N; i++)
    {
      if ( hm == Eq_hkl->h[i] &&  km == Eq_hkl->k[i] &&  lm == Eq_hkl->l[i])
        break;
      if (Eq_hkl->Centric == 0)
        continue;
      if (-hm == Eq_hkl->h[i] && -km == Eq_hkl->k[i] && -lm == Eq_hkl->l[i])
        break;
    }

    if (i == Eq_hkl->N)
    {
      if (Eq_hkl->N >= sizeof Eq_hkl->h / sizeof (*Eq_hkl->h)) {
        SetSgError(IErr_Inc_SymMx);
        return 0;
      }

      Eq_hkl->h[i] = hm;
      Eq_hkl->k[i] = km;
      Eq_hkl->l[i] = lm;

          Eq_hkl->TH[i] = (  lsmx->s.T[0] * h
                           + lsmx->s.T[1] * k
                           + lsmx->s.T[2] * l) % STBF;
      if (Eq_hkl->TH[i] < 0)
          Eq_hkl->TH[i] += STBF;

      Eq_hkl->M++;
      if (Eq_hkl->Centric != 0 && (hm || km || lm))
        Eq_hkl->M++;

      Eq_hkl->N++;
    }
  }

  if (SgInfo->nList % Eq_hkl->N) /* another error trap */ {
    SetSgError(IErr_Inc_SymMx);
    return 0;
  }

  return Eq_hkl->M;
}


int AreSymEquivalent_hkl(const T_SgInfo *SgInfo, int FriedelSym,
                         int h1, int k1, int l1,
                         int h2, int k2, int l2)
{
  int     Centric;
  int     iList, hm, km, lm;
  T_RTMx  *lsmx;


  Centric =    SgInfo->Centric == -1
            || SgInfo->InversionOffOrigin != 0
            || FriedelSym != 0;

  lsmx = SgInfo->ListSeitzMx;

  for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
  {
    HmulR(hm, km, lm, h1, k1, l1, lsmx->s.R);

    if ( h2 == hm &&  k2 == km &&  l2 == lm)
      return  (iList + 1);
    if (Centric == 0)
      continue;
    if (-h2 == hm && -k2 == km && -l2 == lm)
      return -(iList + 1);
  }

  return 0;
}


int SetListMin_hkl(const T_SgInfo *SgInfo, int FriedelSym,
                   int  Maxh, int  Maxk, int  Maxl,
                   int *Minh, int *Mink, int *Minl)
{
  int  Centric;
  int  i, n, u;

  const char xyz[] = "xyz";


  Centric =    SgInfo->Centric == -1
            || SgInfo->InversionOffOrigin != 0
            || FriedelSym != 0;

  *Minh = -Maxh;
  *Mink = -Maxk;
  *Minl = -Maxl;

  for (i = 0; i < 3; i++) {
    if (FindSeitzMx(SgInfo, 3, 0, xyz[i], '=') >= 0) {
      if (Centric) *Minl = 0;
      return 0;
    }
    else if (SgError != NULL)
      return -1;
  }

  n = 0;
  u = 0;

  for (i = 0; i < 3; i++)
    if (FindSeitzMx(SgInfo, 2, 0, xyz[i], '=') >= 0) {
      n++;
      u = i;
    }
    else if (SgError != NULL)
      return -1;

  if (n == 2) {
    SetSgError(IErr_Inc_SymMx);
    return -1;
  }

  if (n == 0)
  {
    if (FriedelSym != 0)
      *Minl = 0;

    return 0;
  }

  if (n == 3)
  {
    if (Centric)
    {
      *Minh = 0;
      *Mink = 0;
      *Minl = 0;

      return 0;
    }

    for (i = 0; i < 3; i++)
      if (FindSeitzMx(SgInfo, 2, 1, xyz[i], '=') >= 0)
        break;
      else if (SgError != NULL)
        return -1;

    if (i == 3) {
      SetSgError(IErr_Inc_SymMx);
      return -1;
    }

    if      (i == 0) {            *Mink = 0; *Minl = 0; }
    else if (i == 1) { *Minh = 0;            *Minl = 0; }
    else             { *Minh = 0; *Mink = 0;            }

    return 0;
  }

  if (Centric)
  {
    if (u == 0) *Minh = 0;
    else        *Mink = 0;

    *Minl = 0;

    return 0;
  }

  if (FindSeitzMx(SgInfo, 2, 1, xyz[u], '=') >= 0)
  {
    if (u == 2) *Mink = 0;
    else        *Minl = 0;

    return 0;
  }
  else if (SgError != NULL) {
    SetSgError(IErr_Inc_SymMx);
    return -1;
  }

  if      (u == 0) *Minh = 0;
  else if (u == 1) *Mink = 0;
  else             *Minl = 0;

  return 0;
}


int CmpEq_hkl(int h1, int k1, int l1,
              int h2, int k2, int l2)
{
  if (l1 >= 0 && l2 <  0) return -1;
  if (l1 <  0 && l2 >= 0) return  1;

  if (k1 >= 0 && k2 <  0) return -1;
  if (k1 <  0 && k2 >= 0) return  1;

  if (h1 >= 0 && h2 <  0) return -1;
  if (h1 <  0 && h2 >= 0) return  1;

  if (abs(l1) < abs(l2)) return -1;
  if (abs(l1) > abs(l2)) return  1;

  if (abs(k1) < abs(k2)) return -1;
  if (abs(k1) > abs(k2)) return  1;

  if (abs(h1) < abs(h2)) return -1;
  if (abs(h1) > abs(h2)) return  1;

  return 0;
}


int IsHidden_hkl(const T_SgInfo *SgInfo, int FriedelSym,
                 int Minh, int Mink, int Minl,
                 int Maxh, int Maxk, int Maxl,
                 int    h, int    k, int    l)
{
  int     iList, Mate, hm, km, lm;
  T_RTMx  *lsmx;


  Mate = 0;

  if (FriedelSym != 0 || SgInfo->Centric == -1)
    Mate = 1;

  lsmx = &SgInfo->ListSeitzMx[0];

  for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
  {
    HmulR(hm, km, lm, h, k, l, lsmx->s.R);

    if (   iList != 0
        && (Minh <=  hm &&  hm <= Maxh)
        && (Mink <=  km &&  km <= Maxk)
        && (Minl <=  lm &&  lm <= Maxl)
        && CmpEq_hkl(h, k, l,  hm,  km,  lm) > 0)
      return  (iList + 1);

    if (   Mate != 0
        && (Minh <= -hm && -hm <= Maxh)
        && (Mink <= -km && -km <= Maxk)
        && (Minl <= -lm && -lm <= Maxl)
        && CmpEq_hkl(h, k, l, -hm, -km, -lm) > 0)
      return -(iList + 1);
  }

  return 0;
}


int Epsilon_hkl(const T_SgInfo *SgInfo, int h, int k, int l)
{
  int     Epsilon;
  int     iList, hm, km, lm;
  T_RTMx  *lsmx;


  if (! (h || k || l))
    return 0;

  Epsilon = 0;

  lsmx = SgInfo->ListSeitzMx;

  for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
  {
    HmulR(hm, km, lm, h, k, l, lsmx->s.R);

    if (   (    hm ==  h && km ==  k && lm ==  l)
        || (    SgInfo->Centric == -1
            &&  hm == -h && km == -k && lm == -l))
      Epsilon++;
  }

  if (Epsilon == 0 || SgInfo->nList % Epsilon != 0) {
    SetSgError(IErr_Inc_SymMx);
    return -1;
  }

  Epsilon *= SgInfo->LatticeInfo->nTrVector;

  return Epsilon;
}


int Mult_hkl(const T_SgInfo *SgInfo, int FriedelSym, int h, int k, int l)
{
  int     Centric;
  int     M, R;
  int     iList, hm, km, lm;
  T_RTMx  *lsmx;


  Centric =    SgInfo->Centric == -1
            || SgInfo->InversionOffOrigin != 0
            || FriedelSym != 0;
  M = 0;
  R = 0;

  lsmx = SgInfo->ListSeitzMx;

  for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
  {
    HmulR(hm, km, lm, h, k, l, lsmx->s.R);

    if      (hm ==  h && km ==  k && lm ==  l)
      M++;
    else if (hm == -h && km == -k && lm == -l)
      R++;
  }

  if (M == 0 || SgInfo->nList % M != 0 || (R != 0 && R != M)) {
    SetSgError(IErr_Inc_SymMx);
    return -1;
  }

  M = SgInfo->nList / M;

  if (Centric != 0 && R == 0 && (h || k || l))
    M *= 2;

  return M;
}
