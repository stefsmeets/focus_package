#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "main.h"
#include "xtal.h"
#include "trialoop.h"
#include "matrix.h"
#include "io.h"


static int H2iH(int H, int NH, int FullRange)
{
  int  m;


  m = (NH - 1) / 2;

  if (FullRange)
  {
    if (-m > H || H > m) return -1;

    if (H < 0)
        H += NH;
  }
  else
    if ( 0 > H || H > m) return -1;

  return H;
}


int iH2H(int iH, int NH, int FullRange)
{
  if (FullRange)
  {
    if (iH < NH / 2) return iH;
    if (iH > NH / 2) return iH - NH;
    if (NH % 2 != 0) return iH;
    InternalError(NULL);
  }

  return iH;
}


static int hkl2ihkl(int   h, int   k, int   l,
                    int *ih, int *ik, int *il)
{
  *ih = H2iH(h, Nx, 1);
  *ik = H2iH(k, Ny, 1);
  *il = H2iH(l, Nz, 0);

  if (*ih < 0 || *ik < 0 || *il < 0) return -1;

  return 0;
}


static void ihkl2hkl(int ih, int ik, int il,
                     int *h, int *k, int *l)
{
  *h = iH2H(ih, Nx, 1);
  *k = iH2H(ik, Ny, 1);
  *l = iH2H(il, Nz, 0);
}


int IndexEq_hkl(T_Eq_hkl *Eq_hkl, int h, int k, int l, int *Conjugate)
{
  int  i;


  for (i = 0; i < Eq_hkl->N; i++)
  {
    if ( h == Eq_hkl->h[i] &&  k == Eq_hkl->k[i] &&  l == Eq_hkl->l[i]) {
      if (Conjugate) *Conjugate =  1;
      return i;
    }

    if (-h == Eq_hkl->h[i] && -k == Eq_hkl->k[i] && -l == Eq_hkl->l[i]) {
      if (Conjugate) *Conjugate = -1;
      return i;
    }
  }

  return -1;
}


int hkl_Compare(int a_h, int a_k, int a_l, int b_h, int b_k, int b_l)
{
  int       val_a, val_b;

  val_a = a_h * a_h + a_k * a_k + a_l * a_l;
  val_b = b_h * b_h + b_k * b_k + b_l * b_l;
  if (val_a < val_b) return -1;
  if (val_a > val_b) return  1;

  val_a = abs(a_h) + abs(a_k) + abs(a_l);
  val_b = abs(b_h) + abs(b_k) + abs(b_l);
  if (val_a < val_b) return -1;
  if (val_a > val_b) return  1;

  if (abs(a_h) > abs(a_k)) val_a = abs(a_h); else val_a = abs(a_k);
  if (abs(a_l) > val_a)    val_a = abs(a_l);
  if (abs(b_h) > abs(b_k)) val_b = abs(b_h); else val_b = abs(b_k);
  if (abs(b_l) > val_b)    val_b = abs(b_l);
  if (val_a < val_b) return -1;
  if (val_a > val_b) return  1;

  if (abs(a_h) < abs(b_h)) return -1;
  if (abs(a_h) > abs(b_h)) return  1;
  if (abs(a_k) < abs(b_k)) return -1;
  if (abs(a_k) > abs(b_k)) return  1;
  if (abs(a_l) < abs(b_l)) return -1;
  if (abs(a_l) > abs(b_l)) return  1;

  if (a_h < b_h) return -1;
  if (a_h > b_h) return  1;
  if (a_k < b_k) return -1;
  if (a_k > b_k) return  1;
  if (a_l < b_l) return -1;
  if (a_l > b_l) return  1;

  return 0;
}


static int Q_FobsRawSortFunction(const T_FobsRaw *a, const T_FobsRaw *b)
{
#define EPS_Q (1.e-6) /* ARBITRARY */

  if (a->Q + EPS_Q < b->Q) return -1;
  if (a->Q - EPS_Q > b->Q) return  1;
  return hkl_Compare(a->h, a->k, a->l, b->h, b->k, b->l);

#undef EPS_Q
}


static int G_FobsRawSortFunction(const T_FobsRaw *a, const T_FobsRaw *b)
{
  if (a->Status < b->Status) return -1;
  if (a->Status > b->Status) return  1;
  return hkl_Compare(a->h, a->k, a->l, b->h, b->k, b->l);
}


void SortFobsRaw(T_FobsRaw *FR, int nFR)
{
  int        iFR, jFR, GroupCount;
  T_FobsRaw  *iFRe, *jFRe;
  T_Eq_hkl   Eq_hkl;


  for (iFR = 0, iFRe = FR; iFR < nFR; iFR++, iFRe++)
    iFRe->Q = Q_hkl(iFRe->h, iFRe->k, iFRe->l, &LatConR);

  if (nFR > 1)
    qsort((void *) FR, nFR, sizeof (*FR),
          (SortFunction) Q_FobsRawSortFunction);

  GroupCount = 0;

  for (iFR = 0, iFRe = FR; iFR < nFR; iFR++, iFRe++)
  {
    if (iFRe->Status < 0)
    {
      iFRe->M = BuildEq_hkl(&SpgrInfo, 1, &Eq_hkl, iFRe->h, iFRe->k, iFRe->l);
      if (SgError != NULL) progerror(SgError);

      iFRe->Status = GroupCount;

      for (jFR = iFR + 1, jFRe = iFRe + 1; jFR < nFR; jFR++, jFRe++)
      {
        if (   jFRe->Status < 0
            && IndexEq_hkl(&Eq_hkl, jFRe->h, jFRe->k, jFRe->l, NULL) >= 0)
        {
          jFRe->M = Eq_hkl.M;
          jFRe->Status = GroupCount;
        }
      }

      GroupCount++;
    }
  }

  if (nFR > 1)
    qsort((void *) FR, nFR, sizeof (*FR),
          (SortFunction) G_FobsRawSortFunction);

  if (Debug0)
  {
    Fprintf(stdout, "After G_FobsRawSortFunction\n");
    for (iFR = 0, iFRe = FR; iFR < nFR; iFR++, iFRe++)
      Fprintf(stdout, " %3d %3d %3d %7.4f %2d %4d\n",
        iFRe->h, iFRe->k, iFRe->l, iFRe->Q, iFRe->M, iFRe->Status);
  }
}


void SetStatus1FobsRaw(T_FobsRaw *FR, int nFR)
{
  int       absent;
  int       PrevGroupCount;


  PrevGroupCount = -1;

  while (nFR--)
  {
    absent = IsSysAbsent_hkl(&SpgrInfo, FR->h, FR->k, FR->l,
                                      &(FR->PhaseRestriction));
    if (SgError != NULL) progerror(SgError);

    if (absent == 0)
      Set_uvw(&SpgrInfo, FR->h, FR->k, FR->l, FR->uvw);

    if (FR->Status == PrevGroupCount)
      FR->Status = FRS_SymEquiv;

    else
    {
      PrevGroupCount = FR->Status;

      if (! (FR->h || FR->k || FR->l))
        FR->Status = FRS_F000;
      else if (absent)
        FR->Status = FRS_SysAbsent;
      else
        FR->Status = FRS_Undefined;
    }

    FR++;
  }
}


void SetStatus2FobsRaw(T_FobsRaw *FR, int nFR)
{
  int       ih, ik, il, i;
  T_Eq_hkl  Eq_hkl;


  nCodeTransl = 0;

  while (nFR--)
  {
    if (FR->Status == FRS_Undefined)
    {
      if (    FR->Fmrg <= AppFabs(FR->SigmaFmrg * SigmaCutOff)
          || (FR->Q > FobsMaxQ && FobsMaxQ > CF0))
        FR->Status = FRS_Omit;

      else
      {
        (void) BuildEq_hkl(&SpgrInfo, 1, &Eq_hkl, FR->h, FR->k, FR->l);
        if (SgError != NULL) progerror(SgError);

        for (i = 0; i < Eq_hkl.N; i++)
        {
          if (   hkl2ihkl( Eq_hkl.h[i],  Eq_hkl.k[i],  Eq_hkl.l[i],
                          &ih, &ik, &il) != 0
              && hkl2ihkl(-Eq_hkl.h[i], -Eq_hkl.k[i], -Eq_hkl.l[i],
                          &ih, &ik, &il) != 0)
          {
            FR->Status = FRS_OffGrid;
            break;
          }
        }

        if (FR->Status == FRS_Undefined)
        {
          if (OverlapAction == OvlAct_Omit && FR->Overlap != 0)
            FR->Status = FRS_Omit;
          else {
            FR->Status = FRS_Active;
            nCodeTransl++;
          }
        }
      }
    }

    FR++;
  }
}


void CalcFmrg(T_FobsRaw *FR, int nFR)
{
  int        nEquiv;
  Fprec      SumF;
  Fprec      SumSF2;
  T_FobsRaw  *FRasym;


  FRasym = NULL;
  SumF   = 0.;
  SumSF2 = 0.;
  nEquiv = 0;

  for (;;)
  {
    if (nFR == 0 || FR->Status != FRS_SymEquiv)
    {
      if (FRasym)
      {
            FRasym->Fmrg      = SumF / nEquiv;
            FRasym->SigmaFmrg = AppSqrt(SumSF2) / nEquiv;

        if (FRasym->Status != FRS_F000)
        {
            FRasym->Fmrg      *= FobsScale;
            FRasym->SigmaFmrg *= FobsScale;
        }
      }

      if (nFR == 0)
        return;

      FRasym = FR;
      SumF   = FR->Fobs;
      SumSF2 = Square(FR->SigmaFobs);
      nEquiv = 1;
    }
    else
    {
      SumF   += FR->Fobs;
      SumSF2 += Square(FR->SigmaFobs);
      nEquiv++;
    }

    nFR--;
     FR++;
  }
}


Fprec TwoThetaDeg(Fprec Q)
{
  Fprec  TT;


  if (Q > 0.)
  {
    TT = LambdaLength * .5 * AppSqrt(Q);
    if (-1. <= TT  && TT <= 1.)
      TT = asin(TT) / PIover180 * 2.;
    else
      TT = 0.;
  }
  else
    TT = 0.;

  return TT;
}


void RedivideAbsent(T_FobsRaw *FR, int nFR)
{
  int        iFR, jFR;
  Fprec      TT0, TTbefore, TTnext;
  Fprec      SumI, I, Ib, In, Io, S, Sb, Sn, So, num, den;
  T_FobsRaw  *FRbefore, *FRnext, *FRj;


  AbsentRedivisionLimit.Frac_I_ignored = 0.;
  AbsentRedivisionLimit.Frac_I_w_moved = 0.;

  SumI = 0.;

  FRbefore = NULL;

  for (iFR = 0; iFR < nFR; iFR++, FR++)
  {
    I = Square(FR->M) * Square(FR->Fmrg);
    S = 2. * FR->M * FR->Fmrg * FR->SigmaFmrg;

    if (   FR->Status != FRS_F000
        && FR->Status != FRS_SymEquiv)
      SumI += I;

    if (FR->Status == FRS_SysAbsent)
    {
      FRnext = NULL;

      FRj = FR + 1;

      for (jFR = iFR + 1; jFR < nFR; jFR++, FRj++)
      {
        if (   FRj->Status != FRS_F000
            && FRj->Status != FRS_SysAbsent
            && FRj->Status != FRS_SymEquiv)
        {
          FRnext = FRj;
          break;
        }
      }

      TT0 = TwoThetaDeg(FR->Q);

      if (FRbefore && AbsentRedivisionLimit.Value > 0.)
      {
        TTbefore = TwoThetaDeg(FRbefore->Q);

        if (AbsentRedivisionLimit.Type == ARLT_Degree2Theta)
        {
          if (TT0 - TTbefore > AbsentRedivisionLimit.Value)
            TTbefore = -1.;
        }
        else
        {
          if (TT0 - TTbefore > FR->FWHM * AbsentRedivisionLimit.Value
                            || FR->FWHM == 0.)
            TTbefore = -1.;
        }
      }
      else
        TTbefore = -1.;

      if (FRnext && AbsentRedivisionLimit.Value > 0.)
      {
        TTnext = TwoThetaDeg(FRnext->Q);

        if (AbsentRedivisionLimit.Type == ARLT_Degree2Theta)
        {
          if (TTnext - TT0 > AbsentRedivisionLimit.Value)
            TTnext = -1.;
        }
        else
        {
          if (TTnext - TT0 > FR->FWHM * AbsentRedivisionLimit.Value
                          || FR->FWHM == 0.)
            TTnext = -1.;
        }
      }
      else
        TTnext = -1.;

      if (TTbefore >= 0. && TTnext >= 0.)
      {
        num = TTnext - TT0;
        den = TTnext - TTbefore;

        if (num >= den)
        {
          Ib = I * .5;
          Sb = S * .5;
          In = Ib;
          Sn = Sb;
        }
        else
        {
          Ib = I * num / den;
          Sb = S * num / den;
          In = I - Ib;
          Sn = S - Sb;
        }

        FR->Fmrg = 0.;
      }
      else if (TTbefore >= 0.)
      {
        Ib = I;
        Sb = S;
        In = 0.;
        Sn = 0.;

        FR->Fmrg = 0.;
      }
      else if (TTnext   >= 0.)
      {
        Ib = 0.;
        Sb = 0.;
        In = I;
        Sn = S;

        FR->Fmrg = 0.;
      }
      else
      {
        Ib = 0.;
        Sb = 0.;
        In = 0.;
        Sn = 0.;

        AbsentRedivisionLimit.Frac_I_ignored += I;
      }

      if (Ib != 0.)
      {
        AbsentRedivisionLimit.Frac_I_w_moved += Ib * (TT0 - TTbefore);

        Io = Square(FRbefore->M) * Square(FRbefore->Fmrg);
        So = 2. * FRbefore->M * FRbefore->Fmrg * FRbefore->SigmaFmrg;

        FRbefore->Fmrg = AppSqrt(Io + Ib) / FRbefore->M;
            den =   2. * AppSqrt(Io + Ib) * FRbefore->M;
        if (den != 0.)
          FRbefore->SigmaFmrg = AppSqrt(Square(So) + Square(Sb)) / den;
      }

      if (In != 0.)
      {
        AbsentRedivisionLimit.Frac_I_w_moved += In * (TTnext - TT0);

        Io = Square(FRnext->M) * Square(FRnext->Fmrg);
        So = 2. * FRnext->M * FRnext->Fmrg * FRnext->SigmaFmrg;

        FRnext->Fmrg = AppSqrt(Io + In) / FRnext->M;
            den = 2. * AppSqrt(Io + In) * FRnext->M;
        if (den != 0.)
          FRnext->SigmaFmrg = AppSqrt(Square(So) + Square(Sn)) / den;
      }
    }
    else if (   FR->Status != FRS_F000
             && FR->Status != FRS_SymEquiv)
      FRbefore = FR;
  }

  if (SumI)
  {
    AbsentRedivisionLimit.Frac_I_ignored /= SumI;
    AbsentRedivisionLimit.Frac_I_w_moved /= SumI;
  }
}


int SetOverlapGroups(T_FobsRaw *FR, int nFR)
{
  int        Group, Ignore;
  Fprec      TT1, TT2, DeltaTT, Critical;
  T_FobsRaw  *PrevValidFR;


  if (OverlapFactor == 0.)
    return 0;

  Ignore = 0;
  Group = 0;
  PrevValidFR = NULL;

  while (nFR--)
  {
    if (Ignore && FR->Status != FRS_SymEquiv)
        Ignore = 0;

    if (Ignore == 0)
    {
      if (   FR->Status == FRS_F000
          || FR->Status == FRS_SysAbsent)
        Ignore = 1;

      else
      {
        if (FR->Status == FRS_SymEquiv)
        {
          if (PrevValidFR)
            FR->Overlap = PrevValidFR->Overlap;
        }
        else
        {
          if (PrevValidFR != NULL)
          {
                      TT2       = TwoThetaDeg(FR->Q);
                            TT1 = TwoThetaDeg(PrevValidFR->Q);
            DeltaTT = TT2 - TT1;

            Critical = (PrevValidFR->FWHM + FR->FWHM) * OverlapFactor;

            if (DeltaTT <= Critical)
            {
              if (PrevValidFR->Overlap == 0)
              {
                Group++;

                do
                {
                  PrevValidFR->Overlap = Group;
                  PrevValidFR++;
                }
                while (PrevValidFR->Status == FRS_SymEquiv);
              }

              FR->Overlap = Group;
            }
          }

          PrevValidFR = FR;
        }
      }
    }

    FR++;
  }

  return Group;
}


void TreatOverlap(T_FobsRaw *FR, int nFR)
{
  int        nGroup;
  int        SumM;
  Fprec      SumMF2, SumM2F2sF2;
  Fprec      X, Y;
  T_FobsRaw  *FRfirst;


  if (   OverlapAction != OvlAct_EqualF2
      && OverlapAction != OvlAct_EqualMF2)
    return;

  FRfirst    = NULL;
  SumM       = 0;
  SumMF2     = 0.;
  SumM2F2sF2 = 0.;
  nGroup     = 0;

  for (;;)
  {
    if (   nFR == 0
        || FR->Overlap && (FRfirst == NULL || FR->Overlap != FRfirst->Overlap))
    {
      if (FRfirst)
      {
        if (OverlapAction == OvlAct_EqualF2)
        {
          if (SumM == 0)
            InternalError("Corrupt SumM");

          X = AppSqrt(SumMF2 / SumM);

          if (SumMF2 != 0.)
            Y = AppSqrt(SumM2F2sF2 / SumMF2 / SumM) / nGroup;
          else
            Y = 0.;

          do
          {
            if (FRfirst->Overlap && FRfirst->Status != FRS_SymEquiv)
            {
                FRfirst->Fmrg = X;
                FRfirst->SigmaFmrg = Y;
            }

            FRfirst++;
          }
          while (FRfirst < FR);
        }
        else
        {
          X = SumMF2 / nGroup;

          if (SumMF2 != 0.)
            Y = SumM2F2sF2 / SumMF2 / nGroup;
          else
            Y = 0.;

          do
          {
            if (FRfirst->Overlap && FRfirst->Status != FRS_SymEquiv)
            {
                FRfirst->Fmrg = AppSqrt(X / FRfirst->M);
                FRfirst->SigmaFmrg = AppSqrt(Y / FRfirst->M);
            }

            FRfirst++;
          }
          while (FRfirst < FR);
        }
      }

      if (nFR == 0)
        return;

      X = FR->M * Square(FR->Fmrg);

      FRfirst     = FR;
      SumM        = FR->M;
      SumMF2      = X;
      SumM2F2sF2  = FR->M * X * Square(FR->SigmaFmrg);
      nGroup      = 1;
    }
    else if (FR->Overlap && FR->Status != FRS_SymEquiv)
    {
      X = FR->M * Square(FR->Fmrg);

      SumM       += FR->M;
      SumMF2     += X;
      SumM2F2sF2 += FR->M * X * Square(FR->SigmaFmrg);
      nGroup++;
    }

    nFR--;
     FR++;
  }
}


void SetF000(void)
{
  int        nFR, i;
  T_FobsRaw  *FR;


  F000cal = 0.;

  for (i = 0; i < nAtomType; i++)
  {
    F000cal +=   AtomType[i].SF_Info.f_stol_0
               * AtomType[i].nPerUnitCell
               * AtomType[i].OccDefault;
  }

  F000mtp = F000cal;
  F000mrg = 0.;
  nFR = nFobsRaw;
   FR =  FobsRaw;

  while (nFR--)
  {
    if (FR->Status == FRS_F000)
    {
      F000mrg = FR->Fmrg;
      F000mtp = FR->Fmrg;
      break;
    }

    FR++;
  }
}


static Fprec CalcSumFmrg(T_CodeTransl *CT, int nCT)
{
  Fprec       Sum = 0.;

  while (nCT--)
  {
    Sum += CT->M_Fmrg;
    CT++;
  }

  return Sum;
}


#define EPS_Fmrg (1.e-6) /* ARBITRARY */

static int CT_U_SortFunction(const T_CodeTransl *a, const T_CodeTransl *b)
{
  if (  (a->FobsRaw->Status == FRS_Undefined)
      > (b->FobsRaw->Status == FRS_Undefined)) return -1;
  if (  (a->FobsRaw->Status == FRS_Undefined)
      < (b->FobsRaw->Status == FRS_Undefined)) return  1;

  if (a->M_Fmrg - EPS_Fmrg > b->M_Fmrg) return -1;
  if (a->M_Fmrg + EPS_Fmrg < b->M_Fmrg) return  1;

  return hkl_Compare(a->FobsRaw->h, a->FobsRaw->k, a->FobsRaw->l,
                     b->FobsRaw->h, b->FobsRaw->k, b->FobsRaw->l);
}


static int CT_F_SortFunction(const T_CodeTransl *a, const T_CodeTransl *b)
{
  if (a->M_Fmrg - EPS_Fmrg > b->M_Fmrg) return -1;
  if (a->M_Fmrg + EPS_Fmrg < b->M_Fmrg) return  1;

  return hkl_Compare(a->FobsRaw->h, a->FobsRaw->k, a->FobsRaw->l,
                     b->FobsRaw->h, b->FobsRaw->k, b->FobsRaw->l);
}

#undef EPS_Fmrg


Fprec InitCodeTransl(T_CodeTransl *CT, int nCT, T_FobsRaw *FR, int nFR)
{
  int           iCT;
  T_CodeTransl  *CT0;
  Fprec         MaxF;


  CT0 = CT;
  iCT = 0;
  MaxF = 0.;

  while (nFR--)
  {
    if (FR->Status == FRS_Undefined || FR->Status == FRS_Active)
    {
      if (iCT >= nCodeTransl)
        InternalError("Corrupt nCodeTransl");

      CT->FobsRaw = FR;
      if (MaxF < CT->FobsRaw->Fmrg) MaxF = CT->FobsRaw->Fmrg;
      CT->M_Fmrg = FR->M * CT->FobsRaw->Fmrg;
      CT->FF = NULL;
      CT->Fcal.r = CT->Fcal.i = 0.;
      CT->FcalAbs = CT->FcalPhi = 0.;
      CT->Phi360Fmrg = CT->Phi360Fcal = 0;
      CT->nCTData = 0;
      CT->CTData = NULL;

      iCT++;
      CT++;
    }

    FR++;
  }

  if (iCT != nCT)
    InternalError("Corrupt nCodeTransl");

  if (nCT > 1)
    qsort((void *) CT0, nCT, sizeof (*CT),
          (SortFunction) CT_U_SortFunction);

  return MaxF;
}


void SetSleeping(T_CodeTransl *CT, int nCT)
{
  int    iCT;
  Fprec  SumDyn, Sum;


  if (nCT > 1)
    qsort((void *) CT, nCT, sizeof (*CT), (SortFunction) CT_F_SortFunction);

  if (ReflectionUsage.Percent == 0)
  {
    if (ReflectionUsage.Value < 0.)
      return;

    iCT  = (int)(ReflectionUsage.Value + .5);
     CT += iCT;
  }
  else
  {
    SumDyn = CalcSumFmrg(CT, nCT);
    Sum = 0.;

    for (iCT = 0; iCT < nCT; iCT++, CT++)
    {
          Sum += CT->M_Fmrg;
      if (Sum / SumDyn > ReflectionUsage.Value * .01)
        break;
    }
  }

  for ( ; iCT < nCT; iCT++, CT++)
    CT->FobsRaw->Status = FRS_Sleeping;
}


static Fprec TabCos360[360];
static Fprec TabSin360[360];


static void CompleteCodeTransl(T_CodeTransl *CT, int nCT)
{
  int       iCT, h, k, l, i;
  T_CTData  *CTD;
  T_Eq_hkl  Eq_hkl;
  Fprec     phi;


  for (i = 0; i < 360; i++)
  {                       phi = PIover180 * i;
    TabCos360[i] = AppCos(phi);
    TabSin360[i] = AppSin(phi);
  }

  for (iCT = 0; iCT < nCT; iCT++)
  {
    if (CT->FobsRaw->Status != FRS_Active)
      InternalError("Corrupt nCT");

    h = CT->FobsRaw->h;
    k = CT->FobsRaw->k;
    l = CT->FobsRaw->l;

    CheckMalloc(CT->FF, uAtomType);

    for (i = 0; i < uAtomType; i++)
      CT->FF[i] = CalcFormFactor(&AtomType[i].SF_Info, h, k, l,
                                  AtomType[i].OccDefault,
                                  AtomType[i].UisoDefault);

    (void) BuildEq_hkl(&SpgrInfo, 1, &Eq_hkl, h, k, l);
    if (SgError != NULL) progerror(SgError);

    CT->nCTData = Eq_hkl.N;
    CheckMalloc(CT->CTData, CT->nCTData);
    CTD = CT->CTData;

    for (i = 0; i < Eq_hkl.N; i++)
    {
        CTD->Sign =  1;

      if (  hkl2ihkl( Eq_hkl.h[i],  Eq_hkl.k[i],  Eq_hkl.l[i],
                      &CTD->ih,     &CTD->ik,     &CTD->il) != 0)
      {
        CTD->Sign = -1;

        if (hkl2ihkl(-Eq_hkl.h[i], -Eq_hkl.k[i], -Eq_hkl.l[i],
                      &CTD->ih,     &CTD->ik,     &CTD->il) != 0)
          InternalError("OffGrid");
      }

      /* phi(Rs.H) = phi(H) - 2.pi.Ts.H */

      CTD->TH = -2 * Eq_hkl.TH[i];
      if (CT->FobsRaw->PhaseRestriction >= 0)
        CTD->TH += CT->FobsRaw->PhaseRestriction;
      CTD->TH = iModPositive(CTD->TH * 180 / STBF, 360);
      CTD++;
    }

    CT++;
  }
}


void DoCodeTranslation(T_PhaseCode *PhaseCode,
                       Fprec *Density, int Mx, int My, int Mz, int Friedel)
{
  int           iCT, iCTD;
  T_CodeTransl  *CT;
  T_CTData      *CTD;
  int           iTab, iD, jD;


  for (iD = 0; iD < Mx * My * Mz; iD++)
    Density[iD] = 0.;

  for (iCT = 0, CT = CodeTransl; iCT < nActivePhase; iCT++, CT++)
  {
    for (iCTD = 0, CTD = CT->CTData; iCTD < CT->nCTData; iCTD++, CTD++)
    {
      iTab = CTD->TH + (*PhaseCode);
      if (iTab >= 360) iTab -= 360;
      if (iCTD == 0) CT->Phi360Fmrg = iTab;
      if (iTab != 0 && CTD->Sign < 0) iTab = 360 - iTab;

      iD = (CTD->ih * My + CTD->ik) * Mz + CTD->il * 2;

#ifndef NO_EXTRA_CHECKS
      if (iTab < 0 || iTab >= 360)
        InternalError("Corrupt iTab");

      if (iD < 0 || iD + 1 >= Mx * My * Mz)
        InternalError(NULL);
#endif

      /* use conjugate complex because library fft does real to complex
         transform with "-1"
       */
      Density[iD    ] =  CT->FobsRaw->Fmrg * TabCos360[iTab];
      Density[iD + 1] = -CT->FobsRaw->Fmrg * TabSin360[iTab];

      if (Friedel == -1 || (Friedel == 1 && CTD->il == 0))
      {
        jD =   (((Nx - CTD->ih) % Nx) * My
             +  ((Ny - CTD->ik) % Ny)) * Mz
             +  ((Nz - CTD->il) % Nz) * 2;

#ifndef NO_EXTRA_CHECKS
        if (jD < 0 || jD + 1 >= Mx * My * Mz)
          InternalError(NULL);
#endif
        Density[jD    ] =  Density[iD    ];
        Density[jD + 1] = -Density[iD + 1];
      }

      if (Debug0) /* Debug: Print DCT */
      {
        int  h, k, l;

        if (iCTD == 0)
          Fprintf(stdout, "DCT %3d %3d %3d\n",
            CT->FobsRaw->h, CT->FobsRaw->k, CT->FobsRaw->l);

        ihkl2hkl(CTD->ih, CTD->ik, CTD->il, &h, &k, &l);

        Fprintf(stdout,
          "    %3d %3d %3d   %7.2f  %7.2f %7.2f  %c(%3d + %3d) = %3d\n",
          h, k, l,
          CT->FobsRaw->Fmrg, Density[iD], Density[iD + 1],
          (CTD->Sign < 0 ? '-' : ' '),
          CTD->TH, *PhaseCode, iTab);
      }
    }

    PhaseCode++;
  }
}


void CalcCCandR(T_CodeTransl *CT, int nCT, Fprec *CC, Fprec *R)
{
  int           iCT;
  T_CodeTransl  *CT0;
  Fprec         Diff, SumDiff, SumObs, SumCalc, MeanObs, MeanCalc;
  Fprec         Soc, So2, Sc2, o, c, denom;
  int           N;


  CT0 = CT;

  SumDiff = 0.;
  SumObs  = 0.;
  SumCalc = 0.;

  N = 0;

  for (iCT = 0, CT = CT0; iCT < nCT; iCT++, CT++)
  {
    Diff = CT->FcalAbs - CT->FobsRaw->Fmrg;
    SumDiff += CT->FobsRaw->M * AppFabs(Diff);

    SumObs  += CT->FobsRaw->M * CT->FobsRaw->Fmrg;
    SumCalc += CT->FobsRaw->M * CT->FcalAbs;

    N += CT->FobsRaw->M;
  }

  if (SumObs == 0.)
    *R = 9.99;
  else
    *R = SumDiff / SumObs;

  MeanObs  = SumObs  / N;
  MeanCalc = SumCalc / N;

  Soc = 0.;
  So2 = 0.;
  Sc2 = 0.;

  for (iCT = 0, CT = CT0; iCT < nCT; iCT++, CT++)
  {
    o = CT->FobsRaw->Fmrg - MeanObs;
    c = CT->FcalAbs       - MeanCalc;

    Soc += CT->FobsRaw->M * (o * c);
    So2 += CT->FobsRaw->M * (o * o);
    Sc2 += CT->FobsRaw->M * (c * c);
  }

      denom = sqrt(So2 * Sc2);
  if (denom == 0.)
    *CC = -2.;
  else
    *CC = Soc / denom;
}


void CheckRestrictedPhase(T_CodeTransl *CT, int *GC)
{
#define Tol 10
  if (   (CT->FcalAbs > (Fprec) 10. * BogFmrg) /* ARBITRARY */
      && (*GC > Tol && *GC < 180 - Tol || *GC > 180 + Tol && *GC < 360 - Tol))
#undef Tol
  {
    Fprintf(stdout, "\n%s %d %d %d  %.8g %.8g  %d\n\n",
      "# WARNING: Resetting illegal code for restricted phase",
      CT->FobsRaw->h, CT->FobsRaw->k, CT->FobsRaw->l,
      CT->FcalAbs, CT->FcalPhi, *GC);

    CountWarnings++;
  }

  if (*GC > 90 && *GC < 270) *GC = 180;
  else                       *GC =   0;
}


static int Find_iInputChemistry(int iAT, int jAT, int Type)
{
  int               iIC, Hit_iIC;
  T_InputChemistry  *IC;


  Hit_iIC = -1;

  IC = InputChemistry;

  for (iIC = 0; iIC < nInputChemistry; iIC++, IC++)
  {
    if (IC->Type == Type)
    {
      if (Type == ICT_MinDistance || Type == ICT_MaxDistance)
      {
        if (   (   AtomType[iAT].Class       == IC->AtomID[0].Class
                && AtomType[iAT].SF_Info.Lbl == IC->AtomID[0].SF_Info.Lbl
                && AtomType[jAT].Class       == IC->AtomID[1].Class
                && AtomType[jAT].SF_Info.Lbl == IC->AtomID[1].SF_Info.Lbl)
            || (   AtomType[iAT].Class       == IC->AtomID[1].Class
                && AtomType[iAT].SF_Info.Lbl == IC->AtomID[1].SF_Info.Lbl
                && AtomType[jAT].Class       == IC->AtomID[0].Class
                && AtomType[jAT].SF_Info.Lbl == IC->AtomID[0].SF_Info.Lbl))
        {
          if (Hit_iIC != -1)
            progerror("Chemistry: duplicate input");

          Hit_iIC = iIC;
        }
      }
      else
        InternalError("Type not implemented");
    }
  }

  return Hit_iIC;
}


static void BuildMinDistance2Mx(void)
{
  int    iAT, jAT, iIC;
  Fprec  *MD2Mx;


  MinOfMD2Mx = MaxLatticeTr2;
  MD2Mx = MinDistance2Mx;

  for (iAT = 0; iAT < uAtomType; iAT++)
  {
    for (jAT = 0; jAT < uAtomType; jAT++)
    {
      if (jAT < iAT)
        *MD2Mx = MinDistance2Mx[jAT * uAtomType + iAT];
      else
      {
            iIC = Find_iInputChemistry(iAT, jAT, ICT_MinDistance);
        if (iIC < 0)
        {
          *MD2Mx  = AtomType[iAT].eListAR->Radius;
          *MD2Mx += AtomType[jAT].eListAR->Radius;
          *MD2Mx *= (Fprec) .9; /* ARBITRARY */
          *MD2Mx *= (*MD2Mx);
        }
        else
        {
          *MD2Mx  = InputChemistry[iIC].Value;
          *MD2Mx *= (*MD2Mx);
        }
      }

      if (MinOfMD2Mx > *MD2Mx)
          MinOfMD2Mx = *MD2Mx;

      MD2Mx++;
    }
  }


  if (Debug1) /* Debug: Print MinDistance2Mx */
  {
    Fprintf(stdout, "MinDistance2Mx\n");

    MD2Mx = MinDistance2Mx;

    for (iAT = 0; iAT < uAtomType; iAT++)
    {
      Fprintf(stdout, "%s  %-4s",
        EchoAtomTypeClass(NULL, AtomType[iAT].Class),
        AtomType[iAT].SF_Info.Lbl);

      for (jAT = 0; jAT <= iAT; jAT++)
        Fprintf(stdout, " %7.3f", AppSqrt(MD2Mx[jAT]));
      putc('\n', stdout);

      MD2Mx += uAtomType;
    }
    putc('\n', stdout);

    Fprintf(stdout, "MinOfMD2Mx = %.3f\n\n", AppSqrt(MinOfMD2Mx));
  }
}


void SetupXtalInternals(void)
{
  CompleteCodeTransl(CodeTransl, nActivePhase);

  if (nFobsRaw > 0)
  {
    SumFmrg = CalcSumFmrg(CodeTransl, nActivePhase);

    if (Debug1) /* Debug: Print SumFmrg */
      Fprintf(stdout, "# SumFmrg = %.0f\n\n", SumFmrg);
  }

  if (F_SiteFrame == 0)
  {
    CheckMalloc(MinDistance2Mx, uAtomType * uAtomType);
    BuildMinDistance2Mx();
  }

  IRho3dVscale = LatConD.v / (Nx * Ny * Nz);
  InitCalcIRius();
}
