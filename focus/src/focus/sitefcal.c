#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "main.h"
#include "xtal.h"
#include "lattice.h"
#include "atominfo.h"
#include "io.h"


typedef struct
  {
    int        nSE;
    T_fVector  *SE;
  }
  T_ListSymEquiv;


typedef struct
  {
    int         h, k, l;
    int         restriction;
    CmplxFprec  Fcal;
    Fprec       Q;
  }
  T_ListFcal;


static T_ListSymEquiv  *ListSE = NULL;

static T_ListFcal  *ListFcal = NULL;
static int         nListFcal = 0;
static Fprec       MaxFcal;


void CompleteSite(T_Site *S, int nS, T_AtomType *AT, int nAT)
{
  int         iS, iAT;
  Fprec       Occ, Uiso;
  char        buf[128];
  T_AtomType  *eAT;


  for (iS = 0; iS < nS; iS++, S++)
  {
    if (S->SF_Info.Lbl)
    {
      if (CompleteSF_Info(&S->SF_Info, 1, 1) != 0)
      {
        Sprintf(buf, "Site  %s  %.60s: Unknown Scattering Factor Label",
          S->Label, S->SF_Info.Lbl);

        progerror(buf);
      }

      AppFree((char *) S->SF_Info.Lbl, strlen(S->SF_Info.Lbl) + 1);
    }
    else
    {
      S->SF_Info.Lbl = S->Label;

      if (CompleteSF_Info(&S->SF_Info, 0, 1) != 0)
      {
        Sprintf(buf, "Site  %s: Unknown Scattering Factor Label", S->Label);
        progerror(buf);
      }
    }

    if (S->SF_Info.SFT)
      S->SF_Info.Lbl = S->SF_Info.SFT->Label;
    else
      S->SF_Info.Lbl = S->SF_Info.CAA->Label;

    if (S->F_Occ == 0 || S->F_Uiso == 0)
    {
      for (iAT = 0, eAT = AT; iAT < nAT; iAT++, eAT++)
        if (eAT->SF_Info.Lbl == S->SF_Info.Lbl) break;

      if (iAT < nAT)
      {
        Occ  = eAT->OccDefault;
        Uiso = eAT->UisoDefault;
      }
      else
      {
        Occ  = Df_OccDefault;
        Uiso = Df_UisoDefault;
      }
    }

    if (S->F_Occ  == 0)
    {
      S->Occ  = Occ;
      S->F_Occ  = -1;
    }

    if (S->F_Uiso == 0)
    {
      S->Uiso = Uiso;
      S->F_Uiso = -1;
    }
  }
}


static void PrepListSE(void)
{
  int             iS, iSE, AtomsPerUnitCell;
  T_Site          *S;
  T_fVector       *SymEquiv;
  int             MaxSymEquiv;
  T_ListSymEquiv  *LSE;
  Fprec              Dist2ConsiderSame;
  Fprec           MaxDist2ConsideredSame;
  Fprec           MinDist2Distinct;


  if (ListSE != NULL)
    InternalError("ListSE != NULL");

  CheckMalloc(ListSE, nSite);

  MaxSymEquiv = SpgrInfo.OrderL;
  CheckMalloc(SymEquiv, MaxSymEquiv);

  Dist2ConsiderSame = .01 * .01; /* ARBITRARY */
  MaxDist2ConsideredSame = -1.;
  MinDist2Distinct = MaxLatticeTr2;

  LSE = ListSE;

  for (iS = 0, S = Site; iS < nSite; iS++, S++)
  {
    LSE->nSE = CalcSymEquiv(S->x, S->y, S->z,
                            SymEquiv, MaxSymEquiv,
                                Dist2ConsiderSame,
                            &MaxDist2ConsideredSame,
                            &MinDist2Distinct,
                            NULL, 0);

    CheckMalloc(LSE->SE, LSE->nSE);

    for (iSE = 0; iSE < LSE->nSE; iSE++)
    {
      LSE->SE[iSE].x = SymEquiv[iSE].x;
      LSE->SE[iSE].y = SymEquiv[iSE].y;
      LSE->SE[iSE].z = SymEquiv[iSE].z;
    }

    LSE++;
  }

  AppFree(SymEquiv, MaxSymEquiv);


  Fprintf(stdout, ">Begin UnitCellContents\n");

  AtomsPerUnitCell = 0;
  LSE = ListSE;

  for (iS = 0, S = Site; iS < nSite; iS++, S++)
  {
    Fprintf(stdout,
      " %4d  %-6s %7.4f  %7.4f  %7.4f  %-4.4s  %5.2f  %7.3f\n",
      LSE->nSE,
      S->Label,
      S->x, S->y, S->z,
      S->SF_Info.Lbl,
      S->Occ, S->Uiso);

    for (iSE = 1; iSE < LSE->nSE; iSE++)
      Fprintf(stdout, "              %7.4f  %7.4f  %7.4f\n",
        LSE->SE[iSE].x, LSE->SE[iSE].y, LSE->SE[iSE].z);

    AtomsPerUnitCell += LSE->nSE;
    LSE++;
  }

  Fprintf(stdout,   "# Atoms / UnitCell = %d\n",
    AtomsPerUnitCell);

  if (MaxDist2ConsideredSame >= 0.)
    Fprintf(stdout, "# MaxDist_ConsideredSame = %.4f\n",
      sqrt(MaxDist2ConsideredSame));

  Fprintf(stdout,   "# MinDist_Distinct       = %.4f\n",
    sqrt(MinDist2Distinct));

  Fprintf(stdout, ">End UnitCellContents\n\n");
}


static void FreeListSE(void)
{
  int             iS;
  T_ListSymEquiv  *LSE;


  LSE = ListSE;

  for (iS = 0; iS < nSite; iS++)
  {
    AppFree(LSE->SE, LSE->nSE);
    LSE++;
  }

  AppFree(ListSE, nSite); ListSE = NULL;
}


static void CalcFcal(CmplxFprec *Fcal, int h, int k, int l)
{
  int        iS, nSE;
  T_Site     *S;
  T_fVector  *SE;
  Fprec      fj, phase;


  Fcal->r =
  Fcal->i = 0.;

  for (iS = 0, S = Site; iS < nSite; iS++, S++)
  {
    fj = CalcFormFactor(&S->SF_Info, h, k, l, S->Occ, S->Uiso);

    nSE = ListSE[iS].nSE;
     SE = ListSE[iS].SE;

    while (nSE--)
    {
      phase = (Fprec) TwoPI * (h * SE->x + k * SE->y + l * SE->z); SE++;
      Fcal->r += fj * AppCos(phase);
      Fcal->i += fj * AppSin(phase);
    }
  }
}


static int ListFcalSortFunction(const T_ListFcal *a, const T_ListFcal *b)
{
#define EPS_Q (1.e-6) /* ARBITRARY */

  if (a->Q + EPS_Q < b->Q) return -1;
  if (a->Q - EPS_Q > b->Q) return  1;
  return hkl_Compare(a->h, a->k, a->l, b->h, b->k, b->l);

#undef EPS_Q
}


static void BuildListFcal(Fprec FcalMaxQ)
{
  int              h, k, l, pass, iListFcal;
  int              Minh, Mink, Minl;
  int              Maxh, Maxk, Maxl;
  int              restriction;
  Fprec            dmin, Qmax, Q;
  Fprec            FcalAbs2;
  CmplxFprec       Fcal;
  T_ListFcal       *lfcal;


  if (FcalMaxQ > 0.)
    dmin = 1. / AppSqrt(FcalMaxQ);
  else
    dmin = 0.;

  CalcMaxhkl(&LatConR, &dmin, LambdaLength, &Maxh, &Maxk, &Maxl);
  Qmax = 1. / Square(dmin);

  (void) SetListMin_hkl(&SpgrInfo, 1, Maxh, Maxk, Maxl, &Minh, &Mink, &Minl);

  MaxFcal = 0;
  nListFcal = iListFcal = 0;

  for (pass = 0; pass < 2; pass++)
  {
    if (pass == 1 && nListFcal > 0)
      CheckMalloc(ListFcal, nListFcal);

    lfcal = ListFcal;

    for (h = Minh; h <= Maxh; h++)
    for (k = Mink; k <= Maxk; k++)
    for (l = Minl; l <= Maxl; l++)
    {
          Q = Q_hkl(h, k, l, &LatConR);
      if (Q > Qmax)
        continue;

      if (IsSysAbsent_hkl(&SpgrInfo, h, k, l, &restriction) != 0)
        continue;
      if (IsHidden_hkl(&SpgrInfo, 1,
                       Minh, Mink, Minl,
                       Maxh, Maxk, Maxl,
                          h,    k,    l) != 0)
        continue;

      if (pass == 0)
        nListFcal++;
      else
      {
        if (iListFcal >= nListFcal)
          InternalError("iListFcal >= nListFcal");

        CalcFcal(&Fcal, h, k, l);

        FcalAbs2 = Fcal.r * Fcal.r + Fcal.i * Fcal.i;
        if (MaxFcal < FcalAbs2) MaxFcal = FcalAbs2;

        lfcal->h = h;
        lfcal->k = k;
        lfcal->l = l;
        lfcal->restriction = restriction;
        lfcal->Fcal.r = Fcal.r;
        lfcal->Fcal.i = Fcal.i;
        lfcal->Q = Q;
        lfcal++;

        iListFcal++;
      }
    }
  }

  MaxFcal = AppSqrt(MaxFcal);

  if (iListFcal != nListFcal)
    InternalError("iListFcal != nListFcal");

  if (nListFcal > 1)
    qsort((void *) ListFcal, nListFcal, sizeof (*ListFcal),
          (SortFunction) ListFcalSortFunction);
}


static void PrintListFcal(void)
{
  int              iListFcal, ss, i;
  T_ListFcal       *lfcal;
  Fprec            d, twotheta;
  Fprec            F, phi;


  Fprintf(stdout, ">Begin SiteFcal\n");

  Fprintf(stdout,
  "#  h   k   l      F         phi        A         B        d     2-Theta\n");

  lfcal = ListFcal;

  for (iListFcal = 0; iListFcal < nListFcal; iListFcal++)
  {
    twotheta = TwoThetaDeg(lfcal->Q);
    d = lfcal->Q; if (d != 0) d = 1. / sqrt(d);

    F = AppSqrt(Square(lfcal->Fcal.r) + Square(lfcal->Fcal.i));
    if (F > MaxFcal * (Fprec) 1.e-5) /* ARBITRARY */
      phi = AppAtan2(lfcal->Fcal.i, lfcal->Fcal.r);
    else
      phi = 0.;

    if      (InInterval(phi,            0., 1.e-5)) /* ARBITRARY */
      phi = 0.;
    else if (InInterval(AppFabs(phi), M_PI, 1.e-5)) /* ARBITRARY */
      phi = M_PI;

    Fprintf(stdout,
      " %3d %3d %3d %10.4f %9.4f %10.4f %9.4f %8.4f %7.3f",
      lfcal->h, lfcal->k, lfcal->l,
      F, phi / PIover180,
      lfcal->Fcal.r, lfcal->Fcal.i,
      d, twotheta);

    ss = Is_ss(&SpgrInfo, lfcal->h, lfcal->k, lfcal->l);
    if (ss == 0) ss = '/';
    else         ss = ':';

    if (lfcal->restriction >= 0)
    {
      i = lfcal->restriction * (180 / STBF);
      Fprintf(stdout, " %3d%c%3d", i, ss, i + 180);
    }
    else
      Fprintf(stdout, "    %c", ss);

    putc('\n', stdout);

    lfcal++;
  }

  Fprintf(stdout, ">End SiteFcal\n");
  putc('\n', stdout);
}


static void PrintList_hkl(void)
{
  int              iListFcal;
  T_ListFcal       *lfcal;
  Fprec            d, twotheta;


  Fprintf(stdout, ">Begin List_hkl\n");
  Fprintf(stdout, "#  h   k   l    d     2-Theta\n");

  lfcal = ListFcal;

  for (iListFcal = 0; iListFcal < nListFcal; iListFcal++)
  {
    if (lfcal->h || lfcal->k || lfcal->l)
    {
      twotheta = TwoThetaDeg(lfcal->Q);
      d = lfcal->Q; if (d != 0) d = 1. / sqrt(d);

      Fprintf(stdout, " %3d %3d %3d %8.4f %7.3f",
        lfcal->h, lfcal->k, lfcal->l, d, twotheta);

      putc('\n', stdout);
    }

    lfcal++;
  }

  Fprintf(stdout, ">End List_hkl\n");
  putc('\n', stdout);
}


static Fprec LPcorr(Fprec TTheta, Fprec Polra) /* See: STEPCO manual */
{
  Fprec  ThetaRad, sinT, cosT, cosTT, LPfac;


  ThetaRad = TTheta * PIover180 / 2.;
  sinT = AppSin(ThetaRad);
  cosT = AppCos(ThetaRad);
  cosTT = AppCos(TTheta * PIover180);

  /* I = K T F**2 LP */

      LPfac = sinT * sinT * cosT * (1. + Polra);
  if (LPfac == 0.)
    return 0.;

  LPfac = (2. * (1. + Polra * cosTT * cosTT)) / LPfac;

  return LPfac;
}


static Fprec CalcIrel(Fprec TTheta, int h, int k, int l, Fprec F2)
{
  int    M;
  Fprec  LP, Irel;


  LP = LPcorr(TTheta, ProfilePOLRA);
  M  = Mult_hkl(&SpgrInfo, 1, h, k, l);

  Irel = LP * (M / 2.) * F2;

  return Irel;
}


static void PrintListCollection(void)
{
  int         iListFcal;
  T_ListFcal  *LFc;
  Fprec       d, twotheta, Irel, RefrIrel, ScaleF;


  RefrIrel = 0.;

  if (ProfileReferenceRefl.Index < 0)
  {
    LFc = ListFcal;

    for (iListFcal = 0; iListFcal < nListFcal; iListFcal++, LFc++)
    {
      twotheta = TwoThetaDeg(LFc->Q);

      Irel = CalcIrel(twotheta, LFc->h, LFc->k, LFc->l,
                      Square(LFc->Fcal.r) + Square(LFc->Fcal.i));

      if (RefrIrel < Irel)
          RefrIrel = Irel;
    }
  }
  else
  {
    LFc = &ListFcal[ProfileReferenceRefl.Index];

    twotheta = TwoThetaDeg(LFc->Q);

    RefrIrel = CalcIrel(twotheta, LFc->h, LFc->k, LFc->l,
                        Square(LFc->Fcal.r) + Square(LFc->Fcal.i));
  }

  if (ProfileReferenceMax && RefrIrel != 0.)
    ScaleF = ProfileReferenceMax / RefrIrel;
  else
    ScaleF = 1.;

  Fprintf(stdout, ">Begin ListCollection\n");
  Fprintf(stdout, "#  h   k   l  2-Theta    d        I(rel)\n");

  LFc = ListFcal;

  for (iListFcal = 0; iListFcal < nListFcal; iListFcal++, LFc++)
  {
    if (! (LFc->h || LFc->k || LFc->l))
      continue;

        twotheta = TwoThetaDeg(LFc->Q);
    if (twotheta + ProfileStep * .5 < ProfileStart)
      continue;
    if (twotheta - ProfileStep * .5 > ProfileEnd)
      break;

    d = LFc->Q; if (d != 0) d = 1. / sqrt(d);

    Irel = CalcIrel(twotheta, LFc->h, LFc->k, LFc->l,
                    Square(LFc->Fcal.r) + Square(LFc->Fcal.i));

    Fprintf(stdout, " %3d %3d %3d  %7.3f %8.4f  %8.2f",
      LFc->h, LFc->k, LFc->l, twotheta, d, Irel * ScaleF);

    putc('\n', stdout);
  }

  Fprintf(stdout, ">End ListCollection\n");
  putc('\n', stdout);
}


static Fprec CalcFWHM(Fprec tanT)
{
  Fprec  FWHM;


  if (ProfileFWHM.Valid)
  {
        FWHM = ProfileFWHM.U + tanT * (ProfileFWHM.V + ProfileFWHM.W * tanT);
    if (FWHM >= 0.)
        FWHM = sqrt(FWHM);
  }
  else
    FWHM = 0.01;

  return FWHM;
}


static int SetReflProfCoef(Fprec Q, int h, int k, int l, Fprec F2,
                           Fprec *Angle, Fprec *Intensity,
                           Fprec *FWHM, Fprec *Asym)
{
  Fprec  tanT;


      *Angle = TwoThetaDeg(Q);
  if (*Angle == 0.)
    return -1;

  *Intensity = CalcIrel(*Angle, h, k, l, F2);

  tanT = AppTan(*Angle * .5 * PIover180);

      *FWHM = CalcFWHM(tanT);
  if (*FWHM < 0.)
  {
    Fprintf(stderr, "%s: Error: Illegal half width function parameters\n",
                    progn);
    AppExit(1);
  }

  if (ProfileAsym.Valid)
    *Asym = ProfileAsym.a1 + (ProfileAsym.a2 + ProfileAsym.a3 / tanT) / tanT;
  else
    *Asym = 0.;

  return 0;
}


static Fprec Q_TTheta(Fprec TTheta, Fprec Lambda)
{
  Fprec  Q = (Fprec) 2. * AppSin(TTheta * (Fprec)(.5 * PIover180)) / Lambda;
  return Q * Q;
}


static Fprec SetHighQ(Fprec TTreference)
{
  int    iPass, InSphere, PrevInSphere, nHigherQ;
  int       h,    k,    l;
  int    Maxh, Maxk, Maxl;
  int    Minh, Mink, Minl;
  Fprec  dmin, Qref, Qmax, Q, HigherQ[2];


      Qref = Q_TTheta(TTreference, LambdaLength);
  if (Qref > 0.)
    dmin = 1. / AppSqrt(Qref);
  else
    dmin = 0.;

  CalcMaxhkl(&LatConR, &dmin, LambdaLength, &Maxh, &Maxk, &Maxl);
  Qref = 1. / Square(dmin);
  Qmax = 4. / Square(LambdaLength);

  PrevInSphere = -1;

  for (iPass = 0;; iPass++)
  {
    if      (iPass % 3 == 0 && iPass != 0)
      Maxh++;
    else if (iPass % 3 == 1)
      Maxk++;
    else if (iPass % 3 == 2)
      Maxl++;

    if (iPass == 0 || iPass % 3 != 0)
      (void) SetListMin_hkl(&SpgrInfo, 1,
                            Maxh, Maxk, Maxl, &Minh, &Mink, &Minl);

    InSphere = 0;
    nHigherQ = 0;

    for (h = Minh; h <= Maxh; h++)
    for (k = Mink; k <= Maxk; k++)
    for (l = Minl; l <= Maxl; l++)
    {
          Q = Q_hkl(h, k, l, &LatConR);
      if (Q > Qmax)
        continue;

      InSphere++;

      if (Q <= Qref)
        continue;

      if (IsSysAbsent_hkl(&SpgrInfo, h, k, l, NULL) != 0)
        continue;
      if (IsHidden_hkl(&SpgrInfo, 1,
                       Minh, Mink, Minl,
                       Maxh, Maxk, Maxl,
                          h,    k,    l) != 0)
        continue;

      if      (nHigherQ == 0)
      {
         HigherQ[0] = Q;
        nHigherQ = 1;
      }
      else if (nHigherQ == 1 || HigherQ[1] > Q)
      {
         HigherQ[1] = Q;
        nHigherQ = 2;

        if (HigherQ[1] < HigherQ[0]) {
            HigherQ[1] = HigherQ[0];
            HigherQ[0] = Q;
        }
      }
    }

    if (nHigherQ == 2)
      return HigherQ[1];

    if (InSphere == PrevInSphere)
      break;

    PrevInSphere = InSphere;
  }

  if (Debug1)
    Fprintf(stdout, "# SetHighQ(): nHigherQ = %d\n\n", nHigherQ);

  return 0.;
}


static Fprec DetFcalMaxQ(int UseRefrRefl, Fprec PeakRange)
{
  int    TT, h, k, l;
  Fprec  TTl, TTe, tanT, FWHM;
  Fprec  Q, FcalMaxQ;


  if (ProfileEnd == 0.)
  {
        TTe = TwoThetaDeg(FobsMaxQ);
    if (TTe == 0.)
      TT = 179;
    else {
          TT = (int) TTe + 1;
      if (TT > 179)
          TT = 179;
    }

    ProfileEnd = (Fprec) TT;
  }

  TTe = ProfileGenEnd = ProfileEnd;

  if (UseRefrRefl)
  {
    if      (ProfileReferenceRefl.Mode == PRRM_hkl)
    {
      h = ProfileReferenceRefl.h;
      k = ProfileReferenceRefl.k;
      l = ProfileReferenceRefl.l;

      if (IsSysAbsent_hkl(&SpgrInfo, h, k, l, NULL) != 0)
        progerror("ProfileReferenceRefl hkl is systematically absent");

      Q = Q_hkl(h, k, l, &LatConR);

      TTe = TwoThetaDeg(Q);
    }
    else if (   ProfileReferenceRefl.Mode == PRRM_TTheta
             && ProfileReferenceRefl.TTheta > 0.)
    {
          Q = SetHighQ(ProfileReferenceRefl.TTheta);
      if (Q == 0.)
        TTe = 179.;
      else
        TTe = TwoThetaDeg(Q);
    }

    if (ProfileGenEnd < TTe)
        ProfileGenEnd = TTe;

    TTe = ProfileGenEnd;
  }

  if (PeakRange > 0.)
  {
    for (TT = (int) ProfileGenEnd; TT < 179; TT++)
    {
                          tanT = AppTan(TT * .5 * PIover180);
          FWHM = CalcFWHM(tanT);
      if (FWHM >= 0.)
      {
            TTl = TT - FWHM * .5 * PeakRange;
        if (TTl <= ProfileGenEnd && TTe < (Fprec)(TT + 1))
          TTe = (Fprec)(TT + 1);
      }
    }
  }

  FcalMaxQ = Q_TTheta(TTe, LambdaLength);

  return FcalMaxQ;
}


static void SetProfileReferenceReflIndex(void)
{
  int         iLFc, iLFcClose, h, k, l;
  Fprec       TT, DeltaTT, CloseDelta;
  T_ListFcal  *LFc;


  iLFcClose  = -1;

  if      (ProfileReferenceRefl.Mode == PRRM_hkl)
  {
    h = ProfileReferenceRefl.h;
    k = ProfileReferenceRefl.k;
    l = ProfileReferenceRefl.l;

    if (IsSysAbsent_hkl(&SpgrInfo, h, k, l, NULL) != 0)
      progerror("ProfileReferenceRefl hkl is systematically absent");

    LFc = ListFcal;

    for (iLFc = 0; iLFc < nListFcal; iLFc++, LFc++)
    {
      if (AreSymEquivalent_hkl(&SpgrInfo, 1,
                               h, k, l,
                               LFc->h, LFc->k, LFc->l)) {
        iLFcClose = iLFc;
        break;
      }
    }
  }
  else if (ProfileReferenceRefl.Mode == PRRM_TTheta)
  {
    CloseDelta = 0.;

    LFc = ListFcal;

    for (iLFc = 0; iLFc < nListFcal; iLFc++, LFc++)
    {
                TT = TwoThetaDeg(LFc->Q);
      DeltaTT = TT - ProfileReferenceRefl.TTheta;
      DeltaTT = AppFabs(DeltaTT);

      if (CloseDelta > DeltaTT || iLFcClose < 0)
      {
          iLFcClose  = iLFc;
          CloseDelta = DeltaTT;
      }
    }
  }

  ProfileReferenceRefl.Index = iLFcClose;

  if (Debug1 && iLFcClose >= 0)
    Fprintf(stdout, "# ProfileReferenceRefl.Index -> %d %d %d  %.6g\n\n",
      ListFcal[iLFcClose].h,
      ListFcal[iLFcClose].k,
      ListFcal[iLFcClose].l,
      TwoThetaDeg(ListFcal[iLFcClose].Q));
}


static Fprec PseudoVoigtAt(Fprec DeltaTT, Fprec FWHM, Fprec FracLorentz)
{
  Fprec  GFac, gau, lor, PV;


  if (FWHM == 0.)
  {
    if (DeltaTT == 0.)
      PV = 1.;
    else
      PV = 0.;
  }
  else
  {
    GFac = .5 / AppSqrt(-log(.5));

    gau = DeltaTT / (GFac * FWHM);
    gau *= gau;
    gau = exp(-gau);

    lor = DeltaTT / (.5 * FWHM);
    lor *= lor;
    lor = 1. / (1. + lor);

    PV = (FracLorentz * lor) + ((1. - FracLorentz) * gau);
  }

  return PV;
}


static void CalcProfile(Fprec *ProfileCounts, int nProfileSteps,
                        const T_StdPeak *StdPeak)
{
  int    iPS, iSPstep, iListFcal, Min_iPS, i;
  Fprec  x, RTT, R, TT, sym, asy;
  Fprec  TTl, TTh;
  Fprec  Angle, Intensity, FWHM, Asym, MinDelta;

  T_StdPeakC  *SPC;
  T_ListFcal  *LFc;

  Fprec  PeakRange;


  if (ProfilePeakShape == PPS_StdPeak) {
    PeakRange = StdPeak->Range;
    SPC = StdPeak->C;
  }
  else {
    PeakRange = PseudoVoigtPeakRange;
    SPC = NULL;
  }

  for (iPS = 0; iPS < nProfileSteps; iPS++)
    ProfileCounts[iPS] = 0.;

  LFc = ListFcal;

  for (iListFcal = 0; iListFcal < nListFcal; iListFcal++, LFc++)
  {
    if (SetReflProfCoef(LFc->Q, LFc->h, LFc->k, LFc->l,
                        Square(LFc->Fcal.r) + Square(LFc->Fcal.i),
                        &Angle, &Intensity, &FWHM, &Asym) != 0)
      continue;

    TTl = Angle - (FWHM * .5 * PeakRange + ProfileStep * .5);
    TTh = Angle + (FWHM * .5 * PeakRange + ProfileStep * .5);

        iPS = (int)(TTl / ProfileStep);
    if (iPS < 0)
        iPS = 0;

    if (FWHM == 0.)
    {
      MinDelta = 0.;
      Min_iPS  = -1;

      for (; iPS < nProfileSteps; iPS++)
      {
            TT = iPS * ProfileStep;
        if (TT > TTh)
          break;

        RTT = Angle - TT;
        RTT = AppFabs(RTT);

        if (MinDelta > RTT || Min_iPS < 0) {
            MinDelta = RTT;
            Min_iPS  = iPS;
        }
      }

          iPS = Min_iPS;
      if (iPS < 0)
          iPS = 0;
    }

    for (; iPS < nProfileSteps; iPS++)
    {
          TT = iPS * ProfileStep;
      if (TT > TTh)
        break;

      if (FWHM == 0.)
      {
        RTT = 0.;
        R   = 0.;
      }
      else
      {
        RTT = TT - Angle;
        R   = AppFabs(RTT) / (FWHM * .5);

        if (R > PeakRange)
          continue;
      }

      if (SPC)
      {
        for (iSPstep = 1; iSPstep < StdPeak->nSteps; iSPstep++)
          if (SPC[iSPstep].R >= R) break;

        if (iSPstep < StdPeak->nSteps)
        {
          i = iSPstep - 1;
          x = (R - SPC[i].R) / (SPC[iSPstep].R - SPC[i].R);
          sym = SPC[i].Sym + x * (SPC[iSPstep].Sym - SPC[i].Sym);
          asy = SPC[i].Asy + x * (SPC[iSPstep].Asy - SPC[i].Asy);
          sym *= Intensity;
          asy *= Asym;
          if (RTT > 0.) x = sym * (1. - asy);
          else          x = sym * (1. + asy);
          ProfileCounts[iPS] += x;
        }
      }
      else
      {
        x = PseudoVoigtAt(RTT, FWHM, PseudoVoigtFracLorentz);
        x *= Intensity;
        ProfileCounts[iPS] += x;
      }

      if (FWHM == 0.)
        break;
    }
  }
}


static Fprec ProfileScaleFactor(const Fprec *ProfileCounts, int nProfileSteps,
                                const T_StdPeak *StdPeak)
{
  int               iPS;
  Fprec             TT, PC_at_hkl, ScaleF;
  const Fprec       *PC, *RefrPC;
  const T_ListFcal  *LFc;


  RefrPC = NULL;

  if (ProfileReferenceRefl.Index >= 0)
  {
    LFc = &ListFcal[ProfileReferenceRefl.Index];

    if      (ProfileReferenceRefl.Mode == PRRM_hkl)
    {
      if (ProfilePeakShape == PPS_StdPeak)
      {
        if (   StdPeak->nSteps < 1
            || StdPeak->C[0].R   != 0.
            || StdPeak->C[0].Asy != 0.)
          InternalError("Corrupt StdPeak");

        ScaleF = StdPeak->C[0].Sym;
      }
      else
        ScaleF = 1.;

      TT = TwoThetaDeg(LFc->Q);

      PC_at_hkl = CalcIrel(TT, LFc->h, LFc->k, LFc->l,
                           Square(LFc->Fcal.r) + Square(LFc->Fcal.i));

      PC_at_hkl *= ScaleF;

      RefrPC = &PC_at_hkl;
    }
    else if (ProfileReferenceRefl.Mode == PRRM_TTheta)
    {
      TT = TwoThetaDeg(LFc->Q);

          iPS = INT_ROUNDED(TT / ProfileStep);
      if (iPS < 0 || iPS >= nProfileSteps)
        InternalError("Corrupt ProfileReferenceRefl.Index");

      RefrPC = &ProfileCounts[iPS];

      if (   iPS - 1 >= 0
          && *RefrPC <  ProfileCounts[iPS - 1])
              RefrPC = &ProfileCounts[iPS - 1];

      if (   iPS + 1 < nProfileSteps
          && *RefrPC <  ProfileCounts[iPS + 1])
              RefrPC = &ProfileCounts[iPS + 1];
    }
  }
  else
  {
    PC = ProfileCounts;

    for (iPS = 0; iPS < nProfileSteps; iPS++, PC++)
    {
          TT = iPS * ProfileStep;
      if (TT + ProfileStep * .5 < ProfileStart)
        continue;
      if (TT - ProfileStep * .5 > ProfileEnd)
        break;

      if (    RefrPC == NULL
          || *RefrPC < *PC)
              RefrPC =  PC;
    }

    if (RefrPC == NULL)
    {
      PC = ProfileCounts;

      for (iPS = 0; iPS < nProfileSteps; iPS++, PC++) {
        if (    RefrPC == NULL
            || *RefrPC < *PC)
                RefrPC =  PC;
      }
    }
  }

  if (ProfileReferenceMax && RefrPC && *RefrPC != 0.)
    ScaleF = (ProfileReferenceMax - ProfileBackground) / (*RefrPC);
  else
    ScaleF = 1.;

  return ScaleF;
}


static void PrintProfile(const Fprec *ProfileCounts, int nProfileSteps,
                         const T_StdPeak *StdPeak)
{
  int          iPS;
  Fprec        TT, Counts, ScaleF, MaxScaledCounts, ESD;


  Fprintf(stdout, ">Begin stepscan\n");

  Fprintf(stdout, "@with g0\n");

  ScaleF = ProfileScaleFactor(ProfileCounts, nProfileSteps, StdPeak);

  MaxScaledCounts = 0.;

  for (iPS = 0; iPS < nProfileSteps; iPS++)
  {
        TT = iPS * ProfileStep;
    if (TT + ProfileStep * .5 < ProfileStart)
      continue;
    if (TT - ProfileStep * .5 > ProfileEnd)
      break;

    Counts = ProfileBackground + ProfileCounts[iPS] * ScaleF;

        ESD = AppSqrt(Counts);
    if (ESD < 0.001)
        ESD = 0.001;

    Fprintf(stdout, "%12.3f %12.2f %12.3f\n",
      TT, Counts, ESD);

    if (MaxScaledCounts < Counts)
        MaxScaledCounts = Counts;
  }

  if (MaxScaledCounts == 0.)
      MaxScaledCounts = ProfileReferenceMax;
  if (MaxScaledCounts == 0.)
      MaxScaledCounts = 1.;

  Fprintf(stdout, "&\n");
  Fprintf(stdout, "@autoscale\n");

  Fprintf(stdout, "@ world xmin %.6g\n", ProfileStart);
  Fprintf(stdout, "@ world xmax %.6g\n", ProfileEnd);
  Fprintf(stdout, "@ world ymin %.6g\n", 0.);
  Fprintf(stdout, "@ world ymax %.6g\n", MaxScaledCounts);

  Fprintf(stdout, ">End stepscan\n");
  putc('\n', stdout);
}


void DoSiteFcal(int OutputMode)
{
  Fprec  FcalMaxQ;


  FcalMaxQ = DetFcalMaxQ(0, 0.);

  PrepListSE();
  BuildListFcal(FcalMaxQ);

  if (OutputMode == 0)
    PrintListFcal();
  else
    PrintList_hkl();

  AppFree(ListFcal, nListFcal);
          ListFcal = NULL;
                    nListFcal = 0;
  FreeListSE();
}


void CollectionLists(int FlagPrintList_hkl,
                     int PowderStepScan, const char *fnStdPeak)
{
  Fprec      FcalMaxQ;
  int        nProfileSteps;
  Fprec      *ProfileCounts, PeakRange;
  T_StdPeak  *StdPeak, StdPeakBuf[1];


  StdPeak = NULL;

  if (PowderStepScan)
  {
    if (ProfilePeakShape == PPS_StdPeak)
    {
                             StdPeak = StdPeakBuf;
      LoadStdPeak(fnStdPeak, StdPeak);

      PeakRange = StdPeak->Range;
    }
    else
      PeakRange = PseudoVoigtPeakRange;

    FcalMaxQ = DetFcalMaxQ(1, PeakRange);
  }
  else
    FcalMaxQ = DetFcalMaxQ(1, 0.);

  EchoProfileSettings(stdout);

  PrepListSE();
  BuildListFcal(FcalMaxQ);

  Fprintf(stdout, "# FcalMaxQ = %.6g => 2-theta %.6g\n",
    FcalMaxQ, TwoThetaDeg(FcalMaxQ));
  Fprintf(stdout, "# nListFcal = %d\n\n", nListFcal);

  SetProfileReferenceReflIndex();

  if (FlagPrintList_hkl)
    PrintListCollection();

  if (PowderStepScan)
  {
    nProfileSteps = INT_ROUNDED(ProfileGenEnd / ProfileStep) + 1;
    CheckMalloc(ProfileCounts, nProfileSteps);

    Fprintf(stdout, "# nProfileSteps = %d\n\n", nProfileSteps);

    CalcProfile(ProfileCounts, nProfileSteps, StdPeak);
    PrintProfile(ProfileCounts, nProfileSteps, StdPeak);

    AppFree(ProfileCounts, nProfileSteps);
  }

  AppFree(ListFcal, nListFcal);
          ListFcal = NULL;
                    nListFcal = 0;
  FreeListSE();

  if (StdPeak) {
    AppFree(StdPeak->Title, strlen(StdPeak->Title) + 1);
    AppFree(StdPeak->C, StdPeak->nSteps);
  }
}


static void SiteCT_Fcal(T_CodeTransl *CT, int nCT)
{
  int         iCT;
  int         h, k, l;
  CmplxFprec  Fcal;


  for (iCT = 0; iCT < nCT; iCT++, CT++)
  {
    h = CT->FobsRaw->h;
    k = CT->FobsRaw->k;
    l = CT->FobsRaw->l;
    Fcal.r =
    Fcal.i = 0.;

    CalcFcal(&Fcal, h, k, l);

    CT->Fcal.r = Fcal.r;
    CT->Fcal.i = Fcal.i;
    CT->FcalAbs = AppSqrt(Fcal.r * Fcal.r + Fcal.i * Fcal.i);
    if (CT->FcalAbs > BogFmrg)
      CT->FcalPhi = AppAtan2(Fcal.i, Fcal.r);
    else
      CT->FcalPhi = 0.;

    if (Debug0) /* Debug: Print h k l Fc */
    {
      Fprintf(stdout, "SiteCT_Fc  %3d %3d %3d  %10.4f  %10.4f\n",
        h, k, l, CT->Fcal.r, CT->Fcal.i);
    }
  }
}


void DoSitePhases(void)
{
  int           iCT, iCTD, code;
  T_CodeTransl  *CT;
  T_CTData      *CTD;
  Fprec         CC, R;


  if (nSite < 1) return;

  PrepListSE();

  SiteCT_Fcal(CodeTransl, nActivePhase);
  CalcCCandR(CodeTransl, nActivePhase, &CC, &R);

  Fprintf(stdout, "# SiteCT_Fcal  CC = %.6g  R = %.3f\n\n", CC, R);

  CT = CodeTransl;

  for (iCT = 0; iCT < nActivePhase; iCT++)
  {
    Phi2Phi360(CT->FcalPhi, CT->Phi360Fcal);

    if (Debug0) /* Debug: Print SitePhases */
    Fprintf(stdout, "SP  %3d %3d %3d  %3d  %8.3f\n",
      CT->FobsRaw->h, CT->FobsRaw->k, CT->FobsRaw->l,
      CT->Phi360Fcal, CT->FcalPhi / PIover180);

    if (CT->nCTData < 1)
      InternalError("CT->nCTData < 1");

    CTD = CT->CTData;

    code = CT->Phi360Fcal - CTD->TH;
    if (code < 0) code += 360;

    if (CT->FobsRaw->PhaseRestriction >= 0)
    {
      if (code != 0 && code != 180)
        CheckRestrictedPhase(CT, &code);
    }

    for (iCTD = 0; iCTD < CT->nCTData; iCTD++)
    {
      CTD->TH += code;
      if (CTD->TH >= 360) CTD->TH -= 360;

      CTD++;
    }

    CT++;
  }

  FreeListSE();
}
