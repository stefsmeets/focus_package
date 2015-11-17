/*
  Space Group Info (c) 1994-97 Ralf W. Grosse-Kunstleve
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define SGCOREDEF__
#include "sginfo.h"


#define Fputs   (void) fputs
#define Fprintf (void) fprintf
#define Fflush  (void) fflush

#define MemCpy(t, s, n)  memcpy((t), (s), (n) * sizeof (*(t)))
#define MemCmp(t, s, n)  memcmp((t), (s), (n) * sizeof (*(t)))


static const char *progn = "sgdev";


static void progerror(const char *message)
{
  Fflush(stdout);
  Fprintf(stderr, "%s: %s\n", progn, message);
  exit(1);
}


#define PMG_Unknown  0
#define PMG_1        1
#define PMG_2        2
#define PMG_222      3
#define PMG_4        4
#define PMG_422      5
#define PMG_3        6
#define PMG_321      7
#define PMG_312      8
#define PMG_6        9
#define PMG_622      10
#define PMG_23       11
#define PMG_432      12


typedef struct
  {
    int  mList;
    int  nList;
    int  *R;
  }
  T_ListRMx;


void ResetListRMx(T_ListRMx *ListRMx)
{
  ListRMx->nList = 0;
}


void InitListRMx(T_ListRMx *ListRMx)
{
  ListRMx->mList = 0;
  ListRMx->R = NULL;
  ResetListRMx(ListRMx);
}


size_t AllocateListRMx(T_ListRMx *ListRMx)
{
  size_t  s;

  if (ListRMx->mList <= 0)
    return 0;
                          s = ListRMx->mList * 9 * sizeof (*ListRMx->R);
      ListRMx->R = malloc(s);
  if (ListRMx->R == NULL)
    return 0;

  ResetListRMx(ListRMx);

  return s;
}


void FreeListRMx(T_ListRMx *ListRMx)
{
  if (ListRMx->R) free(ListRMx->R);
  InitListRMx(ListRMx);
}


static int CoreAdd2ListRMx(T_ListRMx *ListRMx, const int *NewR)
{
  int  iList;
  int  *LR;

  static const char *Err_NonXtalOp =
    "Error: Generators produce non-crystallographic operation";


  LR = ListRMx->R;

  for (iList = 0; iList < ListRMx->nList; iList++, LR += 9)
    if (MemCmp(LR, NewR, 9) == 0)
      return 0; /* matrix is not unique */

  if (SetRotMxInfo(NewR, NULL) == 0) {
    SetSgError(Err_NonXtalOp);
    return -1;
  }

  if (ListRMx->nList >= ListRMx->mList) {
    SetSgError("Internal Error: Allocated space for ListRMx->R too small");
    return -1;
  }

  (void) MemCpy(LR, NewR, 9);

  ListRMx->nList++;

  return 1;
}


int Add2ListRMx(T_ListRMx *ListRMx, const int *NewR)
{
  int   iMult,    jMult;
  int  *iMultR, *jMultR, TrialR[9];


  if (ListRMx->nList == 0)
  {
    /* make sure identity matrix is first in list */

    InitRotMx(TrialR, 1);

    if (CoreAdd2ListRMx(ListRMx, TrialR) < 0)
      return -1;
  }

  (void) MemCpy(TrialR, NewR, 9);

  iMult  =  ListRMx->nList;
  iMultR = &ListRMx->R[iMult * 9];

  jMult  = 1;
  jMultR = &ListRMx->R[jMult * 9]; /* skip first = identity matrix */

  for (;;)
  {
    if (CoreAdd2ListRMx(ListRMx, TrialR) < 0)
      return -1;

    if (jMult > iMult)
    {
      iMult++;
      iMultR += 9;

      jMult  = 1;
      jMultR = &ListRMx->R[jMult * 9]; /* skip first = identity matrix */
    }

    if (iMult == ListRMx->nList)
      break;

    RotMxMultiply(TrialR, jMultR, iMultR);

    jMult++;
    jMultR += 9;
  }

  return 0;
}


static int iRRE_FreeParameters(const int *Mx, int nr, int nc, int *FreeP)
{
  int  pr, pc;
  int  nFreeP;


#define Mx(i, j) Mx[i * nc + j]

  nFreeP = 0;

  pc = 0;

  for (pr = 0; pr < nr && pc < nc; pr++, pc++)
  {
    for (; pc < nc; pc++) {
      if (Mx(pr, pc) != 0) break;
      FreeP[nFreeP] = pc;
            nFreeP++;
    }
  }

#undef Mx

  return nFreeP;
}


static int iRRE_BackSubstitution(const int *Mx, int nr, int nc, const int *V,
                                 const int *FreePi, int nFreeP,
                                 const int *FreePv,
                                 int *SolNom, int *SolDen)
{
  int  pr, pc, ic, i;


#define Mx(i, j) Mx[i * nc + j]

  for (ic = 0; ic < nc; ic++) {
    SolNom[ic] = 0;
    SolDen[ic] = 0;
  }

  for (i = 0; i < nFreeP; i++) {
    SolNom[FreePi[i]] = FreePv[i];
    SolDen[FreePi[i]] = 1;
  }

  for (pr = nr - 1; pr >= 0; pr--)
  {
    for (pc = 0; pc < nc; pc++)
      if (Mx(pr, pc) != 0) goto Substitute;

    continue;

    Substitute:

    SolNom[pc] = V[pr];
    SolDen[pc] = Mx(pr, pc);

    for (ic = pc + 1; ic < nc; ic++)
      SolNom[pc] -= Mx(pr, ic) * SolNom[ic];
  }

  for (ic = 0; ic < nc; ic++)
    if (SolDen[ic] == 0) return -1;

#undef Mx

  return 0;
}


static void AddFracs(const int N1, const int N2,
                     const int D1, const int D2,
                     int *NS,
                     int *DS)
{
  int lcm = iLCM(D1, D2);
  *NS = (N1 * lcm) / D1 + (N2 * lcm) / D2;
  *DS = lcm;
  SimplifyFraction(*NS, *DS, NS, DS);
}


static int iDetFrac(const int *MxN, const int *MxD, int *DetN, int *DetD)
{
  int  N[3], D[3], i;


  AddFracs(MxN[4] * MxN[8], -MxN[5] * MxN[7],
           MxD[4] * MxD[8],  MxD[5] * MxD[7], &N[0], &D[0]);
  N[0] *=  MxN[0];
  D[0] *=  MxD[0];
  AddFracs(MxN[3] * MxN[8], -MxN[5] * MxN[6],
           MxD[3] * MxD[8],  MxD[5] * MxD[6], &N[1], &D[1]);
  N[1] *= -MxN[1];
  D[1] *=  MxD[1];
  AddFracs(MxN[3] * MxN[7], -MxN[4] * MxN[6],
           MxD[3] * MxD[7],  MxD[4] * MxD[6], &N[2], &D[2]);
  N[2] *=  MxN[2];
  D[2] *=  MxD[2];

  for (i = 0; i < 3; i++)
    SimplifyFraction(N[i], D[i], &N[i], &D[i]);

  for (i = 1; i < 3; i++)
    AddFracs(N[0], N[i],
             D[0], D[i], &N[0], &D[0]);

  (*DetN) = N[0];
  (*DetD) = D[0];

  return (*DetN);
}


static int NextLoop_n_m(int mf, int ml, int n, int *ix)
{
  int  p;


  for (p = 0; p < n; p++)
  {
        ix[p]++;
    if (ix[p] <= ml) return 1;
        ix[p] =  mf;
  }

  return 0;
}



static const int *ReferenceRotMx(int Order, const int **PseudoG)
{
  static const int RRMx[10][9] =
  {
    {  1,  0,  0,  /*  1 */
       0,  1,  0,
       0,  0,  1
    },
    { -1,  0,  0,  /* -1 */
       0, -1,  0,
       0,  0, -1
    },
    { -1,  0,  0,  /*  2 */
       0, -1,  0,
       0,  0,  1
    },
    {  1,  0,  0,  /* -2 */
       0,  1,  0,
       0,  0, -1
    },
    {  0, -1,  0,  /*  3 */
       1, -1,  0,
       0,  0,  1
    },
    {  0,  1,  0,  /* -3 */
      -1,  1,  0,
       0,  0, -1
    },
    {  0, -1,  0,  /*  4 */
       1,  0,  0,
       0,  0,  1
    },
    {  0,  1,  0,  /* -4 */
      -1,  0,  0,
       0,  0, -1
    },
    {  1, -1,  0,  /*  6 */
       1,  0,  0,
       0,  0,  1
    },
    { -1,  1,  0,  /* -6 */
      -1,  0,  0,
       0,  0, -1
    }
  };

  int ix;

  ix = Order;
  if (Order < 0) ix *= -1;
  if (ix == 0 || ix == 5 || ix > 6) return NULL;

  if (PseudoG)
  {
    if (ix == 3 || ix == 6)
      *PseudoG = PseudoMetricalMatrix(XS_Hexagonal, '=', 'z', 0);
    else
      *PseudoG = PseudoMetricalMatrix(XS_Cubic,     '=', 'z', 0);
  }

  if (ix > 4) ix = 5;
  ix--;
  ix *= 2;
  if (Order < 0) ix++;

  return RRMx[ix];
}


#define LRBF 12


static int L_CB_SMx(T_RTMx *CSiC,
                    const T_RTMx *CBMx,
                    const T_RTMx *SMx,
                    const T_RTMx *InvCBMx)
{
  int     i;
  T_RTMx  BufMx;


  RTMxMultiply(&BufMx, SMx,  InvCBMx, CTBF / STBF, 0);
  RTMxMultiply(CSiC,   CBMx, &BufMx,  LRBF,        LRBF * CTBF);

  for (i = 0; i < 9; i++)
  {
    if (CSiC->s.R[i] % (LRBF * LRBF)) {
      progerror("Internal Error: Corrupt CBMx/SMx/InvCBMx");
      return -1;
    }

    CSiC->s.R[i] /= (LRBF * LRBF);
  }

  for (i = 0; i < 3; i++)
  {
    if (CSiC->s.T[i] % (LRBF * (CTBF / STBF))) {
      progerror("Internal Error: Out of STBF range");
      return -1;
    }

    CSiC->s.T[i] /= (LRBF * (CTBF / STBF));
  }

  return 0;
}


static int CmpLen(const int *a, const int *b)
{
  if (*a < *b) return -1;
  if (*a > *b) return  1;
  return 0;
}


static void AllCBMx(const int *Mx, int nr, const int *V,
                    const int FreePi[9], int nFreeP,
                    const int Mi[9], const int Mw[9], const int *PseudoGw,
                    int Order)
{
  int        nCBMx, i, j;
  int        FreePv[9], Range;
  int        SolNom[9], SolDen[9];
  int        DetN, DetD, Det;
  T_RTMx     Tiw[1], Twi[1], SMxw[1], SMxi[1];
  int        RtGwR[9], Rt[9], GwR[9], Len[3];
  const int  *PseudoGi;
  int        NewBest, BestTwi[9], BestDet, BestLen[3], MaxBestTwiElement;


  /*
     Mw=Tiw.Mi.Twi
     Gi=Tiwt.Gw.Tiw
   */

  for (i = 0; i < 9; i++)
    SMxi->s.R[i] = Mi[i];

  for (i = 0; i < 3; i++) {
    SMxi->s.T[i] = 0;
     Tiw->s.T[i] = 0;
     Twi->s.T[i] = 0;
  }

  BestDet = 0;
  for (i = 0; i < 3; i++) BestLen[i] = 0;
  for (i = 0; i < 9; i++) BestTwi[i] = 0;

  nCBMx = 0;
  PseudoGi = NULL;
  Range = 2;

  if (nFreeP == 2)
    Range = 4;

  for (i = 0; i < nFreeP; i++) FreePv[i] = -Range;

  do
  {
    if (iRRE_BackSubstitution(Mx, nr, 9, V, FreePi, nFreeP,
                              FreePv, SolNom, SolDen) != 0)
      progerror("Corrupt iReducedRowEchelon()-Mx or iRRE_BackSubstitution()");

    for (i = 0; i < 9; i++) {
      SimplifyFraction(SolNom[i], SolDen[i], &SolNom[i], &SolDen[i]);
      if (SolDen[i] != 1) break;
    }

    if (i == 9 && iDetFrac(SolNom, SolDen, &DetN, &DetD) > 0)
    {
      Fprintf(stdout, "D %s ", FormatFraction(DetN, DetD, 0, NULL, 0));
      MathMx(SolNom, SolDen, 0, 1, 9, 1, 0);
      putc('\n', stdout);
      Fprintf(stdout, "twi=");
      MathMx(SolNom, SolDen, 0, 3, 3, 1, 0);
      putc('\n', stdout);

      for (i = 0; i < 9; i++) {
            Twi->s.R[i] =  SolNom[i] * LRBF;
        if (Twi->s.R[i] %  SolDen[i]) goto NextSolution;
            Twi->s.R[i] /= SolDen[i];
      }

          Det = deterRotMx(Twi->s.R);
      if (Det %  (LRBF * LRBF)) goto NextSolution;
          Det /= (LRBF * LRBF);
      if (Det != DetN * LRBF / DetD) progerror("Corrupt DetN/D");

      iCoFactorMxTp(Twi->s.R, Tiw->s.R);

      Fprintf(stdout, "tiw=");
      MathMx(Tiw->s.R, NULL, LRBF * Det, 3, 3, 1, 0);
      putc('\n', stdout);

      for (i = 0; i < 9; i++) {
        if (Tiw->s.R[i] %  Det) goto NextSolution;
            Tiw->s.R[i] /= Det;
      }

      (void) L_CB_SMx(SMxw, Tiw, SMxi, Twi);

      for (i = 0; i < 9; i++)
        if (SMxw->s.R[i] != Mw[i]) {
          for (i = 0; i < 9; i++)
            Fprintf(stdout, "%2d != %2d\n", SMxw->s.R[i], Mw[i]);
          progerror("Error solving Mi.T=T.Mw");
        }

      if (PseudoGi == NULL)
      {
        if (PseudoGw == NULL) {
          InitRotMx(RtGwR, LRBF * LRBF);
          PseudoGi = RtGwR;
        }
        else
        {
          for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
              Rt[i * 3 + j] = Tiw->s.R[j * 3 + i];

          Fprintf(stdout, "tiwt=");
          MathMx(Rt, NULL, LRBF, 3, 3, 1, 0);
          putc('\n', stdout);

          RotMxMultiply(GwR, PseudoGw, Tiw->s.R);
          RotMxMultiply(RtGwR, Rt, GwR);

          if (RtGwR[0] < 0 || RtGwR[4] < 0 || RtGwR[8] < 0)
            progerror("Error computing PseudoGi");

          PseudoGi = RtGwR;

          Fprintf(stdout, "gw=");
          MathMx(PseudoGw, NULL,           1, 3, 3, 1, 0);
          putc('\n', stdout);
          Fprintf(stdout, "gi=");
          MathMx(PseudoGi, NULL, LRBF * LRBF, 3, 3, 1, 0);
          putc('\n', stdout);
        }
      }

      for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
          Rt[i * 3 + j] = Twi->s.R[j * 3 + i];

      Fprintf(stdout, "twit=");
      MathMx(Rt, NULL, LRBF, 3, 3, 1, 0);
      putc('\n', stdout);

      for (i = 0; i < 3; i++)
        Len[i] = iScalProd(&Rt[i * 3], &Rt[i * 3], PseudoGi);

      MathMx(Len, NULL, LRBF * LRBF * LRBF * LRBF, 1, 3, 1, 0);
      putc('\n', stdout);

      qsort((void *) Len, 3, sizeof (*Len),
            (int (*)(const void *, const void *)) CmpLen);

      MathMx(Len, NULL, LRBF * LRBF * LRBF * LRBF, 1, 3, 1, 0);
      putc('\n', stdout);

      NewBest = 0;

      if      (BestDet > Det || BestDet == 0)
        NewBest = 1;
      else if (BestDet == Det)
      {
        for (i = 0; i < 3; i++)
          if      (BestLen[i] > Len[i]) {
            NewBest = 1;
            break;
          }
          else if (BestLen[i] < Len[i])
            break;
      }

      if (NewBest)
      {
        BestDet = Det;
        for (i = 0; i < 3; i++) BestLen[i] = Len[i];
        for (i = 0; i < 9; i++) BestTwi[i] = Twi->s.R[i];
      }

      nCBMx++;
    }

    NextSolution:;
  }
  while (NextLoop_n_m(-Range, Range, nFreeP, FreePv) != 0);

  if (BestDet)
  {
    MaxBestTwiElement = 0;

    for (i = 0; i < 9; i++)
      if (MaxBestTwiElement < abs(BestTwi[i]))
          MaxBestTwiElement = abs(BestTwi[i]);

    Fprintf(stdout, "Order=%d Best Det=%s",
      Order, FormatFraction(BestDet, LRBF, 0, NULL, 0));
    Fprintf(stdout, " MaxTwiElement=%s Twi=",
      FormatFraction(MaxBestTwiElement, LRBF, 0, NULL, 0));
    MathMx(BestTwi, NULL, LRBF, 3, 3, 1, 0);
    Fprintf(stdout, " Len=");
    MathMx(BestLen, NULL, LRBF * LRBF * LRBF * LRBF, 1, 3, 1, 0);
    putc('\n', stdout);
  }

  Fprintf(stdout, "nCBMx = %d\n", nCBMx);
}


static void SetupMiT_TMw(int *MiT_TMw, int nr,
                         const int Mi[9], const int Mw[9])
{
  int  ib, ir, ic, ix;

  /*
     known: Mi = "is"     matrix
            Mw = "wanted" matrix

     unknown: T = change of basis matrix Xi=T.Xw

     Mi=T.Mw.T^-1 => Mi.T-T.Mw=0

     Mi.T => {{  Mi0,    0,    0,  Mi1,    0,    0,  Mi2,    0,    0},
              {    0,  Mi0,    0,    0,  Mi1,    0,    0,  Mi2,    0},
              {    0,    0,  Mi0,    0,    0,  Mi1,    0,    0,  Mi2},
              {  Mi3,    0,    0,  Mi4,    0,    0,  Mi5,    0,    0},
              {    0,  Mi3,    0,    0,  Mi4,    0,    0,  Mi5,    0},
              {    0,    0,  Mi3,    0,    0,  Mi4,    0,    0,  Mi5},
              {  Mi6,    0,    0,  Mi7,    0,    0,  Mi8,    0,    0},
              {    0,  Mi6,    0,    0,  Mi7,    0,    0,  Mi8,    0},
              {    0,    0,  Mi6,    0,    0,  Mi7,    0,    0,  Mi8}}

    -T.Mw => {{ -Mw0, -Mw3, -Mw6,    0,    0,    0,    0,    0,    0},
              { -Mw1, -Mw4, -Mw7,    0,    0,    0,    0,    0,    0},
              { -Mw2, -Mw5, -Mw8,    0,    0,    0,    0,    0,    0},
              {    0,    0,    0, -Mw0, -Mw3, -Mw6,    0,    0,    0},
              {    0,    0,    0, -Mw1, -Mw4, -Mw7,    0,    0,    0},
              {    0,    0,    0, -Mw2, -Mw5, -Mw8,    0,    0,    0},
              {    0,    0,    0,    0,    0,    0, -Mw0, -Mw3, -Mw6},
              {    0,    0,    0,    0,    0,    0, -Mw1, -Mw4, -Mw7},
              {    0,    0,    0,    0,    0,    0, -Mw2, -Mw5, -Mw8}}
   */

  for (ix = 0; ix < nr * 9; ix++) MiT_TMw[ix] = 0;

  for (ib = 0; ib < 3; ib++)
    for (ir = 0; ir < 3; ir++)
      for (ic = 0; ic < 3; ic++)
        MiT_TMw[ib * 10 + ir * 3 * 9 + ic * 3] = Mi[ir * 3 + ic];

  for (ib = 0; ib < 3; ib++)
    for (ir = 0; ir < 3; ir++)
      for (ic = 0; ic < 3; ic++)
        MiT_TMw[ib * 10 * 3 + ir * 9 + ic] -= Mw[ic * 3 + ir];
}


static void CalcCBMx(const int Mi[9], const T_SMxI *SI)
{
  int        Order, Mul, ic, i;
  int        MiT_TMw[14 * 9], V[14], CumMx[9];
  int        Rank, nFreeP, FreePi[9];
  const int  *Mw, *PseudoGw;


  if ((Order = SI->Order) < 2 || SI->SenseOfRotation < 0) return;

      Mw = ReferenceRotMx(SI->Order, &PseudoGw);
  if (Mw == NULL) progerror("Corrupt ReferenceRotMx");

      Mul = MakeCumRMx(Mi, Order, CumMx);
  if (Mul != Order)
    progerror("Internal Error: CalcCBMx(): Mul != Order");

      Rank = iReducedRowEchelon(CumMx, 3, 3, NULL, 1);
  if (Rank != 1)
    progerror("Internal Error: CalcCBMx(): Rank != 1");

  SetupMiT_TMw(MiT_TMw, 14, Mi, Mw);

  for (i = 0; i < 3; i++) {
    MiT_TMw[ 9 * 9 + 0 + (i * 3)] = CumMx[i];
    MiT_TMw[10 * 9 + 1 + (i * 3)] = CumMx[i];
  }

  for (i = 0; i < 14; i++) V[i] = 0;

  for (i = 0; i < 3; i++) {
    MiT_TMw[(11 + i) * 9 + 2 + (i * 3)] = 1;
           V[11 + i] = SI->Basis[2][i];
  }

  MathMx(MiT_TMw, NULL, 0, 14, 9, 2, 1);
  putc('\n', stdout);
  MathMx(V, NULL, 0, 1, 14, 1, 0);
  putc('\n', stdout);

      Rank = iReducedRowEchelon(MiT_TMw, 14, 9, V, 1);
  if (Rank < 0) progerror(SgError);
  if (Rank > 9)
    progerror("Internal Error: CalcCBMx(): Rank > 9");

  Fprintf(stdout, "Rank = %d\n", Rank);

  MathMx(MiT_TMw, NULL, 0, 14, 9, 2, 1);
  putc('\n', stdout);
  MathMx(V, NULL, 0, 1, 14, 1, 0);
  putc('\n', stdout);

  nFreeP = iRRE_FreeParameters(MiT_TMw, 14, 9, FreePi);

  Fprintf(stdout, "Order = %2d nFreeP = %d (", Order, nFreeP);

  for (ic = 0; ic < nFreeP; ic++)
    Fprintf(stdout, " %d", FreePi[ic]);
  putc(')', stdout);
  putc('\n', stdout);

  if (Rank != 0)
    AllCBMx(MiT_TMw, 9, V, FreePi, nFreeP, Mi, Mw, PseudoGw, Order);
}


void SetSeitzMx(int R0, int R1, int R2,
                int R3, int R4, int R5,
                int R6, int R7, int R8,
                int T0, int T1, int T2,
                T_RTMx *SMx)
{
  SMx->s.R[0] = R0; SMx->s.R[1] = R1; SMx->s.R[2] = R2;
  SMx->s.R[3] = R3; SMx->s.R[4] = R4; SMx->s.R[5] = R5;
  SMx->s.R[6] = R6; SMx->s.R[7] = R7; SMx->s.R[8] = R8;
  SMx->s.T[0] = T0; SMx->s.T[1] = T1; SMx->s.T[2] = T2;
}


static void VerifySeitzMxInfo(const T_RTMx *SMx, const T_SMxI *SI)
{
  int wlD[3], wlT[3], i;


  for (i = 0; i < 3; i++)
    wlD[i] = (SMx->s.T[i] - SI->wg[i]) * (CTBF / STBF);

  RotMx_t_Vector(wlT, SMx->s.R, SI->Tr, 0);
  for (i = 0; i < 3; i++)
    wlT[i] = -wlT[i] + SI->Tr[i];

  for (i = 0; i < 3; i++)
    if (wlT[i] != wlD[i]) {
      Fprintf(stdout,
        "Verify Error: %2d T = %2d %2d %2d wlD = %3d %3d %3d, wlT = %3d %3d %3d\n",
        SI->Order,
        SMx->s.T[0], SMx->s.T[1], SMx->s.T[2],
        wlD[0], wlD[1], wlD[2],
        wlT[0], wlT[1], wlT[2]);
      break;
    }
}


static int EvalShowSeitzMx(const T_RTMx *SMx, int Show)
{
  T_SMxI  SI[1];


  ResetSeitzMxInfo(SI, 0);

  if (SetSeitzMxInfo(SMx, SI) == 0)
  {
    if (SgError != NULL) progerror(SgError);
    return 0;
  }

  if (Show)
  {
    Fprintf(stdout, "%2d ", SI->Order);
    MathMx(SI->Basis[2], NULL, 0, 1, 3, 2, 0);
    putc(' ', stdout);
    if (SI->SenseOfRotation < 0)
      putc('c', stdout);
    else
      putc(' ', stdout);
    putc(' ', stdout);
    MathMx(SMx->s.R, NULL, 0, 3, 3, 2, 0);
    Fprintf(stdout, " wg=");
    MathMx(SI->wg, NULL, STBF, 1, 3, 1, 0);
    Fprintf(stdout, " Tr=");
    MathMx(SI->Tr, NULL, CTBF, 1, 3, 1, 0);
    putc('\n', stdout);
  }

  VerifySeitzMxInfo(SMx, SI);

  return 1;
}


static void VerifyBasis(const int Mi[9], const T_SMxI *SI)
{
  int     Order, Mul, Rank, Det, i, j;
  int     CumMx[9], RRECumMx[9], V000[3];
  T_RTMx  SMx[1], BC_SMx[1], Tzp[1], Tpz[1];


  (void) MemCpy(SMx->s.R, Mi, 9);
  for (i = 0; i < 3; i++) SMx->s.T[i] = 0;

      Order = SI->Order;
  if (Order < 0) {
    for (i = 0; i < 9; i++) SMx->s.R[i] *= -1;
    Order *= -1;
  }

      Mul = MakeCumRMx(SMx->s.R, Order, CumMx);
  if (Mul != Order)
    progerror("Internal Error: VerifyBasis(): Mul != Order");

  if (GetRotMxInfo(SMx->s.R, NULL, NULL) != 0)
    putc('+', stdout);
  else
    putc('-', stdout);

  Fprintf(stdout, "CumMx=");
  MathMx(CumMx, NULL, 0, 3, 3, 1, 0);

  (void) MemCpy(RRECumMx, CumMx, 9);

  Rank = iReducedRowEchelon(RRECumMx, 3, 3, NULL, 1);

  Fprintf(stdout, "=>");
  MathMx(RRECumMx, NULL, 0, 3, 3, 1, 0);

  if (   (Order != 1 && Rank != 1)
      || (Order == 1 && Rank != 3)) {
    putc('\n', stdout);
    progerror("Internal Error: VerifyBasis(): Corrupt Rank");
  }

  Det = deterRotMx((int *) SI->Basis);

  Fprintf(stdout, " Det=%d ", Det);
  MathMx((int *) SI->Basis, NULL, 0, 3, 3, 1, 0);
  putc('\n', stdout);

  if (Det < 1)
    progerror("Internal Error: VerifyBasis(): Det < 1");

  if (Rank == 1)
  {
    for (i = 0; i < 2; i++)
    {
      RotMx_t_Vector(V000, RRECumMx, SI->Basis[i], 0);

      for (j = 0; j < 3; j++)
        if (V000[j] != 0) {
          MathMx(V000, NULL, 0, 1, 3, 1, 0);
          putc('\n', stdout);
          progerror("Internal Error: VerifyBasis(): V000[RRECumMx] != 000");
      }

      RotMx_t_Vector(V000,    CumMx, SI->Basis[i], 0);

      for (j = 0; j < 3; j++)
        if (V000[j] != 0) {
          MathMx(V000, NULL, 0, 1, 3, 1, 0);
          putc('\n', stdout);
          progerror("Internal Error: VerifyBasis(): V000[CumMx] != 000");
      }
    }
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      Tzp->s.R[i * 3 + j] = SI->Basis[j][i] * CRBF;
    }
    Tzp->s.T[i] = 0;
  }

  if (InverseRTMx(Tzp, Tpz, CRBF) == 0)
    progerror("Internal Error: VerifyBasis(): Tzp not invertible");

  Fputs("Tzp=", stdout);
  MathMx(Tzp->s.R, NULL, CRBF, 3, 3, 1, 0);
  putc('\n', stdout);
  Fputs("Tpz=", stdout);
  MathMx(Tpz->s.R, NULL, CRBF, 3, 3, 1, 0);
  putc('\n', stdout);

  if (CB_SMx(BC_SMx, Tpz, SMx, Tzp) != 0) {
    Fprintf(stdout, "%s\n", SgError);
    SgError = NULL;
    CalcCBMx(SMx->s.R, SI);
    progerror("Internal Error: VerifyBasis(): Failure");
  }

  SgError = NULL;
}


int BuildListValidRotMx(int *LVRMx, int mLVRMx,
                        int ProperOnly, int PositiveSenseOnly)
{
  int     R[9], nValid;
  T_SMxI  SI[1];


  nValid = 0;

#define loop(i) for (R[i] = -1; R[i] <= 1; R[i]++)
  loop(0) loop(1) loop(2)
  loop(3) loop(4) loop(5)
  loop(6) loop(7) loop(8)
#undef loop
  {
    if (   SetRotMxInfo(R, SI) != 0
        && (ProperOnly == 0 || SI->Order > 0)
        && (PositiveSenseOnly == 0 || SI->SenseOfRotation >= 0))
    {
      if (nValid == mLVRMx) return -1;

      (void) MemCpy(&LVRMx[nValid * 9], R, 9);
                           nValid++;
    }
  }

  return nValid;
}


/*
  SetSeitzMx( 0,  0, -1,
             -1,  0,  0,
              0,  1,  0,  0,  6,  0, SMx);
  (void) EvalShowSeitzMx(SMx, 1);

  SetSeitzMx( 0, -1,  0,
             -1,  0,  0,
              0,  0,  1,  6,  0,  9, SMx);
  (void) EvalShowSeitzMx(SMx, 1);
 */


static void usage(void)
{
  Fprintf(stderr,
    "usage: %s\n",
    progn);
  exit (1);
}


int main(int argc, char *argv[])
{
  int        i;
  int        mValid, nValid, iValid, jValid;
  int        *ListValidRotMx;
  const int  *Ri, *Rj;
  T_SMxI     SI[1];

  T_ListRMx  ListRMx[1];
  int        nGoodComb, nBadComb;


  for (i = 1; i < argc; i++)
  {
    usage();
  }

  mValid = 960;

      ListValidRotMx = malloc(mValid * 9 * sizeof (*ListValidRotMx));
  if (ListValidRotMx == NULL)
    progerror("Not enough core");

      nValid = BuildListValidRotMx(ListValidRotMx, mValid, 1, 1);
  if (nValid < 0)
    progerror("mValid too small");

  Fprintf(stdout, "nValid = %d\n", nValid);

  mValid = nValid;

      ListValidRotMx
        = realloc(ListValidRotMx, mValid * 9 * sizeof (*ListValidRotMx));
  if (ListValidRotMx == NULL)
    progerror("ListValidRotMx pointer corrupted");

  ListRMx->mList = 24;

  if (AllocateListRMx(ListRMx) == 0)
    progerror("Not enough core");

  nGoodComb = 0;
  nBadComb = 0;

  for (iValid = 0; iValid < nValid - 1; iValid++)
  {
    Ri = &ListValidRotMx[iValid * 9];

    for (jValid = iValid + 1; jValid < nValid; jValid++)
    {
      Rj = &ListValidRotMx[jValid * 9];

      ResetListRMx(ListRMx);

      if (Add2ListRMx(ListRMx, Ri) != 0)
        progerror(SgError);

      if (Add2ListRMx(ListRMx, Rj) != 0)
        nBadComb++;
      else
      {
        nGoodComb++;

        (void) SetRotMxInfo(Ri, SI);
        Fprintf(stdout, "%d ", SI->Order);
        (void) SetRotMxInfo(Rj, SI);
        Fprintf(stdout, "%d ", SI->Order);

        Fprintf(stdout, " %d\n", ListRMx->nList);
      }
    }
  }

  FreeListRMx(ListRMx);

  Fprintf(stdout, "nGoodComb = %d\n", nGoodComb);
  Fprintf(stdout, "nBadComb  = %d\n", nBadComb);

  return 0;
}
