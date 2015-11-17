#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "main.h"
#include "xtal.h"
#include "matrix.h"


void CalcLatticeTr2(T_LatticeConstants *LC, const T_LatticeInfo *LI,
                    Fprec *Min2, Fprec *Max2)
{
  int    i;
  Fprec  G[9], Min_u2, Max_u2, u2, *u, t[3];

  static Fprec List_u[] =
    {
      1., 0., 0.,
      0., 1., 0.,
      0., 0., 1.,
      0., 1., 1.,
      1., 0., 1.,
      1., 1., 0.,
      1., 1., 1.
    };


  Lc2MetricalMx(LC, G);

                u = List_u;
  Min_u2 =
  Max_u2 = DotG(u, G, u); u += 3;

  for (i = 1; i < sizeof List_u / (3 * sizeof (*List_u)); i++)
  {
    u2 = DotG(u, G, u); u += 3;
    if (Min_u2 > u2) Min_u2 = u2;
    if (Max_u2 < u2) Max_u2 = u2;
  }

  /* skip first TrVector which is always 0,0,0
   */
  for (i = 1; i < LI->nTrVector; i++)
  {
    t[0] = LI->TrVector[i * 3 + 0] / fSTBF;
    t[1] = LI->TrVector[i * 3 + 1] / fSTBF;
    t[2] = LI->TrVector[i * 3 + 2] / fSTBF;

    u2 = DotG(t, G, t);
    if (Min_u2 > u2) Min_u2 = u2;
    if (Max_u2 < u2) Max_u2 = u2;
  }

  if (Min2 != NULL) *Min2 = Min_u2;
  if (Max2 != NULL) *Max2 = Max_u2;
}


static void MarkOrbit(T_SgInfo *SgInfo,
                      T_PeakFlags *FlagField,
                      int n_x, int n_y, int n_z,
                      int i_x, int i_y, int i_z,
                      T_WyckoffList *WL_Entry,
                      int *List_e_iff, int *nList_e_iff)
{
  int           iff, nff, d_x, d_y, d_z, f_x, f_y, f_z;
  int           RTx, RTy, RTz, e_x, e_y, e_z, e_iff;
  int           iSymOp, iList;
  const int     *smxr, *smxt;
  const T_RTMx  *lsmx;
  int           nLoopInv, iLoopInv;
  int           iTrV, nTrV;
  const int     *TrV;
  int           nL_e_iff, i;


  iff = (i_x * n_y + i_y) * n_z + i_z;
  nff = n_x * n_y * n_z;
  d_x = 12       * n_y * n_z; f_x = d_x * i_x;
  d_y = 12 * n_x       * n_z; f_y = d_y * i_y;
  d_z = 12 * n_x * n_y      ; f_z = d_z * i_z;

  WL_Entry->nPositions = 0;

  nLoopInv = Sg_nLoopInv(SgInfo);
  nTrV = SgInfo->LatticeInfo->nTrVector;

  nL_e_iff = 0;
  iSymOp = 0;

  lsmx = SgInfo->ListSeitzMx;

  for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
  {
    smxr = lsmx->s.R;
    smxt = lsmx->s.T;

    RTx = smxr[0] * f_x + smxr[1] * f_y + smxr[2] * f_z + nff * smxt[0];
    RTy = smxr[3] * f_x + smxr[4] * f_y + smxr[5] * f_z + nff * smxt[1];
    RTz = smxr[6] * f_x + smxr[7] * f_y + smxr[8] * f_z + nff * smxt[2];

    for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
    {
      WL_Entry->WL_Flag[iSymOp] = WL_F_InActive;

      TrV = SgInfo->LatticeInfo->TrVector;

      for (iTrV = 0; iTrV < nTrV; iTrV++, TrV += 3)
      {
        if (iLoopInv == 0)
        {
          e_x =  RTx + nff * TrV[0];
          e_y =  RTy + nff * TrV[1];
          e_z =  RTz + nff * TrV[2];
        }
        else
        {
          e_x = -RTx + nff * TrV[0];
          e_y = -RTy + nff * TrV[1];
          e_z = -RTz + nff * TrV[2];
        }

        if (e_x % d_x)
          progerror("Inappropriate Grid x");
        if (e_y % d_y)
          progerror("Inappropriate Grid y");
        if (e_z % d_z)
          progerror("Inappropriate Grid z");

        e_x = iModPositive(e_x / d_x, n_x);
        e_y = iModPositive(e_y / d_y, n_y);
        e_z = iModPositive(e_z / d_z, n_z);

        e_iff = (e_x * n_y + e_y) * n_z + e_z;

        if (FlagField != NULL)
        {
          if (FlagField[e_iff] == -1)
          {
            if (iTrV == 0)
              WL_Entry->WL_Flag[iSymOp] = WL_F_Active;
            else
            {
              if (WL_Entry->WL_Flag[iSymOp] != WL_F_Active)
                InternalError("MarkOrbit(): Mismatch: WL_Flag != WL_F_Active");
            }
            FlagField[e_iff] = iff;
            WL_Entry->nPositions++;

            if (List_e_iff != NULL)
              List_e_iff[nL_e_iff++] = e_iff;
          }
          else if (FlagField[e_iff] == iff)
          {
            if (iTrV != 0 && WL_Entry->WL_Flag[iSymOp] == WL_F_Active)
              InternalError("MarkOrbit(): Mismatch: WL_Flag == WL_F_Active");

            if (e_iff == iff)
              WL_Entry->WL_Flag[iSymOp] = WL_F_Identity;
          }
          else
            InternalError("MarkOrbit(): Corrupt ListSeitzMx");
        }
        else if (List_e_iff != NULL)
        {
          for (i = 0; i < nL_e_iff; i++)
            if (List_e_iff[i] == e_iff) break;

          if (i == nL_e_iff)
          {
            if (iTrV == 0)
              WL_Entry->WL_Flag[iSymOp] = WL_F_Active;
            else
            {
              if (WL_Entry->WL_Flag[iSymOp] != WL_F_Active)
                InternalError("MarkOrbit(): Mismatch: WL_Flag != WL_F_Active");
            }

            List_e_iff[nL_e_iff++] = e_iff;
            WL_Entry->nPositions++;
          }
          else
          {
            if (iTrV != 0 && WL_Entry->WL_Flag[iSymOp] == WL_F_Active)
              InternalError("MarkOrbit(): Mismatch: WL_Flag == WL_F_Active");

            if (e_iff == iff)
              WL_Entry->WL_Flag[iSymOp] = WL_F_Identity;
          }
        }
      }

      iSymOp++;
    }
  }

  if (nList_e_iff != NULL) *nList_e_iff = nL_e_iff;
}


static void PrintWL_Flags(const int *WL_Flags)
{
  int  iSymOp;
  int  iList, nLoopInv, iLoopInv;


  nLoopInv = Sg_nLoopInv(&SpgrInfo);

  iSymOp = 0;

  for (iList = 0; iList < SpgrInfo.nList; iList++)
  {
    for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
    {
      switch (WL_Flags[iSymOp++])
      {
        case WL_F_Active:   putc('+', stdout); break;
        case WL_F_InActive: putc('-', stdout); break;
        case WL_F_Identity: putc('I', stdout); break;
        default:
          InternalError("Corrupt WL_Flags");
      }
    }
  }
}


void PrintWyckoffList(void)
{
  int            fld, iSymOp, CN;
  int            iWL, iList, nLoopInv, iLoopInv;
  long           l;
  Fprec          FreeV[3];
  T_WyckoffList  *eWL, **SubWL;


  Fprintf(stdout, ">Begin WyckoffList\n");

  nLoopInv = Sg_nLoopInv(&SpgrInfo);

  eWL = WyckoffList;

  for (iWL = 0; iWL < nWyckoffList; iWL++)
  {
    Fprintf(stdout, "WL %2d: %6d * %3d  ",
            iWL, eWL->Count, eWL->nPositions);

    PrintWL_Flags(eWL->WL_Flag);

    if (eWL->CanBeCN)
    {
      Fprintf(stdout, " => Coordination");
      Fprintf(stdout, " %c", (eWL->CanBeNode ? 'N' : '-'));

      for (CN = 0; CN <= NCNmax; CN++)
        if (eWL->CanBeCN[CN]) Fprintf(stdout, " %d", CN);
    }

    putc('\n', stdout);

    Fprintf(stdout, "                     Free %d", eWL->FreeVect.D);
    if (eWL->FreeVect.D == 1 || eWL->FreeVect.D == 2)
    {
      Tform_cx(eWL->FreeVect.V[0],
               eWL->FreeVect.V[1],
               eWL->FreeVect.V[2],
               FreeV[0],
               FreeV[1],
               FreeV[2]);
      Fprintf(stdout, " => %7.4f %7.4f %7.4f",
        FreeV[0], FreeV[1], FreeV[2]);
    }
    Fprintf(stdout, "\n        ");

    fld = 0;
    iSymOp = 0;

    for (iList = 0; iList < SpgrInfo.nList; iList++)
    {
      for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
      {
        if (eWL->WL_Flag[iSymOp++] == WL_F_Active)
        {
          if (fld++ == 8)
          {
            Fprintf(stdout, "\n        ");
            fld = 1;
          }

          putc(' ', stdout);
          if (iList + 1 < 10) putc(' ', stdout);
          if (iLoopInv == 0) putc(' ', stdout);
          putc('(', stdout);
          if (iLoopInv != 0) putc('-', stdout);
          Fprintf(stdout, "%d)", iList + 1);
        }
      }
    }
    putc('\n', stdout);

    if (eWL->SubWL != NULL && *eWL->SubWL != NULL)
    {
      fld = 0;

      for (SubWL = eWL->SubWL; *SubWL != NULL; SubWL++)
      {
        if (fld == 8)
        {
          putc('\n', stdout);
          fld = 0;
        }

        if (fld++ == 0)
          Fprintf(stdout, "        ");

        l = (*SubWL) - WyckoffList;
        if (l < 10) putc(' ', stdout);
        Fprintf(stdout, " [%ld]", l);
      }
      putc('\n', stdout);
    }

    eWL++;
  }

  Fprintf(stdout, ">End WyckoffList\n\n");
}


int SetToUnitLength(Fprec v[3])
{
  Fprec  l;


      l = DotC(v, v);
      l = AppSqrt(l);
  if (l == 0.)
    return 0;

  v[0] /= l;
  v[1] /= l;
  v[2] /= l;

  return 1;
}


int CombineFreeEigenVects(const T_FreeVect *PreFV,
                                T_FreeVect *SymFV,
                                T_FreeVect *NewFV,
                          int iList, int iLoopInv)
{
  int          Order, Orientation, PosCorrType;
  Fprec        FxN[3];
  Fprec        AbsP2, AbsS2, AbsP_dot_S;
  T_RotMxInfo  *rmxi;


  rmxi = &SpgrInfo.ListRotMxInfo[iList];

  if (iLoopInv == 0) Order =  rmxi->Order;
  else               Order = -rmxi->Order;

  if (Debug0) /* Debug: Print EV Order RefAxis DirCode */
  {
    Fprintf(stdout, "      [%2d %2d %2d]",
                     rmxi->EigenVector[0],
                     rmxi->EigenVector[1],
                     rmxi->EigenVector[2]);
    Fprintf(stdout, " %2d", Order);
    if (rmxi->Inverse != 0) Fprintf(stdout, " -1");
    else                    Fprintf(stdout, "   ");
    if (rmxi->RefAxis) Fprintf(stdout, " '%c'", rmxi->RefAxis);
    if (rmxi->DirCode) Fprintf(stdout, " '%c'", rmxi->DirCode);
    Fprintf(stdout, "  (%d)\n", iList + 1);
  }

  Tform_xc(rmxi->EigenVector[0],
           rmxi->EigenVector[1],
           rmxi->EigenVector[2],
           SymFV->V[0],
           SymFV->V[1],
           SymFV->V[2]);

  if (Order < 0)
  {
    if (Order == -2)
      SymFV->D = 2;
    else
      SymFV->D = 0;
  }
  else
    SymFV->D = 1;

  if (SymFV->D && SetToUnitLength(SymFV->V) == 0)
    InternalError("CombineFreeEigenVects(): Corrupt SymFV->V");

  NewFV->D    = -1;
  PosCorrType = PCT_Unknown;

  if (SymFV->D == 0)
  {
    NewFV->D    = 0;
    PosCorrType = PCT_None;
  }
  else if (PreFV->D == 3)
  {
    (void) memcpy(NewFV, SymFV, sizeof (*SymFV));
    PosCorrType = PCT_None;
  }
  else if (PreFV->D == 1 || PreFV->D == 2)
  {
    AbsP2 = DotC(PreFV->V, PreFV->V);
    AbsS2 = DotC(SymFV->V, SymFV->V);
    AbsP_dot_S = AppSqrt(AbsP2 * AbsS2);
    if (AbsP_dot_S == 0.)
      InternalError("CombineFreeEigenVects(): Corrupt EigenVectors");
    AbsP_dot_S = DotC(PreFV->V, SymFV->V) / AbsP_dot_S;
    AbsP_dot_S = AppFabs(AbsP_dot_S);

    Orientation = EvalAbsDot(AbsP_dot_S);

    if (PreFV->D == 2)
    {
      if (Order == -2) /* Improper + Improper */
      {
        if (Orientation !=  1) /* not parallel */
        {
          CrossC(PreFV->V, SymFV->V, FxN);
          NewFV->V[0] = FxN[0];
          NewFV->V[1] = FxN[1];
          NewFV->V[2] = FxN[2];
          NewFV->D    = 1;
          PosCorrType = PCT_II;
        }
        else                       /* parallel */
        {
          NewFV->V[0] = PreFV->V[0];
          NewFV->V[1] = PreFV->V[1];
          NewFV->V[2] = PreFV->V[2];
          NewFV->D    = 2;
          PosCorrType = PCT_NoChange;
        }
      }
      else /* Improper + Proper */
      {
        if (Orientation != -1) /* not perpendicular */
        {
          NewFV->D    = 0;
          PosCorrType = PCT_IP;
        }
        else                       /* perpendicular */
        {
          NewFV->V[0] = SymFV->V[0];
          NewFV->V[1] = SymFV->V[1];
          NewFV->V[2] = SymFV->V[2];
          NewFV->D    = 1;
          PosCorrType = PCT_None;
        }
      }
    }
    else /* if (PreFV->D == 1) */
    {
      if (Order == -2) /* Proper + Improper */
      {
        if (Orientation != -1) /* not perpendicular */
        {
          NewFV->D    = 0;
          PosCorrType = PCT_PI;
        }
        else                       /* perpendicular */
        {
          NewFV->V[0] = PreFV->V[0];
          NewFV->V[1] = PreFV->V[1];
          NewFV->V[2] = PreFV->V[2];
          NewFV->D    = 1;
          PosCorrType = PCT_None;
        }
      }
      else /* Proper + Proper */
      {
        if (Orientation !=  1) /* not parallel */
        {
          NewFV->D    = 0;
          PosCorrType = PCT_PP;
        }
        else                       /* parallel */
        {
          NewFV->V[0] = PreFV->V[0];
          NewFV->V[1] = PreFV->V[1];
          NewFV->V[2] = PreFV->V[2];
          NewFV->D    = 1;
          PosCorrType = PCT_NoChange;
        }
      }
    }
  }
  else /* if (PreFV->D != 1 && PreFV->D != 2 && PreFV->D != 3) */
    InternalError("CombineFreeEigenVects(): Corrupt PreFV->D");

  if      (NewFV->D == 0)
  {
    NewFV->V[0] =
    NewFV->V[1] =
    NewFV->V[2] = 0.;
  }
  else if (NewFV->D == 1 || NewFV->D == 2)
  {
    if (SetToUnitLength(NewFV->V) == 0)
      InternalError("CombineFreeEigenVects(): Corrupt NewFV->V");
  }
  else if (NewFV->D != 3)
    InternalError("CombineFreeEigenVects(): Corrupt NewFV->D");

  if (PosCorrType == PCT_Unknown)
    InternalError("CombineFreeEigenVects(): Corrupt PosCorrType");

  return PosCorrType;
}


static void SetFreeVect(T_WyckoffList *eWL)
{
  int         iSymOp;
  int         iList, nLoopInv, iLoopInv;
  int         n;
  T_FreeVect  FV[1], SymFV[1], NewFV[1];


  nLoopInv = Sg_nLoopInv(&SpgrInfo);

  FV->V[0] =
  FV->V[1] =
  FV->V[2] = 0.;
  FV->D    = 3;

  iSymOp = 0;

  for (iList = 0; iList < SpgrInfo.nList; iList++)
  {
    for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
    {
      if (eWL->WL_Flag[iSymOp++] == WL_F_Identity)
      {
        (void) CombineFreeEigenVects(FV, SymFV, NewFV, iList, iLoopInv);

        (void) memcpy(FV, NewFV, sizeof (*FV));

        if (FV->D == 0)
          goto NiceFV;
      }
    }
  }

  NiceFV:
                     n = 0;
  if (FV->V[0] < 0.) n++;
  if (FV->V[1] < 0.) n++;
  if (FV->V[2] < 0.) n++;

  if (n > 1)
  {
    FV->V[0] = -FV->V[0];
    FV->V[1] = -FV->V[1];
    FV->V[2] = -FV->V[2];
  }

  (void) memcpy(&eWL->FreeVect, FV, sizeof eWL->FreeVect);

  if (Debug0) /* Debug: Print FV->D FV->V */
  {
    Fprintf(stdout, "    Free %d", eWL->FreeVect.D);
    if (eWL->FreeVect.D == 1 || eWL->FreeVect.D == 2)
    {
      Tform_cx(eWL->FreeVect.V[0],
               eWL->FreeVect.V[1],
               eWL->FreeVect.V[2],
               FV->V[0],
               FV->V[1],
               FV->V[2]);
      Fprintf(stdout, " => %7.4f %7.4f %7.4f",
        FV->V[0], FV->V[1], FV->V[2]);
    }
    putc('\n', stdout);
  }
}


static void AddSubWL(void)
{
  int            iWL, jWL, nBufSubWL, nSymOp, iSymOp, i;
  T_WyckoffList  *i_eWL, *j_eWL;
  T_WyckoffList  *BufSubWL[MaxWyckoffList];


  nSymOp = SpgrInfo.OrderP;

  i_eWL = WyckoffList;

  for (iWL = 0; iWL < nWyckoffList; iWL++)
  {
    nBufSubWL = 0;

    if (i_eWL->FreeVect.D != 0)
    {
             jWL = iWL;
      while (jWL-- > 0)
      {
        j_eWL = &WyckoffList[jWL];

        for (iSymOp = 0; iSymOp < nSymOp; iSymOp++)
          if (   i_eWL->WL_Flag[iSymOp] != WL_F_Active
              && i_eWL->WL_Flag[iSymOp] != j_eWL->WL_Flag[iSymOp]
              && j_eWL->WL_Flag[iSymOp] != WL_F_Identity) break;

        if (iSymOp == nSymOp)
        {
          if (nBufSubWL >= MaxWyckoffList - 1)
            InternalError("AddSubWL()");

          BufSubWL[nBufSubWL++] = j_eWL;
        }
      }
    }

    BufSubWL[nBufSubWL++] = NULL;

    CheckMalloc(i_eWL->SubWL, nBufSubWL);

    for (i = 0; i < nBufSubWL; i++)
      i_eWL->SubWL[i] = BufSubWL[i];

    i_eWL++;
  }
}


static int WL_SortFunction(const T_WyckoffList *a, const T_WyckoffList *b)
{
  int  i, nSymOp;


  if (a->nPositions < b->nPositions) return -1;
  if (a->nPositions > b->nPositions) return  1;

  if (a->FreeVect.D < b->FreeVect.D) return -1;
  if (a->FreeVect.D > b->FreeVect.D) return  1;

  nSymOp = SpgrInfo.OrderP;

  for (i = 0; i < nSymOp; i++)
  {
    if (abs(a->WL_Flag[i]) > abs(b->WL_Flag[i])) return -1;
    if (abs(a->WL_Flag[i]) < abs(b->WL_Flag[i])) return  1;
  }

  for (i = 0; i < nSymOp; i++)
  {
    if (a->WL_Flag[i] > b->WL_Flag[i]) return -1;
    if (a->WL_Flag[i] < b->WL_Flag[i]) return  1;
  }

  return 0;
}


void BuildWyckoffList(T_SgInfo *SgInfo)
{
  int            Max_nPositions, nSymOp, i;
  int            FFDx, FFDy, FFDz, i_x, i_y, i_z, e_x, e_y, e_z;
  int            nff, *ff, *FlagField;
  int            *List_e_iff, nList_e_iff, iList_e_iff;
  int            *BufList_e_iff, nBufList_e_iff;
  int            iWL;
  T_WyckoffList  WL[MaxWyckoffList], Buf_eWL;


  Max_nPositions = SgInfo->OrderL;
  nSymOp         = SgInfo->OrderP;

  if      (   SgInfo->XtalSystem == XS_Trigonal
           || SgInfo->XtalSystem == XS_Hexagonal)
    FFDx = FFDy = FFDz = 12;
  else if (   SgInfo->XtalSystem == XS_Cubic
           && SgInfo->LatticeInfo->nTrVector != 1)
    FFDx = FFDy = FFDz = 16;
  else
    FFDx = FFDy = FFDz =  8;

  nff = FFDx * FFDy * FFDz;

  CheckMalloc(BufList_e_iff, Max_nPositions);
  CheckMalloc(List_e_iff, Max_nPositions);
  CheckMalloc(FlagField, nff);

  nWyckoffList = 0;
  Buf_eWL.WL_Flag = NULL;

  ff = FlagField;
         i = nff;
  while (i--) *ff++ = -1;

  ff = FlagField;

  for (i_x = 0; i_x < FFDx; i_x++)
  for (i_y = 0; i_y < FFDy; i_y++)
  for (i_z = 0; i_z < FFDz; i_z++)
  {
    if (*ff == -1)
    {
      if (Buf_eWL.WL_Flag == NULL)
        CheckMalloc(Buf_eWL.WL_Flag, nSymOp);

      MarkOrbit(SgInfo, FlagField, FFDx, FFDy, FFDz, i_x, i_y, i_z,
                &Buf_eWL,
                List_e_iff, &nList_e_iff);

      if (Debug0) /* Debug: Print MO WL_Image */
      {
        Buf_eWL.Count = 0;
        Fprintf(stdout, "MO %2d %2d %2d  =>  ", i_x, i_y, i_z);
        PrintWL_Flags(Buf_eWL.WL_Flag);
        putc('\n', stdout);
      }

      for (iList_e_iff = 0; iList_e_iff < nList_e_iff; iList_e_iff++)
      {
        if (iList_e_iff != 0)
        {
          if (Buf_eWL.WL_Flag == NULL)
            CheckMalloc(Buf_eWL.WL_Flag, nSymOp);

                e_x = List_e_iff[iList_e_iff];
          e_z = e_x %  FFDz;
                e_x /= FFDz;
          e_y = e_x %  FFDy;
                e_x /= FFDy;

          MarkOrbit(SgInfo, NULL, FFDx, FFDy, FFDz, e_x, e_y, e_z,
                    &Buf_eWL,
                    BufList_e_iff, &nBufList_e_iff);

          if (nBufList_e_iff != nList_e_iff)
            InternalError("BuildWyckoffList(): nBufList_e_iff != nList_e_iff");
        }

        for (iWL = 0; iWL < nWyckoffList; iWL++)
        {
          if (WL[iWL].nPositions == Buf_eWL.nPositions)
          {
            for (i = 0; i < nSymOp; i++)
              if (WL[iWL].WL_Flag[i] != Buf_eWL.WL_Flag[i]) break;

            if (i == nSymOp)
              break;
          }
        }

        if (iWL == nWyckoffList)
        {
          if (nWyckoffList >= MaxWyckoffList)
            InternalError("BuildWyckoffList(): MaxWyckoffList too small");

          WL[nWyckoffList].nPositions = Buf_eWL.nPositions;
          WL[nWyckoffList].WL_Flag    = Buf_eWL.WL_Flag;
          WL[nWyckoffList].CheckCount = 0;
             nWyckoffList++;

          Buf_eWL.WL_Flag = NULL;
        }

        if (iList_e_iff == 0)
        {
          WL[iWL].CheckCount++;

          if (Buf_eWL.nPositions == Max_nPositions)
            break;
        }
      }
    }

    ff++;
  }

  if (Buf_eWL.WL_Flag != NULL)
    AppFree(Buf_eWL.WL_Flag, nSymOp);

  AppFree(FlagField, nff);

  AppFree(List_e_iff, Max_nPositions);
  AppFree(BufList_e_iff, Max_nPositions);

  if (nWyckoffList == 0)
    InternalError("BuildWyckoffList(): nWyckoffList == 0");

  CheckMalloc(WyckoffList, nWyckoffList);

  for (iWL = 0; iWL < nWyckoffList; iWL++)
  {
    WyckoffList[iWL].nPositions  = WL[iWL].nPositions;
    WyckoffList[iWL].WL_Flag     = WL[iWL].WL_Flag;
    SetFreeVect(&WyckoffList[iWL]);
    WyckoffList[iWL].SubWL       = NULL;
    WyckoffList[iWL].CheckCount  = WL[iWL].CheckCount;
    WyckoffList[iWL].Count       = 0;
    WyckoffList[iWL].CanBeNode   = 0;
    WyckoffList[iWL].CanBeCN     = NULL;
  }

  qsort((void *) WyckoffList, nWyckoffList, sizeof (*WyckoffList),
        (SortFunction) WL_SortFunction);

  AddSubWL();

  if (Debug0) /* Debug: */ PrintWyckoffList();

                  iWL = nWyckoffList - 1;
  if (WyckoffList[iWL].nPositions != Max_nPositions)
    InternalError("BuildWyckoffList(): No general position");
}


int FindWL_Entry(T_WyckoffList *WL, T_WyckoffList *Buf_eWL)
{
  int            iWL, nSymOp, iSymOp;


  nSymOp = SpgrInfo.OrderP;

  for (iWL = 0; iWL < nWyckoffList; iWL++)
  {
    if (WL->nPositions == Buf_eWL->nPositions)
    {
      for (iSymOp = 0; iSymOp < nSymOp; iSymOp++)
        if (WL->WL_Flag[iSymOp] != Buf_eWL->WL_Flag[iSymOp]) break;

      if (iSymOp == nSymOp)
      {
        WL->Count++;
        return iWL;
      }
    }

    WL++;
  }

  return -1;
}


void DS_MarkEquiv(T_SgInfo *SgInfo,
                  T_PeakFlags *FlagField, int n_x, int n_y, int n_z)
{
  int            iWL, nSymOp, nAsyGridPoints;
  int            *ff, i_x, i_y, i_z, i;
  T_WyckoffList  Buf_eWL;


  for (iWL = 0; iWL < nWyckoffList; iWL++)
    WyckoffList[iWL].Count = 0;

  nSymOp = SgInfo->OrderP;
  CheckMalloc(Buf_eWL.WL_Flag, nSymOp);

  ff = FlagField;
         i = n_x * n_y * n_z;
  while (i--) *ff++ = -1;

  nAsyGridPoints = 0;

  ff = FlagField;

  for (i_x = 0; i_x < n_x; i_x++)
  for (i_y = 0; i_y < n_y; i_y++)
  for (i_z = 0; i_z < n_z; i_z++)
  {
    if (*ff == -1)
    {
      MarkOrbit(SgInfo, FlagField, n_x, n_y, n_z, i_x, i_y, i_z,
                &Buf_eWL,
                NULL, NULL);

          *ff = FindWL_Entry(WyckoffList, &Buf_eWL);
      if (*ff < 0)
        InternalError("FindWL_Entry() < 0");

      *ff |= PF_HighBit;

      nAsyGridPoints++;
    }

    ff++;
  }

  AppFree(Buf_eWL.WL_Flag, nSymOp);

  if (Debug1) /* Debug: Print nAsyGridPoints */
    Fprintf(stdout, "# Grid-points per asymmetric unit = %d\n\n",
      nAsyGridPoints);

  if (Debug1) /* Debug: */ PrintWyckoffList();

  for (i = 0; i < nWyckoffList; i++)
    if (WyckoffList[i].Count == 0 && WyckoffList[i].CheckCount != 0)
      progerror("Inappropriate Grid_xyz");

  if (Debug0) /* Debug: Print FlagField */
  {
    for (i_z = 0; i_z < n_z; i_z++)
    { for (i_y = 0; i_y < n_y; i_y++)
      { for (i_x = 0; i_x < n_x; i_x++)
        {
          i = (i_x * n_y + i_y) * n_z + i_z;

          if (FlagField[i] & PF_HighBit)
            Fprintf(stdout, " %3d,%d", i, (int)(FlagField[i] & PF_iWL_Mask));
          else
            Fprintf(stdout, " %3d  ", (int) FlagField[i]);
        }
        putc('\n', stdout);
      }
      putc('\n', stdout);
    }
    putc('\n', stdout);
  }
}


static int CanBeCN(const T_SgInfo *SgInfo,
                   const int *WL_Flag, const int *NoSPO)
{
  int  iSymOp, Order;


  /* skip first SymOp which is always the identity (and WL_F_Active)
   */
  WL_Flag++;

  for (iSymOp = 1; iSymOp < SgInfo->OrderP; iSymOp++)
  {
    if (*WL_Flag++ == WL_F_Identity)
    {
      if (SgInfo->Centric == -1)
      {
        Order = SgInfo->ListRotMxInfo[iSymOp / 2].Order;
        if (iSymOp % 2 != 0) Order = -Order;
      }
      else
        Order = SgInfo->ListRotMxInfo[iSymOp    ].Order;

      if (NoSPO[Order + 6]) return 0;
    }
  }

  return 1;
}


void SetCanBeCN(T_SgInfo *SgInfo)
{
  int            IntersectNoSPO[13], Order;
  int            iWL;
  T_WyckoffList  *WL;
  int            iNT, CN;
  T_NodeType     *NT;


  for (Order = -6; Order <= 6; Order++)
    IntersectNoSPO[Order + 6] = 1;

  NT = NodeTypes;

  for (iNT = 0; iNT < nNodeTypes; iNT++, NT++)
    for (Order = -6; Order <= 6; Order++)
      if (NT->NoSPO[Order + 6] == 0)
        IntersectNoSPO[Order + 6] = 0;

  WL = WyckoffList;

  for (iWL = 0; iWL < nWyckoffList; iWL++, WL++)
  {
    WL->CanBeNode = CanBeCN(SgInfo, WL->WL_Flag, IntersectNoSPO);

    CheckMalloc(WL->CanBeCN, NCNmax + 1);

    for (CN = 0; CN <= NCNmax; CN++) WL->CanBeCN[CN] = 0;

    NT = NodeTypes;

    for (iNT = 0; iNT < nNodeTypes; iNT++, NT++)
      WL->CanBeCN[NT->CN] = CanBeCN(SgInfo, WL->WL_Flag, NT->NoSPO);
  }
}


static int MinDist2Shift(T_fVector *Pos1, T_fVector *Pos2,
                         Fprec *MinDist2, T_iVector *m0pShift)
{
  int    Update;
  int    shift_x, shift_y, shift_z;
  Fprec  dx0, dy0, dz0, dx, dy, dz, dxc, dyc, dzc, Dist2;

  static Fprec  fm0p[] = {   -1., 0.,   1. };
  static int    im0p[] = { -STBF, 0,  STBF };


  Update = 0;

  dx0 = Pos2->x - Pos1->x;
  dy0 = Pos2->y - Pos1->y;
  dz0 = Pos2->z - Pos1->z;

  for (shift_x = 0; shift_x < 3; shift_x++)
  {
        dx = dx0 + fm0p[shift_x];

    for (shift_y = 0; shift_y < 3; shift_y++)
    {
        dy = dy0 + fm0p[shift_y];

      for (shift_z = 0; shift_z < 3; shift_z++)
      {
        dz = dz0 + fm0p[shift_z];

        Tform_xc(dx, dy, dz, dxc, dyc, dzc);
        Dist2 = dxc * dxc + dyc * dyc + dzc * dzc;

        if (*MinDist2 > Dist2)
        {
          *MinDist2 = Dist2;

          if (m0pShift != NULL)
          {
            m0pShift->x = im0p[shift_x];
            m0pShift->y = im0p[shift_y];
            m0pShift->z = im0p[shift_z];
          }

          Update = 1;
        }
      }
    }
  }

  return Update;
}


int CalcSymEquiv(Fprec xx, Fprec yx, Fprec zx,
                 T_fVector *SE, int MaxSE,
                 Fprec     Dist2ConsiderSame,
                 Fprec *MaxDist2ConsideredSame,
                 Fprec *MinDist2Distinct,
                 int *Buf_WL_Flag, int MaxBuf_WL_Flag)
{
  int              nSE, iSE, iList, nTrV, iTrV;
  int              nLoopInv, iLoopInv;
  Fprec            RTx, RTy, RTz;
  const T_fRTMx    *lfsmx;
  const T_fVector  *fTrV;
  const Fprec      *fsmxr, *fsmxt;
  T_fVector        Equiv, *eSE;
  int              WL_Flag, iSymOp;
  Fprec            MinDist2;


  nLoopInv = Sg_nLoopInv(&SpgrInfo);
  nTrV = SpgrInfo.LatticeInfo->nTrVector;

  nSE = 0;
  eSE = SE;

  iSymOp = 0;

  lfsmx = List_fSeitzMx;

  for (iList = 0; iList < SpgrInfo.nList; iList++, lfsmx++)
  {
    fsmxr = lfsmx->s.R;
    fsmxt = lfsmx->s.T;

    RTx = fsmxr[0] * xx + fsmxr[1] * yx + fsmxr[2] * zx + fsmxt[0];
    RTy = fsmxr[3] * xx + fsmxr[4] * yx + fsmxr[5] * zx + fsmxt[1];
    RTz = fsmxr[6] * xx + fsmxr[7] * yx + fsmxr[8] * zx + fsmxt[2];

    for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
    {
      WL_Flag = WL_F_InActive;

      for (iTrV = 0, fTrV = fTrVector; iTrV < nTrV; iTrV++, fTrV++)
      {
        if (iLoopInv == 0)
        {
          Equiv.x =  RTx + fTrV->x;
          Equiv.y =  RTy + fTrV->y;
          Equiv.z =  RTz + fTrV->z;
        }
        else
        {
          Equiv.x = -RTx + fTrV->x;
          Equiv.y = -RTy + fTrV->y;
          Equiv.z = -RTz + fTrV->z;
        }

        NormOf(Equiv.x);
        NormOf(Equiv.y);
        NormOf(Equiv.z);

        if (nSE == 0) MinDist2 = MinLatticeTr2;
        else          MinDist2 = MaxLatticeTr2;

                 eSE = SE;
                 iSE = nSE;
        if      (iSE--)
        {
          (void) MinDist2Shift(eSE++, &Equiv, &MinDist2, NULL);

          if (MinDist2 <= Dist2ConsiderSame && WL_Flag == WL_F_InActive)
            WL_Flag = WL_F_Identity;

          while (iSE--)
            (void) MinDist2Shift(eSE++, &Equiv, &MinDist2, NULL);
        }

        if (MinDist2 > Dist2ConsiderSame)
        {
          if (nSE >= MaxSE)
            InternalError("CalcSymEquiv(): nSE >= MaxSE");

          if (*MinDist2Distinct > MinDist2 && nSE != 0)
              *MinDist2Distinct = MinDist2;

          eSE->x = Equiv.x;
          eSE->y = Equiv.y;
          eSE->z = Equiv.z;
          eSE++;
          nSE++;

          if (iTrV == 0)
            WL_Flag = WL_F_Active;
          else if (WL_Flag != WL_F_Active)
            InternalError("CalcSymEquiv(): WL_Flag != WL_F_Active");
        }
        else
        {
          if (*MaxDist2ConsideredSame < MinDist2)
              *MaxDist2ConsideredSame = MinDist2;

          if (iTrV != 0 && WL_Flag == WL_F_Active)
            InternalError("CalcSymEquiv(): WL_Flag == WL_F_Active");
        }
      }

      if (Buf_WL_Flag != NULL)
      {
        if (iSymOp >= MaxBuf_WL_Flag)
          InternalError("CalcSymEquiv(): iSymOp >= MaxBuf_WL_Flag");

        Buf_WL_Flag[iSymOp++] = WL_Flag;
      }
    }
  }

  if (nSE == 1 && *MinDist2Distinct > MinLatticeTr2)
                  *MinDist2Distinct = MinLatticeTr2;

  return nSE;
}


void CheckWL_Entry(T_eD_PeakList *EPL)
{
  static T_fVector  *SymEquiv = NULL;
  static int        MaxSymEquiv = 0;
  static int        *Buf_WL_Flag = NULL;
  static int        nSymOp = 0;
  int               nSE, iSymOp;
  Fprec                Dist2ConsiderSame;
  Fprec             MaxDist2ConsideredSame;
  Fprec             MinDist2Distinct;


  if (SymEquiv == NULL)
  {
    MaxSymEquiv = SpgrInfo.OrderL;
    CheckMalloc(SymEquiv, MaxSymEquiv);

    nSymOp = SpgrInfo.OrderP;
    CheckMalloc(Buf_WL_Flag, nSymOp);
  }

  if (Debug0) /* Debug: Print Check WL_Entry */
    Fprintf(stdout, "Check WL_Entry %2d  for  %7.4f %7.4f %7.4f\n",
      (int)(EPL->WL_Entry - WyckoffList),
      EPL->Position.x, EPL->Position.y, EPL->Position.z);

  Dist2ConsiderSame = .001 * .001; /* ARBITRARY */
  MaxDist2ConsideredSame = -1.;
  MinDist2Distinct = MaxLatticeTr2;

  nSE = CalcSymEquiv(EPL->Position.x, EPL->Position.y, EPL->Position.z,
                     SymEquiv, MaxSymEquiv,
                         Dist2ConsiderSame,
                     &MaxDist2ConsideredSame,
                     &MinDist2Distinct,
                     Buf_WL_Flag, nSymOp);

  if (nSE != EPL->WL_Entry->nPositions)
    InternalError("CheckWL_Entry(): nSE != nPositions");

  for (iSymOp = 0; iSymOp < nSymOp; iSymOp++)
    if (EPL->WL_Entry->WL_Flag[iSymOp] != Buf_WL_Flag[iSymOp])
      InternalError("CheckWL_Entry(): WL_Flag mismatch");

  if (MinDist2Distinct <= CatchDistance2 * (Fprec)(1. - 1.e-5)) /* ARBITRARY */
    InternalError("CheckWL_Entry(): MinDist2Distinct < CatchDistance2");
}


static void FillList_fSeitzMx(const T_SgInfo *SgInfo, T_fRTMx *List_fSMx)
{
  int           iList, i;
  const T_RTMx  *lsmx;


  lsmx = SgInfo->ListSeitzMx;

  for (iList = 0; iList < SgInfo->nList; iList++, lsmx++, List_fSMx++)
  {
    for (i = 0; i < 9; i++) List_fSMx->s.R[i] = lsmx->s.R[i];
    for (i = 0; i < 3; i++) List_fSMx->s.T[i] = lsmx->s.T[i] / fSTBF;
  }
}


void fSymOps(T_SgInfo *SgInfo, T_fRTMx *List_fSMx, T_fVector *fTrV)
{
  int        nTrV;
  const int  *TrV;


  FillList_fSeitzMx(SgInfo, List_fSMx);

  nTrV = SgInfo->LatticeInfo->nTrVector;
   TrV = SgInfo->LatticeInfo->TrVector;

  while (nTrV--)
  {
    fTrV->x = (*TrV++) / fSTBF;
    fTrV->y = (*TrV++) / fSTBF;
    fTrV->z = (*TrV++) / fSTBF;
    fTrV++;
  }
}


void cSymOps(T_SgInfo *SgInfo, T_fRTMx *List_cSMx, T_fVector *cTrV,
             Fprec *Txc, Fprec *Tcx)
{
  int        iList, nTrV;
  const int  *TrV;
  Fprec      BufMx[9];


  FillList_fSeitzMx(SgInfo, List_cSMx);

  /* Definition of Txc: [X]c = Txc . [X]x
       where [X]x is        a vector in terms of the crystallographic basis
         and [X]c is the same vector in terms of the cartesian        basis
     Definition of Tcx: Tcx = Txc^-1, i.e. Tcx is the inverse of Txc

     Transformation of Symmetry Matrices: [A]c = Txc . [A]x . Tcx
   */

  for (iList = 0; iList < SgInfo->nList; iList++, List_cSMx++)
  {
    /* rotation part */
    MxMultiply(BufMx, List_cSMx->s.R, Tcx, 3, 3, 3);
    MxMultiply(List_cSMx->s.R, Txc, BufMx, 3, 3, 3);

    /* translation part */
    MxMultiply(BufMx, Txc, List_cSMx->s.T, 3, 3, 1);
    List_cSMx->s.T[0] = BufMx[0];
    List_cSMx->s.T[1] = BufMx[1];
    List_cSMx->s.T[2] = BufMx[2];
  }

  nTrV = SgInfo->LatticeInfo->nTrVector;
   TrV = SgInfo->LatticeInfo->TrVector;

  while (nTrV--)
  {
    cTrV->x = (Txc[0] * TrV[0] + Txc[1] * TrV[1] + Txc[2] * TrV[2]) / fSTBF;
    cTrV->y = (Txc[3] * TrV[0] + Txc[4] * TrV[1] + Txc[5] * TrV[2]) / fSTBF;
    cTrV->z = (Txc[6] * TrV[0] + Txc[7] * TrV[1] + Txc[8] * TrV[2]) / fSTBF;
    cTrV++;
    TrV += 3;
  }
}
