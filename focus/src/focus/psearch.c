#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "main.h"
#include "xtal.h"
#include "matrix.h"
#include "io.h"
#include "lib.h"


void MapInd(int nMap, Fprec *DensMap, const T_PeakFlags *FlagMap,
            int ChkInDep, int MaxInDep, Fprec *CorrCoeff)
{
  int     iMap, nDep;
  Fprec   minx, maxx, miny, maxy;
  double  x, y, Sx, Sx2, Sy, Sy2, Sxy, LRb, LRm, LRc;

  Fprec              *Dens;
  const T_PeakFlags  *FMap;


  if (CorrCoeff) *CorrCoeff = 0.;

  nDep = 0;

  if (ChkInDep)
  {
    minx = 0.;
    maxx = 0.;
    miny = 0.;
    maxy = 0.;
    Sx  = 0.;
    Sx2 = 0.;
    Sy  = 0.;
    Sy2 = 0.;
    Sxy = 0.;

    Dens = DensMap;
    FMap = FlagMap;

    for (iMap = 0; iMap < nMap; iMap++, Dens++, FMap++)
    {
      if (! (*FMap & PF_HighBit)) {
        x = *Dens;
        y =  DensMap[*FMap];
        if (nDep == 0) {
          minx = x;
          maxx = x;
          miny = y;
          maxy = y;
        }
        else {
          if (minx > x) minx = x;
          if (maxx < x) maxx = x;
          if (miny > y) miny = y;
          if (maxy < y) maxy = y;
        }
        Sx  = Sx  + x;
        Sx2 = Sx2 + x * x;
        Sy  = Sy  + y;
        Sy2 = Sy2 + y * y;
        Sxy = Sxy + x * y;
        nDep++;
      }
    }

    if (nDep > 0)
    {
      Sum2bmc(nDep,
              (double) minx, (double) maxx,
              (double) miny, (double) maxy,
              Sx, Sx2, Sy, Sy2, Sxy,
              &LRb, &LRm, &LRc);

      if (   LRm == 0.
          && LRc == 0.
          && maxy - miny < 1.e-6 * (AppFabs(miny) > AppFabs(maxy) ?
                                    AppFabs(miny) : AppFabs(maxy)))
        LRc = 1.;

      if (CorrCoeff)
         *CorrCoeff = LRc;
      else if (LRc < 0.9999) {
        Fprintf(stdout,
         "\n# Correlation of dependent and independent grid points = %.6g\n\n",
          LRc);
        InternalError(
          "Poor correlation of dependent and independent grid points");
      }
    }
  }

  if (MaxInDep && (! ChkInDep || nDep > 0))
  {
    Dens = DensMap;
    FMap = FlagMap;

    for (iMap = 0; iMap < nMap; iMap++, Dens++, FMap++)
      if (! (*FMap & PF_HighBit))
        if (DensMap[*FMap] < *Dens)
            DensMap[*FMap] = *Dens;
  }

  Dens = DensMap;
  FMap = FlagMap;

  for (iMap = 0; iMap < nMap; iMap++, Dens++, FMap++)
    if (! (*FMap & PF_HighBit))
      *Dens = DensMap[*FMap];
}


static Fprec PeakSearch(Fprec *eDensity, T_PeakFlags *PeakFlags, int Level)
{
  int          nk;
  int  im, jm, km;
  int  i0, j0, k0;
  int  ip, jp, kp;
  int  ibreak, jbreak, kbreak;
  int  nj_nk, ni_nj_nk;

  Fprec        *eD;
  int          iPF;
  T_PeakFlags  *PF, *AsyPF;

#define ResolutionRawPeakHistogram 128
  static int  nBox = ResolutionRawPeakHistogram;
  static int  Box[ResolutionRawPeakHistogram];
  static int  iLastBox = ResolutionRawPeakHistogram - 1;
  int         iBox, cum;
  Fprec       dBox, dDens;
  Fprec       Dyn_eD_CutOff, Fix_eD_CutOff;


  for (iBox = 0; iBox < nBox; iBox++)
    Box[iBox] = 0;

  dBox = (eDmaxRho - eDminRho) / (Fprec) nBox;

        nk =           Nz;
     nj_nk =      Ny * Nz;
  ni_nj_nk = Nx * Ny * Nz;

  PF = PeakFlags;

  for (iPF = 0; iPF < ni_nj_nk; iPF++, PF++)
    if (*PF & PF_HighBit)
      *PF &= (~0 ^ PF_IsPeakBit); /* reset IsPeakBit */

  eD = eDensity;
  PF = PeakFlags;

  im = ni_nj_nk - nj_nk;
  i0 = 0;
  ip = nj_nk;
  ibreak = ni_nj_nk;
  while (ip < ibreak)
  {
    jm = nj_nk - nk;
    j0 = 0;
    jp = nk;
    jbreak = nj_nk;
    while (jp < jbreak)
    {
      km = nk - 1;
      k0 = 0;
      kp = 1;
      kbreak = nk;
      while (kp < kbreak)
      {
        if (*PF & PF_HighBit) AsyPF = PF;
        else                  AsyPF = &PeakFlags[*PF];

        if (! (*AsyPF & PF_IsPeakBit))
        {
          if (Level >= 1) {
            /* m00 0m0 00m
               p00 0p0 00p
             */
            if ((*eD) < eDensity[im + j0 + k0]) goto NextGridPoint;
            if ((*eD) < eDensity[ip + j0 + k0]) goto NextGridPoint;
            if ((*eD) < eDensity[i0 + jm + k0]) goto NextGridPoint;
            if ((*eD) < eDensity[i0 + jp + k0]) goto NextGridPoint;
            if ((*eD) < eDensity[i0 + j0 + km]) goto NextGridPoint;
            if ((*eD) < eDensity[i0 + j0 + kp]) goto NextGridPoint;

            if (Level >= 2) {
              /* mm0 m0m 0mm mp0 m0p 0mp
                 pp0 p0p 0pp pm0 p0m 0pm
               */
              if ((*eD) < eDensity[im + jm + k0]) goto NextGridPoint;
              if ((*eD) < eDensity[ip + jp + k0]) goto NextGridPoint;
              if ((*eD) < eDensity[im + j0 + km]) goto NextGridPoint;
              if ((*eD) < eDensity[ip + j0 + kp]) goto NextGridPoint;
              if ((*eD) < eDensity[i0 + jm + km]) goto NextGridPoint;
              if ((*eD) < eDensity[i0 + jp + kp]) goto NextGridPoint;
              if ((*eD) < eDensity[im + jp + k0]) goto NextGridPoint;
              if ((*eD) < eDensity[ip + jm + k0]) goto NextGridPoint;
              if ((*eD) < eDensity[im + j0 + kp]) goto NextGridPoint;
              if ((*eD) < eDensity[ip + j0 + km]) goto NextGridPoint;
              if ((*eD) < eDensity[i0 + jm + kp]) goto NextGridPoint;
              if ((*eD) < eDensity[i0 + jp + km]) goto NextGridPoint;

              if (Level >= 3) {
                /* mmm mmp mpm mpp
                   ppp ppm pmp pmm
                 */
                if ((*eD) < eDensity[im + jm + km]) goto NextGridPoint;
                if ((*eD) < eDensity[ip + jp + kp]) goto NextGridPoint;
                if ((*eD) < eDensity[im + jm + kp]) goto NextGridPoint;
                if ((*eD) < eDensity[ip + jp + km]) goto NextGridPoint;
                if ((*eD) < eDensity[im + jp + km]) goto NextGridPoint;
                if ((*eD) < eDensity[ip + jm + kp]) goto NextGridPoint;
                if ((*eD) < eDensity[im + jp + kp]) goto NextGridPoint;
                if ((*eD) < eDensity[ip + jm + km]) goto NextGridPoint;
              }
            }
          }

          *AsyPF |= PF_IsPeakBit;

              dDens = (*eD) - eDminRho;
          if (dDens == 0. || dDens < dBox)
            iBox = 0;
          else {
                     iBox = (int)(dDens / dBox);
            if      (iBox > iLastBox) iBox = iLastBox;
            else if (iBox < 0)        iBox = 0;
          }
          Box[iBox] += WyckoffList[*AsyPF & PF_iWL_Mask].nPositions;
        }

        NextGridPoint:

        eD++;
        PF++;

        km = k0;
        k0 = kp;
        kp++;
        if (kp == nk) { kp = 0; kbreak = 1; }
      }
      jm = j0;
      j0 = jp;
      jp += nk;
      if (jp == nj_nk) { jp = 0; jbreak = nk; }
    }
    im = i0;
    i0 = ip;
    ip += nj_nk;
    if (ip == ni_nj_nk) { ip = 0; ibreak = nj_nk; }
  }


  cum = 0;
  iBox = nBox;
  while (iBox--)
  {
    cum += Box[iBox];
    if (cum > MaxRawPeaks) break;
  }
  iBox++;

  Dyn_eD_CutOff = (eDminRho + (Fprec) iBox * dBox)
                 + dBox * (Fprec) 1.e-4; /* compensate round-off errors */
                               /* ARBITRARY */

  Fix_eD_CutOff = eDensityCutOff.Value;

  if (eDensityCutOff.Percent)
    Fix_eD_CutOff *= .01 * eDmaxRho;

  if (Dyn_eD_CutOff < Fix_eD_CutOff)
      Dyn_eD_CutOff = Fix_eD_CutOff;

  if (Debug0) /* Debug: Print Dyn_eD_CutOff */
    Fprintf(stdout, "Dyn_eD_CutOff = %.1f\n", Dyn_eD_CutOff);

  if (F_eD_Histogram) PutHistogram(Box, nBox, NULL);

  return Dyn_eD_CutOff;
}


static void AtA_ni_10(Fprec *A, int ni, Fprec *AtA)
{
  const int  nj = 10;
  int        i, j, k;
  int        knjj; /* k * nj + j */
  Fprec      *A_j, *A_k, *eAtA;


  /* AtA[j][k] = Sum(i)(A[i][j] * A[i][k]); */

  eAtA = AtA;

  for (j = 0; j < nj; j++)
    for (k = 0, knjj = j; k < nj; k++, knjj += nj)
    {
      if (j > k)
        *eAtA = AtA[knjj];
      else
      {
        A_j = A + j;
        A_k = A + k;

        *eAtA = (*A_j) * (*A_k);

        for (i = 1; i < ni; i++)
        {
                     A_j += nj;
                              A_k += nj;
          *eAtA += (*A_j) * (*A_k);
        }
      }

      eAtA++;
    }
}


static void AtV_ni_10(Fprec *A, int ni, Fprec *V, Fprec *AtV)
{
  int    i;
  Fprec  *eAtV;


    eAtV = AtV;
  (*eAtV++) = (*A++) * (*V);
  (*eAtV++) = (*A++) * (*V);
  (*eAtV++) = (*A++) * (*V);
  (*eAtV++) = (*A++) * (*V);
  (*eAtV++) = (*A++) * (*V);
  (*eAtV++) = (*A++) * (*V);
  (*eAtV++) = (*A++) * (*V);
  (*eAtV++) = (*A++) * (*V);
  (*eAtV++) = (*A++) * (*V);
  (*eAtV  ) = (*A++) * (*V);

  for (i = 1; i < ni; i++)
  {
      eAtV = AtV;           V++;
    (*eAtV++) += (*A++) * (*V);
    (*eAtV++) += (*A++) * (*V);
    (*eAtV++) += (*A++) * (*V);
    (*eAtV++) += (*A++) * (*V);
    (*eAtV++) += (*A++) * (*V);
    (*eAtV++) += (*A++) * (*V);
    (*eAtV++) += (*A++) * (*V);
    (*eAtV++) += (*A++) * (*V);
    (*eAtV++) += (*A++) * (*V);
    (*eAtV  ) += (*A++) * (*V);
  }
}


static int ProjectShift(T_eD_PeakList *eEPL,
                        Fprec *xs, Fprec *ys, Fprec *zs)
{
  Fprec  *e, s[3], exs[3], n[3];
  Fprec  e_s, n_n, n_s;
  Fprec  cxs, cys, czs;
  Fprec  ParallelShift, PerpendicularShift;


  e = eEPL->WL_Entry->FreeVect.V;

  Tform_xc(*xs, *ys, *zs, s[0], s[1], s[2]);

  e_s = DotC(e, s);
  CrossC(e, s, exs);
  CrossC(e, exs, n);
  n_n = DotC(n, n);

  if (n_n == CF0)
    n_s = 0.; /* s is exactly parallel e */
  else
    n_s = DotC(n, s);

  if      (eEPL->WL_Entry->FreeVect.D == 1)
  {
    cxs = e_s * e[0];
    cys = e_s * e[1];
    czs = e_s * e[2];

    if (n_n != CF0) n_s /= AppSqrt(n_n);

    ParallelShift      = AppFabs(e_s);
    PerpendicularShift = AppFabs(n_s);
  }
  else if (eEPL->WL_Entry->FreeVect.D == 2)
  {
    if (n_n != CF0) n_s /= n_n;

    cxs = n_s * n[0];
    cys = n_s * n[1];
    czs = n_s * n[2];

    ParallelShift      = cxs * cxs + cys * cys + czs * czs;
    PerpendicularShift = AppFabs(e_s);
  }
  else
    InternalError("Superfluous call");

  if (Debug0) /* Debug: Print updated MaxPerpendicularShift */
  {
    static Fprec  MaxPerpendicularShift = -1.;

    if (MaxPerpendicularShift < PerpendicularShift)
    {
      MaxPerpendicularShift = PerpendicularShift;
        Fprintf(stdout,
          "MaxPerpendicularShift %.7g  %.7g  %8.4f %8.4f %8.4f\n",
          MaxPerpendicularShift, ParallelShift,
          *xs, *ys, *zs);
    }
  }

  Tform_cx(cxs, cys, czs, *xs, *ys, *zs);

  if ((    PerpendicularShift > ParallelShift * (Fprec) .5 /* ARBITRARY */
        && PerpendicularShift > (Fprec) .05)               /* ARBITRARY */
      || PerpendicularShift > (Fprec) .1)                  /* ARBITRARY */
    return -1;

  return 0;
}


#define About_ln_2 ((Fprec) .69)   /* ARBITRARY limits                    */
#define Min_eD_PfI ((Fprec) 1.e-6) /*   established to make interpolation */
                                   /*   numerically stable                */

#define Add_PfI(i_, j_, k_, x, y, z)\
  {\
    i = i_ + j_ + k_;\
    if ((PeakFlags[i] & PF_HighBit) == 0) i = PeakFlags[i];\
    if (eDensity[i] > Min_eD_PfI)\
    {\
      /* ln rho */ *OV++ = AppLog(eDensity[i]);\
           /* a */ *DMx++ = 1.;\
           /* b */ *DMx++ = x;\
           /* c */ *DMx++ = y;\
           /* d */ *DMx++ = z;\
           /* e */ *DMx++ = x * x;\
           /* f */ *DMx++ = y * y;\
           /* g */ *DMx++ = z * z;\
           /* h */ *DMx++ = y * z;\
           /* k */ *DMx++ = x * z;\
           /* l */ *DMx++ = x * y;\
      nPfI++;\
    }\
  }


static void InterpolatePeakPosition(const Fprec *eDensity,
                                    const T_PeakFlags *PeakFlags,
                                    T_eD_PeakList *eEPL,
                                    int ix0, int iy0, int iz0)
{
  int    ixm,    ixp, iym,    iyp, izm,    izp;
  Fprec   xm, x0, xp,  ym, y0, yp,  zm, z0, zp;
  int    nPfI, IllDefined, i0, i;
  Fprec  *DMx, *OV, *CVc;
  Fprec  a, b, c, d, e, f, g, h, k, l;
  Fprec  e_, f_, g_, h_, k_, l_;
  Fprec  ef, g4, h2, l2, lk, k2;
  Fprec  Denom, TwoDenom;
  Fprec  xxOff, yxOff, zxOff, xx, yx, zx;
  Fprec  xs, ys, zs, In, Mx;
  Fprec  Bog_xs, Bog_ys, Bog_zs;

#define MaxPfI (19)
#define nCoeff (10)

         Fprec  ObservationVector[MaxPfI];
  static Fprec  DesignMatrix[MaxPfI * nCoeff];
  static Fprec  LU_NormalMatrix[nCoeff * nCoeff];
         int    LU_Permutation[nCoeff];
  static Fprec  LU_19NormalMatrix[nCoeff * nCoeff];
  static int    LU_19Permutation[nCoeff];
  static int    LU_19Set = 0;
         Fprec  *LU_N;
         int    *LU_P;
         Fprec  CoeffVector[nCoeff];


  IllDefined = 0;

  xxOff = (Fprec) ix0 / Nx;
  yxOff = (Fprec) iy0 / Ny;
  zxOff = (Fprec) iz0 / Nz;

  if (eEPL->WL_Entry->FreeVect.D == 0)
    i0 = (ix0 * Ny + iy0) * Nz + iz0;

  else
  {
    x0 = y0 = z0 = 0.;
    xp = 1. / Nx;
    xm = -xp;
    yp = 1. / Ny;
    ym = -yp;
    zp = 1. / Nz;
    zm = -zp;

    ixp = ix0 + 1;
    if (ix0 == 0)
      ixm = Nx - 1;
    else
    {
      ixm = ix0 - 1;
      if (ixp == Nx) ixp = 0;
    }

    ixm *= Ny * Nz;
    ix0 *= Ny * Nz;
    ixp *= Ny * Nz;

    iyp = iy0 + 1;
    if (iy0 == 0)
      iym = Ny - 1;
    else
    {
      iym = iy0 - 1;
      if (iyp == Ny) iyp = 0;
    }

    iym *= Nz;
    iy0 *= Nz;
    iyp *= Nz;

    izp = iz0 + 1;
    if (iz0 == 0)
      izm = Nz - 1;
    else
    {
      izm = iz0 - 1;
      if (izp == Nz) izp = 0;
    }

    /* Design Matrix layout according to Numerical Recipes:

                  names of coefficients: see
                    J.S. Rollett (ed.) (1965)
                    Computing Methods in Crystallography, pp. 35-37

                  a b c d e f g h k l
             000
             m00
             p00
             0m0
             0p0
             00m
             00p
             mm0
             pp0
             m0m
             p0p
             0mm
             0pp
             mp0
             pm0
             m0p
             p0m
             0mp
             0pm
     */

    OV = ObservationVector;
    DMx = DesignMatrix;
    nPfI = 0;

    Add_PfI(ix0, iy0, iz0,  x0, y0, z0); i0 = i;
    Add_PfI(ixm, iy0, iz0,  xm, y0, z0);
    Add_PfI(ixp, iy0, iz0,  xp, y0, z0);
    Add_PfI(ix0, iym, iz0,  x0, ym, z0);
    Add_PfI(ix0, iyp, iz0,  x0, yp, z0);
    Add_PfI(ix0, iy0, izm,  x0, y0, zm);
    Add_PfI(ix0, iy0, izp,  x0, y0, zp);
    Add_PfI(ixm, iym, iz0,  xm, ym, z0);
    Add_PfI(ixp, iyp, iz0,  xp, yp, z0);
    Add_PfI(ixm, iy0, izm,  xm, y0, zm);
    Add_PfI(ixp, iy0, izp,  xp, y0, zp);
    Add_PfI(ix0, iym, izm,  x0, ym, zm);
    Add_PfI(ix0, iyp, izp,  x0, yp, zp);
    Add_PfI(ixm, iyp, iz0,  xm, yp, z0);
    Add_PfI(ixp, iym, iz0,  xp, ym, z0);
    Add_PfI(ixm, iy0, izp,  xm, y0, zp);
    Add_PfI(ixp, iy0, izm,  xp, y0, zm);
    Add_PfI(ix0, iym, izp,  x0, ym, zp);
    Add_PfI(ix0, iyp, izm,  x0, yp, zm);

    if      (nPfI < MinPfI)
    {
      IllDefined = -1;
      goto RetainPositionOnGrid;
    }
    else if (nPfI == LU_19Set)
    {
      LU_N = LU_19NormalMatrix;
      LU_P = LU_19Permutation;
    }
    else
    {
      LU_N = LU_NormalMatrix;
      LU_P = LU_Permutation;
      AtA_ni_10(DesignMatrix, nPfI, LU_N);
      if (MxLU_Decomposition(LU_N, nCoeff, LU_P) != 0)
      {
        IllDefined = -1;
        goto RetainPositionOnGrid;
      }

      if (nPfI == MaxPfI)
      {
        for (i = 0; i < nCoeff * nCoeff; i++)
          LU_19NormalMatrix[i] = LU_N[i];
        for (i = 0; i < nCoeff; i++)
          LU_19Permutation[i] = LU_P[i];
        LU_19Set = MaxPfI;
      }
    }

    AtV_ni_10(DesignMatrix, nPfI, ObservationVector, CoeffVector);
    MxLU_BackSubstitution(LU_N, nCoeff, LU_P, CoeffVector);

    /* equation to solve d(lnRho)/dx = d(lnRho)/dy = d(lnRho)/dz = 0:

         [2e   l   k] [x]   [-b]
         [ l  2f   h] [y] = [-c]
         [ k   h  2g] [z]   [-d]

       calculation is done with this matrix:

                            [   e    1/2 l  1/2 k ]
                            [                     ]
                      A :=  [ 1/2 l    f    1/2 h ]
                            [                     ]
                            [ 1/2 k  1/2 h    g   ]

       inverse of A is:

                    [            2                                ]
                    [   4 f g - h      2 l g - k h   l h - 2 k f  ]
                    [   ----------   - -----------   -----------  ]
                    [       %1              %1            %1      ]
                    [                                             ]
                    [                           2                 ]
                    [   2 l g - k h    4 e g - k      2 e h - l k ]
              Ai := [ - -----------    ----------   - ----------- ]
                    [        %1            %1              %1     ]
                    [                                             ]
                    [                                          2  ]
                    [  l h - 2 k f     2 e h - l k    4 e f - l   ]
                    [  -----------   - -----------    ----------  ]
                    [       %1              %1            %1      ]


                                   2    2              2
                %1 := 4 e f g - e h  - l  g + l k h - k  f := 4 det(A)


                   [ b ]
               u = [ c ]
                   [ d ]
                                             t  -1
                          3/2   exp(a - 1/4(u  A   u))
       Integral(Peak) = PI     ------------------------
                                     sqrt(det(A))


       From: stewart@cs.umd.edu (G. W. Stewart)
       Newsgroups: sci.math.num-analysis
       Subject: Re: Q: eigenvalues of a symmetric real 3*3 matrix
       Date: 8 Dec 1993 08:03:50 -0500
       Organization:
       U of Maryland, Dept. of Computer Science, Coll. Pk., MD 20742
       Message-ID: <2e4jbm$6s2@thales.cs.umd.edu>

                                          [ e  L  K ]
                                          [         ]
                                     C := [ L  f  H ]
                                          [         ]
                                          [ K  H  g ]

       ... I seem to recall that the person who posted the problem
       originally, later said that what he really needed to know was whether
       all the eigenvalues were negative.  This can be tested by performing
       Gaussian elimination on C.  Specifically, if e, f, and g are not all
       negative then the eigenvalues are not all negative.  Now form the
       matrix

          [ f' H' ]   [ f-LL/e  H-LK/e]
          [       ] = [               ].
          [ *  g' ]   [   *     g-KK/e]

       If f' and g' are not negative, then the eigenvalues are not all
       negative.  Finally form

          g'' = g' - H'H'/f'.

       If g'' is not negative then the eigenvalues are not all negative;
       otherwise they all are negative.  This test is due to Lagrange.

     */

         CVc = CoeffVector;
    a = *CVc++;
    b = *CVc++;
    c = *CVc++;
    d = *CVc++;
    e = *CVc++;
    f = *CVc++;
    g = *CVc++;
    h = *CVc++;
    k = *CVc++;
    l = *CVc;

    /* Specifically, if e, f, and g are not all negative then the
       eigenvalues are not all negative.
     */

    if (   e  >= CF0
        || f  >= CF0
        || g  >= CF0) { IllDefined = -1; goto RetainPositionOnGrid; }

    ef = e * f;
    g4 = g * 4;
    h2 = h * h;
    l2 = l * l;
    lk = l * k;
    k2 = k * k;

    Denom = -ef * g4 + e * h2 + l2 * g - lk * h + k2 * f;

    if (Denom <= CF0) { IllDefined = -1; goto RetainPositionOnGrid; }

    /* Now form the matrix...
     */
                  xx = (Fprec) .25 / e;
    f_ = f - l2 * xx;
    g_ = g - k2 * xx;

    if (   f_ >= CF0
        || g_ >= CF0) { IllDefined = -1; goto RetainPositionOnGrid; }

    h_ = h * (Fprec) .5 - lk * xx;
    g_ = g_ - h_ * h_ / f_;

    if (   g_ >= CF0) { IllDefined = -1; goto RetainPositionOnGrid; }

    /* ...otherwise they all are negative.
     */

    /* division by Denom is omitted here */
    e_ = (f * g4 - h2);
    l_ = (k * h - 2 * l * g);
    k_ = (l * h - 2 * k * f);
    f_ = (e * g4 - k2);
    h_ = (lk - 2 * e * h);
    g_ = (ef * 4 - l2);

    TwoDenom = Denom + Denom;
    xs = (e_ * b + l_ * c + k_ * d) / TwoDenom;
    ys = (l_ * b + f_ * c + h_ * d) / TwoDenom;
    zs = (k_ * b + h_ * c + g_ * d) / TwoDenom;

    if (xm > xs || xs > xp || ym > ys || ys > yp || zm > zs || zs > zp) {
      IllDefined = -1; goto RetainPositionOnGrid;
    }

    Bog_xs = AppFabs(xp) * (Fprec) 1.e-3;
    Bog_ys = AppFabs(yp) * (Fprec) 1.e-3; /* ARBITRARY */
    Bog_zs = AppFabs(zp) * (Fprec) 1.e-3;

    if (AppFabs(xs) < Bog_xs) xs = 0.;
    if (AppFabs(ys) < Bog_ys) ys = 0.;
    if (AppFabs(zs) < Bog_zs) zs = 0.;

                    h_ = (double) b * xs + (double) c * ys + (double) d * zs;
    Mx =   (double) a
         + (double) h_
         + (double) e * xs * xs + (double) f * ys * ys + (double) g * zs * zs
         + (double) h * ys * zs + (double) k * xs * zs + (double) l * xs * ys;

    if (Mx > ObservationVector[0] + About_ln_2) {
      IllDefined = -1; goto RetainPositionOnGrid;
    }

    Mx = AppExp(Mx);

    In = a + h_ * (Fprec) .5;
    In =   (Fprec)(2. * pow(M_PI, 1.5) * LatConD.v)
         * (AppExp(In) / AppSqrt(Denom));

    if (   eEPL->WL_Entry->FreeVect.D == 1
        || eEPL->WL_Entry->FreeVect.D == 2)
      if (ProjectShift(eEPL, &xs, &ys, &zs) != 0) nPfI = -nPfI;

        xx = xxOff + xs;
    if (xx <  CF0) xx += CF1;
    if (xx >= CF1) xx -= CF1;
        yx = yxOff + ys;
    if (yx <  CF0) yx += CF1;
    if (yx >= CF1) yx -= CF1;
        zx = zxOff + zs;
    if (zx <  CF0) zx += CF1;
    if (zx >= CF1) zx -= CF1;

    eEPL->Position.x = xx;
    eEPL->Position.y = yx;
    eEPL->Position.z = zx;
    eEPL->Grid_eD  = eDensity[i0];
    eEPL->Maximum  = Mx;
    eEPL->Integral = In;
    eEPL->nPfI = nPfI;

    /* NormOf inserted 08/31/99. See Reconstruct_fUC_Shift() */
    NormOf(eEPL->Position.x);
    NormOf(eEPL->Position.y);
    NormOf(eEPL->Position.z);

    if (Debug0) /* Debug: IPP old => new */
    {
      /* MxDump(ObservationVector, 1, abs(nPfI), "ObservationVector"); */
      Fprintf(stdout,
        "IPP  %7.4f %7.4f %7.4f  =>  %7.4f %7.4f %7.4f\n",
        xxOff, yxOff, zxOff, xx, yx, zx);
    }

    SucPeakInterpolation++;
    return;
  }

  RetainPositionOnGrid:

  eEPL->Position.x = xxOff;
  eEPL->Position.y = yxOff;
  eEPL->Position.z = zxOff;
  eEPL->Grid_eD  = eDensity[i0];
  eEPL->Maximum  = eDensity[i0];
  eEPL->Integral = eDensity[i0];

  if (IllDefined == 0)
    eEPL->nPfI = 0;
  else
  {
    eEPL->nPfI = -1;
    IllPeakInterpolation++;
  }

  return;
}

#undef About_ln_2
#undef Min_eD_PfI
#undef Add_PfI
#undef MaxPfI
#undef nCoeff


static void InterSectionII(Fprec a0, Fprec a1, Fprec a2, Fprec p[3],
                           Fprec b0, Fprec b1, Fprec b2, Fprec q[3],
                           Fprec *x0, Fprec *x1, Fprec *x2)
{
  Fprec  a[3], b[3], u[3], b_a[3], x[3], ka, kb, ku, aka, akb;
  Fprec  A[9], InvA[9], detA, pxu[3], qxu[3];


  Tform_xc(a0, a1, a2, a[0], a[1], a[2]);
  Tform_xc(b0, b1, b2, b[0], b[1], b[2]);

  CrossC(p, q, u);
  CrossC(p, u, pxu); if (SetToUnitLength(pxu) == 0) InternalError(NULL);
  CrossC(q, u, qxu); if (SetToUnitLength(qxu) == 0) InternalError(NULL);

  A[0] = pxu[0];  A[1] = -qxu[0];  A[2] = u[0];
  A[3] = pxu[1];  A[4] = -qxu[1];  A[5] = u[1];
  A[6] = pxu[2];  A[7] = -qxu[2];  A[8] = u[2];

      detA = deter33fMx(A);
  if (detA == 0.)
    InternalError("InterSectionII(): detA = 0");

  Inverse33fMx(A, InvA);

  b_a[0] = b[0] - a[0];
  b_a[1] = b[1] - a[1];
  b_a[2] = b[2] - a[2];

  MxMultiply(x, InvA, b_a, 3, 3, 1);

  ka = x[0] / detA; aka = AppFabs(ka);
  kb = x[1] / detA; akb = AppFabs(kb);
  ku = x[2] / detA;

  x[0] = a[0];
  x[1] = a[1];
  x[2] = a[2];

  if (aka != 0.)
  {
    ku *= aka / (aka + akb);

    x[0] += ka * pxu[0] + ku * u[0];
    x[1] += ka * pxu[1] + ku * u[1];
    x[2] += ka * pxu[2] + ku * u[2];
  }

  Tform_cx(x[0], x[1], x[2], *x0, *x1, *x2);

  if (Debug0) /* Print I-II a, b, x */
  {
    Fprintf(stdout, "I-II a = %.5f %.5f %.5f [%.6g %.6g %.6g]\n",
      a0, a1, a2, p[0], p[1], p[2]);
    Fprintf(stdout, "I-II b = %.5f %.5f %.5f [%.6g %.6g %.6g]\n",
      b0, b1, b2, q[0], q[1], q[2]);
    Fprintf(stdout, "I-II x = %.5f %.5f %.5f\n",
      *x0, *x1, *x2);
  }
}


static void InterSectionIP(Fprec a0, Fprec a1, Fprec a2, Fprec p[3],
                           Fprec b0, Fprec b1, Fprec b2, Fprec q[3],
                           Fprec *x0, Fprec *x1, Fprec *x2)
{
  /* Plane: p.(x-a)=0
     Line:  b+k*q=x   => p.(b+k*q-a)=0
                      => k*q.p+p.(b-a)=0
                      => k=-(p.(b-a)/(q.p))
   */

  Fprec  qp, x[3], b_a[3], b_ap, k;


      qp = DotC(q, p);
  if (qp == 0.)
    InternalError("InterSectionIP(): qp = 0");

  x[0] = b0 - a0;
  x[1] = b1 - a1;
  x[2] = b2 - a2;

  Tform_xc(x[0], x[1], x[2], b_a[0], b_a[1], b_a[2]);

  b_ap = DotC(b_a, p);

  k = -b_ap / qp;

  Tform_xc(b0, b1, b2, x[0], x[1], x[2]);

  x[0] += k * q[0];
  x[1] += k * q[1];
  x[2] += k * q[2];

  Tform_cx(x[0], x[1], x[2], *x0, *x1, *x2);

  if (Debug0) /* Print I-IP a, b, x */
  {
    Fprintf(stdout, "I-IP a = %.5f %.5f %.5f [%.6g %.6g %.6g]\n",
      a0, a1, a2, p[0], p[1], p[2]);
    Fprintf(stdout, "I-IP b = %.5f %.5f %.5f [%.6g %.6g %.6g]\n",
      b0, b1, b2, q[0], q[1], q[2]);
    Fprintf(stdout, "I-IP x = %.5f %.5f %.5f\n",
      *x0, *x1, *x2);
  }
}


static int  InterSectionPP(Fprec a0, Fprec a1, Fprec a2, Fprec p[3],
                           Fprec b0, Fprec b1, Fprec b2, Fprec q[3],
                           Fprec *x0, Fprec *x1, Fprec *x2)
{
  /*
      Find the point at the middle of the shortest vector between the
      skew lines a+n*p and b+m*q.
   */

  Fprec  pxq[3], A[9], detA, InvA[9], b_a[3], x[3], k, m, d2;


  CrossC(p, q, pxq);

  A[0] = p[0];  A[1] = -q[0];  A[2] = pxq[0];
  A[3] = p[1];  A[4] = -q[1];  A[5] = pxq[1];
  A[6] = p[2];  A[7] = -q[2];  A[8] = pxq[2];

      detA = deter33fMx(A);
  if (detA == 0.)
    InternalError("InterSectionPP(): detA = 0");

  Inverse33fMx(A, InvA);

  x[0] = b0 - a0;
  x[1] = b1 - a1;
  x[2] = b2 - a2;

  Tform_xc(x[0], x[1], x[2], b_a[0], b_a[1], b_a[2]);

  MxMultiply(x, InvA, b_a, 3, 3, 1);

  k = x[0] / detA;
  m = x[2] / detA;

  d2 = m * m * DotC(pxq, pxq);

  if (d2 > (Fprec)(Square(1.e-4))) /* ARBITRARY */
    return 0; /* lines do not intersect */

  Tform_xc(a0, a1, a2, x[0], x[1], x[2]);

                     m *= .5;
  x[0] += k * p[0] + m * pxq[0];
  x[1] += k * p[1] + m * pxq[1];
  x[2] += k * p[2] + m * pxq[2];

  Tform_cx(x[0], x[1], x[2], *x0, *x1, *x2);

  if (Debug0) /* Print I-PP a, b, x */
  {
    Fprintf(stdout, "I-PP a = %.5f %.5f %.5f [%.6g %.6g %.6g]\n",
      a0, a1, a2, p[0], p[1], p[2]);
    Fprintf(stdout, "I-PP b = %.5f %.5f %.5f [%.6g %.6g %.6g]\n",
      b0, b1, b2, q[0], q[1], q[2]);
    Fprintf(stdout, "I-PP x = %.5f %.5f %.5f\n",
      *x0, *x1, *x2);
  }

  return 1;
}


static int MoreSpecial(T_eD_PeakList *EPL,
                       T_FreeVect *CurFV,
                       int iSymOp, int iList, int iLoopInv, int iTrV,
                       T_iVector *UC_Shift0,
                       T_iVector *UC_Shift,
                       const T_iVector *m0pShift)
{
  const T_RTMx  *lsmx;
  const int     *TrV;
  T_RTMx        SMx;
  int           i, Mul;
  int           CumRMx[9], ScGl[3], ShTrMx[9], Tr[3], RmI[9], RmI_Tr[3];
  Fprec         pxx, pyx, pzx, nxx, nyx, nzx, dx, dy, dz, d[3], absd2;
  int           PosCorrType;
  T_FreeVect    NewFV[1], SymFV[1];


  if (CurFV->D == 0)
    return 0;

  if (Debug0) /* Debug: Print Before MoreSpecial */
    Fprintf(stdout, "Before MoreSpecial  %8.4f %8.4f %8.4f\n",
      EPL->Position.x, EPL->Position.y, EPL->Position.z);

  TrV = &SpgrInfo.LatticeInfo->TrVector[iTrV * 3];

  if (Debug0) /* Debug: Print MoreSpecial() internals */
  {
    Fprintf(stdout,
      "WL_Entry %2d:  NumberOfPositions %3d  before MoreSpecial\n",
      (int)(EPL->WL_Entry - WyckoffList),
      EPL->WL_Entry->nPositions);

    Fprintf(stdout, "iSymOp = %d  iList = %d  iLoopInv = %d  iTrV = %d\n",
      iSymOp, iList, iLoopInv, iTrV);

    Fprintf(stdout, "+TrV       = %3d %3d %3d\n",
      TrV[0], TrV[1], TrV[2]);

    Fprintf(stdout, "+UC_Shift  = %3d %3d %3d\n",
      UC_Shift->x,
      UC_Shift->y,
      UC_Shift->z);

    Fprintf(stdout, "-UC_Shift0 = %3d %3d %3d\n",
      UC_Shift0->x,
      UC_Shift0->y,
      UC_Shift0->z);

    Fprintf(stdout, "+m0pShift  = %3d %3d %3d\n",
      m0pShift->x, m0pShift->y, m0pShift->z);
  }

  lsmx = &SpgrInfo.ListSeitzMx[iList];

  if (iLoopInv == 0)
    for (i = 0; i < 12; i++) SMx.a[i] =  lsmx->a[i];
  else
    for (i = 0; i < 12; i++) SMx.a[i] = -lsmx->a[i];

  SMx.s.T[0] += TrV[0] + UC_Shift->x - UC_Shift0->x + m0pShift->x;
  SMx.s.T[1] += TrV[1] + UC_Shift->y - UC_Shift0->y + m0pShift->y;
  SMx.s.T[2] += TrV[2] + UC_Shift->z - UC_Shift0->z + m0pShift->z;

  Mul = BuildCumRMx(&SpgrInfo, iList, iLoopInv, CumRMx);

  RotMx_t_Vector(ScGl, CumRMx, SMx.s.T, 0);

  for (i = 0; i < 3; i++) {
    if (ScGl[i] %  Mul) InternalError("Corrupt ScGl");
        ScGl[i] /= Mul;
        ScGl[i] = iModPositive(ScGl[i], STBF);
    if (ScGl[i] != 0) {
      Fprintf(stdout, "# WARNING: MoreSpecial() <-> ScGl != 0\n");
      CountWarnings++;
      return 0;
    }
  }

  if (GetRotMxInfo(SMx.s.R, NULL, ShTrMx) == 0)
    InternalError(NULL);

  RotMx_t_Vector(Tr, ShTrMx, SMx.s.T, 0);

  for (i = 0; i < 3; i++) {
    if (Tr[i] %  (STBF * STBF / CTBF)) InternalError("Corrupt Tr");
        Tr[i] /= (STBF * STBF / CTBF);
  }

#ifndef NO_EXTRA_CHECKS
  SetRminusI(SMx.s.R, RmI, 0);
  RotMx_t_Vector(RmI_Tr, RmI, Tr, 0);

  for (i = 0; i < 3; i++)
    if (   iModPositive(RmI_Tr[i], CTBF)
        != iModPositive(SMx.s.T[i] * (CTBF / STBF), CTBF))
      InternalError("MoreSpecial() <-> SMx.s.T");
#endif

  SetRminusI(CumRMx, RmI, 0);
  RotMx_t_Vector(RmI_Tr, RmI, Tr, 0);

  pxx = EPL->Position.x;
  pyx = EPL->Position.y;
  pzx = EPL->Position.z;
  nxx = (CumRMx[0] * pxx + CumRMx[1] * pyx + CumRMx[2] * pzx) / Mul;
  nyx = (CumRMx[3] * pxx + CumRMx[4] * pyx + CumRMx[5] * pzx) / Mul;
  nzx = (CumRMx[6] * pxx + CumRMx[7] * pyx + CumRMx[8] * pzx) / Mul;
  nxx += RmI_Tr[0] / fCTBF;
  nyx += RmI_Tr[1] / fCTBF;
  nzx += RmI_Tr[2] / fCTBF;

  PosCorrType = CombineFreeEigenVects(CurFV, SymFV, NewFV, iList, iLoopInv);

  switch (PosCorrType)
  {
    case PCT_None:
      break;
    case PCT_NoChange:
      dx = nxx - pxx;
      dy = nyx - pyx;
      dz = nzx - pzx;
      Tform_xc(dx, dy, dz, d[0], d[1], d[2]);
      absd2 = DotC(d, d);
      if (absd2 > (Fprec)(Square(1.e-4))) { /* ARBITRARY */
        Fprintf(stdout,
          "# WARNING: MoreSpecial() <-> PCT_NoChange vs. absd2\n");
        CountWarnings++;
        return 0;
      }
      break;
    case PCT_II:
      InterSectionII(pxx, pyx, pzx, CurFV->V,
                     nxx, nyx, nzx, SymFV->V,
                     &nxx, &nyx, &nzx);
      CountInterSectionII++;
      break;
    case PCT_IP:
      InterSectionIP(pxx, pyx, pzx, CurFV->V,
                     nxx, nyx, nzx, SymFV->V,
                     &nxx, &nyx, &nzx);
      CountInterSectionIP++;
      break;
    case PCT_PI:
      InterSectionIP(nxx, nyx, nzx, SymFV->V,
                     pxx, pyx, pzx, CurFV->V,
                     &nxx, &nyx, &nzx);
      CountInterSectionPI++;
      break;
    case PCT_PP:
      if (InterSectionPP(pxx, pyx, pzx, CurFV->V,
                         nxx, nyx, nzx, SymFV->V,
                         &nxx, &nyx, &nzx) == 0) {
        Fprintf(stdout,
          "# WARNING: MoreSpecial() <-> InterSectionPP() == 0\n");
        CountWarnings++;
        return 0;
      }
      CountInterSectionPP++;
      break;
    default:
      InternalError("Corrupt PosCorrType");
  }

  (void) memcpy(CurFV, NewFV, sizeof (*CurFV));

  EPL->Position.x = nxx;
  EPL->Position.y = nyx;
  EPL->Position.z = nzx;

  if (Debug0) /* Debug: Print After  MoreSpecial */
    Fprintf(stdout, "After  MoreSpecial  %8.4f %8.4f %8.4f\n",
      EPL->Position.x, EPL->Position.y, EPL->Position.z);

  NormOf(EPL->Position.x);
  NormOf(EPL->Position.y);
  NormOf(EPL->Position.z);

  return 1;
}


static void SetSymEquiv(T_SgInfo *SgInfo, T_eD_PeakList *EPL, int nEPL)
{
  int              iEPL, iEquiv;
  int              iSymOp, iSymTr, iList, nTrV, iTrV, nLoopInv, iLoopInv;
  T_fVector        RT, Eq, *fTrV;
  T_iVector        UC_Shift0, UC_Shift;
  const T_iVector  *m0pShift;
  T_fVector        *xSymEquiv;
  T_fVector        *cSymEquiv;
  const T_fRTMx    *lfsmx;
  int              nSymOp, Buf_WL_Flag[48];
  int              PosUpdated, WL_F_Updated;
  Fprec            AbsF_dot_F;
  T_FreeVect       FreeVect;
  Fprec            dx0, dy0, dz0, dx, dy, dz, Dist2, *NPMx;
  T_WyckoffList    **SubWL;
  int              iSource, iTarget;

#include "mindist2.h"


  nSymOp = SgInfo->OrderP;

  if (nSymOp > sizeof Buf_WL_Flag / sizeof (*Buf_WL_Flag))
    InternalError("Buf_WL_Flag too small");

  nLoopInv = Sg_nLoopInv(SgInfo);
  nTrV = SgInfo->LatticeInfo->nTrVector;

  NPMx = NextPeakMx;

  for (iEPL = 0; iEPL < nEPL; iEPL++)
  {
    for (iSymOp = 0; iSymOp < nSymOp; iSymOp++)
      Buf_WL_Flag[iSymOp] = EPL->WL_Entry->WL_Flag[iSymOp];

    WL_F_Updated = 0;

    FreeVect.V[0] = EPL->WL_Entry->FreeVect.V[0];
    FreeVect.V[1] = EPL->WL_Entry->FreeVect.V[1];
    FreeVect.V[2] = EPL->WL_Entry->FreeVect.V[2];
    FreeVect.D    = EPL->WL_Entry->FreeVect.D;

    PosUpdated = 0;

    Restart: /* after MoreSpecial() */

    *NPMx = MaxLatticeTr2;

    xSymEquiv = EPL->xSymEquiv;
    cSymEquiv = EPL->cSymEquiv;

    iEquiv =
    iSymOp =
    iSymTr = 0;

    for (iList = 0; iList < SpgrInfo.nList; iList++)
    {
      lfsmx = NULL;

      for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
      {
        if (Buf_WL_Flag[iSymOp] != WL_F_Active)
        {
          iSymTr += nTrV;
          iSymOp++;
        }
        else
        {
          if (lfsmx == NULL)
          {                      lfsmx = &List_fSeitzMx[iList];
            fRTMx_t_fVector(&RT, lfsmx, &EPL->Position);
          }

          for (iTrV = 0, fTrV = fTrVector; iTrV < nTrV; iTrV++, fTrV++)
          {
            if (iLoopInv == 0)
            {
              Eq.x =  RT.x + fTrV->x;
              Eq.y =  RT.y + fTrV->y;
              Eq.z =  RT.z + fTrV->z;
            }
            else
            {
              Eq.x = -RT.x + fTrV->x;
              Eq.y = -RT.y + fTrV->y;
              Eq.z = -RT.z + fTrV->z;
            }

            if (iSymTr == 0)
            {
              CountNormOf(Eq.x, UC_Shift0.x); xSymEquiv->x = Eq.x;
              CountNormOf(Eq.y, UC_Shift0.y); xSymEquiv->y = Eq.y;
              CountNormOf(Eq.z, UC_Shift0.z); xSymEquiv->z = Eq.z;
            }
            else
            {
              CountNormOf(Eq.x, UC_Shift.x);  xSymEquiv->x = Eq.x;
              CountNormOf(Eq.y, UC_Shift.y);  xSymEquiv->y = Eq.y;
              CountNormOf(Eq.z, UC_Shift.z);  xSymEquiv->z = Eq.z;
            }

            Tform_xc(Eq.x, Eq.y, Eq.z,
                     cSymEquiv->x, cSymEquiv->y, cSymEquiv->z);

            if (Debug0) /* Debug: Print xSymEquiv cSymEquiv */
            {
              T_iVector  *UCS = (iSymTr == 0 ? &UC_Shift0 : &UC_Shift);
              Fprintf(stdout,
              " %3d %2d  %8.4f %8.4f %8.4f  %8.4f %8.4f %8.4f  %2d %2d %2d\n",
                iEPL, iSymOp,
                xSymEquiv->x, xSymEquiv->y, xSymEquiv->z,
                cSymEquiv->x, cSymEquiv->y, cSymEquiv->z,
                UCS->x / STBF, UCS->y / STBF, UCS->z / STBF);
            }

            PosUpdated = 0;

#define EvalDist2(dx_, dy_, dz_)\
  {\
    Dist2 = dx_ * dx_ + dy_ * dy_ + dz_ * dz_;\
    if (*NPMx > Dist2) *NPMx = Dist2;\
    if (Dist2 <= CatchDistance2)\
    {\
      PosUpdated = MoreSpecial(EPL,\
                               &FreeVect,\
                               iSymOp, iList, iLoopInv, iTrV,\
                               &UC_Shift0, &UC_Shift, m0pShift);\
      Buf_WL_Flag[iSymOp] = WL_F_Identity;\
      WL_F_Updated = 1;\
      if (PosUpdated != 0)\
        goto Restart;\
    }\
    m0pShift++;\
  }
            m0pShift = List_m0pShift;

            dx0 = cSymEquiv->x - EPL->cSymEquiv->x;
            dy0 = cSymEquiv->y - EPL->cSymEquiv->y;
            dz0 = cSymEquiv->z - EPL->cSymEquiv->z;

            if (iSymTr) EvalDist2(dx0, dy0, dz0);

#include "mindist2.c" /* hand-made inline code */
#undef EvalDist2

            iSymTr++;
            xSymEquiv++;
            cSymEquiv++;
            iEquiv++;

            if (iEquiv == EPL->WL_Entry->nPositions)
              goto AllThrough;
          }

          iSymOp++;
        }
      }
    }

    AllThrough:

    if (WL_F_Updated != 0)
    {
      SubWL = EPL->WL_Entry->SubWL;

      while (*SubWL != NULL)
      {
        if ((*SubWL)->FreeVect.D == FreeVect.D)
        {
          for (iSymOp = 0; iSymOp < nSymOp; iSymOp++)
          {
            if (         Buf_WL_Flag[iSymOp] == WL_F_Identity
                && (*SubWL)->WL_Flag[iSymOp] != WL_F_Identity) break;
            if (   (*SubWL)->WL_Flag[iSymOp] == WL_F_Active
                &&       Buf_WL_Flag[iSymOp] != WL_F_Active) break;
          }

          if (iSymOp == nSymOp)
          {
            if (FreeVect.D == 0)
            {
              EPL->WL_Entry = (*SubWL);
              break;
            }
            else
            {
              AbsF_dot_F = DotC(FreeVect.V, (*SubWL)->FreeVect.V);
              AbsF_dot_F = AppFabs(AbsF_dot_F);

              if (AbsF_dot_F > (Fprec) (1. - 1.e-5)) /* ARBITRARY */
              {
                EPL->WL_Entry = (*SubWL);
                break;
              }
            }
          }
        }

        SubWL++;
      }

      if (*SubWL == NULL)
        progerror(
          "Can't find a position with Self-Distance > Catch-Distance\n"
          /* "Check your input." */ );

      if (Debug0) /* Debug: Print Updated WL_Entry */
        Fprintf(stdout,
          "WL_Entry %2d:  NumberOfPositions %3d  Updated\n",
          (int)(EPL->WL_Entry - WyckoffList),
          EPL->WL_Entry->nPositions);

      iSource = iTarget = 0;

      for (iSymOp = 0; iSymOp < nSymOp; iSymOp++)
      {
        if      (EPL->WL_Entry->WL_Flag[iSymOp] == WL_F_Active)
        {
          for (iTrV = 0; iTrV < nTrV; iTrV++)
          {
            if (iSource != iTarget)
            {
              EPL->xSymEquiv[iTarget] = EPL->xSymEquiv[iSource];
              EPL->cSymEquiv[iTarget] = EPL->cSymEquiv[iSource];
            }

            iSource++;
            iTarget++;
          }
        }
        else if (           Buf_WL_Flag[iSymOp] == WL_F_Active)
          iSource += nTrV;
      }

      if (iTarget != EPL->WL_Entry->nPositions)
        InternalError(NULL);

      nCatchRawPeak++;
    }

#ifndef NO_EXTRA_CHECKS
    CheckWL_Entry(EPL);
#endif

    NPMx += (NeD_PeakList + 1);
    EPL++;
  }
}


static void BuildNextPeakMx(void)
{
  int            r_iEPL, c_iEPL, c_nPos;
  T_fVector      *r_cSymEquiv, *c_cSymEquiv;
  T_eD_PeakList  *rEPL, *cEPL;
  Fprec          *NPMx;
  Fprec          dx0, dy0, dz0, dx, dy, dz, Dist2;


  /* compute upper triangle elements of NextPeakMx
     copy to lower triangle
   */

  NPMx = NextPeakMx;
  rEPL = eD_PeakList;

  for (r_iEPL = 0; r_iEPL < NeD_PeakList; r_iEPL++)
  {
    r_cSymEquiv = rEPL->cSymEquiv;
    cEPL = eD_PeakList;

    for (c_iEPL = 0; c_iEPL < NeD_PeakList; c_iEPL++)
    {
      if (r_iEPL > c_iEPL) /* NextPeakMx is symmetric */
        *NPMx = NextPeakMx[c_iEPL * NeD_PeakList + r_iEPL];

      else if (r_iEPL != c_iEPL)
      {
        c_cSymEquiv = cEPL->cSymEquiv;
        c_nPos = cEPL->WL_Entry->nPositions;

        *NPMx = MaxLatticeTr2;

        while (c_nPos--)
        {
#define EvalDist2(dx_, dy_, dz_)\
  {\
    Dist2 = dx_ * dx_ + dy_ * dy_ + dz_ * dz_;\
    if (*NPMx > Dist2) *NPMx = Dist2;\
  }
          dx0 = c_cSymEquiv->x - r_cSymEquiv->x;
          dy0 = c_cSymEquiv->y - r_cSymEquiv->y;
          dz0 = c_cSymEquiv->z - r_cSymEquiv->z;
          EvalDist2(dx0, dy0, dz0);
#include "mindist2.c" /* hand-made inline code */
#undef EvalDist2
          c_cSymEquiv++;
        }
      }

      NPMx++;
      cEPL++;
    }
    rEPL++;
  }


  if (Debug0) /* Debug: Print NextPeakMx */
  {
    int  r, c;

    NPMx = NextPeakMx;

    for (r = 0; r < NeD_PeakList; r++)
    {
      for (c = 0; c < NeD_PeakList; c++)
      {
        Fprintf(stdout, "NPMx[%d][%d] = %8.4f\n", r, c, AppSqrt(*NPMx));
        NPMx++;
      }
    }
  }
}


static void AssignAtoms(void)
{
  int            nPeaksPlus, iPUC, iPUC_Plus, nPos;
  int            iEPL, jEPL, iAT, jAT;
  T_eD_PeakList  *EPL;
  Fprec          *NPMx, *MD2Mx;


  nAsyPeaks = 0;
  nPeaks    = 0;

  for (iAT = 0; iAT < uAtomType; iAT++)
  {
    iPUC = 0;
    EPL = eD_PeakList;

    for (iEPL = 0; iEPL < nCeD_PeakList; iEPL++)
    {
      if (   EPL->iAtomType == -1
          && (   AtomType[iAT].Class != ATC_Node
              || eD_PeakList[iEPL].WL_Entry->CanBeNode))
      {
        nPos = EPL->WL_Entry->nPositions;

            nPeaksPlus = nPeaks + nPos;
        if (nPeaksPlus <= MaxRecycledAtoms)
        {
              iPUC_Plus = iPUC + nPos;
          if (iPUC_Plus <= AtomType[iAT].nPerUnitCell)
          {
            NPMx = NextPeakMx + iEPL * NeD_PeakList;
            MD2Mx = MinDistance2Mx + iAT * uAtomType;

            if (NPMx[iEPL] < MD2Mx[iAT])
            {
              if (NPMx[iEPL] < MinOfMD2Mx)
                EPL->iAtomType = -2;
            }
            else
            {
              for (jEPL = 0; jEPL < nCeD_PeakList; jEPL++)
              {
                    jAT = eD_PeakList[jEPL].iAtomType;
                if (jAT >= 0 &&  NPMx[jEPL] < MD2Mx[jAT])
                {
                  if (NPMx[jEPL] < MinOfMD2Mx)
                    EPL->iAtomType = -2;

                  break;
                }
              }

              if (jEPL == nCeD_PeakList)
              {
                EPL->iAtomType = iAT;
                nAsyPeaks++;
                nPeaks = nPeaksPlus;

                if (nPeaks == MaxRecycledAtoms) return;

                iPUC = iPUC_Plus;
              }
            }
          }
        }
      }

      EPL++;
    }
  }
}


static void AssignLargestFwFragment(void)
{
  int  iAT, jAT, iEPL, nBonds, nBondsFirst;
  int  nPiAT, NewP, Over;


  nAsyPeaks = 0;
  nPeaks    = 0;

  if (LargestFwFragment.nAsyN == 0.)
    return;

  for (iAT = 0; iAT < uAtomType; iAT++)
    if (AtomType[iAT].Class == ATC_Node)
      break;

  if (iAT == uAtomType)
    return;

  if (FwSearchMethod == FwSM_FwTracking)
  {
      jAT = iAT + 1;
    nPiAT = 0;

    for (iEPL =  LargestFwFragment.Low_iEPL;
         iEPL <= LargestFwFragment.Top_iEPL;
         iEPL++)
    {
      if (LargestFwFragment.LnB[iEPL])
      {
        NewP = eD_PeakList[iEPL].WL_Entry->nPositions;
        Over = NewP + nPiAT - AtomType[iAT].nPerUnitCell;

        if (Over > 0)
        {
          for (; jAT < uAtomType; jAT++)
            if (AtomType[jAT].Class == ATC_Node)
              break;

          if (jAT < uAtomType && Over > NewP / 2)
          {
              iAT = jAT;
            nPiAT = 0;
             Over = 0;
          }
        }

        eD_PeakList[iEPL].iAtomType = iAT;
        nAsyPeaks++;
        nPeaks += NewP;
        nPiAT  += NewP;

        if (jAT < uAtomType && Over > 0)
        {
            iAT = jAT;
          nPiAT = 0;
        }
      }
    }
  }
  else /* if (FwSearchMethod == FwSM_AltFwTracking) */
  {
    for (jAT = iAT + 1; jAT < uAtomType; jAT++)
      if (AtomType[jAT].Class == ATC_Node)
        break;

    if (jAT == uAtomType)
      return;

    nBondsFirst = 0;

    for (iEPL =  LargestFwFragment.Low_iEPL;
         iEPL <= LargestFwFragment.Top_iEPL;
         iEPL++)
    {
          nBonds = LargestFwFragment.LnB[iEPL];
      if (nBonds)
      {
        if (nBondsFirst == 0)
            nBondsFirst = nBonds;

        if (   nBonds < 0 && nBondsFirst < 0
            || nBonds > 0 && nBondsFirst > 0)
          eD_PeakList[iEPL].iAtomType = iAT;
        else
          eD_PeakList[iEPL].iAtomType = jAT;

        nAsyPeaks++;
        nPeaks += eD_PeakList[iEPL].WL_Entry->nPositions;
      }
    }
  }
}


static int Internal_eD_PSE = -1;

static int eD_PeakListSortFunction(const T_eD_PeakList *a,
                                   const T_eD_PeakList *b)
{
  switch (Internal_eD_PSE)
  {
    case PSE_Grid_eD:
      if (a->Grid_eD > b->Grid_eD) return -1;
      if (a->Grid_eD < b->Grid_eD) return  1;
      return 0;
    case PSE_Maximum:
      if (a->Maximum > b->Maximum) return -1;
      if (a->Maximum < b->Maximum) return  1;
      return 0;
    case PSE_Integral:
      if (a->Integral > b->Integral) return -1;
      if (a->Integral < b->Integral) return  1;
      return 0;
    case PSE_lLPNNB:
      if (a->lLPNNB > b->lLPNNB) return -1;
      if (a->lLPNNB < b->lLPNNB) return  1;
      if (a->Index < b->Index) return -1;
      if (a->Index > b->Index) return  1;
      return 0;
    case PSE_Index:
      if (a->Index < b->Index) return -1;
      if (a->Index > b->Index) return  1;
      return 0;
    default:
      InternalError("Corrupt Internal_eD_PSE");
      break;
  }
  return 0;
}


void Sort_eD_PeakList(int PeaksSortElement, int nEPL)
{
  Internal_eD_PSE = PeaksSortElement;

  qsort((void *) eD_PeakList, nEPL, sizeof (T_eD_PeakList),
        (SortFunction) eD_PeakListSortFunction);
}


static int Asy_nEPL(int MaxSym)
{
  int  nSymEquiv, nAsy, iEPL;


  nSymEquiv = 0;
  nAsy      = 0;

  for (iEPL = 0; iEPL < NeD_PeakList; iEPL++)
  {
        nSymEquiv += eD_PeakList[iEPL].WL_Entry->nPositions;
    if (nSymEquiv > MaxSym)
      break;

    nAsy++;
  }

  return nAsy;
}


void Free_eD_PeakList(void)
{
  if (eD_PeakList) {
    if (   List_xRawSymEquiv == NULL
        || List_cRawSymEquiv == NULL
        || NextPeakMx == NULL) {
      InternalError(NULL);
    }
    AppFree(NextPeakMx, NeD_PeakList * NeD_PeakList);
            NextPeakMx = NULL;
    AppFree(List_cRawSymEquiv, nList__RawSymEquiv);
            List_cRawSymEquiv = NULL;
    AppFree(List_xRawSymEquiv, nList__RawSymEquiv);
            List_xRawSymEquiv = NULL;
    AppFree(eD_PeakList, NeD_PeakList);
            eD_PeakList = NULL;
    NeD_PeakList = 0;
    nList__RawSymEquiv = 0;
  }
}


int Build_eD_PeakList(Fprec *eDensity, T_PeakFlags *PeakFlags,
                      int UseLargestFwFragment)
{
  int    ix, iy, iz, i;
  int    nSymEquiv, iSymEquiv;
  int    FwFragWasUsed;
  Fprec  Dyn_eD_CutOff;

  Fprec          *eD;
  T_PeakFlags    *PF;
  T_eD_PeakList  *EPL;


  MapInd(Nx * Ny * Nz, eDensity, PeakFlags, 1, 1, NULL);

  Dyn_eD_CutOff = PeakSearch(eDensity, PeakFlags, PeakSearchLevel);

  Free_eD_PeakList();

  nSymEquiv = 0;
  NeD_PeakList = 0;

  PF = PeakFlags;
  eD = eDensity;

  for (i = 0; i < Nx * Ny * Nz; i++, PF++, eD++)
  {
    if (   ((*PF) & PF_HighBit)
        && ((*PF) & PF_IsPeakBit)
        && (*eD) > Dyn_eD_CutOff)
    {
      nSymEquiv += WyckoffList[(*PF) & PF_iWL_Mask].nPositions;
      NeD_PeakList++;
    }
  }

  if (nSymEquiv > MaxRawPeaks)
    InternalError("Corrupt Dyn_eD_CutOff");

  nList__RawSymEquiv = nSymEquiv;
  CheckMalloc(eD_PeakList, NeD_PeakList);
  CheckMalloc(List_xRawSymEquiv, nList__RawSymEquiv);
  CheckMalloc(List_cRawSymEquiv, nList__RawSymEquiv);
  CheckMalloc(NextPeakMx, NeD_PeakList * NeD_PeakList);

  iSymEquiv = 0;
  EPL = eD_PeakList;

  PF = PeakFlags;
  eD = eDensity;

  for (ix = 0; ix < Nx; ix++)
  for (iy = 0; iy < Ny; iy++)
  for (iz = 0; iz < Nz; iz++, PF++, eD++)
  {
    if (   ((*PF) & PF_HighBit)
        && ((*PF) & PF_IsPeakBit)
        && (*eD) > Dyn_eD_CutOff)
    {
      EPL->WL_Entry  = &WyckoffList[(*PF) & PF_iWL_Mask];
      EPL->xSymEquiv = &List_xRawSymEquiv[iSymEquiv];
      EPL->cSymEquiv = &List_cRawSymEquiv[iSymEquiv];
      EPL->Index     = -1;
      EPL->iAtomType = -1;
      EPL->Color     =  0;
      EPL->IsAOS     = -1;
      InterpolatePeakPosition(eDensity, PeakFlags, EPL, ix, iy, iz);

      iSymEquiv += WyckoffList[(*PF) & PF_iWL_Mask].nPositions;
      EPL++;
    }
  }

  if (iSymEquiv != nSymEquiv)
    InternalError(NULL);

  Sort_eD_PeakList(eD_PeaksSortElement, NeD_PeakList);

  for (i = 0; i < NeD_PeakList; i++)
    eD_PeakList[i].Index = i;

  SetSymEquiv(&SpgrInfo, eD_PeakList, NeD_PeakList);
  BuildNextPeakMx();

  LargestFwFragment.nAsyN    = 0;
  LargestFwFragment.nSymN    = 0;
  LargestFwFragment.mBpN     = 0.;
  LargestFwFragment.Low_iEPL = 0;
  LargestFwFragment.Top_iEPL = 0;
  LargestFwFragment.LnB      = NULL;

  if (MaxSymNodes != 0)
  {
    if (UseLargestFwFragment)
    {
      nCeD_PeakList = Asy_nEPL(MaxPeaksFwFragmentSearch);
      CheckMalloc(LargestFwFragment.LnB, nCeD_PeakList);
    }
    else
      nCeD_PeakList = Asy_nEPL(MaxPeaksFwSearch);

    FwSearch(UseLargestFwFragment);
  }

  nAsyPeaks = 0;
  nPeaks    = 0;

  if (UseLargestFwFragment)
    AssignLargestFwFragment();

  if (LargestFwFragment.LnB)
  {
    AppFree(LargestFwFragment.LnB, nCeD_PeakList);
            LargestFwFragment.LnB = NULL;
  }

  FwFragWasUsed = 0;

  if (nAsyPeaks == 0)
  {
    nCeD_PeakList = Asy_nEPL(MaxPotentialAtoms);
    AssignAtoms();
  }
  else
    FwFragWasUsed = 1;

  if (Debug0) /* Debug: Print eD_PeakList */
  {
    Fprintf(stdout, "eD_PeakList\n");
    for (i = 0; i < NeD_PeakList; i++)
    {
      Fprintf(stdout, "Pu%-3d", i);
      Print_eEPL(stdout, eD_PeakList, i, -1, 1);
      putc('\n', stdout);
    }
    putc('\n', stdout);
  }

  return FwFragWasUsed;
}


void SiteFrame(void)
{
  int            MaxSymEquiv, nSymOp, nSymEquiv, iWL, iS, nPos;
  T_fVector      SymEquiv[192];
  T_WyckoffList  Buf_eWL;
  int            WL_Flag[48];
  T_Site         *S;
  Fprec             Dist2ConsiderSame;
  Fprec          MaxDist2ConsideredSame;
  Fprec          MinDist2Distinct;
  T_eD_PeakList  *EPL;


      MaxSymEquiv = SpgrInfo.OrderL;
  if (MaxSymEquiv > sizeof SymEquiv / sizeof (*SymEquiv))
    InternalError("Corrupt SpgrInfo.OrderL");

      nSymOp = SpgrInfo.OrderP;
  if (nSymOp > sizeof WL_Flag / sizeof (*WL_Flag))
    InternalError("Corrupt SpgrInfo.OrderP");

  Buf_eWL.WL_Flag = WL_Flag;

                    /* (must be very low) */
     Dist2ConsiderSame = .001 * .001; /* ARBITRARY */
  MaxDist2ConsideredSame = -1.;
  MinDist2Distinct = MaxLatticeTr2;

  for (iWL = 0; iWL < nWyckoffList; iWL++)
    WyckoffList[iWL].Count = 0;

  NeD_PeakList = nSite;

  CheckMalloc(eD_PeakList, NeD_PeakList);

  MaxRawPeaks = 0;

  for (iS = 0, S = Site; iS < nSite; iS++, S++)
  {
    Buf_eWL.nPositions = CalcSymEquiv(S->x, S->y, S->z,
                                      SymEquiv, MaxSymEquiv,
                                          Dist2ConsiderSame,
                                      &MaxDist2ConsideredSame,
                                      &MinDist2Distinct,
                                      Buf_eWL.WL_Flag, nSymOp);

        iWL = FindWL_Entry(WyckoffList, &Buf_eWL);
    if (iWL < 0)
      InternalError("Corrupt iWL");

    eD_PeakList[iS].WL_Entry = &WyckoffList[iWL];

    MaxRawPeaks += Buf_eWL.nPositions;
  }

  CheckMalloc(List_xRawSymEquiv, MaxRawPeaks);
  CheckMalloc(List_cRawSymEquiv, MaxRawPeaks);
  CheckMalloc(NextPeakMx, NeD_PeakList * NeD_PeakList);

  nSymEquiv = 0;
  EPL = eD_PeakList;

  for (iS = 0, S = Site; iS < nSite; iS++, S++)
  {
    nPos = EPL->WL_Entry->nPositions;

    if (nSymEquiv + nPos > MaxRawPeaks)
      InternalError("Corrupt MaxRawPeaks");

    EPL->xSymEquiv  = &List_xRawSymEquiv[nSymEquiv];
    EPL->cSymEquiv  = &List_cRawSymEquiv[nSymEquiv];
    EPL->Index      = iS;
    EPL->Color      =  0;
    EPL->IsAOS      = -1;
    EPL->iAtomType  = -1;
    EPL->Position.x = S->x; NormOf(EPL->Position.x);
    EPL->Position.y = S->y; NormOf(EPL->Position.y);
    EPL->Position.z = S->z; NormOf(EPL->Position.z);
    EPL->Grid_eD    = 0.;
    EPL->Maximum    = 0.;
    EPL->Integral   = 0.;
    EPL->nPfI       = -1;

    nSymEquiv += nPos;
    EPL++;
  }

  SetSymEquiv(&SpgrInfo, eD_PeakList, NeD_PeakList);

  if (F_nAllAsyPoints >= 0) {
    PrintTicks(stdout, NULL,
      "# Time SiteFrame()->BuildNextPeakMx() Start  ", "\n");
    Fflush(stdout);
  }

  BuildNextPeakMx();

  if (F_nAllAsyPoints >= 0) {
    PrintTicks(stdout, NULL,
      "# Time SiteFrame()->BuildNextPeakMx() End  ", "\n\n");
    Fflush(stdout);
  }

  if (Debug1) /* Debug: Print eD_PeakList */
  {
    Fprintf(stdout, ">Begin SiteFrame_eD_PeakList\n");

    for (iS = 0; iS < NeD_PeakList; iS++)
    {
      Fprintf(stdout, "%-6s", Site[iS].Label);
      Print_eEPL(stdout, eD_PeakList, iS, -1, 1);
      putc('\n', stdout);
    }

    Fprintf(stdout, ">End SiteFrame_eD_PeakList\n\n");
    Fflush(stdout);
  }

  if (MaxSymNodes == 0)
  {
    for (iS = 0; iS < NeD_PeakList; iS++)
      MaxSymNodes += eD_PeakList[iS].WL_Entry->nPositions;
  }

  if (MaxSymNodes != 0)
  {
    MaxPeaksFwSearch =
    MaxPeaksFwFragmentSearch = MaxRawPeaks;

    nCeD_PeakList = NeD_PeakList;

    FwSearch(F_ShowLargestFwFragment);
  }

  AppFree(NextPeakMx, NeD_PeakList * NeD_PeakList);
  AppFree(List_cRawSymEquiv, MaxRawPeaks);
  AppFree(List_xRawSymEquiv, MaxRawPeaks);
  AppFree(eD_PeakList, NeD_PeakList);

  eD_PeakList = NULL;
  List_xRawSymEquiv = NULL;
  List_cRawSymEquiv = NULL;
  NextPeakMx = NULL;
}


#include "ranmar.h"


void RandomSites(int NumberToGenerate)
{
  int            MaxSymEquiv, iS, f;
  long           l;
  T_eD_PeakList  BufEPL[1];
  Fprec          BufNPMx[1];
  long           CII, CIP, CPI, CPP;
  Fprec          px, py, pz;


  eD_PeakList = BufEPL;
  NextPeakMx  = BufNPMx;

  MaxSymEquiv = SpgrInfo.OrderL;

  CheckMalloc(List_xRawSymEquiv, MaxSymEquiv);
  CheckMalloc(List_cRawSymEquiv, MaxSymEquiv);

  eD_PeakList->xSymEquiv  = List_xRawSymEquiv;
  eD_PeakList->cSymEquiv  = List_cRawSymEquiv;
  eD_PeakList->Index      =  0;
  eD_PeakList->Color      =  0;
  eD_PeakList->IsAOS      = -1;
  eD_PeakList->iAtomType  = -1;
  eD_PeakList->Grid_eD    = 0.;
  eD_PeakList->Maximum    = 0.;
  eD_PeakList->Integral   = 0.;
  eD_PeakList->nPfI       = -1;

  if (nWyckoffList < 1)
    InternalError("Corrupt nWyckoffList");

  DoRandomInitialization();

  for (iS = 0; NumberToGenerate == 0 || iS < NumberToGenerate; iS++)
  {
    if (QuitProgram)
      break;

    if (iS && iS % 5000 == 0)
      Fprintf(stdout, "# %d\n", iS);

    eD_PeakList->Position.x = ranmar();
    eD_PeakList->Position.y = ranmar();
    eD_PeakList->Position.z = ranmar();

    eD_PeakList->WL_Entry = &WyckoffList[nWyckoffList - 1];

    CII = CountInterSectionII;
    CIP = CountInterSectionIP;
    CPI = CountInterSectionPI;
    CPP = CountInterSectionPP;

    px = eD_PeakList->Position.x;
    py = eD_PeakList->Position.y;
    pz = eD_PeakList->Position.z;

    SetSymEquiv(&SpgrInfo, eD_PeakList, 1);

    f = 0;

    for (l = CII; l != CountInterSectionII; l++) { Fputs("II ", stdout); f++; }
    for (l = CIP; l != CountInterSectionIP; l++) { Fputs("IP ", stdout); f++; }
    for (l = CPI; l != CountInterSectionPI; l++) { Fputs("PI ", stdout); f++; }
    for (l = CPP; l != CountInterSectionPP; l++) { Fputs("PP ", stdout); f++; }

    if (f)
    {
      Fprintf(stdout, "%d> %7.5f %7.5f %7.5f\n", iS, px, py, pz);

      for (l = CII; l != CountInterSectionII; l++) Fputs("II ", stdout);
      for (l = CIP; l != CountInterSectionIP; l++) Fputs("IP ", stdout);
      for (l = CPI; l != CountInterSectionPI; l++) Fputs("PI ", stdout);
      for (l = CPP; l != CountInterSectionPP; l++) Fputs("PP ", stdout);

      Fprintf(stdout, "%d< %7.5f %7.5f %7.5f\n", iS,
        eD_PeakList->Position.x,
        eD_PeakList->Position.y,
        eD_PeakList->Position.z);

      putc('\n', stdout);
      Fflush(stdout);
    }
  }

  eD_PeakList = NULL;
  NextPeakMx  = NULL;
  AppFree(List_xRawSymEquiv, MaxSymEquiv);
  AppFree(List_cRawSymEquiv, MaxSymEquiv);
}


static void SnapSpecial(T_eD_PeakList *EPL,
                        T_FreeVect *CurFV,
                        int iSymOp)
{
  int              nLoopInv, nTrV;
  int              iList, iLoopInv, iTrV;
  int               ShortiTrV;
  const T_iVector  *Shortm0pSh, *m0pShift;
  Fprec             ShortDist2, Dist2;
  T_fVector        cPos[1], RT, *fTrV, cSymEquiv[1], xSymEquiv[1];
  Fprec            dx0, dy0, dz0, dx, dy, dz;
  T_iVector        UC_Shift0, UC_Shift, ShortUC_Shift;

#include "mindist2.h"


  Tform_xc(EPL->Position.x, EPL->Position.y, EPL->Position.z,
           cPos->x, cPos->y, cPos->z);

  UC_Shift0.x =
  UC_Shift0.y =
  UC_Shift0.z = 0;

  nLoopInv = Sg_nLoopInv(&SpgrInfo);
  nTrV     = SpgrInfo.LatticeInfo->nTrVector;

  iList    = iSymOp / nLoopInv;
  iLoopInv = iSymOp % nLoopInv;

  ShortUC_Shift.x =
  ShortUC_Shift.y =
  ShortUC_Shift.z = 0;
  ShortiTrV       = -1;
  Shortm0pSh      = NULL;
  ShortDist2      = MaxLatticeTr2;

#define EvalDist2(dx_, dy_, dz_)\
  {\
    Dist2 = dx_ * dx_ + dy_ * dy_ + dz_ * dz_;\
    if (ShortDist2 > Dist2)\
    {\
      ShortDist2      = Dist2;\
      Shortm0pSh      = m0pShift;\
      ShortiTrV       = iTrV;\
      ShortUC_Shift.x = UC_Shift.x;\
      ShortUC_Shift.y = UC_Shift.y;\
      ShortUC_Shift.z = UC_Shift.z;\
    }\
    m0pShift++;\
  }

  fRTMx_t_fVector(&RT, &List_fSeitzMx[iList], &EPL->Position);

  for (iTrV = 0, fTrV = fTrVector; iTrV < nTrV; iTrV++, fTrV++)
  {
    if (iLoopInv == 0)
    {
      xSymEquiv->x =  RT.x + fTrV->x;
      xSymEquiv->y =  RT.y + fTrV->y;
      xSymEquiv->z =  RT.z + fTrV->z;
    }
    else
    {
      xSymEquiv->x = -RT.x + fTrV->x;
      xSymEquiv->y = -RT.y + fTrV->y;
      xSymEquiv->z = -RT.z + fTrV->z;
    }

    CountNormOf(xSymEquiv->x, UC_Shift.x);
    CountNormOf(xSymEquiv->y, UC_Shift.y);
    CountNormOf(xSymEquiv->z, UC_Shift.z);

    Tform_xc(xSymEquiv->x, xSymEquiv->y, xSymEquiv->z,
             cSymEquiv->x, cSymEquiv->y, cSymEquiv->z);

    if (Debug0) /* Debug: Print xSymEquiv cSymEquiv */
    {
      T_iVector  *UCS = &UC_Shift;
      Fprintf(stdout,
      " %3d %2d  %8.4f %8.4f %8.4f  %8.4f %8.4f %8.4f  %2d %2d %2d\n",
        iTrV, iSymOp,
        xSymEquiv->x, xSymEquiv->y, xSymEquiv->z,
        cSymEquiv->x, cSymEquiv->y, cSymEquiv->z,
        UCS->x / STBF, UCS->y / STBF, UCS->z / STBF);
    }

    m0pShift = List_m0pShift;

    dx0 = cSymEquiv->x - cPos->x;
    dy0 = cSymEquiv->y - cPos->y;
    dz0 = cSymEquiv->z - cPos->z;
    EvalDist2(dx0, dy0, dz0);

#include "mindist2.c" /* hand-made inline code */
  }

  if (ShortiTrV < 0)
    InternalError(NULL);

  (void) MoreSpecial(EPL, CurFV, iSymOp, iList, iLoopInv, ShortiTrV,
                     &UC_Shift0, &ShortUC_Shift, Shortm0pSh);
#undef EvalDist2
}


void TestMoreSpecial(int NumberToGenerate)
{
  int            iS, iWL;
  int            nSymOp, iSymOp, jSymOp;
  T_eD_PeakList  EPLi[1], EPLj[1];
  T_FreeVect     CurFVi[1], CurFVj[1];
  T_fVector      Pos;


  EPLi->xSymEquiv  = NULL;
  EPLi->cSymEquiv  = NULL;
  EPLi->Index      =  0;
  EPLi->Color      =  0;
  EPLi->IsAOS      = -1;
  EPLi->iAtomType  = -1;
  EPLi->Grid_eD    = 0.;
  EPLi->Maximum    = 0.;
  EPLi->Integral   = 0.;
  EPLi->nPfI       = -1;

  nSymOp   = SpgrInfo.OrderP;

  DoRandomInitialization();

  for (iS = 0; NumberToGenerate == 0 || iS < NumberToGenerate; iS++)
  {
    Pos.x = ranmar();
    Pos.y = ranmar();
    Pos.z = ranmar();

    for (iWL = 0; iWL < nWyckoffList; iWL++)
    {
      for (iSymOp = 1; iSymOp < nSymOp; iSymOp++)
      {
        if (WyckoffList[iWL].WL_Flag[iSymOp] != WL_F_Identity)
          continue;

        EPLi->Position.x = Pos.x;
        EPLi->Position.y = Pos.y;
        EPLi->Position.z = Pos.z;
        EPLi->WL_Entry   = &WyckoffList[iWL];

        CurFVi->V[0] =
        CurFVi->V[1] =
        CurFVi->V[2] = 0.;
        CurFVi->D    = 3;

        if (Debug0) /* Print SnapSpecial iSymOp */
          Fprintf(stdout, "SnapSpecial %d\n", iSymOp);

        SnapSpecial(EPLi, CurFVi, iSymOp);

        if (CurFVi->D)
        {
          for (jSymOp = 1; jSymOp < nSymOp; jSymOp++)
          {
            if (WyckoffList[iWL].WL_Flag[jSymOp] != WL_F_Identity)
              continue;

            (void) memcpy(EPLj, EPLi, sizeof (*EPLj));
            (void) memcpy(CurFVj, CurFVi, sizeof (*CurFVj));

            if (Debug0) /* Print SnapSpecial iSymOp jSymOp */
              Fprintf(stdout, "SnapSpecial %d %d\n", iSymOp, jSymOp);

            SnapSpecial(EPLj, CurFVj, jSymOp);
          }
        }
      }
    }
  }
}
