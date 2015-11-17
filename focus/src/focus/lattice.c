#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef STAND_ALONE
#define A_AppGlobal__
#endif

#include "main.h"
#include "lattice.h"
#include "matrix.h"


double Q_hkl(int h, int k, int l, T_LatticeConstants *rlc)
{
  double    Q;

  Q =   (int)(h * h) * (rlc->a * rlc->a)
      + (int)(k * k) * (rlc->b * rlc->b)
      + (int)(l * l) * (rlc->c * rlc->c)
      + (int)(2 * h * k) * (rlc->a * rlc->b * rlc->cg)
      + (int)(2 * h * l) * (rlc->a * rlc->c * rlc->cb)
      + (int)(2 * k * l) * (rlc->b * rlc->c * rlc->ca);

  if (Q <= 0.) Q = 0.;

  return Q;
}


#define EpsPI (1.e-6) /* ARBITRARY */

static double sinC(double arg)
{
  if (M_PI_2 - EpsPI <= arg && arg <= M_PI_2 + EpsPI)
    return 1.;

  return sin(arg);
}


static double cosC(double arg)
{
  if (M_PI_2 - EpsPI <= arg && arg <= M_PI_2 + EpsPI)
    return 0.;

  return cos(arg);
}


int Lc2RLc(T_LatticeConstants *lc, T_LatticeConstants *rlc)
{
  /* Transformation Lattice Constants -> Reciprocal Lattice Constants
     after Kleber, W., 17. Aufl., Verlag Technik GmbH Berlin 1990, P.352

     V = a * b * c * sqrt(1 - cos(alpha)^2 - cos(beta)^2 - cos(gamma)^2
                            + 2 * cos(alpha) * cos(beta) * cos(gamma));
   */

  double  D;


  lc->sa = sinC(lc->alpha);
  lc->sb = sinC(lc->beta);
  lc->sg = sinC(lc->gamma);
  lc->ca = cosC(lc->alpha);
  lc->cb = cosC(lc->beta);
  lc->cg = cosC(lc->gamma);

  D = 1. - lc->ca * lc->ca - lc->cb * lc->cb - lc->cg * lc->cg
         + 2. * lc->ca * lc->cb * lc->cg;
  if (D < 0.) return -1;

  lc->v = lc->a * lc->b * lc->c * sqrt(D);
  if (lc->v == 0.) return -1;

  if (lc->sa == 0. || lc->sb == 0. || lc->sg == 0.) return -1;

  if (rlc != NULL)
  {
    rlc->a = lc->b * lc->c * lc->sa / lc->v;
    rlc->b = lc->c * lc->a * lc->sb / lc->v;
    rlc->c = lc->a * lc->b * lc->sg / lc->v;
    rlc->ca = (lc->cb * lc->cg - lc->ca) / (lc->sb * lc->sg);
    rlc->cb = (lc->cg * lc->ca - lc->cb) / (lc->sg * lc->sa);
    rlc->cg = (lc->ca * lc->cb - lc->cg) / (lc->sa * lc->sb);
    rlc->alpha = acos(rlc->ca);
    rlc->beta  = acos(rlc->cb);
    rlc->gamma = acos(rlc->cg);
    rlc->sa = sinC(rlc->alpha);
    rlc->sb = sinC(rlc->beta);
    rlc->sg = sinC(rlc->gamma);
    rlc->v = 1. / lc->v;
  }

  return 0;
}


int SetupTrCrystCarte(T_LatticeConstants *lc, T_LatticeConstants *rlc,
                      int Convention,
                      Fprec *TM, int FlagInverse)
{
  /*
    Crystallographic Basis: D = {a,b,c}
    Cartesian Basis:        C = {i,j,k}

    Boisen & Gibbs: Mathematical Crystallography; Reviews in Mineralogy,
    Volume 15, Revised Edition; Mineralogical Society of America;
    Washington D.C. 1990, pp. 72-83.

    k || c
    j || c x a
    i || j x k

    [r]D = position in fractional coordinates
    [r]C = position in cartesian coordinates

    [r]C = A . [r]D
    [r]D = inverse(A) . [r]C

        [ a.sin(beta)  -b.sin(alpha).cos(gamma*)    0 ]
        [                                             ]
    A = [     0         b.sin(alpha).sin(gamma*)    0 ]
        [                                             ]
        [ a.cos(beta)         b.cos(alpha)          c ]

          [ a sb  - b sa cg  0 ]
          [                    ]
        = [   0    b sa sg   0 ]
          [                    ]
          [ a cb     b ca    c ]

    determinant(A) = volume of the unit cell

                 [    1             cg             ]
                 [  ----         -------        0  ]
                 [  a sb         a sb sg           ]
                 [                                 ]
                 [                  1              ]
    inverse(A) = [    0          -------        0  ]
                 [               b sa sg           ]
                 [                                 ]
                 [    cb     sb ca + sa cg cb      ]
                 [ - ----  - ----------------  1/c ]
                 [   sb c       sb sa sg c         ]


    PDB convention:
      i || a
      j is in (a,b) plane
      k = i x j
   */

  T_LatticeConstants  Buf_lc;
  Fprec               s1rca2;


  if      ( lc == NULL)
  {
    if (rlc == NULL) return -1;
    if (Lc2RLc(rlc, &Buf_lc) != 0) return -1;
    lc = &Buf_lc;
  }
  else if (rlc == NULL)
  {
    if ( lc == NULL) return -1;
    if (Lc2RLc( lc, &Buf_lc) != 0) return -1;
    rlc = &Buf_lc;
  }

  if (Convention == 0) /* Boisen & Gibbs */
  {
    if (FlagInverse == 0)
    {
      /* fractional to cartesian */
      TM[0] = lc->a * lc->sb;
      TM[1] = -lc->b * lc->sa * rlc->cg;
      TM[2] = 0.;
      TM[3] = 0.;
      TM[4] = lc->b * lc->sa * rlc->sg;
      TM[5] = 0.;
      TM[6] = lc->a * lc->cb;
      TM[7] = lc->b * lc->ca;
      TM[8] = lc->c;
    }
    else
    {
      if (   lc->a == 0. || lc->b == 0. || lc->c == 0.
          || lc->sa == 0. || lc->sb == 0. || rlc->sg == 0.)
        return -1;

      /* cartesian to fractional */
      TM[0] = (Fprec) 1. / (lc->a * lc->sb);
      TM[1] = rlc->cg / (lc->a * lc->sb * rlc->sg);
      TM[2] = 0.;
      TM[3] = 0.;
      TM[4] = (Fprec) 1. / (lc->b * lc->sa * rlc->sg);
      TM[5] = 0.;
      TM[6] = -lc->cb / (lc->sb * lc->c);
      TM[7] = -(lc->sb * lc->ca + lc->sa * rlc->cg * lc->cb)
              / (lc->sb * lc->sa * rlc->sg * lc->c);
      TM[8] = (Fprec) 1. / lc->c;
    }
  }
  else /* PDB */
  {
    s1rca2 = sqrt(1 - rlc->ca * rlc->ca);

    if (FlagInverse == 0)
    {
      /* fractional to cartesian */
      TM[0] = lc->a;
      TM[1] = lc->cg * lc->b;
      TM[2] = lc->cb * lc->c;
      TM[3] = 0.;
      TM[4] = lc->sg * lc->b;
      TM[5] = -lc->sb * rlc->ca * lc->c;
      TM[6] = 0.;
      TM[7] = 0.;
      TM[8] = lc->sb * lc->c * s1rca2;
    }
    else
    {
      if (   lc->a == 0. || lc->b == 0. || lc->c == 0.
          || lc->sb == 0. || lc->sg == 0. || s1rca2 == 0.)
        return -1;

      /* cartesian to fractional */
      TM[0] = 1. / lc->a;
      TM[1] = -lc->cg / (lc->sg * lc->a);
      TM[2] =   -(lc->cg * lc->sb * rlc->ca + lc->cb * lc->sg)
              / (lc->sb * s1rca2 * lc->sg * lc->a);
      TM[3] = 0.;
      TM[4] = 1. / (lc->sg * lc->b);
      TM[5] = rlc->ca / (s1rca2 * lc->sg * lc->b);
      TM[6] = 0.;
      TM[7] = 0.;
      TM[8] = 1. / (lc->sb * s1rca2 * lc->c);
    }
  }

  return 0;
}


void Lc2MetricalMx(T_LatticeConstants *lc, Fprec *G)
{
  G[0] =        lc->a * lc->a;
  G[1] = G[3] = lc->a * lc->b * lc->cg;
  G[2] = G[6] = lc->a * lc->c * lc->cb;

  G[4] =        lc->b * lc->b;
  G[5] = G[7] = lc->b * lc->c * lc->ca;

  G[8] =        lc->c * lc->c;
}


Fprec DotG(Fprec *u, Fprec *G, Fprec *v)
{
  Fprec      uGv;

  uGv =  u[0] * (G[0] * v[0] + G[1] * v[1] + G[2] * v[2]);
  uGv += u[1] * (G[3] * v[0] + G[4] * v[1] + G[5] * v[2]);
  uGv += u[2] * (G[6] * v[0] + G[7] * v[1] + G[8] * v[2]);

  return uGv;
}


#ifndef STAND_ALONE

void CrossG(Fprec sqrtdetG, Fprec *G, Fprec *r, Fprec *s, Fprec *rxs)
{
  Fprec      Gr[3], Gs[3];

  MxMultiply(Gr, G, r, 3, 3, 1);
  MxMultiply(Gs, G, s, 3, 3, 1);

  rxs[0] = sqrtdetG * (Gr[1] * Gs[2] - Gs[1] * Gr[2]);
  rxs[1] = sqrtdetG * (Gr[2] * Gs[0] - Gs[2] * Gr[0]);
  rxs[2] = sqrtdetG * (Gr[0] * Gs[1] - Gs[0] * Gr[1]);
}


void CalcMaxhkl(T_LatticeConstants *RLc, Fprec *dmin, Fprec Lambda,
                int *h, int *k, int *l)
{
  int      i;
  Fprec    G[9], u[3], v[3], uxv[3], uxv2;
  Fprec    rSphere, foo;
  int      Maxhkl[3];


  /* calculate radius of sphere of reflection */

  if (Lambda <= 0.)
  {
    if (*dmin <= 0.)
    {
      rSphere = 1.;
      *dmin = 1.;
    }
    else
      rSphere = 1. / *dmin;
  }
  else
  {
    rSphere = 2. / Lambda;

    if (*dmin > 0.)
    {
      foo = 1. / *dmin;
      if (foo <= rSphere)
        rSphere = foo;
      else
        *dmin = 1. / rSphere;
    }
    else
      *dmin = 1. / rSphere;
  }

  /* calculation of maximum hkl */

  Lc2MetricalMx(RLc, G);

  for (i = 0; i < 3; i++)
  {
    u[0] = (i == 2 ? 1 : 0);
    u[1] = (i == 0 ? 1 : 0);
    u[2] = (i == 1 ? 1 : 0);

    v[0] = (i == 1 ? 1 : 0);
    v[1] = (i == 2 ? 1 : 0);
    v[2] = (i == 0 ? 1 : 0);

    CrossG(1., G, u, v, uxv); /* Since lenght of uxv is not used
                                 sqrt(det(G)) = v is set to 1 */
    uxv2 = DotG(uxv, G, uxv);
    Maxhkl[i] =  (int) (1.e-6 + uxv[i] / AppSqrt(uxv2) * rSphere);
                     /* ARBITRARY */
  }

  *h = Maxhkl[0];
  *k = Maxhkl[1];
  *l = Maxhkl[2];
}

#endif  /* ! STAND_ALONE */


#ifdef STAND_ALONE


static void usage(void)
{
  fflush(stdout);
  fprintf(stderr, "usage: %s a b c alpha beta gamma\n", progn);
  exit(1);
}


int main(int argc, const char *argv[])
{
  int                 iarg, i;
  double              x;
  char                extrac;
  T_LatticeConstants  Lc[1], RLc[1];
  Fprec               fcTrMx[9], cfTrMx[9];


  progn = "lattice";

  if (argc != 7) usage();

  for (iarg = 1 ; iarg <= 6; iarg++)
  {
    if (sscanf(argv[iarg], "%lf %c", &x, &extrac) != 1 || x <= 0.) usage();

    switch (iarg)
    {
      case 1: Lc->a = x; break;
      case 2: Lc->b = x; break;
      case 3: Lc->c = x; break;
      case 4: Lc->alpha = x * PIover180; break;
      case 5: Lc->beta  = x * PIover180; break;
      case 6: Lc->gamma = x * PIover180; break;
    }
  }

  if (   Lc2RLc(Lc, RLc) != 0
      || SetupTrCrystCarte(Lc, RLc, 1, fcTrMx, 0) != 0
      || SetupTrCrystCarte(Lc, RLc, 1, cfTrMx, 1) != 0) {
    fprintf(stderr, "%s: corrupt lattice parameters\n", progn);
    exit(1);
  }

  printf("D-space unit cell: %.6g %.6g %.6g %.6g %.6g %.6g\n",
    Lc->a, Lc->b, Lc->c,
    Lc->alpha / PIover180, Lc->beta / PIover180, Lc->gamma / PIover180);

  printf("R-space unit cell: %.6g %.6g %.6g %.6g %.6g %.6g\n",
    RLc->a, RLc->b, RLc->c,
    RLc->alpha / PIover180, RLc->beta / PIover180, RLc->gamma / PIover180);

  printf("frac->cart = {");
  for (i = 0; i < 9; i++) {
    if (i % 3 == 0) putchar('{');
    printf("%.6g", fcTrMx[i]);
    if (i % 3 == 2) putchar('}');
    if (i < 8) putchar(',');
  }
  putchar('}');
  putchar('\n');

  printf("cart->frac = {");
  for (i = 0; i < 9; i++) {
    if (i % 3 == 0) putchar('{');
    printf("%.6g", cfTrMx[i]);
    if (i % 3 == 2) putchar('}');
    if (i < 8) putchar(',');
  }
  putchar('}');
  putchar('\n');

  exit(0);
  return 0;
}


#endif  /* STAND_ALONE */
