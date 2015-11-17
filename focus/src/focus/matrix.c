#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#ifdef STAND_ALONE
#define A_AppGlobal__
#endif


#include "main.h"


#define MaxMxInverse_n  3
#define MaxMxLU_n      10


#define SWAP(a, b) { dum = (a); (a) = (b); (b) = dum; }
#define a(i, j) (a + (i) * n)[j]
#define b(i, j) (b + (i) * n)[j]


int MxInverse(Fprec *a, int n, Fprec *b, int m)
{
  int         i, icol, irow, j, k, l, ll;
  Fprec       big, dum, pivinv;

  int  indxc[MaxMxInverse_n];
  int  indxr[MaxMxInverse_n];
  int   ipiv[MaxMxInverse_n];


  if (n > MaxMxInverse_n)
    InternalError("n > MaxMxInverse_n");

  for (j = 0; j < n; j++) ipiv[j] = 0;

  for (i = 0; i < n; i++)
  {
    big = 0.0;

    for (j = 0; j < n; j++)
    {
      if (ipiv[j] != 1)
      {
        for (k = 0; k < n; k++)
        {
          if (ipiv[k] == 0)
          {
            dum = AppFabs(a(j, k));

            if (dum >= big)
            {
              big = dum;
              irow = j;
              icol = k;
            }
          }
          else if (ipiv[k] > 1)
            return -1;
        }
      }
    }

    (ipiv[icol])++;

    if (irow != icol)
    {
      for (l = 0; l < n; l++) SWAP(a(irow, l), a(icol, l))
      for (l = 0; l < m; l++) SWAP(b(l, irow), b(l, icol))
    }

    indxr[i] = irow;
    indxc[i] = icol;

    if (a(icol, icol) == (Fprec) 0.)
      return -2;

    pivinv = (Fprec) 1. / a(icol, icol);
    a(icol, icol) = 1.;

    for (l = 0; l < n; l++) a(icol, l) *= pivinv;
    for (l = 0; l < m; l++) b(l, icol) *= pivinv;

    for (ll = 0; ll < n; ll++)
    {
      if (ll != icol)
      {
        dum = a(ll, icol);
        a(ll, icol) = 0.;
        for (l = 0; l < n; l++) a(ll, l) -= a(icol, l) * dum;
        for (l = 0; l < m; l++) b(l, ll) -= b(l, icol) * dum;
      }
    }
  }

  l = n;
  while (l--)
  {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
        SWAP(a(k, indxr[l]), a(k, indxc[l]));
  }

  return 0;
}

#undef b


int MxLU_Decomposition(Fprec *a, int n, int *perm)
{
  int    imax, i, j, k;
  Fprec  big, dum, sum;

  Fprec  LU_vv[MaxMxLU_n];


  if (n > MaxMxLU_n)
    InternalError("n > MaxMxLU_n");

  for (i = 0; i < n; i++)
  {
    big = 0.;
    for (j = 0; j < n; j++)
    {
      dum = AppFabs(a(i, j));
      if (dum > big) big = dum;
    }
    if (big == (Fprec) 0.) return -1; /* singular matrix */
    LU_vv[i] = (Fprec) 1. / big;
  }

  for (j = 0; j < n; j++)
  {
    for (i = 0; i < j; i++)
    {
      sum = a(i, j);
      for (k = 0; k < i; k++) sum -= a(i, k) * a(k, j);
      a(i, j) = sum;
    }

    big = 0.;
    for ( ; i < n; i++)
    {
      sum = a(i, j);
      for (k = 0; k < j; k++) sum -= a(i, k) * a(k, j);
      a(i, j) = sum;

      dum = LU_vv[i] * AppFabs(sum);
      if (dum >= big)
      {
        big = dum;
        imax = i;
      }
    }

    if (j != imax)
    {
      for (k = 0; k < n; k++) SWAP(a(imax, k), a(j, k));
      LU_vv[imax] = LU_vv[j]; /* no swap, we don't need LU_vv[j] any more */
    }

    perm[j] = imax;

    if (a(j, j) == (Fprec) 0.) return -2; /* singular matrix */

        i = j + 1;
    if (i < n)
    {
      dum = (Fprec) 1. / a(j, j);
      for ( ; i < n; i++)
        a(i, j) *= dum;
    }
  }

  return 0;
}

#undef SWAP


void MxLU_BackSubstitution(Fprec *a, int n, int *perm, Fprec *b)
{
  int    ii, ip, i, j;
  Fprec  sum;


  ii = -1;
  for (i = 0; i < n; i++)
  {
    ip = perm[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii >= 0)
      for (j = ii; j <= i - 1; j++) sum -= a(i, j) * b[j];
    else if (sum)
      ii = i;
    b[i] = sum;
  }

  for (i = n; --i >= 0; )
  {
    sum = b[i];
    for (j = i + 1; j < n; j++) sum -= a(i, j) * b[j];
    b[i] = sum / a(i, i);
  }
}

#undef a


void MxMultiply(Fprec *ab, Fprec *a, Fprec *b, int ma, int na, int nb)
{
  int         i, j, k;
  Fprec       *ai, *aij, *bk, *bkj;

  ai = a;

  for (i = 0; i < ma; i++)
  {
    bk = b;

    for (k = 0; k < nb; k++)
    {
      aij = ai;
      bkj = bk;

      *ab = 0.;

      for (j = 0; j < na; j++)
      {
        *ab += (*aij) * (*bkj);

        aij++;
        bkj += nb;
      }

      ab++;
      bk++;
    }

    ai += na;
  }
}


void MxTranspose(Fprec *a, int ma, int na, Fprec *t)
{
  int         i, j;
  Fprec       *tji;


  for (i = 0; i < ma; i++)
  {
    tji = t;

    for (j = 0; j < na; j++)
    {
      *tji = *a++;
       tji += ma;
    }

    t++;
  }
}


void MxDump(Fprec *a, int ma, int na, char *note)
{
  int       n;

  if (note) Fprintf(stdout, "%s\n", note);

  while (ma--)
  {
    n = na;

    while (n--)
    {
      Fprintf(stdout, "  %5.2f", *a++);
    }
    Fprintf(stdout, "\n");
  }

  Fflush(stdout);
}



Fprec deter33fMx(const Fprec *Mx)
{
  Fprec  det;

  det =  Mx[0] * (Mx[4] * Mx[8] - Mx[5] * Mx[7]);
  det -= Mx[1] * (Mx[3] * Mx[8] - Mx[5] * Mx[6]);
  det += Mx[2] * (Mx[3] * Mx[7] - Mx[4] * Mx[6]);

  return det;
}


void Inverse33fMx(const Fprec *Mx, Fprec *InvMx)
{
  InvMx[0] =   Mx[4] * Mx[8] - Mx[5] * Mx[7];
  InvMx[1] = - Mx[1] * Mx[8] + Mx[2] * Mx[7];
  InvMx[2] =   Mx[1] * Mx[5] - Mx[2] * Mx[4];
  InvMx[3] = - Mx[3] * Mx[8] + Mx[5] * Mx[6];
  InvMx[4] =   Mx[0] * Mx[8] - Mx[2] * Mx[6];
  InvMx[5] = - Mx[0] * Mx[5] + Mx[2] * Mx[3];
  InvMx[6] =   Mx[3] * Mx[7] - Mx[4] * Mx[6];
  InvMx[7] = - Mx[0] * Mx[7] + Mx[1] * Mx[6];
  InvMx[8] =   Mx[0] * Mx[4] - Mx[1] * Mx[3];
}


#ifdef STAND_ALONE


void I_nternalError(const char *message, const char *file, const int line)
{
  Fflush(stdout);
  Fprintf(stderr, "\n%s(%s:%d): Internal Error", progn, file, line);
  if (message) Fprintf(stderr, ": %s", message);
  putc('\n', stderr);
  exit(3);
}


int main(void)
{
  int    Cycle, i;
  int    perm[3];
  Fprec  a[3 * 3], b[3];


  Debug0 = 0;
  Debug1 = 1;


  for (Cycle = 0; Cycle < 10; Cycle++)
  {
    for (i = 0; i < 9; i++) a[i] = (Fprec)(i % 4 ? i + Cycle : i - Cycle);
    MxDump(a, 3, 3, "Original a");

        i = MxLU_Decomposition(a, 3, perm);
    if (i != 0)
      Fprintf(stdout, "MxLU_Decomposition() returns %d\n", i);
    else
    {
      MxDump(a, 3, 3, "LU");

      for (i = 0; i < 3; i++) b[i] = (Fprec)(i % 2 ? i + Cycle : i - Cycle);
      MxDump(b, 3, 1, "Original b");

      MxLU_BackSubstitution(a, 3, perm, b);
      MxDump(b, 3, 1, "x: A.x = b");
    }


    for (i = 0; i < 9; i++) a[i] = (Fprec)(i % 4 ? i + Cycle : i - Cycle);
    MxDump(a, 3, 3, "Original a");

    for (i = 0; i < 3; i++) b[i] = (Fprec)(i % 2 ? i + Cycle : i - Cycle);
    MxDump(b, 3, 1, "Original b");

        i = MxInverse(a, 3, b, 1);
    if (i == 0)
      MxDump(b, 3, 1, "x: A.x = b");
    else
      Fprintf(stdout, "MxInverse() returns %d\n", i);

    putc('\n', stdout);
    putc('\n', stdout);
  }

  exit(0);
  return 0;
}


#endif  /* STAND_ALONE */
