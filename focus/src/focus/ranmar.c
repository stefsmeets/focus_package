#if defined(STAND_ALONE_1) || defined(STAND_ALONE_2)
#define A_AppGlobal__
#endif

#include "ranmar.h"

static double u[98], c, cd, cm;
static int i97, j97;

int rmarin(int ij, int kl)
{
  int i, j, k, l, ii, jj, m;
  double s, t;

  if (ij<0 || ij>31328 || kl<0 || kl>30081)
    return 1;

  i = (ij/177)%177 + 2;
  j = ij%177 + 2;
  k = (kl/169)%178 + 1;
  l = kl%169;

  for (ii=1; ii<=97; ii++) {
    s = 0.0;
    t = 0.5;
    for (jj=1; jj<=24; jj++) {
      m = (((i*j)%179)*k) % 179;
      i = j;
      j = k;
      k = m;
      l = (53*l + 1) % 169;
      if ((l*m)%64 >= 32) s += t;
      t *= 0.5;
    }
    u[ii] = s;
  }

  c = 362436.0 / 16777216.0;
  cd = 7654321.0 / 16777216.0;
  cm = 16777213.0 / 16777216.0;

  i97 = 97;
  j97 = 33;

  nCallRanmar = 0;

  return 0;
}

double ranmar(void)
{
  double uni;

  uni = u[i97] - u[j97];
  if (uni < 0.0) uni += 1.0;
  u[i97] = uni;
  i97--;
  if (i97==0) i97 = 97;
  j97--;
  if (j97==0) j97 = 97;
  c -= cd;
  if (c<0.0) c += cm;
  uni -= c;
  if (uni<0.0) uni += 1.0;

  nCallRanmar++;

  return uni;
}


#ifdef STAND_ALONE_1  /* main function for debugging purposes */

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  long    seed, noutput, nskip;
  int     ij, kl;
  long    i, j, less;
  double  x;


  if (argc != 4 || sscanf(argv[1], "%ld", &seed) != 1
                || sscanf(argv[2], "%ld", &noutput) != 1
                || sscanf(argv[3], "%ld", &nskip) != 1)
  {
    printf("usage: ranmar seed# #output #skip\n");
    return 1;
  }

  ij = abs(seed);

  kl = 1.234 * ij;
  if (ij < 0) ij = -ij;
  if (kl < 0) kl = -kl;
  ij %= 31328;
  kl %= 30081;

  if (rmarin(ij, kl) != 0)
  {
    printf("Internal Error: call rmarin()\n");
    return 1;
  }

  less = 0;

  for (i = 0; i < noutput; i++)
  {
    for (j = 0; j < nskip; j++)
    {
      x = ranmar();
      if (x < .5) less++;
    }

    printf("%15d  %17.15f  %15ld  %6.4f\n",
      nCallRanmar, x, less, (double) less / nCallRanmar);
  }

  return 0;
}

#endif  /* STAND_ALONE_1 */


#ifdef STAND_ALONE_2  /* main function */

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  long    seed, noutput;
  int     ij, kl;
  long    i;
  double  x;


  if (argc != 3 || sscanf(argv[1], "%ld", &seed) != 1
                || sscanf(argv[2], "%ld", &noutput) != 1)
  {
    printf("usage: ranmar seed# #output\n");
    return 1;
  }

  ij = abs(seed);

  kl = 1.234 * ij;
  if (ij < 0) ij = -ij;
  if (kl < 0) kl = -kl;
  ij %= 31328;
  kl %= 30081;

  if (rmarin(ij, kl) != 0)
  {
    printf("Internal Error: call rmarin()\n");
    return 1;
  }

  for (i = 0; i < noutput; i++)
  {
    x = ranmar();
    Fprintf(stdout, "%8.6f\n", x);
  }

  return 0;
}

#endif  /* STAND_ALONE_2 */
