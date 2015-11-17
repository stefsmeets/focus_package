#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>


#include "lib.h"


int str_icmp(const char *s, const char *t)
{
  char     cs, ct;


  while (*s || *t)
  {
    cs = toupper(*s++);
    ct = toupper(*t++);
    if (cs < ct) return -1;
    if (cs > ct) return  1;
  }

  return 0;
}


char *lower(char *s)
{
  char   *p;

  p = s;
  while (*s)
  { switch (*s)
    {
#ifdef __TURBOC__
      case 'ô' : *s = 'î'; break;
      case 'é' : *s = 'Ñ'; break;
      case 'ö' : *s = 'Å'; break;
#endif
      default  : *s = tolower(*s);
    }
    s++;
  }
  return p;
}


char *upper(char *s)
{
  char   *p;

  p = s;
  while (*s)
  { switch (*s)
    {
#ifdef __TURBOC__
      case 'î' : *s = 'ô'; break;
      case 'Ñ' : *s = 'é'; break;
      case 'Å' : *s = 'ö'; break;
#endif
      default  : *s = toupper(*s);
    }
    s++;
  }
  return p;
}


char firstnonblank(const char *s)
{
  while (*s && isspace(*s)) s++;
  return *s;
}


int noblks(char *s)
{
  int         i;
  char        *cp, *cpf, *cpl;

  cp = s;
  while (*cp == ' ') cp++;
  cpf = cpl = cp;
  while (*cp) if (*cp++ != ' ') cpl = cp;

  *cpl = '\0';
  i = (int)(cpl - cpf);

  if (*s == ' ')
  { while (*cpf) *s++ = *cpf++;
    *s = '\0';
  }

  return i;
}


int lnblnk(const char *s)
{
  int    i = 0, lnb = 0;
  char   c;

  while ((c = s[i++]) != '\0')
    if (c != ' ') lnb = i;

  return lnb;
}


void noleadblks(char *s)
{
  char   *p;

  p = s;
  while (*p == ' ') p++;
  if (p != s)
  { while (*p) *s++ = *p++;
    *s = '\0';
  }
}


void notrailblks(char *s)
{
  char   *p;

  for (p = NULL; *s; s++) if (*s != ' ') p = s;
  if (p) *(++p) = '\0';
}


int f_get_line(FILE *fpin, char s[], int size_s)
{
  int         last_s, c, i;


  last_s = size_s - 1;

  i = 0;

  while ((c = getc(fpin)) != EOF && c != '\n')
  {
#ifdef __MSDOS__
    if (c == 0x1A /* CtrlZ */) {
      ungetc(c, fpin);
      c = EOF;
      break;
    }
#endif

    if (i < last_s) s[i++] = (char) c;
  }

  s[i] = '\0';

  if (i == 0 && c == EOF)
    return 0;

  return 1;
}


int fgetIline(FILE *fpin)
{
  int         c, i, notab;


  notab = Imode & I_notabs;
  c = 0;

  do
  {
    i = 0;

    while ((c = getc(fpin)) != EOF && c != '\n')
    {
#ifdef __MSDOS__
      if (c == 0x1A /* CtrlZ */) {
        ungetc(c, fpin);
        c = EOF;
        break;
      }
#endif

      if (i < LNIline) {
        if (c == '\t' && notab) {
          do
            Iline[i++] = ' ';
          while (i < LNIline && i % Itabsize);
        }
        else
          Iline[i++] = (char) c;
      }
    }

    Iline[i] = '\0';
    Ilinec++;

    if (Imode & I_noleadblks) noleadblks(Iline);
    if (Imode & I_notrailblks) notrailblks(Iline);
  }
  while (Imode & I_noemptyline && *Iline == '\0' && c != EOF);

  if (i == 0 && c == EOF)
    return 0;

  if      (Imode & I_lowercase) (void) lower(Iline);
  else if (Imode & I_uppercase) (void) upper(Iline);

  return 1;
}


int str_begin(const char *s1, const char *s2) /* string begin */
{
  while (*s1 && *s2)
  { if      (*s1   < *s2  ) return -1;
    else if (*s1++ > *s2++) return  1;
  }
  if (*s2) return -1;
  return 0;
}


int str_ibegin(const char *s1, const char *s2) /* string ignore-case begin */
{
  char     u1, u2;

  while (*s1 && *s2)
  {
    u1 = toupper(*s1++);
    u2 = toupper(*s2++);
    if      (u1 < u2) return -1;
    else if (u1 > u2) return  1;
  }
  if (*s2) return -1;
  return 0;
}


char *str_xcpy(char *dest, const char *src, size_t max)  /* string max copy */
{                                                 /* working like strncpy() */
  int         i = 0;                              /* but dest always        */
                                                  /* terminates with '\0'   */
  while (i < max && *src) dest[i++] = *src++;
  dest[i] = '\0';

  return dest;
}


void fillblks(char *s, int len)
{
  int         i, fill = 0;

  for (i = 0; i < len; i++)
  { if (fill)
    { s[i] = ' ';
    }
    else if (s[i] == '\0')
    { fill = 1;
      s[i] = ' ';
    }
  }
  s[i] = '\0';
}


int str_arg(char *arg, const char *s, int n)  /* - string argument -       */
{                                        /* return n-th argument of string */
  char  *a = arg;                        /* s in arg.                      */
                                         /* n = 0 for the first argument.  */

    while (*s && isspace(*s) != 0) s++;

  while (n--)
  {
    while (*s && isspace(*s) == 0) s++;
    while (*s && isspace(*s) != 0) s++;
  }

    while (*s && isspace(*s) == 0) *a++ = *s++;

  *a = '\0';

  return (a != arg);
}


long flinec(FILE *fp)
{
  int    c, incompletelastline;
  long   count;

  (void) fseek(fp, 0L, SEEK_SET);
  incompletelastline = 0;
  count = 0;

#ifdef __MSDOS__
  while ((c = getc(fp)) != EOF && c != CtrlZ)
#else
  while ((c = getc(fp)) != EOF)
#endif
  { if (c == '\n') count++, incompletelastline = 0;
    else                    incompletelastline = 1;
  }
  (void) fseek(fp, 0L, SEEK_SET);

  return count + incompletelastline;
}


void Sum2bmc(int Nxy,
             double minx, double maxx,
             double miny, double maxy,
             double Sx, double Sx2,
             double Sy, double Sy2, double Sxy,
             double *b, double *m, double *c)
{
  double       dx, dy, d;
  const double ddyn = 1.e-6;


  *b = *m = *c = 0.;

  if (Nxy < 1) return;

  if (minx == maxx) return;
  if (miny == maxy) {
    *b = miny;
    return;
  }

#define MAX(a, b) ((a) >= (b) ? (a) : (b))
  dx =  MAX(fabs(minx - Sx / Nxy),
            fabs(maxx - Sx / Nxy));
  dy =  MAX(fabs(miny - Sy / Nxy),
            fabs(maxy - Sy / Nxy));
#undef  MAX

  if (dx == 0.) return;
  if (dy == 0.) {
    *b = Sy / Nxy;
    return;
  }

  if (dx < dy * ddyn) return;
  if (dy < dx * ddyn) {
    *b = Sy / Nxy;
    return;
  }

      d = Nxy * Sx2 - Sx * Sx;
  if (d != 0.) {
    *b = (Sx2 * Sy - Sx * Sxy) / d;
    *m = (Nxy * Sxy - Sx * Sy) / d;
  }

  d =   (Sx2 - Sx * Sx / Nxy)
      * (Sy2 - Sy * Sy / Nxy);
  if (d > 0.)
    *c = (Sxy - Sx * Sy / Nxy) / sqrt(d);
}


#include "main.h"


int getnum(char *str, int **num, int max_nnum)
{
  int     i, n, nnum;
  char    buf[128], xtrac, *s;
  double  fn;


  nnum = 0;
  i = 0;

  for (s = str; *s; s++)
  {
    if (*s == ',')
    {
      nnum++;
      i = -1;
    }
    else if (isspace(*s))
    {
      if (i == 1)
      {
        nnum++;
        i = 0;
      }
    }
    else
      i = 1;
  }

  if (i)
    nnum++;

  if (max_nnum && nnum > max_nnum)
    return -1;

  if (nnum)
    CheckMalloc(*num, nnum);

  s = str;

  for (i = 0; i < nnum; i++)
  {
    while (*s && isspace(*s)) s++;

    n = 0;

    while (n < sizeof buf / sizeof (*buf))
    {
      if (*s == '\0' || *s == ',' || isspace(*s))
        break;

      buf[n++] = *s++;
    }

    if (n == sizeof buf / sizeof (*buf))
    {
      AppFree(*num, nnum);
      *num = NULL;
      return -(i + 1);
    }

    buf[n] = '\0';

    if (n == 0)
      fn = 0.;
    else
    {
          n = sscanf(buf, "%lf%c", &fn, &xtrac);
      if (n != 1)
      {
        AppFree(*num, nnum);
        *num = NULL;
        return -(i + 1);
      }

      if (fn < 0.)
          fn -= .5;
      else
          fn += .5;
    }

    (*num)[i] = (int) fn;

    if (*s == ',')
      s++;
  }

  return nnum;
}
