#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>


static const char  *progn = "xrs2scandata";


#define IllegalLine(fnin, lcount) I_llegalLine(fnin, lcount, __LINE__)

static void I_llegalLine(const char *fnin, int lcount, int source_code_line)
{
  fflush(stdout);
  fprintf(stderr, "%s(%d): Illegal line #%d",
    progn, source_code_line, lcount);
  if (fnin != NULL) fprintf(stderr, " in %s", fnin);
  putc('\n', stderr);
  exit(1);
}


static int fgetline(FILE *fpin, char s[], int size_s)
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


static char firstnonblank(char *s)
{
  while (*s == ' ' || *s == '\t') ++s;
  return *s;
}


static int str_ibegin(char *s1, char *s2) /* string ignore case begin */
{
  int      u1, u2;

  while (*s1 && *s2)
  { u1 = toupper(*s1++); u2 = toupper(*s2++);
    if      (u1 < u2) return -1;
    else if (u1 > u2) return  1;
  }
  if (*s2) return -1;
  return 0;
}


static void usage(void)
{
  fprintf(stderr, "usage: %s [-Scale=<factor>] [xrs-File]\n", progn);
  exit(0);
}


int main(int argc, char *argv[])
{
  int     iarg, ifld, idata, n;
  int     lcount;
  FILE    *fpin;
  char    *fnin, buf[256], xtrac;
  double  Scale;
  double  Ttheta, Intens;


  Scale = 1.;

  fnin = NULL;
  fpin = NULL;

  for (iarg = 1; iarg < argc; iarg++)
  {
    if (argv[iarg][0] == '-') {
      if (str_ibegin(argv[iarg], "-Scale=") == 0) {
            n = sscanf(argv[iarg] + 7, "%lf %c", &Scale, &xtrac);
        if (n != 1) usage();
      }
      else
        usage();
    }
    else {
      if (fnin) usage();
      fnin = argv[iarg];
    }
  }

  if (fnin == NULL)
    fpin = stdin;
  else
  {
        fpin = fopen(fnin, "r");
    if (fpin == NULL)
    {
      fprintf(stderr, "%s: Can't open %s\n", progn, fnin);
      exit(1);
    }
  }

  ifld   = 0;
  idata  = 0;
  lcount = 0;

  while (fgetline(fpin, buf, sizeof buf))
  {
    lcount++;

    if (firstnonblank(buf) != '#')
    {
          n = sscanf(buf, "%lf%lf", &Ttheta, &Intens);
      if (n != 2)
        IllegalLine(fnin, lcount);

      Intens *= Scale;

      if (Intens < 1.) Intens = 1.;

      if (ifld == 8)
      {
        putc('\n', stdout);
	ifld = 0;
      }

      if (ifld == 0)
        putc(' ', stdout);

      fprintf(stdout, "%8.0f", Intens);
      ifld++;

          idata++;
      if (idata == 4095) break;
    }
  }

  if (fpin != stdin)
    fclose(fpin);

  if (ifld)
    putc('\n', stdout);

  fprintf(stdout, " %8.0f\n", 0.);

  exit(0);
  return 0;
}
