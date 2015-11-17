#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


static char  *progn = "genseq";


static int is_legal_c_format(char *cfmt)
{
  int   count, fmt_on, fmt_width;
  char  *start_digits;


  count = 0;
  start_digits = NULL;

  for (fmt_on = fmt_width = 0; *cfmt; cfmt++)
  {
    if (isdigit(*cfmt))
    {
      if (start_digits == NULL)
          start_digits = cfmt;
    }

    if (fmt_on == 0)
    {
      if (*cfmt == '%')
      {
        fmt_on = 1;
        fmt_width = 0;
      }
    }
    else
    {
      if      (*cfmt == '%')
      {
        if (fmt_width) return 0;
        fmt_on = 0;
      }
      else if (*cfmt == '$')
      {
        if (start_digits == NULL) return 0;
        if (cfmt - start_digits != 1) return 0;
        if (*start_digits != '1') return 0;
      }
      else if (*cfmt == '*')
        return 0;
      else if (isalpha(*cfmt))
      {
        if (strchr("feEgG", *cfmt) == NULL) return 0;
        fmt_on = 0;
        count++;
      }

      fmt_width++;
    }

    if (isdigit(*cfmt) == 0)
      start_digits = NULL;
  }

  if (fmt_on || count != 1) return 0;

  return 1;
}


static void usage(void)
{
  fprintf(stderr, "usage: %s First Last Step [C-Format]\n", progn);
  exit(1);
}


int main(int argc, char *argv[])
{
  int     n;
  long    i;
  double  First, Last, Step, Val;
  char    *c_format, xtrac;


  if (argc != 4 && argc != 5)
    usage();

      n = sscanf(argv[1], "%lf%c", &First, &xtrac);
  if (n != 1) usage();
      n = sscanf(argv[2], "%lf%c", &Last,  &xtrac);
  if (n != 1) usage();
      n = sscanf(argv[3], "%lf%c", &Step,  &xtrac);
  if (n != 1 || Step == 0.) usage();

  if (argc > 4)
  {                       c_format = argv[4];
    if (is_legal_c_format(c_format) == 0)
    {
      fprintf(stderr, "%s: Illegal \"C-Format\"\n", progn);
      exit(1);
    }
  }
  else
    c_format = "%10g";

  for (i = 0;; i++)
  {
    Val = First + i * Step;

    if (Step < 0.)
    {
      if (Val < Last) break;
    }
    else
    {
      if (Val > Last) break;
    }

        n = fprintf(stdout, c_format, Val); putc('\n', stdout);
    if (n < 1)
      exit(1);
  }
  
  return 0;
}
