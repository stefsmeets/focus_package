#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


static const char *progn = "coseq_reduce";


#define progerror(message) p_rogerror(message, __LINE__)

static void p_rogerror(const char *message, const int source_code_line)
{
  fflush(stdout);
  fprintf(stderr, "%s(%d): %s\n", progn, source_code_line, message);
  exit(1);
}


#define NotEnoughCore() N_otEnoughCore(__LINE__)

static void N_otEnoughCore(const int source_code_line)
{
  p_rogerror("Not enough core", source_code_line);
}


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
  else
    return 1;
}


typedef struct
  {
    int  n;
    int  *Nk;
  }
  T_NodeSeq;


static int NodeSeqSF(const T_NodeSeq *a, const T_NodeSeq *b)
{
  int  i;


  if (a->n < b->n) return -1;
  if (a->n > b->n) return  1;

  for (i = 0; i < a->n; i++)
  {
    if (a->Nk[i] < b->Nk[i]) return -1;
    if (a->Nk[i] > b->Nk[i]) return  1;
  }

  return 0;
}


static void PrintNodeSeq(const char *Lbl, T_NodeSeq *NS, int nNS)
{
  int  iNS, i;


  qsort(NS, nNS, sizeof (*NS),
        (int (*)(const void *, const void *)) NodeSeqSF);

  for (iNS = 0; iNS < nNS; iNS++, NS++)
  {
    fprintf(stdout, "%s", Lbl);
    for (i = 0; i < NS->n; i++) fprintf(stdout, " %d", NS->Nk[i]);
    putc('\n', stdout);
  }
}


static void usage(void)
{
  fprintf(stderr, "usage: %s coseq_file\n", progn);
  exit(1);
}


int main(int argc, char *argv[])
{
  int         i, n;
  char        buf[256];
  int         lcount;
  const char  *fnin;
  FILE        *fpin;

  int  ValidLbl, more;
  char Lbl[256], *fld, *cp;

  int        mNodeSeq, nNS, iNS, mSeq;
  T_NodeSeq  *NodeSeq, *NS;

  fnin = NULL;

  for (i = 1; i < argc; i++)
  {
    if (argv[i][0] == '-')
      usage();
    else if (fnin == NULL)
      fnin = argv[i];
    else
      usage();
  }

  if (fnin)
  {
        fpin = fopen(fnin, "r");
    if (fpin == NULL)
    {
      fprintf(stderr, "%s: Can't open %s\n", progn, fnin);
      exit(1);
    }
  }
  else
    fpin = stdin;

  mNodeSeq = 32;
      mSeq = 64;

      NodeSeq = malloc(mNodeSeq * sizeof (*NodeSeq));
  if (NodeSeq == NULL)
    NotEnoughCore();

  NS = NodeSeq;

  for (iNS = 0; iNS < mNodeSeq; iNS++, NS++)
  {
        NS->Nk = malloc(mSeq * sizeof (*NS->Nk));
    if (NS->Nk == NULL)
      NotEnoughCore();
  }

  ValidLbl = 0;
   NS = NodeSeq;
  nNS = 0;

  lcount = 0;

  while (fgetline(fpin, buf, sizeof buf / sizeof (*buf)))
  {
    lcount++;

    if (strstr(buf, "equal") != NULL) continue;

    fld = buf;

    while (*fld && isspace(*fld)) fld++;
    if (*fld == '\0' || *fld == '#') continue;

    cp = fld; while (*cp && ! isspace(*cp)) cp++;
    if (*cp) { *cp++ = '\0'; more = 1; } else more = 0;

    if (ValidLbl && strcmp(Lbl, fld) != 0)
    {
      PrintNodeSeq(Lbl, NodeSeq, nNS);
      ValidLbl = 0;
       NS = NodeSeq;
      nNS = 0;
    }
    else if (nNS == mNodeSeq)
      progerror("mNodeSeq too small");

    if (ValidLbl == 0) {
      strcpy(Lbl, fld);
      ValidLbl = 1;
    }

    NS->n = 0;

    while (more)
    {
      while (*cp && isspace(*cp)) cp++;
      if (*cp == '\0') break;
      n = sscanf(cp, "%d", &i);
      if (n != 1) IllegalLine(fnin, lcount);
      while (*cp && ! isspace(*cp)) cp++;
      if (NS->n == mSeq) progerror("mSeq too small");
      NS->Nk[NS->n] = i;
             NS->n++;
    }

    for (iNS = 0; iNS < nNS; iNS++)
      if (NodeSeqSF(NS, &NodeSeq[iNS]) == 0) break;

    if (iNS == nNS)
    {
      nNS++;
       NS++;
    }
  }

  if (ValidLbl && strcmp(Lbl, fld) != 0)
    PrintNodeSeq(Lbl, NodeSeq, nNS);

  exit(0);
  return 0;
}
