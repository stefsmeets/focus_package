#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


static const char *progn = "coseq_cmp";


#define BUFLEN 4095


#define progerror(message) p_rogerror(message, __LINE__)

static void p_rogerror(const char *message, int source_code_line)
{
  fflush(stdout);
  fprintf(stderr, "%s(%d): %s\n", progn, source_code_line, message);
  exit(1);
}


#define InternalError(message) I_nternalError(message, __LINE__)

static void I_nternalError(const char *message, int source_code_line)
{
  fflush(stdout);
  fprintf(stderr, "%s(%d): Internal Error: %s\n",
    progn, source_code_line, message);
  exit(1);
}


#define NotEnoughCore() N_otEnoughCore(__LINE__)

static void N_otEnoughCore(int source_code_line)
{
  p_rogerror("Not enough core", source_code_line);
}


#define IllegalLine(fnin, lcount) I_llegalLine(fnin, lcount, __LINE__)

static void I_llegalLine(const char *fnin, long lcount, int source_code_line)
{
  fflush(stdout);
  fprintf(stderr, "%s(%d): Illegal line #%ld",
    progn, source_code_line, lcount);
  if (fnin != NULL) fprintf(stderr, " in %s", fnin);
  putc('\n', stderr);
  exit(1);
}


#ifdef DBG_MALLOC

static void *dbgmalloc(size_t n)
{
  fprintf(stdout, "malloc(%ld)\n", (long) n);
  return malloc(n);
}

#define malloc(n) dbgmalloc(n)

#endif


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


static char *str_dup(const char *s)
{
  char *c = (char *) malloc((strlen(s) + 1) * sizeof (*s));
  if (c != NULL) strcpy(c, s);
  return c;
}


static int str_icmp(const char *s, const char *t)
{
  char     cs, ct;

  while (*s || *t)
  { cs = toupper(*s++);
    ct = toupper(*t++);
    if (cs < ct) return -1;
    if (cs > ct) return  1;
  }
  return 0;
}


static int str_ibegin(const char *s1, const char *s2)
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


static void WrnIfUnsorted(const char *fnin, long lcount,
                          const int *Seq, int aSeq, int bSeq, int nSeq)
{
  int  i, j;


  for (i = aSeq, j = bSeq; i < bSeq && j < nSeq; i++, j++)
  {
    if (Seq[i] > Seq[j])
    {
      fflush(stdout);
      fprintf(stderr,
        "%s: Warning: Sequences not sorted before or at line #%ld",
        progn, lcount);
      if (fnin != NULL) fprintf(stderr, " in %s", fnin);
      putc('\n', stderr);
      break;
    }

    if (Seq[i] < Seq[j])
      break;
  }
}


typedef struct S_ch_str
  {
    struct S_ch_str  *Next;
    char             *str;
  }
  T_ch_str;


typedef struct S_CS_Block
  {
    struct S_CS_Block  *Next;
    T_ch_str           *ch_Lbl;
    int                nSeq;
    int                *Seq;
  }
  T_CS_Block;


static void AddCS_Block(T_CS_Block **pCS_Blocks, char *Lbl, int *Seq, int nSeq)
{
  T_ch_str    **p_ch_Lbl;
  T_CS_Block  *CSB;


  while ((CSB = *pCS_Blocks) != NULL)
  {
    if (    nSeq == CSB->nSeq
        && (nSeq == 0 || memcmp(CSB->Seq, Seq, nSeq * sizeof (*Seq)) == 0))
      break;

    pCS_Blocks = &CSB->Next;
  }

  if (! CSB)
  {
          CSB = *pCS_Blocks = malloc(sizeof (*CSB));
    if (! CSB)
      NotEnoughCore();

    CSB->Next = NULL;

    p_ch_Lbl = &CSB->ch_Lbl;

    CSB->nSeq = nSeq;

    if (nSeq == 0)
      CSB->Seq = NULL;
    else
    {
            CSB->Seq = malloc(nSeq * sizeof (*CSB->Seq));
      if (! CSB->Seq)
        NotEnoughCore();

      memcpy(CSB->Seq, Seq, nSeq * sizeof (*Seq));
    }
  }
  else
  {
    p_ch_Lbl = &CSB->ch_Lbl->Next;
    while (*p_ch_Lbl) p_ch_Lbl = &(*p_ch_Lbl)->Next;
  }

        *p_ch_Lbl = malloc(sizeof (*(*p_ch_Lbl)));
  if (! *p_ch_Lbl)
    NotEnoughCore();

  (*p_ch_Lbl)->Next = NULL;

        (*p_ch_Lbl)->str = str_dup(Lbl);
  if (! (*p_ch_Lbl)->str)
    NotEnoughCore();
}


static T_CS_Block *ReadCSfile(char *fnin, FILE *fpin, int *nInputBlocks)
{
  int   i, n, more;
  long  lcount;
  char  *fld, *cp;

  char  xtrac;

  int   *Seq;
  int   mSeq, nSeq, aSeq, bSeq;

  static char  buf[BUFLEN + 1];
  static char  Lbl[BUFLEN + 1];

  T_CS_Block  *CS_Blocks;
  int         nIB_Buf;


  if (nInputBlocks == NULL)
      nInputBlocks = &nIB_Buf;

                   mSeq = 1000;
      Seq = malloc(mSeq * sizeof (*Seq));
  if (Seq == NULL)
    NotEnoughCore();

  CS_Blocks = NULL;
  *nInputBlocks = 0;

  Lbl[0] = '\0';
  nSeq = 0;
  aSeq = 0;
  bSeq = 0;

  lcount = 0;

  while (fgetline(fpin, buf, sizeof buf / sizeof (*buf)))
  {
    lcount++;

    fld = buf;

    for (;;)
    {
      while (*fld && isspace(*fld)) fld++;

      if (*fld == '\0' || *fld == '#')
        break;

      cp = fld; while (*cp && ! isspace(*cp)) cp++;
      if (*cp) { *cp = '\0'; more = 1; } else more = 0;

      if (! isdigit(*fld))
      {
        if (strcmp(Lbl, fld) != 0)
        {
          if (Lbl[0]) {
            AddCS_Block(&CS_Blocks, Lbl, Seq, nSeq);
            (*nInputBlocks)++;
          }

          strcpy(Lbl, fld);
          nSeq = 0;
          aSeq = 0;
          bSeq = 0;
        }
        else
        {
          WrnIfUnsorted(fnin, lcount, Seq, aSeq, bSeq, nSeq);
          aSeq = bSeq;
          bSeq = nSeq;
        }
      }
      else
      {
        if (Lbl[0] == '\0')
          IllegalLine(fnin, lcount);

            n = sscanf(fld, "%d%c", &i, &xtrac);
        if (n != 1)
          IllegalLine(fnin, lcount);

        if (nSeq == mSeq)
        {
                                 mSeq += 1000;
              Seq = realloc(Seq, mSeq * sizeof (*Seq));
          if (Seq == NULL)
            NotEnoughCore();
        }

        Seq[nSeq++] = i;
      }

      fld = cp; if (more) fld++;
    }
  }

  if (Lbl[0])
  {
    AddCS_Block(&CS_Blocks, Lbl, Seq, nSeq);
    (*nInputBlocks)++;

    WrnIfUnsorted(fnin, lcount, Seq, aSeq, bSeq, nSeq);
  }

  free(Seq);

  return CS_Blocks;
}


typedef  struct
  {
    T_CS_Block  *CSB;
    int         n;
  }
  T_Index_CSB;


static int Index_CSB_SortFunction(T_Index_CSB *a, T_Index_CSB *b)
{
  int       o;
  T_ch_str  *al, *bl;


  if (a->n > b->n) return -1;
  if (a->n < b->n) return  1;

  for (al = a->CSB->ch_Lbl, bl = b->CSB->ch_Lbl;
       al && bl;
       al = al->Next, bl = bl->Next)
  {
        o = strcmp(al->str, bl->str);
    if (o) return o;
  }

  if (al || bl)
    InternalError("al || bl");

  return 0;
}


static T_CS_Block *FindCSB(T_CS_Block *CSB, int *Seq, int nSeq)
{
  for (; CSB != NULL; CSB = CSB->Next)
    if (    nSeq == CSB->nSeq
        && (nSeq == 0 || memcmp(CSB->Seq, Seq, nSeq * sizeof (*Seq)) == 0))
      break;

  return CSB;
}


static void PrintCSB_Labels(T_CS_Block *CSB, int F_all)
{
  int       i;
  T_ch_str  *ch_Lbl;


  i = 0;

  for (ch_Lbl = CSB->ch_Lbl; ch_Lbl; ch_Lbl = ch_Lbl->Next)
  {
    if (i == 0 || F_all)
    {
      if (i) putc(' ', stdout);

      fprintf(stdout, "%s", ch_Lbl->str);
    }
    else
      break;

    i++;
  }
}


static void usage(void)
{
  fprintf(stderr,
    "usage: %s [-v] [-r] [-uniq [-nPerLine=#]] [file [CS-DataBase]]\n",
    progn);

  exit(1);
}


int main(int argc, char *argv[])
{
  int   i, n;
  char  xtrac;
  char  *fnin;
  FILE  *fpin;

  int   F_v, F_r, F_uniq, nPerLine;

  T_ch_str    *ch_Lbl;
  T_CS_Block  *CS_Blocks, *CSB;
  int         nInputBlocks;

  T_Index_CSB  *Index_CSB;
  int          nIndex_CSB;

  char        *fndb;
  FILE        *fpdb;
  T_CS_Block  *DB_CS_Blocks, *DBCSB;


  F_v    = 0;
  F_r    = 0;
  F_uniq = 0;
  nPerLine = 10;
  fnin = NULL;
  fndb = NULL;

  for (i = 1; i < argc; i++)
  {
    if (argv[i][0] == '-')
    {
      if      (str_icmp(argv[i], "-v") == 0)
      {
        if (F_v) usage();
        F_v = 1;
      }
      else if (str_icmp(argv[i], "-r") == 0)
      {
        if (F_r) usage();
        F_r = 1;
      }
      else if (str_icmp(argv[i], "-uniq") == 0)
      {
        if (F_uniq) usage();
        F_uniq = 1;
      }
      else if (str_ibegin(argv[i], "-nPerLine=") == 0)
      {
            n = sscanf(&argv[i][10], "%d %c", &nPerLine, &xtrac);
        if (n != 1 || nPerLine < 1)
          usage();
      }
      else if (! fnin && argv[i][1] == '\0')
        fnin = "";
      else
        usage();
    }
    else
    {
      if      (! fnin) fnin = argv[i];
      else if (! fndb) fndb = argv[i];
      else
        usage();
    }
  }

  if (fnin && fnin[0] == '\0')
    fnin = NULL;

  DB_CS_Blocks = NULL;

  if (fndb == NULL)
      fndb = getenv("COSEQ_DB");
  
  if (fndb == NULL)
      fndb = "/User/smeets/focus/kriber_f/coseq";

  if (fndb)
  {
        fpdb = fopen(fndb, "r");
    if (fpdb == NULL)
      fprintf(stderr, "%s: Warning: Can't open %s\n", progn, fndb);
    else
      DB_CS_Blocks = ReadCSfile(fndb, fpdb, NULL);
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

  CS_Blocks = ReadCSfile(fnin, fpin, &nInputBlocks);

  if (F_uniq == 0)
  {
    nIndex_CSB = 0;

    for (CSB = CS_Blocks; CSB; CSB = CSB->Next)
      nIndex_CSB++;

    if (nIndex_CSB > 0)
    {
            Index_CSB = malloc(nIndex_CSB * sizeof (*Index_CSB));
      if (! Index_CSB)
        NotEnoughCore();

      i = 0;

      for (CSB = CS_Blocks; CSB; CSB = CSB->Next)
      {
        n = 0;

        for (ch_Lbl = CSB->ch_Lbl; ch_Lbl; ch_Lbl = ch_Lbl->Next)
          n++;

        Index_CSB[i].CSB = CSB;
        Index_CSB[i].n   = n;
                  i++;
      }

      if (i != nIndex_CSB)
        InternalError("i != nIndex_CSB");

      qsort(Index_CSB, nIndex_CSB, sizeof (*Index_CSB),
            (int (*)(const void *, const void *)) Index_CSB_SortFunction);

      for (i = 0; i < nIndex_CSB; i++)
      {
        CSB = Index_CSB[i].CSB;

        n = 0;

        for (ch_Lbl = CSB->ch_Lbl; ch_Lbl; ch_Lbl = ch_Lbl->Next)
          n++;

        if (n != Index_CSB[i].n)
          InternalError("n != Index_CSB[i].n");

        if (F_r)
          fprintf(stdout, "%.6g %d : ", (double) n / nInputBlocks, n);

        PrintCSB_Labels(CSB, F_v);

        if (F_r == 0)
          fprintf(stdout, " : %d %.6g", n, (double) n / nInputBlocks);

        if (DB_CS_Blocks)
        {
              DBCSB = FindCSB(DB_CS_Blocks, CSB->Seq, CSB->nSeq);
          if (DBCSB)
          {
            fprintf(stdout, " = ");
            PrintCSB_Labels(DBCSB, 1);
          }
        }

        putc('\n', stdout);
      }
    }
  }
  else
  {
    for (CSB = CS_Blocks; CSB; CSB = CSB->Next)
    {
      fprintf(stdout, "%s", CSB->ch_Lbl->str);

      for (i = 0; i < CSB->nSeq; i++)
      {
        if (i != 0 && i % nPerLine == 0)
          fprintf(stdout, "\n%s", CSB->ch_Lbl->str);

        fprintf(stdout, " %d", CSB->Seq[i]);
      }

      if (i)
        putc('\n', stdout);
    }
  }

  return 0;
}
