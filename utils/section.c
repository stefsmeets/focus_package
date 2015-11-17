#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


#define Fprintf (void) fprintf
#define Fclose  (void) fclose
#define Fflush  (void) fflush


static const char *progn = "section";


#define progerror(message) p_rogerror(message, __LINE__)

static void p_rogerror(char *message, int source_code_line)
{
  Fflush(stdout);
  Fprintf(stderr, "%s(%d): %s\n", progn, source_code_line, message);
  exit(1);
}


#define InternalError() I_nternalError(__LINE__)

static void I_nternalError(const int source_code_line)
{
  p_rogerror("Internal Error", source_code_line);
}


#define NotEnoughCore() N_otEnoughCore(__LINE__)

static void N_otEnoughCore(int source_code_line)
{
  p_rogerror("Not enough core", source_code_line);
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


static char *str_dup(char *s)
{
  char *c = (char *) malloc((strlen(s) + 1) * sizeof (*s));
  if (c != NULL) (void) strcpy(c, s);
  return c;
}


static int str_icmp(const char *s, const char *t)
{
  char  cs, ct;


  while (*s || *t)
  {
    cs = toupper(*s++);
    ct = toupper(*t++);
    if (cs < ct) return -1;
    if (cs > ct) return  1;
  }

  return 0;
}


static int str_ibegin(char *s1, char *s2) /* string ignore-case begin */
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


static int nospaces(char *s)
{
  int   i;
  char  *cp, *cpf, *cpl;

  cp = s; while (*cp && isspace(*cp)) cp++;
  cpf = cpl = cp;
  while (*cp) if (isspace(*cp++) == 0) cpl = cp;

  *cpl = '\0';
  i = cpl - cpf;

  if (*s && isspace(*s))
  {
    while (*cpf) *s++ = *cpf++;
    *s = '\0';
  }

  return i;
}


static int is_legal_c_format(char *cfmt, char *types)
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
        if (strchr(types, *cfmt) == NULL) return 0;
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


static int keycmp(char *s1, char *s2)
{
  while (*s1)
    if (*s1++ != *s2++) return -1;

  if (*s2 && *s2 != ':') return -1;

  return 0;
}


typedef struct
  {
    int           iPos;
    unsigned int  iBlock;
    char          *currkey;
    long          Offset;
  }
  T_BlockPtr;


static int BlockBlockPtrSF(const T_BlockPtr *a, const T_BlockPtr *b)
{
  if (a->iBlock < b->iBlock) return -1;
  if (a->iBlock > b->iBlock) return  1;
  return 0;
}


static int PosBlockPtrSF(const T_BlockPtr *a, const T_BlockPtr *b)
{
  if (a->iPos < b->iPos) return -1;
  if (a->iPos > b->iPos) return  1;
  return 0;
}


static void usage(void)
{
  Fprintf(stderr, "usage: %s %s\n",
    progn,
    "[-v] [-BeginEnd] [-BlockNo[=\"c-format\"]]"
    " [-Take=file|-] [file|- [key [#key]]]");

  exit(1);
}


#define BUFLEN 4095


int main(int argc, char *argv[])
{
  int   i, n;
  char  *fnin;
  FILE  *fpin;

  int   F_v, F_BeginEnd;
  int   cfmt_u, nkey;
  char  *cfmt, *key;
  char  *fntake;
  FILE  *fptake;

  int           nBlockPtr, iBlockPtr, jBlockPtr;
  T_BlockPtr    *BlockPtr;
  int           EchoInput, InBlock, DelayPrint;
  unsigned int  iBlock;

  static char  buf[BUFLEN + 1], currkey[BUFLEN + 1];


  fnin = NULL;
  fpin = NULL;
  F_v = 0;
  F_BeginEnd = 0;
  cfmt_u = 0;
  cfmt = NULL;
  fntake = NULL;
  fptake = NULL;
  nkey = -1;
  key = NULL;

  for (i = 1; i < argc; i++)
  {
    if      (str_icmp(argv[i], "-v") == 0)
      F_v = 1;
    else if (str_icmp(argv[i], "-BeginEnd") == 0)
      F_BeginEnd = 1;
    else if (str_ibegin(argv[i], "-BlockNo") == 0)
    {
      if (cfmt) usage();

           cfmt = &argv[i][8];
      if (*cfmt == '\0')
      {
        cfmt = "%u";
        cfmt_u = 1;
      }
      else if (*cfmt++ == '=')
      {
        if      (is_legal_c_format(cfmt, "di"))
          cfmt_u = 0;
        else if (is_legal_c_format(cfmt, "ouxX"))
          cfmt_u = 1;
        else
          usage();
      }
      else
        usage();
    }
    else if (str_ibegin(argv[i], "-Take=") == 0)
    {
      if (fntake || fptake) usage();

      fntake = argv[i] + 6;

      if (strcmp(fntake, "-") == 0)
      {
        if (fpin) usage();

        fntake = NULL;
        fptake = stdin;
      }
    }
    else if (strcmp(argv[i], "-") == 0)
    {
      if (fnin || fpin || fptake) usage();
      fpin = stdin;
    }
    else if (argv[i][0] == '-')
      usage();
    else
    {
      if      (fnin == NULL && fpin == NULL)
        fnin = argv[i];
      else if (key == NULL)
        key = argv[i];
      else if (nkey == -1)
      {
            n = sscanf(argv[i++], "%d", &nkey);
        if (n != 1 || nkey < 0) usage();
      }
      else
        usage();
    }
  }

  if (F_v)
  {
    putc('#', stdout);

    for (i = 0; i < argc;)
      Fprintf(stdout, " %s", argv[i++]);

    putc('\n', stdout);
  }

  nBlockPtr = 0;
   BlockPtr = NULL;

  if (fntake || fptake)
  {
    if (fntake)
    {
          fptake = fopen(fntake, "r");
      if (fptake == NULL)
      {
        Fprintf(stderr, "%s: Can't open %s\n", progn, fntake);
        exit(1);
      }
    }

       nBlockPtr = 1000;
        BlockPtr = malloc(nBlockPtr * sizeof (*BlockPtr));
    if (BlockPtr == NULL)
      NotEnoughCore();

    iBlockPtr = 0;

    if (nkey >= 0)
    {
      BlockPtr[iBlockPtr].iPos    = iBlockPtr;
      BlockPtr[iBlockPtr].iBlock  = nkey;
      BlockPtr[iBlockPtr].currkey = NULL;
      BlockPtr[iBlockPtr].Offset  = -1;
               iBlockPtr++;
    }

    while (fgetline(fptake, buf, sizeof buf))
    {
      (void) nospaces(buf);
      if (! isdigit(buf[0])) continue;

          n = sscanf(buf, "%d", &i);
      if (n == 1 && i >= 0)
      {
        for (jBlockPtr = 0; jBlockPtr < iBlockPtr; jBlockPtr++)
          if (BlockPtr[jBlockPtr].iBlock == i)
            break;

        if (jBlockPtr == iBlockPtr)
        {
          if (iBlockPtr == nBlockPtr)
          {
               nBlockPtr += 10000;
                BlockPtr = realloc(BlockPtr, nBlockPtr * sizeof (*BlockPtr));
            if (BlockPtr == NULL)
              NotEnoughCore();
          }

          BlockPtr[iBlockPtr].iPos    = iBlockPtr;
          BlockPtr[iBlockPtr].iBlock  = i;
          BlockPtr[iBlockPtr].currkey = NULL;
          BlockPtr[iBlockPtr].Offset  = -1;
                   iBlockPtr++;
        }
      }
    }

    if (fptake != stdin)
      Fclose(fptake);

    nBlockPtr = iBlockPtr;

    if (nBlockPtr == 0)
    {
      if (fntake == NULL)
          fntake = "take file";

      Fprintf(stderr, "%s: No numbers in %s\n", progn, fntake);
      exit(1);
    }

        BlockPtr = realloc(BlockPtr, nBlockPtr * sizeof (*BlockPtr));
    if (BlockPtr == NULL)
      InternalError();

    qsort(BlockPtr, nBlockPtr, sizeof (*BlockPtr),
          (int (*)(const void *, const void *)) BlockBlockPtrSF);

    if (F_v)
    {
      Fprintf(stdout, "# Take %d blocks from %u ... %u\n",
        nBlockPtr, BlockPtr[0].iBlock, BlockPtr[nBlockPtr - 1].iBlock);
    }
  }

  if (fpin == NULL)
  {
    if (fnin == NULL)
      fpin = stdin;
    else
    {
          fpin = fopen(fnin, "r");
      if (fpin == NULL)
      {
        Fprintf(stderr, "%s: Can't open %s\n", progn, fnin);
        exit(1);
      }
    }
  }

  DelayPrint = 0;

  if (nBlockPtr && fpin != stdin)
  {
    for (iBlockPtr = 1; iBlockPtr < nBlockPtr; iBlockPtr++) {
      if (BlockPtr[iBlockPtr - 1].iPos > BlockPtr[iBlockPtr].iPos) {
        DelayPrint = 1;
        break;
      }
    }
  }

      iBlockPtr = 0;
  if (iBlockPtr < nBlockPtr)
    nkey = BlockPtr[iBlockPtr++].iBlock;

  InBlock   = 0;
  iBlock    = 0;
  EchoInput = 0;

  while (fgetline(fpin, buf, sizeof buf))
  {
    if (InBlock == 0)
    {
      if (str_ibegin(buf, ">Begin ") == 0)
      {
        InBlock = 1;

               (void) nospaces(buf + 7);
        (void) strcpy(currkey, buf + 7);

        if (key == NULL || keycmp(key, currkey) == 0)
        {
          if (nkey == -1 || iBlock == nkey)
          {
            if (DelayPrint == 0)
            {
              if (F_BeginEnd)
                Fprintf(stdout, ">Begin %s\n", currkey);
            }
            else
            {
              if (key && strcmp(key, currkey) == 0)
                BlockPtr[iBlockPtr - 1].currkey = key;
              else
                BlockPtr[iBlockPtr - 1].currkey = str_dup(currkey);

              BlockPtr[iBlockPtr - 1].Offset = ftell(fpin);
            }

            EchoInput = 1;
          }
          else
            EchoInput = -1;
        }
      }
    }
    else
    {
      if (str_ibegin(buf, ">End ") == 0)
      {
        (void) nospaces(buf + 5);

        if (keycmp(buf + 5, currkey) == 0)
        {
          if (EchoInput != 0)
          {
            if (F_BeginEnd && EchoInput == 1 && DelayPrint == 0)
              Fprintf(stdout, ">End %s\n\n", currkey);

            iBlock++;

            if (EchoInput == 1)
            {
              if (iBlockPtr < nBlockPtr)
                nkey = BlockPtr[iBlockPtr++].iBlock;
              else if (nkey != -1)
                break;
            }
          }

          InBlock   = 0;
          EchoInput = 0;
        }
      }

      if (EchoInput == 1 && DelayPrint == 0)
      {
        if (cfmt)
        {
          if (cfmt_u)
            Fprintf(stdout, cfmt,       iBlock);
          else
            Fprintf(stdout, cfmt, (int) iBlock);

          if (buf[0]) putc(' ', stdout);
        }

        Fprintf(stdout, "%s\n", buf);
      }
    }
  }

  if (DelayPrint)
  {
    qsort(BlockPtr, nBlockPtr, sizeof (*BlockPtr),
          (int (*)(const void *, const void *)) PosBlockPtrSF);

    for (iBlockPtr = 0; iBlockPtr < nBlockPtr; iBlockPtr++)
    {
      if (BlockPtr[iBlockPtr].Offset < 0) continue;

      if (fseek(fpin, BlockPtr[iBlockPtr].Offset, SEEK_SET) != 0)
      {
        Fprintf(stderr, "%s: Random access error for file \"%s\"\n",
          progn, fnin);
      }

      if (F_BeginEnd)
        Fprintf(stdout, ">Begin %s\n", BlockPtr[iBlockPtr].currkey);

      while (fgetline(fpin, buf, sizeof buf))
      {
        if (str_ibegin(buf, ">End ") == 0)
        {
          (void) nospaces(buf + 5);

          if (keycmp(buf + 5, BlockPtr[iBlockPtr].currkey) == 0)
          {
            if (F_BeginEnd)
              Fprintf(stdout, ">End %s\n\n", BlockPtr[iBlockPtr].currkey);

            break;
          }
        }

        if (cfmt)
        {
          if (cfmt_u)
            Fprintf(stdout, cfmt,       BlockPtr[iBlockPtr].iBlock);
          else
            Fprintf(stdout, cfmt, (int) BlockPtr[iBlockPtr].iBlock);

          if (buf[0]) putc(' ', stdout);
        }

        Fprintf(stdout, "%s\n", buf);
      }
    }
  }

  if (fpin != stdin)
    Fclose(fpin);

  return 0;
}
