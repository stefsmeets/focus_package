/*
  Space Group Info (c) 1994-97 Ralf W. Grosse-Kunstleve
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>


/*
  Macintosh extras (Courtesy Jon Tischler <TischlerJZ@ornl.gov>)
 */
#if defined(__THINK__) || defined(__MWERKS__)
#include <console.h>
#define CONSOLE_LINES   36  /* number of lines to use for console */
#define CONSOLE_COLUMNS 90  /* number of columns to use for console */
#ifdef __MWERKS__
#include <sioux.h>
#endif
#endif


#define SGCOREDEF__
#include "sginfo.h"


/* sgssvfy.c */

int Check_ssVM(T_SgInfo *SgInfo);
int xtal32ss(T_SgInfo *SgInfo);


/* sghklvfy.c */

int Verify_sghkl(const T_SgInfo *SgInfo, int FriedelSym,
                 int Maxh, int Maxk, int Maxl);


/* cnslib.c */

int CNSrecord(FILE *fpout, T_SgInfo *SgInfo, int F_All);
int MakeCNSlib(FILE *fpout, int F_All);


/* sgnorm.c */

void TestAddlG(void);


#define Fputs   (void) fputs
#define Fprintf (void) fprintf
#define Fflush  (void) fflush


static const char *progn = "sginfo";


static void progerror(const char *message)
{
  Fflush(stdout);
  Fprintf(stderr, "%s: %s\n", progn, message);
  exit(1);
}


static void NotEnoughCore(void)
{
  progerror("Not enough core");
}


static void PrintClearSgError(int ClearError, int CertainSgError)
{
  if (CertainSgError && SgError == NULL)
    SetSgError("Internal Error: SgError not set but should be");

  if (SgError)
  {
    Fflush(stdout);
    Fprintf(stderr, "%s: %s\n", progn, SgError);
    if (ClearError == 0) exit(1);
    SgError = NULL;
  }
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


static int str_ibegin(const char *s1, const char *s2) /* string ignore-case */
{                                                     /* begin              */
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


static const char *LegendTabSgName[] =
  {
    "",
    "  Extensions",
    "  ----------",
    "    Monoclinic             unique axis b   unique axis c   unique axis a",
    "                             abc   c-ba      abc   ba-c      abc   -acb",
    "                            ------------    ------------    ------------",
    "             cell choice 1   :b1   :-b1      :c1   :-c1      :a1   :-a1",
    "                         2   :b2   :-b2      :c2   :-c2      :a2   :-a2",
    "                         3   :b3   :-b3      :c3   :-c3      :a3   :-a3",
    "",
    "    Orthorhombic   :ba-c    change of basis abc -> ba-c",
    "                   :1       origin choice 1",
    "                   :2ba-c   origin choice 2, change of basis abc -> ba-c",
    "",
    "    Tetragonal     :1       origin choice 1",
    "           Cubic   :2       origin choice 2",
    "",
    "    Trigonal       :H       hexagonal    axes",
    "                   :R       rhombohedral axes",
    "",
    "  Number   Schoenflies   Hermann-Mauguin             Hall",
    "  ------   -----------   ---------------             ----",
    NULL,
  };


static void ListTabSgName(int WantedSgNumber, int VolLetter, FILE *fpout)
{
  int                i;
  const char         *sgl, *ext, **ltsgn;
  const T_TabSgName  *tsgn, *show, *show_later;


  if (WantedSgNumber == -1)
    for (ltsgn = LegendTabSgName; *ltsgn; ltsgn++)
      Fprintf(fpout, "%s\n", *ltsgn);

  if (VolLetter == '1')
    VolLetter = 'I';
  else
    VolLetter = toupper(VolLetter);

  show = show_later = NULL;

  for (tsgn = TabSgName; tsgn->HallSymbol; tsgn++)
  {
    if (   WantedSgNumber == -1
        || WantedSgNumber == tsgn->SgNumber)
    {
      if (tsgn->SgNumber >= 3 && tsgn->SgNumber < 16)
      {
        if (VolLetter == 'I')
        {
               ext = tsgn->Extension;
          if (*ext == '-')
               ext++;

          if (       tsgn->Extension[0] == 'b'
              && (   tsgn->Extension[1] == '\0'
                  || tsgn->Extension[1] == '1'))
            show_later = tsgn;
          else if (  ext[0] == 'c')
          {
            if (ext[1] == '\0')
              show = tsgn;
            else
            {
              i = 0;
              for (sgl = tsgn->SgLabels; *sgl; sgl++)
                if (*sgl == '=') i++;

              if (i == 2)
                show = tsgn;
            }
          }
        }
        else if (VolLetter == 'A')
        {
          if (   tsgn->Extension[0] != '-'
              && tsgn->Extension[0] != 'a')
            show = tsgn;
        }
        else
          show = tsgn;
      }
      else if (   tsgn->Extension[0] == 'H'
               && VolLetter == 'I')
        show_later = tsgn;
      else if (   VolLetter == 'A'
               || VolLetter == 'I')
      {
        if (   tsgn->Extension[0] == '\0'
            || tsgn->Extension[1] == '\0')
          show = tsgn;
      }
      else
        show = tsgn;

      if (show)
      {
        putc(' ', fpout);
        PrintTabSgNameEntry(show, 1, 0, 0, fpout);
        putc('\n', fpout);
        show = NULL;

        if (show_later)
        {
          putc(' ', fpout);
          PrintTabSgNameEntry(show_later, 1, 0, 0, fpout);
          putc('\n', fpout);
          show_later = NULL;
        }
      }
    }
  }
}


static const char *LegendTabPgName[] =
  {
    "",
    "  Number                 Hermann-Mauguin             Hall",
    "  ------                 ---------------             ----",
    NULL,
  };


static void ListTabPgName(int WantedPgNumber, FILE *fpout)
{
  const char         **ltpgn;
  const T_TabSgName  *tpgn;


  if (WantedPgNumber == -1)
    for (ltpgn = LegendTabPgName; *ltpgn; ltpgn++)
      Fprintf(fpout, "%s\n", *ltpgn);

  for (tpgn = TabPgName; tpgn->HallSymbol; tpgn++)
  {
    if (   WantedPgNumber == -1
        || WantedPgNumber == abs(tpgn->SgNumber))
    {
      putc(' ', fpout);
      PrintTabSgNameEntry(tpgn, 1, 0, 0, fpout);
      putc('\n', fpout);
    }
  }
}


static void ListCIF(FILE *fpout)
{
  int                n;
  const char         **loop, *lbl;
  const T_TabSgName  *tsgn;

  static const char *loop_monoclinic_extensions[] =
    {
  "_monoclinic_extension   # cf. _symmetry_space_group_id",
  "_monoclinic_axis        # cf. IT Vol. A 1983 sec. 2.16.",
  "_monoclinic_setting     # cf. IT Vol. A 1983 tab. 2.16.1.",
  "_monoclinic_cellchoice  # cf. IT Vol. A 1983 sec. 2.16.(i) & fig. 2.6.4.",
  "",
  " b   b  abc   1",
  " b1  b  abc   1",
  " b2  b  abc   2",
  " b3  b  abc   3",
  "-b   b  c-ba  1",
  "-b1  b  c-ba  1",
  "-b2  b  c-ba  2",
  "-b3  b  c-ba  3",
  " c   c  abc   1",
  " c1  c  abc   1",
  " c2  c  abc   2",
  " c3  c  abc   3",
  "-c   c  ba-c  1",
  "-c1  c  ba-c  1",
  "-c2  c  ba-c  2",
  "-c3  c  ba-c  3",
  " a   a  abc   1",
  " a1  a  abc   1",
  " a2  a  abc   2",
  " a3  a  abc   3",
  "-a   a  -acb  1",
  "-a1  a  -acb  1",
  "-a2  a  -acb  2",
  "-a3  a  -acb  3",
      NULL
    };

  static const char *loop_symmetry_space_group[] =
    {
  "_symmetry_space_group_id",
  "_symmetry_space_group_name_sch",
  "_symmetry_space_group_name_h-m   # recognised IUCr CIF data names",
  "_symmetry_space_group_name_hall  # recognised IUCr CIF data names",
      NULL
    };


  Fprintf(fpout, "data_ notation\n\n");

  Fprintf(fpout, "loop_\n");

  for (loop = loop_monoclinic_extensions; *loop; loop++) {
    if ((*loop)[0]) Fprintf(fpout, "    %s", *loop);
    putc('\n', fpout);
  }

  putc('\n', fpout);
  putc('\n', fpout);

  Fprintf(fpout, "loop_\n");

  for (loop = loop_symmetry_space_group; *loop; loop++) {
    if ((*loop)[0]) Fprintf(fpout, "    %s", *loop);
    putc('\n', fpout);
  }

  putc('\n', fpout);

  for (tsgn = TabSgName; tsgn->HallSymbol; tsgn++)
  {
    n = fprintf(fpout, "    %3d", tsgn->SgNumber);

    if (tsgn->Extension[0])
      n += fprintf(fpout, ":%s", tsgn->Extension);

    if (tsgn->SgNumber < 1 || tsgn->SgNumber > 230) {
      SetSgError("Internal Error: ListCIF()");
      return;
    }

    while (n < 14) { putc(' ', fpout); n++; }
    putc(' ', fpout); n++;

    n += fprintf(fpout, "%s", SchoenfliesSymbols[tsgn->SgNumber]);

    while (n < 22) { putc(' ', fpout); n++; }
    putc(' ', fpout); n++;

    n += PrintFullHM_SgName(tsgn, '_', 0, fpout);

    while (n < 36) { putc(' ', fpout); n++; }
    putc(' ', fpout);

    for (lbl = tsgn->HallSymbol; *lbl; lbl++)
    {
      if (*lbl == ' ' && lbl != tsgn->HallSymbol)
        putc('_', fpout);
      else
        putc(*lbl, fpout);
    }

    putc('\n', fpout);
  }
}


void PutAllXYZ(const T_SgInfo *SgInfo, const char *F_CNS, FILE *fpout)
{
  int           iList, f, i;
  int           nTrV, iTrV, nLoopInv, iLoopInv;
  const int     *TrV;
  T_RTMx        SMx;
  const T_RTMx  *lsmx;
  const char    *xyz, *Separator;
  char          buf0[16], buf1[16], buf2[16];


  if (F_CNS) Separator = ",";
  else       Separator = ", ";

  nLoopInv = Sg_nLoopInv(SgInfo);

  nTrV = SgInfo->LatticeInfo->nTrVector;
   TrV = SgInfo->LatticeInfo->TrVector;

  for (iTrV = 0; iTrV < nTrV; iTrV++, TrV += 3)
  {
    for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
    {
      if (iLoopInv == 0) f =  1;
      else               f = -1;

      lsmx = SgInfo->ListSeitzMx;

      if (! F_CNS && (nLoopInv > 1 || nTrV > 1))
      {
        putc('#', fpout);

        if (nTrV > 1)
          Fprintf(fpout, " +(%s %s %s)",
          FormatFraction(TrV[0], STBF, 0, buf0, sizeof buf0 / sizeof (*buf0)),
          FormatFraction(TrV[1], STBF, 0, buf1, sizeof buf1 / sizeof (*buf1)),
          FormatFraction(TrV[2], STBF, 0, buf2, sizeof buf2 / sizeof (*buf2)));

        if (nLoopInv > 1)
          Fprintf(fpout, " Inversion-Flag = %d", iLoopInv);

        putc('\n', fpout);
      }

      for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
      {
        for (i = 0; i < 9; i++)
          SMx.s.R[i] =              f * lsmx->s.R[i];

        for (i = 0; i < 3; i++)
          SMx.s.T[i] = iModPositive(f * lsmx->s.T[i] + TrV[i], STBF);

            xyz = RTMx2XYZ(&SMx, 1, STBF, 0, 0, 1, Separator, NULL, 0);
        if (xyz)
        {
          if (F_CNS)
            Fprintf(fpout, "%ssymm=(%s)\n", F_CNS, xyz);
          else
            Fprintf(fpout, "%s\n", xyz);
        }
        else
        {
          SetSgError("Internal Error: PutAllXYZ()");
          return;
        }
      }
    }
  }
}


static void PutRmI(const T_SgInfo *SgInfo, int Fmt, FILE *fpout)
{
  int           iList, f, i;
  int           nLoopInv, iLoopInv;
  T_RTMx        RmI;
  const T_RTMx  *lsmx;
  int           iMatrix;


  iMatrix = 0;

  nLoopInv = Sg_nLoopInv(SgInfo);

  for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
  {
    if (iLoopInv == 0) f =  1;
    else               f = -1;

    lsmx = SgInfo->ListSeitzMx;

    if (nLoopInv > 1)
      putc('#', fpout);

    if (nLoopInv > 1)
      Fprintf(fpout, " Inversion-Flag = %d", iLoopInv);

    if (nLoopInv > 1)
      putc('\n', fpout);

    for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
    {
      SetRminusI(lsmx->s.R, RmI.s.R, iLoopInv);

      for (i = 0; i < 3; i++)
        RmI.s.T[i] = f * lsmx->s.T[i];

      Fprintf(fpout, "ri%d", ++iMatrix);
      if (Fmt == 0)
        PrintMapleRTMx(&RmI, 1, STBF, NULL, fpout);
      else
        PrintMathMRTMx(&RmI, 1, STBF, NULL, fpout);
    }
  }

  putc('\n', fpout);
}


static void PutCharacterization(const T_SgInfo *SgInfo, FILE *fpout)
{
  int           iList, iSymTr, f, i;
  int           nTrV, iTrV, nLoopInv, iLoopInv;
  const int     *TrV;
  const T_RTMx  *lsmx;
  T_RTMx        SMx[1];
  int           Mul, CumRMx[9], ScGl[3], Sh[3], ShTrMx[9], Tr[3];
  T_RotMxInfo   RMxI[1];
  const char    *ff;
  char          buf0[16], buf1[16], buf2[16];


  iSymTr = 0;

  nLoopInv = Sg_nLoopInv(SgInfo);

  nTrV = SgInfo->LatticeInfo->nTrVector;
   TrV = SgInfo->LatticeInfo->TrVector;

  for (iTrV = 0; iTrV < nTrV; iTrV++, TrV += 3)
  {
    for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
    {
      if (nLoopInv > 1 || nTrV > 1)
        putc('#', fpout);

      if (nTrV > 1)
        Fprintf(fpout, " +(%s %s %s)",
          FormatFraction(TrV[0], STBF, 0, buf0, sizeof buf0 / sizeof (*buf0)),
          FormatFraction(TrV[1], STBF, 0, buf1, sizeof buf1 / sizeof (*buf1)),
          FormatFraction(TrV[2], STBF, 0, buf2, sizeof buf2 / sizeof (*buf2)));

      if (nLoopInv > 1)
        Fprintf(fpout, " Inversion-Flag = %d", iLoopInv);

      if (nLoopInv > 1 || nTrV > 1)
        putc('\n', fpout);

      if (iLoopInv == 0) f =  1;
      else               f = -1;

      lsmx = SgInfo->ListSeitzMx;

      for (iList = 0; iList < SgInfo->nList; iList++, lsmx++, iSymTr++)
      {
        for (i = 0; i < 9; i++)
          SMx->s.R[i] = f * lsmx->s.R[i];

        for (i = 0; i < 3; i++)
          SMx->s.T[i] = f * lsmx->s.T[i] + TrV[i];

        Mul = BuildCumRMx(SgInfo, iList, iLoopInv, CumRMx);

        RotMx_t_Vector(ScGl, CumRMx, SMx->s.T, 0);

        for (i = 0; i < 3; i++) {
          if (ScGl[i] %  Mul) progerror("Corrupt ScGl");
              ScGl[i] /= Mul;
              ScGl[i] = iModPositive(ScGl[i], STBF);
          if (ScGl[i] > STBF / 2) ScGl[i] -= STBF;
        }

        for (i = 0; i < 3; i++)
          Sh[i] = SMx->s.T[i] - ScGl[i];

        if (GetRotMxInfo(SMx->s.R, RMxI, ShTrMx) == 0)
          progerror("Corrupt SMx");

        RotMx_t_Vector(Tr, ShTrMx, Sh, 0);

        for (i = 0; i < 3; i++) {
          if (Tr[i] %  (STBF * STBF / CTBF)) progerror("Corrupt Tr");
              Tr[i] /= (STBF * STBF / CTBF);
        }

#ifndef NoExtraTest
        {
          int  RmI[9], RS[3];

          SetRminusI(SMx->s.R, RmI, 0);
          RotMx_t_Vector(RS, RmI, Tr, 0);

          for (i = 0; i < 3; i++)
            if (   iModPositive(RS[i], CTBF)
                != iModPositive(Sh[i] * (CTBF / STBF), CTBF))
              progerror("Corrupt ShTrMx");
        }
#endif
        for (i = 0; i < 3; i++)
          Tr[i] = -Tr[i];

#ifdef MarkTabXtalMatrices
        if (    RMxI->Inverse == 0
            && (    RMxI->RefAxis == 'z'
                || (RMxI->RefAxis == 'o'
                    && RMxI->EigenVector[0] == 1
                    && RMxI->EigenVector[1] == 1
                    && RMxI->EigenVector[2] == 1))) putc('@', stdout);
#endif

        Fprintf(fpout, "(%d) [%d %d %d] %d",
          iSymTr + 1,
          RMxI->EigenVector[0],
          RMxI->EigenVector[1],
          RMxI->EigenVector[2],
          RMxI->Order);

        if (RMxI->Inverse)
          Fputs("^-1", stdout);

        Fputs(" (", stdout);

        for (i = 0; i < 3; i++) {
              ff = FormatFraction(ScGl[i], STBF, 0, NULL, 0);
          if (ff == NULL)
            progerror("Internal Error");
          if (i)
            putc(' ', stdout);
          Fputs(ff, stdout);
        }

        Fputs(")", stdout);

        for (i = 0; i < 3; i++) {
              ff = FormatFraction(Tr[i], CTBF, 0, NULL, 0);
          if (ff == NULL)
            progerror("Internal Error");
          putc(' ', stdout);
          Fputs(ff, stdout);
        }

        putc('\n', stdout);
      }
    }
  }

  putc('\n', fpout);
}


static void PutHarkerInfo(T_SgInfo *SgInfo, int Grid[3], FILE *fpout)
{
  int           mHI, nHI, iHI;
  T_HarkerInfo  *HarkerInfo;
  T_SgInfo      PgInfo[1];
  const char    *xyz;
  char          buf[20];


  mHI = MaxHarkerInfo;

      HarkerInfo = malloc(mHI * sizeof (*HarkerInfo));
  if (HarkerInfo == NULL) NotEnoughCore();

  PgInfo->MaxList = 192;

  PgInfo->ListSeitzMx
    = malloc(PgInfo->MaxList * sizeof (*PgInfo->ListSeitzMx));
  if (PgInfo->ListSeitzMx == NULL) NotEnoughCore();

  PgInfo->ListRotMxInfo
    = malloc(PgInfo->MaxList * sizeof (*PgInfo->ListRotMxInfo));
  if (PgInfo->ListRotMxInfo == NULL) NotEnoughCore();

  if (SetPgInfo(SgInfo, PgInfo) != 0)
    progerror(SgError);

  if (PgInfo->TabSgName)
  {
    Fprintf(fpout, "Patterson Space Group  ");
    PrintTabSgNameEntry(PgInfo->TabSgName, 0, 0, 0, fpout);
    putc('\n', fpout);
  }
  else
    Fprintf(fpout, "Patterson Hall Symbol  %s\n", PgInfo->HallSymbol);

  putc('\n', fpout);

      nHI = SetHarkerInfo(SgInfo, PgInfo, HarkerInfo, mHI, Grid);
  if (nHI < 0)
    progerror(SgError);

  Fprintf(fpout, "Unique Harker Operators: %d\n", nHI);

  for (iHI = 0; iHI < nHI; iHI++)
  {
    xyz = RTMx2XYZ(&HarkerInfo[iHI].Ox, 1, STBF, 0, 0, 1, ", ", NULL, 0);

    Fprintf(fpout, " %2d  %-26s  %dD",
      HarkerInfo[iHI].m, xyz,
      HarkerInfo[iHI].Dim);

    if (HarkerInfo[iHI].Dim == 1 || HarkerInfo[iHI].Dim == 2)
    {
      (void) sprintf(buf, "[%d %d %d]",
        HarkerInfo[iHI].N[0],
        HarkerInfo[iHI].N[1],
        HarkerInfo[iHI].N[2]);

      Fprintf(fpout, "  %-8s", buf);

      Fprintf(fpout, "  [%s",
        FormatFraction(HarkerInfo[iHI].Ox.s.T[0], STBF, 0, NULL, 0));
      Fprintf(fpout, " %s",
        FormatFraction(HarkerInfo[iHI].Ox.s.T[1], STBF, 0, NULL, 0));
      Fprintf(fpout, " %s]",
        FormatFraction(HarkerInfo[iHI].Ox.s.T[2], STBF, 0, NULL, 0));
    }

    putc('\n', fpout);
  }

  putc('\n', fpout);

  free(PgInfo->ListSeitzMx);
  free(PgInfo->ListRotMxInfo);
  free(HarkerInfo);
}


static void VerifyTabShTrMx(const T_TabXtalRotMx *txrmx,
                            int iNextBasis,
                            int iLoopInv,
                            int iLoopSgn,
                            const int *mmx, int F_Debug)
{
  const int  *ShTrMx, *NBRMx, *InvNBRMx;
  int        ShTrMxBuf[9];
  int        RmI[9], Tr[3], Sh[3], RT[3], RS[3], i;


  ShTrMx = txrmx->ShTrMx[iLoopInv * 2 + iLoopSgn];

  if (iNextBasis)
  {
    (void) memcpy(ShTrMxBuf, ShTrMx, sizeof ShTrMxBuf);

    (void) nNextBasis_of_DirCode(txrmx->DirCode, &NBRMx, &InvNBRMx);

    while (iNextBasis--)
      RotateRotMx(ShTrMxBuf, NBRMx, InvNBRMx);

    ShTrMx = ShTrMxBuf;
  }

  SetRminusI(mmx, RmI, 0);

  for (Tr[0] = -4; Tr[0] <= 4; Tr[0]++)
  for (Tr[1] = -4; Tr[1] <= 4; Tr[1]++)
  for (Tr[2] = -4; Tr[2] <= 4; Tr[2]++)
  {
    RotMx_t_Vector(Sh, RmI,    Tr, 0);
    RotMx_t_Vector(RT, ShTrMx, Sh, 0);
    RotMx_t_Vector(RS, RmI,    RT, 0);

    if (F_Debug)
      Fprintf(stdout, "[%d %d %d] -> [%d %d %d] -> [%d %d %d]\n",
        Tr[0] * STBF, Tr[1] * STBF, Tr[2] * STBF,
        Sh[0], Sh[1], Sh[2],
        RT[0], RT[1], RT[2]);

    for (i = 0; i < 3; i++)
      if (RS[i] != Sh[i] * STBF) progerror("Corrupt ShTrMx");
  }
}


static void PrintMathRotMxInfo(const T_RotMxInfo *RMxI, FILE *fpout)
{
  Fprintf(fpout, "Print[\"[%d %d %d] %d",
    RMxI->EigenVector[0],
    RMxI->EigenVector[1],
    RMxI->EigenVector[2],
    RMxI->Order);
  if (RMxI->Inverse) Fputs("^-1", fpout);
  Fprintf(fpout, " '%c' '%c'\"]\n",
     RMxI->RefAxis,
    (RMxI->DirCode != '"' ? RMxI->DirCode : '+'));
}


static int FindShTrMxWithMathematica(int F_Debug)
{
  int          nNextBasis, iNextBasis, i;
  int          nLoopInv, iLoopInv, iLoopSgn;
  int          MatchMx[9], InvMatchMx[9], SgnMatchMx[9], REV[3], *mmx;
  const int    *NBRMx, *InvNBRMx;
  T_RotMxInfo  RotMxInfo[1];

  const T_TabXtalRotMx  *txrmx;


  Fputs("<<Calculus`VectorAnalysis`\n", stdout);
  Fputs("tomx[pol_] := (s0 = 1; s1 = 0; s2 = 0; r0=pol;\n", stdout);
  Fputs("               s0 = 0; s1 = 1; s2 = 0; r1=pol;\n", stdout);
  Fputs("               s0 = 0; s1 = 0; s2 = 1; r2=pol;\n", stdout);
  Fputs("               Clear[s0,s1,s2];\n", stdout);
  Fprintf(stdout,
        "               Transpose[%d*{r0,r1,r2}])\n", STBF);
  Fputs("s={s0,s1,s2}\n", stdout);
  Fputs("t={t0,t1,t2}\n", stdout);
  Fputs("Ghz={{1,-1/2,0},{-1/2,1,0},{0,0,1}}\n", stdout);
  Fputs("Ghx={{1,0,0},{0,1,-1/2},{0,-1/2,1}}\n", stdout);
  Fputs("Ghy={{1,0,-1/2},{0,1,0},{-1/2,0,1}}\n", stdout);
  Fputs("Print[\"--------------------\"]\n", stdout);

  for (txrmx = TabXtalRotMx; txrmx->Order; txrmx++)
  {
    nNextBasis = nNextBasis_of_DirCode(txrmx->DirCode, &NBRMx, &InvNBRMx);

    if (nNextBasis < 0)
      return -1;

    if (txrmx->Order > 2) nLoopInv = 2;
    else                  nLoopInv = 1;

    for (i = 0; i < 9; i++) MatchMx[i] = txrmx->RMx[i];

    for (iNextBasis = 0; iNextBasis < nNextBasis; iNextBasis++)
    {
      if (iNextBasis)
        RotateRotMx(MatchMx, NBRMx, InvNBRMx);

      mmx = MatchMx;

      for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
      {
        if (iLoopInv)
        {
          iCoFactorMxTp(MatchMx, InvMatchMx);
          mmx = InvMatchMx;
        }

        for (iLoopSgn = 0; iLoopSgn < 2; iLoopSgn++)
        {
          if (iLoopSgn)
          {
            for (i = 0; i < 9; i++) SgnMatchMx[i] = -mmx[i];
            mmx = SgnMatchMx;

            RotMxInfo->Order = -txrmx->Order;
          }
          else
            RotMxInfo->Order =  txrmx->Order;

          RotMxInfo->Inverse = iLoopInv;

          if (nNextBasis == 3)
          {
            switch(iNextBasis)
            {
              case 0: RotMxInfo->RefAxis = 'z'; break;
              case 1: RotMxInfo->RefAxis = 'x'; break;
              case 2: RotMxInfo->RefAxis = 'y'; break;
            }
          }
          else
            RotMxInfo->RefAxis = 'o';

          RotMxInfo->DirCode = txrmx->DirCode;

          for (i = 0; i < 3; i++)
            RotMxInfo->EigenVector[i] = txrmx->EigenVector[i];

          for (i = iNextBasis; i--;)
          {
            RotMx_t_Vector(REV, NBRMx, RotMxInfo->EigenVector, 0);

            if (i-- == 0)
            {
              for (i = 0; i < 3; i++)
                RotMxInfo->EigenVector[i] = REV[i];

              break;
            }

            RotMx_t_Vector(RotMxInfo->EigenVector, NBRMx, REV, 0);
          }

          PrintMathRotMxInfo(RotMxInfo, stdout);

          Fprintf(stdout, "RI={{%d,%d,%d},{%d,%d,%d},{%d,%d,%d}}\n",
            mmx[0] - 1, mmx[1],     mmx[2],
            mmx[3],     mmx[4] - 1, mmx[5],
            mmx[6],     mmx[7],     mmx[8] - 1);

          Fprintf(stdout, "e={%d,%d,%d}\n",
            RotMxInfo->EigenVector[0],
            RotMxInfo->EigenVector[1],
            RotMxInfo->EigenVector[2]);

          if      (RotMxInfo->Order == 1)
            Fputs("Print[\"t = {0,0,0}\"]\n", stdout);
          else if (RotMxInfo->Order == -2)
          {
            if (txrmx->CarHex & 0x1)
            {
              Fputs("Solve[{RI.t==s, CrossProduct[s,e]==0,"
                                   " CrossProduct[t,e]==0}, t]\n",
                    stdout);
              Fputs("t /. %\n", stdout);
              Fputs("Mc=tomx[%[[1]]]\n", stdout);
              PrintMathRotMxInfo(RotMxInfo, stdout);
              Fputs("Print[\"Car\"]\n", stdout);
              Fputs("MatrixForm[Mc]\n", stdout);
              VerifyTabShTrMx(txrmx, iNextBasis, iLoopInv, iLoopSgn, mmx,
                              F_Debug);
            }
            if (txrmx->CarHex & 0x2)
            {
              if (   RotMxInfo->RefAxis != 'z'
                  && RotMxInfo->RefAxis != 'x'
                  && RotMxInfo->RefAxis != 'y') {
                SetSgError("Corrupt RefAxis");
                return -1;
              }
              Fprintf(stdout,
                      "Solve[{RI.t==s, CrossProduct[Gh%c.s,Gh%c.e]==0,"
                      " CrossProduct[Gh%c.t,Gh%c.e]==0}, t]\n",
                      RotMxInfo->RefAxis, RotMxInfo->RefAxis,
                      RotMxInfo->RefAxis, RotMxInfo->RefAxis);
              Fputs("t /. %\n", stdout);
              Fputs("Mh=tomx[%[[1]]]\n", stdout);
              PrintMathRotMxInfo(RotMxInfo, stdout);
              Fputs("Print[\"Hex\"]\n", stdout);
              Fputs("MatrixForm[Mh]\n", stdout);
              VerifyTabShTrMx(txrmx, iNextBasis, iLoopInv, iLoopSgn, mmx,
                              F_Debug);
            }
            if (txrmx->CarHex & 0x1 && txrmx->CarHex & 0x2)
            {
              Fputs("Print[\"Mc-Mh\"]\n", stdout);
              Fputs("Simplify[Mc-Mh]\n", stdout);
            }
          }
          else if (RotMxInfo->Order < 0)
          {
            Fputs("Solve[RI.t==s, t]\n", stdout);
            Fputs("t /. %\n", stdout);
            Fputs("M=tomx[%[[1]]]\n", stdout);
            PrintMathRotMxInfo(RotMxInfo, stdout);
            Fputs("MatrixForm[M]\n", stdout);
            VerifyTabShTrMx(txrmx, iNextBasis, iLoopInv, iLoopSgn, mmx,
                            F_Debug);
          }
          else
          {
            if (txrmx->CarHex & 0x1)
            {
              Fputs("Solve[{RI.t==s, s.e==0, t.e==0}, t]\n", stdout);
              Fputs("t /. %\n", stdout);
              Fputs("Mc=tomx[%[[1]]]\n", stdout);
              PrintMathRotMxInfo(RotMxInfo, stdout);
              Fputs("Print[\"Car\"]\n", stdout);
              Fputs("MatrixForm[Mc]\n", stdout);
              VerifyTabShTrMx(txrmx, iNextBasis, iLoopInv, iLoopSgn, mmx,
                              F_Debug);
            }
            if (txrmx->CarHex & 0x2)
            {
              if (   RotMxInfo->RefAxis != 'z'
                  && RotMxInfo->RefAxis != 'x'
                  && RotMxInfo->RefAxis != 'y') {
                SetSgError("Corrupt RefAxis");
                return -1;
              }
              Fprintf(stdout,
                      "Solve[{RI.t==s, s.Gh%c.e==0, t.Gh%c.e==0}, t]\n",
                      RotMxInfo->RefAxis, RotMxInfo->RefAxis);
              Fputs("t /. %\n", stdout);
              Fputs("Mh=tomx[%[[1]]]\n", stdout);
              PrintMathRotMxInfo(RotMxInfo, stdout);
              Fputs("Print[\"Hex\"]\n", stdout);
              Fputs("MatrixForm[Mh]\n", stdout);
              VerifyTabShTrMx(txrmx, iNextBasis, iLoopInv, iLoopSgn, mmx,
                              F_Debug);
            }
            if (txrmx->CarHex & 0x1 && txrmx->CarHex & 0x2)
            {
              Fputs("Print[\"Mc-Mh\"]\n", stdout);
              Fputs("Simplify[Mc-Mh]\n", stdout);
            }
          }

          Fputs("Print[\"--------------------\"]\n", stdout);
        }
      }
    }
  }

  return 0;
}


static void PutSymbAlgb(const T_SgInfo *SgInfo, int Fmt, FILE *fpout)
{
  int           iList, f, i;
  int           nTrV, iTrV, nLoopInv, iLoopInv;
  const int     *TrV;
  T_RTMx        SMx;
  const T_RTMx  *lsmx;
  int           iMatrix;
  char          buf0[8], buf1[8], buf2[8];


  iMatrix = 0;

  nLoopInv = Sg_nLoopInv(SgInfo);

  nTrV = SgInfo->LatticeInfo->nTrVector;
   TrV = SgInfo->LatticeInfo->TrVector;

  for (iTrV = 0; iTrV < nTrV; iTrV++, TrV += 3)
  {
    for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
    {
      if (iLoopInv == 0) f =  1;
      else               f = -1;

      lsmx = SgInfo->ListSeitzMx;

      if (nLoopInv > 1 || nTrV > 1)
        putc('#', fpout);

      if (nTrV > 1)
        Fprintf(fpout, " +(%s %s %s)",
          FormatFraction(TrV[0], STBF, 0, buf0, sizeof buf0 / sizeof (*buf0)),
          FormatFraction(TrV[1], STBF, 0, buf1, sizeof buf1 / sizeof (*buf1)),
          FormatFraction(TrV[2], STBF, 0, buf2, sizeof buf2 / sizeof (*buf2)));

      if (nLoopInv > 1)
        Fprintf(fpout, " Inversion-Flag = %d", iLoopInv);

      if (nLoopInv > 1 || nTrV > 1)
        putc('\n', fpout);

      for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
      {
        for (i = 0; i < 9; i++)
          SMx.s.R[i] =              f * lsmx->s.R[i];

        for (i = 0; i < 3; i++)
          SMx.s.T[i] = iModPositive(f * lsmx->s.T[i] + TrV[i], STBF);

        Fprintf(fpout, "m%d", ++iMatrix);
        if (Fmt == 0)
          PrintMapleRTMx(&SMx, 1, STBF, NULL, fpout);
        else
          PrintMathMRTMx(&SMx, 1, STBF, NULL, fpout);
      }
    }
  }

  putc('\n', fpout);
}


static void PutSpaceSymFile(const T_SgInfo *SgInfo, FILE *fpout)
{
  unsigned int       SgID;
  int                iList, SuppressMx, f, i;
  int                nTrV, iTrV, nLoopInv, iLoopInv;
  const int          *TrV;
  const T_RTMx       *lsmx;
  const T_TabSgName  *tsgn;


      tsgn = SgInfo->TabSgName;
  if (tsgn && tsgn->SgLabels == NULL) tsgn = NULL;

  SgID = 0;

  if (tsgn != NULL)
    SgID = SgID_Number(tsgn);

  Fprintf(fpout, "%u '", SgID);

  if (tsgn != NULL)
    (void) PrintFullHM_SgName(tsgn, 0, 0, fpout);
  else if (SgInfo->HallSymbol[0])
    Fprintf(fpout, "%s", SgInfo->HallSymbol);
  else
    Fprintf(fpout, "Unknown");

  putc('\'', fpout);
  putc('\n', fpout);

  nLoopInv = Sg_nLoopInv(SgInfo);

  nTrV = SgInfo->LatticeInfo->nTrVector;
   TrV = SgInfo->LatticeInfo->TrVector;

      iList = SgInfo->OrderL;
  if (iList > 1)
  {
    iList--;
    SuppressMx = 1;
  }
  else
    SuppressMx = 0;

  Fprintf(fpout, "%d\n", iList);

  for (iTrV = 0; iTrV < nTrV; iTrV++, TrV += 3)
  {
    for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
    {
      if (iLoopInv == 0) f =  1;
      else               f = -1;

      lsmx = SgInfo->ListSeitzMx;

      for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
      {
        if (SuppressMx == 0)
        {
          for (i = 0; i < 3; i++)
            Fprintf(fpout, " %12.8f %12.8f %12.8f %12.8f\n",
              (double) f * lsmx->s.R[3 * i + 0],
              (double) f * lsmx->s.R[3 * i + 1],
              (double) f * lsmx->s.R[3 * i + 2],
              (double) iModPositive(f * lsmx->s.T[i] + TrV[i], STBF) / STBF);

          putc(':',  fpout);
          putc('\n', fpout);
        }

        SuppressMx = 0;
      }
    }
  }
}


static void PutShelx(const T_SgInfo *SgInfo, FILE *fpout)
{
  int           Latt_N = 0, iList;
  const T_RTMx  *lsmx;
  const char    *xyz;


  if (SgInfo->InversionOffOrigin != 0)
    Fprintf(fpout, "***WARNING***: %s\n",
      "Shelx manual: the origin MUST lie on a center of symmetry");

  switch (SgInfo->LatticeInfo->Code)
  {
    case 'P': Latt_N = 1; break;
    case 'A': Latt_N = 5; break;
    case 'B': Latt_N = 6; break;
    case 'C': Latt_N = 7; break;
    case 'I': Latt_N = 2; break;
    case 'R':
              if (SgInfo->ExtraInfo == EI_Obverse)
                Latt_N = 3; break;
    case 'S':
    case 'T':
              SetSgError("Shelx supports R-obverse only");
              return;
    case 'H': SetSgError("Shelx does not support H-centred cells"); return;
    case 'K': SetSgError("Shelx does not support K-centred cells"); return;
    case 'L': SetSgError("Shelx does not support L-centred cells"); return;
    case 'F': Latt_N = 4; break;
    default:
      goto ReturnError;
  }

  /* N must be made negative if the structure is non-centrosymmetric
   */
  if (SgInfo->Centric != -1)
    Latt_N = -Latt_N;

  Fprintf(fpout, "LATT %2d\n", Latt_N);

  lsmx = &SgInfo->ListSeitzMx[1]; /* skip first = identity matrix */

  for (iList = 1; iList < SgInfo->nList; iList++, lsmx++)
  {
        xyz = RTMx2XYZ(lsmx, 1, STBF, 1, 1, 0, ", ", NULL, 0);
    if (xyz)
      Fprintf(fpout, "SYMM %s\n", xyz);
    else
      goto ReturnError;
  }

  putc('\n', fpout);

  return;

  ReturnError:

  SetSgError("Internal Error: PutShelx()");
  return;
}


static void PutSchakal(const T_SgInfo *SgInfo, FILE *fpout)
{
  int           iList, nMx, i;
  int           nTrV, iTrV;
  const int     *TrV;
  T_RTMx        SMx;
  const T_RTMx  *lsmx;
  const char    *xyz;


  if (Sg_nLoopInv(SgInfo) == 2)
    Fprintf(fpout, "DU -x,-y,-z\n");

  nTrV = SgInfo->LatticeInfo->nTrVector;
   TrV = SgInfo->LatticeInfo->TrVector;

  if (nTrV > 1)
  {
    Fprintf(fpout, "DU");

    InitRotMx(SMx.s.R, 1);

    TrV += 3;

    for (iTrV = 1; iTrV < nTrV; iTrV++, TrV += 3)
    {
      for (i = 0; i < 3; i++)
        SMx.s.T[i] = TrV[i];

          xyz = RTMx2XYZ(&SMx, 1, STBF, 0, 0, 1, ",", NULL, 0);
      if (xyz)
      {
        if (iTrV > 1)
          Fprintf(fpout, " ;");

        Fprintf(fpout, " %s", xyz);
      }
      else
      {
        putc('\n', fpout);
        goto ReturnError;
      }
    }

    putc('\n', fpout);
  }

  nMx = 0;

  lsmx = &SgInfo->ListSeitzMx[1];

  for (iList = 1; iList < SgInfo->nList; iList++, lsmx++)
  {
        xyz = RTMx2XYZ(lsmx, 1, STBF, 0, 0, 1, ",", NULL, 0);
    if (xyz)
    {
      if (nMx % 4 == 0)
      {
        if (nMx) putc('\n', fpout);
        Fprintf(fpout, "SY %s", xyz);
      }
      else
        Fprintf(fpout, " ; %s", xyz);
    }
    else
    {
      putc('\n', fpout);
      goto ReturnError;
    }

    nMx++;
  }

  if (nMx)
    putc('\n', fpout);

  putc('\n', fpout);

  return;

  ReturnError:

  SetSgError("Internal Error: PutSchakal()");
  return;
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


static void IllegalLine(const int lcount, const char *msg)
{
  Fflush(stdout);
  Fprintf(stderr, "%s: Illegal input line #%d", progn, lcount);
  if (msg) Fprintf(stderr, ": %s", msg);
  putc('\n', stderr);
  exit(1);
}


static char firstnonblank(const char *s)
{
  while (*s && isspace(*s)) s++;
  return *s;
}


static void ReadXYZ(FILE *fpin, T_SgInfo *SgInfo)
{
  int     lcount, c, i, n;
  char    buf[256];
  T_RTMx  SMx[1];


  lcount = 0;

  while (fgetline(fpin, buf, sizeof buf))
  {
    lcount++;

        c = firstnonblank(buf);
    if (c == '\0' || strchr("!#$%", c) != NULL)
      continue;

    for (i = 0; buf[i];)
    {
      if (isspace(buf[i])) {
        n = 0;
        while (buf[i] && isspace(buf[i])) { n++; i++; }
        if (n >= 2 && buf[i] && buf[i] != ',') buf[i - 1] = ',';
      }

      if (buf[i] == ',') {
        i++;
        while (buf[i] && isspace(buf[i])) i++;
      }
      else
        i++;
    }

    if (ParseStrXYZ(buf, SMx, 1, STBF) != 0)
      IllegalLine(lcount, NULL);

    Fprintf(stdout, "SymMx  %s\n",
      RTMx2XYZ(SMx, 1, STBF, 0, 0, 1, ", ", NULL, 0));

    if (Add2ListSeitzMx(SgInfo, SMx) < 0)
      IllegalLine(lcount, SgError);
  }

  if (lcount)
    putc('\n', stdout);
}


static void Simple_hklList(T_SgInfo *SgInfo, int FriedelSym,
                           int Maxh, int Maxk, int Maxl,
                           int ListSysAbsent)
{
  int        h, k, l, iList, restriction, M, n, i;
  int        Minh, Mink, Minl;
  int        uvw[3];
  int        CCMx_PL[9], deterCCMx_LP = 0, hP, kP, lP;


  if (SgInfo->LatticeInfo->Code != 'P')
  {
    deterCCMx_LP = deterRotMx(SgInfo->CCMx_LP);
                iCoFactorMxTp(SgInfo->CCMx_LP, CCMx_PL);

    if (deterCCMx_LP < 1)
      goto ReturnError;
  }

  if (SetListMin_hkl(SgInfo, FriedelSym,
                      Maxh,  Maxk,  Maxl,
                     &Minh, &Mink, &Minl) != 0)
    return;

  Fprintf(stdout, ">Begin hklList\n");

  for (h = Minh; h <= Maxh; h++)
  for (k = Mink; k <= Maxk; k++)
  for (l = Minl; l <= Maxl; l++)
  {
    iList = IsSysAbsent_hkl(SgInfo, h, k, l, &restriction);
    if (SgError != NULL)
      return;

    M = Mult_hkl(SgInfo, FriedelSym, h, k, l);
    if (SgError != NULL)
      return;

    if (iList == 0)
    {
      if ((iList = IsHidden_hkl(SgInfo, FriedelSym,
                                Minh, Mink, Minl,
                                Maxh, Maxk, Maxl,
                                   h,    k,    l)) != 0)
        n = fprintf(stdout, "# %3d %3d %3d  %3d  [%d]",
                            h, k, l, M, iList);
      else
        n = fprintf(stdout, "  %3d %3d %3d  %3d",
                            h, k, l, M);

      if (restriction >= 0)
      {
        while (n < 27) { n++; putc(' ', stdout); }
        n += fprintf(stdout, " %2d/%d", restriction, STBF);
      }

      while (n < 34) { n++; putc(' ', stdout); }
      if (Is_ss(SgInfo, h, k, l) == 1)
        n += fprintf(stdout, " s.s.");

      while (n < 41) { n++; putc(' ', stdout); }
      Set_uvw(SgInfo, h, k, l, uvw);
      for (i = 0; i < SgInfo->n_ssVM; i++)
        n += fprintf(stdout, " %3d", uvw[i]);

      if (SgInfo->LatticeInfo->Code != 'P')
      {
        hP = h * CCMx_PL[0] + k * CCMx_PL[3] + l * CCMx_PL[6];
        kP = h * CCMx_PL[1] + k * CCMx_PL[4] + l * CCMx_PL[7];
        lP = h * CCMx_PL[2] + k * CCMx_PL[5] + l * CCMx_PL[8];

        if (hP % deterCCMx_LP || kP % deterCCMx_LP || lP % deterCCMx_LP)
          goto ReturnError;

        hP /= deterCCMx_LP;
        kP /= deterCCMx_LP;
        lP /= deterCCMx_LP;

        while (n < 55) { n++; putc(' ', stdout); }
          n += fprintf(stdout, " P  %3d %3d %3d",
                               hP, kP, lP);
      }

      putc('\n', stdout);
    }
    else if (ListSysAbsent)
      Fprintf(stdout, "# %3d %3d %3d  %3d  (%d)\n",
                      h, k, l, M, iList);
  }

  Fprintf(stdout, ">End hklList\n");

  return;

  ReturnError:

  SetSgError("Internal Error: Simple_hklList()");
  return;
}


/* ****************************************************************************
   some code for the handling of lattice constants
 */


#include <math.h>


typedef struct {
                 double   a, b, c;
                 double   alpha, beta, gamma;
                 double   sa, sb, sg;
                 double   ca, cb, cg;
                 double   v;
                 char     calcs, calcc;
               }
               T_LatticeConstants;


#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2          1.57079632679489661923
#endif

#define PIover180 (M_PI / 180.)

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


static int Lc2RLc(T_LatticeConstants *lc, T_LatticeConstants *rlc)
{
  /* Transformation Lattice Constants -> Reciprocal Lattice Constants
     after Kleber, W., 17. Aufl., Verlag Technik GmbH Berlin 1990, P.352
   */

  double  D;


  if (lc->calcs)
  { lc->sa = sinC(lc->alpha); lc->sb = sinC(lc->beta); lc->sg = sinC(lc->gamma);
    lc->calcs = 0;
  }

  if (lc->calcc)
  { lc->ca = cosC(lc->alpha); lc->cb = cosC(lc->beta); lc->cg = cosC(lc->gamma);
    lc->calcc = 0;
  }

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
    rlc->calcs = 0;
    rlc->calcc = 0;
  }

  return 0;
}


static void Lc2MetricalMx(T_LatticeConstants *lc, double *G)
{
  G[0] =        lc->a * lc->a;
  G[1] = G[3] = lc->a * lc->b * lc->cg;
  G[2] = G[6] = lc->a * lc->c * lc->cb;

  G[4] =        lc->b * lc->b;
  G[5] = G[7] = lc->b * lc->c * lc->ca;

  G[8] =        lc->c * lc->c;
}


static int HarmonizeSgLatCon(T_SgInfo *SgInfo, T_LatticeConstants *lc, int np)
{
  switch(SgInfo->XtalSystem)
  {
    case XS_Triclinic:
      if (np != 6) goto IllUnitCell;
      break;
    case XS_Monoclinic:
      if (np != 4 && np != 6) goto IllUnitCell;
      switch (SgInfo->UniqueRefAxis)
      {
        case 'x': lc->beta  = lc->gamma = 90. * PIover180; break;
        case 'y': if (np != 6) lc->beta  = lc->alpha;
                  lc->alpha = lc->gamma = 90. * PIover180; break;
        case 'z': if (np != 6) lc->gamma = lc->alpha;
                  lc->alpha = lc->beta  = 90. * PIover180; break;
        default:
          goto IntErr;
      }
      break;
    case XS_Orthorhombic:
      if (np != 3 && np != 6) goto IllUnitCell;
      lc->alpha = lc->beta = lc->gamma = 90. * PIover180;
      break;
    case XS_Tetragonal:
      if (np != 2 && np != 6) goto IllUnitCell;
      switch (SgInfo->UniqueRefAxis)
      {
        case 'x': lc->c = lc->b; break;
        case 'y': lc->c = lc->a; break;
        case 'z': if (np != 6) lc->c = lc->b;
                  lc->b = lc->a; break;
        default:
          goto IntErr;
      }
      lc->alpha = lc->beta = lc->gamma = 90. * PIover180;
      break;
    case XS_Trigonal:
      if (np != 2 && np != 6) goto IllUnitCell;
      if (SgInfo->UniqueDirCode == '*')
      {
        if (np != 6) lc->alpha = lc->b * PIover180;
        lc->c = lc->b = lc->a;
        lc->gamma = lc->beta = lc->alpha;
        break;
      }
    case XS_Hexagonal:
      if (np != 2 && np != 6) goto IllUnitCell;
      switch (SgInfo->UniqueRefAxis)
      {
        case 'x': lc->c = lc->b;
                  lc->alpha = 120. * PIover180;
                  lc->beta  = lc->gamma = 90. * PIover180; break;
        case 'y': lc->c = lc->a;
                  lc->beta  = 120. * PIover180;
                  lc->alpha = lc->gamma = 90. * PIover180; break;
        case 'z': if (np != 6) lc->c = lc->b;
                  lc->b = lc->a;
                  lc->gamma = 120. * PIover180;
                  lc->alpha = lc->beta  = 90. * PIover180; break;
        default:
          goto IntErr;
      }
      break;
    case XS_Cubic:
      if (np != 1 && np != 6) goto IllUnitCell;
      lc->c = lc->b = lc->a;
      lc->alpha = lc->beta = lc->gamma = 90. * PIover180;
      break;
    default:
      goto IntErr;
  }

  return  0;

  IntErr: SetSgError("Internal Error: HarmonizeSgLatCon()");
  return -1;

  IllUnitCell: SetSgError("Error: Illegal UnitCell or SpaceGroup");
  return -1;
}


static void MxMultiply(double *ab, double *a, double *b, int ma, int na, int nb)
{
  int     i, j, k;
  double  *ai, *aij, *bk, *bkj;

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


static int TransformLatticeConstants(T_LatticeConstants *LatConA,
                                     int np,
                                     T_LatticeConstants *LatConB,
                                     T_SgInfo *SgInfo,
                                     int *InvCBMxR)
{
  int     i, j;
  double  GA[9], GB[9], GAR[9], R[9], Rt[9];


  if (HarmonizeSgLatCon(SgInfo, LatConA, np) != 0)
    return -1;

  LatConA->calcs = 1;
  LatConA->calcc = 1;

  /* just to check LatConA and to compute sin and cos of angles
   */
  if (Lc2RLc(LatConA, LatConB) != 0) {
    SetSgError("Error: Illegal UnitCell");
    return -1;
  }

  Lc2MetricalMx(LatConA, GA);

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) {
       R[i * 3 + j] = InvCBMxR[i * 3 + j] / (double) CRBF;
      Rt[i * 3 + j] = InvCBMxR[j * 3 + i] / (double) CRBF;
    }

  MxMultiply(GAR, GA, R, 3, 3, 3);
  MxMultiply(GB, Rt, GAR, 3, 3, 3);

  if (GB[0] < 0. || GB[4] < 0. || GB[8] < 0.)
    goto ReturnError;

  LatConB->a = sqrt(GB[0]);
  LatConB->b = sqrt(GB[4]);
  LatConB->c = sqrt(GB[8]);

  LatConB->alpha = GB[5] / LatConB->b / LatConB->c;
  LatConB->beta  = GB[2] / LatConB->c / LatConB->a;
  LatConB->gamma = GB[1] / LatConB->a / LatConB->b;

  if (   LatConB->alpha < -1. || LatConB->alpha > 1.
      || LatConB->beta  < -1. || LatConB->beta  > 1.
      || LatConB->gamma < -1. || LatConB->gamma > 1.)
    goto ReturnError;

  LatConB->alpha = acos(LatConB->alpha);
  LatConB->beta  = acos(LatConB->beta );
  LatConB->gamma = acos(LatConB->gamma);

  LatConB->calcs = 1;
  LatConB->calcc = 1;

  return 0;

  ReturnError:

  SetSgError("InternalError: Corrupt InvCBMxR");
  return -1;
}


/* ****************************************************************************
 */


static void usage(int F_Verbose)
{
  static const char *quick_help[] =
    {
  "-Hall|VolA|VolI|Plane  select conventions",
  "-ListTable[=#]         print [parts of] internal table",
  "-CIF                   print internal table in CIF format",
  "-ReadXYZ               read \"xyz\" symmetry operations from stdin",
  "-XYZ                   print something like \"-x, y+1/2, z\"",
  "-AllXYZ                print all symmetry operations",
  "-Characterize          print translation components and positions",
  "-Harker                print unique Harker operators",
  "-Shelx                 print Shelx LATT & SYMM cards",
  "-CNS                   print CNS symmetry record",
  "-Schakal               print Schakal DU & SY cards",
  "-Space                 print symmetry file for AVS SpaceModule",
  "-Standard              compute transformation to \"standard\" setting",
"-UnitCell=\"a..g\"       unit cell constants a, b, c, alpha, beta, gamma",
  "-v                     be more verbose (also lists more options)",
  NULL
    };

  static const char *more_help[] =
    {
  "-CNSlib                generate CNS spacegroup.lib",
  "-Maple                 print symmetry matrices in Maple format",
  "-Mathematica           print symmetry matrices in Mathematica format",
"-CBMx=\"x,y,z\"          transform space group with CBMx",
  "-ShowPrimitive         show primitive symmetry operations",
  "-hklList[=#]           print simple hkl listing",
  "-Conditions            generate hkl conditions",
  "-FindShTrMx            generate Mathematica input",
  "-Xtal32ss              use xtal32ss to set s.s. vectors",
  "-Verify                verify s.s. vectors & c.o.b. matrices & sghkl",
  "-ClearError            clear errors and continue",
  "-Debug                 debug option",
  NULL
    };

  const char  **hlp;


  if (F_Verbose)
    Fprintf(stderr, "\n%s\n\n", SgInfoVersion());

  Fprintf(stderr,
    "usage: %s [options] [SpaceGroupName_or_# [SpaceGroupName_or_#]]\n",
    progn);

  for (hlp = quick_help; *hlp; hlp++)
    Fprintf(stderr, "  %s\n", *hlp);

  if (F_Verbose)
    for (hlp = more_help; *hlp; hlp++)
      Fprintf(stderr, "  %s\n", *hlp);

  Fprintf(stderr, "\n");
  Fprintf(stderr, "examples: %s 68\n",                          progn);
  Fprintf(stderr, "          %s C2/m:c2 -XYZ\n",                progn);
  Fprintf(stderr, "          %s \"Oh^3\" -Shelx\n",             progn);
  Fprintf(stderr, "          %s -Hall \"-F 4y 2\" -Standard\n", progn);
  Fprintf(stderr, "          %s -VolI 15 -VolA 15\n",           progn);
  Fprintf(stderr, "          %s -ListTable=68\n",               progn);

  if (F_Verbose)
    Fprintf(stderr, "\n");

  exit(1);
}


static void ShowCBMx(T_RTMx *CBMx, T_RTMx *InvCBMx, int F_SymbAlgb)
{
  if (F_SymbAlgb % 2) {
    PrintMapleRTMx(   CBMx, CRBF, CTBF, "   CBMx", stdout);
    PrintMapleRTMx(InvCBMx, CRBF, CTBF, "InvCBMx", stdout);
  }
  if (F_SymbAlgb / 2) {
    PrintMathMRTMx(   CBMx, CRBF, CTBF, "   CBMx", stdout);
    PrintMathMRTMx(InvCBMx, CRBF, CTBF, "InvCBMx", stdout);
  }
  if (F_SymbAlgb == 0) {
    Fprintf(stdout, "   CBMx = %s\n",
      RTMx2XYZ(   CBMx, CRBF, CTBF, 0, 0, 1, ", ", NULL, 0));
    Fprintf(stdout, "InvCBMx = %s\n",
      RTMx2XYZ(InvCBMx, CRBF, CTBF, 0, 0, 1, ", ", NULL, 0));
  }
}


static int ShowTransformedSeitzMatrices(T_SgInfo *SgInfo,
                                        T_RTMx *CBMx, T_RTMx *InvCBMx)
{
  int     i;
  T_RTMx  BufCBMx[1];

  int           iList, f;
  int           nLoopInv, iLoopInv;
  T_RTMx        SMx, BC_SMx;
  const T_RTMx  *lsmx;
  T_SMxI        SI[1];


  if (CBMx == NULL && InvCBMx == NULL) goto ReturnError;

  if      (InvCBMx == NULL) {
           InvCBMx = BufCBMx;
    if (InverseRTMx(CBMx, InvCBMx, CRBF) == 0) goto ReturnError;
  }
  else if (   CBMx == NULL) {
              CBMx = BufCBMx;
    if (InverseRTMx(InvCBMx, CBMx, CRBF) == 0) goto ReturnError;
  }

  ShowCBMx(CBMx, InvCBMx, 0);
  ShowCBMx(CBMx, InvCBMx, 2);
  putc('\n', stdout);

  nLoopInv = Sg_nLoopInv(SgInfo);

  for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
  {
    if (iLoopInv == 0) f =  1;
    else               f = -1;

    lsmx = SgInfo->ListSeitzMx;

    for (iList = 0; iList < SgInfo->nList; iList++, lsmx++)
    {
      for (i = 0; i < 9; i++)
        SMx.s.R[i] = f * lsmx->s.R[i];

      for (i = 0; i < 3; i++)
        SMx.s.T[i] = f * lsmx->s.T[i];

      if (CB_SMx(&BC_SMx, CBMx, &SMx, InvCBMx) != 0)
        goto ReturnError;

      ResetSeitzMxInfo(SI, 0);

      if (SetSeitzMxInfo(&BC_SMx, SI) == 0)
      {
        if (SgError != NULL) return -1;
        goto ReturnError;
      }

      if (GetRotMxInfo(BC_SMx.s.R, NULL, NULL) != 0)
        putc(' ', stdout);
      else
        putc('@', stdout);

      Fprintf(stdout, " %2d", SI->Order);
      Fprintf(stdout, " [ %2d %2d %2d ] ",
        SI->Basis[2][0], SI->Basis[2][1], SI->Basis[2][2]);
      if (SI->SenseOfRotation < 0)
        putc('c', stdout);
      else
        putc(' ', stdout);
      Fprintf(stdout, " wg = ");
      for (i = 0; i < 3; i++)
        Fprintf(stdout, " %6s", FormatFraction(SI->wg[i], STBF, 0, NULL, 0));
      Fprintf(stdout, " Tr = ");
      for (i = 0; i < 3; i++)
        Fprintf(stdout, " %6s", FormatFraction(SI->Tr[i], CTBF, 0, NULL, 0));
      putc('\n', stdout);
#define R BC_SMx.s.R
      Fprintf(stdout, "    {{ %2d, %2d, %2d },\n", R[0], R[1], R[2]);
      Fprintf(stdout, "     { %2d, %2d, %2d },\n", R[3], R[4], R[5]);
      Fprintf(stdout, "     { %2d, %2d, %2d }}\n", R[6], R[7], R[8]);
#undef R
    }
  }

  putc('\n', stdout);

  return 0;

  ReturnError:

  SetSgError("Internal Error: ShowTransformedSeitzMatrices()");
  return -1;
}


static int ShowPrimitiveSeitzMatrices(T_SgInfo *SgInfo)
{
  int     i;
  T_RTMx  CBMx[1];


  if (SgInfo->LatticeInfo->Code == 'P') return 0;

  for (i = 0; i < 9; i++) CBMx->s.R[i] = CRBF * SgInfo->CCMx_LP[i];
  for (i = 0; i < 3; i++) CBMx->s.T[i] = 0;

  return ShowTransformedSeitzMatrices(SgInfo, CBMx, NULL);
}


static void FreeSgInfo(T_SgInfo *SgInfo)
{
#define FreeReset(p) if (p) { free(p); p = NULL; }

  FreeReset(SgInfo->ListSeitzMx)
  FreeReset(SgInfo->ListRotMxInfo)
  FreeReset(SgInfo->ReflCond)
  FreeReset(SgInfo->RestCond)
  FreeReset(SgInfo->SysEnhanced)

#undef FreeReset
}


static int AllocSgInfo(T_SgInfo *SgInfo, int F_Conditions)
{
  SgInfo->MaxList = 192;

  SgInfo->ListSeitzMx = NULL;
  SgInfo->ListRotMxInfo = NULL;
  SgInfo->ReflCond = NULL;
  SgInfo->RestCond = NULL;
  SgInfo->SysEnhanced = NULL;

  SgInfo->ListSeitzMx
    = malloc(SgInfo->MaxList * sizeof (*SgInfo->ListSeitzMx));
  if (SgInfo->ListSeitzMx == NULL) goto ReturnError;

#ifndef No_ListRotMxInfo
  SgInfo->ListRotMxInfo
    = malloc(SgInfo->MaxList * sizeof (*SgInfo->ListRotMxInfo));
  if (SgInfo->ListRotMxInfo == NULL) goto ReturnError;
#endif

  if (F_Conditions)
  {
    SgInfo->ReflCond
      = malloc(SgInfo->MaxList * sizeof (*SgInfo->ReflCond));
    if (SgInfo->ReflCond == NULL) goto ReturnError;

    SgInfo->RestCond
      = malloc(SgInfo->MaxList * sizeof (*SgInfo->RestCond));
    if (SgInfo->RestCond == NULL) goto ReturnError;

    SgInfo->SysEnhanced
      = malloc(SgInfo->MaxList * sizeof (*SgInfo->SysEnhanced));
    if (SgInfo->SysEnhanced == NULL) goto ReturnError;
  }

  return 0;

  ReturnError:
  SetSgError("Not enough core");
  FreeSgInfo(SgInfo);
  return -1;
}


typedef struct
  {
    int                Convention;
    const char         *SgName;
    const T_TabSgName  *InpTSgN;
    const T_TabSgName  *RefTSgN;
    T_RTMx             CBMx, InvCBMx;
  }
  T_SgList;


int main(int argc, char *argv[])
{
  int                 i, n, HaveSpace, pos_hsym;
  int                 F_Convention, Last_F_Convention;
  int                 F_ListTable, F_CIF;
  int                 F_XYZ, F_AllXYZ, F_Maple, F_MathM;
  int                 F_Space, F_Shelx, F_Schakal, F_CNS, F_CNSlib;
  int                 F_hklList, F_Conditions;
  int                 F_Standard, F_UnitCell;
  int                 F_Verbose, F_Verify, F_ClearError, F_Debug;
  int                 F_Characterize, F_FindShTrMx, F_Harker;
  int                 F_Xtal32ss, F_ReadXYZ;
  T_LatticeConstants  LatConA, LatConB;
  char                *cp, xtrac;
  const char          *SgName;
  const T_TabSgName   *tsgn;
  T_SgInfo            SpgrInfo[2], BC_SgInfo, *SgInfo;
  int                 nSgList, iSgList;
  T_SgList             SgList[2];
  T_RTMx              *CBMx, *InvCBMx;
  T_RTMx              CCBMx, CInvCBMx;
  int                 Grid[4][3];
  int                 F_CBMx, F_ShowPrimitive;
  T_RTMx              CmdLnCBMx[2];


/*
  Macintosh extras (Courtesy Jon Tischler <TischlerJZ@ornl.gov>)
 */
#ifdef __THINK__
  console_options.nrows = CONSOLE_LINES;
  console_options.ncols = CONSOLE_COLUMNS;
  console_options.title = "\pSgInfo";
#endif
#ifdef __MWERKS__
  SIOUXSettings.autocloseonquit = FALSE;
  SIOUXSettings.asktosaveonclose = TRUE;
  SIOUXSettings.columns = CONSOLE_COLUMNS;
  SIOUXSettings.rows = CONSOLE_LINES;
#endif
#if defined(__THINK__) || defined(__MWERKS__)
  argc = ccommand(&argv);
#endif


  nSgList = 0;

  F_Convention = 'A'; Last_F_Convention = 0;
  F_ListTable = 0;
  F_CIF = 0;
  F_XYZ = 0;
  F_AllXYZ = 0;
  F_Maple = 0;
  F_MathM = 0;
  F_Space = 0;
  F_Shelx = 0;
  F_Schakal = 0;
  F_CNS = 0;
  F_CNSlib = 0;
  F_hklList = 0;
  F_Conditions = 0;
  F_Standard = 0;
  F_UnitCell = 0;
  F_Verbose = 0;
  F_Verify = 0;
  F_ClearError = 0;
  F_Debug = 0;
  F_Characterize = 0;
  F_FindShTrMx = 0;
  F_Harker = 0;
  F_Xtal32ss = 0;
  F_ReadXYZ = 0;
  F_CBMx = 0;
  F_ShowPrimitive = 0;

  for (i = 1; i < argc; i++)
  {
    if      (str_icmp(argv[i], "-Hall") == 0) {
      F_Convention = 'H';
      Last_F_Convention = 0;
    }
    else if (str_icmp(argv[i], "-VolA") == 0) {
      F_Convention = 'A';
      Last_F_Convention = 'A';
    }
    else if (   str_icmp(argv[i], "-VolI") == 0
             || str_icmp(argv[i], "-Vol1") == 0) {
      F_Convention = 'I';
      Last_F_Convention = 'I';
    }
    else if (   str_icmp(argv[i], "-Plane") == 0) {
      F_Convention = 'P';
      Last_F_Convention = 0;
    }
    else if (str_ibegin(argv[i], "-ListTable") == 0)
    {
                cp = argv[i] + 10;
      if      (*cp == '\0')
        F_ListTable = -1;
      else if (*cp++ == '=')
      {
        n = sscanf(cp, "%d %c", &F_ListTable, &xtrac);
        if (n != 1 || F_ListTable <   1
                   || F_ListTable > 230) usage(0);
      }
      else
        usage(0);
    }
    else if (str_icmp(argv[i], "-CIF") == 0)
      F_CIF = 1;

    else if (str_icmp(argv[i], "-XYZ") == 0)
      F_XYZ = 1;

    else if (str_icmp(argv[i], "-AllXYZ") == 0)
      F_AllXYZ = 1;

    else if (str_icmp(argv[i], "-Maple") == 0)
      F_Maple = 1;

    else if (str_icmp(argv[i], "-Mathematica") == 0)
      F_MathM = 2;

    else if (str_icmp(argv[i], "-Space") == 0)
      F_Space = 1;

    else if (str_icmp(argv[i], "-Shelx") == 0)
      F_Shelx = 1;

    else if (str_icmp(argv[i], "-Schakal") == 0)
      F_Schakal = 1;

    else if (str_icmp(argv[i], "-CNS") == 0)
      F_CNS = 1;

    else if (str_icmp(argv[i], "-CNSlib") == 0)
      F_CNSlib = 1;

    else if (str_ibegin(argv[i], "-hklList") == 0)
    {
           cp = &argv[i][8];
      if (*cp == '\0')
        F_hklList = 4;
      else
      {
        if (*cp++ != '=') usage(1);

            n = sscanf(cp, "%d %c", &F_hklList, &xtrac);
        if (n != 1)
          usage(1);
      }
    }

    else if (str_icmp(argv[i], "-Conditions") == 0)
      F_Conditions = 1;

    else if (str_icmp(argv[i], "-Standard") == 0)
      F_Standard = 1;

    else if (str_ibegin(argv[i], "-UnitCell=") == 0)
    {
      F_UnitCell = sscanf(&argv[i][10], "%lf%lf%lf%lf%lf%lf",
        &LatConA.a,     &LatConA.b,    &LatConA.c,
        &LatConA.alpha, &LatConA.beta, &LatConA.gamma);

      if (F_UnitCell < 1)
        usage(0);

      if (F_UnitCell > 3) LatConA.alpha *= PIover180;
      if (F_UnitCell > 4) LatConA.beta  *= PIover180;
      if (F_UnitCell > 5) LatConA.gamma *= PIover180;
    }
    else if (str_icmp(argv[i], "-v") == 0)
      F_Verbose = 1;

    else if (str_icmp(argv[i], "-Verify") == 0)
      F_Verify = 1;

    else if (str_icmp(argv[i], "-ClearError") == 0)
      F_ClearError = 1;

    else if (str_icmp(argv[i], "-Debug") == 0)
      F_Debug = 1;

    else if (str_icmp(argv[i], "-Characterize") == 0)
      F_Characterize = 1;

    else if (str_icmp(argv[i], "-FindShTrMx") == 0)
      F_FindShTrMx = 1;

    else if (str_icmp(argv[i], "-Harker") == 0)
      F_Harker = 1;

    else if (str_icmp(argv[i], "-Xtal32ss") == 0)
      F_Xtal32ss = 1;

    else if (str_icmp(argv[i], "-ReadXYZ") == 0)
      F_ReadXYZ = 1;

    else if (str_ibegin(argv[i], "-CBMx=") == 0) {
      if (F_CBMx) usage(1);
      if (ParseStrXYZ(&argv[i][6], &CmdLnCBMx[0], CRBF, CTBF) != 0) usage(1);
      F_CBMx = 1;
      if (InverseRTMx(&CmdLnCBMx[0], &CmdLnCBMx[1], CRBF) <= 0) usage(1);
    }

    else if (str_icmp(argv[i], "-ShowPrimitive") == 0)
      F_ShowPrimitive = 1;

    else if (nSgList < 2)
    {
      SgName = argv[i];

      while (*SgName == ' ' || *SgName == '\t') SgName++;

      if (F_Convention == 'H' && isdigit(*SgName))
        SgList[nSgList].Convention = 'A';
      else
        SgList[nSgList].Convention = F_Convention;

      SgList[nSgList].SgName  = SgName;
      SgList[nSgList].InpTSgN = NULL;
      SgList[nSgList].RefTSgN = NULL;

      nSgList++;
    }
    else
      usage(0);
  }

  if (F_Debug && F_Verify) {
    TestAddlG();
    PrintClearSgError(0, 0);
    exit(0);
  }

  if (F_FindShTrMx)
  {
    (void) FindShTrMxWithMathematica(F_Debug);
    PrintClearSgError(1, 0);
    putc('\n', stdout);
  }

  if (F_ListTable)
  {
    if (F_Convention != 'P')
      ListTabSgName(F_ListTable, Last_F_Convention, stdout);
    else if (F_ListTable <= 17)
      ListTabPgName(F_ListTable, stdout);
    else
      usage(0);

    PrintClearSgError(1, 0);
    putc('\n', stdout);
  }

  if (F_CIF && F_Convention != 'P')
  {
    ListCIF(stdout);
    PrintClearSgError(1, 0);
    putc('\n', stdout);
  }

  if (F_CNSlib)
  {
    if (MakeCNSlib(stdout, F_Verbose) != 0) PrintClearSgError(1, 1);
    putc('\n', stdout);
  }

  if (nSgList == 0)
  {
    if (F_ReadXYZ)
    {
      SgList[nSgList].Convention = 'H';
      SgList[nSgList].SgName  = "P 1";
      SgList[nSgList].InpTSgN = NULL;
      SgList[nSgList].RefTSgN = NULL;
             nSgList++;
    }
    else if (   F_FindShTrMx == 0
             && F_ListTable == 0
             && F_CIF == 0
             && F_CNSlib == 0)
      usage(F_Verbose);
    else
      exit(0);
  }

  if (F_Space == 0)
  {
    putc('#', stdout);

    for (i = 0; i < argc; i++)
    {
      putc(' ', stdout);

      HaveSpace = 0;

      if (i) {
        for (n = 0; argv[i][n]; n++) {
          if (isspace(argv[i][n])) {
            HaveSpace = 1;
            break;
          }
        }
      }

      if (HaveSpace == 0)
        Fprintf(stdout, "%s", argv[i]);
      else
      {
        putc('"', stdout);

        for (n = 0; argv[i][n]; n++)
          if (argv[i][n] == '"') putc('+',        stdout);
          else                   putc(argv[i][n], stdout);

        putc('"', stdout);
      }
    }

    putc('\n', stdout);
    putc('\n', stdout);
  }

  if (AllocSgInfo(&SpgrInfo[0], F_Conditions) != 0) PrintClearSgError(0, 1);
  if (nSgList > 1 || F_Standard || F_CBMx)
    if (AllocSgInfo(&SpgrInfo[1], F_Conditions) != 0) PrintClearSgError(0, 1);

  BC_SgInfo.MaxList = 0;

  for (iSgList = 0; iSgList < nSgList; iSgList++)
  {
    if (iSgList) putc('\n', stdout);

    if (nSgList > 1 || F_Standard)
      Fprintf(stdout, "Setting %c:\n\n", "AB"[iSgList]);

    if (F_CBMx && iSgList == 0)
      SgInfo = &SpgrInfo[1];
    else
      SgInfo = &SpgrInfo[iSgList];

    F_Convention = SgList[iSgList].Convention;
    SgName       = SgList[iSgList].SgName;

    tsgn = NULL;

    if (F_Convention == 'A' || F_Convention == 'I' || F_Convention == 'P')
    {
      if (F_Convention == 'P')
        tsgn = FindTabPgNameEntry(SgName);
      else
        tsgn = FindTabSgNameEntry(0, SgName, F_Convention);

      if (tsgn == NULL)
      {
        PrintClearSgError(1, 0);

        if (F_Convention == 'P')
          progerror("Error: Unknown Plane Group Symbol");
        else
          progerror("Error: Unknown Space Group Symbol");
      }

      if (F_Space == 0 && F_CBMx == 0 && (F_ReadXYZ == 0 || iSgList > 0))
      {
        if (F_Convention == 'P') {
          Fprintf(stdout, "Plane Group  ");
          PrintTabSgNameEntry(tsgn, 0, 0, 0, stdout);
        }
        else {
          Fprintf(stdout, "Space Group  ");
          PrintTabSgNameEntry(tsgn, 0, 0, 0, stdout);
        }

        putc('\n', stdout);
      }

      SgName = tsgn->HallSymbol;

      if (F_Convention == 'P' || F_CBMx || (F_ReadXYZ && iSgList == 0))
        tsgn = NULL;
    }

    SgList[iSgList].InpTSgN = tsgn;

    ResetSgInfo(SgInfo);

    SgInfo->TabSgName = tsgn;
    if (tsgn) SgInfo->GenOption = 1;

    pos_hsym = ParseHallSymbol(SgName, SgInfo);

    if (SgError != NULL)
    {
      Fprintf(stdout, "    %s\n", SgName);
      for (i = 0; i < pos_hsym; i++) putc('-', stdout);
      Fprintf(stdout, "---^\n");
      Fprintf(stdout, "%s\n", SgError);
      exit(1);
    }

    if (F_ReadXYZ && iSgList == 0)
      ReadXYZ(stdin, SgInfo);

    if (CompleteSgInfo(SgInfo) != 0)
      PrintClearSgError(F_ClearError, 1);

    if (F_CBMx && iSgList == 0)
    {
      ShowCBMx(&CmdLnCBMx[0], &CmdLnCBMx[1], F_Maple + F_MathM);
      PrintClearSgError(F_ClearError, 0);
      putc('\n', stdout);

      ResetSgInfo(&SpgrInfo[0]);

      if (TransformSgInfo(SgInfo, &CmdLnCBMx[0], &CmdLnCBMx[1],
                          &SpgrInfo[0]) != 0)
        PrintClearSgError(F_ClearError, 1);

      SgInfo = &SpgrInfo[0];

      if (CompleteSgInfo(SgInfo) != 0)
        PrintClearSgError(F_ClearError, 1);
    }

    if (tsgn == NULL && F_Space == 0)
    {
      if (SgInfo->TabSgName)
      {
        Fprintf(stdout, "Space Group  ");
        PrintTabSgNameEntry(SgInfo->TabSgName, 0, 0, 0, stdout);
        putc('\n', stdout);
      }
      else
        Fprintf(stdout, "Hall Symbol  %s\n", SgInfo->HallSymbol);
    }

    PrintClearSgError(F_ClearError, 0);

    if (F_Xtal32ss) {
      if (xtal32ss(SgInfo) < 0)
        PrintClearSgError(F_ClearError, 1);
    }
    else {
      if (Set_ss(SgInfo) < 0)
        PrintClearSgError(F_ClearError, 1);
    }

    if (F_Space == 0) {
      ListSgInfo(SgInfo, F_XYZ, F_Verbose, stdout);
      PrintClearSgError(F_ClearError, 0);
    }

    if (F_ShowPrimitive) {
      (void) ShowPrimitiveSeitzMatrices(SgInfo);
      PrintClearSgError(F_ClearError, 0);
    }

    if (SgInfo->ReflCond)
    {
      if (SetReflCond(SgInfo) < 0)
        PrintClearSgError(F_ClearError, 1);
      else
      {
        ListReflCond(stdout, SgInfo->ReflCond, SgInfo->nReflCond,
                     F_Debug);

        if (SgInfo->RestCond)
        {
          if (SetRestCond(SgInfo) < 0)
            PrintClearSgError(F_ClearError, 1);
          else
            ListRestCond(stdout, SgInfo->RestCond, SgInfo->nRestCond,
                         F_Debug);
        }

        if (SetSysEnhanced(SgInfo) != 0)
          PrintClearSgError(F_ClearError, 1);
        else
          ListSysEnhanced(stdout, SgInfo->ReflCond, SgInfo->SysEnhanced,
                          SgInfo->nReflCond);
      }
    }

    if (F_Verify)
    {
          n = Check_ssVM(SgInfo);
      if (n == 0)
        SetSgError("Verify Error: Corrupt s.s. vectors and moduli");

      if (n <= 0)
        PrintClearSgError(F_ClearError, 1);
      else
        Fprintf(stdout, "Verify O.K.: s.s. vectors and moduli\n\n");
    }

    if (F_Characterize) {
      PutCharacterization(SgInfo, stdout);
      PrintClearSgError(F_ClearError, 0);
    }

    if (F_Harker) {
      PutHarkerInfo(SgInfo, Grid[3], stdout);
      PrintClearSgError(F_ClearError, 0);
    }

    if (F_AllXYZ) {
      PutAllXYZ(SgInfo, NULL, stdout);
      putc('\n', stdout);
      PrintClearSgError(F_ClearError, 0);
    }

    if (SgGrid(SgInfo, Grid, F_Debug) != 0)
      PrintClearSgError(F_ClearError, 1);
    else {
      Fprintf(stdout,
        "Grid Special Positions         %2d %2d %2d\n",
        Grid[0][0], Grid[0][1], Grid[0][2]);
      Fprintf(stdout,
        "Grid Symmetry Operations       %2d %2d %2d\n",
        Grid[1][0], Grid[1][1], Grid[1][2]);
      Fprintf(stdout,
        "Grid Screw/Glide Translations  %2d %2d %2d\n",
        Grid[2][0], Grid[2][1], Grid[2][2]);
      if (F_Harker) Fprintf(stdout,
        "Grid Harker Translations       %2d %2d %2d\n",
        Grid[3][0], Grid[3][1], Grid[3][2]);
      putc('\n', stdout);
    }

    if (F_Maple) {
      PutSymbAlgb(SgInfo, 0, stdout);
      PrintClearSgError(F_ClearError, 0);
      if (F_Verbose) PutRmI(SgInfo, 0, stdout);
      PrintClearSgError(F_ClearError, 0);
    }

    if (F_MathM) {
      PutSymbAlgb(SgInfo, 1, stdout);
      PrintClearSgError(F_ClearError, 0);
      if (F_Verbose) PutRmI(SgInfo, 1, stdout);
      PrintClearSgError(F_ClearError, 0);
    }

    if (F_Space) {
      PutSpaceSymFile(SgInfo, stdout);
      PrintClearSgError(F_ClearError, 0);
    }

    if (F_Shelx) {
      PutShelx(SgInfo, stdout);
      PrintClearSgError(F_ClearError, 0);
    }

    if (F_Schakal) {
      PutSchakal(SgInfo, stdout);
      PrintClearSgError(F_ClearError, 0);
    }

    if (F_CNS) {
      if (CNSrecord(stdout, SgInfo, F_Verbose) != 0)
        PrintClearSgError(F_ClearError, 1);
      putc('\n', stdout);
    }

    if (F_hklList) {
      Simple_hklList(SgInfo, (F_hklList > 0 ? 1 : 0),
                     abs(F_hklList), abs(F_hklList), abs(F_hklList),
                     F_Verbose);
      PrintClearSgError(F_ClearError, 0);
    }

    if (F_Verify) {
                                i = F_hklList; if (i == 0) i = 4;
      i = Verify_sghkl(SgInfo, (i > 0 ? 1 : 0), abs(i), abs(i), abs(i));
      if      (i < 0)
        PrintClearSgError(F_ClearError, 0);
      else if (i == 0)
        Fprintf(stdout, "Verify O.K.: sghkl\n\n");
    }

    if (nSgList > 1 || F_Standard)
    {
         CBMx = &SgList[iSgList].CBMx;
      InvCBMx = &SgList[iSgList].InvCBMx;

      SgList[iSgList].RefTSgN = FindReferenceSpaceGroup(SgInfo,
                                                        CBMx, InvCBMx);
      PrintClearSgError(F_ClearError, 0);

      if (F_Verify)
      {
        T_RTMx  TestInvCBMx[1];

        if (InverseRTMx(CBMx, TestInvCBMx, CRBF) == 0)
          progerror("Internal Error: Corrupt InverseRTMx()");

        for (i = 0; i < 9; i++)
          if (TestInvCBMx->s.R[i] != InvCBMx->s.R[i])
            progerror("Internal Error: Corrupt InverseRTMx()");

        for (i = 0; i < 3; i++)
          if (iModPositive(TestInvCBMx->s.T[i], CTBF) != InvCBMx->s.T[i])
            progerror("Internal Error: Corrupt InverseRTMx()");
      }

      if (SgList[iSgList].RefTSgN)
      {
        if (F_Verbose || F_Verify)
        {
          Fprintf(stdout, "Change of Basis => Reference Setting  ");
          PrintTabSgNameEntry(SgList[iSgList].RefTSgN, 0, 0, 0, stdout);
          putc('\n', stdout);

          ShowCBMx(CBMx, InvCBMx, F_Maple + F_MathM);
          PrintClearSgError(F_ClearError, 0);
        }

        if (F_Verify)
        {
          if (BC_SgInfo.MaxList == 0 && AllocSgInfo(&BC_SgInfo, 0) != 0)
            PrintClearSgError(0, 1);

          ResetSgInfo(&BC_SgInfo);

          if (TransformSgInfo(SgInfo, CBMx, InvCBMx, &BC_SgInfo) == 0)
            (void) CompleteSgInfo(&BC_SgInfo);

          if (SgError)
          {
            PrintClearSgError(F_ClearError, 0);
            SgList[iSgList].RefTSgN = NULL;
          }
          else if (BC_SgInfo.TabSgName != SgList[iSgList].RefTSgN)
          {
            Fprintf(stdout, "Hall Symbol  %s\n", BC_SgInfo.HallSymbol);
            SetSgError("Verify Error: Wrong CBMx/InvCBMx");
            PrintClearSgError(F_ClearError, 0);
            SgList[iSgList].RefTSgN = NULL;
          }
          else
            Fprintf(stdout, "Verify O.K.\n\n");
        }
      }

          tsgn = SgList[iSgList].RefTSgN;
      if (tsgn && F_Standard && nSgList == 1)
      {
        if (Last_F_Convention == 'A' || Last_F_Convention == 'I')
          SgList[nSgList].Convention = Last_F_Convention;
        else
          SgList[nSgList].Convention = 'A';

        SgList[nSgList].SgName  = SchoenfliesSymbols[tsgn->SgNumber];
        SgList[nSgList].InpTSgN = NULL;
        SgList[nSgList].RefTSgN = NULL;

        nSgList++;
      }
    }
  }

  if (   nSgList == 2
      && SgList[0].RefTSgN &&           SgList[1].RefTSgN
      && SgList[0].RefTSgN->SgNumber == SgList[1].RefTSgN->SgNumber)
  {
    putc('\n', stdout);
    Fprintf(stdout, "Change of Basis Setting A -> Setting B:\n");

    RTMxMultiply(   &CCBMx, &SgList[1].InvCBMx, &SgList[0].CBMx,
                 CRBF, CRBF * CTBF);
    RTMxMultiply(&CInvCBMx, &SgList[0].InvCBMx, &SgList[1].CBMx,
                 CRBF, CRBF * CTBF);

    for (i = 0; i < 12; i++)
    {
      if (   CCBMx.a[i] % CRBF) break;
      if (CInvCBMx.a[i] % CRBF) break;

         CCBMx.a[i] /= CRBF;
      CInvCBMx.a[i] /= CRBF;
    }

    if (i < 12)
    {
      SetSgError("Internal Error: Can't combine CBMx's");
      PrintClearSgError(1, 1);
    }
    else
    {
      ShowCBMx(&CCBMx, &CInvCBMx, F_Maple + F_MathM);
      PrintClearSgError(F_ClearError, 0);

      if (F_Verify)
      {
        ResetSgInfo(&BC_SgInfo);

        if (TransformSgInfo(&SpgrInfo[0], &CCBMx, &CInvCBMx, &BC_SgInfo) == 0)
          (void) CompleteSgInfo(&BC_SgInfo);

        if (SgError)
          PrintClearSgError(F_ClearError, 1);

        else if (strcmp(SpgrInfo[1].HallSymbol, BC_SgInfo.HallSymbol) != 0)
        {
          Fprintf(stdout, "Hall Symbol  %s\n", SpgrInfo[1].HallSymbol);
          Fprintf(stdout, "Hall Symbol  %s\n", BC_SgInfo.HallSymbol);
          SetSgError("Verify Error: Wrong CBMx/InvCBMx");
          PrintClearSgError(F_ClearError, 1);
        }
        else
          Fprintf(stdout, "Verify O.K.\n");
      }

      if (F_UnitCell)
      {
        putc('\n', stdout);

        if (TransformLatticeConstants(&LatConA, F_UnitCell,
                                      &LatConB, &SpgrInfo[0],
                                      CInvCBMx.s.R) != 0)
          PrintClearSgError(0, 1);

        Fprintf(stdout,
          "Setting A UnitCell  %.6g %.6g %.6g %.6g %.6g %.6g\n",
          LatConA.a, LatConA.b, LatConA.c,
          LatConA.alpha / PIover180,
          LatConA.beta  / PIover180,
          LatConA.gamma / PIover180);

        Fprintf(stdout,
          "Setting B UnitCell  %.6g %.6g %.6g %.6g %.6g %.6g\n",
          LatConB.a, LatConB.b, LatConB.c,
          LatConB.alpha / PIover180,
          LatConB.beta  / PIover180,
          LatConB.gamma / PIover180);
      }
    }

    putc('\n', stdout);
  }

  exit(0); /* old VAX didn't like "return 0;" */
  return 0;
}
