/*
  Space Group Info (c) 1994-97 Ralf W. Grosse-Kunstleve
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


#include "sginfo.h"


/* sginfo.c */

void PutAllXYZ(const T_SgInfo *SgInfo, const char *F_CNS, FILE *fpout);


#define Fputs   (void) fputs
#define Fprintf (void) fprintf


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


static int AllocSgInfo(T_SgInfo *SgInfo)
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

  SgInfo->ListRotMxInfo
    = malloc(SgInfo->MaxList * sizeof (*SgInfo->ListRotMxInfo));
  if (SgInfo->ListRotMxInfo == NULL) goto ReturnError;

  return 0;

  ReturnError:
  SetSgError("Not enough core");
  FreeSgInfo(SgInfo);
  return -1;
}


static int BuildSgInfo(T_SgInfo *SgInfo, const T_TabSgName *tsgn)
{
  ResetSgInfo(SgInfo);
  SgInfo->GenOption = 1;
  SgInfo->TabSgName = tsgn;

  (void) ParseHallSymbol(tsgn->HallSymbol, SgInfo);
  if (SgError != NULL) return -1;

  if (CompleteSgInfo(SgInfo) != 0) return -1;

  if (Set_ss(SgInfo) < 0) return -1;

  return 0;
}


static void ShowSgLabel(FILE *fpout, const char *Lbl)
{
  int  PrevDigit;


  PrevDigit = 0;

  for (; *Lbl && ! isspace(*Lbl); Lbl++)
  {
    if (isdigit(*Lbl))
    {
      if (PrevDigit) putc('(', fpout);
      putc(*Lbl, fpout);
      if (PrevDigit) putc(')', fpout);
      PrevDigit = 1;
    }
    else
    {
      PrevDigit = 0;
      if (*Lbl != '_') putc(*Lbl, fpout);
    }
  }
}


static void Show_ss(FILE *fpout, const char *LeadString,
                    const T_SgInfo *SgInfo)
{
  int           i_ssVM;
  const T_ssVM  *ssVM;


  ssVM = SgInfo->ssVM;

  for (i_ssVM = 0; i_ssVM < SgInfo->n_ssVM; i_ssVM++, ssVM++)
    Fprintf(fpout, "%seval($ssVM_%d = \"(%2d %2d %2d  %2d)\")\n",
      LeadString, i_ssVM + 1,
      ssVM->V[0], ssVM->V[1], ssVM->V[2], ssVM->M);
}


static int ShowGenAddlG(FILE *fpout, const char *LeadString,
                        const T_SgInfo *SgInfo)
{
  T_RTMx       AddlG[3];
  int         nAddlG;
  const char  *xyz;


           nAddlG = GetEuNoGen(SgInfo->TabSgName->SgNumber,
                               SgInfo->TabSgName->HallSymbol, 1, AddlG);
  if      (nAddlG < 0) {
    Fprintf(fpout, "%s{ eval($GenK2L = \"(unknown)\") }\n", LeadString);
  }
  else if (nAddlG == 1)
  {
        xyz = RTMx2XYZ(&AddlG[0], 1, STBF, 0, 0, 1, ",", NULL, 0);
    if (xyz == NULL) {
      SetSgError("Internal Error: ShowGenAddlG(): Corrupt AddlG");
      return -1;
    }

    Fprintf(fpout, "%seval($GenK2L = \"(%s)\")\n", LeadString, xyz);
  }
  else if (nAddlG != 0) {
    SetSgError("Internal Error: ShowGenAddlG(): Corrupt GetEuNoGen()");
    return -1;
  }

  return 0;
}


static int ShowSymGrid(FILE *fpout, const char *LeadString,
                       T_SgInfo *SgInfo, T_SgInfo *PgInfo)
{
  int           Grid[4][3], SymGrid[3], gcd, i;
  int           i_ssVM;
  const T_ssVM  *ssVM;

  const char *xyz = "xyz";


  if (SgGrid(SgInfo, Grid, 0) != 0) return -1;
  if (SetHarkerInfo(SgInfo, PgInfo, NULL, 0, Grid[3]) != 0) return -1;

  for (i = 0; i < 3; i++) {
    SymGrid[i] = iLCM(Grid[0][i], Grid[2][i]);
    SymGrid[i] = iLCM(SymGrid[i], Grid[3][i]);
  }

  ssVM = SgInfo->ssVM;

  for (i_ssVM = 0; i_ssVM < SgInfo->n_ssVM; i_ssVM++, ssVM++) {
    if (ssVM->M) {
      for (i = 0; i < 3; i++) {
        gcd = iGCD(ssVM->V[i], ssVM->M);
        SymGrid[i] = iLCM(SymGrid[i], ssVM->M / gcd);
      }
    }
  }

  for (i = 0; i < 3; i++)
    Fprintf(fpout, "%seval($SymGrid_%c = %d)\n",
      LeadString, xyz[i], SymGrid[i]);

  return 0;
}


static void ShowAsymUnit(FILE *fpout, const char *LeadString, const char **AU)
{
  int  first;


  if (*AU == NULL) return;

  for (first = 1; *AU; AU++)
  {
    if (first) {
      Fprintf(fpout, "%sasym=(%s", LeadString, *AU);
      first = 0;
    }
    else {
      Fprintf(fpout, "\n%s  and %s", LeadString, *AU);
    }
  }

  Fputs(")\n", fpout);
}


static int CoreCNSrecord(FILE *fpout, const char *LeadString, int F_All,
                         const T_TabSgName *ListTSgN[],
                         T_SgInfo *SgInfo, T_SgInfo *PgInfo)
{
  int         pgix;
  const char  **AU;


  if (SgInfo->TabSgName == NULL) {
    SetSgError(
      "Internal Error: CoreCNSrecord() works only for tabulated settings.");
    return -1;
  }

#ifndef OLDSTYLE
             Fprintf(fpout, "%seval($sg_number = %d)\n",
               LeadString, SgInfo->TabSgName->SgNumber);
#endif
             Fprintf(fpout, "%seval($nsym = %d)\n",
               LeadString, SgInfo->OrderL);
  if (F_All) Fprintf(fpout, "%seval($nsym_p = %d)\n",
               LeadString, SgInfo->OrderP);
  if (F_All) Fprintf(fpout, "%seval($nsym_p_acent = %d)\n",
               LeadString, SgInfo->nList);

  Fprintf(fpout, "%seval($sgname = \"", LeadString);
  ShowSgLabel(fpout, SgInfo->TabSgName->SgLabels);
  Fputs("\")\n", fpout);

  if (F_All) Fprintf(fpout, "%seval($crysys = \"%s\")\n",
               LeadString, XS_Name[SgInfo->XtalSystem]);

  pgix = PG_Index(SgInfo->PointGroup);

  Fprintf(fpout, "%seval($laue_class = \"%s\")\n",
    LeadString, PG_Names[PG_Index(LG_Code_of_PG_Index[pgix])]);

  if (F_All) Fprintf(fpout, "%seval($point_group = \"%s\")\n",
               LeadString, PG_Names[pgix]);

  if (SetPgInfo(SgInfo, PgInfo) != 0) return -1;
  if (PgInfo->TabSgName == NULL) {
    SetSgError("Internal Error: Cannot generate Patterson Space Group Number");
    return -1;
  }
  Fprintf(fpout, "%seval($patt_symm = \"", LeadString);
  ShowSgLabel(fpout, PgInfo->TabSgName->SgLabels);
  Fputs("\")", fpout);
  if (PgInfo->TabSgName != ListTSgN[PgInfo->TabSgName->SgNumber])
    Fprintf(fpout, " { Non-Standard Setting }");
  putc('\n', fpout);

  if (F_All) Fprintf(fpout, "%seval($latt_type = %c)\n",
               LeadString, SgInfo->TabSgName->SgLabels[0]);

  Show_ss(fpout, LeadString, SgInfo);
  if (ShowGenAddlG(fpout, LeadString, SgInfo) != 0) return -1;

  if (ShowSymGrid(fpout, LeadString, SgInfo, PgInfo) != 0) return -1;

  AU = NULL;

  if (SgInfo->TabSgName == ListTSgN[SgInfo->TabSgName->SgNumber])
    AU = StdAsymUnit(SgInfo->TabSgName->SgNumber, NULL);

  if (AU != NULL)
    ShowAsymUnit(fpout, LeadString, AU);
  else
    Fprintf(fpout, "%s{ asym=(unknown) }\n", LeadString);

  PutAllXYZ(SgInfo, LeadString, fpout);

  return 0;
}


static int MakeListTSgN(const T_TabSgName *ListTSgN[231])
{
  int   SgNumber;


  ListTSgN[0] = NULL;

  for (SgNumber = 1; SgNumber <= 230; SgNumber++)
  {
        ListTSgN[SgNumber] = FindTabSgNameEntry(SgNumber, NULL, 'A');
    if (ListTSgN[SgNumber] == NULL) {
      SetSgError("Internal Error: MakeListTSgN()");
      return -1;
    }
  }

  return 0;
}


int CNSrecord(FILE *fpout, T_SgInfo *SgInfo, int F_All)
{
  const T_TabSgName  *ListTSgN[231];
  T_SgInfo           PgInfo[1];


  if (AllocSgInfo(PgInfo) != 0) return -1;

  if (MakeListTSgN(ListTSgN) != 0) goto ReturnError;

  if (CoreCNSrecord(fpout, "", F_All, ListTSgN, SgInfo, PgInfo) != 0)
    goto ReturnError;

  FreeSgInfo(PgInfo);
  return 0;

  ReturnError:

  if (SgError == NULL) SetSgError("Internal Error: CNSrecord()");
  FreeSgInfo(PgInfo);

  return -1;
}


int MakeCNSlib(FILE *fpout, int F_All)
{
  int         SgNumber;
  const char  *NextIf;

  const T_TabSgName  *ListTSgN[231];
  T_SgInfo           SgInfo[1];
  T_SgInfo           PgInfo[1];


  if (AllocSgInfo(SgInfo) != 0) return -1;
  if (AllocSgInfo(PgInfo) != 0) { FreeSgInfo(SgInfo); return -1; }

  if (MakeListTSgN(ListTSgN) != 0) goto ReturnError;

  Fputs("! file  xtallib/spacegroup.lib\n", fpout);
  Fputs("! library of symmetry operators for CNS\n", fpout);
  Fprintf(fpout, "! Automatically generated with %s\n\n",
    SgInfoVersion());

  Fputs("set echo=off end\n", fpout);
  Fputs("set mess=off end\n\n", fpout);

  Fputs("! Reference for $GenK2L:\n", fpout);
  Fputs("!   International Tables Volume A (1987)\n", fpout);
  Fputs("!   Sections 15.2, 15.3 & 15.4.2 and column 6 of Table 15.3.2\n\n",
    fpout);

  Fputs("eval($sg_number = \"VOID\")\n", fpout);
  Fputs("eval($nsym = \"VOID\")\n", fpout);
  Fputs("eval($sgname = \"VOID\")\n", fpout);
  Fputs("eval($laue_class = \"VOID\")\n", fpout);
  Fputs("eval($patt_symm = \"VOID\")\n", fpout);
  Fputs("eval($ssVM_1 = \"VOID\")\n", fpout);
  Fputs("eval($ssVM_2 = \"VOID\")\n", fpout);
  Fputs("eval($ssVM_3 = \"VOID\")\n", fpout);
  Fputs("eval($GenK2L = \"VOID\")\n", fpout);
  Fputs("eval($SymGrid_x = \"VOID\")\n", fpout);
  Fputs("eval($SymGrid_y = \"VOID\")\n", fpout);
  Fputs("eval($SymGrid_z = \"VOID\")\n\n", fpout);

#ifdef OLDSTYLE
  NextIf = "if    ";

  for (SgNumber = 1; SgNumber <= 230; SgNumber++)
  {
    Fprintf(fpout, "%s ($sg = \"", NextIf);
    ShowSgLabel(fpout, ListTSgN[SgNumber]->SgLabels);
    Fputs("\") then\n", fpout);
    Fprintf(fpout, "  eval($sg_number = %3d)\n", SgNumber);

    NextIf = "elseif";
  }

  Fputs("else\n", fpout);
  Fputs("  eval($sg_number=$sg)\n", fpout);
  Fputs("end if\n\n", fpout);
#endif

  Fputs("eval($found = true)\n\n", fpout);

  NextIf = "if    ";

  for (SgNumber = 1; SgNumber <= 230; SgNumber++) {
#ifdef OLDSTYLE
    Fprintf(fpout, "%s ($sg_number = %d) then\n", NextIf, SgNumber);
#else
    Fprintf(fpout, "%s ($sg = \"", NextIf);
    ShowSgLabel(fpout, ListTSgN[SgNumber]->SgLabels);
    Fputs("\") then\n", fpout);
#endif
    if (BuildSgInfo(SgInfo, ListTSgN[SgNumber]) != 0)
      goto ReturnError;
    if (CoreCNSrecord(fpout, "  ", F_All, ListTSgN, SgInfo, PgInfo) != 0)
      goto ReturnError;
    NextIf = "elseif";
  }

  Fputs("else\n", fpout);
  Fputs("  eval($found = false)\n", fpout);
  Fputs("  display %SYMMETRY-ERR: space group $sg", fpout);
  Fputs(" does not exist in the library\n", fpout);
  Fputs("  abort\n", fpout);
  Fputs("end if\n\n", fpout);

  Fputs("if ( $found = true ) then\n", fpout);
  Fputs("  display SYMMETRY: found symmetry operators for", fpout);
  Fputs(" space group $sg in library\n", fpout);
  Fputs("end if\n\n", fpout);

  Fputs("set echo=on end\n", fpout);
  Fputs("set mess=on end\n", fpout);

  FreeSgInfo(SgInfo);
  FreeSgInfo(PgInfo);
  return 0;

  ReturnError:

  if (SgError == NULL) SetSgError("Internal Error: MakeCNSlib()");
  FreeSgInfo(SgInfo);
  FreeSgInfo(PgInfo);

  return -1;
}
