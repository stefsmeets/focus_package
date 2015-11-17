#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "unixstd.h"

#include "main.h"
#include "lib.h"
#include "matrix.h"
#include "xtal.h"
#include "trialoop.h"
#include "atominfo.h"
#include "nodsurf.h"
#include "io.h"




int BuildSpgrInfo(T_SgInfo *SgInfo, const char *SgName)
{
  int                VolLetter;
  const T_TabSgName  *tsgn;


  while (*SgName == ' ' || *SgName == '\t') SgName++;

  VolLetter = -1;

  if      (isdigit(*SgName))
    VolLetter = 'A';
  else if (str_ibegin(SgName, "VolA") == 0)
  {
    VolLetter = 'A';
    SgName += 4;
  }
  else if (   str_ibegin(SgName, "VolI") == 0
           || str_ibegin(SgName, "Vol1") == 0)
  {
    VolLetter = 'I';
    SgName += 4;
  }
  else if (str_ibegin(SgName, "Hall") == 0)
  {
    VolLetter = 0;
    SgName += 4;
  }

  while (*SgName == ' ' || *SgName == '\t') SgName++;

  if (*SgName == '#')
  {
    SgName++;
    if (VolLetter == -1) VolLetter = 'A';
  }

  if (VolLetter == -1)
    VolLetter = 'A';

  tsgn = NULL;

  if (VolLetter)
  {
    tsgn = FindTabSgNameEntry(0, SgName, VolLetter);
    if (tsgn == NULL) return -1;
    SgName = tsgn->HallSymbol;
  }

  SgInfo->MaxList = 192;

  CheckMalloc(SgInfo->ListSeitzMx, SgInfo->MaxList);
  CheckMalloc(SgInfo->ListRotMxInfo, SgInfo->MaxList);

  SgInfo->ReflCond = NULL;
  SgInfo->RestCond = NULL;
  SgInfo->SysEnhanced = NULL;

  ResetSgInfo(SgInfo);
  SgInfo->TabSgName = tsgn;

  (void) ParseHallSymbol(SgName, SgInfo);
  if (SgError != NULL) return -1;

  if (CompleteSgInfo(SgInfo) != 0) return -1;

  if (Set_ss(SgInfo) < 0) return -1;
  if (SgError != NULL) return -1;

  return 0;
}


static int ReadBreakIf(T_BreakIf *BreakIf,
                       char *arg, int pass, int OnlyPhaseDiff)
{
  int        ia, n, repeat;
  double     Value;
  char       extrac;
  T_ValPerW  PhaseDiff, DeltaR, *target;
  int        Branch;


  PhaseDiff.Value = 0.;
  PhaseDiff.Percent = 0;
  PhaseDiff.Weighted = 0;
  DeltaR.Value = 0.;
  DeltaR.Percent = 0;
  DeltaR.Weighted = 0;
  Branch = 0;
  target = NULL;

  ia = 1;
  n = 0;
  repeat = (OnlyPhaseDiff == 0);

  do
  {
    if (! str_arg(arg, Iline, ia++))
      return RCF_MissingParameter;

    if      (str_icmp(arg, "PhaseDiff") == 0)
    {
      if (Branch & BrkIf_Branch_PhaseDiff)
        return RCF_WrongParameter;

      Branch |= BrkIf_Branch_PhaseDiff;
      target = &PhaseDiff;
    }
    else if (OnlyPhaseDiff == 0 && str_icmp(arg, "DeltaR") == 0)
    {
      if (Branch & BrkIf_Branch_DeltaR)
        return RCF_WrongParameter;

      Branch |= BrkIf_Branch_DeltaR;
      target = &DeltaR;
    }
    else if (OnlyPhaseDiff == 0 && str_icmp(arg, "DeltaRw") == 0)
    {
      if (Branch & BrkIf_Branch_DeltaR)
        return RCF_WrongParameter;

      Branch |= BrkIf_Branch_DeltaR;
      target = &DeltaR;
      target->Weighted = 1;
    }
    else
      return RCF_WrongParameter;

    if (! str_arg(arg, Iline, ia++))
      return RCF_MissingParameter;

    if (strcmp(arg, "<") != 0)
      return RCF_WrongParameter;

    if (! str_arg(arg, Iline, ia++))
      return RCF_MissingParameter;

    n = sscanf(arg, "%lf%c", &Value, &extrac);

    if (n == 2)
    {
      if (extrac != '%')
        return RCF_WrongParameter;

      target->Percent = 1;
    }
    else if (n != 1)
      return RCF_WrongParameter;

    target->Value = Value;

    n = str_arg(arg, Iline, ia++);

    if (n && target->Percent == 0)
    {
      if (strcmp(arg, "%") == 0)
      {
        target->Percent = 1;

        n = str_arg(arg, Iline, ia++);
      }
    }

    if (n)
    {
      if (repeat)
      {
        if (str_icmp(arg, "or") == 0)
          Branch |= BrkIf_Branch_Or;
        else if (str_icmp(arg, "and") != 0)
          return RCF_WrongParameter;
      }
      else
        return RCF_WrongParameter;
    }
  }
  while (n && repeat--);

  if (pass == 1)
  {
    if (PhaseDiff.Percent) PhaseDiff.Value *= 0.01;
    if (DeltaR.Percent) DeltaR.Value *= 0.01;

    BreakIf->PhaseDiff.Value = PhaseDiff.Value;
    BreakIf->PhaseDiff.Percent = PhaseDiff.Percent;
    BreakIf->PhaseDiff.Weighted = PhaseDiff.Weighted;
    BreakIf->DeltaR.Value = DeltaR.Value;
    BreakIf->DeltaR.Percent = DeltaR.Percent;
    BreakIf->DeltaR.Weighted = DeltaR.Weighted;
    BreakIf->Branch = Branch;
  }

  return RCF_OK;
}


static int FindAtomTypeClass(char *arg)
{
  if      (strcmp(arg, "*") == 0)
    return ATC_General;
  else if (str_icmp(arg, "Node") == 0)
    return ATC_Node;
  else if (str_icmp(arg, "NodeBridge") == 0)
    return ATC_NodeBridge;

  return   ATC_Unknown;
}


static int ReadValueAndUncertainty(char *s, double *val, double *unc)
{
  int     n;
  char    ob, cb, extrac;
  double  v, u;


      n = sscanf(s, "%lf%c%lf%c%c", &v, &ob, &u, &cb, &extrac);
  if (n != 1 && n != 4)
    return -1;

  if (n == 1)
    u = 0.;
  else
  {
    if (ob != '(' || cb != ')')
      return -1;
  }

  if (val != NULL) *val = v;
  if (unc != NULL) *unc = u;

  if (unc != NULL)
    InternalError("unc not implemented");

  return 0;
}


static int TT_FWHM_2_UVW(Fprec *TT_FWHM, int nTT_FWHM,
                         Fprec *U, Fprec *V, Fprec *W)
{
  int    i, n;
  Fprec  a[9], b[3], *ai, *bi, tanT;


  if (nTT_FWHM % 2 || nTT_FWHM > 6)
    return -1;

  n = nTT_FWHM / 2;

  ai = a;
  bi = b;

  for (i = 0; i < n; i++)
  {
    if (*TT_FWHM < 1.e-3 || *TT_FWHM > 179.)
      return -1;

    tanT = AppTan(*TT_FWHM++ * .5 * PIover180);

      *ai++ = 1.;
    if (n > 1)
      *ai++ = tanT;
    if (n > 2)
      *ai++ = Square(tanT);

    *bi++ = Square(*TT_FWHM);
                    TT_FWHM++;
  }

  if (n && MxInverse(a, n, b, 1) != 0)
    return -1;

  *U = (n > 0 ? b[0] : 0.);
  *V = (n > 1 ? b[1] : 0.);
  *W = (n > 2 ? b[2] : 0.);

  return 0;
}


static int TT_Asym_2_ai(Fprec *TT_Asym, int nTT_Asym,
                        Fprec *a1, Fprec *a2, Fprec *a3)
{
  int    i, n;
  Fprec  a[9], b[3], *ai, *bi, tanT;


  if (nTT_Asym % 2 || nTT_Asym > 6)
    return -1;

  n = nTT_Asym / 2;

  ai = a;
  bi = b;

  for (i = 0; i < n; i++)
  {
        tanT = AppTan(*TT_Asym++ * .5 * PIover180);
    if (tanT == 0.)
      return -1;

      *ai++ = 1.;
    if (n > 1)
      *ai++ = 1. / tanT;
    if (n > 2)
      *ai++ = 1. / Square(tanT);

    *bi++ = *TT_Asym;
             TT_Asym++;
  }

  if (n && MxInverse(a, n, b, 1) != 0)
    return -1;

  *a1 = (n > 0 ? b[0] : 0.);
  *a2 = (n > 1 ? b[1] : 0.);
  *a3 = (n > 2 ? b[2] : 0.);

  return 0;
}


typedef struct { char *text; int class; } T_KeyWord;

#define MaxKeyWords 59

static const T_KeyWord KeyWords[MaxKeyWords + 1] =
  {
    { "Title",                    'I' },
    { "End",                      'I' },
    { "UnitCell",                 'I' },
    { "SpaceGroup",               'I' },

    { "AtomType",                 'I' },
    { "Chemistry",                'I' },
    { "MaxPotentialAtoms",        'I' },
    { "MaxRecycledAtoms",         'I' },

    { "FwSearchMethod",           'I' },
    { "MaxPeaksFwSearch",         'I' },
    { "MaxPeaksFwFragmentSearch", 'I' },
    { "MinNodeDistance",          'I' },
    { "MaxNodeDistance",          'I' },
    { "MinSymNodes",              'I' },
    { "MaxSymNodes",              'I' },
    { "NodeType",                 'I' },
    { "MinLoopSize",              'I' },
    { "MaxLoopSize",              'I' },
    { "EvenLoopSizesOnly",        'I' },
    { "Check3DimConnectivity",    'I' },
    { "IdealT_NodeDistance",      'I' },
    { "IdealTetrahedralEdge",     'I' },
    { "IdealTetrahedralVolume",   'I' },
    { "CheckTetrahedralGeometry", 'I' },

    { "RandomInitialization",     'I' },
    { "RandomBlindCalls",         'I' },
    { "FeedBackCycles",           'I' },
    { "FeedBackBreakIf",          'I' },

    { "Grid_xyz",                 'I' },
    { "PeakSearchLevel",          'I' },
    { "eDensityCutOff",           'I' },
    { "MinPfI",                   'I' },
    { "CatchDistance",            'I' },
    { "eD_PeaksSortElement",      'I' },

    { "Lambda",                   'I' },
    { "FobsMin_d",                'I' },
    { "FobsScale",                'I' },
    { "SigmaCutOff",              'I' },
    { "OverlapFactor",            'I' },
    { "OverlapAction",            'I' },
    { "ReflectionUsage",          'I' },
    { "GenerateFWHM",             'I' },
    { "AbsentRedivisionLimit",    'I' },

    { "tRius",                    'I' },
    { "MinConsecutiveNodes",      'I' },

    { "ProfileStartEndStep",      'I' },
    { "ProfilePOLRA",             'I' },
    { "ProfileFWHM",              'I' },
    { "ProfileAsym",              'I' },
    { "ProfilePeakShape",         'I' },
    { "PseudoVoigtPeakRange",     'I' },
    { "PseudoVoigtFracLorentz",   'I' },
    { "ProfileBackground",        'I' },
    { "ProfileReferenceRefl",     'I' },
    { "ProfileReferenceMax",      'I' },

    { "NormalizeSurface",         'I' },

    { "Site",                     'I' },
    { "ScatteringFactor",         'I' },
    
    { "ScatteringFactorTable",    'I' },

    { NULL,                      '\0' }
  };

enum
  {
    KW_Title = 0,
    KW_End,
    KW_UnitCell,
    KW_SpaceGroup,

    KW_AtomType,
    KW_Chemistry,
    KW_MaxPotentialAtoms,
    KW_MaxRecycledAtoms,

    KW_FwSearchMethod,
    KW_MaxPeaksFwSearch,
    KW_MaxPeaksFwFragmentSearch,
    KW_MinNodeDistance,
    KW_MaxNodeDistance,
    KW_MinSymNodes,
    KW_MaxSymNodes,
    KW_NodeType,
    KW_MinLoopSize,
    KW_MaxLoopSize,
    KW_EvenLoopSizesOnly,
    KW_Check3DimConnectivity,
    KW_IdealT_NodeDistance,
    KW_IdealTetrahedralEdge,
    KW_IdealTetrahedralVolume,
    KW_CheckTetrahedralGeometry,

    KW_RandomInitialization,
    KW_RandomBlindCalls,
    KW_FeedBackCycles,
    KW_FeedBackBreakIf,

    KW_Grid_xyz,
    KW_PeakSearchLevel,
    KW_eDensityCutOff,
    KW_MinPfI,
    KW_CatchDistance,
    KW_eD_PeaksSortElement,

    KW_Lambda,
    KW_FobsMin_d,
    KW_FobsScale,
    KW_SigmaCutOff,
    KW_OverlapFactor,
    KW_OverlapAction,
    KW_ReflectionUsage,
    KW_GenerateFWHM,
    KW_AbsentRedivisionLimit,

    KW_tRius,
    KW_MinConsecutiveNodes,

    KW_ProfileStartEndStep,
    KW_ProfilePOLRA,
    KW_ProfileFWHM,
    KW_ProfileAsym,
    KW_ProfilePeakShape,
    KW_PseudoVoigtPeakRange,
    KW_PseudoVoigtFracLorentz,
    KW_ProfileBackground,
    KW_ProfileReferenceRefl,
    KW_ProfileReferenceMax,

    KW_NormalizeSurface,

    KW_Site,
    KW_ScatteringFactor,
    
    KW_ScatteringFactorTable,

    KW_User
  };

static char *SepCom(char *s)
{
  char        *eol, c;

  eol = s;
  while ((c = *s) != '\0' && c != '#')
  {
    s++;
    if (c != ' ' && c != '\t') eol = s;
  }

  if (*s)
  {
    *eol = '\0';
    s++;
  }
  else
    s = NULL;

  return s;
}

int ReadCmdFile(FILE *fpcmd)
{
  int              pass;
  char             arg[LNLINE + 1];
  /* char          *comment; */
  const T_KeyWord  *kw;
  int              ikw, ContKW;

  char    extrac, *cp;
  int     iarg, i, n;
  double  x;

  int  iRandomBlindCalls, iFobsRaw, iAtomType, iInputChemistry, iSite;
  int  iNodeTypes, iFeedBackCycles;
  int  iSF_Tables, jSF_Tables;


  iSF_Tables =
  jSF_Tables =
  nSF_Tables = 0;

  Itabsize = 8;
  Imode = I_notabs | I_noleadblks | I_notrailblks | I_noemptyline;

  /* nKeyWords = KW_User; */

  for (pass = 0; pass < 2; pass++)
  {
    (void) fseek(fpcmd, 0L, SEEK_SET);

    iRandomBlindCalls =
    iFobsRaw =
    iAtomType =
    iInputChemistry =
    iNodeTypes =
    iFeedBackCycles =
    iSite =
    iSF_Tables = 0;

    ContKW = KW_End;

    Ilinec = 0;

    while (fgetIline(fpcmd))
    {
      if      (*Iline == '&')
      {
        if (ContKW == KW_End)
          return RCF_UnexpectedEOC;

        ContKW = KW_End;
      }
      else if (*Iline != '#' && *Iline != '$')
      {
        /* comment = */ (void) SepCom(Iline);

        (void) str_arg(arg, Iline, 0);

        if (ContKW != KW_End)
        {
          ikw = ContKW;
          iarg = 0;
        }
        else
        {
          for (ikw = 0, kw = KeyWords; kw->class != '\0'; ikw++, kw++)
          {
            if (kw->class == 'I')
            {
              if (str_icmp(arg, kw->text) == 0) break;
            }
            else
              InternalError(NULL);
          }

          if (ikw == KW_End)
            break;

          if (kw->class == '\0' && isalpha(*Iline))
          {
            ikw = KW_Site;
            iarg = 0;
          }
          else
            iarg = 1;
        }

        if (ikw < KW_User)
        {
          switch (ikw)
          {
            case KW_Title:
              if (pass == 1)
              {
                if (TopTitle != NULL) return RCF_WrongParameter;
                AppStrdup(Iline + 5, TopTitle);
                if (TopTitle == NULL) NotEnoughCore();
                (void) noblks(TopTitle);
              }
              break;

            case KW_UnitCell:
            {
              for (iarg = 1 ; iarg <= 6; iarg++)
              {
                if (! str_arg(arg, Iline, iarg))
                  break;

                if (sscanf(arg, "%lf%c", &x, &extrac) != 1 || x <= 0.)
                  return RCF_WrongParameter;

                if (pass == 1) switch (iarg)
                {
                  case 1: LatConD.a = x; break;
                  case 2: LatConD.b = x; break;
                  case 3: LatConD.c = x; break;
                  case 4: LatConD.alpha = x * PIover180; break;
                  case 5: LatConD.beta  = x * PIover180; break;
                  case 6: LatConD.gamma = x * PIover180; break;
                }
              }

              if (iarg == 1)
                return RCF_MissingParameter;

              if (pass == 1)
                nInputLatConD = iarg - 1;

              break;
            }

            case KW_SpaceGroup:
            {
              if (pass == 0)
              {
                        cp = Iline;
                while (*cp == ' ') cp++;
                while (*cp != ' ' && *cp != '\0') cp++;
                if (BuildSpgrInfo(&SpgrInfo, cp) != 0)
                  return RCF_WrongParameter;
              }
              break;
            }

            case KW_AtomType:
            {
              T_AtomType    AT;
              T_AtomRadius  *AR;

              if (! str_arg(arg, Iline, iarg++))
                return RCF_MissingParameter;

              if      (arg[0] == '-') AT.On = 0;
              else if (arg[0] == '+') AT.On = 1;
              else
                return RCF_WrongParameter;

              if (arg[1] == '\0') {
                if (! str_arg(arg, Iline, iarg++))
                  return RCF_MissingParameter;
                cp = &arg[0];
              }
              else
                cp = &arg[1];

                  AT.Class = FindAtomTypeClass(cp);
              if (AT.Class == ATC_Unknown)
                return RCF_WrongParameter;

              if (! str_arg(arg, Iline, iarg++))
                return RCF_MissingParameter;

              if (pass == 1)
              {
                AppStrdup(arg, AT.SF_Info.Lbl);
                if (AT.SF_Info.Lbl == NULL)
                  NotEnoughCore();
              }

              if (! str_arg(arg, Iline, iarg++))
                return RCF_MissingParameter;
              if (sscanf(arg, "%d%c", &i, &extrac) != 1 || i < 0)
                return RCF_WrongParameter;
              AT.nPerUnitCell = i;

              i = str_arg(arg, Iline, iarg++);
              if (i == 0 || strcmp(arg, "*") == 0)
                AT.OccDefault  = Df_OccDefault;
              else
              {
                if (sscanf(arg, "%lf%c", &x, &extrac) != 1 || x <= 0.)
                  return RCF_WrongParameter;
                AT.OccDefault  = x;
              }

              if (i != 0) i = str_arg(arg, Iline, iarg++);
              if (i == 0 || strcmp(arg, "*") == 0)
                AT.UisoDefault = Df_UisoDefault;
              else
              {
                if (sscanf(arg, "%lf%c", &x, &extrac) != 1 || x < 0.)
                  return RCF_WrongParameter;
                AT.UisoDefault = x;
              }

              if (i != 0) i = str_arg(arg, Iline, iarg++);
              if (i == 0 || strcmp(arg, "*") == 0)
                AT.eListAR = NULL;
              else
              {
                if (isalpha(arg[0]))
                {
                  AT.eListAR = FindAtomRadius(arg, 1);
                  if (AT.eListAR == NULL)
                    return RCF_WrongParameter;
                }
                else
                {
                  if (sscanf(arg, "%lf%c", &x, &extrac) != 1)
                    return RCF_WrongParameter;

                  if (pass == 1)
                  {
                    CheckMalloc(AR, 1);
                    AR->Label  = NULL;
                    AR->Radius = x;
                    AT.eListAR = AR;
                  }
                }
              }

              if (pass == 1)
              {
                if (iAtomType >= nAtomType)
                  return RCF_ChangingInput;

                if (AT.eListAR == NULL)
                {
                      AT.eListAR = FindAtomRadius(AT.SF_Info.Lbl, 0);
                  if (AT.eListAR == NULL)
                    return RCF_WrongParameter;
                }

                AtomType[iAtomType].On               = AT.On;
                AtomType[iAtomType].Class            = AT.Class;
                AtomType[iAtomType].SF_Info.Lbl      = AT.SF_Info.Lbl;
                AtomType[iAtomType].SF_Info.SFT      = NULL;
                AtomType[iAtomType].SF_Info.CAA      = NULL;
                AtomType[iAtomType].SF_Info.f_stol_0 = -1.;
                AtomType[iAtomType].nPerUnitCell     = AT.nPerUnitCell;
                AtomType[iAtomType].OccDefault       = AT.OccDefault;
                AtomType[iAtomType].UisoDefault      = AT.UisoDefault;
                AtomType[iAtomType].eListAR          = AT.eListAR;
              }

              iAtomType++;

              break;
            }

            case KW_Chemistry:
            {
              T_InputChemistry  IC;
              T_AtomID          *AtomID;

              if (! str_arg(arg, Iline, iarg++))
                return RCF_MissingParameter;

              if      (str_icmp(arg, "MinDistance") == 0) {
                IC.Type = ICT_MinDistance; n = 2; }
              else if (str_icmp(arg, "MaxDistance") == 0) {
                IC.Type = ICT_MaxDistance; n = 2; }
              else if (str_icmp(arg, "AddBondAtom") == 0) {
                IC.Type = ICT_AddBondAtom; n = 3; }
              else
                return RCF_WrongParameter;

              for (i = 0; i < n; i++)
              {
                if (! str_arg(arg, Iline, iarg++))
                  return RCF_MissingParameter;

                    IC.AtomID[i].Class = FindAtomTypeClass(arg);
                if (IC.AtomID[i].Class == ATC_Unknown)
                  return RCF_WrongParameter;

                if (! str_arg(arg, Iline, iarg++))
                  return RCF_MissingParameter;

                if (pass == 1)
                {
                  AppStrdup(arg, IC.AtomID[i].SF_Info.Lbl);
                  if (NULL ==    IC.AtomID[i].SF_Info.Lbl) NotEnoughCore();
                }
              }

              if (n == 2)
              {
                if (! str_arg(arg, Iline, iarg++))
                  return RCF_MissingParameter;

                if (sscanf(arg, "%lf%c", &x, &extrac) != 1)
                  return RCF_WrongParameter;

                IC.Value = x;
              }

              if (pass == 1)
              {
                InputChemistry[iInputChemistry].Type  = IC.Type;
                InputChemistry[iInputChemistry].Value = IC.Value;

                AtomID = InputChemistry[iInputChemistry].AtomID;

                for (i = 0; i < n; i++)
                {
                  AtomID[i].Class = IC.AtomID[i].Class;

                  AtomID[i].SF_Info.Lbl      = IC.AtomID[i].SF_Info.Lbl;
                  AtomID[i].SF_Info.SFT      = NULL;
                  AtomID[i].SF_Info.CAA      = NULL;
                  AtomID[i].SF_Info.f_stol_0 = -1.;
                }
              }

              iInputChemistry++;

              break;
            }

            case KW_MaxPotentialAtoms:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%d%c", &i, &extrac) != 1 || i < 0)
                return RCF_WrongParameter;
              if (pass == 1) MaxPotentialAtoms = i;
              break;
            }

            case KW_MaxRecycledAtoms:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%d%c", &i, &extrac) != 1 || i < 0)
                return RCF_WrongParameter;
              if (pass == 1) MaxRecycledAtoms = i;
              break;
            }

            case KW_FwSearchMethod:
            {
              int  FwSM;

              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if      (str_icmp(arg, "FwTracking") == 0)
                           FwSM = FwSM_FwTracking;
              else if (str_icmp(arg, "AltFwTracking") == 0)
                           FwSM = FwSM_AltFwTracking;
              else
                return RCF_WrongParameter;

              if (pass == 1)
                FwSearchMethod = FwSM;

              break;
            }

            case KW_MaxPeaksFwSearch:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%d%c", &i, &extrac) != 1 || i < 0)
                return RCF_WrongParameter;
              if (pass == 1) MaxPeaksFwSearch = i;
              break;
            }

            case KW_MaxPeaksFwFragmentSearch:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%d%c", &i, &extrac) != 1 || i < 0)
                return RCF_WrongParameter;
              if (pass == 1) MaxPeaksFwFragmentSearch = i;
              break;
            }

            case KW_MinNodeDistance:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%lf%c", &x, &extrac) != 1)
                return RCF_WrongParameter;
              if (pass == 1) MinNodeDist2 = x * x;
              break;
            }

            case KW_MaxNodeDistance:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%lf%c", &x, &extrac) != 1)
                return RCF_WrongParameter;
              if (pass == 1) MaxNodeDist2 = x * x;
              break;
            }

            case KW_MinSymNodes:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%d%c", &i, &extrac) != 1 || i < 0)
                return RCF_WrongParameter;
              if (pass == 1) MinSymNodes = i;
              break;
            }

            case KW_MaxSymNodes:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (strcmp(arg, "*") == 0)
                i = -1;
              else if (sscanf(arg, "%d%c", &i, &extrac) != 1 || i < 0)
                return RCF_WrongParameter;
              if (pass == 1) MaxSymNodes = i;
              break;
            }

            case KW_NodeType:
            {
              T_NodeType NT;

              if (! str_arg(arg, Iline, iarg++))
                return RCF_MissingParameter;
              if (sscanf(arg, "%d%c", &i, &extrac) != 1 || i < 1)
                return RCF_WrongParameter;

                                        NT.CN           = i;
                                        NT.MaxNpAU      = 0;
              for (i = -6; i <= 6; i++) NT.NoSPO[i + 6] = 0;

              if (str_arg(arg, Iline, iarg++))
              {
                if (strcmp(arg, "*") != 0)
                {
                  if (sscanf(arg, "%d%c", &i, &extrac) != 1 || i < 1)
                    return RCF_WrongParameter;

                  NT.MaxNpAU = i;
                }

                while (str_arg(arg, Iline, iarg++))
                {
                  if (sscanf(arg, "%d%c", &i, &extrac) != 1)
                    return RCF_WrongParameter;

                  if (   i < -6 || i > 6
                      || i == 5 || i == -5
                      || i == 0 || i == 1)
                    return RCF_WrongParameter;

                  NT.NoSPO[i + 6] = 1;
                }
              }

              if (pass == 1)
              {
                for (i = 0; i < iNodeTypes; i++)
                  if (NT.CN == NodeTypes[i].CN)
                    return RCF_WrongParameter;

                (void) memcpy(&NodeTypes[iNodeTypes], &NT, sizeof NT);
              }

              iNodeTypes++;

              break;
            }

            case KW_MinLoopSize:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%d%c", &i, &extrac) != 1 || i < 0
                  || i == 1 || i == 2)
                return RCF_WrongParameter;
              if (pass == 1) MinLoopSize = i;
              break;
            }

            case KW_MaxLoopSize:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%d%c", &i, &extrac) != 1 || i < 3)
                return RCF_WrongParameter;
              if (pass == 1) MaxLoopSize = i;
              break;
            }

            case KW_EvenLoopSizesOnly:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if      (str_icmp(arg, "Off") == 0)
                i = 0;
              else if (str_icmp(arg, "On") == 0)
                i = 1;
              else
                return RCF_WrongParameter;
              if (pass == 1) EvenLoopSizesOnly = i;
              break;
            }

            case KW_Check3DimConnectivity:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if      (str_icmp(arg, "Off") == 0)
                i = 0;
              else if (str_icmp(arg, "On") == 0)
                i = 1;
              else
                return RCF_WrongParameter;
              if (pass == 1) Check3DimConnectivity = i;
              break;
            }

            case KW_IdealT_NodeDistance:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%lf%c", &x, &extrac) != 1)
                return RCF_WrongParameter;
              if (pass == 1) IdealT_NodeDist2 = x * x;
              break;
            }
            case KW_IdealTetrahedralEdge:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%lf%c", &x, &extrac) != 1)
                return RCF_WrongParameter;
              if (pass == 1) IdealTetrahedralEdge2 = x * x;
              break;
            }

            case KW_IdealTetrahedralVolume:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%lf%c", &x, &extrac) != 1 || x <= 0.)
                return RCF_WrongParameter;
              if (pass == 1)
                IdealTetrahedralVol2 = x * x;

              if (str_arg(arg, Iline, 2))
              {
                if (sscanf(arg, "%lf%c", &x, &extrac) != 1 || x <= 0.)
                  return RCF_WrongParameter;
                if (pass == 1)
                  IdealTetrahedralVolFrac = x;
              }
              break;
            }

            case KW_CheckTetrahedralGeometry:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if      (str_icmp(arg, "Off") == 0)
                i = 0;
              else if (str_icmp(arg, "Normal") == 0)
                i = 1;
              else if (str_icmp(arg, "Hard") == 0)
                i = 2;
              else
                return RCF_WrongParameter;
              if (pass == 1) ModeCheckTetrahedra = i;
              break;
            }

            case KW_RandomInitialization:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;

              if (str_icmp(arg, "Time") == 0)
              {
                if (pass == 1)
                {
                  RandomInitialization = 0;
                  FlagTimeRandomInitialization = 1;
                }
              }
              else
              {
                if (sscanf(arg, "%d%c", &i, &extrac) != 1)
                  return RCF_WrongParameter;

                if (pass == 1)
                {
                  RandomInitialization = i;
                  FlagTimeRandomInitialization = 0;
                }
              }
              break;
            }

            case KW_RandomBlindCalls:
            {
              while (str_arg(arg, Iline, iarg))
              {
                iarg++;

                if (sscanf(arg, "%d%c", &i, &extrac) != 1 || i < 0)
                  return RCF_WrongParameter;

                if (pass == 1)
                  RandomBlindCalls[iRandomBlindCalls] = i;

                iRandomBlindCalls++;
              }

              break;
            }

            case KW_FeedBackCycles:
            {
              while (str_arg(arg, Iline, iarg))
              {
                iarg++;

                if (sscanf(arg, "%d%c", &i, &extrac) != 1 || i < 0)
                  return RCF_WrongParameter;

                if (pass == 1)
                  FeedBackCycles[iFeedBackCycles] = i;

                iFeedBackCycles++;
              }

              break;
            }

            case KW_FeedBackBreakIf:
            {
                  n = ReadBreakIf(&FeedBackBreakIf, arg, pass, 0);
              if (n != RCF_OK)
                return n;
              break;
            }

            case KW_Grid_xyz:
            {
              for (iarg = 1; iarg <= 3; iarg++)
              {
                if (! str_arg(arg, Iline, iarg))
                  return RCF_MissingParameter;
                if (sscanf(arg, "%d%c", &n, &extrac) != 1)
                  return RCF_WrongParameter;

                if (pass == 1) switch (iarg)
                {
                  case 1: Nx = n; break;
                  case 2: Ny = n; break;
                  case 3: Nz = n; break;
                }
              }
              break;
            }

            case KW_PeakSearchLevel:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%d%c", &i, &extrac) != 1 || i < 0 || i > 3)
                return RCF_WrongParameter;
              if (pass == 1) PeakSearchLevel = i;
              break;
            }

            case KW_eDensityCutOff:
            {
              char  onec;

              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;

                  n = sscanf(arg, "%lf%c%c", &x, &onec, &extrac);
              if (n < 1 || n > 2 || (n == 2 && onec != '%'))
                return RCF_WrongParameter;
              if (n == 1)
              {
                if (str_arg(arg, Iline, 2))
                {
                  if (sscanf(arg, "%c%c", &onec, &extrac) != 1 || onec != '%')
                    return RCF_WrongParameter;
                  n = 2;
                }
              }
              if (pass == 1)
              {
                eDensityCutOff.Value = x;
                eDensityCutOff.Percent = (n == 2);
              }
              break;
            }

            case KW_MinPfI:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%d%c", &i, &extrac) != 1 || i < 10 || i > 19)
                return RCF_WrongParameter;
              if (pass == 1) MinPfI = i;
              break;
            }

            case KW_CatchDistance:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%lf%c", &x, &extrac) != 1)
                return RCF_WrongParameter;
              if (pass == 1)
                CatchDistance2 = x * x;
              break;
            }

            case KW_eD_PeaksSortElement:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if      (str_icmp(arg, "Grid_eD") == 0)
                i = PSE_Grid_eD;
              else if (str_icmp(arg, "Maximum") == 0)
                i = PSE_Maximum;
              else if (str_icmp(arg, "Integral") == 0)
                i = PSE_Integral;
              else
                return RCF_WrongParameter;
              if (pass == 1)
                eD_PeaksSortElement = i;
              break;
            }

            case KW_Lambda:
            {
              const T_ChXrayWaveLength  *cxw;

              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;

                  cxw = ChXrayWaveLengthOf(arg);
              if (cxw != NULL)
              {
                if (pass == 1)
                {
                  LambdaName   = cxw->Label;
                  LambdaLength = cxw->Length;
                }
              }
              else
              {
                if (sscanf(arg, "%lf%c", &x, &extrac) != 1 || x <= 0.)
                  return RCF_WrongParameter;
                if (pass == 1)
                {
                  LambdaName   = "";
                  LambdaLength = x;
                }
              }
              break;
            }

            case KW_FobsMin_d:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%lf%c", &x, &extrac) != 1)
                return RCF_WrongParameter;
              if (pass == 1)
              {
                if (x > 0.) FobsMaxQ = 1. / (x * x);
                else        FobsMaxQ = 0.;
              }
              break;
            }

            case KW_FobsScale:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%lf%c", &x, &extrac) != 1 || x <= 0.)
                return RCF_WrongParameter;
              if (pass == 1) FobsScale = x;
              break;
            }

            case KW_SigmaCutOff:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%lf%c", &x, &extrac) != 1)
                return RCF_WrongParameter;
              if (pass == 1) SigmaCutOff = x;
              break;
            }

            case KW_OverlapFactor:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%lf%c", &x, &extrac) != 1 || x < 0.)
                return RCF_WrongParameter;
              if (pass == 1) OverlapFactor = x;
              break;
            }

            case KW_OverlapAction:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if      (str_icmp(arg, "NoAction") == 0)
                          i = OvlAct_NoAction;
              else if (str_icmp(arg, "EqualF2") == 0)
                          i = OvlAct_EqualF2;
              else if (str_icmp(arg, "EqualMF2") == 0)
                          i = OvlAct_EqualMF2;
              else if (str_icmp(arg, "Omit") == 0)
                          i = OvlAct_Omit;
              else
                return RCF_WrongParameter;
              if (pass == 1) OverlapAction = i;
              break;
            }

            case KW_ReflectionUsage:
            {
              char  onec;

              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;

                  n = sscanf(arg, "%lf%c%c", &x, &onec, &extrac);
              if (n < 1 || n > 2 || (n == 2 && onec != '%') || x < 0.)
                return RCF_WrongParameter;
              if (n == 1)
              {
                if (str_arg(arg, Iline, 2))
                {
                  if (sscanf(arg, "%c%c", &onec, &extrac) != 1 || onec != '%')
                    return RCF_WrongParameter;
                  n = 2;
                }
              }
              if (n == 1)
                x = (int)(x + .5);
              else if (x > 100.)
                return RCF_WrongParameter;
              if (pass == 1)
              {
                ReflectionUsage.Value = x;
                ReflectionUsage.Percent = (n == 2);
              }
              break;
            }

            case KW_GenerateFWHM:
            {
              int    nval, mode;
              Fprec   val[6];

              nval = 0;
              mode = 0;

              if (str_arg(arg, Iline, iarg++))
              {
                if (str_icmp(arg, "UVW") == 0)
                {
                  mode = 1;

                  if (! str_arg(arg, Iline, iarg++))
                    mode = 0;
                }
                else
                  mode = 2;

                if (mode)
                {
                  while (nval < 6)
                  {
                    if (sscanf(arg, "%lf%c", &x, &extrac) != 1)
                      return RCF_WrongParameter;

                    val[nval++] = x;

                    if (mode == 1 && nval == 3)
                      break;

                    if (! str_arg(arg, Iline, iarg++))
                      break;
                  }

                  if (mode == 2 && nval % 2 != 0)
                    return RCF_MissingParameter;
                }
              }

              if (pass == 1)
              {
                if (mode == 2)
                {
                  GenerateFWHM.Valid = nval / 2;

                  if (TT_FWHM_2_UVW(val, nval, &GenerateFWHM.U,
                                               &GenerateFWHM.V,
                                               &GenerateFWHM.W) != 0)
                    return RCF_WrongParameter;
                }
                else
                {
                  GenerateFWHM.Valid = nval;
                  GenerateFWHM.U = (nval > 0 ? val[0] : 0.);
                  GenerateFWHM.V = (nval > 1 ? val[1] : 0.);
                  GenerateFWHM.W = (nval > 2 ? val[2] : 0.);
                }
              }

              break;
            }

            case KW_AbsentRedivisionLimit:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;

                  n = sscanf(arg, "%lf%c", &x, &extrac);
              if (n != 1 || x < 0.)
                return RCF_WrongParameter;

              if (! str_arg(arg, Iline, 2))
                 n = ARLT_FWHM;
              else
              {
                if (str_icmp(arg, "FWHM") == 0)
                  n = ARLT_FWHM;
                else if (   str_icmp(arg, "2Theta") == 0
                         || str_icmp(arg, "2-Theta") == 0
                         || str_icmp(arg, "Degree2Theta") == 0)
                  n = ARLT_Degree2Theta;
                else
                  return RCF_WrongParameter;
              }

              if (pass == 1)
              {
                AbsentRedivisionLimit.Value = x;
                AbsentRedivisionLimit.Type  = n;
                AbsentRedivisionLimit.Frac_I_ignored = 0.;
                AbsentRedivisionLimit.Frac_I_w_moved = 0.;

              }
              break;
            }

            case KW_tRius:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%lf%c", &x, &extrac) != 1)
                return RCF_WrongParameter;
              if (pass == 1)
                tRius = x;
              break;
            }
            case KW_MinConsecutiveNodes:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (str_icmp(arg, "NoDisplay") == 0)
                i = -1;
              else if (sscanf(arg, "%d%c", &i, &extrac) != 1 || i < 0)
                return RCF_WrongParameter;
              if (pass == 1) MinConsecutiveNodes = i;
              break;
            }

            case KW_ProfileStartEndStep:
            {
              int    nval;
              Fprec   val[3];

              val[0] = ProfileStart;
              val[1] = ProfileEnd;
              val[2] = ProfileStep;

              nval = 0;

              while (str_arg(arg, Iline, iarg++))
              {
                if (nval == 3)
                  return RCF_WrongParameter;

                if (str_icmp(arg, "*") != 0)
                {
                  if (sscanf(arg, "%lf%c", &x, &extrac) != 1 || x < 0.)
                    return RCF_WrongParameter;

                  if (nval < 2)
                  {
                    if (x > 179.)
                    return RCF_WrongParameter;
                  }
                  else if (x == 0.)
                    return RCF_WrongParameter;

                  val[nval] = x;
                }

                nval++;
              }

              if (pass == 1)
              {
                ProfileStart = val[0];
                ProfileEnd   = val[1];
                ProfileStep  = val[2];

                if (ProfileEnd != 0. && ProfileStart > ProfileEnd)
                  ProfileStart = ProfileEnd;
              }

              break;
            }

            case KW_ProfilePOLRA:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%lf%c", &x, &extrac) != 1 || x < 0. || x > 1.)
                return RCF_WrongParameter;
              if (pass == 1) ProfilePOLRA = x;
              break;
            }

            case KW_ProfileFWHM:
            {
              int    nval, mode;
              Fprec   val[6];

              nval = 0;
              mode = 0;

              if (str_arg(arg, Iline, iarg++))
              {
                if (str_icmp(arg, "UVW") == 0)
                {
                  mode = 1;

                  if (! str_arg(arg, Iline, iarg++))
                    mode = 0;
                }
                else
                  mode = 2;

                if (mode)
                {
                  while (nval < 6)
                  {
                    if (sscanf(arg, "%lf%c", &x, &extrac) != 1)
                      return RCF_WrongParameter;

                    val[nval++] = x;

                    if (mode == 1 && nval == 3)
                      break;

                    if (! str_arg(arg, Iline, iarg++))
                      break;
                  }

                  if (mode == 2 && nval % 2 != 0)
                    return RCF_MissingParameter;
                }
              }

              if (pass == 1)
              {
                if (mode == 2)
                {
                  ProfileFWHM.Valid = nval / 2;

                  if (TT_FWHM_2_UVW(val, nval, &ProfileFWHM.U,
                                               &ProfileFWHM.V,
                                               &ProfileFWHM.W) != 0)
                    return RCF_WrongParameter;
                }
                else
                {
                  ProfileFWHM.Valid = nval;
                  ProfileFWHM.U = (nval > 0 ? val[0] : 0.);
                  ProfileFWHM.V = (nval > 1 ? val[1] : 0.);
                  ProfileFWHM.W = (nval > 2 ? val[2] : 0.);
                }
              }

              break;
            }

            case KW_ProfileAsym:
            {
              int    nval, mode;
              Fprec   val[6];

              nval = 0;
              mode = 0;

              if (str_arg(arg, Iline, iarg++))
              {
                if (str_icmp(arg, "a(i)") == 0)
                {
                  mode = 1;

                  if (! str_arg(arg, Iline, iarg++))
                    mode = 0;
                }
                else
                  mode = 2;

                if (mode)
                {
                  while (nval < 6)
                  {
                    if (sscanf(arg, "%lf%c", &x, &extrac) != 1)
                      return RCF_WrongParameter;

                    val[nval++] = x;

                    if (mode == 1 && nval == 3)
                      break;

                    if (! str_arg(arg, Iline, iarg++))
                      break;
                  }

                  if (mode == 2 && nval % 2 != 0)
                    return RCF_MissingParameter;
                }
              }

              if (pass == 1)
              {
                if (mode == 2)
                {
                  ProfileAsym.Valid = nval / 2;

                  if (TT_Asym_2_ai(val, nval, &ProfileAsym.a1,
                                              &ProfileAsym.a2,
                                              &ProfileAsym.a3) != 0)
                    return RCF_WrongParameter;
                }
                else
                {
                  ProfileAsym.Valid = nval;
                  ProfileAsym.a1 = (nval > 0 ? val[0] : 0.);
                  ProfileAsym.a2 = (nval > 1 ? val[1] : 0.);
                  ProfileAsym.a3 = (nval > 2 ? val[2] : 0.);
                }
              }

              break;
            }

            case KW_ProfilePeakShape:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if      (str_icmp(arg, "PseudoVoigt") == 0)
                ProfilePeakShape = PPS_PseudoVoigt;
              else if (str_icmp(arg, "StdPeak") == 0)
                ProfilePeakShape = PPS_StdPeak;
              else
                return RCF_WrongParameter;
              break;
            }

            case KW_PseudoVoigtPeakRange:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%lf%c", &x, &extrac) != 1 || x < 1.)
                return RCF_WrongParameter;
              if (pass == 1) PseudoVoigtPeakRange = x;
              break;
            }

            case KW_PseudoVoigtFracLorentz:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%lf%c", &x, &extrac) != 1 || x < 0. || x > 1.)
                return RCF_WrongParameter;
              if (pass == 1) PseudoVoigtFracLorentz = x;
              break;
            }

            case KW_ProfileBackground:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%lf%c", &x, &extrac) != 1 || x < 0.)
                return RCF_WrongParameter;
              if (pass == 1) ProfileBackground = x;
              break;
            }

            case KW_ProfileReferenceRefl:
            {
              int    nval;
              Fprec   val[4];

              for (nval = 0; nval < 4; nval++)
              {
                if (! str_arg(arg, Iline, iarg++))
                  break;

                if (sscanf(arg, "%lf%c", &x, &extrac) != 1)
                  return RCF_WrongParameter;

                val[nval] = x;
              }

              if      (nval == 1)
              {
                if (val[0] < 0. || val[0] > 179.)
                  return RCF_WrongParameter;

                if (pass == 1) {
                  ProfileReferenceRefl.Mode = PRRM_TTheta;
                  ProfileReferenceRefl.TTheta = val[0];
                }
              }
              else if (nval == 3)
              {
                if (pass == 1) {
                  ProfileReferenceRefl.Mode = PRRM_hkl;
                  ProfileReferenceRefl.h = INT_ROUNDED(val[0]);
                  ProfileReferenceRefl.k = INT_ROUNDED(val[1]);
                  ProfileReferenceRefl.l = INT_ROUNDED(val[2]);
                }
              }
              else if (nval)
                return RCF_WrongParameter;

              break;
            }

            case KW_ProfileReferenceMax:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if (sscanf(arg, "%lf%c", &x, &extrac) != 1 || x < 0.)
                return RCF_WrongParameter;
              if (pass == 1) ProfileReferenceMax = x;
              break;
            }

            case KW_NormalizeSurface:
            {
              int     iv;
              double  v[3];

              for (iv = 0; iv < 3; iv++) v[iv] = 0.;

              for (iv = 0; iv < 3; iv++)
              {
                if (! str_arg(arg, Iline, iv + 1))
                  break;

                if (sscanf(arg, "%lf%c", &v[iv], &extrac) != 1)
                  return RCF_WrongParameter;
              }

              if (   (iv == 1 && v[0] == 0.)
                  || (iv >= 2 && v[1] - v[0] == 0.)
                  || (iv == 3 && v[2] <= 0.))
                return RCF_WrongParameter;

              if (pass == 1) {
                nNormalizeSurface = iv;
                for (iv = 0; iv < 3; iv++)
                  vNormalizeSurface[iv] = v[iv];
              }

              break;
            }

            case KW_Site: /* input Label [ScatFact] x y z [Occ [Uiso]] */
            {
              T_Site  S;

              if (! str_arg(arg, Iline, iarg++))
                return RCF_MissingParameter;

              if (pass == 1)
              {
                AppStrdup(arg, S.Label);
                if (NULL    == S.Label) NotEnoughCore();
              }

              if (! str_arg(arg, Iline, iarg++))
                return RCF_MissingParameter;

              if (isalpha(*arg))
              {
                if (pass == 1)
                {
                  AppStrdup(arg, S.SF_Info.Lbl);
                  if (NULL ==    S.SF_Info.Lbl) NotEnoughCore();
                }

                if (! str_arg(arg, Iline, iarg++))
                  return RCF_MissingParameter;
              }
              else if (strcmp(arg, "*") == 0)
              {
                if (! str_arg(arg, Iline, iarg++))
                  return RCF_MissingParameter;

                S.SF_Info.Lbl = NULL;
              }
              else
                S.SF_Info.Lbl = NULL;

              if (ReadValueAndUncertainty(arg, &x, NULL) != 0)
                return RCF_WrongParameter;
              S.x = x;

              if (! str_arg(arg, Iline, iarg++))
                return RCF_MissingParameter;
              if (ReadValueAndUncertainty(arg, &x, NULL) != 0)
                return RCF_WrongParameter;
              S.y = x;

              if (! str_arg(arg, Iline, iarg++))
                return RCF_MissingParameter;
              if (ReadValueAndUncertainty(arg, &x, NULL) != 0)
                return RCF_WrongParameter;
              S.z = x;

              i = str_arg(arg, Iline, iarg++);
              if (i == 0 || strcmp(arg, "*") == 0)
              {
                S.Occ = 0.;
                S.F_Occ = 0;
              }
              else
              {
                if (ReadValueAndUncertainty(arg, &x, NULL) != 0 || x <= 0.)
                  return RCF_WrongParameter;
                S.Occ = x;
                S.F_Occ = 1;
              }

              if (i != 0) i = str_arg(arg, Iline, iarg++);
              if (i == 0 || strcmp(arg, "*") == 0)
              {
                S.Uiso = 0.;
                S.F_Uiso = 0;
              }
              else
              {
                if (ReadValueAndUncertainty(arg, &x, NULL) != 0 || x < 0.)
                  return RCF_WrongParameter;
                S.Uiso = x;
                S.F_Uiso = 1;
              }

              if (pass == 1)
              {
                if (iSite >= nSite)
                  return RCF_ChangingInput;

                Site[iSite].Label            = S.Label;
                Site[iSite].SF_Info.Lbl      = S.SF_Info.Lbl;
                Site[iSite].SF_Info.SFT      = NULL;
                Site[iSite].SF_Info.CAA      = NULL;
                Site[iSite].SF_Info.f_stol_0 = -1.;
                Site[iSite].x                = S.x;
                Site[iSite].y                = S.y;
                Site[iSite].z                = S.z;
                Site[iSite].Occ              = S.Occ;
                Site[iSite].F_Occ            = S.F_Occ;
                Site[iSite].Uiso             = S.Uiso;
                Site[iSite].F_Uiso           = S.F_Uiso;
              }

              iSite++;

              break;
            }

            case KW_ScatteringFactor:
            {
              if (! str_arg(arg, Iline, iarg++))
                return RCF_MissingParameter;

              if (ContKW == KW_End)
              {
                for (jSF_Tables = 0; jSF_Tables < iSF_Tables; jSF_Tables++)
                  if (strcmp(SF_Tables[jSF_Tables].Label, arg) == 0)
                    break;

                if (jSF_Tables == nSF_Tables)
                {
                  if (pass)
                    return RCF_ChangingInput;
                                                     n = nSF_Tables + 10;
                  CheckRealloc(SF_Tables, SF_Tables, n,  nSF_Tables);

                  nSF_Tables = n;
                }

                if (jSF_Tables == iSF_Tables)
                {
                  if (pass == 0)
                  {
                    AppStrdup(arg, SF_Tables[iSF_Tables].Label);
                      if (SF_Tables[iSF_Tables].Label == NULL)
                        NotEnoughCore();

                    SF_Tables[iSF_Tables].mTab  = 0;
                    SF_Tables[iSF_Tables].nTab  = 0;
                    SF_Tables[iSF_Tables].Tab   = NULL;
                  }

                  iSF_Tables++;
                }

                ContKW = KW_ScatteringFactor;
              }
              else
              {
                T_SF_Tab  Tab;

                if (sscanf(arg, "%lf%c", &x, &extrac) != 1 || x < 0.)
                  return RCF_WrongParameter;

                Tab.stol = x;

                if (! str_arg(arg, Iline, iarg++))
                  return RCF_MissingParameter;

                if (sscanf(arg, "%lf%c", &x, &extrac) != 1 || x < 0.)
                  return RCF_WrongParameter;

                Tab.f = x;

                if (pass == 1)
                {
                      n =  SF_Tables[jSF_Tables].nTab;
                  if (n == SF_Tables[jSF_Tables].mTab)
                    return RCF_ChangingInput;

                  SF_Tables[jSF_Tables].Tab[n].stol = Tab.stol;
                  SF_Tables[jSF_Tables].Tab[n].f    = Tab.f;
                }

                SF_Tables[jSF_Tables].nTab++;
              }

              break;
            }

            case KW_ScatteringFactorTable:
            {
              if (! str_arg(arg, Iline, 1))
                return RCF_MissingParameter;
              if      (str_icmp(arg, "xray") == 0)
                          i = 0;
              //else if (str_icmp(arg, "IT92") == 0)
              //            i = 1;
              else if (str_icmp(arg, "electron") == 0)
                          i = 2;
              else if (str_icmp(arg, "WK95") == 0)
                          i = 0;
              else if (str_icmp(arg, "IT4322") == 0)
                          i = 2;
              else if (str_icmp(arg, "IT4323") == 0)
                          i = 3;
              else
                return RCF_WrongParameter;
              if (pass == 1) 
              {
                ModeScatteringFactorTable = i;
              }
              break;
            }

            default:
              InternalError(NULL);
              break;
          }
        }

        else /* input: h k l Fobs [sigmaFobs [FWHM]] */
        {
          int     h, k, l;
          double  Fobs, Sigma, FWHM;

          if (sscanf(arg, "%d%c", &h, &extrac) != 1)
            return RCF_UnrecognizedLine;
          if (! str_arg(arg, Iline, 1))
            return RCF_MissingParameter;
          if (sscanf(arg, "%d%c", &k, &extrac) != 1)
            return RCF_UnrecognizedLine;
          if (! str_arg(arg, Iline, 2))
            return RCF_MissingParameter;
          if (sscanf(arg, "%d%c", &l, &extrac) != 1)
            return RCF_UnrecognizedLine;

          if (! str_arg(arg, Iline, 3))
            return RCF_MissingParameter;
          if (sscanf(arg, "%lf%c", &Fobs, &extrac) != 1)
            return RCF_WrongParameter;

          Sigma = FWHM = 0.;

          if (str_arg(arg, Iline, 4))
          {
            if (strcmp(arg, "*") != 0)
            {
              if (sscanf(arg, "%lf%c", &Sigma, &extrac) != 1 || Sigma < 0.)
                return RCF_WrongParameter;
            }

            if (str_arg(arg, Iline, 5))
            {
              if (strcmp(arg, "*") != 0)
              {
                if (sscanf(arg, "%lf%c", &FWHM, &extrac) != 1 || FWHM < 0.)
                  return RCF_WrongParameter;
              }
            }
          }

          if (pass == 1)
          {
            if (iFobsRaw >= nFobsRaw)
              return RCF_ChangingInput;

            FobsRaw[iFobsRaw].h = h;
            FobsRaw[iFobsRaw].k = k;
            FobsRaw[iFobsRaw].l = l;
            FobsRaw[iFobsRaw].Fobs = Fobs;
            FobsRaw[iFobsRaw].SigmaFobs = Sigma;
            FobsRaw[iFobsRaw].FWHM = FWHM;
            FobsRaw[iFobsRaw].Overlap = 0;
            FobsRaw[iFobsRaw].Q = 0.;
            FobsRaw[iFobsRaw].M = 0;
            FobsRaw[iFobsRaw].PhaseRestriction = -1;
            FobsRaw[iFobsRaw].uvw[0] = 0;
            FobsRaw[iFobsRaw].uvw[1] = 0;
            FobsRaw[iFobsRaw].uvw[2] = 0;
            FobsRaw[iFobsRaw].Status = -1;
            FobsRaw[iFobsRaw].Fmrg = 0.;
            FobsRaw[iFobsRaw].SigmaFmrg = 0.;
          }

          iFobsRaw++;
        }
      }
    }

    if (pass == 0)
    {
      if (nSF_Tables != iSF_Tables)
        CheckRealloc(SF_Tables, SF_Tables, iSF_Tables, nSF_Tables);

      nRandomBlindCalls = iRandomBlindCalls;
      nAtomType = iAtomType;
      nNodeTypes = iNodeTypes;
      nFeedBackCycles = iFeedBackCycles;
      nInputChemistry = iInputChemistry;
      nSite = iSite;
      nSF_Tables = iSF_Tables;
      nFobsRaw = iFobsRaw;

      if (nRandomBlindCalls > 0)
        CheckMalloc(RandomBlindCalls, nRandomBlindCalls);

      if (nAtomType > 0)
        CheckMalloc(AtomType, nAtomType);

      if (nNodeTypes > 0)
        CheckMalloc(NodeTypes, nNodeTypes);

      if (nFeedBackCycles > 0)
        CheckMalloc(FeedBackCycles, nFeedBackCycles);

      if (nInputChemistry > 0)
        CheckMalloc(InputChemistry, nInputChemistry);

      if (nSite > 0)
        CheckMalloc(Site, nSite);

      n = 0;

      for (iSF_Tables = 0; iSF_Tables < nSF_Tables; iSF_Tables++)
      {
        n += SF_Tables[iSF_Tables].nTab;

             SF_Tables[iSF_Tables].mTab
           = SF_Tables[iSF_Tables].nTab;

             SF_Tables[iSF_Tables].nTab = 0;
      }

      if (n)
        CheckMalloc(SF_Tables[0].Tab, n);

           jSF_Tables = 0;
      for (iSF_Tables = 1; iSF_Tables < nSF_Tables; iSF_Tables++, jSF_Tables++)
            SF_Tables[iSF_Tables].Tab
        =   SF_Tables[jSF_Tables].Tab
          + SF_Tables[jSF_Tables].mTab;

      if (nFobsRaw > 0)
        CheckMalloc(FobsRaw, nFobsRaw);
    }
    else
    {
      if (   iRandomBlindCalls != nRandomBlindCalls
          || iAtomType         != nAtomType
          || iNodeTypes        != nNodeTypes
          || iFeedBackCycles   != nFeedBackCycles
          || iInputChemistry   != nInputChemistry
          || iSite             != nSite
          || iSF_Tables        != nSF_Tables
          || iFobsRaw          != nFobsRaw )
        return RCF_ChangingInput;

      for (iSF_Tables = 0; iSF_Tables < nSF_Tables; iSF_Tables++)
        if (   SF_Tables[iSF_Tables].mTab
            != SF_Tables[iSF_Tables].nTab)
          return RCF_ChangingInput;
    }
  }

  return RCF_OK;
}


static void EchoBreakIf(FILE *fpout, T_BreakIf *BreakIf, char *label)
{
  if (BreakIf->Branch)
  {
    Fprintf(fpout, "%s  ", label);

    if (BreakIf->Branch & BrkIf_Branch_PhaseDiff)
    {
      Fprintf(fpout, "PhaseDiff <");

      if (BreakIf->PhaseDiff.Percent)
        Fprintf(fpout, " %.2f %%",      BreakIf->PhaseDiff.Value * 100.);
      else
        Fprintf(fpout, " %d",     (int) BreakIf->PhaseDiff.Value);
    }
    if (   BreakIf->Branch & BrkIf_Branch_PhaseDiff
        && BreakIf->Branch & BrkIf_Branch_DeltaR)
    {
      if (BreakIf->Branch & BrkIf_Branch_Or)
        Fprintf(fpout, " or ");
      else
        Fprintf(fpout, " and ");
    }
    if (BreakIf->Branch & BrkIf_Branch_DeltaR)
    {
      if (BreakIf->DeltaR.Weighted)
        Fprintf(fpout, "DeltaRw <");
      else
        Fprintf(fpout, "DeltaR <");

      if (BreakIf->DeltaR.Percent)
        Fprintf(fpout, " %.2f %%", BreakIf->DeltaR.Value * 100.);
      else
        Fprintf(fpout, " %.4f",    BreakIf->DeltaR.Value);
    }
    Fprintf(fpout, "\n");
  }
}


char *EchoAtomTypeClass(FILE *fpout, int Class)
{
  char  *s = NULL;


  switch (Class)
  {
    case ATC_General:    s = "*         "; break;
    case ATC_Node:       s = "Node      "; break;
    case ATC_NodeBridge: s = "NodeBridge"; break;
    default:
      InternalError("Corrupt AtomTypeClass");
      break;
  }

  if (fpout)
    Fprintf(fpout, "%s", s);

  return s;
}


void EchoProfileSettings(FILE *fpout)
{
  const char  *s;


  Fprintf(fpout, "ProfileStartEndStep  %.6g", ProfileStart);
  if (ProfileEnd == 0.)
    Fprintf(fpout, " *");
  else
    Fprintf(fpout, " %.6g", ProfileEnd);
  Fprintf(fpout, " %.6g", ProfileStep);
  putc('\n', fpout);

  Fprintf(fpout, "ProfilePOLRA  %.6g\n", ProfilePOLRA);

  Fprintf(fpout, "ProfileFWHM");
  if (ProfileFWHM.Valid)
    Fprintf(fpout, "  UVW  %.6g %.6g %.6g",
      ProfileFWHM.U, ProfileFWHM.V, ProfileFWHM.W);
  putc('\n', fpout);

  Fprintf(fpout, "ProfileAsym");
  if (ProfileAsym.Valid)
    Fprintf(fpout, "  a(i)  %.6g %.6g %.6g",
      ProfileAsym.a1, ProfileAsym.a2, ProfileAsym.a3);
  putc('\n', fpout);

  switch (ProfilePeakShape)
  {
    case PPS_PseudoVoigt: s = "PseudoVoigt"; break;
    case PPS_StdPeak:     s = "StdPeak";     break;
    default:
      InternalError("Corrupt ProfilePeakShape");
      break;
  }
  Fprintf(fpout, "ProfilePeakShape  %s\n", s);

  Fprintf(fpout, "PseudoVoigtPeakRange  %.6g\n", PseudoVoigtPeakRange);
  Fprintf(fpout, "PseudoVoigtFracLorentz  %.6g\n", PseudoVoigtFracLorentz);

  Fprintf(fpout, "ProfileBackground  %.6g\n", ProfileBackground);

  Fprintf(fpout, "ProfileReferenceRefl");
  switch (ProfileReferenceRefl.Mode)
  {
    case PRRM_Free: break;
    case PRRM_hkl:
      Fprintf(fpout, "  %d %d %d",
        ProfileReferenceRefl.h,
        ProfileReferenceRefl.k,
        ProfileReferenceRefl.l);
      break;
    case PRRM_TTheta:
      Fprintf(fpout, "  %.6g", ProfileReferenceRefl.TTheta);
      break;
    default:
      InternalError("Corrupt ProfileReferenceRefl.Mode");
      break;
  }
  putc('\n', fpout);

  Fprintf(fpout, "ProfileReferenceMax  %.6g\n", ProfileReferenceMax);

  putc('\n', fpout);
}


void EchoInput(FILE *fpout, int Full)
{
  int         i, j, n, si, iFR, PrevOverlap;
  const char  *cp;
  Fprec       d;
  T_FobsRaw   *FR, *EquivFR;


  Fprintf(fpout, "Title");
  if (TopTitle != NULL) Fprintf(fpout, " %s", TopTitle);
  Fprintf(fpout, "\n");

  Fprintf(fpout, "SpaceGroup  ");
  if (SpgrInfo.TabSgName)
  {
    (void) PrintFullHM_SgName(SpgrInfo.TabSgName, ' ', 0, fpout);
    if (Full) {
      Fprintf(fpout, "  #  ");
      PrintTabSgNameEntry(SpgrInfo.TabSgName, 0, 0, 0, fpout);
    }
  }
  else if (SpgrInfo.HallSymbol[0])
    Fprintf(fpout, "Hall  %s", SpgrInfo.HallSymbol);
  else
    Fprintf(fpout, "Unknown");

  putc('\n', fpout);

  Fprintf(fpout, "UnitCell  %.6g %.6g %.6g  %.6g %.6g %.6g",
    LatConD.a,
    LatConD.b,
    LatConD.c,
    LatConD.alpha / PIover180,
    LatConD.beta  / PIover180,
    LatConD.gamma / PIover180);
  if (Full) {
    Fprintf(fpout, "  #  %.6g", LatConD.v);
    putc('\n', fpout);
    Fprintf(fpout, "# Min/MaxLatticeTr  %.6g  %.6g",
      AppSqrt(MinLatticeTr2), AppSqrt(MaxLatticeTr2));
  }
  putc('\n', fpout);
  putc('\n', fpout);

  if (nAtomType == 0)
  {
    Fprintf(fpout, "# AtomType Use Class ScatFact #PerUnitCell"
                   " OccDefault UisoDefault Radius\n");
    Fprintf(fpout, "# AtomType  +  Node  Si  8\n");
  }
  else
  {
    Fprintf(fpout, "#AtomTypes  %d\n", nAtomType);

    Fprintf(fpout, "#         Use Class ScatFact #PerUnitCell"
                   " OccDefault UisoDefault Radius\n");

    for (i = 0; i < nAtomType; i++)
    {
      Fprintf(fpout, "AtomType  %c  %s  %-4s  %4d  %5.2f  %7.3f  ",
        (AtomType[i].On ? '+' : '-'),
        EchoAtomTypeClass(NULL, AtomType[i].Class),
        AtomType[i].SF_Info.Lbl,
        AtomType[i].nPerUnitCell,
        AtomType[i].OccDefault,
        AtomType[i].UisoDefault);

      if (AtomType[i].eListAR == NULL)
        Fprintf(fpout, "*     ");
      else
        Fprintf(fpout, "%6.3f", AtomType[i].eListAR->Radius);

      Fprintf(fpout, " # e = %5.4g", AtomType[i].SF_Info.f_stol_0);

      if (AtomType[i].SF_Info.CAA)
      {
        const T_PSE *ePSE = FindInPSE(AtomType[i].SF_Info.CAA->Label, 0);
        if    (ePSE)
          Fprintf(fpout, "  Z = %3d", ePSE->Z);
      }

      putc('\n', fpout);
    }
  }

  if (Full)
  {
    Fprintf(fpout, "# F000cal = %.6g  F000mrg = %.6g  F000mtp = %.6g\n",
      F000cal, F000mrg, F000mtp);
    Fprintf(fpout, "# MaxFmrg scaled = %.6g  unscaled = %.6g\n",
      MaxFmrg, MaxFmrg / FobsScale);
  }

  putc('\n', fpout);

  if (nInputChemistry == 0)
  {
    Fprintf(fpout, "# Chemistry  MinDistance  Node  Si  Node  Si  2.6\n");
  }
  else
  {
    for (i = 0; i < nInputChemistry; i++)
    {
      switch (InputChemistry[i].Type)
      {
        case ICT_MinDistance: cp = "MinDistance"; n = 2; break;
        case ICT_MaxDistance: cp = "MaxDistance"; n = 2; break;
        case ICT_AddBondAtom: cp = "AddBondAtom"; n = 3; break;
        default:
          InternalError("Corrupt InputChemistry Type");
          break;
      }

      Fprintf(fpout, "Chemistry  %s", cp);

      for (j = 0; j < n; j++)
      {
        Fprintf(fpout, "  %s  %-4s",
          EchoAtomTypeClass(NULL, InputChemistry[i].AtomID[j].Class),
          InputChemistry[i].AtomID[j].SF_Info.Lbl);
      }

      if (n == 2)
        Fprintf(fpout, "  %.6g", InputChemistry[i].Value);

      putc('\n', fpout);
    }
  }

  Fprintf(fpout, "MaxPotentialAtoms  %d\n", MaxPotentialAtoms);
  Fprintf(fpout, "MaxRecycledAtoms  %d\n", MaxRecycledAtoms);
  putc('\n', fpout);

  Fprintf(fpout, "FwSearchMethod  ");
  switch (FwSearchMethod)
  {
    case FwSM_FwTracking:
      Fprintf(fpout, "FwTracking\n");
      break;
    case FwSM_AltFwTracking:
      Fprintf(fpout, "AltFwTracking\n");
      break;
    default:
      InternalError("Corrupt FwSearchMethod");
  }

  Fprintf(fpout, "MaxPeaksFwSearch  %d\n", MaxPeaksFwSearch);
  Fprintf(fpout, "MaxPeaksFwFragmentSearch  %d\n", MaxPeaksFwFragmentSearch);

  Fprintf(fpout, "MinNodeDistance  %.6g\n", AppSqrt(MinNodeDist2));
  Fprintf(fpout, "MaxNodeDistance  %.6g\n", AppSqrt(MaxNodeDist2));

  Fprintf(fpout, "MinSymNodes  %d\n", MinSymNodes);
  Fprintf(fpout, "MaxSymNodes  ");
    if (MaxSymNodes < 0) Fprintf(fpout, "*\n");
    else                 Fprintf(fpout, "%d\n", MaxSymNodes);

  for (i = 0; i < nNodeTypes; i++)
  {
    Fprintf(fpout, "NodeType  %d  ", NodeTypes[i].CN);

    if (NodeTypes[i].MaxNpAU)
      Fprintf(fpout, "%d", NodeTypes[i].MaxNpAU);
    else
      putc('*', fpout);

    n = 1;

    for (j = -6; j <= 6; j++) {
      if (NodeTypes[i].NoSPO[j + 6]) {
        if (n) { n = 0; putc(' ', fpout); }
        Fprintf(fpout, " %d", j);
      }
    }

    putc('\n', fpout);
  }

  Fprintf(fpout, "MinLoopSize  %d\n", MinLoopSize);
  Fprintf(fpout, "MaxLoopSize  %d\n", MaxLoopSize);

  switch (EvenLoopSizesOnly)
  {
    case 0: cp = "Off"; break;
    case 1: cp = "On";  break;
    default:
      InternalError("Corrupt EvenLoopSizesOnly");
  }
  Fprintf(fpout, "EvenLoopSizesOnly  %s\n", cp);

  switch (Check3DimConnectivity)
  {
    case 0: cp = "Off"; break;
    case 1: cp = "On";  break;
    default:
      InternalError("Corrupt Check3DimConnectivity");
  }
  Fprintf(fpout, "Check3DimConnectivity  %s\n", cp);

  Fprintf(fpout, "IdealT_NodeDistance  %.6g\n",
          AppSqrt(IdealT_NodeDist2));
  if (Full) {
    Fprintf(fpout, "IdealTetrahedralEdge  %.6g\n",
            AppSqrt(IdealTetrahedralEdge2));
    Fprintf(fpout, "IdealTetrahedralVolume  %.6g  %.6g\n",
            AppSqrt(IdealTetrahedralVol2),
                    IdealTetrahedralVolFrac);
  }

  switch (ModeCheckTetrahedra)
  {
    case 0: cp = "Off";    break;
    case 1: cp = "Normal"; break;
    case 2: cp = "Hard";   break;
    default:
      InternalError("Corrupt ModeCheckTetrahedra");
  }
  Fprintf(fpout, "CheckTetrahedralGeometry  %s\n", cp);

  putc('\n', fpout);

  Fprintf(fpout, "RandomInitialization  ");
  if (FlagTimeRandomInitialization) Fprintf(fpout, "Time");
  if (Full)
    Fprintf(fpout, "  #  %ld\n", RandomInitialization);

  if ((RandomBlindCalls && nRandomBlindCalls > 0) || Full)
    Fprintf(fpout, "RandomBlindCalls");

  if  (RandomBlindCalls && nRandomBlindCalls > 0)
  {
    n = 0;

    putc(' ', fpout);
    for (i = 0; i < nRandomBlindCalls; i++)
    {
      if (n++ == 5) {
        Fprintf(fpout, "\nRandomBlindCalls ");
        n = 1;
      }

      Fprintf(fpout, " %ld", RandomBlindCalls[i]);
    }
  }
  putc('\n', fpout);

  Fprintf(fpout, "FeedBackCycles");
  if (FeedBackCycles && nFeedBackCycles > 0)
  {
    n = 0;

    putc(' ', fpout);
    for (i = 0; i < nFeedBackCycles; i++)
    {
      if (n++ == 20) {
        Fprintf(fpout, "\nFeedBackCycles ");
        n = 1;
      }

      Fprintf(fpout, " %d", FeedBackCycles[i]);
    }
  }
  putc('\n', fpout);

  EchoBreakIf(fpout, &FeedBackBreakIf, "FeedBackBreakIf");

  putc('\n', fpout);

  Fprintf(fpout, "Grid_xyz  %d  %d  %d", Nx, Ny, Nz);
  if (Full)
    Fprintf(fpout, "  # Resolution xyz %.6g %.6g %.6g",
      LatConD.a / Nx,
      LatConD.b / Ny,
      LatConD.c / Nz);
  putc('\n', fpout);

  Fprintf(fpout, "PeakSearchLevel  %d\n", PeakSearchLevel);

  Fprintf(fpout, "eDensityCutOff  %.6g", eDensityCutOff.Value);
  if (eDensityCutOff.Percent) Fprintf(fpout, " %%");
  putc('\n', fpout);

  Fprintf(fpout, "MinPfI  %d\n", MinPfI);

  Fprintf(fpout, "CatchDistance  %.6g\n", AppSqrt(CatchDistance2));

  switch (eD_PeaksSortElement)
  {
    case PSE_Grid_eD:  cp = "Grid_eD"; break;
    case PSE_Maximum:  cp = "Maximum"; break;
    case PSE_Integral: cp = "Integral"; break;
    default:
      InternalError("Corrupt eD_PeaksSortElement");
  }
  Fprintf(fpout, "eD_PeaksSortElement  %s\n", cp);

  putc('\n', fpout);

  Fprintf(fpout, "Lambda");
  if (LambdaName[0] != '\0') Fprintf(fpout, "  %s", LambdaName);
  Fprintf(fpout, "  %.6g\n", LambdaLength);

  Fprintf(fpout, "FobsMin_d  ");
  if (FobsMaxQ > 0.) d = 1. / AppSqrt(FobsMaxQ);
  else               d = 0.;
  Fprintf(fpout, "%.6g\n", d);

  Fprintf(fpout, "FobsScale  %.6g # Fabs = Fobs * FobsScale\n", FobsScale);
  Fprintf(fpout, "SigmaCutOff  %.6g\n", SigmaCutOff);

  Fprintf(fpout, "OverlapFactor  %.6g\n", OverlapFactor);
  switch(OverlapAction)
  {
    case OvlAct_NoAction: cp = "NoAction"; break;
    case OvlAct_EqualF2:  cp = "EqualF2";  break;
    case OvlAct_EqualMF2: cp = "EqualMF2"; break;
    case OvlAct_Omit:     cp = "Omit";     break;
    default:
      InternalError("Corrupt OverlapAction");
      break;
  }
  Fprintf(fpout, "OverlapAction  %s\n", cp);

  Fprintf(fpout, "ReflectionUsage  %.6g", ReflectionUsage.Value);
  if (ReflectionUsage.Percent) Fprintf(fpout, " %%");
  putc('\n', fpout);

  putc('\n', fpout);

  if (ModeScatteringFactorTable == 0)
    Fprintf(fpout, "ScatteringFactorTable xray # WK95\n");
  else if (ModeScatteringFactorTable == 2)
    Fprintf(fpout, "ScatteringFactorTable electron # IT4322\n");
  else if (ModeScatteringFactorTable == 3)
    Fprintf(fpout, "ScatteringFactorTable IT4323 # electron\n");

  putc('\n', fpout);

  Fprintf(fpout, "GenerateFWHM");
  if (GenerateFWHM.Valid)
    Fprintf(fpout, "  UVW  %.6g %.6g %.6g",
      GenerateFWHM.U, GenerateFWHM.V, GenerateFWHM.W);
  putc('\n', fpout);

  Fprintf(fpout, "AbsentRedivisionLimit  ");
  if (AbsentRedivisionLimit.Type == ARLT_FWHM)
    Fprintf(fpout, "%.6g FWHM\n",         AbsentRedivisionLimit.Value);
  else
    Fprintf(fpout, "%.6g Degree2Theta\n", AbsentRedivisionLimit.Value);

  putc('\n', fpout);

  if (Full)
  {
    Fprintf(fpout, "tRius  %.6g\n", tRius);

    Fprintf(fpout, "MinConsecutiveNodes  ");
    if (MinConsecutiveNodes < 0)
      Fprintf(fpout, "NoDisplay");
    else
      Fprintf(fpout, "%d", MinConsecutiveNodes);
    putc('\n', fpout);

    putc('\n', fpout);
  }

  EchoProfileSettings(fpout);

  Fprintf(fpout, "NormalizeSurface");
  for (i = 0; i < nNormalizeSurface; i++)
    Fprintf(fpout, "  %.6g", vNormalizeSurface[i]);
  if (i == 0)
    Fprintf(fpout, "  #  [[Min] Max]");
  putc('\n', fpout);

  putc('\n', fpout);

  if (nSite == 0)
  {
    Fprintf(fpout,
      "#  Site  Label  ScatFact  x  y  z  Occ  Uiso\n");
    Fprintf(fpout,
      "#  Site  Si4  Si4+  %8.5f  %8.5f  %8.5f   %5.2f  %7.3f\n",
      0., 0., 0., Df_OccDefault, Df_UisoDefault);
  }
  else
  {
    Fprintf(fpout, "#Sites  %d\n", nSite);

    Fprintf(fpout, "#     Label  ScatFact  x  y  z  Occ  Uiso\n");

    for (i = 0; i < nSite; i++)
    {
          cp =  Site[i].SF_Info.Lbl;
      if (cp == NULL)
          cp = "";

      Fprintf(fpout, "Site  %-6s %-4s  %8.5f  %8.5f  %8.5f",
        Site[i].Label, cp,
        Site[i].x, Site[i].y, Site[i].z);

      if (Site[i].F_Occ == 0)
        Fprintf(fpout, "     *  ");
      else
        Fprintf(fpout, "   %5.2f", Site[i].Occ);

      if (Site[i].F_Uiso == 0)
        Fprintf(fpout, "     *\n");
      else
        Fprintf(fpout, "  %7.3f\n", Site[i].Uiso);
    }
  }

  putc('\n', fpout);

  for (i = 0; i < nSF_Tables; i++)
  {
    Fprintf(stdout, "ScatteringFactor  %s  #  %d\n",
      SF_Tables[i].Label, SF_Tables[i].nTab);

    for (j = 0; j < SF_Tables[i].nTab; j++)
    {
      Fprintf(stdout, " %8.6f %10.6f\n",
        SF_Tables[i].Tab[j].stol,
        SF_Tables[i].Tab[j].f);
    }

    Fprintf(stdout, "&\n");
  }

  if (nSF_Tables)
    putc('\n', fpout);

  if (Full)
    Fprintf(fpout, "#Fobs  %d\n", nFobsRaw);

  Fprintf(fpout, "#  h   k   l       Fobs      Sigma     FWHM\n");

  PrevOverlap = 0;
  EquivFR = FR = FobsRaw;

  for (iFR = 0; iFR < nFobsRaw; iFR++)
  {
    if (FR->Q > 0.) d = 1. / AppSqrt(FR->Q);
    else            d = 0.;

    if (FR->Overlap == 0)
      i = ' ';
    else if (FR->Overlap == PrevOverlap)
      i = '|';
    else
      i = '_';

    PrevOverlap = FR->Overlap;

    Fprintf(fpout, " %3d %3d %3d %10.4f %10.4f %8.5f # %6.2f%c %7.4f %2d ",
      FR->h, FR->k, FR->l,
      FR->Fobs, FR->SigmaFobs,
      FR->FWHM,
      TwoThetaDeg(FR->Q), i, d, FR->M);

    if (   (FR->Status == FRS_SymEquiv && EquivFR->Status != FRS_SysAbsent)
        || (FR->Status != FRS_SymEquiv &&      FR->Status != FRS_SysAbsent))
    {
      si = Is_ss(&SpgrInfo, FR->h, FR->k, FR->l);
      if (si == 0) si = '/';
      else         si = ':';
    }
    else
      si = ' ';

    if (FR->PhaseRestriction >= 0)
    {
      i = FR->PhaseRestriction * (180 / STBF);
      Fprintf(fpout, "%3d%c%3d ", i, si, i + 180);
    }
    else
      Fprintf(fpout, "   %c    ", si);

    for (i = 0; i < SpgrInfo.n_ssVM; i++)
    {
      if (si != ' ')
        Fprintf(fpout, "%3d ", FR->uvw[i]);
      else
        Fprintf(fpout, "    ");
    }

    switch (FR->Status)
    {
      case FRS_Undefined: cp = "Undefined"; break;
      case FRS_F000:      cp = "F000";      break;
      case FRS_SysAbsent: cp = "SysAbsent"; break;
      case FRS_SymEquiv:  cp = "==";        break;
      case FRS_Omit:      cp = "Omit";      break;
      case FRS_OffGrid:   cp = "OffGrid";   break;
      case FRS_Active:    cp = "Active";    break;
      case FRS_Sleeping:  cp = "Sleeping";  break;
      default:
        InternalError("Illegal Fobs Status");
    }

    Fprintf(fpout, "%s\n", cp);

    if (FR->Status != FRS_SymEquiv)
    {
      if (FR->Status != FRS_F000)
      {
        Fprintf(fpout, "#            %10.4f %10.4f\n",
          FR->Fmrg / FobsScale, FR->SigmaFmrg / FobsScale);
      }

      EquivFR = FR;
    }
    else
    {
      if (FR->Fmrg != 0. || FR->SigmaFmrg != 0.)
        InternalError("Corrupt FR->[Sigma]Fmrg");
    }

    FR++;
  }

  Fprintf(fpout, "End\n");
  Fflush(fpout);
}


void PrintFmrg(FILE *fpout)
{
  int        iFR;
  T_FobsRaw  *FR;


  Fprintf(fpout, ">Begin Fmrg\n");

  Fprintf(fpout, "#  h   k   l       Fmrg      Sigma     FWHM\n");

  for (FR = FobsRaw, iFR = 0; iFR < nFobsRaw; iFR++, FR++)
  {
    if (FR->Status != FRS_SymEquiv && FR->Status != FRS_F000)
    {
      Fprintf(fpout, " %3d %3d %3d %10.4f %10.4f %8.5f\n",
        FR->h, FR->k, FR->l,
        FR->Fmrg / FobsScale, FR->SigmaFmrg / FobsScale,
        FR->FWHM);
    }
  }

  Fprintf(fpout, ">End Fmrg\n\n");
}


void Print_eEPL(FILE *fpout, T_eD_PeakList *EPL, int iEPL,
                int CoordinationNumber, int FlagAll)
{
  Fprec  x, y, z;


  x = EPL[iEPL].Position.x;
  y = EPL[iEPL].Position.y;
  z = EPL[iEPL].Position.z;

  Fprintf(fpout, " %8.5f %8.5f %8.5f", x, y, z);

  if (CoordinationNumber >= 0)
    Fprintf(fpout, " %2d", CoordinationNumber);

  if (FlagAll)
    Fprintf(fpout, " # %2d %8.3f %8.3f %8.3f %3d %2d %2d %c [%d]",
      iEPL,
      EPL[iEPL].Grid_eD,
      EPL[iEPL].Maximum,
      EPL[iEPL].Integral,
      EPL[iEPL].WL_Entry->nPositions,
      EPL[iEPL].nPfI,
      EPL[iEPL].iAtomType,
      EPL[iEPL].WL_Entry->CanBeNode != 0 ? 'N' : '-',
      (int)(EPL[iEPL].WL_Entry - WyckoffList));
}


static void PrintStrudatSeperatorLine(FILE *fpout)
{
  int  i;


  for (i = 0; i < 33; i++)
    putc('-', fpout);

  putc('\n', fpout);
}


void DoPut_strudat(T_eD_PeakList *EPL, int nEPL, int iCycle, char *Label,
                   int FlagAll)
{
  int          iAT, iEPL;
  int          *LbCounter;
  char         buf[80], *fnout;
  static FILE  *fpout = NULL;


  CheckMalloc(LbCounter, uAtomType);

  if (fpout == NULL)
  {
    if (F_Put_strudatFileName == NULL)
      fpout = stdout;
    else
    {
      fnout = F_Put_strudatFileName;

      fpout = fopen(fnout, "w");
      if (fpout == NULL)
      {
        Fprintf(stdout, "DoPut_strudat(): Can't open %s\n", fnout);
        fpout = stdout;
      }
    }
  }

  if (fpout == stdout)
    Fprintf(fpout, "\n>Begin strudat\n");

  Fprintf(fpout, "*c%d\n", iCycle);

  if (TopTitle != NULL) Fprintf(fpout, "%s\n", TopTitle);
  else                  Fprintf(fpout, "Title\n");

  if (Label != NULL) Fprintf(fpout, "%s\n", Label);
  else               Fprintf(fpout, "Label\n");

  if (SpgrInfo.TabSgName)
    (void) PrintFullHM_SgName(SpgrInfo.TabSgName, ' ', 0, fpout);
  else if (SpgrInfo.HallSymbol[0])
    Fprintf(fpout, "Hall %s", SpgrInfo.HallSymbol);
  else
    Fprintf(fpout, "Unknown");

  putc('\n', fpout);

  Fprintf(fpout, "  %9.5f  %9.5f  %9.5f  %7.3f  %7.3f  %7.3f\n",
    LatConD.a,
    LatConD.b,
    LatConD.c,
    LatConD.alpha / PIover180,
    LatConD.beta  / PIover180,
    LatConD.gamma / PIover180);


  for (iAT = 0; iAT < uAtomType; iAT++) LbCounter[iAT] = 0;

  for (iEPL = 0; iEPL < nEPL; iEPL++)
  {
        iAT = EPL[iEPL].iAtomType;
    if (iAT >= 0 || FlagAll)
    {
      if (iAT >= 0)
        Sprintf(buf, "%.40s%d", AtomType[iAT].SF_Info.Lbl, ++LbCounter[iAT]);
      else
        Sprintf(buf, "Pu%d", iEPL);

      Fprintf(fpout, "%-8s", buf);
      Print_eEPL(fpout, EPL, iEPL, -1, FlagAll);
      putc('\n', fpout);
    }
  }

  PrintStrudatSeperatorLine(fpout);

  if (fpout == stdout) Fprintf(fpout, ">End strudat\n\n");

  AppFree(LbCounter, uAtomType);
}


void PutHistogram(int *Box, int nBox, char *legend)
{
  int           iBox, *bx, MaxY;
  static FILE   *fpout = NULL;
  char          *fnout;


  MaxY = 0;
  for (iBox = 0, bx = Box; iBox < nBox; iBox++, bx++)
    if (MaxY < *bx) MaxY = *bx;

  if (fpout == NULL)
  {
    fnout = F_eD_HistogramFileName;

    if (fnout == NULL)
      fpout = stdout;
    else
    {
      fpout = fopen(fnout, "w");
      if (fpout == NULL)
      {
        Fprintf(stdout, "Put_Histogram(): Can't open %s\n", fnout);
        fpout = stdout;
      }
    }
  }

  Fprintf(fpout, "Frame ");
  if (TopTitle != NULL) Fprintf(fpout, " %s", TopTitle);
  if (legend   != NULL) Fprintf(fpout, " %s", legend);
  Fprintf(fpout, "\n");

  Fprintf(fpout, "DimensionX 0 %d\n", nBox - 1);
  Fprintf(fpout, "DimensionY 0 %d\n", MaxY);

  for (iBox = 0, bx = Box; iBox < nBox; iBox++, bx++)
    Fprintf(fpout, "%d %d\n", iBox, *bx);

  Fprintf(fpout, "End\n");

  Fflush(fpout);
}


void Put_eDmap(const Fprec *eD, char *note)
{
  int           ix, iy, iz, ieD;
  static FILE   *fpout = NULL;
  char          *fnout;


  if (fpout == NULL)
  {
    fnout = F_Put_eDmapFileName;

    if (fnout == NULL)
      fpout = stdout;
    else
    {
      fpout = fopen(fnout, "w");
      if (fpout == NULL)
      {
        Fprintf(stdout, "Put_eDmap(): Can't open %s\n", fnout);
        fpout = stdout;
      }
    }
  }

  Fprintf(fpout, ">Begin eDmap\n");

  if (note == NULL) note = TopTitle;
  if (note == NULL) note = "eDmap";
  Fprintf(fpout, "%s\n", note);

  Fprintf(fpout, "# lattice constants\n");

  Fprintf(fpout, "  %9.5f  %9.5f  %9.5f  %7.3f  %7.3f  %7.3f\n",
    LatConD.a,
    LatConD.b,
    LatConD.c,
    LatConD.alpha / PIover180,
    LatConD.beta  / PIover180,
    LatConD.gamma / PIover180);

  Fprintf(fpout, "# Nx Ny Nz\n");
  Fprintf(fpout, "  %d  %d  %d\n", Nx, Ny, Nz);

  Fprintf(fpout, "# x y z start end\n");

  Fprintf(fpout, "  %8.5f  %8.5f\n", 0., (Nx - 1.) / Nx);
  Fprintf(fpout, "  %8.5f  %8.5f\n", 0., (Ny - 1.) / Ny);
  Fprintf(fpout, "  %8.5f  %8.5f\n", 0., (Nz - 1.) / Nz);

  Fprintf(fpout, "# x fast y medium z slow\n");

  for (iz = 0; iz < Nz; iz++)
  for (iy = 0; iy < Ny; iy++)
  for (ix = 0; ix < Nx; ix++)
  {
    ieD = (ix * Ny + iy) * Nz + iz;
    Fprintf(fpout, "%.2f\n", eD[ieD]);
  }

  Fprintf(fpout, ">End eDmap\n");
}


void PrintFramework(int Colored, const int *LnB, int Low_iEPL, int Top_iEPL,
                    const T_LargestFwFragment *LFwF)
{
  int         iEPL;
  int         LblNo[100], N;
  static int  iFf = 0;


  if (LFwF == NULL)
  {
    Fprintf(stdout, ">Begin Framework");
    if (IndicateFw) Fprintf(stdout, ":%4.4d", iFramework);
    putc('\n', stdout);

    Fprintf(stdout, "*Fw");
    if (IndicateFw) Fprintf(stdout, "%4.4d", iFramework);
    putc('\n', stdout);

    if (TopTitle != NULL) Fprintf(stdout, "%s", TopTitle);
    putc('\n', stdout);

    if      (F_nAllAsyPoints >= 0) Fprintf(stdout, "AllAsyPoints");
    else if (F_SiteFrame)          Fprintf(stdout, "SiteFrame");
    else                           Fprintf(stdout, "Framework");

    if (IndicateFw) Fprintf(stdout, " %4.4d", iFramework);

    if (F_nAllAsyPoints < 0 && F_SiteFrame == 0)
      Fprintf(stdout, " %.6g %ld %d",
        CurrCorrelationCoefficient,
        CurrPhaseCode_nCallRanmar,
        Re_i_Cycle);

    putc('\n', stdout);
  }
  else
  {
    LnB      = LFwF->LnB;
    Low_iEPL = LFwF->Low_iEPL;
    Top_iEPL = LFwF->Top_iEPL;

    Fprintf(stdout, ">Begin FwFragment\n");
    Fprintf(stdout, "*Ff%4.4d\n", iFf++);
    if (TopTitle != NULL) Fprintf(stdout, "%s", TopTitle);
    putc('\n', stdout);
    Fprintf(stdout, "FwFragment %d %d %.6g\n",
      LFwF->nAsyN, LFwF->nSymN, LFwF->mBpN);
  }

  if (SpgrInfo.TabSgName)
    (void) PrintFullHM_SgName(SpgrInfo.TabSgName, ' ', 0, stdout);
  else if (SpgrInfo.HallSymbol[0])
    Fprintf(stdout, "Hall %s", SpgrInfo.HallSymbol);
  else
    Fprintf(stdout, "Unknown");
  putc('\n', stdout);

  Fprintf(stdout, " %.6g %.6g %.6g %.6g %.6g %.6g\n",
    LatConD.a,
    LatConD.b,
    LatConD.c,
    LatConD.alpha / PIover180,
    LatConD.beta  / PIover180,
    LatConD.gamma / PIover180);

  for (N = 0; N < 100; N++) LblNo[N] = 0;

  for (iEPL = Low_iEPL; iEPL <= Top_iEPL; iEPL++)
  {
    if (LnB[iEPL])
    {
          N = eD_PeakList[iEPL].Index;
      if (N > 99)
          N = 99;

      while (N && LblNo[N]) N--;

      LblNo[N] = 1;

      if (LnB[iEPL] > 0 && (Colored == 0 || NodeColor(iEPL) == 0))
        Fprintf(stdout, "T%2.2d", N);
      else
        Fprintf(stdout, "U%2.2d", N);

      Print_eEPL(stdout, eD_PeakList, iEPL, abs(LnB[iEPL]), 1);

      putc('\n', stdout);
    }
  }

  PrintStrudatSeperatorLine(stdout);

  if (LFwF == NULL)
    Fprintf(stdout, ">End Framework\n\n");
  else
    Fprintf(stdout, ">End FwFragment\n\n");

  Fflush(stdout);
}


void LoadDensity(Fprec *Density)
{
  int     ix, iy, iz, iD;
  int     c = EOF;
  double  val;

  static FILE  *fpdens = NULL;


  if (fpdens == NULL)
  {
        fpdens = fopen(F_LoadDensityFileName, "r");
    if (fpdens == NULL)
    {
      Fprintf(stderr, "%s: Can't open density file %s\n",
        progn, F_LoadDensityFileName);
      AppExit(1);
    }
  }

  for (iz = 0; iz < Nz; iz++)
  for (iy = 0; iy < Ny; iy++)
  for (ix = 0; ix < Nx; ix++)
  {
    for (;;)
    {
      do
        c = fgetc(fpdens);
      while (c != EOF && isspace(c));

      if (c != '#')
        break;

      do
        c = fgetc(fpdens);
      while (c != EOF && c != '\n');

      if (c == EOF)
        break;
    }

    if (c != EOF)
      c = ungetc(c, fpdens);

    if (c != EOF)
      c = fscanf(fpdens, "%lf", &val);

    if (c == EOF)
    {
      if (ix == 0 && iy == 0 && iz == 0)
      {
        Fprintf(stdout, "\n# End of density file %s\n\n",
          F_LoadDensityFileName);
        AppExit(0);
      }

      Fprintf(stderr, "%s: Not enough data in density file %s\n",
        progn, F_LoadDensityFileName);
      AppExit(1);
    }
    else if (c != 1)
    {
      Fprintf(stderr, "%s: Illegal data in density file %s\n",
        progn, F_LoadDensityFileName);
      AppExit(1);
    }

            iD = (ix * Ny + iy) * Nz + iz;
    Density[iD] = (Fprec) val;
  }
}


static int iModCentre(int i, int m)
{
  if (m > 0)
  {
        i %= m;
    if (i > (m - 1) / 2)
        i -= m;
  }

  return i;
}


long LoadNextPhaseSet(int *Code)
{
  int           mode, n, h, k, l, iCT;
  double        phi;
  T_Eq_hkl      Eq_hkl[1], Eq_FRhkl[1];
  T_CodeTransl  *CT;
  T_FobsRaw     *FR;
  int           phi360, newc, iEq, Conjugate;
  char          buf[128], xtrac;
  static FILE   *fpphas = NULL;
  static long   lcount = 0;


  if (Code == NULL)
    return lcount;

  if (fpphas == NULL)
  {
        fpphas = fopen(F_PhaseSetsFileName, "r");
    if (fpphas == NULL)
    {
      Fprintf(stderr, "%s: Can't open phase sets file %s\n",
        progn, F_PhaseSetsFileName);
      AppExit(1);
    }

    lcount = 0;
  }

  for (iCT = 0; iCT < nActivePhase; iCT++) Code[iCT] = -1;

  mode = 0;

  while (fgetline(fpphas, buf))
  {
    lcount++;

    if (noblks(buf) != 0 && buf[0] != '#' && buf[0] != '$')
    {
      if (mode == 0)
      {
        if (str_ibegin(buf, ">Begin ") == 0)
        {
          (void) noblks(&buf[7]);

          if (str_icmp(&buf[7], "phase_set") == 0)
            mode = 1;
        }
      }
      else
      {
        if (str_ibegin(buf, ">End ") == 0)
        {
          (void) noblks(&buf[5]);

          if (str_icmp(&buf[5], "phase_set") == 0)
          {
            mode = 2;
            break;
          }
        }

            n = sscanf(buf, "%d%d%d%lf%c", &h, &k, &l, &phi, &xtrac);
        if (n != 4 && (n != 5 || isspace(xtrac) == 0))
          IllegalLine(F_PhaseSetsFileName, lcount);

                   phi *= 2. * M_PI;
        Phi2Phi360(phi, phi360);

        (void) BuildEq_hkl(&SpgrInfo, 1, Eq_hkl, h, k, l);
        if (SgError != NULL) progerror(SgError);

        CT = CodeTransl;

        for (iCT = 0; iCT < nActivePhase; iCT++, CT++)
        {
                                  FR = CT->FobsRaw;
          if (IndexEq_hkl(Eq_hkl, FR->h, FR->k, FR->l, NULL) >= 0)
          {
            (void) BuildEq_hkl(&SpgrInfo, 1, Eq_FRhkl, FR->h, FR->k, FR->l);

                iEq = IndexEq_hkl(Eq_FRhkl, h, k, l, &Conjugate);
            if (iEq < 0 || iEq >= CT->nCTData)
              InternalError(NULL);

            newc = Conjugate * phi360 - CT->CTData[iEq].TH;
            newc = iModPositive(newc, 360);

#define Tol 3 /* ARBITRARY */

            if (FR->PhaseRestriction >= 0)
            {
              if      (abs(iModCentre(newc -   0, 360)) <= Tol)
                newc = 0;
              else if (abs(iModCentre(newc - 180, 360)) <= Tol)
                newc = 180;
              else
                IllegalLine(F_PhaseSetsFileName, lcount);

              if      (Code[iCT] == -1)
                       Code[iCT] =  newc;
              else if (Code[iCT] != newc)
                IllegalLine(F_PhaseSetsFileName, lcount);
            }
            else
            {
              if      (Code[iCT] == -1)
                       Code[iCT] = newc;
              else if (abs(iModCentre(newc - Code[iCT], 360)) > Tol)
                IllegalLine(F_PhaseSetsFileName, lcount);
            }
#undef Tol
            break;
          }
        }

        if (iCT == nActivePhase)
          Fprintf(stdout, "# Ignoring hkl %3d %3d %3d\n", h, k, l);
      }
    }
  }

  if (mode != 2) {
    Fclose(fpphas);
    fpphas = NULL;
    lcount = 0;
  }

  if      (mode == 0) {
    Fprintf(stdout, "\n# End of phase sets file %s\n\n",
      F_PhaseSetsFileName);
    return -1;
  }
  else if (mode == 1) {
    Fprintf(stderr, "%s: Unexpected end of phase sets file %s\n",
      progn, F_PhaseSetsFileName);
    AppExit(1);
  }

  n = 0;

  for (iCT = 0; iCT < nActivePhase; iCT++)
    if (Code[iCT] < 0) n++;

  return n;
}


/* StdPeak of Anton Meden's MnPO4 sample
 */
static const char *InternalStdPeak[] =
  {
    "TITLE    STANDARD PEAK MNPO4. WHOLE RANGE",
    "PEAKF    1   51     11.00",
    "              0.0000E+00  0.3345E+00  0.0000E+00  0.0000E+00  0.3995E+00",
    "              0.2200E+00  0.3147E+00  0.9051E-01 -0.1625E+00  0.7089E+00",
    "              0.4400E+00  0.2791E+00  0.2472E+00 -0.2071E+00  0.1107E+01",
    "              0.6600E+00  0.2340E+00  0.4900E+00 -0.2057E+00  0.1414E+01",
    "              0.8800E+00  0.1889E+00  0.8008E+00 -0.1753E+00  0.1687E+01",
    "              0.1100E+01  0.1505E+00  0.1169E+01 -0.1364E+00  0.1589E+01",
    "              0.1320E+01  0.1204E+00  0.1516E+01 -0.1006E+00  0.1224E+01",
    "              0.1540E+01  0.9820E-01  0.1786E+01 -0.7457E-01  0.8947E+00",
    "              0.1760E+01  0.8174E-01  0.1983E+01 -0.5507E-01  0.6504E+00",
    "              0.1980E+01  0.6957E-01  0.2127E+01 -0.4302E-01  0.4828E+00",
    "              0.2200E+01  0.6008E-01  0.2233E+01 -0.3317E-01  0.3023E+00",
    "              0.2420E+01  0.5276E-01  0.2300E+01 -0.2670E-01  0.1875E+00",
    "              0.2640E+01  0.4686E-01  0.2342E+01 -0.2264E-01  0.1818E+00",
    "              0.2860E+01  0.4187E-01  0.2383E+01 -0.1976E-01  0.2152E+00",
    "              0.3080E+01  0.3752E-01  0.2429E+01 -0.1781E-01  0.1653E+00",
    "              0.3300E+01  0.3360E-01  0.2465E+01 -0.1615E-01  0.2492E-01",
    "              0.3520E+01  0.3004E-01  0.2471E+01 -0.1433E-01 -0.7035E-01",
    "              0.3740E+01  0.2689E-01  0.2456E+01 -0.1208E-01 -0.1038E+00",
    "              0.3960E+01  0.2424E-01  0.2433E+01 -0.1009E-01 -0.1185E+00",
    "              0.4180E+01  0.2201E-01  0.2407E+01 -0.9013E-02 -0.1369E+00",
    "              0.4400E+01  0.2002E-01  0.2377E+01 -0.8216E-02 -0.1593E+00",
    "              0.4620E+01  0.1822E-01  0.2342E+01 -0.7237E-02 -0.1748E+00",
    "              0.4840E+01  0.1662E-01  0.2303E+01 -0.6441E-02 -0.1894E+00",
    "              0.5060E+01  0.1520E-01  0.2262E+01 -0.5901E-02 -0.2099E+00",
    "              0.5280E+01  0.1390E-01  0.2215E+01 -0.5625E-02 -0.2430E+00",
    "              0.5500E+01  0.1266E-01  0.2162E+01 -0.5466E-02 -0.2868E+00",
    "              0.5720E+01  0.1146E-01  0.2099E+01 -0.5189E-02 -0.3267E+00",
    "              0.5940E+01  0.1032E-01  0.2027E+01 -0.4589E-02 -0.3528E+00",
    "              0.6160E+01  0.9315E-02  0.1950E+01 -0.3899E-02 -0.3716E+00",
    "              0.6380E+01  0.8454E-02  0.1868E+01 -0.3543E-02 -0.4128E+00",
    "              0.6600E+01  0.7673E-02  0.1777E+01 -0.3433E-02 -0.4615E+00",
    "              0.6820E+01  0.6917E-02  0.1675E+01 -0.3352E-02 -0.5134E+00",
    "              0.7040E+01  0.6180E-02  0.1562E+01 -0.3267E-02 -0.5667E+00",
    "              0.7260E+01  0.5462E-02  0.1438E+01 -0.3106E-02 -0.6155E+00",
    "              0.7480E+01  0.4780E-02  0.1303E+01 -0.2723E-02 -0.6203E+00",
    "              0.7700E+01  0.4182E-02  0.1166E+01 -0.2259E-02 -0.5970E+00",
    "              0.7920E+01  0.3684E-02  0.1035E+01 -0.1935E-02 -0.5739E+00",
    "              0.8140E+01  0.3256E-02  0.9086E+00 -0.1797E-02 -0.5620E+00",
    "              0.8360E+01  0.2861E-02  0.7850E+00 -0.1737E-02 -0.5488E+00",
    "              0.8580E+01  0.2478E-02  0.6643E+00 -0.1681E-02 -0.5273E+00",
    "              0.8800E+01  0.2109E-02  0.5484E+00 -0.1618E-02 -0.4967E+00",
    "              0.9020E+01  0.1753E-02  0.4392E+00 -0.1528E-02 -0.4517E+00",
    "              0.9240E+01  0.1417E-02  0.3399E+00 -0.1417E-02 -0.4013E+00",
    "              0.9460E+01  0.1105E-02  0.2516E+00 -0.1283E-02 -0.3454E+00",
    "              0.9680E+01  0.8232E-03  0.1757E+00 -0.1119E-02 -0.2835E+00",
    "              0.9900E+01  0.5779E-03  0.1135E+00 -0.8545E-03 -0.2027E+00",
    "              0.1012E+02  0.3905E-03  0.6894E-01 -0.5197E-03 -0.1195E+00",
    "              0.1034E+02  0.2750E-03  0.4236E-01 -0.3194E-03 -0.7034E-01",
    "              0.1056E+02  0.2040E-03  0.2670E-01 -0.2083E-03 -0.4284E-01",
    "              0.1078E+02  0.1578E-03  0.1718E-01 -0.1404E-03 -0.2642E-01",
    "              0.1100E+02  0.1275E-03  0.1146E-01  0.0000E+00  0.0000E+00",
    NULL
  };


void LoadStdPeak(const char *fnin, T_StdPeak *StdPeak)
{
  int     n, mode, lcount;
  int     nSteps, iStep;
  double  Range, R, Sym, Asy, dSym, dAsy;
  FILE    *fpin;
  char    buf[128];


  fpin = NULL;

  if (fnin)
  {
        fpin = fopen(fnin, "r");
    if (fpin == NULL)
    {
      Fprintf(stderr, "%s: Can't open standard peak file %s\n",
        progn, fnin);
      AppExit(1);
    }
  }
  else
    fnin = "InternalStdPeak[]";

  mode =
  iStep =
  lcount = 0;

  for (;;)
  {
    if (fpin) {
      if (fgetline(fpin, buf) == 0)
        break;
    }
    else {
      if (InternalStdPeak[lcount] == NULL)
        break;

      (void) str_xcpy(buf, InternalStdPeak[lcount],
                      sizeof buf / sizeof (*buf) - 1);
    }

    lcount++;

        n = firstnonblank(buf);
    if (n != '#' && n != '\0')
    {
      if      (mode == 0)
      {
        if (str_ibegin(buf, "TITLE ") != 0)
          IllegalLine(fnin, lcount);

        (void) noblks(buf + 6);
        AppStrdup(buf + 6, StdPeak->Title);
        if (StdPeak->Title == NULL) NotEnoughCore();

        mode = 1;
      }
      else if (mode == 1)
      {
        if (str_ibegin(buf, "PEAKF ") != 0)
          IllegalLine(fnin, lcount);

            n = sscanf(buf + 6, "%*d%d%lf", &nSteps, &Range);
        if (n != 2 || nSteps < 1 || Range < 1.)
          IllegalLine(fnin, lcount);

        StdPeak->nSteps = nSteps;
        StdPeak->Range  = Range;

        CheckMalloc(StdPeak->C, nSteps);

        mode = 2;
      }
      else if (mode == 2)
      {
        if (iStep >= nSteps)
        {
          Fprintf(stderr, "%s: Error: Too much data lines in file %s\n",
            progn, fnin);
          AppExit(1);
        }

            n = sscanf(buf, "%lf%lf%lf%lf%lf", &R, &Sym, &Asy, &dSym, &dAsy);
        if (n != 5)
          IllegalLine(fnin, lcount);

        StdPeak->C[iStep].R    = R;
        StdPeak->C[iStep].Sym  = Sym;
        StdPeak->C[iStep].Asy  = Asy;
        StdPeak->C[iStep].dSym = dSym;
        StdPeak->C[iStep].dAsy = dAsy;
                 ++iStep;
      }
      else
        InternalError(NULL);
    }
  }

  if (fpin)
    Fclose(fpin);

  if (mode != 2 || iStep < nSteps)
  {
    Fprintf(stderr, "%s: Error: Not enough data lines in file %s\n",
      progn, fnin);
    AppExit(1);
  }

  if (Range > StdPeak->C[nSteps - 1].R)
      Range = StdPeak->C[nSteps - 1].R;

  if (Range < 1.)
  {
    Fprintf(stderr,
      "%s: Error: Illegal peak range: Check file %s\n",
      progn, fnin);
    AppExit(1);
  }
}


static char *fmt_ticks(long ticks, double ticks_per_second)
{
  double       fseconds;
  long         hundrets, seconds, minutes, hours, days;
  static char  buf[80];


                     fseconds = (double) ticks / ticks_per_second;
  hundrets =  (long)(fseconds * 100. + .5);
  days     =  hundrets / (100 * 60 * 60 * 24);
  hundrets -= days     * (100 * 60 * 60 * 24);
  hours    =  hundrets / (100 * 60 * 60);
  hundrets -= hours    * (100 * 60 * 60);
  minutes  =  hundrets / (100 * 60);
  hundrets -= minutes  * (100 * 60);
  seconds  =  hundrets /  100;
  hundrets -= seconds  *  100;

  if      (days)

    Sprintf(buf, "%.2f = %ld:%2.2ld:%2.2ld:%2.2ld.%2.2ld", fseconds,
      days, hours, minutes, seconds, hundrets);

  else if (hours)

    Sprintf(buf, "%.2f = %ld:%2.2ld:%2.2ld.%2.2ld", fseconds,
            hours, minutes, seconds, hundrets);

  else if (minutes >= 10)

    Sprintf(buf, "%.2f = %ld:%2.2ld.%2.2ld", fseconds,
                     minutes, seconds, hundrets);
  else
    Sprintf(buf, "%.2f", fseconds);

  return buf;
}


void CoseqProtocol(int iEPL, int iSphere, int nNew)
{
  static T_Ticks  TickBuf1, TickBuf2;
  static T_Ticks  *CurrTicks = &TickBuf1;
  static T_Ticks  *LastTicks = NULL;

  static FILE  *fpcp = NULL;

  long  ElapsedCPUticks;
  long  ElapsedSeconds;


  (void) GetTicks(CurrTicks);

  if (iEPL >= 0)
  {
    if (LastTicks)
    {
      ElapsedCPUticks =  CurrTicks->cpu_user   - LastTicks->cpu_user;
      ElapsedCPUticks += CurrTicks->cpu_system - LastTicks->cpu_system;

      ElapsedSeconds =   CurrTicks->sec_since_some_day
                       - LastTicks->sec_since_some_day;
    }
    else
    {
      ElapsedCPUticks = 0;
      ElapsedSeconds  = 0;
    }

    if (fpcp == NULL)
    {
          fpcp = fopen(F_CoseqProtocolFileName, "a");
      if (fpcp == NULL)
      {
        Fprintf(stdout,
          "\n\n# WARNING: Cannot append to %s: giving up protocol\n\n",
          F_CoseqProtocolFileName);
        CountWarnings++;

        AppFree(F_CoseqProtocolFileName, strlen(F_CoseqProtocolFileName) + 1);
                F_CoseqProtocolFileName = NULL;

        return;
      }
    }

    Fprintf(fpcp, "T%3.3d(%d) %d [%s] [%ld]",
      iEPL, iSphere + 1, nNew,
      fmt_ticks(ElapsedCPUticks, CurrTicks->ticks_per_second),
      ElapsedSeconds);

    PrintTicks(fpcp, CurrTicks, " [", "]\n");

    if (iSphere >= 20 - 1) /* ARBITRARY */
      Fclose(fpcp), fpcp = NULL;
  }

  if (CurrTicks == &TickBuf1)
  {
    CurrTicks = &TickBuf2;
    LastTicks = &TickBuf1;
  }
  else
  {
    CurrTicks = &TickBuf1;
    LastTicks = &TickBuf2;
  }

  if (iEPL < 0 && fpcp)
    Fclose(fpcp), fpcp = NULL;
}


#define LEN_FNCSVE  255


void AdminCoseqSave(char *fncsve, int ID)
{
  int   i, j;
  char  *cp;

  static int   LastID = 0;
  static int   n_back_fncsve = 0;
  static char  **back_fncsve = NULL;


  if (F_CoseqKeep < 1) return;

  if (fncsve == NULL)
  {
    if (back_fncsve)
    {
      AppFree(back_fncsve[0], F_CoseqKeep * (LEN_FNCSVE + 1));
      AppFree(back_fncsve,    F_CoseqKeep);

      LastID = 0;
      n_back_fncsve = 0;
        back_fncsve = NULL;
    }

    return;
  }

  if (back_fncsve == NULL)
  {
    CheckMalloc(back_fncsve,    F_CoseqKeep);
    CheckMalloc(back_fncsve[0], F_CoseqKeep * (LEN_FNCSVE + 1));

    for (i = 1; i < F_CoseqKeep; i++)
      back_fncsve[i] = back_fncsve[i - 1] + (LEN_FNCSVE + 1);

    n_back_fncsve = 0; LastID = ID;
  }
  else if (ID != LastID)
    n_back_fncsve = 0; LastID = ID;

  if (n_back_fncsve < F_CoseqKeep)
    i = n_back_fncsve++;
  else
  {
    cp = back_fncsve[0];

    for (i = 0, j = 1; j < n_back_fncsve; i = j++)
      back_fncsve[i] = back_fncsve[j];

    back_fncsve[i] = cp;

    (void) unlink(back_fncsve[i]);
  }

  (void) strcpy(back_fncsve[i], fncsve);
}


void CoseqSave(const T_ListSymNodes *LSymN, int FxFyFz,
               int iEPL, int SymN_Offset, int iSphere,
               int *SphereStorage, int nSphereStorage, int iAsyN)
{
  int   i, iSN, iPB, nseq;
  char   fncsve[LEN_FNCSVE + 1];
  FILE  *fpcsve;

  int  OffsCurr, OffsNew;
  int  Min_nPUC, Max_nPUC, Sum_nPUC;

  T_iVector   UC;
  T_SymNodes  *SN;
  T_PUC_Buf   *PBi;
  T_PUC       *PUC;


  Sprintf(fncsve, "%s_%3.3d_%3.3d", F_CoseqSaveFileName, iEPL, iSphere);

  if ((int) strlen(fncsve) > LEN_FNCSVE)
    InternalError("LEN_FNCSVE too small");

      fpcsve = fopen(fncsve, "w");
  if (fpcsve == NULL)
  {
    Fprintf(stdout, "\n\n# WARNING: Cannot write to %s\n\n", fncsve);
    CountWarnings++;

    fpcsve = stdout;
  }

  Fprintf(fpcsve, "# MaxMemoryUsed = %ld\n", (long) MaxMemoryUsed);

  OffsCurr = 0;
  OffsNew  = 0;

  Min_nPUC = 0;
  Max_nPUC = 0;
  Sum_nPUC = 0;

  SN  = LSymN->SymN;
  PBi = SN->PB;

  for (iSN = 0; iSN < LSymN->nSymN; iSN++, SN++)
  {
    for (iPB = 0; iPB < FxFyFz; iPB++, PBi++)
    {
      OffsCurr += PBi->oCurr;
      OffsNew  += PBi->oNew;

      if (iSN || iPB)
      {
        if (Min_nPUC > PBi->oEnd)
            Min_nPUC = PBi->oEnd;
        if (Max_nPUC < PBi->oEnd)
            Max_nPUC = PBi->oEnd;
      }
      else
      {
        Min_nPUC =
        Max_nPUC = PBi->oEnd;
      }

      Sum_nPUC += PBi->oEnd;
    }
  }

  Fprintf(fpcsve, "# nOperators * FxFyFz = %d * %d = %d\n",
    LSymN->nSymN, FxFyFz, LSymN->nSymN * FxFyFz);

  Fprintf(fpcsve, "#  Min_nPUC = %d\n", Min_nPUC);
  Fprintf(fpcsve, "#  Max_nPUC = %d\n", Max_nPUC);
  Fprintf(fpcsve, "# Mean_nPUC = %.1f\n",
    (double)(Sum_nPUC) / (LSymN->nSymN * FxFyFz));
  Fprintf(fpcsve, "#  Sum_nPUC = %d\n", Sum_nPUC);

  Fprintf(fpcsve, ">Begin coseq_save\n");

  Fprintf(fpcsve, "iEPL %d\n", iEPL);
  Fprintf(fpcsve, "OpOffset %d\n", SymN_Offset);
  Fprintf(fpcsve, "iSphere %d\n", iSphere);
  Fprintf(fpcsve, "OffsCurr %d\n", OffsCurr);
  Fprintf(fpcsve, "OffsNew %d\n", OffsNew);
  Fprintf(fpcsve, "nSphereStorage %d\n", nSphereStorage);
  Fprintf(fpcsve, "nAsyT %d\n", LSymN->nAsyN);
  Fprintf(fpcsve, "iAsyT %d\n", iAsyN);

  Fprintf(fpcsve, "SphereStorage:\n");

  nseq = iAsyN * nSphereStorage + MIN(iSphere, nSphereStorage);

  for (i = 0; i < nseq; i++)
    Fprintf(fpcsve, "%d\n", SphereStorage[i]);

  Fprintf(fpcsve, "Pool:\n");

  Fprintf(fpcsve, "# Back\n");

  SN  = LSymN->SymN;
  PBi = SN->PB;

  for (iSN = 0; iSN < LSymN->nSymN; iSN++, SN++)
  {
    for (iPB = 0; iPB < FxFyFz; iPB++, PBi++)
    {
        i = 0;
      PUC = PBi->PUC;

      for (; i < PBi->oCurr; i++, PUC++)
      {
        UnpackUC(*PUC, &UC);
        if (fprintf(fpcsve, "%d %d %d %d\n", iSN, UC.x, UC.y, UC.z) < 0)
          goto GiveUp;
      }
    }
  }

  Fprintf(fpcsve, "# Current\n");

  SN  = LSymN->SymN;
  PBi = SN->PB;

  for (iSN = 0; iSN < LSymN->nSymN; iSN++, SN++)
  {
    for (iPB = 0; iPB < FxFyFz; iPB++, PBi++)
    {
        i = PBi->oCurr;
      PUC = PBi->PUC + i;

      for (; i < PBi->oNew;  i++, PUC++)
      {
        UnpackUC(*PUC, &UC);
        if (fprintf(fpcsve, "%d %d %d %d\n", iSN, UC.x, UC.y, UC.z) < 0)
          goto GiveUp;
      }
    }
  }

  if (fprintf(fpcsve, ">End coseq_save\n") < 0)
    goto GiveUp;

  if (fpcsve != stdout)
  {
            Fclose(fpcsve);
    AdminCoseqSave(fncsve, iEPL);
  }

  return;

  GiveUp:

  if (fpcsve != stdout)
  {
    (void) unlink(fncsve);

    Fprintf(stdout, "\n\n# Error while writing coseq protocol file %s\n\n",
                    fncsve);
  }
}


static int csld_getvar(char *str, char *key, int *var, int min_var)
{
  int   n;
  char  xtrac;


  if (str_ibegin(str, key) == 0)
  {
    n = sscanf(str + strlen(key), "%d%c", var, &xtrac);

    if (   (n == 1 || n == 2 && (isspace(xtrac) || xtrac == '#'))
        && *var >= min_var)
      return 1;
  }

  return 0;
}


int CoseqReload(const T_ListSymNodes *LSymN, int FxFyFz,
                int *iEPL, int *SymN_Offset, int *iSphere,
                int *SphereStorage, int nSphereStorage, int *iAsyN)
{
  int   i, n;
  int   mode, nseq, iseq, jseq, nSS;
  long  lcount;
  char  buf[256], xtrac;
  FILE  *fpcsld;

  int         OffsCurr, OffsNew;
  int         iSN, iPB;
  T_PUC_Buf   *PBi;
  T_PUC       PUC;
  T_iVector   UC;

  int   F_Pool, F_sPool, F_AllocIncrement;
  int   F_iEPL;
  int   F_SymN_Offset, F_iSphere;
  int   F_OffsCurr, F_OffsNew;
  int   F_SphereStorage, F_nSphereStorage, F_nAsyN, F_iAsyN;


  F_Pool = F_sPool = F_AllocIncrement =
  F_iEPL =
  F_SymN_Offset = F_iSphere =
  F_OffsCurr = F_OffsNew =
  F_SphereStorage = F_nSphereStorage = F_nAsyN = F_iAsyN = 0;

      fpcsld = fopen(F_CoseqReloadFileName, "r");
  if (fpcsld == NULL)
  {
    Fprintf(stderr, "%s: Can't open %s\n", progn, F_CoseqReloadFileName);
    AppExit(1);
  }

  lcount = 0;
  mode = 0;

  while (fgetline(fpcsld, buf))
  {
    lcount++;

    if (noblks(buf) > 0 && buf[0] != '#' && buf[0] != '$')
    {
      if (mode == 0)
      {
        if (str_icmp(buf, ">Begin coseq_save") == 0)
          mode = 1;
      }
      else if (str_icmp(buf, ">End coseq_save") == 0)
        break;
      else if (mode == 1)
      {
        if      (   F_sPool == 0
                 && csld_getvar(buf, "sPool ", &i, 1))
                    F_sPool =  1;
        else if (   F_AllocIncrement == 0
                 && csld_getvar(buf, "AllocIncrement ", &i, 1))
                    F_AllocIncrement =  1;
        else if (   F_iEPL == 0
                 && csld_getvar(buf, "iEPL ", iEPL, 0))
                    F_iEPL =  1;
        else if (   F_SymN_Offset == 0
                 && csld_getvar(buf, "OpOffset ", SymN_Offset, 0))
                    F_SymN_Offset =  1;
        else if (   F_iSphere == 0
                 && csld_getvar(buf, "iSphere ", iSphere, 1))
                    F_iSphere =  1;
        else if (   F_OffsCurr == 0
                 && (   csld_getvar(buf, "OffsCurrent ", &OffsCurr, 0)
                     || csld_getvar(buf, "OffsCurr ",    &OffsCurr, 0)))
                    F_OffsCurr =  1;
        else if (   F_OffsNew == 0
                 && csld_getvar(buf, "OffsNew ", &OffsNew, 0))
                    F_OffsNew =  1;
        else if (   F_nSphereStorage == 0
                 && csld_getvar(buf, "nSphereStorage ", &nSS, 0))
        {
                    F_nSphereStorage =  1;
          if (nSS < nSphereStorage)
            IllegalLine(F_CoseqReloadFileName, lcount);
        }
        else if (   F_nAsyN == 0
                 && csld_getvar(buf, "nAsyT ", &i, 1))
        {
                    F_nAsyN =  1;
          if (i != LSymN->nAsyN)
            IllegalLine(F_CoseqReloadFileName, lcount);
        }
        else if (   F_iAsyN == 0
                 && csld_getvar(buf, "iAsyT ", iAsyN, 0))
                    F_iAsyN =  1;
        else if (   F_Pool == 0
                 && str_icmp(buf, "Pool:") == 0)
        {
                    F_Pool =  1;

          if (F_OffsCurr == 0 || F_OffsNew == 0)
            IllegalLine(F_CoseqReloadFileName, lcount);

                  /* ARBITRARY */
              n = (int)(4. * (OffsNew - OffsCurr)
                  / (LSymN->nSymN * FxFyFz) + .5);
          if (n < 10)
              n = 10; /* ARBITRARY */

          PBi = LSymN->SymN->PB;

          for (iPB = 0; iPB < LSymN->nSymN * FxFyFz; iPB++, PBi++)
          {
                                  PBi->sPUC = n;
            CheckMalloc(PBi->PUC, PBi->sPUC);
          }

          mode = 2;
          iseq = 0;
        }
        else if (   F_SphereStorage == 0
                 && str_icmp(buf, "SphereStorage:") == 0)
        {
                    F_SphereStorage  = 1;

          if (F_iAsyN == 0 || F_iSphere == 0 || F_nSphereStorage == 0)
            IllegalLine(F_CoseqReloadFileName, lcount);

              nseq = (*iAsyN) * nSS + MIN(*iSphere, nSS);
          if (nseq)
          {
            mode = 3;
            iseq = 0;
            jseq = 0;
          }
        }
        else
          IllegalLine(F_CoseqReloadFileName, lcount);
      }
      else if (mode == 2) /* Pool: */
      {
        n = sscanf(buf, "%d%d%d%d%c",
                   &iSN, &UC.x, &UC.y, &UC.z, &xtrac);

        if (    n != 4
            && (n != 5 || (isspace(xtrac) == 0 && xtrac != '#'))
            || iSN < 0 || iSN > LSymN->nSymN)
          IllegalLine(F_CoseqReloadFileName, lcount);

        if (CheckUC(&UC))
          progerror("CoseqReload(): |UC| overflow");

               PackUC(PUC, &UC);
        Calc_iPUC_Buf(iPB, &UC);

        PBi = &(LSymN->SymN[iSN].PB[iPB]);

        if (PBi->oEnd == PBi->sPUC)
        {
              n = (int)(PBi->sPUC * 1.3 + .5); /* ARBITRARY */
          if (n - PBi->sPUC < 10) /* ARBITRARY */
              n = PBi->sPUC + 10; /* ARBITRARY */

          CheckRealloc(PBi->PUC, PBi->PUC, n, PBi->sPUC);

          PBi->sPUC = n;
        }

        PBi->PUC[PBi->oEnd] = PUC;
                 PBi->oEnd++;

        if (iseq < OffsCurr)
          PBi->oCurr++;

        PBi->oNew++;

        iseq++;

        if (iseq == OffsNew)
          mode = 1;
      }
      else if (mode == 3) /* SphereStorage: */
      {
               n = sscanf(buf, "%d%c", &i, &xtrac);
        if (   n != 1 && (n != 2 || (isspace(xtrac) == 0 && xtrac != '#'))
            || i < 0)
          IllegalLine(F_CoseqReloadFileName, lcount);

        if (iseq++ % nSS < nSphereStorage)
          SphereStorage[jseq++] = i;

        if (iseq == nseq)
          mode = 1;
      }
      else
        InternalError("Corrupt mode");

      if (   F_OffsCurr && F_OffsNew && OffsCurr > OffsNew
          || F_iAsyN && *iAsyN >= LSymN->nAsyN)
        IllegalLine(F_CoseqReloadFileName, lcount);
    }
  }

  if (! (   F_iEPL
         && F_SymN_Offset && F_iSphere
         && F_OffsCurr && F_OffsNew
         && F_SphereStorage && F_nSphereStorage && F_nAsyN && F_iAsyN)
      || mode != 1)
    IllegalLine(F_CoseqReloadFileName, lcount);

  Fclose(fpcsld);

  return OffsNew - OffsCurr;
}
