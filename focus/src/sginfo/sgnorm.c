/*
  Space Group Info (c) 1994-97 Ralf W. Grosse-Kunstleve
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "sginfo.h"


#define Fprintf (void) fprintf

#define MemCpy(t, s, n)  memcpy((t), (s), (n) * sizeof (*(t)))
#define MemCmp(t, s, n)  memcmp((t), (s), (n) * sizeof (*(t)))


typedef struct
  {                   /* Ref.: Int. Tab. Vol. A Section 15.3 */
    const char  *K2L; /* operations which generate L from K  */
    const char  *L2N; /* operations which generate N from L  */
  }
  T_EuNoGen;


static const T_EuNoGen EuNoGen[] =
  {
              { NULL,   NULL },
    /*   1 */ { "-1",   NULL },
    /*   2 */ { NULL,   NULL },
    /*   3 */ { "-1",   NULL },
    /*   4 */ { "-1",   NULL },
    /*   5 */ { "-1",   NULL },
    /*   6 */ { "-1",   NULL },
    /*   7 */ { "-1",   NULL },
    /*   8 */ { "-1",   NULL },
    /*   9 */ { "-1",   NULL },
    /*  10 */ { NULL,   NULL },
    /*  11 */ { NULL,   NULL },
    /*  12 */ { NULL,   NULL },
    /*  13 */ { NULL,   NULL },
    /*  14 */ { NULL,   NULL },
    /*  15 */ { NULL,   NULL },
    /*  16 */ { "-1",   NULL },
    /*  17 */ { "-1",   NULL },
    /*  18 */ { "-1",   NULL },
    /*  19 */ { "-1",   NULL },
    /*  20 */ { "-1",   NULL },
    /*  21 */ { "-1",   NULL },
    /*  22 */ { "-1",   NULL },
    /*  23 */ { "-1",   NULL },
    /*  24 */ { "-1",   NULL },
    /*  25 */ { "-1",   NULL },
    /*  26 */ { "-1",   NULL },
    /*  27 */ { "-1",   NULL },
    /*  28 */ { "-1",   NULL },
    /*  29 */ { "-1",   NULL },
    /*  30 */ { "-1",   NULL },
    /*  31 */ { "-1",   NULL },
    /*  32 */ { "-1",   NULL },
    /*  33 */ { "-1",   NULL },
    /*  34 */ { "-1",   NULL },
    /*  35 */ { "-1",   NULL },
    /*  36 */ { "-1",   NULL },
    /*  37 */ { "-1",   NULL },
    /*  38 */ { "-1",   NULL },
    /*  39 */ { "-1",   NULL },
    /*  40 */ { "-1",   NULL },
    /*  41 */ { "-1",   NULL },
    /*  42 */ { "-1",   NULL },
    /*  43 */ { "-1uv", NULL },
    /*  44 */ { "-1",   NULL },
    /*  45 */ { "-1",   NULL },
    /*  46 */ { "-1",   NULL },
    /*  47 */ { NULL,   NULL },
    /*  48 */ { NULL,   NULL },
    /*  49 */ { NULL,   NULL },
    /*  50 */ { NULL,   NULL },
    /*  51 */ { NULL,   NULL },
    /*  52 */ { NULL,   NULL },
    /*  53 */ { NULL,   NULL },
    /*  54 */ { NULL,   NULL },
    /*  55 */ { NULL,   NULL },
    /*  56 */ { NULL,   NULL },
    /*  57 */ { NULL,   NULL },
    /*  58 */ { NULL,   NULL },
    /*  59 */ { NULL,   NULL },
    /*  60 */ { NULL,   NULL },
    /*  61 */ { NULL,   NULL },
    /*  62 */ { NULL,   NULL },
    /*  63 */ { NULL,   NULL },
    /*  64 */ { NULL,   NULL },
    /*  65 */ { NULL,   NULL },
    /*  66 */ { NULL,   NULL },
    /*  67 */ { NULL,   NULL },
    /*  68 */ { NULL,   NULL },
    /*  69 */ { NULL,   NULL },
    /*  70 */ { NULL,   NULL },
    /*  71 */ { NULL,   NULL },
    /*  72 */ { NULL,   NULL },
    /*  73 */ { NULL,   NULL },
    /*  74 */ { NULL,   NULL },
    /*  75 */ { "-1",   "-2'" },
    /*  76 */ { NULL,   "2\"" },
    /*  77 */ { "-1",   "-2'" },
    /*  78 */ { NULL,   "2\"" },
    /*  79 */ { "-1",   "-2'" },
    /*  80 */ { "-1a",  "2\"" },
    /*  81 */ { "-1",   "-2'" },
    /*  82 */ { "-1",   "-2'" },
    /*  83 */ { NULL,   "-2'" },
    /*  84 */ { NULL,   "-2'" },
    /*  85 */ { NULL,   "-2'" },
    /*  86 */ { NULL,   "-2'" },
    /*  87 */ { NULL,   "-2'" },
    /*  88 */ { NULL,   "2\"" },
    /*  89 */ { "-1",   NULL },
    /*  90 */ { "-1",   NULL },
    /*  91 */ { NULL,   NULL },
    /*  92 */ { NULL,   NULL },
    /*  93 */ { "-1",   NULL },
    /*  94 */ { "-1",   NULL },
    /*  95 */ { NULL,   NULL },
    /*  96 */ { NULL,   NULL },
    /*  97 */ { "-1",   NULL },
    /*  98 */ { "-1aw", NULL },
    /*  99 */ { "-1",   NULL },
    /* 100 */ { "-1",   NULL },
    /* 101 */ { "-1",   NULL },
    /* 102 */ { "-1",   NULL },
    /* 103 */ { "-1",   NULL },
    /* 104 */ { "-1",   NULL },
    /* 105 */ { "-1",   NULL },
    /* 106 */ { "-1",   NULL },
    /* 107 */ { "-1",   NULL },
    /* 108 */ { "-1",   NULL },
    /* 109 */ { "-1a",  NULL },
    /* 110 */ { "-1a",  NULL },
    /* 111 */ { "-1",   NULL },
    /* 112 */ { "-1",   NULL },
    /* 113 */ { "-1",   NULL },
    /* 114 */ { "-1",   NULL },
    /* 115 */ { "-1",   NULL },
    /* 116 */ { "-1",   NULL },
    /* 117 */ { "-1",   NULL },
    /* 118 */ { "-1",   NULL },
    /* 119 */ { "-1",   NULL },
    /* 120 */ { "-1",   NULL },
    /* 121 */ { "-1",   NULL },
    /* 122 */ { "-1aw", NULL },
    /* 123 */ { NULL,   NULL },
    /* 124 */ { NULL,   NULL },
    /* 125 */ { NULL,   NULL },
    /* 126 */ { NULL,   NULL },
    /* 127 */ { NULL,   NULL },
    /* 128 */ { NULL,   NULL },
    /* 129 */ { NULL,   NULL },
    /* 130 */ { NULL,   NULL },
    /* 131 */ { NULL,   NULL },
    /* 132 */ { NULL,   NULL },
    /* 133 */ { NULL,   NULL },
    /* 134 */ { NULL,   NULL },
    /* 135 */ { NULL,   NULL },
    /* 136 */ { NULL,   NULL },
    /* 137 */ { NULL,   NULL },
    /* 138 */ { NULL,   NULL },
    /* 139 */ { NULL,   NULL },
    /* 140 */ { NULL,   NULL },
    /* 141 */ { NULL,   NULL },
    /* 142 */ { NULL,   NULL },
    /* 143 */ { "-1",   "2 -2'" },
    /* 144 */ { NULL,   "2 2\"" },
    /* 145 */ { NULL,   "2 2\"" },
    /* 146 */ { "-1",   "-2\"" },
    /* 147 */ { NULL,   "2 -2'" },
    /* 148 */ { NULL,   "-2\"" },
    /* 149 */ { "-1",   "2" },
    /* 150 */ { "-1",   "2" },
    /* 151 */ { NULL,   "2" },
    /* 152 */ { NULL,   "2" },
    /* 153 */ { NULL,   "2" },
    /* 154 */ { NULL,   "2" },
    /* 155 */ { "-1",   NULL },
    /* 156 */ { "-1",   "2" },
    /* 157 */ { "-1",   "2" },
    /* 158 */ { "-1",   "2" },
    /* 159 */ { "-1",   "2" },
    /* 160 */ { "-1",   NULL },
    /* 161 */ { "-1",   NULL },
    /* 162 */ { NULL,   "2" },
    /* 163 */ { NULL,   "2" },
    /* 164 */ { NULL,   "2" },
    /* 165 */ { NULL,   "2" },
    /* 166 */ { NULL,   NULL },
    /* 167 */ { NULL,   NULL },
    /* 168 */ { "-1",   "-2'" },
    /* 169 */ { NULL,   "2\"" },
    /* 170 */ { NULL,   "2\"" },
    /* 171 */ { NULL,   "2\"" },
    /* 172 */ { NULL,   "2\"" },
    /* 173 */ { "-1",   "-2'" },
    /* 174 */ { "-1",   "-2'" },
    /* 175 */ { NULL,   "-2'" },
    /* 176 */ { NULL,   "-2'" },
    /* 177 */ { "-1",   NULL },
    /* 178 */ { NULL,   NULL },
    /* 179 */ { NULL,   NULL },
    /* 180 */ { NULL,   NULL },
    /* 181 */ { NULL,   NULL },
    /* 182 */ { "-1",   NULL },
    /* 183 */ { "-1",   NULL },
    /* 184 */ { "-1",   NULL },
    /* 185 */ { "-1",   NULL },
    /* 186 */ { "-1",   NULL },
    /* 187 */ { "-1",   NULL },
    /* 188 */ { "-1",   NULL },
    /* 189 */ { "-1",   NULL },
    /* 190 */ { "-1",   NULL },
    /* 191 */ { NULL,   NULL },
    /* 192 */ { NULL,   NULL },
    /* 193 */ { NULL,   NULL },
    /* 194 */ { NULL,   NULL },
    /* 195 */ { "-1",   "-2'" },
    /* 196 */ { "-1",   "-2'" },
    /* 197 */ { "-1",   "-2'" },
    /* 198 */ { "-1",   "-2'd" },
    /* 199 */ { "-1",   "-2'd" },
    /* 200 */ { NULL,   "-2'" },
    /* 201 */ { NULL,   "-2'" },
    /* 202 */ { NULL,   "-2'" },
    /* 203 */ { NULL,   "-2'" },
    /* 204 */ { NULL,   "-2'" },
    /* 205 */ { NULL,   "-2'd" },
    /* 206 */ { NULL,   NULL },
    /* 207 */ { "-1",   NULL },
    /* 208 */ { "-1",   NULL },
    /* 209 */ { "-1",   NULL },
    /* 210 */ { "-1d",  NULL },
    /* 211 */ { "-1",   NULL },
    /* 212 */ { NULL,   NULL },
    /* 213 */ { NULL,   NULL },
    /* 214 */ { "-1",   NULL },
    /* 215 */ { "-1",   NULL },
    /* 216 */ { "-1",   NULL },
    /* 217 */ { "-1",   NULL },
    /* 218 */ { "-1",   NULL },
    /* 219 */ { "-1",   NULL },
    /* 220 */ { "-1",   NULL },
    /* 221 */ { NULL,   NULL },
    /* 222 */ { NULL,   NULL },
    /* 223 */ { NULL,   NULL },
    /* 224 */ { NULL,   NULL },
    /* 225 */ { NULL,   NULL },
    /* 226 */ { NULL,   NULL },
    /* 227 */ { NULL,   NULL },
    /* 228 */ { NULL,   NULL },
    /* 229 */ { NULL,   NULL },
    /* 230 */ { NULL,   NULL }
  };


typedef struct
  {
    const char *Token;
    const char *xyz;
  }
  T_ListHallTokens;

static const T_ListHallTokens ListHallTokens[] =
  {
    { "-1",   "-x,-y,-z" },
    { "-1a",  "-x+1/2,-y,-z" },
    { "-1aw", "-x+1/2,-y,-z+1/4" },
    { "-1d",  "-x+1/4,-y+1/4,-z+1/4" },
    { "-1uv", "-x+1/4,-y+1/4,-z" },
    { "-2\"", "-y,-x,z" },
    { "-2'",  "y,x,z" },
    { "-2'd", "y+1/4,x+1/4,z+1/4" },
    { "2",    "-x,-y,z" },
    { "2\"",  "y,x,-z" },
    { NULL, NULL }
  };


typedef struct
  {
    int  v[3];
  }
  T_LTr;


#define SgOps_mListLTr (144)
#define SgOps_mListSMx  (24)

typedef struct
  {
    int     nListLTr;
    T_LTr    ListLTr[SgOps_mListLTr];
    int     K[3];
    int     nListSMx;
    T_RTMx   ListSMx[SgOps_mListSMx];
    int     nLSL;
    int     nSSL;
  }
  T_SgOps;


static void ResetSgOps(T_SgOps *SgOps)
{
  int  i;

  SgOps->nListLTr = 1;
  for (i = 0; i < 3; i++) SgOps->ListLTr[0].v[i] = 0;

  for (i = 0; i < 3; i++) SgOps->K[i] = -1;

  SgOps->nListSMx = 1;
  InitSeitzMx(&SgOps->ListSMx[0], 1);

  SgOps->nLSL = 1;
  SgOps->nSSL = 1;
}


static int Add2ListLTr(T_SgOps *SgOps, const int *NewLTr)
{
  int    NLTr[3], i;
  int    iLLTr;
  T_LTr  *LLTr;


  for (i = 0; i < 3; i++)
    NLTr[i] = iModPositive(NewLTr[i], STBF);

  LLTr = SgOps->ListLTr;

  for (iLLTr = 0; iLLTr < SgOps->nListLTr; iLLTr++, LLTr++)
    if (MemCmp(LLTr->v, NLTr, 3) == 0)
      return 0; /* lattice translation already in list */

  if (SgOps->nListLTr >= SgOps_mListLTr) {
    SetSgError("Internal Error: SgOps_mListLTr too small");
    return -1;
  }

  (void) MemCpy(LLTr->v, NLTr, 3);

  SgOps->nListLTr++;

  return 1;
}


static int AddLtrDueToK(T_SgOps *SgOps, const T_RTMx *LSMx)
{
  int  NewLTr[3], i;


  RotMx_t_Vector(NewLTr, LSMx->s.R, SgOps->K, 0);
  for (i = 0; i < 3; i++) NewLTr[i] += 2 * LSMx->s.T[i] - SgOps->K[i];
  return Add2ListLTr(SgOps, NewLTr);
}


static int Add2CentInv(T_SgOps *SgOps, const int *K)
{
  int     NewLTr[3], i;
  int     iLSMx;
  T_RTMx  *LSMx;

  const int  NNN[] = { 0, 0, 0};


  if (! K) K = NNN;

  if (SgOps->K[0] >= 0) /* there is a centre of inversion already */
  {
    for (i = 0; i < 3; i++) NewLTr[i] = SgOps->K[i] - K[i];
    return Add2ListLTr(SgOps, NewLTr);
  }

  for (i = 0; i < 3; i++) SgOps->K[i] = iModPositive(K[i], STBF);

  LSMx = &SgOps->ListSMx[1];

  for (iLSMx = 1; iLSMx < SgOps->nListSMx; iLSMx++, LSMx++)
    if (AddLtrDueToK(SgOps, LSMx) < 0)
      return -1;

  return 1;
}


static int Add2ListSMx(T_SgOps *SgOps, const T_RTMx *NewSMx)
{
  int     mR[9], NewLTr[3], K[3], i;
  int     iLSMx;
  T_RTMx  *LSMx;


  for (i = 0; i < 9; i++) mR[i] = -NewSMx->s.R[i];

  LSMx = SgOps->ListSMx;

  for (iLSMx = 0; iLSMx < SgOps->nListSMx; iLSMx++, LSMx++)
  {
    if (MemCmp(LSMx->s.R, NewSMx->s.R, 9) == 0) {
      for (i = 0; i < 3; i++)   NewLTr[i] = LSMx->s.T[i] - NewSMx->s.T[i];
      return Add2ListLTr(SgOps, NewLTr);
    }

    if (MemCmp(LSMx->s.R, mR, 9) == 0) {
      for (i = 0; i < 3; i++)   K[i]      = LSMx->s.T[i] + NewSMx->s.T[i];
      return Add2CentInv(SgOps, K);
    }
  }

  if (SetRotMxInfo(NewSMx->s.R, NULL) == 0) {
    SetSgError("Error: Non-crystallographic rotation matrix encountered");
    return -1;
  }

  if (SgOps->nListSMx >= SgOps_mListSMx) {
    SetSgError("Internal Error: SgOps_mListSMx too small");
    return -1;
  }

  (void) MemCpy(LSMx->s.R, NewSMx->s.R, 9);
  for (i = 0; i < 3; i++) LSMx->s.T[i] = iModPositive(NewSMx->s.T[i], STBF);

  SgOps->nListSMx++;

  if (SgOps->K[0] >= 0 && AddLtrDueToK(SgOps, LSMx) < 0)
    return -1;

  return 1;
}


static int DoMulSMxLTr(T_SgOps *SgOps, int iLSMx, int iLLTr, int OldOnly)
{
  const T_RTMx  *LSMxi;
  const T_LTr   *LLTri;
  int           NewLTr[3];


  LSMxi = &SgOps->ListSMx[iLSMx];

  for (; iLSMx < SgOps->nListSMx; iLSMx++, LSMxi++)
  {
    LLTri = &SgOps->ListLTr[iLLTr];

    for (; iLLTr < (OldOnly ? SgOps->nLSL : SgOps->nListLTr);
           iLLTr++, LLTri++)
    {
      RotMx_t_Vector(NewLTr, LSMxi->s.R, LLTri->v, 0);
      if (Add2ListLTr(SgOps, NewLTr) < 0)
        return -1;
    }
  }

  return 0;
}


static int XpndListLTr(T_SgOps *SgOps, const int *NewLTr)
{
  int           TrialLTr[3], i;
  int           iLLTr,  jLLTr;
  const T_LTr   *LLTri, *LLTrj;


  if (DoMulSMxLTr(SgOps, SgOps->nSSL, 1, 1) < 0) return -1;
  SgOps->nSSL = SgOps->nListSMx;

  iLLTr  =  SgOps->nLSL;
   LLTri = &SgOps->ListLTr[iLLTr];

  jLLTr  = 1;
   LLTrj = &SgOps->ListLTr[jLLTr];

  for (;;)
  {
    if (NewLTr) {
      if (Add2ListLTr(SgOps, NewLTr) < 0)
        return -1;
    }

    if (DoMulSMxLTr(SgOps, 1, SgOps->nLSL, 0) < 0) return -1;
    SgOps->nLSL = SgOps->nListLTr;

    if (jLLTr > iLLTr)
    {
      iLLTr++;
       LLTri++;

      jLLTr  = 1;
       LLTrj = &SgOps->ListLTr[jLLTr];
    }

    if (iLLTr == SgOps->nListLTr)
      break;

    for (i = 0; i < 3; i++)
      TrialLTr[i] = LLTrj->v[i] + LLTri->v[i];

    NewLTr = TrialLTr;

    jLLTr++;
     LLTrj++;
  }

  return 0;
}


static int XpndCentInv(T_SgOps *SgOps, const int *K)
{
  if (Add2CentInv(SgOps, K) < 0)
    return -1;

  return XpndListLTr(SgOps, NULL);
}


static int XpndListSMx(T_SgOps *SgOps, const T_RTMx *NewSMx)
{
  int           iLSMx,  jLSMx;
  const T_RTMx  *LSMxi, *LSMxj;
  T_RTMx        TrialSMx[1];


  iLSMx  =  SgOps->nListSMx;
   LSMxi = &SgOps->ListSMx[iLSMx];

  jLSMx  = 1;
   LSMxj = &SgOps->ListSMx[jLSMx];

  for (;;)
  {
    if (NewSMx && Add2ListSMx(SgOps, NewSMx) < 0)
      return -1;

    if (jLSMx > iLSMx)
    {
      iLSMx++;
       LSMxi++;

      jLSMx  = 1;
       LSMxj = &SgOps->ListSMx[jLSMx];
    }

    if (iLSMx == SgOps->nListSMx)
      break;

    SeitzMxMultiply(TrialSMx, LSMxj, LSMxi);

    NewSMx = TrialSMx;

    jLSMx++;
     LSMxj++;
  }

  return XpndListLTr(SgOps, NULL);
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


int GetEuNoGen(int SgNumber, const char *HallSymbol, int Which, T_RTMx *AddlG)
{
  const T_TabSgName       *TSgN;
  int                     nAddlG, iType;
  const char              *GHT;
  const T_ListHallTokens  *LHT;
  char                    buf[16], *cp;


  if (Which < 1 || Which > 3) {
    SetSgError(
      "Internal Error: GetEuNoGen(): Illegal value for parameter \"Which\"");
    return -1;
  }

  if (SgNumber < 1 || SgNumber > 230) return -1;

  if (HallSymbol != NULL)
  {
        TSgN = FindTabSgNameEntry(SgNumber, NULL, 'A');
    if (TSgN == NULL) {
      SetSgError("Internal Error: Corrupt FindTabSgNameEntry()");
      return -1;
    }

    if (strcmp(HallSymbol, TSgN->HallSymbol) != 0)
      return -1;
  }

  nAddlG = 0;

  for (iType = 0; iType < 2; iType++)
  {
    GHT = NULL;

    if      (iType == 0 && (Which == 1 || Which == 3))
      GHT = EuNoGen[SgNumber].K2L;
    else if (iType == 1 && (Which == 2 || Which == 3))
      GHT = EuNoGen[SgNumber].L2N;

    while (GHT)
    {
      cp = buf;
      while (*GHT && *GHT != ' ') *cp++ = *GHT++;
      *cp = '\0';

      for (LHT = ListHallTokens; LHT->Token; LHT++)
        if (strcmp(buf, LHT->Token) == 0) break;

      if (   ! LHT
          || nAddlG == 3
          || ParseStrXYZ(LHT->xyz, &AddlG[nAddlG], 1, STBF) != 0) {
        SetSgError("Internal Error: Corrupt ListHallTokens");
        return -1;
      }

      nAddlG++;

      if (! *GHT) break;

      GHT++;
    }
  }

  return nAddlG;
}


static int TestNewGen(T_SgInfo *SgInfo)
{
  int       iList, jList, OrderP;
  T_SgInfo  TstSgInfo[1];
  T_SgOps   SgOps[1];

  int     iSMx, iInv, iLTr, i;
  T_RTMx  SMx[1], TrialSMx[1];


  if (ExpandListSeitzMx(SgInfo) != 0)
    return -1;

  if (AllocSgInfo(TstSgInfo) != 0) return -1;

  for (iList = 1; iList < SgInfo->OrderL - 1; iList++) {
    for (jList = iList + 1; jList < SgInfo->OrderL; jList++)
    {
      Fprintf(stdout, " %2d %2d\n", iList, jList);

      ResetSgInfo(TstSgInfo);

      if (Add2ListSeitzMx(TstSgInfo, &SgInfo->ListSeitzMx[iList]) != 0)
        goto CleanExit;
      if (Add2ListSeitzMx(TstSgInfo, &SgInfo->ListSeitzMx[jList]) != 0)
        goto CleanExit;

      ResetSgOps(SgOps);

      if (XpndListSMx(SgOps, &SgInfo->ListSeitzMx[iList]) < 0)
        goto CleanExit;
      if (XpndListSMx(SgOps, &SgInfo->ListSeitzMx[jList]) < 0)
        goto CleanExit;

      OrderP = SgOps->nListSMx;
      if (SgOps->K[0] >= 0) OrderP *= 2;

      if (TstSgInfo->nList != OrderP * SgOps->nListLTr) {
        SetSgError("Internal Error: OrderL mismatch");
        goto CleanExit;
      }

      for (iSMx = 0; iSMx < SgOps->nListSMx; iSMx++)
      {
        (void) MemCpy(SMx->a, SgOps->ListSMx[iSMx].a, 12);

        for (iInv = 0; iInv < (SgOps->K[0] < 0 ? 1 : 2); iInv++)
        {
          if (iInv) {
            for (i = 0; i < 12; i++) SMx->a[i] *= -1;
            for (i = 0; i < 3; i++) SMx->s.T[i] += SgOps->K[i];
          }

          (void) MemCpy(TrialSMx->s.R, SMx->s.R, 9);

          for (iLTr = 0; iLTr < SgOps->nListLTr; iLTr++)
          {
            for (i = 0; i < 3; i++)
              TrialSMx->s.T[i] = SMx->s.T[i] + SgOps->ListLTr[iLTr].v[i];

            if (Add2ListSeitzMx(TstSgInfo, &SgInfo->ListSeitzMx[jList]) != 0)
              goto CleanExit;

            if (TstSgInfo->nList != OrderP * SgOps->nListLTr) {
              SetSgError("Internal Error: TestNewGen() failure");
              goto CleanExit;
            }
          }
        }
      }
    }
  }

  FreeSgInfo(TstSgInfo);

  return 0;

  CleanExit:

  FreeSgInfo(TstSgInfo);
  return -1;
}


void TestAddlG(void)
{
  int                SgNumber;
  int                iList, iTrV, nTrV, OldK0, OrderP;
  const int          *TrV;
  T_RTMx             AddlG[3];
  T_SgInfo           SgInfo[1];
  T_SgOps            SgOps[1];
  const T_TabSgName  *TSgN;

  int  nAddlG, iAddlG, AddlLtr;
  int  NewGen, OldGen;


  NewGen = 1;
  OldGen = 1;

  if (AllocSgInfo(SgInfo) != 0) return;

  for (SgNumber = 1; SgNumber <= 230; SgNumber++)
  {
    Fprintf(stdout, "%3d\n", SgNumber);

    ResetSgInfo(SgInfo);
    SgInfo->GenOption = -1;

        SgInfo->TabSgName = FindTabSgNameEntry(SgNumber, NULL, 'A');
    if (SgInfo->TabSgName == NULL) {
      SetSgError("Internal Error: TestAddlG(): FindTabSgNameEntry()");
      return;
    }

    (void) ParseHallSymbol(SgInfo->TabSgName->HallSymbol, SgInfo);
    if (SgError != NULL) {
      SetSgError("Internal Error: TestAddlG(): ParseHallSymbol()");
      return;
    }

    if (NewGen)
    {
      ResetSgOps(SgOps);

      if (SgInfo->Centric == -1)
        if (XpndCentInv(SgOps, NULL) < 0) return;

      nTrV =  SgInfo->LatticeInfo->nTrVector;
       TrV = &SgInfo->LatticeInfo->TrVector[3];

      for (iTrV = 1; iTrV < nTrV; iTrV++, TrV += 3)
        if (XpndListLTr(SgOps, TrV) < 0) return;

      for (iList = 1; iList < SgInfo->nList; iList++)
        if (XpndListSMx(SgOps, &SgInfo->ListSeitzMx[iList]) < 0) return;
    }

    if (OldGen)
    {
      TSgN = SgInfo->TabSgName;
      ResetSgInfo(SgInfo);
      SgInfo->GenOption = 1;
      SgInfo->TabSgName = TSgN;

      (void) ParseHallSymbol(SgInfo->TabSgName->HallSymbol, SgInfo);
      if (SgError != NULL) {
        SetSgError("Internal Error: TestAddlG(): ParseHallSymbol()");
        return;
      }

      if (CompleteSgInfo(SgInfo) != 0) return;

      /* if (TestNewGen(SgInfo) < 0) return; */
    }

    if (NewGen && OldGen)
    {
      OrderP = SgOps->nListSMx;
      if (SgOps->K[0] >= 0) OrderP *= 2;

      if (SgInfo->OrderP != OrderP) {
        SetSgError("Internal Error: OrderP mismatch");
        return;
      }

      if (SgInfo->OrderL != OrderP * SgOps->nListLTr) {
        SetSgError("Internal Error: OrderL mismatch");
        return;
      }
    }

        nAddlG = GetEuNoGen(SgNumber, SgInfo->TabSgName->HallSymbol, 3, AddlG);
    if (nAddlG < 0) {
      SetSgError("Internal Error: TestAddlG(): GetEuNoGen()");
      return;
    }

    for (iAddlG = 0; iAddlG < nAddlG; iAddlG++) {
      const char
           *xyz = RTMx2XYZ(&AddlG[iAddlG], 1, STBF, 0, 0, 1, ", ", NULL, 0);
      if (! xyz) return;
      Fprintf(stdout, "%s\n", xyz);
    }

    if (NewGen)
    {
      AddlLtr = 0;

      for (iAddlG = 0; iAddlG < nAddlG; iAddlG++)
      {
        iTrV  = SgOps->nListLTr;
        iList = SgOps->nListSMx;
        OldK0 = SgOps->K[0];

        if (XpndListSMx(SgOps, &AddlG[iAddlG]) < 0) return;

        if (   iTrV  == SgOps->nListLTr
            && iList == SgOps->nListSMx
            && OldK0 == SgOps->K[0]) {
          SetSgError("Internal Error: Generator has no effect");
          return;
        }

        if (iTrV != SgOps->nListLTr)
          AddlLtr++;
      }

      if (AddlLtr)
        Fprintf(stdout, " +%d\n", AddlLtr);
    }
  }

  FreeSgInfo(SgInfo);
}
