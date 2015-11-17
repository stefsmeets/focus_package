#ifndef A_XTAL_H__
#define A_XTAL_H__


#include "atominfo.h"
#include "sginfo.h"
#include "lattice.h"
#include "cputime.h"


#ifndef Ex
#ifdef A_AppGlobal__
#define Ex
#else
#define Ex extern
#endif
#endif



Ex int ModeScatteringFactorTable;


typedef struct { Fprec r, i; } CmplxFprec;


typedef struct
  {
    int    h, k, l;
    Fprec  Fobs, SigmaFobs;
    Fprec  FWHM;
    int    Overlap;
    Fprec  Q;
    int    M;
    int    PhaseRestriction;
    int    uvw[3];
    int    Status;
    Fprec  Fmrg, SigmaFmrg;
  }
  T_FobsRaw;

/* T_FobsRaw Status */

#define FRS_Undefined    0
#define FRS_F000         1
#define FRS_SysAbsent    2
#define FRS_SymEquiv     3
#define FRS_Omit         4
#define FRS_OffGrid      5
#define FRS_Active       6
#define FRS_Sleeping     7


typedef struct
  {
    int  TH;
    int  Sign;
    int  ih, ik, il;
  }
  T_CTData;

typedef struct
  {
    T_FobsRaw   *FobsRaw;
    Fprec       M_Fmrg;
    Fprec       *FF;       /* FormFactor_hkl(iAtomType) */
    CmplxFprec  Fcal;
    Fprec       FcalAbs;
    Fprec       FcalPhi;
    int         Phi360Fmrg;
    int         Phi360Fcal;
    int         nCTData;
    T_CTData    *CTData;
  }
  T_CodeTransl;


#define Phi2Phi360(Phi_, Phi360_)\
  {\
    if ((Phi_) < (Fprec) 0.)\
    {\
      Phi360_ = (int) ((Phi_) / (Fprec) PIover180 - (Fprec) .5) % 360;\
      if ((Phi360_) != 0) Phi360_ += 360;\
    }\
    else\
      Phi360_ = (int) ((Phi_) / (Fprec) PIover180 + (Fprec) .5) % 360;\
  }


typedef struct
  {
    Fprec Value; int Percent, Weighted;
  }
  T_ValPerW;

typedef struct
  {
    T_ValPerW  PhaseDiff;
    T_ValPerW  DeltaR;
    int        Branch;
  }
  T_BreakIf;

#define BrkIf_Branch_PhaseDiff 0x4
#define BrkIf_Branch_Or        0x2
#define BrkIf_Branch_DeltaR    0x1

/*
  000  0  Dumb Rule
  001  1  DeltaR
  010  2  undefined (equiv. Dumb Rule)
  011  3  undefined (equiv. DeltaR)
  100  4  PhaseDiff
  101  5  PhaseDiff and DeltaR
  110  6  undefined (equiv. PhaseDiff)
  111  7  PhaseDiff or DeltaR
*/


typedef struct
  {
    Fprec  stol;
    Fprec  f;
  }
  T_SF_Tab;

typedef struct
  {
    char      *Label;
    int       mTab;
    int       nTab;
    T_SF_Tab  *Tab;
  }
  T_SF_Tables;

typedef struct
  {
    const char           *Lbl;
    T_SF_Tables          *SFT;
    const T_SF_ALL_CAA  *CAA;
    Fprec                f_stol_0;
  }
  T_SF_Info;


typedef struct
  {
    int                 On;
    int                 Class;
    T_SF_Info           SF_Info;
    int                 nPerUnitCell;
    Fprec               OccDefault;
    Fprec               UisoDefault;
    const T_AtomRadius  *eListAR;
  }
  T_AtomType;

#define ATC_Unknown    (-1)
#define ATC_General     (0)
#define ATC_Node        (1)
#define ATC_NodeBridge  (2)

#define Df_OccDefault   (1.)
#define Df_UisoDefault  (0.035)


typedef struct
  {
    int        Class;
    T_SF_Info  SF_Info;
  }
  T_AtomID;


typedef struct
  {
    int       Type;
    T_AtomID  AtomID[3];
    Fprec     Value;
  }
  T_InputChemistry;

#define ICT_MinDistance (0)
#define ICT_MaxDistance (1)
#define ICT_AddBondAtom (2)


typedef struct
  {
    char       *Label;
    T_SF_Info  SF_Info;
    Fprec      x, y, z;
    Fprec      Occ;
    int        F_Occ;
    Fprec      Uiso;
    int        F_Uiso;
  }
  T_Site;


typedef struct
  {
    Fprec  V[3];
    int    D;
  }
  T_FreeVect;


typedef struct T_WyckoffList
  {
    int                   nPositions;
    int                   *WL_Flag;
    T_FreeVect            FreeVect;
    struct T_WyckoffList  **SubWL;
    int                   CheckCount;
    int                   Count;
    int                   CanBeNode;
    int                   *CanBeCN;
  }
  T_WyckoffList;

#define WL_F_Active    (0)
#define WL_F_Identity  (1)
#define WL_F_InActive (-1)

#define MaxWyckoffList 43 /* IT VolA 1983 8.3.2: Max = 27 */


typedef struct
  {
    Fprec  x, y, z;
  }
  T_fVector;


typedef struct
  {
    int  x, y, z;
  }
  T_iVector;


typedef struct
  {
    int        Info;  /* +-(iPos * 27 + i_m0pSh + 1)      */
    T_fVector  Vect;  /* +  => asy                        */
    Fprec      Abs2;  /*  - => sym equiv to previous bond */
  }
  T_PotentialNNBond;


typedef struct T_ListPNNBonds
  {
    int                    iEPL;
    int                    nPNNB;
    T_PotentialNNBond      *PNNB;
    struct T_ListPNNBonds  *ReturnLPNNB;
  }
  T_ListPNNBonds;


typedef struct
  {
    T_fVector       Position;
    Fprec           Grid_eD;
    Fprec           Maximum;
    Fprec           Integral;
    T_WyckoffList   *WL_Entry;
    T_fVector       *xSymEquiv;
    T_fVector       *cSymEquiv;
    int             nPfI;
    int             Index;
    int             iAtomType;
    int             mLPNNB;
    int             nLPNNB;
    int             lLPNNB;
    T_ListPNNBonds  *LPNNB;
    int             Color;
    int             IsAOS; /* Is Allowed Origin Shift */
  }
  T_eD_PeakList;


#define fSTBF ((Fprec) STBF)
#define fCTBF ((Fprec) CTBF)


typedef union
  {
    struct { Fprec R[9], T[3]; } s;
    Fprec                        a[12];
  }
  T_fRTMx;


#define iRTMx_t_fVector(RTMxV_, RTMx_, V_)\
  {\
    (RTMxV_)->x =   (RTMx_)->s.R[0] * (V_)->x\
                  + (RTMx_)->s.R[1] * (V_)->y\
                  + (RTMx_)->s.R[2] * (V_)->z + (RTMx_)->s.T[0] / fSTBF;\
    (RTMxV_)->y =   (RTMx_)->s.R[3] * (V_)->x\
                  + (RTMx_)->s.R[4] * (V_)->y\
                  + (RTMx_)->s.R[5] * (V_)->z + (RTMx_)->s.T[1] / fSTBF;\
    (RTMxV_)->z =   (RTMx_)->s.R[6] * (V_)->x\
                  + (RTMx_)->s.R[7] * (V_)->y\
                  + (RTMx_)->s.R[8] * (V_)->z + (RTMx_)->s.T[2] / fSTBF;\
  }


#define fRTMx_t_fVector(RTMxV_, RTMx_, V_)\
  {\
    (RTMxV_)->x =   (RTMx_)->s.R[0] * (V_)->x\
                  + (RTMx_)->s.R[1] * (V_)->y\
                  + (RTMx_)->s.R[2] * (V_)->z + (RTMx_)->s.T[0];\
    (RTMxV_)->y =   (RTMx_)->s.R[3] * (V_)->x\
                  + (RTMx_)->s.R[4] * (V_)->y\
                  + (RTMx_)->s.R[5] * (V_)->z + (RTMx_)->s.T[1];\
    (RTMxV_)->z =   (RTMx_)->s.R[6] * (V_)->x\
                  + (RTMx_)->s.R[7] * (V_)->y\
                  + (RTMx_)->s.R[8] * (V_)->z + (RTMx_)->s.T[2];\
  }


#define NormOf(x_)\
  {\
    if ((x_) <  (Fprec) 0.)\
      do { (x_) += (Fprec) 1.; } while ((x_) <  (Fprec) 0.);\
    if ((x_) >= (Fprec) 1.)\
      do { (x_) -= (Fprec) 1.; } while ((x_) >= (Fprec) 1.);\
  }

#define CountNormOf(x_, n_)\
  {\
    (n_) = 0;\
    if ((x_) <  (Fprec) 0.)\
      do { (x_) += (Fprec) 1.; (n_) += STBF; } while ((x_) <  (Fprec) 0.);\
    if ((x_) >= (Fprec) 1.)\
      do { (x_) -= (Fprec) 1.; (n_) -= STBF; } while ((x_) >= (Fprec) 1.);\
  }

#define Tform_xc(xx, yx, zx, xc, yc, zc)\
  {\
    /* since we know TMx_xc[2], TMx_xc[3], and TMx_xc[5] always equal 0,\
       we can save some calculations\
     */\
    xc =    TMx_xc[0] * (xx) +    TMx_xc[1] * (yx) /* + TMx_xc[2] * (zx) */;\
    yc = /* TMx_xc[3] * (xx) + */ TMx_xc[4] * (yx) /* + TMx_xc[5] * (zx) */;\
    zc =    TMx_xc[6] * (xx) +    TMx_xc[7] * (yx)    + TMx_xc[8] * (zx)   ;\
  }

#define Tform_cx(xc, yc, zc, xx, yx, zx)\
  {\
    /* since we know TMx_cx[2], TMx_cx[3], and TMx_cx[5] always equal 0,\
       we can save some calculations\
     */\
    xx =    TMx_cx[0] * (xc) +    TMx_cx[1] * (yc) /* + TMx_cx[2] * (zc) */;\
    yx = /* TMx_cx[3] * (xc) + */ TMx_cx[4] * (yc) /* + TMx_cx[5] * (zc) */;\
    zx =    TMx_cx[6] * (xc) +    TMx_cx[7] * (yc)    + TMx_cx[8] * (zc)   ;\
  }

#define COS___1_DEG ( .999847695156391)
#define COS__30_DEG ( .866025403784439)
#define COS__50_DEG ( .642787609686539)
#define COS__60_DEG ( .500000000000000)
#define COS__75_DEG ( .258819045102521)
#define COS__89_DEG ( .017452406437284)
#define COS_120_DEG (-.500000000000000)
#define COS_150_DEG (-.866025403784439)
#define COS_175_DEG (-.996194698091746)

#define EvalAbsDot(AbsDot_)\
     (  (AbsDot_) > (Fprec) COS___1_DEG  /* ARBITRARY */\
        ?  1 /* parallel */\
        :\
       ((AbsDot_) < (Fprec) COS__89_DEG  /* ARBITRARY */\
        ? -1 /* perpendicular */\
        :  0 /* nothing special */\
       )\
     )\


typedef struct
  {
    struct T_SymNodes  *NextSymN;
    T_iVector          UC_Offset;
  }
  T_SymNodeBonds;


typedef unsigned int T_PUC;

typedef struct
  {
    T_PUC  *PUC;
    int    sPUC;
    int    oCurr;
    int    oNew;
    int    oEnd;
#ifdef COSEQ_TRACE
    int    *Ngen;
    int    sNgen;
#endif
  }
  T_PUC_Buf;


typedef struct
  {
    int  Loop;
    int  Count;
  }
  T_LoopCount;


typedef struct
  {
    int          nAngles;
    T_LoopCount  *LC;
  }
  T_LoopConf;


typedef struct T_SymNodes
  {
    int             nBonds;
    T_SymNodeBonds  *NB;
    T_PUC_Buf       *PB;
  }
  T_SymNodes;


typedef struct
  {
    int             Colored;
    const int       *LnB;
    int             Low_iEPL;
    int             Top_iEPL;
    int             nAsyN;
    int             nSymN;
    int             nSymNB;
    T_SymNodes      *SymN;
    T_SymNodeBonds  *SymNB;
    int             MaxLoopSize;
    T_LoopConf      *LCf;
    int             *LS_Flags;
  }
  T_ListSymNodes;


#define NodeColor(i) (eD_PeakList[i].Color)


#define PUC_BPC 10
#define PUC_Mask ((T_PUC) 0x000003FF)
#define PUC_MinC (-512)
#define PUC_MaxC   511

#define CheckUC(UC_)\
  (   (UC_)->x < PUC_MinC || (UC_)->x > PUC_MaxC\
   || (UC_)->y < PUC_MinC || (UC_)->y > PUC_MaxC\
   || (UC_)->z < PUC_MinC || (UC_)->z > PUC_MaxC)

#define PackUC(PUC, UC_)\
  PUC =   (T_PUC)((UC_)->x - PUC_MinC);\
  PUC <<= PUC_BPC;\
  PUC |=  (T_PUC)((UC_)->y - PUC_MinC);\
  PUC <<= PUC_BPC;\
  PUC |=  (T_PUC)((UC_)->z - PUC_MinC)

#define UnpackUC(PUC, UC_)\
  (UC_)->x = (int)(((PUC) >> (2 * PUC_BPC))           ) + PUC_MinC;\
  (UC_)->y = (int)(((PUC) >> (1 * PUC_BPC)) & PUC_Mask) + PUC_MinC;\
  (UC_)->z = (int)(((PUC)                 ) & PUC_Mask) + PUC_MinC

#define Calc_iPUC_Buf(iPB, UC_)\
  iPB =  abs((UC_)->x + 1) % F_CoseqSplit_x;\
  iPB *= F_CoseqSplit_y;\
  iPB += abs((UC_)->y + 1) % F_CoseqSplit_y;\
  iPB *= F_CoseqSplit_z;\
  iPB += abs((UC_)->z + 1) % F_CoseqSplit_z\


/* input values & related */

Ex Fprec  FobsMaxQ;

Ex int        nFobsRaw;
Ex T_FobsRaw  *FobsRaw;

Ex T_SgInfo   SpgrInfo;
Ex Fprec      TMx_xc[9];
Ex Fprec      TMx_cx[9];
Ex T_fVector  *fTrVector;
Ex T_fRTMx    *List_fSeitzMx;

Ex int                 nInputLatConD;
Ex T_LatticeConstants  LatConD;
Ex T_LatticeConstants  LatConR;
Ex Fprec               MinLatticeTr2;
Ex Fprec               MaxLatticeTr2;

Ex const char  *LambdaName;
Ex Fprec        LambdaLength;

Ex int  nPeaks;
Ex int  nAsyPeaks;
Ex int  PeakSearchLevel;
Ex struct { Fprec Value; int Percent; }
        eDensityCutOff;
Ex int  MinPfI;
Ex int  eD_PeaksSortElement;

#define PSE_Grid_eD  0
#define PSE_Maximum  1
#define PSE_Integral 2
#define PSE_lLPNNB   3
#define PSE_Index    4

Ex int         nAtomType;
Ex int         uAtomType;
Ex T_AtomType  *AtomType;

Ex int               nInputChemistry;
Ex T_InputChemistry  *InputChemistry;


Ex int  MinConsecutiveNodes;

#define MaxNCNmax (16) /* we want some variables to be on the stack */

typedef struct
  {
    int  CN;        /* Connectivity Number */
    int  MaxNpAU;   /* Max number of Nodes per asymmetric unit */
    int  NoSPO[13]; /* Not on Special Position with Order (i - 6) */
  }
  T_NodeType;

Ex T_NodeType  *NodeTypes;
Ex int         nNodeTypes;

Ex int  Check3DimConnectivity;

#define FwSM_FwTracking     (0)
#define FwSM_AltFwTracking  (1)

#define PNNB_TooClose  (-1)
#define PNNB_TooMany   (-2)

Ex int    IndicateFw;
Ex int    FwSearchMethod;
Ex int    MinLoopSize;
Ex int    MaxLoopSize;
Ex int    EvenLoopSizesOnly;
Ex int    NCNmin;
Ex int    NCNmax;
Ex Fprec  MinNodeDist2;
Ex Fprec  MaxNodeDist2;
Ex int    MinSymNodes;
Ex int    MaxSymNodes;
Ex Fprec  IdealT_NodeDist2;
Ex Fprec  IdealTetrahedralEdge2;
Ex Fprec  IdealTetrahedralVol2;
Ex Fprec  IdealTetrahedralVolFrac;
Ex int    ModeCheckTetrahedra;

Ex int     nSite;
Ex T_Site  *Site;

Ex int          nSF_Tables;
Ex T_SF_Tables  *SF_Tables;

Ex Fprec  CatchDistance2;
Ex Fprec  FobsScale;
Ex Fprec  SigmaCutOff;
Ex Fprec  OverlapFactor;
Ex int    OverlapAction;
Ex struct { int Valid; Fprec U, V, W; }
          GenerateFWHM;
Ex struct { Fprec Value; int Percent; }
          ReflectionUsage;
Ex struct { Fprec  Value;
            int    Type;
            Fprec  Frac_I_ignored;
            Fprec  Frac_I_w_moved;
          }
          AbsentRedivisionLimit;

#define ARLT_FWHM          0
#define ARLT_Degree2Theta  1

#define OvlAct_NoAction  0
#define OvlAct_EqualF2   1
#define OvlAct_EqualMF2  2
#define OvlAct_Omit      3


Ex Fprec  ProfileStart, ProfileEnd, ProfileStep;
Ex Fprec  ProfileGenEnd;
Ex Fprec  ProfilePOLRA;
Ex struct { int Valid; Fprec U, V, W; }
          ProfileFWHM;
Ex struct { int Valid; Fprec a1, a2, a3; }
          ProfileAsym;
Ex int    ProfilePeakShape;
Ex Fprec  PseudoVoigtPeakRange;
Ex Fprec  PseudoVoigtFracLorentz;
Ex Fprec  ProfileBackground;
Ex struct { int Mode; int h, k, l; Fprec TTheta; int Index; }
          ProfileReferenceRefl;
Ex Fprec  ProfileReferenceMax;

#define PPS_PseudoVoigt  0
#define PPS_StdPeak      1

#define PRRM_Free    0
#define PRRM_hkl     1
#define PRRM_TTheta  2

typedef struct
  {
    Fprec  R, Sym, Asy, dSym, dAsy;
  }
  T_StdPeakC;

typedef struct
  {
    char        *Title;
    int         nSteps;
    Fprec       Range;
    T_StdPeakC  *C;
  }
  T_StdPeak;


typedef struct
  {
    int   nAsyN;
    int   nSymN;
    float mBpN;
    int   Low_iEPL;
    int   Top_iEPL;
    int   *LnB;
  }
  T_LargestFwFragment;

Ex  T_LargestFwFragment  LargestFwFragment;


Ex  int        nFeedBackCycles;
Ex  int        *FeedBackCycles;
Ex  T_BreakIf   FeedBackBreakIf;

Ex Fprec  F000cal;
Ex Fprec  F000mrg;
Ex Fprec  F000mtp;

Ex int  nF000;
Ex int  nSysAbsent;
Ex int  nOmit;
Ex int  nOffGrid;
Ex int  nSymEquivFobs;
Ex int  nActivePhase;
Ex int  nSleeping;

typedef struct
  {
    Fprec  Q;
    int    h, k, l;
  }
  T_MaxQ_MaxH;

Ex T_MaxQ_MaxH  Ac_MaxQ_MaxH;
Ex T_MaxQ_MaxH  AS_MaxQ_MaxH;

Ex int           nCodeTransl;
Ex T_CodeTransl  *CodeTransl;

Ex int    Nx, Ny, Nz;
Ex Fprec  eDminRho;
Ex Fprec  eDmaxRho;
Ex Fprec  IRho3dVscale;
Ex Fprec  IRho3dV;
Ex Fprec  IRius;
Ex Fprec  tRius;

Ex Fprec  SumFmrg;
Ex Fprec  MaxFmrg;
Ex Fprec  BogFmrg;

typedef int T_PeakFlags;
                               /* Attention: assuming 8 bits/byte */
#define PF_HighBit   (0x1u << (sizeof (T_PeakFlags) * 8 - 1))
#define PF_IsPeakBit (0x100)
#define PF_iWL_Mask  (0xff)

Ex T_WyckoffList  *WyckoffList;
Ex int            nWyckoffList;
Ex int            MaxPotentialAtoms;
Ex int            MaxRecycledAtoms;
Ex int            MaxPeaksFwSearch;
Ex int            MaxPeaksFwFragmentSearch;
Ex int            MaxRawPeaks;
Ex int            NeD_PeakList;
Ex T_eD_PeakList  *eD_PeakList;
Ex int            nCeD_PeakList;
Ex int            nList__RawSymEquiv;
Ex T_fVector      *List_xRawSymEquiv;
Ex T_fVector      *List_cRawSymEquiv;
Ex Fprec          *NextPeakMx;
Ex Fprec          *MinDistance2Mx;
Ex Fprec          MinOfMD2Mx;

Ex T_Ticks  FT_TimeStart;
Ex T_Ticks  FT_TimeEnd;
Ex long     FT_SumTicks;

Ex T_Ticks  FwSearchTimeStart;
Ex T_Ticks  FwSearchTimeEnd;
Ex long     FwSearchSumTicks;

Ex int   nFourierTransform;
Ex long  nCatchRawPeak;
Ex long  CountInterSectionII;
Ex long  CountInterSectionIP;
Ex long  CountInterSectionPI;
Ex long  CountInterSectionPP;
Ex long  IllPeakInterpolation;
Ex long  SucPeakInterpolation;
Ex int   Re_i_Cycle;
Ex int   iFramework;
Ex int   nBadTetrahedraFW;
Ex int   nNo3DimConFW;
Ex int   nSmallLoopsFW;
Ex int   nRejectedOddLoopsFW;




#ifndef A_AppGlobal__
extern
const int TetrahedraNoSPO[];
#else
const int TetrahedraNoSPO[] = { 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1 };
                            /* -6    -4 -3 -2 -1     1  2  3  4     6 */
#endif


#undef Ex


#define PCT_Unknown  -1
#define PCT_None      0
#define PCT_NoChange  1
#define PCT_II        2
#define PCT_IP        3
#define PCT_PI        4
#define PCT_PP        5


typedef int T_PhaseCode;


typedef struct
  {
    int    Nx, Ny, Nz;
    int    Mx, My, Mz, Friedel;
    int    mWsp;
    float  *Wsp;
    int    mDensity;
    Fprec  *Density;
    int    SkipFrom_iy;
    int    SkipTo_iy;
    int    Max_iz;
    int    Testing;
    long   TicksFFTPACKgeneral;
    long   TicksFFTPACKcache;
    long   TicksFFTW;
    long   TicksSGImath;
    long   TicksDECdxml;
    int    UseType;
  }
  T_FourierParameters;

#define FPUT_unknown         (0)
#define FPUT_fftpackGeneral  (1)
#define FPUT_fftpackCache    (2)
#define FPUT_fftw            (3)
#define FPUT_sgimath         (4)
#define FPUT_decdxml         (5)
#define FPUT_evaluate        (6)


/* xtal.c */

int iH2H(int iH, int NH, int FullRange);
int IndexEq_hkl(T_Eq_hkl *Eq_hkl, int h, int k, int l, int *Conjugate);
int hkl_Compare(int a_h, int a_k, int a_l, int b_h, int b_k, int b_l);
void SortFobsRaw(T_FobsRaw *FR, int nFR);
void SetStatus1FobsRaw(T_FobsRaw *FR, int nFR);
void SetStatus2FobsRaw(T_FobsRaw *FR, int nFR);
void CalcFmrg(T_FobsRaw *FR, int nFR);
Fprec TwoThetaDeg(Fprec Q);
void RedivideAbsent(T_FobsRaw *FR, int nFR);
int SetOverlapGroups(T_FobsRaw *FR, int nFR);
void TreatOverlap(T_FobsRaw *FR, int nFR);
void SetF000(void);
Fprec InitCodeTransl(T_CodeTransl *CT, int nCT, T_FobsRaw *FR, int nFR);
void SetSleeping(T_CodeTransl *CT, int nCT);
void DoCodeTranslation(T_PhaseCode *PhaseCode,
                       Fprec *Density, int Mx, int My, int Mz, int Friedel);
void CalcCCandR(T_CodeTransl *CT, int nCT, Fprec *CC, Fprec *R);
void CheckRestrictedPhase(T_CodeTransl *CT, int *GC);
void SetupXtalInternals(void);


/* scatfact.c */

void SortSF_Tables(void);
Fprec CalcFormFactor(T_SF_Info *SFI, int h, int k, int l,
                     Fprec Occ, Fprec Uiso);
int CompleteSF_Info(T_SF_Info *SFI, int LblExact, int Set_f_stol_0);


/* fourier.c */

void InitFourierParameters(T_FourierParameters *FP,
                           const T_CodeTransl *CT, int nCT);
void FreeFourierParameters(T_FourierParameters *FP);
void FourierTransform(T_FourierParameters *FP);


/* sitefcal.c */

void CompleteSite(T_Site *S, int nS, T_AtomType *AT, int nAT);
void DoSiteFcal(int OutPutMode);
void CollectionLists(int PrintList_hkl,
                     int PowderStepScan, const char *fnStdPeak);
void DoSitePhases(void);


/* psearch.c */

void MapInd(int nMap, Fprec *DensMap, const T_PeakFlags *FlagMap,
            int ChkInDep, int MaxInDep, Fprec *CorrCoeff);
void Sort_eD_PeakList(int PSE, int nEPL);
void Free_eD_PeakList(void);
int Build_eD_PeakList(Fprec *eDensity, T_PeakFlags *PeakFlags,
                      int UseLargestFwFragment);
void SiteFrame(void);
void RandomSites(int NumberToGenerate);
void TestMoreSpecial(int NumberToGenerate);


/* posnmx.c */

void SeitzMx_of_iPos(int iPos, int *WL_Flag, T_RTMx *SMx, T_fRTMx *fSMx);
int IsUnitTr(Fprec Tr);
int iPos_of_SeitzMx(T_RTMx *SMx, int *WL_Flag, T_fVector *Pos);
void Reconstruct_fUC_Shift(int iEPL, int iPos,
                           T_RTMx *SMx, T_fVector *fUC_Shift_ii);


/* tetrahed.c */

void IdealTetrahedron(void);
int CheckTetrahedra(int Colored, const int *LnB, int Low_iEPL, int Top_iEPL,
                    int *CanNotBeTetrahedron);

/* findpnnb.c */

void BuildLPNNB(void);
void FreeLPNNB(void);
void PrintLPNNB(void);


/* setaos.c */

void SetAOS(void);


/* frame.c */

void FwSearch(int FindLargestFwFragment);


/* evalfw.c */

int CheckFragmentLoops(int Colored, const int *LnB, int Low_iEPL, int Top_iEPL);
void EvalFw(int Colored, const int *LnB, int Low_iEPL, int Top_iEPL);


/* sgextend.c */

void CalcLatticeTr2(T_LatticeConstants *LC, const T_LatticeInfo *LI,
                    Fprec *Min2, Fprec *Max2);
void PrintWyckoffList(void);
int SetToUnitLength(Fprec v[3]);
int CombineFreeEigenVects(const T_FreeVect *PreFV,
                                T_FreeVect *SymFV,
                                T_FreeVect *NewFV,
                          int iList, int iLoopInv);
void BuildWyckoffList(T_SgInfo *SgInfo);
int FindWL_Entry(T_WyckoffList *WL, T_WyckoffList *Buf_eWL);
void DS_MarkEquiv(T_SgInfo *SgInfo, int *FlagField, int n_x, int n_y, int n_z);
void SetCanBeCN(T_SgInfo *SgInfo);
int CalcSymEquiv(Fprec xx, Fprec yx, Fprec zx,
                 T_fVector *SE, int MaxSE,
                 Fprec     Dist2ConsiderSame,
                 Fprec *MaxDist2ConsideredSame,
                 Fprec *MinDist2Distinct,
                 int *Buf_WL_Flag, int MaxBuf_WL_Flag);
void CheckWL_Entry(T_eD_PeakList *EPL);
void fSymOps(T_SgInfo *SgInfo, T_fRTMx *List_fSMx, T_fVector *fTrV);
void cSymOps(T_SgInfo *SgInfo, T_fRTMx *List_cSMx, T_fVector *cTrV,
             Fprec *Txc, Fprec *Tcx);


#endif /* A_XTAL_H__ */
