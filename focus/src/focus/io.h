#ifndef A_IO_H__
#define A_IO_H__

#ifndef Ex
#ifdef A_AppGlobal__
#define Ex
#else
#define Ex extern
#endif
#endif


#define RCF_OK                 (0)
#define RCF_MissingParameter   (1)
#define RCF_WrongParameter     (2)
#define RCF_UnrecognizedLine   (3)
#define RCF_ChangingInput      (4)
#define RCF_UnexpectedEOC      (5)


int BuildSpgrInfo(T_SgInfo *SgInfo, const char *SgName);
int ReadCmdFile(FILE *fpcmd);
char *EchoAtomTypeClass(FILE *fpout, int Class);
void EchoProfileSettings(FILE *fpout);
void EchoInput(FILE *fpout, int Full);
void PrintFmrg(FILE *fpout);
void Print_eEPL(FILE *fpout, T_eD_PeakList *EPL, int iEPL,
                int CoordinationNumber, int FlagAll);
void DoPut_strudat(T_eD_PeakList *EPL, int nEPL, int iCycle, char *Label,
                   int FlagAll);
void PutHistogram(int *Box, int nBox, char *legend);
void Put_eDmap(const Fprec *eD, char *note);
void PrintFramework(int Colored, const int *LnB, int Low_iEPL, int Top_iEPL,
                    const T_LargestFwFragment *LFwF);
void LoadDensity(Fprec *Density);
long LoadNextPhaseSet(int *Code);
void LoadStdPeak(const char *fnin, T_StdPeak *StdPeak);
void CoseqProtocol(int iEPL, int iSphere, int nNew);
void AdminCoseqSave(char *fncsve, int ID);
void CoseqSave(const T_ListSymNodes *LSymN, int FxFyFz,
               int iEPL, int SymN_Offset, int iSphere,
               int *SphereStorage, int nSphereStorage, int iAsyN);
int CoseqReload(const T_ListSymNodes *LSymN, int FxFyFz,
                int *iEPL, int *SymN_Offset, int *iSphere,
                int *SphereStorage, int nSphereStorage, int *iAsyN);
#undef Ex

#endif /* A_IO_H__ */
