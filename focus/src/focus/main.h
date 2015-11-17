#ifndef A_MAIN_H__
#define A_MAIN_H__

#ifndef Ex
#ifdef A_AppGlobal__
#define Ex
#else
#define Ex extern
#endif
#endif


#ifdef __TURBOC__
typedef unsigned long  size_m;
#else
typedef size_t  size_m;
#endif

#if defined(__sgi) || defined(__PARAGON__)
#define USE_FUNCF
#endif


Ex    char           *TopTitle;

Ex    char           *progn;
Ex    size_m         MemoryUsed;
Ex    size_m         MaxMemoryUsed;
Ex    volatile int   QuitProgram;
Ex    int            ReportOnExit;
Ex    int            CountWarnings;

Ex    int            Debug0;
Ex    int            Debug1;


void p_rogerror(const char *message, const char *file, const int line);
void I_nternalError(const char *message, const char *file, const int line);
void N_otEnoughCore(const char *file, const int line);
void I_llegalLine(const char *fnin, long lcount,
                  const char *file, const int line);

#define progerror(m)      p_rogerror(m, __FILE__, __LINE__)
#define InternalError(m)  I_nternalError(m, __FILE__, __LINE__)
#define NotEnoughCore()   N_otEnoughCore(__FILE__, __LINE__)
#define IllegalLine(f, l) I_llegalLine(f, l, __FILE__, __LINE__)

void AppExit(int status);

void AppSignal(void);
void CheckSignalFile(void);

void *App_Malloc(int n, size_m size, char *var_name);
void *App_Realloc(void *ptr, int nn, int n, size_m size,
                  char *nptr_name, char *ptr_name);
void App_Free(void *ptr, size_m n, size_m size, char *var_name);
char *App_Strdup(char *s, char *var_name);

#define AppMalloc(ptr, n) ptr = App_Malloc(n, sizeof (*(ptr)), # ptr)
#define AppRealloc(nptr, ptr, nn, n)\
  nptr = App_Realloc(ptr, nn, n, sizeof (*(ptr)), # nptr, # ptr)
#define AppFree(ptr, n) App_Free(ptr, n, sizeof (*(ptr)), # ptr)
#define AppStrdup(s, sdup) sdup = App_Strdup(s, # sdup)

#define CheckMalloc(ptr, n)\
  {\
    (ptr) = App_Malloc(n, sizeof (*(ptr)), # ptr);\
    if ((ptr) == NULL && (n) != 0) NotEnoughCore();\
  }

#define CheckRealloc(nptr, ptr, nn, n)\
  {\
    (nptr) = App_Realloc(ptr, nn, n, sizeof (*(ptr)), # nptr, # ptr);\
    if ((nptr) == NULL && (nn) != 0) NotEnoughCore();\
  }

void DoRandomInitialization(void);

/* some general definitions */

#define Fputs   (void) fputs
#define Fprintf (void) fprintf
#define Fflush  (void) fflush
#define Fclose  (void) fclose
#define Sprintf (void) sprintf

#ifndef SEEK_SET
#define SEEK_SET 0
#endif

#ifdef MIN
#undef MIN
#endif
#ifdef MAX
#undef MAX
#endif
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2          1.57079632679489661923
#endif
#define PIover180 (M_PI / 180.)
#define TwoPI (2. * M_PI)
#define Square(x) ((x) * (x))
#define Cubic(x) ((x) * (x) * (x))

#ifndef USE_FUNCF
#define fmodf(x, y) ((Fprec) fmod((double) (x), (double) (y)))
#define cosf(x) ((Fprec) cos((double) (x)))
#define sinf(x) ((Fprec) sin((double) (x)))
#define tanf(x) ((Fprec) tan((double) (x)))
#define acosf(x) ((Fprec) acos((double) (x)))
#define asinf(x) ((Fprec) asin((double) (x)))
#define atan2f(y, x) ((Fprec) atan2((double) (y), (double) (x)))
#define sqrtf(x) ((Fprec) sqrt((double) (x)))
#define logf(x) ((Fprec) log((double) (x)))
#define expf(x) ((Fprec) exp((double) (x)))
#endif

#ifndef FP_USE_DOUBLE
#define Fprec float
#define AppFmod(x, y) fmodf(x, y)
#define AppCos(x) cosf(x)
#define AppSin(x) sinf(x)
#define AppTan(x) tanf(x)
#define AppAcos(x) acosf(x)
#define AppAsin(x) asinf(x)
#define AppAtan2(y, x) atan2f(y, x)
#define AppSqrt(x) sqrtf(x)
#define AppLog(x) logf(x)
#define AppExp(x) expf(x)
#else
#define Fprec double
#define AppFmod(x, y) fmod(x, y)
#define AppCos(x) cos(x)
#define AppSin(x) sin(x)
#define AppTan(x) tan(x)
#define AppAcos(x) acos(x)
#define AppAsin(x) asin(x)
#define AppAtan2(y, x) atan2(y, x)
#define AppSqrt(x) sqrt(x)
#define AppLog(x) log(x)
#define AppExp(x) exp(x)
#endif /* FP_USE_DOUBLE */


#define CF0  ((Fprec) 0.)
#define CF1  ((Fprec) 1.)

#define AppFabs(Value) ((Value) < CF0 ? -(Value) : (Value))

#define InInterval(Value, Center, Delta)\
  ((Center) - (Delta) <= (Value) && (Value) <= (Center) + (Delta))

#define INT_ROUNDED(x_) ((int)((x_) + (Fprec)((x_) < CF0 ? -.5 : +.5)))


#define SignalFileCheckInterval 120 /* seconds between SignalFile lookup's */


typedef int (*SortFunction)(const void *, const void *);


Ex   char  *F_SignalFileName;
Ex   int    F_SiteFrame;
Ex   int    F_SiteLabel;
Ex   int    F_AllSitesOn;
Ex   int    F_SitePhases;
Ex   int    F_AllPhaseCodes0;
Ex   int    F_SetAllFmrg;
Ex   Fprec  F_SetAllFmrgValue;
Ex   int    F_CoseqValue;
Ex   int    F_CoseqSelect;
Ex   int    F_CoseqUnitCell;
Ex   char  *F_CoseqSaveFileName;
Ex   int    F_CoseqSaveTime;
Ex   int    F_CoseqKeep;
Ex   char  *F_CoseqReloadFileName;
Ex   char  *F_CoseqProtocolFileName;

Ex   int   *F_CoseqSplit;
Ex   int    F_nCoseqSplit;
Ex   int    F_CoseqSplit_x;
Ex   int    F_CoseqSplit_y;
Ex   int    F_CoseqSplit_z;
Ex   int    F_Put_strudat;
Ex   char  *F_Put_strudatFileName;
Ex   int    F_PutAllPeaks;
Ex   int    F_ShowLargestFwFragment;
Ex   int    F_eD_Histogram;
Ex   char  *F_eD_HistogramFileName;
Ex   int    F_Put_eDmap;
Ex   char  *F_Put_eDmapFileName;
Ex   int    F_nAllAsyPoints;
Ex   int    *F_AllAsyPoints;
Ex   char  *F_LoadDensityFileName;
Ex   char  *F_PhaseSetsFileName;
Ex   char  *F_SurfaceFileName;


#undef Ex

#endif /* A_MAIN_H__ */
