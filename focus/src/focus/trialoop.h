#ifndef A_TRAILOOP_H__
#define A_TRAILOOP_H__

#ifndef Ex
#ifdef A_AppGlobal__
#define Ex
#else
#define Ex extern
#endif
#endif


/* input values */

Ex long  RandomInitialization;
Ex int   nRandomBlindCalls;
Ex long  *RandomBlindCalls;
Ex int   FlagTimeRandomInitialization;


/* internal only */

Ex int    nTrials;
Ex int    nConvergedPhases;
Ex int    nNotConvergedPhases;
Ex Fprec  CurrCorrelationCoefficient;
Ex long   CurrPhaseCode_nCallRanmar;


void TriaLoop(int MaximumTrials, int nTrialsPerPhaseSet);

void InitCalcIRius(void);
void RecyclePhases(T_FourierParameters *FP,
                   T_PeakFlags *PeakFlags, Fprec *Surface,
                   T_PhaseCode *PhaseCode);

#undef Ex

#endif /* A_TRAILOOP_H__ */
