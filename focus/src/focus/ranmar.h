#ifndef A_RANMAR_H__
#define A_RANMAR_H__

#ifndef Ex
#ifdef A_AppGlobal__
#define Ex
#else
#define Ex extern
#endif
#endif

Ex   int          nCallRanmar;

int rmarin(int ij, int kl);
double ranmar(void);

#define SlightlyLessThan(x) ((x) - 1.e-6)

#endif /* A_RANMAR_H__ */
