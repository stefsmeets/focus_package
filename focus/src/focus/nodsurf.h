#ifndef A_NODSURF_H__
#define A_NODSURF_H__

#ifndef Ex
#ifdef A_AppGlobal__
#define Ex
#else
#define Ex extern
#endif
#endif

Ex int    nNormalizeSurface;
Ex Fprec  vNormalizeSurface[3];

#undef Ex

void LoadSurface(Fprec *Surface, int Nx, int Ny, int Nz,
                 const T_PeakFlags *FlagMap);

#endif /* A_NODSURF_H__ */
