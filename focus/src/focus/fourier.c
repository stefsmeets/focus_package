#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if defined(FFTfftw)
#include <fftw3.h>
#endif

#include "main.h"
#include "xtal.h"
#include "io.h"


#define MemCmp(t, s, n)  memcmp((t), (s), (n) * sizeof (*(t)))


#if    defined(__sgi)\
    || (defined(__alpha) && defined(__osf__))\
    || defined(__linux)
#define CALL_FORTRAN_UNDERSCORE
#elif defined(_CRAY)
#define CALL_FORTRAN_UPPERCASE
#endif


#define mWspR(N) (2 * N + 15)
#define mWspC(N) (4 * N + 15)


#if defined(FFTfftpack)

#define MaxSeq  256

#if   defined(CALL_FORTRAN_UNDERSCORE)

#define rffti srffti_
#define rfftf srfftf_
#define rfftb srfftb_
#define cffti scffti_
#define cfftb scfftb_

#elif defined(CALL_FORTRAN_UPPERCASE)

#define rffti SRFFTI
#define rfftf SRFFTF
#define rfftb SRFFTB
#define cffti SCFFTI
#define cfftb SCFFTB

#else

#define rffti srffti
#define rfftf srfftf
#define rfftb srfftb
#define cffti scffti
#define cfftb scfftb

#endif

extern void rffti(int *n, float *wsave);
extern void rfftf(int *n, float *r, float *wsave);
extern void rfftb(int *n, float *r, float *wsave);
extern void cffti(int *n, float *wsave);
extern void cfftb(int *n, float *c, float *wsave);

#endif


#if defined(FFTsgimath)

extern float *sfft3dui(int n1, int n2, int n3, float *wsp);

extern int sfft3du(int job, int n1, int n2, int n3, float *seq,
                   int ld1, int ld2, float *wsp);
#endif


#if defined(FFTdecdxml)

#if   defined(CALL_FORTRAN_UNDERSCORE)

#define sfft_3d sfft_3d_

#elif defined(CALL_FORTRAN_UPPERCASE)

#define sfft_3d SFFT_3D_

#endif

int sfft_3d(char *input_format, char *output_format, char *direction,
            float *in, float *out,
            int *ni, int *nj, int *nk,
            int *lda_i, int *lda_j,
            int *ni_stride, int *nj_stride, int *nk_stride,
            int len_input_format, int len_output_format, int len_direction);

#endif


void trpcmplx(int nx, int ny, int mya, const float *a,
                              int mxb,       float *b);


static void SetSkipAndMaxFP(T_FourierParameters *FP,
                            const T_CodeTransl *CT, int nCT)
{
  int             iCT, iCTD;
  const T_CTData  *CTD;
  int             ik, il, Friedel;


  FP->SkipFrom_iy = 0;
  FP->SkipTo_iy = FP->Ny;
  FP->Max_iz = 0;

  for (iCT = 0; iCT < nCT; iCT++, CT++)
  {
    for (iCTD = 0, CTD = CT->CTData; iCTD < CT->nCTData; iCTD++, CTD++)
    {
      ik = CTD->ik;
      il = CTD->il;

      if (FP->Friedel == -1 || (FP->Friedel == 1 && il == 0))
        Friedel = 1;
      else
        Friedel = 0;

      for (;;)
      {
        if (iH2H(ik, FP->Ny, 1) >= 0)
        {
          if (FP->SkipFrom_iy <= ik)
              FP->SkipFrom_iy =  ik + 1;
        }
        else
        {
          if (FP->SkipTo_iy > ik)
              FP->SkipTo_iy = ik;
        }

        if (FP->Max_iz < il)
            FP->Max_iz = il;

        if (Friedel) {
            Friedel = 0;
          ik = (Ny - CTD->ik) % Ny;
          il = (Nz - CTD->il) % Nz;
        }
        else
          break;
      }
    }
  }

  Fprintf(stdout, "# Ny = %d  SkipFrom_iy = %d  SkipTo_iy = %d\n",
    FP->Ny, FP->SkipFrom_iy, FP->SkipTo_iy);
  Fprintf(stdout, "# Nz = %d  Max_iz = %d\n\n",
    FP->Nz, FP->Max_iz);
}


static void SetNextUseType(T_FourierParameters *FP)
{
  FP->UseType++;

  switch (FP->UseType)
  {
    case FPUT_fftpackGeneral:
#if defined(FFTfftpack)
      break;
#else
      FP->UseType++;
#endif

    case FPUT_fftpackCache:
#if defined(FFTfftpack)
      break;
#else
      FP->UseType++;
#endif

    case FPUT_fftw:
#if defined(FFTfftw)
      break;
#else
      FP->UseType++;
#endif

    case FPUT_sgimath:
#if defined(FFTsgimath)
      break;
#else
      FP->UseType++;
#endif

    case FPUT_decdxml:
#if defined(FFTdecdxml)
      break;
#else
      FP->UseType++;
#endif

    case FPUT_evaluate:
      break;

    default:
      InternalError(NULL);
      break;
  }
}


void InitFourierParameters(T_FourierParameters *FP,
                           const T_CodeTransl *CT, int nCT)
{
  FP->Nx = Nx;
  FP->Ny = Ny;
  FP->Nz = Nz;

  FP->Mx = FP->Nx;
  FP->My = FP->Ny;
  FP->Mz = FP->Nz + 2 - (FP->Nz % 2); /* n+2-(n%2) = (n/2+1)*2 */
  FP->Friedel = 1;

  FP->mWsp = mWspC(FP->Nx) + mWspC(FP->Ny) + mWspR(FP->Nz);
  CheckMalloc(FP->Wsp, FP->mWsp);

  FP->mDensity = FP->Mx * FP->My * FP->Mz;
  CheckMalloc(FP->Density, FP->mDensity);

  SetSkipAndMaxFP(FP, CT, nCT);

  FP->Testing = 1;
  FP->TicksFFTPACKgeneral = -1;
  FP->TicksFFTPACKcache   = -1;
  FP->TicksFFTW           = -1;
  FP->TicksSGImath        = -1;
  FP->TicksDECdxml        = -1;
  FP->UseType = FPUT_unknown;

  SetNextUseType(FP);

  if (FP->UseType == FPUT_evaluate)
    InternalError("Missing FFT library");
}


void FreeFourierParameters(T_FourierParameters *FP)
{
  if (FP->Density) AppFree(FP->Density, FP->mDensity);
      FP->Density = NULL;

  if (FP->Wsp) AppFree(FP->Wsp, FP->mWsp);
      FP->Wsp = NULL;
}


#if defined(FFTfftpack)

static void ccTTfftpack(T_FourierParameters *FP)
{
  int    ix, iy, iz;
  int    idx, idy;
  float  Wspx[mWspC(MaxSeq)];
  float  Wspy[mWspC(MaxSeq)];
  int    mXYplane;
  float  *XYplane, *XYplTrp, *XY;


  if (FP->Nx > MaxSeq || FP->Ny > MaxSeq)
    InternalError(NULL);

  cffti(&FP->Nx, Wspx);
  cffti(&FP->Ny, Wspy);

                           mXYplane = 2 * FP->Nx * FP->Ny;
  CheckMalloc(XYplane, 2 * mXYplane);
  XYplTrp = XYplane + mXYplane;

  for (iz = 0; iz <= FP->Nz; iz += 2)
  {
    if (iz > FP->Max_iz * 2)
      break;

    XY = XYplane;
    idx = iz;
    for (ix = 0; ix < FP->Nx; ix++, idx += FP->My * FP->Mz) {
      idy = idx;
      for (iy = 0; iy < FP->Ny; iy++, idy += FP->Mz) {
        *XY++ = FP->Density[idy    ];
        *XY++ = FP->Density[idy + 1];
      }
    }

    trpcmplx(FP->Nx, FP->Ny, 2 * FP->Ny, XYplane,
                             2 * FP->Nx, XYplTrp);

    XY = XYplTrp;
    for (iy = 0; iy < FP->Ny; iy++, XY += 2 * FP->Nx) {
      if (iy == FP->SkipFrom_iy) {
          iy =  FP->SkipTo_iy;
          XY = XYplTrp + 2 * FP->Nx * iy;
        if (iy == FP->Ny)
          break;
      }
      cfftb(&FP->Nx, XY, Wspx);
    }

    trpcmplx(FP->Ny, FP->Nx, 2 * FP->Nx, XYplTrp,
                             2 * FP->Ny, XYplane);

    XY = XYplane;
    for (ix = 0; ix < FP->Nx; ix++, XY += 2 * FP->Ny)
      cfftb(&FP->Ny, XY, Wspy);

    XY = XYplane;
    idx = iz;
    for (ix = 0; ix < FP->Nx; ix++, idx += FP->My * FP->Mz) {
      idy = idx;
      for (iy = 0; iy < FP->Ny; iy++, idy += FP->Mz) {
        FP->Density[idy    ] = *XY++;
        FP->Density[idy + 1] = *XY++;
      }
    }
  }

  AppFree(XYplane, 2 * mXYplane);
}


static void ccXYfftpack(T_FourierParameters *FP)
{
  int    ix, iy, iz;
  int    id;
  float  Seq[2 * MaxSeq];
  float  Wspx[mWspC(MaxSeq)];
  float  Wspy[mWspC(MaxSeq)];


  if (FP->Nx > MaxSeq || FP->Ny > MaxSeq)
    InternalError(NULL);

  cffti(&FP->Nx, Wspx);
  cffti(&FP->Ny, Wspy);

  for (iz = 0; iz <= FP->Nz; iz += 2)
  {
    if (iz > FP->Max_iz * 2)
      break;

    for (iy = 0; iy < FP->Ny; iy++)
    {
      if (iy == FP->SkipFrom_iy) {
          iy =  FP->SkipTo_iy;
        if (iy == FP->Ny)
          break;
      }

      id = iy * FP->Mz + iz;
      for (ix = 0; ix < FP->Nx; ix++, id += FP->My * FP->Mz) {
        Seq[2 * ix    ] = FP->Density[id    ];
        Seq[2 * ix + 1] = FP->Density[id + 1];
      }

      cfftb(&FP->Nx, Seq, Wspx);

      id = iy * FP->Mz + iz;
      for (ix = 0; ix < FP->Nx; ix++, id += FP->My * FP->Mz) {
        FP->Density[id    ] = Seq[2 * ix    ];
        FP->Density[id + 1] = Seq[2 * ix + 1];
      }
    }

    for (ix = 0; ix < FP->Nx; ix++)
    {
      id = ix * FP->My * FP->Mz + iz;
      for (iy = 0; iy < FP->Ny; iy++, id += FP->Mz) {
        Seq[2 * iy    ] = FP->Density[id    ];
        Seq[2 * iy + 1] = FP->Density[id + 1];
      }

      cfftb(&FP->Ny, Seq, Wspy);

      id = ix * FP->My * FP->Mz + iz;
      for (iy = 0; iy < FP->Ny; iy++, id += FP->Mz) {
        FP->Density[id    ] = Seq[2 * iy    ];
        FP->Density[id + 1] = Seq[2 * iy + 1];
      }
    }
  }
}


static void crZfftpack(T_FourierParameters *FP)
{
  int    id, ix, iy, iz;
  float  *Ds, *Dt;
  float  Wsp[mWspR(MaxSeq)];


  if (FP->Nz > MaxSeq)
    InternalError(NULL);

  rffti(&FP->Nz, Wsp);

  for (ix = 0; ix < FP->Nx; ix++)
  for (iy = 0; iy < FP->Ny; iy++)
  {
                      id = (ix * FP->My + iy) * FP->Mz;
    Dt = &FP->Density[id + 1];
    Ds = &FP->Density[id + 2];

    for (iz = 0; iz < FP->Nz; iz++)
      *Dt++ = *Ds++;

    rfftb(&FP->Nz, &FP->Density[id], Wsp);
  }
}


#if defined(REAL2REAL)

static void DisplayDensity(int nx, int ny, int nz, int my, int mz,
                           const float *dens, const char *label)
{
  int  ix, iy, iz, id;


  Fprintf(stdout, ">Begin %s\n", label);

  for (ix = 0; ix < nx; ix++)
  {
    for (iy = 0; iy < ny; iy++)
    {
      Fprintf(stdout, "%2.2d,%2.2d", ix, iy);

      for (iz = 0; iz < nz; iz++)
      {
        id = (ix * my + iy) * mz + iz;

        Fprintf(stdout, " %7.2f", dens[id]);
      }

      putc('\n', stdout);
    }

    putc('\n', stdout);
  }

  Fprintf(stdout, ">End %s\n\n", label);
  Fflush(stdout);
}


static void rrInitFourierParameters(T_FourierParameters *FP,
                                    int Nx, int Ny, int Nz)
{
  FP->Nx = Nx;
  FP->Ny = Ny;
  FP->Nz = Nz;

  FP->Mx = FP->Nx / 2 + 1;
  FP->My = FP->Ny;
  FP->Mz = FP->Nz + 2 - (FP->Nz % 2);
  FP->Friedel = 1;

  FP->mWsp = 0;
  FP->Wsp  = NULL;

  FP->mDensity = FP->Mx * FP->My * FP->Mz;
  CheckMalloc(FP->Density, FP->mDensity);

  FP->SkipFrom_iy = FP->My;
  FP->SkipTo_iy = FP->My;
  FP->Max_iz = FP->Mz;

  FP->Testing = 1;
  FP->TicksFFTPACKgeneral = -1;
  FP->TicksFFTPACKcache   = -1;
  FP->TicksFFTW           = -1;
  FP->TicksSGImath        = -1;
  FP->TicksDECdxml        = -1;
  FP->UseType = FPUT_unknown;
}


static void rcXfftpack(T_FourierParameters *FP)
{
  int    id, ix, iy, iz;
  float  Seq[MaxSeq];
  float  Wsp[mWspR(MaxSeq)];


  if (FP->Nx > MaxSeq)
    InternalError(NULL);

  rffti(&FP->Nx, Wsp);

  for (iy = 0; iy < FP->Ny; iy++)
  {
    for (iz = 0; iz <= FP->Nz; iz += 2)
    {
      id = iy * FP->Mz + iz;

      for (ix = 0; ix < FP->Nx;)
      {
        Seq[ix] = FP->Density[id];
            ix++;

        if (ix == FP->Nx)
          break;

        Seq[ix] = FP->Density[id + 1];
            ix++;             id += FP->My * FP->Mz;
      }

      rfftf(&FP->Nx, Seq, Wsp);

      id = iy * FP->Mz + iz;

      FP->Density[id    ] = Seq[0];
      FP->Density[id + 1] = 0.F;
                  id += FP->My * FP->Mz;

      ix = 1;

      for (;;)
      {
        if (ix == FP->Nx)
          break;

        FP->Density[id] =  Seq[ix];
                               ix++;
        if (ix == FP->Nx) {
          FP->Density[id + 1] = 0.F;
          break;
        }

        FP->Density[id + 1] = -Seq[ix];
                                   ix++;
                    id += FP->My * FP->Mz;
      }
    }
  }

  FP->Nx = FP->Mx;
}


static void ccYfftpack(T_FourierParameters *FP)
{
  int    id, ix, iy, iz;
  float  Seq[2 * MaxSeq];
  float  Wsp[mWspC(MaxSeq)];


  if (FP->Ny > MaxSeq)
    InternalError(NULL);

  cffti(&FP->Ny, Wsp);

  for (ix = 0; ix <  FP->Nx; ix++)
  for (iz = 0; iz <= FP->Nz; iz += 2)
  {
    if (iz > FP->Max_iz * 2)
      break;

    id = ix * FP->My * FP->Mz + iz;
    for (iy = 0; iy < FP->Ny; iy++, id += FP->Mz) {
      Seq[2 * iy    ] = FP->Density[id    ];
      Seq[2 * iy + 1] = FP->Density[id + 1];
    }

    cfftb(&FP->Ny, Seq, Wsp);

    id = ix * FP->My * FP->Mz + iz;
    for (iy = 0; iy < FP->Ny; iy++, id += FP->Mz) {
      FP->Density[id    ] = Seq[2 * iy    ];
      FP->Density[id + 1] = Seq[2 * iy + 1];
    }
  }
}


static void RealToReal(T_FourierParameters *FP)
{
  int    ix, iy, iz;
  int    ie, io, je, jo;

  T_FourierParameters  rrFP[1];


  rrInitFourierParameters(rrFP, FP->Nx, FP->Ny, FP->Nz);

  for (ix = 0; ix < rrFP->Mx; ix++)
  for (iy = 0; iy < rrFP->Ny; iy++)
  {
    ie = (ix * rrFP->My + iy) * rrFP->Mz;
    io = ie + 1;

    je = (2 * ix * FP->My + iy) * FP->Mz;
    jo = je +  FP->My * FP->Mz;

    for (iz = 0; iz <= rrFP->Nz; iz += 2)
    {
      rrFP->Density[ie] = FP->Density[je];

      if (2 * ix + 1 < FP->Nx)
        rrFP->Density[io] = FP->Density[jo];

      ie += 2;
      io += 2;
      je += 2;
      jo += 2;
    }
  }

  rcXfftpack(rrFP);
  ccYfftpack(rrFP);
  crZfftpack(rrFP);

  DisplayDensity(rrFP->Nx, rrFP->Ny, rrFP->Nz, rrFP->My, rrFP->Mz,
                 rrFP->Density, "rrDensity");

  FreeFourierParameters(rrFP);
}

#endif /* defined(REAL2REAL) */
#endif /* defined(FFTfftpack) */


void FourierTransform(T_FourierParameters *FP)
{
  int    ix, iy, iz;
  int    i, j;
  Fprec  eDscale;
  long   MinTicks;

#if defined(FFTfftw)
  fftwf_plan fftw_plan_;
#endif


  FP->Density[0] = F000mtp;

  if (FP->Testing)
  {
    const char  *what;

    if (FP->Testing < 0) {
        FP->Testing = 0;
        what = "Using";
    }
    else
        what = "Timing";

    switch (FP->UseType)
    {
#if defined(FFTfftpack)
      case FPUT_fftpackGeneral:
        Fprintf(stdout, "# FFT: %s libfftpack(General)\n\n", what);
        Fflush(stdout);
        break;

      case FPUT_fftpackCache:
        Fprintf(stdout, "# FFT: %s libfftpack(Cache)\n\n", what);
        Fflush(stdout);
        break;
#endif

#if defined(FFTfftw)
      case FPUT_fftw:
        Fprintf(stdout, "# FFT: %s fftw\n\n", what);
        Fflush(stdout);
        if (   FP->Mx != FP->Nx
            || FP->My != FP->My
            || FP->Mz != FP->Nz + 2 || FP->Nz % 2 != 0)
          progerror("Improper array dimensions for fftw.");
        break;
#endif

#if defined(FFTsgimath)
      case FPUT_sgimath:
        Fprintf(stdout, "# FFT: %s libcomplib.sgimath\n\n", what);
        Fflush(stdout);
        (void) sfft3dui(FP->Nz, FP->Ny, FP->Nx, FP->Wsp);
        break;
#endif

#if defined(FFTdecdxml)
      case FPUT_decdxml:
        Fprintf(stdout, "# FFT: %s dxml\n\n", what);
        Fflush(stdout);
        if (FP->Nz % 2 != 0)
          progerror("Number of grid points in z-direction must be even.");
        break;
#endif

      default:
        InternalError(NULL);
        break;
    }
  }

  (void) GetTicks(&FT_TimeStart);

  switch (FP->UseType)
  {
#if defined(FFTfftpack)
    case FPUT_fftpackGeneral:
#if defined(REAL2REAL)
      RealToReal(FP);
#endif
      ccXYfftpack(FP);
      crZfftpack(FP);
#if defined(REAL2REAL)
      DisplayDensity(FP->Nx, FP->Ny, FP->Nz, FP->My, FP->Mz,
                     FP->Density, "Density");
#endif
      break;

    case FPUT_fftpackCache:
      ccTTfftpack(FP);
      crZfftpack(FP);
      break;
#endif

#if defined(FFTfftw)
    case FPUT_fftw:
      fftw_plan_ = fftwf_plan_dft_c2r_3d(
        FP->Nx, FP->Ny, FP->Nz,
        (fftw_complex *) FP->Density,
        FP->Density,
        FFTW_ESTIMATE);
      fftwf_execute(fftw_plan_);
      fftwf_destroy_plan(fftw_plan_);
      break;
#endif

#if defined(FFTsgimath)
    case FPUT_sgimath:
      sfft3du(1, FP->Nz, FP->Ny, FP->Nx, FP->Density,
                 FP->Mz, FP->My, FP->Wsp);
      break;
#endif

#if defined(FFTdecdxml)
    case FPUT_decdxml:
      i = 1;
      if (sfft_3d("c", "r", "b",
                  FP->Density, FP->Density,
                  &FP->Nz, &FP->Ny, &FP->Nx, &FP->Mz, &FP->My,
                  &i, &i, &i,
                  1, 1, 1) != 0)
        InternalError("status(sfft_3d) != 0");
      break;
#endif

    default:
      InternalError(NULL);
      break;
  }

  (void) GetTicks(&FT_TimeEnd);

  FT_SumTicks += FT_TimeEnd.cpu_user   - FT_TimeStart.cpu_user;
  FT_SumTicks += FT_TimeEnd.cpu_system - FT_TimeStart.cpu_system;

  FT_TimeStart.second = -1;

  nFourierTransform++;

  if (FP->UseType != FPUT_decdxml)
    eDscale = 1. / LatConD.v;
  else
    eDscale = 1. / LatConD.v * FP->Nx * FP->Ny * FP->Nz;

  i = 0;

  for (ix = 0; ix < FP->Nx; ix++)
  for (iy = 0; iy < FP->Ny; iy++)
  {
    j = (ix * FP->My + iy) * FP->Mz;

    for (iz = 0; iz < FP->Nz; iz++, i++, j++)
      FP->Density[i] = FP->Density[j] * eDscale;
  }

  if (FP->Testing)
  {
    switch (FP->UseType)
    {
#if defined(FFTfftpack)
      case FPUT_fftpackGeneral:
        FP->TicksFFTPACKgeneral = FT_TimeEnd.cpu_user - FT_TimeStart.cpu_user;
        Fprintf(stdout, "# FFT: Timing libfftpack(General): %.6g s\n\n",
          (double) FP->TicksFFTPACKgeneral / FT_TimeStart.ticks_per_second);
        Fflush(stdout);
        break;

      case FPUT_fftpackCache:
        FP->TicksFFTPACKcache = FT_TimeEnd.cpu_user - FT_TimeStart.cpu_user;
        Fprintf(stdout, "# FFT: Timing libfftpack(Cache): %.6g s\n\n",
          (double) FP->TicksFFTPACKcache / FT_TimeStart.ticks_per_second);
        Fflush(stdout);
        break;
#endif

#if defined(FFTfftw)
      case FPUT_fftw:
        FP->TicksFFTW = FT_TimeEnd.cpu_user - FT_TimeStart.cpu_user;
        Fprintf(stdout, "# FFT: Timing fftw: %.6g s\n\n",
          (double) FP->TicksFFTW / FT_TimeStart.ticks_per_second);
        Fflush(stdout);
        break;
#endif

#if defined(FFTsgimath)
      case FPUT_sgimath:
        FP->TicksSGImath = FT_TimeEnd.cpu_user - FT_TimeStart.cpu_user;
        Fprintf(stdout, "# FFT: Timing libcomplib.sgimath: %.6g s\n\n",
          (double) FP->TicksSGImath / FT_TimeStart.ticks_per_second);
        Fflush(stdout);
        break;
#endif

#if defined(FFTdecdxml)
      case FPUT_decdxml:
        FP->TicksDECdxml = FT_TimeEnd.cpu_user - FT_TimeStart.cpu_user;
        Fprintf(stdout, "# FFT: Timing dxml: %.6g s\n\n",
          (double) FP->TicksDECdxml / FT_TimeStart.ticks_per_second);
        Fflush(stdout);
        break;
#endif

      default:
        InternalError(NULL);
        break;
    }

    SetNextUseType(FP);

    if (FP->UseType == FPUT_evaluate)
    {
      MinTicks = -1;

      if ((MinTicks < 0 || MinTicks > FP->TicksFFTPACKgeneral)
                                  &&  FP->TicksFFTPACKgeneral >= 0) {
          MinTicks = FP->TicksFFTPACKgeneral;
          FP->UseType =  FPUT_fftpackGeneral;
      }

      if ((MinTicks < 0 || MinTicks > FP->TicksFFTPACKcache)
                                  &&  FP->TicksFFTPACKcache >= 0) {
          MinTicks = FP->TicksFFTPACKcache;
          FP->UseType =  FPUT_fftpackCache;
      }

      if ((MinTicks < 0 || MinTicks > FP->TicksFFTW)
                                   && FP->TicksFFTW >= 0) {
          MinTicks = FP->TicksFFTW;
          FP->UseType =  FPUT_fftw;
      }

      if ((MinTicks < 0 || MinTicks > FP->TicksSGImath)
                                   && FP->TicksSGImath >= 0) {
          MinTicks = FP->TicksSGImath;
          FP->UseType =  FPUT_sgimath;
      }

      if ((MinTicks < 0 || MinTicks > FP->TicksDECdxml)
                                   && FP->TicksDECdxml >= 0) {
          MinTicks = FP->TicksDECdxml;
          FP->UseType =  FPUT_decdxml;
      }

      if (MinTicks < 0)
        InternalError(NULL);

      FP->Testing = -1;
    }
  }
}


#if defined(CALL_FORTRAN_UNDERSCORE)
#undef CALL_FORTRAN_UNDERSCORE
#endif
#if defined(CALL_FORTRAN_UPPERCASE)
#undef CALL_FORTRAN_UPPERCASE
#endif
