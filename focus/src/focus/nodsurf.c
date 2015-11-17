#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>


#include "main.h"
#include "xtal.h"
#include "nodsurf.h"


static void NormalizeSurface(Fprec *Surface, int Nx, int Ny, int Nz)
{
  int    nD, iD;
  Fprec  densmax, densmin, ddens, f, g, m;


  nD = Nx * Ny * Nz;

          densmin = densmax = 0.;
  if (nD) densmin = densmax = Surface[0];

  for (iD = 1; iD < nD; iD++)
    if      (Surface[iD] < densmin) densmin = Surface[iD];
    else if (Surface[iD] > densmax) densmax = Surface[iD];

  Fprintf (stdout, "#      Input Surface Min/Max  %.6g %.6g\n",
    densmin, densmax);

  if      (nNormalizeSurface == 1)
  {
    if (densmax != 0.)
    {
      f = vNormalizeSurface[0] / densmax;

      for (iD = 0; iD < nD; iD++)
        Surface[iD] *= f;
    }
  }
  else if (nNormalizeSurface == 2 || nNormalizeSurface == 3)
  {
        ddens = densmax - densmin;
    if (ddens != 0.)
    {
      m = vNormalizeSurface[1] - vNormalizeSurface[0];
      f = m / ddens;
      g = (vNormalizeSurface[0] + vNormalizeSurface[1]) * TwoPI / m;

      for (iD = 0; iD < nD; iD++)
      {
        Surface[iD] = (Surface[iD] - densmin) * f + vNormalizeSurface[0];

        if (nNormalizeSurface == 3)
          Surface[iD] = (  AppAtan2(vNormalizeSurface[2]
                                    * (Surface[iD] - g), 1.)
                         + M_PI_2) / M_PI_2;
      }
    }
  }

  if (nNormalizeSurface)
  {
            densmin = densmax = 0.;
    if (nD) densmin = densmax = Surface[0];

    for (iD = 1; iD < nD; iD++)
      if      (Surface[iD] < densmin) densmin = Surface[iD];
      else if (Surface[iD] > densmax) densmax = Surface[iD];

    Fprintf (stdout, "# Normalized Surface Min/Max  %.6g %.6g\n",
      densmin, densmax);
  }

  putc('\n', stdout);
}


void LoadSurface(Fprec *Surface, int Nx, int Ny, int Nz,
                 const T_PeakFlags *FlagMap)
{
  int          ix, iy, iz, iD, dim, NFile[6];
  int          c = EOF;
  double       val;
  T_PeakFlags  *FM;
  Fprec        CorrCoeff;



  static FILE  *fpdens = NULL;


  if (fpdens == NULL)
  {
        fpdens = fopen(F_SurfaceFileName, "r");
    if (fpdens == NULL)
    {
      Fprintf(stderr, "%s: Can't open surface file %s\n",
        progn, F_SurfaceFileName);
      AppExit(1);
    }
  }

/*
Reading of grid dimension
*/

for (dim = 1; dim <= 3; dim++)
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
      c = fscanf(fpdens, "%d", &NFile[dim]);

  }

    Fprintf(stdout, "# Grid_xyz  %d %d %d\n", Nx, Ny, Nz);
    Fprintf(stdout, "# %s xyz  %d %d %d\n",
      F_SurfaceFileName, NFile[1], NFile[2], NFile[3]);

    if ((Nx != NFile[1]) || (Ny != NFile[2]) || (Nz != NFile[3]))
    {
         Fprintf(stdout, "\n# Error: non compatible grid in %s\n\n",
         F_SurfaceFileName);
         AppExit(1);
    }

/*  NFile[1] = NFile[2] = NFile[3] = 0; */

/*
Reading of the iso values of the surface
*/

  for (ix = 0; ix < Nx; ix++)
  for (iy = 0; iy < Ny; iy++)
  for (iz = 0; iz < Nz; iz++)

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

          Fprintf(stdout, "\n# End of surface file %s\n\n",
          F_SurfaceFileName);
          AppExit(0);


      }



      Fprintf(stderr, "%s: Not enough data in surface file %s\n",
      progn, F_SurfaceFileName);
      AppExit(1);
    }
    else if (c != 1)
    {

      Fprintf(stderr, "%s: Illegal data in surface file %s\n",
      progn, F_SurfaceFileName);
      AppExit(1);
    }

            iD = (Nz *iy +iz) + (Ny * Nz * ix);
    Surface[iD] = (Fprec) val;
    /* Fprintf(stdout, " iD  %4d  Surface[iD]  %6.2f \n", iD,  Surface[iD]); */
  }
 /*  c = '8'; */

  (void) fclose (fpdens);

  Fprintf(stdout, "# %s xyz  %d %d %d\n\n",
    F_SurfaceFileName, NFile[1], NFile[2], NFile[3]);

  FM = NULL;

  if (FlagMap == NULL) {
    CheckMalloc(FM, Nx * Ny * Nz);
    DS_MarkEquiv(&SpgrInfo, FM, Nx, Ny, Nz);
    FlagMap = FM;
  }

  MapInd(Nx * Ny * Nz, Surface, FlagMap, 1, 1, &CorrCoeff);

  Fprintf(stdout, "# Symmetry enforced on Surface map.\n");
  Fprintf(stdout, "# Correlation of symmetry dependent grid points\n");
  Fprintf(stdout, "# in input map was %.6g\n\n", CorrCoeff);

  NormalizeSurface(Surface, Nx, Ny, Nz);

  MapInd(Nx * Ny * Nz, Surface, FlagMap, 1, 1, &CorrCoeff);

  if (CorrCoeff < 0.9999)
    InternalError("Corrupt Surface map after NormalizeSurface()");

  if (FM) AppFree(FM, Nx * Ny * Nz);
}
