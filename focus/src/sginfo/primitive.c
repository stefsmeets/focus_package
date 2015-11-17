static int PrimitiveLatticeConstants(T_LatticeConstants *LatConA,
                                     T_LatticeConstants *LatConB,
                                     T_SgInfo *SgInfo)
{
  int     i, j;
  int     CCMx_PL[9], deterCCMx_LP;
  double  GA[9], GB[9], GAR[9], R[9], Rt[9];


  LatConA->calcs = 1;
  LatConA->calcc = 1;

  /* just to check LatConA and to compute sin and cos of angles
   */
  if (Lc2RLc(LatConA, LatConB) != 0) {
    SetSgError("Error: Illegal UnitCell");
    return -1;
  }

  Lc2MetricalMx(LatConA, GA);

  deterCCMx_LP = deterRotMx(SgInfo->CCMx_LP);
              iCoFactorMxTp(SgInfo->CCMx_LP, CCMx_PL);

  if (deterCCMx_LP < 1)
    goto ReturnError;

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) {
       R[i * 3 + j] = CCMx_PL[i * 3 + j] / (double) deterCCMx_LP;
      Rt[i * 3 + j] = CCMx_PL[j * 3 + i] / (double) deterCCMx_LP;
    }

  MxMultiply(GAR, GA, R, 3, 3, 3);
  MxMultiply(GB, Rt, GAR, 3, 3, 3);

  if (GB[0] < 0. || GB[4] < 0. || GB[8] < 0.)
    goto ReturnError;

  LatConB->a = sqrt(GB[0]);
  LatConB->b = sqrt(GB[4]);
  LatConB->c = sqrt(GB[8]);

  LatConB->alpha = GB[5] / LatConB->b / LatConB->c;
  LatConB->beta  = GB[2] / LatConB->c / LatConB->a;
  LatConB->gamma = GB[1] / LatConB->a / LatConB->b;

  if (   LatConB->alpha < -1. || LatConB->alpha > 1.
      || LatConB->beta  < -1. || LatConB->beta  > 1.
      || LatConB->gamma < -1. || LatConB->gamma > 1.)
    goto ReturnError;

  LatConB->alpha = acos(LatConB->alpha);
  LatConB->beta  = acos(LatConB->beta );
  LatConB->gamma = acos(LatConB->gamma);

  LatConB->calcs = 1;
  LatConB->calcc = 1;

  return 0;

  ReturnError:

  SetSgError("InternalError: PrimitiveLatticeConstants()");
  return -1;
}

  if (F_UnitCell)
  {
    putc('\n', stdout);

    if (HarmonizeSgLatCon(&SpgrInfo[0], &LatConA, F_UnitCell) != 0)
      PrintClearSgError(0, 1);

    if (PrimitiveLatticeConstants(&LatConA,
                                  &LatConB, &SpgrInfo[0]) != 0)
    PrintClearSgError(0, 1);

    fprintf(stdout,
      "L UnitCell  %.6g %.6g %.6g %.6g %.6g %.6g\n",
      LatConA.a, LatConA.b, LatConA.c,
      LatConA.alpha / PIover180,
      LatConA.beta  / PIover180,
      LatConA.gamma / PIover180);

    fprintf(stdout,
      "P UnitCell  %.6g %.6g %.6g %.6g %.6g %.6g\n",
      LatConB.a, LatConB.b, LatConB.c,
      LatConB.alpha / PIover180,
      LatConB.beta  / PIover180,
      LatConB.gamma / PIover180);
  }
