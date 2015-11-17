#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "main.h"
#include "xtal.h"


static Fprec V2Tetrahedron(Fprec a2, Fprec b2, Fprec c2,
                           Fprec p2, Fprec q2, Fprec r2)
{
  /* Ref.: Bronstein, Semendjajew; Taschenbuch der Mathematik
           Chapter 2.6.2.3 Polyeder
   */

  Fprec  V2;

  V2 = -r2 * b2 * a2;
  V2 -= q2 * c2 * a2;
  V2 += q2 * p2 * a2;
  V2 += q2 * b2 * a2;
  V2 += r2 * p2 * c2;
  V2 += r2 * b2 * q2;
  V2 += r2 * b2 * c2;
  V2 += r2 * p2 * a2;
  V2 += b2 * q2 * c2;
  V2 += p2 * c2 * a2;
  V2 -= p2 * b2 * c2;
  V2 += r2 * q2 * c2;
  V2 += p2 * b2 * q2;
  V2 += b2 * p2 * a2;
  V2 += r2 * c2 * a2;
  V2 -= r2 * p2 * q2;
  V2 -= r2 * c2 * c2;
  V2 -= p2 * p2 * a2;
  V2 -= b2 * b2 * q2;
  V2 -= r2 * r2 * c2;
  V2 -= b2 * q2 * q2;
  V2 -= p2 * a2 * a2;
  V2 *= (Fprec)(1. / 144.);

  return V2;
}


void IdealTetrahedron(void)
{
  /* IdealTetrahedralAngle = acos(-1. / 3.); */

  if (IdealTetrahedralEdge2 < 0.)
      IdealTetrahedralEdge2 = 8. / 3. * IdealT_NodeDist2;

  if (IdealTetrahedralVol2 < 0.)
      IdealTetrahedralVol2 = V2Tetrahedron(IdealTetrahedralEdge2,
                                           IdealTetrahedralEdge2,
                                           IdealTetrahedralEdge2,
                                           IdealTetrahedralEdge2,
                                           IdealTetrahedralEdge2,
                                           IdealTetrahedralEdge2);
}


static int CheckOneTetrahedron(int iEPL,
                               T_PotentialNNBond **ActiveBonds,
                               int                nActiveBonds)
{
  int    iBond, jBond;
  int    count;
  Fprec  CosAngle[6], *CA, Vol2;


  if (Debug0) /* Debug: Print ActiveBonds */
  {
    for (iBond = 0; iBond < nActiveBonds; iBond++)
    {
      Fprintf(stdout, "AB  %8.4f %8.4f %8.4f # %8.4f %d %d\n",
                ActiveBonds[iBond]->Vect.x,
                ActiveBonds[iBond]->Vect.y,
                ActiveBonds[iBond]->Vect.z,
        AppSqrt(ActiveBonds[iBond]->Abs2),
        iEPL, eD_PeakList[iEPL].Index);
    }
  }

#define DotCV(v1_, v2_)\
  ((v1_).x * (v2_).x + (v1_).y * (v2_).y + (v1_).z * (v2_).z)

#define CrossCV(axb_, a_, b_)\
  {\
    (axb_).x = (a_).y * (b_).z - (b_).y * (a_).z;\
    (axb_).y = (a_).z * (b_).x - (b_).z * (a_).x;\
    (axb_).z = (a_).x * (b_).y - (b_).x * (a_).y;\
  }

  count = 0;
  CA = CosAngle;

  for (iBond = 0; iBond < nActiveBonds; iBond++)
  {
    for (jBond = iBond + 1; jBond < nActiveBonds; jBond++)
    {
      *CA  = DotCV(ActiveBonds[iBond]->Vect, ActiveBonds[jBond]->Vect);
      *CA /= AppSqrt(ActiveBonds[iBond]->Abs2 * ActiveBonds[jBond]->Abs2);

      if (*CA < (Fprec) COS_175_DEG) /* ARBITRARY */
        return -1;

      if (*CA > (Fprec) COS__75_DEG) /* ARBITRARY */
        if (++count > 2)
          return -1;

      if (Debug0) /* Print CosAngle */
        Fprintf(stdout, "TT_CA %d %d %d %d %8.4f\n",
          iEPL, eD_PeakList[iEPL].Index,
          iBond, jBond, acos(*CA) / PIover180);

      CA++;
    }
  }

#define TestAngles(i_, j_, k_)\
  {                         /* ARBITRARY */ count = 0;\
    if (CosAngle[i_] > (Fprec) COS__60_DEG) count++;\
    if (CosAngle[j_] > (Fprec) COS__60_DEG) count++;\
    if (CosAngle[k_] > (Fprec) COS__60_DEG) count++;\
    if (count > 1)\
      return -1;\
  }

  if (nActiveBonds == 4)
  {
    Fprec      a2, b2, c2, p2, q2, r2;
    T_fVector  av, bv, cv, pv, qv, rv;

    TestAngles(0, 1, 2);
    TestAngles(0, 3, 4);
    TestAngles(1, 3, 5);
    TestAngles(2, 4, 5);

#define DistDist2(i_, j_, vij_, v2_)\
    vij_.x = ActiveBonds[j_]->Vect.x - ActiveBonds[i_]->Vect.x;\
    vij_.y = ActiveBonds[j_]->Vect.y - ActiveBonds[i_]->Vect.y;\
    vij_.z = ActiveBonds[j_]->Vect.z - ActiveBonds[i_]->Vect.z;\
    v2_ = vij_.x * vij_.x + vij_.y * vij_.y + vij_.z * vij_.z;

    DistDist2(0, 1, av, a2);
    DistDist2(0, 2, bv, b2);
    DistDist2(0, 3, cv, c2);
    DistDist2(2, 3, pv, p2);
    DistDist2(3, 1, qv, q2);
    DistDist2(1, 2, rv, r2);

#undef DistDist2

        Vol2 = V2Tetrahedron(a2, b2, c2, p2, q2, r2);
    if (Vol2 < IdealTetrahedralVol2 * IdealTetrahedralVolFrac)
      return -1;

    /* idea: at least on of the bonds must be +- perpendicular to the
             plane of the other three vertices
       since tolerance is +-60 degree this is a weak filter
     */
    {
      Fprec      n2, nb;
      T_fVector  nv;

      CrossCV(nv, av, rv);
      n2 = DotCV(nv, nv);
      nb  = DotCV(nv, ActiveBonds[3]->Vect);
      nb /= AppSqrt(n2 * ActiveBonds[3]->Abs2);
      if (nb >= (Fprec) COS__60_DEG || nb <= (Fprec) COS_120_DEG)
        goto P_TestPassed;  /* ARBITRARY */

      CrossCV(nv, cv, qv);
      n2 = DotCV(nv, nv);
      nb  = DotCV(nv, ActiveBonds[2]->Vect);
      nb /= AppSqrt(n2 * ActiveBonds[2]->Abs2);
      if (nb >= (Fprec) COS__60_DEG || nb <= (Fprec) COS_120_DEG)
        goto P_TestPassed;  /* ARBITRARY */

      CrossCV(nv, bv, pv);
      n2 = DotCV(nv, nv);
      nb  = DotCV(nv, ActiveBonds[1]->Vect);
      nb /= AppSqrt(n2 * ActiveBonds[1]->Abs2);
      if (nb >= (Fprec) COS__60_DEG || nb <= (Fprec) COS_120_DEG)
        goto P_TestPassed;  /* ARBITRARY */

      CrossCV(nv, rv, qv);
      n2 = DotCV(nv, nv);
      nb  = DotCV(nv, ActiveBonds[0]->Vect);
      nb /= AppSqrt(n2 * ActiveBonds[0]->Abs2);
      if (nb >= (Fprec) COS__60_DEG || nb <= (Fprec) COS_120_DEG)
        goto P_TestPassed;  /* ARBITRARY */

      return -1;
    }

    P_TestPassed:;

    if (ModeCheckTetrahedra > 1)
    {
      Fprec      ne, nb;
      T_fVector  nv;

#define TestInside(vi_, vj_, vk_, vk2_, i_, s_)\
      CrossCV(nv, vi_, vj_);\
      ne  = s_ * DotCV(nv, vk_);\
      ne /= AppSqrt(vk2_);\
      nb  = -DotCV(nv, ActiveBonds[i_]->Vect);\
      nb /= AppSqrt(ActiveBonds[i_]->Abs2);\
      \
      if (   ne > CF0 && (nb <= CF0 || nb > ne)\
          || ne < CF0 && (nb >= CF0 || nb < ne))\
        return -1 /* outside */

      TestInside(av, bv, cv, c2, 0,  1);
      TestInside(av, cv, bv, b2, 0,  1);
      TestInside(bv, cv, av, a2, 0,  1);
      TestInside(qv, rv, av, a2, 1, -1);

#undef TestInside
    }
  }
  else if (nActiveBonds == 3)
    TestAngles(0, 1, 2);

#undef DotCV
#undef CrossCV
#undef TestAngles

  return 0;
}


int CheckTetrahedra(int Colored, const int *LnB, int Low_iEPL, int Top_iEPL,
                    int *CanNotBeTetrahedron)
{
  int                iEPL, jEPL, iLPNNB, iPNNB;
  int                nPNNB, nActiveBonds;
  T_PotentialNNBond  *PNNB, *ActiveBonds[4];
  int                nBad;


  nBad = 0;

  for (iEPL = Low_iEPL; iEPL <= Top_iEPL; iEPL++)
  {
    if (CanNotBeTetrahedron)
        CanNotBeTetrahedron[iEPL] = 0;

    if (LnB[iEPL] != 4) continue;

    nActiveBonds = 0;

    for (iLPNNB = 0; iLPNNB < eD_PeakList[iEPL].lLPNNB; iLPNNB++)
    {
              jEPL = eD_PeakList[iEPL].LPNNB[iLPNNB].iEPL;
      if (LnB[jEPL] == 0) continue;

      if (Colored == 0 || NodeColor(iEPL) != NodeColor(jEPL))
      {
        nPNNB = eD_PeakList[iEPL].LPNNB[iLPNNB].nPNNB;
         PNNB = eD_PeakList[iEPL].LPNNB[iLPNNB].PNNB;

        for (iPNNB = 0; iPNNB < nPNNB; iPNNB++, PNNB++)
        {
          if (nActiveBonds == 4)
            InternalError("Corrupt LnB");

          ActiveBonds[nActiveBonds] = PNNB;
                      nActiveBonds++;
        }
      }
    }

    if (nActiveBonds != 4)
      InternalError("Corrupt LnB");

    if (CheckOneTetrahedron(iEPL, ActiveBonds, nActiveBonds) != 0)
    {
      if (CanNotBeTetrahedron == NULL)
        return -1;

      CanNotBeTetrahedron[iEPL] = 1;
      nBad++;
    }
  }

  return nBad;
}
