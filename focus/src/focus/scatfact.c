#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>


#include "main.h"
#include "lib.h"
#include "xtal.h"


static Fprec
      Flag_SF_Tab_SortFunction = -1.;

static int SF_Tab_SortFunction(T_SF_Tab *a, T_SF_Tab *b)
{
  if (a->stol < b->stol) return -1;
  if (a->stol > b->stol) return  1;

  if (Flag_SF_Tab_SortFunction == -1. && a->f != b->f)
      Flag_SF_Tab_SortFunction = a->stol;

  return 0;
}


void SortSF_Tables(void)
{
  int  i;


  Flag_SF_Tab_SortFunction = -1.;

  for (i = 0; i < nSF_Tables; i++)
  {
    if (SF_Tables[i].nTab > 1)
    {
      qsort((void *)  SF_Tables[i].Tab, SF_Tables[i].nTab,
            sizeof (*(SF_Tables[i].Tab)),
            (SortFunction) SF_Tab_SortFunction);

      if (Flag_SF_Tab_SortFunction != -1.)
      {
        char  buf[128];

        Sprintf(buf,
          "ScatteringFactor %.60s: duplicate sin(Theta)/Lambda %.6g",
          SF_Tables[i].Label, Flag_SF_Tab_SortFunction);

        progerror(buf);
      }
    }
  }
}


int CompleteSF_Info(T_SF_Info *SFI, int LblExact, int Set_f_stol_0)
{
  int          m;
  const char   *lt, *li;
  int          iSFT;
  T_SF_Tables  *SFT;


  if (SFI->SFT || SFI->CAA || SFI->f_stol_0 != -1.)
    InternalError("CompleteSF_Info(): SFI already set");

  SFT = SF_Tables;

  for (iSFT = 0; iSFT < nSF_Tables; iSFT++, SFT++)
  {
    if (str_icmp(SFT->Label, SFI->Lbl) == 0)
    {
      SFI->SFT = SFT;
      break;
    }
  }

  if (SFI->SFT == NULL && LblExact == 0)
  {
    SFT = SF_Tables;

    for (iSFT = 0; iSFT < nSF_Tables; iSFT++, SFT++)
    {
      m = 0;

      for (lt = SFT->Label, li = SFI->Lbl; isalpha(*lt); lt++, li++)
      {
        if (toupper(*lt) != toupper(*li))
        {
          m = 0;
          break;
        }

        m++;
      }

      if (m)
      {
        SFI->SFT = SFT;
        break;
      }
    }
  }

  if (SFI->SFT)
  {
    if (Set_f_stol_0)
    {
      if (SFI->SFT->Tab[0].stol != 0.)
      {
        char  buf[128];

        Sprintf(buf,
          "ScatteringFactor %.60s: missing entry for sinTovL = 0",
          SFI->SFT->Label);

        progerror(buf);
      }

      SFI->f_stol_0 = SFI->SFT->Tab[0].f;
    }
  }
  else
  {
    if (ModeScatteringFactorTable == 0)
    {
      SFI->CAA = FindSF_WK95_CAA(SFI->Lbl, LblExact);
      // printf("atom: %s -> using xray scattering factors from WK95\n", SFI->Lbl);
      if (NULL == SFI->CAA)
        return -1;
  
      if (Set_f_stol_0)
        SFI->f_stol_0 = CalcSF_WK95_CAA(SFI->CAA, 0.);
    }
 //   else if (ModeScatteringFactorTable == 1)
 //   {
 //     SFI->CAA = FindSF_IT92_CAA(SFI->Lbl, LblExact);
 //     printf("%s -> using scattering factors from IT92\n", SFI->Lbl);
 //     if (NULL == SFI->CAA)
 //       return -1;
 // 
 //     if (Set_f_stol_0)
 //       SFI->f_stol_0 = CalcSF_IT92_CAA(SFI->CAA, 0.);
 //   }
    else if (ModeScatteringFactorTable == 2)
    {
      SFI->CAA = FindSF_IT4322(SFI->Lbl, LblExact); // StefS 
      // printf("atom: %s -> using electron scattering factors from IT 4.3.2.2\n", SFI->Lbl); 
      if (NULL == SFI->CAA)
        return -1;
  
      if (Set_f_stol_0)
        SFI->f_stol_0 = CalcSF_IT4322(SFI->CAA, 0.);
    }
    else if (ModeScatteringFactorTable == 3)
    {
      SFI->CAA = FindSF_IT4323(SFI->Lbl, LblExact); // StefS 
      // printf("atom: %s -> using electron scattering factors from IT 4.3.2.3\n", SFI->Lbl); 
      if (NULL == SFI->CAA)
        return -1;
  
      if (Set_f_stol_0)
        SFI->f_stol_0 = CalcSF_IT4323(SFI->CAA, 0.);
    }
  }

  return 0;
}


static Fprec CalcSF_SFT(T_SF_Tables *SFT, Fprec stol)
{
  int       hit;
  Fprec     f0;
  int       iTab;
  T_SF_Tab  *Tab;


  Tab = SFT->Tab;

  hit = 0;

  for (iTab = 0; iTab < SFT->nTab; iTab++, Tab++)
  {
    if (stol <= Tab->stol)
    {
      hit = 1;
      break;
    }
  }

  if (hit == 0 || stol < SFT->Tab[0].stol)
  {
    char  buf[128];

    Sprintf(buf,
      "ScatteringFactor %.60s: sinTovL = %.6g out of range",
      SFT->Label, stol);

    progerror(buf);
  }

  if (iTab == 0 || stol == Tab->stol)
    f0 = Tab->f;
  else
  {      Tab--;
    f0 = Tab[0].f + (Tab[1].f - Tab[0].f) * (Tab[1].stol -        stol)
                                          / (Tab[1].stol - Tab[0].stol);
  }

  return f0;
}


Fprec CalcFormFactor(T_SF_Info *SFI, int h, int k, int l,
                     Fprec Occ, Fprec Uiso)
{
  double  stol2, f0, B;
  Fprec   FF;


  stol2 = .25 * Q_hkl(h, k, l, &LatConR);

  if (SFI->SFT)
    f0 = CalcSF_SFT(SFI->SFT, sqrt(stol2));
  else if (SFI->CAA)
  {
    if (ModeScatteringFactorTable == 0) 
      f0 = CalcSF_WK95_CAA(SFI->CAA, stol2);
 //   else if (ModeScatteringFactorTable == 1) 
 //     f0 = CalcSF_IT92(SFI->CAA, stol2);
    else if (ModeScatteringFactorTable == 2) 
      f0 = CalcSF_IT4322(SFI->CAA, stol2);
    else if (ModeScatteringFactorTable == 3) 
      f0 = CalcSF_IT4323(SFI->CAA, stol2);
  }
  else
    InternalError("CalcFormFactor(): SFI SFT and CAA not set");

  B = 8. * M_PI * M_PI * Uiso;
  FF = f0 * exp(-B * stol2);
  FF *= Occ;

  if (Debug0) /* Debug: Print h k l sqrt(stol2) Occ Uiso f0 B FF */
  {
    Fprintf(stdout,
      " %3d %3d %3d  %8.4f %8.4f %8.4f %8.4f %8.4f  %8.4f\n",
      h, k, l, sqrt(stol2), Occ, Uiso, f0, B, FF);
  }

  return FF;
}
