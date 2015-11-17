/*
  Space Group Info (c) 1994-97 Ralf W. Grosse-Kunstleve
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>


#include "sginfo.h"


#define Fprintf (void) fprintf
#define Sprintf (void) sprintf


void ListReflCond(FILE *fpout, T_hklCond *cnd, int cndn, int debug)
{
  int   i,j,n,m;
  char  hkl[3]="hkl";
  int   x[3];
  int   S[3][3];
  int   first, HaveCondition;
  char  buf[128];


  Fprintf(fpout, "Reflection conditions\n");

  HaveCondition = 0;

  for(m=0;m<cndn;m++,cnd++) {
    if (!debug && !cnd->Exist ) continue;
    if (debug) Sprintf(buf,"%02d%c ", cnd->iSymTr, cnd->Exist?' ':'-');
    else buf[0]=0;
    (void) memcpy(S,cnd->Mx,sizeof(S));
    (void) memcpy(x,cnd->Gl,sizeof(x));

    for(j=0; j<3; j++){
      for(n=i=0; i<3; i++) if((n=S[i][j]) != 0) break;
      if(i==3)
        Sprintf(&buf[strlen(buf)],"0");
      else {
        if (n==1) Sprintf(&buf[strlen(buf)],"%1c",hkl[i]);
        else if (n==-1) Sprintf(&buf[strlen(buf)],"-%1c",hkl[i]);
        else Sprintf(&buf[strlen(buf)],"%d%1c",n,hkl[i]);
      }
    }

    Sprintf(&buf[strlen(buf)],": ");
    n = cnd->Denom;
    for(i=0,j=0; i<3; i++) j+=abs(x[i] % n);
    if (j) {
      if(!cnd->Exist) Sprintf(&buf[strlen(buf)],"(");
      for(i=0,first=1; i<3; i++) if(x[i] % n) {
        if(x[i] <0) Sprintf(&buf[strlen(buf)],"-");
        else if (!first) Sprintf(&buf[strlen(buf)],"+");
        if(abs(x[i]) != 1) Sprintf(&buf[strlen(buf)],"%1d",abs(x[i]));
        first=0;

        Sprintf(&buf[strlen(buf)],"%1c",hkl[i]);
      }
      Sprintf(&buf[strlen(buf)],"=%dn",n);
      if(!cnd->Exist) Sprintf(&buf[strlen(buf)],")");
    } else {
      Sprintf(&buf[strlen(buf)],"No conditions");
    }

    if (strlen(buf) != 0) {
      Fprintf(fpout, "  %s\n", buf);
      HaveCondition = 1;
    }
  }

  if (HaveCondition == 0)
    Fprintf(fpout, "  hkl: No Condition\n");

  putc('\n', fpout);
}


void ListSysEnhanced(FILE *fpout,
                     const T_hklCond *cnd, const int *SysEnhanced, int cndn)
{
  int   i,j,n,m;
  char  hkl[3]="hkl";
  int   x[3];
  int   S[3][3];
  char  buf[128];


  Fprintf(fpout, "Systematically enhanced reflections\n");

  for(m=0;m<cndn;m++,cnd++) {
    if(SysEnhanced[m]==0) continue;
    (void) memcpy(S,cnd->Mx,sizeof(S));
    (void) memcpy(x,cnd->Gl,sizeof(x));

    buf[0] = 0;
    for(j=0; j<3; j++){
      for(n=i=0; i<3; i++) if((n=S[i][j]) != 0) break;
      if(i==3)
        Sprintf(&buf[strlen(buf)],"0");
      else {
        if (n==1) Sprintf(&buf[strlen(buf)],"%1c",hkl[i]);
        else if (n==-1) Sprintf(&buf[strlen(buf)],"-%1c",hkl[i]);
        else Sprintf(&buf[strlen(buf)],"%d%1c",n,hkl[i]);
      }
    }

    Sprintf(&buf[strlen(buf)],": e=%d",SysEnhanced[m]);

    Fprintf(fpout, "  %s\n", buf);
  }

  putc('\n', fpout);
}


void ListRestCond(FILE *fpout, T_hklCond *rst, int rstn, int debug)
{
  int   i,j,n,m;
  char  hkl[3]="hkl";
  int   x[3];
  int   S[3][3];
  int   first;
  char  buf[128];


  Fprintf(fpout, "Reflections with phase restriction\n");

  for(m=0;m<rstn;m++,rst++) {
    if (!debug && !rst->Exist ) continue;
    if (debug) Sprintf(buf,"%02d%c ", rst->iSymTr, rst->Exist?' ':'-');
    else buf[0]=0;
    (void) memcpy(S,rst->Mx,sizeof(S));
    (void) memcpy(x,rst->Gl,sizeof(x));

    for(j=0; j<3; j++){
      for(n=i=0; i<3; i++) if((n=S[i][j]) != 0) break;
      if(i==3)
        Sprintf(&buf[strlen(buf)],"0");
      else {
        if (n==1) Sprintf(&buf[strlen(buf)],"%1c",hkl[i]);
        else if (n==-1) Sprintf(&buf[strlen(buf)],"-%1c",hkl[i]);
        else Sprintf(&buf[strlen(buf)],"%d%1c",n,hkl[i]);
      }
    }

    n = rst->Denom;
    Sprintf(&buf[strlen(buf)],": ");
    if(!rst->Exist) Sprintf(&buf[strlen(buf)],"<");
    for(i=0,j=0; i<3; i++) j+=abs(x[i] % n);
    if (j) {
      Sprintf(&buf[strlen(buf)],"%d deg (",180/n);
      for(i=0,first=1; i<3; i++) if(x[i] % n) {
        if(x[i] <0) Sprintf(&buf[strlen(buf)],"-");
        else if (!first) Sprintf(&buf[strlen(buf)],"+");
        if(abs(x[i]) != 1) Sprintf(&buf[strlen(buf)],"%1d",abs(x[i]));
        first=0;

        Sprintf(&buf[strlen(buf)],"%1c",hkl[i]);
      }
      Sprintf(&buf[strlen(buf)],")");
    } else {
      Sprintf(&buf[strlen(buf)],"0 deg");
    }
    if(!rst->Exist) Sprintf(&buf[strlen(buf)],">");

    Fprintf(fpout, "  %s\n", buf);
  }

  putc('\n', fpout);
}
