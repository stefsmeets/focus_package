/*
  The source code in this file was contributed by:

    Maxim Larine and Slava Klimkovich
    SoftHard Technology, Ltd.
    E-mail: max@softhard.sk

  Space Group Info (c) 1994-97 Ralf W. Grosse-Kunstleve
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>


#include "sginfo.h"


#define PreCook

static void CookCond(const int *Mx, const int *Gl, int *cMx, int *cGl)
{
  int  i, j, k, p;


  for (i = 0; i < 9; i++) cMx[i] = Mx[i];

  if (cGl != NULL)
    for (i = 0; i < 3; i++) cGl[i] = Gl[i];

  p = 0;

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
      if (cMx[i * 3 + j] != 0) break;

    if (j < 3) continue;

    for (j = 0; j < 3; j++) {
          p = cMx[j * 3 + i];
      if (p == 1 || p == -1) break;
    }

    if (j == 3) continue;

    if (j < i && cMx[j * 3 + j] == 1) continue;

    for (k = 0; k < 3; k++)
      if (k != i && k != j) break;

    if (cMx[j * 3 + k] != 0) continue;

    if (cGl != NULL && (cGl[i] || cGl[j])) continue;

    for (k = 0; k < 3; k++)
      cMx[i * 3 + k] = p * cMx[j * 3 + k];

    for (k = 0; k < 3; k++)
      cMx[j * 3 + k] = 0;

    break;
  }

  for (i = 0; i < 3; i++)
  {
    if (cMx[i * 3 + i] != -1) continue;

    for (j = 0; j < 3; j++)
      cMx[i * 3 + j] *= -1;
  }

  return;
}


#define  AHKL  7
typedef  int  T_AA[2 * AHKL + 1][2 * AHKL + 1][2 * AHKL + 1];


static int next(int n)
{
  return (n+1) % 3;
}


static int primV(int v[3])
{
  int i,n,s;
  for(n=i=1; i<=STBF; i++) if( ((v[0]%i) | (v[1]%i) | (v[2]%i)) ==0) n=i;
  for(s=i=0; i<3; i++) if(v[i]>0) s++; else if(v[i]<0) s--;
  if(s<0) n*= -1;
  else if(s==0) {
    if(v[0]<0) n*=-1;
    else if (v[0]==0 && v[1]<0) n*=-1;
  }
  for(i=0; i<3; i++) v[i] /= n;
  return abs(n);
}


static int primVns(int v[3]) /* Does not change the sign */
{
  int i,n;
  for(n=i=1; i<=STBF; i++) if( ((v[0]%i) | (v[1]%i) | (v[2]%i)) ==0) n=i;
  for(i=0; i<3; i++) v[i] /= n;
  return abs(n);
}


#if defined(COMPILE_UNUSED)
static int primC(int v[3])
{
  int i,n,s;
  for(n=i=1; i<=STBF; i++) if( ((v[0]%i) | (v[3]%i) | (v[6]%i)) ==0) n=i;

  for(s=i=0; i<3; i++) if(v[i*3]>0) s++; else if(v[i*3]<0) s--;
  if(s<0) n*= -1;
  else if(s==0) {
    if(v[0]<0) n*=-1;
    else if (v[0]==0 && v[3]<0) n*=-1;
  }

  for(i=0; i<3; i++) v[i*3] /= n;
  return abs(n);
}
#endif


static int primM(int m[3][3])
{
  int i,n,s,j,x;
  int *v = m[0];
  for(n=i=1; i<=STBF; i++) {
    for(j=0,x=0;j<3*3;j++) x|=v[j]%i;
    if (x==0) n=i;
  }
  for(s=i=j=0; i<3*3; i++) {
    if(v[i]>0) s++;
    if(v[i]!=0 && j==0) j=v[i]/abs(v[i]);
    if(v[i]<0) s--;
  }
  if (j==0) return 1;
  if(s<0) n*= -1;
  else if(s==0) n*=j;
  for(i=0; i<3*3; i++) v[i] /= n;
  return abs(n);
}


static int min1(int m[3][3], int v[3])
{
  int i,j,d;
  v[0]=v[1]=v[2]=0;
  for(i=0; i<3; i++) for(j=2; j>=0; j--) {
    d =  m[i][j];
    if(d) { v[0] = d; v[1] = i; v[2] = j; return 1;}
  }
  return 0;
}


static int min2(int m[3][3], int v[3])
{
  int t[3][2] = {{1,2},{0,2},{0,1}};
  int i,j,d;
  v[0]=v[1]=v[2]=0;
  for(j=0; j<3; j++) for(i=0; i<3; i++) {
    d =  m[t[i][0]][t[j][0]] * m[t[i][1]][t[j][1]] -
         m[t[i][1]][t[j][0]] * m[t[i][0]][t[j][1]];
/*  if(d) { v[0] = d; v[1] = i; v[2] = j; return 1;} */
/* We have to choose a minor with the maximum absoulte value of determinant */
    if(d && abs(d)>abs(v[0])) { v[0] = d; v[1] = i; v[2] = j; }
  }
  if (v[0]) return 1;
  return 0;
}


static int det(int m[3][3])
{
  int i, p, s, d=0;
  for(p=0; p<3; p++) {
    for(s=1, i=0; i<3; i++) s *= m[i][(i+p)%3];    d += s;
    for(s=1, i=0; i<3; i++) s *= m[2-i][(i+p)%3];  d -= s;
  }
  return d;
}


static int rm(int m[3][3], int v[3])
{
  v[0] = det(m);
  if(v[0]) return 3;
  if(min2(m,v)) return 2;
  if(min1(m,v)) return 1;
  return 0;
}


static int ConCmp(const T_hklCond *a, const T_hklCond *b)
{
  int              i, j;


  if (a->Rank < b->Rank) return -1;
  if (a->Rank > b->Rank) return  1;

  for (i = 0; i < 3; i++)
  for (j = 0; j < 3; j++)
  {
    if (a->Mx[i][j] < b->Mx[i][j]) return  1;
    if (a->Mx[i][j] > b->Mx[i][j]) return -1;
  }

  if (a->Denom < b->Denom) return  1;
  if (a->Denom > b->Denom) return -1;

  for (i = 0; i < 3; i++)
  {
    if (a->Gl[i] < b->Gl[i]) return  1;
    if (a->Gl[i] > b->Gl[i]) return -1;
  }

  return 0;
}


static int RestCmp(const T_hklCond *a, const T_hklCond *b)
{
  int              i, j;


  if (a->Rank < b->Rank) return -1;
  if (a->Rank > b->Rank) return  1;

  for (i = 0; i < 3; i++)
  for (j = 0; j < 3; j++)
  {
    if (a->Mx[i][j] < b->Mx[i][j]) return  1;
    if (a->Mx[i][j] > b->Mx[i][j]) return -1;
  }

  if (a->Denom < b->Denom) return -1;
  if (a->Denom > b->Denom) return  1;

  for (i = 0; i < 3; i++)
  {
    if (a->Gl[i] < b->Gl[i]) return  1;
    if (a->Gl[i] > b->Gl[i]) return -1;
  }

  return 0;
}


int CondIsSysAbsent_hkl(const T_hklCond *cond, int ncond, int h, int k, int l)
{
  int    i;
#ifdef PostCook
  int        Mx[9], Gl[3];
#else
  const int  *Mx, *Gl;
#endif


  for(i = 0; i < ncond; i++, cond++)
  {
    if (! cond->Exist)
      continue;

#ifdef PostCook
    CookCond(cond->Mx[0], cond->Gl, Mx, Gl);
#else
    Mx = cond->Mx[0];
    Gl = cond->Gl;
#endif

    if (   (h == Mx[0] * h + Mx[3] * k + Mx[6] * l)
        && (k == Mx[1] * h + Mx[4] * k + Mx[7] * l)
        && (l == Mx[2] * h + Mx[5] * k + Mx[8] * l)
        && (((h * Gl[0] + k * Gl[1] + l * Gl[2]) % cond->Denom) != 0))
      return 1;
  }

  return 0;
}


int Get_hklRestriction(const T_hklCond *cond, int ncond, int h, int k, int l)
{
  int    Restr, i;
#ifdef PostCook
  int        Mx[9], Gl[3];
#else
  const int  *Mx, *Gl;
#endif


  if (! (h || k || l))
    return 0;

  for(i = 0; i < ncond; i++, cond++)
  {
    if (! cond->Exist)
      continue;

#ifdef PostCook
    CookCond(cond->Mx[0], cond->Gl, Mx, Gl);
#else
    Mx = cond->Mx[0];
    Gl = cond->Gl;
#endif

    if (   (h == Mx[0] * h + Mx[3] * k + Mx[6] * l)
        && (k == Mx[1] * h + Mx[4] * k + Mx[7] * l)
        && (l == Mx[2] * h + Mx[5] * k + Mx[8] * l))
    {
      Restr = (h * Gl[0] + k * Gl[1] + l * Gl[2]) % cond->Denom;
      if (Restr < 0) Restr += cond->Denom;
      Restr *= STBF;
      if (Restr % cond->Denom) {
        SetSgError("Internal Error: Get_hklRestriction()");
        return -1;
      }
      Restr /= cond->Denom;
      return Restr;
    }
  }

  return -1;
}


int Get_hklEpsilon(const T_hklCond *cond, int ncond, int *SysEnhanced,
                   int h, int k, int l)
{
  int    i, Epsilon;
#ifdef PostCook
  int        Mx[9];
#else
  const int  *Mx;
#endif


  if (! (h || k || l))
    return 0;

  Epsilon = 0;

  for(i = 0; i < ncond; i++, cond++)
  {
    if (SysEnhanced[i] <= Epsilon)
      continue;

#ifdef PostCook
    CookCond(cond->Mx[0], NULL, Mx, NULL);
#else
    Mx = cond->Mx[0];
#endif

    if (   (h == Mx[0] * h + Mx[3] * k + Mx[6] * l)
        && (k == Mx[1] * h + Mx[4] * k + Mx[7] * l)
        && (l == Mx[2] * h + Mx[5] * k + Mx[8] * l))
      Epsilon = SysEnhanced[i];
  }

  return Epsilon;
}


static int IsInSpace(T_hklCond *cond, int h, int k, int l)
{
  int hm, km, lm;
  int   r=0;
  int *M=cond->Mx[0];

  hm = M[0] * h + M[3] * k + M[6] * l;
  km = M[1] * h + M[4] * k + M[7] * l;
  lm = M[2] * h + M[5] * k + M[8] * l;
  if ((h == hm) && (k == km) && (l == lm)) return 1; /* primitive case */
  if (h && hm) {
    if (hm%h) return 0;
    if (abs(hm/h)>abs(r)) r = hm/h;
  }
  if (k && km) {
    if (km%k) return 0;
    if (abs(km/k)>abs(r)) r = km/k;
  }
  if (l && lm) {
    if (lm%l) return 0;
    if (abs(lm/l)>abs(r)) r = lm/l;
  }
  if (abs(r)<1) return 0;
  if (hm%r || km%r || lm%r) return 0;
  if ((h == hm/r) && (k == km/r) && (l == lm/r)) return 1;
  return 0;
}


int SetRestCond(T_SgInfo *SgInfo)
{
  int     r,i,j,a,b,c,d,f,n,m,h,k,l;
  int      v[3],x[3], t[3]={6,6,6};
  int      R[3][3],S[3][3];
  int     nTrV, iTrV, nLoopInv, iLoopInv, iList, iSymTr, rstn;
  int      *X, good;
  int    *N;
  const int  *TrV;
  T_RTMx    SMx;
  T_RTMx          *lsmx;
  T_hklCond  *rst, *Rest;

  const T_hklCond  *Cond;
  int              CondNum;

  T_AA  *AA;


  SgInfo->nRestCond = -1;

  CondNum = SgInfo->nReflCond;

      Cond = SgInfo->ReflCond;
  if (Cond == NULL) {
    SetSgError("Internal Error: SgInfo->ReflCond == NULL");
    return -1;
  }

      Rest = SgInfo->RestCond;
  if (Rest == NULL) {
    SetSgError("Internal Error: SgInfo->RestCond == NULL");
    return -1;
  }

      AA = malloc(sizeof (*AA));
  if (AA == NULL) {
    SetSgError("Not enough core");
    return -1;
  }

  rst = Rest;

  nLoopInv = Sg_nLoopInv( SgInfo );
  nTrV = SgInfo->LatticeInfo->nTrVector;
   TrV = SgInfo->LatticeInfo->TrVector;
  iSymTr = 0;
  rstn = 0;
  for( iTrV=0; iTrV<nTrV; TrV+=3, iTrV++ ) {
    for ( iLoopInv=0; iLoopInv<nLoopInv; iLoopInv++ ) {
      if (iLoopInv == 0 ) f =  1;
      else                f = -1;
      lsmx = SgInfo->ListSeitzMx;

      for (iList = 0; iList < SgInfo->nList; iList++, lsmx++) {
  iSymTr++;
        for(i=0;i<9;i++)
          SMx.s.R[i] =               f*lsmx->s.R[i];
        for(i=0;i<3;i++)
          SMx.s.T[i] = iModPositive( f*lsmx->s.T[i]+TrV[i], STBF );

#if defined(COMPILE_UNUSED)
        for(i=0;i<3;i++) if (SMx.s.T[i]) break;
        if (i==3) continue;
#endif

        for(i=0; i<3; i++) for(j=0; j<3; j++)  R[j][i] = SMx.a[i*3+j];
        for(i=0; i<3; i++) t[i]=SMx.s.T[i];
        for(i=0; i<3; i++) if(t[i]>STBF/2) t[i]-=STBF;

        R[0][0]+=1;  R[1][1]+=1;  R[2][2]+=1;

        r = rm(R,v);

        (void) memset(S,0,sizeof (S));

        switch (r) {
          case 0:  S[0][0] = S[1][1] = S[2][2] = 1; break;
          case 1:
            i = v[1]; a = v[2]; b=next(a); c = next(b);
            S[b][a] = -R[i][b]; S[b][b]=R[i][a]; S[b][c]=0;
            S[c][a] = -R[i][c]; S[c][c]=R[i][a]; S[c][b]=0;
            break;
          case 2:
            i = next(v[1]); j = next(i);
            a = v[2]; b=next(a); c = next(b); d=v[0];
            S[a][a] = -d;
            S[a][b] = R[i][a]*R[j][c] - R[j][a]*R[i][c];
            S[a][c] = R[i][b]*R[j][a] - R[j][b]*R[i][a];
            break;

          case 3: continue;
        }
        (void) primM(S);
        if (r==2) for(i=0;i<3;i++) (void) primV(S[i]);
        else      for(i=0;i<3;i++) (void) primVns(S[i]);
#if defined(COMPILE_UNUSED)
        for(i=0;i<3;i++) primC(&S[0][i]);
#endif

        (void) memset(x,0,sizeof (x));
        for(i=0; i<3; i++) for(j=0; j<3; j++) x[i] += S[i][j]*t[j];
        for(i=0; i<3; i++) { x[i] %= STBF; if(x[i]>STBF/2) x[i]-=STBF; }
        n = STBF / primVns(x);
        for(i=0; i<3; i++) {
          if (2*abs(x[i])>n) {
            if (x[i]>0) x[i]-=n;
            else x[i]+=n;
          }
        }

        for(i=0,j=0; i<3; i++) j+=abs(x[i] % n);
#if defined(COMPILE_UNUSED)
          if (j) {
#endif
          rst->Exist = 1;
          rst->iSymTr = iSymTr;
          rst->Rank = r;
#ifdef PreCook
          CookCond(&S[0][0], x, &rst->Mx[0][0], rst->Gl);
#else
          memcpy(rst->Mx,S,sizeof (S));
          memcpy(rst->Gl,x,sizeof (x));
#endif
          rst->Denom = n;
          rst++; rstn++;
#if defined(COMPILE_UNUSED)
          }
#endif
      }
    }
  }

/* Now sorting Cond */

  qsort((void *) Rest, rstn, sizeof (*Rest),
        (int (*)(const void *, const void *)) RestCmp);

  for(m=0,n=0;m<rstn;m++) {
#if defined(DelTrivCond)
/* To delete trivials: */
    if(RestCmp(&Rest[n], &Rest[m]) != 0) Rest[++n]=Rest[m];
#else
/* To mark trivials: */
    if(RestCmp(&Rest[n], &Rest[m]) != 0) {
      n=m;
    } else {
      if (m!=n) Rest[m].Exist=0;
    }
#endif
  }

  j=0;
  (void) memset(AA, 0, sizeof (*AA));
  for(m=0;m<rstn;m++) {
    X = Rest[m].Gl;
    for(h=-AHKL;h<AHKL;h++) for(k=-AHKL;k<AHKL;k++) for(l=-AHKL;l<AHKL;l++) {
      if(h==0&&l==0&&k==0) continue;
      if (h==1&&k==-2&&l==0 && m==2)
        j++;
      if (IsInSpace(&Rest[m], h,k,l) &&
          ! CondIsSysAbsent_hkl( Cond, CondNum, h,k,l)) {
        N = &(*AA)[h+AHKL][k+AHKL][l+AHKL];
        i = h*X[0]+k*X[1]+l*X[2];
        i %= Rest[m].Denom;
        if (i<0) i+=Rest[m].Denom;
        i = (i*STBF)/Rest[m].Denom;
        if (i==0) i=STBF;
        if (*N==0) {
          *N = i;
        } else {
          if(*N != i ) {
            if (j==0) {
              SetSgError("Internal Error: CreateRestrictConditions()");
              free(AA);
              return -1;
            }
            j++;
          }
        }
      }
    }
  }

  for(m=0;m<rstn;m++) {
    X = Rest[m].Gl;
    good = 0;
    for(h=-AHKL;h<AHKL;h++) for(k=-AHKL;k<AHKL;k++) for(l=-AHKL;l<AHKL;l++) {
      if(h==0&&l==0&&k==0) continue;
      if (IsInSpace( &Rest[m], h,k,l ) &&
        ! CondIsSysAbsent_hkl( Cond, CondNum, h,k,l)) {
        N = &(*AA)[h+AHKL][k+AHKL][l+AHKL];
        good|=*N;
        *N=0;
      }
    }
    if (good==0) Rest[m].Exist = 0;
  }

  free(AA);

         SgInfo->nRestCond = rstn;
  return SgInfo->nRestCond;
}


int SetReflCond(T_SgInfo *SgInfo)
{
  int     r,i,j,a,b,c,d,f,n,m,h,k,l;
  int      v[3],x[3], t[3]={6,6,6};
  int      R[3][3],S[3][3];
  int     nTrV, iTrV, nLoopInv, iLoopInv, iList, iSymTr, cndn;
  int      *X, good;
  const int  *TrV;
  T_RTMx    SMx;
  T_RTMx     *lsmx;
  T_hklCond  *cnd, *Cond;

  T_AA  *AA;


  SgInfo->nReflCond = -1;

      Cond = SgInfo->ReflCond;
  if (Cond == NULL) {
    SetSgError("Internal Error: SgInfo->ReflCond == NULL");
    return -1;
  }

      AA = malloc(sizeof (*AA));
  if (AA == NULL) {
    SetSgError("Not enough core");
    return -1;
  }

  cnd = Cond;

  nLoopInv = Sg_nLoopInv( SgInfo );
  nTrV = SgInfo->LatticeInfo->nTrVector;
   TrV = SgInfo->LatticeInfo->TrVector;
  iSymTr = 0;
  cndn = 0;
  for( iTrV=0; iTrV<nTrV; TrV+=3, iTrV++ ) {
    for ( iLoopInv=0; iLoopInv<nLoopInv; iLoopInv++ ) {
      if (iLoopInv == 0 ) f =  1;
      else                f = -1;
      lsmx = SgInfo->ListSeitzMx;

      for (iList = 0; iList < SgInfo->nList; iList++, lsmx++) {
  iSymTr++;
        for(i=0;i<9;i++)
          SMx.s.R[i] =               f*lsmx->s.R[i];
        for(i=0;i<3;i++)
          SMx.s.T[i] = iModPositive( f*lsmx->s.T[i]+TrV[i], STBF );

        for(i=0; i<3; i++) for(j=0; j<3; j++)  R[j][i] = SMx.a[i*3+j];
        for(i=0; i<3; i++) t[i]=SMx.s.T[i];
        for(i=0; i<3; i++) if(t[i]>STBF/2) t[i]-=STBF;

        R[0][0]-=1;  R[1][1]-=1;  R[2][2]-=1;

        r = rm(R,v);

        (void) memset(S,0,sizeof (S));

        switch (r) {
          case 0:  S[0][0] = S[1][1] = S[2][2] = 1; break;
          case 1:
            i = v[1]; a = v[2]; b=next(a); c = next(b);
            S[b][a] = -R[i][b]; S[b][b]=R[i][a]; S[b][c]=0;
            S[c][a] = -R[i][c]; S[c][c]=R[i][a]; S[c][b]=0;
            break;
          case 2:
            i = next(v[1]); j = next(i);
            a = v[2]; b=next(a); c = next(b); d=v[0];
            S[a][a] = -d;
            S[a][b] = R[i][a]*R[j][c] - R[j][a]*R[i][c];
            S[a][c] = R[i][b]*R[j][a] - R[j][b]*R[i][a];
            break;

          case 3: continue;
        }
        for(i=0;i<3;i++) (void) primV(S[i]);

        (void) memset(x,0,sizeof (x));
        for(i=0; i<3; i++) for(j=0; j<3; j++) x[i] += S[i][j]*t[j];
        for(i=0; i<3; i++) { x[i] %= STBF; if(x[i]>STBF/2) x[i]-=STBF; }
        n = STBF / primV(x);
        for(i=0; i<3; i++) {
          if (2*abs(x[i])>n) {
            if (x[i]>0) x[i]-=n;
            else x[i]+=n;
          }
        }
        (void) primV(x);

        for(i=0,j=0; i<3; i++) j+=abs(x[i] % n);
#if defined(COMPILE_UNUSED)
          if (j) {
#endif
          cnd->Exist = 1;
          cnd->iSymTr = iSymTr;
          cnd->Rank = r;
#ifdef PreCook
          CookCond(&S[0][0], x, &cnd->Mx[0][0], cnd->Gl);
#else
          memcpy(cnd->Mx,S,sizeof (S));
          memcpy(cnd->Gl,x,sizeof (x));
#endif
          cnd->Denom = n;
          cnd++; cndn++;
#if defined(COMPILE_UNUSED)
          }
#endif
      }
    }
  }

/* Now sorting Cond */

  qsort((void *) Cond, cndn, sizeof (*Cond),
        (int (*)(const void *, const void *)) ConCmp);

  for(m=0,n=0;m<cndn;m++) {
#if defined(DelTrivCond)
/* To delete trivials: */
      if(ConCmp(&Cond[n], &Cond[m]) != 0) Cond[++n]=Cond[m];
#else
/* To mark trivials: */
    if(ConCmp(&Cond[n], &Cond[m]) != 0) {
      n=m;
    } else {
      if (m!=n) Cond[m].Exist=0;
    }
#endif
  }

  /* Reducing spaces for conditions limiting possible regflections */
  (void) memset(AA, 0, sizeof (*AA));
  for(m=0;m<cndn;m++) {
    if (!Cond[m].Exist) continue;
    X = Cond[m].Gl;
    good=0;
    for(h=-AHKL;h<AHKL;h++) for(k=-AHKL;k<AHKL;k++) for(l=-AHKL;l<AHKL;l++) {
      if (IsInSpace( &Cond[m], h,k,l ) &&
          (((h*X[0]+k*X[1]+l*X[2])%Cond[m].Denom) != 0)) {
        if((*AA)[h+AHKL][k+AHKL][l+AHKL]==0) {
           (*AA)[h+AHKL][k+AHKL][l+AHKL]=1;
          good=1;
        }
      }
    }
    if (good==0) Cond[m].Exist=0;
  }

  free(AA);

         SgInfo->nReflCond = cndn;
  return SgInfo->nReflCond;
}


int SetSysEnhanced(T_SgInfo *SgInfo)
{
  int   i,m,h,k,l;
  int   *M, good, M0[9];
  T_AA  *AA;


  if (SgInfo->nReflCond < 0) {
    SetSgError("Internal Error: SgInfo->ReflCond not set");
    return -1;
  }

  if (SgInfo->SysEnhanced == NULL) {
    SetSgError("Internal Error: SgInfo->SysEnhanced == NULL");
    return -1;
  }

      AA = malloc(sizeof (*AA));
  if (AA == NULL) {
    SetSgError("Not enough core");
    return -1;
  }

  /* Epsilon calculations */
  (void) memset(AA, 0, sizeof (*AA));
  for(m=0;m<SgInfo->nReflCond;m++) {
    M = SgInfo->ReflCond[m].Mx[0];
    for(h=-AHKL;h<AHKL;h++) for(k=-AHKL;k<AHKL;k++) for(l=-AHKL;l<AHKL;l++) {
      if (IsInSpace( &SgInfo->ReflCond[m], h,k,l ))
        (*AA)[h+AHKL][k+AHKL][l+AHKL]++;
    }
  }
  (*AA)[AHKL][AHKL][AHKL]=0;  /* i.e. 000 */

  (void) memset(M0,-1,sizeof (M0));
  for(m=SgInfo->nReflCond-1;m>=0;m--) {
    M = SgInfo->ReflCond[m].Mx[0];
    if(M[0]!=M0[0] || M[1]!=M0[1] || M[2]!=M0[2] ||
       M[3]!=M0[3] || M[4]!=M0[4] || M[5]!=M0[5] ||
       M[6]!=M0[6] || M[7]!=M0[7] || M[8]!=M0[8]) {

      for(i=0;i<9;i++) M0[i]=M[i];
      good=0;
      for(h=-AHKL;h<AHKL;h++) for(k=-AHKL;k<AHKL;k++) for(l=-AHKL;l<AHKL;l++) {
        if (IsInSpace( &SgInfo->ReflCond[m], h,k,l ) ) {
          if (good < (*AA)[h+AHKL][k+AHKL][l+AHKL])
              good = (*AA)[h+AHKL][k+AHKL][l+AHKL];

          (*AA)[h+AHKL][k+AHKL][l+AHKL]=0;
        }
      }
      SgInfo->SysEnhanced[m] = good;
    }
  }

  free(AA);

  return 0;
}
