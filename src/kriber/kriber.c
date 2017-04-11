/* Output from p2c, the Pascal-to-C translator */
/* From input file "kriber.pas" */

/* 94-04-01 Ralf Grosse Kunstleve, ETH Zurich
        all struct's are passed by reference now
 */

// Updated by Stef Smeets, Jul 2012
//      - Write cif files
//      - Allow display of atom labels > 99


#include "p2c.h"


/*#### BK ut */
/*#### BK loadatinput */
/*#### BK ew */
/******************************************************************************/
/*                                                                            */
/*   Program  KRIBER              Version 1.0       (January 1991)            */
/*   ===============                                                          */
/*                                                                            */
/*   Author:   Roland Bialek, Institut fuer Kristallografie und Petrografie   */
/*             ETH-Zentrum, CH-8092 Zuerich, Switzerland                      */
/*                                                                            */
/*   an interactive PASCAL program to                                         */
/*    - calculate distances and angles                                        */
/*    - generate input files for the programs DLS-76 and LOADAT (XRS-82)      */
/*    - calculate coordination sequences and loop configurations              */
/*                                                                            */
/*   Program Installation and Modification                                    */
/*   -------------------------------------                                    */
/*  The program is written in standard PASCAL and has been implemented and    */
/*  fully tested on a CYBER 855 running under NOS/VE 1.5.2. It has also been  */
/*  used successfully on a MicroVax II machine under VMS 5.2                  */
/*                                                                            */
/*  There is only one extension to the standard PASCAL used: The identifiers  */
/*  contain the underscore (_) and the dollar sign ($). If your PASCAL        */
/*  does not have this extension, you can remove these two signs with your    */
/*  editor.                                                                   */
/*                                                                            */
/*  Input:  The program runs interactively.                                   */
/*                                                                            */
/*  Output: The output is written to the terminal screen (longest line is     */
/*          80 columns wide).                                                 */
/*                                                                            */
/*  Files:  SYMDAT, STRUDAT, DISTDAT and COSEQ are input files, DLSINPUT,     */
/*          LOADATINPUT, LIST and STRUDATNEW are output files.                */
/*                                                                            */
/*  Arrays: Following constants can be changed in the constant declaration    */
/*          part:                                                             */
/*          maxanzahlatomlagen    :   maximum number of atom positions        */
/*          maxanzahlbindungen    :   maximum number of bonds per atom        */
/*          maxanzahlsymeqkarten  :   maximum number of symeq cards           */
/*          maxindepwinkel        :   maximum number of independent angles    */
/*          maxindepdistanzen     :   maximum number of independent distances */
/*          maxanzahlgegwinkel    :   maximum number of prescribed angles     */
/*          maxanzahlgegdistanzen :   maximum number of prescribed distances  */
/*          maxanzahldistanzen    :   maximum number of calculated distances  */
/*                                                                            */
/*  Modified by Bernt Karasch, modifications marked with '### BK'             */
/*  Maintained by Stef Smeets (2011-)                                                                          */
/*                                                                            */
/******************************************************************************/

#define maxanzahlatomlagen  500
#define maxanzahlbindungen  16
#define maxanzahlsymeqkarten  200 /* 200 is the limit of DLS-76 */
#define maxindepwinkel  2000
#define maxindepdistanzen  500
#define maxanzahlgegwinkel  100
#define maxanzahlgegdistanzen  100
#define maxanzahldistanzen  16


typedef double vektor[3];
typedef long gittervektor[3];
typedef double tensor[3][3];
typedef long rotationsmatrix[3][3];

typedef struct symmetrieoperator {
  rotationsmatrix r;
  vektor t;
} symmetrieoperator;

typedef symmetrieoperator st_symop_array[48];
typedef Char st_kristallsystem[4];
typedef Char st_raumgruppename[30];

typedef struct st_raumgruppe {
  long raumgruppenr;
  st_raumgruppename raumgruppename;
  st_kristallsystem kristallsystem;
  long zaehligkeit;
  symmetrieoperator *symop[192];
  long anz_rotmat, anztranslationen;
  vektor translation[4];
  boolean zentrosymm;
  long multtab[48][48];
} st_raumgruppe;

typedef Char abkbezeichnung[16];
typedef Char string80[80];

typedef struct zellparameter {
  double laenge[3], winkel[3];
} zellparameter;

typedef Char kt_atomlagename[6];

typedef struct kt_atom {
  struct kt_atomart *dessen_atomart;
  gittervektor elementarzelle;
} kt_atom;

typedef struct kt_atomart {
  struct kt_atomlage *deren_atomlage;
  long atomartnr;
  vektor pkoord;
  symmetrieoperator trilinpol;
  kt_atom bindung[maxanzahlbindungen];
} kt_atomart;

typedef Char kt_element[2];

typedef struct kt_atomlage {
  long atomlagenr;
  vektor atomlagekoord;
  symmetrieoperator *symmetrielage;
  kt_atomlagename aname;
  kt_element element;
  long koordinationszahl, anzatome;
  kt_atomart *atom[192];
  struct distanzen_tabellen_eintrag *dist_entry;
  kt_atom grundatom;
} kt_atomlage;

typedef struct kt_strukturdaten {
  string80 name, lit;
  abkbezeichnung abk;
  zellparameter zellpar;
  st_raumgruppename raumgruppe;
  st_kristallsystem kristallsystem;
  long anzatomlagen;
  kt_atomlage *atomlage[maxanzahlatomlagen];
} kt_strukturdaten;

typedef struct distanzen_tabellen_eintrag {
  double laenge;
  kt_atom endatom;
  struct distanzen_tabellen_eintrag *naechste_distanz;
} distanzen_tabellen_eintrag;

typedef long tt_koord_sequenz[10];
typedef enum {
  struktur, distanzen, bindungen
} berechenbare_dinge;
typedef Char mt_befehl[8];


Static _TEXT strudat, symdat;
/*#### BK ew */
Static _TEXT strudatn, distdat;
/*#### BK ut */
Static _TEXT dlsinp;
/*#### BK loadatinput */
Static _TEXT loadainp, list, coseq, cif;
Static mt_befehl befehl;
Static uchar berechnet;
Static boolean strudatnew_offen;
Static tensor st_metrik_tensor;
Static kt_strukturdaten strukturdaten;
Static st_raumgruppe rginfo;


static void InternalError(char *msg)
{
  fflush(stdout);
  fprintf(stderr, "kriber: Internal Error: %s\n", msg);
  exit(1);
}


Static double skalar_produkt(double *a, double *b)
{
  double x;
  long i, j;

  x = 0.0;
  for (i = 0; i <= 2; i++) {
    for (j = 0; j <= 2; j++)
      x += a[i] * st_metrik_tensor[i][j] * b[j];
  }
  return x;
}


Static Void load_zellparameter(double r1, double r2, double r3,
                               double w1, double w2, double w3)
{
  st_metrik_tensor[0][0] = r1 * r1;
  st_metrik_tensor[1][1] = r2 * r2;
  st_metrik_tensor[2][2] = r3 * r3;
  st_metrik_tensor[0][1] = r1 * r2 * cos(w3 * 0.0174533);
  st_metrik_tensor[1][2] = r2 * r3 * cos(w1 * 0.0174533);
  st_metrik_tensor[2][0] = r3 * r1 * cos(w2 * 0.0174533);
  st_metrik_tensor[1][0] = r1 * r2 * cos(w3 * 0.0174533);
  st_metrik_tensor[2][1] = r2 * r3 * cos(w1 * 0.0174533);
  st_metrik_tensor[0][2] = r3 * r1 * cos(w2 * 0.0174533);
}


Static double volumen(double a, double b, double c,
                      double aa, double bb, double cc)
{
  double acos, bcos, ccos, v;

  acos = cos(aa * 0.0174533);
  if (fabs(acos) < 0.00000001)
    acos = 0.0;
  bcos = cos(bb * 0.0174533);
  if (fabs(bcos) < 0.00000001)
    bcos = 0.0;
  ccos = cos(cc * 0.0174533);
  if (fabs(ccos) < 0.00000001)
    ccos = 0.0;
  ccos = cos(cc * 0.0174533);
  if (fabs(ccos) < 0.00000001)
    ccos = 0.0;
  v = 1 - acos * acos - bcos * bcos - ccos * ccos + 2 * acos * bcos * ccos;
  if (fabs(v) < 0.00000001)
    v = 0.0;
  return (a * b * c * sqrt(v));
}


Static boolean gleiche_vektoren(double *v1, double *v2)
{
  return (fabs(v1[0] - v2[0]) < 0.001 && fabs(v1[1] - v2[1]) < 0.001 &&
          fabs(v1[2] - v2[2]) < 0.001);
}


Static boolean gleiche_rotationsmatrizen(long (*x)[3], long (*y)[3])
{
  return (x[0][0] == y[0][0] && x[0][1] == y[0][1] && x[0][2] == y[0][2] &&
          x[1][0] == y[1][0] && x[1][1] == y[1][1] && x[1][2] == y[1][2] &&
          x[2][0] == y[2][0] && x[2][1] == y[2][1] && x[2][2] == y[2][2]);
}


Static boolean gleiche_symmetrieoperatoren(symmetrieoperator *x,
                                           symmetrieoperator *y)
{
  return (gleiche_rotationsmatrizen(x->r, y->r) & gleiche_vektoren(x->t, y->t));
}


#define EPS (1.e-5)

Static Void frac_vektor(double *va, double *vb)
{
  long  i;

  for (i = 0; i <= 2; i++)
  {
    vb[i] = va[i] - (long) va[i];
    for (;;) if (vb[i] >= 1.) vb[i] -= 1.; else break;
    for (;;) if (vb[i] <  0.) vb[i] += 1.; else break;
    if      ( 1. - EPS <= vb[i] && vb[i] <=  1. + EPS) vb[i] = 0.;
    else if (-1. - EPS <= vb[i] && vb[i] <= -1. + EPS) vb[i] = 0.;
    else if (    - EPS <= vb[i] && vb[i] <=       EPS) vb[i] = 0.;
  }
}


Static Void int_vektor(double *v, long *g)
{
  long    i;
  double  vn;

  for (i = 0; i <= 2; i++)
  {
                g[i] = (long) v[i];
    vn = v[i] - g[i];
    for (;;) if (vn >= 1.) { vn -= 1.; g[i]++; } else break;
    for (;;) if (vn <  0.) { vn += 1.; g[i]--; } else break;
    if      ( 1. - EPS <= vn && vn <=  1. + EPS) g[i]++;
    else if (-1. - EPS <= vn && vn <= -1. + EPS) g[i]--;
  }
}

#undef EPS


Static Void write_vektor(_TEXT *f, double *v)
{
  fprintf(f->f, " (%8.5f %8.5f %8.5f)", v[0], v[1], v[2]);
}


Static Void set_einheits_symop(symmetrieoperator *s)
{
  long i, j;

  for (i = 0; i <= 2; i++)
    s->t[i] = 0.0;
  for (i = 0; i <= 2; i++) {
    for (j = 0; j <= 2; j++)
      s->r[i][j] = 0;
  }
  for (i = 0; i <= 2; i++)
    s->r[i][i] = 1;
}


Static Void mult_rotmat_mit_vektor(long (*r)[3], double *v, double *vn)
{
  vn[0] = r[0][0] * v[0] + r[0][1] * v[1] + r[0][2] * v[2];
  vn[1] = r[1][0] * v[0] + r[1][1] * v[1] + r[1][2] * v[2];
  vn[2] = r[2][0] * v[0] + r[2][1] * v[1] + r[2][2] * v[2];
}


Static Void mult_symop_mit_vektor(symmetrieoperator *s, double *v, double *vn)
{
  vn[0] = s->t[0] + s->r[0][0] * v[0] + s->r[0][1] * v[1] + s->r[0][2] * v[2];
  vn[1] = s->t[1] + s->r[1][0] * v[0] + s->r[1][1] * v[1] + s->r[1][2] * v[2];
  vn[2] = s->t[2] + s->r[2][0] * v[0] + s->r[2][1] * v[1] + s->r[2][2] * v[2];
}


Static Void mult_rotationsmatrizen(long (*a)[3], long (*b)[3], long (*c)[3])
{
  long i, j;

  for (i = 0; i <= 2; i++) {
    for (j = 0; j <= 2; j++)
      c[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
  }
}


Static Void mult_symmetrieoperatoren(symmetrieoperator *a,
                                     symmetrieoperator *b,
                                     symmetrieoperator *c)
{
  symmetrieoperator  ab;

  mult_rotationsmatrizen(a->r, b->r, ab.r);
  mult_symop_mit_vektor(a, b->t, ab.t);

  memcpy(c, &ab, sizeof (*c));
}



Static Void subtrahiere_symmetrieoperatoren(symmetrieoperator *x,
                                            symmetrieoperator *y,
                                            symmetrieoperator *z)
{
  long i, j;

  for (i = 0; i <= 2; i++) {
    for (j = 0; j <= 2; j++)
      z->r[i][j] = x->r[i][j] - y->r[i][j];
  }
  for (i = 0; i <= 2; i++)
    z->t[i] = x->t[i] - y->t[i];
}



Static Void write_symmetrieoperator_with_centering(_TEXT *t, symmetrieoperator *x, vektor trans)
{
  long     i, j, nn;
  boolean  ersterterm;

  for (i = 0; i < 3; i++)
  {
    if (i) putc(',', t->f);

    ersterterm = true;

    for (j = 0; j < 3; j++)
    {
      if (x->r[i][j] != 0)
      {
        if (x->r[i][j] < 0)
          putc('-', t->f);
        else if (!ersterterm)
          putc('+', t->f);

        if (labs(x->r[i][j]) != 1)
          fprintf(t->f, "%ld", labs(x->r[i][j]));

        if (j == 0)
          putc('X', t->f);
        if (j == 1)
          putc('Y', t->f);
        if (j == 2)
          putc('Z', t->f);

        ersterterm = false;
      }
    }

    if (fabs(x->t[i]+trans[i]) > 0.001)
    {
      if (x->t[i]+trans[i] < 0.)
        putc('-', t->f);
      else
        putc('+', t->f);

      nn = 0;
      do {
        nn++;
      } while (fabs((x->t[i]+trans[i]) * nn - (long)floor((x->t[i]+trans[i]) * nn + 0.5)) >= 0.001);
      fprintf(t->f, "%ld", (long)floor(fabs((x->t[i]+trans[i]) * nn) + 0.5));
      if (nn > 1)
        fprintf(t->f, "/%ld", nn);

      ersterterm = false;
    }
    else if (ersterterm)
      putc('0', t->f);
  }
}

Static Void write_symmetrieoperator(_TEXT *t, symmetrieoperator *x)
{
  long     i, j, nn;
  boolean  ersterterm;


  for (i = 0; i < 3; i++)
  {
    if (i) putc(',', t->f);

    ersterterm = true;

    for (j = 0; j < 3; j++)
    {
      if (x->r[i][j] != 0)
      {
        if (x->r[i][j] < 0)
          putc('-', t->f);
        else if (!ersterterm)
          putc('+', t->f);

        if (labs(x->r[i][j]) != 1)
          fprintf(t->f, "%ld", labs(x->r[i][j]));

        if (j == 0)
          putc('X', t->f);
        if (j == 1)
          putc('Y', t->f);
        if (j == 2)
          putc('Z', t->f);

        ersterterm = false;
      }
    }

    if (fabs(x->t[i]) > 0.001)
    {
      if (x->t[i] < 0.)
        putc('-', t->f);
      else
        putc('+', t->f);

      nn = 0;
      do {
        nn++;
      } while (fabs(x->t[i] * nn - (long)floor(x->t[i] * nn + 0.5)) >= 0.001);
      fprintf(t->f, "%ld", (long)floor(fabs(x->t[i] * nn) + 0.5));
      if (nn > 1)
        fprintf(t->f, "/%ld", nn);

      ersterterm = false;
    }
    else if (ersterterm)
      putc('0', t->f);
  }
}

               /* till last non blank */
static void print_tlnb(FILE *fpout, const char *s, int n,
                       char *Head, char *Tail)
{
  int  i, j;


  if (Head)
    fprintf(fpout, "%s", Head);

  j = 0;

  for (i = 0; i < n; i++)
  {
    if (s[i] == '\0')
      break;

    if (s[i] != ' ')
      while (j <= i) putc(s[j++], fpout);
  }

  if (Tail)
    fprintf(fpout, "%s", Tail);
}


Static Void write_st_raumgruppe(_TEXT *t, st_raumgruppe *rginfo)
{
  long i, j;

  fprintf(t->f, " Symmetry information from file SYMDAT\n");

  print_tlnb(t->f, rginfo->raumgruppename, 30, " space group ", NULL);
  fprintf(t->f, "  (%ld)\n\n", rginfo->raumgruppenr);

  if (rginfo->anztranslationen > 0) {
    fprintf(t->f, "        ");
    for (i = 0; i < rginfo->anztranslationen; i++) {
      putc('(', t->f);
      for (j = 0; j <= 2; j++) {
        if (fabs(rginfo->translation[i][j]) < 0.0001)
          fprintf(t->f, " 0");
        if (fabs(rginfo->translation[i][j] - 0.5) < 0.0001)
          fprintf(t->f, " 1/2");
        if (fabs(rginfo->translation[i][j] - 1.0 / 3) < 0.0001)
          fprintf(t->f, " 1/3");
        if (fabs(rginfo->translation[i][j] - 2.0 / 3) < 0.0001)
          fprintf(t->f, " 2/3");
        if (j + 1 != 3)
          putc(',', t->f);
      }
      fprintf(t->f, ")+   ");
    }
    fprintf(t->f, "\n\n");
  }
  for (i = 1; i <= rginfo->anz_rotmat; i++) {
    fprintf(t->f, " (%3ld) ", i);
    write_symmetrieoperator(t, rginfo->symop[i - 1]);
    putc('\n', t->f);
  }
  putc('\n', t->f);
}


Static Void read_symmetrieoperator(_TEXT *f, symmetrieoperator *symop)
{
  double i1, i2;
  long vorz, i, j, k;
  double rr;

  while (P_peek(f->f) == ' ')
    getc(f->f);
  for (i = 0; i <= 2; i++) {
    for (j = 0; j <= 2; j++)
      symop->r[i][j] = 0;
  }
  for (i = 0; i <= 2; i++)
    symop->t[i] = 0.0;
  k = 1;
  do {
    vorz = 1;
    if (P_peek(f->f) == ',') {
      k++;
      getc(f->f);
    }
    if (P_peek(f->f) == '-') {
      vorz = -1;
      getc(f->f);
    }
    if (P_peek(f->f) == '+')
      getc(f->f);
    if (isdigit(P_peek(f->f))) {
      fscanf(f->f, "%lg", &i1);
      if (!P_eoln(f->f)) {
        if (P_peek(f->f) == '/') {
          getc(f->f);
          fscanf(f->f, "%lg", &i2);
          rr = i1 / i2;
        } else
          rr = i1;
      } else
        rr = i1;
    } else if (P_peek(f->f) == 'Z' || P_peek(f->f) == 'Y' ||
               P_peek(f->f) == 'X')
      rr = 1.0;
    rr *= vorz;
    if (P_peek(f->f) == 'X')
      symop->r[k - 1][0] = (long)floor(rr + 0.5);
    else if (P_peek(f->f) == 'Y')
      symop->r[k - 1][1] = (long)floor(rr + 0.5);
    else if (P_peek(f->f) == 'Z')
      symop->r[k - 1][2] = (long)floor(rr + 0.5);
    else
      symop->t[k - 1] = rr;
    if (P_peek(f->f) == 'Z' || P_peek(f->f) == 'Y' || P_peek(f->f) == 'X')
      getc(f->f);
    if (P_peek(f->f) == ')')
      getc(f->f);
  } while (!P_eoln(f->f));
}


Static Void readln_raumgruppename(_TEXT *t, Char *rn)
{
  long i;

  for (i = 0; i <= 29; i++)
    rn[i] = ' ';
  while (P_peek(t->f) == ' ')
    getc(t->f);
  i = 0;
  while (!P_eoln(t->f)) {
    i++;
    rn[i - 1] = getc(t->f);
    if (rn[i - 1] == '\n')
      rn[i - 1] = ' ';
    if (islower(rn[i - 1]))
      rn[i - 1] = _toupper(rn[i - 1]);
  }
  fscanf(t->f, "%*[^\n]");
  getc(t->f);
}


Static Void suche_raumgruppe_entry(_TEXT *t, Char *rn_, symmetrieoperator *tra)
{
  st_raumgruppename rn;
  long i, j;
  st_raumgruppename rnd;
  symmetrieoperator tran;
  _TEXT tt;
  boolean b = false;

  memcpy(rn, rn_, sizeof(st_raumgruppename));
  tt.f = NULL;
  *tt.name = '\0';
  j = 0;
  for (i = 1; i <= 30; i++) {
    if (rn[i - 1] == '(')
      j = i;
  }
  if (j > 0) {
    rn[j - 1] = ' ';
    if (*tt.name != '\0') {
      if (tt.f != NULL)
        tt.f = freopen(tt.name, "w", tt.f);
      else
        tt.f = fopen(tt.name, "w");
    } else {
      if (tt.f != NULL)
        rewind(tt.f);
      else
        tt.f = tmpfile();
    }
    if (tt.f == NULL)
      _EscIO(FileNotFound);
    SETUPBUF(tt.f, Char);
    i = j;
    while (rn[i - 1] != ')') {
      i++;
      putc(rn[i - 1], tt.f);
    }
    putc('\n', tt.f);
    if (*tt.name != '\0') {
      if (tt.f != NULL)
        tt.f = freopen(tt.name, "r", tt.f);
      else
        tt.f = fopen(tt.name, "r");
    } else
      rewind(tt.f);
    if (tt.f == NULL)
      _EscIO(FileNotFound);
    RESETBUF(tt.f, Char);
    read_symmetrieoperator(&tt, tra);
    fscanf(tt.f, "%*[^\n]");
    getc(tt.f);
    for (i = j - 1; i <= 29; i++)
      rn[i] = ' ';
  } else
    set_einheits_symop(tra);
  if (*t->name != '\0') {
    if (t->f != NULL)
      t->f = freopen(t->name, "r", t->f);
    else
      t->f = fopen(t->name, "r");
  } else
    rewind(t->f);
  if (t->f == NULL)
    _EscIO(FileNotFound);
  RESETBUF(t->f, Char);
  j = 0;
  for (i = 0; i <= 29; i++) {
    if (rn[i] != ' ') {
      j++;
      rn[j - 1] = rn[i];
    }
  }
  for (i = j; i <= 29; i++)
    rn[i] = ' ';
  for (i = 0; i <= 29; i++)
    rnd[i] = ' ';
  do {
    b = false;
    do {
      if (BUFEOF(t->f))
        b = true;
      else if (P_peek(t->f) == '*')
        b = true;
      else {
        fscanf(t->f, "%*[^\n]");
        getc(t->f);
      }
    } while (!b);
    if (BUFEOF(t->f))
      b = false;
    else {
      getc(t->f);
      for (j = 0; j <= 15; j++) {
        rnd[j] = getc(t->f);
        if (rnd[j] == '\n')
          rnd[j] = ' ';
      }
      if (!strncmp(rnd, "END                           ",
                   sizeof(st_raumgruppename)))
        b = false;
    }
  } while (strncmp(rnd, rn, sizeof(st_raumgruppename)) && b);
  if (!b) {
    printf(" ERROR  space group %.30s not found => P1\n", rn);
    memcpy(rn, "P1                            ", sizeof(st_raumgruppename));
    suche_raumgruppe_entry(t, rn, tra);
  }
  if (P_peek(t->f) == '=') {
    getc(t->f);
    getc(t->f);
    getc(t->f);
    readln_raumgruppename(t, rn);
    suche_raumgruppe_entry(t, rn, &tran);
    mult_symmetrieoperatoren(tra, &tran, tra);
  }
  if (tt.f != NULL)
    fclose(tt.f);
}


Static Void erstelle_multtab(st_raumgruppe *rginfo)
{
  long i, j;
  rotationsmatrix r;
  boolean b;
  long FORLIM, FORLIM1;

  FORLIM = rginfo->anz_rotmat;
  for (i = 0; i < FORLIM; i++) {
    FORLIM1 = rginfo->anz_rotmat;
    for (j = 0; j < FORLIM1; j++) {
      mult_rotationsmatrizen(rginfo->symop[i]->r, rginfo->symop[j]->r, r);
      b = false;
      rginfo->multtab[i][j] = rginfo->anz_rotmat + 1;
      do {
        rginfo->multtab[i][j]--;
        if (rginfo->multtab[i][j] == 0) {
          printf(" ERROR  space group is wrong (not complete?)\n");
          b = true;
        } else if (gleiche_rotationsmatrizen(r,
                     rginfo->symop[rginfo->multtab[i][j] - 1]->r))
          b = true;
      } while (!b);
    }
  }
}


typedef double augm_matrix[4][4];


Local Void austauschschritt(double (*s)[4], long x, long y)
{
  long i, j;

  for (i = 0; i <= 3; i++) {
    for (j = 0; j <= 3; j++) {
      if (i + 1 != x && j + 1 != y)
        s[i][j] -= s[x - 1][j] * s[i][y - 1] / s[x - 1][y - 1];
    }
  }
  for (i = 0; i <= 3; i++) {
    if (i + 1 != x)
      s[i][y - 1] /= s[x - 1][y - 1];
  }
  for (j = 0; j <= 3; j++) {
    if (j + 1 != y)
      s[x - 1][j] = -(s[x - 1][j] / s[x - 1][y - 1]);
  }
  s[x - 1][y - 1] = 1 / s[x - 1][y - 1];
}


Static Void inverser_symmetrieoperator(symmetrieoperator *sa,
                                       symmetrieoperator *sb)
{
  augm_matrix s, ss;
  long i, j, n;
  long z[4], k[4];
  boolean b;

  for (i = 0; i <= 2; i++) {
    for (j = 0; j <= 2; j++)
      s[i][j] = sa->r[i][j];
  }
  for (i = 0; i <= 2; i++)
    s[i][3] = sa->t[i];
  for (j = 0; j <= 2; j++)
    s[3][j] = 0.0;
  s[3][3] = 1.0;
  for (j = 1; j <= 4; j++)
    z[j - 1] = j;
  for (i = 1; i <= 4; i++)
    k[i - 1] = -i;
  b = false;
  do {
    i = 0;
    do {
      i++;
    } while (k[i - 1] >= 0);
    j = 0;
    do {
      j++;
    } while (z[j - 1] <= 0 || s[i - 1][j - 1] == 0);
    austauschschritt(s, i, j);
    n = k[i - 1];
    k[i - 1] = z[j - 1];
    z[j - 1] = n;
    if (z[0] < 0 && z[1] < 0 && z[2] < 0 && z[3] < 0)
      b = true;
  } while (!b);
  for (i = 0; i <= 3; i++) {
    for (j = 0; j <= 3; j++)
      ss[k[i] - 1][-z[j] - 1] = s[i][j];
  }
  for (i = 0; i <= 2; i++) {
    for (j = 0; j <= 2; j++)
      sb->r[i][j] = (long)floor(ss[i][j] + 0.5);
  }
  for (i = 0; i <= 2; i++)
    sb->t[i] = ss[i][3];
}


Local Void berechne_bravaistyp(st_raumgruppe *rginfo)
{
  rginfo->anztranslationen = 0;
  if (rginfo->raumgruppename[0] == 'P') {
    rginfo->anztranslationen = 0;
    return;
  }
  if (rginfo->raumgruppename[0] == 'A') {
    rginfo->anztranslationen = 1;
    rginfo->translation[0][0] = 0.0;
    rginfo->translation[0][1] = 1.0 / 2;
    rginfo->translation[0][2] = 1.0 / 2;
    return;
  }
  if (rginfo->raumgruppename[0] == 'B') {
    rginfo->anztranslationen = 1;
    rginfo->translation[0][0] = 1.0 / 2;
    rginfo->translation[0][1] = 0.0;
    rginfo->translation[0][2] = 1.0 / 2;
    return;
  }
  if (rginfo->raumgruppename[0] == 'C') {
    rginfo->anztranslationen = 1;
    rginfo->translation[0][0] = 1.0 / 2;
    rginfo->translation[0][1] = 1.0 / 2;
    rginfo->translation[0][2] = 0.0;
    return;
  }
  if (rginfo->raumgruppename[0] == 'F') {
    rginfo->anztranslationen = 3;
    rginfo->translation[0][0] = 0.0;
    rginfo->translation[0][1] = 1.0 / 2;
    rginfo->translation[0][2] = 1.0 / 2;
    rginfo->translation[1][0] = 1.0 / 2;
    rginfo->translation[1][1] = 0.0;
    rginfo->translation[1][2] = 1.0 / 2;
    rginfo->translation[2][0] = 1.0 / 2;
    rginfo->translation[2][1] = 1.0 / 2;
    rginfo->translation[2][2] = 0.0;
    return;
  }
  if (rginfo->raumgruppename[0] == 'I') {
    rginfo->anztranslationen = 1;
    rginfo->translation[0][0] = 1.0 / 2;
    rginfo->translation[0][1] = 1.0 / 2;
    rginfo->translation[0][2] = 1.0 / 2;
    return;
  }
  if (rginfo->raumgruppename[0] != 'R') {
    printf(" ERROR  space group %c has unknown centering\n",
           rginfo->raumgruppename[0]);
    return;
  }
  if (!strncmp(rginfo->kristallsystem, "RHO ", sizeof(st_kristallsystem)))
    return;
  rginfo->anztranslationen = 2;
  rginfo->translation[0][0] = 1.0 / 3;
  rginfo->translation[0][1] = 2.0 / 3;
  rginfo->translation[0][2] = 2.0 / 3;
  rginfo->translation[1][0] = 2.0 / 3;
  rginfo->translation[1][1] = 1.0 / 3;
  rginfo->translation[1][2] = 1.0 / 3;
}

Local Void spieglesymopansymzent(st_raumgruppe *rginfo)
{
  long i;
  symmetrieoperator inv;
  long FORLIM;

  set_einheits_symop(&inv);
  for (i = 0; i <= 2; i++)
    inv.r[i][i] = -1;
  FORLIM = rginfo->zaehligkeit;
  for (i = 0; i < FORLIM; i++) {
    rginfo->symop[rginfo->zaehligkeit + i] =
      (symmetrieoperator *)Malloc(sizeof(symmetrieoperator));
    mult_symmetrieoperatoren(&inv, rginfo->symop[i],
                             rginfo->symop[rginfo->zaehligkeit + i]);
  }
  rginfo->zaehligkeit *= 2;
}

Local Void bestimme_setting(st_raumgruppe *rginfo, Char *kristallsystem)
{
  boolean n[3];
  long rr[3];
  long i, j;
  symmetrieoperator *WITH;

  for (i = 0; i <= 2; i++)
    n[i] = false;
  for (i = 0; i < rginfo->zaehligkeit; i++) {
    WITH = rginfo->symop[i];
    for (j = 0; j <= 2; j++)
      rr[j] = WITH->r[0][j] + WITH->r[1][j] + WITH->r[2][j];
    if (rr[0] < 0 && rr[1] < 0 && rr[2] > 0)
      n[2] = true;
    if (rr[0] > 0 && rr[1] > 0 && rr[2] < 0)
      n[2] = true;
    if (rr[0] < 0 && rr[1] > 0 && rr[2] < 0)
      n[1] = true;
    if (rr[0] > 0 && rr[1] < 0 && rr[2] > 0)
      n[1] = true;
    if (rr[0] > 0 && rr[1] < 0 && rr[2] < 0)
      n[0] = true;
    if (rr[0] < 0 && rr[1] > 0 && rr[2] > 0)
      n[0] = true;
  }
  if (n[0] == false && n[1] == false && n[2] == true)
    memcpy(kristallsystem, "MON1", sizeof(st_kristallsystem));
  if (n[0] == false && n[1] == true && n[2] == false)
    memcpy(kristallsystem, "MON2", sizeof(st_kristallsystem));
  if (n[0] == true && n[1] == false && n[2] == false)
    memcpy(kristallsystem, "MON3", sizeof(st_kristallsystem));
}


Static Void read_raumgruppe_entry(_TEXT *symdat, Char *rn,
                                  symmetrieoperator *transf,
                                  st_raumgruppe *rginfo)
{
  long i, j, k;
  symmetrieoperator inv_transf;
  long FORLIM, FORLIM1;
  symmetrieoperator *WITH1;

  memcpy(rginfo->raumgruppename, rn, sizeof(st_raumgruppename));
  fscanf(symdat->f, "%ld", &rginfo->raumgruppenr);
  while (P_peek(symdat->f) == ' ')
    getc(symdat->f);
  for (i = 0; i <= 3; i++) {
    rginfo->kristallsystem[i] = getc(symdat->f);
    if (rginfo->kristallsystem[i] == '\n')
      rginfo->kristallsystem[i] = ' ';
  }
  fscanf(symdat->f, "%ld%*[^\n]", &i);
  getc(symdat->f);
  if (i == -1)
    rginfo->zentrosymm = true;
  else
    rginfo->zentrosymm = false;
  rginfo->anz_rotmat = 0;
  while (P_peek(symdat->f) == '*') {
    fscanf(symdat->f, "%*[^\n]");
    getc(symdat->f);
  }
  while (P_peek(symdat->f) != '=') {
    rginfo->anz_rotmat++;
    rginfo->symop[rginfo->anz_rotmat - 1] =
      (symmetrieoperator *)Malloc(sizeof(symmetrieoperator));
    read_symmetrieoperator(symdat, rginfo->symop[rginfo->anz_rotmat - 1]);
    fscanf(symdat->f, "%*[^\n]");
    getc(symdat->f);
  }
  rginfo->zaehligkeit = rginfo->anz_rotmat;
  if (rginfo->zentrosymm)
    spieglesymopansymzent(rginfo);
  inverser_symmetrieoperator(transf, &inv_transf);
  FORLIM = rginfo->zaehligkeit;
  for (i = 0; i < FORLIM; i++) {
    mult_symmetrieoperatoren(transf, rginfo->symop[i], rginfo->symop[i]);
    mult_symmetrieoperatoren(rginfo->symop[i], &inv_transf, rginfo->symop[i]);
    frac_vektor(rginfo->symop[i]->t, rginfo->symop[i]->t);
  }
  berechne_bravaistyp(rginfo);
  rginfo->anz_rotmat = rginfo->zaehligkeit;
  FORLIM = rginfo->anztranslationen;
  for (i = 1; i <= FORLIM; i++) {
    FORLIM1 = rginfo->zaehligkeit;
    for (j = 0; j < FORLIM1; j++) {
      rginfo->symop[i * rginfo->zaehligkeit + j] =
        (symmetrieoperator *)Malloc(sizeof(symmetrieoperator));
      *rginfo->symop[i * rginfo->zaehligkeit + j] = *rginfo->symop[j];
      WITH1 = rginfo->symop[i * rginfo->zaehligkeit + j];
      for (k = 0; k <= 2; k++)
        WITH1->t[k] += rginfo->translation[i - 1][k];
      WITH1 = rginfo->symop[i * rginfo->zaehligkeit + j];
      frac_vektor(WITH1->t, WITH1->t);
    }
  }
  rginfo->zaehligkeit *= rginfo->anztranslationen + 1;
  erstelle_multtab(rginfo);
  if (!strncmp(rginfo->kristallsystem, "MON ", sizeof(st_kristallsystem)))
    bestimme_setting(rginfo, rginfo->kristallsystem);
}


Static Void read_raumgruppe(_TEXT *t, Char *rn, st_raumgruppe *rg)
{
  symmetrieoperator tra;

  suche_raumgruppe_entry(t, rn, &tra);
  read_raumgruppe_entry(t, rn, &tra, rg);
}


Static boolean polare_richtung(double *v, st_raumgruppe *rginfo)
{
  boolean b;
  long i;
  vektor v2;

  b = true;
  for (i = 0; i < rginfo->zaehligkeit; i++) {
    mult_rotmat_mit_vektor(rginfo->symop[i]->r, v, v2);
    if (!gleiche_vektoren(v2, v))
      b = false;
  }
  return b;
}


Static boolean gleiche_atome(kt_atom *x, kt_atom *y)
{
  return (x->dessen_atomart    == y->dessen_atomart &&
          x->elementarzelle[0] == y->elementarzelle[0] &&
          x->elementarzelle[1] == y->elementarzelle[1] &&
          x->elementarzelle[2] == y->elementarzelle[2]);
}


Static Void berechne_atom_trilinpol(kt_atom *a, symmetrieoperator *ssy)
{
  long i;

  *ssy = a->dessen_atomart->trilinpol;
  for (i = 0; i <= 2; i++)
    ssy->t[i] += a->elementarzelle[i];
}


typedef double sympar[193][4];


typedef long reihe[3];
typedef double matrix[3][4];


Local Void austauschschritt_(double (*s)[4], long x, long y)
{
  long i, j;

  for (i = 1; i <= 192; i++) {
    for (j = 0; j <= 3; j++) {
      if (i != x && j + 1 != y)
        s[i][j] -= s[x][j] * s[i][y - 1] / s[x][y - 1];
    }
  }
  for (j = 0; j <= 3; j++) {
    if (j + 1 != y)
      s[x][j] = -(s[x][j] / s[x][y - 1]);
  }
  for (i = 1; i <= 192; i++)
    s[i][y - 1] = 0.0;
}

Local Void sym(double (*sy)[4], long *n, double (*yy)[4])
{
  long i, j, m;
  long pa[3];

  for (i = 0; i <= 3; i++)
    sy[0][i] = 0.0;
  for (i = 0; i <= 2; i++)
    pa[i] = 0;
  for (i = 0; i <= 2; i++) {
    m = 192;
    do {
      m--;
    } while ((fabs(sy[m][n[i] - 1]) <= 0.00001 || m == pa[0] || m == pa[1] ||
              m == pa[2]) && m != 0);
    if (m != 0) {
      austauschschritt_(sy, m, n[i]);
      pa[n[i] - 1] = m;
    }
  }
  for (i = 0; i <= 2; i++) {
    if (pa[i] != 0) {
      for (j = 0; j <= 3; j++)
        yy[i][j] = sy[pa[i]][j];
    } else {
      for (j = 0; j <= 3; j++)
        yy[i][j] = 0.0;
      yy[i][i] = 1.0;
    }
  }
}

Local Void simplex(double (*sy)[4], symmetrieoperator *y)
{
  matrix yy;
  reihe n;
  sympar ssy;
  long i, j;
  boolean loesung_io;

  memcpy(ssy, sy, sizeof(sympar));
  n[0] = 1;
  n[1] = 2;
  n[2] = 3;
  sym(ssy, n, yy);
  loesung_io = true;
  for (i = 0; i <= 2; i++) {
    for (j = 0; j <= 2; j++) {
      if (fabs((long)floor(yy[i][j] + 0.5) - yy[i][j]) > 0.00001)
        loesung_io = false;
    }
  }
  for (i = 0; i <= 2; i++) {
    if (fabs(yy[i][i]) > 0.00001 && fabs(yy[i][3]) > 0.00001)
      loesung_io = false;
  }
  if (!loesung_io) {
    memcpy(ssy, sy, sizeof(sympar));
    n[0] = 2;
    n[1] = 1;
    n[2] = 3;
    sym(ssy, n, yy);
  }
  for (i = 0; i <= 2; i++) {
    for (j = 0; j <= 2; j++)
      y->r[i][j] = (long)floor(yy[i][j] + 0.5);
  }
  for (i = 0; i <= 2; i++)
    y->t[i] = yy[i][3];
}


Static Void berechne_symmetrielage(st_raumgruppe *rginfo, double *v,
                                   symmetrieoperator *s)
{
  long i, j, k, n;
  vektor w;
  gittervektor g;
  sympar a;
  symmetrieoperator ss;

  n = 0;
  for (i = 1; i <= 192; i++) {
    for (j = 0; j <= 3; j++)
      a[i][j] = 0.0;
  }
  for (i = 1; i < rginfo->zaehligkeit; i++) {
    mult_symop_mit_vektor(rginfo->symop[i], v, w);
    int_vektor(w, g);
    frac_vektor(w, w);
    if (gleiche_vektoren(v, w)) {
      ss = *rginfo->symop[i];
      for (j = 0; j <= 2; j++)
        ss.t[j] -= g[j];
      subtrahiere_symmetrieoperatoren(rginfo->symop[0], &ss, &ss);
      n++;
      for (j = 1; j <= 3; j++) {
        for (k = 0; k <= 2; k++)
          a[n * 3 + j - 3][k] = ss.r[j - 1][k];
        a[n * 3 + j - 3][3] = ss.t[j - 1];
      }
    }
  }
  if (n > 0)
    simplex(a, s);
  else
    *s = *rginfo->symop[0];
}


Static Void berechne_atomlagekoordinaten(st_raumgruppe *rginfo,
                                         kt_atomlage *atomlage)
{
  long i, j;
  gittervektor g;
  symmetrieoperator s;
  long FORLIM1;
  kt_atomart *WITH1;
  long intifa;

  frac_vektor(atomlage->atomlagekoord, atomlage->atomlagekoord);
  atomlage->symmetrielage
    = (symmetrieoperator *)Malloc(sizeof(symmetrieoperator));
  berechne_symmetrielage(rginfo, atomlage->atomlagekoord,
                         atomlage->symmetrielage);
  int_vektor(atomlage->symmetrielage->t, g);
  frac_vektor(atomlage->symmetrielage->t, atomlage->symmetrielage->t);
  for (i = 0; i <= 2; i++)
    atomlage->atomlagekoord[i] -= g[i];
  atomlage->anzatome = 0;
  for (i = 0; i < rginfo->zaehligkeit; i++) {
    mult_symmetrieoperatoren(rginfo->symop[i], atomlage->symmetrielage, &s);
    frac_vektor(s.t, s.t);
    intifa = atomlage->anzatome + 1;
    FORLIM1 = atomlage->anzatome;
    for (j = 1; j <= FORLIM1; j++) {
      if (gleiche_symmetrieoperatoren(&s, &atomlage->atom[j - 1]->trilinpol))
        intifa = j;
    }
    if (intifa > atomlage->anzatome) {
      atomlage->anzatome++;
      atomlage->atom[atomlage->anzatome - 1] =
        (kt_atomart *)Malloc(sizeof(kt_atomart));
      WITH1 = atomlage->atom[atomlage->anzatome - 1];
      WITH1->trilinpol = s;
      WITH1->atomartnr = atomlage->anzatome;
      WITH1->deren_atomlage = atomlage;
      for (j = 0; j < maxanzahlbindungen; j++)
        WITH1->bindung[j].dessen_atomart = NULL;
    }
  }
  FORLIM1 = atomlage->anzatome;
  for (i = 0; i < FORLIM1; i++) {
    WITH1 = atomlage->atom[i];
    mult_symop_mit_vektor(&WITH1->trilinpol, atomlage->atomlagekoord,
                          WITH1->pkoord);
    int_vektor(WITH1->pkoord, g);
    frac_vektor(WITH1->pkoord, WITH1->pkoord);
    for (j = 0; j <= 2; j++)
      WITH1->trilinpol.t[j] -= g[j];
  }
  memcpy(atomlage->atomlagekoord, atomlage->atom[0]->pkoord, sizeof(vektor));
}


Static Void write_atom_trilinpol(_TEXT *t, kt_atom *a)
{
  symmetrieoperator ssy;

  berechne_atom_trilinpol(a, &ssy);
  write_symmetrieoperator(t, &ssy);
}


Static Void alg_symmetrieoperation(symmetrieoperator *s, kt_atom *a, kt_atom *b)
{
  long i;
  kt_atomlage *p;
  symmetrieoperator sa, sr;
  gittervektor ga, gr;

  p = a->dessen_atomart->deren_atomlage;
  berechne_atom_trilinpol(a, &sa);
  mult_symmetrieoperatoren(s, &sa, &sa);
  int_vektor(sa.t, ga);
  frac_vektor(sa.t, sa.t);
  i = 0;
  do {
    i++;
    sr = p->atom[i - 1]->trilinpol;
    int_vektor(sr.t, gr);
    frac_vektor(sr.t, sr.t);
  } while (!gleiche_symmetrieoperatoren(&sr, &sa));
  b->dessen_atomart = p->atom[i - 1];
  for (i = 0; i <= 2; i++)
    b->elementarzelle[i] = ga[i] - gr[i];
}


Static Void berechne_atom_koordinaten(kt_atom *a, double *v)
{
  long i;

  for (i = 0; i <= 2; i++)
    v[i] = a->elementarzelle[i] + a->dessen_atomart->pkoord[i];
}


Static Void write_atomart_name(_TEXT *f, kt_atomart *a)
{
  long i;

  putc(' ', f->f);
  for (i = 0; i <= 4; i++) 
    // StefS 17-09-2012 -> increase from i from 3 to 4 
    // to properly print long atom names (ie Si120)
    putc(a->deren_atomlage->aname[i], f->f);
  fprintf(f->f, "(%3ld)", a->atomartnr);
}


Static Void write_atom_name(_TEXT *f, kt_atom *a)
{
  write_atomart_name(f, a->dessen_atomart);
  fprintf(f->f, "(%2ld%2ld%2ld)",
          a->elementarzelle[0], a->elementarzelle[1], a->elementarzelle[2]);
}


Static Void write_atom_koordinaten(_TEXT *t, kt_atom *a)
{
  vektor v;

  berechne_atom_koordinaten(a, v);
  write_vektor(t, v);
}


Static Void write_strukturdaten(_TEXT *t, kt_strukturdaten *strukturdaten)
{
  long         i;
  double       vol;
  kt_atomlage  *al;


  vol = volumen(strukturdaten->zellpar.laenge[0],
      strukturdaten->zellpar.laenge[1], strukturdaten->zellpar.laenge[2],
      strukturdaten->zellpar.winkel[0], strukturdaten->zellpar.winkel[1],
      strukturdaten->zellpar.winkel[2]);

  print_tlnb(t->f, strukturdaten->name, 80, "\n Name: ", "\n");
  print_tlnb(t->f, strukturdaten->lit,  80,   " Lit: ",  "\n\n");

  fprintf(t->f,
          " Cell parameters:  a = %6.3f  b = %6.3f  c = %6.3f\n",
          strukturdaten->zellpar.laenge[0], strukturdaten->zellpar.laenge[1],
          strukturdaten->zellpar.laenge[2]);
  fprintf(t->f,
          "                   alpha = %6.3f  beta = %6.3f  gamma = %6.3f\n",
          strukturdaten->zellpar.winkel[0], strukturdaten->zellpar.winkel[1],
          strukturdaten->zellpar.winkel[2]);
  fprintf(t->f,
          "                   volume = %.6g\n\n", vol);

  print_tlnb(t->f, strukturdaten->raumgruppe, 30, " Space group: ", NULL);
  fprintf(t->f, "  (%ld)\n\n", rginfo.raumgruppenr);

  fprintf(t->f, " Atoms with coordination number and atoms per unit cell\n");

  for (i = 0; i < strukturdaten->anzatomlagen; i++)
  {
    al = strukturdaten->atomlage[i];

    fprintf(t->f, " (%2ld)   %.6s", i + 1, al->aname);
    write_atom_koordinaten(t, &al->grundatom);
    fprintf(t->f, "%4ld%4ld   ", al->koordinationszahl, al->anzatome);
    write_atom_trilinpol(t, &al->grundatom);
    putc('\n', t->f);
  }
}


Static Void write_cssr(_TEXT *t, kt_strukturdaten *strukturdaten)
{
  long         i, j;
  kt_atomlage  *WITH;
  vektor       v;


  for (i = 0; i < 38; i++) putc(' ', t->f);
  fprintf(t->f, "%8.3f%8.3f%8.3f\n",
          strukturdaten->zellpar.laenge[0],
          strukturdaten->zellpar.laenge[1],
          strukturdaten->zellpar.laenge[2]);

  for (i = 0; i < 21; i++) putc(' ', t->f);
  fprintf(t->f, "%8.3f%8.3f%8.3f    SPGR =%3ld ",
          strukturdaten->zellpar.winkel[0],
          strukturdaten->zellpar.winkel[1],
          strukturdaten->zellpar.winkel[2],
          rginfo.raumgruppenr);
  print_tlnb(t->f, strukturdaten->raumgruppe, 30, NULL, "\n");

  fprintf(t->f, "%4ld%4d ", strukturdaten->anzatomlagen, 0);
  print_tlnb(t->f, strukturdaten->name, 80, NULL, "\n");

  print_tlnb(t->f, strukturdaten->lit, 80, "       ", "\n");

  for (i = 1; i <= strukturdaten->anzatomlagen; i++) {
    WITH = strukturdaten->atomlage[i - 1];
    berechne_atom_koordinaten(&WITH->grundatom, v);
    fprintf(t->f, "%4ld %5.5s %9.5f %9.5f %9.5f ",
            i,
            WITH->aname,
            v[0], v[1], v[2]);

    for (j = 0; j < 8; j++)
      fprintf(t->f, "%4d", 0);

    fprintf(t->f, " %7.3f\n", 0.);
  }
}





Static Void write_cif(_TEXT *t, kt_strukturdaten *strukturdaten, st_raumgruppe *rginfo) // Stefs 05-07-2012
{
  long         i, j;
  kt_atomlage  *WITH;
  vektor       v;
  

  print_tlnb(t->f, strukturdaten->abk, 16, "data_", "\n");
  //fprintf(t->f, "data_kriber");
  fprintf(t->f, "_chemical_name_mineral ??\n");
  fprintf(t->f, "_cell_length_a  %10.5f\n", strukturdaten->zellpar.laenge[0]);
  fprintf(t->f, "_cell_length_b  %10.5f\n", strukturdaten->zellpar.laenge[1]);
  fprintf(t->f, "_cell_length_c  %10.5f\n", strukturdaten->zellpar.laenge[2]);
  fprintf(t->f, "_cell_angle_alpha %10.5f\n", strukturdaten->zellpar.winkel[0]);
  fprintf(t->f, "_cell_angle_beta  %10.5f\n", strukturdaten->zellpar.winkel[1]);
  fprintf(t->f, "_cell_angle_gamma %10.5f\n", strukturdaten->zellpar.winkel[2]);
  fprintf(t->f, "_cell_volume ??\n");
  fprintf(t->f, "_symmetry_space_group_name_H-M '");
  print_tlnb(t->f, strukturdaten->raumgruppe, 30, NULL, "'\n");
  fprintf(t->f, "\n");
  fprintf(t->f, "loop_\n");
  fprintf(t->f, "  _symmetry_equiv_pos_as_xyz\n");
  
  for (i = 1; i <= rginfo->anz_rotmat; i++) {
    write_symmetrieoperator(t, rginfo->symop[i - 1]);
    putc('\n', t->f);
  }

  if (rginfo->anztranslationen > 0) {
    for (i = 0; i < rginfo->anztranslationen; i++) {
      for (j = 0; j < rginfo->anz_rotmat; j++) {
        write_symmetrieoperator_with_centering(t, rginfo->symop[j], rginfo->translation[i]);
        putc('\n', t->f);
      }
    }
  }
  putc('\n', t->f);

  fprintf(t->f, "\n");
  fprintf(t->f, "loop_\n");
  fprintf(t->f, "  _atom_site_label\n");
  //fprintf(t->f, "  _atom_site_type_symbol\n");
  //fprintf(t->f, "  _atom_site_symmetry_multiplicity\n");
  fprintf(t->f, "  _atom_site_fract_x\n");
  fprintf(t->f, "  _atom_site_fract_y\n");
  fprintf(t->f, "  _atom_site_fract_z\n");

  for (i = 1; i <= strukturdaten->anzatomlagen; i++) {
    WITH = strukturdaten->atomlage[i - 1];
    berechne_atom_koordinaten(&WITH->grundatom, v);
    fprintf(t->f, "%5.5s %9.5f %9.5f %9.5f \n",
            WITH->aname,
            v[0], v[1], v[2]);


  }
}


Static Void write_xtl(_TEXT *t, kt_strukturdaten *strukturdaten)
{
  long         i;
  kt_atomlage  *al;
  vektor       v;

  int ii,jj,ll;
  

  print_tlnb(t->f, strukturdaten->name, 80, "TITLE ", "\n");
  print_tlnb(t->f, strukturdaten->lit,  80, "TITLE ", "\n");

  fprintf(t->f, "CELL\n");
  fprintf(t->f, "    %.6g %.6g %.6g %.6g %.6g %.6g\n",
    strukturdaten->zellpar.laenge[0],
    strukturdaten->zellpar.laenge[1],
    strukturdaten->zellpar.laenge[2],
    strukturdaten->zellpar.winkel[0],
    strukturdaten->zellpar.winkel[1],
    strukturdaten->zellpar.winkel[2]);

  fprintf(t->f, "SYMMETRY LABEL ");

  for (i = 0; i < 30; i++)
    if (isspace(strukturdaten->raumgruppe[i]) == 0)
      putc(strukturdaten->raumgruppe[i], t->f);

  putc('\n', t->f);

    /* from here... */
  for(ii=0;ii<=rginfo.zaehligkeit-1;ii++)
  {
    fprintf(t->f,"SYM MAT");
    for(jj=0;jj<=2;jj++)
      for(ll=0;ll<=2;ll++)
        fprintf(t->f,"%5.1lf",(double)(rginfo.symop[ii]->r[jj][ll]));
      for(jj=0;jj<=2;jj++)
        fprintf(t->f,"%7.4lf",(double)(rginfo.symop[ii]->t[jj]));
      fprintf(t->f,"\n");
  }
  fprintf(t->f,"\n");
  /* ...to here -- modified by Christian and Sini$a */


  fprintf(t->f, "ATOMS\n");
  fprintf(t->f, " NAME    X         Y         Z\n");

  for (i = 0; i < strukturdaten->anzatomlagen; i++)
  {
    al = strukturdaten->atomlage[i];

    berechne_atom_koordinaten(&al->grundatom, v);

    fprintf(t->f, "%5.5s %9.5f %9.5f %9.5f\n",
      al->aname, v[0], v[1], v[2]);
  }

  fprintf(t->f, "EOF\n");
}


Static Void write_mvm(_TEXT *t, kt_strukturdaten *strukturdaten, int flag_num)
{
  long         i, sg_name_last_char, sg_number;
  kt_atomlage  *WITH;
  vektor       v;


  sg_name_last_char = ' ';

  for (i = 0; i < 30; i++)
  {
    if (strukturdaten->raumgruppe[i] != ' ')
      sg_name_last_char = strukturdaten->raumgruppe[i];
  }

  if (   sg_name_last_char == 'Z' || sg_name_last_char == 'R'
      || sg_name_last_char == 'z' || sg_name_last_char == 'r')
    sg_number = 20000 + rginfo.raumgruppenr;
  else
    sg_number =         rginfo.raumgruppenr;

  fprintf(t->f, "title: %.80s\n", strukturdaten->name);
  fprintf(t->f, "#      %.80s\n", strukturdaten->lit);
  fprintf(t->f, "type: crystal\n");
  fprintf(t->f, "sg_name: %.30s\n", strukturdaten->raumgruppe);
  fprintf(t->f, "sg_number: %ld\n", sg_number);
  fprintf(t->f, "cell_units: angstroms\n");
  fprintf(t->f, "cell_lengths: %.5f %.5f %.5f\n",
          strukturdaten->zellpar.laenge[0],
          strukturdaten->zellpar.laenge[1],
          strukturdaten->zellpar.laenge[2]);
  fprintf(t->f, "cell_angles: %.3f %.3f %.3f\n",
          strukturdaten->zellpar.winkel[0],
          strukturdaten->zellpar.winkel[1],
          strukturdaten->zellpar.winkel[2]);
  fprintf(t->f, "num_atoms: %ld\n", strukturdaten->anzatomlagen);
  if (flag_num == 0)
    fprintf(t->f, "atom_info: symbol frac\n");
  else
    fprintf(t->f, "atom_info: at_number frac\n");
  fprintf(t->f, "atom_list:\n");

  for (i = 1; i <= strukturdaten->anzatomlagen; i++) {
    WITH = strukturdaten->atomlage[i - 1];
    if (flag_num == 0)
    {
                                   putc(WITH->aname[0], t->f);
      if (isalpha(WITH->aname[1])) putc(WITH->aname[1], t->f);
      else                         putc(' ', t->f);
    }
    else
    {
      fprintf(t->f, "%3ld", (long)(50 + i));
    }
    berechne_atom_koordinaten(&WITH->grundatom, v);
    fprintf(t->f, "  %9.5f %9.5f %9.5f # %5.5s\n",
            v[0], v[1], v[2], WITH->aname);
  }

  fprintf(t->f, "stop:\n");
}


Static Void write_focus(_TEXT *t, kt_strukturdaten *strukturdaten)
{
  long         i, c, i_last_char;
  kt_atomlage  *WITH;
  vektor       v;


  print_tlnb(t->f, strukturdaten->name, 80, "Title ", "\n");
  print_tlnb(t->f, strukturdaten->lit,  80, "#     ", "\n");

  fprintf(t->f, "SpaceGroup  ");

  i_last_char = 0;

  for (i = 0; i < 30; i++)
  {
    if (strukturdaten->raumgruppe[i] != ' ')
      i_last_char = i;
  }

  for (i = 0; i <= i_last_char; i++)
  {
        c = strukturdaten->raumgruppe[i];
    if (c != ' ')
    {
      if      (i != i_last_char)
        putc(c, t->f);
      else if (strchr("Zz", c) != NULL)
      {
        putc(':', t->f);
        putc('2', t->f);
      }
      else if (strchr("Rr", c) != NULL)
      {
        putc(':', t->f);
        putc('R', t->f);
      }
      else if (strchr("SHsh", c) == NULL)
        putc(c, t->f);
    }
  }

  putc('\n', t->f);

  fprintf(t->f, "UnitCell  %.5f %.5f %.5f %.3f %.3f %.3f\n",
          strukturdaten->zellpar.laenge[0],
          strukturdaten->zellpar.laenge[1],
          strukturdaten->zellpar.laenge[2],
          strukturdaten->zellpar.winkel[0],
          strukturdaten->zellpar.winkel[1],
          strukturdaten->zellpar.winkel[2]);

  for (i = 1; i <= strukturdaten->anzatomlagen; i++) {
    WITH = strukturdaten->atomlage[i - 1];
    berechne_atom_koordinaten(&WITH->grundatom, v);
    fprintf(t->f, "%.6s  %9.5f %9.5f %9.5f\n",
            WITH->aname, v[0], v[1], v[2]);
  }

  fprintf(t->f, "End\n");
}


Static Void write_GSAS_crystal_structure(_TEXT *t,
                                         kt_strukturdaten *strukturdaten,
                                         int phase_no)
{
  long         i;
  double       vol;
  kt_atomlage  *WITH;
  vektor       v;


  fprintf(t->f, "CRS%d    PNAM  %.6s\n",
          phase_no, strukturdaten->abk);

  fprintf(t->f, "CRS%d   NATOM%5ld\n",
          phase_no, strukturdaten->anzatomlagen);

  fprintf(t->f, "CRS%d  ABC   %10.6f%10.6f%10.6f    N    0\n",
          phase_no,
          strukturdaten->zellpar.laenge[0],
          strukturdaten->zellpar.laenge[1],
          strukturdaten->zellpar.laenge[2]);

  fprintf(t->f, "CRS%d  ANGLES%10.4f%10.4f%10.4f\n",
          phase_no,
          strukturdaten->zellpar.winkel[0],
          strukturdaten->zellpar.winkel[1],
          strukturdaten->zellpar.winkel[2]);

  for (i = 1; i <= strukturdaten->anzatomlagen; i++)
  {
    WITH = strukturdaten->atomlage[i - 1];
    berechne_atom_koordinaten(&WITH->grundatom, v);

    fprintf(t->f,
      "CRS%d  AT%3ldA  %2.2s      %10.6f%10.6f%10.6f%10.6f%6.6s  %4ld 000\n",
      phase_no, i,
      WITH->element,
      v[0], v[1], v[2], 1.,
      WITH->aname,
      WITH->anzatome);

    fprintf(t->f, "CRS%d  AT%3ldB  0.025000%52.52sI\n", phase_no, i, "");
  }

  vol = volumen(strukturdaten->zellpar.laenge[0],
      strukturdaten->zellpar.laenge[1], strukturdaten->zellpar.laenge[2],
      strukturdaten->zellpar.winkel[0], strukturdaten->zellpar.winkel[1],
      strukturdaten->zellpar.winkel[2]);

  fprintf(t->f, "CRS%d  CELVOL%10.3f\n", phase_no, vol);
  fprintf(t->f, "CRS%d  SG SYM  %20.20s\n",
    phase_no, strukturdaten->raumgruppe);
  fprintf(t->f, "CRS%d  SPAXIS    0    0    1\n", phase_no);
}


Static Void write_GSAS_L_A_I(_TEXT *t,
                             kt_strukturdaten *strukturdaten)
{
  long         i;
  kt_atomlage  *WITH;
  vektor       v;


  for (i = 1; i <= strukturdaten->anzatomlagen; i++)
  {
    WITH = strukturdaten->atomlage[i - 1];
    berechne_atom_koordinaten(&WITH->grundatom, v);

    fprintf(t->f,
      "I N %2.2s %10.6f %10.6f %10.6f %10.6f %6.6s I /\n",
      WITH->element,
      v[0], v[1], v[2], 1.,
      WITH->aname);
  }
}


static void wrie_header(FILE *fpout, kt_strukturdaten *struda,
                        char *PreCode, char *SpaceGroup, int *Fxyz)
{
  if (PreCode == NULL)
      PreCode = "";

  if (SpaceGroup == NULL)
      SpaceGroup = struda->raumgruppe;

  fprintf(fpout, "*%s", PreCode);
  print_tlnb(fpout, struda->abk, sizeof (abkbezeichnung), NULL, "\n");

  print_tlnb(fpout, struda->name, 80, NULL, "\n");
  print_tlnb(fpout, struda->lit,  80, NULL, "\n");

  print_tlnb(fpout, SpaceGroup, 30, NULL, "\n");

  if (Fxyz == NULL)
    fprintf(fpout, " %7.4f %7.4f %7.4f %7.3f %7.3f %7.3f\n",
            struda->zellpar.laenge[0], struda->zellpar.laenge[1],
            struda->zellpar.laenge[2], struda->zellpar.winkel[0],
            struda->zellpar.winkel[1], struda->zellpar.winkel[2]);
  else
    fprintf(fpout, " %7.4f %7.4f %7.4f %7.3f %7.3f %7.3f\n",
            struda->zellpar.laenge[0] * Fxyz[0],
            struda->zellpar.laenge[1] * Fxyz[1],
            struda->zellpar.laenge[2] * Fxyz[2],
            struda->zellpar.winkel[0],
            struda->zellpar.winkel[1],
            struda->zellpar.winkel[2]);
}


Static Void write_strukturdaten_karte(_TEXT *t, kt_strukturdaten *struda)
{
  long i;
  kt_atomlage *al;


  wrie_header(t->f, struda, "", NULL, NULL);

  for (i = 0; i < struda->anzatomlagen; i++) {
    al = struda->atomlage[i];
    fprintf(t->f, "%.6s %8.5f %8.5f %8.5f",
            al->aname,
            al->atomlagekoord[0],
            al->atomlagekoord[1],
            al->atomlagekoord[2]);
    if (al->koordinationszahl > 0)
      fprintf(t->f, " %7ld", al->koordinationszahl);
    putc('\n', t->f);
  }

  for (i = 0; i < 72; i++)
    putc('-', t->f);

  putc('\n', t->f);
}


static void write_P_strukturdaten_karte(_TEXT *t, kt_strukturdaten *struda)
{
  int          i;
  long         ial, ia;
  kt_atomlage  *al;


  wrie_header(t->f, struda, "P_", "P 1", NULL);

  for (ial = 0; ial < struda->anzatomlagen; ial++)
  {
    al = struda->atomlage[ial];

    for (ia = 0; ia < al->anzatome; ia++)
    {
      fprintf(t->f, "%.6s %8.5f %8.5f %8.5f",
              al->aname,
              al->atom[ia]->pkoord[0],
              al->atom[ia]->pkoord[1],
              al->atom[ia]->pkoord[2]);

      if (al->koordinationszahl > 0)
        fprintf(t->f, " %7ld", al->koordinationszahl);

      putc('\n', t->f);
    }
  }

  for (i = 0; i < 72; i++)
    putc('-', t->f);

  putc('\n', t->f);
}


static void write_PR_strukturdaten_karte(_TEXT *t, kt_strukturdaten *struda,
                                         int *Fxyz)
{
  int          i, j;
  long         ial, ia, jal, nel;
  kt_atomlage  *al, *alj;
  char         buf[64];
  int          ixyz[3], FxyzBuf[3];
  double       coor[3];


  if (Fxyz == NULL)
  {
    Fxyz = FxyzBuf;
    Fxyz[0] =
    Fxyz[1] =
    Fxyz[2] = 1;
  }

  wrie_header(t->f, struda, "P_", "P 1", Fxyz);

  for (ial = 0; ial < struda->anzatomlagen; ial++)
  {
    al = struda->atomlage[ial];

    nel = 1;

    for (jal = 0; jal < ial; jal++) {
                               alj = struda->atomlage[jal];
      if (strncmp(al->element, alj->element, 2) == 0)
        nel += alj->anzatome * Fxyz[0] * Fxyz[1] * Fxyz[2];
    }

    for (ia = 0; ia < al->anzatome; ia++)
    {
      for (ixyz[0] = 0; ixyz[0] < Fxyz[0]; ixyz[0]++)
      for (ixyz[1] = 0; ixyz[1] < Fxyz[1]; ixyz[1]++)
      for (ixyz[2] = 0; ixyz[2] < Fxyz[2]; ixyz[2]++)
      {
        for (i = j = 0; i < 2; i++)
          if (al->element[i] != ' ') buf[j++] = al->element[i];

        sprintf(&buf[j], "%d", nel++);

        for (i = 0; i < 3; i++)
          coor[i] = al->atom[ia]->pkoord[i] / Fxyz[i] + (double) ixyz[i] / Fxyz[i];

        fprintf(t->f, "%-6s %8.5f %8.5f %8.5f",
                buf, coor[0], coor[1], coor[2]);

        if (al->koordinationszahl > 0)
          fprintf(t->f, " %7ld", al->koordinationszahl);

        putc('\n', t->f);
      }
    }
  }

  for (i = 0; i < 72; i++)
    putc('-', t->f);

  putc('\n', t->f);
}


Static Void write_atomlage_koordinaten(_TEXT *f, kt_atomlage *atomlage)
{
  long i;
  kt_atomart *WITH;

  fprintf(f->f, "\n Coordinates for atom %.6s\n", atomlage->aname);
  for (i = 0; i < atomlage->anzatome; i++) {
    WITH = atomlage->atom[i];
    write_atomart_name(f, atomlage->atom[i]);
    fprintf(f->f, "  ");
    write_vektor(f, WITH->pkoord);
    fprintf(f->f, "    ");
    write_symmetrieoperator(f, &WITH->trilinpol);
    putc('\n', f->f);
  }
}


Static Void read_string80(_TEXT *t, Char *s)
{
  int  c, i;

  i = 0;

  while ((c = getc(t->f)) != EOF && c != '\n')
    if (i < 80) s[i++] = c;

  while (i < 80)
    s[i++] = ' ';
}


Static Void set_zellparameter(double r1, double r2, double r3,
                              double r4, double r5, double r6,
                              zellparameter *z)
{
  z->laenge[0] = r1;
  z->laenge[1] = r2;
  z->laenge[2] = r3;
  z->winkel[0] = r4;
  z->winkel[1] = r5;
  z->winkel[2] = r6;
  load_zellparameter(r1, r2, r3, r4, r5, r6);
}


Static Void read_zellparameter(_TEXT *t,
                               zellparameter *zellpar, Char *kristallsystem)
{
  long i, ii;
  double r[6];

  if (!(P_peek(t->f) == ' ' || P_peek(t->f) == '+' || P_peek(t->f) == '-' ||
        isdigit(P_peek(t->f)))) {
    for (i = 0; i <= 3; i++) {
      kristallsystem[i] = getc(t->f);
      if (kristallsystem[i] == '\n')
        kristallsystem[i] = ' ';
    }
  }
  ii = 0;
  while (!P_eoln(t->f)) {
    ii++;
    fscanf(t->f, "%lg", &r[ii - 1]);
    while ((!P_eoln(t->f)) & (P_peek(t->f) == ' '))
      getc(t->f);
  }
  fscanf(t->f, "%*[^\n]");
  getc(t->f);
  if (ii == 6) {
    set_zellparameter(r[0], r[1], r[2], r[3], r[4], r[5], zellpar);
    return;
  }
  if (!strncmp(kristallsystem, "CUB ", sizeof(st_kristallsystem)) && ii == 1) {
    set_zellparameter(r[0], r[0], r[0], 90.0, 90.0, 90.0, zellpar);
    return;
  }
  if (!strncmp(kristallsystem, "TET ", sizeof(st_kristallsystem)) && ii == 2) {
    set_zellparameter(r[0], r[0], r[1], 90.0, 90.0, 90.0, zellpar);
    return;
  }
  if (!strncmp(kristallsystem, "ORT ", sizeof(st_kristallsystem)) && ii == 3) {
    set_zellparameter(r[0], r[1], r[2], 90.0, 90.0, 90.0, zellpar);
    return;
  }
  if (!strncmp(kristallsystem, "TRIG", sizeof(st_kristallsystem)) && ii == 2) {
    set_zellparameter(r[0], r[0], r[1], 90.0, 90.0, 120.0, zellpar);
    return;
  }
  if (!strncmp(kristallsystem, "HEX ", sizeof(st_kristallsystem)) && ii == 2) {
    set_zellparameter(r[0], r[0], r[1], 90.0, 90.0, 120.0, zellpar);
    return;
  }
  if (!strncmp(kristallsystem, "MON3", sizeof(st_kristallsystem)) && ii == 4) {
    set_zellparameter(r[0], r[1], r[2], r[3], 90.0, 90.0, zellpar);
    return;
  }
  if (!strncmp(kristallsystem, "MON1", sizeof(st_kristallsystem)) && ii == 4) {
    set_zellparameter(r[0], r[1], r[2], 90.0, 90.0, r[3], zellpar);
    return;
  }
  if (!strncmp(kristallsystem, "MON2", sizeof(st_kristallsystem)) && ii == 4) {
    set_zellparameter(r[0], r[1], r[2], 90.0, r[3], 90.0, zellpar);
    return;
  }
  if (!strncmp(kristallsystem, "RHO ", sizeof(st_kristallsystem)) && ii == 2)
    set_zellparameter(r[0], r[0], r[0], r[1], r[1], r[1], zellpar);
  else
    printf(" ERROR  cannot read the cell parameters\n");
}


Static Void suche_naechsten_eintrag(_TEXT *t, Char *abk, boolean *b)
{
  long i;

  *b = false;
  for (i = 0; i <= 15; i++)
    abk[i] = ' ';
  do {
    if (BUFEOF(t->f))
      *b = true;
    else if (P_peek(t->f) == '*')
      *b = true;
    else {
      fscanf(t->f, "%*[^\n]");
      getc(t->f);
    }
  } while (!*b);
  if (BUFEOF(t->f)) {
    *b = false;
    return;
  }
  getc(t->f);
  i = 0;
  while (!P_eoln(t->f) && i < 16) {
    i++;
    abk[i - 1] = getc(t->f);
    if (abk[i - 1] == '\n')
      abk[i - 1] = ' ';
  }
  fscanf(t->f, "%*[^\n]");
  getc(t->f);
  if (abk[15] != ' ')
    printf(" WARNING  entry name on file STRUDAT too long\n");
  if (!strncmp(abk, "END             ", sizeof(abkbezeichnung)) ||
      !strncmp(abk, "end             ", sizeof(abkbezeichnung)))
    *b = false;
}


Static Void suche_strukturdaten_eintrag(_TEXT *t, Char *gegabk, boolean *b)
{
  abkbezeichnung abk;

  do {
    suche_naechsten_eintrag(t, abk, b);
  } while (strncmp(abk, gegabk, sizeof(abkbezeichnung)) && *b);
  if (*b)
    printf(" ... entry %.16s found\n", abk);
  else
    printf(" ERROR  entry %.16s not found on file STRUDAT\n", gegabk);
}


Static Void erzeuge_atomlage(Char *e, long n, double *koord, long kz,
                             kt_atomlage **a)
{
  kt_atomlage *WITH;

  *a = (kt_atomlage *)Malloc(sizeof(kt_atomlage));
  WITH = *a;
  memcpy(WITH->atomlagekoord, koord, sizeof(vektor));
  memcpy(WITH->element, e, sizeof(kt_element));
  memcpy(WITH->aname, "      ", sizeof(kt_atomlagename)); // could expand this to make longer names?
  WITH->aname[0] = e[0];
  WITH->aname[1] = e[1];
  WITH->koordinationszahl = kz;
  WITH->dist_entry = NULL;

  if (n > 0) {
    if (n < 10)
      WITH->aname[2] = (Char)(n - 1 + '1');
    else if (n < 100) { // modified 04/07/2012 - StefS - added check for 100 atoms and 'else' to support up to 999 for the atom label
      WITH->aname[2] = (Char)(n / 10 - 1 + '1');
      WITH->aname[3] = (Char)(n % 10 - 1 + '1');
/* p2c: kriber.pas, line 1071: Note: Using % for possibly-negative arguments [317]
 * Note: Using % for possibly-negative arguments [317] */
    }
    else {
      WITH->aname[2] = (Char)(n / 100 - 1 + '1');
      WITH->aname[3] = (Char)( (n % 100) / 10  - 1 + '1');
      WITH->aname[4] = (Char)(n % 10  - 1 + '1');
    }


  }
  if (WITH->aname[1] == ' ' && WITH->aname[2] != ' ') {
    WITH->aname[1] = WITH->aname[2];
    WITH->aname[2] = WITH->aname[3];
    WITH->aname[3] = WITH->aname[4];
    WITH->aname[4] = ' ';
  }
  berechne_atomlagekoordinaten(&rginfo, *a);
  (*a)->grundatom.elementarzelle[0] = 0;
  (*a)->grundatom.elementarzelle[1] = 0;
  (*a)->grundatom.elementarzelle[2] = 0;
  (*a)->grundatom.dessen_atomart = (*a)->atom[0];
}


Static Void lade_atomlage(kt_atomlage **neue_atomlage,
                          kt_strukturdaten *struda,
                          boolean *b)
{
  long i, j;
  kt_atomlage *WITH1;
  long FORLIM, FORLIM1;
  _TEXT TEMP;

  if (struda->anzatomlagen < maxanzahlatomlagen) {
    *b = true;
    WITH1 = *neue_atomlage;
    FORLIM = struda->anzatomlagen;
    for (i = 0; i < FORLIM; i++) {
      FORLIM1 = struda->atomlage[i]->anzatome;
      for (j = 0; j < FORLIM1; j++) {
        if (gleiche_vektoren(WITH1->atomlagekoord,
                             struda->atomlage[i]->atom[j]->pkoord)) {
          printf(" ERROR  atom position %.6s has the same coordinates as atom",
                 WITH1->aname);
          TEMP.f = stdout;
          *TEMP.name = '\0';
          write_atomart_name(&TEMP, struda->atomlage[i]->atom[j]);
          printf("\n        and therefore has not been loaded\n");
          *b = false;
        }
      }
      if (!strncmp(WITH1->aname, struda->atomlage[i]->aname,
                   sizeof(kt_atomlagename)))
        printf(
          " WARNING  atom position %.6s has same name as another atom position\n",
          WITH1->aname);
    }
    if (!*b)
      return;
    struda->anzatomlagen++;
    struda->atomlage[struda->anzatomlagen - 1] = *neue_atomlage;
    struda->atomlage[struda->anzatomlagen - 1]->atomlagenr = struda->anzatomlagen;
    return;
  }
  *b = false;
  printf(" ERROR  too many atom positions. New atom position cannot be loaded.\n");
  printf("        The program is compiled for a maximum of %ld atom positions.\n",
         (long)maxanzahlatomlagen);
  printf("        Increase constant \"maxanzahlatomlagen\" and compile again.\n");
}


Static Void read_atomlage(_TEXT *t, kt_atomlage **atomlage)
{
  kt_element element;
  long atom_label, koordinationszahl;
  vektor koord;

  if (isalpha(P_peek(t->f))) {
    *element = getc(t->f);
    if (element[0] == '\n')
      element[0] = ' ';
  } else
    element[0] = ' ';
  if (isalpha(P_peek(t->f))) {
    element[1] = getc(t->f);
    if (element[1] == '\n')
      element[1] = ' ';
  } else
    element[1] = ' ';
  if (P_peek(t->f) == '0' || P_peek(t->f) >= '1' && P_peek(t->f) <= '9')
    fscanf(t->f, "%ld", &atom_label);
  else
    atom_label = 0;
  fscanf(t->f, "%lg%lg%lg", koord, &koord[1], &koord[2]);
  while ((!P_eoln(t->f)) & (P_peek(t->f) == ' '))
    getc(t->f);
  if (!P_eoln(t->f))
    fscanf(t->f, "%ld", &koordinationszahl);
  else {
    koordinationszahl = 0;
    if (!strncmp(element, "SI", sizeof(kt_element)) ||
        !strncmp(element, "Si", sizeof(kt_element)))
      koordinationszahl = 4;
    if (!strncmp(element, "AL", sizeof(kt_element)) ||
        !strncmp(element, "Al", sizeof(kt_element)))
      koordinationszahl = 4;
    if (!strncmp(element, "P ", sizeof(kt_element)))
      koordinationszahl = 4;
    if (!strncmp(element, "GA", sizeof(kt_element)) ||
        !strncmp(element, "Ga", sizeof(kt_element)))
      koordinationszahl = 4;
    if (!strncmp(element, "ZN", sizeof(kt_element)) ||
        !strncmp(element, "Zn", sizeof(kt_element)))
      koordinationszahl = 4;
    if (!strncmp(element, "T ", sizeof(kt_element)))
      koordinationszahl = 4;
    if (!strncmp(element, "O ", sizeof(kt_element)))
      koordinationszahl = 2;
  }
  fscanf(t->f, "%*[^\n]");
  getc(t->f);

  // ATOM LABEL is correct here -> Stefs
  // printf("%d\n", atom_label);

  erzeuge_atomlage(element, atom_label, koord, koordinationszahl, atomlage);
}


Static Void read_strukturdaten(_TEXT *t,
                               kt_strukturdaten *strukturdaten,
                               boolean frage)
{
  kt_atomlage *neue_atomlage;
  boolean b;

  if (frage)
    printf(" Name: ");
  read_string80(t, strukturdaten->name);
  if (frage)
    printf(" Lit:  ");
  read_string80(t, strukturdaten->lit);
  if (frage)
    printf(" Space group: ");
  readln_raumgruppename(t, strukturdaten->raumgruppe);
  read_raumgruppe(&symdat, strukturdaten->raumgruppe, &rginfo);
  memcpy(strukturdaten->kristallsystem, rginfo.kristallsystem,
         sizeof(st_kristallsystem));
  if (frage)
    printf(" Cell parameters: ");
  read_zellparameter(t, &strukturdaten->zellpar,
                     strukturdaten->kristallsystem);
  if (frage) {
    printf(" Add atom position,  end with '*' \n");
    printf(" label (element symbol + number)   coordinates   coordination number\n");
  }
  strukturdaten->anzatomlagen = 0;
  if (frage)
    printf(" Atom%3ld: ", strukturdaten->anzatomlagen + 1);
  while (isalpha(P_peek(t->f))) {
    read_atomlage(t, &neue_atomlage);
    frac_vektor(neue_atomlage->atomlagekoord, neue_atomlage->atomlagekoord);
    lade_atomlage(&neue_atomlage, strukturdaten, &b);
    if (frage)
      printf(" Atom%3ld: ", strukturdaten->anzatomlagen + 1);
  }
  if (frage) {
    fscanf(t->f, "%*[^\n]");
    getc(t->f);
  }
}


Static boolean atom_hat_gleiches_element(kt_atom *a, Char *e)
{
  kt_atomlage *WITH;

  WITH = a->dessen_atomart->deren_atomlage;
  return (WITH->element[0] == e[0] && WITH->element[1] == e[1]);
}


Static boolean symmetrieaequivalente_atome(kt_atom *a, kt_atom *b)
{
  return (   a->dessen_atomart->deren_atomlage
          == b->dessen_atomart->deren_atomlage);
}


Static Void welchsymop(kt_atom *a, kt_atom *b, st_raumgruppe *rginfo,
                       long *anzsymop, symmetrieoperator *sr)
{
  long i, j;
  symmetrieoperator sa, sb, sc;
  gittervektor gb, gc;

  *anzsymop = 0;
  if (!symmetrieaequivalente_atome(a, b))
    return;
  berechne_atom_trilinpol(a, &sa);
  berechne_atom_trilinpol(b, &sb);
  int_vektor(sb.t, gb);
  frac_vektor(sb.t, sb.t);
  for (i = 0; i < rginfo->zaehligkeit; i++) {
    mult_symmetrieoperatoren(rginfo->symop[i], &sa, &sc);
    int_vektor(sc.t, gc);
    frac_vektor(sc.t, sc.t);
    if (gleiche_symmetrieoperatoren(&sc, &sb)) {
      (*anzsymop)++;
      sr[*anzsymop - 1] = *rginfo->symop[i];
      for (j = 0; j <= 2; j++)
        sr[*anzsymop - 1].t[j] += gb[j] - gc[j];
    }
  }
}


Static boolean sym_aequival_distanzen(kt_atom *a1, kt_atom *a2,
                                      kt_atom *a3, kt_atom *a4)
{
  st_symop_array symop;
  long anzsymop, i;
  kt_atom a3n;
  boolean b;

  b = false;
  welchsymop(a2, a4, &rginfo, &anzsymop, symop);
  for (i = 0; i < anzsymop; i++) {
    alg_symmetrieoperation(&symop[i], a1, &a3n);
    if (gleiche_atome(&a3n, a3))
      b = true;
  }
  return b;
}


Static boolean sym_aequival_winkel(kt_atom *a1, kt_atom *a2, kt_atom *a3,
                                   kt_atom *a4, kt_atom *a5, kt_atom *a6)
{
  st_symop_array symop;
  long anzsymop, i;
  kt_atom a4n, a6n;
  boolean b;

  b = false;
  welchsymop(a2, a5, &rginfo, &anzsymop, symop);
  for (i = 0; i < anzsymop; i++) {
    alg_symmetrieoperation(&symop[i], a1, &a4n);
    alg_symmetrieoperation(&symop[i], a3, &a6n);
    if (gleiche_atome(&a4n, a4) & gleiche_atome(&a6n, a6))
      b = true;
  }
  return b;
}


Static Void delete_atomlage(kt_atomlage *atomlage, kt_strukturdaten *struda)
{
  long j, n, FORLIM;

  n = atomlage->atomlagenr;
  atomlage->atomlagenr = 0;
  struda->anzatomlagen--;
  FORLIM = struda->anzatomlagen;
  for (j = n; j <= FORLIM; j++) {
    struda->atomlage[j - 1] = struda->atomlage[j];
    struda->atomlage[j - 1]->atomlagenr = j;
  }
}


Static Void read_atomlage_bez(_TEXT *f, kt_atomlage **pala)
{
  long i;
  Char nn[6];
  long FORLIM;

  *pala = NULL;
  while ((!P_eoln(f->f)) & (P_peek(f->f) == ' '))
    getc(f->f);
  if (!P_eoln(f->f)) {
    for (i = 0; i <= 5; i++)
      nn[i] = ' ';
    i = 0;
    while ((P_peek(f->f) != ' ') & (P_peek(f->f) != '(') & (!P_eoln(f->f))) {
      i++;
      nn[i - 1] = getchar();
      if (nn[i - 1] == '\n')
        nn[i - 1] = ' ';
    }
    FORLIM = strukturdaten.anzatomlagen;
    for (i = 0; i < FORLIM; i++) {
      if (!strncmp(nn, strukturdaten.atomlage[i]->aname, 6)) {
/* p2c: kriber.pas, line 1284: Warning: Incompatible array sizes [164] */
        *pala = strukturdaten.atomlage[i];
      }
    }
  }
  if (*pala == NULL)
    printf(" ERROR  atom position %.6s not found\n", nn);
}


Static Void readln_atom(_TEXT *t, kt_atom *a)
{
  symmetrieoperator ss, ssy;
  long i, j, n;
  kt_atomlage *pala;
  gittervektor v;
  boolean b;
  long FORLIM;
  _TEXT TEMP;

  a->dessen_atomart = NULL;
  read_atomlage_bez(t, &pala);
  if (pala != NULL) {
    a->dessen_atomart = pala->atom[0];
    a->elementarzelle[0] = 0;
    a->elementarzelle[1] = 0;
    a->elementarzelle[2] = 0;
    while ((!P_eoln(t->f)) & (P_peek(t->f) == ' '))
      getc(t->f);
    if (!P_eoln(t->f)) {
      if (P_peek(t->f) == '(') {
        while ((!P_eoln(t->f)) & ((P_peek(t->f) == ' ') | (P_peek(t->f) ==
                                    ')') | (P_peek(t->f) == '(')))
          getc(t->f);
        if (!P_eoln(t->f))
          scanf("%ld", &n);
        if (n > 0 && n <= a->dessen_atomart->deren_atomlage->anzatome) {
          a->dessen_atomart = pala->atom[n - 1];
          while ((!P_eoln(t->f)) & ((P_peek(t->f) == ' ') | (P_peek(t->f) ==
                                      ')') | (P_peek(t->f) == '(')))
            getc(t->f);
         
          if (!P_eoln(t->f))
            fscanf(t->f, "%ld%ld%ld", a->elementarzelle,
                   &a->elementarzelle[1], &a->elementarzelle[2]);
        } else {
          printf(" ERROR  there is no atom ");
          printf("%.6s (%ld)\n", a->dessen_atomart->deren_atomlage->aname, n);
          a->dessen_atomart = NULL;
        }
      } else {
        read_symmetrieoperator(t, &ss);
        int_vektor(ss.t, a->elementarzelle);
        frac_vektor(ss.t, ss.t);
        b = true;
        FORLIM = pala->anzatome;
        for (i = 0; i < FORLIM; i++) {
          ssy = pala->atom[i]->trilinpol;
          int_vektor(ssy.t, v);
          frac_vektor(ssy.t, ssy.t);
          if (gleiche_symmetrieoperatoren(&ss, &ssy) && b) {
            a->dessen_atomart = pala->atom[i];
            for (j = 0; j <= 2; j++)
              a->elementarzelle[j] -= v[j];
            b = false;
          }
        }
        if (b) {
          printf(" ERROR  atom %.6s with symmetry operator ",
                 a->dessen_atomart->deren_atomlage->aname);
          TEMP.f = stdout;
          *TEMP.name = '\0';
          write_symmetrieoperator(&TEMP, &ss);
          printf(" not found\n");
          a->dessen_atomart = NULL;
        }
      }
    }
  }
  if (a->dessen_atomart == NULL)
    printf(" ERROR  atom not found\n");
  fscanf(t->f, "%*[^\n]");
  getc(t->f);
}


Static long koordinationszahl_des_atoms(kt_atom *a)
{
  return (a->dessen_atomart->deren_atomlage->koordinationszahl);
}


Static Void verschiebe_atome_leicht(kt_atomlage *atomlage)
{
  long i;
  double ran;
  long FORLIM;

  ran = 0.23456789;
  ran = fabs(ran + atomlage->atomlagekoord[0]) -
        (long)fabs(ran + atomlage->atomlagekoord[0]);
  atomlage->atomlagekoord[0] += (2 * ran - 1) / 1000;
  ran = fabs(ran + atomlage->atomlagekoord[1]) -
        (long)fabs(ran + atomlage->atomlagekoord[1]);
  atomlage->atomlagekoord[1] += (2 * ran - 1) / 1000;
  ran = fabs(ran + atomlage->atomlagekoord[2]) -
        (long)fabs(ran + atomlage->atomlagekoord[2]);
  atomlage->atomlagekoord[2] += (2 * ran - 1) / 1000;
  mult_symop_mit_vektor(atomlage->symmetrielage, atomlage->atomlagekoord,
                        atomlage->atomlagekoord);
  FORLIM = atomlage->anzatome;
  for (i = 0; i < FORLIM; i++)
    mult_symop_mit_vektor(&atomlage->atom[i]->trilinpol,
                          atomlage->atomlagekoord, atomlage->atom[i]->pkoord);
}


Static Void berechne_distanzen_algebraisch(st_raumgruppe *rginfo)
{
  long i;
  symmetrieoperator s;
  _TEXT TEMP;

  printf("  Atom - Atom      Distance\n");
  for (i = 2; i <= rginfo->zaehligkeit; i++) {
    subtrahiere_symmetrieoperatoren(rginfo->symop[0], rginfo->symop[i - 1], &s);
    printf("   1     %2ld     ", i);
    TEMP.f = stdout;
    *TEMP.name = '\0';
    write_symmetrieoperator(&TEMP, &s);
    putchar('\n');
  }
}


Static Void berechne_distanzen(kt_atom *zatom, kt_strukturdaten *strukturdaten,
                               long max_anzahl, double max_laenge)
{
  long ii, i, j, k, m, n, o, me, ne, oe;
  vektor vv;
  boolean b;
  distanzen_tabellen_eintrag *hh, *neue_distanz, *neue_dist_tab;
  double vol, TEMP;
  kt_atomlage *WITH;
  long FORLIM1;
  kt_atom *WITH1;

  vol = volumen(strukturdaten->zellpar.laenge[0],
      strukturdaten->zellpar.laenge[1], strukturdaten->zellpar.laenge[2],
      strukturdaten->zellpar.winkel[0], strukturdaten->zellpar.winkel[1],
      strukturdaten->zellpar.winkel[2]);
  TEMP = cos(strukturdaten->zellpar.winkel[0] * 0.0174533);
  me = (long)(max_laenge * strukturdaten->zellpar.laenge[1] *
         strukturdaten->zellpar.laenge[2] * sqrt(1 - TEMP * TEMP) / vol) + 1;
  TEMP = cos(strukturdaten->zellpar.winkel[1] * 0.0174533);
  ne = (long)(max_laenge * strukturdaten->zellpar.laenge[0] *
         strukturdaten->zellpar.laenge[2] * sqrt(1 - TEMP * TEMP) / vol) + 1;
  TEMP = cos(strukturdaten->zellpar.winkel[2] * 0.0174533);
  oe = (long)(max_laenge * strukturdaten->zellpar.laenge[0] *
         strukturdaten->zellpar.laenge[1] * sqrt(1 - TEMP * TEMP) / vol) + 1;
  neue_dist_tab = NULL;
  neue_distanz = (distanzen_tabellen_eintrag *)
                 Malloc(sizeof(distanzen_tabellen_eintrag));
  max_laenge *= max_laenge;
  ii = 0;
  for (j = 0; j < strukturdaten->anzatomlagen; j++) {
    WITH = strukturdaten->atomlage[j];
    FORLIM1 = WITH->anzatome;
    for (k = 0; k < FORLIM1; k++) {
      for (m = -me; m <= me; m++) {
        for (n = -ne; n <= ne; n++) {
          for (o = -oe; o <= oe; o++) {
            vv[0] =   WITH->atom[k]->pkoord[0]
                    + m - zatom->dessen_atomart->pkoord[0];
            vv[1] =   WITH->atom[k]->pkoord[1]
                    + n - zatom->dessen_atomart->pkoord[1];
            vv[2] =   WITH->atom[k]->pkoord[2]
                    + o - zatom->dessen_atomart->pkoord[2];
            neue_distanz->laenge = skalar_produkt(vv, vv);
            WITH1 = &neue_distanz->endatom;
            WITH1->elementarzelle[0] = m;
            WITH1->elementarzelle[1] = n;
            WITH1->elementarzelle[2] = o;
            WITH1->dessen_atomart = WITH->atom[k];
            if (!gleiche_atome(zatom, &neue_distanz->endatom) &&
                neue_distanz->laenge <= max_laenge) {
              b = false;
              if (neue_dist_tab != NULL) {
                if (neue_dist_tab->laenge < neue_distanz->laenge) {
                  hh = neue_dist_tab;
                  do {
                    if (hh->naechste_distanz == NULL)
                      b = true;
                    else if (hh->naechste_distanz->laenge > neue_distanz->laenge)
                      b = true;
                    else
                      hh = hh->naechste_distanz;
                  } while (!b);
                  neue_distanz->naechste_distanz = hh->naechste_distanz;
                  hh->naechste_distanz = neue_distanz;
                }
              }
              if (b == false) {
                neue_distanz->naechste_distanz = neue_dist_tab;
                neue_dist_tab = neue_distanz;
              }
              ii++;
              if (ii > max_anzahl) {
                hh = neue_dist_tab;
                for (i = 2; i <= max_anzahl; i++)
                  hh = hh->naechste_distanz;
                neue_distanz = hh->naechste_distanz;
                hh->naechste_distanz = NULL;
              } else
                neue_distanz = (distanzen_tabellen_eintrag *)
                               Malloc(sizeof(distanzen_tabellen_eintrag));
            }
          }
        }
      }
    }
  }
  hh = neue_dist_tab;
  while (hh != NULL) {
    hh->laenge = sqrt(hh->laenge);
    hh = hh->naechste_distanz;
  }
  zatom->dessen_atomart->deren_atomlage->dist_entry = neue_dist_tab;
}


Static Void write_distanzen(_TEXT *t, kt_atom *a)
{
  distanzen_tabellen_eintrag *hh;
  long i, anzs;
  st_symop_array symop;
  kt_atom aa;
  distanzen_tabellen_eintrag *WITH;

  fprintf(t->f, "\n Distances from atom ");
  write_atom_name(t, a);
  fprintf(t->f, "  to\n");
  welchsymop(&a->dessen_atomart->deren_atomlage->grundatom, a, &rginfo, &anzs,
             symop);
  hh = a->dessen_atomart->deren_atomlage->dist_entry;
  i = 0;
  while (hh != NULL) {
    WITH = hh;
    i++;
    fprintf(t->f, " (%ld)   -", i);
    if (anzs > 0)
      alg_symmetrieoperation(&symop[0], &WITH->endatom, &aa);
    else
      aa = WITH->endatom;
    write_atom_name(t, &aa);
    fprintf(t->f, "%11.4f\n", WITH->laenge);
    hh = WITH->naechste_distanz;
  }
}


Static Void write_distanzen_uniq(_TEXT *t, kt_atom *a)
{
  distanzen_tabellen_eintrag *hh;
  long anzs;
  st_symop_array symop;
  kt_atom aa;

  welchsymop(&a->dessen_atomart->deren_atomlage->grundatom, a, &rginfo, &anzs,
             symop);

  hh = a->dessen_atomart->deren_atomlage->dist_entry;

  while (hh != NULL)
  {
    if (anzs > 0)
      alg_symmetrieoperation(&symop[0], &hh->endatom, &aa);
    else
      aa = hh->endatom;

    if (   aa.dessen_atomart->deren_atomlage->atomlagenr
        >= a->dessen_atomart->deren_atomlage->atomlagenr)
      putc(' ', t->f);
    else
      putc('#', t->f);

    write_atom_name(t, a); fprintf(t->f, " -");

    write_atom_name(t, &aa);
    fprintf(t->f, "%11.4f\n", hh->laenge);

    hh = hh->naechste_distanz;
  }
}


Static Void berechne_nte_distanz(kt_atom *a, long n, double *d, kt_atom *b)
{
  distanzen_tabellen_eintrag *hh;
  long i;

  b->dessen_atomart = NULL;
  *d = 0.0;
  hh = a->dessen_atomart->deren_atomlage->dist_entry;
  i = 0;
  while (hh != NULL && i < n) {
    i++;
    *b = hh->endatom;
    *d = hh->laenge;
    hh = hh->naechste_distanz;
  }
  if (i != n)
    b->dessen_atomart = NULL;
  if (b->dessen_atomart == NULL)
    printf(" ERROR  distance no. %ld has not been calculated\n", n);
}


Static double distanz_zwischen_atome(kt_atom *a, kt_atom *b)
{
  vektor v, va;
  long i;

  berechne_atom_koordinaten(a, va);
  berechne_atom_koordinaten(b, v);
  for (i = 0; i <= 2; i++)
    v[i] -= va[i];
  return sqrt(skalar_produkt(v, v));
}


Static double winkel_zwischen_atome(kt_atom *a, kt_atom *b, kt_atom *c)
{
  vektor va, vb, vc;
  long i;
  double cn, cz, cosw, sin2w;

  berechne_atom_koordinaten(a, va);
  berechne_atom_koordinaten(b, vb);
  berechne_atom_koordinaten(c, vc);
  for (i = 0; i <= 2; i++)
    va[i] -= vb[i];
  for (i = 0; i <= 2; i++)
    vc[i] -= vb[i];
  cn = skalar_produkt(va, vc);
  cz = sqrt(skalar_produkt(va, va) * skalar_produkt(vc, vc));
  if (cz == 0)
    cosw = 1.0;
  else
    cosw = cn / cz;
  if (cosw > 1.0)
    cosw = 1.0;
  if (fabs(cosw) < 0.0000000001)
    return 90.0;
  else {
    sin2w = 1 - cosw * cosw;
    if (sin2w < 0)
      sin2w = 0.0;
    if (cosw > 0)
      return (atan(fabs(sqrt(sin2w) / cosw)) * 57.29577951);
    else
      return (180 - atan(fabs(sqrt(sin2w) / cosw)) * 57.29577951);
  }
}


Static Void write_f_distanzen_und_winkel(_TEXT *t, kt_atom *a)
{
  distanzen_tabellen_eintrag *hh1, *hh2;
  long i;
  double w;

  fprintf(t->f, " Distances and angles for atom ");
  write_atom_name(t, a);
  fprintf(t->f,
    "\n                                       (1)      (2)      (3)      (4)      (5)\n");
  i = 0;
  hh1 = a->dessen_atomart->deren_atomlage->dist_entry;
  while (hh1 != NULL && i < 6) {
    i++;
    fprintf(t->f, " (%ld)   -", i);
    write_atom_name(t, &hh1->endatom);
    fprintf(t->f, "%8.4f", hh1->laenge);
    hh2 = a->dessen_atomart->deren_atomlage->dist_entry;
    while (hh2 != hh1) {
      w = winkel_zwischen_atome(&hh1->endatom, a, &hh2->endatom);
      fprintf(t->f, "%9.3f", w);
      hh2 = hh2->naechste_distanz;
    }
    putc('\n', t->f);
    hh1 = hh1->naechste_distanz;
  }
  putc('\n', t->f);
}


Static Void write_u_distanzen_und_winkel(_TEXT *t, kt_atom *a)
{
  double  w;
  int     i;
  distanzen_tabellen_eintrag *hh1, *hh2;


  i = 0;

  for (hh1 = a->dessen_atomart->deren_atomlage->dist_entry;
       hh1 != NULL;
       hh1 = hh1->naechste_distanz)
  {
    if (i++ >= 6) break;

    fprintf(t->f, "%4.4s %4.4s %8.4f\n",
                a->dessen_atomart->deren_atomlage->aname,
      hh1->endatom.dessen_atomart->deren_atomlage->aname,
      hh1->laenge);
  }

  i = 0;

  for (hh1 = a->dessen_atomart->deren_atomlage->dist_entry;
       hh1 != NULL;
       hh1 = hh1->naechste_distanz)
  {
    if (i++ >= 6) break;

    for (hh2 = a->dessen_atomart->deren_atomlage->dist_entry;
         hh2 != hh1;
         hh2 = hh2->naechste_distanz)
    {
      w = winkel_zwischen_atome(&hh1->endatom, a, &hh2->endatom);

      fprintf(t->f, "%4.4s %4.4s %4.4s %9.3f\n",
        hh1->endatom.dessen_atomart->deren_atomlage->aname,
                  a->dessen_atomart->deren_atomlage->aname,
        hh2->endatom.dessen_atomart->deren_atomlage->aname, w);
    }
  }
}


Static Void write_distanzen_und_vektoren(_TEXT *t, kt_atom *a)
{
  distanzen_tabellen_eintrag *hh;
  long i, j;
  vektor v, va;

  fprintf(t->f, " Distances and interatomic vectors for atom ");
  write_atom_name(t, a);
  putc('\n', t->f);
  i = 0;
  hh = a->dessen_atomart->deren_atomlage->dist_entry;
  while (hh != NULL) {
    i++;
    fprintf(t->f, " (%ld)   -", i);
    write_atom_name(t, &hh->endatom);
    fprintf(t->f, "%8.4f", hh->laenge);
    berechne_atom_koordinaten(a, va);
    berechne_atom_koordinaten(&hh->endatom, v);
    for (j = 0; j <= 2; j++)
      v[j] -= va[j];
    fprintf(t->f, "     ");
    write_vektor(t, v);
    putc('\n', t->f);
    hh = hh->naechste_distanz;
  }
  putc('\n', t->f);
}


Static Void berechne_nte_bindung(kt_atom *x, long n, kt_atom *y)
{
  *y = x->dessen_atomart->bindung[n - 1];
  y->elementarzelle[0] += x->elementarzelle[0];
  y->elementarzelle[1] += x->elementarzelle[1];
  y->elementarzelle[2] += x->elementarzelle[2];
}


Static long anzahl_bindungen_des_atoms(kt_atom *a)
{
  long j;
  kt_atomart *WITH;

  WITH = a->dessen_atomart;
  if (WITH->bindung[maxanzahlbindungen - 1].dessen_atomart != NULL) {
    j = maxanzahlbindungen;
    return j;
  }
  j = 0;
  do {
    j++;
  } while (WITH->bindung[j - 1].dessen_atomart != NULL);
  j--;
  return j;
}


Static boolean atome_miteinander_verbunden(kt_atom *a1, kt_atom *a2)
{
  long i;
  boolean b;
  kt_atom a3;
  long FORLIM;

  b = false;
  FORLIM = anzahl_bindungen_des_atoms(a1);
  for (i = 1; i <= FORLIM; i++) {
    berechne_nte_bindung(a1, i, &a3);
    if (gleiche_atome(a2, &a3))
      b = true;
  }
  return b;
}


Static Void erzeuge_bindung(kt_atom *a1, kt_atom *a2)
{
  long i, j, j_a1, j_a2, k;
  kt_atom b1, b2;
  boolean zuwenig_bindungen1, zuwenig_bindungen2;
  _TEXT TEMP;
  long FORLIM;
  kt_atomart *WITH;
  kt_atomlage *WITH1;

  if (atome_miteinander_verbunden(a1, a2)) {
    printf(" WARNING  there is already a bond ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    write_atom_name(&TEMP, a1);
    printf("  -- ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    write_atom_name(&TEMP, a2);
    putchar('\n');
    return;
  }
  zuwenig_bindungen1 = false;
  zuwenig_bindungen2 = false;
  j_a1 = anzahl_bindungen_des_atoms(a1);
  j_a2 = anzahl_bindungen_des_atoms(a2);
  FORLIM = rginfo.zaehligkeit;
  for (i = 0; i < FORLIM; i++) {
    alg_symmetrieoperation(rginfo.symop[i], a1, &b1);
    alg_symmetrieoperation(rginfo.symop[i], a2, &b2);
    for (k = 0; k <= 2; k++)
      b2.elementarzelle[k] -= b1.elementarzelle[k];
    j = j_a1;
    WITH = b1.dessen_atomart;
    do {
      j++;
      if (WITH->bindung[j - 1].dessen_atomart == NULL)
        WITH->bindung[j - 1] = b2;
    } while (!(gleiche_atome(&WITH->bindung[j - 1], &b2) ||
               j == maxanzahlbindungen));
    if (!gleiche_atome(&WITH->bindung[j - 1], &b2))
      zuwenig_bindungen1 = true;
    for (k = 0; k <= 2; k++)
      b1.elementarzelle[k] = -b2.elementarzelle[k];
    j = j_a2;
    WITH = b2.dessen_atomart;
    do {
      j++;
      if (WITH->bindung[j - 1].dessen_atomart == NULL)
        WITH->bindung[j - 1] = b1;
    } while (!(gleiche_atome(&WITH->bindung[j - 1], &b1) ||
               j == maxanzahlbindungen));
    if (!gleiche_atome(&WITH->bindung[j - 1], &b1))
      zuwenig_bindungen2 = true;
  }
  if (zuwenig_bindungen1) {
    printf(" ERROR  too many bonds for atom ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    write_atom_name(&TEMP, a1);
    printf("\n        The program is compiled for a maximum of %ld bonds.\n",
           (long)maxanzahlbindungen);
    printf("        Increase constant \"maxanzahlbindungen\" and compile again.\n");
  }
  if (zuwenig_bindungen2) {
    printf(" ERROR  too many bonds for ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    write_atom_name(&TEMP, a2);
    printf("\n        The program is compiled for a maximum of %ld bonds.\n",
           (long)maxanzahlbindungen);
    printf("        Increase constant \"maxanzahlbindungen\" and compile again.\n");
  }
  if (zuwenig_bindungen1 || zuwenig_bindungen2) {
    printf(" ERROR  bond ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    write_atom_name(&TEMP, a1);
    printf("  -- ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    write_atom_name(&TEMP, a2);
    printf(" not created \n");
    WITH1 = a1->dessen_atomart->deren_atomlage;
    FORLIM = WITH1->anzatome;
    for (j = 0; j < FORLIM; j++) {
      for (k = j_a1; k < maxanzahlbindungen; k++)
        WITH1->atom[j]->bindung[k].dessen_atomart = NULL;
    }
    WITH1 = a2->dessen_atomart->deren_atomlage;
    FORLIM = WITH1->anzatome;
    for (j = 0; j < FORLIM; j++) {
      for (k = j_a2; k < maxanzahlbindungen; k++)
        WITH1->atom[j]->bindung[k].dessen_atomart = NULL;
    }
  } else {
    printf(" ... new bond ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    write_atom_name(&TEMP, a1);
    printf("  -- ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    write_atom_name(&TEMP, a2);
    printf(" created \n");
  }
  if (anzahl_bindungen_des_atoms(a1) > koordinationszahl_des_atoms(a1))
    printf(" WARNING  too many bonds for %.6s\n",
           a1->dessen_atomart->deren_atomlage->aname);
  if (anzahl_bindungen_des_atoms(a2) > koordinationszahl_des_atoms(a2)) {
    if (   a1->dessen_atomart->deren_atomlage
        != a2->dessen_atomart->deren_atomlage)
      printf(" WARNING  too many bonds for %.6s\n",
             a2->dessen_atomart->deren_atomlage->aname);
  }
}


Static Void loesche_bindung(kt_atom *a1, kt_atom *a2)
{
  kt_atom aa;
  boolean b1[maxanzahlbindungen], b2[maxanzahlbindungen];
  long i, j, k, FORLIM;
  kt_atomlage *WITH;
  long FORLIM1;
  kt_atomart *WITH1;
  _TEXT TEMP;

  if (atome_miteinander_verbunden(a1, a2)) {
    for (i = 0; i < maxanzahlbindungen; i++) {
      berechne_nte_bindung(a1, i + 1, &aa);
      b1[i] = false;
      if (aa.dessen_atomart != NULL) {
        if ((!sym_aequival_distanzen(a1, a2, a1, &aa)) &
            (!sym_aequival_distanzen(a2, a1, a1, &aa)))
          b1[i] = true;
      }
      berechne_nte_bindung(a2, i + 1, &aa);
      b2[i] = false;
      if (aa.dessen_atomart != NULL) {
        if ((!sym_aequival_distanzen(a1, a2, a2, &aa)) &
            (!sym_aequival_distanzen(a2, a1, a2, &aa)))
          b2[i] = true;
      }
    }
    k = 0;
    WITH = a1->dessen_atomart->deren_atomlage;
    for (i = 0; i < maxanzahlbindungen; i++) {
      if (b1[i]) {
        k++;
        FORLIM1 = WITH->anzatome;
        for (j = 0; j < FORLIM1; j++) {
          WITH1 = WITH->atom[j];
          WITH1->bindung[k - 1] = WITH1->bindung[i];
        }
      }
    }
    for (i = k; i < maxanzahlbindungen; i++) {
      FORLIM = WITH->anzatome;
      for (j = 0; j < FORLIM; j++) {
        WITH1 = WITH->atom[j];
        WITH1->bindung[i].dessen_atomart = NULL;
      }
    }
    if (   a1->dessen_atomart->deren_atomlage
        != a2->dessen_atomart->deren_atomlage) {
      k = 0;
      WITH = a2->dessen_atomart->deren_atomlage;
      for (i = 0; i < maxanzahlbindungen; i++) {
        if (b2[i]) {
          k++;
          FORLIM = WITH->anzatome;
          for (j = 0; j < FORLIM; j++) {
            WITH1 = WITH->atom[j];
            WITH1->bindung[k - 1] = WITH1->bindung[i];
          }
        }
      }
      for (i = k; i < maxanzahlbindungen; i++) {
        FORLIM = WITH->anzatome;
        for (j = 0; j < FORLIM; j++) {
          WITH1 = WITH->atom[j];
          WITH1->bindung[i].dessen_atomart = NULL;
        }
      }
    }
    printf(" ... old bond ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    write_atom_name(&TEMP, a1);
    printf("  -- ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    write_atom_name(&TEMP, a2);
    printf(" broken \n");
    return;
  }
  printf(" WARNING  there is no bond ");
  TEMP.f = stdout;
  *TEMP.name = '\0';
  write_atom_name(&TEMP, a1);
  printf("  -- ");
  TEMP.f = stdout;
  *TEMP.name = '\0';
  write_atom_name(&TEMP, a2);
  putchar('\n');
}


Static Void berechne_alle_bindungen(kt_strukturdaten *struda,
                                    int OnlyDifferentElements)
{
  long i, j;
  distanzen_tabellen_eintrag *hh;
  long FORLIM;
  kt_atomlage *WITH1;
  long FORLIM1;
  _TEXT TEMP;

  FORLIM = struda->anzatomlagen;
  for (j = 0; j < FORLIM; j++) {
    WITH1 = struda->atomlage[j];
    hh = WITH1->dist_entry;
    FORLIM1 = WITH1->koordinationszahl;
    for (i = 1; i <= FORLIM1; ) {
      if (hh == NULL) {
        printf(" ERROR  not enough distances for atom ");
        TEMP.f = stdout;
        *TEMP.name = '\0';
        write_atom_name(&TEMP, &WITH1->grundatom);
        putchar('\n');
        break;
      } else {
        if (atome_miteinander_verbunden(&WITH1->grundatom, &hh->endatom))
          i++;
        else if (! OnlyDifferentElements || strncmp(WITH1->element,
               hh->endatom.dessen_atomart->deren_atomlage->element,
                                                    sizeof(kt_element)))
        {
          erzeuge_bindung(&WITH1->grundatom, &hh->endatom);
          i++;
        }
        hh = hh->naechste_distanz;
      }
    }
  }
}


Static Void berechne_mitte(kt_atom *a, kt_atom *b, double *v)
{
  long i;

  for (i = 0; i <= 2; i++)
    v[i] = (a->dessen_atomart->pkoord[i] + a->elementarzelle[i] +
            b->dessen_atomart->pkoord[i] + b->elementarzelle[i]) / 2;
}


Static Void baue_o_atome_ein(kt_strukturdaten *strukturdaten,
                             boolean *status)
{
  long i, j, k, kk;
  kt_atom b1, b2, oa;
  vektor v;
  gittervektor g;
  kt_atomlage *neue_atomlage;
  boolean b, bb = false;
  long FORLIM;
  kt_atomlage *WITH1;
  long FORLIM1;

  k = 0;
  FORLIM = strukturdaten->anzatomlagen;
  for (i = 0; i < FORLIM; i++) {
    WITH1 = strukturdaten->atomlage[i];
    if (!strncmp(WITH1->element, "O ", sizeof(kt_element))) {
      kk = WITH1->aname[1] - '0';
      if (WITH1->aname[2] != ' ')
        kk = kk * 10 + WITH1->aname[2] - '0';
      if (kk > k)
        k = kk;
    }
  }
  *status = true;
  FORLIM = strukturdaten->anzatomlagen;
  for (i = 0; i < FORLIM; i++) {
    b1 = strukturdaten->atomlage[i]->grundatom;
    if (!atom_hat_gleiches_element(&b1, "O ")) {
      do {
        bb = true;
        FORLIM1 = anzahl_bindungen_des_atoms(&b1);
        for (j = 1; j <= FORLIM1; j++) {
          berechne_nte_bindung(&b1, j, &b2);
          if (b2.dessen_atomart != NULL) {
            if (!atom_hat_gleiches_element(&b2, "O ")) {
              bb = false;
              berechne_mitte(&b1, &b2, v);
              int_vektor(v, g);
              frac_vektor(v, v);
              k++;
              erzeuge_atomlage("O ", k, v, 2L, &neue_atomlage);
              lade_atomlage(&neue_atomlage, strukturdaten, &b);
              loesche_bindung(&b1, &b2);
              if (b) {
                verschiebe_atome_leicht(neue_atomlage);
                oa.dessen_atomart =
                  strukturdaten->atomlage[strukturdaten->anzatomlagen - 1]->
                  atom[0];
                memcpy(oa.elementarzelle, g, sizeof(gittervektor));
                erzeuge_bindung(&b1, &oa);
                if (!atome_miteinander_verbunden(&b2, &oa))
                  erzeuge_bindung(&b2, &oa);
              }
            }
          }
        }
      } while (!bb);
    }
  }
  FORLIM = strukturdaten->anzatomlagen;
  for (i = 0; i < FORLIM; i++) {
    WITH1 = strukturdaten->atomlage[i];
    if (!strncmp(WITH1->element, "O ", sizeof(kt_element)) &&
        WITH1->koordinationszahl > 2) {
      printf(" ERROR  there are too many bonds for %.6s\n", WITH1->aname);
      *status = false;
    }
  }
}


Static Void write_bindungen(_TEXT *t, kt_atom *atom)
{
  long j, k;
  kt_atom bind[maxanzahlbindungen];
  boolean b;

  fprintf(t->f, "\n Bonds from atom ");
  write_atom_name(t, atom);
  putc('\n', t->f);
  for (k = 1; k <= maxanzahlbindungen; k++) {
    berechne_nte_bindung(atom, k, &bind[k - 1]);
    if (bind[k - 1].dessen_atomart != NULL) {
      fprintf(t->f, " (%ld)   -", k);
      write_atom_name(t, &bind[k - 1]);
      fprintf(t->f, "%10.4f", distanz_zwischen_atome(atom, &bind[k - 1]));
      b = true;
      for (j = 1; j < k; j++) {
        if ((sym_aequival_distanzen(atom, &bind[k - 1], atom, &bind[j - 1]) |
             sym_aequival_distanzen(atom, &bind[k - 1], &bind[j - 1], atom))
             && b) {
          b = false;
          fprintf(t->f, "    symm. equival. to bond (%ld)", j);
        }
      }
      putc('\n', t->f);
    }
  }
}


Static Void write_bindungs_matrix(_TEXT *t, kt_atomlage *atomlage)
{
  long j, k;
  kt_atomart *WITH;

  fprintf(t->f, " Atom %.6s\n", atomlage->aname);
  for (j = 1; j <= atomlage->anzatome; j++) {
    WITH = atomlage->atom[j - 1];
    fprintf(t->f, "%4ld:", j);
    for (k = 0; k < maxanzahlbindungen; k++) {
      if (WITH->bindung[k].dessen_atomart != NULL)
        write_atom_name(t, &WITH->bindung[k]);
    }
    putc('\n', t->f);
  }
  putc('\n', t->f);
}


Static Void verknuepfe_t_atome_direkt(void)
{
  long i, di, j;
  double d, d2;
  kt_atom b1, b2;
  long FORLIM;
  kt_atomlage *WITH;
  long FORLIM1;

  FORLIM = strukturdaten.anzatomlagen;
  for (i = 0; i < FORLIM; i++) {
    WITH = strukturdaten.atomlage[i];
    if (!strncmp(WITH->element, "O ", sizeof(kt_element))) {
      while (WITH->koordinationszahl > 2) {
        printf(
          " WARNING  for atom position %.6s the coordination number %ld is greater then 2\n",
          WITH->aname, WITH->koordinationszahl);
        d = distanz_zwischen_atome(&WITH->grundatom,
                                   &WITH->atom[0]->bindung[0]);
        di = 1;
        FORLIM1 = WITH->koordinationszahl;
        for (j = 2; j <= FORLIM1; j++) {
          d2 = distanz_zwischen_atome(&WITH->grundatom,
                                      &WITH->atom[0]->bindung[j - 1]);
          if (d2 > d) {
            d = d2;
            di = j;
          }
        }
        loesche_bindung(&WITH->grundatom, &WITH->atom[0]->bindung[di - 1]);
      }
    }
  }
  for (i = strukturdaten.anzatomlagen - 1; i >= 0; i--) {
    WITH = strukturdaten.atomlage[i];
    if (!strncmp(WITH->element, "O ", sizeof(kt_element))) {
      if (anzahl_bindungen_des_atoms(&WITH->grundatom) > 0) {
        berechne_nte_bindung(&WITH->grundatom, 1L, &b1);
        if (anzahl_bindungen_des_atoms(&WITH->grundatom) == 2) {
          berechne_nte_bindung(&WITH->grundatom, 2L, &b2);
          loesche_bindung(&WITH->grundatom, &b2);
          if (atome_miteinander_verbunden(&WITH->grundatom, &b1))
            loesche_bindung(&WITH->grundatom, &b1);
          erzeuge_bindung(&b1, &b2);
        }
        if (atome_miteinander_verbunden(&WITH->grundatom, &b1))
          loesche_bindung(&WITH->grundatom, &b1);
        delete_atomlage(strukturdaten.atomlage[i], &strukturdaten);
      }
    }
  }
}


Static Void berechne_koord_seq(kt_atom *atoco, long anz, long *koord_seq)
{
  long i, j, k, l, n;
  kt_atom sphae[2000];
  long z[12];
  long FORLIM, FORLIM1, FORLIM2;

  z[0] = 0;
  z[1] = 1;
  z[2] = koordinationszahl_des_atoms(atoco) + 1;
  sphae[0] = *atoco;
  FORLIM = koordinationszahl_des_atoms(atoco);
  for (i = 1; i <= FORLIM; i++)
    berechne_nte_bindung(atoco, i, &sphae[i]);
  n = z[2] + 1;
  for (i = 3; i <= anz; i++) {
    FORLIM1 = z[i - 1];
    for (j = z[i - 2]; j < FORLIM1; j++) {
      FORLIM2 = koordinationszahl_des_atoms(&sphae[j]);
      for (k = 1; k <= FORLIM2; k++) {
        berechne_nte_bindung(&sphae[j], k, &sphae[n - 1]);
        l = z[i - 3];
        do {
          l++;
        } while (!gleiche_atome(&sphae[l - 1], &sphae[n - 1]));
        if (l == n)
          n++;
      }
    }
    z[i] = n - 1;
  }
  for (i = 1; i <= 10; i++)
    koord_seq[i - 1] = z[i + 1] - z[i];
}


Static Void write_koord_seq(_TEXT *f, long *z)
{
  long i;

  for (i = 0; i <= 9; i++)
    fprintf(f->f, "%4ld", z[i]);
  putc('\n', f->f);
}


Static Void bestimme_topologie(_TEXT *Tfpout)
{
  long i, j, k;
  tt_koord_sequenz z[maxanzahlatomlagen];
  tt_koord_sequenz ze[200];
  long anzze;
  Char zecod[6], zecodneu[6];
  static Char zcod[maxanzahlatomlagen][200][6]; /* 50 = max number of       */
  boolean bz[100], bze[100];                   /*      possible topologies */
  long nzcod;
  boolean b;
  long sdaal, TD10, SumN, SumM;
  kt_atomlage *al;

  verknuepfe_t_atome_direkt();
  fprintf(Tfpout->f, "\n Coordination sequences:\n");
  sdaal = strukturdaten.anzatomlagen;
  for (i = 0; i < sdaal; i++) {
    al = strukturdaten.atomlage[i];
    if (al->koordinationszahl > 0) {
      berechne_koord_seq(&al->grundatom, 11L, z[i]);
      fprintf(Tfpout->f, " %.6s", al->aname);
      write_koord_seq(Tfpout, z[i]);
    }
  }

  fprintf(Tfpout->f, " Mean  ");

  for (j = 0; j <= 9; j++)
  {
    SumM = 0;
    SumN = 0;

    for (i = 0; i < sdaal; i++)
    {
          al = strukturdaten.atomlage[i];
      if (al->koordinationszahl > 0)
      {
        SumM += al->anzatome;
        SumN += al->anzatome * z[i][j];
      }
    }

    fprintf(Tfpout->f, " %.2lf", (double) SumN / SumM);
  }

  putc('\n', Tfpout->f);

  TD10 = 0;
     k = 0;

  for (i = 0; i < sdaal; i++)
  {
        al = strukturdaten.atomlage[i];
    if (al->koordinationszahl > 0)
    {
      SumN = 1;

      for (j = 0; j <= 9; j++)
        SumN += z[i][j];

      TD10 += al->anzatome * SumN;
         k += al->anzatome;
    }
  }

  print_tlnb(Tfpout->f, strukturdaten.abk, sizeof (abkbezeichnung),
    " TD10 ", NULL);
  fprintf(Tfpout->f, " %.6g\n", (double) TD10 / k);

  fprintf(Tfpout->f, " Possible frameworks: \n");
  for (i = 0; i < sdaal; i++) {
    al = strukturdaten.atomlage[i];
    if (al->koordinationszahl > 0) {
      for (j = 0; j <= 9; j++) {
        for (k = 0; k <= 5; k++)
          zcod[i][j][k] = ' ';
      }
    }
  }
  nzcod = 0;
  if (*coseq.name != '\0') {
    if (coseq.f != NULL)
      coseq.f = freopen(coseq.name, "r", coseq.f);
    else
      coseq.f = fopen(coseq.name, "r");
  } else
    rewind(coseq.f);
  if (coseq.f == NULL)
    _EscIO(FileNotFound);
  RESETBUF(coseq.f, Char);
  for (i = 0; i <= 5; i++) {
    zecodneu[i] = getc(coseq.f);
    if (zecodneu[i] == '\n')
      zecodneu[i] = ' ';
  }

  if (zecodneu[5] != ' ' && zecodneu[5] != '\t')
  {
    do
      i = getc(coseq.f);
    while (i != EOF && i != ' ' && i != '\t' && i != '\0');
  }

  while (!BUFEOF(coseq.f)) {
    anzze = 0;
    do {
      memcpy(zecod, zecodneu, 6L);
      anzze++;
      for (i = 0; i <= 9; i++)
        fscanf(coseq.f, "%ld", &ze[anzze - 1][i]);
      fscanf(coseq.f, "%*[^\n]");
      getc(coseq.f);
      if (!BUFEOF(coseq.f)) {
        for (i = 0; i <= 5; i++) {
          zecodneu[i] = getc(coseq.f);
          if (zecodneu[i] == '\n')
            zecodneu[i] = ' ';
        }

        if (zecodneu[5] != ' ' && zecodneu[5] != '\t')
        {
          do
            i = getc(coseq.f);
          while (i != EOF && i != ' ' && i != '\t' && i != '\0');
        }

      } else
        memcpy(zecodneu, "******", 6L);
    } while (!strncmp(zecod, zecodneu, 6));
    for (i = 0; i < sdaal; i++)
      bz[i] = false;
    for (j = 0; j < anzze; j++)
      bze[j] = false;
    for (i = 0; i < sdaal; i++) {
      al = strukturdaten.atomlage[i];
      if (al->koordinationszahl > 0) {
        for (j = 0; j < anzze; j++) {
          b = true;
          for (k = 0; k <= 9; k++) {
            if (ze[j][k] != z[i][k])
              b = false;
          }
          if (b) {
            bz[i] = true;
            bze[j] = true;
          }
        }
      }
    }
    b = true;
    for (j = 0; j < anzze; j++) {
      if (!bze[j])
        b = false;
    }
    if (!b)
      continue;
    nzcod++;
    for (i = 0; i < sdaal; i++) {
      al = strukturdaten.atomlage[i];
      if (al->koordinationszahl > 0) {
        if (bz[i])
          memcpy(zcod[i][nzcod - 1], zecod, 6L);
      }
    }
  }
  for (i = 0; i < sdaal; i++) {
    al = strukturdaten.atomlage[i];
    if (al->koordinationszahl > 0) {
      fprintf(Tfpout->f, " %.6s", al->aname);
      for (j = 0; j < nzcod; j++)
        fprintf(Tfpout->f, "%.6s", zcod[i][j]);
      putc('\n', Tfpout->f);
    }
  }
}


typedef struct wegtyp {
  kt_atom atom[20];
  long anzatome;
} wegtyp;


/* Local variables for berechne_loop_konfig: */
struct LOC_berechne_loop_konfig {
  long anzloesungen;
  wegtyp loesung[50];
} ;

Local Void suche_weg(kt_atom *ort, kt_atom *ziel,
                     long schritte, wegtyp *weg,
                     struct LOC_berechne_loop_konfig *LINK)
{
  long i, j;
  kt_atom neuer_ort;
  boolean b;
  long FORLIM;

  if (schritte > 0) {
    FORLIM = anzahl_bindungen_des_atoms(ort);
    for (i = 1; i <= FORLIM; i++) {
      berechne_nte_bindung(ort, i, &neuer_ort);
      weg->atom[schritte - 1] = neuer_ort;
      b = true;
      for (j = weg->anzatome - 1; j >= schritte; j--) {
        if (gleiche_atome(&weg->atom[j], &neuer_ort))
          b = false;
      }
      if (b)
        suche_weg(&neuer_ort, ziel, schritte - 1, weg, LINK);
    }
    return;
  }
  if (gleiche_atome(ort, ziel)) {
    LINK->anzloesungen++;
    LINK->loesung[LINK->anzloesungen - 1] = *weg;
  }
}


Static Void berechne_loop_konfig(_TEXT *Tfpout, kt_atom *zentralatom)
{
  struct LOC_berechne_loop_konfig V;
  long i, j, k, l;
  wegtyp weg;
  kt_atom nachbaratom[maxanzahlbindungen];
  boolean b;
  long FORLIM;
  long FORLIM2;
  wegtyp *WITH;
  long FORLIM3;

  b = true;
  FORLIM = anzahl_bindungen_des_atoms(zentralatom);
  for (i = 0; i < FORLIM; i++) {
    berechne_nte_bindung(zentralatom, i + 1, &nachbaratom[i]);
    if (gleiche_atome(zentralatom, &nachbaratom[i])) {
      fprintf(Tfpout->f, " ERROR  bond ");
      write_atom_name(Tfpout, zentralatom);
      fprintf(Tfpout->f, "  -- ");
      write_atom_name(Tfpout, &nachbaratom[i]);
      fprintf(Tfpout->f, " not allowed\n");
      b = false;
    }
  }
  if (!b)
    return;
  FORLIM = anzahl_bindungen_des_atoms(zentralatom);
  for (i = 1; i <= FORLIM; i++) {
    for (j = 0; j <= i - 2; j++) {
      weg.anzatome = 1;
      do {
        weg.anzatome++;
        V.anzloesungen = 0;
        weg.atom[weg.anzatome - 1] = *zentralatom;
        weg.atom[weg.anzatome - 2] = nachbaratom[i - 1];
        suche_weg(&nachbaratom[i - 1], &nachbaratom[j],
                  weg.anzatome - 2, &weg, &V);
      } while (V.anzloesungen <= 0);
      fprintf(Tfpout->f, "\n Angle ");
      write_atomart_name(Tfpout, nachbaratom[i - 1].dessen_atomart);
      fprintf(Tfpout->f, "  -- ");
      write_atomart_name(Tfpout, zentralatom->dessen_atomart);
      fprintf(Tfpout->f, "  -- ");
      write_atomart_name(Tfpout, nachbaratom[j].dessen_atomart);
      if (V.anzloesungen == 1)
        fprintf(Tfpout->f, "  has 1 loop  with %ld atoms\n",
          V.loesung[0].anzatome);
      else
        fprintf(Tfpout->f, "  has %ld loops with %ld atoms\n",
               V.anzloesungen, V.loesung[0].anzatome);
      FORLIM2 = V.anzloesungen;
      for (k = 0; k < FORLIM2; k++) {
        WITH = &V.loesung[k];
        fprintf(Tfpout->f, "  -");
        write_atom_name(Tfpout, &WITH->atom[WITH->anzatome - 1]);
        putc(' ', Tfpout->f);
        FORLIM3 = WITH->anzatome;
        for (l = 1; l < FORLIM3; l++) {
          if ((l & 3) == 0 && l > 1)
            fprintf(Tfpout->f, "   ");
          write_atom_name(Tfpout, &WITH->atom[l - 1]);
          putc(' ', Tfpout->f);
          if ((l & 3) == 3)
            putc('\n', Tfpout->f);
        }
        if ((WITH->anzatome & 3) != 0)
          putc('\n', Tfpout->f);
      }
    }
  }
}


Static Void loesche_alle_bindungen(kt_strukturdaten *struda)
{
  long i;
  kt_atom a;
  long FORLIM;
  kt_atomlage *WITH1;

  FORLIM = struda->anzatomlagen;
  for (i = 0; i < FORLIM; i++) {
    WITH1 = struda->atomlage[i];
    while (anzahl_bindungen_des_atoms(&WITH1->grundatom) > 0) {
      berechne_nte_bindung(&WITH1->grundatom, 1L, &a);
      loesche_bindung(&WITH1->grundatom, &a);
    }
  }
}


typedef Char nametyp[6];

typedef struct _REC_symeq {
  kt_atom atona;
  nametyp erstername, zweitername;
} _REC_symeq;

typedef struct dls_symeq_karten {
  _REC_symeq symeq[maxanzahlsymeqkarten][20];
  long anzsymeq[maxanzahlsymeqkarten];
} dls_symeq_karten;

typedef struct _REC_tetcon {
  nametyp n[9];
  kt_atom a[9];
} _REC_tetcon;

typedef struct dls_tetcon_karten {
  _REC_tetcon tetcon[maxanzahlatomlagen];
  long anztetcon;
} dls_tetcon_karten;

typedef struct _REC_distcon {
  kt_atom atom1, atom2;
  nametyp atomname1, atomname2;
} _REC_distcon;

typedef struct dls_distan_karten {
  _REC_distcon distcon[maxindepdistanzen];
  long anzdistcon;
} dls_distan_karten;

typedef struct _REC_anglcon {
  kt_atom atom1, atom2, atomm;
  nametyp atomname1, atomnamem, atomname2;
} _REC_anglcon;

typedef struct dls_angle_karten {
  _REC_anglcon anglcon[maxindepwinkel];
  long anzanglcon;
} dls_angle_karten;

typedef struct _REC_dt_geg_distanz {
  kt_element el1, el2;
  long kz1, kz2;
  double distanz, dvar;
} _REC_dt_geg_distanz;

typedef _REC_dt_geg_distanz dt_geg_distanz[maxanzahlgegdistanzen];

typedef struct _REC_dt_geg_winkel {
  kt_element el1, el2, el3;
  long kz1, kz2, kz3;
  double winkel, wvar;
} _REC_dt_geg_winkel;

typedef _REC_dt_geg_winkel dt_geg_winkel[maxanzahlgegwinkel];


typedef struct _REC_bondis_type {
  kt_element e1, e2, e3;
  kt_atom a1, a2, a3;
} _REC_bondis_type;

typedef _REC_bondis_type bondis_type[20];


/* Local variables for write_input_file: */
struct LOC_write_input_file {
  kt_strukturdaten strukturdaten;
} ;

Local Void read_bindungs_groessen(_TEXT *f,
                                  _REC_dt_geg_distanz *dist,
                                  _REC_dt_geg_winkel *wink)
{
  long i, anz_d, anz_w;
  kt_element e1, e2, e3;
  long k1, k2, k3;
  _REC_dt_geg_distanz *WITH;
  _REC_dt_geg_winkel *WITH1;

  anz_d = 1;
  anz_w = 1;
  if (*f->name != '\0') {
    if (f->f != NULL)
      f->f = freopen(f->name, "r", f->f);
    else
      f->f = fopen(f->name, "r");
  } else
    rewind(f->f);
  if (f->f == NULL)
    _EscIO(FileNotFound);
  RESETBUF(f->f, Char);
  while (!BUFEOF(f->f)) {
    e3[0] = ' ';
    e1[1] = ' ';
    e2[1] = ' ';
    e3[1] = ' ';
    while (P_peek(f->f) == ']' || P_peek(f->f) == '[' || P_peek(f->f) == ' ')
      getc(f->f);
    *e1 = getc(f->f);
    if (e1[0] == '\n')
      e1[0] = ' ';
    if (P_peek(f->f) == ' ' || isalpha(P_peek(f->f))) {
      e1[1] = getc(f->f);
      if (e1[1] == '\n')
        e1[1] = ' ';
    }
    while (P_peek(f->f) == ']' || P_peek(f->f) == '[' || P_peek(f->f) == ' ')
      getc(f->f);
    fscanf(f->f, "%ld", &k1);
    while (P_peek(f->f) == ']' || P_peek(f->f) == '[' || P_peek(f->f) == ' ')
      getc(f->f);
    *e2 = getc(f->f);
    if (e2[0] == '\n')
      e2[0] = ' ';
    if (P_peek(f->f) == ' ' || isalpha(P_peek(f->f))) {
      e2[1] = getc(f->f);
      if (e2[1] == '\n')
        e2[1] = ' ';
    }
    while (P_peek(f->f) == ']' || P_peek(f->f) == '[' || P_peek(f->f) == ' ')
      getc(f->f);
    fscanf(f->f, "%ld", &k2);
    while (P_peek(f->f) == ']' || P_peek(f->f) == '[' || P_peek(f->f) == ' ')
      getc(f->f);
    if (P_peek(f->f) == ' ' || isalpha(P_peek(f->f))) {
      *e3 = getc(f->f);
      if (e3[0] == '\n')
        e3[0] = ' ';
      if (P_peek(f->f) == ' ' || isalpha(P_peek(f->f))) {
        e3[1] = getc(f->f);
        if (e3[1] == '\n')
          e3[1] = ' ';
      }
      while (P_peek(f->f) == ']' || P_peek(f->f) == '[' || P_peek(f->f) == ' ')
        getc(f->f);
      fscanf(f->f, "%ld", &k3);
      while (P_peek(f->f) == ']' || P_peek(f->f) == '[' || P_peek(f->f) == ' ')
        getc(f->f);
    }
    if (e3[0] == ' ') {
      WITH = &dist[anz_d - 1];
      memcpy(WITH->el1, e1, sizeof(kt_element));
      WITH->kz1 = k1;
      memcpy(WITH->el2, e2, sizeof(kt_element));
      WITH->kz2 = k2;
      fscanf(f->f, "%lg%lg%*[^\n]", &WITH->distanz, &WITH->dvar);
      getc(f->f);
      anz_d++;
      continue;
    }
    WITH1 = &wink[anz_w - 1];
    memcpy(WITH1->el1, e1, sizeof(kt_element));
    WITH1->kz1 = k1;
    memcpy(WITH1->el2, e2, sizeof(kt_element));
    WITH1->kz2 = k2;
    memcpy(WITH1->el3, e3, sizeof(kt_element));
    WITH1->kz3 = k3;
    fscanf(f->f, "%lg%lg%*[^\n]", &WITH1->winkel, &WITH1->wvar);
    getc(f->f);
    anz_w++;
  }
  for (i = anz_d - 1; i < maxanzahlgegdistanzen; i++) {
    WITH = &dist[i];
    memcpy(WITH->el1, "  ", sizeof(kt_element));
    memcpy(WITH->el2, "  ", sizeof(kt_element));
    WITH->kz1 = 0;
    WITH->kz2 = 0;
    WITH->distanz = 0.0;
    WITH->dvar = 0.0;
  }
  for (i = anz_w - 1; i < maxanzahlgegwinkel; i++) {
    WITH1 = &wink[i];
    memcpy(WITH1->el1, "  ", sizeof(kt_element));
    memcpy(WITH1->el2, "  ", sizeof(kt_element));
    memcpy(WITH1->el3, "  ", sizeof(kt_element));
    WITH1->kz1 = 0;
    WITH1->kz2 = 0;
    WITH1->kz3 = 0;
    WITH1->winkel = 0.0;
    WITH1->wvar = 0.0;
  }
}

Local Void write_dls_atom_karten(_TEXT *f, kt_strukturdaten *s, int *Polar)
{
  long i;
  vektor v;
  kt_atomlage *WITH;

  int         mode;
  char        *NOREF_aname;
  static char PolarLetter[] = { 'X', 'Y', 'Z' };


  NOREF_aname = NULL;

  for (i = 0; i < s->anzatomlagen; i++) {
    WITH = s->atomlage[i];
    berechne_atom_koordinaten(&WITH->grundatom, v);
    fprintf(f->f, "ATOM   %.6s%8.5f%8.5f%8.5f  %.2s   ",
            WITH->aname, v[0], v[1], v[2], WITH->element);
    write_atom_trilinpol(f, &WITH->grundatom);
    putc('\n', f->f);
    if (NOREF_aname == NULL) NOREF_aname = WITH->aname;
  }

  if (NOREF_aname)
  {
    mode = 0;

    for (i = 0; i < 3; i++)
    {
      if (Polar[i])
      {
        if (mode == 0)
        {
          fprintf(f->f, "NOREF  %.6s", NOREF_aname);
          mode = 1;
        }
        else
          fprintf(f->f, "  ");

        putc(PolarLetter[i], f->f);
      }
    }

    if (mode)
      putc('\n', f->f);
  }
}

Local Void write_loadat_atom_karten(_TEXT *f, kt_strukturdaten *s)
{
  long i;
  kt_atomlage *WITH;

  for (i = 0; i < s->anzatomlagen; i++) {
    WITH = s->atomlage[i];
    fprintf(f->f, "ATOM   %.6s%8.5f%8.5f%8.5f%6.3f\n",
            WITH->aname, WITH->atomlagekoord[0], WITH->atomlagekoord[1],
            WITH->atomlagekoord[2], 0.025);
  }
}

Local Void write_dls_symeq_karten(_TEXT *f,
                                  dls_symeq_karten *dls_symeq,
                                  struct LOC_write_input_file *LINK)
{
  long i, j, anzsymop;
  st_symop_array symop;
  long FORLIM, FORLIM1;
  _REC_symeq *WITH;

  FORLIM = LINK->strukturdaten.anzatomlagen;
  for (i = 0; i < FORLIM; i++) {
    FORLIM1 = dls_symeq->anzsymeq[i];
    for (j = 0; j < FORLIM1; j++) {
      WITH = &dls_symeq->symeq[i][j];
      fprintf(f->f, "SYMEQ  %.6s %.6s ", WITH->erstername, WITH->zweitername);
      welchsymop(&WITH->atona.dessen_atomart->deren_atomlage->grundatom,
                 &WITH->atona, &rginfo, &anzsymop, symop);
      write_symmetrieoperator(f, &symop[0]);
      putc('\n', f->f);
    }
  }
}

Local Void write_dls_tetcon_karten(_TEXT *f, dls_tetcon_karten *dls_tetcon)
{
  long i, j, FORLIM;
  _REC_tetcon *WITH;

  FORLIM = dls_tetcon->anztetcon;
  for (i = 0; i < FORLIM; i++) {
    WITH = &dls_tetcon->tetcon[i];
    fprintf(f->f, "TETCON       ");
    for (j = 0; j <= 8; j++)
      fprintf(f->f, "%.6s ", WITH->n[j]);
    putc('\n', f->f);
  }
}

Local Void bestimme_geg_distanzen(kt_atom *atom1, kt_atom *atom2,
                                  _REC_dt_geg_distanz *geg_dist_,
                                  double *d, double *v, boolean *b)
{
  dt_geg_distanz geg_dist;
  long l, kz_a1, kz_a2;
  kt_element el_a1, el_a2;
  _REC_dt_geg_distanz *WITH;

  memcpy(geg_dist, geg_dist_, sizeof(dt_geg_distanz));
  *b = false;
  *d = 0.0;
  *v = 0.0;
  memcpy(el_a1, atom1->dessen_atomart->deren_atomlage->element,
         sizeof(kt_element));
  memcpy(el_a2, atom2->dessen_atomart->deren_atomlage->element,
         sizeof(kt_element));
  kz_a1 = koordinationszahl_des_atoms(atom1);
  kz_a2 = koordinationszahl_des_atoms(atom2);
  l = 0;
  do {
    l++;
    WITH = &geg_dist[l - 1];
    if ((!strncmp(el_a1, WITH->el1, sizeof(kt_element)) &&
         kz_a1 == WITH->kz1 &&
         !strncmp(el_a2, WITH->el2, sizeof(kt_element)) &&
         kz_a2 == WITH->kz2) ||
        (!strncmp(el_a1, WITH->el2, sizeof(kt_element)) &&
         kz_a1 == WITH->kz2 &&
         !strncmp(el_a2, WITH->el1, sizeof(kt_element)) &&
         kz_a2 == WITH->kz1)) {
      *d = WITH->distanz;
      *v = WITH->dvar;
      *b = true;
    }
  } while (strncmp(geg_dist[l - 1].el1, "  ", sizeof(kt_element)) && !*b);
  if (!*b)
    printf(" WARNING  no information for distance %.2s[%ld] - %.2s[%ld]\n",
           el_a1, kz_a1, el_a2, kz_a2);
}

Local Void bestimme_geg_winkel(kt_atom *atom1, kt_atom *atom2, kt_atom *atom3,
                               _REC_dt_geg_winkel *geg_wink_,
                               double *w, double *wv, boolean *b)
{
  dt_geg_winkel geg_wink;
  long l, kz_a1, kz_a2, kz_a3;
  kt_element el_a1, el_a2, el_a3;
  _REC_dt_geg_winkel *WITH;

  memcpy(geg_wink, geg_wink_, sizeof(dt_geg_winkel));
  *w = 0.0;
  *wv = 0.0;
  *b = false;
  memcpy(el_a1, atom1->dessen_atomart->deren_atomlage->element,
         sizeof(kt_element));
  memcpy(el_a2, atom2->dessen_atomart->deren_atomlage->element,
         sizeof(kt_element));
  memcpy(el_a3, atom3->dessen_atomart->deren_atomlage->element,
         sizeof(kt_element));
  kz_a1 = koordinationszahl_des_atoms(atom1);
  kz_a2 = koordinationszahl_des_atoms(atom2);
  kz_a3 = koordinationszahl_des_atoms(atom3);
  l = 0;
  do {
    l++;
    WITH = &geg_wink[l - 1];
    if ((!strncmp(el_a1, WITH->el1, sizeof(kt_element)) &&
         kz_a1 == WITH->kz1 &&
         !strncmp(el_a2, WITH->el2, sizeof(kt_element)) &&
         kz_a2 == WITH->kz2 &&
         !strncmp(el_a3, WITH->el3, sizeof(kt_element)) &&
         kz_a3 == WITH->kz3) ||
        (!strncmp(el_a1, WITH->el3, sizeof(kt_element)) &&
         kz_a1 == WITH->kz3 &&
         !strncmp(el_a2, WITH->el2, sizeof(kt_element)) &&
         kz_a2 == WITH->kz2 &&
         !strncmp(el_a3, WITH->el1, sizeof(kt_element)) &&
         kz_a3 == WITH->kz1)) {
      *w = WITH->winkel;
      *wv = WITH->wvar;
      *b = true;
    }
  } while (strncmp(geg_wink[l - 1].el1, "  ", sizeof(kt_element)) && !*b);
  if (!*b)
    printf(" WARNING  no information for angle %.2s[%ld] - %.2s[%ld] - %.2s[%ld]\n",
           el_a1, kz_a1, el_a2, kz_a2, el_a3, kz_a3);
}

Local Void write_dls_distan_karten(_TEXT *f,
                                   _REC_dt_geg_distanz *geg_dist,
                                   dls_distan_karten *dls_distan)
{
  long i;
  boolean b;
  double distanz, dvar;
  long FORLIM;
  _REC_distcon *WITH;

  FORLIM = dls_distan->anzdistcon;
  for (i = 0; i < FORLIM; i++) {
    WITH = &dls_distan->distcon[i];
    bestimme_geg_distanzen(&WITH->atom1, &WITH->atom2, geg_dist, &distanz,
                           &dvar, &b);
    fprintf(f->f, "DISTAN %.6s %.6s ", WITH->atomname1, WITH->atomname2);
    if (b)
      fprintf(f->f, "%10.6f%9.4f\n", distanz, 1 / dvar / 50);
    else
      fprintf(f->f, "%10.6f%9.4f\n", 1.0, 0.0);
  }
}

Local Void write_loadat_distan_karten(_TEXT *f,
                                      _REC_dt_geg_distanz *geg_dist,
                                      dls_distan_karten *dls_distan)
{
  long i;
  boolean b;
  double distanz, dvar;
  long FORLIM;
  _REC_distcon *WITH;

  FORLIM = dls_distan->anzdistcon;
  for (i = 0; i < FORLIM; i++) {
    WITH = &dls_distan->distcon[i];
    bestimme_geg_distanzen(&WITH->atom1, &WITH->atom2, geg_dist, &distanz,
                           &dvar, &b);
    fprintf(f->f, "DISTAN %.6s %.6s ", WITH->atomname1, WITH->atomname2);
    if (b)
      fprintf(f->f, "%10.6f%9.4f\n", distanz, dvar);
    else
      fprintf(f->f, "%10.6f%9.4f\n", 1.0, 100.0);
  }
}

Local Void berechne_t_t_distanzen(double d1,  double d2,  double w,
                                  double d1v, double d2v, double wv,
                                  double *d, double *dv)
{
  w /= 57.29577951;
  wv /= 57.29577951;
  *d = sqrt(d1 * d1 + d2 * d2 - 2 * d1 * d2 * cos(w));
  *dv = (d1 * d1v + d2 * d2v - (d1 * d2v + d2 * d1v) * cos(w) +
         d1 * d2 * sin(w) * wv) / *d;
}

Local Void write_dls_angle_karten(_TEXT *f,
                                  _REC_dt_geg_winkel *geg_wink,
                                  _REC_dt_geg_distanz *geg_dist,
                                  dls_angle_karten *dls_angle)
{
  boolean b, bb;
  double d1, d2, d1v, d2v, w, wv, distanz, dvar;
  long i, FORLIM;
  _REC_anglcon *WITH;

  FORLIM = dls_angle->anzanglcon;
  for (i = 0; i < FORLIM; i++) {
    WITH = &dls_angle->anglcon[i];
    bestimme_geg_distanzen(&WITH->atom1, &WITH->atomm,
                           geg_dist, &d1, &d1v, &bb);
    bestimme_geg_distanzen(&WITH->atom2, &WITH->atomm,
                           geg_dist, &d2, &d2v, &b);
    if (!b)
      bb = false;
    bestimme_geg_winkel(&WITH->atom1, &WITH->atomm, &WITH->atom2, geg_wink, &w,
                        &wv, &b);
    if (!b)
      bb = false;
    if (bb) {
      berechne_t_t_distanzen(d1, d2, w, d1v, d2v, wv, &distanz, &dvar);
      fprintf(f->f, "DISTAN %.6s %.6s %10.6f%9.4f\n",
              WITH->atomname1, WITH->atomname2, distanz, 1 / dvar / 50);
    } else
      fprintf(f->f, "DISTAN %.6s %.6s %10.6f%9.4f\n",
              WITH->atomname1, WITH->atomname2, 1.0, 0.0);
  }
}

Local Void write_loadat_angle_karten(_TEXT *f,
                                     _REC_dt_geg_winkel *geg_wink,
                                     dls_angle_karten *dls_angle)
{
  boolean b;
  double winkel, wvar;
  long i, FORLIM;
  _REC_anglcon *WITH;

  FORLIM = dls_angle->anzanglcon;
  for (i = 0; i < FORLIM; i++) {
    WITH = &dls_angle->anglcon[i];
    bestimme_geg_winkel(&WITH->atom1, &WITH->atomm, &WITH->atom2, geg_wink,
                        &winkel, &wvar, &b);
    fprintf(f->f, "ANGLE  %.6s %.6s %.6s ",
            WITH->atomname1, WITH->atomnamem, WITH->atomname2);
    if (b)
      fprintf(f->f, "%10.6f%9.4f\n", winkel, wvar);
    else
      fprintf(f->f, "%10.6f%9.4f\n", 90.0, 0.0);
  }
}

Local Void write_dls_bondis_karten(_TEXT *t,
                                   dls_tetcon_karten *dls_tetcon,
                                   _REC_dt_geg_distanz *geg_dist,
                                   _REC_dt_geg_winkel *geg_wink)
{
  bondis_type bondis;
  long i, j, l, nn, kzz;
  double a, b, c, d, e, f, g, d1, d2, d1v, d2v, w1, w2, w1v, w2v, d11, d12,
         d11v, d12v;
  boolean d1b, d2b, w1b, w2b, bb;
  long FORLIM;
  _REC_tetcon *WITH;
  _REC_bondis_type *WITH1;

  a = 1.0000;
  b = -0.0004;
  c = 0.0000;
  d = 145.00;
  e = 1.0000;
  f = 1.0000;
  g = 1.0000;
  nn = 0;
  FORLIM = dls_tetcon->anztetcon;
  for (i = 0; i < FORLIM; i++) {
    WITH = &dls_tetcon->tetcon[i];
    for (j = 1; j <= 4; j++) {
      WITH1 = &bondis[nn];
      WITH1->a1 = WITH->a[0];
      memcpy(WITH1->e1, WITH1->a1.dessen_atomart->deren_atomlage->element,
             sizeof(kt_element));
      WITH1->a2 = WITH->a[j];
      memcpy(WITH1->e2, WITH1->a2.dessen_atomart->deren_atomlage->element,
             sizeof(kt_element));
      WITH1->a3 = WITH->a[j + 4];
      memcpy(WITH1->e3, WITH1->a3.dessen_atomart->deren_atomlage->element,
             sizeof(kt_element));
      bb = true;
      for (l = 0; l < nn; l++) {
        if (!strncmp(WITH1->e1, bondis[l].e1, sizeof(kt_element)) &&
            !strncmp(WITH1->e2, bondis[l].e2, sizeof(kt_element)) &&
            !strncmp(WITH1->e3, bondis[l].e3, sizeof(kt_element)))
          bb = false;
      }
      if (bb)
        nn++;
    }
  }
  for (i = 0; i < nn; i++) {
    WITH1 = &bondis[i];
    kzz = WITH1->a3.dessen_atomart->deren_atomlage->koordinationszahl;
    WITH1->a3.dessen_atomart->deren_atomlage->koordinationszahl = 4;
    bestimme_geg_distanzen(&WITH1->a1, &WITH1->a2, geg_dist, &d1, &d1v, &d1b);
    bestimme_geg_distanzen(&WITH1->a3, &WITH1->a2, geg_dist, &d2, &d2v, &d2b);
    bestimme_geg_winkel(&WITH1->a1, &WITH1->a2, &WITH1->a3, geg_wink, &w1, &w1v,
                        &w1b);
    bestimme_geg_winkel(&WITH1->a2, &WITH1->a1, &WITH1->a2, geg_wink, &w2, &w2v,
                        &w2b);
    WITH1->a3.dessen_atomart->deren_atomlage->koordinationszahl = kzz;
    if (d1b) {
      a = d1;
      e = 1 / d1v;
    }
    if (w1b)
      d = w1;
    if (d1b && w2b) {
      berechne_t_t_distanzen(d1, d1, w2, d1v, d1v, w2v, &d11, &d11v);
      f = 1 / d11v;
    }
    if (d1b && d2b && w1b) {
      berechne_t_t_distanzen(d1, d2, w1, d1v, d2v, w1v, &d12, &d12v);
      g = 1 / d12v;
    }
    fprintf(t->f,
      "BONDIS %.2s  %.2s  %.2s   %10.5f%10.5f%10.5f%10.5f%5.1f%5.1f%5.1f\n",
      WITH1->e1, WITH1->e2, WITH1->e3, a, b, c, d, e / 50, f / 50, g / 50);
  }
}

Local Void korrigiere_koordinationszahlen(boolean *b,
                                          struct LOC_write_input_file *LINK)
{
  long i, kzn;
  Char c;
  long FORLIM;
  kt_atomlage *WITH;

  *b = true;
  FORLIM = LINK->strukturdaten.anzatomlagen;
  for (i = 0; i < FORLIM; i++) {
    WITH = LINK->strukturdaten.atomlage[i];
    kzn = anzahl_bindungen_des_atoms(&WITH->grundatom);
    if (kzn != WITH->koordinationszahl) {
      if (kzn == 1)
        printf(" ERROR  for atom %.6s there is 1 bond \n", WITH->aname);
      else
        printf(" ERROR  for atom %.6s there are %ld bonds \n",
               WITH->aname, kzn);
      printf("        but the given coordination number is %ld\n",
             WITH->koordinationszahl);
      *b = false;
    }
  }
  if (*b)
    return;
  printf(" change coordination numbers to number of bonds? (y/n):");
  if (P_eoln(stdin))
    return;
  scanf("%c%*[^\n]", &c);
  getchar();
  if (c == '\n')
    c = ' ';
  if (c != 'j' && c != 'J' && c != 'y' && c != 'Y')
    return;
  *b = true;
  FORLIM = LINK->strukturdaten.anzatomlagen;
  for (i = 0; i < FORLIM; i++) {
    WITH = LINK->strukturdaten.atomlage[i];
    kzn = anzahl_bindungen_des_atoms(&WITH->grundatom);
    if (kzn != WITH->koordinationszahl) {
      printf(" ... coordination number for atom position %.6s", WITH->aname);
      printf(" changed to %ld\n", kzn);
      WITH->koordinationszahl = kzn;
    }
  }
}


Local Void berechne_dls_input(dls_tetcon_karten *dls_tetcon,
                              dls_distan_karten *dls_distan,
                              dls_angle_karten *dls_angle,
                              dls_symeq_karten *dls_symeq,
                              boolean *b,
                              struct LOC_write_input_file *LINK)
{
  long i, j, k, l, anzsymop;
  st_symop_array symeq0;
  kt_atom a1, a2, am, atom0;
  _REC_distcon *WITH1;
  long FORLIM, FORLIM1;
  kt_atomart *WITH2;
  _REC_tetcon *WITH3;
  long FORLIM2, FORLIM3;
  _REC_anglcon *WITH4;

  korrigiere_koordinationszahlen(b, LINK);
  if (!*b) {
    printf(" ... generation of input file interrupted\n");
    return;
  }
  for (i = 0; i < maxanzahlsymeqkarten; i++)
    dls_symeq->anzsymeq[i] = 0;
  for (i = 0; i <= 3; i++) {
    for (j = 0; j <= 8; j++)
      memcpy(dls_tetcon->tetcon[i].n[j], "      ", sizeof(nametyp));
  }
  dls_tetcon->anztetcon = 0;
  FORLIM = LINK->strukturdaten.anzatomlagen;
  for (i = 0; i < FORLIM; i++) {
    if (LINK->strukturdaten.atomlage[i]->koordinationszahl == 4) {
      *b = true;
      WITH2 = LINK->strukturdaten.atomlage[i]->atom[0];
      for (j = 0; j <= 3; j++) {
        if (koordinationszahl_des_atoms(&WITH2->bindung[j]) != 2)
          *b = false;
      }
      if (*b) {
        dls_tetcon->anztetcon++;
        WITH3 = &dls_tetcon->tetcon[dls_tetcon->anztetcon - 1];
        WITH3->a[0] = LINK->strukturdaten.atomlage[i]->grundatom;
        for (j = 1; j <= 4; j++) {
          berechne_nte_bindung(WITH3->a, j, &WITH3->a[j]);
          berechne_nte_bindung(&WITH3->a[j], 1L, &WITH3->a[j + 4]);
          if (gleiche_atome(&WITH3->a[j + 4], &WITH3->a[0]))
            berechne_nte_bindung(&WITH3->a[j], 2L, &WITH3->a[j + 4]);
        }
      }
    }
  }
  dls_distan->anzdistcon = 0;
  FORLIM = LINK->strukturdaten.anzatomlagen;
  for (i = 0; i < FORLIM; i++) {
    FORLIM1 = LINK->strukturdaten.atomlage[i]->koordinationszahl;
    for (j = 1; j <= FORLIM1; j++) {
      a1 = LINK->strukturdaten.atomlage[i]->grundatom;
      berechne_nte_bindung(&a1, j, &a2);
      *b = true;
      FORLIM2 = dls_distan->anzdistcon;
      for (k = 0; k < FORLIM2; k++) {
        WITH1 = &dls_distan->distcon[k];
        if (sym_aequival_distanzen(&a1, &a2, &WITH1->atom1, &WITH1->atom2) |
            sym_aequival_distanzen(&a1, &a2, &WITH1->atom2, &WITH1->atom1))
          *b = false;
      }
      if (*b) {
        if (dls_distan->anzdistcon < maxindepdistanzen) {
          dls_distan->anzdistcon++;
          dls_distan->distcon[dls_distan->anzdistcon - 1].atom1 = a1;
          dls_distan->distcon[dls_distan->anzdistcon - 1].atom2 = a2;
        } else {
          printf(" ERROR  too many symmetrically independent distances\n");
          printf(
            "        The program is compiled for a maximum of %ld independent distances.\n",
            (long)maxindepdistanzen);
          printf(
            "        Increase constant \"maxindepdistanzen\" and compile again.\n");
        }
      }
    }
  }
  dls_angle->anzanglcon = 0;
  FORLIM = LINK->strukturdaten.anzatomlagen;
  for (i = 0; i < FORLIM; i++) {
    FORLIM1 = LINK->strukturdaten.atomlage[i]->koordinationszahl;
    for (j = 1; j <= FORLIM1; j++) {
      for (l = 1; l < j; l++) {
        am = LINK->strukturdaten.atomlage[i]->grundatom;
        berechne_nte_bindung(&am, j, &a1);
        berechne_nte_bindung(&am, l, &a2);
        *b = true;
        FORLIM3 = dls_angle->anzanglcon;
        for (k = 0; k < FORLIM3; k++) {
          WITH4 = &dls_angle->anglcon[k];
          if (sym_aequival_winkel(&a1, &am, &a2, &WITH4->atom1, &WITH4->atomm,
                                  &WITH4->atom2) | sym_aequival_winkel(&a1, &am,
                &a2, &WITH4->atom2, &WITH4->atomm, &WITH4->atom1))
            *b = false;
        }
        if (*b) {
          if (dls_angle->anzanglcon < maxindepwinkel) {
            dls_angle->anzanglcon++;
            dls_angle->anglcon[dls_angle->anzanglcon - 1].atom1 = a1;
            dls_angle->anglcon[dls_angle->anzanglcon - 1].atom2 = a2;
            dls_angle->anglcon[dls_angle->anzanglcon - 1].atomm = am;
          } else {
            printf(" ERROR  too many symmetrically independent angles\n");
            printf(
              "        The program is compiled for a maximum of %ld independent angles.\n",
              (long)maxindepwinkel);
            printf(
              "        Increase constant \"maxindepwinkel\" and compile again.\n");
          }
        }
      }
    }
  }

  FORLIM = dls_angle->anzanglcon;
  /* bringe atom1 in Lage x,y,z */
  for (i = 0; i < FORLIM; i++) {
    WITH4 = &dls_angle->anglcon[i];
    atom0 = WITH4->atom1.dessen_atomart->deren_atomlage->grundatom;
    welchsymop(&WITH4->atom1, &atom0, &rginfo, &anzsymop, symeq0);
    alg_symmetrieoperation(&symeq0[0], &WITH4->atom1, &WITH4->atom1);
    alg_symmetrieoperation(&symeq0[0], &WITH4->atomm, &WITH4->atomm);
    alg_symmetrieoperation(&symeq0[0], &WITH4->atom2, &WITH4->atom2);
  }
  *b = true;
}

Local Void erzeuge_atomname(kt_atom *atoco, Char *neuername,
                            dls_symeq_karten *dls_symeq)
{
  long i;
  boolean b;
  kt_atom gruato;
  kt_atomlage *WITH1;
  long FORLIM;
  _REC_symeq *WITH2;

  memcpy(neuername, atoco->dessen_atomart->deren_atomlage->aname,
         sizeof(kt_atomlagename));

  memset(neuername + sizeof(kt_atomlagename), ' ',
         sizeof(nametyp) - sizeof(kt_atomlagename));
  gruato = atoco->dessen_atomart->deren_atomlage->grundatom;
  if (gleiche_atome(atoco, &gruato))
    return;
  WITH1 = atoco->dessen_atomart->deren_atomlage;
  if (WITH1->atomlagenr - 1 >= maxanzahlsymeqkarten)
    InternalError("maxanzahlsymeqkarten exceeded (see source code kriber.c)");
  b = true;
  FORLIM = dls_symeq->anzsymeq[WITH1->atomlagenr - 1];
  for (i = 0; i < FORLIM; i++) {
    if (gleiche_atome(atoco,
                      &dls_symeq->symeq[WITH1->atomlagenr - 1][i].atona)) {
      b = false;
      memcpy(neuername,
             dls_symeq->symeq[WITH1->atomlagenr - 1][i].zweitername,
             sizeof(nametyp));
    }
  }
  if (!b)
    return;
  WITH2 = &dls_symeq->symeq[WITH1->atomlagenr - 1]
    [dls_symeq->anzsymeq[WITH1->atomlagenr - 1]];
  WITH2->atona = *atoco;
  memcpy(WITH2->erstername, neuername, sizeof(nametyp));
  i = 1;
  while (WITH2->erstername[i - 1] != ' ')
    i++;
  dls_symeq->anzsymeq[WITH1->atomlagenr - 1]++;
  neuername[i - 1] = '*';
  neuername[i] = (Char)(dls_symeq->anzsymeq[WITH1->atomlagenr - 1] + 'A' - 1);
  memcpy(WITH2->zweitername, neuername, sizeof(nametyp));
}

Local Void loesche_redundante_karten(dls_tetcon_karten *dls_tetcon,
                                     dls_distan_karten *dls_distan,
                                     dls_angle_karten *dls_angle)
{
  long i, j, k, l, m;
  boolean b;
  long FORLIM, FORLIM1;
  _REC_distcon *WITH2;
  _REC_tetcon *WITH3;
  _REC_anglcon *WITH4;

  k = 0;
  FORLIM = dls_distan->anzdistcon;
  for (i = 0; i < FORLIM; i++) {
    b = true;
    FORLIM1 = dls_tetcon->anztetcon;
    for (j = 0; j < FORLIM1; j++) {
      WITH2 = &dls_distan->distcon[i];
      WITH3 = &dls_tetcon->tetcon[j];
      for (l = 1; l <= 4; l++) {
        if (sym_aequival_distanzen(&WITH2->atom1, &WITH2->atom2, &WITH3->a[0],
              &WITH3->a[l]) | sym_aequival_distanzen(&WITH2->atom1,
              &WITH2->atom2, &WITH3->a[l], &WITH3->a[0]))
          b = false;
      }
    }
    if (b) {
      k++;
      dls_distan->distcon[k - 1] = dls_distan->distcon[i];
    }
  }
  dls_distan->anzdistcon = k;
  k = 0;
  FORLIM = dls_angle->anzanglcon;
  for (i = 0; i < FORLIM; i++) {
    b = true;
    FORLIM1 = dls_tetcon->anztetcon;
    for (j = 0; j < FORLIM1; j++) {
      WITH4 = &dls_angle->anglcon[i];
      WITH3 = &dls_tetcon->tetcon[j];
      for (l = 1; l <= 4; l++) {
        if (sym_aequival_winkel(&WITH4->atom1, &WITH4->atomm, &WITH4->atom2,
              &WITH3->a[0], &WITH3->a[l],
              &WITH3->a[l + 4]) | sym_aequival_winkel(&WITH4->atom1,
              &WITH4->atomm, &WITH4->atom2, &WITH3->a[l + 4], &WITH3->a[l],
              &WITH3->a[0]))
          b = false;
        for (m = 1; m < l; m++) {
          if (sym_aequival_winkel(&WITH4->atom1, &WITH4->atomm, &WITH4->atom2,
              &WITH3->a[l], &WITH3->a[0],
              &WITH3->a[m]) | sym_aequival_winkel(&WITH4->atom1, &WITH4->atomm,
              &WITH4->atom2, &WITH3->a[m], &WITH3->a[0], &WITH3->a[l]))
            b = false;
        }
      }
    }
    if (b) {
      k++;
      dls_angle->anglcon[k - 1] = dls_angle->anglcon[i];
    }
  }
  dls_angle->anzanglcon = k;
}


Static Void write_input_file(_TEXT *f, long num,
                             kt_strukturdaten *strukturdaten_,
                             _TEXT *distdat)
{
  struct LOC_write_input_file V;
  dls_tetcon_karten dls_tetcon;
  dls_distan_karten dls_distan;
  dls_angle_karten dls_angle;
  dls_symeq_karten dls_symeq;
  dt_geg_distanz geg_dist;
  dt_geg_winkel geg_wink;
  long i, j;
  boolean b;
  vektor v;
  long FORLIM;
  _REC_tetcon *WITH;
  _REC_distcon *WITH1;
  _REC_anglcon *WITH2;

  int Polar[3];


  for (i = 0; i < 3; i++)
  {
    v[0] = (i == 0 ? 1 : 0);
    v[1] = (i == 1 ? 1 : 0);
    v[2] = (i == 2 ? 1 : 0);
    if (polare_richtung(v, &rginfo)) Polar[i] = 1;
    else                             Polar[i] = 0;
  }

  V.strukturdaten = *strukturdaten_;
  read_bindungs_groessen(distdat, geg_dist, geg_wink);
  berechne_dls_input(&dls_tetcon, &dls_distan, &dls_angle, &dls_symeq, &b, &V);
  if (b) {
    if (num == 1) {
      loesche_redundante_karten(&dls_tetcon, &dls_distan, &dls_angle);
      FORLIM = dls_tetcon.anztetcon;
      for (i = 0; i < FORLIM; i++) {
        WITH = &dls_tetcon.tetcon[i];
        for (j = 0; j <= 8; j++)
          erzeuge_atomname(&WITH->a[j], WITH->n[j], &dls_symeq);
      }
    }
    FORLIM = dls_distan.anzdistcon;
    for (i = 0; i < FORLIM; i++) {
      WITH1 = &dls_distan.distcon[i];
      erzeuge_atomname(&WITH1->atom1, WITH1->atomname1, &dls_symeq);
      erzeuge_atomname(&WITH1->atom2, WITH1->atomname2, &dls_symeq);
    }
    if (num == 1 || num == 2) {
      FORLIM = dls_angle.anzanglcon;
      for (i = 0; i < FORLIM; i++) {
        WITH2 = &dls_angle.anglcon[i];
        erzeuge_atomname(&WITH2->atom1, WITH2->atomname1, &dls_symeq);
        erzeuge_atomname(&WITH2->atom2, WITH2->atomname2, &dls_symeq);
      }
    }
    if (num == 3) {
      FORLIM = dls_angle.anzanglcon;
      for (i = 0; i < FORLIM; i++) {
        WITH2 = &dls_angle.anglcon[i];
        erzeuge_atomname(&WITH2->atom1, WITH2->atomname1, &dls_symeq);
        erzeuge_atomname(&WITH2->atomm, WITH2->atomnamem, &dls_symeq);
        erzeuge_atomname(&WITH2->atom2, WITH2->atomname2, &dls_symeq);
      }
    }
    if (num == 1 || num == 2) {
      fprintf(f->f, "TITLE %.80s\n", V.strukturdaten.name);
      if (num == 1)
        fprintf(f->f, "DLS-76  -5    30   2 0         2                1\n");
      if (num == 2)
        fprintf(f->f, "DLS-76  -5    30     0         2                1\n");
      if (!strncmp(V.strukturdaten.kristallsystem, "TRIG",
                   sizeof(st_kristallsystem)))
        memcpy(V.strukturdaten.kristallsystem, "HEX ",
               sizeof(st_kristallsystem));
      fprintf(f->f,
              "CELL   %.4s         %10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n",
              V.strukturdaten.kristallsystem,
              V.strukturdaten.zellpar.laenge[0],
              V.strukturdaten.zellpar.laenge[1],
              V.strukturdaten.zellpar.laenge[2],
              V.strukturdaten.zellpar.winkel[0],
              V.strukturdaten.zellpar.winkel[1],
              V.strukturdaten.zellpar.winkel[2]);
      write_dls_atom_karten(f, &V.strukturdaten, Polar);
      write_dls_symeq_karten(f, &dls_symeq, &V);
      if (num == 1)
        write_dls_bondis_karten(f, &dls_tetcon, geg_dist, geg_wink);
      write_dls_distan_karten(f, geg_dist, &dls_distan);
      write_dls_angle_karten(f, geg_wink, geg_dist, &dls_angle);
    }
    if (num == 1)
      write_dls_tetcon_karten(f, &dls_tetcon);
    if (num == 3) {
      if (*f->name != '\0') {
        if (f->f != NULL)
          f->f = freopen(f->name, "w", f->f);
        else
          f->f = fopen(f->name, "w");
      } else {
        if (f->f != NULL)
          rewind(f->f);
        else
          f->f = tmpfile();
      }
      if (f->f == NULL)
        _EscIO(FileNotFound);
      SETUPBUF(f->f, Char);
      write_loadat_atom_karten(f, &V.strukturdaten);
      write_dls_symeq_karten(f, &dls_symeq, &V);
      write_loadat_distan_karten(f, geg_dist, &dls_distan);
      write_loadat_angle_karten(f, geg_wink, &dls_angle);
    }
    fprintf(f->f, "END\n");
    fprintf(f->f, "FINISH\n");
  }
}


static void determine_isym_ilatt_ia_ib_ic(kt_atom *a,
                                          long *isym, long *ilatt,
                                          long *ia, long *ib, long *ic)
{
  int                nt, nz, ns, it, iz, is, i;
  double             dc;
  vektor             gc, ac, sc;
  gittervektor       g;
  symmetrieoperator  inv, invr;


#ifdef JUNK
  {
    _TEXT  t;

    t.f = stdout;
    t.name[0] = '\0';

    write_atom_koordinaten(&t, &a->dessen_atomart->deren_atomlage->grundatom);
    write_atom_koordinaten(&t, a);
    putc(' ', stdout);
  }
#endif

  berechne_atom_koordinaten(&a->dessen_atomart->deren_atomlage->grundatom, gc);
  berechne_atom_koordinaten(a, ac);

  set_einheits_symop(&inv);
  for (i = 0; i < 3; i++)
    inv.r[i][i] = -1;

  nt = rginfo.anztranslationen + 1;

  if (rginfo.zentrosymm)
    nz = 2;
  else
    nz = 1;

  ns = rginfo.anz_rotmat;
  if (rginfo.zentrosymm)
    ns /= 2;

  for (it = 0; it < nt; it++)
  {
    for (iz = 0; iz < nz; iz++)
    {
      for (is = 0; is < ns; is++)
      {
        if (iz)
        {
          mult_symmetrieoperatoren(&inv, rginfo.symop[is], &invr);
          mult_symop_mit_vektor(&invr, gc, sc);
        }
        else
          mult_symop_mit_vektor(rginfo.symop[is], gc, sc);

        if (it)
          for (i = 0; i < 3; i++)
            sc[i] += rginfo.translation[it - 1][i];

        for (i = 0; i < 3; i++)
        {
          dc = sc[i] - ac[i];

                    g[i] = (long) dc;
          dc = dc - g[i];
          for (;;) if (dc >= 1.) { dc -= 1.; g[i]++; } else break;
          for (;;) if (dc <  0.) { dc += 1.; g[i]--; } else break;
                   if (dc > 0.5) { dc -= 1.; g[i]++; }

          if (dc < -1.e-4 || dc > 1.e-4)
            break;
        }

        if (i == 3)
        {
          *isym = is + 1;
          if (iz)
            (*isym) *= -1;

          *ilatt = it;
          *ia    = -g[0];
          *ib    = -g[1];
          *ic    = -g[2];

          return;
        }
      }
    }
  }

  InternalError("determine_isym_ilatt_ia_ib_ic()");
}


Static Void write_GSAS_soft_constraints(_TEXT *t,
                                        kt_strukturdaten *strukturdaten_,
                                        _TEXT *distdat,
                                        int phase_no, int hst_no)
{
  struct LOC_write_input_file V;

  dls_tetcon_karten  dls_tetcon;
  dls_distan_karten  dls_distan;
  dls_angle_karten   dls_angle;
  dls_symeq_karten   dls_symeq;
  dt_geg_distanz     geg_dist;
  dt_geg_winkel      geg_wink;

  long          i;
  long          isym, ilatt, ia, ib, ic;
  boolean       b, bb;
  double        d1, d2, d1v, d2v, w, wv, distanz, dvar;
  _REC_distcon  *dcon;
  _REC_anglcon  *acon;

  int  iHST_hhBND;


  V.strukturdaten = *strukturdaten_;
  read_bindungs_groessen(distdat, geg_dist, geg_wink);
  berechne_dls_input(&dls_tetcon, &dls_distan, &dls_angle, &dls_symeq, &b, &V);

  if (! b) return;

  fprintf(t->f, " EXPR  HTYP%d  PXC  RSN\n", phase_no);
  fprintf(t->f, " EXPR  NHST %5d\n", hst_no);

  fprintf(t->f, "HST %2d NBNDS%5ld\n",
    hst_no, dls_distan.anzdistcon + dls_angle.anzanglcon);

  iHST_hhBND = 0;

  dcon = dls_distan.distcon;

  for (i = 0; i < dls_distan.anzdistcon; i++, dcon++)
  {
    bestimme_geg_distanzen(&dcon->atom1, &dcon->atom2, geg_dist, &distanz,
                           &dvar, &b);

    determine_isym_ilatt_ia_ib_ic(&dcon->atom2, &isym, &ilatt, &ia, &ib, &ic);

    fprintf(t->f, "HST %2dBND%3d %2d%3ld%3ld %2ld %2ld %2ld %2ld %2ld",
      hst_no, ++iHST_hhBND, phase_no,
      dcon->atom1.dessen_atomart->deren_atomlage->atomlagenr,
      dcon->atom2.dessen_atomart->deren_atomlage->atomlagenr,
      isym, ilatt, ia, ib, ic);

    if (! b) { distanz = 1.; dvar = 0.; }

    fprintf(t->f, "%6.3f%6.3f\n", distanz, dvar);
  }

  acon = dls_angle.anglcon;

  for (i = 0; i < dls_angle.anzanglcon; i++, acon++)
  {
    bestimme_geg_distanzen(&acon->atom1, &acon->atomm,
                           geg_dist, &d1, &d1v, &bb);
    bestimme_geg_distanzen(&acon->atom2, &acon->atomm,
                           geg_dist, &d2, &d2v, &b);
    if (! b) bb = false;

    bestimme_geg_winkel(&acon->atom1, &acon->atomm, &acon->atom2,
                        geg_wink, &w, &wv, &b);
    if (! b) bb = false;

    determine_isym_ilatt_ia_ib_ic(&acon->atom2, &isym, &ilatt, &ia, &ib, &ic);

    fprintf(t->f, "HST %2dBND%3d %2d%3ld%3ld %2ld %2ld %2ld %2ld %2ld",
      hst_no, ++iHST_hhBND, phase_no,
      acon->atom1.dessen_atomart->deren_atomlage->atomlagenr,
      acon->atom2.dessen_atomart->deren_atomlage->atomlagenr,
      isym, ilatt, ia, ib, ic);

    if (bb)
      berechne_t_t_distanzen(d1, d2, w, d1v, d2v, wv, &distanz, &dvar);
    else {
      distanz = 1.; dvar = 0.; }

    fprintf(t->f, "%6.3f%6.3f\n", distanz, dvar);
  }
}


Static Void m_initialize(_TEXT *f) // Stefs 05-07-2012 -> initializing necessary to prevent segmentation fault
{
  if (f->name != '\0') {
    if (f->f != NULL)
      f->f = freopen(f->name, "w", f->f);
    else
      f->f = fopen(f->name, "w");
  } else {
    if (f->f != NULL)
      rewind(f->f);
    else
      f->f = tmpfile();
  }
  if (f->f == NULL)
    _EscIO(FileNotFound);
  SETUPBUF(f->f, Char);
}


Static Void m_initialisiere(kt_strukturdaten *strukturdaten)
{
  memcpy(strukturdaten->raumgruppe, "                              ",
         sizeof(st_raumgruppename));
  berechnet &= ~(1 << 0);
  berechnet &= ~(1 << ((long)distanzen - (long)struktur));
  berechnet &= ~(1 << ((long)bindungen - (long)struktur));
  strudatnew_offen = false;

  //// XXX this should be redundant now? /ss
  //if (*list.name != '\0') {
  //  if (list.f != NULL)
  //    list.f = freopen(list.name, "w", list.f);
  //  else
  //    list.f = fopen(list.name, "w");
  //} else {
  //  if (list.f != NULL)
  //    rewind(list.f);
  //  else
  //    list.f = tmpfile();
  //}
  //if (list.f == NULL)
  //  _EscIO(FileNotFound);
  //SETUPBUF(list.f, Char);

  printf(
    " ******************************************************************************\n");
  printf(
    " *                                                                            *\n");
  printf(
    " *   Program  KRIBER        Version 1.2    (January 2013)                     *\n");
  printf(
    " *   ===============                                                          *\n");
  printf(
    " *                                                                            *\n");
  printf(
    " *   an interactive PASCAL program to                                         *\n");
  printf(
    " *    - calculate distances and angles                                        *\n");
  printf(
    " *    - generate input files for the programs DLS-76 and LOADAT (XRS-82)      *\n");
  printf(
    " *    - calculate coordination sequences and loop configurations              *\n");
  printf(
    " *                                                                            *\n");
  printf(
    " *   Author:   Roland Bialek, Institut fuer Kristallografie und Petrografie   *\n");
  printf(
    " *             ETH-Z, CH-8092 Zuerich, Switzerland                            *\n");
  printf(
    " *                                                                            *\n");
  printf(
    " *     ==> type help to get a summary of the commands                         *\n");
  printf(
    " *              quit or exit to leave the program                             *\n");
  printf(
    " *                                                                            *\n");
  printf(
    " ******************************************************************************\n");
}


static void type_return_to_continue(void)
{
  printf("\n");
  printf("        Type <ret> to continue\n");
  scanf("%*[^\n]");
  getchar();
}


Static Void m_write_hilfe(void)
{
  printf("\n *******************************************************************\n");
  printf(" *                                                                 *\n");
  printf(" *          S U M M A R Y    O F   A L L   C O M M A N D S         *\n");
  printf(" *                                                                 *\n");
  printf(" *******************************************************************\n\n");
  printf(" read      crystal structure data                         ==> reacs \n");
  printf(" read      data from first entry in file STRUDAT          ==> reafcs\n");
  printf(" read      data from next entry in file STRUDAT           ==> reancs\n");
  printf(" write     structure data on file STRUDATNEW              ==> wrie  \n");
  printf(" write     primitive structure data on file STRUDATNEW    ==> wriep \n");
  printf(" write     primitive structure with renumbered elements   ==> wriepn\n");
  printf(" write     primitive structure with resized cell          ==> wriepr\n");
  printf(" display   all entries in file STRUDAT                    ==> dise  \n");
  printf("\n");
  printf(" display   crystal structure                              ==> discs \n");
  printf(" write     crystal structure on file  LIST                ==> wrics \n");
  printf(" display   crystal structure in CSSR  format              ==> discssr\n");
  printf(" write     crystal structure in CSSR  format on file LIST ==> wricssr\n");
  printf(" display   crystal structure in XTL   format              ==> disxtl\n");
  printf(" write     crystal structure in XTL   format on file LIST ==> wrixtl\n");
  printf(" display   crystal structure in CIF   format              ==> discif\n");
  printf(" write     crystal structure to structure.cif             ==> wricif\n");
  //type_return_to_continue();
  printf("                                *******                             \n\n");
  printf(" display   crystal structure in Focus format              ==> disfocus\n");
  printf(" write     crystal structure in Focus format on file LIST ==> wrifocus\n");
  printf(" display   crystal structure in MVM   format              ==> dismvm\n");
  printf(" write     crystal structure in MVM   format on file LIST ==> wrimvm\n");
  printf(" display   crystal structure in MVM   format              ==> dismvm#\n");
  printf(" write     crystal structure in MVM   format on file LIST ==> wrimvm#\n");
  printf(" display   crystal structure suitable for GSAS L A I      ==> disglai\n");
  printf(" write     crystal structure suitable for GSAS L A I      ==> wriglai\n");
  printf(" display   GSAS crystal structure records                 ==> disgcs\n");
  printf(" write     GSAS crystal structure records    on file LIST ==> wrigcs\n");
  printf(" display   GSAS soft constraint records                   ==> disgsc\n");
  printf(" write     GSAS soft constraint records      on file LIST ==> wrigsc\n");
  printf("\n");
  printf(" add       atom position                                  ==> adda  \n");
  printf(" delete    atom position                                  ==> dela  \n");
  printf(" display   space group                                    ==> diss  \n");
  printf(" write     space group on file LIST                       ==> wris  \n");
  printf(" display   coordinates                                    ==> disc  \n");
  printf(" write     coordinates on file LIST                       ==> wric  \n");
  printf(" calculate distances                                      ==> cald  \n");
  //type_return_to_continue();
  printf("                                *******                             \n\n");
  printf(" display   all distances                                  ==> disad \n");
  printf(" display   all \"unique\" distances                         ==> disadu\n");
  printf(" display   distances from one atom                        ==> disd  \n");
  printf(" display   distances and angles                           ==> disda \n");
  printf(" display   distances and angles unformatted               ==> disdau\n");
  printf(" display   distances and interatomic vectors              ==> disdv \n");
  printf(" write     all distances on file LIST                     ==> wriad \n");
  printf(" write     all \"unique\" distances on file LIST            ==> wriadu\n");
  printf(" write     distances from one atom on file LIST           ==> wrid  \n");
  printf(" write     distances and angles on file LIST              ==> wrida \n");
  printf(" write     distances and angles unformatted  on file LIST ==> wridau\n");
  printf(" write     distances and interatomic vectors on file LIST ==> wridv \n");
  printf("\n");
  printf(" create    all bonds using shortest distances             ==> creab \n");
  printf(" create    all bonds between different elements only      ==> creabd\n");
  printf(" create    a bond                                         ==> creb  \n");
  printf(" delete    a bond                                         ==> delb  \n");
  printf(" delete    all bonds                                      ==> delab \n");
  printf(" display   bonds from one atom                            ==> disb  \n");
  printf(" display   all bonds                                      ==> disab \n");
  //type_return_to_continue();
  printf("                               *******                              \n\n");
  printf(" write     bonds from one atom on file LIST               ==> wrib  \n");
  printf(" write     all bonds on file LIST                         ==> wriab \n");
  printf("\n");
  printf(" add       O atoms (oxygen bridges)                       ==> addo  \n");
  printf(" delete    oxygen bridges, connect non O-atoms directly   ==> delo  \n");
  printf(" calculate coordination sequences                         ==> calcs \n");
  printf(" write     coordination sequences    on file LIST         ==> wricsq\n");
  printf(" calculate loop configuration                             ==> callc \n");
  printf(" write     loop configuration        on file LIST         ==> wrilc \n");
  printf(" determine topology of the framework                      ==> dett  \n");
  printf(" write     topology of the framework on file LIST         ==> writ  \n");
  printf("\n");
  printf(" write     input file for DLS-76 on file DLSINP           ==> wriid \n");
  printf(" write     DLS input file w/ tetcon cards to file DLSINP  ==> writc \n");
  printf(" write     input file for LOADAT on file LOADATINPUT      ==> wriil \n");
  printf("\n");
  printf(" display   a summary of the commands                      ==> help  \n");
  printf(" quit or exit the program                                 ==> quit  \n");
  printf("                                                       or ==> exit  \n");
  printf("\n");
  printf(                                          
    " ******************************************** updated 31-01-2013 ***\n");
  /* Fuer Kontrollzwecke  cons : kontrolliere Raumgruppe                      */
  /* Fuer Kontrollzwecke  conb : schreibe alle Bindungen                      */
  /* Fuer Kontrollzwecke  calda: berechne Distanzen algebraisch               */
}


static int m_read_strukturdaten(void)
{
  long i;
  boolean eintrag_gefunden;
  _TEXT TEMP;

  if (*strudat.name != '\0') {
    if (strudat.f != NULL)
      strudat.f = freopen(strudat.name, "r", strudat.f);
    else
      strudat.f = fopen(strudat.name, "r");
  } else
    rewind(strudat.f);
  if (strudat.f == NULL)
    _EscIO(FileNotFound);
  RESETBUF(strudat.f, Char);
  printf(" Enter compound identification code to read data from the file\n");
  printf(" STRUDAT or <ret> to enter the data interactively:");
  if (P_eof(stdin)) {
    rewind(stdin);
    TEMP.f = stdin;
    *TEMP.name = '\0';
    read_strukturdaten(&TEMP, &strukturdaten, true);
    memcpy(strukturdaten.abk, "???             ", sizeof(abkbezeichnung));
    berechnet |= 1 << 0;
    berechnet &= ~(1 << ((long)distanzen - (long)struktur));
    berechnet &= ~(1 << ((long)bindungen - (long)struktur));
    return 0;
  }
  for (i = 0; i <= 15; i++)
    strukturdaten.abk[i] = ' ';
  i = 0;
  while (!P_eoln(stdin) && i < 16) {
    i++;
    strukturdaten.abk[i - 1] = getchar();
    if (strukturdaten.abk[i - 1] == '\n')
      strukturdaten.abk[i - 1] = ' ';
  }
  if (strukturdaten.abk[15] != ' ')
    printf(" ERROR  compound identification code too long\n");
  scanf("%*[^\n]");
  getchar();
  if (*strudat.name != '\0') {
    if (strudat.f != NULL)
      strudat.f = freopen(strudat.name, "r", strudat.f);
    else
      strudat.f = fopen(strudat.name, "r");
  } else
    rewind(strudat.f);
  if (strudat.f == NULL)
    _EscIO(FileNotFound);
  RESETBUF(strudat.f, Char);
  suche_strukturdaten_eintrag(&strudat, strukturdaten.abk, &eintrag_gefunden);
  if (!eintrag_gefunden)
    return -1;
  read_strukturdaten(&strudat, &strukturdaten, false);
  berechnet |= 1 << 0;
  berechnet &= ~(1 << ((long)distanzen - (long)struktur));
  berechnet &= ~(1 << ((long)bindungen - (long)struktur));
  return 0;
}


Static Void m_read_erste_strukturdaten(void)
{
  boolean b;

  if (*strudat.name != '\0') {
    if (strudat.f != NULL)
      strudat.f = freopen(strudat.name, "r", strudat.f);
    else
      strudat.f = fopen(strudat.name, "r");
  } else
    rewind(strudat.f);
  if (strudat.f == NULL)
    _EscIO(FileNotFound);
  RESETBUF(strudat.f, Char);
  suche_naechsten_eintrag(&strudat, strukturdaten.abk, &b);
  if (!b) {
    printf(" ERROR  no entry found on file STRUDAT\n");
    return;
  }
  printf(" ... entry %.16s found\n", strukturdaten.abk);
  read_strukturdaten(&strudat, &strukturdaten, false);
  berechnet |= 1 << 0;
  berechnet &= ~(1 << ((long)distanzen - (long)struktur));
  berechnet &= ~(1 << ((long)bindungen - (long)struktur));
}


Static Void m_read_naechste_strukturdaten(void)
{
  boolean b;

  suche_naechsten_eintrag(&strudat, strukturdaten.abk, &b);
  if (!b) {
    printf(" ERROR  no entry found on file STRUDAT\n");
    return;
  }
  printf(" ... entry %.16s found\n", strukturdaten.abk);
  read_strukturdaten(&strudat, &strukturdaten, false);
  berechnet |= 1 << 0;
  berechnet &= ~(1 << ((long)distanzen - (long)struktur));
  berechnet &= ~(1 << ((long)bindungen - (long)struktur));
}


Static Void m_berechne_distanzen(long anzdistanzen)
{
  long i, FORLIM;

  if ((berechnet & (1 << 0)) == 0)
    if (m_read_strukturdaten() != 0) return;
  FORLIM = strukturdaten.anzatomlagen;
  for (i = 0; i < FORLIM; i++)
    berechne_distanzen(&strukturdaten.atomlage[i]->grundatom, &strukturdaten,
                       anzdistanzen, 4.0);
  berechnet |= 1 << ((long)distanzen - (long)struktur);
}


Static Void m_berechne_bindungen(int OnlyDifferentElements)
{
  if ((berechnet & (1 << ((long)distanzen - (long)struktur))) == 0)
    m_berechne_distanzen((long)maxanzahldistanzen);
  berechne_alle_bindungen(&strukturdaten, OnlyDifferentElements);
  berechnet |= 1 << ((long)bindungen - (long)struktur);
}


Static Void m_write_strukturdaten(_TEXT *f, kt_strukturdaten *strukturdaten)
{
  if ((berechnet & (1 << 0)) == 0)
    if (m_read_strukturdaten() != 0) return;
  write_strukturdaten(f, strukturdaten);
}


Static Void m_write_cif(_TEXT *f, kt_strukturdaten *strukturdaten, st_raumgruppe *rginfo)
{
  if ((berechnet & (1 << 0)) == 0)
    if (m_read_strukturdaten() != 0) return;
  write_cif(f, strukturdaten, rginfo);
}


Static Void m_write_cssr(_TEXT *f, kt_strukturdaten *strukturdaten)
{
  if ((berechnet & (1 << 0)) == 0)
    if (m_read_strukturdaten() != 0) return;
  write_cssr(f, strukturdaten);
}


Static Void m_write_xtl(_TEXT *f, kt_strukturdaten *strukturdaten)
{
  if ((berechnet & (1 << 0)) == 0)
    if (m_read_strukturdaten() != 0) return;
  write_xtl(f, strukturdaten);
}


Static Void m_write_mvm(_TEXT *f, kt_strukturdaten *strukturdaten, int flag_num)
{
  if ((berechnet & (1 << 0)) == 0)
    if (m_read_strukturdaten() != 0) return;
  write_mvm(f, strukturdaten, flag_num);
}


Static Void m_write_focus(_TEXT *f, kt_strukturdaten *strukturdaten)
{
  if ((berechnet & (1 << 0)) == 0)
    if (m_read_strukturdaten() != 0) return;
  write_focus(f, strukturdaten);
}


static int get_GSAS_phase_no(void)
{
  int  phase_no, c;


  for (;;)
  {
    printf(" Enter GSAS Phase No: ");

    if (scanf("%d", &phase_no) != 1 || phase_no < 1 || phase_no > 9)
      phase_no = 0;

    while ((c = getchar()) != EOF && c != '\n')
      if (! isspace(c)) phase_no = 0;

    if (phase_no)
      break;
  }

  return phase_no;
}


static int get_GSAS_hst_no(void)
{
  int  hst_no, c;


  for (;;)
  {
    printf(" Enter GSAS HST No: ");

    if (scanf("%d", &hst_no) != 1 || hst_no < 1 || hst_no > 99)
      hst_no = 0;

    while ((c = getchar()) != EOF && c != '\n')
      if (! isspace(c)) hst_no = 0;

    if (hst_no)
      break;
  }

  return hst_no;
}


static void m_write_GSAS_crystal_structure(_TEXT *f,
                                           kt_strukturdaten *strukturdaten,
                                           int phase_no)
{
  if ((berechnet & (1 << 0)) == 0)
    if (m_read_strukturdaten() != 0) return;

  if (phase_no == 0)
    phase_no = get_GSAS_phase_no();

  write_GSAS_crystal_structure(f, strukturdaten, phase_no);
}


static void m_write_GSAS_L_A_I(_TEXT *f,
                               kt_strukturdaten *strukturdaten)
{
  if ((berechnet & (1 << 0)) == 0)
    if (m_read_strukturdaten() != 0) return;

  write_GSAS_L_A_I(f, strukturdaten);
}


static Void m_write_GSAS_soft_constraints(_TEXT *t, int phase_no, int hst_no)
{
  if ((berechnet & (1 << ((long)bindungen - (long)struktur))) == 0)
    m_berechne_bindungen(0);

  if (phase_no == 0)
    phase_no = get_GSAS_phase_no();

  if (hst_no == 0)
    hst_no = get_GSAS_hst_no();

  write_GSAS_soft_constraints(t, &strukturdaten, &distdat, phase_no, hst_no);
}


Static Void m_contr_st_raumgruppe(void)
{
  vektor v;
  _TEXT TEMP;

  TEMP.f = stdin;
  *TEMP.name = '\0';
  readln_raumgruppename(&TEMP, strukturdaten.raumgruppe);
  read_raumgruppe(&symdat, strukturdaten.raumgruppe, &rginfo);
  memcpy(strukturdaten.kristallsystem, rginfo.kristallsystem,
         sizeof(st_kristallsystem));
  TEMP.f = stdout;
  *TEMP.name = '\0';
  write_st_raumgruppe(&TEMP, &rginfo);
  v[0] = 1.0;
  v[1] = 0.0;
  v[2] = 0.0;
  if (polare_richtung(v, &rginfo)) {
    printf(" possible shift in ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    write_vektor(&TEMP, v);
    putchar('\n');
  }
  v[0] = 0.0;
  v[1] = 1.0;
  v[2] = 0.0;
  if (polare_richtung(v, &rginfo)) {
    printf(" possible shift in ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    write_vektor(&TEMP, v);
    putchar('\n');
  }
  v[0] = 0.0;
  v[1] = 0.0;
  v[2] = 1.0;
  if (!polare_richtung(v, &rginfo))
    return;
  printf(" possible shift in ");
  TEMP.f = stdout;
  *TEMP.name = '\0';
  write_vektor(&TEMP, v);
  putchar('\n');
}


Static Void m_write_koordinaten(_TEXT *f)
{
  long i, FORLIM;

  if ((berechnet & (1 << 0)) == 0)
    if (m_read_strukturdaten() != 0) return;
  FORLIM = strukturdaten.anzatomlagen;
  for (i = 0; i < FORLIM; i++)
    write_atomlage_koordinaten(f, strukturdaten.atomlage[i]);
}


Static Void m_write_distanzen(_TEXT *f)
{
  kt_atom a;
  _TEXT TEMP;

  if ((berechnet & (1 << ((long)distanzen - (long)struktur))) == 0)
    m_berechne_distanzen((long)maxanzahldistanzen);
  printf(" display distances from atom:");
  TEMP.f = stdin;
  *TEMP.name = '\0';
  readln_atom(&TEMP, &a);
  if (a.dessen_atomart != NULL)
    write_distanzen(f, &a);
  else
    printf(" ERROR  cannot display distances\n");
}


Static Void m_write_alle_distanzen(_TEXT *f, int uniq)
{
  long i;

  if ((berechnet & (1 << ((long)distanzen - (long)struktur))) == 0)
    m_berechne_distanzen((long)maxanzahldistanzen);

  for (i = 0; i < strukturdaten.anzatomlagen; i++)
  {
    if (uniq == 0)
      write_distanzen(f, &strukturdaten.atomlage[i]->grundatom);
    else
      write_distanzen_uniq(f, &strukturdaten.atomlage[i]->grundatom);
  }
}


Static Void m_write_distanzen_und_winkel(_TEXT *f, int fmt)
{
  long i, FORLIM;

  if ((berechnet & (1 << ((long)distanzen - (long)struktur))) == 0)
    m_berechne_distanzen((long)maxanzahldistanzen);
  FORLIM = strukturdaten.anzatomlagen;
  for (i = 0; i < FORLIM; i++)
  {
    if (fmt)
      write_f_distanzen_und_winkel(f, &strukturdaten.atomlage[i]->grundatom);
    else
      write_u_distanzen_und_winkel(f, &strukturdaten.atomlage[i]->grundatom);
  }
}


Static Void m_write_distanzen_und_vektoren(_TEXT *f)
{
  long i, FORLIM;

  if ((berechnet & (1 << ((long)distanzen - (long)struktur))) == 0)
    m_berechne_distanzen((long)maxanzahldistanzen);
  FORLIM = strukturdaten.anzatomlagen;
  for (i = 0; i < FORLIM; i++)
    write_distanzen_und_vektoren(f, &strukturdaten.atomlage[i]->grundatom);
}


Static Void m_write_bindungen(_TEXT *f)
{
  kt_atom a;
  _TEXT TEMP;

  if ((berechnet & (1 << ((long)bindungen - (long)struktur))) == 0)
    m_berechne_bindungen(0);
  printf(" display bonds from atom:");
  TEMP.f = stdin;
  *TEMP.name = '\0';
  readln_atom(&TEMP, &a);
  if (a.dessen_atomart != NULL)
    write_bindungen(f, &a);
  else
    printf(" ERROR  cannot display bonds\n");
}


Static Void m_write_alle_bindungen(_TEXT *f)
{
  long i, FORLIM;

  if ((berechnet & (1 << ((long)bindungen - (long)struktur))) == 0)
    m_berechne_bindungen(0);
  FORLIM = strukturdaten.anzatomlagen;
  for (i = 0; i < FORLIM; i++)
    write_bindungen(f, &strukturdaten.atomlage[i]->grundatom);
}


Static Void m_write_bindungs_matrix(_TEXT *f)
{
  long i, FORLIM;

  FORLIM = strukturdaten.anzatomlagen;
  for (i = 0; i < FORLIM; i++)
    write_bindungs_matrix(f, strukturdaten.atomlage[i]);
}


Static Void m_DLSinput(void)
{
  Char c;
  long n;

  if ((berechnet & (1 << ((long)bindungen - (long)struktur))) == 0)
        /*#### BK ut */
          m_berechne_bindungen(0);
  if (*dlsinp.name != '\0') {
    if (dlsinp.f != NULL)
      dlsinp.f = freopen(dlsinp.name, "w", dlsinp.f);
    else
      dlsinp.f = fopen(dlsinp.name, "w");
  } else {
    if (dlsinp.f != NULL)
      rewind(dlsinp.f);
    else
      dlsinp.f = tmpfile();
  }
  if (dlsinp.f == NULL)
    _EscIO(FileNotFound);
  SETUPBUF(dlsinp.f, Char);
  n = 2;
  //printf(" DLSinput with TETCON cards? (y/n): => n");
  // if (!P_eoln(stdin)) {   /*#### BK ut */
  //   scanf("%c%*[^\n]", &c);
  //   getchar();
  //   if (c == '\n')
  //     c = ' ';
  //   if (c == 'j' || c == 'J' || c == 'y' || c == 'Y')
  //     n = 1;
  //   if (c == 'n' || c == 'N')
  //     n = 2;
  // }
  write_input_file(&dlsinp, n, &strukturdaten, &distdat);
}

Static Void m_DLSinput_tetcon(void)
{
  Char c;
  long n;

  if ((berechnet & (1 << ((long)bindungen - (long)struktur))) == 0)
        /*#### BK ut */
          m_berechne_bindungen(0);
  if (*dlsinp.name != '\0') {
    if (dlsinp.f != NULL)
      dlsinp.f = freopen(dlsinp.name, "w", dlsinp.f);
    else
      dlsinp.f = fopen(dlsinp.name, "w");
  } else {
    if (dlsinp.f != NULL)
      rewind(dlsinp.f);
    else
      dlsinp.f = tmpfile();
  }
  if (dlsinp.f == NULL)
    _EscIO(FileNotFound);
  SETUPBUF(dlsinp.f, Char);
  n = 1;
  write_input_file(&dlsinp, n, &strukturdaten, &distdat);
}


Static Void m_LOADATinput(void)
{   /*#### BK loadatinput*/
  if (*loadainp.name != '\0') {
    if (loadainp.f != NULL)
      loadainp.f = freopen(loadainp.name, "w", loadainp.f);
    else
      loadainp.f = fopen(loadainp.name, "w");
  } else {
    if (loadainp.f != NULL)
      rewind(loadainp.f);
    else
      loadainp.f = tmpfile();
  }
  if (loadainp.f == NULL)
    _EscIO(FileNotFound);
  SETUPBUF(loadainp.f, Char);
  if ((berechnet & (1 << ((long)bindungen - (long)struktur))) == 0)
        /*#### BK loadatinput */
          m_berechne_bindungen(0);
  write_input_file(&loadainp, 3L, &strukturdaten, &distdat);
}


Static Void m_berechne_o_positionen(void)
{
  boolean status;

  if ((berechnet & (1 << ((long)bindungen - (long)struktur))) == 0)
    m_berechne_bindungen(0);
  baue_o_atome_ein(&strukturdaten, &status);
  berechnet &= ~(1 << ((long)distanzen - (long)struktur));
}


static int fgetline(FILE *fpin, char s[], int size_s)
{
  int         last_s, c, i;


  last_s = size_s - 1;

  i = 0;

  while ((c = getc(fpin)) != EOF && c != '\n')
  {
#ifdef __MSDOS__
    if (c == 0x1A /* CtrlZ */) {
      ungetc(c, fpin);
      c = EOF;
      break;
    }
#endif

    if (i < last_s) s[i++] = c;
  }

  s[i] = '\0';

  if (i == 0 && c == EOF)
    return 0;
  else
    return 1;
}


static int getnum(char *str, int **num, int max_nnum)
{
  int     i, n, nnum;
  char    buf[128], xtrac, *s;
  double  fn;


  nnum = 0;
  i = 0;

  for (s = str; *s; s++)
  {
    if (*s == ',')
    {
      nnum++;
      i = -1;
    }
    else if (isspace(*s))
    {
      if (i == 1)
      {
        nnum++;
        i = 0;
      }
    }
    else
      i = 1;
  }

  if (i)
    nnum++;

  if (max_nnum && nnum > max_nnum)
    return -1;

  if (nnum)
    *num = Malloc(nnum * sizeof (*(*num)));

  s = str;

  for (i = 0; i < nnum; i++)
  {
    while (*s && isspace(*s)) s++;

    n = 0;

    while (n < sizeof buf / sizeof (*buf))
    {
      if (*s == '\0' || *s == ',' || isspace(*s))
        break;

      buf[n++] = *s++;
    }

    if (n == sizeof buf / sizeof (*buf))
    {
      free(*num);
      *num = NULL;
      return -(i + 1);
    }

    buf[n] = '\0';

    if (n == 0)
      fn = 0.;
    else
    {
          n = sscanf(buf, "%lf%c", &fn, &xtrac);
      if (n != 1)
      {
        free(*num);
        *num = NULL;
        return -(i + 1);
      }

      if (fn < 0.)
          fn -= .5;
      else
          fn += .5;
    }

    (*num)[i] = (int) fn;

    if (*s == ',')
      s++;
  }

  return nnum;
}


Static Void m_write_strukturdaten_karte(int Flag)
{
  if (!strudatnew_offen) {   /*#### BK ew */
    if (*strudatn.name != '\0') {
      if (strudatn.f != NULL)
        strudatn.f = freopen(strudatn.name, "w", strudatn.f);
      else
        strudatn.f = fopen(strudatn.name, "w");
    } else {
      if (strudatn.f != NULL)
        rewind(strudatn.f);
      else
        strudatn.f = tmpfile();
    }
    if (strudatn.f == NULL)
      _EscIO(FileNotFound);
    SETUPBUF(strudatn.f, Char);
    strudatnew_offen = true;
  }
  /*#### BK ew */

  if      (Flag == 1)
    write_P_strukturdaten_karte(&strudatn, &strukturdaten);
  else if (Flag == 2)
    write_PR_strukturdaten_karte(&strudatn, &strukturdaten, NULL);
  else if (Flag == 3)
  {
    int   *Fxyz;
    char  buf[80];

    for (;;)
    {
      fprintf(stdout, " Enter factors for a,b,c:\n");
      fgetline(stdin, buf, sizeof buf);

      if (getnum(buf, &Fxyz, 3) == 3)
      {
        if (Fxyz[0] > 0 && Fxyz[1] > 0 && Fxyz[2] > 0)
          break;
      }

      if (Fxyz) {
        free(Fxyz); Fxyz = NULL;
      }
    }

    write_PR_strukturdaten_karte(&strudatn, &strukturdaten, Fxyz);
  }
  else
    write_strukturdaten_karte(&strudatn, &strukturdaten);
}


Static Void m_loesche_atomlage(void)
{
  kt_atomlage *pala;
  _TEXT TEMP;
  kt_atomart *WITH;

  printf(" Atom position to be deleted:");
  if (P_eof(stdin))
    return;
  TEMP.f = stdin;
  *TEMP.name = '\0';
  read_atomlage_bez(&TEMP, &pala);
  scanf("%*[^\n]");
  getchar();
  if (pala == NULL) {
    printf(" ERROR  cannot delete this atom\n");
    return;
  }
  delete_atomlage(pala, &strukturdaten);
  berechnet &= ~(1 << ((long)distanzen - (long)struktur));
  WITH = pala->atom[0];
  while (WITH->bindung[0].dessen_atomart != NULL)
    loesche_bindung(&pala->grundatom, &WITH->bindung[0]);
}


Static Void m_erzeuge_atomlage(void)
{
  kt_atomlage *neue_atomlage;
  boolean b;
  _TEXT TEMP;

  printf(" Add atom position no. %ld\n", strukturdaten.anzatomlagen + 1);
  printf(" label (element symbol + number)   coordinates   coordination number\n");
  neue_atomlage = (kt_atomlage *)Malloc(sizeof(kt_atomlage));
  printf(" (%2ld)  :", strukturdaten.anzatomlagen + 1);
  TEMP.f = stdin;
  *TEMP.name = '\0';
  read_atomlage(&TEMP, &neue_atomlage);
  frac_vektor(neue_atomlage->atomlagekoord, neue_atomlage->atomlagekoord);
  lade_atomlage(&neue_atomlage, &strukturdaten, &b);
  berechnet &= ~(1 << ((long)distanzen - (long)struktur));
}


Static Void m_erzeuge_bindung(void)
{
  kt_atom a1, a2;
  double d;
  long n;
  _TEXT TEMP;

  printf(" Create new bond:             first atom:");
  TEMP.f = stdin;
  *TEMP.name = '\0';
  readln_atom(&TEMP, &a1);
  if (a1.dessen_atomart != NULL) {
    a2.dessen_atomart = NULL;
    printf(" second atom or distance no. (from disd):");
    while ((!P_eoln(stdin)) & (P_peek(stdin) == ' '))
      getc(stdin);
    if (!P_eoln(stdin)) {
      if (P_peek(stdin) == '(')
        getc(stdin);
      if (P_peek(stdin) >= '1' && P_peek(stdin) <= '9' || P_peek(stdin) == '0') {
        scanf("%ld%*[^\n]", &n);
        getchar();
        berechne_nte_distanz(&a1, n, &d, &a2);
      } else {
        TEMP.f = stdin;
        *TEMP.name = '\0';
        readln_atom(&TEMP, &a2);
      }
    }
    if (a2.dessen_atomart != NULL)
      erzeuge_bindung(&a1, &a2);
    else
      printf(" ERROR  cannot create this bond\n");
  } else
    printf(" ERROR  cannot create this bond\n");
  berechnet |= 1 << ((long)bindungen - (long)struktur);
}


Static Void m_loesche_bindung(void)
{
  kt_atom a1, a2;
  long n;
  _TEXT TEMP;

  printf(" Delete bond:             first atom:");
  TEMP.f = stdin;
  *TEMP.name = '\0';
  readln_atom(&TEMP, &a1);
  if (a1.dessen_atomart == NULL) {
    printf(" ERROR  cannot delete this bond\n");
    return;
  }
  a2.dessen_atomart = NULL;
  printf(" second atom or bond no. (from disb):");
  while ((!P_eoln(stdin)) & (P_peek(stdin) == ' '))
    getc(stdin);
  if (!P_eoln(stdin)) {
    if (P_peek(stdin) == '(')
      getc(stdin);
    if (P_peek(stdin) >= '1' && P_peek(stdin) <= '9' || P_peek(stdin) == '0') {
      scanf("%ld%*[^\n]", &n);
      getchar();
      if (n > 0 && n <= maxanzahlbindungen)
        berechne_nte_bindung(&a1, n, &a2);
      else
        a2.dessen_atomart = NULL;
      if (a2.dessen_atomart == NULL) {
        printf(" ERROR  there is no bond no. %ld for atom \n", n);
        TEMP.f = stdout;
        *TEMP.name = '\0';
        write_atom_name(&TEMP, &a1);
      }
    } else {
      TEMP.f = stdin;
      *TEMP.name = '\0';
      readln_atom(&TEMP, &a2);
    }
  }
  if (a2.dessen_atomart != NULL)
    loesche_bindung(&a1, &a2);
  else
    printf(" ERROR  cannot delete this bond\n");
}


Static Void m_loesche_alle_bindung(void)
{
  loesche_alle_bindungen(&strukturdaten);
}


Static Void m_verknuepfe_t_atome_direkt(void)
{
  if ((berechnet & (1 << ((long)bindungen - (long)struktur))) == 0)
    m_berechne_bindungen(0);
  verknuepfe_t_atome_direkt();
}


Static Void m_berechne_koord_seq(_TEXT *Tfpout)
{
  long i;
  tt_koord_sequenz z;

  if ((berechnet & (1 << ((long)bindungen - (long)struktur))) == 0)
    m_berechne_bindungen(0);
  fprintf(Tfpout->f, "\n Coordination sequences of the T-atoms:\n");
  for (i = 0; i < strukturdaten.anzatomlagen; i++) {
    berechne_koord_seq(&strukturdaten.atomlage[i]->grundatom, 11L, z);
    fprintf(Tfpout->f, " %.6s", strukturdaten.atomlage[i]->aname);
    write_koord_seq(Tfpout, z);
  }
}


Static Void m_bestimme_topologie(_TEXT *Tfpout)
{
  if ((berechnet & (1 << ((long)bindungen - (long)struktur))) == 0)
    m_berechne_bindungen(0);
  bestimme_topologie(Tfpout);
}


Static Void m_berechne_loop_konfig(_TEXT *Tfpout)
{
  long i, FORLIM;

  if ((berechnet & (1 << ((long)bindungen - (long)struktur))) == 0)
    m_berechne_bindungen(0);
  FORLIM = strukturdaten.anzatomlagen;
  for (i = 0; i < FORLIM; i++) {
    fprintf(Tfpout->f, "\n Loop configuration for ");
    write_atom_name(Tfpout, &strukturdaten.atomlage[i]->grundatom);
    putc('\n', Tfpout->f);
    berechne_loop_konfig(Tfpout, &strukturdaten.atomlage[i]->grundatom);
  }
}


Static Void m_write_entry_liste(void)
{
  long i;

  if (*strudat.name != '\0') {
    if (strudat.f != NULL)
      strudat.f = freopen(strudat.name, "r", strudat.f);
    else
      strudat.f = fopen(strudat.name, "r");
  } else
    rewind(strudat.f);
  if (strudat.f == NULL)
    _EscIO(FileNotFound);
  RESETBUF(strudat.f, Char);
  printf(" Entries on file STRUDAT\n");
  while (!BUFEOF(strudat.f)) {
    if (P_peek(strudat.f) == '*') {
      i = 0;
      putchar(' ');
      while (!P_eoln(strudat.f)) {
        putchar(P_peek(strudat.f));
        getc(strudat.f);
        i++;
      }
      putchar('\n');
      if (i > 15)
        printf(" ERROR  compound identification code too long\n");
    }
    fscanf(strudat.f, "%*[^\n]");
    getc(strudat.f);
  }
}


Static Void m_read_befehl(Char *befehl, int argc, char *argv[])
{
  int  i, c, mode;

  static int  iarg = 1;


  printf(" ==>");

  mode = 0;

  if (iarg < argc)
  {
    printf("%s\n", argv[iarg]);

    for (i = 0; argv[iarg][i]; i++)
    {
      if (i == 8)
      {
        i = 0;
        mode = -1;
        break;
      }

      befehl[i] = tolower(argv[iarg][i]);
    }

    iarg++;
  }
  else
  {
    i = 0;

    while ((c = getc(stdin)) != EOF && c != '\n')
    {
      if (mode == 0)
      {
        if (isspace(c))
        {
          if (i)
            mode = 1;
        }
        else if (i == 8)
        {
          i = 0;
          mode = -1;
        }
        else
          befehl[i++] = tolower(c);
      }
    }
  }

  while (i < 8) befehl[i++] = ' ';

  if (mode == -1)
    printf(" ERROR  command too long\n");
}


int main(int argc, char *argv[])
{  /* main */
  _TEXT Tstdout;

  PASCAL_MAIN(argc, argv);
  coseq.f = NULL;
  strcpy(coseq.name, "coseq");
  list.f = NULL;
  strcpy(list.name, "list");
  cif.f = NULL;                       // Stefs 05-07-2012
  strcpy(cif.name, "structure.cif"); 
  loadainp.f = NULL;
  strcpy(loadainp.name, "loadainp");
  dlsinp.f = NULL;
  strcpy(dlsinp.name, "dlsinp");
  distdat.f = NULL;
  strcpy(distdat.name, "distdat");
  strudatn.f = NULL;
  strcpy(strudatn.name, "strudatn");
  symdat.f = NULL;
  strcpy(symdat.name, "symdat");
  strudat.f = NULL;
  strcpy(strudat.name, "strudat");
  m_initialisiere(&strukturdaten);

   Tstdout.f    = stdout;
  *Tstdout.name = '\0';

  do {
    m_read_befehl(befehl, argc, argv);
    if (strncmp(befehl, "        ", sizeof(mt_befehl))) {
      if      (!strncmp(befehl, "calcs   ", sizeof(mt_befehl)))
        m_berechne_koord_seq(&Tstdout);
      else if (!strncmp(befehl, "wricsq  ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.csq");
        m_initialize(&list);
        m_berechne_koord_seq(&list);
      }
      else if (!strncmp(befehl, "cald    ", sizeof(mt_befehl)))
        m_berechne_distanzen((long)maxanzahldistanzen);
      else if (!strncmp(befehl, "calda   ", sizeof(mt_befehl)))
        berechne_distanzen_algebraisch(&rginfo);
      else if (!strncmp(befehl, "creab   ", sizeof(mt_befehl)))
        m_berechne_bindungen(0);
      else if (!strncmp(befehl, "creabd  ", sizeof(mt_befehl)))
        m_berechne_bindungen(1);
      else if (!strncmp(befehl, "callc   ", sizeof(mt_befehl)))
        m_berechne_loop_konfig(&Tstdout);
      else if (!strncmp(befehl, "wrilc   ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.lc");
        m_initialize(&list);
        m_berechne_loop_konfig(&list);
      }
      else if (!strncmp(befehl, "addo    ", sizeof(mt_befehl)))
        m_berechne_o_positionen();
      else if (!strncmp(befehl, "adda    ", sizeof(mt_befehl)))
        m_erzeuge_atomlage();
      else if (!strncmp(befehl, "creb    ", sizeof(mt_befehl)))
        m_erzeuge_bindung();
      else if (!strncmp(befehl, "conb    ", sizeof(mt_befehl)))
        m_write_bindungs_matrix(&Tstdout);
      else if (!strncmp(befehl, "cons    ", sizeof(mt_befehl)))
        m_contr_st_raumgruppe();
      else if (!strncmp(befehl, "dela    ", sizeof(mt_befehl)))
        m_loesche_atomlage();
      else if (!strncmp(befehl, "delab   ", sizeof(mt_befehl)))
        m_loesche_alle_bindung();
      else if (!strncmp(befehl, "delb    ", sizeof(mt_befehl)))
        m_loesche_bindung();
      else if (!strncmp(befehl, "delo    ", sizeof(mt_befehl)))
        m_verknuepfe_t_atome_direkt();
      else if (!strncmp(befehl, "dett    ", sizeof(mt_befehl)))
        m_bestimme_topologie(&Tstdout);
      else if (!strncmp(befehl, "writ    ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.t");
        m_initialize(&list);  
        m_bestimme_topologie(&list);
      }
      else if (!strncmp(befehl, "disab   ", sizeof(mt_befehl)))
        m_write_alle_bindungen(&Tstdout);
      else if (!strncmp(befehl, "disb    ", sizeof(mt_befehl)))
        m_write_bindungen(&Tstdout);
      else if (!strncmp(befehl, "discs   ", sizeof(mt_befehl)))
        m_write_strukturdaten(&Tstdout, &strukturdaten);

      else if (!strncmp(befehl, "discif  ", sizeof(mt_befehl))) { // Stefs 05-07-2012
        putc('\n', Tstdout.f);
        m_write_cif(&Tstdout, &strukturdaten, &rginfo);
      }
      else if (!strncmp(befehl, "discssr ", sizeof(mt_befehl))) {
        putc('\n', Tstdout.f);
        m_write_cssr(&Tstdout, &strukturdaten);
      }
      else if (!strncmp(befehl, "disxtl  ", sizeof(mt_befehl))) {
        putc('\n', Tstdout.f);
        m_write_xtl(&Tstdout, &strukturdaten);
      }
      else if (!strncmp(befehl, "dismvm  ", sizeof(mt_befehl))) {
        putc('\n', Tstdout.f);
        m_write_mvm(&Tstdout, &strukturdaten, 0);
      }
      else if (!strncmp(befehl, "dismvm# ", sizeof(mt_befehl))) {
        putc('\n', Tstdout.f);
        m_write_mvm(&Tstdout, &strukturdaten, 1);
      }
      else if (!strncmp(befehl, "disfocus", sizeof(mt_befehl))) {
        putc('\n', Tstdout.f);
        m_write_focus(&Tstdout, &strukturdaten);
      }
      else if (!strncmp(befehl, "disgcs  ", sizeof(mt_befehl))) {
        putc('\n', Tstdout.f);
        m_write_GSAS_crystal_structure(&Tstdout, &strukturdaten, 1);
      }
      else if (!strncmp(befehl, "disglai ", sizeof(mt_befehl))) {
        putc('\n', Tstdout.f);
        m_write_GSAS_L_A_I(&Tstdout, &strukturdaten);
      }
      else if (!strncmp(befehl, "disgsc  ", sizeof(mt_befehl))) {
        putc('\n', Tstdout.f);
        m_write_GSAS_soft_constraints(&Tstdout, 1, 2);
      }

      else if (!strncmp(befehl, "disad   ", sizeof(mt_befehl)))
        m_write_alle_distanzen(&Tstdout, 0);
      else if (!strncmp(befehl, "disadu  ", sizeof(mt_befehl)))
        m_write_alle_distanzen(&Tstdout, 1);

      else if (!strncmp(befehl, "disd    ", sizeof(mt_befehl)))
        m_write_distanzen(&Tstdout);
      else if (!strncmp(befehl, "disda   ", sizeof(mt_befehl)))
        m_write_distanzen_und_winkel(&Tstdout, 1);
      else if (!strncmp(befehl, "disdau  ", sizeof(mt_befehl)))
        m_write_distanzen_und_winkel(&Tstdout, 0);
      else if (!strncmp(befehl, "disdv   ", sizeof(mt_befehl)))
        m_write_distanzen_und_vektoren(&Tstdout);
      else if (!strncmp(befehl, "disc    ", sizeof(mt_befehl)))
        m_write_koordinaten(&Tstdout);
      else if (!strncmp(befehl, "dise    ", sizeof(mt_befehl)))
        m_write_entry_liste();
      else if (!strncmp(befehl, "diss    ", sizeof(mt_befehl)))
        write_st_raumgruppe(&Tstdout, &rginfo);
      else if (!strncmp(befehl, "help    ", sizeof(mt_befehl)))
        m_write_hilfe();
      else if (!strncmp(befehl, "quit    ", sizeof(mt_befehl)))
        putchar('\n');
      else if (!strncmp(befehl, "exit    ", sizeof(mt_befehl)))
        putchar('\n');
      else if (!strncmp(befehl, "reacs   ", sizeof(mt_befehl)))
        m_read_strukturdaten();
      else if (!strncmp(befehl, "reafcs  ", sizeof(mt_befehl)))
        m_read_erste_strukturdaten();
      else if (!strncmp(befehl, "reancs  ", sizeof(mt_befehl)))
        m_read_naechste_strukturdaten();
      else if (!strncmp(befehl, "wrie    ", sizeof(mt_befehl)))
        m_write_strukturdaten_karte(0);
      else if (!strncmp(befehl, "wrics   ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.cs");
        m_initialize(&list);
        m_write_strukturdaten(&list, &strukturdaten);
      }
      else if (!strncmp(befehl, "wricif  ", sizeof(mt_befehl))) // Stefs 05-07-2012
      {
        m_initialize(&cif);
        m_write_cif(&cif, &strukturdaten, &rginfo);
        fclose(cif.f);
        cif.f = NULL;
      }
      else if (!strncmp(befehl, "wricssr ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.cssr");
        m_initialize(&list);
        m_write_cssr(&list, &strukturdaten);
      }
      else if (!strncmp(befehl, "wrixtl  ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.xtl");
        m_initialize(&list);
        m_write_xtl(&list, &strukturdaten);
      }
      else if (!strncmp(befehl, "wrimvm  ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.mvm");
        m_initialize(&list);
        m_write_mvm(&list, &strukturdaten, 0);
      }
      else if (!strncmp(befehl, "wrimvm# ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.mvm");
        m_initialize(&list);
        m_write_mvm(&list, &strukturdaten, 1);
      }
      else if (!strncmp(befehl, "wrifocus", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.foc");
        m_initialize(&list);
        m_write_focus(&list, &strukturdaten);
      }
      else if (!strncmp(befehl, "wrigcs  ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.gcs");
        m_initialize(&list);
        m_write_GSAS_crystal_structure(&list, &strukturdaten, 0);
      }
      else if (!strncmp(befehl, "wriglai ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.glai");
        m_initialize(&list);
        m_write_GSAS_L_A_I(&list, &strukturdaten);
      }
      else if (!strncmp(befehl, "wrigsc  ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.gsc");
        m_initialize(&list);
        m_write_GSAS_soft_constraints(&list, 0, 0);
      }
      else if (!strncmp(befehl, "wriep   ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.ep");
        m_initialize(&list);
        m_write_strukturdaten_karte(1);
      }
      else if (!strncmp(befehl, "wriepn  ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.epn");
        m_initialize(&list);
        m_write_strukturdaten_karte(2);
      }
      else if (!strncmp(befehl, "wriepr  ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.epr");
        m_initialize(&list);
        m_write_strukturdaten_karte(3);
      }
      else if (!strncmp(befehl, "wrida   ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.da");
        m_initialize(&list);
        m_write_distanzen_und_winkel(&list, 1);
      }
      else if (!strncmp(befehl, "wridau  ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.dau");
        m_initialize(&list);
        m_write_distanzen_und_winkel(&list, 0);
      }
      else if (!strncmp(befehl, "wridv   ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.dv");
        m_initialize(&list);
        m_write_distanzen_und_vektoren(&list);
      }
      else if (!strncmp(befehl, "wris    ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.s");
        m_initialize(&list);
        write_st_raumgruppe(&list, &rginfo);
      }
      else if (!strncmp(befehl, "wric    ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.c");
        m_initialize(&list);
        m_write_koordinaten(&list);
      }
      else if (!strncmp(befehl, "wriad   ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.ad");
        m_initialize(&list);
        m_write_alle_distanzen(&list, 0);
      }
      else if (!strncmp(befehl, "wriadu  ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.adu");
        m_initialize(&list);
        m_write_alle_distanzen(&list, 1);
      }
      else if (!strncmp(befehl, "wrid    ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.d");
        m_initialize(&list);
        m_write_distanzen(&list);
      }
      else if (!strncmp(befehl, "wriab   ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.ab");
        m_initialize(&list);
        m_write_alle_bindungen(&list);
      }
      else if (!strncmp(befehl, "wrib    ", sizeof(mt_befehl)))
      {
        strcpy(list.name, "list.b");
        m_initialize(&list);
        m_write_bindungen(&list);
      }
      else if (!strncmp(befehl, "wriid   ", sizeof(mt_befehl)))
      {
        m_DLSinput();
        fclose(dlsinp.f);
      }      
      else if (!strncmp(befehl, "writc   ", sizeof(mt_befehl)))
      {
        m_DLSinput_tetcon();
        fclose(dlsinp.f);
      }
      else if (!strncmp(befehl, "wriil   ", sizeof(mt_befehl)))
      {  
        m_LOADATinput();
        fclose(loadainp.f);
      }
      else
        printf(" ERROR  unknown command %.8s\n", befehl);




    }
    
    // close list file

    if (list.f != NULL)
      fclose(list.f);
    
    putchar('\n');
  } while (strncmp(befehl, "quit    ", sizeof(mt_befehl)) &&
           strncmp(befehl, "exit    ", sizeof(mt_befehl)));
  printf(
    " ******************************************************************************\n");
  printf(
    " *   end of program  KRIBER                                                   *\n");
  printf(
    " ******************************************************************************\n");
  if (strudat.f != NULL)
    fclose(strudat.f);
  if (symdat.f != NULL)
    fclose(symdat.f);
  if (strudatn.f != NULL)
    fclose(strudatn.f);
  if (distdat.f != NULL)
    fclose(distdat.f);
  if (dlsinp.f != NULL)
    fclose(dlsinp.f);
  if (loadainp.f != NULL)
    fclose(loadainp.f);
  if (list.f != NULL)
    fclose(list.f);
  if (cif.f != NULL);   // Stefs 05-07-2012
    fclose(cif.f);
  if (coseq.f != NULL)
    fclose(coseq.f);
  exit(EXIT_SUCCESS);
  return 0;
}



/* End. */
