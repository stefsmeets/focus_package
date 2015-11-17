#ifndef _ATOMINFODECL
#define _ATOMINFODECL

typedef struct
  { 
    char *atom[10];
    int atom_ratio[10];
  }
  T_ChemicalComp;


typedef struct
  {
    const char   *Label;
    double        a[4], b[4], c;
  }
  T_SF_IT92_CAA;


typedef struct
  {
    const char   *Label;
    double        a[5], b[5], c;
  }
  T_SF_WK95_CAA;

typedef struct
  {
    int         Z;
    const char  *Symbol;
    const char  *Name;
    double       Mass;
  }
  T_PSE;

typedef struct
  {
    const char  *Label;
    double       Radius;
  }
  T_AtomRadius;

#endif

