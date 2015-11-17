#ifndef A_LATTICE_H__
#define A_LATTICE_H__


typedef struct {
                 double   a, b, c;
                 double   alpha, beta, gamma;
                 double   sa, sb, sg;
                 double   ca, cb, cg;
                 double   v;
               }
               T_LatticeConstants;


#define DotC(u, v) ((u)[0] * (v)[0] + (u)[1] * (v)[1] + (u)[2] * (v)[2])

#define CrossC(r, s, rxs)\
  {\
    (rxs)[0] = (r)[1] * (s)[2] - (s)[1] * (r)[2];\
    (rxs)[1] = (r)[2] * (s)[0] - (s)[2] * (r)[0];\
    (rxs)[2] = (r)[0] * (s)[1] - (s)[0] * (r)[1];\
  }


double Q_hkl(int h, int k, int l, T_LatticeConstants *rlc);
int Lc2RLc(T_LatticeConstants *lc, T_LatticeConstants *rlc);
int SetupTrCrystCarte(T_LatticeConstants *lc, T_LatticeConstants *rlc,
                      int Convention,
                      Fprec *TM, int FlagInverse);
void Lc2MetricalMx(T_LatticeConstants *lc, Fprec *G);
Fprec DotG(Fprec *u, Fprec *G, Fprec *v);
void CrossG(Fprec sqrtdetG, Fprec *G, Fprec *r, Fprec *s, Fprec *rxs);
void CalcMaxhkl(T_LatticeConstants *RLc, Fprec *dmin, Fprec Lambda,
                int *h, int *k, int *l);


#endif /* A_LATTICE_H__ */
