#ifndef A_MATRIX_H__
#define A_MATRIX_H__


int MxInverse(Fprec *a, int n, Fprec *b, int m);
int MxLU_Decomposition(Fprec *a, int n, int *perm);
void MxLU_BackSubstitution(Fprec *a, int n, int *perm, Fprec *b);
void MxMultiply(Fprec *ab, Fprec *a, Fprec *b, int ma, int na, int nb);
void MxTranspose(Fprec *a, int ma, int na, Fprec *t);
void MxDump(Fprec *a, int ma, int na, char *note);
Fprec deter33fMx(const Fprec *Mx);
void Inverse33fMx(const Fprec *Mx, Fprec *InvMx);


#endif /* A_MATRIX_H__ */
