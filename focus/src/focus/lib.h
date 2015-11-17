#ifndef A_LIB_H__
#define A_LIB_H__

#ifndef Ex
#ifdef A_AppGlobal__
#define Ex
#else
#define Ex extern
#endif
#endif

#ifndef SEEK_SET
#define SEEK_SET 0
#endif

#define LNLINE 4096
#define FNLEN 255

#define LNIline LNLINE

Ex long       Ilinec;
Ex int        Imode;
Ex int        Itabsize;
Ex char       Iline[LNIline + 1];

#define I_raw          0x00
#define I_noemptyline  0x01
#define I_uppercase    0x02
#define I_lowercase    0x04
#define I_noleadblks   0x08
#define I_notrailblks  0x10
#define I_notabs       0x20

#define CtrlZ 0x001A

int str_icmp(const char *s, const char *t);
char *lower(char *s);
char *upper(char *s);
char firstnonblank(const char *s);
int noblks(char *s);
int lnblnk(const char *s);
void noleadblks(char *s);
void notrailblks(char *s);
int f_get_line(FILE *fp, char *s, int size_s);
int fgetIline(FILE *fpin);
int str_begin(const char *s1, const char *s2); /* string begin */
int str_ibegin(const char *s1, const char *s2); /* string ignore-case begin */
char *str_xcpy(char *dest, const char *src, size_t max); /* string max copy */
void fillblks(char *s, int len);
int str_arg(char *arg, const char *s, int n);
long flinec(FILE *fp);
void Sum2bmc(int Nxy,
             double minx, double maxx,
             double miny, double maxy,
             double Sx, double Sx2,
             double Sy, double Sy2, double Sxy,
             double *b, double *m, double *c);
int getnum(char *str, int **num, int max_nnum);


#define fgetline(fpin, s) f_get_line((fpin), (s), sizeof (s) / sizeof (*(s)))


#undef Ex

#endif /* A_LIB_H__ */
