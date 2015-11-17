#ifndef A_CGITBX_H__
#define A_CGITBX_H__


#define RQM_NONE  0
#define RQM_GET   1
#define RQM_POST  2


int putcURL(int c, FILE *fpout);
int fputsURL(const char *s, FILE *fpout);
char *EscapeURL(const char *URL);
char *UnescapeURL(char *eURL);
char *PlusToSpace(char *s);
char *fgetword(FILE *fpin, char sep, int *cl);
char *sgetword(const char *s, char sep, int iword);
char *SeparateWords(char *s, int sep);
void chr2html(FILE *fpout, int c);
void str2html(FILE *fpout, const char *s);
void htmlHead(FILE *fpout,
              const char *base,
              const char *title,
              const char *body_props);
void htmlTail(FILE *fpout);


#endif /* A_CGITBX_H__ */
