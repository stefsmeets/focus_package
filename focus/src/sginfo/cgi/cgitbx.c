#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>


#include "cgitbx.h"


static const char HexCodes[] = "0123456789ABCDEF";


typedef struct
  {
    int         c;
    const char  *m;
  }
  T_EscapeMap;


/*
  http://www.netspace.org/users/dwb/url-guide.html
 */

static const T_EscapeMap URL_EscapeMap[] =
  {
    /* Unsafe */
    { ' ',  "20" },
    { '<',  "3C" },
    { '>',  "3E" },
    { '#',  "23" },
    { '%',  "25" },
    { '{',  "7B" },
    { '}',  "7D" },
    { '|',  "7C" },
    { '\\', "74" },
    { '^',  "5E" },
    { '~',  "7E" },
    { '[',  "5B" },
    { ']',  "5D" },
    { '`',  "60" },
    /* Reserved */
    { ';',  "3B" },
    { '/',  "2F" },
    { '?',  "3F" },
    { ':',  "3A" },
    { '@',  "40" },
    { '=',  "3D" },
    { '&',  "26" },
    /* Added by rwgk */
    { '"',  "22" },
    { '\0', NULL }
  };


static const T_EscapeMap *FindURL_Escape(int c)
{
  const T_EscapeMap  *uem;


  for (uem = URL_EscapeMap; uem->m; uem++)
    if (c == uem->c) return uem;

  return NULL;
}


static int putcURLnOut = 0;


int putcURL(int c, FILE *fpout)
{
  const T_EscapeMap  *uem;


  if (c == ' ')
  {
    c = putc('+', fpout);
    putcURLnOut = 1;
  }
  else
  {
        uem = FindURL_Escape(c);
    if (uem)
    {
                    c = putc('%', fpout);
      if (c != EOF) c = putc(uem->m[0], fpout);
      if (c != EOF) c = putc(uem->m[1], fpout);

      putcURLnOut = 3;
    }
    else
    {
      c = putc(c, fpout);
      putcURLnOut = 1;
    }
  }

  return c;
}


int fputsURL(const char *s, FILE *fpout)
{
  int  n;


  for (n = 0; *s; s++, n += putcURLnOut)
    if (putcURL(*s, fpout) == EOF)
      return EOF;

  return n;
}


char *EscapeURL(const char *URL)
{
  int                l;
  char               *eURL, *e;
  const char         *u;
  const T_EscapeMap  *uem;


  l = 0;

  for (u = URL; *u; u++)
  {
        uem = FindURL_Escape(*u);
    if (uem)
      l += 3;
    else
      l++;
  }

  e = eURL = malloc((l + 1) * sizeof (*eURL));
  if (eURL == NULL)
    return NULL;

  for (u = URL; *u; u++)
  {
        uem = FindURL_Escape(*u);
    if (uem)
    {
      *e++ = '%';
      *e++ = uem->m[0];
      *e++ = uem->m[1];
    }
    else
      *e++ = *u;
  }

  *e = '\0';

  return eURL;
}


static int e2c(char *e)
{
  int  i, j, e0, e1;


  if (e[0] == '\0' || e[1] == '\0')
    return '\0';

  e0 = toupper(e[0]);
  e1 = toupper(e[1]);

  for (i = 0; HexCodes[i]; i++)
    if (HexCodes[i] == e0) break;

  if (HexCodes[i] == '\0')
    return ' ';

  for (j = 0; HexCodes[j]; j++)
    if (HexCodes[j] == e1) break;

  if (HexCodes[j] == '\0')
    return ' ';

  return (i << 4) | j;
}


char *UnescapeURL(char *eURL)
{
  char  *e, *URL;


  e = URL = eURL;

  while (*e)
  {
    if (*e == '%')
    {
      e++;
      *URL++ = e2c(e);
      if (*e++ == '\0') break;
      if (*e++ == '\0') break;
    }
    else
      *URL++ = *e++;
  }

  *URL = '\0';

  return eURL;
}


char *PlusToSpace(char *s)
{
  char  *s0;


  s0 = s;

  for (; *s; s++)
    if (*s == '+') *s = ' ';

  return s0;
}


char *fgetword(FILE *fpin, char sep, int *cl)
{
  int   i, c;
  int   mword;
  char  *word;


      mword = 511;
  if (mword > *cl) mword = *cl;
  if (mword <   0) mword = 0;

      word = malloc((mword + 1) * sizeof (*word));
  if (word == NULL)
    return NULL;

  i = 0;

  for (;;)
  {
    if (*cl <= 0) break;

        c = fgetc(fpin); (*cl)--;
    if (c == EOF) break;
    if (c == sep) break;

    if (i == mword)
    {                           mword += 512;
          word = realloc(word, (mword + 1) * sizeof (*word));
      if (word == NULL)
        return NULL;
    }

    word[i++] = c;
  }

  word[i] = '\0';

  if (i < mword)
    word = realloc(word, (i + 1) * sizeof (*word));

  return word;
}


char *sgetword(const char *s, char sep, int iword)
{
  int         l;
  char        *bw, *w;
  const char  *sw;


  while (*s && iword)
  {
    while (*s && *s != sep) s++;
    if (*s) s++;
    iword--;
  }

  sw = s;

  l = 0;

  while (*s && *s != sep) {
    s++;
    l++;
  }

  s = sw;

  w = bw = malloc((l + 1) * sizeof (*w));

  while (*s && *s != sep)
    *w++ = *s++;

  *w = '\0';

  return bw;
}


char *SeparateWords(char *s, int sep)
{
  for (; *s; s++)
  {
    if (*s == sep)
    {
      *s = '\0';
      return s + 1;
    }
  }

  return s;
}


void chr2html(FILE *fpout, int c)
{
  switch (c)
  {
    case '&': fprintf(fpout, "&amp;");  break;
    case '<': fprintf(fpout, "&lt;");   break;
    case '>': fprintf(fpout, "&gt;");   break;
    case '"': fprintf(fpout, "&quot;"); break;
    default:
      putc(c, fpout); break;
  }
}


void str2html(FILE *fpout, const char *s)
{
  while (*s) chr2html(fpout, *s++);
}


void htmlHead(FILE *fpout,
              const char *base,
              const char *title,
              const char *body_props)
{
  fprintf(fpout, "Content-type: text/html\n\n");

  fprintf(fpout, "<html>\n");
  fprintf(fpout, "<head>\n");

  if (base) fprintf(fpout, "<base href=\"%s\">\n", base);

  fprintf(fpout, "<title>%s</title>\n", title);
  fprintf(fpout, "</head>\n");

  fprintf(fpout, "<body");
  if (body_props) fprintf(fpout, " %s", body_props);
  fprintf(fpout, ">\n");
}


void htmlTail(FILE *fpout)
{
  fprintf(fpout, "</body>\n");
  fprintf(fpout, "</html>\n");
}
