#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#if defined(__sgi)
#include <sys/uio.h>
#include <sys/time.h>
#endif
#include <sys/types.h>
#include <unistd.h>

#if defined(__sgi)
#include <sys/timers.h>
#include <sys/siginfo.h>
#endif
#include <signal.h>


#include "cgitbx.h"


static const char *cmd_name =                        "sginfo";
static const char *SGINFO_CMD_PATH  = "/usr/local/bin/sginfo";

static const char *SGINFO_MAIN_URL  = "/sginfo";
static const char *SGINFO_GIF_URL   = "/sginfo/sginfo.gif";
static const char *SGINFO_QUERY_URL = "/user-cgi-bin/sginfo-query";


#define SetFree(ptr) if (ptr) free(ptr); ptr = NULL

#define fs(s) fputs(s, stdout)


static char *str_dup(const char *s)
{
  char *c = (char *) malloc((strlen(s) + 1) * sizeof (*s));
  if (c != NULL) (void) strcpy(c, s);
  return c;
}


static int str_icmp(const char *s, const char *t)
{
  char  cs, ct;


  while (*s || *t)
  {
    cs = toupper(*s++);
    ct = toupper(*t++);
    if (cs < ct) return -1;
    if (cs > ct) return  1;
  }

  return 0;
}


static char firstnonblank(const char *s)
{
  while (*s && isspace(*s)) s++;
  return *s;
}


static void htmlAddress(FILE *fpout)
{
  fprintf(fpout, "<hr>\n");
  fprintf(fpout, "<address>\n");
  fprintf(fpout, "<a href=\"http://atb.csb.yale.edu/~rwgk/\"\n");
  fprintf(fpout, ">Ralf W. Grosse-Kunstleve</a>\n");
  fprintf(fpout, "&lt;<a href=\"mailto:rwgk@laplace.csb.yale.edu\"\n");
  fprintf(fpout, "                   >rwgk@laplace.csb.yale.edu</a>&gt;\n");
  fprintf(fpout, "</address>\n");
}


static void SendError(const char *msg)
{
  fs("<h3>");
  fs(msg);
  fs("</h3>\n");

  htmlAddress(stdout);
  htmlTail(stdout);

  fclose(stdin);
  fclose(stdout);
  fclose(stderr);

  exit(0);
}


#define SendInternalError(msg) S_endInternalError(msg, __FILE__, __LINE__)

static void S_endInternalError(const char *msg, const char *srcf, int srcl)
{
  fs("<h3>Internal Error:\n");
  str2html(stdout, msg);
  fs("</h3>\n");

  fprintf(stdout, "Source Code File &quot;%s&quot; at line %d\n",
    srcf, srcl);

  htmlAddress(stdout);
  htmlTail(stdout);

  fclose(stdin);
  fclose(stdout);
  fclose(stderr);

  exit(0);
}


static const char *copyarg(const char *s, char **arg)
{
  int         n, b, q;
  char        *a;
  const char  *p;


  if (arg) *arg = NULL;

  a = NULL;

  for (;;)
  {
    p = s;
    n = 0;
    b = 0;
    q = '\0';

    for (;;)
    {
      if (b)
      {
        if (   isspace(*p) == 0
            && *p != '"'
            && *p != '\''
            && *p != '\\')
        {
          if (a) *a++ = '\\';
          n++;
        }

        if (a) *a++ = *p;
        n++;

        b = 0;
      }
      else if (*p == '\\')
        b = 1;

      else if (q)
      {
        if (*p == q)
          q = '\0';
        else
        {
          if (a) *a++ = *p;
          n++;
        }
      }
      else if (   *p == '"'
               || *p == '\'')
        q = *p;

      else if (isspace(*p) == 0)
      {
        if (a) *a++ = *p;
        n++;
      }
      else
        break;

      p++;

      if (*p == '\0')
      {
        if (q)
          SendError("Unmatched quote.");

        if (b)
        {
          if (a) *a++ = '\\';
          n++;
        }

        break;
      }
    }

    if (a) *a = '\0';
    n++;

    if (arg == NULL || *arg) break;

        (*arg) = malloc(n * sizeof (**arg));
    if ((*arg) == NULL)
      SendError("Server Error: Not enough core");

    a = (*arg);
  }

  return p;
}


static char **mkargv(const char *s)
{
  int         argc;
  char        **argv;
  const char  *p;


  argv = NULL;

  for (;;)
  {
    argc = 1;

    p = s;

    for (;;)
    {
      while (*p && isspace(*p)) p++;

      if (*p == '\0') break;

      if (argv == NULL)
        p = copyarg(p, NULL);
      else
        p = copyarg(p, &argv[argc]);

      argc++;
    }

    if (argv) break;

        argv = malloc((argc + 1) * sizeof (*argv));
    if (argv == NULL)
      SendError("Server Error: Not enough core");

    argv[0] = (char *) cmd_name;
  }

  argv[argc] = NULL;

  if (argc > 1 && strcmp(cmd_name, argv[1]) == 0)
  {
    int  i;

    free(argv[1]);

    for (i = 1; i < argc; i++)
      argv[i] = argv[i + 1];

    argc--;

        argv = realloc(argv, (argc + 1) * sizeof (*argv));
    if (argv == NULL)
      SendError("Server Error: Not enough core");
  }

  return argv;
}


static char **addargv(char **argv, const char *s)
{
  int  i;


  for (i = 1; argv[i]; i++)
    if (str_icmp(argv[i], s) == 0)
      return argv;

      argv = realloc(argv, (i + 2) * sizeof (*argv));
  if (argv == NULL)
    SendError("Server Error: Not enough core");

      argv[i] = str_dup(s);
  if (argv[i] == NULL)
    SendError("Server Error: Not enough core");

  i++;

  argv[i] = NULL;

  return argv;
}


static char **cpargv(int argc, char *argv[])
{
  int   i;
  char  **cargv;


      cargv = malloc((argc + 1) * sizeof (*cargv));
  if (cargv == NULL)
    SendError("Server Error: Not enough core");

  cargv[0] = (char *) cmd_name;

  for (i = 1; i < argc; i++)
    cargv[i] = argv[i];

  cargv[i] = NULL;

  return cargv;
}


static void CatchSIGPIPE(int sig)
{
  (void) signal(SIGPIPE, SIG_IGN);

  if (sig != SIGPIPE)
    SendInternalError(NULL);

  SendError("Broken Pipe");
}


static void RunSginfo(char *const *argv, const char *symxyz)
{
  int    lastc, c, i;
  int    pio[2], fdp2cr, fdp2cw, fdc2pr, fdc2pw;
  pid_t  pid;


  if (pipe(pio) != 0)
    SendInternalError("pipe() failed");

  fdp2cr = pio[0];
  fdp2cw = pio[1];

  if (pipe(pio) != 0)
    SendInternalError("pipe() failed");

  fdc2pr = pio[0];
  fdc2pw = pio[1];

  fflush(stdout);
  fflush(stderr);

  pid = fork();

  if (pid < 0)
    SendInternalError("fork() failed");

  if (pid == 0)
  {
    close(fdp2cw);
    close(fdc2pr);

    close(0); i = dup(fdp2cr); close(fdp2cr);

    if (i != 0)
      SendInternalError("dup() failed");

    close(1); i = dup(fdc2pw); close(fdc2pw);

    if (i != 1)
      SendInternalError("dup() failed");

    close(2); i = dup(1);

    if (i != 2)
      SendInternalError("dup() failed");

    (void) execv(SGINFO_CMD_PATH, argv);

    SendInternalError("execv() failed");
  }
  else
  {
    FILE  *fpp2cw;
    FILE  *fpc2pr;

    (void) signal(SIGPIPE, CatchSIGPIPE);

    close(fdp2cr);
    close(fdc2pw);

        fpp2cw = fdopen(fdp2cw, "w");
    if (fpp2cw == NULL)
      SendInternalError("fdopen() failed");

        fpc2pr = fdopen(fdc2pr, "r");
    if (fpc2pr == NULL)
      SendInternalError("fdopen() failed");

    if (symxyz)
    {
      lastc = '\n';

      for (c = *symxyz++; c; c = *symxyz++)
      {
        if (c == '\r' && (*symxyz == '\n' || *symxyz == '\0'))
          continue;

        putc(c, fpp2cw);
        lastc = c;
      }

      if (lastc != '\n')
        putc('\n', fpp2cw);
    }

    fclose(fpp2cw);

    fflush(stdout);

    fs("<pre>\n");

    while ((c = fgetc(fpc2pr)) != EOF)
      chr2html(stdout, c);

    fs("</pre>\n");

    fclose(fpc2pr);
  }
}


static void GetCGI_Variables(int ReqMethod,
                             const char *QueryString, int ContentLength,
                             char **argl, char **symxyz)
{
  int   iword;
  char  *FormItem, *ItemValue;


  FormItem = NULL;

  iword = 0;

  for (;;)
  {
    SetFree(FormItem);

    if     (ReqMethod == RQM_GET)
    {
          FormItem = sgetword(QueryString, '&', iword);
      if (FormItem == NULL)
        SendError("Server Error: Not enough core");

      if (FormItem[0] == '\0')
        break;
    }
    else if (ReqMethod == RQM_POST)
    {
      if (ContentLength == 0)
        break;

          FormItem = fgetword(stdin, '&', &ContentLength);
      if (FormItem == NULL)
        SendError("Server Error: Not enough core");
    }
    else
      break;

    PlusToSpace(FormItem);
    UnescapeURL(FormItem);

    ItemValue = SeparateWords(FormItem, '=');

    if      (str_icmp(FormItem, "argl") == 0)
    {
      if (*argl)
        SendError("Corrupt Input Form: Multiple \"argl\"");

      if (firstnonblank(ItemValue) != '\0')
      {
            (*argl) = str_dup(ItemValue);
        if ((*argl) == NULL)
          SendError("Server Error: Not enough core");
      }
    }
    else if (str_icmp(FormItem, "symxyz") == 0)
    {
      if (*symxyz)
        SendError("Corrupt Input Form: Multiple \"symxyz\"");

      if (firstnonblank(ItemValue) != '\0')
      {
            (*symxyz) = str_dup(ItemValue);
        if ((*symxyz) == NULL)
          SendError("Server Error: Not enough core");
      }
    }
    else
      SendError("Corrupt Input Form: Unknown item name");

    iword++;
  }

  SetFree(FormItem);
}


int main(int argc, char *argv[])
{
  int         ReqMethod, ContentLength;
  const char  *QueryString;

  char  *env;
  char  *argl, *symxyz, **run_cmd_argv;


#define GetEnv(s)\
        env = getenv( # s );\
    if (env)\
      s = env

  GetEnv(SGINFO_CMD_PATH);
  GetEnv(SGINFO_MAIN_URL);
  GetEnv(SGINFO_GIF_URL);
  GetEnv(SGINFO_QUERY_URL);

#undef GetEnv

  htmlHead(stdout, NULL, "sginfo gateway", NULL);
  fs("<h3>SgInfo - Space group Info</h3>\n");
  fs("<hr>\n");

  ReqMethod = RQM_NONE;
  ContentLength = 0;
  QueryString = NULL;

      env = getenv("REQUEST_METHOD");
  if (env)
  {
    if      (str_icmp(env, "GET") == 0)
    {
      ReqMethod = RQM_GET;

          env = getenv("QUERY_STRING");
      if (env == NULL)
        SendError("Unknown QUERY_STRING");

      QueryString = env;
    }
    else if (str_icmp(env, "POST") == 0)
    {
      ReqMethod = RQM_POST;

          env = getenv("CONTENT_TYPE");
      if (env == NULL || strcmp(env,"application/x-www-form-urlencoded"))
        SendError("Unknown CONTENT_TYPE");

          env = getenv("CONTENT_LENGTH");
      if (env == NULL)
        SendError("Unknown CONTENT_LENGTH");

          ContentLength = atoi(env);
      if (ContentLength < 0)
        SendError("Corrupt CONTENT_LENGTH");
    }
    else
      SendError("Unknown REQUEST_METHOD");
  }

  argl   = NULL;
  symxyz = NULL;

  GetCGI_Variables(ReqMethod, QueryString, ContentLength,
                   &argl, &symxyz);

  if (argl)
  {
    run_cmd_argv = mkargv(argl);

    if (symxyz)
      run_cmd_argv = addargv(run_cmd_argv, "-ReadXYZ");
  }
  else if (symxyz)
    run_cmd_argv = mkargv("-ReadXYZ");
  else
    run_cmd_argv = cpargv(argc, argv);

  if (argl || symxyz)
  {
    RunSginfo(run_cmd_argv, symxyz);
    fs("<hr>\n");
  }
  else if (ReqMethod == RQM_NONE || ReqMethod == RQM_GET)
  {
    fs("<img alt=\"SgInfo|ofnIgS\"");
    printf(" src=\"%s\">\n", SGINFO_GIF_URL);
    fs("<hr>\n");
  }

  fs("<form method=\"post\"");
  printf(" action=\"%s\">\n", SGINFO_QUERY_URL);
  fs("<code><strong>sginfo </strong></code><input\n");
  fs("name=\"argl\" size=60\n");
  fs("value=\"");

  if (argl && ! symxyz) str2html(stdout, argl);

  fs("\"><p>\n");
  fs("<input type=\"submit\" value=\"Run sginfo\">\n");
  fs("<input type=\"reset\"  value=\"Reset\">\n");
  printf("<strong><a href=\"%s\">", SGINFO_MAIN_URL);
  fs("On-line Documentation</a></strong>\n");

  if ((argl || symxyz) && (ReqMethod == RQM_NONE || ReqMethod == RQM_GET))
  {
    fs("<p>Running sginfo without arguments will\n");
    fs("   produce a little help screen.\n");
  }

  fs("<p>With option <code>-ReadXYZ</code>,\n");
  fs("   enter your symmetry operations below:\n");
  fs("<p>\n");
  fs("<textarea name=\"symxyz\" rows=8 cols=40></textarea>\n");
  fs("</form>\n");

  if (! (argl || symxyz))
  {
    fs("<hr>\n");
    RunSginfo(run_cmd_argv, symxyz);
  }

  htmlAddress(stdout);
  htmlTail(stdout);

  exit(0);
  return 0;
}
