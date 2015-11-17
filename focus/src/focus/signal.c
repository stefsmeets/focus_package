#include <stdio.h>
#include <stdlib.h>
#include <signal.h>


#ifndef SIGHUP
#define SIGHUP -1
#else
static int  NoSIGHUP = 0;
#endif
#ifndef  SIGINT
#define SIGINT -2
#else
static int  NoSIGINT = 0;
#endif
#ifndef  SIGQUIT
#define SIGQUIT -3
#else
static int  NoSIGQUIT = 0;
#endif
#ifndef  SIGALRM
#define SIGALRM -14
#else
static int  NoSIGALRM = 0;
#endif
#ifndef  SIGTERM
#define SIGTERM -15
#else
static int  NoSIGTERM = 0;
#endif


#include "main.h"


static char *NameOfSignal(int sig)
{
  switch (sig)
  { case SIGHUP:    return "hangup";
    case SIGINT:    return "interrupt";
    case SIGQUIT:   return "quit";
    case SIGALRM:   return "alarm clock";
    case SIGTERM:   return "software termination signal";
    default:        break;
  }

  return "unknown signal number";
}

static void CatchSignal(int sig)
{
#if SIGINT >= 0
                (void) signal(SIGINT, SIG_IGN);
#endif
#if SIGHUP >= 0
                (void) signal(SIGHUP, SIG_IGN);
#endif
#if SIGQUIT >= 0
                (void) signal(SIGQUIT, SIG_IGN);
#endif
#if SIGALRM >= 0
                (void) signal(SIGALRM, SIG_IGN);
#endif
#if SIGTERM >= 0
                (void) signal(SIGTERM, SIG_IGN);
#endif

  switch (sig)
  {
    case SIGHUP:
    case SIGALRM:
      Fflush(stdout);
      AppSignal();
      break;

    case SIGINT:
    case SIGQUIT:
    case SIGTERM:
    default:
      Fprintf(stdout,
        "\n\n# Signal #%d: %s\n\n", sig, NameOfSignal(sig));

      if (sig == SIGINT)
      {
        QuitProgram = 1;
        NoSIGINT = 1;
        AppSignal();
      }
      else
        AppExit(2);

      break;
  }
}

void AppSignal(void)
{
#if SIGINT >= 0
             if (NoSIGINT)
                (void) signal(SIGINT, SIG_IGN);
             else
                (void) signal(SIGINT, CatchSignal);
#endif
#if SIGHUP >= 0
             if (NoSIGHUP)
                (void) signal(SIGHUP, SIG_IGN);
             else
                (void) signal(SIGHUP, CatchSignal);
#endif
#if SIGQUIT >= 0
             if (NoSIGQUIT)
                (void) signal(SIGQUIT, SIG_IGN);
             else
                (void) signal(SIGQUIT, CatchSignal);
#endif
#if SIGALRM >= 0
             if (NoSIGALRM)
                (void) signal(SIGALRM, SIG_IGN);
             else
                (void) signal(SIGALRM, CatchSignal);
#endif
#if SIGTERM >= 0
             if (NoSIGTERM)
                (void) signal(SIGTERM, SIG_IGN);
             else
                (void) signal(SIGTERM, CatchSignal);
#endif
}


void CheckSignalFile(void)
{
  int   c;
  FILE  *fpsig;


  if (F_SignalFileName)
  {
        fpsig = fopen(F_SignalFileName, "r");
    if (fpsig != NULL)
    {
      c = fgetc(fpsig);
      Fclose(fpsig);

      Fprintf(stdout,
        "\n\n# Signal file \"%s\" detected\n\n", F_SignalFileName);

      QuitProgram = 2;

      if (c == '!') AppExit(2);
    }
  }
}
