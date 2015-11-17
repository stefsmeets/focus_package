#include <stdio.h>
#include <stdlib.h>

#include "cputime.h"


#if defined(__TURBOC__)

#include <time.h>

void InitTicks(void)
{
}

long GetTicks(T_Ticks *Ticks)
{
  time_t     sec_since_some_day;
  struct tm  *stm;


  sec_since_some_day = time(NULL);

  if (Ticks)
  {
    Ticks->sec_since_some_day = (long) sec_since_some_day;

    Ticks->wallclock  =
    Ticks->cpu_user   = clock();
    Ticks->cpu_system = 0;
    Ticks->ticks_per_second = CLK_TCK;

    stm = localtime(&sec_since_some_day);

    Ticks->year   = stm->tm_year;
    Ticks->month  = stm->tm_mon + 1;
    Ticks->day    = stm->tm_mday;
    Ticks->hour   = stm->tm_hour;
    Ticks->minute = stm->tm_min;
    Ticks->second = stm->tm_sec;
  }

  return (long) sec_since_some_day;
}


#elif (defined(__ALPHA) && defined(__VMS))

#include <time.h>

static time_t TicksRef;

void InitTicks(void)
{
  TicksRef = time(NULL);
}

long GetTicks(T_Ticks *Ticks)
{
  time_t     sec_since_some_day;
  tbuffer_t  stms;
  struct tm  *stm;


  sec_since_some_day = time(NULL);

  if (Ticks)
  {
    times(&stms);
    Ticks->sec_since_some_day = (long) sec_since_some_day;

    Ticks->wallclock  = (sec_since_some_day - TicksRef) * CLK_TCK;
    Ticks->cpu_user   = stms.proc_user_time;
    Ticks->cpu_system = stms.proc_system_time;
    Ticks->ticks_per_second = CLK_TCK;

    stm = localtime(&sec_since_some_day);

    Ticks->year   = stm->tm_year;
    Ticks->month  = stm->tm_mon + 1;
    Ticks->day    = stm->tm_mday;
    Ticks->hour   = stm->tm_hour;
    Ticks->minute = stm->tm_min;
    Ticks->second = stm->tm_sec;
  }

  return (long) sec_since_some_day;
}


#elif     defined(_CRAY)\
      || (defined(__alpha) && defined(__osf__))\
      ||  defined(__hpux)\
      ||  defined(_AIX)

#include <sys/types.h>
#include <sys/times.h>
#include <time.h>

static time_t TicksRef;

void InitTicks(void)
{
  struct tms  stms;

  TicksRef = times(&stms);
}

long GetTicks(T_Ticks *Ticks)
{
  time_t      sec_since_some_day;
  struct tms  stms;
  struct tm   *stm;


  sec_since_some_day = time(NULL);

  if (Ticks)
  {
    Ticks->sec_since_some_day = (long) sec_since_some_day;

    Ticks->wallclock  = times(&stms) - TicksRef;
    Ticks->cpu_user   = stms.tms_utime;
    Ticks->cpu_system = stms.tms_stime;
    Ticks->ticks_per_second = CLK_TCK;

    stm = localtime(&sec_since_some_day);

    Ticks->year   = stm->tm_year;
    Ticks->month  = stm->tm_mon + 1;
    Ticks->day    = stm->tm_mday;
    Ticks->hour   = stm->tm_hour;
    Ticks->minute = stm->tm_min;
    Ticks->second = stm->tm_sec;
  }

  return (long) sec_since_some_day;
}


#elif defined(__linux__) || defined(__APPLE__)

#include <sys/times.h>
#include <unistd.h>
#include <time.h>

static time_t TicksRef;
static long   TicksCLK_TCK;

void InitTicks(void)
{
  struct tms  stms;

  TicksRef     = times(&stms);
  TicksCLK_TCK = sysconf(_SC_CLK_TCK);
}

long GetTicks(T_Ticks *Ticks)
{
  time_t      sec_since_some_day;
  struct tms  stms;
  struct tm   *stm;


  sec_since_some_day = time(NULL);

  if (Ticks)
  {
    Ticks->sec_since_some_day = (long) sec_since_some_day;

    Ticks->wallclock  = times(&stms) - TicksRef;
    Ticks->cpu_user   = stms.tms_utime;
    Ticks->cpu_system = stms.tms_stime;
    Ticks->ticks_per_second = TicksCLK_TCK;

    stm = localtime(&sec_since_some_day);

    Ticks->year   = stm->tm_year;
    Ticks->month  = stm->tm_mon + 1;
    Ticks->day    = stm->tm_mday;
    Ticks->hour   = stm->tm_hour;
    Ticks->minute = stm->tm_min;
    Ticks->second = stm->tm_sec;
  }

  return (long) sec_since_some_day;
}


#elif defined(__sgi) || defined(__PARAGON__)

#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>
#include <time.h>

static time_t TicksRef;

void InitTicks(void)
{
  struct tms  stms;

  TicksRef = times(&stms);
}

long GetTicks(T_Ticks *Ticks)
{
  time_t      sec_since_some_day;
  struct tms  stms;
  struct tm   *stm;


  sec_since_some_day = time(NULL);

  if (Ticks)
  {
    Ticks->sec_since_some_day = (long) sec_since_some_day;

    Ticks->wallclock  = times(&stms) - TicksRef;
    Ticks->cpu_user   = stms.tms_utime;
    Ticks->cpu_system = stms.tms_stime;
    Ticks->ticks_per_second = HZ;

    stm = localtime(&sec_since_some_day);

    Ticks->year   = stm->tm_year;
    Ticks->month  = stm->tm_mon + 1;
    Ticks->day    = stm->tm_mday;
    Ticks->hour   = stm->tm_hour;
    Ticks->minute = stm->tm_min;
    Ticks->second = stm->tm_sec;
  }

  return (long) sec_since_some_day;
}


#else /* vanilla definition */

void InitTicks(void)
{
}

long GetTicks(T_Ticks *Ticks)
{
  if (Ticks)
  {
    Ticks->sec_since_some_day = 0;
    Ticks->wallclock  = 0;
    Ticks->cpu_user   = 0;
    Ticks->cpu_system = 0;
    Ticks->ticks_per_second = 1.;
    Ticks->year   = 0;
    Ticks->month  = 0;
    Ticks->day    = 0;
    Ticks->hour   = 0;
    Ticks->minute = 0;
    Ticks->second = 0;
  }

  return 0;
}


#endif /* end of machine dependent routines */


void PrintTicks(FILE *fpout, T_Ticks *Ticks, char *Head, char *Tail)
{
  T_Ticks     LocalTicks, *T;

  if (Ticks == NULL)
  {
    T = &LocalTicks;
    (void) GetTicks(T);
  }
  else
    T = Ticks;

  if (Head != NULL)
  {
    (void) fprintf(fpout, "%s", Head);
  }

  (void) fprintf(fpout, "%2.2d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d",
    T->year, T->month, T->day,
    T->hour, T->minute, T->second);

  (void) fprintf(fpout, "  %.2f  %.2f  %.2f  %.2f",
    (double) T->wallclock / T->ticks_per_second,
    (double) (T->cpu_user + T->cpu_system) / T->ticks_per_second,
    (double) T->cpu_user / T->ticks_per_second,
    (double) T->cpu_system / T->ticks_per_second);

  if (Tail != NULL)
  {
    (void) fprintf(fpout, "%s", Tail);
  }
}


#ifdef STAND_ALONE  /* main function for debugging purposes */

#include <math.h>

int main(int argc, char *argv[])
{
  int           i;
  long          l, n;
  double        d;
  T_Ticks       Ticks;


  InitTicks();

  if (argc != 2 || sscanf(argv[1], "%ld", &n) != 1)
  {
    (void) printf("usage: cputime #cycles\n");
  }
  else
  {
    (void) GetTicks(&Ticks);
    (void) fprintf(stdout, "ticks_per_second = %g\n", Ticks.ticks_per_second);
    PrintTicks(stdout, &Ticks, NULL, NULL);
    putc('\n', stdout);

    do
    {
      d = 1.;
      for (l = 0; l < n; l++) d = cos(d + l);

      (void) GetTicks(&Ticks);
      PrintTicks(stdout, &Ticks, NULL, NULL);
      putc('\n', stdout);

      i = getchar();
    }
    while (i != 'x' && i != 'X' && i != 'q' && i != 'Q');
  }

  return 0;
}

#endif  /* STAND_ALONE */
