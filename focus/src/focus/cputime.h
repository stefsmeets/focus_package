#ifndef A_CPUTIME_H__
#define A_CPUTIME_H__


typedef struct
  {
    long    sec_since_some_day;
    long    wallclock;
    long    cpu_user;
    long    cpu_system;
    double  ticks_per_second;
    int     year;
    int     month;
    int     day;
    int     hour;
    int     minute;
    int     second;
  }
  T_Ticks;


void InitTicks(void);
long GetTicks(T_Ticks *Ticks);
void PrintTicks(FILE *fpout, T_Ticks *Ticks, char *Head, char *Tail);


#endif /* A_CPUTIME_H__ */
