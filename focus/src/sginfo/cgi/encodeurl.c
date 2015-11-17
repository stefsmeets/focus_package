#include <stdio.h>
#include <stdlib.h>


#include "cgitbx.h"


int main(void)
{
  int  c;


  while ((c = getc(stdin)) != EOF)
    putcURL(c, stdout);
  
  exit(0);
  return 0;
}
