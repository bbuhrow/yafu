/*
  Copyright (C) 2002 Jens Franke, T. Kleinjung.
  This file is part of gnfs4linux, distributed under the terms of the
  GNU General Public Licence and WITHOUT ANY WARRANTY.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/


#include <stdio.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <sys/timeb.h>
#include <stdlib.h>
#include <string.h>
#include "siever-config.h"
#include "../if.h"
#include "zeit.h"

clock_t *zeitcounter;
u64_t *asmzeitcounter;
double *zeitsum;
size_t zeitcounteranz;

u64_t asmgetclock(void);

void zeita(size_t i)
{
  zeitcounter[i]=clock();
  asmzeitcounter[i]-=asmgetclock();
}

#if 0
void zeitA(size_t i)
{
  asmzeitcounter[i]-=asmgetclock();
}
#endif

void zeitb(size_t i)
{
  zeitsum[i]+=((double)(clock()-zeitcounter[i]));
  asmzeitcounter[i]+=asmgetclock();
}

#if 0
void zeitB(size_t i)
{
  asmzeitcounter[i]+=asmgetclock();
}
#endif

void initzeit(size_t i)
{
  zeitsum=(double *)xmalloc(i*sizeof(double));
  memset(zeitsum,0,i*sizeof(double));
  zeitcounter=(clock_t *)xmalloc(i*sizeof(clock_t));
  memset(zeitcounter,0,i*sizeof(clock_t));
  asmzeitcounter=(u64_t *)xmalloc(i*sizeof(u64_t));
  memset(asmzeitcounter,0,i*sizeof(u64_t));
  zeitcounteranz=i;
}


#if 0
/* For use with mpqst.c */
void printzeit(void)
{
  long i;

  printf("\nZeit: ");
  for (i=0; i<zeitcounteranz; i++)
    printf("%ld: %.3fs  ",i,(double)(zeitcounter[i])/CLOCKS_PER_SEC);
  printf("\n");
#if 1
  printf("\nCycles: ");
  for (i=0; i<zeitcounteranz; i++)
    printf("%ld: %llu  ",i,asmzeitcounter[i]);
  printf("\n");
#endif
}
#else
void printzeit(size_t i)
{
  if(i>=zeitcounteranz) {
    fprintf(stderr,"Attempt to print time %zd of %zd\n",i,zeitcounteranz);
    abort();
  }
  printf("%zu: %.3fs (%llu) ",i,zeitsum[i]/CLOCKS_PER_SEC,asmzeitcounter[i]);
}

#endif
