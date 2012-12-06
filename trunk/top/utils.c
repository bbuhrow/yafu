/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

Some parts of the code (and also this header), included in this 
distribution have been reused from other sources. In particular I 
have benefitted greatly from the work of Jason Papadopoulos's msieve @ 
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St. 
Denis Tom's Fast Math library.  Many thanks to their kind donation of 
code to the public domain.
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/

#ifdef _MSC_VER
#include <intrin.h>
#pragma intrinsic(__rdtsc)
#endif

#include "yafu.h"
#include "common.h"
#include "util.h"
#include "arith.h"
#include "yafu_string.h"

const char* szFeatures[] =
{
    "x87 FPU On Chip",
    "Virtual-8086 Mode Enhancement",
    "Debugging Extensions",
    "Page Size Extensions",
    "Time Stamp Counter",
    "RDMSR and WRMSR Support",
    "Physical Address Extensions",
    "Machine Check Exception",
    "CMPXCHG8B Instruction",
    "APIC On Chip",
    "Unknown1",
    "SYSENTER and SYSEXIT",
    "Memory Type Range Registers",
    "PTE Global Bit",
    "Machine Check Architecture",
    "Conditional Move/Compare Instruction",
    "Page Attribute Table",
    "36-bit Page Size Extension",
    "Processor Serial Number",
    "CFLUSH Extension",
    "Unknown2",
    "Debug Store",
    "Thermal Monitor and Clock Ctrl",
    "MMX Technology",
    "FXSAVE/FXRSTOR",
    "SSE Extensions",
    "SSE2 Extensions",
    "Self Snoop",
    "Multithreading Technology",
    "Thermal Monitor",
    "Unknown4",
    "Pending Break Enable"
};

fp_digit spRand(fp_digit lower, fp_digit upper)
{
	// advance the state of the LCG and return the appropriate result

#if BITS_PER_DIGIT == 64
	//we need to do some gymnastics to prevent the potentially 64 bit value
	//of (upper - lower) from being truncated to a 53 bit double

	uint32 n = spBits(upper-lower);

	if (n > 32)
	{
		fp_digit boundary = 4294967296ULL;
		fp_digit l,u;
		l = spRand(lower,boundary-1);
		u = spRand(0,upper>>32);
		return  l + (u << 32);
	}
	else
	{
		LCGSTATE = 6364136223846793005ULL * LCGSTATE + 1442695040888963407ULL;
		return lower + (fp_digit)(
			(double)(upper - lower) * (double)(LCGSTATE >> 32) / 4294967296.0);
	}

	
#else

	LCGSTATE = 6364136223846793005ULL * LCGSTATE + 1442695040888963407ULL;
	return lower + (fp_digit)(
			(double)(upper - lower) * (double)(LCGSTATE >> 32) / 4294967296.0);
#endif

}

void zRandb(z *n, int bits)
{
	//generate random number of 'bits' bits or less
	int full_words = bits/BITS_PER_DIGIT;
	int last_bits = bits%BITS_PER_DIGIT;
	fp_digit mask;
	int i;

	if (bits == 0)
	{
		zCopy(&zZero, n);
		return;
	}

	if (n->alloc < full_words + 1)
		zGrow(n,full_words+1);

	for (i=0;i<full_words;i++)
		n->val[i] = spRand(0,MAX_DIGIT);

	mask = HIBITMASK;
	if (last_bits == 0)
	{
		//then bits_per_digit divides bits, and we
		//need to ensure the top bit is set
		n->val[full_words-1] |= mask;
	}
	else
	{
		fp_digit last_word = spRand(0,MAX_DIGIT);
		//need to generate a last word
		mask = 1;
		for (i=0; i<last_bits; i++)
			mask |= mask << 1;

		n->val[full_words] = mask & last_word;
	}

	n->size = full_words + (last_bits > 0);

	//truncate leading 0 digits
	for(i = n->size - 1; i >= 0; i--)
	{
		if (n->val[i] != 0)
			break;
	}
	n->size = i+1;

	if (n->size == 0)
	{
		n->val[0] = 0;
		n->size = 1;
	}

	return;
}

void zRand(z *n, uint32 ndigits)
{
	int full_words = ndigits/DEC_DIGIT_PER_WORD;
	int rem = ndigits%DEC_DIGIT_PER_WORD;
	int i;
	
	if (ndigits == 0)
	{ 
		zCopy(&zZero, n);
		return;
	}

	if (full_words+1 > n->alloc)
		zGrow(n,full_words + LIMB_BLKSZ);	

	for (i=0;i<full_words;i++)
		n->val[i] = spRand(0,MAX_DIGIT);

	if (rem)
	{
		n->val[i] = spRand(0,(fp_digit)pow(10,rem));
		i++;
	}

	for(;i>=0;i--)
	{
		if (n->val[i] != 0)
			break;
	}
	n->size = i+1;

	if (n->size == 0)
	{
		n->val[0] = 0;
		n->size = 1;
	}

	return;
}


void sInit(str_t *s)
{
	s->s = (char *)malloc(GSTR_MAXSIZE * sizeof(char));
	if (s->s == NULL)
	{
		printf("couldn't allocate str_t in sInit\n");
		exit(-1);
	}
	s->s[0] = '\0';
	s->nchars = 1;
	s->alloc = GSTR_MAXSIZE;
	return;
}

void sFree(str_t *s)
{
	free(s->s);
	return;
}

void sClear(str_t *s)
{
	s->s[0] = '\0';
	s->nchars = 1;
	return;
}

void toStr(char *src, str_t *dest)
{
	if ((int)strlen(src) > dest->alloc)
	{
		sGrow(dest,strlen(src) + 10);
		dest->alloc = strlen(src) + 10;
	}
	memcpy(dest->s,src,strlen(src) * sizeof(char));
	dest->s[strlen(src)] = '\0';
	dest->nchars = strlen(src) + 1;

	return;
}

void sGrow(str_t *str, int size)
{
	//printf("growing str_t size...\n");
	str->s = (char *)realloc(str->s,size * sizeof(char));
	if (str->s == NULL)
	{
		printf("unable to reallocate string in sGrow\n");
		exit(-1);
	}
	str->alloc = size;

	return;
}

void sAppend(const char *src, str_t *dest)
{
	if (((int)strlen(src) + dest->nchars) >= dest->alloc)
	{
		sGrow(dest,strlen(src) + dest->nchars + 10);
		dest->alloc = strlen(src) + dest->nchars + 10;
	}

	memcpy(dest->s + dest->nchars - 1,src,strlen(src) * sizeof(char));
	dest->nchars += strlen(src);	//already has a null char accounted for
	dest->s[dest->nchars - 1] = '\0';

	return;
}

void sCopy(str_t *dest, str_t *src)
{
	if (dest->alloc < src->nchars + 2)
	{
		//free(dest->s);
		//dest->s = (char *)malloc((src->nchars+2) * sizeof(char));
		//printf("growing str in sCopy...\n");
		dest->s = (char *)realloc(dest->s, (src->nchars+2) * sizeof(char));
		dest->alloc = src->nchars+2;
	}
	memcpy(dest->s,src->s,src->nchars * sizeof(char));
	dest->nchars = src->nchars;
	return;
}

void logprint(FILE *infile, char *args, ...)
{
	va_list ap;
	time_t curtime;
	struct tm *loctime;
	char *s;

	if (!LOGFLAG)
		return;

	s = (char *)malloc(GSTR_MAXSIZE*sizeof(char));
	curtime = time(NULL);
	loctime = localtime(&curtime);

	if (infile == NULL)
		infile = stdout;

	//args has at least one argument, which may be NULL
	va_start(ap, args);

	//print the date and version stamp
	strftime(s,256,"%m/%d/%y %H:%M:%S",loctime);
	//fprintf(infile,"%s v%d.%02d @ %s, ",s,MAJOR_VER,MINOR_VER,sysname);
	fprintf(infile,"%s v%s @ %s, ",s,VERSION_STRING,sysname);

	vfprintf(infile,args,ap);

	va_end(ap);
	free(s);
	return;
}

void logprint_oc(const char *name, const char *method, char *args, ...)
{
	//print formatted string to logfile, given the file name
	//and access method
	va_list ap;
	time_t curtime;
	struct tm *loctime;
	char *s;
	FILE *logfile;

	if (!LOGFLAG)
		return;

	logfile = fopen(name, method);
	if (logfile == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("could not open logfile for method %c\n", method[0]);
		logfile = stdout;
	}

	s = (char *)malloc(GSTR_MAXSIZE*sizeof(char));
	curtime = time(NULL);
	loctime = localtime(&curtime);

	//args has at least one argument, which may be NULL
	va_start(ap, args);

	//print the date and version stamp
	strftime(s,256,"%m/%d/%y %H:%M:%S",loctime);
	//fprintf(infile,"%s v%d.%02d @ %s, ",s,MAJOR_VER,MINOR_VER,sysname);
	fprintf(logfile,"%s v%s @ %s, ",s,VERSION_STRING,sysname);

	vfprintf(logfile,args,ap);

	va_end(ap);
	free(s);

	if (logfile != NULL)
		fclose(logfile);

	return;
}

char *gettimever(char *s)
{
	time_t curtime;
	struct tm *loctime;

	curtime = time(NULL);
	loctime = localtime(&curtime);

	strcpy(s,"");
	//print the date stamp
	strftime(s,256,"%a %b %d %Y %H:%M:%S",loctime);
	//sprintf(s,"%s v%d.%d, ",s,MAJOR_VER,MINOR_VER);
	sprintf(s,"%s v%s, ",s,VERSION_STRING);

	return s;
}


void zNextPrime(z *n, z *p, int dir)
{
	//return the next prime after n, in the direction indicated by dir
	//if dir is positive, next higher else, next lower prime
	//n may get altered
	int szn = abs(n->size);

	if (szn > p->alloc)
		zGrow(p,szn + LIMB_BLKSZ);

	if (zCompare(n,&zThree) <= 0)
	{
		// handle special cases of really small inputs
		if (dir > 0)
		{
			if (n->val[0] == 3)
				sp2z(5,p);
			else if (n->val[0] == 2)
				sp2z(3,p);
			else
				sp2z(2,p);
		}
		else
		{
			sp2z(2,p);
		}
		return;
	}

	//make sure start is odd
	if (!(n->val[0] & 1))
	{
		if (dir>0)
			zShortAdd(n,1,p);
		else
			zShortSub(n,1,p);
		if (isPrime(p))
			return;
	}
	else
		zCopy(n,p);

	while (1)
	{			
		if (dir > 0)
		{
			zShortAdd(p,2,p);
			if (zShortMod(p,5) == 0)
				zShortAdd(p,2,p);
		}
		else
		{
			zShortSub(p,2,p);
			if (zShortMod(p,5) == 0 && zCompare(p,&zFive) != 0)
				zShortSub(p,2,p);
		}

		if (isPrime(p))
			return;
	}	

	return;
}

void zNextPrime_1(fp_digit n, fp_digit *p, z *work, int dir)
{
	//return the next prime after n, in the direction indicated by dir
	//if dir is positive, next higher else, next lower prime

	if (n <= 3)
	{
		// handle special cases of really small inputs
		if (dir > 0)
		{
			if (n == 3)
				*p = 5;
			else if (n == 2)
				*p = 3;
			else
				*p = 2;
		}
		else
		{
			*p = 2;
		}
		return;
	}

	//make sure start is odd
	if (!(n & 1))
	{
		if (dir>0)
			*p = n+1;
		else
			*p = n-1;
	}
	else
		*p = n;

	while (1)
	{			
		if (dir > 0)
		{
			*p += 2;
			if (((*p % 5) == 0))
				*p += 2;
		}
		else
		{
			*p -= 2;
			if (((*p % 5) == 0) && (*p != 5))
				*p -= 2;
		}

		sp2z(*p,work);
		work->type = UNKNOWN;
		if (isPrime(work))
			break;
	}	

	return;
}

void generate_pseudoprime_list(int num, int bits)
{
	//generate a list of 'num' pseudoprimes, each of size 'bits'
	//save to pseudoprimes.dat

	FILE *out;
	int i;
	z tmp1,tmp2,tmp3;

	zInit(&tmp1);
	zInit(&tmp2);
	zInit(&tmp3);

	out = fopen("pseudoprimes.dat","w");
	if (out == NULL)
	{
		printf("couldn't open pseudoprimes.dat for writing\n");
		return;
	}
	
	for (i=0; i<num; i++)
	{
		zRandb(&tmp3,bits/2);
		zNextPrime(&tmp3,&tmp1,1);
		zRandb(&tmp3,bits/2);
		zNextPrime(&tmp3,&tmp2,1);
		zMul(&tmp1,&tmp2,&tmp3);
		fprintf(out,"%s,%s,%s\n",
			z2decstr(&tmp3,&gstr1),
			z2decstr(&tmp1,&gstr2),
			z2decstr(&tmp2,&gstr3));
	}
	fclose(out);

	printf("generated %d pseudoprimes in file pseudoprimes.dat\n",num);
	zFree(&tmp1);
	zFree(&tmp2);
	zFree(&tmp3);
	return;
}

void gordon(int bits, z *p)
{
	//find a random strong prime of size 'bits'
	//follows the Handbook of applied cryptography
	/*
	SUMMARY: a strong prime p is generated.
	1. Generate two large random primes s and t of roughly equal bitlength (see Note 4.54).
	2. Select an integer i0. Find the first prime in the sequence 2it + 1, for i = i0; i0 +
		1; i0 + 2; : : : (see Note 4.54). Denote this prime by r = 2it+ 1.
	3. Compute p0 = 2(sr-2 mod r)s - 1.
	4. Select an integer j0. Find the first prime in the sequence p0 +2jrs, for j = j0; j0 +
		1; j0 + 2; : : : (see Note 4.54). Denote this prime by p = p0 + 2jrs.
	5. Return(p).

  4.54 Note (implementing Gordon’s algorithm)
	(i) The primes s and t required in step 1 can be probable primes generated by Algorithm
	4.44. TheMiller-Rabin test (Algorithm 4.24) can be used to test each candidate
	for primality in steps 2 and 4, after ruling out candidates that are divisible by a small
	prime less than some boundB. See Note 4.45 for guidance on selecting B. Since the
	Miller-Rabin test is a probabilistic primality test, the output of this implementation
	of Gordon’s algorithm is a probable prime.
	(ii) By carefully choosing the sizes of primes s, t and parameters i0, j0, one can control
	the exact bitlength of the resulting prime p. Note that the bitlengths of r and s will
	be about half that of p, while the bitlength of t will be slightly less than that of r.
	*/
	
	int i,j,s_len,n_words;
	z s,t,r,tmp,tmp2,p0;

	zInit(&s);
	zInit(&t);
	zInit(&r);
	zInit(&tmp);
	zInit(&tmp2);
	zInit(&p0);

	//need to check allocation of tmp vars.  how big do they get?

	//1. s and t should be about half the bitlength of p
	//random s of bitlength bits/2
	s_len = bits/2 - 4;
	n_words = s_len/BITS_PER_DIGIT;
	s.size = n_words + (s_len % BITS_PER_DIGIT != 0);
	for (i=0; i<s_len; i++)
		s.val[i/BITS_PER_DIGIT] |= (fp_digit)(rand() & 0x1) << (fp_digit)(i%BITS_PER_DIGIT);
	
	//force the most significant bit = 1, to keep the bitlength about right
	s.val[s.size - 1] |= (fp_digit)0x1 << (fp_digit)((s_len-1)%BITS_PER_DIGIT);

	//random t of bitlength bits/2
	t.size = n_words + (s_len % BITS_PER_DIGIT != 0);
	for (i=0; i<s_len; i++)
		t.val[i/BITS_PER_DIGIT] |= (fp_digit)(rand() & 0x1) << (fp_digit)(i%BITS_PER_DIGIT);
	
	//force the most significant bit = 1, to keep the bitlength about right
	t.val[t.size - 1] |= (fp_digit)0x1 << (fp_digit)((s_len-1)%BITS_PER_DIGIT);

	//2. Select an integer i0. Find the first prime in the sequence 2it + 1, for i = i0; i0 +
	//1; i0 + 2; : : : (see Note 4.54). Denote this prime by r = 2it+ 1.
	i=1;
	zShiftLeft(&r,&t,1);
	zShortAdd(&r,1,&r);
	while (!isPrime(&r))
	{
		i++;
		zShiftLeft(&r,&t,1);
		zShortMul(&r,i,&r);
		zShortAdd(&r,1,&r);
	}

	//3. Compute p0 = 2(sr-2 mod r)s - 1.
	zMul(&s,&r,&tmp);
	zShortSub(&tmp,2,&p0);
	zDiv(&p0,&r,&tmp,&tmp2);
	zMul(&tmp2,&s,&p0);
	zShiftLeft(&p0,&p0,1);
	zShortSub(&p0,1,&p0);

	//4. Select an integer j0. Find the first prime in the sequence p0 +2jrs, for j = j0; j0 +
	//1; j0 + 2; : : : (see Note 4.54). Denote this prime by p = p0 + 2jrs.
	j=1;
	zMul(&r,&s,&tmp);
	zShiftLeft(&tmp,&tmp,1);
	zAdd(&p0,&tmp,p);
	while (!isPrime(p))
	{
		j++;
		zMul(&r,&s,&tmp);
		zShiftLeft(&tmp,&tmp,1);
		zShortMul(&tmp,j,&tmp);
		zAdd(&p0,&tmp,p);
	}

	zFree(&s);
	zFree(&t);
	zFree(&r);
	zFree(&tmp);
	zFree(&tmp2);
	zFree(&p0);
	return;
}

void build_RSA(int bits, z *n)
{
	z p,q;
	int i;
	int words,subwords;

	if (bits < 65)
	{
		printf("bitlength to small\n");
		return;
	}

	zInit(&p);
	zInit(&q);

	words = bits/BITS_PER_DIGIT + (bits%BITS_PER_DIGIT != 0);
	subwords = (int)ceil((double)words/2.0);

	if (subwords > p.alloc)
		zGrow(&p,subwords + LIMB_BLKSZ);

	if (subwords > q.alloc)
		zGrow(&q,subwords + LIMB_BLKSZ);

	if (words > n->alloc)
		zGrow(n,words + LIMB_BLKSZ);

	zClear(n);
	i=0;
	while (zBits(n) != bits)
	{
		gordon(bits/2,&p);
		gordon(bits/2,&q);
		zMul(&p,&q,n);
		//printf("trial %d has %d bits\n",i,zBits(n));
		i++;
	}

	//printf("took %d trials\n",i);
	zFree(&p);
	zFree(&q);
	return;
}

//user dimis:
//http://cboard.cprogramming.com/cplusplus-programming/
//101085-how-measure-time-multi-core-machines-pthreads.html
//
TIME_DIFF * my_difftime (struct timeval * start, struct timeval * end)
{
	TIME_DIFF * diff = (TIME_DIFF *) malloc ( sizeof (TIME_DIFF) );

	if (start->tv_sec == end->tv_sec) {
		diff->secs = 0;
		diff->usecs = end->tv_usec - start->tv_usec;
	}
	else {
		diff->usecs = 1000000 - start->tv_usec;
		diff->secs = end->tv_sec - (start->tv_sec + 1);
		diff->usecs += end->tv_usec;
		if (diff->usecs >= 1000000) {
			diff->usecs -= 1000000;
			diff->secs += 1;
		}
	}
	
	return diff;
}

//http://www.openasthra.com/c-tidbits/gettimeofday-function-for-windows/
#if defined (_MSC_VER)
	int gettimeofday(struct timeval *tv, struct timezone *tz)
	{
	  FILETIME ft;
	  unsigned __int64 tmpres = 0;
	  static int tzflag;
	 
	  if (NULL != tv)
	  {
		GetSystemTimeAsFileTime(&ft);
	 
		tmpres |= ft.dwHighDateTime;
		tmpres <<= 32;
		tmpres |= ft.dwLowDateTime;
	 
		/*converting file time to unix epoch*/
		tmpres /= 10;  /*convert into microseconds*/
		tmpres -= DELTA_EPOCH_IN_MICROSECS; 
		tv->tv_sec = (long)(tmpres / 1000000UL);
		tv->tv_usec = (long)(tmpres % 1000000UL);
	  }
	 
	  if (NULL != tz)
	  {
		if (!tzflag)
		{
		  _tzset();
		  tzflag++;
		}
		tz->tz_minuteswest = _timezone / 60;
		tz->tz_dsttime = _daylight;
	  }
	 
	  return 0;
	}
#endif

int qcomp_uint16(const void *x, const void *y)
{
	uint16 *xx = (uint16 *)x;
	uint16 *yy = (uint16 *)y;
	
	if (*xx > *yy)
		return 1;
	else if (*xx == *yy)
		return 0;
	else
		return -1;
}

int qcomp_uint32(const void *x, const void *y)
{
	uint32 *xx = (uint32 *)x;
	uint32 *yy = (uint32 *)y;
	
	if (*xx > *yy)
		return 1;
	else if (*xx == *yy)
		return 0;
	else
		return -1;
}

int qcomp_uint64(const void *x, const void *y)
{
	uint64 *xx = (uint64 *)x;
	uint64 *yy = (uint64 *)y;
	
	if (*xx > *yy)
		return 1;
	else if (*xx == *yy)
		return 0;
	else
		return -1;
}

int qcomp_int(const void *x, const void *y)
{
	int *xx = (int *)x;
	int *yy = (int *)y;
	
	if (*xx > *yy)
		return 1;
	else if (*xx == *yy)
		return 0;
	else
		return -1;
}

#ifdef _MSC_VER

	/* Core aware timing on Windows, courtesy of Brian Gladman */

	#if defined( _WIN64 )

		#define current_processor_number GetCurrentProcessorNumber

	#else

		unsigned long current_processor_number(void)
		{
			__asm
			{
				mov     eax,1
				cpuid
				shr     ebx,24
				mov     eax, ebx
			}
		}

	#endif

	int lock_thread_to_core(void)
	{   DWORD_PTR afp, afs;

		if(GetProcessAffinityMask(GetCurrentProcess(), &afp, &afs))
		{
			afp &= (DWORD_PTR)(1 << current_processor_number());
			if(SetThreadAffinityMask(GetCurrentThread(), afp))
				return EXIT_SUCCESS;
		}
		return EXIT_FAILURE;
	}

	int unlock_thread_from_core(void)
	{   DWORD_PTR afp, afs;

        if(GetProcessAffinityMask(GetCurrentProcess(), &afp, &afs))
		{
			if(SetThreadAffinityMask(GetCurrentThread(), afp))
				return EXIT_SUCCESS;
		}
		return EXIT_FAILURE;
	}

    double cycles_per_second = 0.0;
    double ticks_per_second = 0.0;
    double cycles_per_tick = 0.0;

    uint64 measure_processor_speed(void)
    {   unsigned long long cycles;

        lock_thread_to_core();
        cycles = __rdtsc();
        Sleep(100);
        cycles = __rdtsc() - cycles;
        unlock_thread_from_core();
        cycles_per_second = 10.0 * (double)cycles;

        if(ticks_per_second == 0.0)
        {   LARGE_INTEGER ll;
            QueryPerformanceFrequency(&ll);
            ticks_per_second = (double)ll.QuadPart;
            cycles_per_tick = cycles_per_second / ticks_per_second;
        }
        return cycles;
    }

    double get_tsc_time(void)
    {
        if(cycles_per_second == 0.0)
            measure_processor_speed(); 
        return __rdtsc() / cycles_per_second;
    }

    double get_pfc_time(void)
    {   LARGE_INTEGER ll;

        if(ticks_per_second == 0.0)
            measure_processor_speed(); 
        QueryPerformanceCounter(&ll);
        return ll.QuadPart / ticks_per_second;
    }

#else

	uint64 measure_processor_speed(void)
	{   
		uint64 cycles;
		struct timeval start, stop;
		double t_time;
		TIME_DIFF *	difference;

		gettimeofday(&start,NULL);

		cycles = yafu_read_clock(); 
		do
		{
			gettimeofday (&stop, NULL);
			difference = my_difftime (&start, &stop);
			t_time = ((double)difference->secs + 
				(double)difference->usecs / 1000000);
			free(difference);
		}
		while (t_time < 0.1);
		cycles = yafu_read_clock() - cycles;

		return cycles;                  /* return cycles per second  */
	}

#endif


//functions below here courtesy of Jason Papadopoulos

uint64 yafu_read_clock(void) 
{
#if defined(GCC_ASM32X) || defined(GCC_ASM64X) 
	uint32 lo, hi;
	ASM_G("rdtsc":"=d"(hi),"=a"(lo));
	return (uint64)hi << 32 | lo;

#elif defined(_MSC_VER)
    LARGE_INTEGER ll;
    QueryPerformanceCounter(&ll);
	return (uint64)(ll.QuadPart * cycles_per_tick);
#else
	struct timeval thistime;   
	gettimeofday(&thistime, NULL);
	return (uint64)(cycles_per_second * 
        (thistime.tv_sec + thistime.tv_usec / 1000000.0));
#endif
}

/* macro to execute the x86 CPUID instruction. Note that
   this is more verbose than it needs to be; Intel Macs reserve
   the EBX or RBX register for the PIC base address, and so
   this register cannot get clobbered by inline assembly */

typedef union {
	uint32 data;

	struct {
		uint32 cache_type : 5;
		uint32 cache_level : 3;
		uint32 i_dont_care : 24;
	} s;
} cache_type_t;

typedef union {
	uint32 data;

	struct {
		uint32 line_size : 12;
		uint32 num_lines : 10;
		uint32 ways : 10;
	} s;
} cache_size_t;

#if defined(GCC_ASM32X)
	#define HAS_CPUID
	#define CPUID(code, a, b, c, d) 			\
		ASM_G volatile(					\
			"movl %%ebx, %%esi   \n\t"		\
			"cpuid               \n\t"		\
			"movl %%ebx, %1      \n\t"		\
			"movl %%esi, %%ebx   \n\t"		\
			:"=a"(a), "=m"(b), "=c"(c), "=d"(d) 	\
			:"0"(code) : "%esi")
	#define CPUID2(code1, code2, a, b, c, d) 			\
		ASM_G volatile(					\
			"movl %%ebx, %%esi   \n\t"		\
			"cpuid               \n\t"		\
			"movl %%ebx, %1      \n\t"		\
			"movl %%esi, %%ebx   \n\t"		\
			:"=a"(a), "=m"(b), "=c"(c), "=d"(d) 	\
			:"0"(code1), "2"(code2) : "%esi")

#elif defined(GCC_ASM64X)
	#define HAS_CPUID
	#define CPUID(code, a, b, c, d) 			\
		ASM_G volatile(					\
			"movq %%rbx, %%rsi   \n\t"		\
			"cpuid               \n\t"		\
			"movl %%ebx, %1      \n\t"		\
			"movq %%rsi, %%rbx   \n\t"		\
			:"=a"(a), "=m"(b), "=c"(c), "=d"(d) 	\
			:"0"(code) : "%rsi")
	#define CPUID2(code1, code2, a, b, c, d)		\
		ASM_G volatile(					\
			"movq %%rbx, %%rsi   \n\t"		\
			"cpuid               \n\t"		\
			"movl %%ebx, %1      \n\t"		\
			"movq %%rsi, %%rbx   \n\t"		\
			:"=a"(a), "=m"(b), "=c"(c), "=d"(d) 	\
			:"0"(code1), "2"(code2) : "%rsi")

#elif defined(_MSC_VER)
	#include <intrin.h>
	#define HAS_CPUID
	#define CPUID(code, a, b, c, d)	\
	{	uint32 _z[4]; \
		__cpuid(_z, code); \
		a = _z[0]; \
		b = _z[1]; \
		c = _z[2]; \
		d = _z[3]; \
	}
	#define CPUID2(code1, code2, a, b, c, d) \
	{	uint32 _z[4]; \
		__cpuidex(_z, code1, code2); \
		a = _z[0]; \
		b = _z[1]; \
		c = _z[2]; \
		d = _z[3]; \
	}
#endif


void yafu_get_cache_sizes(uint32 *level1_size_out,
			uint32 *level2_size_out) {

	/* attempt to automatically detect the size of
	   the L2 cache; this helps tune the choice of
	   parameters or algorithms used in the sieve-
	   based methods. It should guess right for most 
	   PCs and Macs when using gcc.

	   Otherwise, you have the source so just fill in
	   the correct number. */

	uint32 cache_size1 = DEFAULT_L1_CACHE_SIZE; 
	uint32 cache_size2 = DEFAULT_L2_CACHE_SIZE; 

#if defined(HAS_CPUID) && !defined(__APPLE__)

	/* reading the CPU-specific features of x86
	   processors is a simple 57-step process.
	   The following should be able to retrieve
	   the L1/L2/L3 cache size of any Intel or AMD
	   processor made after ~1995 */

	uint32 a, b, c, d;
	uint8 is_intel, is_amd;

	CPUID(0, a, b, c, d);
	is_intel = ((b & 0xff) == 'G');		/* "GenuineIntel" */
	is_amd = ((b & 0xff) == 'A');		/* "AuthenticAMD" */

	if (is_intel && a >= 2) {

		uint32 i; 
		uint8 features[15];
		uint32 max_special;
		uint32 j1 = 0;
		uint32 j2 = 0;

		/* handle newer Intel */

		if (a >= 4) {
			for (i = 0; i < 100; i++) {
				uint32 num_sets;
				cache_type_t type;
				cache_size_t size;

				CPUID2(4, i, type.data, size.data, num_sets, d);

				/* must be data cache or unified cache */

				if (type.s.cache_type == 0)
					break;
				else if (type.s.cache_type != 1 &&
					 type.s.cache_type != 3)
					continue;

				d = (size.s.line_size + 1) *
				    (size.s.num_lines + 1) *
				    (size.s.ways + 1) *
				    (num_sets + 1);

				if (type.s.cache_level == 1)
					j1 = MAX(j1, d);
				else
					j2 = MAX(j2, d);
			}
		}

		CPUID(0x80000000, max_special, b, c, d);
		if (max_special >= 0x80000006) {
			CPUID(0x80000006, a, b, c, d);
			j2 = MAX(j2, 1024 * (c >> 16));
		}

		/* handle older Intel, possibly overriding the above */

		CPUID(2, a, b, c, d);

		features[0] = (a >> 8);
		features[1] = (a >> 16);
		features[2] = (a >> 24);
		features[3] = b;
		features[4] = (b >> 8);
		features[5] = (b >> 16);
		features[6] = (b >> 24);
		features[7] = c;
		features[8] = (c >> 8);
		features[9] = (c >> 16);
		features[10] = (c >> 24);
		features[11] = d;
		features[12] = (d >> 8);
		features[13] = (d >> 16);
		features[14] = (d >> 24);

		/* use the maximum of the (known) L2 and L3 cache sizes */

		for (i = 0; i < sizeof(features); i++) {
			switch (features[i]) {
			/* level 1 cache codes */
			case 0x06:
			case 0x0a:
			case 0x66:
				j1 = MAX(j1, 8*1024); break;
			case 0x08:
			case 0x0c:
			case 0x0d:
			case 0x60:
			case 0x67:
				j1 = MAX(j1, 16*1024); break;
			case 0x0e:
				j1 = MAX(j1, 24*1024); break;
			case 0x09:
			case 0x2c:
			case 0x30:
			case 0x68:
				j1 = MAX(j1, 32*1024); break;

			/* level 2 and level 3 cache codes */
			case 0x41:
			case 0x79:
				j2 = MAX(j2, 128*1024); break;
			case 0x21:
			case 0x42:
			case 0x7a:
			case 0x82:
				j2 = MAX(j2, 256*1024); break;
			case 0x22:
			case 0x43:
			case 0x7b:
			case 0x7f:
			case 0x80:
			case 0x83:
			case 0x86:
				j2 = MAX(j2, 512*1024); break;
			case 0x23:
			case 0x44:
			case 0x78:
			case 0x7c:
			case 0x84:
			case 0x87:
				j2 = MAX(j2, 1*1024*1024); break;
			case 0x25:
			case 0x45:
			case 0x7d:
			case 0x85:
				j2 = MAX(j2, 2*1024*1024); break;
			case 0x48:
				j2 = MAX(j2, 3*1024*1024); break;
			case 0x29:
			case 0x46:
			case 0x49:
				j2 = MAX(j2, 4*1024*1024); break;
			case 0x4a:
			case 0x4e:
				j2 = MAX(j2, 6*1024*1024); break;
			case 0x47:
			case 0x4b:
			case 0xe4:
				j2 = MAX(j2, 8*1024*1024); break;
			case 0x4c:
			case 0xea:
				j2 = MAX(j2, 12*1024*1024); break;
			case 0x4d:
				j2 = MAX(j2, 16*1024*1024); break;
			case 0xeb:
				j2 = MAX(j2, 18*1024*1024); break;
			case 0xec:
				j2 = MAX(j2, 24*1024*1024); break;
			}
		}
		if (j1 > 0)
			cache_size1 = j1;
		if (j2 > 0)
			cache_size2 = j2;
	}
	else if (is_amd) {

		uint32 max_special;
		CPUID(0x80000000, max_special, b, c, d);

		if (max_special >= 0x80000005) {
			CPUID(0x80000005, a, b, c, d);
			cache_size1 = 1024 * (c >> 24);

			if (max_special >= 0x80000006) {
				CPUID(0x80000006, a, b, c, d);
				cache_size2 = MAX(1024 * (c >> 16),
						  512 * 1024 * (d >> 18));
			}
		}
	}
#endif

	*level1_size_out = cache_size1;
	*level2_size_out = cache_size2;
}


/*--------------------------------------------------------------------*/
enum cpu_type yafu_get_cpu_type(void) {

	enum cpu_type cpu = cpu_generic;

#if defined(HAS_CPUID)
	uint32 a, b, c, d;

	CPUID(0, a, b, c, d);
	if ((b & 0xff) == 'G') {	/* "GenuineIntel" */

		uint8 family, model;

		switch (a) {
		case 1:
			cpu = cpu_pentium;
			break;
		case 2:
			CPUID(1, a, b, c, d);
			family = (a >> 8) & 0xf;
			model = (a >> 4) & 0xf;
			if (family == 6) {
				if (model == 9 || model == 13)
					cpu = cpu_pentium_m;
				else
					cpu = cpu_pentium2;
			}
			else if (family == 15) {
				cpu = cpu_pentium4;
			}
			break;
		case 3:
			cpu = cpu_pentium3;
			break;
		case 5:
		case 6:
			cpu = cpu_pentium4;
			break;
		default:  
			/* a = 10+; some subspecies of core or core2 */
			CPUID(1, a, b, c, d);
			family = (a >> 8) & 0xf;
			model = (a >> 4) & 0xf;

			if (model >= 10)
				cpu = cpu_core; //cpu_nehalem;
			else
				cpu = cpu_core;

			break;
		}
	}
	else if ((b & 0xff) == 'A') {		/* "AuthenticAMD" */

		uint8 family, model;

		CPUID(1, a, b, c, d);
		family = (a >> 8) & 0xf;
		model = (a >> 4) & 0xf;
		if (family == 15)
			cpu = cpu_opteron;
		else if (family == 6) {
			CPUID(0x80000001, a, b, c, d);
			if (d & 0x1000000)		/* full SSE */
				cpu = cpu_athlon_xp;
			else				/* partial SSE */
				cpu = cpu_athlon;
		}
	}
#endif

	return cpu;
}

// http://msdn.microsoft.com/en-us/library/hskdteyh.aspx
// cpuid.cpp 
// processor: x86, x64
// Use the __cpuid intrinsic to get information about a CPU
// modified for c compliers and use of CPUID macros 
//		- brb, 10/26/10



int extended_cpuid(char *idstr, int *cachelinesize, char *bSSE41Extensions, int do_print)
{
    char CPUString[0x20];
    char CPUBrandString[0x40];
    int CPUInfo[4] = {-1};
    int nSteppingID = 0;
    int nModel = 0;
    int nFamily = 0;
    int nProcessorType = 0;
    int nExtendedmodel = 0;
    int nExtendedfamily = 0;
    int nBrandIndex = 0;
    int nCLFLUSHcachelinesize = 0;
    int nLogicalProcessors = 0;
    int nAPICPhysicalID = 0;
    int nFeatureInfo = 0;
    int nCacheLineSize = 0;
    int nL2Associativity = 0;
    int nCacheSizeK = 0;
    int nPhysicalAddress = 0;
    int nVirtualAddress = 0;
    int nRet = 0;

    int nCores = 0;
    int nCacheType = 0;
    int nCacheLevel = 0;
    int nMaxThread = 0;
    int nSysLineSize = 0;
    int nPhysicalLinePartitions = 0;
    int nWaysAssociativity = 0;
    int nNumberSets = 0;

    unsigned    nIds, nExIds, i;

    char    bSSE3Instructions = 0;
    char    bMONITOR_MWAIT = 0;
    char    bCPLQualifiedDebugStore = 0;
    char    bVirtualMachineExtensions = 0;
    char    bEnhancedIntelSpeedStepTechnology = 0;
    char    bThermalMonitor2 = 0;
    char    bSupplementalSSE3 = 0;
    char    bL1ContextID = 0;
    char    bCMPXCHG16B = 0;
    char    bxTPRUpdateControl = 0;
    char    bPerfDebugCapabilityMSR = 0;
    //char    bSSE41Extensions = 0;
    char    bSSE42Extensions = 0;
    char    bPOPCNT = 0;

    char    bMultithreading = 0;

    char    bLAHF_SAHFAvailable = 0;
    char    bCmpLegacy = 0;
    char    bSVM = 0;
    char    bExtApicSpace = 0;
    char    bAltMovCr8 = 0;
    char    bLZCNT = 0;
    char    bSSE4A = 0;
    char    bMisalignedSSE = 0;
    char    bPREFETCH = 0;
    char    bSKINITandDEV = 0;
    char    bSYSCALL_SYSRETAvailable = 0;
    char    bExecuteDisableBitAvailable = 0;
    char    bMMXExtensions = 0;
    char    bFFXSR = 0;
    char    b1GBSupport = 0;
    char    bRDTSCP = 0;
    char    b64Available = 0;
    char    b3DNowExt = 0;
    char    b3DNow = 0;
    char    bNestedPaging = 0;
    char    bLBRVisualization = 0;
    char    bFP128 = 0;
    char    bMOVOptimization = 0;

    char    bSelfInit = 0;
    char    bFullyAssociative = 0;

	*bSSE41Extensions = 0;

    // __cpuid with an InfoType argument of 0 returns the number of
    // valid Ids in CPUInfo[0] and the CPU identification string in
    // the other three array elements. The CPU identification string is
    // not in linear order. The code below arranges the information 
    // in a human readable form.
	CPUID(0,CPUInfo[0],CPUInfo[1],CPUInfo[2],CPUInfo[3]);
    //__cpuid(CPUInfo, 0);
    nIds = CPUInfo[0];
    memset(CPUString, 0, sizeof(CPUString));
    *((int*)CPUString) = CPUInfo[1];
    *((int*)(CPUString+4)) = CPUInfo[3];
    *((int*)(CPUString+8)) = CPUInfo[2];

    // Get the information associated with each valid Id
    for (i=0; i<=nIds; ++i)
    {
		CPUID(i,CPUInfo[0],CPUInfo[1],CPUInfo[2],CPUInfo[3]);
        //__cpuid(CPUInfo, i);

		if (do_print)
		{
			printf("\nFor InfoType %d\n", i); 
			printf("CPUInfo[0] = 0x%x\n", CPUInfo[0]);
			printf("CPUInfo[1] = 0x%x\n", CPUInfo[1]);
			printf("CPUInfo[2] = 0x%x\n", CPUInfo[2]);
			printf("CPUInfo[3] = 0x%x\n", CPUInfo[3]);
		}

        // Interpret CPU feature information.
        if  (i == 1)
        {
            nSteppingID = CPUInfo[0] & 0xf;
            nModel = (CPUInfo[0] >> 4) & 0xf;
            nFamily = (CPUInfo[0] >> 8) & 0xf;
            nProcessorType = (CPUInfo[0] >> 12) & 0x3;
            nExtendedmodel = (CPUInfo[0] >> 16) & 0xf;
            nExtendedfamily = (CPUInfo[0] >> 20) & 0xff;
            nBrandIndex = CPUInfo[1] & 0xff;
            *cachelinesize = nCLFLUSHcachelinesize = ((CPUInfo[1] >> 8) & 0xff) * 8;
            nLogicalProcessors = ((CPUInfo[1] >> 16) & 0xff);
            nAPICPhysicalID = (CPUInfo[1] >> 24) & 0xff;
            bSSE3Instructions = (CPUInfo[2] & 0x1) || 0;
            bMONITOR_MWAIT = (CPUInfo[2] & 0x8) || 0;
            bCPLQualifiedDebugStore = (CPUInfo[2] & 0x10) || 0;
            bVirtualMachineExtensions = (CPUInfo[2] & 0x20) || 0;
            bEnhancedIntelSpeedStepTechnology = (CPUInfo[2] & 0x80) || 0;
            bThermalMonitor2 = (CPUInfo[2] & 0x100) || 0;
            bSupplementalSSE3 = (CPUInfo[2] & 0x200) || 0;
            bL1ContextID = (CPUInfo[2] & 0x300) || 0;
            bCMPXCHG16B= (CPUInfo[2] & 0x2000) || 0;
            bxTPRUpdateControl = (CPUInfo[2] & 0x4000) || 0;
            bPerfDebugCapabilityMSR = (CPUInfo[2] & 0x8000) || 0;
            *bSSE41Extensions = (CPUInfo[2] & 0x80000) || 0;
            bSSE42Extensions = (CPUInfo[2] & 0x100000) || 0;
            bPOPCNT= (CPUInfo[2] & 0x800000) || 0;
            nFeatureInfo = CPUInfo[3];
            bMultithreading = (nFeatureInfo & (1 << 28)) || 0;
        }
    }

    // Calling __cpuid with 0x80000000 as the InfoType argument
    // gets the number of valid extended IDs.
    CPUID(0x80000000,CPUInfo[0],CPUInfo[1],CPUInfo[2],CPUInfo[3]);
	//__cpuid(CPUInfo, 0x80000000);
    nExIds = CPUInfo[0];
    memset(CPUBrandString, 0, sizeof(CPUBrandString));

    // Get the information associated with each extended ID.
    for (i=0x80000000; i<=nExIds; ++i)
    {
        CPUID(i,CPUInfo[0],CPUInfo[1],CPUInfo[2],CPUInfo[3]);
		//__cpuid(CPUInfo, i);
		if (do_print)
		{
			printf("\nFor InfoType %x\n", i); 
			printf("CPUInfo[0] = 0x%x\n", CPUInfo[0]);
			printf("CPUInfo[1] = 0x%x\n", CPUInfo[1]);
			printf("CPUInfo[2] = 0x%x\n", CPUInfo[2]);
			printf("CPUInfo[3] = 0x%x\n", CPUInfo[3]);
		}

        if  (i == 0x80000001)
        {
            bLAHF_SAHFAvailable = (CPUInfo[2] & 0x1) || 0;
            bCmpLegacy = (CPUInfo[2] & 0x2) || 0;
            bSVM = (CPUInfo[2] & 0x4) || 0;
            bExtApicSpace = (CPUInfo[2] & 0x8) || 0;
            bAltMovCr8 = (CPUInfo[2] & 0x10) || 0;
            bLZCNT = (CPUInfo[2] & 0x20) || 0;
            bSSE4A = (CPUInfo[2] & 0x40) || 0;
            bMisalignedSSE = (CPUInfo[2] & 0x80) || 0;
            bPREFETCH = (CPUInfo[2] & 0x100) || 0;
            bSKINITandDEV = (CPUInfo[2] & 0x1000) || 0;
            bSYSCALL_SYSRETAvailable = (CPUInfo[3] & 0x800) || 0;
            bExecuteDisableBitAvailable = (CPUInfo[3] & 0x10000) || 0;
            bMMXExtensions = (CPUInfo[3] & 0x40000) || 0;
            bFFXSR = (CPUInfo[3] & 0x200000) || 0;
            b1GBSupport = (CPUInfo[3] & 0x400000) || 0;
            bRDTSCP = (CPUInfo[3] & 0x8000000) || 0;
            b64Available = (CPUInfo[3] & 0x20000000) || 0;
            b3DNowExt = (CPUInfo[3] & 0x40000000) || 0;
            b3DNow = (CPUInfo[3] & 0x80000000) || 0;
        }

        // Interpret CPU brand string and cache information.
        if  (i == 0x80000002)
            memcpy(CPUBrandString, CPUInfo, sizeof(CPUInfo));
        else if  (i == 0x80000003)
            memcpy(CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
        else if  (i == 0x80000004)
            memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));
        else if  (i == 0x80000006)
        {
            nCacheLineSize = CPUInfo[2] & 0xff;
            nL2Associativity = (CPUInfo[2] >> 12) & 0xf;
            nCacheSizeK = (CPUInfo[2] >> 16) & 0xffff;
        }
        else if  (i == 0x80000008)
        {
           nPhysicalAddress = CPUInfo[0] & 0xff;
           nVirtualAddress = (CPUInfo[0] >> 8) & 0xff;
        }
        else if  (i == 0x8000000A)
        {
            bNestedPaging = (CPUInfo[3] & 0x1) || 0;
            bLBRVisualization = (CPUInfo[3] & 0x2) || 0;
        }
        else if  (i == 0x8000001A)
        {
            bFP128 = (CPUInfo[0] & 0x1) || 0;
            bMOVOptimization = (CPUInfo[0] & 0x2) || 0;
        }
    }
	strcpy(idstr,CPUBrandString);
    // Display all the information in user-friendly format.
	if (do_print)
		printf("\n\nCPU String: %s\n", CPUString);

    if  (nIds >= 1)
    {
		if (do_print)
		{
			if  (nSteppingID)
				printf("Stepping ID = %d\n", nSteppingID);
			if  (nModel)
				printf("Model = %d\n", nModel);
			if  (nFamily)
				printf("Family = %d\n", nFamily);
			if  (nProcessorType)
				printf("Processor Type = %d\n", nProcessorType);
			if  (nExtendedmodel)
				printf("Extended model = %d\n", nExtendedmodel);
			if  (nExtendedfamily)
				printf("Extended family = %d\n", nExtendedfamily);
			if  (nBrandIndex)
				printf("Brand Index = %d\n", nBrandIndex);
			if  (nCLFLUSHcachelinesize)
				printf("CLFLUSH cache line size = %d\n",
						 nCLFLUSHcachelinesize);
			if (bMultithreading && (nLogicalProcessors > 0))
			   printf("Logical Processor Count = %d\n", nLogicalProcessors);
			if  (nAPICPhysicalID)
				printf("APIC Physical ID = %d\n", nAPICPhysicalID);

			if  (nFeatureInfo || bSSE3Instructions ||
				 bMONITOR_MWAIT || bCPLQualifiedDebugStore ||
				 bVirtualMachineExtensions || bEnhancedIntelSpeedStepTechnology ||
				 bThermalMonitor2 || bSupplementalSSE3 || bL1ContextID || 
				 bCMPXCHG16B || bxTPRUpdateControl || bPerfDebugCapabilityMSR || 
				 *bSSE41Extensions || bSSE42Extensions || bPOPCNT || 
				 bLAHF_SAHFAvailable || bCmpLegacy || bSVM ||
				 bExtApicSpace || bAltMovCr8 ||
				 bLZCNT || bSSE4A || bMisalignedSSE ||
				 bPREFETCH || bSKINITandDEV || bSYSCALL_SYSRETAvailable || 
				 bExecuteDisableBitAvailable || bMMXExtensions || bFFXSR || b1GBSupport ||
				 bRDTSCP || b64Available || b3DNowExt || b3DNow || bNestedPaging || 
				 bLBRVisualization || bFP128 || bMOVOptimization )
			{
				printf("\nThe following features are supported:\n");

				if  (bSSE3Instructions)
					printf("\tSSE3\n");
				if  (bMONITOR_MWAIT)
					printf("\tMONITOR/MWAIT\n");
				if  (bCPLQualifiedDebugStore)
					printf("\tCPL Qualified Debug Store\n");
				if  (bVirtualMachineExtensions)
					printf("\tVirtual Machine Extensions\n");
				if  (bEnhancedIntelSpeedStepTechnology)
					printf("\tEnhanced Intel SpeedStep Technology\n");
				if  (bThermalMonitor2)
					printf("\tThermal Monitor 2\n");
				if  (bSupplementalSSE3)
					printf("\tSupplemental Streaming SIMD Extensions 3\n");
				if  (bL1ContextID)
					printf("\tL1 Context ID\n");
				if  (bCMPXCHG16B)
					printf("\tCMPXCHG16B Instruction\n");
				if  (bxTPRUpdateControl)
					printf("\txTPR Update Control\n");
				if  (bPerfDebugCapabilityMSR)
					printf("\tPerf\\Debug Capability MSR\n");
				if  (*bSSE41Extensions)
					printf("\tSSE4.1 Extensions\n");
				if  (bSSE42Extensions)
					printf("\tSSE4.2 Extensions\n");
				if  (bPOPCNT)
					printf("\tPPOPCNT Instruction\n");

				i = 0;
				nIds = 1;
				while (i < (sizeof(szFeatures)/sizeof(const char*)))
				{
					if  (nFeatureInfo & nIds)
					{
						printf("\t");
						printf("%s", szFeatures[i]);
						printf("\n");
					}

					nIds <<= 1;
					++i;
				}
				if (bLAHF_SAHFAvailable)
					printf("\tLAHF/SAHF in 64-bit mode\n");
				if (bCmpLegacy)
					printf("\tCore multi-processing legacy mode\n");
				if (bSVM)
					printf("\tSecure Virtual Machine\n");
				if (bExtApicSpace)
					printf("\tExtended APIC Register Space\n");
				if (bAltMovCr8)
					printf("\tAltMovCr8\n");
				if (bLZCNT)
					printf("\tLZCNT instruction\n");
				if (bSSE4A)
					printf("\tSSE4A (EXTRQ, INSERTQ, MOVNTSD, MOVNTSS)\n");
				if (bMisalignedSSE)
					printf("\tMisaligned SSE mode\n");
				if (bPREFETCH)
					printf("\tPREFETCH and PREFETCHW Instructions\n");
				if (bSKINITandDEV)
					printf("\tSKINIT and DEV support\n");
				if (bSYSCALL_SYSRETAvailable)
					printf("\tSYSCALL/SYSRET in 64-bit mode\n");
				if (bExecuteDisableBitAvailable)
					printf("\tExecute Disable Bit\n");
				if (bMMXExtensions)
					printf("\tExtensions to MMX Instructions\n");
				if (bFFXSR)
					printf("\tFFXSR\n");
				if (b1GBSupport)
					printf("\t1GB page support\n");
				if (bRDTSCP)
					printf("\tRDTSCP instruction\n");
				if (b64Available)
					printf("\t64 bit Technology\n");
				if (b3DNowExt)
					printf("\t3Dnow Ext\n");
				if (b3DNow)
					printf("\t3Dnow! instructions\n");
				if (bNestedPaging)
					printf("\tNested Paging\n");
				if (bLBRVisualization)
					printf("\tLBR Visualization\n");
				if (bFP128)
					printf("\tFP128 optimization\n");
				if (bMOVOptimization)
					printf("\tMOVU Optimization\n");
			}
		}		 
    }

    if  (nExIds >= 0x80000004 && do_print)
        printf("\nCPU Brand String: %s\n", CPUBrandString);

    if  (nExIds >= 0x80000006 && do_print)
    {
        printf("Cache Line Size = %d\n", nCacheLineSize);
        printf("L2 Associativity = %d\n", nL2Associativity);
        printf("Cache Size = %dK\n", nCacheSizeK);
    }


    for (i=0;;i++)
    {
		CPUID2(0x4,i,CPUInfo[0],CPUInfo[1],CPUInfo[2],CPUInfo[3]);
        //__cpuidex(CPUInfo, 0x4, i);
        if(!(CPUInfo[0] & 0xf0)) break;

        if(i == 0)
        {
            nCores = CPUInfo[0] >> 26;
			if (do_print)
				printf("\n\nNumber of Cores = %d\n", nCores + 1);
        }

        nCacheType = (CPUInfo[0] & 0x1f);
        nCacheLevel = (CPUInfo[0] & 0xe0) >> 5;
        bSelfInit = (CPUInfo[0] & 0x100) >> 8;
        bFullyAssociative = (CPUInfo[0] & 0x200) >> 9;
        nMaxThread = (CPUInfo[0] & 0x03ffc000) >> 14;
        nSysLineSize = (CPUInfo[1] & 0x0fff);
        nPhysicalLinePartitions = (CPUInfo[1] & 0x03ff000) >> 12;
        nWaysAssociativity = (CPUInfo[1]) >> 22;
        nNumberSets = CPUInfo[2];

		if (do_print)
		{
			printf("\n");

			printf("ECX Index %d\n", i);
			switch (nCacheType)
			{
				case 0:
					printf("   Type: Null\n");
					break;
				case 1:
					printf("   Type: Data Cache\n");
					break;
				case 2:
					printf("   Type: Instruction Cache\n");
					break;
				case 3:
					printf("   Type: Unified Cache\n");
					break;
				default:
					 printf("   Type: Unknown\n");
			}

			printf("   Level = %d\n", nCacheLevel + 1); 
			if (bSelfInit)
			{
				printf("   Self Initializing\n");
			}
			else
			{
				printf("   Not Self Initializing\n");
			}
			if (bFullyAssociative)
			{
				printf("   Is Fully Associatve\n");
			}
			else
			{
				printf("   Is Not Fully Associatve\n");
			}
			printf("   Max Threads = %d\n", 
				nMaxThread+1);
			printf("   System Line Size = %d\n", 
				nSysLineSize+1);
			printf("   Physical Line Partions = %d\n", 
				nPhysicalLinePartitions+1);
			printf("   Ways of Associativity = %d\n", 
				nWaysAssociativity+1);
			printf("   Number of Sets = %d\n", 
				nNumberSets+1);
		}
    }

    return  nRet;
}

int bin_search_uint32(int idp, int idm, uint32 q, uint32 *input)
{
	int next = (idp + idm) / 2;
	
	while ((idp - idm) > 10)
	{
		if (input[next] > q)
		{
			idp = next;
			next = (next + idm) / 2;							
		}
		else					
		{
			idm = next;
			next = (idp + next) / 2;							
		}
	}

	for (next = idm; next < idm + 10; next++)
		if (input[next] == q)
			return next;

	if (input[next] != q)
		next = -1;

	return next;
}


