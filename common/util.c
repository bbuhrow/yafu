/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: util.c 964 2014-05-03 03:30:03Z jasonp_sf $
--------------------------------------------------------------------*/

#include <util.h>

/*---------------------------------------------------------------------*/
void *
aligned_malloc(size_t len, uint32 align) {

	/* this code internally calls malloc() with a size > len,
	   rounds the returned pointer up to the next "align"-byte
	   boundary, stores the original pointer in the bytes just
	   before that address and returns a pointer to aligned
	   memory. aligned_free() reverses the process so that the
	   padded address returned here is automatically accounted for.
	
	   The arithmetic is messy because there's no guarantee that
	   type "void *" is the same size as type "unsigned long". For
	   example, on the Alpha you can specify a long to be 4 bytes
	   but pointers are all 8 bytes in size. */
      
	void *ptr, *aligned_ptr;
	unsigned long addr;

	ptr = xmalloc(len+align);

	 /* offset to next ALIGN-byte boundary */

	addr = (unsigned long)ptr;				
	addr = align - (addr & (align-1));
	aligned_ptr = (void *)((uint8 *)ptr + addr);

	*( (void **)aligned_ptr - 1 ) = ptr;
	return aligned_ptr;
}

/*---------------------------------------------------------------------*/
void
aligned_free(void *newptr) {

	void *ptr;

	if (newptr == NULL) 
		return;
	ptr = *( (void **)newptr - 1 );
	free(ptr);
}

/*------------------------------------------------------------------*/
uint64
read_clock(void) {

#if defined(GCC_ASM32X) || defined(GCC_ASM64X) 
	uint32 lo, hi;
	ASM_G("rdtsc":"=d"(hi),"=a"(lo));
	return (uint64)hi << 32 | lo;

#elif defined(_MSC_VER)
	LARGE_INTEGER ret;
	QueryPerformanceCounter(&ret);
	return ret.QuadPart;

#else
	struct timeval thistime;   
	gettimeofday(&thistime, NULL);
	return thistime.tv_sec * 1000000 + thistime.tv_usec;
#endif
}

/*------------------------------------------------------------------*/
double
get_cpu_time(void) {

#if defined(WIN32) || defined(_WIN64)
	FILETIME create_time = {0, 0};
	FILETIME exit_time = {0, 0};
	FILETIME kernel_time = {0, 0};
	FILETIME user_time = {0, 0};

	GetThreadTimes(GetCurrentThread(),
			&create_time,
			&exit_time,
			&kernel_time,
			&user_time);

	return ((uint64)user_time.dwHighDateTime << 32 | 
	               user_time.dwLowDateTime) / 10000000.0;
#else
	struct rusage r_usage;

	#if 0 /* use for linux 2.6.26+ */
	getrusage(RUSAGE_THREAD, &r_usage);
	#else
	getrusage(RUSAGE_SELF, &r_usage);
	#endif

	return ((uint64)r_usage.ru_utime.tv_sec * 1000000 +
	               r_usage.ru_utime.tv_usec) / 1000000.0;
#endif
}


#if defined(_MSC_VER)

#if defined(__clang__)
int msieve_gettimeofday(struct timeval* tv, struct timezone* tz)
{
	struct timespec ts;
	timespec_get(&ts, TIME_UTC);

	//printf("timespec_get returned sec = %lu, nsec = %lu\n", ts.tv_sec, ts.tv_nsec);

	tv->tv_sec = ts.tv_sec;
	tv->tv_usec = ts.tv_nsec / 1000;

	return 0;
}
#else
int msieve_gettimeofday(struct timeval* tv, struct timezone* tz)
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
#endif

double msieve_difftime(struct timeval* start, struct timeval* end)
{
	double secs;
	double usecs;

	if (start->tv_sec == end->tv_sec) {
		secs = 0;
		usecs = end->tv_usec - start->tv_usec;
	}
	else {
		usecs = 1000000 - start->tv_usec;
		secs = end->tv_sec - (start->tv_sec + 1);
		usecs += end->tv_usec;
		if (usecs >= 1000000) {
			usecs -= 1000000;
			secs += 1;
		}
	}

	return secs + usecs / 1000000.;
}

/*--------------------------------------------------------------------*/
void set_idle_priority(void) {

#if (defined(WIN32) || defined(_WIN64))
	SetPriorityClass(GetCurrentProcess(),
			IDLE_PRIORITY_CLASS);
#else
	nice(100);
#endif
}

/*--------------------------------------------------------------------*/
/* safe default values */
#define DEFAULT_L1_CACHE_SIZE (32 * 1024)
#define DEFAULT_L2_CACHE_SIZE (512 * 1024)

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

/* macro to execute the x86 CPUID instruction. Note that
   this is more verbose than it needs to be; Intel Macs reserve
   the EBX or RBX register for the PIC base address, and so
   this register cannot get clobbered by inline assembly */

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

#elif defined(_MSC_VER) && defined(__clang__)
#include <x86intrin.h>
#define HAS_CPUID
#define CPUID(__leaf, __eax, __ebx, __ecx, __edx) \
    __asm("cpuid" : "=a"(__eax), "=b" (__ebx), "=c"(__ecx), "=d"(__edx) \
                  : "0"(__leaf))
#define CPUID2(code1, code2, a, b, c, d) \
	__asm("cpuid" : "=a"(a), "=b" (b), "=c"(c), "=d"(d) \
                  : "0"(code1), "2"(code2))

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

void get_cache_sizes(uint32 *level1_size_out,
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

#if defined(HAS_CPUID)

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
enum cpu_type get_cpu_type(void) {

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

/*--------------------------------------------------------------------*/
uint64 get_file_size(char *name) {

#if defined(WIN32) || defined(_WIN64)
	struct _stati64 tmp;

	if (_stati64(name, &tmp) != 0) {
		char name_gz[256];
		sprintf(name_gz, "%s.gz", name);
		if (_stati64(name_gz, &tmp) != 0) 
			return 0;
		return (tmp.st_size / 11) * 20;
	}
#else
	struct stat tmp;

	if (stat(name, &tmp) != 0) {
		char name_gz[256];
		sprintf(name_gz, "%s.gz", name);
		if (stat(name_gz, &tmp) != 0) 
			return 0;
		return (tmp.st_size / 11) * 20;
	}
#endif

	return tmp.st_size;
}

/*--------------------------------------------------------------------*/
uint64 get_ram_size(void) {

#if defined(WIN32)
	MEMORYSTATUS tmp;

	tmp.dwLength = sizeof(MEMORYSTATUS);
	GlobalMemoryStatus(&tmp);

	return tmp.dwTotalPhys;

#elif defined(_WIN64)
	MEMORYSTATUSEX tmp;

	tmp.dwLength = sizeof(MEMORYSTATUSEX);
	if (GlobalMemoryStatusEx(&tmp) == FALSE)
		return 0;

	return tmp.ullTotalPhys;

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)
	int page_size = sysconf(_SC_PAGESIZE);
	int num_pages = sysconf(_SC_PHYS_PAGES);

	if (page_size < 0 || num_pages < 0)
		return 0;
	return (uint64)page_size * (uint64)num_pages;
#else

	return 0;
#endif
}

/*--------------------------------------------------------------------*/
libhandle_t load_dynamic_lib(const char *libname)
{
#if defined(WIN32) || defined(_WIN64)
	HMODULE h = LoadLibraryA((LPCSTR)libname);

	if (h == NULL)
		printf("cannot load library '%s', error %u\n", 
				libname, (uint32)GetLastError());
#else
	void * h = dlopen(libname, RTLD_LAZY);

	if (h == NULL)
		printf("cannot load library '%s': %s\n", 
				libname, dlerror());
#endif
	return h;
}

/*--------------------------------------------------------------------*/
void unload_dynamic_lib(libhandle_t h)
{
#if defined(WIN32) || defined(_WIN64)
	FreeLibrary(h);
#else
	dlclose(h);
#endif
}

/*--------------------------------------------------------------------*/
void * get_lib_symbol(libhandle_t h, const char *symbol_name)
{
#if defined(WIN32) || defined(_WIN64)
	void *s = GetProcAddress(h, (LPCSTR)symbol_name);

	if (s == NULL)
		printf("cannot load symbol '%s', error %u\n", 
				symbol_name, (uint32)GetLastError());
#else
	void * s = dlsym(h, symbol_name);

	if (s == NULL)
		printf("cannot load symbol '%s': %s\n", 
				symbol_name, dlerror());
#endif
	return s;
}
