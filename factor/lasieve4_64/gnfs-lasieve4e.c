/* gnfs-lasieve4e.c
  By Jens Franke.
  6/13/04: Hacked up for use in GGNFS by Chris Monico.
  9/30/22: Vector AVX512 code contributed by Ben Buhrow

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/


// #include "lasieve.h"

#include <stdio.h> 
#include <sys/types.h> 
#include <math.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <limits.h> 
#include <string.h> 
#include <time.h> 
#include <ctype.h>
#include <immintrin.h>

#include "../gmp-ecm/microecm.h"


#ifdef LINUX
#include <endian.h> 
#endif
#include <gmp.h> 
#include <signal.h> 
#include <setjmp.h> 
#if defined( _MSC_VER ) && defined( _WIN64 )
#include <crtdbg.h>
#endif
#include "asm/siever-config.h"
#ifndef TDS_MPQS
#define TDS_MPQS TDS_SPECIAL_Q
#endif
#ifndef TDS_PRIMALITY_TEST
#define TDS_PRIMALITY_TEST TDS_IMMEDIATELY
#endif

#include "if.h"
#include "primgen32.h"
#include "asm/32bit.h"
#include "redu2.h"
#include "recurrence6.h"
#include "fbgen.h"
#include "real-poly-aux.h"
#include "gmp-aux.h"
#include "lasieve-prepn.h"
#include <immintrin.h>

#define TDS_IMMEDIATELY 0
#define TDS_BIGSS 1
#define TDS_ODDNESS_CLASS 2
#define TDS_SPECIAL_Q 3

#define GCD_SIEVE_BOUND 10
#include "asm/siever-config.c"

#include "asm/lasched.h"
#include "asm/medsched.h"

#define L1_SIZE (1UL<<L1_BITS)

#if 0
#define ZSS_STAT
u32_t nss= 0,nzss[3]= {0,0,0};
#endif

//#define AVX512_SIEVEA
//#define AVX512_SIEVE4
//#define AVX512_SIEVE3
//#define AVX512_SIEVE2
//#define AVX512_SIEVE1

#ifdef _MSC_VER
#define AVX512_SIEVE1
#endif

#ifdef AVX512_SIEVE4
uint64_t*** rtab_s0;
uint64_t*** rtab_s1;

static int presieve_bound_p = 31;
static int rtab_lookup[32] = { 0,0,0,0,0,1,0,2,0,0,0,3,0,4,0,0,0,5,0,6,0,0,0,7,0,0,0,0,0,8,0,9};
static int steps[32] = { 0,0,0,63,0,60,0,63,0,0,0,55,0,52,0,0,0,51,0,57,0,0,0,46,0,0,0,0,0,58,0,62};
static int initialized_rtab_s1[10] = {0,0,0,0,0,0,0,0,0,0};
static int initialized_rtab_s0[10] = {0,0,0,0,0,0,0,0,0,0};
static int all_initialized_s[2] = { 0,0 };
#endif

static float FB_bound[2],sieve_report_multiplier[2];
static u16_t sieve_min[2],max_primebits[2],max_factorbits[2];
static u32_t*(FB[2]),*(proots[2]),FBsize[2];

static double*(tpoly_f[2]);
#define CANDIDATE_SEARCH_STEPS 128
static unsigned char**(sieve_report_bounds[2]);
static i32_t n_srb_i,n_srb_j;

static u32_t b,first_spq,first_spq1,first_root,last_spq,sieve_count;

static mpz_t m,N,aux1,aux2,aux3,sr_a,sr_b;

static mpz_t*(poly[2]);

// for tinyecm/microecm
static uint64_t pran;
static mpz_t factor1, factor2, factor3;

double*(poly_f[2]),poly_norm[2];

i32_t poldeg[2],poldeg_max;

u32_t  keep_factorbase;
u32_t  g_resume;
#define MAX_LPFACTORS 3
static mpz_t rational_rest,algebraic_rest;
mpz_t factors[MAX_LPFACTORS];
static u32_t yield= 0,n_mpqsfail[2]= {0,0},n_mpqsvain[2]= {0,0};
static i64_t mpqs_clock= 0;

static i64_t sieve_clock= 0, lsetup_clock=0, sch_clock= 0,td_clock= 0,tdi_clock= 0;
static i64_t cs_clock[2]= {0,0},Schedule_clock= 0,medsched_clock= 0;
static i64_t si_clock[2]= {0,0},s1_clock[2]= {0,0};
static i64_t s2_clock[2]= {0,0},s3_clock[2]= {0,0};
static i64_t tdsi_clock[2]= {0,0},tds1_clock[2]= {0,0},tds2_clock[2]= {0,0};
static i64_t tds3_clock[2]= {0,0},tds4_clock[2]= {0,0};

char *base_name;
char *input_line= NULL;
size_t input_line_alloc= 0;


static u32_t ncand;
static u16_t*cand;
static unsigned char*fss_sv;
u32_t process_no;
char *sysload_cmd;
double sieveStartTime;

static int tdcand_cmp(const void*x,const void*y)
{
  return(int)(*((u16_t*)x))-(int)(*((u16_t*)y));
}

typedef struct xFBstruct{
  u32_t p,pp,q,qq,r,l;
}*xFBptr;
static volatile xFBptr xFB[2];
static volatile u32_t xFBs[2];

static void xFBtranslate(u16_t*rop,xFBptr op);
static int xFBcmp(const void*,const void*);

static u32_t add_primepowers2xaFB(size_t*aFB_alloc_ptr,u32_t pp_bound,u32_t side,u32_t p,u32_t r);



int getfactor_tecm(mpz_t n, mpz_t f, int target_bits, uint64_t* pran);




i32_t a0,a1,b0,b1;
#if 0
u32_t I_bits;
#endif
u32_t J_bits,i_shift,n_I,n_J;
u32_t root_no;
float sigma;

static u32_t oddness_type;
static u32_t n_i,n_j,i_bits,j_bits;

u32_t spq_i,spq_j,spq_x;

u32_t fbi1[2];

u32_t fbis[2];

#if I_bits<=L1_BITS
static u32_t j_per_strip,jps_bits;
#else
#define j_per_strip 1
#define jps_bits    0
#endif
static u32_t n_strips;
static struct schedule_struct{
  u16_t***schedule;
  u32_t*fbi_bounds;
  u32_t n_pieces;
  unsigned char*schedlogs;
  u16_t n_strips,current_strip;
  size_t alloc,alloc1;
  u32_t*ri;
}*(schedules[2]);
u32_t n_schedules[2];

static u32_t*(LPri[2]);
#define RI_SIZE 2

static u32_t*(current_ij[2]);

static size_t sched_alloc[2];
#define SE_SIZE 2
#define SCHEDFBI_MAXSTEP 0x10000

#define USE_MEDSCHED
#ifdef USE_MEDSCHED
static u16_t**(med_sched[2]);
static u32_t*(medsched_fbi_bounds[2]);
static unsigned char*(medsched_logs[2]);
static size_t medsched_alloc[2];
static u16_t n_medsched_pieces[2];
#endif

static unsigned char*sieve_interval= NULL,*(FB_logs[2]);
static unsigned char*tiny_sieve_buffer;
#define TINY_SIEVEBUFFER_SIZE 420
#define TINY_SIEVE_MIN 8
static double sieve_multiplier[2],FB_maxlog[2];
static u32_t j_offset;

void do_scheduling(struct schedule_struct*,u32_t,u32_t,u32_t);

static u16_t*(smallsieve_aux[2]),*(smallsieve_auxbound[2][5]);
static u16_t*(smallsieve_tinybound[2]);

static u16_t*(smallsieve_aux1[2]),*(smallsieve_aux1_ub_odd[2]);
static u16_t*(smallsieve_aux1_ub[2]),*(smallsieve_tinybound1[2]);

static u16_t*(smallsieve_aux2[2]),*(smallsieve_aux2_ub[2]);

static u16_t*(smallpsieve_aux[2]),*(smallpsieve_aux_ub_pow1[2]);
static u16_t*(smallpsieve_aux_ub_odd[2]),*(smallpsieve_aux_ub[2]);
static unsigned char*horizontal_sievesums;

static u16_t*(x2FB[2]),x2FBs[2];

static u16_t*tinysieve_curpos;
#ifndef MMX_TD
static u16_t**(smalltdsieve_aux[2]);
#ifdef PREINVERT
static u32_t*(smalltd_pi[2]);
#endif
#endif

#ifdef GCD_SIEVE_BOUND
static u32_t np_gcd_sieve;
static unsigned char*gcd_sieve_buffer;
static void gcd_sieve(void);
#endif

u16_t**schedbuf;

static void store_candidate(u16_t,u16_t,unsigned char);

void trial_divide(void);

#ifndef SCHED_TDS_BUFSIZE
#define SCHED_TDS_BUFSIZE 1024
#endif
u16_t*(sched_tds_buffer[SCHED_TDS_BUFSIZE]);

u32_t*mpz_trialdiv(mpz_t N,u32_t*pbuf,u32_t ncp,char*errmsg);

static void output_tdsurvivor(u32_t*,u32_t*,u32_t*,u32_t*,mpz_t,mpz_t);
static void store_tdsurvivor(u32_t*,u32_t*,u32_t*,u32_t*,mpz_t,mpz_t);
static int primality_tests(void);
static void primality_tests_all(void);
static void output_all_tdsurvivors(void);
static u32_t*tds_fbp_buffer;
static i64_t*tds_ab;
static mpz_t*tds_lp;
static size_t max_tds= 0,*tds_fbp,tds_fbp_alloc= 0,total_ntds= 0;
#define MAX_TDS_INCREMENT 1024
#define TDS_FBP_ALLOC_INCREMENT 8192

#if 0
#define OFMT_CWI
#endif
#ifdef OFMT_CWI
static char u32_t2cwi(u32_t);
#endif

void dumpsieve(u32_t j_offset,u32_t side);

u32_t*(td_buf[2]),**td_buf1;
size_t td_buf_alloc[2]= {1024,1024};

static unsigned char tds_coll[UCHAR_MAX];
u32_t**tds_fbi= NULL;
u32_t**tds_fbi_curpos= NULL;
#ifndef TDFBI_ALLOC
#define TDFBI_ALLOC 256
static size_t tds_fbi_alloc= TDFBI_ALLOC;
#endif

static mpz_t td_rests[L1_SIZE];
static mpz_t large_factors[2],*(large_primes[2]);
static mpz_t FBb_sq[2];
static mpz_t FBb_cu[2];


void Usage()
{
  complain("Usage");
}

static u32_t n_prereports= 0,n_reports= 0,n_rep1= 0,n_rep2= 0;
static u32_t n_tdsurvivors[2]= {0,0};
static FILE*g_ofile;
static char*g_ofile_name;

#ifdef STC_DEBUG
FILE*debugfile;
#endif

static u16_t special_q_side,first_td_side,first_sieve_side;
static u16_t append_output,exitval;
static u16_t cmdline_first_sieve_side= USHRT_MAX;
static u16_t cmdline_first_td_side= USHRT_MAX;
static u16_t cmdline_first_mpqs_side= USHRT_MAX;
static u16_t cmdline_first_psp_side= USHRT_MAX;
#define ALGEBRAIC_SIDE 0
#define RATIONAL_SIDE 1
#define NO_SIDE 2

static pr32_struct special_q_ps;
u32_t special_q;
double special_q_log;

#define USER_INTERRUPT 1
#define SCHED_PATHOLOGY 2

#define USER_INTERRUPT_EXITVAL 2
#define LOADTEST_EXITVAL 3

jmp_buf termination_jb;

static void
terminate_sieving(int signo)
{
  exitval= USER_INTERRUPT_EXITVAL;
  longjmp(termination_jb,USER_INTERRUPT);
}

static clock_t last_clock;

#ifdef MMX_TDBENCH
extern u64_t MMX_TdNloop;
#endif

/*******************************************************/
double sTime()
     /*******************************************************/
#if 0 && !defined (_MSC_VER) && !defined (__MINGW32__) && !defined (MINGW32)
{
  static struct  timeval  this_tv;
  static struct  timezone dumbTZ;
  double t;
  
  gettimeofday(&this_tv, &dumbTZ);
  t = this_tv.tv_sec + 0.000001*this_tv.tv_usec;
  return t;
}
#else
{
  return clock() / (double)CLOCKS_PER_SEC;
}
#endif

/**************************************************/
void logTotalTime()
     /**************************************************/
{
  double t=sTime()-sieveStartTime;
  FILE *fp=fopen("ggnfs.log", "a");
  
  fprintf(fp, "\tLatSieveTime: %ld\n", (long)t);
  fclose(fp);
}

/**************************************************/
int parse_q_from_line(char *buf) {
  /**************************************************/
  char *p, *tmp, *next_field;
  u32_t q, q0, i, side;
  static int first=0;
  
  for(p=tmp=buf; *p && isspace(*p); p++);
  if(!*p) return 0; /* empty line, skip */
  
  side = (special_q_side == RATIONAL_SIDE) ? 0 : 1;
  for(i=0; *p; p++) {
    if(*p == ':') {
      if(i++ == side) tmp = p; /* we will only scan this section for a q0 */
    } else if(!(*p == '-' || *p == ',' || isspace(*p) || isxdigit(*p))) {
      if(first++ == 0) printf(" Warning! some corrupt lines in the original file\n");
      return -1;
    }
  }
  if(i!=2) {
    printf(" Warning: an incomplete line in the original file; if just a few, it's ok, they will be skipped\n");
    return -1;           /* must have two ':' some ',' and hexdigits */
  }
  
  q0 = first_spq;
  do {
    q = strtoul(tmp + 1, &next_field, 16);
    if(q >= first_spq && q < first_spq+sieve_count)
      q0 = q;
    tmp = next_field;
  } while(tmp[0] == ',' && isxdigit(tmp[1]));
  
  /* I've seen cases when q0 is not the last reported in the comma-separated list */
  /* However, the closer it is to the end of the line the more likely it was the true q0 */
  /* In 99% cases it is the last value, but we don't want to depend on that */
  
  if(q0 > first_spq && q0 < first_spq+sieve_count) {
    sieve_count -= (q0 - first_spq);
    first_spq = q0;
  }
  return 1;
}  

/**************************************************/
int main(int argc, char **argv)
     /**************************************************/
{
  u16_t zip_output,force_aFBcalc;
  u16_t catch_signals;
  u32_t all_spq_done;
  u32_t n_spq, n_spq_discard;
  double tStart, tNow, lastReport;
  
#if defined (_MSC_VER) && defined (_DEBUG)
  int tmpDbgFlag;
  tmpDbgFlag = _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG);
  /*
    tmpDbgFlag |= _CRTDBG_CHECK_ALWAYS_DF;
    tmpDbgFlag |= _CRTDBG_CHECK_CRT_DF;
    tmpDbgFlag |= _CRTDBG_DELAY_FREE_MEM_DF;
  */
  tmpDbgFlag |= _CRTDBG_LEAK_CHECK_DF;
  _CrtSetDbgFlag(tmpDbgFlag);
#endif
  
  n_spq = 0;
  n_spq_discard = 0;
  
  sieveStartTime = sTime();
  
#ifdef STC_DEBUG
  debugfile= fopen("rtdsdebug","w");
#endif


  {
	  i32_t option;
	  FILE* input_data;
	  u32_t i;

	  sieve_count = U32_MAX;
	  g_ofile_name = NULL;
	  zip_output = 0;
	  special_q_side = NO_SIDE;
	  sigma = 0;
	  keep_factorbase = 0;
	  g_resume = 0;
	  base_name = NULL;
	  first_spq = 0;
	  sieve_count = 1;
	  force_aFBcalc = 0;
	  sysload_cmd = NULL;
	  process_no = 0;
	  catch_signals = 0;

	  J_bits = U32_MAX;

#define NumRead(x) if(sscanf(optarg,"%u",&x)!=1) Usage()
#define NumRead16(x) if(sscanf(optarg,"%hu",&x)!=1) Usage()

	  // parse command line arguments
	  append_output = 0;
	  while ((option = getopt(argc, argv, "FJ:L:M:N:P:RS:ab:c:f:i:kn:o:q:rt:vz")) != -1) {
		  switch (option) {
		  case 'R':
			  g_resume = 1; break;
		  case'F':
			  force_aFBcalc = 1;
			  break;
		  case'J':
			  NumRead(J_bits);
			  break;
		  case'L':
			  sysload_cmd = optarg;
			  break;
		  case'M':
			  if (sscanf(optarg, "%hu", &cmdline_first_mpqs_side) != 1)
				  complain("-M %s ???\n", optarg);
			  break;
		  case'P':
			  if (sscanf(optarg, "%hu", &cmdline_first_psp_side) != 1)
				  complain("-P %s ???\n", optarg);
			  break;
		  case'S':
			  if (sscanf(optarg, "%f", &sigma) != 1) {
				  errprintf("Cannot read floating point number %s\n", optarg);
				  Usage();
			  }
			  break;
		  case'a':
			  if (special_q_side != NO_SIDE) {
				  errprintf("Ignoring -a\n");
				  break;
			  }
			  special_q_side = ALGEBRAIC_SIDE;
			  break;
		  case'b':
			  if (base_name != NULL)errprintf("Ignoring -b %s\n", base_name);
			  else base_name = optarg;
			  break;
		  case'c':
			  NumRead(sieve_count);
			  break;
		  case'f':
			  if (sscanf(optarg, "%u:%u:%u", &first_spq, &first_spq1, &first_root) != 3) {
				  if (sscanf(optarg, "%u", &first_spq) == 1) {
					  first_spq1 = first_spq;
					  first_root = 0;
				  }
				  else Usage();
			  }
			  else append_output = 1;
			  break;
		  case'i':
			  if (sscanf(optarg, "%hu", &cmdline_first_sieve_side) != 1)
				  complain("-i %s ???\n", optarg);
			  break;
		  case'k':
			  keep_factorbase = 1;
			  break;
		  case'n':
			  catch_signals = 1;
		  case'N':
			  NumRead(process_no);
			  break;
		  case'o':
			  g_ofile_name = optarg;
			  break;
		  case'q':
			  NumRead16(special_q_side);
			  break;
		  case'r':
			  if (special_q_side != NO_SIDE) {
				  errprintf("Ignoring -r\n");
				  break;
			  }
			  special_q_side = RATIONAL_SIDE;
			  break;
		  case't':
			  if (sscanf(optarg, "%hu", &cmdline_first_td_side) != 1)
				  complain("-t %s ???\n", optarg);
			  break;
		  case'v':
			  verbose++;
			  break;
		  case'z':
			  zip_output = 1;
			  break;
		  }
	  }

	  char features[1024];
	  sprintf(features, "with asm64");
#ifdef AVX512_TD
	  sprintf(features, "%s,avx-512 mmx-td", features);
#endif
#ifdef AVX512_LASIEVE_SETUP
	  sprintf(features, "%s,avx-512 lasetup", features);
#endif
#ifdef AVX512_LASCHED
	  sprintf(features, "%s,avx-512 lasched", features);
#endif
#ifdef AVX512_SIEVE1
	  sprintf(features, "%s,avx-512 sieve1", features);
#endif
#ifdef AVX512_ECM
	  sprintf(features, "%s,avx-512 ecm", features);
#endif
	  if (verbose) { /* first rudimentary test of automatic $Rev reporting */
		  fprintf(stderr, "gnfs-lasieve4I%de (%s): L1_BITS=%d\n", 
			  I_bits, features, L1_BITS);
	  }

#define LINE_BUF_SIZE 300

	  if (g_resume != 0) {
		  char buf[LINE_BUF_SIZE];
		  int ret;

		  if (zip_output != 0)
			  complain("Cannot resume gzipped file. gunzip, and retry without -z\n");
		  if (g_ofile_name == NULL)
			  complain("Cannot resume without the file name\n");
		  if (strcmp(g_ofile_name, "-") == 0)
			  complain("Cannot resume using stdout\n");
		  if ((g_ofile = fopen(g_ofile_name, "ab+")) == NULL)
			  complain("Cannot open %s for append: %m\n", g_ofile_name);
		  while (fgets(buf, LINE_BUF_SIZE, g_ofile)) {
			  ret = parse_q_from_line(buf);
		  }
		  if (ret < 0) fprintf(g_ofile, "\n"); /* encapsulating the last incomplete line */
		  printf(" Resuming with -f %d -c %d\n", first_spq, sieve_count);
		  first_spq1 = first_spq;
	  }

	  if (J_bits == U32_MAX)
		  J_bits = I_bits - 1;
	  if (cmdline_first_psp_side == USHRT_MAX)
		  cmdline_first_psp_side = cmdline_first_mpqs_side;
#ifndef I_bits
#error Must #define I_bits
#endif

	  if (optind < argc && base_name == NULL) {
		  base_name = argv[optind];
		  optind++;
	  }
	  if (optind < argc)fprintf(stderr, "Ignoring %u trailing command line args\n",
		  argc - optind);
	  if (base_name == NULL)base_name = "gnfs";
	  if ((input_data = fopen(base_name, "r")) == NULL) {
		  complain("Cannot open %s for input of nfs polynomials: %m\n", base_name);
	  }
	  mpz_init(N);
	  mpz_init(m);
	  mpz_init(aux1);
	  mpz_init(aux2);
	  mpz_init(aux3);
	  mpz_init(sr_a);
	  mpz_init(sr_b);
	  mpz_ull_init();
	  mpz_init(rational_rest);
	  mpz_init(algebraic_rest);

	  // for tinyecm
	  mpz_init(factor1);
	  mpz_init(factor2);
	  mpz_init(factor3);
	  pran = 42;

#ifdef AVX512_SIEVE4
	  // for the small-prime sieve
	  rtab_s0 = (uint64_t***)malloc(10 * sizeof(uint64_t**));
	  rtab_s1 = (uint64_t***)malloc(10 * sizeof(uint64_t**));
	  for (i = 0; i < 10; i++) {
		  rtab_s0[i] = (uint64_t**)malloc(32 * sizeof(uint64_t*));
		  rtab_s1[i] = (uint64_t**)malloc(32 * sizeof(uint64_t*));
		  int j;
		  for (j = 0; j < 32; j++) {
			  rtab_s0[i][j] = (uint64_t*)malloc(8 * sizeof(uint64_t));
			  rtab_s1[i][j] = (uint64_t*)malloc(8 * sizeof(uint64_t));
		  }
	  }
#endif

	  //input_poly(N,poly,poldeg,poly+1,poldeg+1,m,input_data);
#if 0
	  if (poldeg[1] > 1) {
		  if (poldeg[0] == 1) {
			  mpz_t* X;
			  poldeg[0] = poldeg[1];
			  poldeg[1] = 1;
			  X = poly[0];
			  poly[0] = poly[1];
			  poly[1] = X;
		  }
		  else {
			  complain("Degrees >1 on both sides not implemented\n");
		  }
	  }
#endif
	  //skip_blanks_comments(&input_line,&input_line_alloc,input_data);
	  fclose(input_data);
	  // parse the input poly file
	  {
		  FILE* fp;
		  char token[256], value[512], thisLine[1024];


		  sieve_min[0] = sieve_min[1] = 0;

		  if (!(fp = fopen(base_name, "rb"))) {
			  printf("Error opening %s for read!\n", base_name);
			  return -1;
		  }
		  input_poly(N, poly, poldeg, poly + 1, poldeg + 1, m, fp);
		  rewind(fp);
		  while (!(feof(fp))) {
			  thisLine[0] = 0;
			  fgets(thisLine, 1023, fp);
			  /* Special case: If there's a polynomial, handle it seperately: */
			  if (strncmp(thisLine, "START_POLY", 10) == 0) {
				  while (!(feof(fp)) && strncmp(thisLine, "END_POLY", 8))
					  fgets(thisLine, 1023, fp);
			  }
			  else  if ((sscanf(thisLine, "%255s %511s", token, value) == 2) &&
				  (thisLine[0] != '#')) {

				  token[sizeof(token) - 1] = 0;
				  if (strncmp(token, "skew:", 5) == 0) {
					  sigma = (float)atof(value);
				  }
				  else if (strncmp(token, "q0:", 3) == 0) {
					  first_spq = atol(value);
					  first_spq1 = first_spq;
					  first_root = 0;

				  }
				  else if (strncmp(token, "qintsize:", 9) == 0) {
					  sieve_count = atol(value);
				  }
				  else if ((strncmp(token, "skip0:", 6) == 0) ||
					  (strncmp(token, "askip:", 6) == 0)) {
					  sieve_min[0] = atol(value);
				  }
				  else if ((strncmp(token, "skip1:", 6) == 0) ||
					  (strncmp(token, "rskip:", 6) == 0)) {
					  sieve_min[1] = atol(value);
				  }
				  else if ((strncmp(token, "lim0:", 5) == 0) ||
					  (strncmp(token, "alim:", 5) == 0)) {
					  FB_bound[0] = (float)atol(value);
				  }
				  else if ((strncmp(token, "lim1:", 5) == 0) ||
					  (strncmp(token, "rlim:", 5) == 0)) {
					  FB_bound[1] = (float)atof(value);
				  }
				  else if ((strncmp(token, "lpb0:", 5) == 0) ||
					  (strncmp(token, "lpba:", 5) == 0)) {
					  max_primebits[0] = atoi(value);
				  }
				  else if ((strncmp(token, "lpb1:", 5) == 0) ||
					  (strncmp(token, "lpbr:", 5) == 0)) {
					  max_primebits[1] = atoi(value);
				  }
				  else if ((strncmp(token, "mfb0:", 5) == 0) ||
					  (strncmp(token, "mfba:", 5) == 0)) {
					  max_factorbits[0] = atoi(value);
				  }
				  else if ((strncmp(token, "mfb1:", 5) == 0) ||
					  (strncmp(token, "mfbr:", 5) == 0)) {
					  value[sizeof(value) - 1] = 0;
					  max_factorbits[1] = atoi(value);
				  }
				  else if ((strncmp(token, "lambda0:", 8) == 0) ||
					  (strncmp(token, "alambda:", 8) == 0)) {
					  sieve_report_multiplier[0] = (float)atof(value);
				  }
				  else if ((strncmp(token, "lambda1:", 8) == 0) ||
					  (strncmp(token, "rlambda:", 8) == 0)) {
					  sieve_report_multiplier[1] = (float)atof(value);
				  }
#ifdef _NO
				  else {
					  printf("Warning: Ignoring input line:\n%s\n", thisLine);
				  }
#endif
			  }
		  }
		  fclose(fp);
	  }



	  last_spq = first_spq + sieve_count;
	  if (last_spq >= I32_MAX / 2) {
		  complain("Cannot handle special q >= %d\n", I32_MAX / 2);
	  }
	  for (i = 0; i < 2; i++) {
		  if (FB_bound[i] < 4 || sieve_report_multiplier[i] <= 0) {
			  complain("Please set all bounds to reasonable values!\n");
		  }
#if 1
		  if (max_primebits[i] > 33) {
			  complain("Only large primes up to 33 bits are allowed.\n");
		  }
#endif
	  }
	  if (sieve_count != 0) {
		  if (sigma == 0)complain("Please set a skewness\n");
		  if (special_q_side == NO_SIDE) {
			  errprintf("Please use -a or -r\n");
			  Usage();
		  }
		  if (FB_bound[special_q_side] > first_spq) {
			  FB_bound[special_q_side] = (float)first_spq - 1;
			  if (verbose) printf("Warning:  lowering FB_bound to %u.\n", first_spq - 1);
			  //complain("Special q lower bound %u below rFB bound %g\n",
			  //first_spq,FB_bound[special_q_side]);
		  }
	  }
	  //fclose(input_data);
	  if (poldeg[0] < poldeg[1])poldeg_max = poldeg[1];
	  else poldeg_max = poldeg[0];



	  i_shift = 1 << (I_bits - 1);
	  n_I = 1 << I_bits;
	  n_i = n_I / 2;
	  i_bits = I_bits - 1;
	  n_J = 1 << J_bits;
	  n_j = n_J / 2;
	  j_bits = J_bits - 1;
	  {
		  u32_t i, j;
		  double x, y, z;

		  x = sqrt(first_spq * sigma) * n_I;
		  y = x / sigma;
		  for (j = 0; j < 2; j++) {
			  poly_f[j] = xmalloc((poldeg[j] + 1) * sizeof(*poly_f[j]));

			  for (i = 0, z = 1, poly_norm[j] = 0;
				  i <= poldeg[j]; i++) {
				  poly_f[j][i] = mpz_get_d(poly[j][i]);
				  poly_norm[j] = poly_norm[j] * y + fabs(poly_f[j][i]) * z;
				  z *= x;
			  }
		  }
	  }

  }
  
  siever_init();
  // handle opening the output file
  if (sieve_count != 0) {
	  if (g_ofile_name == NULL) {
		  if (zip_output == 0) {
			  asprintf(&g_ofile_name, "%s.lasieve-%u.%u-%u", base_name,
				  special_q_side, first_spq, last_spq);
		  }
		  else {
			  asprintf(&g_ofile_name,
				  append_output == 0 ?
				  "gzip --best --stdout > %s.lasieve-%u.%u-%u.gz" :
				  "gzip --best --stdout >> %s.lasieve-%u.%u-%u.gz",
				  base_name, special_q_side, first_spq, last_spq);
		  }
	  }
	  else {
		  if (strcmp(g_ofile_name, "-") == 0) {
			  if (zip_output == 0) {
				  g_ofile = stdout;
				  g_ofile_name = "to stdout";
				  goto done_opening_output;
			  }
			  else g_ofile_name = "gzip --best --stdout";
		  }
		  else {
			  if (fnmatch("*.gz", g_ofile_name, 0) == 0) {
				  char* on1;

				  zip_output = 1;
				  on1 = strdup(g_ofile_name);
				  asprintf(&g_ofile_name, "gzip --best --stdout > %s", on1);
				  free(on1);
			  }
			  else zip_output = 0;
		  }
	  }
	  if (zip_output == 0) {
		  if (g_resume != 0) {
			  goto done_opening_output;
		  }
		  if (append_output > 0) {
			  g_ofile = fopen(g_ofile_name, "a");
		  }
		  else {
			  if ((g_ofile = fopen(g_ofile_name, "r")) != NULL)
				  complain(" Will not overwrite existing file %s for output; rename it, move it away, or use -R option (resume)\n", g_ofile_name);
			  g_ofile = fopen(g_ofile_name, "w");
		  }
		  if (g_ofile == NULL)complain("Cannot open %s for output: %m\n", g_ofile_name);
	  }
	  else {
		  if ((g_ofile = popen(g_ofile_name, "w")) == NULL)
			  complain("Cannot exec %s for output: %m\n", g_ofile_name);
	  }
  done_opening_output:
	  fprintf(g_ofile, "");
  }
  
  // read or build the rational and algebraic factor base
  {
	  size_t FBS_alloc = 4096;
	  u32_t prime;
	  pr32_struct ps;
	  char* afbname;
	  FILE* afbfile;
	  u32_t side;

	  initprime32(&ps);

	  for (side = 0; side < 2; side++) {
		  if (poldeg[side] == 1) {
			  FB[side] = xmalloc(FBS_alloc * sizeof(u32_t));
			  proots[side] = xmalloc(FBS_alloc * sizeof(u32_t));
			  prime = firstprime32(&ps);
			  for (prime = nextprime32(&ps), fbi1[side] = 0, FBsize[side] = 0;
				  prime < FB_bound[side]; prime = nextprime32(&ps)) {
				  u32_t x;
				  x = mpz_fdiv_ui(poly[side][1], prime);
				  if (x > 0) {
					  modulo32 = prime;
					  x = modmul32(modinv32(x), mpz_fdiv_ui(poly[side][0], prime));
					  x = x > 0 ? prime - x : 0;
				  }
				  else x = prime;
				  if (prime < L1_SIZE)fbi1[side] = FBsize[side];
				  if (prime < n_i)fbis[side] = FBsize[side];
				  if (FBsize[side] == FBS_alloc) {
					  FBS_alloc *= 2;
					  FB[side] = xrealloc(FB[side], FBS_alloc * sizeof(u32_t));
					  proots[side] = xrealloc(proots[side], FBS_alloc * sizeof(u32_t));
				  }
				  proots[side][FBsize[side]] = x;
				  FB[side][FBsize[side]++] = prime;
			  }
			  proots[side] = xrealloc(proots[side], FBsize[side] * sizeof(u32_t));
			  FB[side] = xrealloc(FB[side], FBsize[side] * sizeof(u32_t));
			  fbi1[side]++;
			  fbis[side]++;
			  if (fbi1[side] < fbis[side])fbi1[side] = fbis[side];
		  }
		  else {
			  u32_t j, k, l;
			  asprintf(&afbname, "%s.afb.%u", base_name, side);
			  if (force_aFBcalc > 0 || (afbfile = fopen(afbname, "r")) == NULL) {
				  u32_t* root_buffer;
				  size_t aFB_alloc;

				  root_buffer = xmalloc(poldeg[side] * sizeof(*root_buffer));
				  aFB_alloc = 4096;
				  FB[side] = xmalloc(aFB_alloc * sizeof(**FB));
				  proots[side] = xmalloc(aFB_alloc * sizeof(**proots));
				  for (prime = firstprime32(&ps), FBsize[side] = 0;
					  prime < FB_bound[side]; prime = nextprime32(&ps)) {
					  u32_t i, nr;

					  nr = root_finder(root_buffer, poly[side], poldeg[side], prime);
					  for (i = 0; i < nr; i++) {
						  if (aFB_alloc <= FBsize[side]) {
							  aFB_alloc *= 2;
							  FB[side] = xrealloc(FB[side], aFB_alloc * sizeof(**FB));
							  proots[side] = xrealloc(proots[side], aFB_alloc * sizeof(**proots));
						  }
						  FB[side][FBsize[side]] = prime;
						  proots[side][FBsize[side]] = root_buffer[i];
						  if (prime > 2)FBsize[side]++;
					  }
				  }
				  FB[side] = xrealloc(FB[side], FBsize[side] * sizeof(**FB));
				  proots[side] = xrealloc(proots[side], FBsize[side] * sizeof(**proots));
				  free(root_buffer);

				  if (keep_factorbase > 0) {
					  if ((afbfile = fopen(afbname, "w")) == NULL) {
						  complain("Cannot open %s for output of aFB: %m\n", afbname);
					  }
					  if (write_u32(afbfile, &(FBsize[side]), 1) != 1) {
						  complain("Cannot write aFBsize to %s: %m\n", afbname);
					  }
					  if (write_u32(afbfile, FB[side], FBsize[side]) != FBsize[side] ||
						  write_u32(afbfile, proots[side], FBsize[side]) != FBsize[side]) {
						  complain("Cannot write aFB to %s: %m\n", afbname);
					  }
					  if (write_u32(afbfile, &xFBs[side], 1) != 1) {
						  complain("Cannot write aFBsize to %s: %m\n", afbname);
					  }
					  fclose(afbfile);
				  }

			  }
			  else {
				  if (read_u32(afbfile, &(FBsize[side]), 1) != 1) {
					  complain("Cannot read aFB size from %s: %m\n", afbname);
				  }
				  FB[side] = xmalloc(FBsize[side] * sizeof(u32_t));
				  proots[side] = xmalloc(FBsize[side] * sizeof(u32_t));
				  if (read_u32(afbfile, FB[side], FBsize[side]) != FBsize[side] ||
					  read_u32(afbfile, proots[side], FBsize[side]) != FBsize[side]) {
					  complain("Cannot read aFB from %s: %m\n", afbname);
				  }
				  if (read_u32(afbfile, &xFBs[side], 1) != 1) {
					  complain("%s: Cannot read xFBsize\n", afbname);
				  }
				  fclose(afbfile);

			  }


			  for (j = 0, k = 0, l = 0; j < FBsize[side]; j++) {
				  if (FB[side][j] < L1_SIZE)k = j;
				  if (FB[side][j] < n_i)l = j;
				  if (FB[side][j] > L1_SIZE && FB[side][j] > n_I)break;
			  }
			  if (FBsize[side] > 0) {
				  if (k < l)k = l;
				  fbis[side] = l + 1;
				  fbi1[side] = k + 1;
			  }
			  else {
				  fbis[side] = 0;
				  fbi1[side] = 0;
			  }
		  }
	  }

	  {
		  u32_t i, srfbs, safbs;

		  for (i = 0, srfbs = 0; i < xFBs[1]; i++) {
			  if (xFB[1][i].p == xFB[1][i].pp)srfbs++;
		  }
		  for (i = 0, safbs = 0; i < xFBs[0]; i++) {
			  if (xFB[0][i].p == xFB[0][i].pp)safbs++;
		  }
		  logbook(0, "FBsize %u+%u (deg %u), %u+%u (deg %u)\n",
			  FBsize[0], safbs, poldeg[0], FBsize[1], srfbs, poldeg[1]);
	  }
	  free(afbname);

	  {
		  u32_t i;
		  size_t si, sj;

		  n_srb_i = 2 * ((n_i + 2 * CANDIDATE_SEARCH_STEPS - 1) / (2 * CANDIDATE_SEARCH_STEPS));
		  n_srb_j = (n_J + 2 * CANDIDATE_SEARCH_STEPS - 1) / (2 * CANDIDATE_SEARCH_STEPS);
		  sj = n_srb_j * sizeof(*(sieve_report_bounds[0]));
		  si = n_srb_i * sizeof(**(sieve_report_bounds[0]));
		  for (i = 0; i < 2; i++) {
			  u32_t j;

			  tpoly_f[i] = xmalloc((1 + poldeg[i]) * sizeof(**tpoly_f));
			  sieve_report_bounds[i] = xmalloc(sj);
			  for (j = 0; j < n_srb_j; j++)
				  sieve_report_bounds[i][j] = xmalloc(si);
		  }
	  }

  }
  
  if(sieve_count == 0)
	  exit(0);

  // more work on the factor bases
  {
	  u32_t side, i;

	  for (side = 0; side < 2; side++) {
		  u32_t prime, nr;
		  struct xFBstruct* s;
		  u32_t* root_buffer;
		  size_t xaFB_alloc = 0;
		  FB_logs[side] = xmalloc(FBsize[side]);
		  sieve_multiplier[side] = (UCHAR_MAX - 50) / log(poly_norm[side]);

		  root_buffer = xmalloc(poldeg[side] * sizeof(*root_buffer));
		  prime = 2;
		  nr = root_finder(root_buffer, poly[side], poldeg[side], prime);

		  for (i = 0; i < nr; i++) {
			  adjust_bufsize((void**)&(xFB[side]), &xaFB_alloc, 1 + xFBs[side],
				  16, sizeof(**xFB));
			  s = xFB[side] + xFBs[side];
			  s->p = prime;
			  s->pp = prime;
			  if (root_buffer[i] == prime) {
				  s->qq = prime;
				  s->q = 1;
				  s->r = 1;
			  }
			  else {
				  s->qq = 1;
				  s->q = prime;
				  s->r = root_buffer[i];
			  }
			  xFBs[side]++;
			  add_primepowers2xaFB(&xaFB_alloc, n_I, side, 0, 0);
		  }
		  free(root_buffer);
		  for (i = 0; i < FBsize[side]; i++) {
			  double l;
			  u32_t l1;

			  prime = FB[side][i];
			  if (prime > n_I / prime)break;
			  l = log(prime);
			  l1 = add_primepowers2xaFB(&xaFB_alloc, n_I, side, prime, proots[side][i]);
			  FB_logs[side][i] = rint(l1 * l * sieve_multiplier[side]);
		  }
		  while (i < FBsize[side]) {
			  double l;

			  l = log(FB[side][i]);
			  if (l > FB_maxlog[side])FB_maxlog[side] = l;
			  FB_logs[side][i++] = rint(sieve_multiplier[side] * l);
		  }
		  FB_maxlog[side] *= sieve_multiplier[side];
		  qsort(xFB[side], xFBs[side], sizeof(*(xFB[side])), xFBcmp);
	  }
  }
  
  sieve_interval= xvalloc(L1_SIZE);
  cand= xvalloc(L1_SIZE*sizeof(*cand));
  fss_sv= xvalloc(L1_SIZE);
  tiny_sieve_buffer= xmalloc(TINY_SIEVEBUFFER_SIZE);
  if(n_i> L1_SIZE)
    complain("Strip length %u exceeds L1 size %u\n",n_i,L1_SIZE);
#if I_bits<=L1_BITS
  j_per_strip= L1_SIZE/n_i;
  jps_bits= L1_BITS-i_bits;
#endif
  if(j_per_strip!=1<<jps_bits)
    Schlendrian("Expected %u j per strip, calculated %u\n",
		j_per_strip,1<<jps_bits);
  n_strips= n_j>>(L1_BITS-i_bits);
  rec_info_init(n_i,n_j);

  // memory allocation
  {
	  u32_t s;
#define MAX_TINY_2POW 4

	  if (poldeg[0] < poldeg[1])s = poldeg[1];
	  else s = poldeg[0];
	  tinysieve_curpos = xmalloc(TINY_SIEVE_MIN * s * sizeof(*tinysieve_curpos));
	  horizontal_sievesums = xmalloc(j_per_strip * sizeof(*horizontal_sievesums));
	  for (s = 0; s < 2; s++) {
		  u32_t fbi;
		  size_t maxent;

		  smallsieve_aux[s] = xmalloc(4 * fbis[s] * sizeof(*(smallsieve_aux[s])));
#ifndef MMX_TD
#ifdef PREINVERT
		  smalltd_pi[s] = xmalloc(fbis[s] * sizeof(*(smalltd_pi[s])));
#endif
		  smalltdsieve_aux[s] = xmalloc(j_per_strip * sizeof(*(smalltdsieve_aux[s])));
		  for (fbi = 0; fbi < j_per_strip; fbi++)
			  smalltdsieve_aux[s][fbi] =
			  xmalloc(fbis[s] * sizeof(**(smalltdsieve_aux[s])));
#else

		  MMX_TdAllocate(j_per_strip, fbis[0], fbis[1]);
#endif
		  smallsieve_aux1[s] = xmalloc(6 * xFBs[s] * sizeof(*(smallsieve_aux1[s])));


		  maxent = fbis[s];
		  maxent += xFBs[s];
		  smallpsieve_aux[s] = xmalloc(3 * maxent * sizeof(*(smallpsieve_aux[s])));
		  maxent = 0;
		  for (fbi = 0; fbi < xFBs[s]; fbi++) {
			  if (xFB[s][fbi].p == 2)
				  maxent++;
		  }
		  smallsieve_aux2[s] = xmalloc(4 * maxent * sizeof(*(smallsieve_aux2[s])));
		  x2FB[s] = xmalloc(maxent * 6 * sizeof(*(x2FB[s])));
	  }
  }
  
#ifdef GCD_SIEVE_BOUND
  {
    u32_t p,i;
    
    firstprime32(&special_q_ps);
    np_gcd_sieve= 0;
    for(p= nextprime32(&special_q_ps);p<GCD_SIEVE_BOUND;
	p= nextprime32(&special_q_ps))np_gcd_sieve++;
    gcd_sieve_buffer= xmalloc(2*np_gcd_sieve*sizeof(*gcd_sieve_buffer));
    
    firstprime32(&special_q_ps);
    i= 0;
    for(p= nextprime32(&special_q_ps);p<GCD_SIEVE_BOUND;
	p= nextprime32(&special_q_ps))gcd_sieve_buffer[2*i++]= p;
  }
#endif
  
  {
    u32_t s;
    for(s= 0;s<2;s++) {
      if(sieve_min[s]<TINY_SIEVE_MIN&&sieve_min[s]!=0) {
	errprintf("Sieving with all primes on side %u since\n",s);
	errprintf("tiny sieve procedure is being used\n");
	sieve_min[s]= 0;
      }
      current_ij[s]= xmalloc(FBsize[s]*sizeof(*current_ij[s]));
      LPri[s]= xmalloc(FBsize[s]*sizeof(**LPri)*RI_SIZE);
    }
  }
  
  // more memory allocation and bucket sieve setup.
  {
	  u32_t s;
	  size_t total_alloc;
	  u16_t* sched_buf;
	  double pvl_max[2];

	  total_alloc = 0;
	  for (s = 0; s < 2; s++) {
		  u32_t i, fbi_lb;

		  if (sigma >= 1)pvl_max[s] = poldeg[s] * log(last_spq * sqrt(sigma));
		  else pvl_max[s] = poldeg[s] * log(last_spq / sqrt(sigma));
		  pvl_max[s] += log(poly_norm[s]);
		  if (fbi1[s] >= FBsize[s] || i_bits + j_bits <= L1_BITS) {
			  n_schedules[s] = 0;
			  continue;
		  }
		  for (i = 0; i < N_PRIMEBOUNDS; i++)
			  if (FB_bound[s] <= schedule_primebounds[i] ||
				  i_bits + j_bits <= schedule_sizebits[i])
				  break;
		  n_schedules[s] = i + 1;
		  schedules[s] = xmalloc(n_schedules[s] * sizeof(**schedules));
		  fbi_lb = fbi1[s];
		  //printf("%d schedules on side %d\n", n_schedules[s], s);
		  for (i = 0; i < n_schedules[s]; i++) {
			  u32_t fbp_lb, fbp_ub;
			  u32_t fbi, fbi_ub;
			  u32_t sp_i;
			  u32_t n, sl_i;
			  u32_t ns;
			  size_t allocate, all1;

			  if (i == n_schedules[s] - 1)fbp_ub = FB_bound[s];
			  else fbp_ub = schedule_primebounds[i];
			  if (i == 0)fbp_lb = FB[s][fbi1[s]];
			  else fbp_lb = schedule_primebounds[i - 1];

			  if (i_bits + j_bits < schedule_sizebits[i])ns = 1 << (i_bits + j_bits - L1_BITS);
			  else ns = 1 << (schedule_sizebits[i] - L1_BITS);
			  schedules[s][i].n_strips = ns;



#ifndef SCHED_TOL
#ifndef NO_SCHEDTOL
			  /* these values are experimental; report SCHED_PATHOLOGY to http://mersenneforum.org/showthread.php?t=11430 */
#define SCHED_PAD 48
#if I_bits<15
#define SCHED_TOL 2
#else
#define SCHED_TOL 1.2
#endif
#endif
#endif
#ifdef SCHED_TOL
			  allocate = rint(SCHED_PAD + SCHED_TOL * n_i * j_per_strip * log(log(fbp_ub) / log(fbp_lb)));
#else
			  allocate = rint(sched_tol[i] * n_i * j_per_strip * log(log(fbp_ub) / log(fbp_lb)));
#endif
			  allocate *= SE_SIZE;

			  all1 = allocate + n_i * ceil(pvl_max[s] / log(fbp_lb)) * SE_SIZE;
			  schedules[s][i].alloc = allocate;
			  schedules[s][i].alloc1 = all1;
			  //printf("### [%d][%d] alloc = %ld %ld %ld %ld\n",s,i,allocate,all1,n_i,j_per_strip);

			  for (n = 0, fbi = fbi_lb; fbi < FBsize[s];) {
				  u32_t fbi_ub1;
				  fbi_ub1 = fbi + SCHEDFBI_MAXSTEP;
				  if (fbi_ub1 >= FBsize[s])fbi_ub1 = FBsize[s];
				  else {
					  if (FB[s][fbi_ub1] > fbp_ub) {
						  while (FB[s][fbi_ub1] > fbp_ub)
							  fbi_ub1--;
						  fbi_ub1++;
					  }
				  }
				  if (FB_logs[s][fbi] == FB_logs[s][fbi_ub1 - 1]) {
					  n++;
					  fbi = fbi_ub1;
				  }
				  else {
					  u32_t l;
					  n += FB_logs[s][fbi_ub1 - 1] - FB_logs[s][fbi];
					  fbi = fbi_ub1 - 1;
					  l = FB_logs[s][fbi];
					  while (FB_logs[s][fbi] == l)fbi--;
					  fbi++;
				  }
				  if (fbi >= FBsize[s] || FB[s][fbi] > fbp_ub)break;
			  }
			  fbi_ub = fbi;
			  schedules[s][i].n_pieces = n;
			  n++;
			  schedules[s][i].schedule = xmalloc(n * sizeof(*(schedules[s][i].schedule)));
			  for (sl_i = 0; sl_i < n; sl_i++)
				  schedules[s][i].schedule[sl_i] =
				  xmalloc(ns * sizeof(**(schedules[s][i].schedule)));
			  schedules[s][i].schedule[0][0] = (u16_t*)total_alloc;
			  total_alloc += all1;
			  for (sp_i = 1; sp_i < ns; sp_i++) {
				  schedules[s][i].schedule[0][sp_i] = (u16_t*)total_alloc;
				  total_alloc += allocate;
			  }
			  schedules[s][i].fbi_bounds =
				  xmalloc(n * sizeof(*(schedules[s][i].fbi_bounds)));
			  schedules[s][i].schedlogs = xmalloc(n);
			  for (n = 0, fbi = fbi_lb; fbi < fbi_ub;) {
				  u32_t fbi_ub1;
				  fbi_ub1 = fbi + SCHEDFBI_MAXSTEP;
				  if (fbi_ub1 > fbi_ub)fbi_ub1 = fbi_ub;
				  if (FB_logs[s][fbi] == FB_logs[s][fbi_ub1 - 1]) {
					  schedules[s][i].fbi_bounds[n++] = fbi;
					  fbi = fbi_ub1;
				  }
				  else {
					  u32_t l, lmax;

					  lmax = FB_logs[s][fbi_ub1 - 1];
					  for (l = FB_logs[s][fbi]; l < lmax; l++) {
						  schedules[s][i].fbi_bounds[n++] = fbi;
						  while (fbi < fbi_ub && FB_logs[s][fbi] == l)fbi++;
					  }
				  }
			  }
			  if (n != schedules[s][i].n_pieces)
				  Schlendrian("Expected %u schedule pieces on side %u, have %u\n",
					  schedules[s][i].n_pieces, s, n);
			  schedules[s][i].fbi_bounds[n++] = fbi;
			  for (n = 0; n < schedules[s][i].n_pieces; n++)
				  schedules[s][i].schedlogs[n] = FB_logs[s][schedules[s][i].fbi_bounds[n]];
#ifdef CONTIGUOUS_RI
			  schedules[s][i].ri =
				  LPri[s] + (schedules[s][i].fbi_bounds[0] - fbis[s]);// *RI_SIZE;
#else
			  schedules[s][i].ri =
				  LPri[s] + (schedules[s][i].fbi_bounds[0] - fbis[s]) * RI_SIZE;
#endif
			 // printf("\tschedule offset %d = %u, (%u elements)", 
		//		  i, (schedules[s][i].fbi_bounds[0] - fbis[s]), 
		//		  (schedules[s][i].fbi_bounds[0] - fbis[s]));// *RI_SIZE);
		//	  if (i > 0)
		//		  printf(" diff = %u\n", (schedules[s][i].fbi_bounds[0] - fbis[s]) -
		//			  (schedules[s][i - 1].fbi_bounds[0] - fbis[s]));
		//	  else
		//		  printf("\n");
			  fbi_lb = fbi_ub;
		  }
	  }

	  sched_buf = xmalloc((total_alloc + 65536 * SE_SIZE * j_per_strip) *
		  sizeof(**((**schedules).schedule)));
	  for (s = 0; s < 2; s++) {
		  u32_t i;
		  for (i = 0; i < n_schedules[s]; i++) {
			  u32_t sp_i;

			  for (sp_i = 0; sp_i < schedules[s][i].n_strips; sp_i++)
				  schedules[s][i].schedule[0][sp_i] =
				  sched_buf + (size_t)(schedules[s][i].schedule[0][sp_i]);
		  }
	  }

#ifdef USE_MEDSCHED
	  {
		  u32_t s;

		  for (s = 0; s < 2; s++) {
			  if (fbis[s] < fbi1[s]) {
				  u32_t fbi;
				  u32_t n;
				  unsigned char oldlog;


				  medsched_alloc[s] = j_per_strip * (fbi1[s] - fbis[s]) * SE_SIZE;

				  medsched_alloc[s] += n_i * ceil(pvl_max[s] / log(n_i)) * SE_SIZE;
				  n_medsched_pieces[s] = 1 + FB_logs[s][fbi1[s] - 1] - FB_logs[s][fbis[s]];
				  med_sched[s] = xmalloc((1 + n_medsched_pieces[s]) * sizeof(**med_sched));
				  med_sched[s][0] = xmalloc(medsched_alloc[s] * sizeof(***med_sched));

				  medsched_fbi_bounds[s] =
					  xmalloc((1 + n_medsched_pieces[s]) * sizeof(**medsched_fbi_bounds));
				  medsched_logs[s] = xmalloc(n_medsched_pieces[s]);

				  for (n = 0, fbi = fbis[s], oldlog = UCHAR_MAX; fbi < fbi1[s]; fbi++) {
					  if (FB_logs[s][fbi] != oldlog) {
						  medsched_fbi_bounds[s][n] = fbi;
						  oldlog = FB_logs[s][fbi];
						  medsched_logs[s][n++] = oldlog;
					  }
				  }
				  if (n != n_medsched_pieces[s])
					  Schlendrian("Expected %u medium schedule pieces on side %u, have %u\n",
						  n_medsched_pieces[s], s, n);
				  medsched_fbi_bounds[s][n] = fbi;
			  }
			  else {

				  n_medsched_pieces[s] = 0;
			  }
		  }
	  }
#endif

  }
  
  {
	  u32_t s;
	  size_t schedbuf_alloc;

	  for (s = 0, schedbuf_alloc = 0; s < 2; s++) {
		  u32_t i;

		  for (i = 0; i < n_schedules[s]; i++)
			  if (schedules[s][i].n_pieces > schedbuf_alloc)
				  schedbuf_alloc = schedules[s][i].n_pieces;
	  }
	  schedbuf = xmalloc((1 + schedbuf_alloc) * sizeof(*schedbuf));
  }
  
  td_buf1= xmalloc((1+L1_SIZE)*sizeof(*td_buf1));
  td_buf[0]= xmalloc(td_buf_alloc[0]*sizeof(**td_buf));
  td_buf[1]= xmalloc(td_buf_alloc[1]*sizeof(**td_buf));
  
  {
    u32_t i;
    if(tds_fbi == NULL) {
      tds_fbi= xmalloc(UCHAR_MAX*sizeof(*tds_fbi));
      tds_fbi_curpos= xmalloc(UCHAR_MAX*sizeof(*tds_fbi));
      for(i= 0;i<UCHAR_MAX;i++)
	tds_fbi[i]= xmalloc(tds_fbi_alloc*sizeof(**tds_fbi));
    }
  }
  
  {
    u32_t s,i;
    for(i= 0;i<L1_SIZE;i++) {
      mpz_init(td_rests[i]);
    }
    for(s= 0;s<2;s++) {
      mpz_init(large_factors[s]);
      large_primes[s]= xmalloc(max_factorbits[s]*sizeof(*(large_primes[s])));
      for(i= 0;i<max_factorbits[s];i++) {
	mpz_init(large_primes[s][i]);
      }
      mpz_init_set_d(FBb_sq[s],FB_bound[s]);
      mpz_mul(FBb_sq[s],FBb_sq[s],FBb_sq[s]);
      mpz_init_set_d(FBb_cu[s],FB_bound[s]);
      mpz_pow_ui(FBb_cu[s],FBb_cu[s],3);
    }
  }
  
  
  all_spq_done= 1;
  // process all special_q
  {
	  u32_t* r;
	  initprime32(&special_q_ps);
	  last_clock = clock();
	  n_spq = 0;
	  n_spq_discard = 0;
	  r = xmalloc(poldeg_max * sizeof(*r));
	  special_q = pr32_seek(&special_q_ps, first_spq1);
	  if (catch_signals != 0) {
		  signal(SIGTERM, terminate_sieving);
		  signal(SIGINT, terminate_sieving);
	  }
	  tStart = lastReport = sTime();
	  for (; special_q < last_spq && special_q != 0; 
		  special_q = nextprime32(&special_q_ps), first_root = 0) 
	  {
		  u32_t nr;

		  special_q_log = log(special_q);
		  if (cmdline_first_sieve_side == USHRT_MAX) {
#if 1
			  double nn[2];
			  u32_t s;
#if 0
			  for (s = 0; s < 2; s++) {
				  nn[s] = log(poly_norm[s] * (special_q_side == s ? 1 : special_q));
				  nn[s] = nn[s] / log(FB_bound[s]) - sieve_report_multiplier[s];
		  }
#else
			  for (s = 0; s < 2; s++) {
				  nn[s] = log(poly_norm[s] * (special_q_side == s ? 1 : special_q));
				  nn[s] = nn[s] / (sieve_report_multiplier[s] * log(FB_bound[s]));
			  }
#endif
			  if (nn[0] < nn[1])first_sieve_side = 1;
			  else first_sieve_side = 0;
#else
			  if (poly_norm[0] * (special_q_side == 0 ? 1 : special_q)
				  < poly_norm[1] * (special_q_side == 1 ? 1 : special_q)) {
				  first_sieve_side = 1;
			  }
			  else {
				  first_sieve_side = 0;
	  }
#endif
		  }
		  else
		  {
			  first_sieve_side = cmdline_first_sieve_side;
			  if (first_sieve_side >= 2)complain("First sieve side must not be %u\n",
				  (u32_t)first_sieve_side);
		  }

		  logbook(1,"First sieve side: %u\n",(u32_t)first_sieve_side);
		  if(cmdline_first_td_side!=USHRT_MAX)first_td_side= cmdline_first_td_side;
		  else first_td_side= first_sieve_side;
		  if (poldeg[special_q_side] > 1) {
			  nr = root_finder(r, poly[special_q_side], poldeg[special_q_side], special_q);
			  if (nr == 0)continue;
			  if (r[nr - 1] == special_q) {


				  nr--;
			  }
		  }
		  else {
			  u32_t x = mpz_fdiv_ui(poly[special_q_side][1], special_q);
			  if (x == 0) {
				  n_spq_discard++;
				  continue;
			  }
			  modulo32 = special_q;
			  x = modmul32(modinv32(x), mpz_fdiv_ui(poly[special_q_side][0], special_q));
			  r[0] = x == 0 ? 0 : special_q - x;
			  nr = 1;
		  }
      
		  for(root_no= 0;root_no<nr;root_no++) {
			u32_t termination_condition;
	
			if(r[root_no]<first_root)continue;
			if((termination_condition= setjmp(termination_jb))!=0) {
			  if(termination_condition == USER_INTERRUPT)
				{
				  //char*hn,*ofn;
				  //FILE*of;
	      
				  //hn= xmalloc(100);
				  //if(gethostname(hn,99) == 0)
				  //asprintf(&ofn,"%s.%s.last_spq%d",base_name,hn,process_no);
				  //else asprintf(&ofn,"%s.unknown_host.last_spq%d",base_name,process_no);
				  //free(hn);
	      
				  //if((of= fopen(ofn,"w"))!=0) {
				  //fprintf(of,"%u\n",special_q);
				  //fclose(of);
				  //}
				  //free(ofn);
				  all_spq_done= 0;
				  break;
				}
	  
			  else{
	#if 0
			char *cmd;

			asprintf(&cmd,"touch badsched.%s.%u.%u.%u",base_name,
						special_q_side,special_q,r[root_no]);
			system(cmd);
			free(cmd);
	#endif
				continue;
			  }
			}

			if(sysload_cmd!=NULL) {
			  if(system(sysload_cmd)!=0) {
				exitval= LOADTEST_EXITVAL;
				longjmp(termination_jb,USER_INTERRUPT);
			  }
			}
			n_spq++;

			reduce2(&a0,&b0,&a1,&b1,(i32_t)special_q,0,(i32_t)r[root_no],1,sigma*sigma);

			{
				if (b0 % ((i32_t)special_q) == 0 && b1 % ((i32_t)special_q) == 0) {
					i32_t x;

					x = a0 % ((i32_t)special_q);
					if (x < 0)x += (i32_t)special_q;
					spq_i = x;
					x = a1 % ((i32_t)special_q);
					if (x < 0)x += (i32_t)special_q;
					spq_j = x;
				}
				else {
					i32_t x;

					x = b0 % ((i32_t)special_q);
					if (x < 0)x += (i32_t)special_q;
					spq_i = x;
					x = b1 % ((i32_t)special_q);
					if (x < 0)x += (i32_t)special_q;
					spq_j = x;
				}
				modulo32 = special_q;
				spq_x = modmul32(spq_i, i_shift);
			}
	
			//fprintf(g_ofile,"# Start %u %u (%d,%d) (%d,%d)\n",
			//special_q,r[root_no],a0,b0,a1,b1);

			// do all of the work (scheduling/setup ?) for this root of this special_q (?)
			{
			  u32_t subsieve_nr;
	  
			  {
				u32_t absa0,absa1,absb0,absb1;
				char a0s,a1s;
				clock_t new_clock;
	    
		#define GET_ABSSIG(abs,sig,arg) if(arg> 0) { abs= (u32_t)arg; sig= '+';} \
				else { abs= (u32_t)(-arg); sig= '-'; }

				GET_ABSSIG(absa0,a0s,a0);
				GET_ABSSIG(absa1,a1s,a1);
				absb0= b0;
				absb1= b1;

				{
					u32_t s;

					// setup packed sieve structures and bounds.
					// possible vectorization candidate?  No, not enough 
					// time is spent here (Sieve-Change).
					// The sieve and tds routines we want to change to 32x
					// processing are associated with smallsieve_aux.
					// instead of packing them, set them at intervals of fbis[s];
					for (s = 0; s < 2; s++) {
						u32_t fbi;
						u16_t* abuf;
						u16_t* ibuf;

						abuf = smallsieve_aux[s];
						ibuf = smallpsieve_aux[s];
						for (fbi = 0; fbi < fbis[s]; fbi++) {
							u32_t aa, bb;
							modulo32 = FB[s][fbi];

							aa = absa0 % FB[s][fbi];
							if (a0s == '-' && aa != 0)aa = FB[s][fbi] - aa;
							bb = absb0 % FB[s][fbi];
							if (proots[s][fbi] != FB[s][fbi]) {
								u32_t x;
								x = modsub32(aa, modmul32(proots[s][fbi], bb));
								if (x != 0) {
									aa = absa1 % FB[s][fbi];
									if (a1s == '-' && aa != 0)aa = FB[s][fbi] - aa;
									bb = absb1 % FB[s][fbi];
									//x = modmul32(asm_modinv32(x), modsub32(modmul32(proots[s][fbi], bb), aa));
									x = modmul32(modinv32(x), modsub32(modmul32(proots[s][fbi], bb), aa));
									abuf[0] = (u16_t)(FB[s][fbi]);
									abuf[1] = (u16_t)x;
									abuf[2] = (u16_t)(FB_logs[s][fbi]);
									abuf += 4;
								}
								else {
									ibuf[0] = (u16_t)(FB[s][fbi]);
									ibuf[1] = (u16_t)(FB_logs[s][fbi]);
									ibuf += 3;
								}
							}
							else {

								if (bb != 0) {
									u32_t x;
									x = modulo32 - bb;
									bb = absb1 % FB[s][fbi];
									abuf[0] = (u16_t)(FB[s][fbi]);
									abuf[1] = (u16_t)(modmul32(modinv32(x), bb));
									//abuf[1] = (u16_t)(modmul32(asm_modinv32(x), bb));
									abuf[2] = (u16_t)(FB_logs[s][fbi]);
									abuf += 4;
								}
								else {
									ibuf[0] = (u16_t)(FB[s][fbi]);
									ibuf[1] = (u16_t)(FB_logs[s][fbi]);
									ibuf += 3;
								}
							}
						}
						smallsieve_auxbound[s][0] = abuf;
						smallpsieve_aux_ub_pow1[s] = ibuf;
					}
				}
	    
				{
					u32_t s;

					for (s = 0; s < 2; s++) {
						u32_t i;
						u16_t* buf;
						u16_t* buf2;
						u16_t* ibuf;

						buf = smallsieve_aux1[s];
						buf2 = x2FB[s];
						ibuf = smallpsieve_aux_ub_pow1[s];
						for (i = 0; i < xFBs[s]; i++) {
							if (xFB[s][i].p == 2) {
								xFBtranslate(buf2, xFB[s] + i);
								buf2 += 4;
							}
							else {
								xFBtranslate(buf, xFB[s] + i);
								if (buf[0] == 1) {
									ibuf[1] = xFB[s][i].l;
									ibuf[0] = xFB[s][i].pp;
									ibuf += 3;
								}
								else buf += 6;
							}
						}
						x2FBs[s] = (buf2 - x2FB[s]) / 4;
						smallpsieve_aux_ub_odd[s] = ibuf;
						smallsieve_aux1_ub_odd[s] = buf;
					}
				}
	    
				{
					u32_t s;

#ifndef MMX_TD
					for (s = 0; s < 2; s++) {
						u32_t i;
						u16_t* x;

#if 1
						for (i = 0, x = smallsieve_aux[s]; x < smallsieve_auxbound[s][0]-32; i+=8, x += 32) {
							u32_t k;

							__m512i xv = _mm512_load_si512(x);
							__m512i vp = _mm512_and_epi64(xv, _mm512_set1_epi64(0xffff));  // isolate primes
							__m512i vr = _mm512_srli_epi64(xv, 16);	// align roots
							__m128i vp128 = _mm512_cvtepi64_epi16(vp);
							__m128i vr128 = _mm512_cvtepi64_epi16(vr);
							__m128i vpr128 = vr128;

							for (k = 0; k < j_per_strip; k++) {
								//smalltdsieve_aux[s][k][i] = r;
								//r = modadd32(r, pr);
								_mm_storeu_epi16(&smalltdsieve_aux[s][k][i], vr128);
								vr128 = _mm_add_epi16(vr128, vpr128);
								__mmask8 m = _mm_cmpge_epu16_mask(vr128, vp128);
								vr128 = _mm_mask_sub_epi16(vr128, m, vr128, vp128);
							}
#ifdef PREINVERT
							{
								//u32_t pinv;
								__m512i vpi = vp;

								//pinv = modulo32;
								//pinv = 2 * pinv - pinv * pinv * modulo32;
								//pinv = 2 * pinv - pinv * pinv * modulo32;
								//pinv = 2 * pinv - pinv * pinv * modulo32;

								__m512i t2 = _mm512_mullo_epi32(vpi, vpi);
								__m512i t1 = _mm512_slli_epi32(vpi, 1);
								__m512i t3 = _mm512_mullo_epi32(t2, vp);
								vpi = _mm512_sub_epi32(t1, t3);

								
								t2 = _mm512_mullo_epi32(vpi, vpi);
								t1 = _mm512_slli_epi32(vpi, 1);
								t3 = _mm512_mullo_epi32(t2, vp);
								vpi = _mm512_sub_epi32(t1, t3);

								
								t2 = _mm512_mullo_epi32(vpi, vpi);
								t1 = _mm512_slli_epi32(vpi, 1);
								t3 = _mm512_mullo_epi32(t2, vp);
								vpi = _mm512_sub_epi32(t1, t3);

								
								t2 = _mm512_mullo_epi32(vpi, vpi);
								t1 = _mm512_slli_epi32(vpi, 1);
								t3 = _mm512_mullo_epi32(t2, vp);
								vpi = _mm512_sub_epi32(t1, t3);

								__m256i vpi256 = _mm512_cvtepi64_epi32(vpi);
								_mm256_storeu_epi32(&smalltd_pi[s][i], vpi256);

								//smalltd_pi[s][i] = 2 * pinv - pinv * pinv * modulo32;
							}

#endif
						}

						for ( ; x < smallsieve_auxbound[s][0]; i++, x += 4) {
							u32_t k, r, pr;

							modulo32 = *x;
							r = x[1];
							pr = r;
							for (k = 0; k < j_per_strip; k++) {
								smalltdsieve_aux[s][k][i] = r;
								r = modadd32(r, pr);
							}
#ifdef PREINVERT
							{
								u32_t pinv;

								pinv = modulo32;
								pinv = 2 * pinv - pinv * pinv * modulo32;
								pinv = 2 * pinv - pinv * pinv * modulo32;
								pinv = 2 * pinv - pinv * pinv * modulo32;
#if 0
								pinv = 2 * pinv - pinv * pinv * modulo32;
#endif
								smalltd_pi[s][i] = 2 * pinv - pinv * pinv * modulo32;
							}

#endif

						}

#else
						for (i = 0, x = smallsieve_aux[s]; x < smallsieve_auxbound[s][0]; i++, x += 4) {
							u32_t k, r, pr;

							modulo32 = *x;
							r = x[1];
							pr = r;
							for (k = 0; k < j_per_strip; k++) {
								smalltdsieve_aux[s][k][i] = r;
								r = modadd32(r, pr);
							}
#ifdef PREINVERT
							{
								u32_t pinv;

								pinv = modulo32;
								pinv = 2 * pinv - pinv * pinv * modulo32;
								pinv = 2 * pinv - pinv * pinv * modulo32;
								pinv = 2 * pinv - pinv * pinv * modulo32;
#if 0
								pinv = 2 * pinv - pinv * pinv * modulo32;
#endif
								smalltd_pi[s][i] = 2 * pinv - pinv * pinv * modulo32;
							}

#endif

						}
#endif

					}
#endif
				}
	    
				{
					u32_t s;

					for (s = 0; s < 2; s++) {
						u16_t* x, * xx, k, pbound, copy_buf[6];

						k = 0;
						pbound = TINY_SIEVE_MIN;
						for (x = smallsieve_aux[s]; x < smallsieve_auxbound[s][0]; x += 4) {
							if (*x > pbound) {
								if (k == 0)smallsieve_tinybound[s] = x;
								else smallsieve_auxbound[s][5 - k] = x;
								k++;
								if (k < 5)pbound = n_i / (5 - k);
								else break;
							}
						}
						while (k < 5)smallsieve_auxbound[s][5 - (k++)] = x;
						for (x = (xx = smallsieve_aux1[s]); x < smallsieve_aux1_ub_odd[s]; x += 6) {
							if (x[0] < TINY_SIEVE_MIN) {
								if (x != xx) {
									memcpy(copy_buf, x, 6 * sizeof(*x));
									memcpy(x, xx, 6 * sizeof(*x));
									memcpy(xx, copy_buf, 6 * sizeof(*x));
								}
								xx += 6;
							}
						}
						smallsieve_tinybound1[s] = xx;
					}
				}

				new_clock = clock();
				sch_clock += new_clock - last_clock;
				last_clock = new_clock;
	    
				{
					u32_t s;

					for (s = 0; s < 2; s++) {
#ifdef SCHEDULING_FUNCTION_CALCULATES_RI
						lasieve_setup(FB[s] + fbis[s], proots[s] + fbis[s], fbi1[s] - fbis[s],
							a0, a1, b0, b1, LPri[s]);
#else
						lasieve_setup(FB[s] + fbis[s], proots[s] + fbis[s], FBsize[s] - fbis[s],
							a0, a1, b0, b1, LPri[s], FBsize[s]);
						//printf("setup complete on side %d, size %u\n", s, FBsize[s]-fbis[s]);
#endif
					}
				}
	    
				{
					u32_t i, k;
					for (i = 0; i < 2; i++) {
						double large_primes_summand;
						tpol(tpoly_f[i], poly_f[i], poldeg[i], a0, a1, b0, b1);
						large_primes_summand = sieve_report_multiplier[i] * FB_maxlog[i];
						if (i == special_q_side)
							large_primes_summand += sieve_multiplier[i] * log(special_q);
						get_sieve_report_bounds(sieve_report_bounds[i], tpoly_f[i], poldeg[i],
							n_srb_i, n_srb_j, 2 * CANDIDATE_SEARCH_STEPS,
							sieve_multiplier[i], large_primes_summand);
					}
				}
	    
				new_clock= clock();
				lsetup_clock+= new_clock-last_clock;
				last_clock= new_clock;
			  }
	  
			  // do all of the work (sieving/td ?) for all oddness types
			  for(oddness_type= 1;oddness_type<4;oddness_type++) 
			  {
				  {
					  u32_t s;
					  for (s = 0; s < 2; s++) {
						  switch (oddness_type) {
							  u16_t* x;
						  case 1:
							  for (x = smallsieve_aux[s]; x < smallsieve_auxbound[s][0]; x += 4) {
								  u32_t p;

								  p = x[0];
								  x[3] = ((i_shift + p) / 2) % p;
							  }

							  for (x = smallsieve_aux1[s]; x < smallsieve_aux1_ub_odd[s]; x += 6) {
								  u32_t p;

								  p = x[0];

								  x[4] = ((i_shift + p) / 2) % p;
								  x[5] = 0;
							  }

							  for (x = smallpsieve_aux[s]; x < smallpsieve_aux_ub_odd[s]; x += 3)
								  x[2] = 0;

							  {
								  u16_t* x, * y, * z;
								  u32_t i;

								  x = smallsieve_aux1_ub_odd[s];
								  y = smallpsieve_aux_ub_odd[s];
								  z = smallsieve_aux2[s];
								  for (i = 0; i < 4 * x2FBs[s]; i += 4) {
									  u32_t p, pr, d, l;
									  u16_t** a;

									  d = x2FB[s][i + 1];
									  if (d == 1)continue;
									  p = x2FB[s][i];
									  pr = x2FB[s][i + 2];
									  l = x2FB[s][i + 3];
									  if (p < 4) {
										  if (p == 1) {
											  *y = d / 2;
											  *(y + 2) = 0;
										  }
										  else {
											  *y = d;
											  *(y + 2) = d / 2;
										  }
										  *(y + 1) = l;
										  y += 3;
										  continue;
									  }
									  p = p / 2;
									  if (p <= MAX_TINY_2POW)a = &z;
									  else a = &x;
									  **a = p;
									  *(1 + *a) = d;
									  *(2 + *a) = pr % p;
									  *(3 + *a) = l;
									  *(4 + *a) = ((i_shift + pr) / 2) % p;
									  *(5 + *a) = d / 2;
									  *a += 6;
								  }
								  smallsieve_aux1_ub[s] = x;
								  smallpsieve_aux_ub[s] = y;
								  smallsieve_aux2_ub[s] = z;
							  }

							  break;
						  case 2:
							  for (x = smallsieve_aux[s]; x < smallsieve_auxbound[s][0]; x += 4) {
								  u32_t p, pr;

								  p = x[0];
								  pr = x[1];
								  x[3] = (pr % 2 == 0 ? ((i_shift + pr) / 2) % p : ((i_shift + pr + p) / 2) % p);
							  }

							  for (x = smallsieve_aux1[s]; x < smallsieve_aux1_ub_odd[s]; x += 6) {
								  u32_t p, d, pr;

								  p = x[0];
								  d = x[1];
								  pr = x[2];

								  x[4] = (pr % 2 == 0 ? ((i_shift + pr) / 2) % p : ((i_shift + pr + p) / 2) % p);
								  x[5] = d / 2;
							  }

							  for (x = smallpsieve_aux[s]; x < smallpsieve_aux_ub_odd[s]; x += 3)
								  x[2] = (x[0]) / 2;

							  {
								  u16_t* x, * y, * z;
								  u32_t i;

								  x = smallsieve_aux1_ub_odd[s];
								  y = smallpsieve_aux_ub_odd[s];
								  z = smallsieve_aux2[s];
								  for (i = 0; i < 4 * x2FBs[s]; i += 4) {
									  u32_t p, pr, d, l;
									  u16_t** a;

									  d = x2FB[s][i + 1];
									  if (d != 1)continue;
									  pr = x2FB[s][i + 2];
									  if (pr % 2 != 0)continue;
									  p = x2FB[s][i];
									  l = x2FB[s][i + 3];
									  if (p < 4) {

										  if (p == 1) {
											  Schlendrian("Use 1=2^0 for sieving?\n");
										  }
										  *y = d;
										  *(y + 1) = l;
										  *(y + 2) = 0;
										  y += 3;
										  continue;
									  }
									  p = p / 2;
									  if (p <= MAX_TINY_2POW)a = &z;
									  else a = &x;
									  **a = p;
									  *(1 + *a) = d;
									  *(2 + *a) = pr % p;
									  *(3 + *a) = l;
									  *(4 + *a) = ((i_shift + pr) / 2) % p;
									  *(5 + *a) = 0;
									  *a += 6;
								  }
								  smallsieve_aux1_ub[s] = x;
								  smallpsieve_aux_ub[s] = y;
								  smallsieve_aux2_ub[s] = z;
							  }

							  break;
						  case 3:
							  for (x = smallsieve_aux[s]; x < smallsieve_auxbound[s][0]; x += 4) {
								  u32_t p, pr;

								  p = x[0];
								  pr = x[1];
								  x[3] = (pr % 2 == 1 ? ((i_shift + pr) / 2) % p : ((i_shift + pr + p) / 2) % p);
							  }

							  for (x = smallsieve_aux1[s]; x < smallsieve_aux1_ub_odd[s]; x += 6) {
								  u32_t p, d, pr;

								  p = x[0];
								  d = x[1];
								  pr = x[2];

								  x[4] = (pr % 2 == 1 ? ((i_shift + pr) / 2) % p : ((i_shift + pr + p) / 2) % p);
								  x[5] = d / 2;
							  }

							  for (x = smallpsieve_aux[s]; x < smallpsieve_aux_ub_odd[s]; x += 3)
								  x[2] = (x[0]) / 2;

							  {
								  u16_t* x, * y, * z;
								  u32_t i;

								  x = smallsieve_aux1_ub_odd[s];
								  y = smallpsieve_aux_ub_odd[s];
								  z = smallsieve_aux2[s];
								  for (i = 0; i < 4 * x2FBs[s]; i += 4) {
									  u32_t p, pr, d, l;
									  u16_t** a;

									  d = x2FB[s][i + 1];
									  if (d != 1)continue;
									  pr = x2FB[s][i + 2];
									  if (pr % 2 != 1)continue;
									  p = x2FB[s][i];
									  l = x2FB[s][i + 3];
									  if (p < 4) {

										  if (p == 1) {
											  Schlendrian("Use 1=2^0 for sieving?\n");
										  }
										  *y = d;
										  *(y + 1) = l;
										  *(y + 2) = 0;
										  y += 3;
										  continue;
									  }
									  p = p / 2;
									  if (p <= MAX_TINY_2POW)a = &z;
									  else a = &x;
									  **a = p;
									  *(1 + *a) = d;
									  *(2 + *a) = pr % p;
									  *(3 + *a) = l;
									  *(4 + *a) = ((i_shift + pr) / 2) % p;
									  *(5 + *a) = 0;
									  *a += 6;
								  }
								  smallsieve_aux1_ub[s] = x;
								  smallpsieve_aux_ub[s] = y;
								  smallsieve_aux2_ub[s] = z;
							  }

							  break;
						  }
					  }
				  }
	    
	    
		#ifdef GCD_SIEVE_BOUND
				  {
					  u32_t i;

					  for (i = 0; i < np_gcd_sieve; i++) {
						  gcd_sieve_buffer[2 * i + 1] = (oddness_type / 2) * (gcd_sieve_buffer[2 * i] / 2);
					  }
				  }
		#endif
	    
				j_offset= 0;
		#ifndef NOSCHED
				{
					u32_t s;
					clock_t new_clock;

					for (s = 0; s < 2; s++) {
						u32_t i;

						for (i = 0; i < n_schedules[s]; i++) {
							u32_t ns;

							ns = schedules[s][i].n_strips;
							if (ns > n_strips)ns = n_strips;
							do_scheduling(schedules[s] + i, ns, oddness_type, s);
							schedules[s][i].current_strip = 0;
						}
					}
		#ifdef GATHER_STAT
					new_clock = clock();
					Schedule_clock += new_clock - last_clock;
					last_clock = new_clock;
		#endif
				}
		#else 
		#define BADSCHED
		#endif
	    
	    
		#ifdef ZSS_STAT
				nss+= n_strips;
		#endif

				// do all of the work (for this oddness type?)
				for (subsieve_nr = 0; subsieve_nr < n_strips;
					subsieve_nr++, j_offset += j_per_strip)
				{
					u16_t s, stepno;
		#ifdef USE_MEDSCHED
		#ifndef NOSCHED
					for (s = 0; s < 2; s++) {
						u32_t ll, * sched, * ri;

						if (n_medsched_pieces[s] == 0)continue;
						for (ll = 0, sched = (u32_t*)med_sched[s][0], ri = LPri[s];
							ll < n_medsched_pieces[s]; ll++) {
							//printf("medsched %d from ij %u to %u\n", ll, 
							//	medsched_fbi_bounds[s][ll], medsched_fbi_bounds[s][ll+1]);
							ri = medsched(ri, current_ij[s] + medsched_fbi_bounds[s][ll],
								current_ij[s] + medsched_fbi_bounds[s][ll + 1], &sched,
								medsched_fbi_bounds[s][ll], j_offset == 0 ? oddness_type : 0, FBsize[s]);
							med_sched[s][ll + 1] = (u16_t*)sched;
						}
					}
		#endif


					{
						clock_t new_clock;
						new_clock = clock();
						medsched_clock += new_clock - last_clock;
						last_clock = new_clock;
					}
		#endif
					for (s = first_sieve_side, stepno = 0; stepno < 2; stepno++, s = 1 - s)
					{
						clock_t new_clock, clock_diff;

						{
							u32_t j;
							u16_t* x;

							// do all sieving with the tinyest primes into a small buffer,
							// then copy that buffer over the sieve interval.
							for (x = smallsieve_aux[s], j = 0; x < smallsieve_tinybound[s]; x += 4, j++) {
								tinysieve_curpos[j] = x[3];
							}
							for (j = 0; j < j_per_strip; j++) {
								unsigned char* si_ub;
								bzero(tiny_sieve_buffer, TINY_SIEVEBUFFER_SIZE);
								si_ub = tiny_sieve_buffer + TINY_SIEVEBUFFER_SIZE;
								{
									u16_t* x;

									for (x = smallsieve_aux[s]; x < smallsieve_tinybound[s]; x += 4) {
										u32_t p, r, pr;
										unsigned char l, * si;

										p = x[0];
										pr = x[1];
										l = x[2];
										r = x[3];
										si = tiny_sieve_buffer + r;
										while (si < si_ub) {
											*si += l;
											si += p;
										}
										r = r + pr;
										if (r >= p)r = r - p;
										x[3] = r;
									}
								}

								{
									u16_t* x;

									for (x = smallsieve_aux2[s]; x < smallsieve_aux2_ub[s]; x += 6) {
										u32_t p, r, pr, d, d0;
										unsigned char l, * si;

										p = x[0];
										d = x[1];
										pr = x[2];
										l = x[3];
										r = x[4];

										d0 = x[5];
										if (d0 > 0) {
											x[5]--;
											continue;
										}
										si = tiny_sieve_buffer + r;
										while (si < si_ub) {
											*si += l;
											si += p;
										}
										r = r + pr;
										if (r >= p)r = r - p;
										x[4] = r;
										x[5] = d - 1;
									}
								}

								{
									u16_t* x;

									for (x = smallsieve_aux1[s]; x < smallsieve_tinybound1[s]; x += 6) {
										u32_t p, r, pr, d, d0;
										unsigned char l, * si;

										p = x[0];
										d = x[1];
										pr = x[2];
										l = x[3];
										r = x[4];

										d0 = x[5];
										if (d0 > 0) {
											x[5]--;
											continue;
										}
										si = tiny_sieve_buffer + r;
										while (si < si_ub) {
											*si += l;
											si += p;
										}
										r = r + pr;
										if (r >= p)r = r - p;
										x[4] = r;
										x[5] = d - 1;
									}
								}

								{
									unsigned char* si;

									si = sieve_interval + (j << i_bits);
									si_ub = sieve_interval + ((j + 1) << i_bits);
									while (si + TINY_SIEVEBUFFER_SIZE < si_ub) {
										memcpy(si, tiny_sieve_buffer, TINY_SIEVEBUFFER_SIZE);
										si += TINY_SIEVEBUFFER_SIZE;
									}
									memcpy(si, tiny_sieve_buffer, si_ub - si);
								}

							}
							for (x = smallsieve_aux[s], j = 0; x < smallsieve_tinybound[s]; x += 4, j++) {
								x[3] = tinysieve_curpos[j];
							}
						}

		#ifdef ZSS_STAT
						if (s == 1 && ncand == 0)
							nzss[0]++;
		#endif
						new_clock = clock();
						clock_diff = new_clock - last_clock;
						si_clock[s] += clock_diff;
						sieve_clock += clock_diff;
						last_clock = new_clock;


				


		#if defined( ASM_LINESIEVER) && !defined(AVX512_SIEVE4)
						slinie(smallsieve_tinybound[s], smallsieve_auxbound[s][4], sieve_interval);
		#else
						{
							u16_t* x;

#ifdef AVX512_SIEVE4
							// Sadly, this isn't any faster... it seemed like
							// a good idea and helps in the SoE :(
							if (!all_initialized_s[s]) {
								for (x = smallsieve_tinybound[s]; x < smallsieve_auxbound[s][4]; x += 4) {
									u32_t p;
									unsigned char l;
									uint64_t** rtab;

									p = x[0];	// prime
									l = x[2];	// log

									// there isn't any small-prime variation in the siever
									// that I can see and as a result, we sieve with some
									// very small primes.  The idea in this loop is to 
									// do an initial sieve in a byte-mask, and then stream
									// out writes of rotations of this mask, like
									// in SoE presieving.
									if (p > presieve_bound_p) break;

									if ((s == 0) && !initialized_rtab_s0[rtab_lookup[p]])
									{
										int j, k;
										rtab = rtab_s0[rtab_lookup[p]];

										//printf("uint64_t rtab_s0[%d][%d][8] = {\n", rtab_lookup[p], p);
										for (j = 0; j < p; j++)
										{
											for (k = 0; k < 8; k++)
												rtab[j][k] = 0;

											//printf("{");
											// make the 64-byte vector for the prime/offset
											for (k = j; k < 64; k += p)
											{
												rtab[j][k / 8] |= ((uint64_t)l << (8 * (k & 7)));
											}
											//for (k = 0; k < 7; k++)
											//	printf("0x%016lxULL, ", rtab[j][k]);
											//printf("0x%016lxULL},\n", rtab[j][k]);
										}
										//printf("};\n");
										initialized_rtab_s0[rtab_lookup[p]] = 1;
									}
									else if ((s == 1) && !initialized_rtab_s1[rtab_lookup[p]])
									{
										int j, k;
										rtab = rtab_s1[rtab_lookup[p]];
										//printf("uint64_t rtab_s1[%d][%d][8] = {\n", rtab_lookup[p], p);
										for (j = 0; j < p; j++)
										{
											for (k = 0; k < 8; k++)
												rtab[j][k] = 0;

											//printf("{");
											// make the 64-byte vector for the prime/offset
											for (k = j; k < 64; k += p)
											{
												rtab[j][k / 8] |= ((uint64_t)l << (8 * (k & 7)));
											}
											//for (k = 0; k < 7; k++)
											//	printf("0x%016lxULL, ", rtab[j][k]);
											//printf("0x%016lxULL},\n", rtab[j][k]);
										}
										//printf("};\n");
										initialized_rtab_s1[rtab_lookup[p]] = 1;
									}

								}
								all_initialized_s[s] = 1;
							}

							for (x = smallsieve_tinybound[s]; x < smallsieve_auxbound[s][4]; x += 16) {
								unsigned char *y;
								uint64_t** rtab1;
								uint64_t** rtab2;
								uint64_t** rtab3;
								uint64_t** rtab4;

								if (x[12] > presieve_bound_p) break;

								if (s == 0)
								{
									rtab1 = rtab_s0[rtab_lookup[x[0]]];
									rtab2 = rtab_s0[rtab_lookup[x[4]]];
									rtab3 = rtab_s0[rtab_lookup[x[8]]];
									rtab4 = rtab_s0[rtab_lookup[x[12]]];
								}
								else
								{
									rtab1 = rtab_s1[rtab_lookup[x[0]]];
									rtab2 = rtab_s1[rtab_lookup[x[4]]];
									rtab3 = rtab_s1[rtab_lookup[x[8]]];
									rtab4 = rtab_s1[rtab_lookup[x[12]]];
								}
								int step1 = steps[x[ 0]];
								int step2 = steps[x[ 4]];
								int step3 = steps[x[ 8]];
								int step4 = steps[x[12]];
								u32_t r1 = x[3];
								u32_t r2 = x[7];
								u32_t r3 = x[11];
								u32_t r4 = x[15];
								for (y = sieve_interval; y < sieve_interval + L1_SIZE; y += n_i) {
									unsigned char* yy;

									int rr1 = r1;
									int rr2 = r2;
									int rr3 = r3;
									int rr4 = r4;
									for (yy = y; yy < (y + n_i); yy += 64)
									{
										// load all of the presieved 64B chuncks for these
										// four primes as well as this chunck of the sieve interval.
										__m512i vl1 = _mm512_load_si512(rtab1[rr1]);
										__m512i vl2 = _mm512_load_si512(rtab2[rr2]);
										__m512i vl3 = _mm512_load_si512(rtab3[rr3]);
										__m512i vl4 = _mm512_load_si512(rtab4[rr4]);
										__m512i vy = _mm512_load_si512(yy);
										// combine the presieved chuncks and write it out.
										vl1 = _mm512_add_epi8(vl3, vl1);
										vl2 = _mm512_add_epi8(vl4, vl2);
										vy = _mm512_add_epi8(vy, vl1);
										vy = _mm512_add_epi8(vy, vl2);
										_mm512_store_si512(yy, vy);

										// next rotation
										rr1 += step1;
										rr2 += step2;
										rr3 += step3;
										rr4 += step4;
										rr1 = (rr1 >= 64) ? rr1 - 64 : rr1 + x[ 0] - 64;
										rr2 = (rr2 >= 64) ? rr2 - 64 : rr2 + x[ 4] - 64;
										rr3 = (rr3 >= 64) ? rr3 - 64 : rr3 + x[ 8] - 64;
										rr4 = (rr4 >= 64) ? rr4 - 64 : rr4 + x[12] - 64;
									}

									r1 = r1 + x[1];
									r2 = r2 + x[5];
									r3 = r3 + x[9];
									r4 = r4 + x[13];
									if (r1 >= x[ 0])r1 = r1 - x[ 0];
									if (r2 >= x[ 4])r2 = r2 - x[ 4];
									if (r3 >= x[ 8])r3 = r3 - x[ 8];
									if (r4 >= x[12])r4 = r4 - x[12];
								}
							}

							_mm256_zeroupper();

							// primes here are less than n_i/4, we unroll the sieve 4x 
							// over those chunks.
							//x = smallsieve_tinybound[s]
							for ( ; x < smallsieve_auxbound[s][4]; x += 4) {
#else
							for (x = smallsieve_tinybound[s]; x < smallsieve_auxbound[s][4]; x += 4) {
#endif
								u32_t p, r, pr;
								unsigned char l, * y;

								p = x[0];	// prime
								pr = x[1];	// prime root/recurrence?
								l = x[2];	// log
								r = x[3];	// starting root/recurrence?

								for (y = sieve_interval; y < sieve_interval + L1_SIZE; y += n_i) {
									unsigned char* yy, * yy_ub;

									yy_ub = y + n_i - 3 * p;
									for (yy = y + r; yy < yy_ub; yy = yy + 4 * p) {
										*(yy) += l;
										*(yy + p) += l;
										*(yy + 2 * p) += l;
										*(yy + 3 * p) += l;
									}
									while (yy < y + n_i) {
										*(yy) += l;
										yy += p;
									}
									r = r + pr;
									if (r >= p)r = r - p;
								}
							}
						}
		#endif


		#if defined( ASM_LINESIEVER3)  && !defined(AVX512_SIEVE3)
						slinie3(smallsieve_auxbound[s][4], smallsieve_auxbound[s][3], sieve_interval);
		#else
						{
							u16_t* x;

#if defined(AVX512_SIEVE3)
							// in this interval primes hit once; after that we have to check.
							// do 8 at a time.
							__m512i vni = _mm512_set1_epi16(n_i);
							for (x = smallsieve_auxbound[s][4]; x < smallsieve_auxbound[s][3] - 32; x += 32)
							{
								// for reference, sieve info is packed into x thus:
								//p = x[0];
								//pr = x[1];
								//l = x[2];
								//r = x[3];
								unsigned char* y, l = x[30];	// assume these 8 primes have the same log
								__m512i xv = _mm512_load_si512(x);
								__m512i vp = _mm512_slli_epi64(xv, 48);			// align these with r
								__m512i vpr = _mm512_slli_epi64(xv, 32);		// align these with r
								vpr = _mm512_and_epi64(vpr, _mm512_set1_epi64(0xffff000000000000ULL));
								uint16_t xm[32];
								__mmask32 m;
								for (y = sieve_interval; y < sieve_interval + L1_SIZE; y += n_i) {

									_mm512_store_epi64(xm, xv);

									*(y + xm[3]) += l;
									*(y + xm[7]) += l;
									*(y + xm[11]) += l;
									*(y + xm[15]) += l;
									*(y + xm[19]) += l;
									*(y + xm[23]) += l;
									*(y + xm[27]) += l;
									*(y + xm[31]) += l;

									__m512i xv2 = _mm512_add_epi16(xv, vp);

									*(y + xm[0] + xm[3]) += l;
									*(y + xm[4] + xm[7]) += l;
									*(y + xm[8] + xm[11]) += l;
									*(y + xm[12] + xm[15]) += l;
									*(y + xm[16] + xm[19]) += l;
									*(y + xm[20] + xm[23]) += l;
									*(y + xm[24] + xm[27]) += l;
									*(y + xm[28] + xm[31]) += l;

									xv2 = _mm512_add_epi16(xv2, vp);

									*(y + 2*xm[0] + xm[3]) += l;
									*(y + 2*xm[4] + xm[7]) += l;
									*(y + 2*xm[8] + xm[11]) += l;
									*(y + 2* xm[12] + xm[15]) += l;
									*(y + 2* xm[16] + xm[19]) += l;
									*(y + 2* xm[20] + xm[23]) += l;
									*(y + 2* xm[24] + xm[27]) += l;
									*(y + 2* xm[28] + xm[31]) += l;

									xv2 = _mm512_add_epi16(xv2, vp);
									m = _mm512_mask_cmplt_epu16_mask(0x88888888, xv2, vni);

									// as primes get larger, fewer of them will hit twice,
									// so it's faster to go directly to the set bits.
									while (m > 0)
									{
										int id = _tzcnt_u32(m) / 4;
										*(y + 3*xm[id * 4] + xm[id * 4 + 3]) += l;
										m = _blsr_u32(m);
									}

									// now update r
									xv = _mm512_mask_add_epi16(xv, 0x88888888, xv, vpr);
									m = _mm512_mask_cmpge_epu16_mask(0x88888888, xv, vp);
									xv = _mm512_mask_sub_epi16(xv, m, xv, vp);
								}
							}

							//x = smallsieve_auxbound[s][4]
							for (; x < smallsieve_auxbound[s][3]; x += 4) {

#else
							for (x = smallsieve_auxbound[s][4]; x < smallsieve_auxbound[s][3]; x += 4) {
#endif
								u32_t p, r, pr;
								unsigned char l, * y;

								p = x[0];
								pr = x[1];
								l = x[2];
								r = x[3];
								for (y = sieve_interval; y < sieve_interval + L1_SIZE; y += n_i) {
									unsigned char* yy;

									yy = y + r;
									*(yy) += l;
									*(yy + p) += l;
									*(yy + 2 * p) += l;
									yy += 3 * p;
									if (yy < y + n_i)*(yy) += l;
									r = r + pr;
									if (r >= p)r = r - p;
								}
		#if 0
								x[3] = r;
		#endif
							}
						}
		#endif


		#if defined( ASM_LINESIEVER2)  && !defined(AVX512_SIEVE2)
						slinie2(smallsieve_auxbound[s][3], smallsieve_auxbound[s][2], sieve_interval);
		#else
						{
							u16_t* x;


#if defined(AVX512_SIEVE2)
							// in this interval primes hit 3x; after that we have to check.
							// do 8 at a time.
							__m512i vni = _mm512_set1_epi16(n_i);
							for (x = smallsieve_auxbound[s][3]; x < smallsieve_auxbound[s][2] - 32; x += 32)
							{
								// for reference, sieve info is packed into x thus:
								unsigned char* y, l = x[30];	// assume these 8 primes have the same log
								__m512i xv = _mm512_load_si512(x);
								__m512i vp = _mm512_slli_epi64(xv, 48);			// align these with r
								__m512i vpr = _mm512_slli_epi64(xv, 32);		// align these with r
								vpr = _mm512_and_epi64(vpr, _mm512_set1_epi64(0xffff000000000000ULL));
								uint16_t xm[32];
								__mmask32 m;
								for (y = sieve_interval; y < sieve_interval + L1_SIZE; y += n_i) {

									_mm512_store_epi64(xm, xv);

									*(y + xm[3]) += l;
									*(y + xm[7]) += l;
									*(y + xm[11]) += l;
									*(y + xm[15]) += l;
									*(y + xm[19]) += l;
									*(y + xm[23]) += l;
									*(y + xm[27]) += l;
									*(y + xm[31]) += l;

									__m512i xv2 = _mm512_add_epi16(xv, vp);
									_mm512_store_epi64(xm, xv2);

									*(y + xm[3]) += l;
									*(y + xm[7]) += l;
									*(y + xm[11]) += l;
									*(y +  xm[15]) += l;
									*(y +  xm[19]) += l;
									*(y +  xm[23]) += l;
									*(y +  xm[27]) += l;
									*(y +  xm[31]) += l;

									xv2 = _mm512_add_epi16(xv2, vp);
									_mm512_store_epi64(xm, xv2);

									m = _mm512_mask_cmplt_epu16_mask(0x88888888, xv2, vni);

									// as primes get larger, fewer of them will hit twice,
									// so it's faster to go directly to the set bits.
									while (m > 0)
									{
										int id = _tzcnt_u32(m) / 4;
										*(y + xm[id * 4 + 3]) += l;
										m = _blsr_u32(m);
									}

									// now update r
									xv = _mm512_mask_add_epi16(xv, 0x88888888, xv, vpr);
									m = _mm512_mask_cmpge_epu16_mask(0x88888888, xv, vp);
									xv = _mm512_mask_sub_epi16(xv, m, xv, vp);
								}
							}

							//x = smallsieve_auxbound[s][3]
							for (; x < smallsieve_auxbound[s][2]; x += 4) {

#else
							for (x = smallsieve_auxbound[s][3]; x < smallsieve_auxbound[s][2]; x += 4) {
#endif
								u32_t p, r, pr;
								unsigned char l, * y;

								p = x[0];
								pr = x[1];
								l = x[2];
								r = x[3];

								for (y = sieve_interval; y < sieve_interval + L1_SIZE; y += n_i) {
									unsigned char* yy;

									yy = y + r;
									*(yy) += l;
									*(yy + p) += l;
									yy += 2 * p;
									if (yy < y + n_i)*(yy) += l;
									r = r + pr;
									if (r >= p)r = r - p;
								}
		#if 0
								x[3] = r;
		#endif
							}
						}
		#endif


		#if defined( ASM_LINESIEVER1)  && !defined(AVX512_SIEVE1)
						slinie1(smallsieve_auxbound[s][2], smallsieve_auxbound[s][1], sieve_interval);
		#else
						{
							u16_t* x;

#if defined(AVX512_SIEVE1)
							// in this interval primes hit once; after that we have to check.
							// do 8 at a time.
							// possible further improvement: restructure x for contiguous
							// p's, pr's, r's, l's, so we can do 32x at a time.  This touches
							// a lot of the code...
							__m512i vni = _mm512_set1_epi16(n_i);
							for (x = smallsieve_auxbound[s][2]; x < smallsieve_auxbound[s][1] - 32; x += 32)
							{
								// for reference, sieve info is packed into x thus:
								//p = x[0];
								//pr = x[1];
								//l = x[2];
								//r = x[3];
								unsigned char* y, l = x[30];	// assume these 8 primes have the same log
								__m512i xv = _mm512_load_si512(x);
								__m512i vp = _mm512_slli_epi64(xv, 48);			// align these with r
								__m512i vpr = _mm512_slli_epi64(xv, 32);		// align these with r
								vpr = _mm512_and_epi64(vpr, _mm512_set1_epi64(0xffff000000000000ULL));
								uint16_t xm[32];
								__mmask32 m;
								for (y = sieve_interval; y < sieve_interval + L1_SIZE; y += n_i) {

									_mm512_store_epi64(xm, xv);

									*(y + xm[ 3]) += l;
									*(y + xm[ 7]) += l;
									*(y + xm[11]) += l;

									__m512i xv2 = _mm512_add_epi16(xv, vp);

									*(y + xm[15]) += l;
									*(y + xm[19]) += l;
									*(y + xm[23]) += l;

									m = _mm512_mask_cmplt_epu16_mask(0x88888888, xv2, vni);

									*(y + xm[27]) += l;
									*(y + xm[31]) += l;

									// as primes get larger, fewer of them will hit twice,
									// so it's faster to go directly to the set bits.
									while (m > 0)
									{
										int id = _tzcnt_u32(m) / 4;
										*(y + xm[id*4] + xm[id*4+3]) += l;
										m = _blsr_u32(m);
									}

									// now update r
									xv = _mm512_mask_add_epi16(xv, 0x88888888, xv, vpr);
									m = _mm512_mask_cmpge_epu16_mask(0x88888888, xv, vp);
									xv = _mm512_mask_sub_epi16(xv, m, xv, vp);
								}
							}

							//x = smallsieve_auxbound[s][2]
							for (; x < smallsieve_auxbound[s][1]; x += 4) {

#else
							for (x = smallsieve_auxbound[s][2]; x < smallsieve_auxbound[s][1]; x += 4) {
#endif
								u32_t p, r, pr;
								unsigned char l, * y;

								p = x[0];
								pr = x[1];
								l = x[2];
								r = x[3];

								for (y = sieve_interval; y < sieve_interval + L1_SIZE; y += n_i) {
									unsigned char* yy;

									yy = y + r;
									*(yy) += l;
									yy += p;
									if (yy < y + n_i)*(yy) += l;
									r = r + pr;
									if (r >= p)r = r - p;
								}
		#if 0
								x[3] = r;
		#endif
							}
						}
		#endif

		#if 0
						{
							u16_t* x;

							for (x = smallsieve_auxbound[s][1]; x < smallsieve_auxbound[s][0]; x += 4) {
								u32_t p, r, pr;
								unsigned char l, * y;

								p = x[0];
								pr = x[1];
								l = x[2];
								r = x[3];
								for (y = sieve_interval; y < sieve_interval + L1_SIZE; y += n_i) {
									if (r < n_i)*(y + r) += l;
									r = r + pr;
									if (r >= p)r = r - p;
								}
		#if 0
								x[3] = r;
		#endif
							}
						}
		#endif

		#if 1
						{
							u16_t* x;

							// these are the small prime powers...
							for (x = smallsieve_tinybound1[s]; x < smallsieve_aux1_ub[s]; x += 6) {
								u32_t p, r, pr, d, d0;
								unsigned char l;

								p = x[0];
								d = x[1];
								pr = x[2];
								l = x[3];
								r = x[4];

								for (d0 = x[5]; d0 < j_per_strip; d0 += d) {
									unsigned char* y, * yy, * yy_ub;

									y = sieve_interval + (d0 << i_bits);
									yy_ub = y + n_i - 3 * p;
									for (yy = y + r; yy < yy_ub; yy = yy + 4 * p) {
										*(yy) += l;
										*(yy + p) += l;
										*(yy + 2 * p) += l;
										*(yy + 3 * p) += l;
									}
									while (yy < y + n_i) {
										*(yy) += l;
										yy += p;
									}
									r = r + pr;
									if (r >= p)r = r - p;
								}
								x[4] = r;
								x[5] = d0 - j_per_strip;
							}
						}
		#endif

		#if 1
						{
							u16_t* x;

							bzero(horizontal_sievesums, j_per_strip);
							for (x = smallpsieve_aux[s]; x < smallpsieve_aux_ub[s]; x += 3) {
								u32_t p, d;
								unsigned char l;

								p = x[0];
								l = x[1];
								d = x[2];
		#if I_bits==L1_BITS
								// j_per_strip = 2 here, and p is rarely 0 (which is a bug)
								if (d < 2) horizontal_sievesums[d] += l;
								if (p == 0) x[0] = USHRT_MAX - 1; // otherwise will crash in trial_divide()
								else if ((d += p) < 2) horizontal_sievesums[d] += l;
		#else
		#if I_bits<L1_BITS
								while (d < j_per_strip) {
									horizontal_sievesums[d] += l;
									d += p;
								}
		#else
								// j_per_strip = 1 here, and p is rarely 0 (which is a bug) 
								if (d == 0)
									*horizontal_sievesums += l;
		#endif
		#endif
		#if 0
								x[2] = d - j_per_strip;
		#endif
							}
						}
		#else
						bzero(horizontal_sievesums, j_per_strip);
		#endif

						new_clock = clock();
						clock_diff = new_clock - last_clock;
						s1_clock[s] += clock_diff;
						sieve_clock += clock_diff;
						last_clock = new_clock;

		#ifdef BADSCHED
						ncand = 0;
						continue;
		#endif
		#ifndef MEDSCHE_SI_OFFS
		#ifdef BIGENDIAN
		#define MEDSCHED_SI_OFFS 1
		#else
		#define MEDSCHED_SI_OFFS 0
		#endif
		#endif
		#ifdef ASM_SCHEDSIEVE1
						schedsieve(medsched_logs[s], n_medsched_pieces[s], med_sched[s], sieve_interval);
		#else
						{
							u32_t l;

							for (l = 0; l < n_medsched_pieces[s]; l++) {
								unsigned char x;
								u16_t* schedule_ptr;

								x = medsched_logs[s][l];
		#ifdef ASM_SCHEDSIEVE
								schedsieve(x, sieve_interval, med_sched[s][l], med_sched[s][l + 1]);
		#else
								//for reference: #define SE_SIZE 2
								for (schedule_ptr = med_sched[s][l] + MEDSCHED_SI_OFFS;
									schedule_ptr + 3 * SE_SIZE < med_sched[s][l + 1];
									schedule_ptr += 4 * SE_SIZE) {
									sieve_interval[*schedule_ptr] += x;
									sieve_interval[*(schedule_ptr + SE_SIZE)] += x;
									sieve_interval[*(schedule_ptr + 2 * SE_SIZE)] += x;
									sieve_interval[*(schedule_ptr + 3 * SE_SIZE)] += x;
								}
								for (;
									schedule_ptr < med_sched[s][l + 1]; schedule_ptr += SE_SIZE)
									sieve_interval[*schedule_ptr] += x;
		#endif
							}
						}
		#endif

						new_clock = clock();
						clock_diff = new_clock - last_clock;
						s2_clock[s] += clock_diff;
						sieve_clock += clock_diff;
						last_clock = new_clock;
		#ifndef SCHED_SI_OFFS
		#ifdef BIGENDIAN
		#define SCHED_SI_OFFS 1
		#else
		#define SCHED_SI_OFFS 0
		#endif
		#endif
						{
							u32_t j;

							for (j = 0; j < n_schedules[s]; j++) {
								if (schedules[s][j].current_strip == schedules[s][j].n_strips) {
									u32_t ns;

									ns = schedules[s][j].n_strips;
									if (ns > n_strips - subsieve_nr)ns = n_strips - subsieve_nr;
									do_scheduling(schedules[s] + j, ns, 0, s);
									schedules[s][j].current_strip = 0;
								}
							}

		#ifdef GATHER_STAT
							new_clock = clock();
							Schedule_clock += new_clock - last_clock;
							last_clock = new_clock;
		#endif

							for (j = 0; j < n_schedules[s]; j++) {
		#ifdef ASM_SCHEDSIEVE1
								u32_t i, k;

								k = schedules[s][j].current_strip;
								for (i = 0; i <= schedules[s][j].n_pieces; i++) {
									schedbuf[i] = schedules[s][j].schedule[i][k];
								}
								schedsieve(schedules[s][j].schedlogs, schedules[s][j].n_pieces,
									schedbuf, sieve_interval);
		#else
								u32_t l, k;

								k = schedules[s][j].current_strip;
								l = 0;
								while (l < schedules[s][j].n_pieces) {
									unsigned char x;
									u16_t* schedule_ptr, * sptr_ub;

									x = schedules[s][j].schedlogs[l];
									schedule_ptr = schedules[s][j].schedule[l][k] + SCHED_SI_OFFS;
									while (l < schedules[s][j].n_pieces)
										if (schedules[s][j].schedlogs[++l] != x)break;
									sptr_ub = schedules[s][j].schedule[l][k];

		#ifdef ASM_SCHEDSIEVE
									// this is faster, but not by much
									schedsieve(x, sieve_interval, schedule_ptr, sptr_ub);

		#else
									//for reference: #define SE_SIZE 2
									//sieve_interval is a byte array
									//schedule_ptr is a pointer to a list of 
									//random offsets within the sieve_interval.
									//we need to increment each of those offsets by x.
									//so this is a gather/scatter and probably not
									// a good candidate for vectorization.
									// e.g.:
									// sptr_ub=4101207624
									//  * sched_ptr = 9155, +2 = 8861, +4 = 25560, +6 = 10868
									//  * sched_ptr = 31581, +2 = 24146, +4 = 17921, +6 = 1956
									//  * sched_ptr = 23119, +2 = 29339, +4 = 14577, +6 = 7741
									//  * sched_ptr = 22200, +2 = 20647, +4 = 32225, +6 = 27748
									//  * sched_ptr = 12575, +2 = 19747, +4 = 3901, +6 = 12502
									//  * sched_ptr = 17874, +2 = 1026, +4 = 31833, +6 = 4227
									//  sptr_ub = 0x7ffff4737648
									//  sched_ptr = 0x7ffff4737010, +2 = 0x7ffff4737014, +4 = 0x7ffff4737018, +6 = 0x7ffff473701c
									//  sched_ptr = 0x7ffff4737020, +2 = 0x7ffff4737024, +4 = 0x7ffff4737028, +6 = 0x7ffff473702c
									//  sched_ptr = 0x7ffff4737030, +2 = 0x7ffff4737034, +4 = 0x7ffff4737038, +6 = 0x7ffff473703c
									//  sched_ptr = 0x7ffff4737040, +2 = 0x7ffff4737044, +4 = 0x7ffff4737048, +6 = 0x7ffff473704c
									//  sched_ptr = 0x7ffff4737050, +2 = 0x7ffff4737054, +4 = 0x7ffff4737058, +6 = 0x7ffff473705c
									//  sched_ptr = 0x7ffff4737060, +2 = 0x7ffff4737064, +4 = 0x7ffff4737068, +6 = 0x7ffff473706c
									//  sched_ptr = 0x7ffff4737070, +2 = 0x7ffff4737074, +4 = 0x7ffff4737078, +6 = 0x7ffff473707c
									//
									//printf("sptr_ub=%p, sched_ptr=%p\n", sptr_ub, schedule_ptr);
									while (schedule_ptr + 3 * SE_SIZE < sptr_ub) {
										sieve_interval[*schedule_ptr] += x;
										sieve_interval[*(schedule_ptr + SE_SIZE)] += x;
										sieve_interval[*(schedule_ptr + 2 * SE_SIZE)] += x;
										sieve_interval[*(schedule_ptr + 3 * SE_SIZE)] += x;
										//printf("sched_ptr=%p, +2=%p, +4=%p, +6=%p\n",
										  //  schedule_ptr, (schedule_ptr + SE_SIZE),
										  //  (schedule_ptr + 2 * SE_SIZE), (schedule_ptr + 3 * SE_SIZE));
										schedule_ptr += 4 * SE_SIZE;
									}
									while (schedule_ptr < sptr_ub) {
										sieve_interval[*schedule_ptr] += x;
										schedule_ptr += SE_SIZE;
									}
									//exit(1);
		#endif
								}
		#endif
							}
						}

		#if 0
						dumpsieve(j_offset, s);
		#endif
						new_clock = clock();
						clock_diff = new_clock - last_clock;
						sieve_clock += clock_diff;
						s3_clock[s] += clock_diff;
						last_clock = new_clock;

						if (s == first_sieve_side) {
		#ifdef GCD_SIEVE_BOUND
							gcd_sieve();
		#endif
		#ifdef ASM_SEARCH0
							{
								unsigned char* srbs;
								u32_t i;
								srbs = sieve_report_bounds[s][j_offset / CANDIDATE_SEARCH_STEPS];
								ncand = lasieve_search0(sieve_interval, horizontal_sievesums,
									horizontal_sievesums + j_per_strip,
									srbs, srbs + n_i / CANDIDATE_SEARCH_STEPS,
									cand, fss_sv);
								for (i = 0; i < ncand; i++)fss_sv[i] += horizontal_sievesums[cand[i] >> i_bits];
							}
		#else
							{
								unsigned char* srbs;
								u32_t i;

								srbs = sieve_report_bounds[s][j_offset / CANDIDATE_SEARCH_STEPS];
								ncand = 0;
								for (i = 0; i < n_i; i += CANDIDATE_SEARCH_STEPS) {
									unsigned char st;
									u32_t j;

									st = *(srbs++);
									for (j = 0; j < j_per_strip; j++) {
										unsigned char* i_o, * i_max, st1;

										i_o = sieve_interval + (j << i_bits) + i;
										i_max = i_o + CANDIDATE_SEARCH_STEPS;
										if (st <= horizontal_sievesums[j]) {
											while (i_o < i_max) {
												cand[ncand] = i_o - sieve_interval;
												fss_sv[ncand++] = *(i_o++) + horizontal_sievesums[j];
											}
											continue;
										}
										st1 = st - horizontal_sievesums[j];
		#ifndef HAVE_SSIMD
		#ifdef GNFS_CS32 
		#define bc_t unsigned long
		#define BC_MASK 0x80808080
		#else
		#define bc_t unsigned long long
		#define BC_MASK 0x8080808080808080
		#endif
										{
											if (st1 < 0x80) {
												bc_t bc, * i_oo;

												bc = st1;
												bc = (bc << 8) | bc;
												bc = (bc << 16) | bc;
		#ifndef GNFS_CS32
												bc = (bc << 32) | bc;
		#endif
												bc = BC_MASK - bc;
												for (i_oo = (bc_t*)i_o; i_oo < (bc_t*)i_max; i_oo++) {
													bc_t v = *i_oo;
													if (((v & BC_MASK) | ((v + bc) & BC_MASK)) == 0)continue;
													for (i_o = (unsigned char*)i_oo; i_o < (unsigned char*)(i_oo + 1); i_o++) {
														if (*i_o >= st1) {
															cand[ncand] = i_o - sieve_interval;
															fss_sv[ncand++] = *i_o + horizontal_sievesums[j];

														}
													}
												}
											}
											else {
												bc_t* i_oo;

												for (i_oo = (bc_t*)i_o; i_oo < (bc_t*)i_max; i_oo++) {
													if ((*i_oo & BC_MASK) == 0)continue;
													for (i_o = (unsigned char*)i_oo; i_o < (unsigned char*)(i_oo + 1); i_o++) {
														if (*i_o >= st1) {
															cand[ncand] = i_o - sieve_interval;
															fss_sv[ncand++] = *i_o + horizontal_sievesums[j];

														}
													}
												}
											}
										}
		#else
										{
											unsigned long long x;

											x = st1 - 1;
											x |= x << 8;
											x |= x << 16;
											x |= x << 32;
											while (i_o < i_max) {
												asm volatile("movq (%%eax),%%mm7\n"
													"1:\n"
													"movq (%%esi),%%mm1\n"
													"movq 8(%%esi),%%mm0\n"
													"pmaxub 16(%%esi),%%mm1\n"
													"pmaxub 24(%%esi),%%mm0\n"
													"pmaxub %%mm7,%%mm1\n"
													"pmaxub %%mm1,%%mm0\n"
													"pcmpeqb %%mm7,%%mm0\n"
													"pmovmskb %%mm0,%%eax\n"
													"cmpl $255,%%eax\n"
													"jnz 2f\n"
													"leal 32(%%esi),%%esi\n"
													"cmpl %%esi,%%edi\n"
													"ja 1b\n"
													"2:\n"
													"emms":"=S"(i_o) : "a"(&x), "S"(i_o),
													"D"(i_max));
												if (i_o < i_max) {
													unsigned char* i_max2 = i_o + 32;
													while (i_o < i_max2) {
														if (*i_o >= st1) {
															cand[ncand] = i_o - sieve_interval;
															fss_sv[ncand++] = *i_o + horizontal_sievesums[j];

														}
														i_o++;
													}
												}
											}
										}
		#endif

									}
								}
							}
		#endif
		#if 0
							{
								char* ofn;
								FILE* of;
								asprintf(&ofn, "cdump.%u.%u.j%u.ot%u", special_q, r[root_no],
									j_offset, oddness_type);
								if ((of = fopen(ofn, "w")) != NULL) {
									u32_t i;
									fprintf(of, "%u candidates\n", ncand);
									for (i = 0; i < ncand; i++)
										fprintf(of, "%u %u\n", cand[i], fss_sv[i]);
									fclose(of);
								}
								else errprintf("Cannot open debug file %s: %m\n", ofn);
								free(ofn);
							}
		#endif

						}
						else
						{
							u32_t i, nc1;
							unsigned char* srbs;
							static u32_t bad_pvl = 0;

							srbs = sieve_report_bounds[s][j_offset / CANDIDATE_SEARCH_STEPS];
							n_prereports += ncand;
							for (i = 0, nc1 = 0; i < ncand; i++) {
								u16_t st_i, t_j, ii, jj, j;
								double pvl;

								j = cand[i] >> i_bits;
		#ifndef DEBUG_SIEVE_REPORT_BOUNDS
								if (sieve_interval[cand[i]] + horizontal_sievesums[j] <
									srbs[(cand[i] & (n_i - 1)) / CANDIDATE_SEARCH_STEPS])
									continue;
		#endif
								jj = j_offset + j;
								ii = cand[i] & (n_i - 1);
								st_i = 2 * ii + (oddness_type == 2 ? 0 : 1);
								t_j = 2 * jj + (oddness_type == 1 ? 0 : 1);
								pvl = log(fabs(rpol_eval(tpoly_f[s], poldeg[s],
									(double)st_i - (double)i_shift, (double)t_j)));
								if (special_q_side == s)
									pvl -= special_q_log;
								pvl *= sieve_multiplier[s];
								pvl -= sieve_report_multiplier[s] * FB_maxlog[s];
								if ((double)(sieve_interval[cand[i]] + horizontal_sievesums[j]) >= pvl) {
		#ifdef DEBUG_SIEVE_REPORT_BOUNDS
									if (sieve_interval[cand[i]] + horizontal_sievesums[j] <
										srbs[(cand[i] & (n_i - 1)) / CANDIDATE_SEARCH_STEPS]) {
										double pvl1;

										pvl = fabs(rpol_eval(tpoly_f[s], poldeg[s],
											(double)st_i - (double)i_shift, (double)t_j));
										fprintf(stderr, "Bad pvl min %u at (%f,%f),spq=%u\npvl: %.5g->",
											bad_pvl++, (double)st_i - (double)i_shift, (double)t_j,
											special_q, pvl);
										pvl = log(pvl);
										fprintf(stderr, "%.3f->", pvl);
										pvl = sieve_multiplier[s] * pvl;
										fprintf(stderr, "%.3f->", pvl);
										if (special_q_side == s)pvl -= sieve_multiplier[s] * special_q_log;
										fprintf(stderr, "%.3f->", pvl);
										pvl -= sieve_report_multiplier[s] * FB_maxlog[s];
										fprintf(stderr, "%.3f\nLower bound was %u sv was %u=%u+%u\n", pvl,
											(u32_t)srbs[(cand[i] & (n_i - 1)) / CANDIDATE_SEARCH_STEPS],
											(u32_t)sieve_interval[cand[i]] + (u32_t)horizontal_sievesums[j],
											(u32_t)sieve_interval[cand[i]], (u32_t)horizontal_sievesums[j]);
									}
		#endif
									fss_sv[nc1] = fss_sv[i];
									cand[nc1++] = cand[i];
								}
							}
							ncand = nc1;
						}

						new_clock = clock();
						clock_diff = new_clock - last_clock;
						sieve_clock += clock_diff;
						cs_clock[s] += clock_diff;
						last_clock = new_clock;
					}


		#ifndef BADSCHED
					trial_divide();
		#endif
					{
						clock_t new_clock;
						new_clock = clock();
						td_clock += new_clock - last_clock;
						last_clock = new_clock;
					}

		#if TDS_MPQS == TDS_BIGSS
		#error "MPQS at BIGSS not yet for serial siever"
				  output_all_tdsurvivors();
		#else
		#if TDS_PRIMALITY_TEST == TDS_BIGSS
		#error "MPQS at BIGSS not yet for serial siever"
				  primality_tests_all();
		#endif
		#endif
				  
				}

		#if TDS_MPQS == TDS_ODDNESS_CLASS
				output_all_tdsurvivors();
		#else
		#if TDS_PRIMALITY_TEST == TDS_ODDNESS_CLASS
				primality_tests_all();
		#endif
		#endif

		#if TDS_MPQS == TDS_ODDNESS_CLASS || TDS_PRIMALITY_TEST == TDS_ODDNESS_CLASS
				{
				  clock_t new_clock;
				  new_clock= clock();
				  td_clock+= new_clock-last_clock;
				  last_clock= new_clock;
				}
		#endif

			  }

		#if TDS_MPQS == TDS_SPECIAL_Q

			  output_all_tdsurvivors();

		#else
		#if TDS_PRIMALITY_TEST == TDS_SPECIAL_Q
			  primality_tests_all();
		#endif
		#endif
		#if TDS_MPQS == TDS_SPECIAL_Q || TDS_PRIMALITY_TEST == TDS_SPECIAL_Q
			  {
				clock_t new_clock;
				new_clock= clock();
				td_clock+= new_clock-last_clock;
				last_clock= new_clock;
			  }
		#endif

			  
			  
			}
	
			// done with this root
			//fprintf(g_ofile,"# Done %u %u (%d,%d) (%d,%d)\n",
			//special_q,r[root_no],a0,b0,a1,b1);
		  }

		  if (root_no < nr) {
			  break;
		  }
		  tNow = sTime();
		  if (tNow > lastReport + 5.0) {
			  lastReport = tNow;
			  if (verbose) {
				  int eta = (int)(((double)last_spq - special_q) * (tNow - tStart) / ((double)special_q - first_spq + 1) / 60);
				  fprintf(stderr, "\rtotal yield: %u, q=%u (%1.5lf sec/rel; ETA %dh%02dm)  ",
					  (unsigned int)yield, (unsigned int)special_q, (tNow - tStart) / yield,
					  eta / 60, eta % 60);
				  fflush(stderr);
			  }
		  }

		  // done with this special_q
	  }


    fprintf(stderr, "\rtotal yield: %u, q=%u (%1.5lf sec/rel) \n", 
	    (unsigned int)yield, (unsigned int)special_q, (sTime() - tStart)/yield);
    free(r);
  }
  
  if(sieve_count!=0) {
    if(zip_output!=0)pclose(g_ofile);
    else fclose(g_ofile);
  }
  logbook(0,"%u Special q, %u reduction iterations\n",n_spq,n_iter);
  
  if(n_spq_discard> 0)
	  logbook(0,"%u Special q discarded\n",n_spq_discard);

  // report timings
  {
	  u32_t side;
	  logbook(0, "reports: %u->%u->%u->%u->%u->%u\n",
		  n_prereports, n_reports, n_rep1, n_rep2,
		  n_tdsurvivors[first_td_side], n_tdsurvivors[1 - first_td_side]);
	  logbook(0, "Number of relations with k rational and l algebraic primes for (k,l)=:\n");

	  sieve_clock = rint((1000.0 * sieve_clock) / CLOCKS_PER_SEC);
	  sch_clock = rint((1000.0 * sch_clock) / CLOCKS_PER_SEC);
	  lsetup_clock = rint((1000.0 * lsetup_clock) / CLOCKS_PER_SEC);
	  td_clock = rint((1000.0 * td_clock) / CLOCKS_PER_SEC);
	  tdi_clock = rint((1000.0 * tdi_clock) / CLOCKS_PER_SEC);
	  Schedule_clock = rint((1000.0 * Schedule_clock) / CLOCKS_PER_SEC);
	  medsched_clock = rint((1000.0 * medsched_clock) / CLOCKS_PER_SEC);
	  mpqs_clock = rint((1000.0 * mpqs_clock) / CLOCKS_PER_SEC);

	  for (side = 0; side < 2; side++) {
		  cs_clock[side] = rint((1000.0 * cs_clock[side]) / CLOCKS_PER_SEC);
		  si_clock[side] = rint((1000.0 * si_clock[side]) / CLOCKS_PER_SEC);
		  s1_clock[side] = rint((1000.0 * s1_clock[side]) / CLOCKS_PER_SEC);
		  s2_clock[side] = rint((1000.0 * s2_clock[side]) / CLOCKS_PER_SEC);
		  s3_clock[side] = rint((1000.0 * s3_clock[side]) / CLOCKS_PER_SEC);
		  tdsi_clock[side] = rint((1000.0 * tdsi_clock[side]) / CLOCKS_PER_SEC);
		  tds1_clock[side] = rint((1000.0 * tds1_clock[side]) / CLOCKS_PER_SEC);
		  tds2_clock[side] = rint((1000.0 * tds2_clock[side]) / CLOCKS_PER_SEC);
		  tds3_clock[side] = rint((1000.0 * tds3_clock[side]) / CLOCKS_PER_SEC);
		  tds4_clock[side] = rint((1000.0 * tds4_clock[side]) / CLOCKS_PER_SEC);
	  }

	  logbook(0, "\nTotal yield: %u\n", yield);
	  if (n_mpqsfail[0] != 0 || n_mpqsfail[1] != 0 ||
		  n_mpqsvain[0] != 0 || n_mpqsvain[1] != 0) {
		  logbook(0, "%u/%u mpqs failures, %u/%u vain mpqs\n", n_mpqsfail[0],
			  n_mpqsfail[1], n_mpqsvain[0], n_mpqsvain[1]);
	  }
	  logbook(0, "milliseconds total: Sieve %d Sched %d medsched %d\n",
		  (int)sieve_clock, (int)Schedule_clock, (int)medsched_clock);
	  logbook(0, "TD %d (Init %d, MPQS %d) Sieve-Change %d, lasieve_setup %d\n",
		  (int)td_clock, (int)tdi_clock, (int)mpqs_clock, (int)sch_clock, (int)lsetup_clock);
	  for (side = 0; side < 2; side++) {
		  logbook(0, "TD side %d: init/small/medium/large/search: %d %d %d %d %d\n",
			  (int)side, (int)tdsi_clock[side], (int)tds1_clock[side],
			  (int)tds2_clock[side], (int)tds3_clock[side], (int)tds4_clock[side]);
		  logbook(0, "sieve: init/small/medium/large/search: %d %d %d %d %d\n",
			  (int)si_clock[side], (int)s1_clock[side], (int)s2_clock[side],
			  (int)s3_clock[side], (int)cs_clock[side]);
	  }
#ifdef MMX_TDBENCH
	  fprintf(stderr, "MMX-Loops: %qu\n", MMX_TdNloop);
#endif
#ifdef ZSS_STAT
	  fprintf(stderr,
		  "%u subsieves, zero: %u first sieve, %u second sieve %u first td\n",
		  nss, nzss[0], nzss[1], nzss[2]);
#endif
  }

	logTotalTime();
	if (special_q >= last_spq && all_spq_done != 0)exit(0);
	if (exitval == 0)exitval = 1;
	exit(exitval);
}

#ifndef NOSCHED
void do_scheduling(struct schedule_struct* sched, u32_t ns, u32_t ot, u32_t s)
{
	u32_t ll, n1_j, * ri;;

	n1_j = ns << (L1_BITS - i_bits);
	for (ll = 0, ri = sched->ri; ll < sched->n_pieces; ll++) {
		u32_t fbi_lb, fbi_ub, fbio;
		memcpy(sched->schedule[ll + 1], sched->schedule[ll], ns * sizeof(u16_t**));
		fbio = sched->fbi_bounds[ll];
		fbi_lb = fbio;
		fbi_ub = sched->fbi_bounds[ll + 1];

#ifdef SCHEDULING_FUNCTION_CALCULATES_RI
		if (ot == 1)
			lasieve_setup(FB[s] + fbi_lb, proots[s] + fbi_lb, fbi_ub - fbi_lb,
				a0, a1, b0, b1, LPri[s] + (fbi_lb - fbis[s]) * RI_SIZE);
#endif
		//printf("lasched %d on side %d from ij %u to %u\n", ll, s, fbi_lb, fbi_ub);
		ri = lasched(ri, current_ij[s] + fbi_lb, current_ij[s] + fbi_ub,
			n1_j, (u32_t**)(sched->schedule[ll + 1]), fbi_lb - fbio, ot, FBsize[s]);

		{
			u32_t k;
			for (k = 0; k < ns; k++)
				if (sched->schedule[ll + 1][k] >= sched->schedule[0][k] + sched->alloc) {
					if (k == 0 && sched->schedule[ll + 1][k] < sched->schedule[0][k] + sched->alloc1)
						continue;
					/* report SCHED_PATHOLOGY to http://mersenneforum.org/showthread.php?t=11430 */
					fprintf(stderr, "\rSCHED_PATHOLOGY q0=%u k=%d excess=%d                      \n",
						(unsigned int)special_q, k, sched->schedule[ll + 1][k] - (sched->schedule[0][k] + sched->alloc));
					longjmp(termination_jb, SCHED_PATHOLOGY);
				}
		}

	}
}
#endif

#ifdef GCD_SIEVE_BOUND
static void
gcd_sieve()
{
  u32_t i;
  
  for(i= 0;i<np_gcd_sieve;i++) {
    u32_t x,p;
    
    x= gcd_sieve_buffer[2*i+1];
    p= gcd_sieve_buffer[2*i];
    while(x<j_per_strip) {
      unsigned char*z,*z_ub;
      
      z= sieve_interval+(x<<i_bits);
      z_ub= z+n_i-3*p;
      z+= oddness_type == 2?(n_i/2)%p:((n_i+p-1)/2)%p;
      while(z<z_ub) {
	*z= 0;
	*(z+p)= 0;
	z+= 2*p;
	*z= 0;
	*(z+p)= 0;
	z+= 2*p;
      }
      z_ub+= 3*p;
      while(z<z_ub) {
	*z= 0;
	z+= p;
      }
      x= x+p;
    }
    gcd_sieve_buffer[2*i+1]= x-j_per_strip;
  }
}
#endif

static void
xFBtranslate(u16_t*rop,xFBptr op)
{
  u32_t x,y,am,bm,rqq;
  
  modulo32= op->pp;
  rop[3]= op->l;
  am= a1> 0?((u32_t)a1)%modulo32:modulo32-((u32_t)(-a1))%modulo32;
  if(am == modulo32)am= 0;
  bm= b1> 0?((u32_t)b1)%modulo32:modulo32-((u32_t)(-b1))%modulo32;
  if(bm == modulo32)bm= 0;
  x= modsub32(modmul32(op->qq,am),modmul32(op->r,bm));
  am= a0> 0?((u32_t)a0)%modulo32:modulo32-((u32_t)(-a0))%modulo32;
  if(am == modulo32)am= 0;
  bm= b0> 0?((u32_t)b0)%modulo32:modulo32-((u32_t)(-b0))%modulo32;
  if(bm == modulo32)bm= 0;
  y= modsub32(modmul32(op->r,bm),modmul32(op->qq,am));
  rqq= 1;
  if(y!=0) {
    while(y%(op->p) == 0) {
      y= y/(op->p);
      rqq*= op->p;
    }
  } else {
    rqq= op->pp;
  }
  modulo32= modulo32/rqq;
  rop[0]= modulo32;
  rop[1]= rqq;
  if(modulo32> 1)
    rop[2]= modmul32(modinv32(y),x);
  else rop[2]= 0;
  rop[4]= op->l;
}

static int
xFBcmp(const void*opA,const void*opB)
{
  xFBptr op1,op2;
  op1= (xFBptr)opA;
  op2= (xFBptr)opB;
  if(op1->pp<op2->pp)return-1;
  if(op1->pp == op2->pp)return 0;
  return 1;
}

static u32_t
add_primepowers2xaFB(size_t*xaFB_alloc_ptr,u32_t pp_bound,
		     u32_t s,u32_t p,u32_t r)
{
  u32_t a,b,q,qo,*rbuf,nr,*Ar,exponent,init_xFB;
  size_t rbuf_alloc;
  if(xFBs[s] == 0&&p == 0)Schlendrian("add_primepowers2xaFB on empty xaFB\n");
  
  rbuf_alloc= 0;
  Ar= xmalloc((1+poldeg[s])*sizeof(*Ar));
  
  if(p!=0) {
    init_xFB= 0;
    q= p;
    if(r == p) {
      a= 1;
      b= p;
    } else {
      a= r;
      b= 1;
    }
  } else {
    init_xFB= 1;
    q= xFB[s][xFBs[s]-1].pp;
    p= xFB[s][xFBs[s]-1].p;
    a= xFB[s][xFBs[s]-1].r;
    b= xFB[s][xFBs[s]-1].qq;
  }
  
  qo= q;
  exponent= 1;
  for(;;) {
    u32_t j,r;
    if(q> pp_bound/p)break;
    modulo32= p*q;
    for(j= 0;j<=poldeg[s];j++)
      Ar[j]= mpz_fdiv_ui(poly[s][j],modulo32);
    if(b == 1) {
      for(r= a,nr= 0;r<modulo32;r+= qo) {
	u32_t pv;
	for(j= 1,pv= Ar[poldeg[s]];j<=poldeg[s];j++) {
	  pv= modadd32(Ar[poldeg[s]-j],modmul32(pv,r));
	}
	if(pv == 0) {
	  adjust_bufsize((void**)&rbuf,&rbuf_alloc,1+nr,4,sizeof(*rbuf));
	  rbuf[nr++]= r;
	} else if(pv%q!=0)Schlendrian("xFBgen: %u not a root mod %u\n", r,q);
      }
    } else {
      for(r= (modmul32(b,modinv32(a)))%qo,nr= 0;r<modulo32;r+= qo) {
	u32_t pv;
	for(j= 1,pv= Ar[0];j<=poldeg[s];j++) {
	  pv= modadd32(Ar[j],modmul32(pv,r));
	}
	if(pv == 0) {
	  adjust_bufsize((void**)&rbuf,&rbuf_alloc,1+nr,4,sizeof(*rbuf));
	  rbuf[nr++]= r;
	} else if(pv%q!=0)Schlendrian("xFBgen: %u^{-1} not a root mod %u\n", r,q);
      }
    }
    
    if(qo*nr!=modulo32)break;
    q= modulo32;
    exponent++;
  }
  if(init_xFB!=0)
    xFB[s][xFBs[s]-1].l= 
      rint(sieve_multiplier[s]*log(q))-rint(sieve_multiplier[s]*log(qo/p));
  if(q<=pp_bound/p) {
    u32_t j;
    for(j= 0;j<nr;j++) {
      xFBptr f;
      
      adjust_bufsize((void**)&(xFB[s]),xaFB_alloc_ptr,1+xFBs[s],16,sizeof(**xFB));
      f= xFB[s]+xFBs[s];
      f->p= p;
      f->pp= q*p;
      if(b == 1) {
	f->qq= 1;
	f->r= rbuf[j];
	f->q= f->pp;
      } else {
	modulo32= (q*p)/b;
	rbuf[j]= rbuf[j]/b;
	if(rbuf[j] == 0) {
	  f->qq= f->pp;
	  f->q= 1;
	  f->r= 1;
	} else {
	  while(rbuf[j]%p == 0) {
	    rbuf[j]= rbuf[j]/p;
	    modulo32= modulo32/p;
	  }
	  f->qq= (f->pp)/modulo32;
	  f->q= modulo32;
	  f->r= modinv32(rbuf[j]);
	}
      }
      
      xFBs[s]++;
      add_primepowers2xaFB(xaFB_alloc_ptr,pp_bound,s,0,0);
    }
  }
  if(rbuf_alloc> 0)free(rbuf);
  free(Ar);
  return exponent;
}

void
trial_divide()
{
  u32_t ci;
  u32_t nc1;
  u16_t side,tdstep;
  clock_t last_tdclock,newclock;
  
#ifdef NO_TDCODE
  return;
#endif
  
  for(ci= 0,nc1= 0;ci<ncand;ci++) {
    u16_t strip_i,strip_j;
    u16_t st_i,true_j;
    u16_t s;
    double pvl;
    
    {
      u16_t jj;
      
      strip_j= cand[ci]>>i_bits;
      jj= j_offset+strip_j;
      strip_i= cand[ci]&(n_i-1);
      st_i= 2*strip_i+(oddness_type == 2?0:1);
      true_j= 2*jj+(oddness_type == 1?0:1);
    }
    
    n_reports++;
    s= first_sieve_side;
#ifdef STC_DEBUG
    fprintf(debugfile,"%hu %hu\n",st_i,true_j);
#endif
    if(gcd32(st_i<i_shift?i_shift-st_i:st_i-i_shift,true_j)!=1)continue;
    n_rep1++;
    pvl= log(fabs(rpol_eval(tpoly_f[s],poldeg[s],
			    (double)st_i-(double)i_shift,(double)true_j)));
    if(special_q_side == s)pvl-= special_q_log;
    pvl*= sieve_multiplier[s];
    if((double)fss_sv[ci]+sieve_report_multiplier[s]*FB_maxlog[s]<pvl)continue;
    n_rep2++;
    
    
    modulo32= special_q;
    if(modadd32(modmul32(st_i,spq_i),modmul32(true_j,spq_j)) == spq_x)continue;
    cand[nc1++]= cand[ci];
  }
  
#ifdef ZSS_STAT
  if(ncand == 0)
    nzss[1]++;
#endif
  last_tdclock= clock();
  tdi_clock+= last_tdclock-last_clock;
  ncand= nc1;
  qsort(cand,ncand,sizeof(*cand),tdcand_cmp);
  td_buf1[0]= td_buf[first_td_side];
  for(side= first_td_side,tdstep= 0;tdstep<2;side= 1-side,tdstep++) {

#ifdef ZSS_STAT
    if(tdstep == 1&&ncand == 0)
      nzss[2]++;
#endif

    {
      u32_t nfbp;
      
      u32_t p_bound;
      
      u16_t last_j,strip_i,strip_j;
      u16_t*smalltdsieve_auxbound;
      
      nfbp= 0;
      
      
      
      
      
#ifndef SET_TDS_PBOUND
      
	  if(ncand> 0)p_bound= (2*n_i*j_per_strip)/(5*ncand);
      else p_bound= U32_MAX;

#else
      p_bound= SET_TDS_PBOUND(n_i,j_per_strip,ncand);
#endif
      
	  {
		  unsigned char ht, allcoll;

		  bzero(sieve_interval, L1_SIZE);
		  bzero(tds_coll, UCHAR_MAX - 1);
		  for (ci = 0, ht = 1, allcoll = 0; ci < ncand; ci++) {
			  unsigned char cht;

			  cht = sieve_interval[cand[ci]];
			  if (cht == 0) {
				  cht = ht;
				  if (ht < UCHAR_MAX)ht++;
				  else {
					  ht = 1;
					  allcoll = 1;
				  }
				  tds_coll[cht - 1] = allcoll;
				  sieve_interval[cand[ci]] = cht;
			  }
			  else {
				  tds_coll[cht - 1] = 1;
			  }
			  fss_sv[ci] = cht - 1;
		  }
	  }
      
#ifdef MMX_TD
      smalltdsieve_auxbound= MMX_TdInit(side,smallsieve_aux[side],
					smallsieve_auxbound[side][0],
					&p_bound,j_offset == 0&&oddness_type == 1);
#else
	  {
		  u16_t* x, * z;

		  x = smallsieve_aux[side];
		  z = smallsieve_auxbound[side][0];
		  if (*x > p_bound)
			  smalltdsieve_auxbound = x;
		  else {

			  while (x + 4 < z) {
				  u16_t* y;

				  y = x + 4 * ((z - x) / 8);
				  if (y == smallsieve_auxbound[side][0] || *y > p_bound)z = y;
				  else x = y;
			  }
			  smalltdsieve_auxbound = z;
		  }
	  }
#endif
      
      newclock= clock();
      tdsi_clock[side]+= newclock-last_tdclock;
      last_tdclock= newclock;
      memcpy(tds_fbi_curpos,tds_fbi,UCHAR_MAX*sizeof(*tds_fbi));
      
#ifdef ASM_SCHEDTDSIEVE
      {
	u32_t x,*(y[2]);
	
	x= 0;
	y[0]= med_sched[side][0];
	y[1]= med_sched[side][n_medsched_pieces[side]];
	schedtdsieve(&x,1,y,sieve_interval,tds_fbi_curpos);
      }
#else
	  {
		  u32_t l;

		  for (l = 0; l < n_medsched_pieces[side]; l++) {
			  u16_t* x, * x_ub;

			  x_ub = med_sched[side][l + 1];

			  for (x = med_sched[side][l] + MEDSCHED_SI_OFFS; x + 6 < x_ub; x += 8) {
				  unsigned char z;
				  if ((sieve_interval[*x] | sieve_interval[*(x + 2)] |
					  sieve_interval[*(x + 4)] | sieve_interval[*(x + 6)]) == 0) {
					  continue;
				  }
				  if ((z = sieve_interval[*x]) != 0)
					  *(tds_fbi_curpos[z - 1]++) = *(x + 1 - 2 * MEDSCHED_SI_OFFS);
				  if ((z = sieve_interval[*(x + 2)]) != 0)
					  *(tds_fbi_curpos[z - 1]++) = *(x + 3 - 2 * MEDSCHED_SI_OFFS);
				  if ((z = sieve_interval[*(x + 4)]) != 0)
					  *(tds_fbi_curpos[z - 1]++) = *(x + 5 - 2 * MEDSCHED_SI_OFFS);
				  if ((z = sieve_interval[*(x + 6)]) != 0)
					  *(tds_fbi_curpos[z - 1]++) = *(x + 7 - 2 * MEDSCHED_SI_OFFS);
			  }
			  while (x < x_ub) {
				  unsigned char z;

				  if ((z = sieve_interval[*x]) != 0)
					  *(tds_fbi_curpos[z - 1]++) = *(x + 1 - 2 * MEDSCHED_SI_OFFS);
				  x += 2;
			  }
		  }
	  }
#endif
      newclock= clock();
      tds2_clock[side]+= newclock-last_tdclock;
      last_tdclock= newclock;
      
	  {
		  u32_t j;

		  for (j = 0; j < n_schedules[side]; j++) {
#ifdef ASM_SCHEDTDSIEVE
			  u32_t i, k;
			  k = schedules[side][j].current_strip++;
			  for (i = 0; i <= schedules[side][j].n_pieces; i++) {
				  schedbuf[i] = schedules[side][j].schedule[i][k];
			  }
			  schedtdsieve(schedules[side][j].fbi_bounds, schedules[side][j].n_pieces,
				  schedbuf, sieve_interval, tds_fbi_curpos);
#else
#if 1
			  u32_t k, l, fbi_offset;
			  u16_t* x, * x_ub;
			  k = schedules[side][j].current_strip++;
			  x = schedules[side][j].schedule[0][k] + SCHED_SI_OFFS;
			  x_ub = schedules[side][j].schedule[schedules[side][j].n_pieces][k];
			  l = 0;
			  fbi_offset = schedules[side][j].fbi_bounds[l];
			  while (x < x_ub) {
				  u16_t** b0, ** b1, ** b0_ub;
#ifdef ASM_SCHEDTDSIEVE2
				  b0 = tdsieve_sched2buf(&x, x_ub, sieve_interval, sched_tds_buffer,
					  sched_tds_buffer + SCHED_TDS_BUFSIZE - 4);
#else

#if 1
				  // very similar speed to ASM_SCHEDTDSIEVE2... gathers are slow.
				  // it might be faster for really big inputs that use large factor
				  // bases.
				  __m512i zero = _mm512_setzero_si512();
				  b0 = sched_tds_buffer;
				  b0_ub = b0 + SCHED_TDS_BUFSIZE;
				  for (; x + 94 < x_ub; x = x + 96) {

					  __m512i vx1 = _mm512_mask_loadu_epi16(zero, 0x55555555, x);
					  __m512i vx2 = _mm512_mask_loadu_epi16(zero, 0x55555555, x + 32);
					  __m512i vx3 = _mm512_mask_loadu_epi16(zero, 0x55555555, x + 64);

					  // better CPI with 3 independent gathers
					  __m512i vsi1 = _mm512_i32gather_epi32(vx1, sieve_interval, 1);
					  __m512i vsi2 = _mm512_i32gather_epi32(vx2, sieve_interval, 1);
					  __m512i vsi3 = _mm512_i32gather_epi32(vx3, sieve_interval, 1);

					  __mmask64 m1 = _mm512_mask_cmpneq_epi8_mask(0x1111111111111111ull, vsi1, zero);
					  __mmask64 m2 = _mm512_mask_cmpneq_epi8_mask(0x1111111111111111ull, vsi2, zero);
					  __mmask64 m3 = _mm512_mask_cmpneq_epi8_mask(0x1111111111111111ull, vsi3, zero);

					  if ((m1 | m2 | m3) == 0)
						  continue;

					  while (m1 > 0)
					  {
						  int id = _tzcnt_u64(m1) / 4;
						  *(b0++) = x + 2*id;
						  m1 = _blsr_u64(m1);
					  }
					  while (m2 > 0)
					  {
						  int id = _tzcnt_u64(m2) / 4;
						  *(b0++) = x + 32 + 2*id;
						  m2 = _blsr_u64(m2);
					  }
					  while (m3 > 0)
					  {
						  int id = _tzcnt_u64(m3) / 4;
						  *(b0++) = x + 64 + 2*id;
						  m3 = _blsr_u32(m3);
					  }

					  if (b0 + 48 > b0_ub)goto sched_tds1;
				  }
				  for (; x < x_ub; x += 2) {
					  if (sieve_interval[*x] != 0)*(b0++) = x;
				  }

#else
				  // could try a gather-load here
				  b0 = sched_tds_buffer;
				  b0_ub = b0 + SCHED_TDS_BUFSIZE;
				  for (; x + 6 < x_ub; x = x + 8) {
					  if ((sieve_interval[*x] | sieve_interval[*(x + 2)] |
						  sieve_interval[*(x + 4)] | sieve_interval[*(x + 6)]) == 0)
						  continue;
					  if (sieve_interval[x[0]] != 0)*(b0++) = x;
					  if (sieve_interval[x[2]] != 0)*(b0++) = x + 2;
					  if (sieve_interval[x[4]] != 0)*(b0++) = x + 4;
					  if (sieve_interval[x[6]] != 0)*(b0++) = x + 6;
					  if (b0 + 4 > b0_ub)goto sched_tds1;
				  }
				  for (; x < x_ub; x += 2) {
					  if (sieve_interval[*x] != 0)*(b0++) = x;
				  }

#endif

#endif
			  sched_tds1:
				  for (b1 = sched_tds_buffer; b1 < b0; b1++) {
					  u16_t* y;
					  u32_t fbi;
					  y = *(b1);
					  if (schedules[side][j].schedule[l + 1][k] <= y) {
						  do {
							  l++;
							  if (l >= schedules[side][j].n_pieces)Schlendrian("XXX\n");
						  } while (schedules[side][j].schedule[l + 1][k] <= y);
						  fbi_offset = schedules[side][j].fbi_bounds[l];
					  }
					  fbi = fbi_offset + *(y + 1 - 2 * SCHED_SI_OFFS);
					  *(tds_fbi_curpos[sieve_interval[*y] - 1]++) = fbi;

#ifdef TDS_FB_PREFETCH
					  TDS_FB_PREFETCH(FB[side] + fbi);
#endif

				  }
			  }
#else
			  u32_t l, k;

			  k = schedules[side][j].current_strip++;
			  for (l = 0; l < schedules[side][j].n_pieces; l++) {
				  u16_t* x, * x_ub;
				  u32_t fbi_offset;

				  x_ub = schedules[side][j].schedule[l + 1][k];
				  fbi_offset = schedules[side][j].fbi_bounds[l];
				  for (x = schedules[side][j].schedule[l][k] + SCHED_SI_OFFS; x + 6 < x_ub; x += 8) {
					  unsigned char z;
					  if ((sieve_interval[*x] | sieve_interval[*(x + 2)] |
						  sieve_interval[*(x + 4)] | sieve_interval[*(x + 6)]) == 0) {
						  continue;
					  }
					  if ((z = sieve_interval[*x]) != 0)
						  *(tds_fbi_curpos[z - 1]++) = fbi_offset + *(x + 1 - 2 * SCHED_SI_OFFS);
					  if ((z = sieve_interval[*(x + 2)]) != 0)
						  *(tds_fbi_curpos[z - 1]++) = fbi_offset + *(x + 3 - 2 * SCHED_SI_OFFS);
					  if ((z = sieve_interval[*(x + 4)]) != 0)
						  *(tds_fbi_curpos[z - 1]++) = fbi_offset + *(x + 5 - 2 * SCHED_SI_OFFS);
					  if ((z = sieve_interval[*(x + 6)]) != 0)
						  *(tds_fbi_curpos[z - 1]++) = fbi_offset + *(x + 7 - 2 * SCHED_SI_OFFS);
				  }
				  while (x < x_ub) {
					  unsigned char z;

					  if ((z = sieve_interval[*x]) != 0)
						  *(tds_fbi_curpos[z - 1]++) = fbi_offset + *(x + 1 - 2 * SCHED_SI_OFFS);
					  x += 2;
				  }
			  }
#endif
#endif
		  }
	  }
      newclock= clock();
      tds3_clock[side]+= newclock-last_tdclock;
      last_tdclock= newclock;
      
	  {
		  u32_t i;

		  for (i = 0; i < UCHAR_MAX && i < ncand; i++) {
			  u32_t* p;

			  for (p = tds_fbi[i]; p < tds_fbi_curpos[i]; p++)
				  *p = FB[side][*p];
		  }
	  }
      
	  {
		  u16_t* x;

//#define AVX512_TDS
		  // tested to be about the same speed.  restructuring for 32x processing
	      // instead of 8x would probably be needed for a win.

#ifdef ASM_TDSLINIE
		  x = smalltdsieve_auxbound;
		  if (x < smallsieve_auxbound[side][4]) {
			  tdslinie(x, smallsieve_auxbound[side][4], sieve_interval, tds_fbi_curpos);
			  x = smallsieve_auxbound[side][4];
		  }
#else
		  for (x = smalltdsieve_auxbound;
			  x < smallsieve_auxbound[side][4]; x = x + 4) {
			  u32_t p, r, pr;
			  unsigned char* y;

			  p = x[0];
			  pr = x[1];
			  r = x[3];
			  modulo32 = p;
			  for (y = sieve_interval; y < sieve_interval + L1_SIZE; y += n_i) {
				  unsigned char* yy, * yy_ub;

				  // basically we are checking for non-zero sieve locations
				  // in this prime's progression and if we find any, we
				  // append the primes to the end of a list for that
				  // sieve location.  
				  // So we basically re-do all of the sieve, except we are
				  // not writing/modifiying the sieve.
				  // yafu's approach in the QS is different for sufficiently large p.
				  // there, we find all of the non-zero locations, then for each
				  // location, determine which prime progressions hit it by
				  // vector-stepping many primes at once through their progressions
				  // and checking for equality to the index in question.
				  // for this loop, that strategy doesn't make much sense because
				  // the primes are too small... too many steps.  But later loops
				  // might be an alternative option to try.
				  yy_ub = y + n_i - 3 * p;
				  yy = y + r;
				  while (yy < yy_ub) {
					  unsigned char o;
					  o = (*yy) | (*(yy + p));
					  yy += 2 * p;
					  if ((o | (*yy) | (*(yy + p))) != 0) {
						  yy = yy - 2 * p;
						  if (*yy != 0)
							  *(tds_fbi_curpos[*yy - 1]++) = p;
						  if (*(yy + p) != 0)
							  *(tds_fbi_curpos[*(yy + p) - 1]++) = p;
						  yy += 2 * p;
						  if (*yy != 0)
							  *(tds_fbi_curpos[*yy - 1]++) = p;
						  if (*(yy + p) != 0)
							  *(tds_fbi_curpos[*(yy + p) - 1]++) = p;
					  }
					  yy += 2 * p;
				  }
				  yy_ub += 2 * p;
				  if (yy < yy_ub) {
					  if (((*yy) | (*(yy + p))) != 0) {
						  if (*yy != 0)
							  *(tds_fbi_curpos[*yy - 1]++) = p;
						  if (*(yy + p) != 0)
							  *(tds_fbi_curpos[*(yy + p) - 1]++) = p;
					  }
					  yy += 2 * p;
				  }
				  yy_ub += p;
				  if (yy < yy_ub) {
					  if (*yy != 0)
						  *(tds_fbi_curpos[*yy - 1]++) = p;
				  }
				  r = modadd32(r, pr);
			  }
			  x[3] = r;
		  }
#endif

#if !defined(AVX512_TDS) && defined( ASM_TDSLINIE3)
		  if (x < smallsieve_auxbound[side][3]) {
			  tdslinie3(x, smallsieve_auxbound[side][3], sieve_interval, tds_fbi_curpos);
			  x = smallsieve_auxbound[side][3];
		  }
#else

#if defined(AVX512_TDS)
		  unsigned char* y;
		  uint32_t nzlocs[16][128];		// 16 segments max, 128 locations max.
		  uint32_t numnzlocs[16];
		  uint32_t thisseg = 0;
		  uint32_t thisloc = 0;
		  uint32_t numlocs = 0;
		  uint32_t numsegs = 0;
		  if (1) {
			  // first find the non-zero locations for each segment of the interval,
			  // relative to the start of each segment
			  __m512i zero = _mm512_setzero_si512();

			  for (y = sieve_interval; y < sieve_interval + L1_SIZE; y += n_i) {

				  int idx = 0;
				  numlocs = 0;

				  for (idx = 0; idx < n_i; idx += 64)
				  {
					  __m512i vy = _mm512_load_si512(y + idx);
					  __mmask64 m64 = _mm512_cmpneq_epi8_mask(vy, zero);

					  while (m64 > 0)
					  {
						  uint32_t id = _tzcnt_u64(m64);
						  nzlocs[thisseg][numlocs++] = idx + id;
						  m64 = _blsr_u64(m64);
					  }
				  }

				  numnzlocs[thisseg] = numlocs;
				  thisseg++;
			  }
			  numsegs = thisseg;
		  }

		  __m512i vni = _mm512_set1_epi16(n_i);
		  uint16_t xm[32];
		  u16_t* xstart = x;

		  for (x = xstart; x < smallsieve_auxbound[side][3] - 32; x = x + 32) {

			  __mmask32 m;
			  __mmask32 mn_i;
			  __m512i xv = _mm512_load_si512(x);
			  __m512i vp = _mm512_slli_epi64(xv, 48);			// align these with r
			  __m512i vpr = _mm512_slli_epi64(xv, 32);		// align these with r
			  vpr = _mm512_and_epi64(vpr, _mm512_set1_epi64(0xffff000000000000ULL));
			  __m512i xv2 = _mm512_mask_add_epi16(xv, 0x88888888, xv, vp);
			  __m512i xv3 = _mm512_mask_add_epi16(xv2, 0x88888888, xv2, vp);
			  __m512i xv4 = _mm512_mask_add_epi16(xv3, 0x88888888, xv3, vp);

			  for (thisseg = 0, y = sieve_interval; y < sieve_interval + L1_SIZE; thisseg++, y += n_i) {

				  //for (thisloc = 0; thisloc < numnzlocs[thisseg]; thisloc++)
				 // {
				  thisloc = 0;
				  while (thisloc < numnzlocs[thisseg])
				  {
					  if ((thisloc + 3) < numnzlocs[thisseg])
					  {
						  uint32_t loc = nzlocs[thisseg][thisloc];

						  __m512i vloc = _mm512_set1_epi16(nzlocs[thisseg][thisloc]);
						  __m512i vloc2 = _mm512_set1_epi16(nzlocs[thisseg][thisloc + 1]);
						  __m512i vloc3 = _mm512_set1_epi16(nzlocs[thisseg][thisloc + 2]);
						  __m512i vloc4 = _mm512_set1_epi16(nzlocs[thisseg][thisloc + 3]);

						  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv) |
							  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv) |
							  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc3, xv) |
							  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc4, xv);

						  if (m > 0)
						  {
							  _mm512_store_si512(xm, xv);
							  loc = nzlocs[thisseg][thisloc];
							  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
							  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
							  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
							  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
							  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
							  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
							  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
							  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

							  loc = nzlocs[thisseg][thisloc + 1];

							  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
							  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
							  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
							  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
							  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
							  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
							  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
							  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

							  loc = nzlocs[thisseg][thisloc + 2];

							  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
							  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
							  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
							  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
							  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
							  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
							  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
							  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

							  loc = nzlocs[thisseg][thisloc + 3];

							  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
							  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
							  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
							  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
							  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
							  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
							  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
							  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
						  }

						  mn_i = _mm512_mask_cmplt_epu16_mask(0x88888888, xv2, vni);

						  if (mn_i > 0)
						  {
							  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv2) |
								  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv2) |
								  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc3, xv2) |
								  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc4, xv2);
							  if (m > 0)
							  {
								  _mm512_store_si512(xm, xv2);
								  loc = nzlocs[thisseg][thisloc];
								  if ((xm[3]) == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if ((xm[7]) == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if ((xm[11]) == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if ((xm[15]) == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if ((xm[19]) == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if ((xm[23]) == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if ((xm[27]) == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if ((xm[31]) == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
								  loc = nzlocs[thisseg][thisloc + 1];

								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

								  loc = nzlocs[thisseg][thisloc + 2];

								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

								  loc = nzlocs[thisseg][thisloc + 3];

								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
							  }

							  mn_i = _mm512_mask_cmplt_epu16_mask(0x88888888, xv3, vni);

							  if (mn_i > 0)
							  {
								  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv3) |
									  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv3) |
									  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc3, xv3) |
									  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc4, xv3);
								  if (m > 0)
								  {
									  _mm512_store_si512(xm, xv3);
									  loc = nzlocs[thisseg][thisloc];
									  if ((xm[3]) == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if ((xm[7]) == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if ((xm[11]) == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if ((xm[15]) == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if ((xm[19]) == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if ((xm[23]) == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if ((xm[27]) == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if ((xm[31]) == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
									  loc = nzlocs[thisseg][thisloc + 1];

									  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

									  loc = nzlocs[thisseg][thisloc + 2];

									  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

									  loc = nzlocs[thisseg][thisloc + 3];

									  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
								  }

								  mn_i = _mm512_mask_cmplt_epu16_mask(0x88888888, xv4, vni);

								  if (mn_i > 0)
								  {
									  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv4) |
										  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv4) |
										  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc3, xv4) |
										  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc4, xv4);
									  if (m > 0)
									  {
										  _mm512_store_si512(xm, xv4);
										  loc = nzlocs[thisseg][thisloc];
										  if ((xm[3]) == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
										  if ((xm[7]) == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
										  if ((xm[11]) == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
										  if ((xm[15]) == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
										  if ((xm[19]) == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
										  if ((xm[23]) == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
										  if ((xm[27]) == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
										  if ((xm[31]) == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
										  loc = nzlocs[thisseg][thisloc + 1];

										  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
										  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
										  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
										  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
										  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
										  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
										  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
										  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

										  loc = nzlocs[thisseg][thisloc + 2];

										  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
										  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
										  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
										  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
										  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
										  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
										  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
										  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

										  loc = nzlocs[thisseg][thisloc + 3];

										  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
										  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
										  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
										  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
										  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
										  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
										  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
										  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
									  }
								  }
							  }
						  }

						  thisloc += 4;
		  }
					  else if ((thisloc + 1) < numnzlocs[thisseg])
					  {
						  uint32_t loc = nzlocs[thisseg][thisloc];

						  __m512i vloc = _mm512_set1_epi16(nzlocs[thisseg][thisloc]);
						  __m512i vloc2 = _mm512_set1_epi16(nzlocs[thisseg][thisloc + 1]);

						  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv) |
							  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv);

						  if (m > 0)
						  {
							  _mm512_store_si512(xm, xv);
							  loc = nzlocs[thisseg][thisloc];
							  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
							  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
							  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
							  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
							  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
							  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
							  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
							  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

							  loc = nzlocs[thisseg][thisloc + 1];

							  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
							  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
							  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
							  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
							  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
							  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
							  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
							  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
						  }

						  mn_i = _mm512_mask_cmplt_epu16_mask(0x88888888, xv2, vni);

						  if (mn_i > 0)
						  {
							  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv2) |
								  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv2);
							  if (m > 0)
							  {
								  _mm512_store_si512(xm, xv2);
								  loc = nzlocs[thisseg][thisloc];
								  if ((xm[3]) == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if ((xm[7]) == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if ((xm[11]) == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if ((xm[15]) == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if ((xm[19]) == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if ((xm[23]) == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if ((xm[27]) == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if ((xm[31]) == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
								  loc = nzlocs[thisseg][thisloc + 1];

								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
							  }

							  mn_i = _mm512_mask_cmplt_epu16_mask(0x88888888, xv3, vni);

							  if (mn_i > 0)
							  {
								  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv3) |
									  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv3);
								  if (m > 0)
								  {
									  _mm512_store_si512(xm, xv3);
									  loc = nzlocs[thisseg][thisloc];
									  if ((xm[3]) == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if ((xm[7]) == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if ((xm[11]) == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if ((xm[15]) == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if ((xm[19]) == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if ((xm[23]) == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if ((xm[27]) == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if ((xm[31]) == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
									  loc = nzlocs[thisseg][thisloc + 1];

									  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
								  }

								  mn_i = _mm512_mask_cmplt_epu16_mask(0x88888888, xv4, vni);

								  if (mn_i > 0)
								  {
									  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv4) |
										  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv4);
									  if (m > 0)
									  {
										  _mm512_store_si512(xm, xv4);
										  loc = nzlocs[thisseg][thisloc];
										  if ((xm[3]) == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
										  if ((xm[7]) == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
										  if ((xm[11]) == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
										  if ((xm[15]) == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
										  if ((xm[19]) == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
										  if ((xm[23]) == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
										  if ((xm[27]) == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
										  if ((xm[31]) == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
										  loc = nzlocs[thisseg][thisloc + 1];

										  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
										  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
										  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
										  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
										  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
										  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
										  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
										  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
									  }
								  }
							  }
						  }

						  thisloc += 2;
					  }
					  else
					  {
						  uint32_t loc = nzlocs[thisseg][thisloc];

						  __m512i vloc = _mm512_set1_epi16(loc);
						  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv);

						  if (m > 0)
						  {
							  _mm512_store_si512(xm, xv);
							  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
							  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
							  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
							  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
							  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
							  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
							  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
							  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
						  }

						  mn_i = _mm512_mask_cmplt_epu16_mask(0x88888888, xv2, vni);

						  if (mn_i > 0)
						  {
							  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv2);
							  if (m > 0)
							  {
								  _mm512_store_si512(xm, xv2);
								  if ((xm[3]) == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if ((xm[7]) == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if ((xm[11]) == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if ((xm[15]) == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if ((xm[19]) == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if ((xm[23]) == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if ((xm[27]) == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if ((xm[31]) == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
							  }

							  mn_i = _mm512_mask_cmplt_epu16_mask(0x88888888, xv3, vni);

							  if (mn_i > 0)
							  {
								  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv3);
								  if (m > 0)
								  {
									  _mm512_store_si512(xm, xv3);
									  if ((xm[3]) == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if ((xm[7]) == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if ((xm[11]) == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if ((xm[15]) == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if ((xm[19]) == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if ((xm[23]) == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if ((xm[27]) == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if ((xm[31]) == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
								  }

								  mn_i = _mm512_mask_cmplt_epu16_mask(0x88888888, xv4, vni);

								  if (mn_i > 0)
								  {
									  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv4);
									  if (m > 0)
									  {
										  _mm512_store_si512(xm, xv4);
										  if ((xm[3]) == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
										  if ((xm[7]) == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
										  if ((xm[11]) == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
										  if ((xm[15]) == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
										  if ((xm[19]) == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
										  if ((xm[23]) == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
										  if ((xm[27]) == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
										  if ((xm[31]) == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
									  }
								  }
							  }
						  }

						  thisloc++;

					  }

					  }

				  // now update r
				  xv = _mm512_mask_add_epi16(xv, 0x88888888, xv, vpr);
				  m = _mm512_mask_cmpge_epu16_mask(0x88888888, xv, vp);
				  xv = _mm512_mask_sub_epi16(xv, m, xv, vp);

				  xv2 = _mm512_mask_add_epi16(xv, 0x88888888, xv, vp);
				  xv3 = _mm512_mask_add_epi16(xv2, 0x88888888, xv2, vp);
				  xv4 = _mm512_mask_add_epi16(xv3, 0x88888888, xv3, vp);

				  }
			  _mm512_store_si512(x, xv);
			  }

		  for (; x < smallsieve_auxbound[side][3]; x = x + 4) {
			  u32_t p, r, pr;

			  p = x[0];
			  pr = x[1];
			  r = x[3];
			  modulo32 = p;

			  for (thisseg = 0, y = sieve_interval; y < sieve_interval + L1_SIZE; thisseg++, y += n_i) {

				  for (thisloc = 0; thisloc < numnzlocs[thisseg]; thisloc++)
				  {
					  uint32_t loc = nzlocs[thisseg][thisloc];

					  if (r == loc) *(tds_fbi_curpos[y[r] - 1]++) = p;
					  if (((r + p) == loc) && ((r + p) < n_i)) *(tds_fbi_curpos[y[r + p] - 1]++) = p;
					  if (((r + 2 * p) == loc) && ((r + 2 * p) < n_i)) *(tds_fbi_curpos[y[r + 2 * p] - 1]++) = p;
					  if (((r + 3 * p) == loc) && ((r + 3 * p) < n_i)) *(tds_fbi_curpos[y[r + 3 * p] - 1]++) = p;
				  }
				  r = modadd32(r, pr);
			  }
			  x[3] = r;
		  }

#else
		  for (; x < smallsieve_auxbound[side][3]; x = x + 4) {
			  u32_t p, r, pr;
			  unsigned char* y;

			  p = x[0];
			  pr = x[1];
			  r = x[3];
			  modulo32 = p;
			  for (y = sieve_interval; y < sieve_interval + L1_SIZE; y += n_i) {
				  unsigned char* yy, * yy_ub;

				  yy_ub = y + n_i;
				  yy = y + r;
				  if (((*yy) | (*(yy + p)) | (*(yy + 2 * p))) != 0) {
					  if (*yy != 0)*(tds_fbi_curpos[*yy - 1]++) = p;
					  if (*(yy + p) != 0)*(tds_fbi_curpos[*(yy + p) - 1]++) = p;
					  if (*(yy + 2 * p) != 0)*(tds_fbi_curpos[*(yy + 2 * p) - 1]++) = p;
				  }
				  yy += 3 * p;
				  if (yy < yy_ub) {
					  if (*yy != 0)
						  *(tds_fbi_curpos[*yy - 1]++) = p;
				  }
				  r = modadd32(r, pr);
			  }
			  x[3] = r;
		  }
#endif

#endif

#if !defined(AVX512_TDS) && defined( ASM_TDSLINIE2)
		  if (x < smallsieve_auxbound[side][2]) {
			  tdslinie2(x, smallsieve_auxbound[side][2], sieve_interval, tds_fbi_curpos);
			  x = smallsieve_auxbound[side][2];
		  }
#else

#if defined(AVX512_TDS)
		  for (; x < smallsieve_auxbound[side][2] - 32; x = x + 32) {

			  __mmask32 m;
			  __mmask32 mn_i;
			  __m512i xv = _mm512_load_si512(x);
			  __m512i vp = _mm512_slli_epi64(xv, 48);			// align these with r
			  __m512i vpr = _mm512_slli_epi64(xv, 32);		// align these with r
			  vpr = _mm512_and_epi64(vpr, _mm512_set1_epi64(0xffff000000000000ULL));
			  __m512i xv2 = _mm512_mask_add_epi16(xv, 0x88888888, xv, vp);
			  __m512i xv3 = _mm512_mask_add_epi16(xv2, 0x88888888, xv2, vp);
			  __m512i xv4 = _mm512_mask_add_epi16(xv3, 0x88888888, xv3, vp);

			  for (thisseg = 0, y = sieve_interval; y < sieve_interval + L1_SIZE; thisseg++, y += n_i) {

				  //for (thisloc = 0; thisloc < numnzlocs[thisseg]; thisloc++)
				 // {
				  thisloc = 0;
				  while (thisloc < numnzlocs[thisseg])
				  {
					  if ((thisloc + 3) < numnzlocs[thisseg])
					  {
						  uint32_t loc = nzlocs[thisseg][thisloc];

						  __m512i vloc = _mm512_set1_epi16(nzlocs[thisseg][thisloc]);
						  __m512i vloc2 = _mm512_set1_epi16(nzlocs[thisseg][thisloc + 1]);
						  __m512i vloc3 = _mm512_set1_epi16(nzlocs[thisseg][thisloc + 2]);
						  __m512i vloc4 = _mm512_set1_epi16(nzlocs[thisseg][thisloc + 3]);

						  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv) |
							  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv) |
							  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc3, xv) |
							  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc4, xv);

						  if (m > 0)
						  {
							  _mm512_store_si512(xm, xv);
							  loc = nzlocs[thisseg][thisloc];
							  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
							  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
							  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
							  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
							  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
							  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
							  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
							  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

							  loc = nzlocs[thisseg][thisloc + 1];

							  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
							  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
							  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
							  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
							  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
							  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
							  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
							  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

							  loc = nzlocs[thisseg][thisloc + 2];

							  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
							  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
							  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
							  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
							  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
							  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
							  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
							  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

							  loc = nzlocs[thisseg][thisloc + 3];

							  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
							  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
							  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
							  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
							  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
							  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
							  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
							  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
						  }

						  mn_i = _mm512_mask_cmplt_epu16_mask(0x88888888, xv2, vni);

						  if (mn_i > 0)
						  {
							  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv2) |
								  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv2) |
								  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc3, xv2) |
								  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc4, xv2);
							  if (m > 0)
							  {
								  _mm512_store_si512(xm, xv2);
								  loc = nzlocs[thisseg][thisloc];
								  if ((xm[3]) == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if ((xm[7]) == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if ((xm[11]) == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if ((xm[15]) == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if ((xm[19]) == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if ((xm[23]) == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if ((xm[27]) == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if ((xm[31]) == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
								  loc = nzlocs[thisseg][thisloc + 1];

								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

								  loc = nzlocs[thisseg][thisloc + 2];

								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

								  loc = nzlocs[thisseg][thisloc + 3];

								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
							  }

							  mn_i = _mm512_mask_cmplt_epu16_mask(0x88888888, xv3, vni);

							  if (mn_i > 0)
							  {
								  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv3) |
									  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv3) |
									  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc3, xv3) |
									  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc4, xv3);
								  if (m > 0)
								  {
									  _mm512_store_si512(xm, xv3);
									  loc = nzlocs[thisseg][thisloc];
									  if ((xm[3]) == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if ((xm[7]) == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if ((xm[11]) == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if ((xm[15]) == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if ((xm[19]) == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if ((xm[23]) == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if ((xm[27]) == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if ((xm[31]) == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
									  loc = nzlocs[thisseg][thisloc + 1];

									  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

									  loc = nzlocs[thisseg][thisloc + 2];

									  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

									  loc = nzlocs[thisseg][thisloc + 3];

									  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
								  }

								  
							  }
						  }

						  thisloc += 4;
		  }
					  else if ((thisloc + 1) < numnzlocs[thisseg])
					  {
						  uint32_t loc = nzlocs[thisseg][thisloc];

						  __m512i vloc = _mm512_set1_epi16(nzlocs[thisseg][thisloc]);
						  __m512i vloc2 = _mm512_set1_epi16(nzlocs[thisseg][thisloc + 1]);

						  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv) |
							  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv);

						  if (m > 0)
						  {
							  _mm512_store_si512(xm, xv);
							  loc = nzlocs[thisseg][thisloc];
							  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
							  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
							  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
							  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
							  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
							  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
							  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
							  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

							  loc = nzlocs[thisseg][thisloc + 1];

							  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
							  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
							  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
							  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
							  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
							  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
							  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
							  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
						  }

						  mn_i = _mm512_mask_cmplt_epu16_mask(0x88888888, xv2, vni);

						  if (mn_i > 0)
						  {
							  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv2) |
								  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv2);
							  if (m > 0)
							  {
								  _mm512_store_si512(xm, xv2);
								  loc = nzlocs[thisseg][thisloc];
								  if ((xm[3]) == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if ((xm[7]) == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if ((xm[11]) == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if ((xm[15]) == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if ((xm[19]) == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if ((xm[23]) == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if ((xm[27]) == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if ((xm[31]) == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
								  loc = nzlocs[thisseg][thisloc + 1];

								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
							  }

							  mn_i = _mm512_mask_cmplt_epu16_mask(0x88888888, xv3, vni);

							  if (mn_i > 0)
							  {
								  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv3) |
									  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv3);
								  if (m > 0)
								  {
									  _mm512_store_si512(xm, xv3);
									  loc = nzlocs[thisseg][thisloc];
									  if ((xm[3]) == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if ((xm[7]) == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if ((xm[11]) == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if ((xm[15]) == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if ((xm[19]) == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if ((xm[23]) == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if ((xm[27]) == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if ((xm[31]) == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
									  loc = nzlocs[thisseg][thisloc + 1];

									  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
								  }
							  }
						  }

						  thisloc += 2;
					  }
					  else
					  {
						  uint32_t loc = nzlocs[thisseg][thisloc];

						  __m512i vloc = _mm512_set1_epi16(loc);
						  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv);

						  if (m > 0)
						  {
							  _mm512_store_si512(xm, xv);
							  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
							  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
							  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
							  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
							  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
							  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
							  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
							  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
						  }

						  mn_i = _mm512_mask_cmplt_epu16_mask(0x88888888, xv2, vni);

						  if (mn_i > 0)
						  {
							  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv2);
							  if (m > 0)
							  {
								  _mm512_store_si512(xm, xv2);
								  if ((xm[3]) == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if ((xm[7]) == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if ((xm[11]) == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if ((xm[15]) == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if ((xm[19]) == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if ((xm[23]) == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if ((xm[27]) == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if ((xm[31]) == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
							  }

							  mn_i = _mm512_mask_cmplt_epu16_mask(0x88888888, xv3, vni);

							  if (mn_i > 0)
							  {
								  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv3);
								  if (m > 0)
								  {
									  _mm512_store_si512(xm, xv3);
									  if ((xm[3]) == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if ((xm[7]) == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if ((xm[11]) == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if ((xm[15]) == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if ((xm[19]) == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if ((xm[23]) == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if ((xm[27]) == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if ((xm[31]) == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
								  }

							  }
						  }

						  thisloc++;

					  }

					  }

				  // now update r
				  xv = _mm512_mask_add_epi16(xv, 0x88888888, xv, vpr);
				  m = _mm512_mask_cmpge_epu16_mask(0x88888888, xv, vp);
				  xv = _mm512_mask_sub_epi16(xv, m, xv, vp);

				  xv2 = _mm512_mask_add_epi16(xv, 0x88888888, xv, vp);
				  xv3 = _mm512_mask_add_epi16(xv2, 0x88888888, xv2, vp);

				  }
			  _mm512_store_si512(x, xv);
			  }

		  for (; x < smallsieve_auxbound[side][2]; x = x + 4) {
			  u32_t p, r, pr;

			  p = x[0];
			  pr = x[1];
			  r = x[3];
			  modulo32 = p;

			  for (thisseg = 0, y = sieve_interval; y < sieve_interval + L1_SIZE; thisseg++, y += n_i) {

				  for (thisloc = 0; thisloc < numnzlocs[thisseg]; thisloc++)
				  {
					  uint32_t loc = nzlocs[thisseg][thisloc];

					  if (r == loc) *(tds_fbi_curpos[y[r] - 1]++) = p;
					  if (((r + p) == loc) && ((r + p) < n_i)) *(tds_fbi_curpos[y[r + p] - 1]++) = p;
					  if (((r + 2 * p) == loc) && ((r + 2 * p) < n_i)) *(tds_fbi_curpos[y[r + 2 * p] - 1]++) = p;
				  }
				  r = modadd32(r, pr);
			  }
			  x[3] = r;
		  }


#else
		  for (; x < smallsieve_auxbound[side][2]; x = x + 4) {
			  u32_t p, r, pr;
			  unsigned char* y;

			  p = x[0];
			  pr = x[1];
			  r = x[3];
			  modulo32 = p;
			  for (y = sieve_interval; y < sieve_interval + L1_SIZE; y += n_i) {
				  unsigned char* yy, * yy_ub;

				  yy_ub = y + n_i;
				  yy = y + r;
				  if (((*yy) | (*(yy + p))) != 0) {
					  if (*yy != 0)*(tds_fbi_curpos[*yy - 1]++) = p;
					  if (*(yy + p) != 0)*(tds_fbi_curpos[*(yy + p) - 1]++) = p;
				  }
				  yy += 2 * p;
				  if (yy < yy_ub) {
					  if (*yy != 0)
						  *(tds_fbi_curpos[*yy - 1]++) = p;
				  }
				  r = modadd32(r, pr);
			  }
			  x[3] = r;
		  }
#endif

#endif

		  //#define AVX512_TDS1

#if !defined(AVX512_TDS) && defined( ASM_TDSLINIE1)
		  if (x < smallsieve_auxbound[side][1]) {
			  tdslinie1(x, smallsieve_auxbound[side][1], sieve_interval, tds_fbi_curpos);
			  x = smallsieve_auxbound[side][1];
		  }
#else

#if 0 //def AVX512_TDS1
		  // basically we are checking for non-zero sieve locations
			// in this prime's progression and if we find any, we
			// append the primes to the end of a list for that
			// sieve location.  
			// So we basically re-do all of the sieve, except we are
			// not writing/modifiying the sieve.
			// yafu's approach in the QS is different for sufficiently large p.
			// there, we find all of the non-zero locations, then for each
			// location, determine which prime progressions hit it by
			// vector-stepping many primes at once through their progressions
			// and checking for equality to the index in question.
			// for the small-prime loop, that strategy doesn't make much sense because
			// the primes are too small... too many steps.  But later loops
			// might be an alternative option to try.
		  {
			  unsigned char* y;

			  uint32_t nzlocs[16][128];		// 16 segments max, 128 locations max.
			  uint32_t numnzlocs[16];
			  uint32_t thisseg = 0;
			  uint32_t thisloc = 0;
			  uint32_t numlocs = 0;
			  uint32_t numsegs = 0;
			  if (1) {
				  // first find the non-zero locations for each segment of the interval,
				  // relative to the start of each segment
				  __m512i zero = _mm512_setzero_si512();

				  for (y = sieve_interval; y < sieve_interval + L1_SIZE; y += n_i) {

					  int idx = 0;
					  numlocs = 0;

					  for (idx = 0; idx < n_i; idx += 64)
					  {
						  __m512i vy = _mm512_load_si512(y + idx);
						  __mmask64 m64 = _mm512_cmpneq_epi8_mask(vy, zero);

						  while (m64 > 0)
						  {
							  uint32_t id = _tzcnt_u64(m64);
							  nzlocs[thisseg][numlocs++] = idx + id;
							  m64 = _blsr_u64(m64);
						  }
					  }

					  numnzlocs[thisseg] = numlocs;
					  thisseg++;
				  }
				  numsegs = thisseg;
			  }

			  __m512i vni = _mm512_set1_epi16(n_i);
			  uint16_t xm[32];

			  for ( ; x < smallsieve_auxbound[side][1] - 32; x = x + 32) {

				  __mmask32 m;
				  __mmask32 mn_i;
				  __m512i xv = _mm512_load_si512(x);
				  __m512i vp = _mm512_slli_epi64(xv, 48);			// align these with r
				  __m512i vpr = _mm512_slli_epi64(xv, 32);		// align these with r
				  vpr = _mm512_and_epi64(vpr, _mm512_set1_epi64(0xffff000000000000ULL));
				  __m512i xv2 = _mm512_mask_add_epi16(xv, 0x88888888, xv, vp);

				  for (thisseg = 0, y = sieve_interval; y < sieve_interval + L1_SIZE; thisseg++, y += n_i) {

					  thisloc = 0;
					  while (thisloc < numnzlocs[thisseg])
					  {
#if 0
						  if ((thisloc + 3) < numnzlocs[thisseg])
						  {
							  uint32_t loc = nzlocs[thisseg][thisloc];

							  __m512i vloc = _mm512_set1_epi16(nzlocs[thisseg][thisloc]);
							  __m512i vloc2 = _mm512_set1_epi16(nzlocs[thisseg][thisloc + 1]);
							  __m512i vloc3 = _mm512_set1_epi16(nzlocs[thisseg][thisloc + 2]);
							  __m512i vloc4 = _mm512_set1_epi16(nzlocs[thisseg][thisloc + 3]);

							  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv) |
								  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv) |
								  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc3, xv) |
								  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc4, xv);

							  if (m > 0)
							  {
								  _mm512_store_si512(xm, xv);
								  loc = nzlocs[thisseg][thisloc];
								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

								  loc = nzlocs[thisseg][thisloc + 1];

								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

								  loc = nzlocs[thisseg][thisloc + 2];

								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

								  loc = nzlocs[thisseg][thisloc + 3];

								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
							  }

							  mn_i = _mm512_mask_cmplt_epu16_mask(0x88888888, xv2, vni);

							  if (mn_i > 0)
							  {
								  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv2) |
									  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv2) |
									  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc3, xv2) |
									  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc4, xv2);
								  if (m > 0)
								  {
									  _mm512_store_si512(xm, xv2);
									  loc = nzlocs[thisseg][thisloc];
									  if ((xm[3]) == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if ((xm[7]) == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if ((xm[11]) == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if ((xm[15]) == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if ((xm[19]) == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if ((xm[23]) == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if ((xm[27]) == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if ((xm[31]) == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
									  loc = nzlocs[thisseg][thisloc + 1];

									  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

									  loc = nzlocs[thisseg][thisloc + 2];

									  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

									  loc = nzlocs[thisseg][thisloc + 3];

									  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
								  }

							  }

							  thisloc += 4;
						  }
						  else if ((thisloc + 1) < numnzlocs[thisseg])
						  {
							  uint32_t loc = nzlocs[thisseg][thisloc];

							  __m512i vloc = _mm512_set1_epi16(nzlocs[thisseg][thisloc]);
							  __m512i vloc2 = _mm512_set1_epi16(nzlocs[thisseg][thisloc+1]);

							  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv) |
								  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv);

							  if (m > 0)
							  {
								  _mm512_store_si512(xm, xv);
								  loc = nzlocs[thisseg][thisloc];
								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

								  loc = nzlocs[thisseg][thisloc + 1];

								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
							  }

							  mn_i = _mm512_mask_cmplt_epu16_mask(0x88888888, xv2, vni);

							  if (mn_i > 0)
							  {
								  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv2) |
									  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv2);
								  if (m > 0)
								  {
									  _mm512_store_si512(xm, xv2);
									  loc = nzlocs[thisseg][thisloc];
									  if ((xm[3]) == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if ((xm[7]) == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if ((xm[11]) == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if ((xm[15]) == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if ((xm[19]) == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if ((xm[23]) == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if ((xm[27]) == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if ((xm[31]) == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
									  loc = nzlocs[thisseg][thisloc + 1];

									  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
								  }

							  }

							  thisloc += 2;
						  }
						  else
#endif
						  {
							  uint32_t loc = nzlocs[thisseg][thisloc];

							  __m512i vloc = _mm512_set1_epi16(loc);
							  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv);

							  if (m > 0)
							  {
								  _mm512_store_si512(xm, xv);
								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
							  }

							  mn_i = _mm512_mask_cmplt_epu16_mask(0x88888888, xv2, vni);

							  if (mn_i > 0)
							  {
								  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv2);
								  if (m > 0)
								  {
									  _mm512_store_si512(xm, xv2);
									  if ((xm[3]) == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
									  if ((xm[7]) == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
									  if ((xm[11]) == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
									  if ((xm[15]) == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
									  if ((xm[19]) == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
									  if ((xm[23]) == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
									  if ((xm[27]) == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
									  if ((xm[31]) == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
								  }

							  }

							  thisloc++;

						  }

					  }

					  // now update r
					  xv = _mm512_mask_add_epi16(xv, 0x88888888, xv, vpr);
					  m = _mm512_mask_cmpge_epu16_mask(0x88888888, xv, vp);
					  xv = _mm512_mask_sub_epi16(xv, m, xv, vp);

					  xv2 = _mm512_mask_add_epi16(xv, 0x88888888, xv, vp);
				  }
				  _mm512_store_si512(x, xv);
			  }

			  for ( ; x < smallsieve_auxbound[side][1]; x = x + 4) {
			  	  u32_t p, r, pr;
			  	  
			  	  p = x[0];
			  	  pr = x[1];
			  	  r = x[3];
			  	  modulo32 = p;

				  for (thisseg = 0, y = sieve_interval; y < sieve_interval + L1_SIZE; thisseg++, y += n_i) {

			  		  for (thisloc = 0; thisloc < numnzlocs[thisseg]; thisloc++)
			  		  {
			  	  		uint32_t loc = nzlocs[thisseg][thisloc];
			  	  
			  	  		if (r == loc) *(tds_fbi_curpos[y[r] - 1]++) = p;
			  	  		if (((r + p) == loc) && ((r + p) < n_i)) *(tds_fbi_curpos[y[r + p] - 1]++) = p;
			  		  }
			  		  r = modadd32(r, pr);
				  }
				  x[3] = r;
			  }
		  }


#else
		  for (; x < smallsieve_auxbound[side][1]; x = x + 4) {
			  u32_t p, r, pr;
			  unsigned char* y;

			  p = x[0];
			  pr = x[1];
			  r = x[3];
			  modulo32 = p;
			  for (y = sieve_interval; y < sieve_interval + L1_SIZE; y += n_i) {
				  unsigned char* yy, * yy_ub;

				  yy_ub = y + n_i;
				  yy = y + r;
				  if (*yy != 0)*(tds_fbi_curpos[*yy - 1]++) = p;
				  yy += p;
				  if (yy < yy_ub) {
					  if (*yy != 0)
						  *(tds_fbi_curpos[*yy - 1]++) = p;
				  }
				  r = modadd32(r, pr);
			  }
			  x[3] = r;
		  }
#endif

#endif

#if !defined(AVX512_TDS) && defined( ASM_TDSLINIE0)
		  if (x < smallsieve_auxbound[side][0]) {
			  tdslinie0(x, smallsieve_auxbound[side][0], sieve_interval, tds_fbi_curpos);
			  x = smallsieve_auxbound[side][0];
		  }

#else


#ifdef AVX512_TDS
		  // basically we are checking for non-zero sieve locations
			// in this prime's progression and if we find any, we
			// append the primes to the end of a list for that
			// sieve location.  
			// So we basically re-do all of the sieve, except we are
			// not writing/modifiying the sieve.
			// yafu's approach in the QS is different for sufficiently large p.
			// there, we find all of the non-zero locations, then for each
			// location, determine which prime progressions hit it by
			// vector-stepping many primes at once through their progressions
			// and checking for equality to the index in question.
			// for the small-prime loop, that strategy doesn't make much sense because
			// the primes are too small... too many steps.  But later loops
			// might be an alternative option to try.
		  {
			  unsigned char* y;

			  for (; x < smallsieve_auxbound[side][0] - 32; x = x + 32) {

				  __mmask32 m;
				  __mmask32 mn_i;
				  __m512i xv = _mm512_load_si512(x);
				  __m512i vp = _mm512_slli_epi64(xv, 48);			// align these with r
				  __m512i vpr = _mm512_slli_epi64(xv, 32);		// align these with r
				  vpr = _mm512_and_epi64(vpr, _mm512_set1_epi64(0xffff000000000000ULL));

				  for (thisseg = 0, y = sieve_interval; y < sieve_interval + L1_SIZE; thisseg++, y += n_i) {

					  //for (thisloc = 0; thisloc < numnzlocs[thisseg]; thisloc++)
					 // {
					  thisloc = 0;
					  while (thisloc < numnzlocs[thisseg])
					  {
						  if ((thisloc + 3) < numnzlocs[thisseg])
						  {
							  uint32_t loc = nzlocs[thisseg][thisloc];

							  __m512i vloc = _mm512_set1_epi16(nzlocs[thisseg][thisloc]);
							  __m512i vloc2 = _mm512_set1_epi16(nzlocs[thisseg][thisloc + 1]);
							  __m512i vloc3 = _mm512_set1_epi16(nzlocs[thisseg][thisloc + 2]);
							  __m512i vloc4 = _mm512_set1_epi16(nzlocs[thisseg][thisloc + 3]);

							  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv) |
								  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv) |
								  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc3, xv) |
								  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc4, xv);

							  if (m > 0)
							  {
								  _mm512_store_si512(xm, xv);
								  loc = nzlocs[thisseg][thisloc];
								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

								  loc = nzlocs[thisseg][thisloc + 1];

								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

								  loc = nzlocs[thisseg][thisloc + 2];

								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

								  loc = nzlocs[thisseg][thisloc + 3];

								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
							  }

							  thisloc += 4;
						  }
						  else if ((thisloc + 1) < numnzlocs[thisseg])
						  {
							  uint32_t loc = nzlocs[thisseg][thisloc];

							  __m512i vloc = _mm512_set1_epi16(nzlocs[thisseg][thisloc]);
							  __m512i vloc2 = _mm512_set1_epi16(nzlocs[thisseg][thisloc + 1]);

							  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv) |
								  _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc2, xv);

							  if (m > 0)
							  {
								  _mm512_store_si512(xm, xv);
								  loc = nzlocs[thisseg][thisloc];
								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];

								  loc = nzlocs[thisseg][thisloc + 1];

								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
							  }

							  thisloc += 2;
						  }
						  else
						  {
							  uint32_t loc = nzlocs[thisseg][thisloc];

							  __m512i vloc = _mm512_set1_epi16(loc);
							  m = _mm512_mask_cmpeq_epu16_mask(0x88888888, vloc, xv);

							  if (m > 0)
							  {
								  _mm512_store_si512(xm, xv);
								  if (xm[3] == loc) *(tds_fbi_curpos[y[xm[3]] - 1]++) = x[0];
								  if (xm[7] == loc) *(tds_fbi_curpos[y[xm[7]] - 1]++) = x[4];
								  if (xm[11] == loc) *(tds_fbi_curpos[y[xm[11]] - 1]++) = x[8];
								  if (xm[15] == loc) *(tds_fbi_curpos[y[xm[15]] - 1]++) = x[12];
								  if (xm[19] == loc) *(tds_fbi_curpos[y[xm[19]] - 1]++) = x[16];
								  if (xm[23] == loc) *(tds_fbi_curpos[y[xm[23]] - 1]++) = x[20];
								  if (xm[27] == loc) *(tds_fbi_curpos[y[xm[27]] - 1]++) = x[24];
								  if (xm[31] == loc) *(tds_fbi_curpos[y[xm[31]] - 1]++) = x[28];
							  }

							  thisloc++;

						  }

					  }

					  // now update r
					  xv = _mm512_mask_add_epi16(xv, 0x88888888, xv, vpr);
					  m = _mm512_mask_cmpge_epu16_mask(0x88888888, xv, vp);
					  xv = _mm512_mask_sub_epi16(xv, m, xv, vp);
				  }
				  _mm512_store_si512(x, xv);
			  }

			  for (; x < smallsieve_auxbound[side][0]; x = x + 4) {
				  u32_t p, r, pr;

				  p = x[0];
				  pr = x[1];
				  r = x[3];
				  modulo32 = p;

				  for (thisseg = 0, y = sieve_interval; y < sieve_interval + L1_SIZE; thisseg++, y += n_i) {

					  for (thisloc = 0; thisloc < numnzlocs[thisseg]; thisloc++)
					  {
						  uint32_t loc = nzlocs[thisseg][thisloc];

						  if (((r) == loc) && ((r) < n_i)) *(tds_fbi_curpos[y[r] - 1]++) = p;
					  }
					  r = modadd32(r, pr);
				  }
				  x[3] = r;
			  }
		  }


#else

#if 0

			// crucial!  we are just coming out of the legacy tds assembly.
			_mm256_zeroupper();

			__m512i vni = _mm512_set1_epi16(n_i);
			for (; x < smallsieve_auxbound[side][0] - 32; x = x + 32) {
				unsigned char* y;

				__m512i xv = _mm512_load_si512(x);
				__m512i vp = _mm512_slli_epi64(xv, 48);			// align these with r
				__m512i vpr = _mm512_slli_epi64(xv, 32);		// align these with r
				vpr = _mm512_and_epi64(vpr, _mm512_set1_epi64(0xffff000000000000ULL));

				for (y = sieve_interval; y < sieve_interval + L1_SIZE; y += n_i) {

					__mmask8 m = _mm512_mask_cmplt_epu16_mask(0x88888888, xv, vni);

					while (m > 0)
					{
						int id = _tzcnt_u32(m) / 4;
						unsigned char* yy = y + x[4 * id + 3];
						u32_t p = x[4 * id];
						if (*yy != 0){
							*(tds_fbi_curpos[*yy - 1]++) = p;
							printf("adding prime %u to tds_fbi_curpos[%u-1]\n", p, *yy);
						}
						m = _blsr_u32(m);
					}

					//yy_ub = y + n_i;
					//yy = y + r;
					//if (yy < yy_ub) {
					//	if (*yy != 0)
					//		*(tds_fbi_curpos[*yy - 1]++) = p;
					//}
					//r = modadd32(r, pr);

					xv = _mm512_mask_add_epi16(xv, 0x88888888, xv, vpr);
					m = _mm512_mask_cmpge_epu16_mask(0x88888888, xv, vp);
					xv = _mm512_mask_sub_epi16(xv, 0x88888888, xv, vp);

				}
				//x[3] = r;
				_mm512_store_si512(x, xv);
			}

			for (; x < smallsieve_auxbound[side][0]; x = x + 4) {
				u32_t p, r, pr;
				unsigned char* y;

				p = x[0];
				pr = x[1];
				r = x[3];
				modulo32 = p;
				for (y = sieve_interval; y < sieve_interval + L1_SIZE; y += n_i) {
					unsigned char* yy, * yy_ub;

					yy_ub = y + n_i;
					yy = y + r;
					if (yy < yy_ub) {
						if (*yy != 0)
							*(tds_fbi_curpos[*yy - 1]++) = p;
					}
					r = modadd32(r, pr);
				}
				x[3] = r;
			}

#else

		  for (; x < smallsieve_auxbound[side][0]; x = x + 4) {
			  u32_t p, r, pr;
			  unsigned char* y;

			  p = x[0];
			  pr = x[1];
			  r = x[3];
			  modulo32 = p;
			  for (y = sieve_interval; y < sieve_interval + L1_SIZE; y += n_i) {
				  unsigned char* yy, * yy_ub;

				  yy_ub = y + n_i;
				  yy = y + r;
				  if (yy < yy_ub) {
					  if (*yy != 0)
						  *(tds_fbi_curpos[*yy - 1]++) = p;
				  }
				  r = modadd32(r, pr);
			  }
			  x[3] = r;
		  }
		  
#endif

#endif

#endif
		  newclock = clock();
		  tds1_clock[side] += newclock - last_tdclock;
		  last_tdclock = newclock;
	  }
      
      last_j= 0;
	  for (ci = 0, nc1 = 0; ci < ncand; ci++) {
		  u32_t* fbp_buf;
		  u32_t* fbp_ptr;
		  u16_t st_i, true_j;
		  i32_t true_i;

		  u32_t coll;
		  {
			  u16_t jj;

			  strip_j = cand[ci] >> i_bits;
			  jj = j_offset + strip_j;
			  strip_i = cand[ci] & (n_i - 1);
			  st_i = 2 * strip_i + (oddness_type == 2 ? 0 : 1);
			  true_j = 2 * jj + (oddness_type == 1 ? 0 : 1);
		  }

		  if (strip_j != last_j) {
			  u16_t j_step;
			  if (strip_j <= last_j)
				  Schlendrian("TD: Not sorted\n");
			  j_step = strip_j - last_j;
			  last_j = strip_j;


//#define AVX512_MMX_TD

#ifdef MMX_TD
			  MMX_TdUpdate(side, j_step);
#else
			  {
				  u32_t i;
				  u16_t* x, * y;

				  y = smalltdsieve_aux[side][j_step - 1];
				  for (i = 0, x = smallsieve_aux[side]; x < smallsieve_auxbound[side][0]; i++, x += 4) {
					  modulo32 = x[0];
					  if (modulo32 > p_bound)break;
					  x[3] = modadd32((u32_t)x[3], (u32_t)y[i]);
				  }
			  }
#endif
			  {
				  u16_t* x;
				  for (x = smallpsieve_aux[side]; x < smallpsieve_aux_ub[side]; x += 3) {
					  modulo32 = x[0];
					  x[2] = modsub32(x[2], (j_step) % modulo32);
				  }
			  }

		  }
		  true_i = (i32_t)st_i - (i32_t)i_shift;
		  mpz_set_si(aux1, true_i);
		  mpz_mul_si(aux1, aux1, a0);
		  mpz_set_si(aux2, a1);
		  mpz_mul_ui(aux2, aux2, (u32_t)true_j);
		  mpz_add(sr_a, aux1, aux2);

		  mpz_set_si(aux1, true_i);
		  mpz_mul_si(aux1, aux1, b0);
		  mpz_set_si(aux2, b1);
		  mpz_mul_ui(aux2, aux2, (u32_t)true_j);
		  mpz_add(sr_b, aux1, aux2);

		  if (mpz_sgn(sr_b) < 0) {
			  mpz_neg(sr_b, sr_b);
			  mpz_neg(sr_a, sr_a);
		  }

		  {
			  u32_t i;
			  i = 1;
			  mpz_set(aux2, sr_a);
			  mpz_set(aux1, poly[side][0]);
			  for (;;) {

				  mpz_mul(aux1, aux1, sr_b);
				  mpz_mul(aux3, aux2, poly[side][i]);
				  mpz_add(aux1, aux1, aux3);
				  if (++i > poldeg[side])break;

				  mpz_mul(aux2, aux2, sr_a);
			  }
		  }

		  if (td_buf_alloc[side] < nfbp + mpz_sizeinbase(aux1, 2)) {

			  td_buf_alloc[side] += 1024;
			  while (td_buf_alloc[side] < nfbp + mpz_sizeinbase(aux1, 2)) {
				  td_buf_alloc[side] += 1024;
			  }
			  td_buf[side] = xrealloc(td_buf[side], td_buf_alloc[side] * sizeof(**td_buf));
			  if (side == first_td_side) {
				  u32_t i, * oldptr;

				  oldptr = td_buf1[0];
				  for (i = 0; i <= nc1; i++)
					  td_buf1[i] = td_buf[side] + (td_buf1[i] - oldptr);
			  }
		  }
		  if (side == first_td_side)fbp_buf = td_buf1[nc1];
		  else fbp_buf = td_buf[side];
		  fbp_ptr = fbp_buf;

		  if (side == special_q_side) {
			  *(fbp_ptr++) = special_q;
		  }

		  {
			  int np, x;

			  x = fss_sv[ci];
			  np = tds_fbi_curpos[x] - tds_fbi[x];
			  memcpy(fbp_ptr, tds_fbi[x], np * sizeof(*fbp_ptr));
			  fbp_ptr += np;
		  }

		  {
			  u16_t* x;

#ifndef MMX_TD

#ifdef PREINVERT
			  {
				  u32_t* p_inv;

				  p_inv = smalltd_pi[side];

#if AVX512_MMX_TD
				  __m512i vstrip_i = _mm512_set1_epi64(strip_i);
				  __m512i zero = _mm512_setzero_si512();
				  for (x = smallsieve_aux[side];
					  (x < smallsieve_auxbound[side][0]-32) && x[28] <= p_bound; x += 32, p_inv+=8) {

					  __m512i xv = _mm512_load_si512(x);
					  __m512i vp = _mm512_and_epi64(xv, _mm512_set1_epi64(0xffff));  // isolate primes
					  __m256i iv256 = _mm256_loadu_epi32(p_inv);
					  __m512i vr = _mm512_srli_epi64(xv, 48);	// align roots
					  __m512i iv = _mm512_cvtepu32_epi64(iv256);

					  __mmask8 m = _mm512_cmpgt_epu64_mask(vr, vstrip_i);
					  vr = _mm512_sub_epi64(vstrip_i, vr);
					  vr = _mm512_mask_add_epi64(vr, m, vr, vp);

					  iv = _mm512_mul_epu32(iv, vr);				  
					  m = _mm512_cmpeq_epu64_mask(_mm512_and_epi64(iv, _mm512_set1_epi64(0xffff0000)), zero);
					  uint32_t m32 = m & 0xff;
					  while (m32 > 0)
					  {
						  uint32_t id = _tzcnt_u32(m32);
						  *(fbp_ptr++) = x[4*id];
						  m32 = _blsr_u32(m32);
					  }
				  }

				  for (;
					  x < smallsieve_auxbound[side][0] && *x <= p_bound; x += 4, p_inv++) {
					  modulo32 = *x;

					  if (((modsub32((u32_t)strip_i, (u32_t)(x[3])) * (*p_inv)) & 0xffff0000) == 0) {
						  *(fbp_ptr++) = *x;
					  }
				  }
#else
				  for (x = smallsieve_aux[side];
					  x < smallsieve_auxbound[side][0] && *x <= p_bound; x += 4, p_inv++) {
					  modulo32 = *x;

					  if (((modsub32((u32_t)strip_i, (u32_t)(x[3])) * (*p_inv)) & 0xffff0000) == 0) {
						  *(fbp_ptr++) = *x;
					  }
				  }
#endif
			  }

#else

//#define AVX512_MMX_TD

#if 0 //def AVX512_MMX_TD

			  // needs a faster vector divide
			  __m512i vstrip_i = _mm512_set1_epi64(strip_i);
			  __m512i zero = _mm512_setzero_si512();
			  for (x = smallsieve_aux[side];
				  (x < (smallsieve_auxbound[side][0] - 32)) && (x[28] <= p_bound); x += 32) {

				  __m512i xv = _mm512_load_si512(x);
				  __m512i vp = _mm512_and_epi64(xv, _mm512_set1_epi64(0xffff));  // isolate primes
				  __m512i vr = _mm512_srli_epi64(xv, 48);	// align roots

				  __m512i vt = _mm512_rem_epu64(vstrip_i, vp);
				  __mmask8 m = _mm512_cmpeq_epu64_mask(vt, vr);

				  while (m > 0)
				  {
					  uint32_t id = _tzcnt_u32(m);
					  *(fbp_ptr++) = x[4 * id];
					  m = _blsr_u32(m);
				  }
			  }

			  for ( ;
				  x < smallsieve_auxbound[side][0] && *x <= p_bound; x += 4) {
				  u32_t p;

				  p = *x;
				  if (strip_i % p == x[3])
					  *(fbp_ptr++) = p;
			  }

#else
			  for (x = smallsieve_aux[side];
				  x < smallsieve_auxbound[side][0] && *x <= p_bound; x += 4) {
				  u32_t p;

				  p = *x;
				  if (strip_i % p == x[3])
					  *(fbp_ptr++) = p;
			  }
#endif

#endif
#else
			  fbp_ptr = MMX_Td(fbp_ptr, side, strip_i);
#endif
			  for (x = smallpsieve_aux[side]; x < smallpsieve_aux_ub_pow1[side]; x += 3) {
				  if (x[2] == 0) {
					  *(fbp_ptr++) = *x;
				  }
			  }
		  }

		  fbp_ptr = mpz_trialdiv(aux1, fbp_buf, fbp_ptr - fbp_buf,
			  tds_coll[fss_sv[ci]] == 0 ? "td error" : NULL);

		  if (mpz_sizeinbase(aux1, 2) <= max_factorbits[side]) {
			  n_tdsurvivors[side]++;
			  if (side == first_td_side) {
				  if (mpz_sgn(aux1) > 0)
					  mpz_set(td_rests[nc1], aux1);
				  else
					  mpz_neg(td_rests[nc1], aux1);
				  cand[nc1++] = cand[ci];
				  td_buf1[nc1] = fbp_ptr;
				  nfbp = fbp_ptr - td_buf[side];
				  continue;
			  }
			  if (mpz_sgn(aux1) < 0)mpz_neg(aux1, aux1);
#if TDS_MPQS == TDS_IMMEDIATELY
			  output_tdsurvivor(td_buf1[ci], td_buf1[ci + 1], fbp_buf, fbp_ptr,
				  td_rests[ci], aux1);
#else
#if TDS_PRIMALITY_TEST == TDS_IMMEDIATELY
			  mpz_set(large_factors[first_td_side], td_rests[ci]);
			  mpz_set(large_factors[1 - first_td_side], aux1);
			  if (primality_tests() == 1) {
				  store_tdsurvivor(td_buf1[ci], td_buf1[ci + 1], fbp_buf, fbp_ptr,
					  large_factors[first_td_side],
					  large_factors[1 - first_td_side]);
			  }
#else
			  store_tdsurvivor(td_buf1[ci], td_buf1[ci + 1], fbp_buf, fbp_ptr,
				  td_rests[ci], aux1);
#endif 
#endif 

		  }
		  else continue;

	  }
#ifndef MMX_TD
      {
	u16_t j_step;
	
	j_step= j_per_strip-last_j;
#ifdef MMX_TD
	MMX_TdUpdate(side,j_step);
#else
	{
	  u32_t i;
	  u16_t*x,*y;
	  
	  y= smalltdsieve_aux[side][j_step-1];
	  for(i= 0,x= smallsieve_aux[side];x<smallsieve_auxbound[side][0];i++,x+= 4) {
	    modulo32= x[0];
	    if(modulo32> p_bound)break;
	    x[3]= modadd32((u32_t)x[3],(u32_t)y[i]);
	  }
	}
#endif
	{
	  u16_t*x;
	  for(x= smallpsieve_aux[side];x<smallpsieve_aux_ub[side];x+= 3) {
	    modulo32= x[0];
	    x[2]= modsub32(x[2],(j_step)%modulo32);
	  }
	}
	
      }
#else
	  {
		  u16_t* x, j_step;
		  j_step = j_per_strip - last_j;
		  for (x = smallpsieve_aux[side]; x < smallpsieve_aux_ub[side]; x += 3) {
			  if ((modulo32 = x[0]))
				  x[2] = modsub32(x[2], (j_step) % modulo32);
		  }
	  }
#endif
      newclock= clock();
      tds4_clock[side]+= newclock-last_tdclock;
      last_tdclock= newclock;
      ncand= nc1;
    }
    
  }
}

#ifndef ASM_MPZ_TD

static mpz_t mpz_td_aux;
static u32_t initialized= 0;

u32_t* mpz_trialdiv(mpz_t N,u32_t*pbuf,u32_t ncp,char*errmsg)
{
  u32_t np,np1,i,e2;
  
  if(initialized == 0) {
    mpz_init(mpz_td_aux);
    initialized= 1;
  }
  e2= 0;
  while((mpz_get_ui(N)%2) == 0) {
    mpz_fdiv_q_2exp(N,N,1);
    e2++;
  }
  if(errmsg!=NULL) {
    for(i= 0,np= 0;i<ncp;i++) {
      if(mpz_fdiv_q_ui(N,N,pbuf[i])!=0)
	Schlendrian("%s : %u does not divide\n",errmsg,pbuf[i]);
      pbuf[np++]= pbuf[i];
    }
  } else {
    for(i= 0,np= 0;i<ncp;i++) {
      if(mpz_fdiv_q_ui(mpz_td_aux,N,pbuf[i]) == 0) {
	mpz_set(N,mpz_td_aux);
	pbuf[np++]= pbuf[i];
      }
    }
  }
  np1= np;
  for(i= 0;i<np1;i++) {
    while(mpz_fdiv_q_ui(mpz_td_aux,N,pbuf[i]) == 0) {
      mpz_set(N,mpz_td_aux);
      pbuf[np++]= pbuf[i];
    }
  }
  for(i= 0;i<e2;i++)
    pbuf[np++]= 2;
  return pbuf+np;
}
#endif

static void
store_tdsurvivor(fbp_buf0,fbp_buf0_ub,fbp_buf1,fbp_buf1_ub,lf0,lf1)
     u32_t*fbp_buf0,*fbp_buf1,*fbp_buf0_ub,*fbp_buf1_ub;
     mpz_t lf0,lf1;
{
  size_t n0,n1,n;
  
  if(total_ntds>=max_tds) {
    size_t i;
    if(max_tds == 0) {
      tds_fbp= xmalloc((2*MAX_TDS_INCREMENT+1)*sizeof(*tds_fbp));
      tds_fbp[0]= 0;
      tds_ab= xmalloc(2*MAX_TDS_INCREMENT*sizeof(*tds_ab));
      tds_lp= xmalloc(2*MAX_TDS_INCREMENT*sizeof(*tds_lp));
    } else {
      tds_fbp= xrealloc(tds_fbp,
			(2*(MAX_TDS_INCREMENT+max_tds)+1)*sizeof(*tds_fbp));
      tds_ab= xrealloc(tds_ab,2*(MAX_TDS_INCREMENT+max_tds)*sizeof(*tds_ab));
      tds_lp= xrealloc(tds_lp,2*(MAX_TDS_INCREMENT+max_tds)*sizeof(*tds_lp));
    }
    for(i= 2*total_ntds;i<2*(MAX_TDS_INCREMENT+max_tds);i++)mpz_init(tds_lp[i]);
    max_tds+= MAX_TDS_INCREMENT;
  }
  
  if(mpz_sizeinbase(lf0,2)> max_factorbits[first_td_side]||
     mpz_sizeinbase(lf1,2)> max_factorbits[1-first_td_side]) {
    fprintf(stderr,"large lp in store_tdsurvivor\n");
    return;
  }
  mpz_set(tds_lp[2*total_ntds],lf0);
  mpz_set(tds_lp[2*total_ntds+1],lf1);
  n0= fbp_buf0_ub-fbp_buf0;
  n1= fbp_buf1_ub-fbp_buf1;
  n= tds_fbp[2*total_ntds];
  if(n+n0+n1> tds_fbp_alloc) {
    size_t a;
    
    a= tds_fbp_alloc;
    while(a<n+n0+n1)a+= TDS_FBP_ALLOC_INCREMENT;
    if(tds_fbp_alloc == 0)tds_fbp_buffer= xmalloc(a*sizeof(*tds_fbp_buffer));
    else tds_fbp_buffer= xrealloc(tds_fbp_buffer,a*sizeof(*tds_fbp_buffer));
    tds_fbp_alloc= a;
  }
  
  
  memcpy(tds_fbp_buffer+n,fbp_buf0,n0*sizeof(*fbp_buf0));
  n+= n0;
  tds_fbp[2*total_ntds+1]= n;
  memcpy(tds_fbp_buffer+n,fbp_buf1,n1*sizeof(*fbp_buf1));
  tds_fbp[2*total_ntds+2]= n+n1;
  tds_ab[2*total_ntds]= mpz_get_sll(sr_a);
  tds_ab[2*total_ntds+1]= mpz_get_sll(sr_b);
  total_ntds++;
}

static int
primality_tests()
{
  int s;
  u16_t first_psp_side= cmdline_first_psp_side;
  int need_test[2];

  if(first_psp_side>1)
    first_psp_side= (mpz_cmp(large_factors[0],large_factors[1])<0) ? 0 : 1;
  
  for(s= 0;s<2;s++) {
    size_t nb;
    need_test[s]= 0;
    
    nb= mpz_sizeinbase(large_factors[s],2);
    if(mpz_cmp(large_factors[s],FBb_sq[s])<0) {
      if(nb<=max_primebits[s])continue;
      return 0;
    }
    if(mpz_cmp(large_factors[s],FBb_cu[s])<0) {
      if(nb<=max_primebits[s]*2) {
        need_test[s]= 1;
        continue;
      }
      return 0;
    }
    need_test[s]= 1;
  }
  for(s= 0;s<2;s++) {
    u16_t s1;
    s1= s^first_psp_side;
    if(!need_test[s1])continue;
    if(psp(large_factors[s1],1)==1)return 0;
    mpz_neg(large_factors[s1],large_factors[s1]);
  }
  return 1;
}

#if (TDS_PRIMALITY_TEST != TDS_IMMEDIATELY) && (TDS_PRIMALITY_TEST != TDS_MPQS)
static void
primality_tests_all()
{
  size_t i,j;
  
  // possibly a candidate for a future vectorized single-limb prime-test
  // routine.  As with splitting cofactors, we'd need a fair number
  // of them per call to keep that efficient.
  for(i= 0,j= 0;i<total_ntds;i++) {
    mpz_set(large_factors[first_td_side],tds_lp[2*i]);
    mpz_set(large_factors[1-first_td_side],tds_lp[2*i+1]);
    if(primality_tests() == 0)continue;
    mpz_set(tds_lp[2*j],large_factors[first_td_side]);
    mpz_set(tds_lp[2*j+1],large_factors[1-first_td_side]);
    tds_fbp[2*j+1]= tds_fbp[2*i+1];
    tds_fbp[2*j+2]= tds_fbp[2*i+2];
    tds_ab[2*j]= tds_ab[2*i];
    tds_ab[2*j+1]= tds_ab[2*i+1];
    j++;
  }
  total_ntds= j;
}
#endif

#if TDS_MPQS != TDS_IMMEDIATELY
static void
output_all_tdsurvivors()
{
  size_t i;
  
  // depending on how big total_ntds gets, this loop is a candidate
  // for use of the vectorized uecm (and/or a future vectorized tecm).
  // There needs to be a largish number for it to be a win, and with
  // uecm already pretty fast, it could be more work than it's worth.
  for(i= 0;i<total_ntds;i++) {
    mpz_set_sll(sr_a,tds_ab[2*i]);
    mpz_set_sll(sr_b,tds_ab[2*i+1]);
    output_tdsurvivor(tds_fbp_buffer+tds_fbp[2*i],
		      tds_fbp_buffer+tds_fbp[2*i+1],
		      tds_fbp_buffer+tds_fbp[2*i+1],
		      tds_fbp_buffer+tds_fbp[2*i+2],
		      tds_lp[2*i],tds_lp[2*i+1]);
  }
  total_ntds= 0;
}
#endif

static void
output_tdsurvivor(u32_t* fbp_buf0, u32_t* fbp_buf0_ub, u32_t* fbp_buf1, u32_t* fbp_buf1_ub,
	mpz_t lf0, mpz_t lf1)
{
  u32_t s,*(fbp_buffers[2]),*(fbp_buffers_ub[2]);
  u32_t nlp[2];
  clock_t cl;
  u16_t first_mpqs_side= cmdline_first_mpqs_side;
  
  s= first_td_side;
  fbp_buffers[s]= fbp_buf0;
  fbp_buffers_ub[s]= fbp_buf0_ub;
  fbp_buffers[1-s]= fbp_buf1;
  fbp_buffers_ub[1-s]= fbp_buf1_ub;
  mpz_set(large_factors[s],lf0);
  mpz_set(large_factors[1-s],lf1);
  
#if TDS_PRIMALITY_TEST == TDS_MPQS
  if(primality_tests() == 0)return;
#endif
  
  if(first_mpqs_side>1)
    first_mpqs_side = (mpz_cmp(large_factors[0],large_factors[1])>0) ? 0 : 1;
  /* Note: here, the large_factors[] values are negative, so the logic is reversed */

  cl= clock();

  for(s= 0;s<2;s++) {
    u16_t s1;
    i32_t i,nf;
    mpz_t*mf;
    
    s1= s^first_mpqs_side;
	if (mpz_sgn(large_factors[s1]) > 0) {
		if (mpz_cmp_ui(large_factors[s1], 1) == 0)
			nlp[s1] = 0;
		else {
			nlp[s1] = 1;
			mpz_set(large_primes[s1][0], large_factors[s1]);
		}
		continue;
	}
    
    mpz_neg(large_factors[s1],large_factors[s1]);

#ifdef GGNFS_MPQS
	// original mpqs code
    if((nf= mpqs_factor(large_factors[s1],max_primebits[s1],&mf))<0) {
      /* did it fail on a square? */
      mpz_sqrtrem(large_primes[s1][0],large_primes[s1][1],large_factors[s1]);
      if(mpz_sgn(large_primes[s1][1]) == 0) { /* remainder == 0? */
	mpz_set(large_primes[s1][1],large_primes[s1][0]);
	nlp[s1]= 2;
#if 0 /* this is now tested well enough, no need for a message */
	if(verbose > 1) {
	  fprintf(stderr," mpqs on a prime square ");
	  mpz_out_str(stderr,10,large_primes[s1][0]);
	  fprintf(stderr,"^2  ");
	}
#endif
	continue;
      }
      if(verbose > 1) {
      fprintf(stderr,"mpqs failed for ");
      mpz_out_str(stderr,10,large_factors[s1]);
      fprintf(stderr,"(a,b): ");
      mpz_out_str(stderr,10,sr_a);
      fprintf(stderr," ");
      mpz_out_str(stderr,10,sr_b);
      fprintf(stderr,"\n");
      }
      n_mpqsfail[s1]++;
      break;
    }
	if(nf == 0) {
		n_mpqsvain[s1]++;
		break;
	}
	gmp_printf("mpqs: %Zd = ", large_factors[s1]);
	for (i = 0; i < nf; i++)
	{
		mpz_set(large_primes[s1][i], mf[i]);
		gmp_printf("%Zd ", mf[i]);
	}
	printf("\n");
	nlp[s1]= nf;

#else

	if (mpz_sizeinbase(large_factors[s1], 2) <= 64) {
		uint64_t n64 = mpz_get_ui(large_factors[s1]);
		uint64_t f = getfactor_uecm(n64, 0, &pran);
		//printf("uecm: %lu = %lu * %lu\n", n64, f, n64 / f);
		if (f > 1)
		{
			mpz_set_ui(factor1, f);
			mpz_tdiv_q_ui(factor2, large_factors[s1], f);
			if (mpz_sizeinbase(factor1, 2) > max_primebits[s1]) {
				n_mpqsvain[s1]++;
				break;
			}
			if (mpz_sizeinbase(factor2, 2) > max_primebits[s1]) {
				n_mpqsvain[s1]++;
				break;
			}
			if (mpz_probab_prime_p(factor1, 1) == 0)
			{
				n_mpqsvain[s1]++;
				break;
			}
			if (mpz_probab_prime_p(factor2, 1) == 0)
			{
				n_mpqsvain[s1]++;
				break;
			}
			mpz_set(large_primes[s1][0], factor1);
			mpz_set(large_primes[s1][1], factor2);
			nlp[s1] = 2;
		}
		else
		{
			break;
		}
	}
	else
	{
		//printf("max_primebits[s1]=%d,sizeinbase(n)=%d\n", 
		//	max_primebits[s1], mpz_sizeinbase(large_factors[s1],2));
		if (getfactor_tecm(large_factors[s1], factor1, 
			mpz_sizeinbase(large_factors[s1], 2) / 3 - 2, &pran) > 0)
		{
			if (mpz_sizeinbase(factor1, 2) < max_primebits[s1])
			{
				mpz_tdiv_q(factor2, large_factors[s1], factor1);

				// if the remaining residue is obviously too big, we're done.
				if (mpz_sizeinbase(factor2, 2) > (max_primebits[s1] * 2))
				{
					break;
				}

				// check if the residue is prime.  could again use
				// a cheaper method.
				if (mpz_probab_prime_p(factor2, 1) > 0)
				{
					break;
				}

				// ok, so we have extracted one suitable factor, and the 
				// cofactor is not prime.  Do more work to split the cofactor,
				// which is now <= 64 bits in size.
				uint64_t q64 = mpz_get_ui(factor2);

				// todo: target this better based on expected factor size.
				uint64_t f64 = getfactor_uecm(q64, 0, &pran);
				
				if (f64 > 1)
				{
					mpz_set_ui(factor3, f64);
					mpz_tdiv_q_ui(factor2, factor2, f64);

					if (mpz_sizeinbase(factor2, 2) > max_primebits[s1]) {
						n_mpqsvain[s1]++;
						break;
					}
					if (mpz_sizeinbase(factor3, 2) > max_primebits[s1]) {
						n_mpqsvain[s1]++;
						break;
					}
					if (mpz_probab_prime_p(factor1, 1) == 0)
					{
						n_mpqsvain[s1]++;
						break;
					}
					if (mpz_probab_prime_p(factor2, 1) == 0)
					{
						n_mpqsvain[s1]++;
						break;
					}
					if (mpz_probab_prime_p(factor3, 1) == 0)
					{
						n_mpqsvain[s1]++;
						break;
					}
					mpz_set(large_primes[s1][0], factor1);
					mpz_set(large_primes[s1][1], factor2);
					mpz_set(large_primes[s1][2], factor3);
					nlp[s1] = 3;
				}
				else
				{
					break;
				}
			}
			else
			{
				// check if the factor is prime.  could again use
				// a cheaper method.
				if (mpz_probab_prime_p(factor1, 1) > 0)
				{
					// if the factor is obviously too big, give up.  This isn't a
					// failure since we haven't expended much effort yet.
					break;
				}
				else
				{
					// tecm found a composite first factor.
					// if it is obviously too big, we're done.
					if (mpz_sizeinbase(factor1, 2) > (max_primebits[s1] * 2))
					{
						break;
					}

					//gmp_printf("composite first factor %Zd\n", factor1);

					// isolate the 2nd smaller factor, and check its size.
					mpz_tdiv_q(factor2, large_factors[s1], factor1);

					if (mpz_sizeinbase(factor2, 2) > (max_primebits[s1]))
					{
						break;
					}

					// this assumes max_primebits is 32 or less...
					uint64_t q64 = mpz_get_ui(factor1);

					// todo: target this better based on expected factor size.
					uint64_t f64 = getfactor_uecm(q64, 0, &pran);

					if (f64 > 1)
					{
						mpz_set_ui(factor3, f64);
						mpz_tdiv_q_ui(factor1, factor1, f64);

						if (mpz_sizeinbase(factor1, 2) > max_primebits[s1]) {
							n_mpqsvain[s1]++;
							break;
						}
						if (mpz_sizeinbase(factor3, 2) > max_primebits[s1]) {
							n_mpqsvain[s1]++;
							break;
						}
						if (mpz_probab_prime_p(factor1, 1) == 0)
						{
							n_mpqsvain[s1]++;
							break;
						}
						if (mpz_probab_prime_p(factor2, 1) == 0)
						{
							n_mpqsvain[s1]++;
							break;
						}
						if (mpz_probab_prime_p(factor3, 1) == 0)
						{
							n_mpqsvain[s1]++;
							break;
						}
						mpz_set(large_primes[s1][0], factor1);
						mpz_set(large_primes[s1][1], factor2);
						mpz_set(large_primes[s1][2], factor3);
						nlp[s1] = 3;
					}
					else
					{
						break;
					}
				}
			}
		}
		else
		{
			// if ecm can't find a factor, give up.  This isn't a
			// failure since we haven't expended much effort yet.
			break;
			
		}
	}

	_mm256_zeroupper();
	 
#endif

  }
  mpqs_clock+= clock()-cl;
  if(s!=2)return;
  yield++;
  
  mpz_out_str(g_ofile, 10, sr_a);
  fprintf(g_ofile, ",");
  mpz_out_str(g_ofile, 10, sr_b);
  
#define OBASE -16
  /* -16 means "use upper-case" hex in gmplib >= 4.0 */
  for(s= 0;s<2;s++) {
    int num = 0;
    u32_t *x = fbp_buffers_ub[1-s];

    fprintf(g_ofile, ":");
    while (num < nlp[1-s]) {
      if (num>0) fprintf(g_ofile, ",");
      mpz_out_str(g_ofile, OBASE, large_primes[1-s][num]);
      num++;
    }
    while (x-- != fbp_buffers[1-s]) {
      if ((unsigned int)*x <1000) continue;
      if (num>0) fprintf(g_ofile, ",");
      fprintf(g_ofile, "%X", (unsigned int)*x);
      num++;
    }
  }
  fprintf(g_ofile, "\n");
}

/* 
   #ifdef OFMT_CWI
   #define CWI_LPB 0x100000
   #define OBASE 10
   {
   u32_t nlp_char[2];
   
   for(s= 0;s<2;s++) {
   u32_t*x,nlp1;
   
   for(x= fbp_buffers[s],nlp1= nlp[s];x<fbp_buffers_ub[s];x++)
   if(*x> CWI_LPB)
   nlp1++;
   if((nlp_char[s]= u32_t2cwi(nlp1)) == '\0')break;
   }
   if(s == 0) {
   errprintf("Conversion to CWI format failed\n");
   continue;
   }
   #ifdef OFMT_CWI_REVERSE
   fprintf(g_ofile,"01%c%c ",nlp_char[1],nlp_char[0]);
   #else
   fprintf(g_ofile,"01%c%c ",nlp_char[0],nlp_char[1]);
   #endif
   }
   #else
   fprintf(g_ofile,"W ");
   #define OBASE 16
   #endif
   mpz_out_str(g_ofile,OBASE,sr_a);
   fprintf(g_ofile," ");
   mpz_out_str(g_ofile,OBASE,sr_b);
   
   
   #ifndef OFMT_CWI_REVERSE
   for(s= 0;s<2;s++) {
   u32_t i,*x;
   #ifndef OFMT_CWI
   fprintf(g_ofile,"\n%c",'X'+s);
   #endif
   for(i= 0;i<nlp[s];i++) {
   fprintf(g_ofile," ");
   mpz_out_str(g_ofile,OBASE,large_primes[s][i]);
   }
   for(x= fbp_buffers[s];x<fbp_buffers_ub[s];x++) {
   #ifndef OFMT_CWI
   fprintf(g_ofile," %x",*x);
   #else
   if(*x> CWI_LPB)
   fprintf(g_ofile," %d",*x);
   #endif
   }
   }
   #else
   for(s= 0;s<2;s++) {
   u32_t i,*x;
   for(i= 0;i<nlp[1-s];i++) {
   fprintf(g_ofile," ");
   mpz_out_str(g_ofile,OBASE,large_primes[1-s][i]);
   }
   for(x= fbp_buffers[1-s];x<fbp_buffers_ub[1-s];x++) {
   if(*x> CWI_LPB)
   fprintf(g_ofile," %d",*x);
   }
   }
   #endif
   #ifndef OFMT_CWI
   fprintf(g_ofile,"\n");
   #else
   fprintf(g_ofile,";\n");
   #endif
   }
*/

#ifdef OFMT_CWI
static char u32_t2cwi(u32_t n)
{
  if(n<10)return'0'+n;
  n= n-10;
  if(n<26)return'A'+n;
  n= n-26;
  if(n<26)return'a'+n;
  return'\0';
}
#endif

#ifdef DEBUG
int mpout(mpz_t X)
{
  mpz_out_str(stdout,10,X);
  puts("");
  return 1;
}
#endif

void dumpsieve(u32_t j_offset,u32_t side)
{
  FILE*g_ofile;
  char*ofn;
  asprintf(&ofn,"sdump4e.ot%u.j%u.s%u",oddness_type,j_offset,side);
  if((g_ofile= fopen(ofn,"w")) == NULL) {
    free(ofn);
    return;
  }
  fwrite(sieve_interval,1,L1_SIZE,g_ofile);
  fclose(g_ofile);
  free(ofn);
  asprintf(&ofn,"hzsdump4e.ot%u.j%u.s%u",oddness_type,j_offset,side);
  if((g_ofile= fopen(ofn,"w")) == NULL) {
    free(ofn);
    return;
  }
  fwrite(horizontal_sievesums,1,j_per_strip,g_ofile);
  fclose(g_ofile);
  free(ofn);
}
