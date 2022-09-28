/*3:*/
// // #line 36 "gnfs-lasieve4e.w"

#include <stdio.h> 
#include <sys/types.h> 
#include <math.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <limits.h> 
#include <string.h> 
#include <time.h> 
#ifdef LINUX
#include <endian.h> 
#endif
#include <gmp.h> 
#include <signal.h> 
#include <setjmp.h> 

#include "asm/siever-config.h"
#ifndef TDS_MPQS
#define TDS_MPQS TDS_SPECIAL_Q
#endif
#ifndef TDS_PRIMALITY_TEST
#define TDS_PRIMALITY_TEST TDS_IMMEDIATELY
#endif

/*:3*//*4:*/
// // #line 61 "gnfs-lasieve4e.w"

#include "if.h"
#include "primgen32.h"
#include "asm/32bit.h"
#include "redu2.h"
#include "recurrence6.h"
#include "fbgen.h"
#include "real-poly-aux.h"
#include "gmp-aux.h"
#include "lasieve-prepn.h"

/*:4*//*5:*/
// #line 75 "gnfs-lasieve4e.w"

#define TDS_IMMEDIATELY 0
#define TDS_BIGSS 1
#define TDS_ODDNESS_CLASS 2
#define TDS_SPECIAL_Q 3

/*:5*//*6:*/
// #line 82 "gnfs-lasieve4e.w"

#define GCD_SIEVE_BOUND 10
#include "asm/siever-config.c"

#include "asm/lasched.h"
#include "asm/medsched.h"

#define L1_SIZE (1UL<<L1_BITS)

#if 0
#define ZSS_STAT
u32_t nss= 0,nzss[3]= {0,0,0};
#endif

static float
FB_bound[2],sieve_report_multiplier[2];
static u16_t sieve_min[2],max_primebits[2],max_factorbits[2];
static u32_t*(FB[2]),*(proots[2]),FBsize[2];


/*44:*/
// #line 1610 "gnfs-lasieve4e.w"

static double*(tpoly_f[2]);
#define CANDIDATE_SEARCH_STEPS 128
static unsigned char**(sieve_report_bounds[2]);
static i32_t n_srb_i,n_srb_j;

/*:44*/
// #line 102 "gnfs-lasieve4e.w"

static u32_t b,first_spq,first_spq1,first_root,last_spq,sieve_count;

static mpz_t m,N,aux1,aux2,aux3,sr_a,sr_b;

static mpz_t*(poly[2]);


double*(poly_f[2]),poly_norm[2];

i32_t poldeg[2],poldeg_max;

u32_t keep_factorbase;
#define MAX_LPFACTORS 3
static mpz_t rational_rest,algebraic_rest;
mpz_t factors[MAX_LPFACTORS];
static u32_t yield= 0,n_mpqsfail[2]= {0,0},n_mpqsvain[2]= {0,0};
static i64_t mpqs_clock= 0;

static i64_t sieve_clock= 0,sch_clock= 0,td_clock= 0,tdi_clock= 0;
static i64_t cs_clock[2]= {0,0},Schedule_clock= 0,medsched_clock= 0;
static i64_t si_clock[2]= {0,0},s1_clock[2]= {0,0};
static i64_t s2_clock[2]= {0,0},s3_clock[2]= {0,0};
static i64_t tdsi_clock[2]= {0,0},tds1_clock[2]= {0,0},tds2_clock[2]= {0,0};
static i64_t tds3_clock[2]= {0,0},tds4_clock[2]= {0,0};

char*basename;
char*input_line= NULL;
size_t input_line_alloc= 0;

/*:6*//*7:*/
// #line 133 "gnfs-lasieve4e.w"

static u32_t ncand;
static u16_t*cand;
static unsigned char*fss_sv;

/*:7*//*8:*/
// #line 139 "gnfs-lasieve4e.w"

static int tdcand_cmp(const void*x,const void*y)
{
	return(int)(*((u16_t*)x))-(int)(*((u16_t*)y));
}

/*:8*//*9:*/
// #line 155 "gnfs-lasieve4e.w"

typedef struct xFBstruct{
	u32_t p,pp,q,qq,r,l;
}*xFBptr;
static volatile xFBptr xFB[2];
static volatile u32_t xFBs[2];

/*:9*//*10:*/
// #line 181 "gnfs-lasieve4e.w"

static void xFBtranslate(u16_t*rop,xFBptr op);
static int xFBcmp(const void*,const void*);

/*:10*//*12:*/
// #line 196 "gnfs-lasieve4e.w"

static u32_t add_primepowers2xaFB(size_t*aFB_alloc_ptr,
								  u32_t pp_bound,u32_t side,u32_t p,u32_t r);

/*:12*//*13:*/
// #line 214 "gnfs-lasieve4e.w"

i32_t a0,a1,b0,b1;
#if 0
u32_t I_bits;
#endif
u32_t J_bits,i_shift,n_I,n_J;
u32_t root_no;
float sigma;

/*:13*//*14:*/
// #line 234 "gnfs-lasieve4e.w"

static u32_t oddness_type;
static u32_t n_i,n_j,i_bits,j_bits;

/*:14*//*15:*/
// #line 239 "gnfs-lasieve4e.w"

/*18:*/
// #line 461 "gnfs-lasieve4e.w"

u32_t spq_i,spq_j,spq_x;

/*:18*//*28:*/
// #line 995 "gnfs-lasieve4e.w"

u32_t fbi1[2];

/*:28*//*29:*/
// #line 1000 "gnfs-lasieve4e.w"

u32_t fbis[2];

/*:29*//*31:*/
// #line 1101 "gnfs-lasieve4e.w"

static u32_t j_per_strip,jps_bits,jps_mask,n_strips;
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

/*:31*//*32:*/
// #line 1115 "gnfs-lasieve4e.w"

static u32_t*(LPri[2]);
#define RI_SIZE 2

/*:32*//*33:*/
// #line 1120 "gnfs-lasieve4e.w"

static u32_t*(current_ij[2]);

/*:33*//*34:*/
// #line 1127 "gnfs-lasieve4e.w"

static size_t sched_alloc[2];
#define SE_SIZE 2
#define SCHEDFBI_MAXSTEP 0x10000

/*:34*//*38:*/
// #line 1338 "gnfs-lasieve4e.w"

#define USE_MEDSCHED
#ifdef USE_MEDSCHED
static u16_t**(med_sched[2]);
static u32_t*(medsched_fbi_bounds[2]);
static unsigned char*(medsched_logs[2]);
static size_t medsched_alloc[2];
static u16_t n_medsched_pieces[2];
#endif

/*:38*//*40:*/
// #line 1392 "gnfs-lasieve4e.w"

static unsigned char*sieve_interval= NULL,*(FB_logs[2]);
static unsigned char*tiny_sieve_buffer;
#define TINY_SIEVEBUFFER_SIZE 420
#define TINY_SIEVE_MIN 8
static double sieve_multiplier[2],FB_maxlog[2];
static u32_t j_offset;

/*:40*//*48:*/
// #line 1682 "gnfs-lasieve4e.w"

void do_scheduling(struct schedule_struct*,u32_t,u32_t,u32_t);

/*:48*//*51:*/
// #line 1729 "gnfs-lasieve4e.w"

static u16_t*(smallsieve_aux[2]),*(smallsieve_auxbound[2][5]);
static u16_t*(smallsieve_tinybound[2]);

/*:51*//*52:*/
// #line 1737 "gnfs-lasieve4e.w"

static u16_t*(smallsieve_aux1[2]),*(smallsieve_aux1_ub_odd[2]);
static u16_t*(smallsieve_aux1_ub[2]),*(smallsieve_tinybound1[2]);

/*:52*//*53:*/
// #line 1744 "gnfs-lasieve4e.w"

static u16_t*(smallsieve_aux2[2]),*(smallsieve_aux2_ub[2]);

/*:53*//*54:*/
// #line 1761 "gnfs-lasieve4e.w"

static u16_t*(smallpsieve_aux[2]),*(smallpsieve_aux_ub_pow1[2]);
static u16_t*(smallpsieve_aux_ub_odd[2]),*(smallpsieve_aux_ub[2]);
static unsigned char*horizontal_sievesums;

/*:54*//*55:*/
// #line 1768 "gnfs-lasieve4e.w"

static u16_t*(x2FB[2]),x2FBs[2];

/*:55*//*56:*/
// #line 1776 "gnfs-lasieve4e.w"

static u16_t*tinysieve_curpos;
#ifndef MMX_TD
static u16_t**(smalltdsieve_aux[2]);
#ifdef PREINVERT
static u32_t*(smalltd_pi[2]);
#endif
#endif

/*:56*//*57:*/
// #line 1788 "gnfs-lasieve4e.w"

#ifdef GCD_SIEVE_BOUND
static u32_t np_gcd_sieve;
static unsigned char*gcd_sieve_buffer;
static void gcd_sieve(void);
#endif

/*:57*//*95:*/
// #line 2777 "gnfs-lasieve4e.w"

u16_t**schedbuf;

/*:95*//*103:*/
// #line 2872 "gnfs-lasieve4e.w"

static void store_candidate(u16_t,u16_t,unsigned char);

/*:103*//*110:*/
// #line 3050 "gnfs-lasieve4e.w"

void trial_divide(void);

/*:110*//*124:*/
// #line 3486 "gnfs-lasieve4e.w"

#ifndef SCHED_TDS_BUFSIZE
#define SCHED_TDS_BUFSIZE 1024
#endif
u16_t*(sched_tds_buffer[SCHED_TDS_BUFSIZE]);

/*:124*//*140:*/
// #line 3896 "gnfs-lasieve4e.w"

u32_t*mpz_trialdiv(mpz_t N,u32_t*pbuf,u32_t ncp,char*errmsg);

/*:140*//*142:*/
// #line 3949 "gnfs-lasieve4e.w"

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

/*:142*//*150:*/
// #line 4229 "gnfs-lasieve4e.w"

#if 0
#define OFMT_CWI
#endif
#ifdef OFMT_CWI
static char u32_t2cwi(u32_t);
#endif

/*:150*//*153:*/
// #line 4263 "gnfs-lasieve4e.w"

void dumpsieve(u32_t j_offset,u32_t side);

/*:153*/
// #line 240 "gnfs-lasieve4e.w"

/*114:*/
// #line 3135 "gnfs-lasieve4e.w"

u32_t*(td_buf[2]),**td_buf1;
size_t td_buf_alloc[2]= {1024,1024};

/*:114*//*119:*/
// #line 3306 "gnfs-lasieve4e.w"

static unsigned char tds_coll[UCHAR_MAX];
u32_t**tds_fbi= NULL;
u32_t**tds_fbi_curpos= NULL;
#ifndef TDFBI_ALLOC
#define TDFBI_ALLOC 256
static size_t tds_fbi_alloc= TDFBI_ALLOC;
#endif

/*:119*//*138:*/
// #line 3871 "gnfs-lasieve4e.w"

static mpz_t td_rests[L1_SIZE];
static mpz_t large_factors[2],*(large_primes[2]);
static mpz_t FBb_sq[2];

/*:138*/
// #line 241 "gnfs-lasieve4e.w"


void Usage()
{
	complain("Usage");
}

static u32_t n_prereports= 0,n_reports= 0,n_rep1= 0,n_rep2= 0;
static u32_t n_tdsurvivors[2]= {0,0};
static FILE*ofile;
static char*ofile_name;

#ifdef STC_DEBUG
FILE*debugfile;
#endif

static u16_t special_q_side,first_td_side,first_sieve_side;
static u16_t first_psp_side,first_mpqs_side,append_output,exitval;
static u16_t cmdline_first_sieve_side= USHRT_MAX;
static u16_t cmdline_first_td_side= USHRT_MAX;
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

main(int argc,char**argv)
{
	u16_t zip_output,force_aFBcalc;
	u16_t catch_signals;
	u32_t all_spq_done;
	u32_t n_spq,n_spq_discard;
	double tStart, tNow, lastReport;
	char*sysload_cmd;
	u32_t process_no;
#ifdef STC_DEBUG
	debugfile= fopen("rtdsdebug","w");
#endif
	/*21:*/
	// #line 562 "gnfs-lasieve4e.w"
	
	{
		i32_t option;
		FILE*input_data;
		u32_t i;
		
		sieve_count= U32_MAX;
		ofile_name= NULL;
		zip_output= 0;
		special_q_side= NO_SIDE;
		sigma= 0;
		keep_factorbase= 0;
		basename= NULL;
		first_spq= 0;
		sieve_count= 1;
		force_aFBcalc= 0;
		sysload_cmd= NULL;
		process_no= 0;
		catch_signals= 0;
		
		first_psp_side= 2;
		first_mpqs_side= 0;
		J_bits= U32_MAX;
		
#define NumRead(x) if(sscanf(optarg,"%u",&x)!=1) Usage()
#define NumRead16(x) if(sscanf(optarg,"%hu",&x)!=1) Usage()
		
		append_output= 0;
		while((option= getopt(argc,argv,"FJ:L:M:N:P:S:ab:c:f:i:kn:o:q:rt:vz"))!=-1){
			switch(option){
			case'F':
				force_aFBcalc= 1;
				break;
			case'J':
				NumRead(J_bits);
				break;
			case'L':
				sysload_cmd= optarg;
				break;
			case'M':
				NumRead16(first_mpqs_side);
				break;
			case'P':
				NumRead16(first_psp_side);
				break;
			case'S':
				if(sscanf(optarg,"%f",&sigma)!=1){
					errprintf("Cannot read floating point number %s\n",optarg);
					Usage();
				}
				break;
			case'a':
				if(special_q_side!=NO_SIDE){
					errprintf("Ignoring -a\n");
					break;
				}
				special_q_side= ALGEBRAIC_SIDE;
				break;
			case'b':
				if(basename!=NULL)errprintf("Ignoring -b %s\n",basename);
				else basename= optarg;
				break;
			case'c':
				NumRead(sieve_count);
				break;
			case'f':
				if(sscanf(optarg,"%u:%u:%u",&first_spq,&first_spq1,&first_root)!=3){
					if(sscanf(optarg,"%u",&first_spq)==1){
						first_spq1= first_spq;
						first_root= 0;
					}else Usage();
				}else append_output= 1;
				break;
			case'i':
				if(sscanf(optarg,"%hu",&cmdline_first_sieve_side)!=1)
					complain("-i %s ???\n",optarg);
				break;
			case'k':
				keep_factorbase= 1;
				break;
			case'n':
				catch_signals= 1;
			case'N':
				NumRead(process_no);
				break;
			case'o':
				ofile_name= optarg;
				break;
			case'q':
				NumRead16(special_q_side);
				break;
			case'r':
				if(special_q_side!=NO_SIDE){
					errprintf("Ignoring -r\n");
					break;
				}
				special_q_side= RATIONAL_SIDE;
				break;
			case't':
				if(sscanf(optarg,"%hu",&cmdline_first_td_side)!=1)
					complain("-t %s ???\n",optarg);
				break;
			case'v':
				verbose++;
				break;
			case'z':
				zip_output= 1;
				break;
			}
		}
		if(J_bits==U32_MAX)J_bits= I_bits-1;
		if(first_psp_side==2)first_psp_side= first_mpqs_side;
#ifndef I_bits
#error Must #define I_bits
#endif
		last_spq= first_spq+sieve_count;
		if(last_spq>=I32_MAX/2){
			
			
			complain("Cannot handle special q >= %d\n",I32_MAX/2);
		}
		if(optind<argc&&basename==NULL){
			basename= argv[optind];
			optind++;
		}
		if(optind<argc)fprintf(stderr,"Ignoring %u trailing command line args\n",
			argc-optind);
		if(basename==NULL)basename= "gnfs";
		if((input_data= fopen(basename,"r"))==NULL){
			complain("Cannot open %s for input of nfs polynomials: %m\n",basename);
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
		// input_poly(N,poly,poldeg,poly+1,poldeg+1,m,input_data);
#if 0
		if(poldeg[1]> 1){
			if(poldeg[0]==1){
				mpz_t*X;
				poldeg[0]= poldeg[1];
				poldeg[1]= 1;
				X= poly[0];
				poly[0]= poly[1];
				poly[1]= X;
			}else{
				complain("Degrees >1 on both sides not implemented\n");
			}
		}
#endif
		//skip_blanks_comments(&input_line,&input_line_alloc,input_data);
		close(input_data);
		{ FILE *fp;
		char token[256], value[256], thisLine[1024];
		
		
		sieve_min[0] = sieve_min[1]=0;
		
		if (!(fp = fopen(basename, "rb"))) {
			printf("Error opening %s for read!\n", basename);
			return -1;
		}
		input_poly(N, poly, poldeg, poly + 1, poldeg + 1, m, fp);
		rewind(fp);
		while (!(feof(fp))) {
			thisLine[0] = 0;
			fgets(thisLine, 1023, fp);
			/* Special case: If there's a polynomial, handle it seperately: */
			if (strncmp(thisLine, "START_POLY", 10)==0) {
				while (!(feof(fp)) && strncmp(thisLine, "END_POLY", 8)) 
					fgets(thisLine, 1023, fp);
			} else  if ((sscanf(thisLine, "%255s %255s", token, value)==2) && 
                (thisLine[0] != '#')) {
				
				token[sizeof(token)-1] = 0;
				if (strncmp(token, "skew:", 5)==0) {
					sigma = (float)atof(value);
				} else if (strncmp(token, "q0:", 3)==0) {
					first_spq = atol(value);
				} else if (strncmp(token, "qintsize:", 9)==0) {
					sieve_count = atol(value);
				} else if ((strncmp(token, "skip0:", 6)==0) ||
					(strncmp(token, "askip:", 6)==0)) {
					sieve_min[0] = atol(value);
				} else if ((strncmp(token, "skip1:", 6)==0) ||
					(strncmp(token, "rskip:", 6)==0)) {
					sieve_min[1] = atol(value);
				} else if ((strncmp(token, "lim0:", 5)==0) ||
					(strncmp(token, "alim:", 5)==0)) {
					FB_bound[0] = (float)atol(value);
				} else if ((strncmp(token, "lim1:", 5)==0)||
					(strncmp(token, "rlim:", 5)==0)) {
					FB_bound[1] = (float)atof(value);
				} else if ((strncmp(token, "lpb0:", 5)==0) ||
					(strncmp(token, "lpba:", 5)==0)) {
					max_primebits[0] = atoi(value);
				} else if ((strncmp(token, "lpb1:", 5)==0) ||
					(strncmp(token, "lpbr:", 5)==0)) {
					max_primebits[1] = atoi(value);
				} else if ((strncmp(token, "mfb0:", 5)==0) ||
					(strncmp(token, "mfba:", 5)==0)) {
					max_factorbits[0] = atoi(value);
				} else if ((strncmp(token, "mfb1:", 5)==0) ||
					(strncmp(token, "mfbr:", 5)==0)) {
					value[sizeof(value)-1] = 0;
					max_factorbits[1] = atoi(value);
				} else if ((strncmp(token, "lambda0:", 8)==0) ||
					(strncmp(token, "alambda:", 8)==0)) {
					sieve_report_multiplier[0] = (float)atof(value);
				} else if ((strncmp(token, "lambda1:", 8)==0) ||
					(strncmp(token, "rlambda:", 8)==0)) {
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
		for(i= 0;i<2;i++){
			if(FB_bound[i]<4||sieve_report_multiplier[i]<=0){
				complain("Please set all bounds to reasonable values!\n");
			}
#if 1
			if(max_primebits[i]> 33){
				complain("Only large primes up to 33 bits are allowed.\n");
			}
#endif
		}
		if(sieve_count!=0){
			if(sigma==0)complain("Please set a skewness\n");
			if(special_q_side==NO_SIDE){
				errprintf("Please use -a or -r\n");
				Usage();
			}
			if(FB_bound[special_q_side]> first_spq){
				complain("Special q lower bound %u below rFB bound %g\n",
					first_spq,FB_bound[special_q_side]);
			}
		}
		fclose(input_data);
		if(poldeg[0]<poldeg[1])poldeg_max= poldeg[1];
		else poldeg_max= poldeg[0];
		
		
		
		i_shift= 1<<(I_bits-1);
		n_I= 1<<I_bits;
		n_i= n_I/2;
		i_bits= I_bits-1;
		n_J= 1<<J_bits;
		n_j= n_J/2;
		j_bits= J_bits-1;
		/*22:*/
		// #line 775 "gnfs-lasieve4e.w"
		
		{
			u32_t i,j;
			double x,y,z;
			
			x= sqrt(first_spq*sigma)*n_I;
			y= x/sigma;
			for(j= 0;j<2;j++){
				poly_f[j]= xmalloc((poldeg[j]+1)*sizeof(*poly_f[j]));
				
				for(i= 0,z= 1,poly_norm[j]= 0;
				i<=poldeg[j];i++){
					poly_f[j][i]= mpz_get_d(poly[j][i]);
					poly_norm[j]= poly_norm[j]*y+fabs(poly_f[j][i])*z;
					z*= x;
				}
			}
		}
		
		/*:22*/
		// #line 771 "gnfs-lasieve4e.w"
		
}

/*:21*/
// #line 301 "gnfs-lasieve4e.w"

siever_init();
/*23:*/
// #line 795 "gnfs-lasieve4e.w"

if(sieve_count!=0){
	
	if(ofile_name==NULL){
		if(zip_output==0){
			asprintf(&ofile_name,"%s.lasieve-%u.%u-%u",basename,
				special_q_side,first_spq,last_spq);
		}else{
			asprintf(&ofile_name,
				append_output==0?
				"gzip --best --stdout > %s.lasieve-%u.%u-%u.gz":
			"gzip --best --stdout >> %s.lasieve-%u.%u-%u.gz",
				basename,special_q_side,first_spq,last_spq);
		}
	}else{
		if(strcmp(ofile_name,"-")==0){
			if(zip_output==0){
				ofile= stdout;
				ofile_name= "to stdout";
				goto done_opening_output;
			}else ofile_name= "gzip --best --stdout";
		}else{
			if(fnmatch("*.gz",ofile_name,0)==0){
				char*on1;
				
				zip_output= 1;
				on1= strdup(ofile_name);
				asprintf(&ofile_name,"gzip --best --stdout > %s",on1);
				free(on1);
			}else zip_output= 0;
		}
	}
	if(zip_output==0){
		if(append_output> 0){
			ofile= fopen(ofile_name,"a");
		}else ofile= fopen(ofile_name,"w");
		if(ofile==NULL)complain("Cannot open %s for output: %m\n",ofile_name);
	}else{
		if((ofile= popen(ofile_name,"w"))==NULL)
			complain("Cannot exec %s for output: %m\n",ofile_name);
	}
done_opening_output:
	//fprintf(ofile,"F 0 X %u 1\n",poldeg[0]);
	fprintf(ofile,"");
}

/*:23*/
// #line 303 "gnfs-lasieve4e.w"

/*24:*/
// #line 845 "gnfs-lasieve4e.w"

{
	size_t FBS_alloc= 4096;
	u32_t prime;
	pr32_struct ps;
	char*afbname;
	FILE*afbfile;
	u32_t side;
	
	initprime32(&ps);
	
	for(side= 0;side<2;side++){
		if(poldeg[side]==1){
			FB[side]= xmalloc(FBS_alloc*sizeof(u32_t));
			proots[side]= xmalloc(FBS_alloc*sizeof(u32_t));
			prime= firstprime32(&ps);
			for(prime= nextprime32(&ps),fbi1[side]= 0,FBsize[side]= 0;
			prime<FB_bound[side];prime= nextprime32(&ps)){
				u32_t x;
				x= mpz_fdiv_ui(poly[side][1],prime);
				if(x> 0){
					modulo32= prime;
					x= modmul32(modinv32(x),mpz_fdiv_ui(poly[side][0],prime));
					x= x> 0?prime-x:0;
				}else x= prime;
				if(prime<L1_SIZE)fbi1[side]= FBsize[side];
				if(prime<n_i)fbis[side]= FBsize[side];
				if(FBsize[side]==FBS_alloc){
					FBS_alloc*= 2;
					FB[side]= xrealloc(FB[side],FBS_alloc*sizeof(u32_t));
					proots[side]= xrealloc(proots[side],FBS_alloc*sizeof(u32_t));
				}
				proots[side][FBsize[side]]= x;
				FB[side][FBsize[side]++]= prime;
			}
			proots[side]= xrealloc(proots[side],FBsize[side]*sizeof(u32_t));
			FB[side]= xrealloc(FB[side],FBsize[side]*sizeof(u32_t));
			fbi1[side]++;
			fbis[side]++;
			if(fbi1[side]<fbis[side])fbi1[side]= fbis[side];
		}else{
			u32_t j,k,l;
			asprintf(&afbname,"%s.afb.%u",basename,side);
			if(force_aFBcalc> 0||(afbfile= fopen(afbname,"r"))==NULL){
				/*25:*/
				// #line 930 "gnfs-lasieve4e.w"
				
				u32_t*root_buffer;
				size_t aFB_alloc;
				
				root_buffer= xmalloc(poldeg[side]*sizeof(*root_buffer));
				aFB_alloc= 4096;
				FB[side]= xmalloc(aFB_alloc*sizeof(**FB));
				proots[side]= xmalloc(aFB_alloc*sizeof(**proots));
				for(prime= firstprime32(&ps),FBsize[side]= 0;
				prime<FB_bound[side];prime= nextprime32(&ps)){
					u32_t i,nr;
					
					nr= root_finder(root_buffer,poly[side],poldeg[side],prime);
					for(i= 0;i<nr;i++){
						if(aFB_alloc<=FBsize[side]){
							aFB_alloc*= 2;
							FB[side]= xrealloc(FB[side],aFB_alloc*sizeof(**FB));
							proots[side]= xrealloc(proots[side],aFB_alloc*sizeof(**proots));
						}
						FB[side][FBsize[side]]= prime;
						proots[side][FBsize[side]]= root_buffer[i];
						if(prime> 2)FBsize[side]++;
					}
				}
				FB[side]= xrealloc(FB[side],FBsize[side]*sizeof(**FB));
				proots[side]= xrealloc(proots[side],FBsize[side]*sizeof(**proots));
				free(root_buffer);
				
				/*:25*/
				// #line 889 "gnfs-lasieve4e.w"
				
				if(keep_factorbase> 0)/*27:*/
					// #line 975 "gnfs-lasieve4e.w"
					
				{
					if((afbfile= fopen(afbname,"w"))==NULL){
						complain("Cannot open %s for output of aFB: %m\n",afbname);
					}
					if(write_u32(afbfile,&(FBsize[side]),1)!=1){
						complain("Cannot write aFBsize to %s: %m\n",afbname);
					}
					if(write_u32(afbfile,FB[side],FBsize[side])!=FBsize[side]||
						write_u32(afbfile,proots[side],FBsize[side])!=FBsize[side]){
						complain("Cannot write aFB to %s: %m\n",afbname);
					}
					if(write_u32(afbfile,&xFBs[side],1)!=1){
						complain("Cannot write aFBsize to %s: %m\n",afbname);
					}
					fclose(afbfile);
				}
				
				/*:27*/
				// #line 890 "gnfs-lasieve4e.w"
				
			}else{
				/*26:*/
				// #line 959 "gnfs-lasieve4e.w"
				
				if(read_u32(afbfile,&(FBsize[side]),1)!=1){
					complain("Cannot read aFB size from %s: %m\n",afbname);
				}
				FB[side]= xmalloc(FBsize[side]*sizeof(u32_t));
				proots[side]= xmalloc(FBsize[side]*sizeof(u32_t));
				if(read_u32(afbfile,FB[side],FBsize[side])!=FBsize[side]||
					read_u32(afbfile,proots[side],FBsize[side])!=FBsize[side]){
					complain("Cannot read aFB from %s: %m\n",afbname);
				}
				if(read_u32(afbfile,&xFBs[side],1)!=1){
					complain("%s: Cannot read xFBsize\n",afbname);
				}
				fclose(afbfile);
				
				/*:26*/
				// #line 892 "gnfs-lasieve4e.w"
				
			}
			
			
			for(j= 0,k= 0,l= 0;j<FBsize[side];j++){
				if(FB[side][j]<L1_SIZE)k= j;
				if(FB[side][j]<n_i)l= j;
				if(FB[side][j]> L1_SIZE&&FB[side][j]> n_I)break;
			}
			if(FBsize[side]> 0){
				if(k<l)k= l;
				fbis[side]= l+1;
				fbi1[side]= k+1;
			}else{
				fbis[side]= 0;
				fbi1[side]= 0;
			}
		}
}

{
	u32_t i,srfbs,safbs;
	
	for(i= 0,srfbs= 0;i<xFBs[1];i++){
		if(xFB[1][i].p==xFB[1][i].pp)srfbs++;
	}
	for(i= 0,safbs= 0;i<xFBs[0];i++){
		if(xFB[0][i].p==xFB[0][i].pp)safbs++;
	}
	logbook(0,"FBsize %u+%u (deg %u), %u+%u (deg %u)\n",
		FBsize[0],safbs,poldeg[0],FBsize[1],srfbs,poldeg[1]);
}
free(afbname);

/*45:*/
// #line 1617 "gnfs-lasieve4e.w"

{
	u32_t i;
	size_t si,sj;
	
	n_srb_i= 2*((n_i+2*CANDIDATE_SEARCH_STEPS-1)/(2*CANDIDATE_SEARCH_STEPS));
	n_srb_j= (n_J+2*CANDIDATE_SEARCH_STEPS-1)/(2*CANDIDATE_SEARCH_STEPS);
	sj= n_srb_j*sizeof(*(sieve_report_bounds[0]));
	si= n_srb_i*sizeof(**(sieve_report_bounds[0]));
	for(i= 0;i<2;i++){
		u32_t j;
		
		tpoly_f[i]= xmalloc((1+poldeg[i])*sizeof(**tpoly_f));
		sieve_report_bounds[i]= xmalloc(sj);
		for(j= 0;j<n_srb_j;j++)
			sieve_report_bounds[i][j]= xmalloc(si);
	}
}

/*:45*/
// #line 926 "gnfs-lasieve4e.w"

}

/*:24*/
// #line 304 "gnfs-lasieve4e.w"

if(sieve_count==0)exit(0);
/*30:*/
// #line 1004 "gnfs-lasieve4e.w"

{
	u32_t side,i;
	
	for(side= 0;side<2;side++){
		u32_t prime,nr;
		struct xFBstruct*s;
		u32_t*root_buffer;
		size_t xaFB_alloc= 0;
		FB_logs[side]= xmalloc(FBsize[side]);
		sieve_multiplier[side]= (UCHAR_MAX-50)/log(poly_norm[side]);
		
		root_buffer= xmalloc(poldeg[side]*sizeof(*root_buffer));
		prime= 2;
		nr= root_finder(root_buffer,poly[side],poldeg[side],prime);
		
		for(i= 0;i<nr;i++){
			adjust_bufsize((void**)&(xFB[side]),&xaFB_alloc,1+xFBs[side],
				16,sizeof(**xFB));
			s= xFB[side]+xFBs[side];
			s->p= prime;
			s->pp= prime;
			if(root_buffer[i]==prime){
				s->qq= prime;
				s->q= 1;
				s->r= 1;
			}else{
				s->qq= 1;
				s->q= prime;
				s->r= root_buffer[i];
			}
			xFBs[side]++;
			add_primepowers2xaFB(&xaFB_alloc,n_I,side,0,0);
		}
		free(root_buffer);
		for(i= 0;i<FBsize[side];i++){
			double l;
			u32_t l1;
			
			prime= FB[side][i];
			if(prime> n_I/prime)break;
			l= log(prime);
			l1= add_primepowers2xaFB(&xaFB_alloc,n_I,side,prime,proots[side][i]);
			FB_logs[side][i]= rint(l1*l*sieve_multiplier[side]);
		}
		while(i<FBsize[side]){
			double l;
			
			l= log(FB[side][i]);
			if(l> FB_maxlog[side])FB_maxlog[side]= l;
			FB_logs[side][i++]= rint(sieve_multiplier[side]*l);
		}
		FB_maxlog[side]*= sieve_multiplier[side];
		qsort(xFB[side],xFBs[side],sizeof(*(xFB[side])),xFBcmp);
	}
}

/*:30*/
// #line 306 "gnfs-lasieve4e.w"

/*35:*/
// #line 1133 "gnfs-lasieve4e.w"

sieve_interval= xvalloc(L1_SIZE);
cand= xvalloc(L1_SIZE*sizeof(*cand));
fss_sv= xvalloc(L1_SIZE);
tiny_sieve_buffer= xmalloc(TINY_SIEVEBUFFER_SIZE);
if(n_i> L1_SIZE)
complain("Strip length %u exceeds L1 size %u\n",n_i,L1_SIZE);
j_per_strip= L1_SIZE/n_i;
jps_bits= L1_BITS-i_bits;
jps_mask= j_per_strip-1;
if(j_per_strip!=1<<jps_bits)
Schlendrian("Expected %u j per strip, calculated %u\n",
			j_per_strip,1<<jps_bits);
n_strips= n_j>>(L1_BITS-i_bits);
rec_info_init(n_i,n_j);
/*58:*/
// #line 1796 "gnfs-lasieve4e.w"

{
	u32_t s;
#define MAX_TINY_2POW 4
	
	if(poldeg[0]<poldeg[1])s= poldeg[1];
	else s= poldeg[0];
	tinysieve_curpos= xmalloc(TINY_SIEVE_MIN*s*sizeof(*tinysieve_curpos));
	horizontal_sievesums= xmalloc(j_per_strip*sizeof(*horizontal_sievesums));
	for(s= 0;s<2;s++){
		u32_t fbi;
		size_t maxent;
		
		smallsieve_aux[s]= xmalloc(4*fbis[s]*sizeof(*(smallsieve_aux[s])));
#ifndef MMX_TD
#ifdef PREINVERT
		smalltd_pi[s]= xmalloc(fbis[s]*sizeof(*(smalltd_pi[s])));
#endif
		smalltdsieve_aux[s]= xmalloc(j_per_strip*sizeof(*(smalltdsieve_aux[s])));
		for(fbi= 0;fbi<j_per_strip;fbi++)
			smalltdsieve_aux[s][fbi]= 
			xmalloc(fbis[s]*sizeof(**(smalltdsieve_aux[s])));
#else
		
		MMX_TdAllocate(j_per_strip,fbis[0],fbis[1]);
#endif
		smallsieve_aux1[s]= xmalloc(6*xFBs[s]*sizeof(*(smallsieve_aux1[s])));
		
		
		maxent= fbis[s];
		maxent+= xFBs[s];
		smallpsieve_aux[s]= xmalloc(3*maxent*sizeof(*(smallpsieve_aux[s])));
		maxent= 0;
		for(fbi= 0;fbi<xFBs[s];fbi++){
			if(xFB[s][fbi].p==2)
				maxent++;
		}
		smallsieve_aux2[s]= xmalloc(4*maxent*sizeof(*(smallsieve_aux2[s])));
		x2FB[s]= xmalloc(maxent*6*sizeof(*(x2FB[s])));
	}
}

/*:58*//*59:*/
// #line 1839 "gnfs-lasieve4e.w"

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

/*:59*/
// #line 1148 "gnfs-lasieve4e.w"

{
	u32_t s;
	for(s= 0;s<2;s++){
		if(sieve_min[s]<TINY_SIEVE_MIN&&sieve_min[s]!=0){
			errprintf("Sieving with all primes on side %u since\n",s);
			errprintf("tiny sieve procedure is being used\n");
			sieve_min[s]= 0;
		}
		current_ij[s]= xmalloc(FBsize[s]*sizeof(*current_ij[s]));
		LPri[s]= xmalloc(FBsize[s]*sizeof(**LPri)*RI_SIZE);
	}
}

/*:35*//*36:*/
// #line 1175 "gnfs-lasieve4e.w"

{
	u32_t s;
	size_t total_alloc;
	u16_t*sched_buf;
	double pvl_max[2];
	
	total_alloc= 0;
	for(s= 0;s<2;s++){
		u32_t i,fbi_lb;
		
		if(sigma>=1)pvl_max[s]= poldeg[s]*log(last_spq*sqrt(sigma));
		else pvl_max[s]= poldeg[s]*log(last_spq/sqrt(sigma));
		pvl_max[s]+= log(poly_norm[s]);
		if(fbi1[s]>=FBsize[s]||i_bits+j_bits<=L1_BITS){
			n_schedules[s]= 0;
			continue;
		}
		for(i= 0;i<N_PRIMEBOUNDS;i++)
			if(FB_bound[s]<=schedule_primebounds[i]||
				i_bits+j_bits<=schedule_sizebits[i])
				break;
			n_schedules[s]= i+1;
			schedules[s]= xmalloc(n_schedules[s]*sizeof(**schedules));
			fbi_lb= fbi1[s];
			for(i= 0;i<n_schedules[s];i++){
				u32_t fbp_lb,fbp_ub;
				u32_t fbi,fbi_ub;
				u32_t sp_i;
				u32_t n,sl_i;
				u32_t ns;
				size_t allocate,all1;
				
				if(i==n_schedules[s]-1)fbp_ub= FB_bound[s];
				else fbp_ub= schedule_primebounds[i];
				if(i==0)fbp_lb= FB[s][fbi1[s]];
				else fbp_lb= schedule_primebounds[i-1];
				
				if(i_bits+j_bits<schedule_sizebits[i])ns= 1<<(i_bits+j_bits-L1_BITS);
				else ns= 1<<(schedule_sizebits[i]-L1_BITS);
				schedules[s][i].n_strips= ns;
				
				
				
#ifndef SCHED_TOL
#ifndef NO_SCHEDTOL
#define SCHED_TOL 2
#endif
#endif
#ifdef SCHED_TOL
				allocate= rint(SCHED_TOL*n_i*j_per_strip*log(log(fbp_ub)/log(fbp_lb)));
#else
				allocate= rint(sched_tol[i]*n_i*j_per_strip*log(log(fbp_ub)/log(fbp_lb)));
#endif
				allocate*= SE_SIZE;
				
				
				
				all1= allocate+n_i*ceil(pvl_max[s]/log(fbp_lb))*SE_SIZE;
				schedules[s][i].alloc= allocate;
				schedules[s][i].alloc1= all1;
				
				for(n= 0,fbi= fbi_lb;fbi<FBsize[s];){
					u32_t fbi_ub1;
					fbi_ub1= fbi+SCHEDFBI_MAXSTEP;
					if(fbi_ub1>=FBsize[s])fbi_ub1= FBsize[s];
					else{
						if(FB[s][fbi_ub1]> fbp_ub){
							while(FB[s][fbi_ub1]> fbp_ub)
								fbi_ub1--;
							fbi_ub1++;
						}
					}
					if(FB_logs[s][fbi]==FB_logs[s][fbi_ub1-1]){
						n++;
						fbi= fbi_ub1;
					}else{
						u32_t l;
						n+= FB_logs[s][fbi_ub1-1]-FB_logs[s][fbi];
						fbi= fbi_ub1-1;
						l= FB_logs[s][fbi];
						while(FB_logs[s][fbi]==l)fbi--;
						fbi++;
					}
					if(fbi>=FBsize[s]||FB[s][fbi]> fbp_ub)break;
				}
				fbi_ub= fbi;
				schedules[s][i].n_pieces= n;
				n++;
				schedules[s][i].schedule= xmalloc(n*sizeof(*(schedules[s][i].schedule)));
				for(sl_i= 0;sl_i<n;sl_i++)
					schedules[s][i].schedule[sl_i]= 
					xmalloc(ns*sizeof(**(schedules[s][i].schedule)));
				schedules[s][i].schedule[0][0]= (u16_t*)total_alloc;
				total_alloc+= all1;
				for(sp_i= 1;sp_i<ns;sp_i++){
					schedules[s][i].schedule[0][sp_i]= (u16_t*)total_alloc;
					total_alloc+= allocate;
				}
				schedules[s][i].fbi_bounds= 
					xmalloc(n*sizeof(*(schedules[s][i].fbi_bounds)));
				schedules[s][i].schedlogs= xmalloc(n);
				for(n= 0,fbi= fbi_lb;fbi<fbi_ub;){
					u32_t fbi_ub1;
					fbi_ub1= fbi+SCHEDFBI_MAXSTEP;
					if(fbi_ub1> fbi_ub)fbi_ub1= fbi_ub;
					if(FB_logs[s][fbi]==FB_logs[s][fbi_ub1-1]){
						schedules[s][i].fbi_bounds[n++]= fbi;
						fbi= fbi_ub1;
					}else{
						u32_t l,lmax;
						
						lmax= FB_logs[s][fbi_ub1-1];
						for(l= FB_logs[s][fbi];l<lmax;l++){
							schedules[s][i].fbi_bounds[n++]= fbi;
							while(fbi<fbi_ub&&FB_logs[s][fbi]==l)fbi++;
						}
					}
				}
				if(n!=schedules[s][i].n_pieces)
					Schlendrian("Expected %u schedule pieces on side %u, have %u\n",
					schedules[s][i].n_pieces,s,n);
				schedules[s][i].fbi_bounds[n++]= fbi;
				for(n= 0;n<schedules[s][i].n_pieces;n++)
					schedules[s][i].schedlogs[n]= FB_logs[s][schedules[s][i].fbi_bounds[n]];
				schedules[s][i].ri= 
					LPri[s]+(schedules[s][i].fbi_bounds[0]-fbis[s])*RI_SIZE;
				fbi_lb= fbi_ub;
}
}
/*37:*/
// #line 1323 "gnfs-lasieve4e.w"

sched_buf= xmalloc((total_alloc+65536*SE_SIZE*j_per_strip)*
				   sizeof(**((**schedules).schedule)));
for(s= 0;s<2;s++){
	u32_t i;
	for(i= 0;i<n_schedules[s];i++){
		u32_t sp_i;
		
		for(sp_i= 0;sp_i<schedules[s][i].n_strips;sp_i++)
			schedules[s][i].schedule[0][sp_i]= 
			sched_buf+(size_t)(schedules[s][i].schedule[0][sp_i]);
	}
}

/*:37*/
// #line 1305 "gnfs-lasieve4e.w"

/*39:*/
// #line 1349 "gnfs-lasieve4e.w"

#ifdef USE_MEDSCHED
{
	u32_t s;
	
	for(s= 0;s<2;s++){
		if(fbis[s]<fbi1[s]){
			u32_t fbi;
			u32_t n;
			unsigned char oldlog;
			
			
			medsched_alloc[s]= j_per_strip*(fbi1[s]-fbis[s])*SE_SIZE;
			
			medsched_alloc[s]+= n_i*ceil(pvl_max[s]/log(n_i))*SE_SIZE;
			n_medsched_pieces[s]= 1+FB_logs[s][fbi1[s]-1]-FB_logs[s][fbis[s]];
			med_sched[s]= xmalloc((1+n_medsched_pieces[s])*sizeof(**med_sched));
			med_sched[s][0]= xmalloc(medsched_alloc[s]*sizeof(***med_sched));
			
			medsched_fbi_bounds[s]= 
				xmalloc((1+n_medsched_pieces[s])*sizeof(**medsched_fbi_bounds));
			medsched_logs[s]= xmalloc(n_medsched_pieces[s]);
			
			for(n= 0,fbi= fbis[s],oldlog= UCHAR_MAX;fbi<fbi1[s];fbi++){
				if(FB_logs[s][fbi]!=oldlog){
					medsched_fbi_bounds[s][n]= fbi;
					oldlog= FB_logs[s][fbi];
					medsched_logs[s][n++]= oldlog;
				}
			}
			if(n!=n_medsched_pieces[s])
				Schlendrian("Expected %u medium schedule pieces on side %u, have %u\n",
				n_medsched_pieces[s],s,n);
			medsched_fbi_bounds[s][n]= fbi;
		}else{
			
			n_medsched_pieces[s]= 0;
		}
	}
}
#endif

/*:39*/
// #line 1306 "gnfs-lasieve4e.w"

}

/*:36*//*96:*/
// #line 2781 "gnfs-lasieve4e.w"

{
	u32_t s;
	size_t schedbuf_alloc;
	
	for(s= 0,schedbuf_alloc= 0;s<2;s++){
		u32_t i;
		
		for(i= 0;i<n_schedules[s];i++)
			if(schedules[s][i].n_pieces> schedbuf_alloc)
				schedbuf_alloc= schedules[s][i].n_pieces;
	}
	schedbuf= xmalloc((1+schedbuf_alloc)*sizeof(*schedbuf));
}

/*:96*/
// #line 307 "gnfs-lasieve4e.w"

/*115:*/
// #line 3140 "gnfs-lasieve4e.w"

td_buf1= xmalloc((1+L1_SIZE)*sizeof(*td_buf1));
td_buf[0]= xmalloc(td_buf_alloc[0]*sizeof(**td_buf));
td_buf[1]= xmalloc(td_buf_alloc[1]*sizeof(**td_buf));

/*:115*//*120:*/
// #line 3316 "gnfs-lasieve4e.w"

{
	u32_t i;
	if(tds_fbi==NULL){
		tds_fbi= xmalloc(UCHAR_MAX*sizeof(*tds_fbi));
		tds_fbi_curpos= xmalloc(UCHAR_MAX*sizeof(*tds_fbi));
		for(i= 0;i<UCHAR_MAX;i++)
			tds_fbi[i]= xmalloc(tds_fbi_alloc*sizeof(**tds_fbi));
	}
}

/*:120*//*139:*/
// #line 3877 "gnfs-lasieve4e.w"

{
	u32_t s,i;
	for(i= 0;i<L1_SIZE;i++){
		mpz_init(td_rests[i]);
	}
	for(s= 0;s<2;s++){
		mpz_init(large_factors[s]);
		large_primes[s]= xmalloc(max_factorbits[s]*sizeof(*(large_primes[s])));
		for(i= 0;i<max_factorbits[s];i++){
			mpz_init(large_primes[s][i]);
		}
		mpz_init_set_d(FBb_sq[s],FB_bound[s]);
		mpz_mul(FBb_sq[s],FBb_sq[s],FBb_sq[s]);
	}
}


/*:139*/
// #line 308 "gnfs-lasieve4e.w"

all_spq_done= 1;
/*16:*/
// #line 325 "gnfs-lasieve4e.w"

{
	u32_t*r;
	initprime32(&special_q_ps);
	last_clock= clock();
	n_spq= 0;
	n_spq_discard= 0;
	r= xmalloc(poldeg_max*sizeof(*r));
	special_q= pr32_seek(&special_q_ps,first_spq1);
	if(catch_signals!=0){
		signal(SIGTERM,terminate_sieving);
		signal(SIGINT,terminate_sieving);
	}
	for(;
	special_q<last_spq&&special_q!=0;
	special_q= nextprime32(&special_q_ps),first_root= 0){
		u32_t nr;
		
		special_q_log= log(special_q);
		if(cmdline_first_sieve_side==USHRT_MAX){
#if 1
			double nn[2];
			u32_t s;
#if 0
			for(s= 0;s<2;s++){
				nn[s]= log(poly_norm[s]*(special_q_side==s?1:special_q));
				nn[s]= nn[s]/log(FB_bound[s])-sieve_report_multiplier[s];
			}
#else
			for(s= 0;s<2;s++){
				nn[s]= log(poly_norm[s]*(special_q_side==s?1:special_q));
				nn[s]= nn[s]/(sieve_report_multiplier[s]*log(FB_bound[s]));
			}
#endif
			if(nn[0]<nn[1])first_sieve_side= 1;
			else first_sieve_side= 0;
#else
			if(poly_norm[0]*(special_q_side==0?1:special_q)
				<poly_norm[1]*(special_q_side==1?1:special_q)){
				first_sieve_side= 1;
			}else{
				first_sieve_side= 0;
			}
#endif
		}else{
			first_sieve_side= cmdline_first_sieve_side;
			if(first_sieve_side>=2)complain("First sieve side must not be %u\n",
				(u32_t)first_sieve_side);
		}
		logbook(1,"First sieve side: %u\n",(u32_t)first_sieve_side);
		if(cmdline_first_td_side!=USHRT_MAX)first_td_side= cmdline_first_td_side;
		else first_td_side= first_sieve_side;
		if(poldeg[special_q_side]> 1){
			nr= root_finder(r,poly[special_q_side],poldeg[special_q_side],special_q);
			if(nr==0)continue;
			if(r[nr-1]==special_q){
				
				
				nr--;
			}
		}else{
			u32_t x= mpz_fdiv_ui(poly[special_q_side][1],special_q);
			if(x==0){
				n_spq_discard++;
				continue;
			}
			modulo32= special_q;
			x= modmul32(modinv32(x),mpz_fdiv_ui(poly[special_q_side][0],special_q));
			r[0]= x==0?0:special_q-x;
			nr= 1;
		}
		
		for(root_no= 0;root_no<nr;root_no++){
			u32_t termination_condition;
			
			if(r[root_no]<first_root)continue;
			if((termination_condition= setjmp(termination_jb))!=0){
				if(termination_condition==USER_INTERRUPT)
					/*17:*/
					// #line 440 "gnfs-lasieve4e.w"
					
				{
					char*hn,*ofn;
					FILE*of;
					
					hn= xmalloc(100);
					if(gethostname(hn,99)==0)
						asprintf(&ofn,"%s.%s.last_spq%d",basename,hn,process_no);
					else asprintf(&ofn,"%s.unknown_host.last_spq%d",basename,process_no);
					free(hn);
					
					if((of= fopen(ofn,"w"))!=0){
						fprintf(of,"%u\n",special_q);
						fclose(of);
					}
					free(ofn);
					all_spq_done= 0;
					break;
				}
				
				/*:17*/
				// #line 403 "gnfs-lasieve4e.w"
				
				else{
					char*cmd;
					asprintf(&cmd,"touch badsched.%s.%u.%u.%u",basename,
						special_q_side,special_q,r[root_no]);
					system(cmd);
					free(cmd);
					continue;
				}
			}
			if(sysload_cmd!=NULL){
				
				if(system(sysload_cmd)!=0){
					exitval= LOADTEST_EXITVAL;
					longjmp(termination_jb,USER_INTERRUPT);
				}
			}
			n_spq++;
			reduce2(&a0,&b0,&a1,&b1,(i32_t)special_q,0,(i32_t)r[root_no],1,sigma*sigma);
			/*19:*/
			// #line 478 "gnfs-lasieve4e.w"
			
			{
				if(b0%((i32_t)special_q)==0&&b1%((i32_t)special_q)==0){
					i32_t x;
					
					x= a0%((i32_t)special_q);
					if(x<0)x+= (i32_t)special_q;
					spq_i= x;
					x= a1%((i32_t)special_q);
					if(x<0)x+= (i32_t)special_q;
					spq_j= x;
				}else{
					i32_t x;
					
					x= b0%((i32_t)special_q);
					if(x<0)x+= (i32_t)special_q;
					spq_i= x;
					x= b1%((i32_t)special_q);
					if(x<0)x+= (i32_t)special_q;
					spq_j= x;
				}
				modulo32= special_q;
				spq_x= modmul32(spq_i,i_shift);
			}
			
			/*:19*/
			// #line 422 "gnfs-lasieve4e.w"
			
			//fprintf(ofile,"# Start %u %u (%d,%d) (%d,%d)\n",
			//special_q,r[root_no],a0,b0,a1,b1);
			/*41:*/
			// #line 1401 "gnfs-lasieve4e.w"
			
			{
				u32_t subsieve_nr;
				
				/*42:*/
				// #line 1551 "gnfs-lasieve4e.w"
				
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
					/*60:*/
					// #line 1858 "gnfs-lasieve4e.w"
					
					{
						u32_t s;
						
						for(s= 0;s<2;s++){
							u32_t fbi;
							u16_t*abuf;
							u16_t*ibuf;
							
							abuf= smallsieve_aux[s];
							ibuf= smallpsieve_aux[s];
							for(fbi= 0;fbi<fbis[s];fbi++){
								u32_t aa,bb;
								modulo32= FB[s][fbi];
								
								aa= absa0%FB[s][fbi];
								if(a0s=='-'&&aa!=0)aa= FB[s][fbi]-aa;
								bb= absb0%FB[s][fbi];
								if(proots[s][fbi]!=FB[s][fbi]){
									u32_t x;
									x= modsub32(aa,modmul32(proots[s][fbi],bb));
									if(x!=0){
										aa= absa1%FB[s][fbi];
										if(a1s=='-'&&aa!=0)aa= FB[s][fbi]-aa;
										bb= absb1%FB[s][fbi];
										x= modmul32(asm_modinv32(x),modsub32(modmul32(proots[s][fbi],bb),aa));
										abuf[0]= (u16_t)(FB[s][fbi]);
										abuf[1]= (u16_t)x;
										abuf[2]= (u16_t)(FB_logs[s][fbi]);
										abuf+= 4;
									}else{
										ibuf[0]= (u16_t)(FB[s][fbi]);
										ibuf[1]= (u16_t)(FB_logs[s][fbi]);
										ibuf+= 3;
									}
								}else{
									
									if(bb!=0){
										u32_t x;
										x= modulo32-bb;
										bb= absb1%FB[s][fbi];
										abuf[0]= (u16_t)(FB[s][fbi]);
										abuf[1]= (u16_t)(modmul32(asm_modinv32(x),bb));
										abuf[2]= (u16_t)(FB_logs[s][fbi]);
										abuf+= 4;
									}else{
										ibuf[0]= (u16_t)(FB[s][fbi]);
										ibuf[1]= (u16_t)(FB_logs[s][fbi]);
										ibuf+= 3;
									}
								}
							}
							smallsieve_auxbound[s][0]= abuf;
							smallpsieve_aux_ub_pow1[s]= ibuf;
						}
					}
					
					/*:60*//*61:*/
					// #line 1916 "gnfs-lasieve4e.w"
					
					{
						u32_t s;
						
						for(s= 0;s<2;s++){
							u32_t i;
							u16_t*buf;
							u16_t*buf2;
							u16_t*ibuf;
							
							buf= smallsieve_aux1[s];
							buf2= x2FB[s];
							ibuf= smallpsieve_aux_ub_pow1[s];
							for(i= 0;i<xFBs[s];i++){
								if(xFB[s][i].p==2){
									xFBtranslate(buf2,xFB[s]+i);
									buf2+= 4;
								}else{
									xFBtranslate(buf,xFB[s]+i);
									if(buf[0]==1){
										ibuf[1]= xFB[s][i].l;
										ibuf[0]= xFB[s][i].pp;
										ibuf+= 3;
									}else buf+= 6;
								}
							}
							x2FBs[s]= (buf2-x2FB[s])/4;
							smallpsieve_aux_ub_odd[s]= ibuf;
							smallsieve_aux1_ub_odd[s]= buf;
						}
					}
					
					/*:61*//*62:*/
					// #line 1949 "gnfs-lasieve4e.w"
					
					{
						u32_t s;
						
#ifndef MMX_TD
						for(s= 0;s<2;s++){
							u32_t i;
							u16_t*x;
							
							for(i= 0,x= smallsieve_aux[s];x<smallsieve_auxbound[s][0];i++,x+= 4){
								u32_t k,r,pr;
								
								modulo32= *x;
								r= x[1];
								pr= r;
								for(k= 0;k<j_per_strip;k++){
									smalltdsieve_aux[s][k][i]= r;
									r= modadd32(r,pr);
								}
#ifdef PREINVERT
								/*63:*/
								// #line 1980 "gnfs-lasieve4e.w"
								
								{
									u32_t pinv;
									
									pinv= modulo32;
									pinv= 2*pinv-pinv*pinv*modulo32;
									pinv= 2*pinv-pinv*pinv*modulo32;
									pinv= 2*pinv-pinv*pinv*modulo32;
#if 0
									pinv= 2*pinv-pinv*pinv*modulo32;
#endif
									smalltd_pi[s][i]= 2*pinv-pinv*pinv*modulo32;
								}
								
								/*:63*/
								// #line 1969 "gnfs-lasieve4e.w"
								
#endif
							}
						}
#endif
					}
					
					/*:62*//*64:*/
					// #line 1998 "gnfs-lasieve4e.w"
					
					{
						u32_t s;
						
						for(s= 0;s<2;s++){
							u16_t*x,*xx,k,pbound,copy_buf[6];
							
							k= 0;
							pbound= TINY_SIEVE_MIN;
							for(x= smallsieve_aux[s];x<smallsieve_auxbound[s][0];x+= 4){
								if(*x> pbound){
									if(k==0)smallsieve_tinybound[s]= x;
									else smallsieve_auxbound[s][5-k]= x;
									k++;
									if(k<5)pbound= n_i/(5-k);
									else break;
								}
							}
							while(k<5)smallsieve_auxbound[s][5-(k++)]= x;
							for(x= (xx= smallsieve_aux1[s]);x<smallsieve_aux1_ub_odd[s];x+= 6){
								if(x[0]<TINY_SIEVE_MIN){
									if(x!=xx){
										memcpy(copy_buf,x,6*sizeof(*x));
										memcpy(x,xx,6*sizeof(*x));
										memcpy(xx,copy_buf,6*sizeof(*x));
									}
									xx+= 6;
								}
							}
							smallsieve_tinybound1[s]= xx;
						}
					}
					
					/*:64*/
					// #line 1563 "gnfs-lasieve4e.w"
					
					/*43:*/
					// #line 1572 "gnfs-lasieve4e.w"
					
					{
						u32_t s;
						
						for(s= 0;s<2;s++){
#ifdef SCHEDULING_FUNCTION_CALCULATES_RI
							lasieve_setup(FB[s]+fbis[s],proots[s]+fbis[s],fbi1[s]-fbis[s],
								a0,a1,b0,b1,LPri[s]);
#else
							lasieve_setup(FB[s]+fbis[s],proots[s]+fbis[s],FBsize[s]-fbis[s],
								a0,a1,b0,b1,LPri[s]);
#endif
						}
					}
					
					/*:43*/
					// #line 1564 "gnfs-lasieve4e.w"
					
					/*46:*/
					// #line 1637 "gnfs-lasieve4e.w"
					
					{
						u32_t i,k;
						for(i= 0;i<2;i++){
							double large_primes_summand;
							tpol(tpoly_f[i],poly_f[i],poldeg[i],a0,a1,b0,b1);
							large_primes_summand= sieve_report_multiplier[i]*FB_maxlog[i];
							if(i==special_q_side)
								large_primes_summand+= sieve_multiplier[i]*log(special_q);
							get_sieve_report_bounds(sieve_report_bounds[i],tpoly_f[i],poldeg[i],
								n_srb_i,n_srb_j,2*CANDIDATE_SEARCH_STEPS,
								sieve_multiplier[i],large_primes_summand);
						}
					}
					
					/*:46*/
					// #line 1565 "gnfs-lasieve4e.w"
					
					new_clock= clock();
					sch_clock+= new_clock-last_clock;
					last_clock= new_clock;
}

/*:42*/
// #line 1405 "gnfs-lasieve4e.w"


for(oddness_type= 1;oddness_type<4;oddness_type++){
	/*65:*/
	// #line 2032 "gnfs-lasieve4e.w"
	
	{
		u32_t s;
		for(s= 0;s<2;s++){
			switch(oddness_type){
				u16_t*x;
			case 1:
				/*68:*/
				// #line 2106 "gnfs-lasieve4e.w"
				
				for(x= smallsieve_aux[s];x<smallsieve_auxbound[s][0];x+= 4){
					u32_t p;
					
					p= x[0];
					x[3]= ((i_shift+p)/2)%p;
				}
				
				/*:68*//*71:*/
				// #line 2136 "gnfs-lasieve4e.w"
				
				for(x= smallsieve_aux1[s];x<smallsieve_aux1_ub_odd[s];x+= 6){
					u32_t p;
					
					p= x[0];
					
					x[4]= ((i_shift+p)/2)%p;
					x[5]= 0;
				}
				
				/*:71*//*74:*/
				// #line 2174 "gnfs-lasieve4e.w"
				
				for(x= smallpsieve_aux[s];x<smallpsieve_aux_ub_odd[s];x+= 3)
					x[2]= 0;
				
				/*:74*//*78:*/
				// #line 2212 "gnfs-lasieve4e.w"
				
				{
					u16_t*x,*y,*z;
					u32_t i;
					
					x= smallsieve_aux1_ub_odd[s];
					y= smallpsieve_aux_ub_odd[s];
					z= smallsieve_aux2[s];
					for(i= 0;i<4*x2FBs[s];i+= 4){
						u32_t p,pr,d,l;
						u16_t**a;
						
						d= x2FB[s][i+1];
						if(d==1)continue;
						p= x2FB[s][i];
						pr= x2FB[s][i+2];
						l= x2FB[s][i+3];
						if(p<4){
							if(p==1){
								*y= d/2;
								*(y+2)= 0;
							}else{
								*y= d;
								*(y+2)= d/2;
							}
							*(y+1)= l;
							y+= 3;
							continue;
						}
						p= p/2;
						if(p<=MAX_TINY_2POW)a= &z;
						else a= &x;
						**a= p;
						*(1+*a)= d;
						*(2+*a)= pr%p;
						*(3+*a)= l;
						*(4+*a)= ((i_shift+pr)/2)%p;
						*(5+*a)= d/2;
						*a+= 6;
					}
					smallsieve_aux1_ub[s]= x;
					smallpsieve_aux_ub[s]= y;
					smallsieve_aux2_ub[s]= z;
				}
				
				/*:78*/
				// #line 2039 "gnfs-lasieve4e.w"
				
				break;
			case 2:
				/*69:*/
				// #line 2115 "gnfs-lasieve4e.w"
				
				for(x= smallsieve_aux[s];x<smallsieve_auxbound[s][0];x+= 4){
					u32_t p,pr;
					
					p= x[0];
					pr= x[1];
					x[3]= (pr%2==0?((i_shift+pr)/2)%p:((i_shift+pr+p)/2)%p);
				}
				
				/*:69*//*72:*/
				// #line 2147 "gnfs-lasieve4e.w"
				
				for(x= smallsieve_aux1[s];x<smallsieve_aux1_ub_odd[s];x+= 6){
					u32_t p,d,pr;
					
					p= x[0];
					d= x[1];
					pr= x[2];
					
					x[4]= (pr%2==0?((i_shift+pr)/2)%p:((i_shift+pr+p)/2)%p);
					x[5]= d/2;
				}
				
				/*:72*//*75:*/
				// #line 2179 "gnfs-lasieve4e.w"
				
				for(x= smallpsieve_aux[s];x<smallpsieve_aux_ub_odd[s];x+= 3)
					x[2]= (x[0])/2;
				
				/*:75*//*79:*/
				// #line 2258 "gnfs-lasieve4e.w"
				
				{
					u16_t*x,*y,*z;
					u32_t i;
					
					x= smallsieve_aux1_ub_odd[s];
					y= smallpsieve_aux_ub_odd[s];
					z= smallsieve_aux2[s];
					for(i= 0;i<4*x2FBs[s];i+= 4){
						u32_t p,pr,d,l;
						u16_t**a;
						
						d= x2FB[s][i+1];
						if(d!=1)continue;
						pr= x2FB[s][i+2];
						if(pr%2!=0)continue;
						p= x2FB[s][i];
						l= x2FB[s][i+3];
						if(p<4){
							
							if(p==1){
								Schlendrian("Use 1=2^0 for sieving?\n");
							}
							*y= d;
							*(y+1)= l;
							*(y+2)= 0;
							y+= 3;
							continue;
						}
						p= p/2;
						if(p<=MAX_TINY_2POW)a= &z;
						else a= &x;
						**a= p;
						*(1+*a)= d;
						*(2+*a)= pr%p;
						*(3+*a)= l;
						*(4+*a)= ((i_shift+pr)/2)%p;
						*(5+*a)= 0;
						*a+= 6;
					}
					smallsieve_aux1_ub[s]= x;
					smallpsieve_aux_ub[s]= y;
					smallsieve_aux2_ub[s]= z;
				}
				
				/*:79*/
				// #line 2042 "gnfs-lasieve4e.w"
				
				break;
			case 3:
				/*70:*/
				// #line 2125 "gnfs-lasieve4e.w"
				
				for(x= smallsieve_aux[s];x<smallsieve_auxbound[s][0];x+= 4){
					u32_t p,pr;
					
					p= x[0];
					pr= x[1];
					x[3]= (pr%2==1?((i_shift+pr)/2)%p:((i_shift+pr+p)/2)%p);
				}
				
				/*:70*//*73:*/
				// #line 2160 "gnfs-lasieve4e.w"
				
				for(x= smallsieve_aux1[s];x<smallsieve_aux1_ub_odd[s];x+= 6){
					u32_t p,d,pr;
					
					p= x[0];
					d= x[1];
					pr= x[2];
					
					x[4]= (pr%2==1?((i_shift+pr)/2)%p:((i_shift+pr+p)/2)%p);
					x[5]= d/2;
				}
				
				/*:73*//*76:*/
				// #line 2184 "gnfs-lasieve4e.w"
				
				for(x= smallpsieve_aux[s];x<smallpsieve_aux_ub_odd[s];x+= 3)
					x[2]= (x[0])/2;
				
				/*:76*//*80:*/
				// #line 2304 "gnfs-lasieve4e.w"
				
				{
					u16_t*x,*y,*z;
					u32_t i;
					
					x= smallsieve_aux1_ub_odd[s];
					y= smallpsieve_aux_ub_odd[s];
					z= smallsieve_aux2[s];
					for(i= 0;i<4*x2FBs[s];i+= 4){
						u32_t p,pr,d,l;
						u16_t**a;
						
						d= x2FB[s][i+1];
						if(d!=1)continue;
						pr= x2FB[s][i+2];
						if(pr%2!=1)continue;
						p= x2FB[s][i];
						l= x2FB[s][i+3];
						if(p<4){
							
							if(p==1){
								Schlendrian("Use 1=2^0 for sieving?\n");
							}
							*y= d;
							*(y+1)= l;
							*(y+2)= 0;
							y+= 3;
							continue;
						}
						p= p/2;
						if(p<=MAX_TINY_2POW)a= &z;
						else a= &x;
						**a= p;
						*(1+*a)= d;
						*(2+*a)= pr%p;
						*(3+*a)= l;
						*(4+*a)= ((i_shift+pr)/2)%p;
						*(5+*a)= 0;
						*a+= 6;
					}
					smallsieve_aux1_ub[s]= x;
					smallpsieve_aux_ub[s]= y;
					smallsieve_aux2_ub[s]= z;
				}
				
				/*:80*/
				// #line 2045 "gnfs-lasieve4e.w"
				
				break;
}
}
}


/*:65*//*66:*/
// #line 2054 "gnfs-lasieve4e.w"

#ifdef GCD_SIEVE_BOUND
{
	u32_t i;
	
	for(i= 0;i<np_gcd_sieve;i++){
		gcd_sieve_buffer[2*i+1]= (oddness_type/2)*(gcd_sieve_buffer[2*i]/2);
	}
}
#endif

/*:66*/
// #line 1408 "gnfs-lasieve4e.w"

j_offset= 0;
/*47:*/
// #line 1653 "gnfs-lasieve4e.w"

#ifndef NOSCHED
{
	u32_t s;
	clock_t new_clock;
	
	for(s= 0;s<2;s++){
		u32_t i;
		
		for(i= 0;i<n_schedules[s];i++){
			u32_t ns;
			
			ns= schedules[s][i].n_strips;
			if(ns> n_strips)ns= n_strips;
			do_scheduling(schedules[s]+i,ns,oddness_type,s);
			schedules[s][i].current_strip= 0;
		}
	}
#ifdef GATHER_STAT
	new_clock= clock();
	Schedule_clock+= new_clock-last_clock;
	last_clock= new_clock;
#endif
}
#else 
#define BADSCHED
#endif

/*:47*/
// #line 1410 "gnfs-lasieve4e.w"


#ifdef ZSS_STAT
nss+= n_strips;
#endif
for(subsieve_nr= 0;subsieve_nr<n_strips;
subsieve_nr++,j_offset+= j_per_strip){
	u16_t s,stepno;
#ifdef USE_MEDSCHED
	/*94:*/
	// #line 2759 "gnfs-lasieve4e.w"
	
#ifndef NOSCHED
	for(s= 0;s<2;s++){
		u32_t ll,*sched,*ri;
		
		if(n_medsched_pieces[s]==0)continue;
		for(ll= 0,sched= (u32_t*)med_sched[s][0],ri= LPri[s];
		ll<n_medsched_pieces[s];ll++){
			ri= medsched(ri,current_ij[s]+medsched_fbi_bounds[s][ll],
				current_ij[s]+medsched_fbi_bounds[s][ll+1],&sched,
				medsched_fbi_bounds[s][ll],j_offset==0?oddness_type:0);
			med_sched[s][ll+1]= (u16_t*)sched;
		}
	}
#endif
	
	/*:94*/
	// #line 1419 "gnfs-lasieve4e.w"
	;
	{
		clock_t new_clock;
		new_clock= clock();
		medsched_clock+= new_clock-last_clock;
		last_clock= new_clock;
	}
#endif
	for(s= first_sieve_side,stepno= 0;stepno<2;stepno++,s= 1-s){
		clock_t new_clock,clock_diff;
		
		/*81:*/
		// #line 2350 "gnfs-lasieve4e.w"
		
		{
			u32_t j;
			u16_t*x;
			
			
			for(x= smallsieve_aux[s],j= 0;x<smallsieve_tinybound[s];x+= 4,j++){
				tinysieve_curpos[j]= x[3];
			}
			for(j= 0;j<j_per_strip;j++){
				unsigned char*si_ub;
				bzero(tiny_sieve_buffer,TINY_SIEVEBUFFER_SIZE);
				si_ub= tiny_sieve_buffer+TINY_SIEVEBUFFER_SIZE;
				/*82:*/
				// #line 2372 "gnfs-lasieve4e.w"
				
				{
					u16_t*x;
					
					for(x= smallsieve_aux[s];x<smallsieve_tinybound[s];x+= 4){
						u32_t p,r,pr;
						unsigned char l,*si;
						
						p= x[0];
						pr= x[1];
						l= x[2];
						r= x[3];
						si= tiny_sieve_buffer+r;
						while(si<si_ub){
							*si+= l;
							si+= p;
						}
						r= r+pr;
						if(r>=p)r= r-p;
						x[3]= r;
					}
				}
				
				/*:82*//*83:*/
				// #line 2396 "gnfs-lasieve4e.w"
				
				{
					u16_t*x;
					
					for(x= smallsieve_aux2[s];x<smallsieve_aux2_ub[s];x+= 6){
						u32_t p,r,pr,d,d0;
						unsigned char l,*si;
						
						p= x[0];
						d= x[1];
						pr= x[2];
						l= x[3];
						r= x[4];
						
						d0= x[5];
						if(d0> 0){
							x[5]--;
							continue;
						}
						si= tiny_sieve_buffer+r;
						while(si<si_ub){
							*si+= l;
							si+= p;
						}
						r= r+pr;
						if(r>=p)r= r-p;
						x[4]= r;
						x[5]= d-1;
					}
				}
				
				/*:83*//*84:*/
				// #line 2428 "gnfs-lasieve4e.w"
				
				{
					u16_t*x;
					
					for(x= smallsieve_aux1[s];x<smallsieve_tinybound1[s];x+= 6){
						u32_t p,r,pr,d,d0;
						unsigned char l,*si;
						
						p= x[0];
						d= x[1];
						pr= x[2];
						l= x[3];
						r= x[4];
						
						d0= x[5];
						if(d0> 0){
							x[5]--;
							continue;
						}
						si= tiny_sieve_buffer+r;
						while(si<si_ub){
							*si+= l;
							si+= p;
						}
						r= r+pr;
						if(r>=p)r= r-p;
						x[4]= r;
						x[5]= d-1;
					}
				}
				
				/*:84*/
				// #line 2363 "gnfs-lasieve4e.w"
				
				/*85:*/
				// #line 2460 "gnfs-lasieve4e.w"
				
				{
					unsigned char*si;
					
					si= sieve_interval+(j<<i_bits);
					si_ub= sieve_interval+((j+1)<<i_bits);
					while(si+TINY_SIEVEBUFFER_SIZE<si_ub){
						memcpy(si,tiny_sieve_buffer,TINY_SIEVEBUFFER_SIZE);
						si+= TINY_SIEVEBUFFER_SIZE;
					}
					memcpy(si,tiny_sieve_buffer,si_ub-si);
				}
				
				/*:85*/
				// #line 2364 "gnfs-lasieve4e.w"
				
}
for(x= smallsieve_aux[s],j= 0;x<smallsieve_tinybound[s];x+= 4,j++){
	x[3]= tinysieve_curpos[j];
}
}

/*:81*/
// #line 1430 "gnfs-lasieve4e.w"

#ifdef ZSS_STAT
if(s==1&&ncand==0)
nzss[0]++;
#endif
new_clock= clock();
clock_diff= new_clock-last_clock;
si_clock[s]+= clock_diff;
sieve_clock+= clock_diff;
last_clock= new_clock;
/*86:*/
// #line 2474 "gnfs-lasieve4e.w"

#ifdef ASM_LINESIEVER
slinie(smallsieve_tinybound[s],smallsieve_auxbound[s][4],sieve_interval);
#else
{
	u16_t*x;
	
	for(x= smallsieve_tinybound[s];x<smallsieve_auxbound[s][4];x+= 4){
		u32_t p,r,pr;
		unsigned char l,*y;
		
		p= x[0];
		pr= x[1];
		l= x[2];
		r= x[3];
		for(y= sieve_interval;y<sieve_interval+L1_SIZE;y+= n_i){
			unsigned char*yy,*yy_ub;
			
			yy_ub= y+n_i-3*p;
			for(yy= y+r;yy<yy_ub;yy= yy+4*p){
				*(yy)+= l;
				*(yy+p)+= l;
				*(yy+2*p)+= l;
				*(yy+3*p)+= l;
			}
			while(yy<y+n_i){
				*(yy)+= l;
				yy+= p;
			}
			r= r+pr;
			if(r>=p)r= r-p;
		}
#if 0
		x[3]= r;
#endif
	}
}
#endif

/*:86*//*87:*/
// #line 2514 "gnfs-lasieve4e.w"

#if 1
#ifdef ASM_LINESIEVER3
slinie3(smallsieve_auxbound[s][4],smallsieve_auxbound[s][3],sieve_interval);
#else
{
	u16_t*x;
	
	for(x= smallsieve_auxbound[s][4];x<smallsieve_auxbound[s][3];x+= 4){
		u32_t p,r,pr;
		unsigned char l,*y;
		
		p= x[0];
		pr= x[1];
		l= x[2];
		r= x[3];
		for(y= sieve_interval;y<sieve_interval+L1_SIZE;y+= n_i){
			unsigned char*yy;
			
			yy= y+r;
			*(yy)+= l;
			*(yy+p)+= l;
			*(yy+2*p)+= l;
			yy+= 3*p;
			if(yy<y+n_i)*(yy)+= l;
			r= r+pr;
			if(r>=p)r= r-p;
		}
#if 0
		x[3]= r;
#endif
	}
}
#endif
#endif

/*:87*//*88:*/
// #line 2551 "gnfs-lasieve4e.w"

#if 1
#ifdef ASM_LINESIEVER2
slinie2(smallsieve_auxbound[s][3],smallsieve_auxbound[s][2],sieve_interval);
#else
{
	u16_t*x;
	
	for(x= smallsieve_auxbound[s][3];x<smallsieve_auxbound[s][2];x+= 4){
		u32_t p,r,pr;
		unsigned char l,*y;
		
		p= x[0];
		pr= x[1];
		l= x[2];
		r= x[3];
		for(y= sieve_interval;y<sieve_interval+L1_SIZE;y+= n_i){
			unsigned char*yy;
			
			yy= y+r;
			*(yy)+= l;
			*(yy+p)+= l;
			yy+= 2*p;
			if(yy<y+n_i)*(yy)+= l;
			r= r+pr;
			if(r>=p)r= r-p;
		}
#if 0
		x[3]= r;
#endif
	}
}
#endif
#endif

/*:88*//*89:*/
// #line 2587 "gnfs-lasieve4e.w"

#if 1
#ifdef ASM_LINESIEVER1
slinie1(smallsieve_auxbound[s][2],smallsieve_auxbound[s][1],sieve_interval);
#else
{
	u16_t*x;
	
	for(x= smallsieve_auxbound[s][2];x<smallsieve_auxbound[s][1];x+= 4){
		u32_t p,r,pr;
		unsigned char l,*y;
		
		p= x[0];
		pr= x[1];
		l= x[2];
		r= x[3];
		for(y= sieve_interval;y<sieve_interval+L1_SIZE;y+= n_i){
			unsigned char*yy;
			
			yy= y+r;
			*(yy)+= l;
			yy+= p;
			if(yy<y+n_i)*(yy)+= l;
			r= r+pr;
			if(r>=p)r= r-p;
		}
#if 0
		x[3]= r;
#endif
	}
}
#endif
#endif

/*:89*//*90:*/
// #line 2622 "gnfs-lasieve4e.w"

#if 0
{
	u16_t*x;
	
	for(x= smallsieve_auxbound[s][1];x<smallsieve_auxbound[s][0];x+= 4){
		u32_t p,r,pr;
		unsigned char l,*y;
		
		p= x[0];
		pr= x[1];
		l= x[2];
		r= x[3];
		for(y= sieve_interval;y<sieve_interval+L1_SIZE;y+= n_i){
			if(r<n_i)*(y+r)+= l;
			r= r+pr;
			if(r>=p)r= r-p;
		}
#if 0
		x[3]= r;
#endif
	}
}
#endif

/*:90*//*91:*/
// #line 2648 "gnfs-lasieve4e.w"

#if 1
{
	u16_t*x;
	
	for(x= smallsieve_tinybound1[s];x<smallsieve_aux1_ub[s];x+= 6){
		u32_t p,r,pr,d,d0;
		unsigned char l;
		
		p= x[0];
		d= x[1];
		pr= x[2];
		l= x[3];
		r= x[4];
		
		for(d0= x[5];d0<j_per_strip;d0+= d){
			unsigned char*y,*yy,*yy_ub;
			
			y= sieve_interval+(d0<<i_bits);
			yy_ub= y+n_i-3*p;
			for(yy= y+r;yy<yy_ub;yy= yy+4*p){
				*(yy)+= l;
				*(yy+p)+= l;
				*(yy+2*p)+= l;
				*(yy+3*p)+= l;
			}
			while(yy<y+n_i){
				*(yy)+= l;
				yy+= p;
			}
			r= r+pr;
			if(r>=p)r= r-p;
		}
		x[4]= r;
		x[5]= d0-j_per_strip;
	}
}
#endif

/*:91*//*92:*/
// #line 2689 "gnfs-lasieve4e.w"

#if 1
{
	u16_t*x;
	
	bzero(horizontal_sievesums,j_per_strip);
	for(x= smallpsieve_aux[s];x<smallpsieve_aux_ub[s];x+= 3){
		u32_t p,d;
		unsigned char l;
		
		p= x[0];
		l= x[1];
		d= x[2];
		while(d<j_per_strip){
			horizontal_sievesums[d]+= l;
			d+= p;
		}
#if 0
		x[2]= d-j_per_strip;
#endif
	}
}
#else
bzero(horizontal_sievesums,j_per_strip);
#endif

/*:92*/
// #line 1440 "gnfs-lasieve4e.w"

new_clock= clock();
clock_diff= new_clock-last_clock;
s1_clock[s]+= clock_diff;
sieve_clock+= clock_diff;
last_clock= new_clock;
#ifdef BADSCHED
ncand= 0;
continue;
#endif
/*93:*/
// #line 2720 "gnfs-lasieve4e.w"

#ifndef MEDSCHE_SI_OFFS
#ifdef BIGENDIAN
#define MEDSCHED_SI_OFFS 1
#else
#define MEDSCHED_SI_OFFS 0
#endif
#endif
#ifdef ASM_SCHEDSIEVE1
schedsieve(medsched_logs[s],n_medsched_pieces[s],med_sched[s],sieve_interval);
#else
{
	u32_t l;
	
	for(l= 0;l<n_medsched_pieces[s];l++){
		unsigned char x;
		u16_t*schedule_ptr;
		
		x= medsched_logs[s][l];
#ifdef ASM_SCHEDSIEVE
		schedsieve(x,sieve_interval,med_sched[s][l],med_sched[s][l+1]);
#else
		for(schedule_ptr= med_sched[s][l]+MEDSCHED_SI_OFFS;
		schedule_ptr+3*SE_SIZE<med_sched[s][l+1];
		schedule_ptr+= 4*SE_SIZE){
			sieve_interval[*schedule_ptr]+= x;
			sieve_interval[*(schedule_ptr+SE_SIZE)]+= x;
			sieve_interval[*(schedule_ptr+2*SE_SIZE)]+= x;
			sieve_interval[*(schedule_ptr+3*SE_SIZE)]+= x;
		}
		for(;
		schedule_ptr<med_sched[s][l+1];schedule_ptr+= SE_SIZE)
			sieve_interval[*schedule_ptr]+= x;
#endif
	}
}
#endif

/*:93*/
// #line 1450 "gnfs-lasieve4e.w"

new_clock= clock();
clock_diff= new_clock-last_clock;
s2_clock[s]+= clock_diff;
sieve_clock+= clock_diff;
last_clock= new_clock;
/*97:*/
// #line 2797 "gnfs-lasieve4e.w"

#ifndef SCHED_SI_OFFS
#ifdef BIGENDIAN
#define SCHED_SI_OFFS 1
#else
#define SCHED_SI_OFFS 0
#endif
#endif
{
	u32_t j;
	
	for(j= 0;j<n_schedules[s];j++){
		if(schedules[s][j].current_strip==schedules[s][j].n_strips){
			u32_t ns;
			
			ns= schedules[s][j].n_strips;
			if(ns> n_strips-subsieve_nr)ns= n_strips-subsieve_nr;
			do_scheduling(schedules[s]+j,ns,0,s);
			schedules[s][j].current_strip= 0;
		}
	}
#ifdef GATHER_STAT
	new_clock= clock();
	Schedule_clock+= new_clock-last_clock;
	last_clock= new_clock;
#endif
	
	for(j= 0;j<n_schedules[s];j++){
#ifdef ASM_SCHEDSIEVE1
		u32_t i,k;
		
		k= schedules[s][j].current_strip;
		for(i= 0;i<=schedules[s][j].n_pieces;i++){
			schedbuf[i]= schedules[s][j].schedule[i][k];
		}
		schedsieve(schedules[s][j].schedlogs,schedules[s][j].n_pieces,
			schedbuf,sieve_interval);
#else
		u32_t l,k;
		
		k= schedules[s][j].current_strip;
		l= 0;
		while(l<schedules[s][j].n_pieces){
			unsigned char x;
			u16_t*schedule_ptr,*sptr_ub;
			
			x= schedules[s][j].schedlogs[l];
			schedule_ptr= schedules[s][j].schedule[l][k]+SCHED_SI_OFFS;
			while(l<schedules[s][j].n_pieces)
				if(schedules[s][j].schedlogs[++l]!=x)break;
				sptr_ub= schedules[s][j].schedule[l][k];
				
#ifdef ASM_SCHEDSIEVE
				schedsieve(x,sieve_interval,schedule_ptr,sptr_ub);
#else
				while(schedule_ptr+3*SE_SIZE<sptr_ub){
					sieve_interval[*schedule_ptr]+= x;
					sieve_interval[*(schedule_ptr+SE_SIZE)]+= x;
					sieve_interval[*(schedule_ptr+2*SE_SIZE)]+= x;
					sieve_interval[*(schedule_ptr+3*SE_SIZE)]+= x;
					schedule_ptr+= 4*SE_SIZE;
				}
				while(schedule_ptr<sptr_ub){
					sieve_interval[*schedule_ptr]+= x;
					schedule_ptr+= SE_SIZE;
				}
#endif
		}
#endif
	}
}

// #line 1 "la-cs.w"
/*:97*/
// #line 1456 "gnfs-lasieve4e.w"

#if 0
dumpsieve(j_offset,s);
#endif
new_clock= clock();
clock_diff= new_clock-last_clock;
sieve_clock+= clock_diff;
s3_clock[s]+= clock_diff;
last_clock= new_clock;

if(s==first_sieve_side){
#ifdef GCD_SIEVE_BOUND
	gcd_sieve();
#endif
	/*98:*/
	// #line 11 "la-cs.w"
	
#ifdef ASM_SEARCH0
	{
		unsigned char*srbs;
		u32_t i;
		srbs= sieve_report_bounds[s][j_offset/CANDIDATE_SEARCH_STEPS];
		ncand= lasieve_search0(sieve_interval,horizontal_sievesums,
			horizontal_sievesums+j_per_strip,
			srbs,srbs+n_i/CANDIDATE_SEARCH_STEPS,
			cand,fss_sv);
		for(i= 0;i<ncand;i++)fss_sv[i]+= horizontal_sievesums[cand[i]>>i_bits];
	}
#else
	{
		unsigned char*srbs;
		u32_t i;
		
		srbs= sieve_report_bounds[s][j_offset/CANDIDATE_SEARCH_STEPS];
		ncand= 0;
		for(i= 0;i<n_i;i+= CANDIDATE_SEARCH_STEPS){
			unsigned char st;
			u32_t j;
			
			st= *(srbs++);
			for(j= 0;j<j_per_strip;j++){
				unsigned char*i_o,*i_max,st1;
				
				i_o= sieve_interval+(j<<i_bits)+i;
				i_max= i_o+CANDIDATE_SEARCH_STEPS;
				if(st<=horizontal_sievesums[j]){
					while(i_o<i_max){
						cand[ncand]= i_o-sieve_interval;
						fss_sv[ncand++]= *(i_o++)+horizontal_sievesums[j];
					}
					continue;
				}
				st1= st-horizontal_sievesums[j];
				/*99:*/
				// #line 71 "la-cs.w"
				
#ifndef HAVE_SSIMD
#ifdef GNFS_CS32 
#define bc_t unsigned long
#define BC_MASK 0x80808080
#else
#define bc_t unsigned long long
#define BC_MASK 0x8080808080808080
#endif
				{
					if(st1<0x80){
						bc_t bc,*i_oo;
						
						bc= st1;
						bc= (bc<<8)|bc;
						bc= (bc<<16)|bc;
#ifndef GNFS_CS32
						bc= (bc<<32)|bc;
#endif
						bc= BC_MASK-bc;
						for(i_oo= (bc_t*)i_o;i_oo<(bc_t*)i_max;i_oo++){
							bc_t v= *i_oo;
							if(((v&BC_MASK)|((v+bc)&BC_MASK))==0)continue;
							for(i_o= (unsigned char*)i_oo;i_o<(unsigned char*)(i_oo+1);i_o++){
								if(*i_o>=st1){
									/*100:*/
									// #line 154 "la-cs.w"
									
									cand[ncand]= i_o-sieve_interval;
									fss_sv[ncand++]= *i_o+horizontal_sievesums[j];
									
									/*:100*/
									// #line 96 "la-cs.w"
									
								}
							}
						}
					}else{
						bc_t*i_oo;
						
						for(i_oo= (bc_t*)i_o;i_oo<(bc_t*)i_max;i_oo++){
							if((*i_oo&BC_MASK)==0)continue;
							for(i_o= (unsigned char*)i_oo;i_o<(unsigned char*)(i_oo+1);i_o++){
								if(*i_o>=st1){
									/*100:*/
									// #line 154 "la-cs.w"
									
									cand[ncand]= i_o-sieve_interval;
									fss_sv[ncand++]= *i_o+horizontal_sievesums[j];
									
									/*:100*/
									// #line 107 "la-cs.w"
									
								}
							}
						}
					}
				}
#else
				{
					unsigned long long x;
					
					x= st1-1;
					x|= x<<8;
					x|= x<<16;
					x|= x<<32;
					while(i_o<i_max){
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
							"emms":"=S"(i_o):"a"(&x),"S"(i_o),
							"D"(i_max));
						if(i_o<i_max){
							unsigned char*i_max2= i_o+32;
							while(i_o<i_max2){
								if(*i_o>=st1){
									/*100:*/
									// #line 154 "la-cs.w"
									
									cand[ncand]= i_o-sieve_interval;
									fss_sv[ncand++]= *i_o+horizontal_sievesums[j];
									
									/*:100*/
									// #line 144 "la-cs.w"
									
								}
								i_o++;
							}
						}
					}
				}
#endif
				
				/*:99*/
				// #line 48 "la-cs.w"
				
}
}
}
#endif
#if 0
{
	char*ofn;
	FILE*of;
	asprintf(&ofn,"cdump.%u.%u.j%u.ot%u",special_q,r[root_no],
		j_offset,oddness_type);
	if((of= fopen(ofn,"w"))!=NULL){
		u32_t i;
		fprintf(of,"%u candidates\n",ncand);
		for(i= 0;i<ncand;i++)
			fprintf(of,"%u %u\n",cand[i],fss_sv[i]);
		fclose(of);
	}else errprintf("Cannot open debug file %s: %m\n",ofn);
	free(ofn);
}
#endif

/*:98*/
// #line 1470 "gnfs-lasieve4e.w"

}
else
/*101:*/
// #line 159 "la-cs.w"

{
	u32_t i,nc1;
	unsigned char*srbs;
	static u32_t bad_pvl= 0;
	
	srbs= sieve_report_bounds[s][j_offset/CANDIDATE_SEARCH_STEPS];
	n_prereports+= ncand;
	for(i= 0,nc1= 0;i<ncand;i++){
		u16_t st_i,t_j,ii,jj,j;
		double pvl;
		
		j= cand[i]>>i_bits;
#ifndef DEBUG_SIEVE_REPORT_BOUNDS
		if(sieve_interval[cand[i]]+horizontal_sievesums[j]<
			srbs[(cand[i]&(n_i-1))/CANDIDATE_SEARCH_STEPS])
			continue;
#endif
		jj= j_offset+j;
		ii= cand[i]&(n_i-1);
		st_i= 2*ii+(oddness_type==2?0:1);
		t_j= 2*jj+(oddness_type==1?0:1);
		pvl= log(fabs(rpol_eval(tpoly_f[s],poldeg[s],
			(double)st_i-(double)i_shift,(double)t_j)));
		if(special_q_side==s)
			pvl-= special_q_log;
		pvl*= sieve_multiplier[s];
		pvl-= sieve_report_multiplier[s]*FB_maxlog[s];
		if((double)(sieve_interval[cand[i]]+horizontal_sievesums[j])>=pvl){
#ifdef DEBUG_SIEVE_REPORT_BOUNDS
			/*102:*/
			// #line 199 "la-cs.w"
			
			if(sieve_interval[cand[i]]+horizontal_sievesums[j]<
				srbs[(cand[i]&(n_i-1))/CANDIDATE_SEARCH_STEPS]){
				double pvl1;
				
				pvl= fabs(rpol_eval(tpoly_f[s],poldeg[s],
					(double)st_i-(double)i_shift,(double)t_j));
				fprintf(stderr,"Bad pvl min %u at (%f,%f),spq=%u\npvl: %.5g->",
					bad_pvl++,(double)st_i-(double)i_shift,(double)t_j,
					special_q,pvl);
				pvl= log(pvl);
				fprintf(stderr,"%.3f->",pvl);
				pvl= sieve_multiplier[s]*pvl;
				fprintf(stderr,"%.3f->",pvl);
				if(special_q_side==s)pvl-= sieve_multiplier[s]*special_q_log;
				fprintf(stderr,"%.3f->",pvl);
				pvl-= sieve_report_multiplier[s]*FB_maxlog[s];
				fprintf(stderr,"%.3f\nLower bound was %u sv was %u=%u+%u\n",pvl,
					(u32_t)srbs[(cand[i]&(n_i-1))/CANDIDATE_SEARCH_STEPS],
					(u32_t)sieve_interval[cand[i]]+(u32_t)horizontal_sievesums[j],
					(u32_t)sieve_interval[cand[i]],(u32_t)horizontal_sievesums[j]);
			}
			// #line 2870 "gnfs-lasieve4e.w"
			
			/*:102*/
			// #line 189 "la-cs.w"
			
#endif
			fss_sv[nc1]= fss_sv[i];
			cand[nc1++]= cand[i];
		}
	}
	ncand= nc1;
}

/*:101*/
// #line 1473 "gnfs-lasieve4e.w"

new_clock= clock();
clock_diff= new_clock-last_clock;
sieve_clock+= clock_diff;
cs_clock[s]+= clock_diff;
last_clock= new_clock;
}
#ifndef BADSCHED
trial_divide();
#endif
{
	clock_t new_clock;
	new_clock= clock();
	td_clock+= new_clock-last_clock;
	last_clock= new_clock;
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

/*:41*/
// #line 425 "gnfs-lasieve4e.w"

//fprintf(ofile,"# Done %u %u (%d,%d) (%d,%d)\n",
//special_q,r[root_no],a0,b0,a1,b1);
}
if(root_no<nr){
	
	
	
	break;
}

tNow = sTime();
if (tNow > lastReport + 5.0) {
	lastReport = sTime();
	fprintf(stderr, "\rtotal yield: %u, q=%u (%1.5lf sec/rel)", 
		(unsigned int)yield, (unsigned int)special_q, (tNow - tStart)/yield);
	fflush(stderr);
	
}}
fprintf(stderr, "\rtotal yield: %u, q=%u (%1.5lf sec/rel)\n", 
		(unsigned int)yield, (unsigned int)special_q, ((clock() / (double)CLOCKS_PER_SEC) - tStart)/yield);
free(r);
}

/*:16*/
// #line 310 "gnfs-lasieve4e.w"

if(sieve_count!=0){
	if(zip_output!=0)pclose(ofile);
	else fclose(ofile);
}
logbook(0,"%u Special q, %u reduction iterations\n",n_spq,n_iter);

if(n_spq_discard> 0)logbook(0,"%u Special q discarded\n",n_spq_discard);
/*20:*/
// #line 504 "gnfs-lasieve4e.w"

{
	u32_t side;
	logbook(0,"reports: %u->%u->%u->%u->%u->%u\n",
		n_prereports,n_reports,n_rep1,n_rep2,
		n_tdsurvivors[first_td_side],n_tdsurvivors[1-first_td_side]);
	logbook(0,"Number of relations with k rational and l algebraic primes for (k,l)=:\n");
	
	sieve_clock= rint((1000.0*sieve_clock)/CLOCKS_PER_SEC);
	sch_clock= rint((1000.0*sch_clock)/CLOCKS_PER_SEC);
	td_clock= rint((1000.0*td_clock)/CLOCKS_PER_SEC);
	tdi_clock= rint((1000.0*tdi_clock)/CLOCKS_PER_SEC);
	Schedule_clock= rint((1000.0*Schedule_clock)/CLOCKS_PER_SEC);
	medsched_clock= rint((1000.0*medsched_clock)/CLOCKS_PER_SEC);
	mpqs_clock= rint((1000.0*mpqs_clock)/CLOCKS_PER_SEC);
	
	for(side= 0;side<2;side++){
		cs_clock[side]= rint((1000.0*cs_clock[side])/CLOCKS_PER_SEC);
		si_clock[side]= rint((1000.0*si_clock[side])/CLOCKS_PER_SEC);
		s1_clock[side]= rint((1000.0*s1_clock[side])/CLOCKS_PER_SEC);
		s2_clock[side]= rint((1000.0*s2_clock[side])/CLOCKS_PER_SEC);
		s3_clock[side]= rint((1000.0*s3_clock[side])/CLOCKS_PER_SEC);
		tdsi_clock[side]= rint((1000.0*tdsi_clock[side])/CLOCKS_PER_SEC);
		tds1_clock[side]= rint((1000.0*tds1_clock[side])/CLOCKS_PER_SEC);
		tds2_clock[side]= rint((1000.0*tds2_clock[side])/CLOCKS_PER_SEC);
		tds3_clock[side]= rint((1000.0*tds3_clock[side])/CLOCKS_PER_SEC);
		tds4_clock[side]= rint((1000.0*tds4_clock[side])/CLOCKS_PER_SEC);
	}
	
	logbook(0,"\nTotal yield: %u\n",yield);
	if(n_mpqsfail[0]!=0||n_mpqsfail[1]!=0||
		n_mpqsvain[0]!=0||n_mpqsvain[1]!=0){
		logbook(0,"%u/%u mpqs failures, %u/%u vain mpqs\n",n_mpqsfail[0],
			n_mpqsfail[1],n_mpqsvain[0],n_mpqsvain[1]);
	}
	logbook(0,"milliseconds total: Sieve %d Sched %d medsched %d\n",
		(int)sieve_clock,(int)Schedule_clock,(int)medsched_clock);
	logbook(0,"TD %d (Init %d, MPQS %d) Sieve-Change %d\n",
		(int)td_clock,(int)tdi_clock,(int)mpqs_clock,(int)sch_clock);
	for(side= 0;side<2;side++){
		logbook(0,"TD side %d: init/small/medium/large/search: %d %d %d %d %d\n",
			(int)side,(int)tdsi_clock[side],(int)tds1_clock[side],
			(int)tds2_clock[side],(int)tds3_clock[side],(int)tds4_clock[side]);
		logbook(0,"sieve: init/small/medium/large/search: %d %d %d %d %d\n",
			(int)si_clock[side],(int)s1_clock[side],(int)s2_clock[side],
			(int)s3_clock[side],(int)cs_clock[side]);
	}
#ifdef MMX_TDBENCH
	fprintf(stderr,"MMX-Loops: %qu\n",MMX_TdNloop);
#endif
#ifdef ZSS_STAT
	fprintf(stderr,
		"%u subsieves, zero: %u first sieve, %u second sieve %u first td\n",
		nss,nzss[0],nzss[1],nzss[2]);
#endif
}

/*:20*/
// #line 318 "gnfs-lasieve4e.w"

if(special_q>=last_spq&&all_spq_done!=0)exit(0);
if(exitval==0)exitval= 1;
exit(exitval);
}

/*:15*//*49:*/
// #line 1686 "gnfs-lasieve4e.w"

#ifndef NOSCHED
void do_scheduling(struct schedule_struct*sched,u32_t ns,u32_t ot,u32_t s)
{
	u32_t ll,n1_j,*ri;;
	
	n1_j= ns<<(L1_BITS-i_bits);
	for(ll= 0,ri= sched->ri;ll<sched->n_pieces;ll++){
		u32_t fbi_lb,fbi_ub,fbio;
		memcpy(sched->schedule[ll+1],sched->schedule[ll],ns*sizeof(u16_t**));
		fbio= sched->fbi_bounds[ll];
		fbi_lb= fbio;
		fbi_ub= sched->fbi_bounds[ll+1];
#ifdef SCHEDULING_FUNCTION_CALCULATES_RI
		if(ot==1)
			lasieve_setup(FB[s]+fbi_lb,proots[s]+fbi_lb,fbi_ub-fbi_lb,
			a0,a1,b0,b1,LPri[s]+(fbi_lb-fbis[s])*RI_SIZE);
#endif
		ri= lasched(ri,current_ij[s]+fbi_lb,current_ij[s]+fbi_ub,
			n1_j,(u32_t**)(sched->schedule[ll+1]),fbi_lb-fbio,ot);
		/*50:*/
		// #line 1712 "gnfs-lasieve4e.w"
		
		{
			u32_t k;
			for(k= 0;k<ns;k++)
				if(sched->schedule[ll+1][k]>=sched->schedule[0][k]+sched->alloc){
					if(k==0&&sched->schedule[ll+1][k]<sched->schedule[0][k]+sched->alloc1)
						continue;
					longjmp(termination_jb,SCHED_PATHOLOGY);
				}
		}
		
		/*:50*/
		// #line 1706 "gnfs-lasieve4e.w"
		
	}
}
#endif

/*:49*//*67:*/
// #line 2066 "gnfs-lasieve4e.w"

#ifdef GCD_SIEVE_BOUND
static void
gcd_sieve()
{
	u32_t i;
	
	for(i= 0;i<np_gcd_sieve;i++){
		u32_t x,p;
		
		x= gcd_sieve_buffer[2*i+1];
		p= gcd_sieve_buffer[2*i];
		while(x<j_per_strip){
			unsigned char*z,*z_ub;
			
			z= sieve_interval+(x<<i_bits);
			z_ub= z+n_i-3*p;
			z+= oddness_type==2?(n_i/2)%p:((n_i+p-1)/2)%p;
			while(z<z_ub){
				*z= 0;
				*(z+p)= 0;
				z+= 2*p;
				*z= 0;
				*(z+p)= 0;
				z+= 2*p;
			}
			z_ub+= 3*p;
			while(z<z_ub){
				*z= 0;
				z+= p;
			}
			x= x+p;
		}
		gcd_sieve_buffer[2*i+1]= x-j_per_strip;
	}
}
#endif

/*:67*//*104:*/
// #line 2876 "gnfs-lasieve4e.w"

static void
xFBtranslate(u16_t*rop,xFBptr op)
{
	u32_t x,y,am,bm,rqq;
	
	modulo32= op->pp;
	rop[3]= op->l;
	am= a1> 0?((u32_t)a1)%modulo32:modulo32-((u32_t)(-a1))%modulo32;
	if(am==modulo32)am= 0;
	bm= b1> 0?((u32_t)b1)%modulo32:modulo32-((u32_t)(-b1))%modulo32;
	if(bm==modulo32)bm= 0;
	x= modsub32(modmul32(op->qq,am),modmul32(op->r,bm));
	am= a0> 0?((u32_t)a0)%modulo32:modulo32-((u32_t)(-a0))%modulo32;
	if(am==modulo32)am= 0;
	bm= b0> 0?((u32_t)b0)%modulo32:modulo32-((u32_t)(-b0))%modulo32;
	if(bm==modulo32)bm= 0;
y= modsub32(modmul32(op->r,bm),modmul32(op->qq,am));
rqq= 1;
if(y!=0){
while(y%(op->p)==0){
y= y/(op->p);
rqq*= op->p;
}
}else{
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

/*:104*//*105:*/
// #line 2913 "gnfs-lasieve4e.w"

static int
xFBcmp(const void*opA,const void*opB)
{
xFBptr op1,op2;
op1= (xFBptr)opA;
op2= (xFBptr)opB;
if(op1->pp<op2->pp)return-1;
if(op1->pp==op2->pp)return 0;
return 1;
}

/*:105*//*106:*/
// #line 2926 "gnfs-lasieve4e.w"

static u32_t
add_primepowers2xaFB(size_t*xaFB_alloc_ptr,u32_t pp_bound,
					 u32_t s,u32_t p,u32_t r)
{
u32_t a,b,q,qo,*rbuf,nr,*Ar,exponent,init_xFB;
size_t rbuf_alloc;
if(xFBs[s]==0&&p==0)Schlendrian("add_primepowers2xaFB on empty xaFB\n");

rbuf_alloc= 0;
Ar= xmalloc((1+poldeg[s])*sizeof(*Ar));

if(p!=0){
init_xFB= 0;
q= p;
if(r==p){
a= 1;
b= p;
}else{
a= r;
b= 1;
}
}else{
init_xFB= 1;
q= xFB[s][xFBs[s]-1].pp;
p= xFB[s][xFBs[s]-1].p;
a= xFB[s][xFBs[s]-1].r;
b= xFB[s][xFBs[s]-1].qq;
}

qo= q;
exponent= 1;
for(;;){
u32_t j,r;
if(q> pp_bound/p)break;
modulo32= p*q;
for(j= 0;j<=poldeg[s];j++)
Ar[j]= mpz_fdiv_ui(poly[s][j],modulo32);
if(b==1)/*107:*/
// #line 2987 "gnfs-lasieve4e.w"

{
for(r= a,nr= 0;r<modulo32;r+= qo){
u32_t pv;
for(j= 1,pv= Ar[poldeg[s]];j<=poldeg[s];j++){
pv= modadd32(Ar[poldeg[s]-j],modmul32(pv,r));
}
if(pv==0){
adjust_bufsize((void**)&rbuf,&rbuf_alloc,1+nr,4,sizeof(*rbuf));
rbuf[nr++]= r;
}else if(pv%q!=0)Schlendrian("xFBgen: %u not a root mod %u\n",
r,q);
}
}

/*:107*/
// #line 2964 "gnfs-lasieve4e.w"

else/*108:*/
// #line 3003 "gnfs-lasieve4e.w"

{
for(r= (modmul32(b,modinv32(a)))%qo,nr= 0;r<modulo32;r+= qo){
u32_t pv;
for(j= 1,pv= Ar[0];j<=poldeg[s];j++){
pv= modadd32(Ar[j],modmul32(pv,r));
}
if(pv==0){
adjust_bufsize((void**)&rbuf,&rbuf_alloc,1+nr,4,sizeof(*rbuf));
rbuf[nr++]= r;
}else if(pv%q!=0)Schlendrian("xFBgen: %u^{-1} not a root mod %u\n",
r,q);
}
}

/*:108*/
// #line 2965 "gnfs-lasieve4e.w"

if(qo*nr!=modulo32)break;
q= modulo32;
exponent++;
}
if(init_xFB!=0)
xFB[s][xFBs[s]-1].l= 
rint(sieve_multiplier[s]*log(q))-rint(sieve_multiplier[s]*log(qo/p));
if(q<=pp_bound/p){
u32_t j;
for(j= 0;j<nr;j++){
/*109:*/
// #line 3019 "gnfs-lasieve4e.w"

xFBptr f;

adjust_bufsize((void**)&(xFB[s]),xaFB_alloc_ptr,1+xFBs[s],16,sizeof(**xFB));
f= xFB[s]+xFBs[s];
f->p= p;
f->pp= q*p;
if(b==1){
f->qq= 1;
f->r= rbuf[j];
f->q= f->pp;
}else{
modulo32= (q*p)/b;
rbuf[j]= rbuf[j]/b;
if(rbuf[j]==0){
f->qq= f->pp;
f->q= 1;
f->r= 1;
}else{
while(rbuf[j]%p==0){
rbuf[j]= rbuf[j]/p;
modulo32= modulo32/p;
}
f->qq= (f->pp)/modulo32;
f->q= modulo32;
f->r= modinv32(rbuf[j]);
}
}

/*:109*/
// #line 2976 "gnfs-lasieve4e.w"

xFBs[s]++;
add_primepowers2xaFB(xaFB_alloc_ptr,pp_bound,s,0,0);
}
}
if(rbuf_alloc> 0)free(rbuf);
free(Ar);
return exponent;
}

/*:106*//*111:*/
// #line 3054 "gnfs-lasieve4e.w"

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

/*112:*/
// #line 3087 "gnfs-lasieve4e.w"

for(ci= 0,nc1= 0;ci<ncand;ci++){
u16_t strip_i,strip_j;
u16_t st_i,true_j;
u16_t s;
double pvl;

/*113:*/
// #line 3116 "gnfs-lasieve4e.w"

{
u16_t jj;

strip_j= cand[ci]>>i_bits;
jj= j_offset+strip_j;
strip_i= cand[ci]&(n_i-1);
st_i= 2*strip_i+(oddness_type==2?0:1);
true_j= 2*jj+(oddness_type==1?0:1);
}

/*:113*/
// #line 3094 "gnfs-lasieve4e.w"

n_reports++;
s= first_sieve_side;
#ifdef STC_DEBUG
fprintf(debugfile,"%hu %hu\n",st_i,true_j);
#endif
if(gcd32(st_i<i_shift?i_shift-st_i:st_i-i_shift,true_j)!=1)continue;
n_rep1++;
pvl= log(fabs(rpol_eval(tpoly_f[s],poldeg[s],
(double)st_i-(double)i_shift,(double)true_j)));
if(special_q_side==s)pvl-= special_q_log;
pvl*= sieve_multiplier[s];
if((double)fss_sv[ci]+sieve_report_multiplier[s]*FB_maxlog[s]<pvl)continue;
n_rep2++;


modulo32= special_q;
if(modadd32(modmul32(st_i,spq_i),modmul32(true_j,spq_j))==spq_x)continue;
cand[nc1++]= cand[ci];
}

/*:112*/
// #line 3067 "gnfs-lasieve4e.w"

#ifdef ZSS_STAT
if(ncand==0)
nzss[1]++;
#endif
last_tdclock= clock();
tdi_clock+= last_tdclock-last_clock;
ncand= nc1;
qsort(cand,ncand,sizeof(*cand),tdcand_cmp);
td_buf1[0]= td_buf[first_td_side];
for(side= first_td_side,tdstep= 0;tdstep<2;side= 1-side,tdstep++){
#ifdef ZSS_STAT
if(tdstep==1&&ncand==0)
nzss[2]++;
#endif
/*116:*/
// #line 3151 "gnfs-lasieve4e.w"

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

/*117:*/
// #line 3246 "gnfs-lasieve4e.w"

{
unsigned char ht,allcoll;

bzero(sieve_interval,L1_SIZE);
bzero(tds_coll,UCHAR_MAX-1);
for(ci= 0,ht= 1,allcoll= 0;ci<ncand;ci++){
unsigned char cht;

cht= sieve_interval[cand[ci]];
if(cht==0){
cht= ht;
if(ht<UCHAR_MAX)ht++;
else{
ht= 1;
allcoll= 1;
}
tds_coll[cht-1]= allcoll;
sieve_interval[cand[ci]]= cht;
}else{
tds_coll[cht-1]= 1;
}
fss_sv[ci]= cht-1;
}
}

/*:117*//*118:*/
// #line 3280 "gnfs-lasieve4e.w"

#ifdef MMX_TD
smalltdsieve_auxbound= MMX_TdInit(side,smallsieve_aux[side],
smallsieve_auxbound[side][0],
&p_bound,j_offset==0&&oddness_type==1);
#else
{
u16_t*x,*z;

x= smallsieve_aux[side];
z= smallsieve_auxbound[side][0];
if(*x> p_bound)smalltdsieve_auxbound= x;
else{
while(x+4<z){
u16_t*y;

y= x+4*((z-x)/8);
if(y==smallsieve_auxbound[side][0]||*y> p_bound)z= y;
else x= y;
}
smalltdsieve_auxbound= z;
}
}
#endif

/*:118*/
// #line 3173 "gnfs-lasieve4e.w"

newclock= clock();
tdsi_clock[side]+= newclock-last_tdclock;
last_tdclock= newclock;
/*121:*/
// #line 3328 "gnfs-lasieve4e.w"

memcpy(tds_fbi_curpos,tds_fbi,UCHAR_MAX*sizeof(*tds_fbi));

/*:121*//*122:*/
// #line 3332 "gnfs-lasieve4e.w"

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

for(l= 0;l<n_medsched_pieces[side];l++){
u16_t*x,*x_ub;

x_ub= med_sched[side][l+1];

for(x= med_sched[side][l]+MEDSCHED_SI_OFFS;x+6<x_ub;x+= 8){
unsigned char z;
if((sieve_interval[*x]|sieve_interval[*(x+2)]|
sieve_interval[*(x+4)]|sieve_interval[*(x+6)])==0){
continue;
}
if((z= sieve_interval[*x])!=0)
*(tds_fbi_curpos[z-1]++)= *(x+1-2*MEDSCHED_SI_OFFS);
if((z= sieve_interval[*(x+2)])!=0)
*(tds_fbi_curpos[z-1]++)= *(x+3-2*MEDSCHED_SI_OFFS);
if((z= sieve_interval[*(x+4)])!=0)
*(tds_fbi_curpos[z-1]++)= *(x+5-2*MEDSCHED_SI_OFFS);
if((z= sieve_interval[*(x+6)])!=0)
*(tds_fbi_curpos[z-1]++)= *(x+7-2*MEDSCHED_SI_OFFS);
}
while(x<x_ub){
unsigned char z;

if((z= sieve_interval[*x])!=0)
*(tds_fbi_curpos[z-1]++)= *(x+1-2*MEDSCHED_SI_OFFS);
x+= 2;
}
}
}
#endif
newclock= clock();
tds2_clock[side]+= newclock-last_tdclock;
last_tdclock= newclock;

/*:122*//*123:*/
// #line 3381 "gnfs-lasieve4e.w"

{
u32_t j;

for(j= 0;j<n_schedules[side];j++){
#ifdef ASM_SCHEDTDSIEVE
u32_t i,k;
k= schedules[side][j].current_strip++;
for(i= 0;i<=schedules[side][j].n_pieces;i++){
schedbuf[i]= schedules[side][j].schedule[i][k];
}
schedtdsieve(schedules[side][j].fbi_bounds,schedules[side][j].n_pieces,
schedbuf,sieve_interval,tds_fbi_curpos);
#else
#if 1
u32_t k,l,fbi_offset;
u16_t*x,*x_ub;
k= schedules[side][j].current_strip++;
x= schedules[side][j].schedule[0][k]+SCHED_SI_OFFS;
x_ub= schedules[side][j].schedule[schedules[side][j].n_pieces][k];
l= 0;
fbi_offset= schedules[side][j].fbi_bounds[l];
while(x<x_ub){
u16_t**b0,**b1,**b0_ub;
#ifdef ASM_SCHEDTDSIEVE2
b0= tdsieve_sched2buf(&x,x_ub,sieve_interval,sched_tds_buffer,
sched_tds_buffer+SCHED_TDS_BUFSIZE-4);
#else
b0= sched_tds_buffer;
b0_ub= b0+SCHED_TDS_BUFSIZE;
for(;x+6<x_ub;x= x+8){
if((sieve_interval[*x]|sieve_interval[*(x+2)]|
sieve_interval[*(x+4)]|sieve_interval[*(x+6)])==0)
continue;
if(sieve_interval[x[0]]!=0)*(b0++)= x;
if(sieve_interval[x[2]]!=0)*(b0++)= x+2;
if(sieve_interval[x[4]]!=0)*(b0++)= x+4;
if(sieve_interval[x[6]]!=0)*(b0++)= x+6;
if(b0+4> b0_ub)goto sched_tds1;
}
for(;x<x_ub;x+= 2){
if(sieve_interval[*x]!=0)*(b0++)= x;
}
#endif
sched_tds1:
for(b1= sched_tds_buffer;b1<b0;b1++){
u16_t*y;
u32_t fbi;
y= *(b1);
if(schedules[side][j].schedule[l+1][k]<=y){
do{
l++;
if(l>=schedules[side][j].n_pieces)Schlendrian("XXX\n");
}while(schedules[side][j].schedule[l+1][k]<=y);
fbi_offset= schedules[side][j].fbi_bounds[l];
}
fbi= fbi_offset+*(y+1-2*SCHED_SI_OFFS);
*(tds_fbi_curpos[sieve_interval[*y]-1]++)= fbi;
#ifdef TDS_FB_PREFETCH
TDS_FB_PREFETCH(FB[side]+fbi);
#endif
}
}
#else
u32_t l,k;

k= schedules[side][j].current_strip++;
for(l= 0;l<schedules[side][j].n_pieces;l++){
u16_t*x,*x_ub;
u32_t fbi_offset;

x_ub= schedules[side][j].schedule[l+1][k];
fbi_offset= schedules[side][j].fbi_bounds[l];
for(x= schedules[side][j].schedule[l][k]+SCHED_SI_OFFS;x+6<x_ub;x+= 8){
unsigned char z;
if((sieve_interval[*x]|sieve_interval[*(x+2)]|
sieve_interval[*(x+4)]|sieve_interval[*(x+6)])==0){
continue;
}
if((z= sieve_interval[*x])!=0)
*(tds_fbi_curpos[z-1]++)= fbi_offset+*(x+1-2*SCHED_SI_OFFS);
if((z= sieve_interval[*(x+2)])!=0)
*(tds_fbi_curpos[z-1]++)= fbi_offset+*(x+3-2*SCHED_SI_OFFS);
if((z= sieve_interval[*(x+4)])!=0)
*(tds_fbi_curpos[z-1]++)= fbi_offset+*(x+5-2*SCHED_SI_OFFS);
if((z= sieve_interval[*(x+6)])!=0)
*(tds_fbi_curpos[z-1]++)= fbi_offset+*(x+7-2*SCHED_SI_OFFS);
}
while(x<x_ub){
unsigned char z;

if((z= sieve_interval[*x])!=0)
*(tds_fbi_curpos[z-1]++)= fbi_offset+*(x+1-2*SCHED_SI_OFFS);
x+= 2;
}
}
#endif
#endif
}
}
newclock= clock();
tds3_clock[side]+= newclock-last_tdclock;
last_tdclock= newclock;

/*:123*//*125:*/
// #line 3494 "gnfs-lasieve4e.w"

{
u32_t i;

for(i= 0;i<UCHAR_MAX&&i<ncand;i++){
u32_t*p;

for(p= tds_fbi[i];p<tds_fbi_curpos[i];p++)
*p= FB[side][*p];
}
}

/*:125*//*126:*/
// #line 3507 "gnfs-lasieve4e.w"

{
u16_t*x;

#ifdef ASM_TDSLINIE
x= smalltdsieve_auxbound;
if(x<smallsieve_auxbound[side][4]){
tdslinie(x,smallsieve_auxbound[side][4],sieve_interval,tds_fbi_curpos);
x= smallsieve_auxbound[side][4];
}
#else
for(x= smalltdsieve_auxbound;
x<smallsieve_auxbound[side][4];x= x+4){
u32_t p,r,pr;
unsigned char*y;

p= x[0];
pr= x[1];
r= x[3];
modulo32= p;
for(y= sieve_interval;y<sieve_interval+L1_SIZE;y+= n_i){
unsigned char*yy,*yy_ub;

yy_ub= y+n_i-3*p;
yy= y+r;
while(yy<yy_ub){
unsigned char o;
o= (*yy)|(*(yy+p));
yy+= 2*p;
if((o|(*yy)|(*(yy+p)))!=0){
yy= yy-2*p;
if(*yy!=0)
*(tds_fbi_curpos[*yy-1]++)= p;
if(*(yy+p)!=0)
*(tds_fbi_curpos[*(yy+p)-1]++)= p;
yy+= 2*p;
if(*yy!=0)
*(tds_fbi_curpos[*yy-1]++)= p;
if(*(yy+p)!=0)
*(tds_fbi_curpos[*(yy+p)-1]++)= p;
}
yy+= 2*p;
}
yy_ub+= 2*p;
if(yy<yy_ub){
if(((*yy)|(*(yy+p)))!=0){
if(*yy!=0)
*(tds_fbi_curpos[*yy-1]++)= p;
if(*(yy+p)!=0)
*(tds_fbi_curpos[*(yy+p)-1]++)= p;
}
yy+= 2*p;
}
yy_ub+= p;
if(yy<yy_ub){
if(*yy!=0)
*(tds_fbi_curpos[*yy-1]++)= p;
}
r= modadd32(r,pr);
}
x[3]= r;
}
#endif
#ifdef ASM_TDSLINIE3
if(x<smallsieve_auxbound[side][3]){
tdslinie3(x,smallsieve_auxbound[side][3],sieve_interval,tds_fbi_curpos);
x= smallsieve_auxbound[side][3];
}
#else
for(;x<smallsieve_auxbound[side][3];x= x+4){
u32_t p,r,pr;
unsigned char*y;

p= x[0];
pr= x[1];
r= x[3];
modulo32= p;
for(y= sieve_interval;y<sieve_interval+L1_SIZE;y+= n_i){
unsigned char*yy,*yy_ub;

yy_ub= y+n_i;
yy= y+r;
if(((*yy)|(*(yy+p))|(*(yy+2*p)))!=0){
if(*yy!=0)*(tds_fbi_curpos[*yy-1]++)= p;
if(*(yy+p)!=0)*(tds_fbi_curpos[*(yy+p)-1]++)= p;
if(*(yy+2*p)!=0)*(tds_fbi_curpos[*(yy+2*p)-1]++)= p;
}
yy+= 3*p;
if(yy<yy_ub){
if(*yy!=0)
*(tds_fbi_curpos[*yy-1]++)= p;
}
r= modadd32(r,pr);
}
x[3]= r;
}
#endif
#ifdef ASM_TDSLINIE2
if(x<smallsieve_auxbound[side][2]){
tdslinie2(x,smallsieve_auxbound[side][2],sieve_interval,tds_fbi_curpos);
x= smallsieve_auxbound[side][2];
}
#else
for(;x<smallsieve_auxbound[side][2];x= x+4){
u32_t p,r,pr;
unsigned char*y;

p= x[0];
pr= x[1];
r= x[3];
modulo32= p;
for(y= sieve_interval;y<sieve_interval+L1_SIZE;y+= n_i){
unsigned char*yy,*yy_ub;

yy_ub= y+n_i;
yy= y+r;
if(((*yy)|(*(yy+p)))!=0){
if(*yy!=0)*(tds_fbi_curpos[*yy-1]++)= p;
if(*(yy+p)!=0)*(tds_fbi_curpos[*(yy+p)-1]++)= p;
}
yy+= 2*p;
if(yy<yy_ub){
if(*yy!=0)
*(tds_fbi_curpos[*yy-1]++)= p;
}
r= modadd32(r,pr);
}
x[3]= r;
}
#endif
#ifdef ASM_TDSLINIE1
if(x<smallsieve_auxbound[side][1]){
tdslinie1(x,smallsieve_auxbound[side][1],sieve_interval,tds_fbi_curpos);
x= smallsieve_auxbound[side][1];
}
#else
for(;x<smallsieve_auxbound[side][1];x= x+4){
u32_t p,r,pr;
unsigned char*y;

p= x[0];
pr= x[1];
r= x[3];
modulo32= p;
for(y= sieve_interval;y<sieve_interval+L1_SIZE;y+= n_i){
unsigned char*yy,*yy_ub;

yy_ub= y+n_i;
yy= y+r;
if(*yy!=0)*(tds_fbi_curpos[*yy-1]++)= p;
yy+= p;
if(yy<yy_ub){
if(*yy!=0)
*(tds_fbi_curpos[*yy-1]++)= p;
}
r= modadd32(r,pr);
}
x[3]= r;
}
#endif
#ifdef ASM_TDSLINIE0
if(x<smallsieve_auxbound[side][0]){
tdslinie0(x,smallsieve_auxbound[side][0],sieve_interval,tds_fbi_curpos);
x= smallsieve_auxbound[side][0];
}
#else
for(;x<smallsieve_auxbound[side][0];x= x+4){
u32_t p,r,pr;
unsigned char*y;

p= x[0];
pr= x[1];
r= x[3];
modulo32= p;
for(y= sieve_interval;y<sieve_interval+L1_SIZE;y+= n_i){
unsigned char*yy,*yy_ub;

yy_ub= y+n_i;
yy= y+r;
if(yy<yy_ub){
if(*yy!=0)
*(tds_fbi_curpos[*yy-1]++)= p;
}
r= modadd32(r,pr);
}
x[3]= r;
}
#endif
newclock= clock();
tds1_clock[side]+= newclock-last_tdclock;
last_tdclock= newclock;
}

/*:126*/
// #line 3177 "gnfs-lasieve4e.w"

last_j= 0;
for(ci= 0,nc1= 0;ci<ncand;ci++){
u32_t*fbp_buf;
u32_t*fbp_ptr;
u16_t st_i,true_j;
i32_t true_i;

u32_t coll;
/*113:*/
// #line 3116 "gnfs-lasieve4e.w"

{
u16_t jj;

strip_j= cand[ci]>>i_bits;
jj= j_offset+strip_j;
strip_i= cand[ci]&(n_i-1);
st_i= 2*strip_i+(oddness_type==2?0:1);
true_j= 2*jj+(oddness_type==1?0:1);
}

/*:113*/
// #line 3186 "gnfs-lasieve4e.w"

if(strip_j!=last_j){
u16_t j_step;
if(strip_j<=last_j)
Schlendrian("TD: Not sorted\n");
j_step= strip_j-last_j;
last_j= strip_j;
/*127:*/
// #line 3702 "gnfs-lasieve4e.w"

#ifdef MMX_TD
MMX_TdUpdate(side,j_step);
#else
{
u32_t i;
u16_t*x,*y;

y= smalltdsieve_aux[side][j_step-1];
for(i= 0,x= smallsieve_aux[side];x<smallsieve_auxbound[side][0];i++,x+= 4){
modulo32= x[0];
if(modulo32> p_bound)break;
x[3]= modadd32((u32_t)x[3],(u32_t)y[i]);
}
}
#endif
{
u16_t*x;
for(x= smallpsieve_aux[side];x<smallpsieve_aux_ub[side];x+= 3){
modulo32= x[0];
x[2]= modsub32(x[2],(j_step)%modulo32);
}
}

/*:127*/
// #line 3193 "gnfs-lasieve4e.w"

}
true_i= (i32_t)st_i-(i32_t)i_shift;
/*128:*/
// #line 3727 "gnfs-lasieve4e.w"

mpz_set_si(aux1,true_i);
mpz_mul_si(aux1,aux1,a0);
mpz_set_si(aux2,a1);
mpz_mul_ui(aux2,aux2,(u32_t)true_j);
mpz_add(sr_a,aux1,aux2);

/*:128*//*129:*/
// #line 3736 "gnfs-lasieve4e.w"

mpz_set_si(aux1,true_i);
mpz_mul_si(aux1,aux1,b0);
mpz_set_si(aux2,b1);
mpz_mul_ui(aux2,aux2,(u32_t)true_j);
mpz_add(sr_b,aux1,aux2);

/*:129*//*130:*/
// #line 3744 "gnfs-lasieve4e.w"

if(mpz_sgn(sr_b)<0){
mpz_neg(sr_b,sr_b);
mpz_neg(sr_a,sr_a);
}

/*:130*/
// #line 3196 "gnfs-lasieve4e.w"

/*131:*/
// #line 3751 "gnfs-lasieve4e.w"

{
u32_t i;
i= 1;
mpz_set(aux2,sr_a);
mpz_set(aux1,poly[side][0]);
for(;;){

mpz_mul(aux1,aux1,sr_b);
mpz_mul(aux3,aux2,poly[side][i]);
mpz_add(aux1,aux1,aux3);
if(++i> poldeg[side])break;

mpz_mul(aux2,aux2,sr_a);
}
}

/*:131*/
// #line 3197 "gnfs-lasieve4e.w"

if(td_buf_alloc[side]<nfbp+mpz_sizeinbase(aux1,2)){

td_buf_alloc[side]+= 1024;
while(td_buf_alloc[side]<nfbp+mpz_sizeinbase(aux1,2)){
td_buf_alloc[side]+= 1024;
}
td_buf[side]= xrealloc(td_buf[side],td_buf_alloc[side]*sizeof(**td_buf));
if(side==first_td_side){
u32_t i,*oldptr;

oldptr= td_buf1[0];
for(i= 0;i<=nc1;i++)
td_buf1[i]= td_buf[side]+(td_buf1[i]-oldptr);
}
}
if(side==first_td_side)fbp_buf= td_buf1[nc1];
else fbp_buf= td_buf[side];
fbp_ptr= fbp_buf;
/*132:*/
// #line 3769 "gnfs-lasieve4e.w"

{
int np,x;

x= fss_sv[ci];
np= tds_fbi_curpos[x]-tds_fbi[x];
memcpy(fbp_ptr,tds_fbi[x],np*sizeof(*fbp_ptr));
fbp_ptr+= np;
}

/*:132*/
// #line 3216 "gnfs-lasieve4e.w"

/*133:*/
// #line 3780 "gnfs-lasieve4e.w"

{
u16_t*x;

#ifndef MMX_TD
#ifdef PREINVERT
/*134:*/
// #line 3808 "gnfs-lasieve4e.w"

{
u32_t*p_inv;

p_inv= smalltd_pi[side];
for(x= smallsieve_aux[side];
x<smallsieve_auxbound[side][0]&&*x<=p_bound;x+= 4,p_inv++){
modulo32= *x;

if(((modsub32((u32_t)strip_i,(u32_t)(x[3]))*(*p_inv))&0xffff0000)==0){
*(fbp_ptr++)= *x;
}
}
}

/*:134*/
// #line 3786 "gnfs-lasieve4e.w"

#else
for(x= smallsieve_aux[side];
x<smallsieve_auxbound[side][0]&&*x<=p_bound;x+= 4){
u32_t p;

p= *x;
if(strip_i%p==x[3])
*(fbp_ptr++)= p;
}
#endif
#else
fbp_ptr= MMX_Td(fbp_ptr,side,strip_i);
#endif
for(x= smallpsieve_aux[side];x<smallpsieve_aux_ub_pow1[side];x+= 3){
if(x[2]==0){
*(fbp_ptr++)= *x;
}
}
}

/*:133*/
// #line 3217 "gnfs-lasieve4e.w"

/*135:*/
// #line 3824 "gnfs-lasieve4e.w"

if(side==special_q_side){
*(fbp_ptr++)= special_q;
}

/*:135*/
// #line 3218 "gnfs-lasieve4e.w"

/*136:*/
// #line 3830 "gnfs-lasieve4e.w"

fbp_ptr= mpz_trialdiv(aux1,fbp_buf,fbp_ptr-fbp_buf,
tds_coll[fss_sv[ci]]==0?"td error":NULL);

/*:136*/
// #line 3219 "gnfs-lasieve4e.w"

/*137:*/
// #line 3836 "gnfs-lasieve4e.w"

if(mpz_sizeinbase(aux1,2)<=max_factorbits[side]){
n_tdsurvivors[side]++;
if(side==first_td_side){
if(mpz_sgn(aux1)> 0)
mpz_set(td_rests[nc1],aux1);
else
mpz_neg(td_rests[nc1],aux1);
cand[nc1++]= cand[ci];
td_buf1[nc1]= fbp_ptr;
nfbp= fbp_ptr-td_buf[side];
continue;
}
if(mpz_sgn(aux1)<0)mpz_neg(aux1,aux1);
#if TDS_MPQS == TDS_IMMEDIATELY
output_tdsurvivor(td_buf1[ci],td_buf1[ci+1],fbp_buf,fbp_ptr,
td_rests[ci],aux1);
#else
#if TDS_PRIMALITY_TEST == TDS_IMMEDIATELY
mpz_set(large_factors[first_td_side],td_rests[ci]);
mpz_set(large_factors[1-first_td_side],aux1);
if(primality_tests()==1){
store_tdsurvivor(td_buf1[ci],td_buf1[ci+1],fbp_buf,fbp_ptr,
large_factors[first_td_side],
large_factors[1-first_td_side]);
}
#else
store_tdsurvivor(td_buf1[ci],td_buf1[ci+1],fbp_buf,fbp_ptr,
td_rests[ci],aux1);
#endif 
#endif 

}else continue;

/*:137*/
// #line 3220 "gnfs-lasieve4e.w"

}
#ifndef MMX_TD
{
u16_t j_step;

j_step= j_per_strip-last_j;
/*127:*/
// #line 3702 "gnfs-lasieve4e.w"

#ifdef MMX_TD
MMX_TdUpdate(side,j_step);
#else
{
u32_t i;
u16_t*x,*y;

y= smalltdsieve_aux[side][j_step-1];
for(i= 0,x= smallsieve_aux[side];x<smallsieve_auxbound[side][0];i++,x+= 4){
modulo32= x[0];
if(modulo32> p_bound)break;
x[3]= modadd32((u32_t)x[3],(u32_t)y[i]);
}
}
#endif
{
u16_t*x;
for(x= smallpsieve_aux[side];x<smallpsieve_aux_ub[side];x+= 3){
modulo32= x[0];
x[2]= modsub32(x[2],(j_step)%modulo32);
}
}

/*:127*/
// #line 3227 "gnfs-lasieve4e.w"

}
#else
{
u16_t*x,j_step;
j_step= j_per_strip-last_j;
for(x= smallpsieve_aux[side];x<smallpsieve_aux_ub[side];x+= 3){
modulo32= x[0];
x[2]= modsub32(x[2],(j_step)%modulo32);
}
}
#endif
newclock= clock();
tds4_clock[side]+= newclock-last_tdclock;
last_tdclock= newclock;
ncand= nc1;
}

/*:116*/
// #line 3082 "gnfs-lasieve4e.w"

}
}

/*:111*//*141:*/
// #line 3900 "gnfs-lasieve4e.w"

#ifndef ASM_MPZ_TD

static mpz_t mpz_td_aux;
static u32_t initialized= 0;

u32_t*
mpz_trialdiv(mpz_t N,u32_t*pbuf,u32_t ncp,char*errmsg)
{
u32_t np,np1,i,e2;

if(initialized==0){
mpz_init(mpz_td_aux);
initialized= 1;
}
e2= 0;
while((mpz_get_ui(N)%2)==0){
mpz_fdiv_q_2exp(N,N,1);
e2++;
}
if(errmsg!=NULL){
for(i= 0,np= 0;i<ncp;i++){
if(mpz_fdiv_q_ui(N,N,pbuf[i])!=0)
Schlendrian("%s : %u does not divide\n",errmsg,pbuf[i]);
pbuf[np++]= pbuf[i];
}
}else{
for(i= 0,np= 0;i<ncp;i++){
if(mpz_fdiv_q_ui(mpz_td_aux,N,pbuf[i])==0){
mpz_set(N,mpz_td_aux);
pbuf[np++]= pbuf[i];
}
}
}
np1= np;
for(i= 0;i<np1;i++){
while(mpz_fdiv_q_ui(mpz_td_aux,N,pbuf[i])==0){
mpz_set(N,mpz_td_aux);
pbuf[np++]= pbuf[i];
}
}
for(i= 0;i<e2;i++)
pbuf[np++]= 2;
return pbuf+np;
}
#endif

/*:141*//*143:*/
// #line 3963 "gnfs-lasieve4e.w"

static void
store_tdsurvivor(fbp_buf0,fbp_buf0_ub,fbp_buf1,fbp_buf1_ub,lf0,lf1)
u32_t*fbp_buf0,*fbp_buf1,*fbp_buf0_ub,*fbp_buf1_ub;
mpz_t lf0,lf1;
{
size_t n0,n1,n;

/*144:*/
// #line 3994 "gnfs-lasieve4e.w"

if(total_ntds>=max_tds){
size_t i;
if(max_tds==0){
tds_fbp= xmalloc((2*MAX_TDS_INCREMENT+1)*sizeof(*tds_fbp));
tds_fbp[0]= 0;
tds_ab= xmalloc(2*MAX_TDS_INCREMENT*sizeof(*tds_ab));
tds_lp= xmalloc(2*MAX_TDS_INCREMENT*sizeof(*tds_lp));
}else{
tds_fbp= xrealloc(tds_fbp,
(2*(MAX_TDS_INCREMENT+max_tds)+1)*sizeof(*tds_fbp));
tds_ab= xrealloc(tds_ab,2*(MAX_TDS_INCREMENT+max_tds)*sizeof(*tds_ab));
tds_lp= xrealloc(tds_lp,2*(MAX_TDS_INCREMENT+max_tds)*sizeof(*tds_lp));
}
for(i= 2*total_ntds;i<2*(MAX_TDS_INCREMENT+max_tds);i++)mpz_init(tds_lp[i]);
max_tds+= MAX_TDS_INCREMENT;
}

/*:144*/
// #line 3971 "gnfs-lasieve4e.w"

if(mpz_sizeinbase(lf0,2)> max_factorbits[first_td_side]||
mpz_sizeinbase(lf1,2)> max_factorbits[1-first_td_side]){
fprintf(stderr,"large lp in store_tdsurvivor\n");
return;
}
mpz_set(tds_lp[2*total_ntds],lf0);
mpz_set(tds_lp[2*total_ntds+1],lf1);
n0= fbp_buf0_ub-fbp_buf0;
n1= fbp_buf1_ub-fbp_buf1;
n= tds_fbp[2*total_ntds];
/*145:*/
// #line 4013 "gnfs-lasieve4e.w"

if(n+n0+n1> tds_fbp_alloc){
size_t a;

a= tds_fbp_alloc;
while(a<n+n0+n1)a+= TDS_FBP_ALLOC_INCREMENT;
if(tds_fbp_alloc==0)tds_fbp_buffer= xmalloc(a*sizeof(*tds_fbp_buffer));
else tds_fbp_buffer= xrealloc(tds_fbp_buffer,a*sizeof(*tds_fbp_buffer));
tds_fbp_alloc= a;
}

/*:145*/
// #line 3982 "gnfs-lasieve4e.w"
;
memcpy(tds_fbp_buffer+n,fbp_buf0,n0*sizeof(*fbp_buf0));
n+= n0;
tds_fbp[2*total_ntds+1]= n;
memcpy(tds_fbp_buffer+n,fbp_buf1,n1*sizeof(*fbp_buf1));
tds_fbp[2*total_ntds+2]= n+n1;
tds_ab[2*total_ntds]= mpz_get_sll(sr_a);
tds_ab[2*total_ntds+1]= mpz_get_sll(sr_b);
total_ntds++;
}

/*:143*//*146:*/
// #line 4025 "gnfs-lasieve4e.w"

static int
primality_tests()
{
int s;

for(s= 0;s<2;s++){
i16_t is_prime;
u16_t s1;

s1= s^first_psp_side;
if(mpz_cmp_ui(large_factors[s1],1)==0)continue;
if(mpz_cmp(large_factors[s1],FBb_sq[s1])<0)is_prime= 1;
else is_prime= psp(large_factors[s1],1);
if(is_prime==1){
if(mpz_sizeinbase(large_factors[s1],2)> max_primebits[s1])return 0;
}else mpz_neg(large_factors[s1],large_factors[s1]);
}
return 1;
}

/*:146*//*147:*/
// #line 4047 "gnfs-lasieve4e.w"

#if (TDS_PRIMALITY_TEST != TDS_IMMEDIATELY) && (TDS_PRIMALITY_TEST != TDS_MPQS)
static void
primality_tests_all()
{
size_t i,j;

for(i= 0,j= 0;i<total_ntds;i++){
mpz_set(large_factors[first_td_side],tds_lp[2*i]);
mpz_set(large_factors[1-first_td_side],tds_lp[2*i+1]);
if(primality_tests()==0)continue;
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

/*:147*//*148:*/
// #line 4071 "gnfs-lasieve4e.w"

#if TDS_MPQS != TDS_IMMEDIATELY
static void
output_all_tdsurvivors()
{
size_t i;

for(i= 0;i<total_ntds;i++){
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

/*:148*//*149:*/
// #line 4092 "gnfs-lasieve4e.w"

static void
output_tdsurvivor(fbp_buf0,fbp_buf0_ub,fbp_buf1,fbp_buf1_ub,lf0,lf1)
u32_t*fbp_buf0,*fbp_buf1,*fbp_buf0_ub,*fbp_buf1_ub;
mpz_t lf0,lf1;
{
u32_t s,*(fbp_buffers[2]),*(fbp_buffers_ub[2]);
u32_t nlp[2];
clock_t cl;

s= first_td_side;
fbp_buffers[s]= fbp_buf0;
fbp_buffers_ub[s]= fbp_buf0_ub;
fbp_buffers[1-s]= fbp_buf1;
fbp_buffers_ub[1-s]= fbp_buf1_ub;
mpz_set(large_factors[s],lf0);
mpz_set(large_factors[1-s],lf1);

#if TDS_PRIMALITY_TEST == TDS_MPQS
if(primality_tests()==0)return;
#endif

cl= clock();
for(s= 0;s<2;s++){
u16_t s1;
i32_t i,nf;
mpz_t*mf;

s1= s^first_mpqs_side;
if(mpz_sgn(large_factors[s1])> 0){
if(mpz_cmp_ui(large_factors[s1],1)==0)
nlp[s1]= 0;
else{
nlp[s1]= 1;
mpz_set(large_primes[s1][0],large_factors[s1]);
}
continue;
}

mpz_neg(large_factors[s1],large_factors[s1]);
if((nf= mpqs_factor(large_factors[s1],max_primebits[s1],&mf))<0){
fprintf(stderr,"mpqs failed for ");
mpz_out_str(stderr,10,large_factors[s1]);
fprintf(stderr,"(a,b): ");
mpz_out_str(stderr,10,sr_a);
fprintf(stderr," ");
mpz_out_str(stderr,10,sr_b);
fprintf(stderr,"\n");
n_mpqsfail[s1]++;
break;
}
if(nf==0){

n_mpqsvain[s1]++;
break;
}
for(i= 0;i<nf;i++)
mpz_set(large_primes[s1][i],mf[i]);
nlp[s1]= nf;
}
mpqs_clock+= clock()-cl;
if(s!=2)return;

yield++;
#define OBASE 16
mpz_out_str(ofile, 10, sr_a);
fprintf(ofile, ",");
mpz_out_str(ofile, 10, sr_b);{ int numR=0;
{ u32_t *x; int i;
    
                          fprintf(ofile, ":");
                          for (i = 0; i < nlp[1]; i++) { /* rational first. */
                            if (i>0) fprintf(ofile, ",");
                            mpz_out_str(ofile, OBASE, large_primes[1][i]);
                            numR++;
                          }
                          for (x = fbp_buffers[1]; x < fbp_buffers_ub[1];x++) {
							  if ((unsigned int)*x > 1000) {
								  if (numR>0) fprintf(ofile, ",%X", (unsigned int)*x);
                                  else { fprintf(ofile, "%X", (unsigned int)*x); numR++;}
							  }
                          }
                        }
                        { int numA=0;
                          u32_t *x; int i;
    
                          fprintf(ofile, ":");
                          for (i = 0; i < nlp[0]; i++) { /* algebraic next. */
                            if (i>0) fprintf(ofile, ",");
                            mpz_out_str(ofile, OBASE, large_primes[0][i]);
                            numA++;
                          }
                          for (x = fbp_buffers[0]; x < fbp_buffers_ub[0];x++) {
                            if ((unsigned int)*x > 1000) {
								if (numA>0) fprintf(ofile, ",%X", (unsigned int)*x);
                                else { fprintf(ofile, "%X", (unsigned int)*x); numA++;}
							}
                          }
                        }
fprintf(ofile, "\n");
}}
/* 
#ifdef OFMT_CWI
#define CWI_LPB 0x100000
#define OBASE 10
{
u32_t nlp_char[2];

for(s= 0;s<2;s++){
u32_t*x,nlp1;

for(x= fbp_buffers[s],nlp1= nlp[s];x<fbp_buffers_ub[s];x++)
if(*x> CWI_LPB)
nlp1++;
if((nlp_char[s]= u32_t2cwi(nlp1))=='\0')break;
}
if(s==0){
errprintf("Conversion to CWI format failed\n");
continue;
}
#ifdef OFMT_CWI_REVERSE
fprintf(ofile,"01%c%c ",nlp_char[1],nlp_char[0]);
#else
fprintf(ofile,"01%c%c ",nlp_char[0],nlp_char[1]);
#endif
}
#else
fprintf(ofile,"W ");
#define OBASE 16
#endif
mpz_out_str(ofile,OBASE,sr_a);
fprintf(ofile," ");
mpz_out_str(ofile,OBASE,sr_b);


#ifndef OFMT_CWI_REVERSE
for(s= 0;s<2;s++){
u32_t i,*x;
#ifndef OFMT_CWI
fprintf(ofile,"\n%c",'X'+s);
#endif
for(i= 0;i<nlp[s];i++){
fprintf(ofile," ");
mpz_out_str(ofile,OBASE,large_primes[s][i]);
}
for(x= fbp_buffers[s];x<fbp_buffers_ub[s];x++){
#ifndef OFMT_CWI
fprintf(ofile," %X",*x);
#else
if(*x> CWI_LPB)
fprintf(ofile," %d",*x);
#endif
}
}
#else
for(s= 0;s<2;s++){
u32_t i,*x;
for(i= 0;i<nlp[1-s];i++){
fprintf(ofile," ");
mpz_out_str(ofile,OBASE,large_primes[1-s][i]);
}
for(x= fbp_buffers[1-s];x<fbp_buffers_ub[1-s];x++){
if(*x> CWI_LPB)
fprintf(ofile," %d",*x);
}
}
#endif
#ifndef OFMT_CWI
fprintf(ofile,"\n");
#else
fprintf(ofile,";\n");
#endif
}
*/
/*:149*//*151:*/
/* 
// #line 4238 "gnfs-lasieve4e.w"

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
*/
/*:151*//*152:*/
// #line 4252 "gnfs-lasieve4e.w"

#ifdef DEBUG
int mpout(mpz_t X)
{
mpz_out_str(stdout,10,X);
puts("");
return 1;
}
#endif

/*:152*//*154:*/
// #line 4267 "gnfs-lasieve4e.w"

void
dumpsieve(u32_t j_offset,u32_t side)
{
FILE*ofile;
char*ofn;
asprintf(&ofn,"sdump4e.ot%u.j%u.s%u",oddness_type,j_offset,side);
if((ofile= fopen(ofn,"w"))==NULL){
free(ofn);
return;
}
fwrite(sieve_interval,1,L1_SIZE,ofile);
fclose(ofile);
free(ofn);
asprintf(&ofn,"hzsdump4e.ot%u.j%u.s%u",oddness_type,j_offset,side);
if((ofile= fopen(ofn,"w"))==NULL){
free(ofn);
return;
}
fwrite(horizontal_sievesums,1,j_per_strip,ofile);
fclose(ofile);
free(ofn);
}/*:154*/
