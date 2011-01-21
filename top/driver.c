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
       				   --bbuhrow@gmail.com 7/28/10
----------------------------------------------------------------------*/

#include "yafu.h"
#include "soe.h"
#include "calc.h"
#include "yafu_string.h"
#include "util.h"
#include "factor.h"
#include "tfm.h"

#ifdef HAVE_GMP
	#include "gmp.h"
#endif

#if defined(HAVE_GMP_ECM) && defined(WIN32) 
	#include <gmp-ecm\config.h>
#elif defined(HAVE_GMP_ECM)
	#include "config.h"
#endif

// the number of recognized command line options
#define NUMOPTIONS 38
// maximum length of command line option strings
#define MAXOPTIONLEN 20

// command line options visible to driver.c
char OptionArray[NUMOPTIONS][MAXOPTIONLEN] = {
	"B1pm1",
	"B1pp1",
	"B1ecm",
	"rhomax",
	"B2pm1",
	"B2pp1",
	"B2ecm",
	"qssave",
	"siqsB",
	"siqsTF",
	"siqsR",
	"siqsT",
	"siqsNB",
	"siqsM",
	"logfile",
	"batchfile",
	"seed",
	"sigma",
	"session",
	"threads",
	"v",
	"silent",
	"pfile",
	"pscreen",
	"forceDLP",
	"fmtmax",
	"noopt",
	"vproc",
	"noecm",
	"ggnfs_dir",
	"tune_info",
	"ecm_qs_ratio",
	"ecm_gnfs_ratio",
	"one",
	"op",
	"of",
	"ou",
	"plan"};


// indication of whether or not an option needs a corresponding argument
int needsArg[NUMOPTIONS] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
	1,1,1,1,1,0,0,0,0,0,1,0,0,0,1,1,1,1,0,1,1,1,1};

// function to read the .ini file and populate options
void readINI(void);
void apply_tuneinfo(char *arg);

// functions to populate the global options with default values, and to free
// those which allocate memory
void set_default_globals(void);
void free_globals(void);

// function containing system commands to get the computer name
void get_computer_info(char *idstr);

// function to print the splash screen to file/screen
void print_splash(int is_cmdline_run, FILE *logfile, char *idstr);

// functions to make a batchfile ready to execute, and to process batchfile lines
void prepare_batchfile(char *input_exp);
int process_batchline(char *input_exp, char *indup);
void finalize_batchline();

// function to process all incoming arguments
int process_arguments(int argc, char **argv, char *input_exp);

int main(int argc, char *argv[])
{
	uint32 i=0,insize = GSTR_MAXSIZE;
	char *input_exp, *ptr, *indup;
	str_t str;
	z tmp;
	int nooutput,offset,slog,is_cmdline_run=0;
	FILE *logfile;

	//the input expression
	input_exp = (char *)malloc(GSTR_MAXSIZE*sizeof(char));
	indup = (char *)malloc(GSTR_MAXSIZE*sizeof(char));
	strcpy(input_exp,"");
	sInit(&str);

	//set defaults for various things
	set_default_globals();

	//get the computer name, cache sizes, etc.  store in globals
	get_computer_info(CPU_ID_STR);	

	//now check for an .ini file, which will override these defaults
	//command line arguments will override the .ini file
	readINI();	

	//check/process input arguments
	is_cmdline_run = process_arguments(argc, argv, input_exp);
	
	// get the batchfile ready, if requested
	if (USEBATCHFILE)
	{
		prepare_batchfile(input_exp);		
		
		//batchfile jobs are command line in nature
		is_cmdline_run = 1;

		//remember the input expression
		strcpy(indup,input_exp);
	}

	//never run silently when run interactively, else the results of
	//calculations will never be displayed.
	if (!is_cmdline_run && VFLAG < 0)
		VFLAG = 0;

	//session log
	logfile = fopen(sessionname,"a");
	if (logfile == NULL)
	{
		printf("couldn't open %s for appending\n",sessionname);
		slog = 0;
	}
	else
		slog = 1;	
		
	// print the splash screen, to the logfile and depending on options, to the screen
	print_splash(is_cmdline_run, logfile, CPU_ID_STR);

	// bigint used in this routine
	zInit(&tmp);
	
	//start the calculator
	//right now this just allocates room for user variables
	calc_init();				
		
	if (USERSEED)
	{
		logprint(logfile,"User random seed:  %u\n\n",g_rand.low);
	}
	else
	{
		get_random_seeds(&g_rand);
		logprint(logfile,"New random seeds: %u, %u\n\n",g_rand.hi,g_rand.low);
	}
	fflush(logfile);

//	printf("WARNING: constant seed is set\n");
//	g_rand.hi = 123;
//	g_rand.low = 123;
	srand(g_rand.low);
#if BITS_PER_DIGIT == 64
	LCGSTATE = (uint64)g_rand.hi << 32 | (uint64)g_rand.low;
#else
	LCGSTATE = g_rand.low;
#endif

	//command line
	while (1)
	{
		//handle a batch file, if passed in.
		if (USEBATCHFILE)
		{
			i = process_batchline(input_exp, indup);
			if (i == 1)
			{
				finalize_batchline();
				break;
			}
			else if (i == 2)
				continue;
		}
		else if (!is_cmdline_run)
		{
			// get command from user
			fgets(input_exp,GSTR_MAXSIZE,stdin);
			while (1)
			{
				if (input_exp[strlen(input_exp) - 1] == 13 || input_exp[strlen(input_exp) - 1] == 10)
				{
					//replace with a null char and continue
					printf("\n");
					fflush(stdout);
					input_exp[strlen(input_exp) - 1] = '\0';
					break;
				}
				else
				{
					//last char is not a carriage return means
					//the input is longer than allocated.
					//reallocate and get another chunk
					insize += GSTR_MAXSIZE;
					input_exp = (char *)realloc(input_exp,insize*sizeof(char));
					if (input_exp == NULL)
					{
						printf("couldn't reallocate string when parsing\n");
						exit(-1);
					}
					fgets(input_exp+strlen(input_exp),GSTR_MAXSIZE,stdin);
				}
			}	
		}
		else
		{
			// input expression already read in.  nothing to do here.

		}
		
		//search for substring help in input
		ptr = strstr(input_exp,"help");
		if (ptr != NULL)
			helpfunc(input_exp);
		else if ((strcmp(input_exp,"quit") == 0) || (strcmp(input_exp,"exit") == 0))
			break;
		else
		{
			logprint(logfile,"Processing expression: %s\n\n",input_exp);
			toStr(input_exp,&str);

			preprocess(&str);
			strcpy(input_exp,str.s);

			//detect an '=' operation, and verify the destination of the = is valid
			//pass on everything to the right of the = to the calculator
			if ((ptr = strchr(str.s,'=')) != NULL)
			{
				*ptr = '\0';
				if (invalid_dest(str.s))
				{
					printf("invalid destination %s\n",str.s);
					offset = ptr-str.s+1;
					sFree(&str);
					sInit(&str);
					memcpy(str.s,input_exp+offset,(GSTR_MAXSIZE-offset) * sizeof(char));
					strcpy(input_exp,"ans");
					str.nchars = strlen(str.s) + 1;
				}
				else
				{
					offset = ptr-str.s+1;
					sFree(&str);
					sInit(&str);
					memcpy(str.s,input_exp+offset,(GSTR_MAXSIZE-offset) * sizeof(char));

					input_exp[offset-1] = '\0';
					str.nchars = strlen(str.s) + 1;
				}
			}
			else
				strcpy(input_exp,"ans");

			//look for a trailing semicolon
			if (str.s[str.nchars-2] == ';')
			{
				nooutput = 1;
				str.s[str.nchars-2] = '\0';
				str.nchars--;
			}
			else
				nooutput = 0;

			if (!calc(&str))
			{
				if (strcmp(str.s,"") != 0)
				{
					str2hexz(str.s,&tmp);
					if (set_uvar(input_exp,&tmp))
						new_uvar(input_exp,&tmp);
					if (nooutput == 0)
					{
						if (OBASE == DEC)
						{
							if (VFLAG >= 0)
								printf("\n%s = %s\n\n",input_exp,z2decstr(&tmp,&gstr1));
						}
						else if (OBASE == HEX)
						{
							if (strstr(str.s,"0x") != NULL)
							{
								if (VFLAG >= 0)
									printf("\n%s = %s\n\n",input_exp,str.s);
							}
							else
							{
								z2hexstr(&tmp,&str);
								if (VFLAG >= 0)
									printf("\n%s = %s\n\n",input_exp,str.s);
							}
						}
					}
				}
			}
			
			//free_factor_list();
		}

#ifdef WIN32
		fflush(stdin);	//important!  otherwise scanf will process printf's output
		
#else
		if (!is_cmdline_run)
		{
			fflush(stdin);	//important!  otherwise scanf will process printf's output
			fflush(stdout);
		}
#endif

		input_exp = (char *)realloc(input_exp,GSTR_MAXSIZE*sizeof(char));
		if (input_exp == NULL)
		{
			printf("couldn't reallocate string during cleanup\n");
			exit(-1);
		}
		input_exp[0] = '\0';

		if (is_cmdline_run)
		{
			if (USEBATCHFILE)
			{
				// the line from the batchfile finished.  make the temporary file
				// created in processs_batchline the new batchfile, with the line
				// we just finished removed.
				finalize_batchline();
			}
			else
				break;
		}
		else
			printf(">> ");

	}

	if (slog)
		fclose(logfile);

	calc_finalize();
	free_globals();
	sFree(&str);
	free(input_exp);
	free(indup);	
	zFree(&tmp);

	return 1;
}

void readINI(void)
{
	FILE *doc;
	char *str;
	char *key;
	char *value;
	int len;

	doc = fopen("yafu.ini","r");

	if (doc == NULL)
		return;

	str = (char *)malloc(1024*sizeof(char));
	while (fgets(str,1024,doc) != NULL)
	{
		//if first character is a % sign, skip this line
		if (str[0] == '%')
			continue;

		//if last character of line is newline, remove it
		do 
		{
			len = strlen(str);
			if (str[len - 1] == 10)
				str[len - 1] = '\0';
			else if (str[len - 1] == 13)
				str[len - 1] = '\0';
			else
				break;
		} while (len > 0);

		//read keyword by looking for an equal sign
		key = strtok(str,"=");

		if (key == NULL)
		{
			printf("Invalid line in yafu.ini, use Keyword=Value pairs"
				"See docfile.txt for valid keywords");
			continue;
		}

		//read value
		value = strtok((char *)0,"=");

		if (value == NULL)
		{
			printf("expected value after keyword %s\n",key);
			continue;
		}

		//apply the option... same routine command line options use
		applyOpt(key,value);
	}

	fclose(doc);
	free(str);

	return;
}

void helpfunc(char *s)
{
	//search the docfile for the right entry
	//just search for the heading, and print everything until
	//the next heading is found
	FILE *doc;
	char *func;
	char *str;
	int printtopic = 0;
	int j;

	j=4;
	if (s[j] == '\0')
		j = 0;
	else
		while (isspace(s[j])) j++;		//skip white space
	func = s + j;

	//func now points to a string with the desired help topic
	//open the doc file and search matching topics
	doc = fopen("docfile.txt","r");
	if (doc == NULL)
	{
		printf("documentation file not found\n");
		return;
	}

	printf("searching for help on '%s'\n",func);
	str = (char *)malloc(1024*sizeof(char));
	while (!feof(doc))
	{
		//read a line
		fgets(str,1024,doc);
		//is this a header?
		//printf("(%d) %s",strlen(str),str);
		if ((str[0] == '[') && (str[strlen(str) - 2] == ']'))
		{
			//printf("in printtopic if\n");
			if (printtopic == 1)
				break;
			printtopic = 0;
			//does it match our topic?
			str[strlen(str) - 2] = '\0';
			if (strstr(func,str+1) != NULL)
				printtopic = 1;
		}
		else
		{
			if (printtopic)
				printf("%s",str);
		}
	}
	fclose(doc);
	free(str);
	return;
}

int invalid_dest(char *dest)
{
	//return 1 if invalid, 0 otherwise
	int i;

	if (getFunc(dest,&i) >= 0)
		return 1;	//is a function name

	//global vars are ok
	if (strcmp(dest,"POLLARD_STG1_MAX") == 0) {
		return 0;}
	else if (strcmp(dest,"POLLARD_STG2_MAX") == 0) {
		return 0;}
	else if (strcmp(dest,"WILL_STG1_MAX") == 0) {
		return 0;}
	else if (strcmp(dest,"WILL_STG2_MAX") == 0) {
		return 0;}
	else if (strcmp(dest,"ECM_STG1_MAX") == 0) {
		return 0;}
	else if (strcmp(dest,"ECM_STG2_MAX") == 0) {
		return 0;}
	else if (strcmp(dest,"BRENT_MAX_IT") == 0) {
		return 0;}
	else if (strcmp(dest,"IBASE") == 0) {
		return 0;}
	else if (strcmp(dest,"OBASE") == 0) {
		return 0;}
	else if (strcmp(dest,"QS_DUMP_CUTOFF") == 0) {
		return 0;}
	else if (strcmp(dest,"NUM_WITNESSES") == 0) {
		return 0;}
	else if (strcmp(dest,"LOGFLAG") == 0) {
		return 0;}
	else if (strcmp(dest,"VFLAG") == 0) {
		return 0;}
	else if (strcmp(dest,"PRIMES_TO_FILE") == 0) {
		return 0;}
	else if (strcmp(dest,"PRIMES_TO_SCREEN") == 0) {
		return 0;}

	//check starting char not lower case letter or _ or `
	if (dest[0] < 95 || dest[0] > 122) return 1;

	return 0;
}

int invalid_num(char *num)
{
	//check that num consists of only numeric or alphanumeric characters
	int i=0;
	int nchars = strlen(num);
	
	if (nchars == 0) return 1;

	if (num[0] == '-')
		i++;
	
	//check for 0x, 0d, 0b, or 0o.  nchars must be > 3-i in this case
	if (num[i] == '0' && num[i+1] == 'x' && ((nchars-i) > 2))
	{
		//num is hex, and can have lower or upper case alpha characters
		i += 2;
		for (;i<nchars;i++)
		{ 
			if (num[i] > 102)	//102 == f
				return 1;
			else if (num[i] < 48) 
				return 1;
			else if (num[i] > 57 && num[i] < 65)
				return 1;
			else if (num[i] > 70 && num[i] < 97)	//97 == a
				return 1;
		}
	}
	else if (num[i] == '0' && num[i+1] == 'd' && ((nchars-i) > 2))
	{
		//num is dec, and can have only digits 0 - 9
		i += 2;
		for (;i<nchars;i++)
		{ 
			if (num[i] < 48 || num[i] > 57) 
				return 1;
		}
	}
	else if (num[i] == '0' && num[i+1] == 'b' && ((nchars-i) > 2))
	{
		//num is bin, and can have only digits 0 - 1
		i += 2;
		for (;i<nchars;i++)
		{ 
			if (num[i] < 48 || num[i] > 49) 
				return 1;
		}
	}
	else if (num[i] == '0' && num[i+1] == 'o' && ((nchars-i) > 2))
	{
		//num is oct, and can have only digits 0 - 7
		i += 2;
		for (;i<nchars;i++)
		{ 
			if (num[i] < 48 || num[i] > 55) 
				return 1;
		}
	}
	else
	{
		//no base designator, go by IBASE
		if (IBASE == HEX)
		{
			//num is hex, and can have only upper case alpha characters
			for (;i<nchars;i++)
			{ 
				if (num[i] < 48) 
					return 1;
				else if (num[i] > 57 && num[i] < 65)
					return 1;
				else if (num[i] > 70)	//70 == F
					return 1;
			}
		}
		else if (IBASE == DEC)
		{
			//num is dec, and can have only digits 0 - 9
			for (;i<nchars;i++)
			{ 
				if (num[i] < 48 || num[i] > 57) 
					return 1;
			}
		}
		else if (IBASE == BIN)
		{
			//num is bin, and can have only digits 0 - 1
			for (;i<nchars;i++)
			{ 
				if (num[i] < 48 || num[i] > 49) 
					return 1;
			}
		}
		else if (IBASE == OCT)
		{
			//num is oct, and can have only digits 0 - 7
			for (;i<nchars;i++)
			{ 
				if (num[i] < 48 || num[i] > 55) 
					return 1;
			}
		}
	}

	return 0;
}

void prepare_batchfile(char *input_exp)
{
	char *ptr;
	
	//look for @ symbol in input expression
	ptr = strchr(input_exp,'@');
	if (ptr == NULL)
	{
		printf("missing variable indicator (@) in input expression\n");
		printf("ignoring any input expression: interpreting batchfile lines as input expressions\n");
		sprintf(input_exp,"@");
	}

	return;
}

int process_arguments(int argc, char **argv, char *input_exp)
{
	int is_cmdline_run=0;
	FILE *in = stdin;

	//now check for and handle any incoming arguments, whatever
	//their source.  
	if (argc > 1)
	{
		//process arguments

		if (argv[1][0] == '-')
		{
			//then there are just flags, no expression.  start up normally
			//after processing the flags
			process_flags(argc-1, &argv[1]);
		}
		else
		{
			//assume its a command.  execute once, after processing any
			//flags.  an invalid command will be handled by calc.
			is_cmdline_run = 1;
			strcpy(input_exp,argv[1]);
			if (argc > 2)
				process_flags(argc-2, &argv[2]);
		}		
	}
	else
	{
		//else, need to check for input from a redirect or pipe.
		//this is done differently on unix like systems vs. windows
#ifdef WIN32	//not complete, but ok for now
		fseek(in,-1,SEEK_END);
		if (ftell(in) >= 0)
		{
			rewind(in);
			fgets(input_exp,1024,in);
			is_cmdline_run = 1;
		}

#else
		if (isatty(fileno(in)) == 0)
		{			
			fgets(input_exp,1024,in);
			is_cmdline_run = 1;
		}
#endif
	}

	return is_cmdline_run;

}

void print_splash(int is_cmdline_run, FILE *logfile, char *idstr)
{
	if (VFLAG >= 0)
		printf("\n\n");

	if (VFLAG > 0 || !is_cmdline_run)
	{	
		logprint(NULL,"System/Build Info: \n");
	}
	logprint(logfile,"System/Build Info: \n");
	fflush(stdout);

#if !defined(HAVE_GMP) || !defined(HAVE_GMP_ECM)
	fflush(stdout);
#else
	if (VFLAG > 0 || !is_cmdline_run)
#ifdef _MSC_MPIR_VERSION
		printf("Using GMP-ECM %s, Powered by MPIR %s\n", VERSION,
_MSC_MPIR_VERSION);
#else
	#ifdef VERSION
		printf("Using GMP-ECM %s, Powered by GMP %d.%d.%d\n", VERSION, 
			__GNU_MP_VERSION,__GNU_MP_VERSION_MINOR,__GNU_MP_VERSION_PATCHLEVEL);
	#else
		fprintf(logfile,"Using GMP-ECM, Powered by GMP\n");
	#endif

#endif

#ifdef _MSC_MPIR_VERSION
	fprintf(logfile,"Using GMP-ECM %s, Powered by MPIR %s\n", VERSION,
_MSC_MPIR_VERSION);
#else
#ifdef VERSION
	fprintf(logfile,"Using GMP-ECM %s, Powered by GMP %d.%d.%d\n", VERSION,
		__GNU_MP_VERSION,__GNU_MP_VERSION_MINOR,__GNU_MP_VERSION_PATCHLEVEL);
#else
	fprintf(logfile,"Using GMP-ECM, Powered by GMP\n");
#endif
#endif
	fflush(stdout);
#endif

	fprintf(logfile,"cached %lu primes. pmax = %lu\n",szSOEp,spSOEprimes[szSOEp-1]);
	fprintf(logfile,"detected %s\ndetected L1 = %d bytes, L2 = %d bytes, CL = %d bytes\n",
		idstr,L1CACHE,L2CACHE,CLSIZE);
	fprintf(logfile,"measured cpu frequency ~= %f\n\n",
		MEAS_CPU_FREQUENCY);

	fflush(logfile);

	if (VFLAG > 0 || !is_cmdline_run)
	{		
		printf("detected %s\ndetected L1 = %d bytes, L2 = %d bytes, CL = %d bytes\n",
			idstr,L1CACHE,L2CACHE,CLSIZE);
		printf("measured cpu frequency ~= %f\n\n",
			MEAS_CPU_FREQUENCY);

		printf("===============================================================\n");
		printf("======= Welcome to YAFU (Yet Another Factoring Utility) =======\n");
		printf("=======             bbuhrow@gmail.com                   =======\n");
		printf("=======     Type help at any time, or quit to quit      =======\n");
		printf("===============================================================\n");
		printf("cached %lu primes. pmax = %lu\n\n",szSOEp,spSOEprimes[szSOEp-1]);
		printf("\n>> ");

	}

	return;
}

void get_computer_info(char *idstr)
{

	//read cache sizes
	yafu_get_cache_sizes(&L1CACHE,&L2CACHE);

	//figure out cpu freq in order to scale qs time estimations
	//0.1 seconds won't be very accurate, but hopefully close
	//enough for the rough scaling we'll be doing anyway.
    MEAS_CPU_FREQUENCY = measure_processor_speed() / 1.0e5;

	// run an extended cpuid command to get the cache line size, and
	// optionally print a bunch of info to the screen
	extended_cpuid(idstr, &CLSIZE, VERBOSE_PROC_INFO);

#ifdef WIN32

	sysname_sz = MAX_COMPUTERNAME_LENGTH + 1;
	GetComputerName(sysname,&sysname_sz);
	
#else

	sysname_sz = 255;
	gethostname(sysname,sysname_sz);
	sysname_sz = strlen(sysname);
	
#endif
	return;
}

void set_default_globals(void)
{
	uint64 limit, i;
	uint32 p_est;

	BRENT_MAX_IT=1000;
	FMTMAX=1000000;
	WILL_STG1_MAX=20000;
	WILL_STG2_MAX = 50 * (uint64)WILL_STG1_MAX;
	PP1_STG2_ISDEFAULT = 1;
	POLLARD_STG1_MAX=100000;
	POLLARD_STG2_MAX = 100 * (uint64)POLLARD_STG1_MAX;
	PM1_STG2_ISDEFAULT = 1;
	ECM_STG1_MAX=11000;
	ECM_STG2_MAX = 25 * (uint64)ECM_STG1_MAX;
	ECM_STG2_ISDEFAULT = 1;
	VFLAG = 0;
	VERBOSE_PROC_INFO = 0;
	LOGFLAG = 1;
	QS_DUMP_CUTOFF = 2048;

	QS_EXPONENT = 0;
	QS_MULTIPLIER = 0;
	QS_TUNE_FREQ = 0;

	GNFS_EXPONENT = 0;
	GNFS_MULTIPLIER = 0;
	GNFS_TUNE_FREQ = 0;

	QS_GNFS_XOVER = 0;

	MIN_NFS_DIGITS = 85;

	NUM_WITNESSES = 20;
	AUTO_FACTOR=0;
	WANT_ONLY_1_FACTOR = 0;
	PRIME_THRESHOLD = 100000000;
	PRIMES_TO_FILE = 0;
	PRIMES_TO_SCREEN = 0;
	SCALE = 40;
	NO_ECM = 0;
	USEBATCHFILE = 0;
	USERSEED = 0;
	SIGMA = 0;
	THREADS = 1;
	TARGET_ECM_QS_RATIO = 0.25;
	TARGET_ECM_GNFS_RATIO = 0.25;
	TARGET_ECM_SNFS_RATIO = 0.20;
	gbl_override_B_flag = 0;
	gbl_override_blocks_flag = 0;
	gbl_override_lpmult_flag = 0;
	gbl_force_DLP = 0;
	strcpy(siqs_savefile,"siqs.dat");
	strcpy(flogname,"factor.log");
	strcpy(sessionname,"session.log");
#if defined(_WIN64)
	strcpy(ggnfs_dir,"..\\..\\..\\..\\ggnfs-bin\\x64\\");
#elif defined(WIN32)
	strcpy(ggnfs_dir,"..\\..\\..\\..\\ggnfs-bin\\Win32\\");
#else
	strcpy(ggnfs_dir,"../ggnfs-bin/");
#endif
	szSOEp = 1000100;	
	
	//set some useful globals
	zInit(&zZero);
	zInit(&zOne);
	zInit(&zTwo);
	zInit(&zThree);
	zOne.val[0] = 1;
	zTwo.val[0] = 2;
	zThree.val[0] = 3;

	//whether we want to output certain info to their own files...
	WANT_OUTPUT_PRIMES = 0;
	WANT_OUTPUT_FACTORS = 0;
	WANT_OUTPUT_UNFACTORED = 0;

	//pretesting plan used by factor()
	yafu_pretest_plan = PRETEST_NORMAL;
	strcpy(plan_str,"normal");

	//global strings, used mostly for logprint stuff
	sInit(&gstr1);
	sInit(&gstr2);
	sInit(&gstr3);

	//global i/o base
	IBASE = DEC;
	OBASE = DEC;

	//find, and hold globally, primes less than some N
	p_est = (uint32)((double)szSOEp/log((double)szSOEp)*1.2);
	limit=szSOEp;
	spSOEprimes = (uint64 *)malloc((size_t) (p_est*sizeof(uint64)));
	szSOEp = spSOE(spSOEprimes,0,&limit,0);
	spSOEprimes = (uint64 *)realloc(spSOEprimes,(size_t) (szSOEp*sizeof(uint64)));

	//resident chunk of primes
	PRIMES = (uint64 *)malloc((size_t) (szSOEp*sizeof(uint64)));
	for (i=0;i<szSOEp;i++)
		PRIMES[i] = spSOEprimes[i];

	NUM_P=szSOEp;
	P_MIN=0; 
	P_MAX=spSOEprimes[(uint32)NUM_P-1];

	return;
}

void free_globals()
{
	zFree(&zZero);
	zFree(&zOne);
	zFree(&zTwo);
	zFree(&zThree);
	free(spSOEprimes);
	free(PRIMES);
	sFree(&gstr1);
	sFree(&gstr2);
	sFree(&gstr3);

	return;
}

void finalize_batchline()
{
	rename(batchfilename,"_bkup");
	rename("__tmpbatchfile",batchfilename);
	remove("_bkup");
	remove("__tmpbatchfile");

	return;
}

int process_batchline(char *input_exp, char *indup)
{
	int nChars, j, i;
	char line[GSTR_MAXSIZE], tmpline[GSTR_MAXSIZE], *ptr, *ptr2;
	FILE *batchfile, *tmpfile;

	//try to open the file
	batchfile = fopen(batchfilename,"r");	

	if (batchfile == NULL)
	{
		printf("couldn't open %s for reading\n",batchfilename);
		exit(-1);
	}	

	//load the next line of the batch file and get the expression
	//ready for processing
	strcpy(line,"");
	strcpy(input_exp,"");

	// read a line
	ptr = fgets(line,GSTR_MAXSIZE,batchfile);	

	// copy everything in the file after the line we just read to
	// a temporary file.  if the expression we just read finishes, 
	// the temporary file will become the batch file (effectively 
	// eliminating the expression from the batch file).
	tmpfile = fopen("__tmpbatchfile", "w");

	while (!feof(batchfile))
	{
		ptr2 = fgets(tmpline,GSTR_MAXSIZE,batchfile);
		if (ptr2 == NULL)
			break;
		fputs(tmpline,tmpfile);
	}
	fclose(tmpfile);

	// close the batchfile
	fclose(batchfile);

	// check the line we read
	if (ptr == NULL)
	{
		printf("fgets returned null; done processing batchfile\n");		
		return 1;
	}

	// remove LF an CRs from line
	nChars = 0;
	for (j=0; j<strlen(line); j++)
	{
		switch (line[j])
		{
		case 13:
		case 10:
			break;
		default:
			line[nChars++] = line[j];
			break;
		}
	}
	line[nChars++] = '\0';

	//ignore blank lines
	if (strlen(line) == 0)
	{
		printf("skipping blank line\n");
		return 2;
	}

	//substitute the batchfile line into the '@' symbol in the input expression
	nChars = 0;
	for (i=0; i<strlen(indup); i++)
	{
		if (indup[i] == '@')
		{
			for (j=0; j<strlen(line); j++)
				input_exp[nChars++] = line[j];
		}
		else				
			input_exp[nChars++] = indup[i];

		if (nChars > GSTR_MAXSIZE)
			input_exp = (char *)realloc(input_exp,strlen(indup) + strlen(line) + 2);
	}
	input_exp[nChars++] = '\0';

	printf("=== Starting work on batchfile expression ===\n");
	printf("%s\n",input_exp);
	printf("=============================================\n");

	return 0;
}

unsigned process_flags(int argc, char **argv)
{
    int ch = 0, i,j,valid;
	char optbuf[MAXOPTIONLEN];
	char argbuf[80];

    //argument loop
	i = 0;
	while (i < argc)
	{
		//read in the option
		ch = argv[i][0];
		if (ch != '-')
		{
			printf("no switch detected\n");
			exit(1);
		}
		strcpy(optbuf,argv[i]);

		//check if its valid
		valid = 0;
		for (j=0; j<NUMOPTIONS;j++)
		{
			if (strcmp(OptionArray[j],optbuf+1) == 0)
			{
				valid = 1;
				break;
			}
		}
		if (valid == 0)
		{
			printf("invalid option %s\n",optbuf);
			exit(1);
		}

		//check to see if this option requires an argument
		if (needsArg[j])
		{
			i++;
			if (i == argc)
			{
				printf("argument expected for %s\n",optbuf);
				exit(1);
			}
			strcpy(argbuf,argv[i]);

			//now apply -option argument
			applyOpt(optbuf+1,argbuf);
		}
		else
		{
			//apply -option
			//now apply -option argument
			applyOpt(optbuf+1,NULL);

		}
		i++;
	}

    return 1;
}

void applyOpt(char *opt, char *arg)
{
	char **ptr;
	int i;
	z tmp;

	zInit(&tmp);

	ptr = NULL;
	if (strcmp(opt,OptionArray[0]) == 0)
	{
		//argument should be all numeric
		for (i=0;i<(int)strlen(arg);i++)
		{
			if (!isdigit(arg[i]))
			{
				printf("expected numeric input for option %s\n",opt);
				exit(1);
			}
		}

		POLLARD_STG1_MAX = strtoul(arg,ptr,10);
		if (PM1_STG2_ISDEFAULT)
		{
			//stg2 hasn't been specified yet, so set it to the default value
			POLLARD_STG2_MAX = 100 * (uint64)POLLARD_STG1_MAX;
			//else, we have already specified a B2, so don't overwrite it with
			//the default
		}
	}
	else if (strcmp(opt,OptionArray[1]) == 0)
	{
		//argument should be all numeric
		for (i=0;i<(int)strlen(arg);i++)
		{
			if (!isdigit(arg[i]))
			{
				printf("expected numeric input for option %s\n",opt);
				exit(1);
			}
		}

		WILL_STG1_MAX = strtoul(arg,ptr,10);
		if (PP1_STG2_ISDEFAULT)
		{
			//stg2 hasn't been specified yet, so set it to the default value
			WILL_STG2_MAX = 50 * (uint64)WILL_STG1_MAX;
			//else, we have already specified a B2, so don't overwrite it with
			//the default
		}
	}
	else if (strcmp(opt,OptionArray[2]) == 0)
	{
		//argument should be all numeric
		for (i=0;i<(int)strlen(arg);i++)
		{
			if (!isdigit(arg[i]))
			{
				printf("expected numeric input for option %s\n",opt);
				exit(1);
			}
		}

		ECM_STG1_MAX = strtoul(arg,ptr,10);
		if (ECM_STG2_ISDEFAULT)
		{
			//stg2 hasn't been specified yet, so set it to the default value
			ECM_STG2_MAX = 25 * (uint64)ECM_STG1_MAX;
			//else, we have already specified a B2, so don't overwrite it with
			//the default
		}
	}
	else if (strcmp(opt,OptionArray[3]) == 0)
	{
		//argument should be all numeric
		for (i=0;i<(int)strlen(arg);i++)
		{
			if (!isdigit(arg[i]))
			{
				printf("expected numeric input for option %s\n",opt);
				exit(1);
			}
		}

		BRENT_MAX_IT = strtoul(arg,ptr,10);
	}
	else if (strcmp(opt,OptionArray[4]) == 0)
	{
		//argument should be all numeric
		for (i=0;i<(int)strlen(arg);i++)
		{
			if (!isdigit(arg[i]))
			{
				printf("expected numeric input for option %s\n",opt);
				exit(1);
			}
		}

		str2hexz(arg,&tmp);
		POLLARD_STG2_MAX = z264(&tmp);
		PM1_STG2_ISDEFAULT = 0;
	}
	else if (strcmp(opt,OptionArray[5]) == 0)
	{
		//argument should be all numeric
		for (i=0;i<(int)strlen(arg);i++)
		{
			if (!isdigit(arg[i]))
			{
				printf("expected numeric input for option %s\n",opt);
				exit(1);
			}
		}

		//WILL_STG2_MAX = strtoul(arg,ptr,10);
		str2hexz(arg,&tmp);
		WILL_STG2_MAX = z264(&tmp);
		PP1_STG2_ISDEFAULT = 0;
	}
	else if (strcmp(opt,OptionArray[6]) == 0)
	{
		//argument should be all numeric
		for (i=0;i<(int)strlen(arg);i++)
		{
			if (!isdigit(arg[i]))
			{
				printf("expected numeric input for option %s\n",opt);
				exit(1);
			}
		}

		//ECM_STG2_MAX = strtoul(arg,ptr,10);
		str2hexz(arg,&tmp);
		ECM_STG2_MAX = z264(&tmp);
		ECM_STG2_ISDEFAULT = 0;
	}
	else if (strcmp(opt,OptionArray[7]) == 0)
	{
		//argument is a string
	
		if (strlen(arg) < 1024)
			strcpy(siqs_savefile,arg);
		else
			printf("*** argument to savefile too long, ignoring ***\n");
	}
	else if (strcmp(opt,OptionArray[8]) == 0)
	{
		//argument should be all numeric
		for (i=0;i<(int)strlen(arg);i++)
		{
			if (!isdigit(arg[i]))
			{
				printf("expected numeric input for option %s\n",opt);
				exit(1);
			}
		}

		gbl_override_B = strtoul(arg,ptr,10);
		gbl_override_B_flag = 1;
	}
	else if (strcmp(opt,OptionArray[9]) == 0)
	{
		//argument should be all numeric
		for (i=0;i<(int)strlen(arg);i++)
		{
			if (!isdigit(arg[i]))
			{
				printf("expected numeric input for option %s\n",opt);
				exit(1);
			}
		}

		gbl_override_tf = strtoul(arg,ptr,10);
		gbl_override_tf_flag = 1;
	}
	else if (strcmp(opt,OptionArray[10]) == 0)
	{
		//argument should be all numeric
		for (i=0;i<(int)strlen(arg);i++)
		{
			if (!isdigit(arg[i]))
			{
				printf("expected numeric input for option %s\n",opt);
				exit(1);
			}
		}

		gbl_override_rel = strtoul(arg,ptr,10);
		gbl_override_rel_flag = 1;
	}
	else if (strcmp(opt,OptionArray[11]) == 0)
	{
		//argument should be all numeric
		for (i=0;i<(int)strlen(arg);i++)
		{
			if (!isdigit(arg[i]))
			{
				printf("expected numeric input for option %s\n",opt);
				exit(1);
			}
		}

		gbl_override_time = strtoul(arg,ptr,10);
		gbl_override_time_flag = 1;
	}
	else if (strcmp(opt,OptionArray[12]) == 0)
	{
		//argument should be all numeric
		for (i=0;i<(int)strlen(arg);i++)
		{
			if (!isdigit(arg[i]))
			{
				printf("expected numeric input for option %s\n",opt);
				exit(1);
			}
		}

		gbl_override_blocks = strtoul(arg,ptr,10);
		gbl_override_blocks_flag = 1;
	}
	else if (strcmp(opt,OptionArray[13]) == 0)
	{
		//argument should be all numeric
		for (i=0;i<(int)strlen(arg);i++)
		{
			if (!isdigit(arg[i]))
			{
				printf("expected numeric input for option %s\n",opt);
				exit(1);
			}
		}

		gbl_override_lpmult = strtoul(arg,ptr,10);
		gbl_override_lpmult_flag = 1;
	}
	else if (strcmp(opt,OptionArray[14]) == 0)
	{
		//argument is a string
		if (strlen(arg) < 1024)
			strcpy(flogname,arg);
		else
			printf("*** argument to logfile too long, ignoring ***\n");
	}
	else if (strcmp(opt,OptionArray[15]) == 0)
	{
		//argument is a string
		if (strlen(arg) < 80)
		{
			strcpy(batchfilename,arg);
			USEBATCHFILE = 1;
		}
		else
			printf("*** argument to batchfile too long, ignoring ***\n");
	}
	else if (strcmp(opt,OptionArray[16]) == 0)
	{
		USERSEED = 1;
		sscanf(arg,"%u,%u",&g_rand.hi,&g_rand.low);
	}
	else if (strcmp(opt,OptionArray[17]) == 0)
	{
		//argument should be all numeric
		for (i=0;i<(int)strlen(arg);i++)
		{
			if (!isdigit(arg[i]))
			{
				printf("expected numeric input for option %s\n",opt);
				exit(1);
			}
		}

		SIGMA = strtoul(arg,NULL,10);
	}
	else if (strcmp(opt,OptionArray[18]) == 0)
	{
		//argument is a string
		if (strlen(arg) < 1024)
		{
			strcpy(sessionname,arg);
		}
		else
			printf("*** argument to sessionname too long, ignoring ***\n");
	}
	else if (strcmp(opt,OptionArray[19]) == 0)
	{
		//argument should be all numeric
		for (i=0;i<(int)strlen(arg);i++)
		{
			if (!isdigit(arg[i]))
			{
				printf("expected numeric input for option %s\n",opt);
				exit(1);
			}
		}

		THREADS = strtoul(arg,NULL,10);
	}
	else if (strcmp(opt,OptionArray[20]) == 0)
	{
		VFLAG++;
	}
	else if (strcmp(opt,OptionArray[21]) == 0)
	{
		VFLAG = -1;
	}
	else if (strcmp(opt,OptionArray[22]) == 0)
	{
		PRIMES_TO_FILE = 1;
	}
	else if (strcmp(opt,OptionArray[23]) == 0)
	{
		PRIMES_TO_SCREEN = 1;
	}
	else if (strcmp(opt,OptionArray[24]) == 0)
	{
		gbl_force_DLP = 1;
	}
	else if (strcmp(opt,OptionArray[25]) == 0)
	{
		//argument should be all numeric
		for (i=0;i<(int)strlen(arg);i++)
		{
			if (!isdigit(arg[i]))
			{
				printf("expected numeric input for option %s\n",opt);
				exit(1);
			}
		}

		FMTMAX = strtoul(arg,NULL,10);
	}
	else if (strcmp(opt,OptionArray[26]) == 0)
	{
		NO_SIQS_OPT = 1;
	}
	else if (strcmp(opt,OptionArray[27]) == 0)
	{
		VERBOSE_PROC_INFO++;
	}
	else if (strcmp(opt,OptionArray[28]) == 0)
	{
		yafu_pretest_plan = PRETEST_NOECM;
	}
	else if (strcmp(opt,OptionArray[29]) == 0)
	{
		//argument is a string
		if (strlen(arg) < 1024)
			strcpy(ggnfs_dir,arg);
		else
			printf("*** argument to ggnfs_dir too long, ignoring ***\n");
	}
	else if (strcmp(opt,OptionArray[30]) == 0)
	{
		//parse the tune_info string and if it matches the current OS and CPU, 
		//set the appropriate globals
		apply_tuneinfo(arg);
	}
	else if (strcmp(opt,OptionArray[31]) == 0)
	{
		sscanf(arg, "%lf", &TARGET_ECM_QS_RATIO);
	}
	else if (strcmp(opt,OptionArray[32]) == 0)
	{
		sscanf(arg, "%lf", &TARGET_ECM_GNFS_RATIO);
	}
	else if (strcmp(opt,OptionArray[33]) == 0)
	{
		WANT_ONLY_1_FACTOR = 1;
	}
	else if (strcmp(opt,OptionArray[34]) == 0)
	{
		//argument is a string
		if (strlen(arg) < 1024)
		{
			strcpy(op_str,arg);
			WANT_OUTPUT_PRIMES = 1;
		}
		else
			printf("*** argument to -op too long, ignoring ***\n");
	}
	else if (strcmp(opt,OptionArray[35]) == 0)
	{
		//argument is a string
		if (strlen(arg) < 1024)
		{
			strcpy(of_str,arg);
			WANT_OUTPUT_FACTORS = 1;
		}
		else
			printf("*** argument to -of too long, ignoring ***\n");
	}
	else if (strcmp(opt,OptionArray[36]) == 0)
	{
		//argument is a string
		if (strlen(arg) < 1024)
		{
			strcpy(ou_str,arg);
			WANT_OUTPUT_UNFACTORED = 1;
		}
		else
			printf("*** argument to -ou too long, ignoring ***\n");
	}
	else if (strcmp(opt,OptionArray[37]) == 0)
	{
		//argument is a string
		if (strlen(arg) < 1024)
		{
			strcpy(plan_str,arg);

			// test for recognized options.  for now, ignore "custom" even though
			// it is in the enum, because it is too complicated to implement right now.
			// that's a task for another day.
			if (strcmp(plan_str, "none") == 0)
				yafu_pretest_plan = PRETEST_NONE;
			else if (strcmp(plan_str, "noecm") == 0)
				yafu_pretest_plan = PRETEST_NOECM;
			else if (strcmp(plan_str, "light") == 0)
				yafu_pretest_plan = PRETEST_LIGHT;
			else if (strcmp(plan_str, "deep") == 0)
				yafu_pretest_plan = PRETEST_DEEP;
			else if (strcmp(plan_str, "normal") == 0)
				yafu_pretest_plan = PRETEST_NORMAL;
			else			
			{
				printf("*** unknown plan option, ignoring ***\n");
				strcpy(plan_str,"normal");
			}
		}
		else
			printf("*** argument to -plan too long, ignoring ***\n");
	}
	else
	{
		printf("invalid option %s\n",opt);
		exit(1);
	}

	zFree(&tmp);
	return;
}

void apply_tuneinfo(char *arg)
{
	int i,j;
	char cpustr[80], osstr[80];

	//read up to the first comma - this is the cpu id string
	j=0;
	for (i=0; i<strlen(arg); i++)
	{
		if (arg[i] == 10) break;
		if (arg[i] == 13) break;
		if (arg[i] == ',') break;
		cpustr[j++] = arg[i];
	}
	cpustr[j] = '\0';
	i++;

	//read up to the next comma - this is the OS string
	j=0;
	for ( ; i<strlen(arg); i++)
	{
		if (arg[i] == 10) break;
		if (arg[i] == 13) break;
		if (arg[i] == ',') break;
		osstr[j++] = arg[i];
	}
	osstr[j] = '\0';

	//printf("found OS = %s and CPU = %s in tune_info field\n",osstr, cpustr);


#if defined(_WIN64)
	if ((strcmp(cpustr,CPU_ID_STR) == 0) && (strcmp(osstr, "WIN64") == 0))
	{
		printf("Applying tune_info entry for %s - %s\n",osstr,cpustr);
		
		sscanf(arg + i + 1, "%lg, %lg, %lg, %lg, %lg, %lg",
			&QS_MULTIPLIER, &QS_EXPONENT,
			&GNFS_MULTIPLIER, &GNFS_EXPONENT, 
			&QS_GNFS_XOVER, &GNFS_TUNE_FREQ);
		QS_TUNE_FREQ = GNFS_TUNE_FREQ;

	}
#elif defined(WIN32)
	if ((strcmp(cpustr,CPU_ID_STR) == 0) && (strcmp(osstr, "WIN32") == 0))
	{
		printf("Applying tune_info entry for %s - %s\n",osstr,cpustr);
		sscanf(arg + i + 1, "%lg, %lg, %lg, %lg, %lg, %lg",
			&QS_MULTIPLIER, &QS_EXPONENT,
			&GNFS_MULTIPLIER, &GNFS_EXPONENT, 
			&QS_GNFS_XOVER, &GNFS_TUNE_FREQ);
		QS_TUNE_FREQ = GNFS_TUNE_FREQ;
	}
#elif BITS_PER_DIGIT == 64
	if ((strcmp(cpustr,CPU_ID_STR) == 0) && (strcmp(osstr, "LINUX64") == 0))
	{
		printf("Applying tune_info entry for %s - %s\n",osstr,cpustr);
		
		sscanf(arg + i + 1, "%lg, %lg, %lg, %lg, %lg, %lg",
			&QS_MULTIPLIER, &QS_EXPONENT,
			&GNFS_MULTIPLIER, &GNFS_EXPONENT, 
			&QS_GNFS_XOVER, &GNFS_TUNE_FREQ);
		QS_TUNE_FREQ = GNFS_TUNE_FREQ;
	}
#else 
	if ((strcmp(cpustr,CPU_ID_STR) == 0) && (strcmp(osstr, "LINUX32") == 0))
	{
		printf("Applying tune_info entry for %s - %s\n",osstr,cpustr);
		
		sscanf(arg + i + 1, "%lg, %lg, %lg, %lg, %lg, %lg",
			&QS_MULTIPLIER, &QS_EXPONENT,
			&GNFS_MULTIPLIER, &GNFS_EXPONENT, 
			&QS_GNFS_XOVER, &GNFS_TUNE_FREQ);
		QS_TUNE_FREQ = GNFS_TUNE_FREQ;
	}
#endif	

	//printf("QS_MULTIPLIER = %lg, QS_EXPONENT = %lg\nNFS_MULTIPLIER = %lg, NFS_EXPONENT = %lg\nXOVER = %lg, TUNE_FREQ = %lg\n",
	//	QS_MULTIPLIER, QS_EXPONENT,
	//	GNFS_MULTIPLIER, GNFS_EXPONENT, 
	//	QS_GNFS_XOVER, GNFS_TUNE_FREQ);

	return;
}

//function get_random_seeds courtesy of Jason Papadopoulos
void get_random_seeds(rand_t *r) {

	uint32 tmp_seed1, tmp_seed2;

	/* In a multithreaded program, every msieve object
	   should have two unique, non-correlated seeds
	   chosen for it */

	//in YAFU, make them available everywhere, by putting them in
	//a global structure that holds them.

#ifndef WIN32

	FILE *rand_device = fopen("/dev/urandom", "r");

	if (rand_device != NULL) {

		/* Yay! Cryptographic-quality nondeterministic randomness! */

		fread(&tmp_seed1, sizeof(uint32), (size_t)1, rand_device);
		fread(&tmp_seed2, sizeof(uint32), (size_t)1, rand_device);
		fclose(rand_device);
	}
	else

#endif
	{
		/* <Shrug> For everyone else, sample the current time,
		   the high-res timer (hopefully not correlated to the
		   current time), and the process ID. Multithreaded
		   applications should fold in the thread ID too */

		uint64 high_res_time = yafu_read_clock();
		tmp_seed1 = ((uint32)(high_res_time >> 32) ^
			     (uint32)time(NULL)) * 
			    (uint32)getpid();
		tmp_seed2 = (uint32)high_res_time;
	}

	/* The final seeds are the result of a multiplicative
	   hash of the initial seeds */

	r->low = tmp_seed1 * ((uint32)40499 * 65543);
	r->hi = tmp_seed2 * ((uint32)40499 * 65543);
}

