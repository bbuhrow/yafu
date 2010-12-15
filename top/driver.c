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

char OptionArray[28][20] = {
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
	"vproc"};
int needsArg[28] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,0,0};
int numOptions = 28;
void readINI(void);
void gen_masks();
void gen_maps();

/*
void gen_masks()
{
	//generate masks for the 8 primes 7,11,13,17,19,23,29,31 mod 64
	//for all possible offsets mod p
	int p[9] = {5,7,11,13,17,19,23,29,31};
	int i,offset,k;
	uint64 one = 1;

	for (i=0; i<9; i++)
	{
		//each prime
		int prime = p[i];

		printf("64 mod %d = %d\n",prime,64 % prime);
		for (offset=0;offset<prime;offset++)
		{
			uint64 mask = 0;
			uint64 byte[8];

			//each offset
			for (k=offset; k<64; k+=prime)
				mask |= (uint64)(one << (uint64)k);

			printf("0x%lx,\n",~mask);
		}

	}


	return;
}

void gen_maps()
{
	//generate a map which does the following:
	//given a starting value, n,  in (0,30], and another value, p, taking on values mod 30,
	//return the number of steps of p starting from n which equals 1
	int i,j,k,t,val;
	int steps30[8] = {1,7,11,13,17,19,23,29};
	int steps6[2] = {1,5};
	int steps210[48] = {1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 
		101, 103, 107, 109, 113, 121, 127, 131, 137, 139, 143, 149, 151, 157, 163, 167, 
		169, 173, 179, 181, 187, 191, 193, 197, 199, 209};


	//30 starting values
	for (i=0; i < 30; i++)
	{		
		printf("}\n{");
		//8 step values
		for (k=0; k < 8; k++)
		{
			//30 increments
			val = i;
			for (j=0; j < 30; j++)
			{
				//check for a residue of 1.  other residues can be generated from this table.
				if (val == 1)
				{
					printf("%d,",j);
					break;
				}
				else
				{
					val += steps30[k];
					if (val > 30)
						val -= 30;
				}
			}
		}
	}

	//6 starting values
	for (i=0; i < 6; i++)
	{		
		printf("}\n{");
		//2 step values
		for (k=0; k < 2; k++)
		{
			//6 increments
			val = i;
			for (j=0; j < 6; j++)
			{
				//check for a residue of 1.  other residues can be generated from this table.
				if (val == 1)
				{
					printf("%d,",j);
					break;
				}
				else
				{
					val += steps6[k];
					if (val > 6)
						val -= 6;
				}
			}
		}
	}

	//210 starting values
	for (i=0; i < 210; i++)
	{		
		printf("}\n{");
		//48 step values
		for (k=0; k < 48; k++)
		{
			//210 increments
			val = i;
			for (j=0; j < 210; j++)
			{
				//check for a residue of 1.  other residues can be generated from this table.
				if (val == 1)
				{
					printf("%d,",j);
					break;
				}
				else
				{
					val += steps210[k];
					if (val > 210)
						val -= 210;
				}
			}
		}
	}
	

	return;
}

*/

int main(int argc, char *argv[])
{
	uint32 i=0,insize = GSTR_MAXSIZE;
	uint64 limit;
	char *s, *ptr, *line, *indup, idstr[64];
	str_t str;
	z tmp;
	enum cpu_type cpu;
	int nooutput,offset,slog,do_once=0,cachelinesize;
	FILE *logfile, *in = stdin, *batchfile;

	//the input expression
	s = (char *)malloc(GSTR_MAXSIZE*sizeof(char));
	indup = (char *)malloc(GSTR_MAXSIZE*sizeof(char));
	strcpy(s,"");
	sInit(&str);

	//for reading batch files
	line = (char *)malloc(GSTR_MAXSIZE*sizeof(char));
	strcpy(line,"");

	//set defaults for various things
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
	NUM_WITNESSES = 20;
	AUTO_FACTOR=0;
	PRIME_THRESHOLD = 100000000;
	PRIMES_TO_FILE = 0;
	PRIMES_TO_SCREEN = 0;
	SCALE = 40;
	USEBATCHFILE = 0;
	USERSEED = 0;
	SIGMA = 0;
	THREADS = 1;
	gbl_override_B_flag = 0;
	gbl_override_blocks_flag = 0;
	gbl_override_lpmult_flag = 0;
	gbl_force_DLP = 0;
	strcpy(siqs_savefile,"siqs.dat");
	strcpy(flogname,"factor.log");
	strcpy(sessionname,"session.log");

	//now check for an .ini file, which will override these defaults
	//command line arguments will override the .ini file
	readINI();

	//get the computer name
#ifdef WIN32

	sysname_sz = MAX_COMPUTERNAME_LENGTH + 1;
	GetComputerName(sysname,&sysname_sz);
	
#else

	sysname_sz = 255;
	gethostname(sysname,sysname_sz);
	sysname_sz = strlen(sysname);
	
#endif

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
			do_once = 1;
			strcpy(s,argv[1]);
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
			fgets(s,1024,in);
			do_once = 1;
		}

#else
		if (isatty(fileno(in)) == 0)
		{
			
			fgets(s,1024,in);
			printf("processing redirect %s\n",s);
			do_once = 1;
		}
#endif
	}

	//don't get as many primes if this is a cmd line job.
	//this many should be unnoticable to the user, and
	//routines that need more know how to get more.
	if (do_once)
		szSOEp = 1000100;
	else
		szSOEp = 10000100;

	//never run silently when run interactively, else the results of
	//calculations will never be displayed.
	if (!do_once && VFLAG < 0)
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

	if (VFLAG >= 0)
		printf("\n\n");


	if (VFLAG > 0 || !do_once)
	{	
		logprint(NULL,"System/Build Info: \n");
	}
	logprint(logfile,"System/Build Info: \n");
	fflush(stdout);

#if !defined(HAVE_GMP) || !defined(HAVE_GMP_ECM)
	fflush(stdout);
#else
	if (VFLAG > 0 || !do_once)
#ifdef _MSC_MPIR_VERSION
		printf("Using GMP-ECM %s, Powered by MPIR %s\n", VERSION,
_MSC_MPIR_VERSION);
#else
		printf("Using GMP-ECM %s, Powered by GMP %d.%d.%d\n", VERSION, 
			__GNU_MP_VERSION,__GNU_MP_VERSION_MINOR,__GNU_MP_VERSION_PATCHLEVEL);
#endif
#ifdef _MSC_MPIR_VERSION
	fprintf(logfile,"Using GMP-ECM %s, Powered by MPIR %s\n", VERSION,
_MSC_MPIR_VERSION);
#else
	fprintf(logfile,"Using GMP-ECM %s, Powered by GMP %d.%d.%d\n", VERSION,
		__GNU_MP_VERSION,__GNU_MP_VERSION_MINOR,__GNU_MP_VERSION_PATCHLEVEL);
#endif
	fflush(stdout);
#endif


	//find, and hold globally, primes less than some N
	i = (uint32)((double)szSOEp/log((double)szSOEp)*1.2);
	limit=szSOEp;
	spSOEprimes = (uint64 *)malloc((size_t) (i*sizeof(uint64)));
	szSOEp = spSOE(spSOEprimes,0,&limit,0);
	spSOEprimes = (uint64 *)realloc(spSOEprimes,(size_t) (szSOEp*sizeof(uint64)));

	//resident chunk of primes
	PRIMES = (uint64 *)malloc((size_t) (szSOEp*sizeof(uint64)));
	for (i=0;i<szSOEp;i++)
		PRIMES[i] = spSOEprimes[i];

	NUM_P=szSOEp;
	P_MIN=0; 
	P_MAX=spSOEprimes[(uint32)NUM_P-1];

	get_cache_sizes(&L1CACHE,&L2CACHE);
	cpu = get_cpu_type();

	//figure out cpu freq in order to scale qs time estimations
	//0.1 seconds won't be very accurate, but hopefully close
	//enough for the rough scaling we'll be doing anyway.
    MEAS_CPU_FREQUENCY = measure_processor_speed() / 1.0e5;

	extended_cpuid(idstr, &cachelinesize, VERBOSE_PROC_INFO);
		
	fprintf(logfile,"cached %lu primes. pmax = %lu\n",szSOEp,spSOEprimes[szSOEp-1]);
	fprintf(logfile,"detected %s\ndetected L1 = %d bytes, L2 = %d bytes, CL = %d bytes\n",
		idstr,L1CACHE,L2CACHE,cachelinesize);
	fprintf(logfile,"measured cpu frequency ~= %f\n\n",
		MEAS_CPU_FREQUENCY);

	fflush(logfile);

	if (VFLAG > 0 || !do_once)
	{		
		printf("detected %s\ndetected L1 = %d bytes, L2 = %d bytes, CL = %d bytes\n",
			idstr,L1CACHE,L2CACHE,cachelinesize);
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

	//set some useful globals
	zInit(&zZero);
	zInit(&zOne);
	zInit(&zTwo);
	zInit(&zThree);
	zOne.val[0] = 1;
	zTwo.val[0] = 2;
	zThree.val[0] = 3;
	zInit(&tmp);

	//global strings, used mostly for logprint stuff
	sInit(&gstr1);
	sInit(&gstr2);
	sInit(&gstr3);

	//global i/o base
	IBASE = DEC;
	OBASE = DEC;

	//start the calculator
	//right now this just allocates room for user variables
	calc_init();

	//global list of factors
	MAX_GLOBAL_FACTORS = 256;
	global_factors = (factor_t *)malloc(MAX_GLOBAL_FACTORS * sizeof(factor_t));
	NUM_GLOBAL_FACTORS = 0;

	MSC_ASM = 0;
		
	batchfile = NULL;
	if (USEBATCHFILE)
	{
		//get the first line
		batchfile = fopen(batchfilename,"r");
		if (batchfile == NULL)
		{
			printf("couldn't open %s for reading\n",batchfilename);
			exit(-1);
		}
		fgets(line,1024,batchfile);
		if (strlen(line) == 0)
		{
			printf("nothing to read in %s\n",batchfilename);
			exit(-1);
		}
		ptr = strchr(line,10);
		if (ptr == NULL)
		{
			printf("newline missing in %s\n",batchfilename);
			exit(-1);
		}
		//remember the input expression
		strcpy(indup,s);
	}
		
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

	//gen_masks();
	//gen_maps();

	//command line
	while (1)
	{
		//handle a batch file, if passed in.
		if (USEBATCHFILE)
		{
			if (do_once)
			{
				//if this is a command line job, then we need to do two things
				//1.) preserve the do_once nature of the command line job
				//2.) loop over the inputs in the batchfile
				//do this by loading in ith line of batchfile into the @
				//symbol of the input expression, and looping over the number
				//of lines in the batch file.

				//need special check for @ symbol in input expression
				ptr = strchr(s,'@');
				if (ptr == NULL)
				{
					printf("missing variable indicator (@) input expression\n");
					exit(-1);
				}
				if (strlen(s) + strlen(line) > GSTR_MAXSIZE)
				{
					printf("allocating input expression\n");
					s = realloc(s,strlen(s) + strlen(line) + 2);
				}
				//replace @ with line
				//gcc fgets says the newline is two bytes
				//msvd says the newline is one byte
#if defined(_MSC_VER)
				memcpy(ptr + strlen(line)-1,ptr + 1,s + strlen(s) - ptr);
				memcpy(ptr,line,strlen(line)-1);
#else
				memcpy(ptr + strlen(line)-2,ptr + 1,s + strlen(s) - ptr);
				memcpy(ptr,line,strlen(line)-2);
#endif


				printf("=== Starting work on batchfile expression ===\n");
				printf("%s\n",s);
				printf("=============================================\n");

			}
			else
			{
				//this is a normal execution.  make the lines of the batchfile
				//available to the user in N separate user variable, and
				//announce their names.  if the inputs are expressions, will
				//need to evaulate them
	

			}
		}

		if (!do_once)
		{
			//scanf("%[^\n]",s);	//read everything up to carriage return
			fgets(s,GSTR_MAXSIZE,stdin);
			while (1)
			{
				if (s[strlen(s) - 1] == 13 || s[strlen(s) - 1] == 10)
				{
					//replace with a null char and continue
					printf("\n");
					fflush(stdout);
					//get_random_seeds(&g_rand);
					//logprint(logfile,"New random seeds: %u, %u\n\n",g_rand.hi,g_rand.low);
					//srand(g_rand.low);
					s[strlen(s) - 1] = '\0';
					break;
				}
				else
				{
					//last char is not a carriage return means
					//the input is longer than allocated.
					//reallocate and get another chunk
					insize += GSTR_MAXSIZE;
					s = (char *)realloc(s,insize*sizeof(char));
					if (s == NULL)
					{
						printf("couldn't reallocate string when parsing\n");
						exit(-1);
					}
					fgets(s+strlen(s),GSTR_MAXSIZE,stdin);
				}
			}	
		}
		
		/*
		ptr = strstr(s,"algebraic");
		if (ptr != NULL)
		{
			strcpy(str.s,s);
			str.nchars = strlen(s)+1;
			algebraic(&str);
			strcpy(s,"");
			continue;
		}
		*/

		//search for substring help in input
		ptr = strstr(s,"help");
		if (ptr != NULL)
			helpfunc(s);
		else if ((strcmp(s,"quit") == 0) || (strcmp(s,"exit") == 0))
			break;
		else
		{
			//printf("Processing expression: %s\n\n",s);
			logprint(logfile,"Processing expression: %s\n\n",s);
			toStr(s,&str);

			preprocess(&str);
			strcpy(s,str.s);

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
					memcpy(str.s,s+offset,(GSTR_MAXSIZE-offset) * sizeof(char));
					strcpy(s,"ans");
					str.nchars = strlen(str.s) + 1;
				}
				else
				{
					offset = ptr-str.s+1;
					sFree(&str);
					sInit(&str);
					memcpy(str.s,s+offset,(GSTR_MAXSIZE-offset) * sizeof(char));

					s[offset-1] = '\0';
					str.nchars = strlen(str.s) + 1;
				}
			}
			else
				strcpy(s,"ans");

			//look for a trailing semicolon
			if (str.s[str.nchars-2] == ';')
			{
				nooutput = 1;
				str.s[str.nchars-2] = '\0';
				str.nchars--;
			}
			else
				nooutput = 0;

			//if (slog)
			//	logprint(logfile,"%s\n",str.s);

			if (!calc(&str))
			{
				if (strcmp(str.s,"") != 0)
				{
					//printf("calc returned string: %s\n",str.s);
					str2hexz(str.s,&tmp);
					if (set_uvar(s,&tmp))
						new_uvar(s,&tmp);
					if (nooutput == 0)
					{
						if (OBASE == DEC)
						{
							if (VFLAG >= 0)
								printf("\n%s = %s\n\n",s,z2decstr(&tmp,&gstr1));
						}
						else if (OBASE == HEX)
						{
							if (strstr(str.s,"0x") != NULL)
							{
								if (VFLAG >= 0)
									printf("\n%s = %s\n\n",s,str.s);
							}
							else
							{
								z2hexstr(&tmp,&str);
								if (VFLAG >= 0)
									printf("\n%s = %s\n\n",s,str.s);
							}
						}
					}
				}
			}
			
			free_factor_list();
		}

#ifdef WIN32
		fflush(stdin);	//important!  otherwise scanf will process printf's output
		
#else
		if (!do_once)
		{
			//scanf("%s",s);  //read carriage return;
			fflush(stdin);	//important!  otherwise scanf will process printf's output
			fflush(stdout);
		}
#endif

		s = (char *)realloc(s,GSTR_MAXSIZE*sizeof(char));
		if (s == NULL)
		{
			printf("couldn't reallocate string during cleanup\n");
			exit(-1);
		}
		s[0] = '\0';

		if (do_once)
		{
			if (USEBATCHFILE)
			{
				//get the next line, exit if there is nothing more to get
				strcpy(line,"");
				if (feof(batchfile))
					break;

				fgets(line,GSTR_MAXSIZE,batchfile);

				if (strlen(line) == 0)
					break;

				ptr = strchr(line,10);
				if (ptr == NULL)
				{
					printf("warning: newline missing in %s\n",batchfilename);
					toStr(line,&str);
#if defined(_MSC_VER)
					sAppend("\n",&str);
#else
					sAppend("\n\n",&str);
#endif
					strcpy(line,str.s);
				}

				//restore the input expression
				strcpy(s,indup);
			}
			else
			{
				//we are running from the command line, just execute once
				break;
			}
		}
		else
			printf(">> ");
	}

	if (slog)
		fclose(logfile);

	if (USEBATCHFILE)
		fclose(batchfile);

	calc_finalize();
	free(global_factors);
	free(spSOEprimes);
	free(PRIMES);
	sFree(&gstr1);
	sFree(&gstr2);
	sFree(&gstr3);
	sFree(&str);
	free(s);
	free(line);
	free(indup);
	zFree(&zZero);
	zFree(&zOne);
	zFree(&zTwo);
	zFree(&zThree);
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

unsigned process_flags(int argc, char **argv)
{
    int ch = 0, i,j,valid;
	char optbuf[20];
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
		for (j=0; j<numOptions;j++)
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
	else
	{
		printf("invalid option %s\n",opt);
		exit(1);
	}

	zFree(&tmp);
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

		uint64 high_res_time = read_clock();
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

