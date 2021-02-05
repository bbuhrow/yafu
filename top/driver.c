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
#include "ytools.h"
#include "factor.h"
#include "gmp.h"
#include <ecm.h>

#if defined(__unix__)
#include <termios.h>
#endif

#include "cmdOptions.h"

// function to read the .ini file and populate options
void apply_tuneinfo(yafu_obj_t* yobj, fact_obj_t *fobj, char *arg);

// function to print the splash screen to file/screen
void print_splash(info_t* comp_info, int is_cmdline_run, FILE *logfile, 
    int VFLAG, double freq, int numwit);

// functions to make a batchfile ready to execute, and to process batchfile lines
void prepare_batchfile(char *input_exp);
char * process_batchline(yafu_obj_t* yobj, char *input_exp, char *indup, int *code);
void finalize_batchline(yafu_obj_t* yobj);
int exp_is_open(char *line, int firstline);

// functions to process all incoming arguments
int check_expression(options_t *options);
char * get_input(char *input_exp, uint32 *insize);

#if defined(__unix__)
#define CMDHIST_SIZE 16
static char **CMDHIST;
static int CMDHIST_HEAD = 0;
static int CMDHIST_TAIL = 0;
#endif

int main(int argc, char *argv[])
{
	uint32 insize = GSTR_MAXSIZE;
	char *input_exp, *ptr, *indup, *input_line;
    str_t input_str;
    
	int slog,is_cmdline_run=0;
	FILE *logfile;
    FILE *scriptfile = NULL;
	fact_obj_t *fobj;
    int firstline = 1;
    options_t *options;
    meta_t calc_metadata;
    yafu_obj_t yafu_obj;
    soe_staticdata_t *sdata;
    info_t comp_info;
    int i;

#if defined(__unix__)

    static struct termios oldtio, newtio;
    tcgetattr(0, &oldtio);
    newtio = oldtio;
    newtio.c_lflag &= ~(ICANON | ECHO | ECHOE | ECHOK);

    CMDHIST = (char **)malloc(CMDHIST_SIZE * sizeof(char *));
    for (i = 0; i < CMDHIST_SIZE; i++)
    {
        CMDHIST[i] = (char *)malloc(GSTR_MAXSIZE*sizeof(char));
    }

#endif

	//the input expression
	input_exp = (char *)malloc(GSTR_MAXSIZE*sizeof(char));
	indup = (char *)malloc(GSTR_MAXSIZE*sizeof(char));
    input_line = (char *)malloc(GSTR_MAXSIZE*sizeof(char));
    sInit(&input_str);
	strcpy(input_exp,"");
    strcpy(input_line, "");
	
	// set defaults for various things and read the .ini file, if any.
	yafu_init(&yafu_obj);
    options = initOpt();
    readINI("yafu.ini", options);

    // then process the command line, overriding any .ini settings.
    processOpts(argc, argv, options);

    // some things go into globals, but this is being phased out
    yafu_obj.VFLAG = options->verbosity;
    yafu_obj.THREADS = options->threads;
    yafu_obj.CMD_LINE_REPEAT = options->repeat;
    yafu_obj.VERBOSE_PROC_INFO = options->vproc;
    yafu_obj.LATHREADS = options->lathreads;
    if (strlen(options->factorlog) == 0)
    {
        yafu_obj.LOGFLAG = 0;
    }
    if (strlen(options->batchfile) > 0)
    {
        yafu_obj.USEBATCHFILE = 1;
    }
    if (options->rand_seed == 0)
    {
        uint32 seed1, seed2;
        get_random_seeds(&seed1, &seed2);
        options->rand_seed = ((uint64)seed2 << 32) | (uint64)seed1;
    }
    else
    {
        yafu_obj.USERSEED = 1;
    }

    sdata = soe_init(options->verbosity, options->threads, options->soe_blocksize);

    // initial cache of primes.
    szSOEp = 10000000; 
    PRIMES = soe_wrapper(sdata, 0, szSOEp, 0, &NUM_P, 0, 0);
    szSOEp = NUM_P;

    // save a batch of sieve primes too.
    spSOEprimes = (uint32*)malloc((size_t)(NUM_P * sizeof(uint32)));
    for (i = 0; i < NUM_P; i++)
    {
        spSOEprimes[i] = (uint32)PRIMES[i];
    }

    P_MIN = 0;
    P_MAX = PRIMES[(uint32)NUM_P - 1];


#if !defined( TARGET_KNC )
    // get the computer name, cache sizes, etc.  store in globals
    // we need to have the cpu id string before calling apply_tuneinfo so that
    // any tune_info lines are applied correctly.
    get_computer_info(&comp_info, 0);
#endif

#if !defined( TARGET_KNC )
    // now that we've processed arguments, spit out vproc info if requested
#ifndef __APPLE__
    // 
    if (options->vproc)
    {
        extended_cpuid(comp_info.idstr, &comp_info.cachelinesize,
            &comp_info.bSSE41Extensions, &comp_info.AVX,
            &comp_info.AVX2, options->vproc);
    }
#endif
#endif

	// a factorization object that gets passed around to any factorization routine
	// called out in the input expression.  if no factorization routine is specified,
	// this is not used.  initialize and pass in all of the options.
	fobj = (fact_obj_t *)malloc(sizeof(fact_obj_t));
	init_factobj(fobj, options);
    if (strlen(options->tune_info) > 0)
    {
        printf("tune_info: %s\n", options->tune_info);
        apply_tuneinfo(&yafu_obj, fobj, options->tune_info);
    }
    
    // figure out cpu freq in order to scale qs time estimations
    // 0.1 seconds won't be very accurate, but hopefully close
    // enough for the rough scaling we'll be doing anyway.
    if (options->no_clk_test == 0)
        fobj->MEAS_CPU_FREQUENCY = measure_processor_speed() / 1.0e5;
    else
        fobj->MEAS_CPU_FREQUENCY = 42;

    strcpy(fobj->CPU_ID_STR, comp_info.idstr);
    fobj->HAS_AVX2 = comp_info.AVX2;
    fobj->HAS_AVX = comp_info.AVX;
    fobj->HAS_SSE41 = comp_info.bSSE41Extensions;
    fobj->NUM_WITNESSES = options->num_prp_witnesses;
    fobj->cache_size1 = comp_info.L1cache;
    fobj->cache_size2 = comp_info.L2cache;
    fobj->LOGFLAG = yafu_obj.LOGFLAG;
    fobj->THREADS = yafu_obj.THREADS;

#if BITS_PER_DIGIT == 64
    fobj->lcg_state = options->rand_seed;
#else
    fobj->lcg_state = (uint32)options->rand_seed;
#endif	

    calc_metadata.fobj = fobj;
    calc_metadata.options = options;
    calc_metadata.sdata = sdata;

	// check/process input arguments
	is_cmdline_run = check_expression(options);
    if (is_cmdline_run == 3)
    {
        // a default function applied to text that has no other function.
        int len = strlen(options->inputExpr) + 9;
        options->inputExpr = (char*)xrealloc(options->inputExpr, len);
        input_exp = (char*)xrealloc(input_exp, len);
        sprintf(input_exp, "factor(%s)", options->inputExpr);
        strcpy(options->inputExpr, input_exp);
        strcpy(input_line, options->inputExpr);
        is_cmdline_run = 1;
    }
    else
    {
        strcpy(input_line, options->inputExpr);
        strcpy(input_exp, options->inputExpr);
    }

	if (is_cmdline_run == 2)
	{
		// batchfile from stdin
        yafu_obj.USEBATCHFILE = 2;
	}

	// get the batchfile ready, if requested
	if (yafu_obj.USEBATCHFILE)
	{
		prepare_batchfile(input_exp);		
		
		// batchfile jobs are command line in nature
		is_cmdline_run = 1;		
	}

    if (strlen(yafu_obj.scriptname) > 0)
    {
        scriptfile = fopen(yafu_obj.scriptname, "r");
        if (scriptfile == NULL)
        {
            printf("could not find %s\n", yafu_obj.scriptname);
            exit(1);
        }
        else
        {
            is_cmdline_run = 1;
        }
    }

    if (yafu_obj.USEBATCHFILE || (yafu_obj.CMD_LINE_REPEAT > 0))
    {
        strcpy(indup, input_exp);	//remember the input expression
    }

	//never run silently when run interactively, else the results of
	//calculations will never be displayed.
    if (!is_cmdline_run && yafu_obj.VFLAG < 0)
    {
        yafu_obj.VFLAG = 0;
    }

	//session log
    if (yafu_obj.LOGFLAG)
    {
        logfile = fopen(yafu_obj.sessionname, "a");
        if (logfile == NULL)
        {
            printf("fopen error: %s\n", strerror(errno));
            printf("couldn't open %s for appending\n", yafu_obj.sessionname);
            slog = 0;
        }
        else
            slog = 1;
    }
    else
    {
        logfile = NULL;
        slog = 0;
    }
		
	// print the splash screen, to the logfile and depending on options, to the screen
	print_splash(&comp_info, is_cmdline_run, logfile, yafu_obj.VFLAG, 
        yafu_obj.MEAS_CPU_FREQUENCY, yafu_obj.NUM_WITNESSES);
	
	// start the calculator
	// right now this just allocates room for user variables
	calc_init(options->rand_seed);
		
	logprint(logfile,"Random seed: %" PRIu64 "\n", options->rand_seed);
	fflush(logfile);

	//printf("WARNING: constant seed is set\n");
	//g_rand.hi = 123;
	//g_rand.low = 123;
    srand((unsigned int)options->rand_seed);

	// command line
	while (1)
	{		
        // running interactively, reset the fobj every line.
        reset_factobj(fobj);

		// handle a batch file, if passed in.
		if (yafu_obj.USEBATCHFILE)
		{
			int code;
            input_line = process_batchline(&yafu_obj, input_line, indup, &code);
			if (code == 1)
			{
				finalize_batchline(&yafu_obj);
				break;
			}
			else if (code == 2)
				continue;
		}
        else if (strlen(yafu_obj.scriptname) > 0)
        {
            if (scriptfile != NULL)
            {
                if (fgets(input_line, GSTR_MAXSIZE, scriptfile) == NULL)
                {
                    //    break;
                }
                
            }
        }
		else if (!is_cmdline_run)
		{
#if defined(__unix__)
            tcsetattr(0, TCSANOW, &newtio);
#endif
            input_line = get_input(input_line, &insize);
#if defined(__unix__)
            tcsetattr(0, TCSANOW, &oldtio);
#endif
		}
		else
		{
			// input expression already read in.  nothing to do here.

		}
		
		// help, exit, or execute the current expression...
        ptr = strstr(input_line, "help");
		if (ptr != NULL)
            helpfunc(input_line);
        else if ((strcmp(input_line, "quit") == 0) || (strcmp(input_line, "exit") == 0))
			break;
		else
		{
            char* result;
            sAppend(input_line, &input_str);
            if (exp_is_open(input_line, firstline))
            {
                if (strlen(input_line) > 0)
                    sAppend(",", &input_str);
                firstline = 0;
                continue;
            }

            firstline = 1;
            reset_preprocessor();
            logprint(logfile, "Processing: %s\n", input_str.s);
            result = process_expression(input_str.s, &calc_metadata, 0, 0);
            logprint(logfile, "Result    : %s\n", result);
            sClear(&input_str);
            if (result != NULL)
            {
                free(result);
            }
		}

#if defined(WIN32) && !defined(__MINGW32__)
		fflush(stdin);	//important!  otherwise scanf will process printf's output
		
#else
		if (!is_cmdline_run)
		{
			fflush(stdin);	//important!  otherwise scanf will process printf's output
			fflush(stdout);
		}
#endif

		// get the next expression, if running a batchfile, or
		// re-display the command prompt
		if (yafu_obj.CMD_LINE_REPEAT == 0)
		{
            input_line = (char *)realloc(input_line, GSTR_MAXSIZE*sizeof(char));
            if (input_line == NULL)
			{
				printf("couldn't reallocate string during cleanup\n");
				exit(-1);
			}

            input_line[0] = '\0';
		}

		if (is_cmdline_run)
		{
			if (yafu_obj.USEBATCHFILE)
			{
				// the line from the batchfile finished.  make the temporary file
				// created in processs_batchline the new batchfile, with the line
				// we just finished removed.
				finalize_batchline(&yafu_obj);
			}
            else if (scriptfile != NULL)
            {
                if (feof(scriptfile))
                {
                    if (yafu_obj.CMD_LINE_REPEAT > 0)
                    {
                        yafu_obj.CMD_LINE_REPEAT--;
                        fclose(scriptfile);
                        scriptfile = fopen(yafu_obj.scriptname, "r");
                        if (scriptfile == NULL)
                        {
                            printf("could not find %s\n", yafu_obj.scriptname);
                            exit(1);
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }
			else if (yafu_obj.CMD_LINE_REPEAT > 0)
			{
                yafu_obj.CMD_LINE_REPEAT--;
                strcpy(input_line, indup);
			}
			else
				break;
		}
        else
        {
            printf(">> ");
            fflush(stdout);
        }

	}

    if (scriptfile != NULL)
    {
        fclose(scriptfile);
    }

    if (slog)
    {
        fclose(logfile);
    }

	calc_finalize();
	yafu_finalize(&yafu_obj);	
	free(input_exp);
    free(input_line);
	free(indup);	
	free_factobj(fobj);
	free(fobj);      
    sFree(&input_str);
    free(options->inputExpr);
    soe_finalize(sdata);
    free(sdata);

#if defined(__unix__)
    for (i = 0; i < CMDHIST_SIZE; i++)
    {
        free(CMDHIST[i]);
    }
    free(CMDHIST);
#endif

	return 0;
}

int exp_is_open(char *line, int firstline)
{
    int i;
    static int openp, closedp, openb, closedb;

    if (firstline)
    {
        openp = openb = closedp = closedb = 0;
    }

    for (i = 0; i < strlen(line); i++)
    {
        if (line[i] == '(') openp++;
        if (line[i] == ')') closedp++;
        if (line[i] == '{') openb++;
        if (line[i] == '}') closedb++;
    }
    if ((openp == closedp) && (openb == closedb))
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

char * get_input(char *input_exp, uint32 *insize)
{
#if !defined(__unix__)

    // get command from user
    fgets(input_exp, GSTR_MAXSIZE, stdin);

    while (1)
    {
        if (input_exp[strlen(input_exp) - 1] == 13 || input_exp[strlen(input_exp) - 1] == 10)
        {
            // replace with a null char and continue
            printf("\n");
            fflush(stdout);
            input_exp[strlen(input_exp) - 1] = '\0';
            break;
        }
        else
        {
            // last char is not a carriage return means
            // the input is longer than allocated.
            // reallocate and get another chunk
            *insize += GSTR_MAXSIZE;
            input_exp = (char *)realloc(input_exp, *insize * sizeof(char));
            if (input_exp == NULL)
            {
                printf("couldn't reallocate string when parsing\n");
                exit(-1);
            }
            fgets(input_exp + strlen(input_exp), GSTR_MAXSIZE, stdin);
        }
    }

#else
    
    int n = 0;
    int p = CMDHIST_HEAD;
    strcpy(CMDHIST[p], "");

    while (1)
    {        
        int c = getc(stdin);

        // is this an escape sequence?
        if (c == 27) {
            // "throw away" next two characters which specify escape sequence
            int c1 = getc(stdin);
            int c2 = getc(stdin);
            int i;

            if ((c1 == 91) && (c2 == 65))
            {
                // clear the current screen contents
                for (i = 0; i < n; i++)
                    printf("\b");

                for (i = 0; i < n; i++)
                    printf(" ");

                for (i = 0; i < n; i++)
                    printf("\b");

                // save whatever is currently entered                
                if (p == CMDHIST_HEAD)
                {
                    input_exp[n] = '\0';
                    memcpy(CMDHIST[CMDHIST_HEAD], input_exp, GSTR_MAXSIZE * sizeof(char));
                }

                // uparrow     
                if (CMDHIST_HEAD >= CMDHIST_TAIL)
                {
                    p--;

                    // wrap
                    if (p < 0)
                    {
                        // but not past tail
                        p = 0;
                    }
                }
                else
                {
                    p--;

                    // wrap
                    if (p < 0)
                    {
                        p = CMDHIST_SIZE - 1;
                    }
                }                                

                // and print the previous one
                printf("%s", CMDHIST[p]);
                strcpy(input_exp, CMDHIST[p]);
                n = strlen(input_exp);
            }
            else if ((c1 == 91) && (c2 == 66))
            {
                // downarrow
                // clear the current screen contents
                for (i = 0; i < strlen(CMDHIST[p]); i++)
                    printf("\b");

                for (i = 0; i < strlen(CMDHIST[p]); i++)
                    printf(" ");

                for (i = 0; i < strlen(CMDHIST[p]); i++)
                    printf("\b");

                if (p != CMDHIST_HEAD)
                {
                    p++;

                    // wrap
                    if (p == CMDHIST_SIZE)
                    {
                        p = 0;
                    }
                }                            

                // and print the next one
                printf("%s", CMDHIST[p]);
                strcpy(input_exp, CMDHIST[p]);
                n = strlen(input_exp);
            }
            else if ((c1 == 91) && (c2 == 67))
            {
                // rightarrow
            }
            else if ((c1 == 91) && (c2 == 68))
            {
                // leftarrow
            }
            else
            {
                printf("unknown escape sequence %d %d\n", c1, c2);
            }

            continue;
        }

        // if backspace
        if (c == 0x7f)
        {
            //fprintf(stderr,"saw a backspace\n"); fflush(stderr);
            if (n > 0)
            {
                // go one char left
                printf("\b");
                // overwrite the char with whitespace
                printf(" ");
                // go back to "now removed char position"
                printf("\b");
                n--;
            }
            continue;
        }

        if (c == EOF)
        {
            printf("\n");
            exit(0);
        }

        if ((c == 13) || (c == 10))
        {
            input_exp[n++] = '\0';
            break;
        }

        putc(c, stdout);
        input_exp[n++] = (char)c;

        if (n >= *insize)
        {
            *insize += GSTR_MAXSIZE;
            input_exp = (char *)realloc(input_exp, *insize * sizeof(char));
            if (input_exp == NULL)
            {
                printf("couldn't reallocate string when parsing\n");
                exit(-1);
            }
        }
    }

    printf("\n");
    fflush(stdout);

    if (strlen(input_exp) > 0)
    {
        memcpy(CMDHIST[CMDHIST_HEAD++], input_exp, GSTR_MAXSIZE * sizeof(char));

        if (CMDHIST_TAIL > 0)
        {
            CMDHIST_TAIL++;
            if (CMDHIST_TAIL == CMDHIST_SIZE)
                CMDHIST_TAIL = 0;
        }

        if (CMDHIST_HEAD == CMDHIST_SIZE)
        {
            CMDHIST_HEAD = 0;
            if (CMDHIST_TAIL == 0)
                CMDHIST_TAIL = 1;
        }
    }

#endif

    return input_exp;
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
		while (isspace((int)s[j])) j++;		//skip white space
	func = s + j;

	//func now points to a string with the desired help topic
	//open the doc file and search matching topics
	doc = fopen("docfile.txt","r");
	if (doc == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
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

void prepare_batchfile(char *input_exp)
{
	char *ptr;
	
	//look for @ symbol in input expression
	ptr = strchr(input_exp,'@');
	if (ptr == NULL)
	{
		printf("no variable indicator (@): interpreting batchfile lines as input expressions\n");
		sprintf(input_exp,"@");
	}

	return;
}

int check_expression(options_t* options)
{
    int is_cmdline_run = 0;

    if (strlen(options->inputExpr) > 0)
    {
        // there was an expression on the command line.
        is_cmdline_run = 1;
    }

    // now check for incoming pipes or redirects.  If we see one, ignore the
    // command line expression and process the pipe/redirect.

    // Different ways to detect depending on the environment.  Here's an attempt
    // to cover them; msys2 is difficult.
    // if running in a MSYS2 fake console:
    // from https://github.com/nodejs/node/issues/3006
    // When running node in a fake console, it's useful to imagine that you're running it 
    // with input and output redirected to files(e.g.node foobar.js <infile.txt >outfile.txt).
    // That's basically what node thinks is going on (a pipe and a file look more-or-less 
    // the same to node). Things that would work with file redirections 
    // will work in the fake console.
    //
    // not sure how to sort it out in the case of msys2.  recommended running in normal
    // windows cmd terminal once it is built.
    // detect if stdin is a pipe
    // http://stackoverflow.com/questions/1312922/detect-if-stdin-is-a-terminal-or-pipe-in-c-c-qt
    // But this doesn't work if we are running in a msys console because of how
    // they interface with stdin/out/err through pipes, so there will always
    // be a pipe.
    // https://github.com/msys2/msys2/wiki/Porting
#if defined(__MINGW32__)
    // I'm not sure how to detect at runtime if this is an msys shell.
    // So unfortunately if we compile with mingw we basically have to remove 
    // the ability to process from pipes or redirects.  should be able to use 
    // batchfiles via command line switch still.
    if (0)
    {

#elif defined(WIN32) 
    if (_isatty(_fileno(stdin)) == 0)
    {
        fseek(stdin, -1, SEEK_END);
        if (ftell(stdin) >= 0)
        {
            rewind(stdin);
#else
        if (isatty(fileno(stdin)) == 0)
        {
#endif

            // ok, we also have incoming data.  This is just
            // batchfile mode with the batchfile = stdin.
            is_cmdline_run = 2;
        }
#if defined(WIN32) && !defined(__MINGW32__)		//not complete, but ok for now
    }
#endif
    else
    {
        // no incoming data, just execute the provided expression, if any,
        // or start up an interactive session.
        // special check: if there is no function call, insert a 
        // default function call
        if ((is_cmdline_run == 1) && (strstr(options->inputExpr, "(") == NULL))
        {
            // this indicates we have an input expression with no function call
            is_cmdline_run = 3;
        }
    }

#if defined(__MINGW32__)
    }
#endif

	return is_cmdline_run;

}

void print_splash(info_t *comp_info, int is_cmdline_run, FILE *logfile, 
    int VFLAG, double freq, int numwit)
{
	if (VFLAG >= 0)
		printf("\n\n");

	if ((VFLAG > 0) || !is_cmdline_run)
	{	
		logprint(NULL,"System/Build Info: \n");
        fflush(stdout);
	}
    logprint(logfile, "=====================================\n");
	logprint(logfile,"System/Build Info: \n");

    if ((VFLAG > 0) || !is_cmdline_run)
    {
#ifdef _MSC_MPIR_VERSION
#ifdef ECM_VERSION
        printf("Using GMP-ECM %s, Powered by MPIR %s\n", ECM_VERSION,
            _MSC_MPIR_VERSION);
        logprint(logfile, "Using GMP-ECM %s, Powered by MPIR %s\n", ECM_VERSION,
            _MSC_MPIR_VERSION);
#elif defined(VERSION)

        printf("Using GMP-ECM %s, Powered by MPIR %s\n", VERSION,
            _MSC_MPIR_VERSION);
        logprint(logfile,"Using GMP-ECM %s, Powered by MPIR %s\n", VERSION,
            _MSC_MPIR_VERSION);

#else
        printf("Using GMP-ECM <unknown>, Powered by MPIR %s\n", 
            _MSC_MPIR_VERSION);
        logprint(logfile,"Using GMP-ECM <unknown>, Powered by MPIR %s\n",
            _MSC_MPIR_VERSION);


#endif
#else
#ifdef ECM_VERSION
        printf("Using GMP-ECM %s, Powered by GMP %d.%d.%d\n", ECM_VERSION, 
            __GNU_MP_VERSION,__GNU_MP_VERSION_MINOR,__GNU_MP_VERSION_PATCHLEVEL);
        logprint(logfile,"Using GMP-ECM %s, Powered by GMP %d.%d.%d\n", ECM_VERSION,
            __GNU_MP_VERSION,__GNU_MP_VERSION_MINOR,__GNU_MP_VERSION_PATCHLEVEL);
#else
        printf("Using GMP-ECM, Powered by GMP\n");
        logprint(logfile,"Using GMP-ECM, Powered by GMP\n");
#endif

#endif

        fflush(stdout);
    }

    logprint(logfile,"cached %u primes. pmax = %u\n", szSOEp, spSOEprimes[szSOEp-1]);
    logprint(logfile,"detected %s\ndetected L1 = %d bytes, L2 = %d bytes, CL = %d bytes\n",
		comp_info->idstr, comp_info->L1cache, comp_info->L2cache, comp_info->cachelinesize);
    if (freq > 100.0)
        logprint(logfile,"measured cpu frequency ~= %f\n", freq);
    if (numwit == 1)
        logprint(logfile,"using %u random witness for Rabin-Miller PRP checks\n", numwit);
    else
        logprint(logfile, "using %u random witnesses for Rabin-Miller PRP checks\n", numwit);

	if (VFLAG > 0 || !is_cmdline_run)
	{		
		printf("detected %s\ndetected L1 = %d bytes, L2 = %d bytes, CL = %d bytes\n",
            comp_info->idstr, comp_info->L1cache, comp_info->L2cache, comp_info->cachelinesize);
        if (freq > 100.0)
		    printf("measured cpu frequency ~= %f\n", freq);
        if (numwit == 1)
		    printf("using %u random witness for Rabin-Miller PRP checks\n\n", numwit);
        else
            printf("using %u random witnesses for Rabin-Miller PRP checks\n\n", numwit);

		printf("===============================================================\n");
		printf("======= Welcome to YAFU (Yet Another Factoring Utility) =======\n");
		printf("=======             bbuhrow@gmail.com                   =======\n");
		printf("=======     Type help at any time, or quit to quit      =======\n");
		printf("===============================================================\n");
		printf("cached %u primes. pmax = %u\n\n", 
            szSOEp, spSOEprimes[szSOEp-1]);
		printf("\n>> ");
        fflush(stdout);
	}

	return;
}

void yafu_set_idle_priority(void) {

#if defined(WIN32) || defined(_WIN64)
	SetPriorityClass(GetCurrentProcess(),
			IDLE_PRIORITY_CLASS);
#else
	nice(100);
#endif
}

void yafu_init(yafu_obj_t* yobj)
{
    yobj->VFLAG = 0;
    yobj->VERBOSE_PROC_INFO = 0;
    yobj->LOGFLAG = 1;
    yobj->NUM_WITNESSES = 1;
    yobj->NO_CLK_TEST = 0;
    yobj->USEBATCHFILE = 0;
    yobj->USERSEED = 0;
    yobj->THREADS = 1;
    yobj->LATHREADS = 0;
    yobj->CMD_LINE_REPEAT = 0;
    yobj->MEAS_CPU_FREQUENCY = 0.0;

	strcpy(yobj->sessionname,"session.log");
    strcpy(yobj->scriptname, "");

	return;
}

void yafu_finalize(yafu_obj_t* yobj)
{
	free(spSOEprimes);
	free(PRIMES);

	return;
}

void finalize_batchline(yafu_obj_t* yobj)
{
	if (yobj->USEBATCHFILE == 1)
	{
		rename(yobj->batchfilename,"_bkup");
		rename("__tmpbatchfile", yobj->batchfilename);
		remove("_bkup");
		remove("__tmpbatchfile");
	}

	return;
}

char * process_batchline(yafu_obj_t* yobj, char *input_exp, char *indup, int *code)
{
	int nChars, j, i;
	char *line, tmpline[GSTR_MAXSIZE], *ptr, *ptr2;
	FILE *batchfile, *tmpfile;

	//try to open the file
	if (yobj->USEBATCHFILE == 2)
		batchfile = stdin;
	else
		batchfile = fopen(yobj->batchfilename,"r");

	if (batchfile == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("couldn't open %s for reading\n", yobj->batchfilename);
		exit(-1);
	}	

	//load the next line of the batch file and get the expression
	//ready for processing
	line = (char *)malloc(GSTR_MAXSIZE * sizeof(char));
	strcpy(line,"");
	strcpy(input_exp,"");

	// read a line - skipping blank lines
	do
	{
		while (1)
		{
			ptr = fgets(tmpline,GSTR_MAXSIZE,batchfile);
			strcpy(line + strlen(line), tmpline);
			
			// stop if we didn't read anything
			if (feof(batchfile))
			{
				printf("eof; done processing batchfile\n");
				fclose(batchfile);
				*code = 1;
				free(line);
				return input_exp;
			}

			if (ptr == NULL)
			{
				printf("fgets returned null; done processing batchfile\n");		
				fclose(batchfile);
				*code = 1;
				free(line);
				return input_exp;
			}

			// if we got the end of the line, stop reading
			if ((line[strlen(line)-1] == 0xa) ||
				(line[strlen(line)-1] == 0xd))
				break;

			// else reallocate the buffer and get some more
			line = (char *)realloc(line, (strlen(line) + GSTR_MAXSIZE) * sizeof(char));
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

	} while (strlen(line) == 0);	

	// this only applies for non-stdin batchfiles
	if (yobj->USEBATCHFILE == 1)
	{
		// copy everything in the file after the line we just read to
		// a temporary file.  if the expression we just read finishes, 
		// the temporary file will become the batch file (effectively 
		// eliminating the expression from the batch file).
		tmpfile = fopen("__tmpbatchfile", "w");
	
		if (tmpfile == NULL)
		{
			printf("fopen error: %s\n", strerror(errno));
			printf("couldn't open __tmpbatchfile for reading\n");
			exit(-1);
		}	

		while (!feof(batchfile))
		{
			ptr2 = fgets(tmpline,GSTR_MAXSIZE,batchfile);
			if (ptr2 == NULL)
				break;

			if (strlen(tmpline) == 0)
				continue;

			fputs(tmpline,tmpfile);
		}
		fclose(tmpfile);

		// close the batchfile
		fclose(batchfile);		

	}

	//ignore blank lines
	if (strlen(line) == 0)
	{
		*code = 2;
		free(line);
		return input_exp;
	}

	//ignore comment lines
	if (((line[0] == '/') && (line[1] == '/')) || (line[0] == '%'))
	{
		*code = 2;
		free(line);
		return input_exp;
	}

	//substitute the batchfile line into the '@' symbol in the input expression
	nChars = 0;
	if ((strlen(indup) + strlen(line)) >= strlen(input_exp))
		input_exp = (char *)realloc(input_exp, strlen(indup) + strlen(line) + 2);

	for (i=0; i<strlen(indup); i++)
	{
		if (indup[i] == '@')
		{
			for (j=0; j<strlen(line); j++)
				input_exp[nChars++] = line[j];
		}
		else				
			input_exp[nChars++] = indup[i];
	}
	input_exp[nChars++] = '\0';

	if (yobj->VFLAG >= 0)
	{
		printf("=== Starting work on batchfile expression ===\n");
		printf("%s\n",input_exp);
		printf("=============================================\n");
		fflush(stdout);
	}

	free(line);
	*code = 0;
	return input_exp;;
}

void apply_tuneinfo(yafu_obj_t* yobj, fact_obj_t *fobj, char *arg)
{
	int i,j;
	char cpustr[80], osstr[80];
    int xover;

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

	//printf("found OS = %s and CPU = %s in tune_info field, this cpu is %s\n",
    //    osstr, cpustr, CPU_ID_STR);

    // "xover" trumps tune info.  I.e., if a specific crossover has been
    // specified, prefer this to whatever may be in the tune_info string.
    xover = fobj->autofact_obj.qs_gnfs_xover;

#if defined(_WIN64)
	if ((strcmp(cpustr, yobj->CPU_ID_STR) == 0) && (strcmp(osstr, "WIN64") == 0))
	{
		//printf("Applying tune_info entry for %s - %s\n",osstr,cpustr);
		
		sscanf(arg + i + 1, "%lg, %lg, %lg, %lg, %lg, %lg",
			&fobj->qs_obj.qs_multiplier, &fobj->qs_obj.qs_exponent,
			&fobj->nfs_obj.gnfs_multiplier, &fobj->nfs_obj.gnfs_exponent, 
			&fobj->autofact_obj.qs_gnfs_xover, &fobj->nfs_obj.gnfs_tune_freq);
		fobj->qs_obj.qs_tune_freq = fobj->nfs_obj.gnfs_tune_freq;

	}
#elif defined(WIN32)
	if ((strcmp(cpustr, yobj->CPU_ID_STR) == 0) && (strcmp(osstr, "WIN32") == 0))
	{
		//printf("Applying tune_info entry for %s - %s\n",osstr,cpustr);
		sscanf(arg + i + 1, "%lg, %lg, %lg, %lg, %lg, %lg",
			&fobj->qs_obj.qs_multiplier, &fobj->qs_obj.qs_exponent,
			&fobj->nfs_obj.gnfs_multiplier, &fobj->nfs_obj.gnfs_exponent, 
			&fobj->autofact_obj.qs_gnfs_xover, &fobj->nfs_obj.gnfs_tune_freq);
		fobj->qs_obj.qs_tune_freq = fobj->nfs_obj.gnfs_tune_freq;
	}
#elif BITS_PER_DIGIT == 64
	if ((strcmp(cpustr, yobj->CPU_ID_STR) == 0) && (strcmp(osstr, "LINUX64") == 0))
	{
		//printf("Applying tune_info entry for %s - %s\n",osstr,cpustr);
		
		sscanf(arg + i + 1, "%lg, %lg, %lg, %lg, %lg, %lg",
			&fobj->qs_obj.qs_multiplier, &fobj->qs_obj.qs_exponent,
			&fobj->nfs_obj.gnfs_multiplier, &fobj->nfs_obj.gnfs_exponent, 
			&fobj->autofact_obj.qs_gnfs_xover, &fobj->nfs_obj.gnfs_tune_freq);
		fobj->qs_obj.qs_tune_freq = fobj->nfs_obj.gnfs_tune_freq;
	}
#else 
	if ((strcmp(cpustr, yobj->CPU_ID_STR) == 0) && (strcmp(osstr, "LINUX32") == 0))
	{
		//printf("Applying tune_info entry for %s - %s\n",osstr,cpustr);
		
		sscanf(arg + i + 1, "%lg, %lg, %lg, %lg, %lg, %lg",
			&fobj->qs_obj.qs_multiplier, &fobj->qs_obj.qs_exponent,
			&fobj->nfs_obj.gnfs_multiplier, &fobj->nfs_obj.gnfs_exponent, 
			&fobj->autofact_obj.qs_gnfs_xover, &fobj->nfs_obj.gnfs_tune_freq);
		fobj->qs_obj.qs_tune_freq = fobj->nfs_obj.gnfs_tune_freq;
	}
#endif	

	//printf("QS_MULTIPLIER = %lg, QS_EXPONENT = %lg\nNFS_MULTIPLIER = %lg, NFS_EXPONENT = %lg\nXOVER = %lg, TUNE_FREQ = %lg\n",
	//	fobj->qs_obj.qs_multiplier, fobj->qs_obj.qs_exponent,
	//	fobj->nfs_obj.gnfs_multiplier, fobj->nfs_obj.gnfs_exponent, 
	//	fobj->autofact_obj.qs_gnfs_xover, fobj->qs_obj.qs_tune_freq);

    // restore the user's xover if preferred.
    if (fobj->autofact_obj.prefer_xover)
    {
        fobj->autofact_obj.qs_gnfs_xover = xover;
    }

	return;
}

