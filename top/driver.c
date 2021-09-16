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
#include "autofactor.h"
#include "gmp.h"
#include <ecm.h>
#include <immintrin.h>

#if defined(__unix__)
#include <termios.h>
#elif defined(_MSC_VER)
#include <io.h>
#endif

#include "cmdOptions.h"

// function to read the .ini file and populate options
void apply_tuneinfo(yafu_obj_t* yobj, fact_obj_t *fobj, char *arg);

// function to print the splash screen to file/screen
void print_splash(fact_obj_t* fobj, info_t* comp_info, int is_cmdline_run,
    FILE* logfile, int VFLAG, double freq, int numwit);
void helpfunc(char* s);

// functions to make a batchfile ready to execute, and to process batchfile lines
void prepare_batchfile(char *input_exp);
char * process_batchline(yafu_obj_t* yobj, char *input_exp, char *indup, int *code);
void finalize_batchline(yafu_obj_t* yobj);
int exp_is_open(char *line, int firstline);

// functions to process all incoming arguments
int check_expression(options_t *options);
char * get_input(char *input_exp, uint32_t *insize);
void options_to_factobj(fact_obj_t* fobj, options_t* options);

#if defined(__unix__)
#define CMDHIST_SIZE 16
static char **CMDHIST;
static int CMDHIST_HEAD = 0;
static int CMDHIST_TAIL = 0;
#endif

int main(int argc, char *argv[])
{
	uint32_t insize = GSTR_MAXSIZE;
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
    //soe_staticdata_t* sdata;
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
        uint32_t seed1, seed2;
        get_random_seeds(&seed1, &seed2);
        options->rand_seed = ((uint64_t)seed2 << 32) | (uint64_t)seed1;
    }
    else
    {
        printf("input user seed %" PRIu64" detected\n", options->rand_seed);
        yafu_obj.USERSEED = 1;
    }

#if !defined(__APPLE__)
    // get the computer name, cache sizes, etc.  store in globals
    // we need to have the cpu id string before calling apply_tuneinfo so that
    // any tune_info lines are applied correctly.
    ytools_get_computer_info(&comp_info, options->vproc);
    strncpy(yafu_obj.CPU_ID_STR, comp_info.idstr, 255);
#endif

	// a factorization object that gets passed around to any factorization routine
	// called out in the input expression.  if no factorization routine is specified,
	// this is not used.  initialize and pass in all of the options.
	fobj = (fact_obj_t *)malloc(sizeof(fact_obj_t));
	init_factobj(fobj);
    
    options_to_factobj(fobj, options);
    for (i = 0; i < options->num_tune_info; i++)
    {
        apply_tuneinfo(&yafu_obj, fobj, options->tune_info[i]);
    }

    fobj->MEAS_CPU_FREQUENCY = 42;  // not used anymore
    strcpy(fobj->CPU_ID_STR, comp_info.idstr);
    fobj->HAS_AVX2 = comp_info.AVX2;
    fobj->HAS_AVX = comp_info.AVX;
    fobj->HAS_SSE41 = comp_info.bSSE41Extensions;
    fobj->NUM_WITNESSES = yafu_obj.NUM_WITNESSES = options->num_prp_witnesses;
    fobj->cache_size1 = fobj->L1CACHE = comp_info.L1cache;
    fobj->cache_size2 = fobj->L2CACHE = comp_info.L2cache;
    fobj->LOGFLAG = yafu_obj.LOGFLAG;
    fobj->THREADS = yafu_obj.THREADS;
    fobj->HAS_BMI2 = comp_info.BMI2;

#if defined(__INTEL_COMPILER)
    if (_may_i_use_cpu_feature(_FEATURE_AVX512F))
#elif defined(__GNUC__)
    if (__builtin_cpu_supports("avx512f"))
#else
    if (1)
#endif
    {
        fobj->HAS_AVX512F = comp_info.AVX512F;
    }

#ifdef __INTEL_COMPILER
    if (_may_i_use_cpu_feature(_FEATURE_AVX512BW))
#elif defined(__GNUC__)
    if (__builtin_cpu_supports("avx512bw"))
#else
    if (1)
#endif
    {
        fobj->HAS_AVX512BW = comp_info.AVX512BW;
    }

#if BITS_PER_DIGIT == 64
    fobj->lcg_state = options->rand_seed;
#else
    fobj->lcg_state = (uint32_t)options->rand_seed;
#endif	

    calc_metadata.fobj = fobj;
    calc_metadata.pscreen = options->pscreen;
    calc_metadata.pfile = options->pfile;

	// check/process input arguments
	is_cmdline_run = check_expression(options);
    if (is_cmdline_run == 3)
    {
        // a default function applied to text that has no other function.
        int len = (int)strlen(options->inputExpr) + 9;
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
	print_splash(fobj, &comp_info, is_cmdline_run, logfile, yafu_obj.VFLAG, 
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

        // don't need to re-annouce anything every time the user
        // hits the return key
        options_to_factobj(fobj, options);
        int verbose_level = fobj->VFLAG;
        fobj->VFLAG = yafu_obj.VFLAG = -1;
        for (i = 0; i < options->num_tune_info; i++)
        {
            apply_tuneinfo(&yafu_obj, fobj, options->tune_info[i]);
        }
        fobj->VFLAG = yafu_obj.VFLAG = verbose_level;

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
            result = process_expression(input_str.s, &calc_metadata, fobj->VFLAG < 0, 0);
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
            insize = GSTR_MAXSIZE;
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
    //soe_finalize(sdata);
    //free(sdata);

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

char * get_input(char *input_exp, uint32_t*insize)
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

        if (n >= *insize)
        {
            *insize += GSTR_MAXSIZE;
            input_exp = (char *)xrealloc(input_exp, *insize * sizeof(char));
            if (input_exp == NULL)
            {
                printf("couldn't reallocate string when parsing\n");
                exit(-1);
            }
        }

        input_exp[n++] = (char)c;
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

#ifdef NO_PIPE

#else
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

//#if defined(__MINGW32__)
//    }
//#endif
#endif

	return is_cmdline_run;

}

void print_splash(fact_obj_t *fobj, info_t *comp_info, int is_cmdline_run, 
    FILE* logfile, int VFLAG, double freq, int numwit)
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
        printf("YAFU Version %s\n", YAFU_VERSION_STRING);
        logprint(logfile, "YAFU Version %s\n", YAFU_VERSION_STRING);
#ifdef _MSC_VER
        printf("Built with Microsoft Visual Studio %d\n", _MSC_VER);
        logprint(logfile, "Built with Microsoft Visual Studio % d\n", _MSC_VER);
#elif defined (__INTEL_COMPILER)
        printf("Built with Intel Compiler %d\n", __INTEL_COMPILER);
        logprint(logfile, "Built with Intel Compiler %d\n", __INTEL_COMPILER);
#elif defined (__GNUC__)
        printf("Built with GCC %d\n", __GNUC__);
        logprint(logfile, "Built with GCC %d\n", __GNUC__);
#else 
        printf("Built with undefined compiler\n");
        logprint(logfile, "Built with undefined compiler\n");
#endif

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

    logprint(logfile,"detected %s\ndetected L1 = %d bytes, L2 = %d bytes, CL = %d bytes\n",
		comp_info->idstr, comp_info->L1cache, comp_info->L2cache, comp_info->cachelinesize);
    if (freq > 100.0)
        logprint(logfile,"measured cpu frequency ~= %f\n", freq);
    if (numwit == 1)
        logprint(logfile,"using %u random witness for Rabin-Miller PRP checks\n", numwit);
    else
        logprint(logfile, "using %u random witnesses for Rabin-Miller PRP checks\n", numwit);
    logprint(logfile, "Cached %lu primes: max prime is %lu\n\n", fobj->num_p, fobj->max_p);

	if (VFLAG > 0 || !is_cmdline_run)
	{		
		printf("Detected %s\nDetected L1 = %d bytes, L2 = %d bytes, CL = %d bytes\n",
            comp_info->idstr, comp_info->L1cache, comp_info->L2cache, comp_info->cachelinesize);
        if (freq > 100.0)
		    printf("Measured cpu frequency ~= %f\n", freq);
        if (numwit == 1)
		    printf("Using %u random witness for Rabin-Miller PRP checks\n", numwit);
        else
            printf("Using %u random witnesses for Rabin-Miller PRP checks\n", numwit);

        printf("Cached %lu primes; max prime is %lu\n\n", fobj->num_p, fobj->max_p);

		printf("===============================================================\n");
		printf("======= Welcome to YAFU (Yet Another Factoring Utility) =======\n");
		printf("=======             bbuhrow@gmail.com                   =======\n");
		printf("=======     Type help at any time, or quit to quit      =======\n");
		printf("===============================================================\n");
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
    double xover;

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
    //    osstr, cpustr, yobj->CPU_ID_STR);

    // "xover" trumps tune info.  I.e., if a specific crossover has been
    // specified, prefer this to whatever may be in the tune_info string.
    xover = fobj->autofact_obj.qs_gnfs_xover;

#if defined(_WIN64)
	if ((strcmp(cpustr, yobj->CPU_ID_STR) == 0) && (strcmp(osstr, "WIN64") == 0))
	{
        if (yobj->VFLAG > 0)
		    printf("Applying tune_info entry for %s - %s\n",osstr,cpustr);
		
		sscanf(arg + i + 1, "%lg, %lg, %lg, %lg, %lg, %lg",
			&fobj->qs_obj.qs_multiplier, &fobj->qs_obj.qs_exponent,
			&fobj->nfs_obj.gnfs_multiplier, &fobj->nfs_obj.gnfs_exponent, 
			&fobj->autofact_obj.qs_gnfs_xover, &fobj->nfs_obj.gnfs_tune_freq);
		fobj->qs_obj.qs_tune_freq = fobj->nfs_obj.gnfs_tune_freq;

	}
#elif defined(WIN32)
	if ((strcmp(cpustr, yobj->CPU_ID_STR) == 0) && (strcmp(osstr, "WIN32") == 0))
	{
        if (yobj->VFLAG > 0)
		    printf("Applying tune_info entry for %s - %s\n",osstr,cpustr);
		sscanf(arg + i + 1, "%lg, %lg, %lg, %lg, %lg, %lg",
			&fobj->qs_obj.qs_multiplier, &fobj->qs_obj.qs_exponent,
			&fobj->nfs_obj.gnfs_multiplier, &fobj->nfs_obj.gnfs_exponent, 
			&fobj->autofact_obj.qs_gnfs_xover, &fobj->nfs_obj.gnfs_tune_freq);
		fobj->qs_obj.qs_tune_freq = fobj->nfs_obj.gnfs_tune_freq;
	}
#else 

    //printf("cpustr: %s\nyobj->CPU_ID_STR: %s\nosstr: %s\n", cpustr, yobj->CPU_ID_STR, osstr);
	if ((strcmp(cpustr, yobj->CPU_ID_STR) == 0) && (strcmp(osstr, "LINUX64") == 0))
	{
        if (yobj->VFLAG > 0)
            printf("Applying tune_info entry for %s - %s\n",osstr,cpustr);
		
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

void options_to_factobj(fact_obj_t* fobj, options_t* options)
{
    // set any parameters of fobj changed by user options
    uint32_t seed1, seed2;
    int i;

    fobj->num_threads = options->threads;
    strcpy(fobj->flogname, options->factorlog);

    // initialize global stuff in fobj
    if (options->rand_seed != 0)
    {
        fobj->seed1 = (uint32_t)options->rand_seed & 0xffffffff;
        fobj->seed2 = (uint32_t)(options->rand_seed >> 32);
        fobj->lcg_state = options->rand_seed;
    }

    fobj->NUM_WITNESSES = options->num_prp_witnesses;

    // initialize stuff for rho	
    fobj->rho_obj.iterations = options->rhomax;
    fobj->rho_obj.curr_poly = 0;

    // initialize stuff for pm1	
    fobj->pm1_obj.B1 = options->B1pm1;
    fobj->pm1_obj.B2 = options->B2pm1;
    fobj->pm1_obj.stg2_is_default = (options->B2pm1 == 0);
    fobj->pm1_obj.pm1_exponent = 0;
    fobj->pm1_obj.pm1_multiplier = 0;
    fobj->pm1_obj.pm1_tune_freq = 0;
    strncpy(fobj->pm1_obj.vpm1_work_file, options->vpm1_work_file, 256);
    strncpy(fobj->pm1_obj.resume_file, options->resume_file, 256);

    // initialize stuff for pp1	
    fobj->pp1_obj.B1 = options->B1pp1;
    fobj->pp1_obj.B2 = options->B2pp1;
    fobj->pp1_obj.stg2_is_default = (options->B2pp1 == 0);
    fobj->pp1_obj.pp1_exponent = 0;
    fobj->pp1_obj.pp1_multiplier = 0;
    fobj->pp1_obj.pp1_tune_freq = 0;
    fobj->pp1_obj.lcg_state = hash64(lcg_rand_64(&fobj->lcg_state));
    strncpy(fobj->pp1_obj.vpp1_work_file, options->vpp1_work_file, 256);
    strncpy(fobj->pp1_obj.resume_file, options->resume_file, 256);

    // initialize stuff for ecm	
    fobj->ecm_obj.B1 = options->B1ecm;
    fobj->ecm_obj.B2 = options->B2ecm;
    fobj->ecm_obj.stg2_is_default = (options->B2ecm == 0);
    fobj->ecm_obj.sigma = options->sigma;
    fobj->ecm_obj.num_curves = 90;
    fobj->ecm_obj.curves_run = 0;
    fobj->ecm_obj.ecm_exponent = 0;
    fobj->ecm_obj.ecm_multiplier = 0;
    fobj->ecm_obj.ecm_tune_freq = 0;
    fobj->ecm_obj.bail_on_factor = 1;
    fobj->ecm_obj.save_b1 = options->saveB1;
    fobj->ecm_obj.rand_seed1 = fobj->seed1;
    fobj->ecm_obj.rand_seed2 = fobj->seed2;
    fobj->ecm_obj.lcg_state = (uint64_t*)xrealloc(fobj->ecm_obj.lcg_state,
        options->threads * sizeof(uint64_t));
    for (i = 0; i < (int)options->threads; i++)
    {
        fobj->ecm_obj.lcg_state[i] =
            hash64(lcg_rand_64(&fobj->lcg_state));
    }

    // unlike ggnfs, ecm does not *require* external binaries.  
    // an empty string indicates the use of the built-in GMP-ECM hooks, while
    // a non-empty string (filled in by the user) will indicate the use of
    // an external binary
    strcpy(fobj->ecm_obj.ecm_path, options->ecm_path);
    fobj->ecm_obj.use_external = 0;
    fobj->ecm_obj.prefer_gmpecm = options->prefer_gmpecm;
    fobj->ecm_obj.prefer_gmpecm_stg2 = options->prefer_gmpecm_stg2;
    fobj->ecm_obj.ecm_ext_xover = options->ext_ecm_xover;

    // initialize stuff for squfof
    fobj->squfof_obj.num_factors = 0;

    // initialize stuff for qs	
    fobj->qs_obj.gbl_override_B_flag = (options->siqsB > 0);
    fobj->qs_obj.gbl_override_B = options->siqsB;
    fobj->qs_obj.gbl_override_small_cutoff_flag = (options->siqsTFSm > 0);
    fobj->qs_obj.gbl_override_small_cutoff = options->siqsTFSm;
    fobj->qs_obj.gbl_override_blocks_flag = (options->siqsNB > 0);
    fobj->qs_obj.gbl_override_blocks = options->siqsNB;
    fobj->qs_obj.gbl_override_lpmult_flag = (options->siqsM > 0);
    fobj->qs_obj.gbl_override_lpmult = options->siqsM;
    fobj->qs_obj.gbl_override_rel_flag = (options->siqsR > 0);
    fobj->qs_obj.gbl_override_rel = options->siqsR;
    fobj->qs_obj.gbl_override_tf_flag = (options->siqsTF > 0);
    fobj->qs_obj.gbl_override_tf = options->siqsTF;
    fobj->qs_obj.gbl_override_time_flag = (options->siqsT > 0);
    fobj->qs_obj.gbl_override_time = options->siqsT;
    fobj->qs_obj.inmem_cutoff = options->inmem_cutoff;
    if ((options->siqsMFBD < 1.0) || (options->siqsMFBD > 2.0))
    {
        if (options->siqsMFBD < 1.0)
            options->siqsMFBD = 1.0;
        if (options->siqsMFBD > 2.0)
            options->siqsMFBD = 2.0;
        printf("DLP Tdiv exponent should be between 1.0 and 2.0; setting MFBD = %1.2f\n",
            options->siqsMFBD);
    }
    if ((options->siqsMFBT < 2.0) || (options->siqsMFBT > 3.0))
    {
        if (options->siqsMFBT < 2.0)
            options->siqsMFBT = 2.0;
        if (options->siqsMFBT > 3.0)
            options->siqsMFBT = 3.0;
        printf("TLP Tdiv exponent should be between 2.0 and 3.0; setting MFBT = %1.2f\n",
            options->siqsMFBT);
    }
    fobj->qs_obj.gbl_override_mfbd = options->siqsMFBD;
    fobj->qs_obj.gbl_override_mfbt = options->siqsMFBT;
    fobj->qs_obj.gbl_override_lpb = options->siqsLPB;
    if (options->siqsBDiv < 1.0)
    {
        printf("Batch divisor cannot be less than 1.0; setting BDiv = 1.0\n");
        options->siqsBDiv = 1.0;
    }
    fobj->qs_obj.gbl_override_bdiv = (float)options->siqsBDiv;
    fobj->qs_obj.gbl_override_3lp_bat = options->siqsNobat;
    fobj->qs_obj.gbl_btarget = options->siqsBT;
    fobj->qs_obj.flags = 0;
    fobj->qs_obj.gbl_force_DLP = options->siqsForceDLP;
    fobj->qs_obj.gbl_force_TLP = options->siqsForceTLP;
    fobj->qs_obj.qs_exponent = 0;
    fobj->qs_obj.qs_multiplier = 0;
    fobj->qs_obj.qs_tune_freq = 0;
    fobj->qs_obj.no_small_cutoff_opt = options->no_opt;
    strcpy(fobj->qs_obj.siqs_savefile, options->qssave);

    init_lehman();

    // initialize stuff for trial division	
    fobj->div_obj.print = 0;
    fobj->div_obj.limit = 10000;
    fobj->div_obj.fmtlimit = options->fermat_max;

    //initialize stuff for nfs
    fobj->nfs_obj.snfs = 0;
    fobj->nfs_obj.gnfs = options->force_gnfs;
    fobj->nfs_obj.gnfs_exponent = 0;
    fobj->nfs_obj.gnfs_multiplier = 0;
    fobj->nfs_obj.gnfs_tune_freq = 0;
    fobj->nfs_obj.min_digits = 85;
    // raise min_rels bounds by a percentage
    // on unsuccessful filtering
    fobj->nfs_obj.filter_min_rels_nudge = 1.0 + options->filt_bump / 100.0;

    fobj->nfs_obj.siever = options->ggnfs_siever;
    fobj->nfs_obj.startq = options->sieveQstart;
    fobj->nfs_obj.rangeq = options->sieveQstop;
    if (fobj->nfs_obj.rangeq < fobj->nfs_obj.startq)
    {
        printf("ignoring sieve stop < sieve start\n");
    }
    else
    {
        fobj->nfs_obj.rangeq = fobj->nfs_obj.rangeq - fobj->nfs_obj.startq;
    }
    if (options->sieveQstart > 0)
    {
        if ((options->sieveQstart == 1) && (options->sieveQstop == 1))
        {
            fobj->nfs_obj.startq = fobj->nfs_obj.rangeq = 0;
        }
        fobj->nfs_obj.nfs_phases |= NFS_PHASE_SIEVE;
    }

    fobj->nfs_obj.polystart = options->polystart;
    fobj->nfs_obj.polyrange = options->polystop;
    if (fobj->nfs_obj.rangeq < fobj->nfs_obj.startq)
    {
        printf("ignoring poly stop < poly start\n");
    }
    else
    {
        fobj->nfs_obj.polyrange = fobj->nfs_obj.polyrange - fobj->nfs_obj.polystart;
    }
    if (options->polystart > 0)
    {
        if ((options->polystart == 1) && (options->polystop == 1))
        {
            fobj->nfs_obj.polystart = fobj->nfs_obj.polyrange = 0;
        }
        fobj->nfs_obj.nfs_phases |= NFS_PHASE_POLY;
    }

    if (options->nc)
    {
        fobj->nfs_obj.nfs_phases |= NFS_PHASE_FILTER;
        fobj->nfs_obj.nfs_phases |= NFS_PHASE_LA;
        fobj->nfs_obj.nfs_phases |= NFS_PHASE_SQRT;
    }
    if (options->ncr)
    {
        fobj->nfs_obj.nfs_phases |= NFS_PHASE_LA_RESUME;
    }
    if (options->nc1)
    {
        fobj->nfs_obj.nfs_phases |= NFS_PHASE_FILTER;
    }
    if (options->nc2)
    {
        fobj->nfs_obj.nfs_phases |= NFS_PHASE_LA;
    }
    if (options->nc3)
    {
        fobj->nfs_obj.nfs_phases |= NFS_PHASE_SQRT;
    }

    char tmp[MAXARGLEN];
    char* cptr;
    strcpy(fobj->nfs_obj.outputfile, options->nfs_outfile);
    strcpy(tmp, fobj->nfs_obj.outputfile);
    cptr = strchr(tmp, 46);
    if (cptr == NULL)
    {
        //no . in provided filename
        sprintf(fobj->nfs_obj.outputfile, "%s.dat", fobj->nfs_obj.outputfile);
        sprintf(fobj->nfs_obj.logfile, "%s.log", fobj->nfs_obj.outputfile);
        sprintf(fobj->nfs_obj.fbfile, "%s.fb", fobj->nfs_obj.outputfile);
    }
    else
    {
        cptr[0] = '\0';
        sprintf(fobj->nfs_obj.logfile, "%s.log", tmp);
        sprintf(fobj->nfs_obj.fbfile, "%s.fb", tmp);
    }

    if (options->rat_side)
        fobj->nfs_obj.sq_side = -1;
    else if (options->alg_side)
        fobj->nfs_obj.sq_side = 1;
    else
        fobj->nfs_obj.sq_side = 0;					// default: choose what makes sense

    fobj->nfs_obj.timeout = options->nfs_timeout;
    strcpy(fobj->nfs_obj.job_infile, options->nfs_jobfile);

    // default = fast search
    fobj->nfs_obj.poly_option = 4;
    if (strcmp(options->poly_method, "wide") == 0)
        fobj->nfs_obj.poly_option = 1;
    else if (strcmp(options->poly_method, "deep") == 0)
        fobj->nfs_obj.poly_option = 2;
    else if (strcmp(options->poly_method, "fast") == 0)
        fobj->nfs_obj.poly_option = 0;
    else if (strcmp(options->poly_method, "min") == 0)
        fobj->nfs_obj.poly_option = 3;
    else if (strcmp(options->poly_method, "avg") == 0)
        fobj->nfs_obj.poly_option = 4;
    else if (strcmp(options->poly_method, "good") == 0)
        fobj->nfs_obj.poly_option = 5;

    fobj->nfs_obj.restart_flag = options->nfs_resume;
    fobj->nfs_obj.snfs_testsieve_threshold = options->snfs_testsieve_threshold;
    *fobj->nfs_obj.filearg = '\0';

    fobj->nfs_obj.polybatch = options->poly_batch;
    strcpy(fobj->nfs_obj.ggnfs_dir, options->ggnfs_dir);

    // initialize autofactor object
    // whether we want to output certain info to their own files...
    fobj->autofact_obj.want_output_primes = 0;
    fobj->autofact_obj.want_output_factors = 0;
    fobj->autofact_obj.want_output_unfactored = 0;
    fobj->autofact_obj.want_output_expressions = options->want_output_expr;
    if (strlen(options->opfile) > 0)
    {
        fobj->autofact_obj.want_output_primes = 1;
        strcpy(fobj->autofact_obj.op_str, options->opfile);
    }
    if (strlen(options->offile) > 0)
    {
        fobj->autofact_obj.want_output_factors = 1;
        strcpy(fobj->autofact_obj.of_str, options->offile);
    }
    if (strlen(options->oufile) > 0)
    {
        fobj->autofact_obj.want_output_unfactored = 1;
        strcpy(fobj->autofact_obj.ou_str, options->oufile);
    }

    fobj->autofact_obj.qs_gnfs_xover = options->xover;
    fobj->autofact_obj.qs_snfs_xover = options->qs_snfs_xover;
    // use xover even when timing info is available
    fobj->autofact_obj.prefer_xover = (options->xover != 95);
    fobj->autofact_obj.want_only_1_factor = options->one_factor;
    fobj->autofact_obj.no_ecm = options->no_ecm;
    fobj->autofact_obj.target_pretest_ratio = options->pretest_ratio;
    fobj->autofact_obj.initial_work = options->work;
    fobj->autofact_obj.has_snfs_form = -1;		// not checked yet

    //pretesting plan used by factor()
    fobj->autofact_obj.yafu_pretest_plan = PRETEST_NORMAL;
    strcpy(fobj->autofact_obj.plan_str, options->fact_plan);

    if (strcmp(options->fact_plan, "none") == 0)
        fobj->autofact_obj.yafu_pretest_plan = PRETEST_NONE;
    else if (strcmp(options->fact_plan, "noecm") == 0)
        fobj->autofact_obj.yafu_pretest_plan = PRETEST_NOECM;
    else if (strcmp(options->fact_plan, "light") == 0)
        fobj->autofact_obj.yafu_pretest_plan = PRETEST_LIGHT;
    else if (strcmp(options->fact_plan, "normal") == 0)
        fobj->autofact_obj.yafu_pretest_plan = PRETEST_NORMAL;
    else if (strcmp(options->fact_plan, "deep") == 0)
        fobj->autofact_obj.yafu_pretest_plan = PRETEST_DEEP;
    else if (strcmp(options->fact_plan, "custom") == 0)
        fobj->autofact_obj.yafu_pretest_plan = PRETEST_CUSTOM;
    fobj->autofact_obj.only_pretest = options->pretest;
    fobj->autofact_obj.autofact_active = 0;

    // if a number is <= aprcl_prove_cutoff, we will prove it prime or composite
    fobj->factors->aprcl_prove_cutoff = options->aprcl_p;
    // if a number is >= aprcl_display_cutoff, we will show the APRCL progress
    fobj->factors->aprcl_display_cutoff = options->aprcl_d;

    fobj->VFLAG = options->verbosity;
    if (strlen(options->factorlog) == 0)
    {
        fobj->LOGFLAG = 0;
    }
    else
    {
        fobj->LOGFLAG = 1;
    }

    return;
}


