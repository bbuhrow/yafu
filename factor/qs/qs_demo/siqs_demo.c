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

#include "siqs_demo.h"
#include "soe.h"
#include "calc.h"
#include "ytools.h"
#include "factor.h"
#include <stdlib.h>
#include <stdio.h>
#include "cmdOptions.h"

// function to print the splash screen to file/screen
void print_splash(info_t* comp_info, int is_cmdline_run, FILE *logfile, 
    int VFLAG, int numwit);

// functions to make a batchfile ready to execute, and to process batchfile lines
void prepare_batchfile(char *input_exp);
char * process_batchline(siqs_obj_t* yobj, char *input_exp, char *indup, int *code);
void finalize_batchline(siqs_obj_t* yobj);
int exp_is_open(char *line, int firstline);

// functions to process all incoming arguments
int check_expression(options_t *options);
void options_to_factobj(fact_obj_t* fobj, options_t* options);

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
    siqs_obj_t obj;
    soe_staticdata_t* sdata;
    info_t comp_info;
    int run_type;
    int i;

	//the input expression
	input_exp = (char *)malloc(GSTR_MAXSIZE*sizeof(char));
	indup = (char *)malloc(GSTR_MAXSIZE*sizeof(char));
    input_line = (char *)malloc(GSTR_MAXSIZE*sizeof(char));
    sInit(&input_str);
	strcpy(input_exp,"");
    strcpy(input_line, "");
	
	// set defaults for various things and read the .ini file, if any.
	siqs_init(&obj);
    options = initOpt();
    readINI("siqs.ini", options);

    // then process the command line, overriding any .ini settings.
    processOpts(argc, argv, options);

    // some things go into globals, but this is being phased out
    obj.VFLAG = options->verbosity;
    obj.THREADS = options->threads;
    obj.CMD_LINE_REPEAT = options->repeat;
    obj.VERBOSE_PROC_INFO = options->vproc;

    if (strlen(options->factorlog) == 0)
    {
        obj.LOGFLAG = 0;
    }
    if (strlen(options->batchfile) > 0)
    {
        obj.USEBATCHFILE = 1;
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
        obj.USERSEED = 1;
    }

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
	init_factobj(fobj);

    options_to_factobj(fobj, options);
    strcpy(fobj->CPU_ID_STR, comp_info.idstr);
    fobj->HAS_AVX2 = comp_info.AVX2;
    fobj->HAS_AVX = comp_info.AVX;
    fobj->HAS_SSE41 = comp_info.bSSE41Extensions;
    fobj->NUM_WITNESSES = options->num_prp_witnesses;
    fobj->cache_size1 = fobj->L1CACHE = comp_info.L1cache;
    fobj->cache_size2 = fobj->L2CACHE = comp_info.L2cache;
    fobj->LOGFLAG = obj.LOGFLAG;
    fobj->THREADS = obj.THREADS;

    // put a list of primes in the fobj; many algorithms use it
    sdata = soe_init(0, 1, 32768);
    fobj->primes = soe_wrapper(sdata, 0, 100000000, 0, &fobj->num_p, 0, 0);
    fobj->min_p = 2;
    fobj->max_p = fobj->primes[fobj->num_p - 1];

#if BITS_PER_DIGIT == 64
    fobj->lcg_state = options->rand_seed;
#else
    fobj->lcg_state = (uint32_t)options->rand_seed;
#endif	

    calc_metadata.fobj = fobj;

	// check/process input arguments
	run_type = check_expression(options);
    if (run_type == 3)
    {
        // a default function applied to text that has no other function.
        int len = (int)strlen(options->inputExpr) + 9;
        options->inputExpr = (char*)xrealloc(options->inputExpr, len);
        input_exp = (char*)xrealloc(input_exp, len);
        sprintf(input_exp, "siqs(%s)", options->inputExpr);
        strcpy(options->inputExpr, input_exp);
        strcpy(input_line, options->inputExpr);
        run_type = 1;
    }
    else
    {
        strcpy(input_line, options->inputExpr);
        strcpy(input_exp, options->inputExpr);
    }

	if (run_type == 2)
	{
		// batchfile from stdin
        obj.USEBATCHFILE = 2;
	}

	// get the batchfile ready, if requested
	if (obj.USEBATCHFILE)
	{
		prepare_batchfile(input_exp);		
		
		// batchfile jobs are command line in nature
        run_type = 1;
	}

    if (strlen(obj.scriptname) > 0)
    {
        scriptfile = fopen(obj.scriptname, "r");
        if (scriptfile == NULL)
        {
            printf("could not find %s\n", obj.scriptname);
            exit(1);
        }
        else
        {
            run_type = 1;
        }
    }

    if (obj.USEBATCHFILE || (obj.CMD_LINE_REPEAT > 0))
    {
        strcpy(indup, input_exp);	//remember the input expression
    }

    if (run_type == 0)
    {
        printf("An input expression is required\n");
        exit(0);
    }

	//session log
    if (obj.LOGFLAG)
    {
        logfile = fopen(obj.sessionname, "a");
        if (logfile == NULL)
        {
            printf("fopen error: %s\n", strerror(errno));
            printf("couldn't open %s for appending\n", obj.sessionname);
            slog = 0;
        }
        else
        {
            slog = 1;
        }
    }
    else
    {
        logfile = NULL;
        slog = 0;
    }
		
	// print the splash screen, to the logfile and depending on options, to the screen
	print_splash(&comp_info, is_cmdline_run, logfile, obj.VFLAG, 
        obj.NUM_WITNESSES);
	
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
		if (obj.USEBATCHFILE)
		{
			int code;
            input_line = process_batchline(&obj, input_line, indup, &code);
			if (code == 1)
			{
				finalize_batchline(&obj);
				break;
			}
			else if (code == 2)
				continue;
		}
        else if (strlen(obj.scriptname) > 0)
        {
            if (scriptfile != NULL)
            {
                if (fgets(input_line, GSTR_MAXSIZE, scriptfile) == NULL)
                {
                    //    break;
                }
                
            }
        }
		else
		{
			// input expression already read in.  nothing to do here.

		}
		
		// exit, or execute the current expression...
        if ((strcmp(input_line, "quit") == 0) || (strcmp(input_line, "exit") == 0))
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

		// get the next expression, if running a batchfile, or
		// re-display the command prompt
		if (obj.CMD_LINE_REPEAT == 0)
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
			if (obj.USEBATCHFILE)
			{
				// the line from the batchfile finished.  make the temporary file
				// created in processs_batchline the new batchfile, with the line
				// we just finished removed.
				finalize_batchline(&obj);
			}
            else if (scriptfile != NULL)
            {
                if (feof(scriptfile))
                {
                    if (obj.CMD_LINE_REPEAT > 0)
                    {
                        obj.CMD_LINE_REPEAT--;
                        fclose(scriptfile);
                        scriptfile = fopen(obj.scriptname, "r");
                        if (scriptfile == NULL)
                        {
                            printf("could not find %s\n", obj.scriptname);
                            exit(1);
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }
			else if (obj.CMD_LINE_REPEAT > 0)
			{
                obj.CMD_LINE_REPEAT--;
                strcpy(input_line, indup);
			}
            else
            {
                break;
            }
		}
        else
        {
            break;
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
	free(input_exp);
    free(input_line);
	free(indup);	
	free_factobj(fobj);
	free(fobj);      
    sFree(&input_str);
    free(options->inputExpr);
    soe_finalize(sdata);
    free(sdata);

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
    int VFLAG, int numwit)
{
	if (VFLAG >= 0)
		printf("\n\n");

	if (VFLAG > 0)
	{	
		logprint(NULL,"System/Build Info: \n");
        fflush(stdout);
	}
    logprint(logfile, "=====================================\n");
	logprint(logfile,"System/Build Info: \n");

    if (VFLAG > 0)
    {
#ifdef _MSC_VER
        printf("Built with Microsoft Visual Studio %d\n", _MSC_VER);
#elif defined (__INTEL_COMPILER)
        printf("Built with Intel Compiler %d\n", __INTEL_COMPILER);
#elif defined (__GNUC__)
        printf("Built with GCC %d\n", __GNUC__);
#else 
        printf("Built with undefined compiler\n");
#endif

#ifdef _MSC_MPIR_VERSION
        printf("Powered by MPIR %s\n", _MSC_MPIR_VERSION);
        logprint(logfile,"Powered by MPIR %s\n",_MSC_MPIR_VERSION);
#else
        printf("Powered by GMP\n");
        logprint(logfile,"Powered by GMP\n");
#endif

        fflush(stdout);
    }

    logprint(logfile,"detected %s\ndetected L1 = %d bytes, L2 = %d bytes, CL = %d bytes\n",
		comp_info->idstr, comp_info->L1cache, comp_info->L2cache, comp_info->cachelinesize);
    if (numwit == 1)
        logprint(logfile,"using %u random witness for Rabin-Miller PRP checks\n", numwit);
    else
        logprint(logfile, "using %u random witnesses for Rabin-Miller PRP checks\n", numwit);

	if (VFLAG > 0)
	{		
		printf("Detected %s\nDetected L1 = %d bytes, L2 = %d bytes, CL = %d bytes\n",
            comp_info->idstr, comp_info->L1cache, comp_info->L2cache, comp_info->cachelinesize);
        if (numwit == 1)
		    printf("Using %u random witness for Rabin-Miller PRP checks\n\n", numwit);
        else
            printf("Using %u random witnesses for Rabin-Miller PRP checks\n\n", numwit);
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

void siqs_init(siqs_obj_t* yobj)
{
    yobj->VFLAG = 0;
    yobj->VERBOSE_PROC_INFO = 0;
    yobj->LOGFLAG = 1;
    yobj->NUM_WITNESSES = 1;
    yobj->USEBATCHFILE = 0;
    yobj->USERSEED = 0;
    yobj->THREADS = 1;
    yobj->CMD_LINE_REPEAT = 0;
    yobj->MEAS_CPU_FREQUENCY = 0.0;

	strcpy(yobj->sessionname,"session.log");
    strcpy(yobj->scriptname, "");

	return;
}

void finalize_batchline(siqs_obj_t* yobj)
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

char * process_batchline(siqs_obj_t* yobj, char *input_exp, char *indup, int *code)
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

