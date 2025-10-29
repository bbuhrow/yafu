/*
MIT License

Copyright (c) 2021 Ben Buhrow

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

// ============================================================================
// a simple command line parser.
// supports defining required and options arguments to a program.
// supports defining command line switches with optional arguments.
// ============================================================================


#include "cmdOptions.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "gmp.h"

void printUsage(options_t* options);

// ========================================================================
// In this section, define the options, any long form aliases for them,
// whether they have an argument or not, and help strings to display
// if the user tries -h

// define help text for usage of the program.  
char usageHelp[MAXHELPLEN] = "num-to-factor B1 [B2] [options]";

// command line options, specified by '-'
char OptionArray[NUMOPTIONS][MAXOPTIONLEN] = { 
    "t", "s", "c", "v", "f", 
    "p", "q", "g", "d", "e",
    "x"};

// command line option aliases, specified by '--'
// need the same number of strings here, even if
// some of them are blank (i.e., have no long form alias).
char LongOptionAliases[NUMOPTIONS][MAXOPTIONLEN] = {
    "threads", "sigma", "count", "verbose", "save",
    "prefer_gmpecm", "prefer_gmpecm_stg2", "gpu", "gpudev", "external",
    "xover"};

// indication of whether or not an option needs a corresponding argument.
// needs to be the same length as the above two arrays.
// 0 = no argument
// 1 = argument required
// 2 = argument optional
int needsArg[NUMOPTIONS] = {
    1,1,1,0,0,
    0,0,0,1,1,
    1};

// help strings displayed with -h
// needs to be the same length as the above arrays, even if 
// some are left blank
char OptionHelp[NUMOPTIONS][MAXHELPLEN] = {
    "(int)  Number of threads",
    "(int)  Input curve sigma (suyama parameterization)",
    "(int)  Curve count (will get rounded up to multiple of threads * 8 when using avx-ecm",
    "       Increase verbosity - can specify multiple times",
    "       Save B1 results to file save_b1.txt",
    "(int)  Force gmp-ecm for stage 1 and stage 2",
    "(int)  Force gmp-ecm for stage 2",
    "(int)  Use gpu-ecm with cgbn",
    "(int)  Specify gpu device number",
    "(char) Path to external gmp-ecm executable file",
    "(int)  Crossover to switch to external gmp-ecm"};
// ========================================================================

// ========================================================================
// modify this function to handle input arguments
// ========================================================================
void applyArg(char* arg, int argNum, options_t* options)
{
    if (argNum == 0)
    {
        strncpy(options->inputN, arg, MAXARGLEN);
    }
    else if (argNum == 1)
    {
        strncpy(options->B1, arg, MAXARGLEN);
    }
    else if (argNum == 2)
    {
        strncpy(options->B2, arg, MAXARGLEN);
    }

    return;
}

// ========================================================================
// modify this function to handle each option case, assigning values
// or arguments to the members of the options_t structure.
// ========================================================================
void applyOpt(char* opt, char* arg, options_t* options)
{
    if (strcmp(opt, options->OptionArray[0]) == 0)
    {
        options->threads = atoi(arg);
    }
    else if (strcmp(opt, options->OptionArray[1]) == 0)
    {
        options->sigma = strtoull(arg, NULL, 10);
    }
    else if (strcmp(opt, options->OptionArray[2]) == 0)
    {
        options->numcurves = strtoul(arg, NULL, 10);
    }
    else if (strcmp(opt, options->OptionArray[3]) == 0)
    {
        options->verbosity++;
    }
    else if (strcmp(opt, options->OptionArray[4]) == 0)
    {
        options->save_b1 = 1;
    }
    else if (strcmp(opt, options->OptionArray[5]) == 0)
    {
        options->prefer_gmpecm = 1;
    }
    else if (strcmp(opt, options->OptionArray[6]) == 0)
    {
        options->prefer_gmpecm_stg2 = 1;
    }
    else if (strcmp(opt, options->OptionArray[7]) == 0)
    {
        options->gpu = 1;
    }
    else if (strcmp(opt, options->OptionArray[8]) == 0)
    {
        options->gpudev = atoi(arg);
    }
    else if (strcmp(opt, options->OptionArray[9]) == 0)
    {
        strncpy(options->extecm_path, arg, MAXARGLEN);
    }
    else if (strcmp(opt, options->OptionArray[10]) == 0)
    {
        options->extecm_xover = strtoull(arg, NULL, 10);
    }
    else
    {
        int i;

        printf("invalid option %s\n", opt);
        printf("supported options are: ");
        for (i = 0; i < NUMOPTIONS; i++)
        {
            if (i % 5 == 0) printf("\n");
            printf("%s ", options->OptionArray[i]);
        }
        exit(0);
    }

    return;
}

// ========================================================================
// modify this function to handle assign a default value to each 
// member of the options_t structure, if desired.
// ========================================================================
options_t* initOpt(void)
{
    options_t* options = (options_t*)malloc(sizeof(options_t));
    int i;

    if (options == NULL)
    {
        printf("could not allocate options structure\n");
        exit(1);
    }

    for (i = 0; i < NUMOPTIONS; i++)
    {
        strcpy(options->OptionArray[i], OptionArray[i]);
        strcpy(options->OptionHelp[i], OptionHelp[i]);
        strcpy(options->LongOptionAliases[i], LongOptionAliases[i]);
        options->needsArg[i] = needsArg[i];
    }

    // ========================================================================
    // define how many required and optional arguments there are
    options->numArguments = 3;
    options->numRequired = 2;
    // ========================================================================


    // ========================================================================
    // default values assigned to optional arguments here
    strncpy(options->inputN, "0", MAXARGLEN);
    strncpy(options->B1, "11000", MAXARGLEN);
    strncpy(options->B2, "0", MAXARGLEN);
    // ========================================================================


    // ========================================================================
    // default values assigned to switches here:
    options->verbosity = 0;
    strcpy(options->extecm_path, "");
    options->save_b1 = 0;
    options->extecm_xover = 50000;
    options->threads = 1;
    options->sigma = 0;
    options->prefer_gmpecm = 0;
    options->prefer_gmpecm_stg2 = 0;
    options->gpu = 0;
    options->gpudev = 0;
    options->numcurves = 1;
    // ========================================================================

    return options;
}




// ========================================================================
// this function should not need to be changed
// ========================================================================
int processOpts(int argc, char** argv, options_t* options)
{
    // this function takes the command line arguments and returns
    // an options_t, which will have any valid options specified
    // on the command line filled in.  Return the number of options
    // processed.
    int i, j, k = 0, valid, longswitch = 0;
    char optbuf[MAXOPTIONLEN+1];
    char argbuf[MAXARGLEN+1];
    int numOpt = 0;

    //argument loop
    i = 1;
    while (i < argc)
    {
        longswitch = 0;

        // read in the option
        if ((strlen(argv[i]) > 1) && (argv[i][0] == '-') && (argv[i][1] == '-'))
        {
            longswitch = 1;
        }
        else if (argv[i][0] == '-')
        {
            longswitch = 0;
        }
        else
        {
            if (numOpt >= (options->numArguments))
            {
                printf("Too many inputs\n\n");
                printUsage(options);
                exit(0);
            }
            else
            {
                // process argument
                applyArg(argv[i], numOpt, options);
                numOpt++;
                i++;
                continue;
            }
        }

        if (numOpt < options->numRequired)
        {
            if (options->numRequired == 1)
            {
                printf("Missing %d required argument \n\n", options->numRequired);
            }
            else
            {
                printf("Missing %d required arguments \n\n", options->numRequired);
            }
            printUsage(options);
            exit(0);
        }

        // check if the options is valid
        valid = 0;
        if (longswitch)
        {
            for (j = 0; j < NUMOPTIONS; j++)
            {
                if (strncmp(options->LongOptionAliases[j], &argv[i][2], MAXOPTIONLEN) == 0)
                {
                    valid = 1;
                    strncpy(optbuf, options->OptionArray[j], MAXOPTIONLEN);
                    break;
                }
            }
        }
        else
        {
            for (j = 0; j < NUMOPTIONS; j++)
            {
                if (strncmp(options->OptionArray[j], &argv[i][1], MAXOPTIONLEN) == 0)
                {
                    valid = 1;
                    strncpy(optbuf, &argv[i][1], MAXOPTIONLEN);
                    break;
                }
            }
        }

        if (valid == 0)
        {
            printUsage(options);
            exit(0);
        }

        //check to see if this option requires an argument
        if (options->needsArg[j] == 1)
        {
            i++;
            if ((i == argc) || argv[i][0] == '-')
            {
                printf("argument expected for %s%s\n\n", longswitch ? "--" : "-", optbuf);
                printUsage(options);
                exit(0);
            }
            strncpy(argbuf, argv[i], MAXARGLEN);

            //now apply -option argument
            applyOpt(optbuf, argbuf, options);
            k++;
        }
        else if (options->needsArg[j] == 2)
        {
            // check to see if an argument was supplied
            if (((i + 1) == argc) || argv[i + 1][0] == '-')
            {
                // no option supplied.  use default option
                applyOpt(optbuf, NULL, options);
                k++;
            }
            else
            {
                i++;
                // an option was supplied, pass it on
                strncpy(argbuf, argv[i], MAXARGLEN);

                //now apply -option argument
                applyOpt(optbuf, argbuf, options);
                k++;
            }
        }
        else
        {
            //apply -option
            //now apply -option argument
            applyOpt(optbuf, NULL, options);
            k++;
        }
        i++;
    }

    if (numOpt < options->numRequired)
    {
        if (options->numRequired == 1)
        {
            printf("Missing %d required argument \n\n", options->numRequired);
        }
        else
        {
            printf("Missing %d required arguments \n\n", options->numRequired);
        }
        printUsage(options);
        exit(0);
    }

    return k;
}

void printUsage(options_t* options)
{
    int j;
    printf("Usage: %s\n", usageHelp);
    printf("supported options are: \n");
    for (j = 0; j < NUMOPTIONS; j++)
    {
        if (strlen(options->LongOptionAliases[j]) > 0)
        {
            if (options->needsArg[j])
            {
                printf("%s <value>: %s (alias --%s)\n",
                    options->OptionArray[j], options->OptionHelp[j],
                    options->LongOptionAliases[j]);
            }
            else
            {
                printf("%s        : %s (alias --%s)\n",
                    options->OptionArray[j], options->OptionHelp[j],
                    options->LongOptionAliases[j]);
            }
        }
        else
        {
            if (options->needsArg[j])
            {
                printf("%s <value>: %s\n", options->OptionArray[j], options->OptionHelp[j]);
            }
            else
            {
                printf("%s        : %s\n", options->OptionArray[j], options->OptionHelp[j]);
            }
        }
    }
    return;
}
