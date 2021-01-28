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

#ifndef CMD_PARSE_H
#define CMD_PARSE_H

#ifdef __cplusplus
extern "C" {
#endif

    // the number of recognized command line options
#define NUMOPTIONS 86
// maximum length of command line option strings
#define MAXOPTIONLEN 20
// maximum length of help string for each option
#define MAXHELPLEN 128
// maximum length of an argument to an option
#define MAXARGLEN 256

typedef struct
{
    char usageHelp[MAXHELPLEN];

    // command line options (names to be preceeded by a -)
    char OptionArray[NUMOPTIONS][MAXOPTIONLEN];
    char OptionHelp[NUMOPTIONS][MAXHELPLEN];
    char LongOptionAliases[NUMOPTIONS][MAXOPTIONLEN];

    // indication of whether or not an option needs a corresponding argument
    // 0 = no argument
    // 1 = argument required
    // 2 = argument optional
    int needsArg[NUMOPTIONS];

    // ========================================================================
    // This defines how many required and optional arguments
    // there are.  These are things sent into the program that define its
    // behavior and don't have "-" or "-- in front of them.
    int numArguments;
    int numRequired;
    // ========================================================================


    // ========================================================================
    // These variables define the required and optional arguments
    char inputArg[MAXARGLEN];
    // ========================================================================


    // ========================================================================
    // These variables define things that change program behavior and can be
    // set via command line switches.  I.e., these need a "-" or "--" in
    // front of them on the command line... add items here as necessary
    int optionInt;
    char optionStr[MAXARGLEN];
    double optionDbl;
    int optionNoarg;
    // ========================================================================

} options_t;

extern options_t* initOpt(void);
extern void applyOpt(char* opt, char* arg, options_t* options);
extern int processOpts(int argc, char** argv, options_t* options);



#ifdef __cplusplus
}
#endif

#endif /* #ifndef CMD_PARSE_H */

