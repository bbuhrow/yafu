#include <check.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/wait.h>

#define GZ_SUFFIX ".gz"
#define MAX_NAME_LEN 1024

START_TEST(test_filename_buffer_overflow_protection)
{
    // Invariant: Buffer overflow via strcpy/strcat on filename must not corrupt memory
    // Test by spawning minigzip with adversarial filenames and checking for crash/segfault
    
    const char *payloads[] = {
        "a",                                    // Valid: short filename
        "normalfile.txt",                       // Valid: typical filename
        "x" "x" "x" "x" "x" "x" "x" "x"        // Boundary: 8 chars
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x" "x"
        "x" "x" "x" "x" "x" "x" "x"