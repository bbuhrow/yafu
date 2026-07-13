#include <check.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/*
 * The vulnerable code in zlib/minigzip.c uses sprintf into a fixed buffer
 * (typically 1024 bytes) without bounds checking. We test that calling
 * the minigzip binary with oversized filenames does not cause a crash
 * (buffer overflow). We invoke it via subprocess to exercise the real code.
 */

START_TEST(test_minigzip_buffer_overflow)
{
    /* Invariant: Buffer reads/writes never exceed declared length;
       oversized inputs must be rejected or truncated, not overflow. */
    const char *payloads[] = {
        /* Exact exploit: path that exceeds typical 1024-byte buffer */
        NULL, /* will be generated dynamically */
        NULL, /* 2x buffer: 2048 chars */
        "valid_short_name.txt" /* valid input */
    };

    /* Generate oversized payloads */
    char payload_10x[10240 + 4];
    memset(payload_10x, 'A', 10240);
    strcpy(payload_10x + 10240, ".gz");

    char payload_2x[2048 + 4];
    memset(payload_2x, 'B', 2048);
    strcpy(payload_2x + 2048, ".gz");

    payloads[0] = payload_10x;
    payloads[1] = payload_2x;

    int num_payloads = 3;

    for (int i = 0; i < num_payloads; i++) {
        char cmd[16384];
        snprintf(cmd, sizeof(cmd), "./minigzip -d \"%s\" 2>/dev/null", payloads[i]);
        int ret = system(cmd);
        /* The program should exit gracefully (not crash with signal).
           On UNIX, if killed by signal, WIFSIGNALED is true. */
        if (ret != -1) {
            ck_assert_msg(!WIFSIGNALED(ret) ||
                          (WTERMSIG(ret) != 11 && WTERMSIG(ret) != 6),
                          "minigzip crashed with signal %d on payload %d",
                          WTERMSIG(ret), i);
        }
    }
}
END_TEST

Suite *security_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("Security");
    tc_core = tcase_create("Core");

    tcase_add_test(tc_core, test_minigzip_buffer_overflow);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = security_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}