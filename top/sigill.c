#if !defined(_WIN32)
#define _POSIX_C_SOURCE 200809L
#endif

#include "sigill.h"

#include <signal.h>

#if defined(_WIN32)
#include <io.h>
#include <process.h>
#else
#include <unistd.h>
#endif

static int sigill_log_fd = -1;

#if defined(_WIN32)
typedef int sigill_write_result_t;
#else
typedef ssize_t sigill_write_result_t;
#endif

static void write_sigill_message(int fd)
{
    static const char message[] =
        "Unsupported SIMD instruction set detected. This CPU does not support "
        "the instructions required by this build.\n";
    const char *cursor = message;
    size_t remaining = sizeof(message) - 1;

    while (remaining > 0)
    {
#if defined(_WIN32)
        sigill_write_result_t written =
            _write(fd, cursor, (unsigned int)remaining);
#else
        sigill_write_result_t written = write(fd, cursor, remaining);
#endif
        size_t written_size;

        if (written <= 0)
            break;
        written_size = (size_t)written;
        cursor += written_size;
        remaining -= written_size;
    }
}

static void handle_sigill(int signal_number)
{
    (void)signal_number;

    if (sigill_log_fd >= 0)
        write_sigill_message(sigill_log_fd);

#if defined(_WIN32)
    write_sigill_message(2);
    _exit(1);
#else
    write_sigill_message(STDERR_FILENO);
    _exit(1);
#endif
}

int install_sigill_handler(FILE *logfile)
{
    if (logfile != NULL)
    {
        if (fflush(logfile) != 0)
            return 0;
#if defined(_WIN32)
        sigill_log_fd = _fileno(logfile);
#else
        sigill_log_fd = fileno(logfile);
#endif
        if (sigill_log_fd < 0)
            return 0;
    }

    return signal(SIGILL, handle_sigill) != SIG_ERR;
}
