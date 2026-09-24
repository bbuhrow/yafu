/*
 * arith/mod64.c
 *
 * Portable C replacement for mod64.asm.
 * Any optimising compiler (MSVC, ClangCL, GCC) will emit a single DIV
 * instruction for this with /O2 or -O2, identical to the hand-written asm.
 *
 * To migrate:
 *   - Replace arith/mod64.asm with this file
 *   - In the VS project: remove the YASM item for mod64.asm, add this as ClCompile
 *   - Remove the vsyasm.props and vsyasm.targets imports from yafu.vcxproj
 */

#include <stdint.h>

uint64_t mod_64(uint64_t a, uint64_t b)
{
    return a % b;
}
