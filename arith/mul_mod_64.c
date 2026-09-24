/*
 * arith/mul_mod_64.c
 *
 * Portable C replacement for x64_support_masm.asm.
 *
 *   uint64_t mul_mod_64(uint64_t x, uint64_t y, uint64_t m)
 *   Returns (x * y) % m.
 *
 * The multiply x*y produces a 128-bit intermediate result which is then
 * reduced modulo m using a single 128÷64 DIV instruction.  Any optimising
 * compiler (MSVC, ClangCL, GCC) will emit exactly:
 *   MUL  rdx        ; rdx:rax = x * y  (128-bit product)
 *   DIV  r8         ; rax = quotient, rdx = remainder
 * with /O2 or -O2, identical to the hand-written assembly.
 *
 * Precondition (same as the original asm): m > 0 and the quotient of
 * (x*y) / m must fit in 64 bits, i.e. x < m or y < m (or both).
 * Violating this causes a #DE fault on x86-64, same as the original.
 *
 * To migrate:
 *   - Replace arith/x64_support_masm.asm with this file
 *   - In the VS project: remove the YASM/MASM item, add this as ClCompile
 *   - Remove vsyasm.props / vsyasm.targets imports if mod64.asm was the
 *     last remaining asm file (this file completes that migration)
 */

#include <stdint.h>

#if defined(__GNUC__) || defined(__clang__)

uint64_t mul_mod_64(uint64_t x, uint64_t y, uint64_t m)
{
    return (uint64_t)((__uint128_t)x * y % m);
}

#else  /* MSVC */

#include <intrin.h>  /* _umul128, _udiv128 — VS 2019 16.8+ */

uint64_t mul_mod_64(uint64_t x, uint64_t y, uint64_t m)
{
    uint64_t hi, lo, rem;
    lo = _umul128(x, y, &hi);   /* rdx:rax = x * y */
    _udiv128(hi, lo, m, &rem);  /* rem    = (hi:lo) % m */
    return rem;
}

#endif
