dnl Copyright (C) 2002 Jens Franke, T.Kleinjung
dnl This file is part of gnfs4linux, distributed under the terms of the
dnl GNU General Public Licence and WITHOUT ANY WARRANTY.
dnl
dnl You should have received a copy of the GNU General Public License along
dnl with this program; see the file COPYING.  If not, write to the Free
dnl Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.
dnl
dnl argumente: rdi, rsi, rdx, rcx
dnl funktionen:
dnl asm_zero64(a): a=0
dnl asm_copy64(b,a): b=a
dnl asm_sub64_3(c,a,b): c=a-b mod N
dnl asm_half64(a): a/=2 mod N
dnl asm_sub_n64(b,a): b-=a  mod 2^64
dnl asm_diff64(c,a,b): c=|a-b|
dnl asm_add64(b,a): b+=a mod N
dnl asm_add64_ui(b,a): b+=a mod N, a is ulong
dnl asm_mulm64(): prod=f1*f2 mod N
dnl function_head(asm_sqm64)
dnl asm_inv64(res,b)
dnl

#include "underscore.h"

.comm montgomery_modulo_n,8
.comm montgomery_inv_n,8
dnl
define(rop,%rdi)dnl
define(op,%rsi)dnl
define(op1,%rsi)dnl
define(op2,%rdx)dnl
dnl asm_zero64(a): a=0
function_head(asm_zero64)
	xorq %rax,%rax
	movq %rax,(rop)
	ret
dnl asm_copy64(b,a): b=a
function_head(asm_copy64)
	movq (op),%rax
	movq %rax,(rop)
	ret
dnl asm_sub64_3(c,a,b): c=a-b mod N
function_head(asm_sub64_3)
	movq (op1),%rax
dnl smjs	movq montgomery_modulo_n,%r9
	movq montgomery_modulo_n(%rip),%r9

	xorq %r8,%r8
	subq (op2),%rax
	cmovcq (%r9),%r8
	addq %r8,%rax
	movq %rax,(rop)
	ret
dnl asm_half64(a): a/=2 mod N
function_head(asm_half64)
	movq (rop),%rax
dnl smjs	movq montgomery_modulo_n(%rip),%r8
	movq montgomery_modulo_n(%rip),%r8

	movq $0,%r10
	testq $1,%rax
	movq (%r8),%r9
	cmovnzq %r9,%r10
	addq %r10,%rax
	rcrq $1,%rax
	movq %rax,(rop)
	ret
dnl asm_sub_n64(b,a): b-=a  mod 2^64
dnl smjs For windows calling want either all implemented or non for sub_n, so unimplement here to match 128 and 192
dnl function_head(asm_sub_n64)
function_head(xasm_sub_n64)
	movq (rop),%rax
	subq (op),%rax
	movq %rax,(rop)
	ret
dnl asm_add64(b,a): b+=a mod N
function_head(asm_add64)
	movq (rop),%rax
dnl smjs	movq montgomery_modulo_n,%r10
	movq montgomery_modulo_n(%rip),%r10

	xorq %r8,%r8
	movq (%r10),%r9
	subq %r9,%rax
	addq (op),%rax
	cmovncq %r9,%r8
	addq %r8,%rax
	movq %rax,(rop)
	ret
dnl asm_mulm64(): prod=f1*f2 mod N
function_head(asm_mulm64)
	movq (op2),%rax
	mulq (op1)
	movq montgomery_inv_n(%rip),%r8
	movq montgomery_modulo_n(%rip),%r9
	movq %rax,%r10
	movq %rdx,%r11
	mulq %r8
	movq (%r9),%r9
	xorq %r8,%r8   # useless but seems to save 0.1 cycles
	subq %r9,%r11
	mulq %r9
	addq %r10,%rax
	adcq %r11,%rdx
	cmovncq %r9,%rax
	addq %rax,%rdx
	movq %rdx,(rop)
	ret
undefine(`rop')dnl
undefine(`op')dnl
undefine(`op1')dnl
undefine(`op2')dnl

.if 0
dnl static void ecm_duplicate(ulong *x1, ulong *z1)
dnl {
dnl   asm_copy(mm_w,x1); asm_add2(mm_w,z1);
dnl   asm_squmod(mm_u,mm_w); /* u=(x1+z1)^2 */
dnl   asm_sub(mm_w,x1,z1);
dnl   asm_squmod(mm_v,mm_w); /* v=(x1-z1)^2 */
dnl   asm_mulmod(x1,mm_u,mm_v); /* x2=(u*v) */
dnl   asm_sub(mm_w,mm_u,mm_v); /* w=u-v=4*x1*z1 */
dnl   asm_mulmod(mm_u,mm_b,mm_w);
dnl   asm_add2(mm_u,mm_v); /* u=(v+b*w) */
dnl   asm_mulmod(z1,mm_w,mm_u); /* z2=(w*u) */
dnl }

dnl asm_duplicate()
function_head(asm_duplicate)
	pushq %rbx
	pushq %r12
	pushq %r13
	pushq %r14
	pushq %r15
	movq montgomery_modulo_n,%r10
	movq (%r10),%r15   # n

        movq (%rdi),%r10
        movq (%rsi),%r11
	movq %r10,%rax
        xorq %r8,%r8
        subq %r15,%rax
        addq %r11,%rax
        cmovncq %r15,%r8
        addq %r8,%rax
        movq %rax,%r14

        movq %r10,%rax
        xorq %r8,%r8
        subq %r11,%rax
        cmovcq %r15,%r8
        addq %r8,%rax
        movq %rax,%rbx


        movq %r14,%rax
        mulq %r14
        movq montgomery_inv_n(%rip),%r8
        movq %rax,%r10
        movq %rdx,%r11
        mulq %r8
        movq %r15,%r9
        xorq %r8,%r8   # useless but seems to save 0.1 cycles
        subq %r9,%r11
        mulq %r9
        addq %r10,%rax
        adcq %r11,%rdx
        cmovncq %r9,%rax
        addq %rax,%rdx
        movq %rdx,%r14

        movq %rbx,%rax
        mulq %rbx
        movq montgomery_inv_n(%rip),%r8
        movq %rax,%r10
        movq %rdx,%r11
        mulq %r8
        movq %r15,%r9
        xorq %r8,%r8   # useless but seems to save 0.1 cycles
        subq %r9,%r11
        mulq %r9
        addq %r10,%rax
        adcq %r11,%rdx
        cmovncq %r9,%rax
        addq %rax,%rdx
        movq %rdx,%rbx

        movq %rbx,%rax
        mulq %r14
        movq montgomery_inv_n(%rip),%r8
        movq %rax,%r10
        movq %rdx,%r11
        mulq %r8
        movq %r15,%r9
        xorq %r8,%r8   # useless but seems to save 0.1 cycles
        subq %r9,%r11
        mulq %r9
        addq %r10,%rax
        adcq %r11,%rdx
        cmovncq %r9,%rax
        addq %rax,%rdx
        movq %rdx,(%rdi)

        movq %r14,%rax
        xorq %r8,%r8
        subq %rbx,%rax
        cmovcq %r15,%r8
        addq %r8,%rax
        movq %rax,%r14

        movl    $mm_b, %ecx
        mulq (%ecx)
        movq montgomery_inv_n(%rip),%r8
        movq %rax,%r10
        movq %rdx,%r11
        mulq %r8
        movq %r15,%r9
        xorq %r8,%r8   # useless but seems to save 0.1 cycles
        subq %r9,%r11
        mulq %r9
        addq %r10,%rax
        adcq %r11,%rdx
        cmovncq %r9,%rax
        addq %rax,%rdx

        movq %rdx,%rax
        xorq %r8,%r8
        subq %r15,%rax
        addq %rbx,%rax
        cmovncq %r15,%r8
        addq %r8,%rax

        mulq %r14
        movq montgomery_inv_n(%rip),%r8
        movq %rax,%r10
        movq %rdx,%r11
        mulq %r8
        movq %r15,%r9
        xorq %r8,%r8   # useless but seems to save 0.1 cycles
        subq %r9,%r11
        mulq %r9
        addq %r10,%rax
        adcq %r11,%rdx
        cmovncq %r9,%rax
        addq %rax,%rdx
        movq %rdx,(%esi)

	popq %r15
	popq %r14
	popq %r13
	popq %r12
	popq %rbx
	ret
.endif
