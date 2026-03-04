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
dnl asm_zero192(a): a=0
dnl asm_copy192(b,a): b=a
dnl asm_sub192_3(c,a,b): c=a-b mod N
dnl asm_half192(a): a/=2 mod N
dnl asm_sub_n192(b,a): b-=a  mod 2^192
dnl asm_diff192(c,a,b): c=|a-b|
dnl asm_add192(b,a): b+=a mod N
dnl asm_add192_ui(b,a): b+=a mod N, a is ulong
dnl asm_mulm192(): prod=f1*f2 mod N
dnl function_head(asm_sqm192)
dnl

#include "underscore.h"

.comm montgomery_modulo_n,8
.comm montgomery_inv_n,8
dnl
define(rop,%rdi)dnl
define(op,%rsi)dnl
define(op1,%rsi)dnl
define(op2,%rdx)dnl
dnl asm_zero192(a): a=0
function_head(asm_zero192)
	xorq %rax,%rax
	movq %rax,(rop)
	movq %rax,8(rop)
	movq %rax,16(rop)
	ret
dnl asm_copy192(b,a): b=a
function_head(asm_copy192)
	movq (op),%rax
	movq 8(op),%rdx
	movq %rax,(rop)
	movq %rdx,8(rop)
	movq 16(op),%rax
	movq %rax,16(rop)
	ret
dnl asm_sub192_3(c,a,b): c=a-b mod N
function_head(asm_sub192_3)
	movq (op1),%rax
	movq 8(op1),%rcx
	xorq %r8,%r8
	subq (op2),%rax
dnl smjs	movq montgomery_modulo_n,%r11
	movq montgomery_modulo_n(%rip),%r11

	movq $0,%r9
	sbbq 8(op2),%rcx
	movq 16(op1),%r10
	sbbq 16(op2),%r10
	movq $0,%rdx
	cmovcq (%r11),%r8
	cmovcq 8(%r11),%r9
	cmovcq 16(%r11),%rdx
	addq %r8,%rax
	movq %rax,(rop)
	adcq %r9,%rcx
	movq %rcx,8(rop)
	adcq %r10,%rdx
	movq %rdx,16(rop)
	ret
dnl asm_half192(a): a/=2 mod N
function_head(asm_half192)
	movq (rop),%rax
	testq $1,%rax
	movq 8(rop),%rcx
	movq 16(rop),%rdx
	jnz half_odd192
	shrq $1,%rdx
	rcrq $1,%rcx
	rcrq $1,%rax
	movq %rdx,16(rop)
	movq %rcx,8(rop)
	movq %rax,(rop)
	ret
dnl a is odd, compute (a+N)/2
half_odd192:
dnl smjs	movq montgomery_modulo_n,%r8
	movq montgomery_modulo_n(%rip),%r8

	addq (%r8),%rax
	adcq 8(%r8),%rcx
	adcq 16(%r8),%rdx
	rcrq $1,%rdx
        rcrq $1,%rcx
	rcrq $1,%rax
	movq %rdx,16(rop)
	movq %rcx,8(rop)
	movq %rax,(rop)
	ret
dnl asm_sub_n192(b,a): b-=a  mod 2^192???????????????????
function_head(xasm_sub_n192)
	movq (rop),%rax
	subq (op),%rax
	movq %rax,(rop)
	ret
dnl asm_add192(b,a): b+=a mod N
function_head(asm_add192)
	movq (rop),%rax
	movq 8(rop),%rcx
dnl smjs	movq montgomery_modulo_n,%r11
	movq montgomery_modulo_n(%rip),%r11

	movq (%r11),%r8
	movq 8(%r11),%r9
	subq %r8,%rax
	sbbq %r9,%rcx
	movq 16(rop),%rdx
	movq 16(%r11),%r10
	sbbq %r10,%rdx
	addq (op),%rax
	adcq 8(op),%rcx
	adcq 16(op),%rdx
	movq $0,%r11
	cmovcq %r11,%r8
	cmovcq %r11,%r9
	cmovcq %r11,%r10
	addq %r8,%rax
	adcq %r9,%rcx
	adcq %r10,%rdx
	movq %rax,(rop)
	movq %rcx,8(rop)
	movq %rdx,16(rop)
	ret
define(res0,%r9)dnl
define(res1,%r10)dnl
define(res2,%r14)dnl
define(res3,%r11)dnl
define(res4,%r9)dnl
define(res5,%r10)dnl
define(op2x,%rcx)dnl
define(aux,%r13)dnl
define(h,%r8)dnl
define(n,%rbx)dnl
define(n0,%r12)dnl
define(n1,%r15)dnl
define(n2,%rbx)dnl  same as n
dnl asm_mulm192(): prod=f1*f2 mod N
function_head(asm_mulm192)
dnl first multiplication
	movq (op2),h
	movq (op1),%rax
	movq op2,op2x
	mulq h

	movq %rbx,-16(%rsp)
	movq %r12,-24(%rsp)
	movq %r13,-32(%rsp)

	movq %rax,res0
	movq %rdx,res1
	movq 8(op1),%rax
	mulq h
dnl smjs	movq montgomery_inv_n,aux
	movq montgomery_inv_n(%rip),aux

	movq %r14,-40(%rsp)
	movq %r15,-48(%rsp)
	movq $0,res2
	movq $0,res3
	addq %rax,res1
	adcq %rdx,res2
	movq 16(op1),%rax
	mulq h
	imulq res0,aux
dnl smjs	movq montgomery_modulo_n,n
	movq montgomery_modulo_n(%rip),n

	movq (n),n0
	movq 8(n),n1
	movq 16(n),n2
	addq %rax,res2
	adcq %rdx,res3
dnl first reduction
	movq n0,%rax
	mulq aux
	addq %rax,res0
	adcq %rdx,res1
	adcq $0,res2
	adcq $0,res3
	adcq $0,res4  # res4=res0=0
	movq n1,%rax
	mulq aux
	addq %rax,res1
	adcq %rdx,res2
	adcq $0,res3
	adcq $0,res4
	movq n2,%rax
	mulq aux
	movq 8(op2x),h
	addq %rax,res2
	adcq %rdx,res3
	adcq $0,res4
dnl second multiplication
	movq (op1),%rax
	mulq h
	addq %rax,res1
	adcq %rdx,res2
	adcq $0,res3
	adcq $0,res4
	movq 8(op1),%rax
	mulq h
	addq %rax,res2
	adcq %rdx,res3
	adcq $0,res4
	movq 16(op1),%rax
	mulq h
dnl smjs	movq montgomery_inv_n,aux
	movq montgomery_inv_n(%rip),aux

	imulq res1,aux
	addq %rax,res3
	adcq %rdx,res4
dnl second reduction
	movq n0,%rax
	mulq aux
	addq %rax,res1
	adcq %rdx,res2
	adcq $0,res3
	adcq $0,res4
	adcq $0,res5  # res5=res1=0
	movq n1,%rax
	mulq aux
	addq %rax,res2
	adcq %rdx,res3
	adcq $0,res4
	adcq $0,res5
	movq n2,%rax
	mulq aux
	movq 16(op2x),h
	addq %rax,res3
	adcq %rdx,res4
	adcq $0,res5
dnl third multiplication
	movq (op1),%rax
	mulq h
	addq %rax,res2
	adcq %rdx,res3
	adcq $0,res4
	adcq $0,res5
	movq 8(op1),%rax
	mulq h
	addq %rax,res3
	adcq %rdx,res4
	adcq $0,res5
	movq 16(op1),%rax
	mulq h
dnl smjs	movq montgomery_inv_n,aux
	movq montgomery_inv_n(%rip),aux
	imulq res2,aux
	subq n0,res3
	sbbq n1,res4
	sbbq n2,res5
	addq %rax,res4
	adcq %rdx,res5
dnl third reduction
	movq n0,%rax
	mulq aux
	addq %rax,res2
	adcq %rdx,res3
	adcq $0,res4
	adcq $0,res5
	movq n1,%rax
	mulq aux
	addq %rax,res3
	adcq %rdx,res4
	adcq $0,res5
	movq n2,%rax
	mulq aux
	addq %rax,res4
	adcq %rdx,res5
dnl result-N in (res3,res4,res5); add N if no carry
	cmovcq res2,n0
	cmovcq res2,n1
	cmovcq res2,n2
	movq -40(%rsp),%r14
	movq -32(%rsp),%r13
	addq n0,res3
	adcq n1,res4
	adcq n2,res5
	movq -48(%rsp),%r15
	movq -24(%rsp),%r12
	movq -16(%rsp),%rbx

	movq res3,(rop)
	movq res4,8(rop)
	movq res5,16(rop)
	ret
undefine(`rop')dnl
undefine(`op')dnl
undefine(`op1')dnl
undefine(`op2')dnl
undefine(`res0')dnl
undefine(`res1')dnl
undefine(`res2')dnl
undefine(`res3')dnl
undefine(`res4')dnl
undefine(`res5')dnl
undefine(`op2x')dnl
undefine(`aux')dnl
undefine(`h')dnl
undefine(`n')dnl
undefine(`n0')dnl
undefine(`n1')dnl
undefine(`n2')dnl
