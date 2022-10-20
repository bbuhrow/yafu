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
dnl asm_zero128(a): a=0
dnl asm_copy128(b,a): b=a
dnl asm_sub128_3(c,a,b): c=a-b mod N
dnl asm_half128(a): a/=2 mod N
dnl asm_sub_n128(b,a): b-=a  mod 2^128
dnl asm_diff128(c,a,b): c=|a-b|
dnl asm_add128(b,a): b+=a mod N
dnl asm_add128_ui(b,a): b+=a mod N, a is ulong
dnl asm_mulm128(): prod=f1*f2 mod N
dnl function_head(asm_sqm128)
dnl

#include "underscore.h"

.comm montgomery_modulo_n,8
.comm montgomery_inv_n,8
dnl
define(rop,%rdi)dnl
define(op,%rsi)dnl
define(op1,%rsi)dnl
define(op2,%rdx)dnl
dnl asm_zero128(a): a=0
function_head(asm_zero128)
	xorq %rax,%rax
	movq %rax,(rop)
	movq %rax,8(rop)
	ret
dnl asm_copy128(b,a): b=a
function_head(asm_copy128)
	movq (op),%rax
	movq 8(op),%rdx
	movq %rax,(rop)
	movq %rdx,8(rop)
	ret
dnl asm_sub128_3(c,a,b): c=a-b mod N
function_head(asm_sub128_3)
	movq (op1),%rax
	movq 8(op1),%rcx
	xorq %r8,%r8
	subq (op2),%rax

dnl smjs	movq montgomery_modulo_n,%r10
	movq montgomery_modulo_n(%rip),%r10

	movq $0,%r9
	sbbq 8(op2),%rcx
	cmovcq (%r10),%r8
	cmovcq 8(%r10),%r9
	addq %r8,%rax
	movq %rax,(rop)
	adcq %r9,%rcx
	movq %rcx,8(rop)
	ret
dnl asm_half128(a): a/=2 mod N
function_head(asm_half128)
	movq (rop),%rax
	testq $1,%rax
	movq 8(rop),%rdx
	jnz half_odd128
	shrq $1,%rdx
	rcrq $1,%rax
	movq %rdx,8(rop)
	movq %rax,(rop)
	ret
dnl a is odd, compute (a+N)/2
half_odd128:
dnl smjs	movq montgomery_modulo_n,%r8
	movq montgomery_modulo_n(%rip),%r8
	addq (%r8),%rax
	adcq 8(%r8),%rdx
	rcrq $1,%rdx
	rcrq $1,%rax
	movq %rdx,8(rop)
	movq %rax,(rop)
	ret
dnl asm_sub_n128(b,a): b-=a  mod 2^128???????????????????
function_head(xasm_sub_n128)
	movq (rop),%rax
	subq (op),%rax
	movq %rax,(rop)
	ret
dnl asm_add128(b,a): b+=a mod N
function_head(asm_add128)
	movq (rop),%rax
	movq 8(rop),%rdx
dnl smjs	movq montgomery_modulo_n,%r11
	movq montgomery_modulo_n(%rip),%r11

	movq (%r11),%r8
	movq 8(%r11),%r9
	subq %r8,%rax
	sbbq %r9,%rdx
	addq (op),%rax
	adcq 8(op),%rdx
	movq $0,%r10
	movq $0,%r11
	cmovncq %r8,%r10
	cmovncq %r9,%r11
	addq %r10,%rax
	adcq %r11,%rdx
	movq %rax,(rop)
	movq %rdx,8(rop)
	ret
define(res0,%r11)dnl
define(res1,%r13)dnl
define(res2,%r12)dnl
define(res3,%r11)dnl
define(aux,%r14)dnl
define(h,%r8)dnl
define(h2,%r9)dnl
define(n,%r10)dnl
define(n0,%rcx)dnl
define(n1,%r10)dnl
dnl asm_mulm128(): prod=f1*f2 mod N
function_head(asm_mulm128)
dnl first multiplication
	movq (op2),h
	movq 8(op2),h2
	movq (op1),%rax
	mulq h
	movq %r13,-16(%rsp)
	movq %r14,-24(%rsp)
dnl smjs	movq montgomery_modulo_n,n
	movq montgomery_modulo_n(%rip),n
dnl smjs	movq montgomery_inv_n,aux
	movq montgomery_inv_n(%rip),aux

	movq %rax,res0
	movq %rdx,res1
	movq 8(op1),%rax
	mulq h
	imulq res0,aux
	movq %r12,-8(%rsp)
	movq $0,res2
	movq (n),n0
	movq 8(n),n1
	addq %rax,res1
	adcq %rdx,res2
dnl first reduction
	movq n0,%rax
	mulq aux
	addq %rax,res0
	adcq %rdx,res1
	adcq $0,res2
	adcq $0,res3   # res3=res0=0
	movq n1,%rax
	mulq aux
	addq %rax,res1
	adcq %rdx,res2
	adcq $0,res3
dnl second multiplication
	movq (op1),%rax
	mulq h2
dnl smjs	movq montgomery_inv_n,aux
	movq montgomery_inv_n(%rip),aux
	addq %rax,res1
	adcq %rdx,res2
	adcq $0,res3
	movq 8(op1),%rax
	mulq h2
	imulq res1,aux
	subq n0,res2
	sbbq n1,res3
	addq %rax,res2
	adcq %rdx,res3
dnl second reduction
	movq n0,%rax
	mulq aux
	addq %rax,res1
	adcq %rdx,res2
	adcq $0,res3
	movq n1,%rax
	mulq aux
	movq -24(%rsp),%r14
	addq %rax,res2
	adcq %rdx,res3
dnl result-N in (res2,res3); add N if no carry
	cmovcq res1,n0
	cmovcq res1,n1
	addq res2,n0
	adcq res3,n1
	movq -16(%rsp),%r13
	movq -8(%rsp),%r12
	movq n0,(rop)
	movq n1,8(rop)

	ret
undefine(`rop')dnl
undefine(`op')dnl
undefine(`op1')dnl
undefine(`op2')dnl
undefine(`res0')dnl
undefine(`res1')dnl
undefine(`res2')dnl
undefine(`res3')dnl
undefine(`op2x')dnl
undefine(`h')dnl
undefine(`h2')dnl
undefine(`n0')dnl
undefine(`n1')dnl
