# Copyright (C) 2004 Jens Franke, T.Kleinjung
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# Written by T. Kleinjung
# Modifications by J. Franke

define(relptr,%rdi)dnl
define(minus,%rsi)dnl
define(qx_arg,%rdx)dnl
define(qx,%r8)dnl
define(qxd,%r8d)dnl
define(nr,%r9)dnl
define(nrw,%r9w)dnl
define(nr1,%rsi)dnl
define(nr1w,%si)dnl
define(aux1,%r10)dnl
define(aux2,%r11)dnl
define(aux2d,%r11d)dnl
define(aux3,%rbx)dnl
define(aux3d,%ebx)dnl
define(aux4,%r12)dnl
define(aux4d,%r12d)dnl
define(aux5,%r13)dnl
define(aux5w,%r13w)dnl

	.align 16
.globl asm_td
	.type    asm_td,@function
# mm0: p,p
# mm1: inv,inv
# mm2: s1,s2
# mm3: ind,ind
# mm4: computing
# mm5: 0

asm_td:
	pushq %rbx
	pushq %r12
	pushq %r13
	movq qx_arg,qx
	movzwq (relptr),%rax      # ind
	movq %rax,%rcx
	shlq $16,%rax
	orq %rcx,%rax            # ind,ind
	movq %rax,%rcx
	shlq $32,%rax
	orq %rcx,%rax            # ind,ind,ind,ind
	movd %rax,%xmm5
	movd %rax,%xmm3
	pslldq $8,%xmm5
	paddw %xmm5,%xmm3        # ind,ind,ind,ind,ind,ind,ind,ind
	pxor %xmm5,%xmm5

	movzwq 8(relptr),nr      # nr

dnl Calculate product of primes found by sieving in %rax
#prod:
#	movq $mpqs_FB,aux1
#	movzwq mpqs_nFBk_1(%rip),%rdx
#	addq %rdx,%rdx
#	addq %rdx,%rdx
#	subq %rdx,aux1
#	movq $0,%rcx
#	movq $1,%rax
#prodloop:
#	cmpq nr,%rcx
#	jnc prodend
#	movzwq 10(relptr,%rcx,2),%rdx
#	movzwq (aux1,%rdx,4),%rdx
#	mulq %rdx
#	incq %rcx
#	jmp prodloop
#	
#prodend:

	movzwq mpqs_td_begin(%rip),%rcx

	movq $mpqs_FB_inv_info,aux3
	movq $mpqs_FB_start,aux2

	movw mpqs_nFBk_1(%rip),aux5w
	addw $4,aux5w

dnl Measure time for MMX Loop in zeit0
dnl	movq %rdi,aux5
dnl	xorq %rdi,%rdi
dnl	call zeitA
dnl	movq aux5,%rdi

dnl first consider primes nr 1, 2 and 3
	movaps (aux3),%xmm0
	movaps 16(aux3),%xmm1
	movaps (aux2),%xmm2
	movaps %xmm0,%xmm4
	psubw %xmm2,%xmm4
	paddw %xmm3,%xmm4        # ind+p-s1,ind+p-s2,...
	pmullw %xmm1,%xmm4
	movaps 16(aux2),%xmm2
	leaq 32(aux3),aux3
	leaq 16(aux2),aux2
	pmulhuw %xmm0,%xmm4
	movaps (aux3),%xmm0
	movaps 16(aux3),%xmm1
	psubw %xmm3,%xmm2
	pcmpeqw %xmm5,%xmm4
	pmovmskb %xmm4,aux4d
	andl $0xfff0,aux4d   # only considering primes nr 1,2,3
	movaps %xmm0,%xmm4
	psubw %xmm2,%xmm4
	jz loop2

dnl found divisor
	movl aux4d,%eax
	subw $3,aux5w       # aux5w=mpqs_nFBk_1+1
	shrl $2,%eax
	orl %eax,aux4d

	shrl $5,aux4d
	movw aux5w,10(relptr,nr,2)
	adcq $0,nr
	incw aux5w          # aux5w=mpqs_nFBk_1+2
	shrl $4,aux4d
	movw aux5w,10(relptr,nr,2)
	adcq $0,nr
	incw aux5w          # aux5w=mpqs_nFBk_1+3
	shrl $4,aux4d
	movw aux5w,10(relptr,nr,2)
	adcq $0,nr
	incw aux5w          # aux5w=mpqs_nFBk_1+4

loop2:
	addw $4,aux5w
	subq $4,%rcx
	leaq 32(aux3),aux3
	leaq 16(aux2),aux2
	jz prod

	pmullw %xmm1,%xmm4
	movaps (aux2),%xmm2
	psubw %xmm3,%xmm2
	pmulhuw %xmm0,%xmm4
	movaps (aux3),%xmm0
	movaps 16(aux3),%xmm1
	pcmpeqw %xmm5,%xmm4
	pmovmskb %xmm4,aux4d
	testl aux4d,aux4d
	movaps %xmm0,%xmm4
	psubw %xmm2,%xmm4
	jz loop2
dnl found divisor
	movl aux4d,%eax
	subw $4,aux5w
	shrl $2,%eax
	orl %eax,aux4d  # significant bits at position 0,4,8,12

	shrl $1,aux4d
	movw aux5w,10(relptr,nr,2)
	adcq $0,nr
	incw aux5w
	shrl $4,aux4d
	movw aux5w,10(relptr,nr,2)
	adcq $0,nr
	incw aux5w
	shrl $4,aux4d
	movw aux5w,10(relptr,nr,2)
	adcq $0,nr
	incw aux5w
	shrl $4,aux4d
	movw aux5w,10(relptr,nr,2)
	adcq $0,nr
	incw aux5w

	jmp loop2

prod:
	movq $mpqs_FB,aux1
	movzwq mpqs_nFBk_1(%rip),%rdx
	addq %rdx,%rdx
	addq %rdx,%rdx
	subq %rdx,aux1
	movq $0,%rcx
	movq $1,%rax
prodloop:
	cmpq nr,%rcx
	jnc prodend
	movzwq 10(relptr,%rcx,2),%rdx
	movzwq (aux1,%rdx,4),%rdx
	mulq %rdx
	incq %rcx
	jmp prodloop
prodend:
	movq %rax,aux3

dnl End of measuring MMX Loop in zeit0
dnl	movq %rdi,aux5
dnl	xorq %rdi,%rdi
dnl	call zeitB
dnl	movq aux5,%rdi


	testq $1,minus
dnl nr1 and minus are the same register
	movq nr,nr1
	jz posloop
	movw $0,10(relptr,nr,2)
	incq nr

posloop:
	testq $0x00000001,qx
	jnz odd
	incq nr
	cmpq $27,nr
	jnc gotonext
	movw mpqs_nFBk_1(%rip),%cx
	movw %cx,8(relptr,nr,2)
	shrq $1,qx
	jmp posloop

odd:

dnl	cmpl $0,16(%esp)
dnl	jnz division
dnl	cmpl 12(%esp),%edx
dnl	jnc gotonext

dnl	movq %rdi,aux5
dnl	xorq %rdi,%rdi
dnl	call zeitA
dnl	movq aux5,%rdi

division:
	movq aux3,aux2
	movq aux3,%rdx
	movq aux3,aux5
	movq aux3,%rcx
	andq $0xff,%rdx
	shrq $32,aux5
	shlq $32,%rcx
	shrq $1,%rdx
	incq %rcx
	movl aux3d,%eax
	cmpq %rcx,qx
	movzbl mpqs_256_inv_table(%rdx),aux2d
	adcq $0,aux5
	mull aux2d
	testq aux5,aux5
	jz gotonext
	andl $0xff00,%eax
	mull aux2d
	subl %eax,aux2d
	movl aux3d,%eax
	mull aux2d
	andl $0xffff0000,%eax
	mull aux2d
	subl %eax,aux2d
	movl qxd,%eax
	mull aux2d
	movl %eax,qxd

dnl	movq %rdi,aux5
dnl	xorq %rdi,%rdi
dnl	call zeitB
dnl	movq aux5,%rdi

# trial divison of sieved primes
	movq $mpqs_FB_inv,aux2
	movzwq mpqs_nFBk_1,%rcx
	addq %rcx,%rcx
	addq %rcx,%rcx
	subq %rcx,aux2
	movq $0,%rcx
tdloop:
	testq nr1,nr1
	movl qxd,%eax
	jz tdend
	movzwq 8(relptr,nr1,2),aux5  # bx: ii
	decq nr1
dnl Value of aux1 is still mpqs_FB-4*mpqs_nFBk_1
	movzwl (aux1,aux5,4),%ecx   # cx: p
	movl (aux2,aux5,4),aux3d  # edx: inv
divloop:
	mull aux3d
	movl %eax,aux4d
	mull %ecx
	testl %edx,%edx
	jnz tdloop
	cmpw $27,nrw
	jnc gotonext
	movl aux4d,%eax
	movw aux5w,10(relptr,nr,2)
	incq nr
	movl aux4d,qxd
	jmp divloop


tdend:
# trial division of mpqs_FBk-primes
	cmpl $1,qxd
	jz end

	movq $mpqs_FBk_inv,aux2
	movq $mpqs_FBk,aux1
	movzwq mpqs_nFBk(%rip),nr1
	incq nr1
tdloopk:
	decq nr1
	movl qxd,%eax
	cmpq $0,nr1
	jz tdendk
	movzwl -2(aux1,nr1,2),%ecx   # cx: p
	movl -4(aux2,nr1,4),aux3d  # aux3: inv
divloopk:
	mull aux3d
	movl %eax,aux4d        # rr
	mull %ecx
	testl %edx,%edx
	jnz tdloopk
	cmpw $27,nrw
	jnc gotonext
	movl aux4d,%eax
	movw nr1w,10(relptr,nr,2)
	incq nr
	movl aux4d,qxd
	jmp divloopk

tdendk:
# trial division of mpqs_FB_Adiv-primes
	cmpl $1,%eax
	jz end

	movq $mpqs_FB_A_inv,aux2
	movq $mpqs_Adiv_all,aux1
        movw mpqs_nFB(%rip),aux5w
        addw mpqs_nFBk(%rip),aux5w
	movzwq mpqs_nAdiv_total,nr1
	incq nr1
tdloopa:
	decq nr1
	movl qxd,%eax
	cmpq $0,nr1
	jz tdenda
	movzwl -2(aux1,nr1,2),%ecx   # cx: p
	movl -4(aux2,nr1,4),aux3d  # aux3: inv
divloopa:
	mull aux3d
	movl %eax,aux4d        # rr
	mull %ecx
	testl %edx,%edx
	jnz tdloopa
	movl aux4d,%eax
	cmpw $27,nrw
	jnc gotonext
	addw nr1w,aux5w
	movw aux5w,10(relptr,nr,2)
	incq nr
	subw nr1w,aux5w
	movl aux4d,qxd
	jmp divloopa

tdenda:
	

end:
	movw nrw,8(relptr)
	movl qxd,%eax
	emms	
	popq %r13
	popq %r12
	popq %rbx
	ret

gotonext:
	xorq %rax,%rax
	emms
	popq %r13
	popq %r12
	popq %rbx
	ret

Schlendrian:
	call abort
