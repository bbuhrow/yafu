# Copyright (C) 2004 Jens Franke, T.Kleinjung
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# Written by T. Kleinjung
# Modifications by J. Franke

#include "underscore.h"

dnl n=asm3_tdsieve(fb,fbs,&(buf[0]),mpqs3_nFBk_1+i-1);

define(FB,%rdi)dnl
define(FBs,%rsi)dnl
define(buf,%rdx)dnl
define(disp,%rcx)dnl
define(disp_word,%cx)dnl
define(pend,%rbp)dnl
define(aux1,%rbx)dnl
define(aux2,%rax)dnl
define(sloc1,%r8)dnl
define(sloc2,%r9)dnl
define(prime,%r10)dnl
define(prime2,%r11)dnl
define(sieve_interval,%r12)dnl
define(sieveend,%r13)dnl
define(sievebound,%r14)dnl
define(sl,%r15)dnl
define(acc,%rbx)dnl
define(accb,%bl)dnl

function_head(asm3_tdsieve)
	pushq %rbx
	pushq %rbp
	pushq %r12
	pushq %r13
	pushq %r14
	pushq %r15

dnl smjs	movzwq mpqs3_sievelen(%rip),sl
	movzwq mpqs3_sievelen(%rip),sl

	movq sl,pend
	shrq $2,pend
dnl smjs	movq mpqs3_sievearray(%rip),sieve_interval
	movq mpqs3_sievearray(%rip),sieve_interval

	movq sieve_interval,sieveend
	addq sl,sieveend
tds4loop:
	movzwq (FB),prime
	cmpq pend,prime
	jnc tds4loopend

	leaq 4(FB),FB
	incq disp
	movzwq (FBs),sloc1
	movzwq 2(FBs),sloc2
	leaq 4(FBs),FBs

	cmpq sloc1,sloc2
	movq sloc1,aux1
	movq sloc2,aux2
	cmovcq aux1,sloc2
	cmovcq aux2,sloc1    # sloc1<sloc2

	addq sieve_interval,sloc1
	addq sieve_interval,sloc2

	movq sieveend,sievebound
	leaq (prime,prime),prime2
	subq prime,sievebound
	subq prime2,sievebound
	xorq acc,acc

tds4innerloop:
	movzbq (sloc1),%rax
	addq %rax,acc
	shlq $8,acc
	movzbq (sloc2),%rax
	addq %rax,acc
	shlq $8,acc
	movzbq (sloc1,prime),%rax
	addq %rax,acc
	shlq $8,acc
	movzbq (sloc2,prime),%rax
	addq %rax,acc
	shlq $8,acc
	addq prime2,sloc1
	addq prime2,sloc2
	movzbq (sloc1),%rax
	addq %rax,acc
	shlq $8,acc
	movzbq (sloc2),%rax
	addq %rax,acc
	shlq $8,acc
	movzbq (sloc1,prime),%rax
	addq %rax,acc
	shlq $8,acc
	movzbq (sloc2,prime),%rax
	addq %rax,acc
	leaq (sloc1,prime2),sloc1
	leaq (sloc2,prime2),sloc2

	jnz tds4store
tds4store_return:
	cmpq sievebound,sloc2
	jc tds4innerloop

	addq prime2,sievebound
	cmpq sievebound,sloc2
	jnc tds4check

	movzbq (sloc1),%rax
	addq %rax,acc
	shlq $8,acc
	movzbq (sloc2),%rax
	addq %rax,acc
	shlq $8,acc
	movzbq (sloc1,prime),%rax
	addq %rax,acc
	shlq $8,acc
	movzbq (sloc2,prime),%rax
	addq %rax,acc
	shlq $8,acc
	
	addq prime2,sloc1
	addq prime2,sloc2
tds4check:
	cmpq sieveend,sloc1
	cmovncq sieveend,sloc1
	cmpq sieveend,sloc2
	cmovncq sieveend,sloc2
	movzbq (sloc1),%rax
	addq %rax,acc
	shlq $8,acc
	addq prime,sloc1
	movzbq (sloc2),%rax
	addq %rax,acc
	shlq $8,acc
	cmpq sieveend,sloc1
	cmovncq sieveend,sloc1
	movzbq (sloc1),%rax
	addq %rax,acc
	jz tds4loop

tds4store0loop:
	movq $0xff,%rax
	andq acc,%rax
	jz tds4nostore0
	leaq (buf,%rax,8),prime2
	movq (prime2),%rax
	movw disp_word,(%rax)
	addq $2,%rax
	movq %rax,(prime2)
tds4nostore0:
	shrq $8,acc
	jnz tds4store0loop

	jmp tds4loop

tds4store:
# process entries
.if 0
 pushq %rdi
 pushq %rdx
 movq $22,%rdi
 call zeitA
 popq %rdx
 popq %rdi
.endif

tds4storeloop:
	movq $0xff,%rax
	andq acc,%rax
	jz tds4nostore
	leaq (buf,%rax,8),prime2
	movq (prime2),%rax
	movw disp_word,(%rax)
	addq $2,%rax
	movq %rax,(prime2)
tds4nostore:
	shrq $8,acc
	jnz tds4storeloop

	leaq (prime,prime),prime2
.if 0
 pushq %rdi
 pushq %rdx
 movq $22,%rdi
 call zeitB
 popq %rdx
 popq %rdi
.endif
	jmp tds4store_return

tds4loopend:


.if 0
	movzbl 3(%esp,sv),%eax
	decl sv
	leal (p,%eax,4),svlen
	movl (svlen),%eax
	movw s1w,(%eax)
	leal 2(%eax),%eax
	movl %eax,(svlen)
	jnz tds4store

        movl fbaux,s1    # retrieve p and s1
        movzwl -4(s1),p
	movl s1aux,s1
	jmp tds4innerloop

tds4innerloopend:
	leal (iend,p,2),iend
	xorl %eax,%eax
	movl %eax,byte1
	movl %eax,byte5
	cmpl iend,s3
	jnc tds4check2

	movb (sv,s1),%al
	addl p,s1
	movb %al,byte1
	movb %al,%ah
	movb (sv,s1),%al
	addl p,s1
	movb %al,byte2
	orb %al,%ah

	movb (sv,s3),%al
	addl p,s3
	movb %al,byte5
	orb %al,%ah
	movb (sv,s3),%al
	addl p,s3
	movb %al,byte6
	orb %al,%ah
tds4check2:
	cmpl svlen,s3
	jnc tds4check1

	movb (sv,s1),%al
	addl p,s1
	movb %al,byte3
	orb %al,%ah

	movb (sv,s3),%al
	addl p,s3
	movb %al,byte7
	orb %al,%ah
tds4check1:
	cmpl svlen,s1
	cmovncl svlen,s1
	movb (sv,s1),%al
	addl p,s1
	movb %al,byte4
	orb %al,%ah
	cmpl svlen,s3
	cmovncl svlen,s3
	movb (sv,s3),%al
	addl p,s3
	movb %al,byte8
	orb %al,%ah

	jz tds4loop
# compute mpqs3_nFBk_1+i:
	movl fbaux,s1
	subl fbptr,s1
	shrl $2,s1
	addl nr,s1
	movl buffer,p
# skip zeros in byte1-8:
	xorl sv,sv
	movb byte1,%al
	addb $0xff,%al
	adcl $0,sv
	movb byte2,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv
	movb byte3,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv
	movb byte4,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv
	movb byte5,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv
	movb byte6,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv
	movb byte7,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv
	movb byte8,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv

# process entries
tds4storex:
	movzbl 3(%esp,sv),%eax
	decl sv
	leal (p,%eax,4),svlen
	movl (svlen),%eax
	movw s1w,(%eax)
	leal 2(%eax),%eax
	movl %eax,(svlen)
	jnz tds4storex

	movl mpqs3_sievearray,sv
	movl mpqs3_sievelen,svlen
	jmp tds4loop

tds4loopend:
	movl fbsaux,fbs   # fb is ok

	movl svlen,%eax
	xorl %edx,%edx
	movl $3,s1
	divl s1
	movl %eax,pend
	movl mpqs3_sievearray,sv   # sv is edx
tds3loop:
	movzwl (fb),p
	cmpl pend,p
	jnc tds3loopend

	movzwl (fbs),s1
	leal 4(fb),fb
	movb (sv,s1),%al
	addl p,s1
	movb %al,byte1
	movb %al,%ah
	movb (sv,s1),%al
	addl p,s1
	movb %al,byte2
	orb %al,%ah
	movb (sv,s1),%al
	addl p,s1
	movb %al,byte3
	cmpl svlen,s1
	cmovncl svlen,s1
	orb %al,%ah
	movb (sv,s1),%al
	movzwl 2(fbs),s2
	movb %al,byte4
	orb %al,%ah

	movb (sv,s2),%al
	addl p,s2
	movb %al,byte5
	orb %al,%ah
	movb (sv,s2),%al
	addl p,s2
	movb %al,byte6
	orb %al,%ah
	movb (sv,s2),%al
	addl p,s2
	movb %al,byte7
	orb %al,%ah
	cmpl svlen,s2
	cmovncl svlen,s2
	movb (sv,s2),%al
	leal 4(fbs),fbs
	orb %al,%ah
	movb %al,byte8

	jz tds3loop
# compute mpqs3_nFBk_1+i:
	movl fb,s1
	subl fbptr,s1
	shrl $2,s1
	addl nr,s1
	movl buffer,p
# skip zeros in byte1-8:
	xorl sv,sv
	movb byte1,%al
	addb $0xff,%al
	adcl $0,sv
	movb byte2,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv
	movb byte3,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv
	movb byte4,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv
	movb byte5,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv
	movb byte6,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv
	movb byte7,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv
	movb byte8,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv

# process entries
tds3store:
	movzbl 3(%esp,sv),%eax
	decl sv
	leal (p,%eax,4),svlen
	movl (svlen),%eax
	movw s1w,(%eax)
	leal 2(%eax),%eax
	movl %eax,(svlen)
	jnz tds3store

	movl mpqs3_sievearray,sv
	movl mpqs3_sievelen,svlen
	jmp tds3loop

tds3loopend:
	movl svlen,p
	shrl $1,p
	movl p,pend
tds2loop:
	movzwl (fb),p
	cmpl pend,p
	jnc tds2loopend

	movzwl (fbs),s1
	leal 4(fb),fb
	movb (sv,s1),%al
	addl p,s1
	movb %al,byte1
	movb %al,%ah
	movb (sv,s1),%al
	addl p,s1
	movb %al,byte2
	cmpl svlen,s1
	cmovncl svlen,s1
	orb %al,%ah
	movb (sv,s1),%al
	movzwl 2(fbs),s2
	movb %al,byte3
	orb %al,%ah

	movb (sv,s2),%al
	addl p,s2
	movb %al,byte4
	orb %al,%ah
	movb (sv,s2),%al
	addl p,s2
	movb %al,byte5
	orb %al,%ah
	cmpl svlen,s2
	cmovncl svlen,s2
	movb (sv,s2),%al
	leal 4(fbs),fbs
	orb %al,%ah
	movb %al,byte6

	jz tds2loop
# compute mpqs3_nFBk_1+i:
	movl fb,s1
	subl fbptr,s1
	shrl $2,s1
	addl nr,s1
	movl buffer,p
# skip zeros in byte1-6:
	xorl sv,sv
	movb byte1,%al
	addb $0xff,%al
	adcl $0,sv
	movb byte2,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv
	movb byte3,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv
	movb byte4,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv
	movb byte5,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv
	movb byte6,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv

# process entries
tds2store:
	movzbl 3(%esp,sv),%eax
	decl sv
	leal (p,%eax,4),svlen
	movl (svlen),%eax
	movw s1w,(%eax)
	leal 2(%eax),%eax
	movl %eax,(svlen)
	jnz tds2store

	movl mpqs3_sievearray,sv
	movl mpqs3_sievelen,svlen
	jmp tds2loop

tds2loopend:

tds1loop:
	movzwl (fb),p
	cmpl svlen,p
	jnc tds1loopend

	movzwl (fbs),s1
	leal 4(fb),fb
	movb (sv,s1),%al
	addl p,s1
	movb %al,byte1
	cmpl svlen,s1
	cmovncl svlen,s1
	movb %al,%ah
	movb (sv,s1),%al
	movzwl 2(fbs),s2
	movb %al,byte2
	orb %al,%ah

	movb (sv,s2),%al
	addl p,s2
	movb %al,byte3
	orb %al,%ah
	cmpl svlen,s2
	cmovncl svlen,s2
	movb (sv,s2),%al
	leal 4(fbs),fbs
	orb %al,%ah
	movb %al,byte4

	jz tds1loop
# compute mpqs3_nFBk_1+i:
	movl fb,s1
	subl fbptr,s1
	shrl $2,s1
	addl nr,s1
	movl buffer,p
# skip zeros in byte1-4:
	xorl sv,sv
	movb byte1,%al
	addb $0xff,%al
	adcl $0,sv
	movb byte2,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv
	movb byte3,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv
	movb byte4,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv

# process entries
tds1store:
	movzbl 3(%esp,sv),%eax
	decl sv
	leal (p,%eax,4),svlen
	movl (svlen),%eax
	movw s1w,(%eax)
	leal 2(%eax),%eax
	movl %eax,(svlen)
	jnz tds1store

	movl mpqs3_sievearray,sv
	movl mpqs3_sievelen,svlen
	jmp tds1loop

tds1loopend:

tds0loop:
	movzwl (fb),p
	cmpl $0xffff,p
	jz tds0loopend

	movzwl (fbs),s1
	leal 4(fb),fb
	cmpl svlen,s1
	cmovncl svlen,s1
	movb (sv,s1),%al
	movzwl 2(fbs),s2
	movb %al,byte1
	movb %al,%ah

	cmpl svlen,s2
	cmovncl svlen,s2
	movb (sv,s2),%al
	leal 4(fbs),fbs
	orb %al,%ah
	movb %al,byte2

	jz tds0loop
# compute mpqs3_nFBk_1+i:
	movl fb,s1
	subl fbptr,s1
	shrl $2,s1
	addl nr,s1
	movl buffer,p
# skip zeros in byte1-2:
	xorl sv,sv
	movb byte1,%al
	addb $0xff,%al
	adcl $0,sv
	movb byte2,%al
	movb %al,4(%esp,sv)
	addb $0xff,%al
	adcl $0,sv

# process entries
tds0store:
	movzbl 3(%esp,sv),%eax
	decl sv
	leal (p,%eax,4),svlen
	movl (svlen),%eax
	movw s1w,(%eax)
	leal 2(%eax),%eax
	movl %eax,(svlen)
	jnz tds0store

	movl mpqs3_sievearray,sv
	movl mpqs3_sievelen,svlen
	jmp tds0loop

tds0loopend:
	movl fbptr,fbs
	movl fb,%eax
	subl fbs,%eax
	shrl $1,%eax

	addl $28,%esp
	popl %ebp
	popl %ebx
	popl %esi
	popl %edi
	ret
.endif
	xorq %rax,%rax
	emms
	popq %r15
	popq %r14
	popq %r13
	popq %r12
	popq %rbp
	popq %rbx
	ret

define(FB,%rdi)dnl
define(FBs,%rsi)dnl
define(buf,%rdx)dnl
define(disp,%rcx)dnl
define(disp_word,%cx)dnl
define(pend,%rbp)dnl
define(aux1,%rbx)dnl
define(aux2,%rax)dnl
define(sloc1,%r8)dnl
define(sloc2,%r9)dnl
define(prime,%r10)dnl
define(prime2,%r11)dnl
define(sieve_interval,%r12)dnl
define(sieveend,%r13)dnl
define(sievebound,%r14)dnl
define(sl,%r15)dnl
define(acc,%rbx)dnl

undefine(`FB')dnl
undefine(`FBs')dnl
undefine(`buf')dnl
undefine(`disp')dnl
undefine(`disp_word')dnl
undefine(`pend')dnl
undefine(`aux1')dnl
undefine(`aux2')dnl
undefine(`sloc1')dnl
undefine(`sloc2')dnl
undefine(`prime')dnl
undefine(`prime2')dnl
undefine(`sieve_interval')dnl
undefine(`sieveend')dnl
undefine(`sievebound')dnl
undefine(`sl')dnl
undefine(`acc')dnl

define(relptr,%rdi)dnl
define(relptrw,%rdi)dnl
define(minus,%rsi)dnl
define(qx_arg,%rdx)dnl
define(qx0,%r8)dnl
define(qxd,%r8d)dnl
define(qx1,%r14)dnl
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
define(prod0,%r15)dnl
define(prod1,%rbp)dnl
define(NMAXDIV,$``''25)dnl

dnl mm0: p,p
dnl mm1: inv,inv
dnl mm2: s1,s2
dnl mm3: ind,ind
dnl mm4: computing
dnl mm5: 0

	.align 16
ifelse(windows,`1',
`.globl asm3_td
        .def    asm3_td; .scl    2;      .type   32;     .endef
asm3_td:'
,osx,`1',
`.globl _asm3_td
_asm3_td:'
,linux,`1',
`.globl asm3_td
	.type    asm3_td,@function
asm3_td:'
,)dnl
	pushq %rbx
	pushq %rbp
	pushq %r12
	pushq %r13
	pushq %r14
	pushq %r15
	movq 4(qx_arg),qx1
	shrq $32,qx1
	movq (qx_arg),qx0
	movzwq (relptr),%rax     # ind
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

	movzwq 12(relptr),nr      # nr

dnl Measure time for MMX Loop in zeit20
dnl	movq %rdi,aux5
dnl	movq $20,%rdi
dnl	call zeitA
dnl	movq aux5,%rdi

dnl smjs	movzwq mpqs3_td_begin(%rip),%rcx
	movzwq mpqs3_td_begin(%rip),%rcx

dnl smjs	movq $mpqs3_FB_inv_info,aux3
dnl smjs	movq $mpqs3_FB_start,aux2
	leaq mpqs3_FB_inv_info(%rip),aux3
	leaq mpqs3_FB_start(%rip),aux2

dnl smjs	movw mpqs3_nFBk_1(%rip),aux5w
	movw mpqs3_nFBk_1(%rip),aux5w

	addw $4,aux5w
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
	subw $3,aux5w       # aux5w=mpqs3_nFBk_1+1
	shrl $2,%eax
	orl %eax,aux4d

	shrl $5,aux4d
	movw aux5w,14(relptrw,nr,2)
	adcq $0,nr
	incw aux5w          # aux5w=mpqs3_nFBk_1+2
	shrl $4,aux4d
	movw aux5w,14(relptrw,nr,2)
	adcq $0,nr
	incw aux5w          # aux5w=mpqs3_nFBk_1+3
	shrl $4,aux4d
	movw aux5w,14(relptrw,nr,2)
	adcq $0,nr
	incw aux5w          # aux5w=mpqs3_nFBk_1+4

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
	movw aux5w,14(relptrw,nr,2)
	adcq $0,nr
	incw aux5w
	shrl $4,aux4d
	movw aux5w,14(relptrw,nr,2)
	adcq $0,nr
	incw aux5w
	shrl $4,aux4d
	movw aux5w,14(relptrw,nr,2)
	adcq $0,nr
	incw aux5w
	shrl $4,aux4d
	movw aux5w,14(relptrw,nr,2)
	adcq $0,nr
	incw aux5w

	jmp loop2

dnl Calculate product of primes found by sieving in (prod0,prod1)
dnl Calculate as long as possible in %rax
prod:
dnl smjs	movq $mpqs3_FB,aux1
	leaq mpqs3_FB(%rip),aux1

dnl smjs	movzwq mpqs3_nFBk_1(%rip),%rdx
	movzwq mpqs3_nFBk_1(%rip),%rdx

	addq %rdx,%rdx
	addq %rdx,%rdx
	subq %rdx,aux1

	movq $0,%rcx
	movq $1,%rax
	movq $0,%rdx       # important if nr=0
prodloop1:
	cmpq nr,%rcx
	jnc prodend
	movzwq 14(relptr,%rcx,2),aux5
	incq %rcx
	movzwq (aux1,aux5,4),aux5
	mulq aux5
	jnc prodloop1

prodloop2:
	cmpq nr,%rcx
	jnc prodend
	movzwq 14(relptr,%rcx,2),aux5
	movzwq (aux1,aux5,4),aux5
	movq %rax,prod0
	movq %rdx,%rax
	mulq aux5
	movq %rax,prod1
	movq prod0,%rax
	mulq aux5
	addq prod1,%rdx
	incq %rcx
	jmp prodloop2
	
prodend:
	movq %rax,prod0
	movq %rdx,prod1

dnl End of measuring MMX Loop in zeit20
dnl	movq %rdi,aux5
dnl	movq $20,%rdi
dnl	call zeitB
dnl	movq aux5,%rdi


	testq $1,minus
dnl nr1 and minus are the same register
	movq nr,nr1
	jz posloop
	movw $0,14(relptrw,nr,2)
	incq nr

posloop:
	testq $0x00000001,qx0
	jnz odd
	incq nr
	cmpq NMAXDIV,nr
	jnc gotonext
dnl smjs	movw mpqs3_nFBk_1(%rip),%cx
	movw mpqs3_nFBk_1(%rip),%cx

	movw %cx,12(relptrw,nr,2)
	shrq $1,qx1
	rcrq $1,qx0
	jmp posloop

odd:

division:
	movq prod0,%rax
	movq prod1,aux5
	shrq $32,%rax
	shlq $32,aux5
	movq prod0,aux2
	addq %rax,aux5
	shlq $32,aux2
	subq qx0,aux2
	sbbq qx1,aux5
	jc gotonext

	movq prod0,%rdx
	andq $0xff,%rdx
	shrq $1,%rdx

dnl smjs	movzbl mpqs3_256_inv_table(%rdx),aux2d
dnl smjs THINK ?? I can use aux3 here temporarily and put back after from eax
	leaq mpqs3_256_inv_table(%rip),aux3
	movzbl (aux3, %rdx),aux2d

	movq prod0,aux3
	movl aux3d,%eax

	mull aux2d
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

dnl trial divison of sieved primes
dnl smjs	movq $mpqs3_FB_inv,aux2
	leaq mpqs3_FB_inv(%rip),aux2

dnl smjs	movzwq mpqs3_nFBk_1,%rcx
	movzwq mpqs3_nFBk_1(%rip),%rcx
	addq %rcx,%rcx
	addq %rcx,%rcx
	subq %rcx,aux2
	movq $0,%rcx
tdloop:
	testq nr1,nr1
	movl qxd,%eax
	jz tdend
	movzwq 12(relptr,nr1,2),aux5  # bx: ii
	decq nr1
dnl Value of aux1 is still mpqs3_FB-4*mpqs3_nFBk_1
	movzwl (aux1,aux5,4),%ecx   # cx: p
	movl (aux2,aux5,4),aux3d  # edx: inv
divloop:
	mull aux3d
	movl %eax,aux4d
	mull %ecx
	testl %edx,%edx
	jnz tdloop
	cmpw NMAXDIV,nrw
	jnc gotonext
	movl aux4d,%eax
	movw aux5w,14(relptrw,nr,2)
	incq nr
	movl aux4d,qxd
	jmp divloop


tdend:
dnl trial division of mpqs3_FBk-primes
	cmpl $1,qxd
	jz end

dnl smjs	movq $mpqs3_FBk_inv,aux2
dnl smjs	movq $mpqs3_FBk,aux1
	leaq mpqs3_FBk_inv(%rip),aux2
	leaq mpqs3_FBk(%rip),aux1

dnl smjs	movzwq mpqs3_nFBk(%rip),nr1
	movzwq mpqs3_nFBk(%rip),nr1
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
	cmpw NMAXDIV,nrw
	jnc gotonext
	movl aux4d,%eax
	movw nr1w,14(relptrw,nr,2)
	incq nr
	movl aux4d,qxd
	jmp divloopk

tdendk:
dnl trial division of mpqs3_FB_Adiv-primes
	cmpl $1,%eax
	jz end

dnl smjs	movq $mpqs3_FB_A_inv,aux2
dnl smjs	movq $mpqs3_Adiv_all,aux1
dnl smjs 	movw mpqs3_nFB(%rip),aux5w
dnl smjs	addw mpqs3_nFBk(%rip),aux5w
dnl smjs	movzwq mpqs3_nAdiv_total,nr1
	leaq mpqs3_FB_A_inv(%rip),aux2
	leaq mpqs3_Adiv_all(%rip),aux1
        movw mpqs3_nFB(%rip),aux5w
        addw mpqs3_nFBk(%rip),aux5w
	movzwq mpqs3_nAdiv_total(%rip),nr1

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
	cmpw NMAXDIV,nrw
	jnc gotonext
	addw nr1w,aux5w
	movw aux5w,14(relptrw,nr,2)
	incq nr
	subw nr1w,aux5w
	movl aux4d,qxd
	jmp divloopa

tdenda:
	

end:
	movw nrw,12(relptr)
	movl qxd,%eax
	emms	
	popq %r15
	popq %r14
	popq %r13
	popq %r12
	popq %rbp
	popq %rbx
	ret

gotonext:
	xorq %rax,%rax
	emms
	popq %r15
	popq %r14
	popq %r13
	popq %r12
	popq %rbp
	popq %rbx
	ret

Schlendrian:
	call abort
