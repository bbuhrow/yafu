dnl Copyright (C) 2004 Jens Franke, T.Kleinjung
dnl This file is part of gnfs4linux, distributed under the terms of the
dnl GNU General Public Licence and WITHOUT ANY WARRANTY.

dnl You should have received a copy of the GNU General Public License along
dnl with this program; see the file COPYING.  If not, write to the Free
dnl Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.


define(FB_src,%rdi)dnl
define(proots_src,%rsi)dnl
define(fbsz_src,%rdx)dnl
define(a0_src,%rcx)dnl
define(a1_src,%r8)dnl
define(b0_src,%r9)dnl
define(b1_src,56(%rsp))dnl
define(ri_ptr_src,64(%rsp))dnl
define(FB,%r12)dnl
define(ri_ptr,%r13)dnl
define(FB_ub,%r14)dnl
define(p_reg,%rbx)dnl
define(p_regd,%ebx)dnl
define(q_reg,%r15)dnl
define(q_regd,%r15d)dnl
define(proots,%rbp)dnl
dnl Hilfsvariablen
define(aux0,%rax)dnl
define(aux1,%rdx)dnl

define(pq,%xmm0)dnl
dnl The pair of factor base primes
define(rs,%xmm1)dnl
dnl The pair of corresponding projective roots
define(invpq,%xmm2)dnl
dnl Inverse of the pair of factor base primes
define(a0,%xmm3)dnl
dnl hat a0 im ersten und dritten dword zu stehen
define(a1,%xmm4)dnl
dnl Same for a1, etc
define(b0,%xmm5)dnl
define(b1,%xmm6)dnl
define(x,%xmm7)dnl
define(x1,%xmm8)dnl
define(y,%xmm9)dnl
define(y1,%xmm10)dnl
define(xmm_aux0,%xmm11)dnl
define(xmm_aux1,%xmm12)dnl
define(xmm_aux2,%xmm13)dnl
define(xmm_aux3,%xmm14)dnl

dnl Use this comment for lines contributing to the calculation of x
define(Xcal,)dnl
define(X1cal,)dnl
define(Ycal,)dnl
define(Y1cal,)dnl
define(Ical,)dnl
define(nIcal,)dnl
define(InvX,)dnl
dnl %eax and %edx are used in the division.

dnl In general we have to calculate
dnl  ( +- a_1 + rb_1 ) / (+- a_0 - rb_0 ) modulo p.
dnl The signs of a_0 and a_1 depend on the lattice; for each choice we have
dnl a function asm_lasieve_mm_setup'i' , i=0,1,2,3.
dnl
dnl The calculation is done as follows:
dnl Let p be fixed and R(a)=a/2^32 mod p.
dnl We first compute num=R(+- a_1 + rb_1) and den=R(R(+- a_0 - rb_0)).
dnl For den!=0 the result is R(num*den^-1).
dnl If two successive prime ideals of the factor base lie over the same prime p
dnl we try to save one inversion using a trick of Montgomery (rarely one of the
dnl denominators is zero; in this case we do the two calculations seperately).
dnl R(a) is calculated as in Montgomery multiplication.


dnl case a0>=0, a1>=0:	i=0
dnl case a0>=0, a1<0:	i=1
dnl case a0<0, a1>=0:	i=2
dnl case a0<0, a1<0:	i=3
dnl asm_lasieve_mm_setup`'i`'(FB,proots,fbsz,absa0,b0_ul,absa1,b1_ul,ri_ptr)

forloop(i,0,3,`
function_head(asm_lasieve_mm_setup2`'i)
	pushq %rbx
	pushq %r12
	pushq %r13
	pushq %r14
	pushq %r15
	pushq %rbp

	pxor a0,a0
	shrq $1,fbsz_src
	movq FB_src,FB
	pxor b0,b0
	pxor a1,a1
	adcq $0,fbsz_src
	pxor b1,b1
	movd a0_src,a0
	movq proots_src,proots
	movd b0_src,b0
	movd a1_src,a1
	leaq (FB,fbsz_src,8),FB_ub
	movd b1_src,b1
	pshufd $0x88,a0,a0
	pshufd $0x88,b0,b0
	movq ri_ptr_src,ri_ptr
	pshufd $0x88,a1,a1
	pshufd $0x88,b1,b1

	movq (FB),aux0
	movq aux0,aux1
	movd aux0,pq
	movq (proots),rs
	shrq $32,aux1
	andq $0x7fffffff,aux0
	shrq $1,aux0
	shrq $1,aux1
	andq $0x7f,aux0
	andq $0x7f,aux1
	movzbq mpqs_256_inv_table(aux0),aux0
	movzbq mpqs_256_inv_table(aux1),aux1
	pxor invpq,invpq
	pinsrw $0,aux0,invpq
	pinsrw $4,aux1,invpq
	pshufd $0x98,pq,pq
Ical	movdqa invpq,xmm_aux0
Ical	pmuludq pq,invpq
Ical	pmuludq xmm_aux0,invpq
Ical	pslld $1,xmm_aux0
Ical	psubd invpq,xmm_aux0
Ical	movdqa xmm_aux0,invpq
Ical	pmuludq pq,xmm_aux0
Ical	pmuludq invpq,xmm_aux0
	pshufd $0x98,rs,xmm_aux1
Ical	pslld $1,invpq



dnl Begin of long long long loop
loop2`'i:



	movdqa b0,x1
	movl (FB),p_regd
	movl 4(FB),q_regd

X1cal	movdqa xmm_aux1,rs
	leaq 8(FB),FB
X1cal	pcmpeqd pq,xmm_aux1
ifelse(eval((i-1)*(i-2)),0,`
X1cal	pcmpeqd xmm_aux2,xmm_aux2')

	leaq 8(proots),proots
dnl Following is also for Y1
X1cal	pshufd $0xa0,xmm_aux1,xmm_aux1

Ical	psubd xmm_aux0,invpq

ifelse(eval((i-1)*(i-2)),0,`
X1cal	paddd pq,xmm_aux2')

dnl p-r or r, for x
ifelse(eval(i*(i-1)),0,`
Xcal	movdqa pq,x',`
Xcal	movdqa rs,x')

ifelse(eval((i-1)*(i-2)),0,`
X1cal	pmuludq xmm_aux2,x1')

InvX	pxor xmm_aux2,xmm_aux2
X1cal	pand xmm_aux1,x1

ifelse(eval(i*(i-1)),0,`
Xcal	psubd rs,x')

Xcal	pmuludq b0,x
X1cal	movdqa xmm_aux1,xmm_aux0
ifelse(eval(i*(i-2)),0,`
Ycal	movdqa pq,y',`
Ycal	movdqa rs,y')
Xcal	paddq a0,x
Y1cal	movdqa b1,y1
Xcal	pandn x,xmm_aux0

ifelse(eval(i*(i-2)),0,`
Ycal	psubd  rs,y')
dnl rs is no longer needed. Load its value for next round
	movq (proots),rs

X1cal	por   x1,xmm_aux0

dnl case 0,1: x[0,1]=(p-r)*b0+a0
dnl case 2,3: x[0,1]=r*b0+a0
dnl At this point, the case r==p needs no further
dnl seperate treatment. x1 may be used as aux variable

Xcal	movdqa xmm_aux0,x1
Xcal	pmuludq invpq,xmm_aux0
Xcal	pxor x,x
Y1cal	pand xmm_aux1,y1
Xcal	pmuludq pq,xmm_aux0
Xcal	psubq xmm_aux0,x1
Ycal	pmuludq b1,y
Xcal	psrldq $4,x1
Xcal	pcmpgtd x1,x
Xcal	pand pq,x
Xcal	paddd x,x1

Xcal	movdqa x1,x
Xcal	pmuludq invpq,x1
dnl Following line is for later use: figure out if any of the two components
dnl of x is zero. If this line is postponed until after the final calculation
dnl of x, xmm_aux0 mus be loaded with pq instead of 0.
Ycal	paddq a1,y
InvX	pcmpeqd x,xmm_aux2
Xcal	pmuludq pq,x1
Ycal	pandn y,xmm_aux1
Xcal	psubq x1,x
Xcal	psrldq $4,x
Ycal	por y1,xmm_aux1
Xcal	paddd pq,x

dnl At this point, xmm_aux1 holds intermediate
dnl value for y:
dnl case 0,2 xmm_aux1[0,1]=(p-r)*b1+a1
dnl case 1,3 xmm_aux1[0,1]=r*b1+a1

dnl call to asm_modinv
dnl this version returns 0 if its first arg is 0 or p.
dnl Second arg is p.

cmi32B`'i:

dnl If a component of x vanishes mod p, replace it by one. The
dnl information about these components is in xmm_aux2 (and will
dnl be needed one more time).
InvX	psubd  xmm_aux2,x
InvX	movdqa x,x1
dnl Now x1 is the same as x, but each component of x which is zero
dnl has been replaced by 1.
InvX	pshufd $0xfe,x,xmm_aux0
InvX	pmuludq x1,xmm_aux0
Ycal	movdqa xmm_aux1,y
InvX	movdqa xmm_aux0,x
InvX	pmuludq invpq,xmm_aux0
InvX	pxor xmm_aux3,xmm_aux3
InvX	pmuludq pq,xmm_aux0
InvX	psubq xmm_aux0,x
InvX	psrldq $4,x
InvX	pcmpgtd x,xmm_aux3
InvX	pand pq,xmm_aux3
InvX	paddq xmm_aux3,x
	
InvX	movd x,%edi
InvX	pshufd $0x4e,x1,x1
InvX	movl p_regd,%esi

Ycal	pmuludq invpq,xmm_aux1
InvX	call asm_modinv32b

ifelse(eval(i*(i-3)),0,`
InvX	subl %eax,q_regd')
ifelse(eval(i*(i-3)),0,`
InvX	movd q_regd,x',`
InvX	movd %eax,x')

Ycal	pmuludq pq,xmm_aux1
dnl Prepare calculation of inverse of next prime number pair.
	movq (FB),aux0
InvX	pshufd $0x44,x,x
Ycal	pxor y1,y1
InvX	pmuludq x,x1
Ycal	psubq xmm_aux1,y
InvX	movdqa x1,x
InvX	pmuludq invpq,x1
Ycal	psrldq $4,y
	movq aux0,aux1
InvX	pxor xmm_aux3,xmm_aux3
InvX	pmuludq pq,x1
Ycal	pcmpgtd y,y1
InvX	psubq x1,x
Ycal	pand pq,y1
InvX	psrldq $4,x
InvX	pcmpgtd x,xmm_aux3
Ycal	paddd y1,y
InvX	pand pq,xmm_aux3
InvX	paddq xmm_aux3,x

dnl content of xmm_aux0[0]:
dnl case 0 -(2^64*((p-r)*b0+a0))^-1 mod p
dnl case 1  (2^64*((p-r)*b0+a0))^-1 mod p
dnl case 2  (2^64*(r*b0+a0))^-1 mod p
dnl case 3 -(2^64*(r*b0+a0))^-1 mod p
dnl moreover, if the final result is the infinity point of
dnl the projective line, then the corresponding entry in xmm_aux2
dnl equals 0xffffffff

dnl x1 is no longer needed. Load it with next pq pair
dnl to start calculation of inverse
	movd aux0,x1

dnl y1 is no longer needed. Load it with 0.
	pxor y1,y1

	shrq $33,aux1
	shrq $1,aux0
dnl case 0,2 y[0]=2^32*((p-r)*b1+a1) mod p
dnl case 1,3 y[0]=2^32*(r*b1+a1) mod p

	pmuludq y,x
	andq $0x7f,aux0
	andq $0x7f,aux1
	movdqa x,xmm_aux0
	pshufd $0x98,x1,x1
dnl	Also, prepare call to get_recurrence_info
	movq ri_ptr,%rdi
	pmuludq invpq,xmm_aux0
	pxor invpq,invpq
	movzbq mpqs_256_inv_table(aux0),aux0
	movzbq mpqs_256_inv_table(aux1),aux1
	pmuludq pq,xmm_aux0
	pinsrw $0,aux0,invpq
	pinsrw $4,aux1,invpq
	psubq xmm_aux0,x
dnl	prepare call to get_recurrence_info
	movl p_regd,%esi
nIcal	movdqa invpq,xmm_aux0
	psrldq $4,x
nIcal	pmuludq x1,invpq
dnl recall that x1 is the next pq and that y1 is zero.
	pcmpgtd x,y1
	pand pq,y1
	paddd y1,x

dnl In addition, make sure that we do the right thing if we have the
dnl infinity element of P^1(F_p).

nIcal	pmuludq xmm_aux0,invpq
nIcal	pslld $1,xmm_aux0
	pand xmm_aux2,pq
	pandn x,xmm_aux2
nIcal	psubd invpq,xmm_aux0
	por xmm_aux2,pq
	movdqa pq,x

gri2`'i:
	movd x,%edx
	psrldq $8,x
nIcal	movdqa xmm_aux0,invpq
nIcal	pmuludq x1,xmm_aux0
	call get_recurrence_info

	movl p_regd,%esi
	leaq 8(ri_ptr),%rdi
	movd x,%edx
Ical	pmuludq invpq,xmm_aux0
	movdqa x1,pq
	call get_recurrence_info
	cmpq FB,FB_ub
nIcal	pslld $1,invpq
	pshufd $0x98,rs,xmm_aux1
	leaq 16(ri_ptr),ri_ptr
	ja loop2`'i
	popq %rbp
	popq %r15
	popq %r14
	popq %r13
	popq %r12
	popq %rbx
	ret
')
