dnl Copyright (C) 2002 Jens Franke, T.Kleinjung
dnl This file is part of gnfs4linux, distributed under the terms of the
dnl GNU General Public Licence and WITHOUT ANY WARRANTY.

dnl You should have received a copy of the GNU General Public License along
dnl with this program; see the file COPYING.  If not, write to the Free
dnl Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.

#include "underscore.h" 

define(FB,%rcx)dnl current entry of mpqs_FB
define(FBs,%r8)dnl current entry of mps_FB_start
define(FBl,%r9)dnl current entry of mpqs_FB_log
define(sl4,%r10)dnl mpqs_sievelen/4
define(sl,%r11)dnl mpqs_sievelen
define(sloc1,%rsi)dnl
define(sloc2,%rdi)dnl
define(sievebound,%rbx)dnl
define(aux1,%r12)dnl
define(sieve_interval,%r13)dnl
define(prime,%r14)dnl
define(prime2,%rdx)dnl
define(sieveend,%r15)dnl
function_head(asm3_sieve)
	pushq %rbx
	pushq %r12
	pushq %r13
	pushq %r14
	pushq %r15
dnl smjs	movzwq mpqs3_sievebegin(%rip),%rax
dnl smjs	movzwq mpqs3_sievelen(%rip),sl4
	movzwq mpqs3_sievebegin(%rip),%rax
	movzwq mpqs3_sievelen(%rip),sl4

dnl smjs	movq $mpqs3_FB,FB
dnl smjs	movq $mpqs3_FB_start,FBs
dnl smjs	movq $mpqs3_FB_log,FBl
	leaq mpqs3_FB(%rip),FB
	leaq mpqs3_FB_start(%rip),FBs
	leaq mpqs3_FB_log(%rip),FBl


	leaq (FB,%rax,4),FB
        leaq (FBs,%rax,4),FBs
	leaq (FBl,%rax),FBl
	movq sl4,sl
	shrq $3,sl4
dnl smjs	movq mpqs3_sievearray(%rip),sieve_interval
	movq mpqs3_sievearray(%rip),sieve_interval

	movq sieve_interval,sieveend
	addq sl,sieveend

.if 0
 movq $40,%rdi
 call zeitA
.endif

.align 16
mainloop:
	movzwq (FB),prime   # p
	cmpq prime,sl4
	leaq 4(FB),FB
	jc loop7begin

	movzwq (FBs),sloc1
	movzwq 2(FBs),sloc2

        cmpq sloc1,sloc2   # sort sloc1,sloc2:
        movq sloc1,aux1
        movq sloc2,prime2
        cmovcq aux1,sloc2
        cmovcq prime2,sloc1

	leaq (prime,prime,2),aux1
	leaq (prime,prime),prime2
	movb (FBl),%al
	addq $4,FBs
	addq sieve_interval,sloc1
	addq sieve_interval,sloc2

loop4:
	addb %al,(sloc1)
	addb %al,(sloc1,prime)
	addb %al,(sloc2)
	addb %al,(sloc2,prime)
	addq prime2,sloc1
	addq prime2,sloc2
	addb %al,(sloc1)
	addb %al,(sloc1,prime)
	addb %al,(sloc2)
	addb %al,(sloc2,prime)
	addq prime2,sloc1
	addq prime2,sloc2
        cmpq sieveend,sloc1
	jc loop4
	incq FBl

	jmp mainloop


loop7begin:
.if 0
 movq $40,%rdi
 call zeitB

 movq $41,%rdi
 call zeitA
.endif
	movq sl,%rax
	xorq %rdx,%rdx  # prime2 is not used at the moment
	movq $7,sl4
	divq sl4
	movq %rax,sl4   # sl4=sl/7 in this part
	leaq -4(FB),FB

mainloop7:
	movzwq (FB),prime   # p
	cmpq prime,sl4
	leaq 4(FB),FB
	jc loop6begin

	movzwq (FBs),sloc1
	movzwq 2(FBs),sloc2

	addq sieve_interval,sloc1
	addq sieve_interval,sloc2
	movb (FBl),%al
	leaq (prime,prime),prime2
	addq $4,FBs

	addb %al,(sloc1)
	addb %al,(sloc2)
	addb %al,(sloc1,prime)
	addb %al,(sloc2,prime)
	addq prime2,sloc1
	addq prime2,sloc2
	addb %al,(sloc1)
	addb %al,(sloc2)
	addb %al,(sloc1,prime)
	addb %al,(sloc2,prime)
	addq prime2,sloc1
	addq prime2,sloc2
	addb %al,(sloc1)
	addb %al,(sloc2)
	addb %al,(sloc1,prime)
	addb %al,(sloc2,prime)
	addq prime2,sloc1
	addq prime2,sloc2
	addb %al,(sloc1)
	addb %al,(sloc2)
	addb %al,(sloc1,prime)
	addb %al,(sloc2,prime)
	incq FBl

	jmp mainloop7

loop6begin:
	movq sl,%rax
	xorq %rdx,%rdx  # prime2 is not used at the moment
	movq $6,sl4
	divq sl4
	movq %rax,sl4   # sl4=sl/6 in this part
	leaq -4(FB),FB

mainloop6:
	movzwq (FB),prime   # p
	cmpq prime,sl4
	leaq 4(FB),FB
	jc loop5begin

	movzwq (FBs),sloc1
	movzwq 2(FBs),sloc2

	addq sieve_interval,sloc1
	addq sieve_interval,sloc2
	movb (FBl),%al
	leaq (prime,prime),prime2
	addq $4,FBs

	addb %al,(sloc1)
	addb %al,(sloc2)
	addb %al,(sloc1,prime)
	addb %al,(sloc2,prime)
	addq prime2,sloc1
	addq prime2,sloc2
	addb %al,(sloc1)
	addb %al,(sloc2)
	addb %al,(sloc1,prime)
	addb %al,(sloc2,prime)
	addq prime2,sloc1
	addq prime2,sloc2
	addb %al,(sloc1)
	addb %al,(sloc2)
	addb %al,(sloc1,prime)
	addb %al,(sloc2,prime)
	addb %al,(sloc1,prime2)
	addb %al,(sloc2,prime2)
	incq FBl

	jmp mainloop6

loop5begin:
	movq sl,%rax
	xorq %rdx,%rdx  # prime2 is not used at the moment
	movq $5,sl4
	divq sl4
	movq %rax,sl4   # sl4=sl/5 in this part
	leaq -4(FB),FB

mainloop5:
	movzwq (FB),prime   # p
	cmpq prime,sl4
	leaq 4(FB),FB
	jc loop4begin

	movzwq (FBs),sloc1
	movzwq 2(FBs),sloc2

	addq sieve_interval,sloc1
	addq sieve_interval,sloc2
	movb (FBl),%al
	leaq (prime,prime),prime2
	addq $4,FBs

	addb %al,(sloc1)
	addb %al,(sloc2)
	addb %al,(sloc1,prime)
	addb %al,(sloc2,prime)
	addq prime2,sloc1
	addq prime2,sloc2
	addb %al,(sloc1)
	addb %al,(sloc2)
	addb %al,(sloc1,prime)
	addb %al,(sloc2,prime)
	addq prime2,sloc1
	addq prime2,sloc2
	addb %al,(sloc1)
	addb %al,(sloc2)
	addb %al,(sloc1,prime)
	addb %al,(sloc2,prime)
	incq FBl

	jmp mainloop5

loop4begin:
	movq sl,%rax
	xorq %rdx,%rdx  # prime2 is not used at the moment
	movq $4,sl4
	divq sl4
	movq %rax,sl4   # sl4=sl/4 in this part
	leaq -4(FB),FB

mainloop4:
	movzwq (FB),prime   # p
	cmpq prime,sl4
	leaq 4(FB),FB
	jc loop3begin

	movzwq (FBs),sloc1
	movzwq 2(FBs),sloc2

	addq sieve_interval,sloc1
	addq sieve_interval,sloc2
	movb (FBl),%al
	leaq (prime,prime),prime2
	addq $4,FBs

	addb %al,(sloc1)
	addb %al,(sloc2)
	addb %al,(sloc1,prime)
	addb %al,(sloc2,prime)
	addq prime2,sloc1
	addq prime2,sloc2
	addb %al,(sloc1)
	addb %al,(sloc2)
	addb %al,(sloc1,prime)
	addb %al,(sloc2,prime)
	addb %al,(sloc1,prime2)
	addb %al,(sloc2,prime2)
	incq FBl

	jmp mainloop4



loop3begin:
.if 0
 movq $41,%rdi
 call zeitB

 movq $42,%rdi
 call zeitA
.endif

	movq sl,%rax
	xorq %rdx,%rdx  # prime2 is not used at the moment
	movq $3,sl4
	divq sl4
	movq %rax,sl4   # sl4=sl/3 in this part
	leaq -4(FB),FB

mainloop3:
	movzwq (FB),prime   # p
	cmpq prime,sl4
	leaq 4(FB),FB
	jc loop2begin

	movzwq (FBs),sloc1
	movzwq 2(FBs),sloc2

	addq sieve_interval,sloc1
	addq sieve_interval,sloc2
	movb (FBl),%al
	addq $4,FBs

	addb %al,(sloc1)
	addb %al,(sloc2)
	addb %al,(sloc1,prime)
	addb %al,(sloc2,prime)
	leaq (sloc1,prime,2),sloc1
	leaq (sloc2,prime,2),sloc2
	addb %al,(sloc1)
	addb %al,(sloc2)
	addb %al,(sloc1,prime)
	addb %al,(sloc2,prime)
	incq FBl

	jmp mainloop3

loop2begin:
.if 0
 movq $42,%rdi
 call zeitB

 movq $43,%rdi
 call zeitA
.endif

	movq sl,sl4
	shrq $1,sl4
	leaq -4(FB),FB

mainloop2:
	movzwq (FB),prime   # p
	cmpq prime,sl4
	leaq 4(FB),FB
	jc loop1begin

	movzwq (FBs),sloc1
	movzwq 2(FBs),sloc2

	addq sieve_interval,sloc1
	addq sieve_interval,sloc2
	movb (FBl),%al
	addq $4,FBs

        addb %al,(sloc1,prime)
        addb %al,(sloc2,prime)
	addq prime,prime
        addb %al,(sloc1)
        addb %al,(sloc2)
        addb %al,(sloc1,prime)
        addb %al,(sloc2,prime)
	incq FBl

	jmp mainloop2

loop1begin:
.if 0
 movq $43,%rdi
 call zeitB

 movq $44,%rdi
 call zeitA
.endif

	movq sl,sl4
	subq $1,sl4
	leaq -4(FB),FB

mainloop1:
	movzwq (FB),prime   # p
	cmpq prime,sl4
	leaq 4(FB),FB
	jc end

	movzwq (FBs),sloc1
	movzwq 2(FBs),sloc2

	addq sieve_interval,sloc1
	addq sieve_interval,sloc2
	movb (FBl),%al
	addq $4,FBs

        addb %al,(sloc1)
        addb %al,(sloc2)
        addb %al,(sloc1,prime)
        addb %al,(sloc2,prime)
	incq FBl

	jmp mainloop1
end:
.if 0
 movq $44,%rdi
 call zeitB
.endif
	popq %r15
	popq %r14
	popq %r13
	popq %r12
	popq %rbx
	ret

.if 1
undefine(`FB')dnl
undefine(`FBs')dnl
undefine(`FBl')dnl
undefine(`sl4')dnl
undefine(`sl')dnl
undefine(`sloc1')dnl
undefine(`sloc2')dnl
undefine(`sievebound')dnl
undefine(`aux1')dnl
undefine(`sieve_interval')dnl
undefine(`prime')dnl
undefine(`prime2')dnl
undefine(`sieveend')dnl
define(FB0,%rcx)dnl current entry of mpqs_FB0
define(FBs,%r8)dnl current entry of mps_FB_start
define(sl,%r11)dnl mpqs_sievelen
define(sloc1,%rsi)dnl
define(sloc2,%rdi)dnl
define(sievebound,%rbx)dnl
define(aux1,%r12)dnl
define(sieve_interval,%r13)dnl
define(prime,%r14)dnl
define(prime_byte,%r14b)dnl
define(prime2,%rdx)dnl
define(sieveend,%r15)dnl
function_head(asm3_sievea)
	pushq %rbx
	pushq %r12
	pushq %r13
	pushq %r14
	pushq %r15
dnl smjs	movzwq mpqs3_sievebegin(%rip),%rax
dnl smjs	movzwq mpqs3_sievelen(%rip),sl
	movzwq mpqs3_sievebegin(%rip),%rax
	movzwq mpqs3_sievelen(%rip),sl

dnl smjs	movq $mpqs3_FB0,FB0
dnl smjs	movq $mpqs3_FB_start,FBs
	leaq mpqs3_FB0(%rip),FB0
	leaq mpqs3_FB_start(%rip),FBs

        leaq (FBs,%rax,4),FBs

dnl smjs	movq mpqs3_sievearray(%rip),sieve_interval
	movq mpqs3_sievearray(%rip),sieve_interval

	movq sieve_interval,sieveend
	addq sl,sieveend

.if 0
 movq $21,%rdi
 call zeitA
.endif

.align 16
mainloop8a:
	movzwq (FB0),prime   # p
	testq prime,prime
	leaq 2(FB0),FB0
	jz update8a
return8a:
	movzwq (FBs),sloc1
	movzwq 2(FBs),sloc2
	leaq 4(FBs),FBs

	cmpq sloc1,sloc2   # sort sloc1,sloc2:
	movq sloc1,aux1
	movq sloc2,prime2
	cmovcq aux1,sloc2
	cmovcq prime2,sloc1

	addq sieve_interval,sloc1
	addq sieve_interval,sloc2

	movq sieveend,sievebound
	leaq (prime,prime,2),aux1
	leaq (prime,prime),prime2
	subq aux1,sievebound

loop8a:
        addb %al,(sloc2)
        addb %al,(sloc2,prime)
        addb %al,(sloc1)
        addb %al,(sloc1,prime)
        addb %al,(sloc2,prime2)
        addb %al,(sloc2,aux1)
        addb %al,(sloc1,prime2)
        addb %al,(sloc1,aux1)
        leaq (sloc2,prime2,2),sloc2
        leaq (sloc1,prime2,2),sloc1
        cmpq sievebound,sloc2
	jc loop8a
dnl for the first test a jump seems to be faster than 4 cmov's
	addq prime2,sievebound
	cmpq sievebound,sloc2
	jnc check8a

	addb %al,(sloc1)
	addb %al,(sloc1,prime)
	addb %al,(sloc2)
	addb %al,(sloc2,prime)
	addq prime2,sloc1
	addq prime2,sloc2

check8a:
dnl for the remaining 3 tests cmov' s seem to be faster
        cmpq sieveend,sloc1
        cmovncq sieveend,sloc1
        cmpq sieveend,sloc2
        cmovncq sieveend,sloc2
        addb %al,(sloc1)
        addq prime,sloc1
        addb %al,(sloc2)
        cmpq sieveend,sloc1
        cmovncq sieveend,sloc1
        addb %al,(sloc1)

	jmp mainloop8a

update8a:
	movzwq (FB0),prime
	testq prime,prime
	jz loop4begina
	movb prime_byte,%al
	movzwq 2(FB0),prime     # p
	leaq 4(FB0),FB0
	jmp return8a

loop4begina:
	leaq 2(FB0),FB0

.if 0
 movq %rax,sloc1
 movq $21,%rdi
 call zeitB

 movq $22,%rdi
 call zeitA
 movq sloc1,%rax
.endif

mainloop4a:
	movzwq (FB0),prime   # p
	testq prime,prime
	leaq 2(FB0),FB0
	jz update4a
return4a:
	movzwq (FBs),sloc1
	movzwq 2(FBs),sloc2
	leaq 4(FBs),FBs

	addq sieve_interval,sloc1
	addq sieve_interval,sloc2

	movq sieveend,sievebound
	leaq (prime,prime,2),aux1
	leaq (prime,prime),prime2
	subq prime,sievebound

        addb %al,(sloc2)
        addb %al,(sloc2,prime)
        addb %al,(sloc1)
        addb %al,(sloc1,prime)
        addb %al,(sloc2,prime2)
        addb %al,(sloc2,aux1)
        addb %al,(sloc1,prime2)
        addb %al,(sloc1,aux1)
        leaq (sloc2,prime2,2),sloc2
        leaq (sloc1,prime2,2),sloc1
dnl for the first test a jump seems to be faster than 4 cmov's
	cmpq sievebound,sloc2
	jnc check4a

	addb %al,(sloc1)
	addb %al,(sloc1,prime)
	addb %al,(sloc2)
	addb %al,(sloc2,prime)
	addq prime2,sloc1
	addq prime2,sloc2

check4a:
dnl for the remaining 4 tests cmov' s seem to be faster
        cmpq sieveend,sloc1
        cmovncq sieveend,sloc1
        cmpq sieveend,sloc2
        cmovncq sieveend,sloc2
        addb %al,(sloc1)
        addq prime,sloc1
        addb %al,(sloc2)
        addq prime,sloc2
        cmpq sieveend,sloc1
        cmovncq sieveend,sloc1
        cmpq sieveend,sloc2
        cmovncq sieveend,sloc2
        addb %al,(sloc1)
        addb %al,(sloc2)

	jmp mainloop4a

update4a:
	movzwq (FB0),prime
	testq prime,prime
	jz loop3begina
	movb prime_byte,%al
	movzwq 2(FB0),prime     # p
	leaq 4(FB0),FB0
	jmp return4a


loop3begina:
	leaq 2(FB0),FB0
.if 0
 movq sloc1,%rax
 movq $22,%rdi
 call zeitB

 movq $23,%rdi
 call zeitA
 movq sloc1,%rax
.endif
mainloop3a:
	movzwq (FB0),prime   # p
	testq prime,prime
	leaq 2(FB0),FB0
	jz update3a
return3a:
	movzwq (FBs),sloc1
	movzwq 2(FBs),sloc2
	leaq 4(FBs),FBs

	addq sieve_interval,sloc1
	addq sieve_interval,sloc2

	leaq (prime,prime,2),aux1
	leaq (prime,prime),prime2

        addb %al,(sloc1)
        addb %al,(sloc1,prime)
        addb %al,(sloc2)
        addb %al,(sloc2,prime)
        addb %al,(sloc1,prime2)
        addb %al,(sloc2,prime2)
        addq aux1,sloc1
        addq aux1,sloc2

        cmpq sieveend,sloc1
        cmovncq sieveend,sloc1
        cmpq sieveend,sloc2
        cmovncq sieveend,sloc2
        addb %al,(sloc1)
        addb %al,(sloc2)

	jmp mainloop3a
update3a:
	movzwq (FB0),prime
	testq prime,prime
	jz loop2begina
	movb prime_byte,%al
	movzwq 2(FB0),prime     # p
	leaq 4(FB0),FB0
	jmp return3a


loop2begina:
	leaq 2(FB0),FB0

.if 0
 movq sloc1,%rax
 movq $23,%rdi
 call zeitB

 movq $24,%rdi
 call zeitA
 movq sloc1,%rax
.endif

mainloop2a:
	movzwq (FB0),prime   # p
	testq prime,prime
	leaq 2(FB0),FB0
	jz update2a
return2a:
	movzwq (FBs),sloc1
	movzwq 2(FBs),sloc2
	leaq 4(FBs),FBs

	addq sieve_interval,sloc1
	addq sieve_interval,sloc2

	leaq (prime,prime),prime2

        addb %al,(sloc1)
        addb %al,(sloc1,prime)
        addb %al,(sloc2)
        addb %al,(sloc2,prime)
        addq prime2,sloc1
        addq prime2,sloc2

        cmpq sieveend,sloc1
        cmovncq sieveend,sloc1
        cmpq sieveend,sloc2
        cmovncq sieveend,sloc2
        addb %al,(sloc1)
        addb %al,(sloc2)

	jmp mainloop2a
update2a:
	movzwq (FB0),prime
	testq prime,prime
	jz loop1begina
	movb prime_byte,%al
	movzwq 2(FB0),prime     # p
	leaq 4(FB0),FB0
	jmp return2a


loop1begina:
	leaq 2(FB0),FB0

.if 0
 movq sloc1,%rax
 movq $24,%rdi
 call zeitB

 movq $25,%rdi
 call zeitA
 movq sloc1,%rax
.endif

mainloop1a:
	movzwq (FB0),prime   # p
	testq prime,prime
	leaq 2(FB0),FB0
	jz update1a
return1a:
	movzwq (FBs),sloc1
	movzwq 2(FBs),sloc2
	leaq 4(FBs),FBs

	addq sieve_interval,sloc1
	addq sieve_interval,sloc2

        addb %al,(sloc1)
        addb %al,(sloc2)
        addq prime,sloc1
        addq prime,sloc2

        cmpq sieveend,sloc1
        cmovncq sieveend,sloc1
        cmpq sieveend,sloc2
        cmovncq sieveend,sloc2
        addb %al,(sloc1)
        addb %al,(sloc2)

	jmp mainloop1a
update1a:
	movzwq (FB0),prime
	testq prime,prime
	jz loop0begina
	movb prime_byte,%al
	movzwq 2(FB0),prime     # p
	leaq 4(FB0),FB0
	jmp return1a

loop0begina:
	leaq 2(FB0),FB0

.if 0
 movq sloc1,%rax
 movq $25,%rdi
 call zeitB
 movq sloc1,%rax
.endif

	popq %r15
	popq %r14
	popq %r13
	popq %r12
	popq %rbx
	ret
.endif

