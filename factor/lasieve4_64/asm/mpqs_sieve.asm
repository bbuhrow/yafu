define(FB,%rcx)dnl current entry of mpqs_FB
define(FBs,%r8)dnl current entry of mps_FB_start
define(FBl,%r9)dnl current entry of mpqs_FB_log
define(sl4,%r10)dnl mpqs_sievelen/4
define(sl,%r11)dnl mpqs_sievelen, used only at the beginning
define(sb,%r11)dnl
define(sloc1,%rsi)dnl
define(sloc2,%rdi)dnl
define(sievebound,%rbx)dnl
define(aux1,%r12)dnl
define(sieve_interval,%r13)dnl
define(prime,%r14)dnl
define(prime2,%rdx)dnl
define(sieveend,%r15)dnl
function_head(asm_sieve)
	pushq %rbx
	pushq %r12
	pushq %r13
	pushq %r14
	pushq %r15
	movzwq mpqs_sievebegin(%rip),%rax
	movzwq mpqs_sievelen(%rip),sl4
	movq $mpqs_FB,FB
	movq $mpqs_FB_start,FBs
	movq $mpqs_FB_log,FBl
	leaq (FB,%rax,4),FB
	leaq (FBs,%rax,4),FBs
	leaq (FBl,%rax),FBl
	movq sl4,sl
	shrq $2,sl4
	movq mpqs_sievearray(%rip),sieve_interval
	movq sieve_interval,sieveend
	addq sl,sieveend

.align 16
mainloop:
	movzwq (FB),prime   # p
	cmpq prime,sl4
	leaq 4(FB),FB
	jc loop3begin

	movzwq (FBs),sloc1
	movzwq 2(FBs),sloc2
	leaq 4(FBs),FBs

	cmpq sloc1,sloc2   # sort sloc1,sloc2:
	movq sloc1,aux1
	movq sloc2,prime2
	movb (FBl),%al
	leaq 1(FBl),FBl
	cmovcq aux1,sloc2
	cmovcq prime2,sloc1

	addq sieve_interval,sloc1
	addq sieve_interval,sloc2

	movq sieveend,sievebound
	leaq (prime,prime,2),aux1
	leaq (prime,prime),prime2
	subq aux1,sievebound

#	movq sievebound,sb
#	subq prime2,sb
#	subq prime2,sb

loop4:
	addb %al,(sloc2)
	addb %al,(sloc2,prime)
	addb %al,(sloc1)
	addb %al,(sloc1,prime)
	addb %al,(sloc2,prime2)
	addb %al,(sloc2,aux1)
	addb %al,(sloc1,prime2)
	addb %al,(sloc1,aux1)
#	cmpq sb,sloc2
	leaq (sloc2,prime2,2),sloc2
	leaq (sloc1,prime2,2),sloc1
	cmpq sievebound,sloc2
	jc loop4
dnl for the first test a jump seems to be faster than 4 cmovs
	addq prime2,sievebound
	cmpq sievebound,sloc2
	jnc check

	addb %al,(sloc1)
	addb %al,(sloc1,prime)
	addb %al,(sloc2)
	addb %al,(sloc2,prime)
	addq prime2,sloc1
	addq prime2,sloc2

check:
dnl for the remaining 3 tests cmovs seem to be faster
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

	jmp mainloop

loop3begin:
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
	leaq 4(FBs),FBs

	addq sieve_interval,sloc1
	addq sieve_interval,sloc2
	movb (FBl),%al
	leaq 1(FBl),FBl

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

	jmp mainloop3

loop2begin:
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
	leaq 4(FBs),FBs

	addq sieve_interval,sloc1
	addq sieve_interval,sloc2
	movb (FBl),%al
	leaq 1(FBl),FBl

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

	jmp mainloop2

loop1begin:
	movq sl,sl4
	subq $1,sl4
	leaq -4(FB),FB

mainloop1:
	movzwq (FB),prime   # p
	cmpq prime,sl4
	leaq 4(FB),FB
	jc loop0begin

	movzwq (FBs),sloc1
	movzwq 2(FBs),sloc2
	leaq 4(FBs),FBs

	addq sieve_interval,sloc1
	addq sieve_interval,sloc2
	movb (FBl),%al
	leaq 1(FBl),FBl

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

	jmp mainloop1

loop0begin:
	movq $0xffff,sl4
	leaq -4(FB),FB

mainloop0:
	movzwq (FB),prime   # p
	cmpq prime,sl4
	leaq 4(FB),FB
	jz end

	movzwq (FBs),sloc1
	movzwq 2(FBs),sloc2
	leaq 4(FBs),FBs

	addq sieve_interval,sloc1
	addq sieve_interval,sloc2
	movb (FBl),%al
	leaq 1(FBl),FBl

        cmpq sieveend,sloc1
        cmovncq sieveend,sloc1
        cmpq sieveend,sloc2
        cmovncq sieveend,sloc2
        addb %al,(sloc1)
        addb %al,(sloc2)

	jmp mainloop0
end:
	popq %r15
	popq %r14
	popq %r13
	popq %r12
	popq %rbx
	ret
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

