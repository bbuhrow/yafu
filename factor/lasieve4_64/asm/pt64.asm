define(N,%rdi)dnl
define(one,%rsi)dnl
define(auxreg,%rcx)dnl
dnl auxreg will be destroyed by squaring.
dnl other auxregs which may be used outside the squaring
dnl are %rax and %rdx
dnl auxreg2 is preserved by squaring.
define(mm64_aux,%r8)dnl
define(expon,%r9)dnl
define(b,%r10)dnl
define(t,%r11)dnl

	.align 8
	.type modsq64,@function
modsq64:
	movq t,%rax
	mulq t
	movq %rdx,t
	mulq mm64_aux
	mulq N
	xorq %rax,%rax
	subq %rdx,t
	cmovbq N,%rax
	addq %rax,t
	ret

function_head(pt64)
dnl Calculate one
	movq $1,%rdx
	xorq %rax,%rax
	divq %rdi
	movq N,auxreg
	movq N,expon
	andq $255,auxreg
	subq $1,expon
	shrq $1,auxreg
	movzbq mpqs_256_inv_table(auxreg),mm64_aux
	movq %rdx,one

dnl	movq %rdx,%rax

pt64_auxcalc:
forloop(`i',1,3,`
	movq mm64_aux,%rax
	mulq mm64_aux
	shlq $1,mm64_aux
	mulq N
	subq %rax,mm64_aux')
	
dnl	movq mm64_aux,%rax

pt64_expon_bcalc:
	bsfq expon,auxreg
	shrq %cl,expon
	movq auxreg,b

dnl	movq expon,%rax
dnl	movq b,%rax

pt64_mcalc:
	bsrq expon,auxreg
	movq $1,%rax
	shlq %cl,%rax
	movq %rax,auxreg

	movq one,t
	jmp pt64_modpow_dupt
pt64_modpow:
	call modsq64
	testq auxreg,expon
	jz pt64_modpow_next
pt64_modpow_dupt:
	xorq %rax,%rax
	shlq $1,t
	cmovcq N,%rax
	cmpq N,t
	cmovaeq N,%rax
	subq %rax,t
pt64_modpow_next:
	shrq $1,auxreg
	jnz pt64_modpow

dnl	movq t,%rax

	cmpq t,one
	movq N,auxreg
	je pt64_prime
	subq one,auxreg
pt64_RMtests:
	cmpq t,auxreg
	je pt64_prime
	decq b
	jz pt64_comp
	call modsq64
	jmp pt64_RMtests

	.align 8
pt64_comp:
	xorq %rax,%rax
	ret

	.align 8
pt64_prime:
	movq $1,%rax
	ret
