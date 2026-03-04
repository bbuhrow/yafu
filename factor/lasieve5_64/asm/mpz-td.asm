dnl
dnl Now, define our arguments. We have 4 saved regs and the return address.
dnl Therefore, our first argument is:
dnl
define(prime,%rdi)dnl
define(modular_inverse,%rsi)dnl
define(limb_ptr_arg,%rdx)dnl
define(nlimbs,%rcx)dnl
define(limb_ptr_ub,%r8)dnl
define(carry,%r9)dnl
define(limb_ptr,%r10)dnl
dnl
function_head(mpz_asm_td)
	movq limb_ptr_arg,limb_ptr
	xorq carry,carry
	leaq (limb_ptr_arg,nlimbs,8),limb_ptr_ub
	movq (limb_ptr_arg),%rax
	cmpq $2,nlimbs
	leaq -8(limb_ptr_ub),limb_ptr_ub
	jb mpz_asm_td_last_limb
mpz_asm_td_loop:
	mulq modular_inverse
	movq %rax,(limb_ptr)
	mulq prime
	leaq 8(limb_ptr),limb_ptr
	movq (limb_ptr),%rax
	addq carry,%rdx
	xorq carry,carry
	subq %rdx,%rax
	adcq carry,carry
	cmpq limb_ptr,limb_ptr_ub
	ja mpz_asm_td_loop
mpz_asm_td_last_limb:
	mulq modular_inverse
	movq %rax,(limb_ptr)
	mulq prime
	addq carry,%rdx
	movq %rdx,%rax
	ret
