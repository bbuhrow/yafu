dnl tdslinie(aux_ptr,aux_ptr_ub,sieve_interval,tds_buffer)
dnl We save %ebx,%esi,%edi,%ebp and also have one auto variable and the
dnl return address on the stack. Therefore, the stack offset of the
dnl first arg is 24.
define(aux_ptr,%rdi)dnl
define(aux_ptr_ub,%rsi)dnl
define(sieve_interval,%rdx)dnl
define(tds_buffer,%rcx)dnl
dnl Now, the registers which we are going to use
define(sieve_ptr,%r8)dnl
define(sieve_ptr_ub,%r9)dnl
define(root,%r10)dnl
define(rootw,%r10w)dnl
define(prime,%r11)dnl
define(prime32,%r11d)dnl
define(sv0,%al)dnl
dnl The ax-register may also be used for auxilliary 32-bit values if sv1
dnl is not used
define(auxreg,%r15)dnl
dnl Offset of the various things from this pointer
define(prime_src,(aux_ptr))dnl
define(proot_src,2(aux_ptr))dnl
define(root_src,6(aux_ptr))dnl
dnl We store the int difference projective_root-prime here:
define(proot,%rbx)dnl
define(sieve_ptr2,%r12)dnl
define(sieve_ptr_ub2,%r13)dnl
define(tds_buffer2,%r14)dnl
dnl This macro is taken from the GNU info documentation of m4.
function_head(tdslinie)
	pushq %rbx
	pushq %r12
	pushq %r13
	pushq %r14
	pushq %r15
	cmpq aux_ptr,aux_ptr_ub
	jbe tdslinie_ende
	subq $8,aux_ptr_ub
tdslinie_fbi_loop:
	movzwq proot_src,auxreg
	movzwq prime_src,prime
	movzwq root_src,root
	subq prime,auxreg
	movq auxreg,proot	
	movq sieve_interval,sieve_ptr_ub
forloop(`i',1,j_per_strip,`
	movq root,sieve_ptr
	xorq auxreg,auxreg
	addq proot,root
	leaq (sieve_ptr_ub,sieve_ptr),sieve_ptr
	cmovncq prime,auxreg
	addq $n_i,sieve_ptr_ub
	addq auxreg,root
	leaq (prime,prime,4),auxreg
	subq auxreg,sieve_ptr_ub
tdslinie_loop`'i:
	movb (sieve_ptr),sv0
	orb (sieve_ptr,prime),sv0
	leaq (sieve_ptr,prime,2),sieve_ptr
	orb (sieve_ptr),sv0
	orb (sieve_ptr,prime),sv0
	testb sv0,sv0
	jz tdslinie_looptest`'i
	leaq (sieve_ptr,prime,2),sieve_ptr_ub2
	negq prime
	leaq (sieve_ptr,prime,2),sieve_ptr2
	negq prime
tdslinie_tloop`'i`'a:
	movzbq (sieve_ptr2),auxreg
	testq auxreg,auxreg
	leaq (sieve_ptr2,prime),sieve_ptr2
	jz tdslinie_tloop`'i`'b
	decq auxreg
	leaq (tds_buffer,auxreg,8),auxreg
	movq (auxreg),tds_buffer2
	movl prime32,(tds_buffer2)
	leaq 4(tds_buffer2),tds_buffer2
	movq tds_buffer2,(auxreg)
tdslinie_tloop`'i`'b:
	cmpq sieve_ptr2,sieve_ptr_ub2
	ja tdslinie_tloop`'i`'a
tdslinie_looptest`'i:
	cmpq sieve_ptr,sieve_ptr_ub
	leaq (sieve_ptr,prime,2),sieve_ptr
	ja tdslinie_loop`'i
	leaq (sieve_ptr_ub,prime,4),sieve_ptr_ub
	xorb sv0,sv0
	cmpq sieve_ptr,sieve_ptr_ub
	movq sieve_ptr,sieve_ptr2
	jbe tdslinie_lasttesta`'i
	movb (sieve_ptr),sv0
	orb (sieve_ptr,prime),sv0
	leaq (sieve_ptr,prime,2),sieve_ptr
tdslinie_lasttesta`'i:
	leaq (sieve_ptr_ub,prime),sieve_ptr_ub
	cmpq sieve_ptr,sieve_ptr_ub
	jbe tdslinie_lasttestb`'i
	orb (sieve_ptr),sv0
tdslinie_lasttestb`'i:
	testb sv0,sv0
	jz tdslinie_next_j`'i
tdslinie_tloop`'i`'c:
	movzbq (sieve_ptr2),auxreg
	testq auxreg,auxreg
	leaq (sieve_ptr2,prime),sieve_ptr2
	jz tdslinie_tloop`'i`'d
	decq auxreg
	leaq (tds_buffer,auxreg,8),auxreg
	movq (auxreg),tds_buffer2
	movl prime32,(tds_buffer2)
	leaq 4(tds_buffer2),tds_buffer2
	movq tds_buffer2,(auxreg)
tdslinie_tloop`'i`'d:
	cmpq sieve_ptr2,sieve_ptr_ub
	ja tdslinie_tloop`'i`'c
tdslinie_next_j`'i:
')
	cmpq aux_ptr,aux_ptr_ub
	movw rootw,root_src
	leaq 8(aux_ptr),aux_ptr
	ja tdslinie_fbi_loop
tdslinie_ende:
	popq %r15
	popq %r14
	popq %r13
	popq %r12
	popq %rbx
	ret
