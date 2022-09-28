dnl tdslinie3(aux_ptr,aux_ptr_ub,sieve_interval,tds_buffer)
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
dnl The bx-register may also be used for auxilliary 32-bit values if sv1
dnl is not used
define(auxreg,%r12)dnl
define(auxb,%r12b)dnl
dnl Offset of the various things from this pointer
define(prime_src,(aux_ptr))dnl
define(proot_src,2(aux_ptr))dnl
define(root_src,6(aux_ptr))dnl
dnl We store the int difference projective_root-prime here:
define(proot,%rbx)dnl
define(tds_buffer2,%rax)
function_head(tdslinie3)
	cmpq aux_ptr,aux_ptr_ub
	pushq %rbx
	pushq %r12
	jbe tdslinie3_ende
	subq $8,aux_ptr_ub
tdslinie3_fbi_loop:
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
	movb (sieve_ptr),sv0
	subq prime,sieve_ptr_ub
	orb (sieve_ptr,prime),sv0
	leaq (sieve_ptr,prime,2),sieve_ptr
	orb (sieve_ptr),sv0
	cmpq sieve_ptr,sieve_ptr_ub
	jbe tdslinie3_t2_`'i
	orb (sieve_ptr,prime),sv0
tdslinie3_t2_`'i:
	testb sv0,sv0
	leaq (sieve_ptr_ub,prime),sieve_ptr_ub
	jz tdslinie3_next_j`'i
	subq prime,sieve_ptr
	subq prime,sieve_ptr
tdslinie3_tloop`'i:
	movzbq (sieve_ptr),auxreg
	testq auxreg,auxreg
	leaq (sieve_ptr,prime),sieve_ptr
	jz tdslinie3_tdloop`'i`'a
	decq auxreg
	leaq (tds_buffer,auxreg,8),auxreg
	movq (auxreg),tds_buffer2
	movl prime32,(tds_buffer2)
	leaq 4(tds_buffer2),tds_buffer2
	movq tds_buffer2,(auxreg)
tdslinie3_tdloop`'i`'a:
	cmpq sieve_ptr,sieve_ptr_ub
	ja tdslinie3_tloop`'i
tdslinie3_next_j`'i:
')
	cmpq aux_ptr,aux_ptr_ub
	movw rootw,root_src
	leaq 8(aux_ptr),aux_ptr
	ja tdslinie3_fbi_loop
tdslinie3_ende:
	popq %r12
	popq %rbx
	ret
