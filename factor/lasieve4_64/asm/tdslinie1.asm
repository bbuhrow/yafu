dnl tdslinie1(aux_ptr,aux_ptr_ub,sieve_interval,tds_buffer)
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
define(auxreg,%r12)dnl
define(auxb,%r12b)dnl
dnl Offset of the various things from this pointer
define(prime_src,(aux_ptr))dnl
define(proot_src,2(aux_ptr))dnl
define(root_src,6(aux_ptr))dnl
dnl We store the int difference projective_root-prime here:
define(proot,%rbx)dnl
define(tds_buffer2,%rax)dnl
dnl This macro is taken from the GNU info documentation of m4.
function_head(tdslinie1)
	cmpq aux_ptr,aux_ptr_ub
	pushq %rbx
	pushq %r12
	jbe tdslinie1_ende
	subq $8,aux_ptr_ub
	leaq -8(tds_buffer),tds_buffer
tdslinie1_fbi_loop:
	movzwq proot_src,auxreg
	movzwq prime_src,prime
	movzwq root_src,root
	subq prime,auxreg
	movq auxreg,proot	
	movq sieve_interval,sieve_ptr_ub
	xorb sv0,sv0
forloop(`i',1,j_per_strip,`
	movq root,sieve_ptr
	xorq auxreg,auxreg
	addq proot,root
	leaq (sieve_ptr_ub,sieve_ptr),sieve_ptr
	cmovncq prime,auxreg
	addq $n_i,sieve_ptr_ub
	orb (sieve_ptr),sv0
	addq auxreg,root
	subq prime,sieve_ptr_ub
	xorq auxreg,auxreg
	cmpq sieve_ptr,sieve_ptr_ub
	leaq (sieve_ptr_ub,prime),sieve_ptr_ub
	cmovaq (sieve_ptr,prime),auxreg
	orb auxb,sv0
')
	testb sv0,sv0
	jnz tdslinie1_suche
tdslinie1_nextfbi:
	cmpq aux_ptr,aux_ptr_ub
	movw rootw,root_src
	leaq 8(aux_ptr),aux_ptr
	ja tdslinie1_fbi_loop
tdslinie1_ende:
	popq %r12
	popq %rbx
	ret

	.align 8
tdslinie1_suche:
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
	subq prime,sieve_ptr_ub
	movb (sieve_ptr),sv0
	cmpq sieve_ptr,sieve_ptr_ub
	leaq (sieve_ptr_ub,prime),sieve_ptr_ub
	jbe tdsline1s_t2_`'i
	orb (sieve_ptr,prime),sv0
tdsline1s_t2_`'i:
	testb sv0,sv0
	jz tdsline1s_next_j`'i
	movzbq (sieve_ptr),auxreg
	testq auxreg,auxreg
	leaq (sieve_ptr,prime),sieve_ptr
	jz tdsline1s_s1_`'i
	movq (tds_buffer,auxreg,8),tds_buffer2
	movl prime32,(tds_buffer2)
	leaq 4(tds_buffer2),tds_buffer2
	movq tds_buffer2,(tds_buffer,auxreg,8)
tdsline1s_s1_`'i:
	cmpq sieve_ptr,sieve_ptr_ub
	jbe tdsline1s_next_j`'i
        movzbq (sieve_ptr),auxreg
        testq auxreg,auxreg
        jz tdsline1s_next_j`'i
        movq (tds_buffer,auxreg,8),sieve_ptr
        movl prime32,(sieve_ptr)
        leaq 4(sieve_ptr),sieve_ptr
        movq sieve_ptr,(tds_buffer,auxreg,8)
tdsline1s_next_j`'i:
')
	jmp tdslinie1_nextfbi
