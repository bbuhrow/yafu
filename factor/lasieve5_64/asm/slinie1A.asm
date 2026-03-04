dnl slinie1(aux_ptr,aux_ptr_ub,sieve_interval)
define(aux_ptr,%rdi)dnl
define(aux_ptr_ub,%rsi)dnl
define(sieve_interval,%rdx)dnl
dnl Now, the registers which we are going to use
define(sieve_ptr,%rcx)dnl
define(sieve_ptr_ub,%rax)dnl
define(root,%r8)dnl
define(prime,%r9)dnl
define(sieve_log,%r10b)dnl
define(sv0,%r11b)dnl
define(sv1,%bl)dnl
dnl The bx-register may also be used for auxilliary 32-bit values if sv1
dnl is not used
define(auxreg,%rbx)dnl
dnl Offset of the various things from this pointer
define(prime_src,(aux_ptr))dnl
define(proot_src,2(aux_ptr))dnl
define(log_src,4(aux_ptr))dnl
define(root_src,6(aux_ptr))dnl
dnl We store the int difference projective_root-prime here:
define(proot,%r12)dnl
ifelse(j_per_strip,1,`divert(-1)define(j_per_strip_minus1,1)',
       `define(j_per_strip_minus1,eval(j_per_strip-1))')dnl
function_head(slinie1)
	cmpq aux_ptr,aux_ptr_ub
	pushq %rbx
	pushq %r12
	jbe slinie1_ende
	subq $8,aux_ptr_ub
slinie1_fbi_loop:
	movzwq proot_src,auxreg
	movzwq prime_src,prime
	movzwq root_src,root
	subq prime,auxreg
	movb log_src,sieve_log
	movq auxreg,proot	
	movq sieve_interval,sieve_ptr_ub
forloop(`i',1,j_per_strip_minus1,`
	movq root,sieve_ptr
	xorq auxreg,auxreg
	addq proot,root
	leaq (sieve_ptr_ub,sieve_ptr),sieve_ptr
	cmovncq prime,auxreg
	addq $n_i,sieve_ptr_ub
	addq auxreg,root
	subq prime,sieve_ptr_ub
	movb (sieve_ptr),sv0
	addb sieve_log,sv0
	cmpq sieve_ptr,sieve_ptr_ub
	movb sv0,(sieve_ptr)
	leaq (sieve_ptr_ub,prime),sieve_ptr_ub
	jbe slinie1_next_j`'i
	addb sieve_log,(sieve_ptr,prime)
slinie1_next_j`'i:
')
	movq root,sieve_ptr
	leaq (sieve_ptr_ub,sieve_ptr),sieve_ptr
	addq $n_i,sieve_ptr_ub
	subq prime,sieve_ptr_ub
	movb (sieve_ptr),sv0
	addb sieve_log,sv0
	cmpq sieve_ptr,sieve_ptr_ub
	movb sv0,(sieve_ptr)
	leaq (sieve_ptr_ub,prime),sieve_ptr_ub
	jbe slinie1_next_fbi
	addb sieve_log,(sieve_ptr,prime)
slinie1_next_fbi:
	cmpq aux_ptr,aux_ptr_ub
#	movq root,root_src
	leaq 8(aux_ptr),aux_ptr
	ja slinie1_fbi_loop
slinie1_ende:
	popq %r12
	popq %rbx
	ret
ifelse(j_per_strip,1,`divert`'',`divert(-1)')dnl
function_head(slinie1)
	cmpq aux_ptr,aux_ptr_ub
	pushq %rbx
	pushq %r12
	jbe slinie1_ende
	subq $8,aux_ptr_ub
slinie1_fbi_loop:
	movzwq proot_src,auxreg
	movzwq prime_src,prime
	movzwq root_src,root
	subq prime,auxreg
	movb log_src,sieve_log
	movq auxreg,proot	
	movq sieve_interval,sieve_ptr_ub
	movq root,sieve_ptr
	leaq (sieve_ptr_ub,sieve_ptr),sieve_ptr
	addq $n_i,sieve_ptr_ub
	subq prime,sieve_ptr_ub
	movb (sieve_ptr),sv0
	addb sieve_log,sv0
	cmpq sieve_ptr,sieve_ptr_ub
	movb sv0,(sieve_ptr)
	leaq (sieve_ptr_ub,prime),sieve_ptr_ub
	jbe slinie1_next_fbi
	addb sieve_log,(sieve_ptr,prime)
slinie1_next_fbi:
	cmpq aux_ptr,aux_ptr_ub
#	movq root,root_src
	leaq 8(aux_ptr),aux_ptr
	ja slinie1_fbi_loop
slinie1_ende:
	popq %r12
	popq %rbx
	ret
