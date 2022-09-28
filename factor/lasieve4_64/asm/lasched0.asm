define(ri,%rdi)dnl
define(ij_ptr,%rsi)dnl
define(ij_ptr_ub_arg,%rdx)dnl
define(ij_ptr_ub,(%rsp))
dnl define(n1_j,%rcx)
define(ij_ub,%rcx)dnl
define(ij_ubd,%ecx)dnl
define(sched,%r8)dnl
define(fbi_offs,%r9)dnl
define(fbi_offsd,%r9d)dnl
define(ij,%rax)dnl
define(ijd,%eax)dnl
define(aux0,%rbx)dnl
define(aux0d,%ebx)dnl
define(aux1,%r10)dnl
define(aux1d,%r10d)dnl
define(aux2,%r11)dnl
define(aux2d,%r11d)dnl
define(aux3,%r12)dnl
define(aux3d,%r12d)dnl
define(aux4,%r13)dnl
define(aux4d,%r13d)dnl
define(a,%r14)dnl
define(ad,%r14d)dnl
define(b,%r15)dnl
define(bd,%r15d)dnl
define(v0d,%edx)
define(v1d,%ebp)
ifelse(nt_sched,1,
	`function_head(lasched0nt)',
	`function_head(lasched0)')
	pushq %rbx
	shlq $n_i_bits,ij_ub
	pushq %r12
	movl (ij_ptr),ijd
	pushq %r13
	cmpq ij_ptr,ij_ptr_ub_arg
	leaq -4(ij_ptr_ub_arg),ij_ptr_ub_arg
	pushq %r14
	jbe lasched0_fbi_loop_end
	shrq $16,fbi_offs
	pushq %r15
	pushq %rbp
	pushq ij_ptr_ub_arg
lasched_fbi_loop:
	movl (ri),ad
	movl 4(ri),bd
	cmpl ijd,ij_ubd
	movl ijd,aux0d
	movl ijd,aux1d
	prefetch 128(ri)
	jbe lasched0_ijloop_end
	movl ad,v0d
	negl ad
	movl bd,v1d
	negl bd
	andl $n_i_mask,ad
	andl $n_i_mask,bd
lasched0_ij_loop:
	movl ijd,aux2d
	andq $n_i_mask,aux0
	shrq $l1_bits,aux1
	andl $l1_mask,aux2d
	cmpl aux0d,ad
	movl $0,aux3d
	movl $0,aux4d
	cmovbel v1d(ri),aux3d
	orl fbi_offsd,aux2d
	cmpq aux0,b
	movq (sched,aux1,8),aux0
	cmoval v0d,aux4d
	addl aux3d,ijd
ifelse(nt_sched,1,`
	movnti aux2d,(aux0)',`
	movl aux2d,(aux0)')
	leaq 4(aux0),aux0
	addl aux4d,ijd
	movq aux0,(sched,aux1,8)
	cmpl ijd,ij_ubd
	movl ijd,aux0d
	movl ijd,aux1d
	ja lasched0_ij_loop
lasched0_ijloop_end:
	prefetch 64(ij_ptr)
	subl ij_ubd,ijd
	leaq 8(ri),ri
	movl ijd,(ij_ptr)
	addl $65536,fbi_offsd
	cmpq ij_ptr,ij_ptr_ub
	movl 4(ij_ptr),ijd
	leaq 4(ij_ptr),ij_ptr
	ja lasched_fbi_loop
lasched0_fbi_loop_end:
	movq ri,%rax
	popq %rbp
	popq %rbp
	popq %r15
	popq %r14
	popq %r13
	popq %r12
	popq %rbx
	ret
