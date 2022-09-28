define(ri,%rdi)dnl
define(ij_ptr,%rsi)dnl
define(ij_ptr_ub_arg,%rdx)dnl
define(ij_ptr_ub,(%rsp))dnl
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
define(ot_mask,eval(n_i|1))dnl
define(ot_tester1,eval((ot&1)|((ot&2)*(2**(n_i_bits-1)))))dnl
define(ot_tester2,eval(n_i^ot_tester1))dnl
define(v0d,%edx)dnl
define(v1d,%ebp)dnl
ifelse(nt_sched,1,
	`function_head(lasched`'ot`'nt)',
	`function_head(lasched`'ot)')
	pushq %rbx
	shlq $n_i_bits,ij_ub
	pushq %r12
	pushq %r13
	pushq %r14
	shrq $16,fbi_offs
	pushq %r15
	cmpq ij_ptr,ij_ptr_ub_arg
	pushq %rbp
	leaq -4(ij_ptr_ub_arg),ij_ptr_ub_arg
	pushq ij_ptr_ub_arg
	movl (ri),ad
	movl 4(ri),bd
	jbe lasched`'ot`'_fbi_loop_end
lasched`'ot`'_fbi_loop:
	movl ad,aux0d
	movl bd,aux1d
	andl $ot_mask,aux0d
	andl $ot_mask,aux1d
	cmpl $ot_tester1,aux0d
	movl $0,aux2d
	movl ad,v0d
	cmovnel bd,aux2d
	cmpl $ot_tester2,aux1d
	movl $0,ijd
	movl bd,v1d
	cmovnel ad,ijd
	negl ad
	negl bd
	addl aux2d,ijd
	andl $n_i_mask,ad
	ifelse(ot,1,`addl $n_i,ijd')
	andl $n_i_mask,bd
	shrl $1,ijd
	ifelse(ot,2,`cmpl ad,bd
	jbe lasched2_exception
lasched2_afterexception:')
	cmpl ijd,ij_ubd
	movl ijd,aux0d
	movl ijd,aux1d
	prefetch 128(ri)
	jbe lasched`'ot`'_ijloop_end
lasched`'ot`'_ij_loop:
	movl ijd,aux2d
	andq $n_i_mask,aux0
	shrq $l1_bits,aux1
	andl $l1_mask,aux2d
	cmpl aux0d,ad
	movl $0,aux3d
	movl $0,aux4d
	cmovbel v1d,aux3d
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
	ja lasched`'ot`'_ij_loop
lasched`'ot`'_ijloop_end:
	prefetchw 64(ij_ptr)
	subl ij_ubd,ijd
	leaq 8(ri),ri
	movl ijd,(ij_ptr)
	addl $65536,fbi_offsd
	cmpq ij_ptr,ij_ptr_ub
	leaq 4(ij_ptr),ij_ptr
	movl (ri),ad
	movl 4(ri),bd
	ja lasched`'ot`'_fbi_loop
lasched`'ot`'_fbi_loop_end:
	movq ri,%rax
	popq %rbp
	popq %rbp
	popq %r15
	popq %r14
	popq %r13
	popq %r12
	popq %rbx
	ret

ifelse(ot,2,`
lasched2_exception:
	movl (ri),aux0d
	cmpl aux0d,4(ri)
	jb lasched2_afterexception
	cmpl ad,bd
	movl $n_i,ijd
	movl $0,aux1d
	cmovel 4(ri),ijd
	cmovel aux0d,aux1d
	subl aux1d,ijd
	shrl $1,ijd
	jmp lasched2_afterexception')
