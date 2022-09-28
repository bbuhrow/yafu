define(ri,%rdi)dnl
define(ij_ptr,%rsi)dnl
define(ij_ptr_ub,%rdx)dnl
define(sched_ptr,%rcx)dnl
define(fbi_offs,%r8d)dnl
define(a,%eax)dnl
define(b,%ebx)dnl
define(ij,%r9d)dnl
define(aux1,%r10d)dnl
define(aux2,%r11d)dnl
define(sched,%r12)dnl

function_head(medsched0)
	pushq %rbx
	pushq %r12

	shll $16,fbi_offs
	cmpq ij_ptr,ij_ptr_ub
	leaq -4(ij_ptr_ub),ij_ptr_ub
	movq (sched_ptr),sched
	jbe medsched_fbiloop_end
medsched_fbiloop:
	movl (ij_ptr),ij
	movl (ri),a
	prefetcht0 128(ri)
	cmpl $l1_size,ij
	movl 4(ri),b
	jae  medsched_ij_loop_end
	negl a
	negl b
	movl ij,aux1
	xorl aux2,aux2
	andl $n_i_mask,a
	andl $n_i_mask,b
medsched_ij_loop:
	andl $n_i_mask,aux1
	orl  fbi_offs,ij
	cmpl aux1,a
	movl ij,(sched)
	leaq 4(sched),sched
	cmovbel 4(ri),aux2
	andl $l1_mask,ij
	cmpl aux1,b
	movl $0,aux1
	cmoval  (ri),aux1
	addl aux2,ij
	addl aux1,ij
	xorl aux2,aux2
	cmpl $l1_size,ij
	movl ij,aux1
	jb medsched_ij_loop
medsched_ij_loop_end:
	prefetcht0 128(ij_ptr)
	subl $l1_size,ij
	addl $65536,fbi_offs
	cmpq ij_ptr,ij_ptr_ub
	movl ij,(ij_ptr)
	leaq 8(ri),ri
	leaq 4(ij_ptr),ij_ptr
	ja medsched_fbiloop
medsched_fbiloop_end:
	movq sched,(sched_ptr)
	popq %r12
	movq ri,%rax
	popq %rbx
	ret
