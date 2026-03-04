define(aux0,%esi)dnl
define(aux0q,%rsi)dnl
define(aux1,%r8d)dnl
define(aux1q,%r8)dnl
define(aux2,%r9d)dnl
define(prime_buffer,%rdi)dnl
define(strip_i_arg,aux0)dnl
define(aux_ptr,%rdx)dnl
define(aux_ptr_ub_arg,%rcx)dnl
dnl %ecx is used directly as it has special relation to shift operations
define(aux_ptr_ub,%rax)dnl
function_head(asm_MMX_Td8)
dnl	movl strip_i_arg,aux0
	movl strip_i_arg,aux1
	movdqa (aux_ptr),%xmm2
	shll $16,aux1
	movdqa 16(aux_ptr),%xmm0
	orl aux1,aux0
	movdqa 32(aux_ptr),%xmm1
	leaq 48(aux_ptr),aux_ptr
	movd %rsi,%xmm3
	movq aux_ptr_ub_arg,aux_ptr_ub
	leaq 24(aux_ptr_ub),aux_ptr_ub
	pshufd $0,%xmm3,%xmm3
	pxor %xmm4,%xmm4
	paddw %xmm3,%xmm2
	movd %xmm3,%r10d
MMX_TdLoop:
dnl	prefetcht0 96(aux_ptr)
	cmpq aux_ptr,aux_ptr_ub
	pmullw %xmm2,%xmm1
	movdqa (aux_ptr),%xmm2
	jbe MMX_TdEnde
	pmulhw %xmm1,%xmm0
	movdqa 32(aux_ptr),%xmm1
	pcmpeqw %xmm4,%xmm0
	paddw %xmm3,%xmm2
	pmovmskb %xmm0,aux0
	movdqa 16(aux_ptr),%xmm0
	leaq 48(aux_ptr),aux_ptr
	testl aux0,aux0
	jz MMX_TdLoop
	andl $0x5555,aux0
	xorq aux1q,aux1q
MMX_TdLoop1:
	bsfq aux0q,%rcx
	shrl %cl,aux0
	shrq $1,%rcx
	addq %rcx,aux1q
	movzwl -80(aux_ptr,aux1q,2),aux2
	shrl $2,aux0
	movl aux2,(prime_buffer)
	incl aux1
	testl aux0,aux0
	leaq 4(prime_buffer),prime_buffer
	jz MMX_TdLoop
	jmp MMX_TdLoop1
.align 4
MMX_TdEnde:
	movq prime_buffer,%rax
	ret

define(`auxptr_arg',%rdi)dnl
define(`auxptr_ub_arg',%rsi)dnl
define(`uptr_arg',%rdx)dnl
define(`auxptr',%rax)dnl
define(`auxptr_ub',%rcx)dnl
define(`uptr',%r8)dnl
function_head(asm_TdUpdate8)
	movq auxptr_arg,auxptr
	movq auxptr_ub_arg,auxptr_ub
	cmpq auxptr,auxptr_ub_arg
	movq uptr_arg,uptr
	leaq -48(auxptr_ub),auxptr_ub
	movdqa (auxptr),%xmm2
	movdqa (uptr),%xmm0
	leaq 16(uptr),uptr
	movdqa 16(auxptr),%xmm1
	movdqa %xmm0,%xmm3
	jbe TdUpdateEnde
TdUpdateLoop:
dnl	prefetcht0 96(aux_ptr)
	pcmpgtw %xmm2,%xmm3
	psubw %xmm0,%xmm2
	cmpq auxptr,auxptr_ub
	movdqa (uptr),%xmm0
	pand %xmm1,%xmm3
	movdqa 64(auxptr),%xmm1
dnl	prefetcht0 32(uptr)
	leaq 16(uptr),uptr
	paddw %xmm2,%xmm3
	movdqa 48(auxptr),%xmm2
	movdqa %xmm3,(auxptr)
	movdqa %xmm0,%xmm3
	leaq 48(auxptr),auxptr
	ja TdUpdateLoop
TdUpdateEnde:
	ret
