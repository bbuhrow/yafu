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
function_head(asm_MMX_Td4)
dnl	movl strip_i_arg,aux0
	movl strip_i_arg,aux1
	movq (aux_ptr),%mm2
	shll $16,aux1
	movq 8(aux_ptr),%mm0
	orl aux1,aux0
	movq 16(aux_ptr),%mm1
	leaq 24(aux_ptr),aux_ptr
	movd aux0,%mm3
	movq aux_ptr_ub_arg,aux_ptr_ub
	psllq $32,%mm3
	movd aux0,%mm4
	leaq 24(aux_ptr_ub),aux_ptr_ub
	por %mm4,%mm3
	pxor %mm4,%mm4
	paddw %mm3,%mm2
MMX_TdLoop:
	cmpq aux_ptr,aux_ptr_ub
	pmullw %mm2,%mm1
	movq (aux_ptr),%mm2
	jbe MMX_TdEnde
	pmulhw %mm1,%mm0
	movq 16(aux_ptr),%mm1
	pcmpeqw %mm4,%mm0
	paddw %mm3,%mm2
	pmovmskb %mm0,aux0
	movq 8(aux_ptr),%mm0
	leaq 24(aux_ptr),aux_ptr
	testl aux0,aux0
	jz MMX_TdLoop
	andl $0x55,aux0
	xorq aux1q,aux1q
MMX_TdLoop1:
	bsfq aux0q,%rcx
	shrl %cl,aux0
	shrq $1,%rcx
	addq %rcx,aux1q
	movzwl -40(aux_ptr,aux1q,2),aux2
	shrl $2,aux0
	movl aux2,(prime_buffer)
	incl aux1
	testl aux0,aux0
	leaq 4(prime_buffer),prime_buffer
	jz MMX_TdLoop
	jmp MMX_TdLoop1
.align 4
MMX_TdEnde:
	emms
	movq prime_buffer,%rax
	ret

define(`auxptr_arg',%rdi)dnl
define(`auxptr_ub_arg',%rsi)dnl
define(`uptr_arg',%rdx)dnl
define(`auxptr',%rax)dnl
define(`auxptr_ub',%rcx)dnl
define(`uptr',%r8)dnl
function_head(asm_TdUpdate4)
	movq auxptr_arg,auxptr
	movq auxptr_ub_arg,auxptr_ub
	cmpq auxptr,auxptr_ub_arg
	movq uptr_arg,uptr
	leaq -24(auxptr_ub),auxptr_ub
	movq (auxptr),%mm2
	movq (uptr),%mm0
	leaq 8(uptr),uptr
	movq 8(auxptr),%mm1
	movq %mm0,%mm3
	jbe TdUpdateEnde
TdUpdateLoop:
	pcmpgtw %mm2,%mm3
	psubw %mm0,%mm2
	cmpq auxptr,auxptr_ub
	movq (uptr),%mm0
	pand %mm1,%mm3
	movq 32(auxptr),%mm1
	leaq 8(uptr),uptr
	paddw %mm2,%mm3
	movq 24(auxptr),%mm2
	movq %mm3,(auxptr)
	movq %mm0,%mm3
	leaq 24(auxptr),auxptr
	ja TdUpdateLoop
TdUpdateEnde:
	emms
	ret
