define(sieve_interval,%rdi)dnl
define(hzs_ptr,%rsi)dnl
define(hzs_ptr_ub,%rdx)dnl
define(srb_ptr_arg,%rcx)dnl
dnl %rcx is needed for shift operations
define(srb_ptr,%r11)dnl
define(srb_ptr_ub,%r8)dnl
define(cand_array,%r9)dnl
define(csv_arg,48(%rsp))dnl
define(si_ptr,%r10)dnl
define(si_ptrl,%r10d)dnl
define(si_ptrw,%r10w)dnl
define(srb,%al)dnl
define(srb1,%bl)dnl
define(auxreg1,%rax)dnl
define(auxreg2_byte,%cl)dnl
define(auxreg2,%rcx)dnl
define(auxreg,%rbx)dnl
define(si_index,%r12)dnl
define(si_ptr1,%r13)dnl
define(csv,%r14)dnl
define(si_ptr_ub,%r15)
define(cs_steps,128)dnl
define(cs_2steps,eval(2*cs_steps))dnl
define(n_i_minus_cssteps,eval(n_i-cs_steps))dnl
function_head(lasieve_search0)
	pushq %rbx
	pushq %r12
	pushq %r13
	pushq %r14
	pushq %r15
	movq sieve_interval,si_ptr
	leaq cs_steps`'(sieve_interval),si_ptr_ub
	movq csv_arg,csv
	movq srb_ptr_arg,srb_ptr
	movq si_ptr,si_ptr1
search_j_loop:
	movb (srb_ptr),srb
	movb (hzs_ptr),srb1
	incb srb1
	leaq 1(hzs_ptr),hzs_ptr
	subb srb1,srb
	jb store_many_sievereports
	movzbq srb,auxreg2
	movd auxreg2,%xmm7
	punpcklbw %xmm7,%xmm7
	xorq auxreg,auxreg
	punpcklwd %xmm7,%xmm7
	pshufd $0,%xmm7,%xmm7
search_i_loop:
	movdqa   (si_ptr),%xmm0
	movdqa 16(si_ptr),%xmm1
	movdqa %xmm7,%xmm4
	movdqa 32(si_ptr),%xmm2
	pmaxub %xmm0,%xmm4
	movdqa 48(si_ptr),%xmm3
	movdqa %xmm2,%xmm5
	pmaxub %xmm1,%xmm4
	pmaxub %xmm3,%xmm5
	leaq 64(si_ptr),si_ptr
	pmaxub %xmm4,%xmm5
	pcmpeqb %xmm7,%xmm5
	pmovmskb %xmm5,auxreg
	cmpq $65535,auxreg
	jne search_sievereport
i_loop_entry2:
	cmpq si_ptr,si_ptr_ub
	ja search_i_loop
i_loop_ende:
	cmpq hzs_ptr,hzs_ptr_ub
	leaq n_i`'(si_ptr_ub),si_ptr_ub
	leaq n_i_minus_cssteps`'(si_ptr),si_ptr
	ja search_j_loop
cs_ret2l0:
	leaq -j_per_strip`'(hzs_ptr_ub),hzs_ptr
	leaq 1(srb_ptr),srb_ptr
	movq si_ptr1,si_ptr
	cmpq srb_ptr,srb_ptr_ub
	leaq cs_2steps`'(si_ptr),si_ptr_ub
	leaq cs_steps`'(si_ptr),si_ptr
	movq si_ptr,si_ptr1
	ja search_j_loop
lss_ende:
	movq csv,%rax
	subq csv_arg,%rax
	popq %r15
	popq %r14
	emms
	popq %r13
	popq %r12
	popq %rbx
	ret

.align 4
search_sievereport:
	pmaxub %xmm7,%xmm0
	pmaxub %xmm7,%xmm1
	pmaxub %xmm7,%xmm2
	pmaxub %xmm7,%xmm3
	pcmpeqb %xmm7,%xmm0
	pcmpeqb %xmm7,%xmm1
	pcmpeqb %xmm7,%xmm2
	pcmpeqb %xmm7,%xmm3
	pmovmskb %xmm0,auxreg1
	pmovmskb %xmm2,auxreg2
	salq $32,auxreg2
	pmovmskb %xmm1,auxreg
	orq auxreg2,auxreg1
	pmovmskb %xmm3,auxreg2
	salq $16,auxreg
	salq $48,auxreg2
	orq auxreg,auxreg1
	xorq auxreg,auxreg
	orq auxreg2,auxreg1
	leaq -64(si_ptr),si_ptr
	xorq $0xffffffffffffffff,auxreg1
	subq sieve_interval,si_ptr
loop1:	
	bsfq auxreg1,auxreg2
	addq auxreg2,auxreg
	addq auxreg,si_ptr
	shrq auxreg2_byte,auxreg1
	movw si_ptrw,(cand_array)
	leaq 2(cand_array),cand_array 
	movb (sieve_interval,si_ptr),auxreg2_byte
	subq auxreg,si_ptr
	incq auxreg
	shrq $1,auxreg1
	testq auxreg1,auxreg1
	movb auxreg2_byte,(csv)
	leaq 1(csv),csv
	jnz loop1
	leaq 64(si_ptr),si_ptr
	leaq (sieve_interval,si_ptr),si_ptr
	jmp i_loop_entry2

.align 4
store_many_sievereports:
	subq sieve_interval,si_ptr
	movq si_ptr,auxreg
	shlq $16,auxreg
	addq si_ptr,auxreg
	addq $0x00010000,auxreg
	movd auxreg,%xmm0
	addq $0x00020002,auxreg
	movd auxreg,%xmm1
	addq $0x00020002,auxreg
	movq $0x00080008,auxreg1
	psllq $32,%xmm1
	movd auxreg,%xmm2
	addq $0x00020002,auxreg
	movd auxreg1,%xmm3
	por %xmm1,%xmm0
	addq sieve_interval,si_ptr
	movd auxreg,%xmm1
	movq %xmm3,%xmm4
	psllq $32,%xmm1
	psllq $32,%xmm4
	por %xmm2,%xmm1
	por %xmm4,%xmm3
many_sr_loop:
	movq (si_ptr),%xmm2
	movq 8(si_ptr),%xmm4
	movq %xmm0,(cand_array)
	leaq 16(si_ptr),si_ptr
	paddw %xmm3,%xmm0
	movq %xmm1,8(cand_array)
	cmpq si_ptr,si_ptr_ub
	paddw %xmm3,%xmm1
	movq %xmm2,(csv)
	movq %xmm0,16(cand_array)
	paddw %xmm3,%xmm0
	movq %xmm4,8(csv)
	leaq 16(csv),csv
	movq %xmm1,24(cand_array)
	leaq 32(cand_array),cand_array
	paddw %xmm3,%xmm1
	ja many_sr_loop
	jmp i_loop_ende
