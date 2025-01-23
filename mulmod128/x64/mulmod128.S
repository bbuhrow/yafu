	.global mulmod128
	.type mulmod128, @function
	.align 4

mulmod128:
	push %rbp
	mov %rsp, %rbp
	push %rdx

	mov $0x00000000ffffffff, %r8
	mov %rdi, %r9
	and %r8, %r9
	mov %rdi, %r10
	shr $32, %r10
	mov %r9, %r14
        mov %r10, %r15
	mov %rsi, %r13
	and %r8, %rsi
	shr $32, %r13

	mov %r14, %rax
	mul %rsi
	mov %rax, %r9
	and %r8, %r9
	mov %rax, %r11
	shr $32, %r11

	mov %r14, %rax
	mul %r13
	mov %rax, %r14
	shr $32, %r14
	mov %rax, %r12
	and %r8, %r12

        mov %r15, %rax
	mul %rsi
        mov %rax, %r10
	add %r11, %r10
	add %r12, %r10
        mov %r10, %r11
	shr $32, %r11
	and %r8, %r10

	mov %r15, %rax
	mul %r13
        mov %rax, %r13
	add %r13, %r11
	add %r14, %r11
	mov %r11, %r12
	shr $32, %r12
	and %r8, %r11

	pop %rdx

	mov %rdx, %r15
	mov %rdx, %r14
	shr $32, %r15
        and %r8, %r14
first1:
	cmp $0, %r12
	je next1
	mov %r14, %rax
	mul %r12
	add %rax, %r10
	mov %r10, %r13
	and %r8, %r10
	shr $32, %r13
	add %r13, %r11
	mov %r15, %rax
	mul %r12
	add %rax, %r11
	mov %r11, %r12
	and %r8, %r11
	shr $32, %r12
	jmp first1
next1:
	cmp $0, %r11
	je last
	mov %r14, %rax
	mul %r11
	add %rax, %r9
	mov %r9, %r13
	and %r8, %r9

	shr $32, %r13
	add %r13, %r10
	mov %r15, %rax
	mul %r11
	add %rax, %r10
	mov %r10, %r11
	and %r8, %r10
	shr $32, %r11
	jmp next1
last:
	shl $32, %r10
	add %r10, %r9

	mov %r9, %rax
	mov %rbp, %rsp
	pop %rbp
	ret

