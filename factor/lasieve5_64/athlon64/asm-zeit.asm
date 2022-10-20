function_head(asmgetclock)
	rdtsc
	shlq $32,%rdx
	orq %rdx,%rax
	ret

function_head(zeitA)
	rdtsc
	shlq $32,%rdx
	orq %rdx,%rax
	movq asmzeitcounter(%rip), %rdx
	subq %rax,(%rdx,%rdi,8)
	ret

function_head(zeitB)
	rdtsc
	shlq $32,%rdx
	orq %rdx,%rax
	movq asmzeitcounter(%rip), %rdx
	addq %rax,(%rdx,%rdi,8)
	ret
