    
    BITS 64
    
	global mod_64

mod_64:
    mov     rax, rcx
    mov     rcx, rdx
    xor      rdx, rdx
    div      rcx
    mov     rax, rdx
    ret
