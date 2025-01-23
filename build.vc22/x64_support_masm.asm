
;   uint64 mul_mod(uint64 x, uint64 y, uint64 m)
;   return (x * y) % m

.CODE
     
mul_mod_64 PROC
    mov     rax, rcx
    mul     rdx
    div     r8
    mov     rax, rdx
    ret

mul_mod_64 ENDP

END

