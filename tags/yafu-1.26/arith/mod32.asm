
    text
    bits 32
    global mod_32

mod_32:
    mov     eax, [esp+4]
    xor      edx, edx
    div      dword [esp+8]
    mov     eax, edx
    ret
