	.global mulmod128
	.type mulmod128, @function
	.align 4

mulmod128:
	mov x7, 0x00000000ffffffff
	and x3, x0, x7
	mov x4, x0, lsr #32

	and x9, x3, x7
        and x10, x4, x7

	mov x8, x1, lsr #32
	and x1, x1, x7

	mul x0, x1, x9
	mov x5, x0, lsr #32
	and x3, x0, x7

	mul x0, x8, x9
	mov x9, x0, lsr #32
	and x6, x0, x7

	mul x0, x10, x1
	add x0, x0, x5
	add x0, x0, x6
	mov x5, x0, lsr #32
	and x4, x0, x7

	mul x0, x8, x10
	add x0, x0, x5
	add x0, x0, x9
	mov x6, x0, lsr #32
	and x5, x0, x7

	mov x10, x2, lsr #32
	and x9, x2, x7

first1:
	cmp x0, x6
	beq next1
	mul x0, x6, x9
	add x0, x0, x4
	mov x8, x0, lsr #32
	and x4, x0, x7
	add x5, x8, x5
	mul x0, x6, x10
	add x5, x5, x0
	mov x6, x5, lsr #32
	and x5, x5, x7
	b first1
next1:
	cmp x0, x5
	beq last
	mul x0, x5, x9
	add x0, x0, x3
	mov x8, x0, lsr #32 
	and x3, x0, x7
	add x4, x8, x4
	mul x0, x5, x10
	add x4, x4, x0
	mov x5, x4, lsr #32
	and x4, x4, x7
	b next1
last:
	mov x4, x4, lsl #32
	add x0, x3, x4
	ret

