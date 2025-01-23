	.global mulmod128
	.type mulmod128, @function
	.align 4

mulmod128:
	li t4, 0x00000000ffffffff
	and t0, a0, t4
	srli t1, a0, 32

	add t6, x0, t0
        add a3, x0, t1

	srli t5, a1, 32
	and a1, a1, t4

	mul a0, a1, t6
	srli t2, a0, 32
	and t0, a0, t4

	mul a0, t5, t6
	srli t6, a0, 32
	and t3, a0, t4

	mul a0, a3, a1
	add a0, a0, t2
	add a0, a0, t3
	srli t2, a0, 32
	and t1, a0, t4

	mul a0, t5, a3
	add a0, a0, t2
	add a0, a0, t6
	srli t3, a0, 32
	and t2, a0, t4

	srli a3, a2, 32
	and t6, a2, t4

first1:
	beq x0, t3, next1
	mul a0, t3, t6
	add a0, a0, t1
	srli t5, a0, 32
	and t1, a0, t4
	add t2, t5, t2
	mul a0, t3, a3
	add t2, t2, a0
	srli t3, t2, 32
	and t2, t2, t4
	j first1
next1:
	beq x0, t2, last
	mul a0, t2, t6
	add a0, a0, t0
	srli t5, a0, 32 
	and t0, a0, t4
	add t1, t5, t1
	mul a0, t2, a3
	add t1, t1, a0
	srli t2, t1, 32
	and t1, t1, t4
	j next1
last:
	slli t1, t1, 32
	add a0, t0, t1
	ret
