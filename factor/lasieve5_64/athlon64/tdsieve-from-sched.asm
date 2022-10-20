define(sched,%r9)dnl
define(si0,%r10)dnl
define(si1,%r11)dnl
define(svo,%al)dnl
define(sched_arg,%rdi)dnl
define(sched_ub,%rsi)dnl
define(sieve_interval,%rdx)dnl
define(sched_tds_buffer,%rcx)dnl
define(sched_tds_buf_ub,%r8)dnl
function_head(tdsieve_sched2buf)
	movq (sched_arg),sched
	leaq -12(sched_ub),sched_ub
	cmpq sched,sched_ub
	jbe sched_tds_loop_ende
sched_tds_loop:
	prefetcht0 128(sched)
	movzwq (sched),si0
	movzwq 4(sched),si1
	movb (sieve_interval,si0),svo
	movzwq 8(sched),si0
	orb (sieve_interval,si1),svo
	movzwq 12(sched),si1
	orb (sieve_interval,si0),svo
	orb (sieve_interval,si1),svo
	testb svo,svo
	jnz stds_search_divisor
sched_tds_loop_entry1:
	leaq 16(sched),sched
	cmpq sched,sched_ub
	ja sched_tds_loop
sched_tds_loop_ende:
	leaq 8(sched_ub),sched_ub
	cmpq sched,sched_ub
	leaq 4(sched_ub),sched_ub
	movzwq (sched),si0
	jbe sched_tds_last_test
	movb (sieve_interval,si0),svo
	movzwq 4(sched),si1
	testb svo,svo
	movb (sieve_interval,si1),svo
	jz sched_tds_2last_test
	movq sched,(sched_tds_buffer)
	leaq 8(sched_tds_buffer),sched_tds_buffer
sched_tds_2last_test:
	testb svo,svo
	movzwq 8(sched),si0
	leaq 4(sched),si1
	leaq 8(sched),sched
	jz sched_tds_last_test
	movq si1,(sched_tds_buffer)
	leaq 8(sched_tds_buffer),sched_tds_buffer
sched_tds_last_test:
	cmpq sched,sched_ub
	jbe sched_tds_return
	movb (sieve_interval,si0),svo
	testb svo,svo
	jz sched_tds_ende
	movq sched,(sched_tds_buffer)
	leaq 8(sched_tds_buffer),sched_tds_buffer
sched_tds_ende:
	leaq 4(sched),sched
sched_tds_return:
dnl Note sched_tds_buffer is in the register containing the return value
	movq sched_tds_buffer,%rax
	movq sched,(sched_arg)
	ret

.align 4
stds_search_divisor:
	movzwq (sched),si0
	movb (sieve_interval,si0),svo
	movzwq 4(sched),si1
	testb svo,svo
	jz stds_search2
	movq sched,(sched_tds_buffer)
	leaq 8(sched_tds_buffer),sched_tds_buffer
stds_search2:
	movb (sieve_interval,si1),svo
	movzwq 8(sched),si0
	testb svo,svo
	leaq 4(sched),si1
	jz tds_search3
	movq si1,(sched_tds_buffer)
	leaq 8(sched_tds_buffer),sched_tds_buffer
tds_search3:
	movb (sieve_interval,si0),svo
	leaq 8(sched),si0
	movzwq 12(sched),si1
	testb svo,svo
        jz tds_search4
	movq si0,(sched_tds_buffer)
	leaq 8(sched_tds_buffer),sched_tds_buffer
tds_search4:
	movb (sieve_interval,si1),svo
	leaq 12(sched),si1
	testb svo,svo
	jz sched_tds_test_buffer
	movq si1,(sched_tds_buffer)
        leaq 8(sched_tds_buffer),sched_tds_buffer
sched_tds_test_buffer:
	cmpq sched_tds_buffer,sched_tds_buf_ub
	ja sched_tds_loop_entry1
	leaq 16(sched),sched
	jmp sched_tds_return
