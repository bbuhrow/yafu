# Copyright (C) 2002 Jens Franke, T.Kleinjung
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

define(tmp_j,%r10d)dnl
define(k32,%r11d)dnl

function_head(asm_gauss)
whileloop:
	movl mpqs_gauss_k(%rip),%eax
	decl %eax
	movl $1,%r8d
	movl %eax,mpqs_gauss_k(%rip)
	jl end
	movl %eax,%ecx
	andl $31,%ecx
	shll %cl,%r8d
	movl %eax,%ecx
	shrl $5,%ecx      # k32
	movl %ecx,k32

# find j:
	movl mpqs_gauss_j(%rip),%ecx
	movq $mpqs_gauss_c,%rdi
 
	movq mpqs_gauss_row(%rip),%rsi
	movl k32,%eax
	movl mpqs_gauss_n32(%rip),%edx
	movq (%rsi,%rcx,8),%rsi
	subq %rdx,%rax
	decl %ecx
	leaq (%rsi,%rax,4),%rsi
	movl mpqs_gauss_n32(%rip),%eax

.align 16 
loopj:        # ? cycles per iteration
	incl %ecx
entryj:
        cmpl mpqs_gauss_m(%rip),%ecx
        jnc  end_of_loop
	leaq (%rsi,%rax,4),%rsi
	movl (%rsi),%edx
	andl %r8d,%edx
	jz loopj

endj:       # j in %rcx

	cmpl mpqs_gauss_j(%rip),%ecx
	jz noxch_j
# exchange rows j and mpqs_gauss_j: uses %rax, %rsi, %r9, %rdx

	movq mpqs_gauss_row(%rip),%rsi
	movq (%rsi,%rcx,8),%rax    # mpqs_gauss_row[j]
	movl mpqs_gauss_j(%rip),%edx
	movq (%rsi,%rdx,8),%r9    # mpqs_gauss_row[mpqs_gauss_j]
	xorq %rdx,%rdx
xch_loop:
	movq (%rax,%rdx,4),%mm0
	movq (%r9,%rdx,4),%mm1
	movq %mm0,(%r9,%rdx,4)
	movq %mm1,(%rax,%rdx,4)
	addq $2,%rdx
	cmpl mpqs_gauss_n32(%rip),%edx
	jc xch_loop

xch_loop_end:
	movl mpqs_gauss_j(%rip),%ecx
noxch_j:
	movq $mpqs_gauss_d,%rsi
	movl mpqs_gauss_k(%rip),%eax
	movw %cx,(%rsi,%rax,2)
	movw %ax,(%rdi,%rcx,2)

# update mpqs_gauss_j:
	incl mpqs_gauss_j(%rip)
	movl %ecx,tmp_j

dnl	movq tmp_j,%rcx
	movq mpqs_gauss_mat(%rip),%rsi
	movq $mpqs_gauss_col,%rdi
	movl k32,%edx
	testl %ecx,%ecx
	movl mpqs_gauss_n32(%rip),%eax
	leaq (%rsi,%rdx,4),%rsi
	movq $2,%rdx
	jz searchloopAend
	xorl %ecx,%ecx
.align 16
searchloopA:
	movl (%rsi),%r9d
	andl %r8d,%r9d
	movq $0,%r9
	movw %cx,(%rdi)
	cmovnzq %rdx,%r9
	incl %ecx
	addq %r9,%rdi
	cmpl tmp_j,%ecx
	leaq (%rsi,%rax,4),%rsi
	jc searchloopA

searchloopAend:
	incl %ecx
	cmpl mpqs_gauss_m(%rip),%ecx
	leaq (%rsi,%rax,4),%rsi
	jnc searchloopBend
.align 16
searchloopB:
	movl (%rsi),%r9d
	andl %r8d,%r9d
	movq $0,%r9
	movw %cx,(%rdi)
	cmovnzq %rdx,%r9
	incq %rcx
	addq %r9,%rdi
	cmpl mpqs_gauss_m(%rip),%ecx
	leaq (%rsi,%rax,4),%rsi
	jc searchloopB

searchloopBend:
	movq $mpqs_gauss_col,%r8
	cmpq %r8,%rdi
	jz whileloop
	movl tmp_j,%ecx
	movq mpqs_gauss_row(%rip),%rsi
	movq (%rsi,%rcx,8),%rsi

	movl k32,%edx
	cmpl $2,%edx
	jc entry1
	cmpl $4,%edx
	jc entry2
	cmpl $6,%edx
	jc entry3
# generic code for >192 bit
.align 16
outerloop0:
	movzwq (%r8),%r9
	movq mpqs_gauss_row(%rip),%rax
	movq (%rax,%r9,8),%rax
	xorl %ecx,%ecx
innerloop0:
	movl (%rsi,%rcx,4),%edx
	xorl %edx,(%rax,%rcx,4)
	incl %ecx
	cmpl %ecx,k32
	jnc innerloop0

	leaq 2(%r8),%r8
	cmpq %rdi,%r8
	jc outerloop0


	jmp whileloop



end:
	movl tmp_j,%eax
	emms
	ret

end_of_loop:
	movq $mpqs_gauss_d,%rdi
	movl mpqs_gauss_k(%rip),%eax
	movw $-1,%r9w
	movw %r9w,(%rdi,%rax,2)
	movl %ecx,tmp_j
	jmp end

# code for length <=64 bit
entry1:
	movzwq (%r8),%r9
	movq (%rsi),%mm0
	movq mpqs_gauss_row(%rip),%rsi
	movq (%rsi,%r9,8),%rax
outerloop1:
	movq (%rax),%mm1
	leaq 2(%r8),%r8
	movzwq (%r8),%r9
	pxor %mm0,%mm1
	cmpq %rdi,%r8
	movq %mm1,(%rax)
	movq (%rsi,%r9,8),%rax
	jc outerloop1
	jmp whileloop

# code for length <=128 bit
entry2:
	movq (%rsi),%mm0
	movzwq (%r8),%r9
	movq 8(%rsi),%mm2
	movq mpqs_gauss_row(%rip),%rsi
	movq (%rsi,%r9,8),%rax
outerloop2:
	movq (%rax),%mm1
	leaq 2(%r8),%r8
	movq 8(%rax),%mm3
	pxor %mm0,%mm1
	movzwq (%r8),%r9
	pxor %mm2,%mm3
	cmpq %rdi,%r8
	movq %mm1,(%rax)
	movq %mm3,8(%rax)
	movq (%rsi,%r9,8),%rax
	jc outerloop2
	jmp whileloop

# code for length <=192 bit
entry3:
	movq (%rsi),%mm0
	movzwq (%r8),%r9
	movq 8(%rsi),%mm2
	movq 16(%rsi),%mm4
	movq mpqs_gauss_row(%rip),%rsi
	movq (%rsi,%r9,8),%rax
outerloop3:
	movq (%rax),%mm1
	leaq 2(%r8),%r8
	movq 8(%rax),%mm3
	pxor %mm0,%mm1
	movq 16(%rax),%mm5
	movzwq (%r8),%r9
	pxor %mm2,%mm3
	movq %mm1,(%rax)
	pxor %mm4,%mm5
	movq %mm3,8(%rax)
	cmpq %rdi,%r8
	movq %mm5,16(%rax)
	movq (%rsi,%r9,8),%rax
	jc outerloop3
	jmp whileloop
