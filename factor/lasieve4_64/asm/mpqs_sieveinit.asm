# Copyright (C) 2002 Jens Franke, T.Kleinjung
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.



define(sv,%rdi)dnl
define(len,%rsi)dnl
define(tab_arg,%rdx)dnl
define(maskptr,%rcx)dnl
define(mask1,%rcx)dnl
define(tinyptr,%r8)dnl
define(tinylen,%r9)dnl
define(tinylen4,%r8)dnl
define(ilen,%r10)dnl
define(tptr,%r11)dnl
define(tptrend,%r12)dnl
define(tab,%rbx)dnl
define(mult,%rax)dnl
define(llen,%rdx)dnl
define(mask,%mm2)dnl

function_head(asm_sieve_init)
	pushq %r12
	pushq %rbx
	movq tab_arg,tab
	movq (maskptr),mask1
	movq tinyptr,tptr
	movq tinyptr,tptrend
	movq tinylen,tinylen4
	shlq $4,tinylen4
	addq tinylen4,tptrend  # tinyend
	movzwq (tab),llen     # this is 0
outerloop:
	movzwq 4(tab),ilen
	subq llen,ilen       # length in byte
	shrq $4,ilen         # length in 32 byte
	movzwq 2(tab),mult
	leaq 4(tab),tab
	mulq mask1
	movd %rax,mask
.align 16
innerloop:
	movq (tptr),%mm0
	movq 8(tptr),%mm1
	leaq 16(tptr),tptr
	xorq %rax,%rax
	cmpq tptrend,tptr
	paddb mask,%mm0
	paddb mask,%mm1
	cmovncq tinylen4,%rax
	subq %rax,tptr
	decq ilen            # length>0
	movq %mm0,(sv)
	movq %mm1,8(sv)
	leaq 16(sv),sv

	jnz innerloop

	movzwq (tab),llen
	cmpq len,llen
	jnz outerloop

	emms
	popq %rbx
	popq %r12
	ret


function_head(asm_sieve_init16)
	pushq %r12
	pushq %rbx
	movq tab_arg,tab
	movq (maskptr),mask1
	movq tinyptr,tptr
	movq tinyptr,tptrend
	movq tinylen,tinylen4
	shlq $4,tinylen4
	addq tinylen4,tptrend  # tinyend
	movzwq (tab),llen     # this is 0
outerloop16:
	movzwq 4(tab),ilen
	subq llen,ilen       # length in byte
	shrq $5,ilen         # length in 32 byte
	movzwq 2(tab),mult
	leaq 4(tab),tab
	mulq mask1
	movd %rax,mask
.align 16
innerloop16:
	movq (tptr),%mm0
	movq 8(tptr),%mm1
	movq 16(tptr),%mm4
	movq 24(tptr),%mm5
	leaq 32(tptr),tptr
	paddb mask,%mm0
	paddb mask,%mm1
	cmpq tptrend,tptr
	movq $0,%rax
	cmovncq tinylen4,%rax
	paddb mask,%mm4
	paddb mask,%mm5
	subq %rax,tptr
	decq ilen            # length>0
	movq %mm0,(sv)
	movq %mm1,8(sv)
	movq %mm4,16(sv)
	movq %mm5,24(sv)
	leaq 32(sv),sv

	jnz innerloop16

	movzwq (tab),llen
	cmpq len,llen
	jnz outerloop16

	emms
	popq %rbx
	popq %r12
	ret


