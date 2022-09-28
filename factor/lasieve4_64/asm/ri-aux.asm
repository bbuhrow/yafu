# Copyright (C) 2002 Jens Franke, T.Kleinjung
dnl This file is part of gnfs4linux, distributed under the terms of the
dnl GNU General Public Licence and WITHOUT ANY WARRANTY.

dnl You should have received a copy of the GNU General Public License along
dnl with this program; see the file COPYING.  If not, write to the Free
dnl Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.

dnl The function we want to write.
dnl ulong asm_getbc(x)
dnl Modular inverse of x modulo modulo32
dnl and y satisfies 0<y<x.

define(x,%edi)dnl
define(y,%esi)dnl
define(Aarg,%edx)dnl
define(at,%rcx)dnl
define(bt,%r8)dnl
define(ct,%r9)dnl
define(dt,16(%rsp))dnl
define(xc,%r10d)dnl
define(yc,%r11d)dnl
define(A,%ebx)dnl

dnl Number of trial subtractions before doing a division
define(nts,15)dnl


.text
	.align 4
.globl asm_getbc
	.type	 asm_getbc,@function	
asm_getbc:
	pushq %rbx
	xorl xc,xc
	movl Aarg,A
	xorl yc,yc
	incl xc
	cmpl x,A
	ja   have_bs
divide:

forloop(`i',1,nts,`
	subl x,y
	addl xc,yc
	cmpl x,y
	jb test_y'
)
	movl y,%eax
	xorl %edx,%edx
	divl x
	movl %edx,y
	mull xc
	addl %eax,yc

test_y:
	cmpl y,A
	ja have_ct

forloop(`i',1,nts,`
	subl y,x
	addl yc,xc
	cmpl y,x
	jb test_x'
)

        movl x,%eax
	xorl %edx,%edx
	divl y
	movl %edx,x
	mull yc
	addl %eax,xc
test_x:
	cmpl x,A
	jbe divide
have_bs:
forloop(`i',1,nts,`
	subl x,y
	addl xc,yc
	cmpl y,A
	ja have_bsct'
)
	movl y,%eax
	xorl %edx,%edx
	subl A,%eax
	divl x
	incl %eax
	movl %eax,A
	mull x
	subl %eax,y
	movl A,%eax
	mull xc
	addl %eax,yc
	jmp have_bsct
have_ct:
forloop(`i',1,nts,`
	subl y,x
	addl yc,xc
	cmpl x,A
	ja have_bsct'
)
	movl x,%eax
	xorl %edx,%edx
	subl A,%eax
	divl y
	incl %eax
	movl %eax,A
	mull y
	subl %eax,x
	movl A,%eax
	mull yc
	addl %eax,xc
have_bsct:
	movq dt,%rax
	movl x,(at)
	movl xc,(bt)
	movl y,(ct)
	movl yc,(%rax)
	popq %rbx
	ret
