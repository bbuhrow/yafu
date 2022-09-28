dnl Copyright (C) 2004 Jens Franke, T.Kleinjung
dnl This file is part of gnfs4linux, distributed under the terms of the
dnl GNU General Public Licence and WITHOUT ANY WARRANTY.

dnl You should have received a copy of the GNU General Public License along
dnl with this program; see the file COPYING.  If not, write to the Free
dnl Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.

dnl The function we want to write.
dnl ulong asm_modinv32(x)
dnl Modular inverse of x modulo modulo32
dnl and y satisfies 0<y<x.

define(x,%edi)dnl
define(y,%esi)dnl
define(xc,%ecx)dnl
define(yc,%r8d)dnl

dnl Number of trial subtractions before doing a division
define(nts,15)dnl

function_head(asm_modinv32)
	movl modulo32(%rip),y
	testl x,x
	jz badargs
	cmpl x,y
	jbe badargs
dnl Set xc and yc to their initial values.
	xorl yc,yc
	xorl xc,xc
	incl xc
	cmpl $1,x
	jbe have_inverse2
divide:

forloop(`i',1,nts,`
	subl x,y
	addl xc,yc
	cmpl x,y
	jb xlarger'
)
	movl y,%eax
	xorl %edx,%edx
	divl x
	movl %edx,y
	mull xc
	addl %eax,yc

xlarger:
	cmpl $1,y
	jbe have_inverse1

forloop(`i',1,nts,`
	subl y,x
	addl yc,xc
	cmpl y,x
	jb ylarger'
)

        movl x,%eax
	xorl %edx,%edx
	divl y
	movl %edx,x
	mull yc
	addl %eax,xc
ylarger:
	cmpl $1,x
	ja divide
have_inverse2:
	jne badargs
	movl xc,%eax
	ret

have_inverse1:
	jne badargs
	movl modulo32(%rip),%eax
	subl yc,%eax
	ret
badargs:
	call abort
	.END
