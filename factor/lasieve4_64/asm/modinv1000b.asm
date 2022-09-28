dnl Copyright (C) 2004 Jens Franke, T.Kleinjung
dnl This file is part of gnfs4linux, distributed under the terms of the
dnl GNU General Public Licence and WITHOUT ANY WARRANTY.

dnl You should have received a copy of the GNU General Public License along
dnl with this program; see the file COPYING.  If not, write to the Free
dnl Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.

dnl The function we want to write.
dnl ulong asm_modinv32b(x,y)
dnl Modular inverse of x modulo y
dnl when x satisfies 0<x<y.
dnl If x==0 or x==y, 0 is returned
dnl All args outside these bounds,
dnl and all other non-coprime args, cause abort

define(x,%edi)dnl
define(y,%esi)dnl
define(xc,%ecx)dnl
define(yc,%r8d)dnl
define(p,%r9d)dnl
dnl Number of trial subtractions before doing a division
define(nts,15)dnl

function_head(asm_modinv32b)
	testl x,x
	movl y,p
	jz Modinv_arg0
	cmpl x,y
	jbe BadArgs1
dnl Set xc and yc to their initial values.
	xorl yc,yc
	xorl xc,xc
	incl xc
	cmpl $1,x
	jbe have_Inverse2
divide:

forloop(`i',1,nts,`
	subl x,y
	addl xc,yc
	cmpl x,y
	jb xLarger'
)
	movl y,%eax
	xorl %edx,%edx
	divl x
	movl %edx,y
	mull xc
	addl %eax,yc

xLarger:
	cmpl $1,y
	jbe have_Inverse1

forloop(`i',1,nts,`
	subl y,x
	addl yc,xc
	cmpl y,x
	jb yLarger'
)

        movl x,%eax
	xorl %edx,%edx
	divl y
	movl %edx,x
	mull yc
	addl %eax,xc
yLarger:
	cmpl $1,x
	ja divide
have_Inverse2:
	jne BadArgs
	movl xc,%eax
	ret

have_Inverse1:
	jne BadArgs
	movl p,%eax
	subl yc,%eax
	ret
BadArgs1:
	je Modinv_arg0
BadArgs:
	call abort
Modinv_arg0:
	xorl %eax,%eax
	ret
	.END
