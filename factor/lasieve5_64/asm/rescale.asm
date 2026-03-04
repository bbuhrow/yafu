dnl Copyright (C) 2002,2004 Jens Franke, T.Kleinjung
dnl This file is part of gnfs4linux, distributed under the terms of the
dnl GNU General Public Licence and WITHOUT ANY WARRANTY.

dnl You should have received a copy of the GNU General Public License along
dnl with this program; see the file COPYING.  If not, write to the Free
dnl Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.

dnl rescale_interval1(array,length):
dnl replaces array[i] by (array[i]+1)/2, 0<=i<length
dnl
dnl rescale_interval2(array,length):
dnl replaces array[i] by (array[i]+3)/4, 0<=i<length

define(array,%rdi)dnl
define(length,%rsi)dnl
define(len16,%rcx)dnl
define(x0,%xmm0)dnl
define(x1,%xmm1)dnl

function_head(rescale_interval1)
	movq length,len16
	shrq $4,len16
	testq len16,len16
	jz loop1_1
	xorps x0,x0
loop1_16:
	movaps (array),x1
	pavgb x0,x1
	movaps x1,(array)
	leaq 16(array),array
	decq len16
	jnz loop1_16
	andq $15,length
loop1_1:
	testq length,length
	jz loopend1
	movzbl (array),%eax
	incl %eax
	shrl $1,%eax 
	movb %al,(array)
	incq array
	decq length
	jmp loop1_1
loopend1:
	emms
	ret

function_head(rescale_interval2)
	movq length,len16
	shrq $4,len16
	testq len16,len16
	jz loop2_1
	xorps x0,x0
loop2_16:
	movaps (array),x1
	pavgb x0,x1
	pavgb x0,x1
	movaps x1,(array)
	leaq 16(array),array
	decq len16
	jnz loop2_16
	andq $15,length
loop2_1:
	testq length,length
	jz loopend2
	movzbl (array),%eax
	incl %eax
	shrl $1,%eax 
	movb %al,(array)
	incq array
	decq length
	jmp loop2_1
loopend2:
	emms
	ret

function_head(rescale_interval1_noxmm)
	testq length,length
	jz loopend1n
loop1n:
	movzbl (array),%eax
	incl %eax
	shrl $1,%eax 
	movb %al,(array)
	incq array
	decq length
	jnz loop1n
loopend1n:
	ret

function_head(rescale_interval2_noxmm)
	testq length,length
	jz loopend2n
loop2n:
	movzbl (array),%eax
	addl $3,%eax
	shrl $2,%eax 
	movb %al,(array)
	incq array
	decq length
	jnz loop2n
loopend2n:
	ret

