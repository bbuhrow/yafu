dnl Copyright (C) 2002,2004 Jens Franke, T.Kleinjung
dnl This file is part of gnfs4linux, distributed under the terms of the
dnl GNU General Public Licence and WITHOUT ANY WARRANTY.

dnl You should have received a copy of the GNU General Public License along
dnl with this program; see the file COPYING.  If not, write to the Free
dnl Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.

dnl translate_sieve_interval(array,length,table):
dnl replaces array[i] by table[array[i]] for 0<=i<length

define(array,%rdi)dnl
define(length,%rsi)dnl
define(table,%rdx)dnl
define(i,%rcx)dnl
define(tab,%r8)dnl
define(len8,%r9)dnl
define(aux0,%rcx)dnl
define(aux0b,%cl)dnl
define(aux1,%r11)dnl
define(aux2,%rax)dnl
define(b,%rdx)dnl

function_head(translate_sieve_interval1)
	testq length,length
	jz loopend
loop:
	movzbq (array),%rax
	movb (table,%rax,1),%cl
	movb %cl,(array)
	incq array
	decq length
	jnz loop
loopend:
	ret

function_head(translate_sieve_interval2)
	testq length,length
	jz loopend2
	leaq (array,length),array
	negq length
loop2:
	movzbq (array,length),%rax
	movb (table,%rax,1),%cl
	movb %cl,(array,length)
	incq length
	jnz loop2
loopend2:
	ret

function_head(translate_sieve_interval2l)
	xorl %ecx,%ecx
loop2l:
	mov %ecx,%r8d
	incl %ecx
	cmpl %esi,%ecx
	movzbq (%r8,array),%rax
	movzbl (%rax,table),%r9d
	movb %r9b,(%r8,array)
	jbe loop2l
loopend2l:
	ret

function_head(translate_sieve_interval3)
	pushq %rbx
	testq length,length
	jz loopend3
	movq table,%rbx
	leaq (array,length,1),array
	negq length
loop3:
	movb (array,length,1),%al
	xlatb
	movb %al,(array,length,1)
	incq length
	jnz loop3
loopend3:
	popq %rbx
	ret

function_head(translate_sieve_interval4)
	movq length,len8
	shrq $3,len8
	testq len8,len8
	movq table,tab
	jz loop1b
loop8:
	movq (array),aux0
	movq aux0,b
	andq $0xff,b
	shrq $8,aux0
	movzbq (tab,b),aux1
	movq aux0,b
	andq $0xff,b
	shrq $8,aux0
	movzbq (tab,b),aux2
	shlq $8,aux2
	addq aux2,aux1
	movq aux0,b
	andq $0xff,b
	shrq $8,aux0
	movzbq (tab,b),aux2
	shlq $16,aux2
	addq aux2,aux1
	movq aux0,b
	andq $0xff,b
	shrq $8,aux0
	movzbq (tab,b),aux2
	shlq $24,aux2
	addq aux2,aux1
	movq aux0,b
	andq $0xff,b
	shrq $8,aux0
	movzbq (tab,b),aux2
	shlq $32,aux2
	addq aux2,aux1
	movq aux0,b
	andq $0xff,b
	shrq $8,aux0
	movzbq (tab,b),aux2
	shlq $40,aux2
	addq aux2,aux1
	movq aux0,b
	andq $0xff,b
	shrq $8,aux0
	movzbq (tab,b),aux2
	shlq $48,aux2
	addq aux2,aux1
	movq aux0,b
	andq $0xff,b
	shrq $8,aux0
	movzbq (tab,b),aux2
	shlq $56,aux2
	addq aux2,aux1
	movq aux1,(array)
.if 0
	movzbq (array),%rax
	movb (tab,%rax),%cl
	movb %cl,(array)
	movzbq 1(array),%rax
	movb (tab,%rax),%dl
	movb %dl,1(array)
	movzbq 2(array),%rax
	movb (tab,%rax),%cl
	movb %cl,2(array)
	movzbq 3(array),%rax
	movb (tab,%rax),%dl
	movb %dl,3(array)
	movzbq 4(array),%rax
	movb (tab,%rax),%cl
	movb %cl,4(array)
	movzbq 5(array),%rax
	movb (tab,%rax),%dl
	movb %dl,5(array)
	movzbq 6(array),%rax
	movb (tab,%rax),%cl
	movb %cl,6(array)
	movzbq 7(array),%rax
	movb (tab,%rax),%dl
	movb %dl,7(array)
.endif
	decq len8
	leaq 8(array),array
	jnz loop8
	andq $0x7,length
loop1b:
	testq length,length
	jz loopend4
loop1:
	movzbq (array),%rax
	movb (tab,%rax),%cl
	movb %cl,(array)
	decq length
	leaq 8(array),array
	jnz loop1
loopend4:
	ret

