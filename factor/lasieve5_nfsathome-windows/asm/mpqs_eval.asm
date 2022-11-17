dnl Copyright (C) 2002,2004 Jens Franke, T.Kleinjung
dnl This file is part of gnfs4linux, distributed under the terms of the
dnl GNU General Public Licence and WITHOUT ANY WARRANTY.

dnl You should have received a copy of the GNU General Public License along
dnl with this program; see the file COPYING.  If not, write to the Free
dnl Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.

dnl Written by T. Kleinjung
dnl Modified by J. Franke
dnl Athlon64 version by J. Franke

dnl asm_evaluate(sievebegin,sieveend,buffer,nmax):
dnl scans sievearray between sievebegin and sieveend for entries
dnl >127 and stores them in buffer (2 Bytes), stores at most nmax
dnl entries, returns number of stored entries

dnl edx: counts entries found so far
dnl esi: points to location of array we are investigating
dnl edi: points at end of array (=sieveend)
dnl mm7: 0
dnl mm0-3:
dnl ebx, ecx:

define(sv_arg,%rdi)dnl
define(svend,%rsi)dnl
define(buffer_arg,%rdx)dnl
define(nmax,%rcx)dnl
define(sv,%r8)dnl
define(svw,%r8w)dnl
define(buffer,%rax)dnl
define(buffer_ub,%r9)dnl

function_head(asm_evaluate0)
	movq sv_arg,sv
	movq buffer_arg,buffer
	pxor %mm7,%mm7
	leaq (buffer,nmax,2),buffer_ub
	jmp entry320

loop320:
	leaq 32(sv),sv
	movq %mm7,-32(sv)
	movq %mm7,-24(sv)
	movq %mm7,-16(sv)
	movq %mm7,-8(sv)

entry320:
	cmpq svend,sv
	jz end0

	movq (sv),%mm0
	movq 8(sv),%mm1
	movq 16(sv),%mm2
	movq 24(sv),%mm3

	por %mm0,%mm1
	por %mm2,%mm3
	por %mm1,%mm3
	pmovmskb %mm3,%r11
	testq %r11,%r11
	jz loop320

	movq 8(sv),%mm1
	movq 24(sv),%mm3
	pmovmskb %mm0,%r10
	pmovmskb %mm2,%r11
	salq $16,%r11
	pmovmskb %mm1,%rcx
	orq %r11,%r10
	pmovmskb %mm3,%r11
	salq $8,%rcx
	salq $24,%r11
	orq %rcx,%r10
	subq sv_arg,sv
	orq %r11,%r10
	xorq %r11,%r11
loop10:
	bsfq %r10,%rcx
	addq %rcx,%r11
	incq %rcx
	addq %r11,sv
	shrq %cl,%r10
	movw svw,(buffer)
	leaq 2(buffer),buffer 
	subq %r11,sv
	incq %r11
	cmpq buffer,buffer_ub
	jbe buffer_full0
	testq %r10,%r10
	jnz loop10
	addq sv_arg,sv
	jmp loop320

buffer_full0:
	addq sv_arg,sv
loop0320:
	cmpq svend,sv
	jz end0

	leaq 32(sv),sv
	movq %mm7,-32(sv)
	movq %mm7,-24(sv)
	movq %mm7,-16(sv)
	movq %mm7,-8(sv)
	jmp loop0320

end0:
	subq buffer_arg,buffer
	emms
	shrq $1,buffer
dnl buffer happens to be %rax
	ret


function_head(asm_evaluate)
	movq sv_arg,sv
	movq buffer_arg,buffer
	leaq (buffer,nmax,2),buffer_ub
	subq $32,sv

loop32:
	leaq 32(sv),sv
	cmpq svend,sv
	jz end

	movq (sv),%mm0
	movq 8(sv),%mm1
	movq 16(sv),%mm2
	movq 24(sv),%mm3

	por %mm0,%mm1
	por %mm2,%mm3
	por %mm1,%mm3
	pmovmskb %mm3,%r11
	testq %r11,%r11
	jz loop32

	movq 8(sv),%mm1
	movq 24(sv),%mm3
	pmovmskb %mm0,%r10
	pmovmskb %mm2,%r11
	salq $16,%r11
	pmovmskb %mm1,%rcx
	orq %r11,%r10
	pmovmskb %mm3,%r11
	salq $8,%rcx
	salq $24,%r11
	orq %rcx,%r10
	subq sv_arg,sv
	orq %r11,%r10
	xorq %r11,%r11
loop1:	
	bsfq %r10,%rcx
	addq %rcx,%r11
	incq %rcx
	addq %r11,sv
	shrq %cl,%r10
	movw svw,(buffer)
	leaq 2(buffer),buffer 
	subq %r11,sv
	incq %r11
	cmpq buffer,buffer_ub
	jbe end
	testq %r10,%r10
	jnz loop1
	addq sv_arg,sv
	jmp loop32

end:
	subq buffer_arg,buffer
	emms
	shrq $1,buffer
dnl buffer happens to be %rax
	ret


function_head(asm_evaluate0_xmm)
	movq sv_arg,sv
	movq buffer_arg,buffer
	xorps %xmm7,%xmm7
	leaq (buffer,nmax,2),buffer_ub
	jmp entry32x0

loop32x0:
	leaq 32(sv),sv
entry32x0:
	cmpq svend,sv
	jz endx0

	movaps (sv),%xmm0
	movaps %xmm7,(sv)
	movaps 16(sv),%xmm1
	movaps %xmm7,16(sv)
	movaps %xmm1,%xmm2
	orps %xmm0,%xmm1
	pmovmskb %xmm1,%r11
	testq %r11,%r11
	jz loop32x0

	pmovmskb %xmm0,%r10
	pmovmskb %xmm2,%r11
	shlq $16,%r11
	orq %r11,%r10
	subq sv_arg,sv
	xorq %r11,%r11
loop1x0:
	bsfq %r10,%rcx
	addq %rcx,%r11
	incq %rcx
	addq %r11,sv
	shrq %cl,%r10
	movw svw,(buffer)
	leaq 2(buffer),buffer 
	subq %r11,sv
	incq %r11
	cmpq buffer,buffer_ub
	jbe buffer_fullx0
	testq %r10,%r10
	jnz loop1x0
	addq sv_arg,sv
	jmp loop32x0

buffer_fullx0:
	addq sv_arg,sv
loop032x0:
	cmpq svend,sv
	jz endx0

	leaq 32(sv),sv
	movaps %xmm7,-32(sv)
	movaps %xmm7,-16(sv)
	jmp loop032x0

endx0:
	subq buffer_arg,buffer
	emms
	shrq $1,buffer
dnl buffer happens to be %rax
	ret


function_head(asm_evaluate_xmm)
	movq sv_arg,sv
	movq buffer_arg,buffer
	leaq (buffer,nmax,2),buffer_ub
	subq $32,sv
loop32x:
	leaq 32(sv),sv
	cmpq svend,sv
	jz endx

	movaps (sv),%xmm0
	movaps 16(sv),%xmm1
	movaps %xmm0,%xmm2
	orps %xmm1,%xmm0
	pmovmskb %xmm0,%r11
	testq %r11,%r11
	jz loop32x

	pmovmskb %xmm1,%r11
	pmovmskb %xmm2,%r10
	shlq $16,%r11
	orq %r11,%r10
	subq sv_arg,sv
	xorq %r11,%r11
loop1x:
	bsfq %r10,%rcx
	addq %rcx,%r11
	incq %rcx
	addq %r11,sv
	shrq %cl,%r10
	movw svw,(buffer)
	leaq 2(buffer),buffer 
	subq %r11,sv
	incq %r11
	cmpq buffer,buffer_ub
	jbe endx
	testq %r10,%r10
	jnz loop1x
	addq sv_arg,sv
	jmp loop32x

endx:
	subq buffer_arg,buffer
	emms
	shrq $1,buffer
dnl buffer happens to be %rax
	ret
