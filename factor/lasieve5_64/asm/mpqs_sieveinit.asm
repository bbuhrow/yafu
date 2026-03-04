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
undefine(`sv')dnl
undefine(`len')dnl
undefine(`tab_arg')dnl
undefine(`tab')dnl
undefine(`maskptr')dnl
undefine(`tinyptr')dnl
undefine(`tinylen')dnl
undefine(`tinylen4')dnl
undefine(`mask')dnl
undefine(`mask1')dnl
undefine(`ilen')dnl
undefine(`tptr')dnl
undefine(`tptrend')dnl
undefine(`mult')dnl
undefine(`llen')dnl
define(ta_arg,%rdi)dnl
define(mask16,%rsi)dnl
define(len,%rdx)dnl
define(ms,%mm0)dnl
define(mt,%mm1)dnl
define(ms0,%mm2)dnl
define(msz,%mm3)dnl
define(mtz,%mm4)dnl
define(ms0z,%mm5)dnl
define(src,%r8)dnl
define(targ,%r9)dnl
define(len8,%r10)dnl
define(sh,%mm6)dnl
define(ish,%mm7)dnl
function_head(asm_sieve_init0)
#  memcpy(mpqs3_tinyarray+mpqs3_tiny_prod,mpqs3_tinyarray,mpqs3_tiny_prod);

	movq ta_arg,src
	movq len,len8
	movq len8,%rax
	shrq $3,len8
	leaq (src,len8,8),targ
	andq $7,%rax     # shift in byte
	shlq $3,%rax     # shift in bit
	movd %rax,sh
	negq %rax
	addq $64,%rax
	movd %rax,ish

	movq (src),ms
	movq (targ),mt
	movq ms,ms0
	psllq ish,mt
	psrlq ish,mt
	psllq sh,ms
	pxor ms,mt
	movq mt,(targ)
	testq $1,len8
	jz mcloop1begin
	decq len8
	movq 8(src),ms
	movq ms0,mt
	addq $8,src
	addq $8,targ
	psrlq ish,mt
	movq ms,ms0
	psllq sh,ms
	pxor ms,mt
	movq mt,(targ)

mcloop1begin:
	testq len8,len8
	leaq 8(src),src
	leaq 8(targ),targ
	jz mcloop1end
	shrq $1,len8

mcloop1:
	movq (src),ms
	movq 8(src),msz
	movq ms0,mt
	movq ms,mtz
	psrlq ish,mt
	psrlq ish,mtz
	movq msz,ms0
	psllq sh,ms
	psllq sh,msz
	pxor ms,mt
	pxor msz,mtz
	movq mt,(targ)
	movq mtz,8(targ)
	decq len8
	leaq 16(src),src
	leaq 16(targ),targ
	jnz mcloop1

mcloop1end:
	psrlq ish,ms0
	movq ms0,(targ)

#  memcpy(mpqs3_tinyarray+2*mpqs3_tiny_prod,mpqs3_tinyarray,2*mpqs3_tiny_prod);

	movq ta_arg,src
	movq len,len8
	movq len8,%rax
	shrq $2,len8
	leaq (src,len8,8),targ
	andq $3,%rax     # shift in words
	shlq $4,%rax     # shift in bit
	movd %rax,sh
	negq %rax
	addq $64,%rax
	movd %rax,ish

	movq (src),ms
	movq (targ),mt
	movq ms,ms0
	psllq ish,mt
	psrlq ish,mt
	psllq sh,ms
	pxor ms,mt
	movq mt,(targ)
	testq $1,len8
	jz mcloop2begin
	decq len8
	movq 8(src),ms
	movq ms0,mt
	addq $8,src
	addq $8,targ
	psrlq ish,mt
	movq ms,ms0
	psllq sh,ms
	pxor ms,mt
	movq mt,(targ)

mcloop2begin:
	testq len8,len8
	leaq 8(src),src
	leaq 8(targ),targ
	jz mcloop2end
	shrq $1,len8

mcloop2:
	movq (src),ms
	movq 8(src),msz
	movq ms0,mt
	movq ms,mtz
	psrlq ish,mt
	psrlq ish,mtz
	movq msz,ms0
	psllq sh,ms
	psllq sh,msz
	pxor ms,mt
	pxor msz,mtz
	movq mt,(targ)
	movq mtz,8(targ)
	decq len8
	leaq 16(src),src
	leaq 16(targ),targ
	jnz mcloop2

mcloop2end:
	psrlq ish,ms0
	movq ms0,(targ)

#  memcpy(mpqs3_tinyarray+4*mpqs3_tiny_prod,mpqs3_tinyarray,4*mpqs3_tiny_prod);

	movq ta_arg,src
	movq len,len8
	movq len8,%rax
	shrq $1,len8
	leaq (src,len8,8),targ
	andq $1,%rax     # shift in ints
	shlq $5,%rax     # shift in bit
	movd %rax,sh
	negq %rax
	addq $64,%rax
	movd %rax,ish

	movq (src),ms
	movq (targ),mt
	movq ms,ms0
	psllq ish,mt
	psrlq ish,mt
	psllq sh,ms
	pxor ms,mt
	movq mt,(targ)
	testq $1,len8
	jz mcloop3begin
	decq len8
	movq 8(src),ms
	movq ms0,mt
	addq $8,src
	addq $8,targ
	psrlq ish,mt
	movq ms,ms0
	psllq sh,ms
	pxor ms,mt
	movq mt,(targ)

mcloop3begin:
	testq len8,len8
	leaq 8(src),src
	leaq 8(targ),targ
	jz mcloop3end
	shrq $1,len8

mcloop3:
	movq (src),ms
	movq 8(src),msz
	movq ms0,mt
	movq ms,mtz
	psrlq ish,mt
	psrlq ish,mtz
	movq msz,ms0
	psllq sh,ms
	psllq sh,msz
	pxor ms,mt
	pxor msz,mtz
	movq mt,(targ)
	movq mtz,8(targ)
	decq len8
	leaq 16(src),src
	leaq 16(targ),targ
	jnz mcloop3

mcloop3end:
	psrlq ish,ms0
	movq ms0,(targ)

#  memcpy(mpqs3_tinyarray+8*mpqs3_tiny_prod,mpqs3_tinyarray,8*mpqs3_tiny_prod);
#  memcpy(mpqs3_tinyarray+16*mpqs3_tiny_prod,mpqs3_tinyarray,16);
#  mask0=m64[0]; mask1=m64[1];
#  ullsv=(u64_t *)mpqs3_tinyarray;
#  ullsvend=ullsv+2*mpqs3_tiny_prod+2;
#  while (ullsv<ullsvend) {
#    *ullsv+++=mask0;
#    *ullsv+++=mask1;
#  }

	movq ta_arg,src
	movq len,len8
	leaq (src,len8,8),targ
	movq (mask16),mt
	movq 8(mask16),mtz
# len8 must be odd
	movq (src),ms
	movq ms,ms0
	paddb mt,ms
	paddb mtz,ms0
	movq ms,(src)
	movq ms0,(targ)

	shrq $1,len8
	testq len8,len8
	leaq 8(src),src
	leaq 8(targ),targ
	jz mcloop4end

mcloop4:
	movq (src),ms
	movq 8(src),msz
	movq ms,ms0
	movq msz,ms0z
	paddb mtz,ms
	paddb mt,msz
	paddb mt,ms0
	paddb mtz,ms0z
	movq ms,(src)
	movq msz,8(src)
	movq ms0,(targ)
	movq ms0z,8(targ)
	decq len8
	leaq 16(src),src
	leaq 16(targ),targ
	jnz mcloop4

mcloop4end:
	movq ta_arg,src
	movq (src),ms
	movq 8(src),msz
	movq ms,(targ)
	movq msz,8(targ)

	emms
	ret


