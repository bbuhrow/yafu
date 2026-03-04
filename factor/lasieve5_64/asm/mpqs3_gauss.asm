# Copyright (C) 2002 Jens Franke, T.Kleinjung
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.


# //  asm_re_strip(rowptr,m,dptr,ucmptr):
# //
# //  tab[0]=0ULL;
# //  for (j=0,zz=1; j<4; j++,zz+=zz) {
# //    if (dptr[j]==-1) tab[zz]=0ULL;
# //    else tab[zz]=rowptr[dptr[j]];
# //    for (k=1; k<zz; k++) tab[zz+k]=tab[k]^tab[zz];
# //  }
# //  for (t=0; t<m; t++)
# //    rowptr[t]^=tab[ucmptr[t]];

define(`tab',%r8)dnl
define(`rowptr',%rdi)dnl
define(`cnt',%rsi)dnl
define(`dptr',%rdx)dnl
define(`ucmptr',%rcx)dnl
define(`d',%r9)dnl
define(`h',%r10)dnl
define(`h1',%r11)dnl
function_head(asm_re_strip)
	subq $140,%rsp

	leaq 12(%rsp),tab
	movq tab,%rax
	andq $7,%rax
	subq %rax,tab   # 8 | tab
	pxor %mm0,%mm0
# fill tab:
	movq %mm0,(tab)

	movswq (dptr),d
	leaq (rowptr,d,8),h
	addq $1,d
	cmovcq tab,h
	movq (h),%mm1   # row64[d[c0]+m*l] or 0
	movq %mm1,8(tab)

	movswq 2(dptr),d
	leaq (rowptr,d,8),h
	addq $1,d
	cmovcq tab,h
	movq (h),%mm2   # row64[d[c0+1]+m*l] or 0
	movq %mm2,16(tab)
	movq %mm2,%mm3
	pxor %mm1,%mm3
	movq %mm3,24(tab)

	movswq 4(dptr),d
	leaq (rowptr,d,8),h
	addq $1,d
	cmovcq tab,h
	movq (h),%mm4   # row64[d[c0+2]+m*l] or 0
	movq %mm4,32(tab)
	movq %mm4,%mm5
	movq %mm4,%mm6
	movq %mm4,%mm7
	pxor %mm1,%mm5
	pxor %mm2,%mm6
	pxor %mm3,%mm7
	movq %mm5,40(tab)
	movq %mm6,48(tab)
	movq %mm7,56(tab)

	movswq 6(dptr),d
	leaq (rowptr,d,8),h
	addq $1,d
	cmovcq tab,h
	testq $1,cnt    # flag for jump below
	movq (h),%mm0   # row64[d[c0+3]+m*l] or 0
	pxor %mm0,%mm1
	movq %mm0,64(tab)
	movq %mm1,72(tab)
	pxor %mm0,%mm2
	pxor %mm0,%mm3
	movq %mm2,80(tab)
	movq %mm3,88(tab)
	pxor %mm0,%mm4
	pxor %mm0,%mm5
	movq %mm4,96(tab)
	movq %mm5,104(tab)
	pxor %mm0,%mm6
	pxor %mm0,%mm7
	movq %mm6,112(tab)
	movq %mm7,120(tab)
# apply to strip
	jz odd_end
	movzbq (ucmptr),h
	movq (rowptr),%mm7
	addq $1,ucmptr
	pxor (tab,h,8),%mm7
	movq %mm7,(rowptr)
	leaq 8(rowptr),rowptr
odd_end:
	shrq $1,cnt
re_loop:
	decq cnt
	movzbq (ucmptr),h
	movzbq 1(ucmptr),h1
	movq (rowptr),%mm6
	movq 8(rowptr),%mm7
	pxor (tab,h,8),%mm6
	pxor (tab,h1,8),%mm7
	movq %mm6,(rowptr)
	movq %mm7,8(rowptr)
	leaq 2(ucmptr),ucmptr
	leaq 16(rowptr),rowptr
	jnz re_loop


	emms
	addq $140,%rsp
	ret

