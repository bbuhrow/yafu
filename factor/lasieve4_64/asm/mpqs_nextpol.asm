dnl Copyright (C) 2002 Jens Franke, T.Kleinjung
dnl This file is part of gnfs4linux, distributed under the terms of the
dnl GNU General Public Licence and WITHOUT ANY WARRANTY.

dnl You should have received a copy of the GNU General Public License along
dnl with this program; see the file COPYING.  If not, write to the Free
dnl Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.


define(len,%rcx)dnl
define(disp_ptr,%r8)dnl
define(p_ptr,%r9)dnl
define(inv_ptr,%r10)dnl
define(start_ptr,%r11)dnl

dnl for (i=0; i<mpqs_nFB; i++) {
dnl   p=fb[2*i];
dnl   mmi=mpqs_FB_mm_inv[i];
dnl   pi=fbs[2*i]; bbb=fbs[2*i+1];
dnl   cc=bbb;
dnl   if (cc&1) cc+=p; cc>>=1;
dnl   cc=mpqs_FB_disp[i]+(p-cc); if (cc>=p) cc-=p;
dnl   cc1=fb[2*i+1];
dnl   h32=cc1*pi;
dnl   MMREDUCE; cc1=h32;
dnl   if (cc1>=p) cc1-=p;
dnl   cc2=p-cc1;
dnl   cc1+=cc; if (cc1>=p) cc1-=p;
dnl   cc2+=cc; if (cc2>=p) cc2-=p;
dnl   fbs[2*i]=(ushort)cc1; fbs[2*i+1]=(ushort)cc2;
dnl }

# asm_next_pol11_xmm(len)
function_head(asm_next_pol11_xmm)
	movq %rdi,len    # len
	movq $0x00010001,%rax
	movq %rax,%r9
	shlq $32,%rax
	orq %r9,%rax
	movd %rax,%xmm7
	pslldq $8,%xmm7
	movd %rax,%xmm6
	paddw %xmm6,%xmm7    # 0x00010001000100010001000100010001
	movq $mpqs_FB_disp,disp_ptr
	movq $mpqs_FB_np_px,p_ptr
	movq $mpqs_FB_mm_inv,inv_ptr
	movq $mpqs_FB_start,start_ptr
np11_mainloop:
	movaps (p_ptr),%xmm0     # p
	movaps 16(p_ptr),%xmm1    # sqrt
	movaps (inv_ptr),%xmm2     # mm_inv
	movaps (start_ptr),%xmm3     # pi
	movaps 16(start_ptr),%xmm4    # cc

	movaps %xmm3,%xmm6
	pmullw %xmm1,%xmm6     # low sqrt*cc1
	pmulhuw %xmm1,%xmm3    # high sqrt*cc1
	pmullw %xmm2,%xmm6     # h=low(sqrt*cc1)*mm_inv
	pxor %xmm5,%xmm5
	pcmpeqw %xmm6,%xmm5    # 0xffff iff h=0
	pand %xmm7,%xmm5       # 0x0001 iff h=0
	pmulhuw %xmm0,%xmm6    # high h*p
	paddw %xmm7,%xmm3
	psubw %xmm5,%xmm3      # carry
	paddw %xmm6,%xmm3      # res
	movaps %xmm3,%xmm6       # if >=p subtract p
	paddw %xmm7,%xmm6      # res+1
	pcmpgtw %xmm0,%xmm6    # 0xffff iff res>=p
	pand %xmm0,%xmm6       # p iff res>=p
	psubw %xmm6,%xmm3      # res mod p = cc1

	movaps %xmm4,%xmm5
	pand %xmm7,%xmm5       # 0x0001 iff cc odd
	pcmpeqw %xmm7,%xmm5    # 0xffff iff cc odd
	pand %xmm0,%xmm5       # p iff cc odd
	paddw %xmm4,%xmm5
	psrlw $1,%xmm5        # cc/2 mod p = cc

	movaps (disp_ptr),%xmm1     # disp
	movaps %xmm0,%xmm4
	psubw %xmm5,%xmm4      # p-cc
	paddw %xmm1,%xmm4
	movaps %xmm4,%xmm6       # if >=p subtract p
	paddw %xmm7,%xmm6      # res+1
	pcmpgtw %xmm0,%xmm6    # 0xffff iff res>=p
	pand %xmm0,%xmm6       # p iff res>=p
	psubw %xmm6,%xmm4      # res mod p = cc

	movaps %xmm0,%xmm2
	psubw %xmm3,%xmm2      # p-cc1 = cc2

	paddw %xmm4,%xmm3
	movaps %xmm3,%xmm6       # if >=p subtract p
	paddw %xmm7,%xmm6      # res+1
	pcmpgtw %xmm0,%xmm6    # 0xffff iff res>=p
	pand %xmm0,%xmm6       # p iff res>=p
	psubw %xmm6,%xmm3      # res mod p = cc1

	paddw %xmm4,%xmm2
	movaps %xmm2,%xmm5       # if >=p subtract p
	paddw %xmm7,%xmm5      # res+1
	pcmpgtw %xmm0,%xmm5    # 0xffff iff res>=p
	pand %xmm0,%xmm5       # p iff res>=p
	psubw %xmm5,%xmm2      # res mod p = cc2

	movaps %xmm3,%xmm4
	punpcklwd %xmm2,%xmm3
	punpckhwd %xmm2,%xmm4
	movaps %xmm3,(start_ptr)
	movaps %xmm4,16(start_ptr)

	decq len 
	leaq 32(start_ptr),start_ptr
	leaq 16(disp_ptr),disp_ptr
	leaq 32(p_ptr),p_ptr
	leaq 16(inv_ptr),inv_ptr
	jnz np11_mainloop

	emms
	ret

undefine(`len')dnl
undefine(`disp_ptr')dnl
undefine(`p_ptr')dnl
undefine(`inv_ptr')dnl
undefine(`start_ptr')dnl
define(len,%rcx)dnl
define(invtabptr,%rdx)dnl
define(ropptr,%r8)dnl
define(p_ptr,%r9)dnl
define(inv_ptr,%r10)dnl
define(start_ptr,%r11)dnl


dnl for (i=1; i<mpqs_nFB; i++) {
dnl   p=fb[2*i];
dnl   mmi=mpqs_FB_mm_inv[i];
dnl   cc=invptr[i];
dnl   cc*=bimul; h32=cc;
dnl   MMREDUCE; cc=h32; if (cc>=p) cc-=p;
dnl   ropptr[i]=(ushort)cc;
dnl   pi=fbs[2*i];
dnl   pi*=invptr[i]; h32=pi;
dnl   MMREDUCE; fbs[2*i]=(ushort)h32;
dnl   bbb=fbs[2*i+1]+cc; if (bbb>=p) bbb-=p; fbs[2*i+1]=bbb;
dnl }

# asm_next_pol10_xmm(len,*invptr,*ropptr,bimul)
function_head(asm_next_pol10_xmm)
	movq %rdx,ropptr

	movq $0x00010001,%rax
	movq %rax,%r9
	shlq $32,%rax
	orq %r9,%rax
	movd %rax,%xmm7
	pslldq $8,%xmm7
	movd %rax,%xmm6
	paddw %xmm6,%xmm7    # 0x00010001000100010001000100010001
	mulq %rcx            # bimul
	movd %rax,%xmm1
	pslldq $8,%xmm1
	movd %rax,%xmm5
	paddw %xmm5,%xmm1    # bimul

	movq %rdi,len    # len
	movq $mpqs_FB_np_px,p_ptr
	movq $mpqs_FB_mm_inv,inv_ptr
	movq $mpqs_FB_start,start_ptr
	movq %rsi,invtabptr    # invptr

np10_mainloop:
	movaps (p_ptr),%xmm0     # p
	movaps (inv_ptr),%xmm2     # mm_inv
	movaps (start_ptr),%xmm3     # pi
	movaps (invtabptr),%xmm4     # cc=invptr[i]

	movaps %xmm3,%xmm6
	pmullw %xmm4,%xmm6     # low pi*cc
	pmulhuw %xmm4,%xmm3    # high pi*cc
	pmullw %xmm2,%xmm6     # h=low(pi*cc)*mm_inv
	pxor %xmm5,%xmm5
	pcmpeqw %xmm6,%xmm5    # 0xffff iff h=0
	pand %xmm7,%xmm5       # 0x0001 iff h=0
	pmulhuw %xmm0,%xmm6    # high h*p
	paddw %xmm7,%xmm3
	psubw %xmm5,%xmm3      # carry
	paddw %xmm6,%xmm3      # res
	movaps %xmm3,(start_ptr)

	movaps %xmm4,%xmm6
	pmullw %xmm1,%xmm6     # low pi*cc
	pmulhuw %xmm1,%xmm4    # high pi*cc
	pmullw %xmm2,%xmm6     # h=low(pi*cc)*mm_inv
	pxor %xmm5,%xmm5
	pcmpeqw %xmm6,%xmm5    # 0xffff iff h=0
	pand %xmm7,%xmm5       # 0x0001 iff h=0
	pmulhuw %xmm0,%xmm6    # high h*p
	paddw %xmm7,%xmm4
	psubw %xmm5,%xmm4      # carry
	paddw %xmm6,%xmm4      # res
	movaps %xmm4,%xmm6       # if >=p subtract p
	paddw %xmm7,%xmm6      # res+1
	pcmpgtw %xmm0,%xmm6    # 0xffff iff res>=p
	pand %xmm0,%xmm6       # p iff res>=p
	psubw %xmm6,%xmm4      # res mod p
	movaps %xmm4,(ropptr)

	movaps 16(start_ptr),%xmm3    # fbs[2*i+1]
	paddw %xmm4,%xmm3      # fbs[2*i+1]+cc
	movaps %xmm3,%xmm6       # if >=p subtract p
	paddw %xmm7,%xmm6      # res+1
	pcmpgtw %xmm0,%xmm6    # 0xffff iff res>=p
	pand %xmm0,%xmm6       # p iff res>=p
	psubw %xmm6,%xmm3      # res mod p
	movaps %xmm3,16(start_ptr)

	decq len
	leaq 32(start_ptr),start_ptr
	leaq 16(ropptr),ropptr
	leaq 16(invtabptr),invtabptr
	leaq 32(p_ptr),p_ptr
	leaq 16(inv_ptr),inv_ptr
	jnz np10_mainloop

	emms
	ret

undefine(`len')dnl
undefine(`invtabptr')dnl
undefine(`ropptr')dnl
undefine(`p_ptr')dnl
undefine(`inv_ptr')dnl
undefine(`start_ptr')dnl


dnl asm_next_pol3plus(len,*SI_add)
function_head(asm_next_pol3plus)
	movq $mpqs_FB_start,%rcx
	movq $mpqs_FB_np_p,%rdx
	movq $0x00010001,%rax
	movd %rax,%mm7
	psllq $32,%mm7
	movd %rax,%mm6
	paddw %mm6,%mm7    # 0x0001000100010001

np3plus_mainloop:
	movq (%rdx),%mm2
	movq (%rsi),%mm4
	movq %mm2,%mm0
	movq %mm2,%mm1
	punpcklwd %mm2,%mm0   # p0,p0,p1,p1
	punpckhwd %mm2,%mm1   # p2,p2,p3,p3
	movq %mm4,%mm2
	punpcklwd %mm4,%mm2   # a0,a0,a1,a1
	movq %mm4,%mm3
	punpckhwd %mm4,%mm3   # a2,a2,a3,a3

	movq (%rcx),%mm4
	paddw %mm2,%mm4
	movq %mm7,%mm2
	paddw %mm4,%mm2
	pcmpgtw %mm0,%mm2
	pand %mm0,%mm2
	psubw %mm2,%mm4
	movq %mm4,(%rcx)

	movq 8(%rcx),%mm5
	paddw %mm3,%mm5
	movq %mm7,%mm3
	paddw %mm5,%mm3
	pcmpgtw %mm1,%mm3
	pand %mm1,%mm3
	psubw %mm3,%mm5
	movq %mm5,8(%rcx)

	leaq 16(%rdx),%rdx
	decq %rdi
	leaq 8(%rsi),%rsi
	leaq 16(%rcx),%rcx
	jnz np3plus_mainloop

	emms
	ret

dnl asm_next_pol3plus_xmm(len,*SI_add)
function_head(asm_next_pol3plus_xmm)
	movq $mpqs_FB_start,%rcx
	movq $mpqs_FB_np_px,%rdx
	movq $0x00010001,%rax
	movq %rax,%r9
	shlq $32,%rax
	orq %r9,%rax
	movd %rax,%xmm7
	pslldq $8,%xmm7
	movd %rax,%xmm6
	paddw %xmm6,%xmm7    # 0x00010001000100010001000100010001

np3plus_mainloopx:
	movaps (%rdx),%xmm2
	movaps (%rsi),%xmm4
	movaps %xmm2,%xmm0
	movaps %xmm2,%xmm1
	punpcklwd %xmm2,%xmm0   # p0,p0,p1,p1
	punpckhwd %xmm2,%xmm1   # p2,p2,p3,p3
	movaps %xmm4,%xmm2
	punpcklwd %xmm4,%xmm2   # a0,a0,a1,a1
	movaps %xmm4,%xmm3
	punpckhwd %xmm4,%xmm3   # a2,a2,a3,a3

	movaps (%rcx),%xmm4
	paddw %xmm2,%xmm4
	movaps %xmm7,%xmm2
	paddw %xmm4,%xmm2
	pcmpgtw %xmm0,%xmm2
	pand %xmm0,%xmm2
	psubw %xmm2,%xmm4
	movaps %xmm4,(%rcx)

	movaps 16(%rcx),%xmm5
	paddw %xmm3,%xmm5
	movaps %xmm7,%xmm3
	paddw %xmm5,%xmm3
	pcmpgtw %xmm1,%xmm3
	pand %xmm1,%xmm3
	psubw %xmm3,%xmm5
	movaps %xmm5,16(%rcx)

	leaq 32(%rdx),%rdx
	decq %rdi
	leaq 16(%rsi),%rsi
	leaq 32(%rcx),%rcx
	jnz np3plus_mainloopx

	emms
	ret


dnl asm_next_pol3minus(len,*SI_add)
function_head(asm_next_pol3minus)
	movq $mpqs_FB_start,%rcx
	movq $mpqs_FB_np_p,%rdx
	movq $0x00010001,%rax
	movd %rax,%mm7
	psllq $32,%mm7
	movd %rax,%mm6
	paddw %mm6,%mm7    # 0x0001000100010001

np3minus_mainloop:
	movq (%rdx),%mm2
	movq (%rsi),%mm4
	movq %mm2,%mm0
	movq %mm2,%mm1
	punpcklwd %mm2,%mm0   # p0,p0,p1,p1
	punpckhwd %mm2,%mm1   # p2,p2,p3,p3
	movq %mm4,%mm2
	punpcklwd %mm4,%mm2   # a0,a0,a1,a1
	movq %mm4,%mm3
	punpckhwd %mm4,%mm3   # a2,a2,a3,a3

	movq (%rcx),%mm4
	psubw %mm0,%mm2
	psubw %mm2,%mm4
	movq %mm7,%mm2
	paddw %mm4,%mm2
	pcmpgtw %mm0,%mm2
	pand %mm0,%mm2
	psubw %mm2,%mm4
	movq %mm4,(%rcx)

	movq 8(%rcx),%mm5
	psubw %mm1,%mm3
	psubw %mm3,%mm5
	movq %mm7,%mm3
	paddw %mm5,%mm3
	pcmpgtw %mm1,%mm3
	pand %mm1,%mm3
	psubw %mm3,%mm5
	movq %mm5,8(%rcx)

	leaq 16(%rdx),%rdx
	decq %rdi
	leaq 8(%rsi),%rsi
	leaq 16(%rcx),%rcx
	jnz np3minus_mainloop

	emms
	ret

dnl asm_next_pol3minus_xmm(len,*SI_add)
function_head(asm_next_pol3minus_xmm)
	movq $mpqs_FB_start,%rcx
	movq $mpqs_FB_np_px,%rdx
	movq $0x00010001,%rax
	movq %rax,%r9
	shlq $32,%rax
	orq %r9,%rax
	movd %rax,%xmm7
	pslldq $8,%xmm7
	movd %rax,%xmm6
	paddw %xmm6,%xmm7    # 0x00010001000100010001000100010001

np3minus_mainloopx:
	movaps (%rdx),%xmm2
	movaps (%rsi),%xmm4
	movaps %xmm2,%xmm0
	movaps %xmm2,%xmm1
	punpcklwd %xmm2,%xmm0   # p0,p0,p1,p1
	punpckhwd %xmm2,%xmm1   # p2,p2,p3,p3
	movaps %xmm4,%xmm2
	punpcklwd %xmm4,%xmm2   # a0,a0,a1,a1
	movaps %xmm4,%xmm3
	punpckhwd %xmm4,%xmm3   # a2,a2,a3,a3

	movaps (%rcx),%xmm4
	psubw %xmm0,%xmm2
	psubw %xmm2,%xmm4
	movaps %xmm7,%xmm2
	paddw %xmm4,%xmm2
	pcmpgtw %xmm0,%xmm2
	pand %xmm0,%xmm2
	psubw %xmm2,%xmm4
	movaps %xmm4,(%rcx)

	movaps 16(%rcx),%xmm5
	psubw %xmm1,%xmm3
	psubw %xmm3,%xmm5
	movaps %xmm7,%xmm3
	paddw %xmm5,%xmm3
	pcmpgtw %xmm1,%xmm3
	pand %xmm1,%xmm3
	psubw %xmm3,%xmm5
	movaps %xmm5,16(%rcx)

	leaq 32(%rdx),%rdx
	decq %rdi
	leaq 16(%rsi),%rsi
	leaq 32(%rcx),%rcx
	jnz np3minus_mainloopx

	emms
	ret

