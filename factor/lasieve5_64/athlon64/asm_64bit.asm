dnl Copyright (C) 2006 Jens Franke, T.Kleinjung
dnl This file is part of gnfs4linux, distributed under the terms of the
dnl GNU General Public Licence and WITHOUT ANY WARRANTY.

dnl You should have received a copy of the GNU General Public License along
dnl with this program; see the file COPYING.  If not, write to the Free
dnl Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.

#include "underscore.h"

function_head(asm_modadd64)
dnl smjs	movq modulo64(%rip),%rdx
	movq modulo64(%rip),%rdx
	movq %rdi,%rax
	subq %rdx,%rax
	addq %rsi,%rax
	movq $0,%rcx
	cmovcq %rcx,%rdx
	addq %rdx,%rax
	ret

function_head(asm_modmul64)
	movq %rdi,%rax
	mulq %rsi
dnl smjs	movq modulo64(%rip),%rcx
	movq modulo64(%rip),%rcx
	divq %rcx
	movq %rdx,%rax
	ret
