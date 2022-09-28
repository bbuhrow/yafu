	.globl schedsieve
	.type schedsieve,@function
/*
Copyright (C) 2001 Jens Franke, T. Kleinjung.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.
*/

function_head(schedsieve)
	movl %edi,%eax
	subq $12,%rcx
	cmpq %rdx,%rcx
	movzwq (%rdx),%r8
	movzwq 4(%rdx),%r9
	movzwq 8(%rdx),%r10
	movzwq 12(%rdx),%r11
        leaq 16(%rdx),%rdx
        jbe fat_loop_end
fat_loop:
        prefetch 128(%rdx)
        addb %al,(%rsi,%r8)
        movzwq (%rdx),%r8
        addb %al,(%rsi,%r9)
        movzwq 4(%rdx),%r9
        addb %al,(%rsi,%r10)
        movzwq 8(%rdx),%r10
        addb %al,(%rsi,%r11)
        cmpq %rdx,%rcx
        movzwq 12(%rdx),%r11
        leaq 16(%rdx),%rdx
        ja fat_loop
fat_loop_end:
	addq $28,%rcx
	cmpq %rdx,%rcx
	leaq 4(%rdx),%rdx
	jbe schedsieve_end
        addb %al,(%rsi,%r8)
	cmpq %rdx,%rcx
	leaq 4(%rdx),%rdx
	jbe schedsieve_end
        addb %al,(%rsi,%r9)
	cmpq %rdx,%rcx
	jbe schedsieve_end
        addb %al,(%rsi,%r10)
schedsieve_end:
	ret
