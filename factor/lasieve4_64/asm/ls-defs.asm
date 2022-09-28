define(l1_bits,15)dnl
define(n_i,eval(2**n_i_bits))dnl
define(n_i_mask,eval(n_i-1))dnl
define(l1_size,eval(2**l1_bits))dnl
define(l1_mask,eval(l1_size-1))dnl
define(j_per_strip,eval(2**(l1_bits-n_i_bits)))dnl
define(`forloop',
            `pushdef(`$1', `$2')_forloop(`$1', `$2', `$3', `$4')popdef(`$1')')
     define(`_forloop',
            `$4`'ifelse($1, `$3', ,
     		   `define(`$1', incr($1))_forloop(`$1', `$2', `$3', `$4')')')
define(function_head,.text
	.p2align 4`,,'15
.globl $1
	.type	$1`,' @function
$1:)dnl
