dnl smjs define one of the following (windows, osx, linux)
define(windows,1)dnl
dnl define(osx,1)dnl
dnl define(linux,1)dnl

define(`error_exit',
       `errprint(__file__:__line__`: ERROR: $*
')m4exit(`1')')dnl

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

ifelse(windows,`1', 
`define(function_head,.text
	.p2align 4`,,'15
.globl $1
        .def    $1;     .scl    2;      .type   32;     .endef
$1:)',dnl
osx,`1', 
`define(function_head,.text
	.p2align 4`,,'15
.globl _$1
_$1:)',dnl
linux,`1', 
`define(function_head,.text
	.p2align 4`,,'15
.globl $1
	.type	$1`,' @function
$1:)',dnl
`error_exit(`  You need to edit `asm/ls-defs.asm' to define which OS you are running on (linux, osx or windows)')'
)
