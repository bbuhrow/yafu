\def\bfZ{{\bf{Z}}}
\def\abs#1{{|#1|}}
@f u32_t int
@f u16_t int
@* Obtaining recurrence information for sieving event.

@*3 Copying.
Copyright (C) 2001 Jens Franke.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.

@*3 Mathematical introduction to the problem.

In order to efficiently implement lattice sieving for the number field sieve,
one has to solve the following problem. Assume that sieving for a fixed
special $q$ has been prepared, such that the lattice of $(a,b)$-pairs
which are divisible by this special $q$ can be parametrized by pairs
$(i,j)$. Let the admissible range for $i$ be denoted by $I'$, which is
assumed to be an interval of length $A$. Assume
moreover that we want to sieve with a prime ideal lying above $p$ whose
intersection with the $(i,j)$-lattice is given by the congruence
$i\equiv rj\bmod p$. If $p<A$, sieving is easy to organize: Every $i$-line
contains at least one sieving event. Otherwise, every $i$-line contains at
most one sieving event, and things are more complicated. In order to avoid
waste of computer time, it is necessary to predict, given an $i$-line and
a sieving event located on the intersection of $I'$ with this $i$-line,
the next $i$-line (belonging to the smallest $j$ larger than the $j$ for
the given $i$-line) which also contains a sieving event inside $I'$.

The problem is more conveniently formulated in terms of $I$, the set of
non-negative integers less than $A$. We view $I$ as a subset of $\bfZ/p\bfZ$
and consider that dynamical system $T(x)=x+r$ on $\bfZ/p\bfZ$. It is
then necessary to determine, for each element of $I$, the time of its
return to $I$. It is easy to see that if $x\in I$ returns to $I$ at
time $l$, then one of the endpoints $0$ or $A-1$ of $I$ must also returns to
$I$ at time $l$.

@ The point $0$ recurs to $I$ at time $k$ if and only if there exists
an $l$ such that
$$ 0\le rl-kp<A.\eqno(1)$$
Similarly, the point $A-1$ returns to $I$ at time $k$ if and only if theres
exists $l$ such that
$$ 0\le pk-rl<A.\eqno(2)$$
We have $\abs{lx-k}<A/p$ with $x=r/p$ if and only if one of these inequalities
holds. The $(k,l)$-pair with the smallest $l$ such that (1) or (2) holds is a
Diophantine approximation to $r/p$, hence is given by an intermediate result
$k/l$ of the continued fraction expansion $[0,a_1,\ldots,a_n]$ of $r/p$.
Therefore, it easy to use the continued fraction expansion to $k/l$ to
determine the earliest time $l>0$ at which some point (which can be chosen
to be an endpoint) of $I$ returns to $I$. The next intermediate fraction of
the continued fraction expansion then gives a time at which the other endpoint
returns to the interval. Unfortunately, that this is the smallest return time
can be guaranteed only if $k/l=[0,\ldots,a_j]$ with $a_{j+1}=1$. The general
case is described by the following assertion.

@ {\bf Proposition.}
Let $k/l=[0,a_1,\ldots,a_n]$, $k_-/l_-=[0,a_1,\ldots,a_{n-1}]$
and $k_+/l_+=[0,a_1,\ldots,a_{n+1}]$ be three successive continued fraction
approximations to the real number $x$ satisfying $0<x<1$. Then all pairs
$\tilde k,\tilde l$ of coprime positive integers with
$\abs{\tilde k-x\tilde l}<\abs{k_-x-l_-}$ and $l<\tilde l\le l_+$
% and for which
% $\tilde k-x\tilde l$ and $k_--xl_-$ have the same sign
are given by $\tilde k/\tilde l=[0,a_1,\ldots,a_n,a]$ with $1\le a\le a_{n+1}$.

@ The significance of this fact for the problem posed at the beginning is
obvious. In order to see this, apply the proposition to $x=r/p$ where
$l$ is the earliest time at which some endpoint of $I$ returns to $I$.
We know that the other endpoint returns to $I$ at time $l_+$, for
$k_+/l_+$ is a Diophantine approximation to $x$ from the opposite side and
of a better quality than $k/l$. If the other endpoint of $I$ does not yet
return to $I$ at time $l$, then its return time must be given by an
approximation $\tilde k/\tilde l$ to $r$ with $l<\tilde l\le l_+$. Since
$\abs{k_--xl_-}$ is too large to allow the return of any point of $I$ to
$I$ at time $l_-$, we must have $\abs{\tilde k-x\tilde l}<\abs{k_--xl_-}$ and
the proposition can be applied. Let $n$ have the same meaning as in the
proposition: $k/l=[0,a_1,\dots,a_n]$. Let $r/p=[0,\dots,a_N]$.

It is easy to see that the return time of $y\in I$ coincides
with the earliest time at which $0$ returns to $I'=\{0;\ldots;A-1-y\}$
or $y$ returns to $\{0;\ldots;y\}$. Therefore, it must be the denominator
of a continued fraction $[0,a_0,\dots,a_m,a]$ with $n\le m<N$ and
$0<a\le a_{m+1}$. Therefore, the proposition contains
the complete solution to our problem. 
In addition, let $(\hat k,\hat l)=(k,l)+(\tilde k,\tilde l)$. If
$\tilde k/tilde l=[0,a_0,\dots,a_n,a]$ then
$\hat k/\hat l=[0,a_0,\dots,a_n,a+1]$, and in the case $a=a_{n+1}$ this
also coincides with $[0,a_0,\dots,a_{n+1},1]$. This means that $\hat l$
is the smallest time $>\tilde l$ which is a possible candidate for the
return time of $y$. Since $\hat d\abs{\hat kp-\hat lr}$ is equal to the
absolute value of $\tilde d-d$, where $\tilde d=\abs{\tilde kp-\tilde lr}$
and $d=\abs{kp-lr}$, we have $\hat d+\tilde d\le d\le A$ and also
$\hat d+d\le\tilde d\le A$. Since $\tilde kp-\tilde lr$ and $kp-lr$ have
opposite signs, this is easily seen to imply that $y$ returns to $I$ at
time $\hat l$ if it does not already return at time $l$ or $\tilde l$. This
means that there are only three continued fractions which are relevant to our
problem.

@ {\bf Proof of the Proposition.}
We make use of the following well-known facts about continued fractions:
\itemitem{(a)} We have $kl_--k_-l=\pm1$
\itemitem{(b)} The numbers $k-lx$ and $k_--l_-x$ have opposite signs.
\itemitem{(c)} We have $l_+=a_{n+1}l+l_-$ and $k_+=a_{n+1}k+k_-$.

By the first of these facts, it is possible to find
integers $a$ and $b$ solving the two equations
$$ \tilde k=ak+bk_-\qquad\qquad\tilde l=al+bl_-.$$
Since $\tilde k/\tilde l$ must be different from $k/l$ and $k_-/l_-$, we
have $a\not=0$ and $b\not=0$. By the second of these facts, $a$ and $b$ have
the same sign, or $\abs{\tilde k-x\tilde l}\ge\abs{k-xl}+\abs{k_--xl_-}$ in
contradiction to our assumptions. Since $\tilde l$ is positive, $a$ and
$b$ are both positive. If $a>a_{n+1}$, then $\tilde l>l_+=a_{n+1}l+l_-$ by
fact (c). If $b>1$ and $a\le a_{n+1}$, then
$$\abs{\tilde k-x\tilde l}=\abs{(b-1)(k_--xl_-)+(k_+-xl_+)+(a-a_{n+1})(k-xl)}
\ge\abs{k_--xl_-},$$
where the fact (b) was used two more times, one time with $n$ replaced
by $n+1$. We arrive at $b=1$, hence $\tilde k/\tilde l=[0,a_1,\dots,a_n,a]$
with $1\le a\le a_{n+1}$ as claimed.

@ Our function takes the form
@(recurrence6.h@>=
#include "asm/siever-config.h"
void rec_info_init(u32_t A,u32_t ub);
u32_t ASM_ATTR get_recurrence_info(u32_t *res_ptr,u32_t p,u32_t r);
// SMJS How to remove underscore issue - u32_t 
//      get_recurrence_info(u32_t *res_ptr,u32_t p,u32_t r) asm ("get_recurrence_info");

@ The first argument to |rec_info_init| has the same meaning as above, that
is, it specifies the length of the interval $I'$ of admissible values for
the $i$-coordinate of the sieving lattice. The second argument |ub| does
almost the same thing for |j|. In other words, it gives the strict upper bound
on the recurrence times we are interested in. The meaning of the
second to fourth arguments will be explained below.

The fifth and sixth argument to |get_recurrence_infoA| and
|get_recurrence_infoB| have the meaning described above. The first argument to
this function is a pointer to an array
of |short int| integers to which the recurrence information is written.
The meaning of this information is as follows: Let
$x=$(|res_ptr[0],res_ptr[1]|) and $y=$(|res_ptr[2],res_ptr[3]|). Then
for arbitray $i\in I$, the first recurrence of $i$ to $I$ and its time
are given by $(i,0)+x$, $(i,0)+y$ or $(i,0)+x+y$, depending on the cases
$i<a$, $i\ge b$ or $a\le i<b$. If |get_recurrence_infoB| is used, then
$a$ is |res_ptr[4]| and $b$ is |res_ptr[5]|. Otherwise, they are given
by $a=$|-res_ptr[0]| and $b=$|A-res_ptr[3]|.

The functions |get_recurrence_infoC| and |get_recurrence_intoD| are
similar. However, instead of storing |x| at |res_ptr|, the number
|*((u32_t*)res_ptr)| will equal $x_1*A+x_0$. Similarly, |*(1+(u32_t*)res_ptr)|
will equal $y_1*A+x_1$. The function |get_recurrence_infoC| writes the shorter
format, |get_recurrence_infoD| also stores $a$ and $b$ after $x$ and $y$.

It is a fatal error to call |get_recurrence_info| with
|r>p| or with |p<=ub| or |p<=A| or |r!=0 && r!=p && gcd(p,r)!=1| or
with odd |A|.

All functions returns the number of |short int|-items written
to |res_ptr|.

@* Implementation of the function |get_recurrence_info|.

It is necessary to store auxilliary information in the following
variables:
@c
#include <sys/types.h>
#include <limits.h>
#include "asm/siever-config.h"
#include "if.h"
#include "recurrence6.h"

static u32_t A,A_bits,ub;

void rec_info_init(u32_t A1,u32_t ub1)
{
  u32_t aa;
  if(A1>USHRT_MAX) complain("Recurrence init: A=%u exceeds %u\n",A1,USHRT_MAX);
  A=A1;
  if(A%2 != 0) complain("rec_info_init with odd range %u\n",A1);
  if(ub1>USHRT_MAX/2+2) /* changed for I16!!!, is this ok??? */
    complain("Recurrence init: ub=%u exceeds %u\n",ub1,USHRT_MAX/4);
  if(ub1<2)
    complain("Recurrence time %u too small\n",ub1);
  ub=ub1;
  for(aa=1,A_bits=0;A_bits<=CHAR_BIT*sizeof(u32_t);A_bits++,aa*=2)
    if(aa==A) break;
  if(A_bits>CHAR_BIT*sizeof(u32_t))
    complain("rec_info_init: A=%u not a power of two\n",A);
}

@
@c
u32_t ASM_ATTR
get_recurrence_info(u32_t *res_ptr,u32_t p,u32_t r)
{
  u32_t b,c,s,t;

  if(r==0) {
    @<Treat root zero.@>@;
    goto done;
  }
  if(r>=p) {
    @<Treat root infinity.@>@;
    goto done;
  }

  @<Determine |b|, |c|, |s|, |t|@>@;
  @<Modify them to satisfy our size restrictions@>@;
 done: /* I promise that no code outside this module jumps to this point. */
  @<Write the recurrence info@>@;
#if 0
  @<Write the first sieving events@>@;
#endif
  return 2;
}

@ In principle, we want the lattice $\Gamma'$ generated by (|b|,|s|) and
(-|c|,|t|) be be equal to the lattice $\Gamma$ of all $(i,j)$ with
$i\cong rj\pmod p$. However, the values for |b| and |c| must fit into a
|short int| and the values for |s| and |t| must fit into a
|short unsigned int|. It is easy to see that the following code achieves
this, but we only have $\Gamma\cap\Omega=\Gamma'\cap\Omega|, where
Omega consists of all coprime $(i,j)$ with $-A\le i<A$ and $b<=2*ub$.

@<Treat root zero.@>=
b=1;
s=2*ub+1;
c=A-3;
t=2*ub;  

@ In this case, $\Gamma\cap\Omega$ contains at most the points $(p,k)$,
$-A\le k<A$, $k\not=0$.
@<Treat root infinity.@>=
b=1;
s=0;
c=A-1;
t=2*ub+1;

@ These numbers have the following meaning: $s$ is the smallest
positive integer with $b=s*r\bmod p<A$. Similarly, $t$ is the smallest positive
integer with $c=-t*r\bmod p<A$.

This essentially is a truncated semi-extended Euclidean Algorithm.

@<Determine |b|, |c|, |s|, |t|@>=
#ifdef HAVE_ASM_GETBC
asm_getbc(r,p,A,&b,&s,&c,&t);
#else
{
  b=r;
  s=1;
  c=p;
  t=0;
  for(;;) {
    u32_t k;

    if(b<A) {
    /* Get the smallest |k| with |c-k*b<A|. This is: */
      k=(c-A)/b+1;
      t+=k*s;
      c-=k*b;
      break;
    }
    k=c/b;
    c=c%b;
    t+=k*s;
    if(c<A) {
      /* Now, calculate the smallest |k| satisfying |b-k*c<A| */
      k=(b-A)/c+1;
      s+=k*t;
      b-=k*c;
      break;
    }
    k=b/c;
    b=b%c;
    s+=k*t;
  }
}
#endif

@ You may want to use this if you want to test an asm version of the
above.
@<Determine |b|, |c|, |s|, |t|@>=
#if 0
{
  u32_t b1,c1,s1,t1;

  asm_getbc(r,p,A,&b1,&s1,&c1,&t1);
  if(b1 != b || c1 != c || s1 != s || t1 != t) Schlendrian("Asm ri\n");
}
#endif

@ The same remark as above applies: At the moment, we have a basis
$(b,s)$ and $(-c,t)$ of our sieving lattice $Gamma$, which possibly violates
our size restrictions. We will obtain a quadruple satisfying our size
restrictions and generating a sublattice $\Gamma'$ with
$\Gamma\cap\Omega=\Gamma'\cat\Omega$. This modified quadruple
will still be linearly independent over $\mathop{GF}\nolimits_2$.

In most cases, the elements of $\Omega\cap\Gamma$ will be linear combinations
of $x=(b,s)$ and $y=(-c,t)$ with non-negative coefficients. There is, however,
an exceptional case which must not be ignored because it would cause serious
trouble once in a few million factor base elements. If we denote by
$\tilde\Omega$ the set of all $(i,j)\in\Omega$ with $i>-a$ (instead of just
$i\ge-a$), then it is easy to see that the elements $\tilde\Omega\cap\Gamma$
are non-negative linear combinations of $x$ and $y$. In the limiting
case $A=b+c$ and if $t>s$, then we also have $y-x=(-A,t-s)\in\Omega\cap\Gamma$.

@<Modify them to satisfy our size restrictions@>=
{
  i32_t d;

  d=t-s;
  s=( s<=2*ub ? s : 2*ub+2-(s&1));
  if(t>2*ub) {
    if(d<0 || d>2*ub) d=2*ub+2+(s&1)-(t&1);
    t=s+d;
  }
}

@
@<Write the recurrence info@>=
{
  res_ptr[0]=(s<<A_bits)+b;
  res_ptr[1]=(t<<A_bits)-c;
}

@ Recall that $x=(b,s)$ and $y=(-c,t)$. Define the oddness type of a pair of
integers as |(z[0]%2)+2*(z[1]%2)|. We have to calculate the first sieve
events (i.e., with minimal |z[1]|) for each of the three oddness types.

@<Write the first sieving events@>=
{
  u32_t z[2];
  u32_t have_it[3]={0,0,0};

  z[0]=A+b;
  z[1]=s;
  @<Store |z|@>@;
  if(b+c<=A && s<=t) {
    if(b+c<A) {
      z[0]=A;
      z[1]=1;
    } else {
      z[0]=0;
      z[1]=t-s;
    }
  } else {
    z[0]=A+b-c;
    z[1]=t+s;
  }
#if 0
  z[1]=( z[1]<=2*ub ? z[1] : 2*ub+2-(z[1]&1) );
#else
  z[1]=( z[1]<=USHRT_MAX ? z[1] : 2*ub+2-(z[1]&1) );
#endif
  @<Store |z|@>@;
  z[0]=A-c;
  z[1]=t;
  @<Store |z|@>@;
  if(have_it[0]==0) Schlendrian("???");
}

@
@<Store |z|@>=
{
  u32_t ot; /* Oddness type. */
#ifdef ULONG_RI
  u32_t *xx;
#else
  u16_t *xx;
#endif

  ot=(z[0]&1)+2*(z[1]&1);
  if(ot==1) xx=x1;
  else if(ot==2) xx=x2;
  else xx=x3;
  have_it[ot-1]=1;
#ifdef ULONG_RI
  *xx=((z[1]/2)<<A_bits)|(z[0]/2);
#else
  xx[0]=z[0]/2;
  xx[1]=z[1]/2;
#endif
}
