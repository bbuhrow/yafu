/*
kmpqs: Kleinjung's small-number MPQS from lasieve5 (gnfs4linux), ported
into yafu.

The algorithms are Copyright (C) 2001 Jens Franke, T. Kleinjung, and this
module is distributed under the terms of the GNU General Public Licence and
WITHOUT ANY WARRANTY.  See factor/lasieve5_64/COPYING.

Kernels shared by kmpqs1.c and kmpqs3.c, replacing lasieve's asm.  Each
kernel computes exactly what lasieve's C code computes (same sieve bytes,
same survivors in the same order, same prime lists), so the relation
counters do not change.  The ISA is chosen at compile time from the
compiler's target macros: AVX-512BW, AVX2, SSE4.1/SSE2, else scalar.

Data layout (both variants):
  fb[2*j], fb[2*j+1]      prime p_j and sqrt(kN) mod p_j (uint16)
  fbs[2*j], fbs[2*j+1]    the two sieve roots of p_j, each < p_j
  p2[2*j] = p2[2*j+1]     p_j, duplicated to line up with fbs
  inv2[2*j] = inv2[2*j+1] p_j^-1 mod 2^16, duplicated
Prime index 0 is 2 and is never touched here.
*/

#ifndef KMPQS_KERNELS_H
#define KMPQS_KERNELS_H

#include <stdint.h>
#include <string.h>

/* -DKMPQS_NO_SIMD forces the scalar kernels (for testing, or other CPUs). */
#if !defined(KMPQS_NO_SIMD)
#if defined(__AVX512BW__) && defined(__AVX512F__) && defined(__AVX512VL__)
#define KMPQS_AVX512 1
#endif
#if defined(__AVX2__)
#define KMPQS_AVX2 1
#endif
#if defined(__SSE4_1__)
#define KMPQS_SSE41 1
#endif
#if defined(__SSE2__) || defined(_M_X64)
#define KMPQS_SSE2 1
#endif
#endif

#if defined(KMPQS_AVX512) || defined(KMPQS_AVX2) || defined(KMPQS_SSE2)
#include <immintrin.h>
#endif

#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>
static __inline int kmpqs_ctz32(uint32_t x) { unsigned long r; _BitScanForward(&r, x); return (int)r; }
static __inline int kmpqs_ctz64(uint64_t x) { unsigned long r; _BitScanForward64(&r, x); return (int)r; }
#else
#define kmpqs_ctz32(x) __builtin_ctz(x)
#define kmpqs_ctz64(x) __builtin_ctzll(x)
#endif

/* Padding (in uint16 elements) after arrays that vector loops read or
   write in whole vectors. */
#define KMPQS_VPAD 32

/* Jacobi symbol (a/n) for odd n and 0 <= a < n (binary algorithm); 0 if
   gcd(a, n) > 1.  Replaces lasieve's subtractive version, which loops
   forever on a = 0. */
static inline int kmpqs_jacobi(uint32_t a, uint32_t n)
{
	uint32_t t = 0, tmp;

	while (a) {
		int z = kmpqs_ctz32(a);

		a >>= z;
		if (z & 1)			/* (2/n) = -1 for n = 3, 5 mod 8 */
			t ^= ((n >> 1) ^ (n >> 2)) & 1;
		if (a < n) {			/* reciprocity: flip if both = 3 mod 4 */
			t ^= (a & n & 2) >> 1;
			tmp = a; a = n; n = tmp;
		}
		a -= n;
	}
	return (n == 1) ? (t ? -1 : 1) : 0;
}

/* x mod p, for x given as nw 32-bit words (least significant first) and
   p < 2^31 (so (r << 32) | w fits in 64 bits). */
static inline uint32_t kmpqs_mod_words(const uint32_t *w, int nw, uint32_t p)
{
	uint64_t r = 0;
	int i;

	for (i = nw - 1; i >= 0; i--)
		r = ((r << 32) | w[i]) % p;
	return (uint32_t)r;
}

/* Name of the kernel set this file was compiled with. */
#if defined(KMPQS_AVX512)
#define KMPQS_KERNEL_ISA "avx512bw"
#elif defined(KMPQS_AVX2)
#define KMPQS_KERNEL_ISA "avx2"
#elif defined(KMPQS_SSE2)
#define KMPQS_KERNEL_ISA "sse2"
#else
#define KMPQS_KERNEL_ISA "scalar"
#endif

/* ---------------------------------------------------------------------------
   Sieve: add log p at every x in [0,L) with x == s1 or s2 (mod p), for the
   primes fb[2*i], i >= i0, until the 0xffff terminator.  Primes ascend and
   s1, s2 < p.  sv[L] must exist; it absorbs the branch-free boundary hits
   (as lasieve's asm does), so its value is garbage afterwards.
   --------------------------------------------------------------------------- */
/* Pointers (not 32-bit indices) so the adds use base+stride addressing.
   Positions are compared as integers: a+3p may lie past the array. */
#define KMPQS_CLAMP(q, end) ((uintptr_t)(q) < (uintptr_t)(end) ? (q) : (end))

static inline void kmpqs_sieve_fb(uint8_t *sv, uint32_t L, const uint16_t *fb,
	const uint16_t *fbs, const uint8_t *fblog, uint32_t i)
{
	uint8_t *const end = sv + L;
	uint8_t *a, *b;
	size_t p, r1, r2;
	uint8_t lo;

	/* p <= L/4: unrolled by 4 per root */
	for (;; i++) {
		size_t p2, p3, p4;
		uint8_t *bound;

		p = fb[2 * i];
		if (p > (L >> 2)) break;
		lo = fblog[i];
		r1 = fbs[2 * i]; r2 = fbs[2 * i + 1];
		a = sv + (r1 < r2 ? r1 : r2);
		b = sv + (r1 < r2 ? r2 : r1);
		p2 = 2 * p; p3 = 3 * p; p4 = 4 * p;
		bound = end - p3;
		while (b < bound) {
			a[0] += lo; a[p] += lo; a[p2] += lo; a[p3] += lo;
			b[0] += lo; b[p] += lo; b[p2] += lo; b[p3] += lo;
			a += p4; b += p4;
		}
		/* b+3p >= end: b has at most 3 hits left, a (< b) at most 4 */
		*KMPQS_CLAMP(a, end) += lo;
		*KMPQS_CLAMP(a + p, end) += lo;
		*KMPQS_CLAMP(a + p2, end) += lo;
		*KMPQS_CLAMP(a + p3, end) += lo;
		*KMPQS_CLAMP(b, end) += lo;
		*KMPQS_CLAMP(b + p, end) += lo;
		*KMPQS_CLAMP(b + p2, end) += lo;
	}
	/* L/4 < p <= L/3: 3 hits per root, a 4th maybe */
	for (;; i++) {
		p = fb[2 * i];
		if (p > L / 3) break;
		lo = fblog[i];
		a = sv + fbs[2 * i]; b = sv + fbs[2 * i + 1];
		a[0] += lo; a[p] += lo; a[2 * p] += lo;
		b[0] += lo; b[p] += lo; b[2 * p] += lo;
		*KMPQS_CLAMP(a + 3 * p, end) += lo;
		*KMPQS_CLAMP(b + 3 * p, end) += lo;
	}
	/* L/3 < p <= L/2: 2 hits per root, a 3rd maybe */
	for (;; i++) {
		p = fb[2 * i];
		if (p > (L >> 1)) break;
		lo = fblog[i];
		a = sv + fbs[2 * i]; b = sv + fbs[2 * i + 1];
		a[0] += lo; a[p] += lo;
		b[0] += lo; b[p] += lo;
		*KMPQS_CLAMP(a + 2 * p, end) += lo;
		*KMPQS_CLAMP(b + 2 * p, end) += lo;
	}
	/* L/2 < p <= L: 1 hit per root, a 2nd maybe */
	for (;; i++) {
		p = fb[2 * i];
		if (p > L) break;
		lo = fblog[i];
		a = sv + fbs[2 * i]; b = sv + fbs[2 * i + 1];
		a[0] += lo; b[0] += lo;
		*KMPQS_CLAMP(a + p, end) += lo;
		*KMPQS_CLAMP(b + p, end) += lo;
	}
	/* p > L: at most 1 hit per root */
	for (;; i++) {
		p = fb[2 * i];
		if (p == 0xffff) break;
		lo = fblog[i];
		*KMPQS_CLAMP(sv + fbs[2 * i], end) += lo;
		*KMPQS_CLAMP(sv + fbs[2 * i + 1], end) += lo;
	}
}

/* Same for one root per prime (the divisors of the multiplier k). */
static inline void kmpqs_sieve_one(uint8_t *sv, uint32_t L, size_t p,
	uint32_t r, uint8_t lo)
{
	uint8_t *const end = sv + L;
	uint8_t *a = sv + r;

	while (a + 3 * p < end) {
		a[0] += lo; a[p] += lo; a[2 * p] += lo; a[3 * p] += lo;
		a += 4 * p;
	}
	*KMPQS_CLAMP(a, end) += lo;
	*KMPQS_CLAMP(a + p, end) += lo;
	*KMPQS_CLAMP(a + 2 * p, end) += lo;
}

/* ---------------------------------------------------------------------------
   Sieve initialisation: sv = (periodic tiny-prime pattern) + (log |Q(x)|
   approximation), in 64-bit words as lasieve does (so byte overflows carry
   exactly as there).  Segment i of the interval, from sinit_tab[2i] to
   sinit_tab[2i+2], gets 0x0101..01 * (lo * sinit_tab[2i+1]) added; the
   table ends with an entry equal to L.  All boundaries are multiples of 64.
   pat[0..period) is the pattern (period a multiple of 8); this function
   first extends it by 64 bytes, so pat needs period + 64 bytes.
   --------------------------------------------------------------------------- */
static inline void kmpqs_sieve_init(uint8_t *sv, uint32_t L, const uint16_t *sinit_tab,
	uint8_t *pat, uint32_t period, uint8_t lo)
{
	uint32_t i, k, x = 0, o = 0;

	for (k = 0; k < 64; k++)
		pat[period + k] = pat[k % period];

	for (i = 0; sinit_tab[2 * i] != L; i++) {
		uint32_t end = sinit_tab[2 * i + 2];
		uint64_t mask = 0x0101010101010101ULL * (uint64_t)(lo * sinit_tab[2 * i + 1]);
#if defined(KMPQS_AVX512)
		__m512i vm = _mm512_set1_epi64((long long)mask);
		for (; x < end; x += 64) {
			_mm512_storeu_si512((void *)(sv + x), _mm512_add_epi64(
				_mm512_loadu_si512((const void *)(pat + o)), vm));
			o += 64; while (o >= period) o -= period;
		}
#elif defined(KMPQS_AVX2)
		__m256i vm = _mm256_set1_epi64x((long long)mask);
		for (; x < end; x += 32) {
			_mm256_storeu_si256((__m256i *)(sv + x), _mm256_add_epi64(
				_mm256_loadu_si256((const __m256i *)(pat + o)), vm));
			o += 32; while (o >= period) o -= period;
		}
#elif defined(KMPQS_SSE2)
		__m128i vm = _mm_set1_epi64x((long long)mask);
		for (; x < end; x += 16) {
			_mm_storeu_si128((__m128i *)(sv + x), _mm_add_epi64(
				_mm_loadu_si128((const __m128i *)(pat + o)), vm));
			o += 16; if (o >= period) o -= period;
		}
#else
		for (; x < end; x += 8) {
			uint64_t w;
			memcpy(&w, pat + o, 8);
			w += mask;
			memcpy(sv + x, &w, 8);
			o += 8; if (o >= period) o -= period;
		}
#endif
	}
}

/* ---------------------------------------------------------------------------
   Survivor scan: every byte of sv[0..L) with the top bit set is a survivor.
   Appends survivor positions to buf (in order) and zeroes sv[0..L).
   Like lasieve, after each 4-byte word that held survivors it stops
   collecting once more than nmax-8 are found.  L is a multiple of 64.
   --------------------------------------------------------------------------- */
static inline uint32_t kmpqs_eval_scan(uint8_t *sv, uint32_t L, uint16_t *buf,
	uint32_t nmax)
{
	uint32_t n = 0, base;

	for (base = 0; base < L; base += 64) {
		uint64_t m;
#if defined(KMPQS_AVX512)
		__m512i v = _mm512_loadu_si512((const void *)(sv + base));
		m = (uint64_t)_mm512_movepi8_mask(v);
		_mm512_storeu_si512((void *)(sv + base), _mm512_setzero_si512());
#elif defined(KMPQS_AVX2)
		__m256i v0 = _mm256_loadu_si256((const __m256i *)(sv + base));
		__m256i v1 = _mm256_loadu_si256((const __m256i *)(sv + base + 32));
		m = (uint64_t)(uint32_t)_mm256_movemask_epi8(v0) |
			((uint64_t)(uint32_t)_mm256_movemask_epi8(v1) << 32);
		_mm256_storeu_si256((__m256i *)(sv + base), _mm256_setzero_si256());
		_mm256_storeu_si256((__m256i *)(sv + base + 32), _mm256_setzero_si256());
#elif defined(KMPQS_SSE2)
		__m128i v0 = _mm_loadu_si128((const __m128i *)(sv + base));
		__m128i v1 = _mm_loadu_si128((const __m128i *)(sv + base + 16));
		__m128i v2 = _mm_loadu_si128((const __m128i *)(sv + base + 32));
		__m128i v3 = _mm_loadu_si128((const __m128i *)(sv + base + 48));
		m = (uint64_t)(uint32_t)_mm_movemask_epi8(v0) |
			((uint64_t)(uint32_t)_mm_movemask_epi8(v1) << 16) |
			((uint64_t)(uint32_t)_mm_movemask_epi8(v2) << 32) |
			((uint64_t)(uint32_t)_mm_movemask_epi8(v3) << 48);
		_mm_storeu_si128((__m128i *)(sv + base), _mm_setzero_si128());
		_mm_storeu_si128((__m128i *)(sv + base + 16), _mm_setzero_si128());
		_mm_storeu_si128((__m128i *)(sv + base + 32), _mm_setzero_si128());
		_mm_storeu_si128((__m128i *)(sv + base + 48), _mm_setzero_si128());
#else
		uint32_t k, b;
		m = 0;
		for (k = 0; k < 64; k += 8) {
			uint64_t w;
			memcpy(&w, sv + base + k, 8);
			if (w & 0x8080808080808080ULL)
				for (b = 0; b < 8; b++)
					m |= (uint64_t)(sv[base + k + b] >> 7) << (k + b);
		}
		memset(sv + base, 0, 64);
#endif
		while (m) {
			/* all survivors of the 4-byte word holding the lowest one */
			int w = kmpqs_ctz64(m) & ~3;
			uint64_t wm = m & (0xfULL << w);

			m ^= wm;
			while (wm) {
				buf[n++] = (uint16_t)(base + kmpqs_ctz64(wm));
				wm &= wm - 1;
			}
			if (n > nmax - 8) {
				if (base + 64 < L)
					memset(sv + base + 64, 0, L - base - 64);
				return n;
			}
		}
	}
	return n;
}

/* ---------------------------------------------------------------------------
   Root update for the next polynomial (Gray code step):
     plus:  s = (s + add) mod p      minus:  s = (s + p - add) mod p
   for the elements e in [2, ne) of fbs (both roots of primes 1..ne/2-1).
   Values are < 2^15, so s + add fits in 16 bits.
   --------------------------------------------------------------------------- */
static inline void kmpqs_np_update(uint16_t *fbs, const uint16_t *p2,
	const uint16_t *add2, uint32_t ne, int minus)
{
	uint32_t e = 2;

#if defined(KMPQS_AVX512)
	for (; e < ne; e += 32) {
		__mmask32 k = (ne - e >= 32) ? (__mmask32)0xffffffff :
			(__mmask32)((1u << (ne - e)) - 1);
		__m512i s = _mm512_maskz_loadu_epi16(k, fbs + e);
		__m512i p = _mm512_maskz_loadu_epi16(k, p2 + e);
		__m512i a = _mm512_maskz_loadu_epi16(k, add2 + e);
		if (minus) a = _mm512_sub_epi16(p, a);
		s = _mm512_add_epi16(s, a);
		s = _mm512_min_epu16(s, _mm512_sub_epi16(s, p));
		_mm512_mask_storeu_epi16(fbs + e, k, s);
	}
#else
#if defined(KMPQS_AVX2)
	for (; e + 16 <= ne; e += 16) {
		__m256i s = _mm256_loadu_si256((const __m256i *)(fbs + e));
		__m256i p = _mm256_loadu_si256((const __m256i *)(p2 + e));
		__m256i a = _mm256_loadu_si256((const __m256i *)(add2 + e));
		if (minus) a = _mm256_sub_epi16(p, a);
		s = _mm256_add_epi16(s, a);
		s = _mm256_min_epu16(s, _mm256_sub_epi16(s, p));
		_mm256_storeu_si256((__m256i *)(fbs + e), s);
	}
#endif
#if defined(KMPQS_SSE41)
	for (; e + 8 <= ne; e += 8) {
		__m128i s = _mm_loadu_si128((const __m128i *)(fbs + e));
		__m128i p = _mm_loadu_si128((const __m128i *)(p2 + e));
		__m128i a = _mm_loadu_si128((const __m128i *)(add2 + e));
		if (minus) a = _mm_sub_epi16(p, a);
		s = _mm_add_epi16(s, a);
		s = _mm_min_epu16(s, _mm_sub_epi16(s, p));
		_mm_storeu_si128((__m128i *)(fbs + e), s);
	}
#endif
	for (; e < ne; e++) {
		uint32_t p = p2[e];
		uint32_t s = (uint32_t)fbs[e] + (minus ? p - add2[e] : add2[e]);
		fbs[e] = (uint16_t)(s >= p ? s - p : s);
	}
#endif
}

/* ---------------------------------------------------------------------------
   New polynomial A (lasieve's asm_next_pol10/11): for primes i in [1, nfb),
   with 16-bit Montgomery reduction  red(h) = (h + ((h*mmi) mod 2^16)*p) / 2^16:
     add_j = red(inv_j[i] * bmul[j]) mod p           for j < nadiv
     pi    = red(...red(pi0 * inv_0[i])... * inv_{nadiv-1}[i])
     c     = (disp[i] - (sum_j add_j)/2) mod p        (halving mod p)
     r     = red(sqrt_i * pi) mod p
     roots = (c + r) mod p, (c - r) mod p
   writes add_j to add[j][i] and to add2[j][2i], add2[j][2i+1], and the roots
   to fbs[2i], fbs[2i+1].  fb holds (p, sqrt) pairs; mmi = -p^-1 mod 2^16.
   All intermediate values stay below 2^31 for p < 2^15.
   --------------------------------------------------------------------------- */
#define KMPQS_MMRED(h, mmi, p) (((h) + ((((h) & 0xffff) * (mmi)) & 0xffff) * (p)) >> 16)

static inline void kmpqs_np_newA_scalar(uint32_t i, uint32_t nfb, const uint16_t *fb,
	const uint16_t *mm_inv, const uint16_t *disp, uint32_t nadiv,
	const uint16_t *const *inv, uint16_t *const *add, uint16_t *const *add2,
	const uint32_t *bmul, uint32_t pi0, uint16_t *fbs)
{
	for (; i < nfb; i++) {
		uint32_t p = fb[2 * i], mmi = mm_inv[i], c = 0, pi = pi0, cc, j;

		for (j = 0; j < nadiv; j++) {
			cc = KMPQS_MMRED((uint32_t)inv[j][i] * bmul[j], mmi, p);
			if (cc >= p) cc -= p;
			add[j][i] = (uint16_t)cc;
			add2[j][2 * i] = add2[j][2 * i + 1] = (uint16_t)cc;
			c += cc; if (c >= p) c -= p;
			pi = KMPQS_MMRED(pi * inv[j][i], mmi, p);
		}
		if (c & 1) c += p;
		c >>= 1;
		c = disp[i] + (p - c); if (c >= p) c -= p;
		cc = KMPQS_MMRED((uint32_t)fb[2 * i + 1] * pi, mmi, p);
		if (cc >= p) cc -= p;
		pi = p - cc;		/* second root offset */
		cc += c; if (cc >= p) cc -= p;
		pi += c; if (pi >= p) pi -= p;
		fbs[2 * i] = (uint16_t)cc;
		fbs[2 * i + 1] = (uint16_t)pi;
	}
}

#if defined(KMPQS_AVX512)
static inline __m512i kmpqs_mmred512(__m512i h, __m512i mmi, __m512i p)
{
	__m512i m16 = _mm512_set1_epi32(0xffff);
	__m512i t = _mm512_and_si512(_mm512_mullo_epi32(_mm512_and_si512(h, m16), mmi), m16);
	return _mm512_srli_epi32(_mm512_add_epi32(h, _mm512_mullo_epi32(t, p)), 16);
}
/* x >= p ? x-p : x, for x < 2^31 */
static inline __m512i kmpqs_csub512(__m512i x, __m512i p)
{
	return _mm512_min_epu32(x, _mm512_sub_epi32(x, p));
}
#endif
#if defined(KMPQS_AVX2)
static inline __m256i kmpqs_mmred256(__m256i h, __m256i mmi, __m256i p)
{
	__m256i m16 = _mm256_set1_epi32(0xffff);
	__m256i t = _mm256_and_si256(_mm256_mullo_epi32(_mm256_and_si256(h, m16), mmi), m16);
	return _mm256_srli_epi32(_mm256_add_epi32(h, _mm256_mullo_epi32(t, p)), 16);
}
static inline __m256i kmpqs_csub256(__m256i x, __m256i p)
{
	return _mm256_min_epu32(x, _mm256_sub_epi32(x, p));
}
#endif

static inline void kmpqs_np_newA(uint32_t nfb, const uint16_t *fb,
	const uint16_t *mm_inv, const uint16_t *disp, uint32_t nadiv,
	const uint16_t *const *inv, uint16_t *const *add, uint16_t *const *add2,
	const uint32_t *bmul, uint32_t pi0, uint16_t *fbs)
{
	uint32_t i = 1;

#if defined(KMPQS_AVX512)
	{
		const __m512i lo16 = _mm512_set1_epi32(0xffff);
		for (; i < nfb; i += 16) {
			__mmask16 k = (nfb - i >= 16) ? (__mmask16)0xffff : (__mmask16)((1u << (nfb - i)) - 1);
			__m512i ps = _mm512_maskz_loadu_epi32(k, (const uint32_t *)fb + i);	/* p | sqrt<<16 */
			__m512i p = _mm512_and_si512(ps, lo16);
			__m512i sq = _mm512_srli_epi32(ps, 16);
			__m512i mmi = _mm512_cvtepu16_epi32(_mm256_maskz_loadu_epi16(k, mm_inv + i));
			__m512i ds = _mm512_cvtepu16_epi32(_mm256_maskz_loadu_epi16(k, disp + i));
			__m512i c = _mm512_setzero_si512(), pi = _mm512_set1_epi32((int)pi0), cc, r1, r2;
			uint32_t j;

			for (j = 0; j < nadiv; j++) {
				__m512i iv = _mm512_cvtepu16_epi32(_mm256_maskz_loadu_epi16(k, inv[j] + i));
				cc = kmpqs_mmred512(_mm512_mullo_epi32(iv, _mm512_set1_epi32((int)bmul[j])), mmi, p);
				cc = kmpqs_csub512(cc, p);
				_mm256_mask_storeu_epi16(add[j] + i, k, _mm512_cvtepi32_epi16(cc));
				_mm512_mask_storeu_epi32((uint32_t *)add2[j] + i, k,
					_mm512_or_si512(cc, _mm512_slli_epi32(cc, 16)));
				c = kmpqs_csub512(_mm512_add_epi32(c, cc), p);
				pi = kmpqs_mmred512(_mm512_mullo_epi32(pi, iv), mmi, p);
			}
			/* halve mod p: odd c gets +p first */
			c = _mm512_add_epi32(c, _mm512_and_si512(p, _mm512_sub_epi32(_mm512_setzero_si512(),
				_mm512_and_si512(c, _mm512_set1_epi32(1)))));
			c = _mm512_srli_epi32(c, 1);
			c = kmpqs_csub512(_mm512_sub_epi32(_mm512_add_epi32(ds, p), c), p);
			r1 = kmpqs_csub512(kmpqs_mmred512(_mm512_mullo_epi32(sq, pi), mmi, p), p);
			r2 = _mm512_sub_epi32(p, r1);
			r1 = kmpqs_csub512(_mm512_add_epi32(r1, c), p);
			r2 = kmpqs_csub512(_mm512_add_epi32(r2, c), p);
			_mm512_mask_storeu_epi32((uint32_t *)fbs + i, k,
				_mm512_or_si512(r1, _mm512_slli_epi32(r2, 16)));
		}
	}
#elif defined(KMPQS_AVX2)
	{
		const __m256i lo16 = _mm256_set1_epi32(0xffff);
		for (; i + 8 <= nfb; i += 8) {
			__m256i ps = _mm256_loadu_si256((const __m256i *)((const uint32_t *)fb + i));
			__m256i p = _mm256_and_si256(ps, lo16);
			__m256i sq = _mm256_srli_epi32(ps, 16);
			__m256i mmi = _mm256_cvtepu16_epi32(_mm_loadu_si128((const __m128i *)(mm_inv + i)));
			__m256i ds = _mm256_cvtepu16_epi32(_mm_loadu_si128((const __m128i *)(disp + i)));
			__m256i c = _mm256_setzero_si256(), pi = _mm256_set1_epi32((int)pi0), cc, r1, r2;
			uint32_t j;

			for (j = 0; j < nadiv; j++) {
				__m256i iv = _mm256_cvtepu16_epi32(_mm_loadu_si128((const __m128i *)(inv[j] + i)));
				__m256i dup;
				cc = kmpqs_mmred256(_mm256_mullo_epi32(iv, _mm256_set1_epi32((int)bmul[j])), mmi, p);
				cc = kmpqs_csub256(cc, p);
				/* 8 x 16-bit: pack (values < 2^16, so unsigned saturation is exact) */
				_mm_storeu_si128((__m128i *)(add[j] + i), _mm_packus_epi32(
					_mm256_castsi256_si128(cc), _mm256_extracti128_si256(cc, 1)));
				dup = _mm256_or_si256(cc, _mm256_slli_epi32(cc, 16));
				_mm256_storeu_si256((__m256i *)((uint32_t *)add2[j] + i), dup);
				c = kmpqs_csub256(_mm256_add_epi32(c, cc), p);
				pi = kmpqs_mmred256(_mm256_mullo_epi32(pi, iv), mmi, p);
			}
			c = _mm256_add_epi32(c, _mm256_and_si256(p, _mm256_sub_epi32(_mm256_setzero_si256(),
				_mm256_and_si256(c, _mm256_set1_epi32(1)))));
			c = _mm256_srli_epi32(c, 1);
			c = kmpqs_csub256(_mm256_sub_epi32(_mm256_add_epi32(ds, p), c), p);
			r1 = kmpqs_csub256(kmpqs_mmred256(_mm256_mullo_epi32(sq, pi), mmi, p), p);
			r2 = _mm256_sub_epi32(p, r1);
			r1 = kmpqs_csub256(_mm256_add_epi32(r1, c), p);
			r2 = kmpqs_csub256(_mm256_add_epi32(r2, c), p);
			_mm256_storeu_si256((__m256i *)((uint32_t *)fbs + i),
				_mm256_or_si256(r1, _mm256_slli_epi32(r2, 16)));
		}
	}
#endif
	kmpqs_np_newA_scalar(i, nfb, fb, mm_inv, disp, nadiv, inv, add, add2, bmul, pi0, fbs);
}

/* ---------------------------------------------------------------------------
   Gaussian elimination over GF(2), lasieve's 4-columns-at-a-time method.
   The m x n matrix is stored word-major: bit c of row r is bit (c & 63) of
   M[r + m*(c >> 6)], so for a fixed word the rows are contiguous.  On
   return d[c] is the pivot row of column c, or -1 for a free column
   (d[] must be -1 on entry).  m <= KMPQS_GAUSS_MAXROWS.
   --------------------------------------------------------------------------- */
#define KMPQS_GAUSS_MAXROWS 2048

static inline void kmpqs_rowechelon(uint64_t *M, uint32_t m, uint32_t n,
	uint32_t n64, short *d)
{
	uint32_t i, j, k, l, t, k64, col, row, r0, c0, zz, sh;
	uint64_t tmp, tab[16];
	uint8_t ucm[KMPQS_GAUSS_MAXROWS + 64], mm[4], mk, z;

	col = 0; row = 0;
	while (1) {
		if (row >= m) break;
		if (col >= n) break;

		k64 = col >> 6;
		r0 = row; c0 = col; sh = c0 & 63;
		for (i = 0; i < 4; i++, col++) {
			for (j = row; j < m; j++) {
				mk = (uint8_t)(0xf & (M[j + m * k64] >> sh));
				if (!mk) continue;
				for (k = 0, z = 1; k < i; k++, z += z) if (mk & z) mk ^= mm[k];
				if (mk & z) break;
			}
			if (j >= m) {
				d[col] = -1;
				continue;
			}
			if (j > row) {
				for (l = k64; l < n64; l++) {
					tmp = M[j + m * l];
					M[j + m * l] = M[row + m * l];
					M[row + m * l] = tmp;
				}
			}
			mk = (uint8_t)(M[row + m * k64] >> sh);
			for (k = 0, z = 1; k < i; k++, z += z) {
				if (mk & z) {
					mk ^= mm[k];
					for (l = k64; l < n64; l++)
						M[row + m * l] ^= M[d[c0 + k] + m * l];
				}
			}
			mm[i] = mk;
			for (k = 0; k < i; k++) {
				if (d[c0 + k] == -1) continue;
				if (mm[k] & z) {
					mm[k] ^= mm[i];
					for (l = k64; l < n64; l++)
						M[d[c0 + k] + m * l] ^= M[row + m * l];
				}
			}
			d[col] = (short)row;
			row++;
			if (row >= m) break;
		}

		/* every row gets the combination of the (up to 4) pivot rows its
		   bits c0..c0+3 select; the pivot rows themselves are skipped */
		j = 0;
#if defined(KMPQS_AVX512)
		{
			const __m512i vsh = _mm512_set1_epi64(sh), m4 = _mm512_set1_epi64(0xf);
			for (; j + 8 <= m; j += 8) {
				__m512i w = _mm512_loadu_si512((const void *)(M + j + m * k64));
				w = _mm512_and_si512(_mm512_srlv_epi64(w, vsh), m4);
				_mm_storel_epi64((__m128i *)(ucm + j), _mm512_cvtepi64_epi8(w));
			}
		}
#endif
		for (; j < m; j++)
			ucm[j] = (uint8_t)((M[j + m * k64] >> sh) & 0xf);
		for (j = r0; j < row; j++) ucm[j] = 0;
		for (l = k64; l < n64; l++) {
			uint64_t *Ml = M + m * l;

			tab[0] = 0ULL;
			for (j = 0, zz = 1; j < 4; j++, zz += zz) {
				if (d[c0 + j] == -1) tab[zz] = 0ULL;
				else tab[zz] = Ml[d[c0 + j]];
				for (k = 1; k < zz; k++) tab[zz + k] = tab[k] ^ tab[zz];
			}
			t = 0;
#if defined(KMPQS_AVX512)
			{
				const __m512i t0 = _mm512_loadu_si512((const void *)tab);
				const __m512i t1 = _mm512_loadu_si512((const void *)(tab + 8));
				for (; t + 8 <= m; t += 8) {
					__m512i idx = _mm512_cvtepu8_epi64(_mm_loadl_epi64((const __m128i *)(ucm + t)));
					__m512i x = _mm512_permutex2var_epi64(t0, idx, t1);
					_mm512_storeu_si512((void *)(Ml + t),
						_mm512_xor_si512(_mm512_loadu_si512((const void *)(Ml + t)), x));
				}
			}
#endif
			for (; t < m; t++)
				Ml[t] ^= tab[ucm[t]];
		}
	}
}

/* ---------------------------------------------------------------------------
   Trial division by roots: for the sieve position ind, write the indices j
   in [1, jend) with ind == s1 or s2 (mod p_j) to out, ascending; returns
   how many.  Exact when ind + p_j < 2^16 (the caller checks
   sievelen + pmax < 2^16): with x = ind + p - s < 2^16 and y = x*inv mod
   2^16, p | x  <=>  y*p < 2^16  <=>  mulhi(y, p) == 0.
   --------------------------------------------------------------------------- */
static inline uint32_t kmpqs_td_roots(uint32_t ind, const uint16_t *fbs,
	const uint16_t *p2, const uint16_t *inv2, uint32_t jend, uint16_t *out)
{
	uint32_t n = 0, e = 2, ne = 2 * jend;

#if defined(KMPQS_AVX512)
	__m512i vind = _mm512_set1_epi16((short)ind);
#if defined(__AVX512VBMI2__)
	/* lane l of a block holds root l&1 of prime e/2 + l/2 */
	static const uint16_t lane_prime[32] = {
		0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7,
		8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15 };
	const __m512i vlp = _mm512_loadu_si512((const void *)lane_prime);
#endif
	for (; e < ne; e += 32) {
		__mmask32 k = (ne - e >= 32) ? (__mmask32)0xffffffff :
			(__mmask32)((1u << (ne - e)) - 1);
		__m512i s = _mm512_maskz_loadu_epi16(k, fbs + e);
		__m512i p = _mm512_maskz_loadu_epi16(k, p2 + e);
		__m512i iv = _mm512_maskz_loadu_epi16(k, inv2 + e);
		__m512i x = _mm512_sub_epi16(_mm512_add_epi16(vind, p), s);
		__m512i h = _mm512_mulhi_epu16(_mm512_mullo_epi16(x, iv), p);
		uint32_t m = (uint32_t)_mm512_mask_cmpeq_epi16_mask(k, h, _mm512_setzero_si512());
		/* one bit per prime: either root */
		m = (m | (m >> 1)) & 0x55555555u;
#if defined(__AVX512VBMI2__)
		/* write all hits at once; out needs 32 entries of slack */
		_mm512_storeu_si512((void *)(out + n), _mm512_maskz_compress_epi16((__mmask32)m,
			_mm512_add_epi16(vlp, _mm512_set1_epi16((short)(e >> 1)))));
		n += (uint32_t)_mm_popcnt_u32(m);
#else
		while (m) {
			out[n++] = (uint16_t)((e >> 1) + (kmpqs_ctz32(m) >> 1));
			m &= m - 1;
		}
#endif
	}
#else
#if defined(KMPQS_AVX2)
	__m256i vind = _mm256_set1_epi16((short)ind);
	for (; e + 16 <= ne; e += 16) {
		__m256i s = _mm256_loadu_si256((const __m256i *)(fbs + e));
		__m256i p = _mm256_loadu_si256((const __m256i *)(p2 + e));
		__m256i iv = _mm256_loadu_si256((const __m256i *)(inv2 + e));
		__m256i x = _mm256_sub_epi16(_mm256_add_epi16(vind, p), s);
		__m256i h = _mm256_mulhi_epu16(_mm256_mullo_epi16(x, iv), p);
		uint32_t m = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi16(h, _mm256_setzero_si256()));
		/* 2 mask bits per root, 4 per prime */
		m |= m >> 2;
		m &= 0x11111111u;
		while (m) {
			out[n++] = (uint16_t)((e >> 1) + (kmpqs_ctz32(m) >> 2));
			m &= m - 1;
		}
	}
#endif
#if defined(KMPQS_SSE2)
	{
		__m128i vind = _mm_set1_epi16((short)ind);
		for (; e + 8 <= ne; e += 8) {
			__m128i s = _mm_loadu_si128((const __m128i *)(fbs + e));
			__m128i p = _mm_loadu_si128((const __m128i *)(p2 + e));
			__m128i iv = _mm_loadu_si128((const __m128i *)(inv2 + e));
			__m128i x = _mm_sub_epi16(_mm_add_epi16(vind, p), s);
			__m128i h = _mm_mulhi_epu16(_mm_mullo_epi16(x, iv), p);
			uint32_t m = (uint32_t)_mm_movemask_epi8(_mm_cmpeq_epi16(h, _mm_setzero_si128()));
			m |= m >> 2;
			m &= 0x1111u;
			while (m) {
				out[n++] = (uint16_t)((e >> 1) + (kmpqs_ctz32(m) >> 2));
				m &= m - 1;
			}
		}
	}
#endif
	for (; e < ne; e += 2) {
		uint32_t p = p2[e], iv = inv2[e];
		uint32_t x1 = (ind + p - fbs[e]) & 0xffff, x2 = (ind + p - fbs[e + 1]) & 0xffff;
		uint32_t h1 = (((x1 * iv) & 0xffff) * p) >> 16;
		uint32_t h2 = (((x2 * iv) & 0xffff) * p) >> 16;
		if (h1 == 0 || h2 == 0)
			out[n++] = (uint16_t)(e >> 1);
	}
#endif
	return n;
}

#endif
