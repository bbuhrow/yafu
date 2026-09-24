/* byteio.h -- portable little-endian read/write helpers.
 *
 * Never relies on host endianness or struct layout/padding: every
 * multi-byte field is assembled/disassembled one byte at a time, so
 * the on-disk format is identical on big- and little-endian hosts,
 * and under gcc, clang or MSVC.
 */
#ifndef COLLHARNESS_BYTEIO_H
#define COLLHARNESS_BYTEIO_H

#include <stdio.h>
#include <stdint.h>

static inline void
wr_u32(FILE *f, uint32_t v)
{
	unsigned char b[4];
	b[0] = (unsigned char)(v & 0xFF);
	b[1] = (unsigned char)((v >> 8) & 0xFF);
	b[2] = (unsigned char)((v >> 16) & 0xFF);
	b[3] = (unsigned char)((v >> 24) & 0xFF);
	fwrite(b, 1, 4, f);
}

static inline void
wr_u64(FILE *f, uint64_t v)
{
	wr_u32(f, (uint32_t)(v & 0xFFFFFFFFu));
	wr_u32(f, (uint32_t)((v >> 32) & 0xFFFFFFFFu));
}

static inline void
wr_i64(FILE *f, int64_t v)
{
	wr_u64(f, (uint64_t)v);
}

static inline void
wr_bytes(FILE *f, const void *p, size_t n)
{
	fwrite(p, 1, n, f);
}

/* Return 1 on success, 0 on short read / EOF. */
static inline int
rd_u32(FILE *f, uint32_t *out)
{
	unsigned char b[4];
	if (fread(b, 1, 4, f) != 4)
		return 0;
	*out = (uint32_t)b[0] | ((uint32_t)b[1] << 8) |
			((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
	return 1;
}

static inline int
rd_u64(FILE *f, uint64_t *out)
{
	uint32_t lo, hi;
	if (!rd_u32(f, &lo) || !rd_u32(f, &hi))
		return 0;
	*out = (uint64_t)lo | ((uint64_t)hi << 32);
	return 1;
}

static inline int
rd_i64(FILE *f, int64_t *out)
{
	uint64_t u;
	if (!rd_u64(f, &u))
		return 0;
	*out = (int64_t)u;
	return 1;
}

static inline int
rd_bytes(FILE *f, void *p, size_t n)
{
	return fread(p, 1, n, f) == n;
}

#endif /* COLLHARNESS_BYTEIO_H */
