/*  Written in 2016 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */
#ifndef XOROSHIRO64STAR_H
#define XOROSHIRO64STAR_H
#include <stdint.h>

/* This is xoroshiro64* 1.0, our best and fastest 32-bit small-state
   generator for 32-bit floating-point numbers. We suggest to use its
   upper bits for floating-point generation, as it is slightly faster than
   xoroshiro64**. It passes all tests we are aware of except for linearity
   tests, as the lowest six bits have low linear complexity, so if low
   linear complexity is not considered an issue (as it is usually the
   case) it can be used to generate 32-bit outputs, too.

   We suggest to use a sign test to extract a random Boolean value, and
   right shifts to extract subsets of bits.

   The state must be seeded so that it is not everywhere zero. */


static inline uint32_t rotl32(const uint32_t x, int k) {
	return (x << k) | (x >> (32 - k));
}


static uint32_t seed[2];

static uint8_t short_rand;

static uint32_t long_rand;

uint32_t Next(void) {
	const uint32_t s0 = seed[0];
	uint32_t s1 = seed[1];
	long_rand = s0 * 0x9E3779BB;

	s1 ^= s0;
	seed[0] = rotl32(s0, 26) ^ s1 ^ (s1 << 9); // a, b
	seed[1] = rotl32(s1, 13); // c

	return long_rand;
}

static inline uint8_t NextShort(void) {
	const uint32_t s0 = seed[0];
	uint32_t s1 = seed[1];
	short_rand = s0 * 0x9E3779BB;

	s1 ^= s0;
	seed[0] = rotl32(s0, 26) ^ s1 ^ (s1 << 9); // a, b
	seed[1] = rotl32(s1, 13); // c

	return short_rand;
}
#endif