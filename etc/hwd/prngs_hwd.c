/* Can be used by hwd. */
static inline uint64_t rotate64(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

static inline uint32_t rotate32(const uint32_t x, int k) {
	return (x << k) | (x >> (32 - k));
}

/*
static inline uint32_t rotate32(uint32_t value, unsigned int rot)
{
#if __clang__ && (__x86_64__  || __i386__)
    asm ("roll   %%cl, %0" : "=r" (value) : "0" (value), "c" (rot));
    return value;
#else
    return (value << rot) | (value >> (32 - rot));
#endif
}
*/


//minimover64
//static uint64_t x = UINT64_C(0x9E3779B9);
//static inline uint64_t next()
//{
//    return (x = rotate64(x, 21) * UINT64_C(0x9E3779B9)) * UINT64_C(0x41C64E6D);
//}

//movercounter64
static uint64_t a = UINT64_C(0x9E3779B9), b = UINT64_C(1);
static inline uint64_t next()
{
    return (a = rotate64(a, 21) * (b += UINT64_C(0x9E3779B97F4A7AF6))) * UINT64_C(0x41C64E6B);
}
//uint64_t z = (x = (x ^ UINT64_C(0x6C8E9CF570932BD5)) * UINT64_C(0x41C64E6B));
//z -= z >> 28;
//return z ^ z >> 26;

//    const uint64_t z = x,
//                   s = (x = (z ^ UINT64_C(0x6C8E9CF570932BD5)) * UINT64_C(0x41C64E6B)) - (z >> 31);
//    return s ^ s >> 26;

// thrustalt
//static uint64_t x = UINT64_C(1);
//static inline uint64_t next()
//{
//    const uint64_t s = (x += UINT64_C(0x6C8E9CF570932BD5));
//    const uint64_t z = (s ^ (s >> 25)) * (s | UINT64_C(0xA529));
//    return z ^ (z >> 22);
//}
// linnorm
/*
static uint64_t x = 1ULL;
static inline uint64_t next()
{
    uint64_t z = (x = x * UINT64_C(1103515245) + UINT64_C(1));
	z = (x ^ x >> 27) * UINT64_C(0xAEF17502108EF2D9);
    return z ^ z >> 25;
}
*/
/*
static uint64_t x = 1ULL;
static inline uint64_t next()
{
    uint64_t z = (((x += UINT64_C(0x632BE59BD9B4E019)) ^ UINT64_C(0x9E3779B97F4A7C15))); // * UINT64_C(0xC6BC279692B5CC83)
    z = (z ^ z >> 27) * UINT64_C(0xAEF17502108EF2D9);
    return z ^ z >> 25;
}
*/
// QuixoticRNG
/*
static uint64_t x = 1ULL;
static inline uint64_t next()
{
    uint64_t z = (x = (x ^ UINT64_C(0x6C8E9CF570932BD5)) * UINT64_C(0x2E4DD58B));//UINT64_C(0x42042E4DD58B));//0x5DA942042E4DD58B
	z = (z ^ z >> 27u) + UINT64_C(0xAEF17502108EF2D9);
	return z ^ z >> 25u;
}
*/
/*
// mcg 128
#define __GLIBCXX_BITSIZE_INT_N_0 128
#define __GLIBCXX_TYPE_INT_N_0 __int128

static __uint128_t state = 1;   // can be seeded to any odd number

static inline uint64_t next()
{
    return (state *= 0xda942042e4dd58b5ULL) >> 64;
}
*/
/*
static uint32_t x = 1ULL, y = 1ULL;
static uint32_t inline next() {
	const uint32_t s0 = x;
	uint32_t s1 = y;
	const uint32_t result = s0 + s1; //s0 * 0xA5CB3 ^ s1

	s1 ^= s0;
	x = rotate32(s0, 26) ^ s1 ^ (s1 << 9); // a, b
	y = rotate32(s1, 13); // c

	return result;
}
*/
/*
static uint32_t x = 0xFFFFFFU, y = 0xFFFFFFFU;
static uint32_t inline next() {
        const uint32_t s0 = x;
        uint32_t s1 = y;
        const uint32_t result = s0 + s1;
        s1 ^= s0;
        //x = (s0 << 13 | s0 >> 19) ^ s1 ^ (s1 << 5);
        //y = (s1 << 28 | s1 >> 4);
        //return (result << 10 | result >> 22) + s0;
        
        x = rotate32(s0, 26u) ^ s1 ^ (s1 << 9);
        y = rotate32(s1, 13u);
        return rotate32(result, 10u) + s0;
}*/
/*
        x = rotate32(s0, 26u) ^ s1 ^ (s1 << 9);
        y = rotate32(s1, 13u);
        return rotate32(result, 10u) + s0;


        x = rotate32(s0, 13u) ^ s1 ^ (s1 << 5u);
        y = rotate32(s1, 28u);
        return rotate32(result, 10u) + s0;
*/
/*
        x = (s0 << 13 | s0 >> 19) ^ s1 ^ (s1 << 5);
        y = (s1 << 28 | s1 >> 4);
        return (result << 10 | result >> 22) + s0;
*/

// xoroshiro128+
/*
static uint64_t s[2] = {1, 2};

static inline uint64_t next(void) {
	const uint64_t s0 = s[0];
	uint64_t s1 = s[1];
	const uint64_t result = s0 + s1;

	s1 ^= s0;
	s[0] = rotate64(s0, 24) ^ s1 ^ (s1 << 16); // a, b
	s[1] = rotate64(s1, 37); // c

	return result;
	//return (result ^ result >> 23) + (s[2] += UINT64_C(0x632BE59BD9B4E019));
}
*/

//xoshiro256+
/*
static uint64_t s[4] = {1, 2, 3, 4};

uint64_t next(void) {
	const uint64_t result = s[0] + s[3];
	const uint64_t t = s[1] << 17;
	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];
	s[2] ^= t;
	s[3] = rotate64(s[3], 45);
	return result;
}
*/
/*
// Xoshiro Star Phi
static uint32_t s[4] = {1, 2, 3, 4};

uint32_t next(void) {
	const uint32_t result = s[1] * UINT32_C(31);
	const uint32_t t = s[1] << 9;
	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];
	s[2] ^= t;
	s[3] = rotate32(s[3], 11);
	return rotate32(result, 23) + UINT32_C(0x9E3779B9);
}
*/
/* // Mover64, first formulation
static uint64_t state[2] = {1, 1};
uint64_t next(void) {
    const uint64_t a = state[0] * UINT64_C(0x41C64E6B);
    const uint64_t b = state[1] * UINT64_C(0x9E3779B9);
    return (state[0] = rotate64(a, 28)) ^ (state[1] = rotate64(b, 37));
}
*/
/* Xorshift128 with Lerosu2 scrambler
static uint64_t state[2] = {1, 1};
static inline uint64_t next(void) {
	uint64_t s1 = state[0];
	const uint64_t s0 = state[1];
	const uint64_t result = (s1 << 6) - rotate64(s1, 4); 
	state[0] = s0;
	s1 ^= s1 << 23; // a
	state[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5); // b, c
	return (result << 11) - rotate64(result, 8);
}
*/
/*
//SeaSlater64, xoroshiro128 with Lerosu2 scrambler
static uint64_t s[2] = {1, 1};

static inline uint64_t next(void) {
	const uint64_t s0 = s[0];
	uint64_t s1 = s[1];
	const uint64_t result = (s0 << 7) - rotate64(s0, 5);

	s1 ^= s0;
	s[0] = rotate64(s0, 24) ^ s1 ^ (s1 << 16); // a, b
	s[1] = rotate64(s1, 37); // c

	return (result << 12) - rotate64(result, 9);
}
*/
/* //overdrive
static uint64_t s[2] = {1, 1};

static inline uint64_t next(void) {
    const uint64_t s0 = s[0];
    const uint64_t s1 = s[1];
    const uint64_t result = s0 ^ s1;
    s[0] = rotate64(s0, 28) * UINT64_C(0x41C64E6B);
    s[1] = UINT64_C(0xC6BC279692B5CC8B) - rotate64(s1, 35);
    return result;
}
*/