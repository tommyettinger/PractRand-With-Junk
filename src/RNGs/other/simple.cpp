#include <string>
#include <sstream>
#include "PractRand/config.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include "PractRand/rng_internals.h"

#include "PractRand/RNGs/other/simple.h"

namespace PractRand {
	using namespace Internals;
	namespace RNGs {
		namespace Polymorphic {
			namespace NotRecommended {
				Uint16 xsalta16x3::raw16() {//slightly more complex output function
					Uint16 tmp, old;
					tmp = a + c + ((b >> 11) | (b << 5));
					old = a;
					a = b ^ (b >> 3);
					b = c ^ (c << 7);
					c = c ^ (c >> 8) ^ old;
					return tmp;
				}
				std::string xsalta16x3::get_name() const { return "xsalta16x3"; }
				void xsalta16x3::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c);
					if (!a && !b && !c) a = 1;
				}
				Uint16 xsaltb16x3::raw16() {//radically different output function
					Uint16 tmp, old;
					tmp = (a & b) | (b & c) | (c & a);
					old = a;
					a = b ^ (b >> 3);
					b = c ^ (c << 7);
					c = c ^ (c >> 8) ^ old;
					return c + ((tmp << 5) | (tmp >> 11));
				}
				std::string xsaltb16x3::get_name() const { return "xsaltb16x3"; }
				void xsaltb16x3::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c);
					if (!a && !b && !c) a = 1;
				}
				Uint16 xsaltc16x3::raw16() {//deviating from the standard LFSR state function
					Uint16 old;
					old = a;
					a = b ^ (b >> 5);
					b = c + (c << 3);
					c = c ^ (c >> 7) ^ old;
					return a + b;
				}
				std::string xsaltc16x3::get_name() const { return "xsaltc16x3"; }
				void xsaltc16x3::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c);
					if (!a && !b && !c) a = 1;
				}

				Uint32 xorshift32::raw32() {
					a ^= a << 13;
					a ^= a >> 17;
					a ^= a << 5;
					return a;
				}
				std::string xorshift32::get_name() const { return "xorshift32"; }
				void xorshift32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					if (!a) a = 1;
				}
				Uint64 xorshift64::raw64() {
					a ^= a << 13;
					a ^= a >> 7;
					a ^= a << 17;
					return a;
				}
				std::string xorshift64::get_name() const { return "xorshift64"; }
				void xorshift64::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					if (!a) a = 1;
				}

				Uint64 xorshift128plus::raw64() {
					Uint64 s1 = seed0;
					Uint64 s0 = seed1;
					seed0 = s0;
					s1 ^= s1 << 23;
					seed1 = (s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26));
					//// RandomXS128; fails quickly on mostly the low 8 bits.
					return seed1 + s0;
					//// does fine to 1TB, at least; no anomalies.
					//return rotate64(seed1 * 5, 29) * 257;
					// xoshiro256** scrambler, seems to be fine.
					//return rotate64(seed1 * 5, 7) * 9;
				}
				std::string xorshift128plus::get_name() const { return "xorshift128plus"; }
				void xorshift128plus::walk_state(StateWalkingObject *walker) {
					walker->handle(seed0);
					walker->handle(seed1);
					if (!(seed0|seed1)) seed0 = 1;
				}

				void xorshift64of128::xrs(int bits) {
					if (bits < 64) {
						low ^= low >> bits;
						low ^= high << (64 - bits);
						high ^= high >> bits;
					}
					else if (bits > 64) low ^= high >> (bits - 64);
					else low ^= high;
				}
				void xorshift64of128::xls(int bits) {
					if (bits < 64) {
						high ^= high << bits;
						high ^= low >> (64 - bits);
						low ^= low << bits;
					}
					else if (bits > 64) low ^= high >> (bits - 64);
					else low ^= high;
				}
				Uint64 xorshift64of128::raw64() {
					xls(16);
					xrs(53);
					xls(47);
					return low;
				}
				std::string xorshift64of128::get_name() const { return "xorshift64of128"; }
				void xorshift64of128::walk_state(StateWalkingObject *walker) {
					walker->handle(low); walker->handle(high);
					if (!high && !low) low = 1;
				}

				std::string xorshift32of128::get_name() const { return "xorshift32of128"; }
				void xorshift32of128::walk_state(StateWalkingObject *walker) { impl.walk_state(walker); }
				std::string xorshift32of64::get_name() const { return "xorshift32of64"; }
				void xorshift32of64::walk_state(StateWalkingObject *walker) { impl.walk_state(walker); }
				std::string xorshift16of32::get_name() const { return "xorshift16of32"; }
				void xorshift16of32::walk_state(StateWalkingObject *walker) { impl.walk_state(walker); }

				Uint32 xorshift32x4::raw32() {
					Uint32 tmp = x ^ (x << 15);
					x = y;
					y = z;
					z = w;
					tmp ^= tmp >> 4;
					w ^= (w >> 21) ^ tmp;
					return w;
				}
				std::string xorshift32x4::get_name() const { return "xorshift32x4"; }
				void xorshift32x4::walk_state(StateWalkingObject *walker) {
					walker->handle(x);
					walker->handle(y);
					walker->handle(z);
					walker->handle(w);
					if (!(x || y || z || w)) x = 1;
				}

				Uint32 xorwow32of96::raw32() {
					a += 362437;
					return a + impl.raw32();
				}
				std::string xorwow32of96::get_name() const { return "xorwow32of96"; }
				void xorwow32of96::walk_state(StateWalkingObject *walker) {
					impl.walk_state(walker);
					walker->handle(a);
				}
				Uint32 xorwow32x6::raw32() {
//					//kotlin's version
//					Uint32 t = x;
//					t ^= t >> 2;
//					x = y;
//					y = z;
//					z = w;
//					Uint32 v0 = v;
//					w = v0;
//					t ^= t << 1 ^ v0 ^ v0 << 4;
//					v = t;
//					d += 362437;
//					return t + d;
					////original
					Uint32 tmp = x;
					x = y;
					y = z;
					z = w ^ (w << 1);
					w = v ^ (v >> 7);
					v ^= (v << 4) ^ tmp;
					d += 362437;
					return v + d;
				}
				std::string xorwow32x6::get_name() const { return "xorwow32x6"; }
				void xorwow32x6::walk_state(StateWalkingObject *walker) {
					walker->handle(x);
					walker->handle(y);
					walker->handle(z);
					walker->handle(w);
					walker->handle(v);
					walker->handle(d);
					if (!(x || y || z || w || v)) x = 1;
				}
				Uint64 xoroshiro128plus::raw64() {
					//const uint64_t s0 = state0;
					//uint64_t s1 = state1;
					//const uint64_t result = s0 + s1;
					//s1 ^= s0;
					//state0 = rotate64(s0, 24) ^ s1 ^ (s1 << 16); // a, b
					//state1 = rotate64(s1, 37); // c
					//return result;
					//const uint64_t s0 = state0;
					//const uint64_t s1 = state1 ^ s0;
					//const uint64_t result = (s0 << 5) - s0;
					//state0 = rotate64(s0, 24) ^ s1 ^ (s1 << 16); // a, b
					//state1 = rotate64(s1, 37); // c
					//return rotate64(result, 60) + UINT64_C(0x9E3779B97F4A7C15);

					const uint64_t s0 = state0;
					uint64_t s1 = state1;
					const uint64_t result = s0 * 5;
					s1 ^= s0;
					state0 = rotate64(s0, 24) ^ s1 ^ (s1 << 16); // a, b
					state1 = rotate64(s1, 37); // c
					return rotate64(result, 59) * 17;
				}
				std::string xoroshiro128plus::get_name() const { return "xoroshiro128plus"; }
				void xoroshiro128plus::walk_state(StateWalkingObject *walker) {
					walker->handle(state0);
					walker->handle(state1);
				}
				Uint64 xoroshiro128plus_2p64::raw64() {
					Uint64 result = state0 + state1;
					static const Uint64 JUMP[] = { 0xbeac0467eba5facb, 0xd86b048b86aa9922 };

					Uint64 s0 = 0;
					Uint64 s1 = 0;
					for (int i = 0; i < sizeof JUMP / sizeof *JUMP; i++) {
						for (int b = 0; b < 64; b++) {
							if (JUMP[i] & 1ULL << b) {
								s0 ^= state0;
								s1 ^= state1;
							}
							Uint64 tmp = state0 ^ state1;
							state0 = ((state0 << 55) | (state0 >> (64 - 55))) ^ tmp ^ (tmp << 14);
							state1 = ((tmp << 36) | (tmp >> (64 - 36)));
						}
					}
					state0 = s0;
					state1 = s1;
					return result;
				}
				std::string xoroshiro128plus_2p64::get_name() const { return "xoroshiro128plus_2p64"; }
				void xoroshiro128plus_2p64::walk_state(StateWalkingObject *walker) {
					walker->handle(state0);
					walker->handle(state1);
				}


				Uint32 sapparot::raw32() {
					Uint32 tmp;
					tmp = a + 0x9e3779b9;
					tmp = (tmp << 7) | (tmp >> 25);
					a = b ^ (~tmp) ^ (tmp << 3);
					a = (a << 7) | (a >> 25);
					b = tmp;
					return a ^ b;
				}
				std::string sapparot::get_name() const { return "sapparot"; }
				void sapparot::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
				}

				Uint16 sap16of48::raw16() {
					Uint16 tmp;
					tmp = a + 0x79b9 + c;
					tmp = (tmp << 5) | (tmp >> 11);
					a = b ^ (~tmp) ^ (tmp << 3);
					a = (a << 5) | (a >> 11);
					b = tmp;
					c = (c + a) ^ b;
					return b;
				}
				std::string sap16of48::get_name() const { return "sap16of48"; }
				void sap16of48::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				Uint32 sap32of96::raw32() {
					Uint32 tmp;
					tmp = a + 0x9e3779b9 + c;
					tmp = (tmp << 7) | (tmp >> 25);
					a = b ^ (~tmp) ^ (tmp << 3);
					a = (a << 7) | (a >> 25);
					b = tmp;
					c = (c + a) ^ b;
					return b;
				}
				std::string sap32of96::get_name() const { return "sap32of96"; }
				void sap32of96::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}

				Uint32 flea32x1::raw32() {
					enum { SHIFT1 = 15, SHIFT2 = 27 };
					Uint32 e = a[d % SIZE];
					a[d % SIZE] = ((b << SHIFT1) | (b >> (32 - SHIFT1)));
					b = c + ((d << SHIFT2) | (d >> (32 - SHIFT2)));
					c = d + a[i++ % SIZE];
					d = e + c;
					return b;
				}
				std::string flea32x1::get_name() const { return "flea32x1"; }
				void flea32x1::walk_state(StateWalkingObject *walker) {
					for (int z = 0; z < SIZE; z++) walker->handle(a[z]);
					walker->handle(b);
					walker->handle(c);
					walker->handle(d);
					walker->handle(i);
				}

				/*
				sfc version 1:
				tmp = a ^ counter++;
				a = b + (b << SHIFT1);
				b = ((b << SHIFT2) | (b >> (WORD_BITS - SHIFT2))) + tmp;
				steps:
				load a			load b			load counter
				*	b<<SHIFT1		b<<<SHIFT2		counter + 1		a ^ counter
				+ b				+ a ^ counter	store counter
				store a			store b
				overall:
				quality suffers, but it's fast
				note that (b + (b << SHIFT1)) can be a single LEA on x86 for small shifts
				constants:
				16: 2,5 ?
				32: 5,12 ?
				64: 7,41 ?
				sfc version 2:
				code:
				tmp = a ^ b;
				a = b + (b << SHIFT1);
				b = ((b << SHIFT2) | (b >> (WORD_BITS - SHIFT2))) + tmp + counter++;
				steps:
				load a			load b				load counter
				*	b << SHIFT1		b <<< SHIFT2		a ^ b			counter + 1
				+ b				+((a^b)+counter)***
				store a			store b
				suggested values for SHIFT1,SHIFT2 by output avalanche:
				16 bit: *2,5 - 2,6 - 3,5 - 3,6
				32 bit: 4,9 - 5,9 - *5,12 - 6,10 - 9,13
				64 bit: 6,40 - 7,17 - *7,41 - 10,24
				overall:
				does okay @ 32 & 64 bit, but poorly at 16 bit
				note that (b + (b << SHIFT1)) can be a single LEA on x86 for small shifts
				sfc version 3:
				code
				tmp = a + b + counter++;
				a = b ^ (b >> SHIFT1);
				b = ((b << SHIFT2) | (b >> (WORD_BITS - SHIFT2))) + tmp;
				steps:
				load a			load b			load counter
				b >> SHIFT1		b <<< SHIFT2	a+b+counter***		counter + 1
				*	^ b				+(a+b+counter)	store counter
				store a			store b
				overall:
				good statistical properties
				slower than earlier versions on my CPU, but still fast
				refence values for SHIFT1,SHIFT2:
				16 bit: 2,5
				32 bit: 5,12
				64 bit: 7,41
				sfc version 4:
				code
				Word old = a + b + counter++;
				a = b ^ (b >> SHIFT2);
				b = c + (c << SHIFT3);
				c = old + rotate(c,SHIFT1);
				return old;
				steps:
				load a			load b			load c				load counter
				*	b >> SHIFT2		b << SHIFT3		c <<< SHIFT1		counter + 1		a+b+counter***
				^ b				+ b				+(a+b+counter)		store counter
				store a			store b			store c
				overall:
				very good behavior on statistical tests
				uses an extra word / register - not as nice for inlining
				adequate speed
				refence values for SHIFT1,SHIFT2,SHIFT3:
				8 bit:  3,2,1
				16 bit: 7,5,2
				32 bit: 25,8,3
				64 bit: 25,12,3
				*/
				Uint16 sfc_v1_16::raw16() {
					Uint16 tmp = a ^ counter++;
					a = b + (b << 2);
					enum { BARREL_SHIFT = 5 };
					b = ((b << BARREL_SHIFT) | (b >> (16 - BARREL_SHIFT))) + tmp;
					return a;
				}
				std::string sfc_v1_16::get_name() const { return "sfc_v1_16"; }
				void sfc_v1_16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(counter);
				}
				Uint32 sfc_v1_32::raw32() {
					Uint32 tmp = a ^ counter++;
					a = b + (b << 5);
					enum { BARREL_SHIFT = 12 };
					b = ((b << BARREL_SHIFT) | (b >> (32 - BARREL_SHIFT))) + tmp;
					return a;
				}
				std::string sfc_v1_32::get_name() const { return "sfc_v1_32"; }
				void sfc_v1_32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(counter);
				}
				Uint16 sfc_v2_16::raw16() {
					Uint16 tmp = a ^ b;
					a = b + (b << 2);
					enum { BARREL_SHIFT = 5 };
					b = ((b << BARREL_SHIFT) | (b >> (16 - BARREL_SHIFT))) + tmp + counter++;
					return tmp;
				}
				std::string sfc_v2_16::get_name() const { return "sfc_v2_16"; }
				void sfc_v2_16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(counter);
				}
				Uint32 sfc_v2_32::raw32() {
					Uint32 tmp = a ^ b;
					a = b + (b << 5);
					enum { BARREL_SHIFT = 12 };
					b = ((b << BARREL_SHIFT) | (b >> (32 - BARREL_SHIFT))) + tmp + counter++;
					return tmp;
				}
				std::string sfc_v2_32::get_name() const { return "sfc_v2_32"; }
				void sfc_v2_32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(counter);
				}
				Uint16 sfc_v3_16::raw16() {
					Uint16 tmp = a + b + counter++;
					a = b ^ (b >> 2);
					enum { BARREL_SHIFT = 5 };
					b = ((b << BARREL_SHIFT) | (b >> (16 - BARREL_SHIFT))) + tmp;
					return tmp;
				}
				std::string sfc_v3_16::get_name() const { return "sfc_v3_16"; }
				void sfc_v3_16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(counter);
				}
				Uint32 sfc_v3_32::raw32() {
					Uint32 tmp = a + b + counter++;
					a = b ^ (b >> 5);
					enum { BARREL_SHIFT = 12 };
					b = ((b << BARREL_SHIFT) | (b >> (32 - BARREL_SHIFT))) + tmp;
					return tmp;
				}
				std::string sfc_v3_32::get_name() const { return "sfc_v3_32"; }
				void sfc_v3_32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(counter);
				}
				// A copy of O'Neill's JSF8, I think exactly.
				// I can't reproduce her results with this version of PractRand.
				// Running this command,
				// /RNG_test jsf8 -te 1 -tlmin 10 -tlmax 50 -multithreaded -seed 0
				// Fails at the 256KB mark.
				Uint8 jsf8::raw8() {
					Uint8 e = a - rotate8(b, 1);
					a = b ^ rotate8(c, 4);
					b = c + d;
					c = d + e;
					d = e + a;
					return d;
				}
				std::string jsf8::get_name() const { return "jsf8"; }
				void jsf8::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(d);
					if (!(a | b) && !(c | d)) d++;
				}
// Fails at 32GB (due to mod3n(5):(1,9-0) ) with the command:
// ./RNG_test jsf16 -tf 2 -tlmin 10 -tlmax 50 -multithreaded -seed 0
// with -te 1 instead of -tf2, this fails at 1GB on the same test.

				Uint16 jsf16::raw16() {
					Uint16 e = a - ((b << 13) | (b >> 3));
					a = b ^ ((c << 9) | (c >> 7));
					b = c + d;
					c = d + e;
					d = e + a;
					return d;
				}
				std::string jsf16::get_name() const { return "jsf16"; }
				/*seed:
					a = Uint16(s);
					b = Uint16(s >> 16);
					c = Uint16(s >> 32);
					d = Uint16(s >> 48);
					if (!(a|b) && !(c|d)) d++;
					for (int i = 0; i < 20; i++) raw16();
					*/
				void jsf16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(d);
					if (!(a | b) && !(c | d)) d++;
				}

				Uint32 simpleA::raw32() {
					enum { BARREL_SHIFT1 = 19 };
					Uint32 tmp = b ^ ((a << BARREL_SHIFT1) | (a >> (32 - BARREL_SHIFT1)));
					a = ~b + c;
					b = c;
					c += tmp;
					return b;
				}
				std::string simpleA::get_name() const { return "simpleA"; }
				void simpleA::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				Uint16 simpleB::raw16() {
					enum { BARREL_SHIFT1 = 3, BARREL_SHIFT2 = 5 };
					Uint16 tmp = ((a << BARREL_SHIFT1) | (a >> (16 - BARREL_SHIFT1))) ^ ~b;
					a = b + c;
					b = c;
					c = tmp + ((c << BARREL_SHIFT2) | (c >> (16 - BARREL_SHIFT2)));
					return tmp;
				}
				std::string simpleB::get_name() const { return "simpleB"; }
				void simpleB::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				Uint16 simpleC::raw16() {
					enum { BARREL_SHIFT1 = 3, BARREL_SHIFT2 = 5 };
					Uint16 tmp = ((a << BARREL_SHIFT1) | (a >> (16 - BARREL_SHIFT1))) ^ ~b;
					a = b + ((c << BARREL_SHIFT2) | (c >> (16 - BARREL_SHIFT2)));
					b = c ^ (c >> 2);
					c += tmp;
					return tmp;
				}
				std::string simpleC::get_name() const { return "simpleC"; }
				void simpleC::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				Uint32 simpleD::raw32() {
					enum { BARREL_SHIFT1 = 19 };
					Uint32 old = a;
					Uint32 tmp = b ^ ((a << BARREL_SHIFT1) | (a >> (32 - BARREL_SHIFT1)));
					a = b + c;
					b = c ^ old;
					c = old + tmp;
					return tmp;
				}
				std::string simpleD::get_name() const { return "simpleD"; }
				void simpleD::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				Uint32 simpleE::raw32() {
					Uint32 old = a + b;
					a = b ^ c;
					b = c + old;
					c = old + ((c << 13) | (c >> 19));
					return old;
				}
				std::string simpleE::get_name() const { return "simpleE"; }
				void simpleE::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				Uint16 simpleF::raw16() {
					Uint16 old = a ^ b;
					a = b ^ (c & d);
					b = c + old;
					c = ~d;
					d = old + ((d << 5) | (d >> 11));
					return c;
				}
				std::string simpleF::get_name() const { return "simpleF"; }
				void simpleF::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}

				Uint32 simpleG::raw32() {
					////Passes 32TB with no anomalies. WHY??? I dunno.
//					Uint32 t = a + b ^ d++;
//					a ^= rotate32(b, 13);
//					b += rotate32(c, 3);
//					c += rotate32(t, 25);
//					return t;

//					Uint32 t = (a + b ^ (d += 0x6C8E9CF5U));
					Uint32 t = (a + b ^ d++);
					a ^= rotate32(b, 13);
					b += c;
					c -= rotate32(t, 25); 
					return t;
					////original simpleG
//					Uint32 old = a ^ (b >> 7);
//					a = b + c + d;
//					b = c ^ d;
//					c = d + old;
//					d = old;
//					return a + c;
				}
				std::string simpleG::get_name() const { return "simpleG"; }
				void simpleG::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(d);
				}

				static Uint64 shift_array64(Uint64 vec[2], unsigned long bits) {
					bits -= 64;
					if (!(bits % 64)) return vec[bits / 64];
					return (vec[bits / 64] << (bits & 63)) | (vec[1 + bits / 64] >> (64 - (bits & 63)));
				}
				Uint32 trivium_weakenedA::raw32() {
					Uint32 tmp_a = Uint32(b >> 2) ^ Uint32(b >> 17);
					Uint32 tmp_b = Uint32(a >> 2) ^ Uint32(a >> 23);
					Uint32 new_a = tmp_a ^ Uint32(a >> 16) ^ (Uint32(b >> 16) & Uint32(b >> 18));
					Uint32 new_b = tmp_b ^ Uint32(b >> 13) ^ (Uint32(a >> 22) & Uint32(a >> 24));
					a <<= 32; a |= new_a;
					b <<= 32; b |= new_b;

					return tmp_a ^ tmp_b;
				}
				std::string trivium_weakenedA::get_name() const { return "trivium_weakenedA"; }
				void trivium_weakenedA::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
				}
				Uint16 trivium_weakenedB::raw16() {
					Uint16 tmp_a = Uint16(c >> 2) ^ Uint16(c >> 15);
					Uint16 tmp_b = Uint16(a >> 2) ^ Uint16(a >> 13);
					Uint16 tmp_c = Uint16(b >> 5) ^ Uint16(b >> 10);
					Uint16 new_a = tmp_a ^ Uint16(a >> 9) ^ (Uint16(c >> 14) & Uint16(c >> 16));
					Uint16 new_b = tmp_b ^ Uint16(b >> 7) ^ (Uint16(a >> 12) & Uint16(a >> 14));
					Uint16 new_c = tmp_c ^ Uint16(c >> 11) ^ (Uint16(b >> 9) & Uint16(b >> 11));
					a <<= 16; a |= new_a;
					b <<= 16; b |= new_b;
					c <<= 16; c |= new_c;

					return tmp_a ^ tmp_b ^ tmp_c;
				}
				std::string trivium_weakenedB::get_name() const { return "trivium_weakenedB"; }
				void trivium_weakenedB::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}

				Uint32 mo_Lesr32::raw32() {
					state = (state << 7) - state; state = rotate32(state, 23);
					return state;
				}
				std::string mo_Lesr32::get_name() const { return "mo_Lesr32"; }
				void mo_Lesr32::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint32 mo_ResrRers32::raw32() {
					a = rotate32(a, 21) - a; a = rotate32(a, 26);
					b = rotate32(b, 20) - rotate32(b, 9);
					return a ^ b;
				}
				std::string mo_ResrRers32::get_name() const { return "mo_ResrRers32"; }
				void mo_ResrRers32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
				}
				Uint32 mo_Rers32of64::raw32() {
					state = rotate64(state, 8) - rotate64(state, 29);
					return Uint32(state);
				}
				std::string mo_Rers32of64::get_name() const { return "mo_Rers32of64"; }
				void mo_Rers32of64::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint32 mo_Resr32of64::raw32() {
					state = rotate64(state, 21) - state; state = rotate64(state, 20);
					return Uint32(state);
				}
				std::string mo_Resr32of64::get_name() const { return "mo_Resr32of64"; }
				void mo_Resr32of64::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint32 mo_Resdra32of64::raw32() {
					state = rotate64(state, 42) - state;
					state += rotate64(state, 14);
					return Uint32(state);
				}
				std::string mo_Resdra32of64::get_name() const { return "mo_Resdra32of64"; }
				void mo_Resdra32of64::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint32 murmlacish::raw32() {
					Uint32 tmp = state1 + (state1 << 3);
					tmp ^= tmp >> 8;
					state1 = rotate32(state1, 11) + state2;
					state2 += state3 ^ (state3 >> 7) ^ tmp;
					state3 += tmp + (tmp << 3);
					return state1;
				}
				std::string murmlacish::get_name() const {
					std::ostringstream str;
					str << "murmlacish";
					return str.str();
				}
				void murmlacish::walk_state(StateWalkingObject *walker) {
					walker->handle(state1);
					walker->handle(state2);
					walker->handle(state3);
				}

				Uint16 gjishA::raw16() {
					b += a; c = rotate16(c, 4); a ^= b;
					c += b; a = rotate16(a, 7); b ^= c;
					a += c; b = rotate16(b, 11); c ^= a;
					return a;
				}
				std::string gjishA::get_name() const { return "gjishA"; }
				void gjishA::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				Uint16 gjishB::raw16() {
					b += a; c = rotate16(c, 4); a ^= b;
					c += b; a = rotate16(a, 7); b ^= c;
					c += counter++;
					a += c; b = rotate16(b, 11); c ^= a;
					return a;
				}
				std::string gjishB::get_name() const { return "gjishB"; }
				void gjishB::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(counter);
				}
				Uint32 gjishC::raw32() {
					b += a; c = rotate32(c, 21); a ^= b;
					c += b; a = rotate32(a, 13); b ^= c;
					a += c; b = rotate32(b, 0); c ^= a;
					return a;
				}
				std::string gjishC::get_name() const { return "gjishC"; }
				void gjishC::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				Uint32 gjishD::raw32() {
					//			4,10,17		4,10,19		5,7,11		5,7,12		5,7,13		5,7,14		5,7,15		5,8,10		5,8,11		5,8,12		5,8,13		5,8,14		5,8,15
					//--big		3			5			1			1			1			1			4			6			2			4			3			1			4
					//--huge	78						32			10			31			23									67						79			62

					//			6,8,15		6,8,17		6,10,13		6,11,20		6,11,18		6,11,15		6,11,14		6,11,19		6,15,29		6,9,21		6,9,17		
					//--big		1			4			7			10			6			2			4			5			7			7			10			
					//--huge	41															81																		

					//			7,9,13		7,9,14	(7,14,9)	7,9,15	(7,15,9)	7,9,16		7,10,15		7,10,14		7,11,17		
					//--big		3			1		1			1		3			7			4			4			3			
					//--huge	50			43		35			40		52												76			

					//			9,13,25		8,13,29		
					//--big		7			2			
					//--huge				65			

					//			21,13,0		21,13,3		21,13,5		19,11,6		18,11,5		23,13,7		23,13,6		21,13,7		12,19,7		12,19,5		14,21,6		
					//--big		~1			~9			~9			~5			6			4			7			1?			4			6			11
					//--huge	61																					82									
					b += a; c = rotate32(c, 5); a ^= b;
					c += b; a = rotate32(a, 8); b ^= c;
					a += c; b = rotate32(b, 16); c ^= a;
					return a;
				}
				std::string gjishD::get_name() const { return "gjishD"; }
				void gjishD::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}

				Uint16 ara16::raw16() {
					a += rotate16(b + c, 3);
					b += rotate16(c + a, 5);
					c += rotate16(a + b, 7);
					return a;
				}
				std::string ara16::get_name() const { return "ara16"; }
				void ara16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				Uint32 ara32::raw32() {
					a += rotate32(b + c, 7);
					b += rotate32(c + a, 11);
					c += rotate32(a + b, 15);
					return a;
				}
				std::string ara32::get_name() const { return "ara32"; }
				void ara32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				Uint16 arx16::raw16() {
					a ^= rotate16(b + c, 3);
					b ^= rotate16(c + a, 5);
					c ^= rotate16(a + b, 7);
					return a;
				}
				std::string arx16::get_name() const { return "arx16"; }
				void arx16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				Uint32 arx32::raw32() {
					a ^= rotate32(b + c, 7);
					b ^= rotate32(c + a, 11);
					c ^= rotate32(a + b, 15);
					return a;
				}
				std::string arx32::get_name() const { return "arx32"; }
				void arx32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				Uint16 hara16::raw16() {
					a += rotate16(b + c, 3);
					b ^= rotate16(c + a, 5);
					c += rotate16(a + b, 7);
					return a;
				}
				std::string hara16::get_name() const { return "hara16"; }
				void hara16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				Uint16 harx16::raw16() {
					a ^= rotate16(b + c, 3);
					b += rotate16(c + a, 5);
					c ^= rotate16(a + b, 7);
					return a;
				}
				std::string harx16::get_name() const { return "harx16"; }
				void harx16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}

				Uint16 learx16::raw16() {
					a ^= rotate16(b + c, 3);
					b ^= rotate16(c + (c << 3), 5);
					c ^= rotate16(a + (a << 3), 7);
					return a;
				}
				std::string learx16::get_name() const { return "learx16"; }
				void learx16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				Uint16 hlearx16::raw16() {
					a ^= rotate16(b + c, 3);
					b ^= rotate16(c + (c << 3), 5);
					c += rotate16(a + (a << 3), 7);
					return a;
				}
				std::string hlearx16::get_name() const { return "hlearx16"; }
				void hlearx16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}

				Uint16 alearx16::raw16() {
					a ^= rotate16(b + c, 3);
					b ^= rotate16(c + (c << 3), 5);
					c ^= rotate16(a + (a << 3), 7) + b;
					return a;
				}
				std::string alearx16::get_name() const { return "alearx16"; }
				void alearx16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}

				Uint16 arac16::raw16() {
					a += rotate16(b + c, 3) + counter++;
					b += rotate16(c + a, 5);
					c += rotate16(a + b, 7);
					return a;
				}
				std::string arac16::get_name() const { return "arac16"; }
				void arac16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(counter);
				}
				Uint16 arxc16::raw16() {
					a ^= rotate16(b + c, 3) + counter++;
					b ^= rotate16(c + a, 5);
					c ^= rotate16(a + b, 7);
					return a;
				}
				std::string arxc16::get_name() const { return "arxc16"; }
				void arxc16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(counter);
				}

				Uint16 rarns16::raw16() {					//		378		131		124		111		F31		F11		F12		F21		FF1		FE1		FD1		F1A
					Uint16 old1 = xs1, old2 = xs2, old3 = xs3;//							*				*						*		*				
					Uint16 old = xs1;
					xs1 = xs2 ^ (xs2 >> S1);
					xs2 = xs3 ^ (xs3 << S2);
					xs3 = xs3 ^ (xs3 >> S3) ^ old;
					Uint16 rv = old1 + rotate16(old1, 5);
					return (rv) ^ rotate16(xs1 + rv, 7);  //		16 TB?	512 GB	1+ TB	1 GB	128 GB	256 MB	4 GB	16 GB	64 MB	128 MB	2 GB	2 GB	27
					//return rv^rotate16(rv+xs1,9)^rotate16(rv+xs2,3);  //					64 GB			256 GB					1 TB	128 GB					27
					//return (rv ^ xs2) + rotate16(rv + xs1, 3);//			1+ TB			8 GB			16 GB					2 GB	16 GB					31
					//Uint16 rv = old1 + rotate16(old1, 3);
					//return (rv ^ xs2) + rotate16(xs1 + rv, 7);//							16 GB			16 GB					8 GB	1 GB					30
					//Uint16 rv = old2+old1 + rotate16(old1, 3);
					//return (rv ^ xs2) + rotate16(xs1 + rv, 7);//			4+ TB			64 GB	4+ TB	2+ TB					32 GB	16 GB					34
					//Uint16 rv = old1 + rotate16(old1+old2, 3);
					//return (rv ^ xs2) + rotate16(xs1 + rv, 7);//							32 GB			~1+ TB					~1+ TB	512 GB					35
				}
				std::string rarns16::get_name() const { return "rarns16"; }
				void rarns16::walk_state(StateWalkingObject *walker) {
					walker->handle(xs1);
					walker->handle(xs2);
					walker->handle(xs3);
					walker->handle(history);
					if (!(xs1 | xs2 | xs3)) xs1 = 1;
				}
				void rarns16::seed(Uint64 s) {
					xs1 = Uint16(s >> 0);
					xs2 = Uint16(s >> 16);
					xs3 = Uint16(s >> 32);
					if (!(xs1 | xs2 | xs3)) {
						xs1 = Uint16(s >> 48);
						xs2 = 1;
					}
					raw16();
					return;
				}
				Uint64 thrust64a::raw64() {
					const Uint64 s0 = state0;
					Uint64 s1 = state1;
					const Uint64 result = s0 + (state1 += 0x9E3779B97F4A7AF5ULL);
					s1 ^= s0;
					state0 = ((s0 << 55) | (s0 >> 9)) ^ s1 ^ ((s1 << 28) | (s1 >> 36));
					return result;
				}
				std::string thrust64a::get_name() const { return "thrust64a"; }
				void thrust64a::walk_state(StateWalkingObject *walker) {
					walker->handle(state0);
					walker->handle(state1);
				}



				Uint32 xoroshiro32_32::raw32() {
					// rotate32(state1 + 0x41C64E6D, 17) + 0x9E3779B9
					const uint32_t s0 = state0;
					const uint32_t s1 = s0 ^ state1;
					//<< 5 , rotl 3
					//<< 6 , rotl 4
					//<< 9 , rotl 5
					//<< 9 , rotl 7
					//<< 10, rotl 6
					//<< 11, rotl 5
					//<< 11, rotl 6
					//<< 11, rotl 9
					//<< 12, rotl 6
					//<< 12, rotl 8
					//<< 12, rotl 10
					//<< 15, rotl 11
					//<< 17, rotl 13
					//const Uint32 result = ((s0 << 6) - rotate32(s0, 4));// +0xE32BE5AB;// * 31;// +0x41C64E6D;// (s0 ^ 0xAEF1F2D9) + 0x41C64E6D;

					//Uint32 s2 = (state2 += UINT32_C(0x9E3779BD));
					//Uint32 s2 = (state2 = rotate32(state2 + 0xC68E9CB7U ^ 0xB5402ED7U, 1));

					// +UINT32_C(0x41C64E6D);// + (state2 += UINT32_C(0x9E3779BD));// +0xE32BE5AB;// * 31;// +0x41C64E6D;// (s0 ^ 0xAEF1F2D9) + 0x41C64E6D;
					//const Uint32 result = (s0 << 5) - rotate32(s0, 3) + UINT32_C(0x9E3779BD);
					//state0 = rotate32(s0, 13) ^ s1 ^ (s1 << 5); // a=13, b=5
					//state1 = rotate32(s1, 28); // c=28
					//return ((result << 10) - rotate32(result, 7));

					////// starfish32
					const uint32_t result = s0 * 31;
					state0 = rotate32(s0, 26) ^ s1 ^ (s1 << 9);
					state1 = rotate32(s1, 13);
					return (rotate32(result, 28)) + UINT32_C(0x9E3779BD);

					//const uint32_t result = (s0 << 5) - s0;
					//state0 = rotate32(s0, 26) ^ s1 ^ (s1 << 9);
					//state1 = rotate32(s1, 13);
					//return (result << 11) - rotate32(result, 4);


					//const uint32_t result = s0 * (s1 | 0xA529u);
					//state0 = rotate32(s0, 26) ^ s1 ^ (s1 << 9);
					//state1 = rotate32(s1, 13);
					//return rotate32(result, 21) + s1;

					//return ((result << 6) ^ (result << 28 | result >> 4)) + 0x9E3779BD;
					// (^ 0xAEF1F2D9) + 0x41C64E6D 23 + 0x9E3779B9 // 0xC74EAD55
					//with result=s0*31; passes 32TB no anomalies

					//return (rotate32(result, 28) * 31) + 0x9E3779BD;
					//return (rotate32(result, 17) ^ 0xC74EAD55) + 0x9E3779BD;
					// UINT32_C(0x9E3779BD);
					//s1 += s2;// (s2 ^ s2 >> 15);
					//return (s2 ^ s2 >> 11) + s0;
					//Uint32 s2 = (state2 = (state2 ^ UINT32_C(0x9E3779BD)) * UINT32_C(3));
					//Uint32 s2 = (state2 += UINT32_C(0x9E3779BD));
					//s2 ^= s2 >> 16;
					//s2 ^= s2 >> 16;
					//return (s2 ^ s2 >> 15) + s1;// rotate32(s1, 23);
					//0x632BE5AB
					//rotate32(s1, state2 >> 27)
					//s1 *= state2 | 1;
					//s1 ^= s1 >> 15;
					//s1 ^= state2;
					//s1 += state2;
					//s1 ^= rotate32(s1, 23) ^ rotate32(s1, 14);
					//c1 += rotate32(c0, 19);
					//s1 ^= rotate32(s1, 19) ^ rotate32(s1, 23);
					//return s1;// ^ s1 >> 12;
					//const Uint32 s0 = state0;
					//Uint32 s1 = state1;
					//const Uint32 result = s0 + s1;

					//s1 ^= s0;
					//state0 = rotate32(s0, 13) ^ s1 ^ (s1 << 5); // a, b
					//state1 = rotate32(s1, 28); // c
					//oriole modified
					//return (result ^ result >> 11) + (state2 += 0x632BE5AB);
					//lathe32
					//return rotate32(result, 10) + s0;
					//modified lathe32
					//return (result ^ result >> 10) + s0;

					//return (rotate32(result, 22)) ^ (s1 - s0);
					//return (rotate32(result, 23) ^ rotate32(result, 11) ^ result) + s0;

					//state0 = rotate32(s0, 13) ^ s1 ^ (s1 << 9); // a, b
					//state1 = rotate32(s1, 26); // c
					//return rotate32(result, 29) + (state2 += 0x632BE5C5);
					// 0xC74EAD5D; //0x632BE5AB passes 32TB

					// usually overkill, Weyl works fine with a rotate
					//state2 = (s2 ^ 0xC74EAD55) * 0x95CB3;
					//return result;
					// passes with a,b,c 13,9,26 up to 16TB
					//return (result ^ result >> 17) * 0xE5CB3;

					//return (result ^ result >> 18) * 0xA5CB3; // ^ rotate32(result, 10) ^ rotate32(result, 20)

					//return (result ^ result >> 19) * (s0 | 1); // ^ rotate32(result, 10) ^ rotate32(result, 20)
					//Uint32 tmp = state;
					//state ^= state << 13;
					//state ^= state >> 17;
					//state ^= state << 5;
					//tmp += (tmp << 18) | (tmp >> 14);//barrel shift and add
					//tmp += ((tmp << 7) | (tmp >> 25)) ^ ((tmp << 11) | (tmp >> 21));//two barrel shifts, an xor, and an add
					//return tmp ^ state;
				}
				std::string xoroshiro32_32::get_name() const { return "xoroshiro32_32"; }
				void xoroshiro32_32::walk_state(StateWalkingObject *walker) {
					walker->handle(state0);
					walker->handle(state1);
					if (state0 == 0 && state1 == 0)
						state0 = 1;
					walker->handle(state2);
					//state2 |= 1;
				}
				Uint16 xoroshiro16_16plus::raw16() {
					//Uint16 result = state0 + state1;
					//Uint16 tmp = state0 ^ state1;
					//state0 = ((state0 << 13) | (state0 >> (64 - 13))) ^ tmp ^ (tmp << 5);
					//state1 = ((tmp << 10) | (tmp >> (64 - 10)));
					//return result;
					const Uint16 s0 = state0;
					Uint16 s1 = state1;
					const Uint16 result = s0 + s1;

					s1 ^= s0;
					state0 = rotate16(s0, 13) ^ s1 ^ (s1 << 5); // a, b
					state1 = rotate16(s1, 10); // c
					return rotate16(result, 5) + s0 ^ (state1 - state0);

				}
				std::string xoroshiro16_16plus::get_name() const { return "xoroshiro16_16plus"; }
				void xoroshiro16_16plus::walk_state(StateWalkingObject *walker) {
					walker->handle(state0);
					walker->handle(state1);
					if (state0 == 0 && state1 == 0)
						state0 = 1;
				}

				uint64_t jump(uint64_t state){
				    const uint64_t poly = 0x5556837749D9A17FUL;
				    long val = 0UL, b = 1UL;
				    for (int i = 0; i < 63; i++, b <<= 1) {
				        if(poly & b) val ^= state;
				        state ^= state << 7;
				        state ^= state >> 9;
				    }
				    return val;
				}

				Uint64 oriole64::raw64() {
					//Uint64 result = state0 + state1;
					//Uint64 tmp = state0 ^ state1;
					//state0 = rotate64(state0, 55) ^ tmp ^ (tmp << 14);
					//state1 = rotate64(tmp, 36);
					//return result;
					// const Uint64 s0 = state0;
					// Uint64 s1 = state1;
					// const Uint64 result = s0 + s1;// +(state2 = state2 * Uint64(0x41C64E6D) + Uint64(0x9E3779B97F4A7C15));
					// s1 ^= s0;
					// state0 = rotate64(s0, 47) ^ s1 ^ (s1 << 16); // a, b
					// state1 = rotate64(s1, 37); // c
					// return rotate64(result, 9) + (state2 += Uint64(1));//0x6C8E9CF570932BD9
					// 												   //state0 = rotate64(s0, 24) ^ s1 ^ (s1 << 16); // a, b
					// 												   //state1 = rotate64(s1, 37); // c
					//// BlubberRNG, passes 64TB no anomalies.
					//// based on Konadare192Px++, roughly.
					//// Period is at least 2 to the 64 because state0 is a Weyl sequence.
					//// It could be as high as 2 to the 192, but is likely to be smaller and have sub-cycles.
					//const uint64_t a0 = (state0 += 0xC6BC279692B5C323UL);
					//const uint64_t b0 = state1;
					//const uint64_t c0 = state2;
					//state1 = rotate64(c0 + b0, 39) ^ a0;
					//state2 = rotate64(a0 ^ c0, 23) - b0;
					//return a0 - b0 ^ c0;

					//const uint64_t a0 = (state0 += 0xC6BC279692B5C323UL);
					//const uint64_t b0 = state1;
					//const uint64_t c0 = state2;
					//state1 = rotate64(c0 + b0, 39) ^ a0;
					//state2 = rotate64(a0 ^ c0, 23) - b0;
					//return c0;

					//// GrouchoRNG, passes 64TB no anomalies.
					//// also based on the Konadare family.
					//// Minimal data dependency; state0 does not depend on its previous value, nor does state1 or state2.
					//// Effectively an ARX algorithm, since the ~ it uses is equivalent to "^ 0xFFFFFFFFFFFFFFFFUL"
					//// and the "- n" is equivalent to "+ 1 + (n ^ 0xFFFFFFFFFFFFFFFFUL)"
					//// (Get it, Groucho ARX?)
					//const uint64_t a0 = state0;
					//const uint64_t b0 = state1;
					//const uint64_t c0 = state2;
					//state0 = b0 + ~c0;
					//state1 = rotate64(a0, 39) ^ c0;
					//state2 = rotate64(b0, 23) - a0;
					//return a0 - b0 ^ c0;

					//// heartbreaker. Gets a "very suspicious" right at 64TB:
					////  BCFN(2+2,13-0,T)                  R= +17.2  p =  9.4e-9   very suspicious
					//// this is significantly faster than GrouchoRNG, though, so some aspect of this may be useful.
					//const uint64_t a0 = state0;
					//const uint64_t b0 = state1;
					//const uint64_t c0 = state2;
					//state0 = b0 + ~c0;
					//state1 = rotate64(a0, 45) ^ c0;
					//state2 = rotate64(b0, 23) - a0;
					//return a0;

					//// HarpoRNG, passes 64TB with no anomalies.
					//// Like a pared-down version of GrouchoRNG, it is a tiny bit faster than that and maybe even faster
					//// than RomuTrio, at least when implemented in Java (where both are very quick, especially on Java 16).
					//// Like GrouchoRNG, it has minimal data dependency, and is an ARX-type generator. It does use the
					//// ~ and - operators, but those are constant-time and equivalent to add and xor operations.
					//const uint64_t a0 = state0;
					//const uint64_t b0 = state1;
					//const uint64_t c0 = state2;
					//state0 = b0 + ~c0;
					//state1 = rotate64(a0, 46) ^ c0;
					//state2 = rotate64(b0, 23) - a0;
					//return a0 + b0;

					//// ChicoRNG, passes 64TB with no anomalies.
					//// This is faster than HarpoRNG, due to returning a0 as-is.
					//// Each of the subsequent states can be computed in parallel instructions, without
					//// affecting the current return. This is also an ARX-type generator, though it uses a subtraction.
					//const uint64_t a0 = state0;
					//const uint64_t b0 = state1;
					//const uint64_t c0 = state2;
					//state0 = b0 ^ c0 + 0xC6BC279692B5C323UL;
					//state1 = rotate64(a0, 46) + c0;
					//state2 = rotate64(b0, 23) - a0;
					//return a0;

					//// "Chico2", passes 64TB with no anomalies.
					//// hwd results unknown so far; testing ongoing.
					//// (I was testing a different generator with hwd at the same time as this with PractRand).
					//// (I'm an old man, I'm confused!)
					//// Exhaustive searches for the pairs of rotation amounts with the best avalanche properties
					//// yielded these values, 27 and 42, as outliers with among the best quality.
					//const uint64_t a0 = state0;
					//const uint64_t b0 = state1;
					//const uint64_t c0 = state2;
					//state0 = b0 ^ c0 + 0xC6BC279692B5C323UL;
					//state1 = rotate64(a0, 27) + c0;
					//state2 = rotate64(b0, 42) ^ a0;
					//return b0;

					//// passes 64TB with 3 anomalies, all on the low 4 bits of a subsection.
					//// length= 8 gigabytes (2^33 bytes), time= 33.6 seconds
				  	//// Test Name                         Raw       Processed     Evaluation
				  	//// [Low4/16]DC6-9x1Bytes-1           R=  +6.7  p =  1.0e-3   unusual
					//// length= 256 gigabytes (2^38 bytes), time= 1026 seconds
					//// Test Name                         Raw       Processed     Evaluation
					//// [Low4/16]FPF-14+6/16:all          R=  -5.9  p =1-2.3e-5   unusual
					//// length= 2 terabytes (2^41 bytes), time= 8167 seconds
					//// Test Name                         Raw       Processed     Evaluation
					//// [Low4/32]BCFN(2+2,13-0,T)         R=  +9.1  p =  2.0e-4   unusual
					//// Attempting to pare down RomuTrio for speed.
					////
					//// Also note that this passes 2.5 PB of hwd testing without any issues.

					//const uint64_t fa = state0;
					//const uint64_t fb = state1;
					//const uint64_t fc = state2;
					//state0 = 0xD1342543DE82EF95UL * fc;
					//state1 = fa ^ fb ^ fc;
					//state2 = rotate64(fb, 42) + 0xC6BC279692B5C323UL;
					//return fa;

					//// Experimenting with better bijections.
					//// This currently fails BCFN very early.
					//const uint64_t a0 = state2;
					//const uint64_t b0 = state2 ^ state0;
					//const uint64_t c0 = state0 - state1;
					//state0 = 0xC6BC279692B5C323UL + a0;
					//state1 = rotate64(b0, 23);
					//state2 = rotate64(c0, 56);
					//return c0;

					//const uint64_t a0 = state0;
					//const uint64_t b0 = state1;
					//const uint64_t c0 = state2;
					//state0 = rotate64(b0, 20) ^ c0; 
					//state1 = rotate64(c0, 57) + a0;
					//state2 = b0 + 0xC6BC279692B5C323L;
					//return a0;

					//// MargeRNG, a slight tweak on an earlier generator that lets it pass 64TB with no anomalies.
					//// It's also quite fast, and doesn't enter a known 1-length cycle with an all-0 state.
					//// This was also tested with left-rotations of 42 and 21; 42 had an unusually high amount of anomalies,
					//// while 21 had serious issues in hwd testing at about 200 TB. Using 41 still passes 64TB with no anomalies,
					//// and is still undergoing hwd testing (but has passed 400TB without issue there).
					//const uint64_t fa = state0;
					//const uint64_t fb = state1;
					//const uint64_t fc = state2;
					//state0 = 0xD1342543DE82EF95UL * fc;
					//state1 = fa ^ fb ^ fc;
					//state2 = rotate64(fb, 41) + 0xC6BC279692B5C323UL;
					//return fa;

// + 0x9E3779B97F4A7C15UL;
// + 0xC6BC279692B5C323UL;
// + 0x71AF0AB5118A6E6DUL;

					//state0 = 0x71AF0AB5118A6E6DUL + b0;
					// state0 = 0xC6BC279692B5C323UL + a0 ^ b0;

//	const uint64_t fa = state0;
//	const uint64_t fb = state1;
//	const uint64_t fc = state2;
//	state0 = fb ^ fb << 7;
//	state1 = fa ^ fa >> 9;
//	state2 = fa + rotate64(fc, 41);
//	return fc + fb;
    // surprisingly, not bad... passes at least 32GB without anomalies, could do more?
	// has a fixpoint at all-zero.
//	const uint64_t fa = state0;
//	const uint64_t fb = state1;
//	const uint64_t fc = state2;
//	state0 = fc * 0xF1357AEA2E62A9C5UL;
//	state1 = rotate64(fa, 44);
//	state2 = fa + fb;
//	return fc;
	const uint64_t fa = state0;
	const uint64_t fb = state1;
	const uint64_t fc = state2;
	//// works well, passes 64TB with one anomaly at 32TB ("unusual")
//	state0 = fc * 0xF1357AEA2E62A9C5UL;
//	state1 = rotate64(fa, 44);
//	state2 = fb ^ fa + 0x9E3779B97F4A7C15UL;
//	return fc;
    // Duck192 (faster version)
	// Passes 64TB with no anomalies (seed 0).
	// With seed 1, has an "unusual" anomaly at 512MB, then passes the remaining 64TB without more anomalies.
	//  [Low4/32]DC6-9x1Bytes-1           R=  +7.6  p =  3.1e-4   unusual
	// Passes a huge amount, at least 2^57.32 bytes, of ReMort without suspect results.
	state0 = fc * 0xF1357AEA2E62A9C5UL;
	state1 = rotate64(fa, 44) + 0x9E3779B97F4A7C15UL;
	state2 = fb ^ fa;
	return state2;
				}
				std::string oriole64::get_name() const { return "duck192"; }
				void oriole64::walk_state(StateWalkingObject *walker) {
					// walker->handle(state0);
					// state0 |= (state0 == 0U);
					// state1 = jump(state0);
					// walker->handle(state2);
					walker->handle(state0);
					walker->handle(state1);
					// state0 |= ((state0 | state1) == 0);
					walker->handle(state2);
				}

				Uint64 xoshiro256starstar::raw64() {
//					//xoshiro256+
//					// const uint64_t result = state0 + state3;
//					// xoshiro256++
//					const uint64_t result = rotate64(state0 + state3, 23) + state0;
//					//xoshiro256**
//					// const uint64_t result = rotate64(state1 * 5, 7) * 9;
//					// used with the random rotate option below
//					//const uint64_t result = state0 + state1;
//
//					const uint64_t t = state1 << 17;
//
//					state2 ^= state0;
//					state3 ^= state1;
//					state1 ^= state2;
//					state0 ^= state3;
//
//					state2 ^= t;
//
//					state3 = rotate64(state3, 45);
//					return result;

// This has 64-bit output, but uses the 32-bit xoshiro128 state transition and the 32-bit ++ scrambler for its upper 32 bits.
// That scrambler uses states 0 and 3, while the lower 32 bits scramble states 2 and 1 similarly (using subtraction and a different rotation).
// It passes 64TB of PractRand without any anomalies. It passes ReMort testing through 2 to the 57.32 bytes. The ReMort results aren't especially
// strong, but were never considered suspect over the course of testing.
					const uint64_t result = (uint64_t)(rotate32(state0 + state3, 7) + state0) << 32 ^ (rotate32(state2 - state1, 13) + state2);
					const uint32_t t = state1 << 9;

					state2 ^= state0;
					state3 ^= state1;
					state1 ^= state2;
					state0 ^= state3;

					state2 ^= t;

					state3 = rotate32(state3, 11);
					return result;


					//return rotate64(result, state3 & 63);
					//					return result ^ result >> 28;
				}
				// std::string xoshiro256starstar::get_name() const { return "xoshiro256starstar"; }
				std::string xoshiro256starstar::get_name() const { return "xoshiro128weird"; }
				void xoshiro256starstar::walk_state(StateWalkingObject *walker) {
					walker->handle(state0);
					walker->handle(state1);
					walker->handle(state2);
					walker->handle(state3);
					state0 |= ((state0 | state1 | state2 | state3) == 0);
				}
				Uint64 xoshiro256scramjet::raw64() {
					//// This is an alternate scrambler for xoshiro256.
					//// It keeps the 4D equidistribution of xoshiro256**, but is closer to
					//// the mixers in LXM than the simple ** mixer.
					//// It passes 64TB without anomalies.
						//uint64_t result = (state[1] * 0xF1357AEA2E62A9C5UL); // inverse is 0x781494A55DAAED0DUL
						////// Either of the below two lines passes 64TB, no anomalies.
						//result ^= result >> 31;
						////result ^= result >> 44 ^ result >> 23;
						//result *= 0xF1357AEA2E62A9C5UL;
						//const uint64_t t = state[1] << 17;
						//state[2] ^= state[0];
						//state[3] ^= state[1];
						//state[1] ^= state[2];
						//state[0] ^= state[3];
						//state[2] ^= t;
						//state[3] = rotate64(state[3], 45);
						//return (result ^ result >> 29);
						
						//// Passes 64TB with one anomaly, "unusual" at 512GB:
						//// [Low4/64]DC6-9x1Bytes-1           R=  +5.3  p =  1.6e-3   unusual
//						uint64_t result = state[1] + 0xF1357AEA2E62A9C5UL; // inverse is 0x781494A55DAAED0DUL
//						result ^= result >> 31;
//						result += 0x9E3779B97F4A7C15UL;
//						result ^= result >> 29;
//						const uint64_t t = state[1] << 17;
//						state[2] ^= state[0];
//						state[3] ^= state[1];
//						state[1] ^= state[2];
//						state[0] ^= state[3];
//						state[2] ^= t;
//						state[3] = rotate64(state[3], 45);
//						return result;
						
						const uint64_t result = (state[1] + 0x9E3779B97F4A7C15UL) * 0xD1342543DE82EF95UL; // inverse is 0x781494A55DAAED0DUL
						
						const uint64_t t = state[1] << 17;
						state[2] ^= state[0];
						state[3] ^= state[1];
						state[1] ^= state[2];
						state[0] ^= state[3];
						state[2] ^= t;
						state[3] = rotate64(state[3], 45);
						return result ^ result >> 25 ^ result >> 47;
						
					// if (++current & 3)
					// {
					// 	const uint64_t result = state[current] + first[current];
					// 	return rotate64(current, 33) + second[current];
					// }
					// else
					// {
					// 	current = 0;
					// 	const uint64_t result = state[0] + 0x9E3779B97F4A7C15ULL;

					// 	const uint64_t t = state[1] << 17;

					// 	state[2] ^= state[0];
					// 	state[3] ^= state[1];
					// 	state[1] ^= state[2];
					// 	state[0] ^= state[3];

					// 	state[2] ^= t;

					// 	state[3] = rotate64(state[3], 45);
					// 	return rotate64(result, 33) + 0xDB4F0B9175AE2165ULL;
					// }
				}
				std::string xoshiro256scramjet::get_name() const { return "xoshiro256scramjet"; }
				void xoshiro256scramjet::walk_state(StateWalkingObject *walker) {
					walker->handle(state[0]);
					walker->handle(state[1]);
					walker->handle(state[2]);
					walker->handle(state[3]);
					state[0] |= ((state[0] | state[1] | state[2] | state[3]) == 0);
				}

				Uint32 xoshiro4x32::raw32() {
					//xoshiro256+
					//const uint64_t result = state0 + state3;
					//xoshiro256**
					//const uint32_t result = rotate32(state1 * 5, 7) * 9;
					//const uint32_t result = (state1 ^ UINT32_C(0x41C64E6D)) + UINT32_C(0x9E3779BD);
					// (state1 + UINT32_C(0x9E3779BD));
					// used with the random rotate option below
					//const uint64_t result = state0 + state1;
					// passes 32TB, gets only one anomaly at 32TB: "unusual" Gap-16:B
					// (^ 0xAEF1F2D9) (+ 0x41C64E6D) (rotate by 23) (+ 0x9E3779B9)
					// passes 32TB, one anomaly at 256GB: "unusual" [Low8/32]DC6-9x1Bytes-1
					// rotate32(state1 + 0x41C64E6D, 17) + 0x9E3779B9
					//const uint32_t result = (state1 ^ 0xAEF1F2D9) + 0x41C64E6D;
					//const uint32_t result = state1 + 0x9E3779B9;// 0x41C64E6D;//0xAEF1F2D9 // 0xC74EAD55 //rotate32(result, 23) + 0x9E3779B9;// 0x9E3779BDU;// 0x9E3779BB;
					// passes 32TB with only one anomaly, "unusual" at 4TB: FPF
					//const uint32_t result = (rotate32(state1 * 31U, 23) + UINT32_C(0x9E3779BD));
					// *3U;//(state0 << 8) + state0;// +UINT32_C(0x9E3779BD);
					// StarPhi32, passes 32TB with seed 0x6dcb1257
//					const uint32_t result = state1 * UINT32_C(31);
//					return (rotate32(result, 28)) + UINT32_C(0x9E3779BD);
//                  // Sashimi32, passes 32TB with seed 0; one unusual anomaly at 4TB (Low4/16 BCFN)
//					const uint32_t result = state1 + 0x9E3779BDu;
//					return (result << 7) - rotate32(result, 4);
//                  // Sushi32, passes 32TB with seed 0, no anomalies
//					const uint32_t result = state1 + 0x9E3779B9u;
//					return (result << 7) - rotate32(result, 3);


					const uint32_t result = state1 + 0x41C64E6Du;// + 0x9E3779B9u;
				    //const uint64_t result = (uint64_t)(rol32(s->s[0] + s->s[3], 7) + s->s[0]) << 32 ^ (rol32(s->s[2] - s->s[1], 13) + s->s[2]);
					// const uint64_t result = (uint64_t)(rotate32(state0 + state3, 7) + state0) << 32 ^ (rotate32(state2 - state1, 13) + state2);
					const uint32_t t = state1 << 9;

					state2 ^= state0;
					state3 ^= state1;
					state1 ^= state2;
					state0 ^= state3;

					state2 ^= t;

					state3 = rotate32(state3, 11);
					return result;
					//return rotate32(state3, 23) + (state0 ^ 0x41C64E6Du) * 0x9E3779BBu;
					//return (result << 7) - rotate32(result, 3);
					// return rotate32(result, 17) + 0x9E3779B9u;

					//return result ^ result >> 11;
					//return ((result << 11) - rotate32(result, 9));
					//return (result << 6) - rotate32(result, 4);
					//return rotate32(result, 21) + 0x9E3779B9;// *0xaaaaaaab;// *0xfbffdfff;// (result ^ result >> 17) * 3U;

					//return (result ^ result >> 15) * 3U;


					//return (result ^ rotate32(result, 6) ^ rotate32(result, 23)) * 3U;
					//return (rotate32(result, 23) * 3U);
					//return rotate64(result, state3 & 63);
					//					return result ^ result >> 28;
				}
				std::string xoshiro4x32::get_name() const { return "xoshiro4x32"; }
				void xoshiro4x32::walk_state(StateWalkingObject *walker) {
					walker->handle(state0);
					walker->handle(state1);
					walker->handle(state2);
					walker->handle(state3);
					state0 |= ((state0 | state1 | state2 | state3) == 0);
				}
				Uint32 mover32::raw32() {
					//a = rotate32(a * 0x9E377, 10);
					//b = rotate32(b * 0x64E6D, 22);
					// passes 32TB
					//a = rotate32(a * 0x9E37, 17);
					//b = rotate32(b * 0x4E6D, 14);
					// fails all over
					//a = rotate32(a + 0x9E3779B9, 2);
					//b = rotate32(b + 0x6C8E9CF7, 7);
					// bad at 4TB, very suspicious
					//a = rotate32(a * 0x9E37, 17);
					//b = rotate32(b + 0x9E3779B9, 2);

					//a = rotate32(a + 0x9E3779B9, 2);
					//b = (b ^ 0xC3564E95) * 0x9E37B;
					//return (a ^ a >> 12) + (b ^ rotate32(b, 11) ^ rotate32(b, 17));

					// works well up to at least 8TB
					//int z = (b += 0xC3564E95);
					//z = (z ^ z >> 15 ^ (a = rotate32(a + 0x9E3779B9, 2))) * 0x6C8E9;
					//return z ^ (z >> 15) + (a ^ a >> 12);

					//z = (z ^ z >> 14 ^ 0x6C8C9CF5) * 0x9E37B;
					//z = (z ^ z >> 13 ^ 0x867B8085) * 0x64E73;

					//Uint32 y = (a = rotate32(a + 0x9E3779B9, 2)), z = (b += 0xC3564E95);// (b = (b ^ 0xC3564E95) * 0x9E37B);
					//y ^= z ^ (z >> 15);
					//z ^= (y ^ y >> 12 ^ 0xC3564E95U);
					//return z ^ (z >> 14) + (y ^ y >> 13);

					//Uint32 z = (a ^ (a << 5) - 0x632BE5AD) + (b ^ 0xC3564E95 - (b << 8));
					//z ^= 0xC3564E95 - (z << 7);
					//return z ^ (z >> 12);
					//a = rotate32(a + 0x9E3779B9, 2);
					//b = rotate32(~(b * 0xACED), 28);
					//Uint32 z = (a ^ b ^ 0xC3564E95) * 0x9E37B;
					//return z ^ z >> 15;


					//int y = (a = rotate32(a + 0x9E3779B9, 2)) ^ 0x632BE5AD;
					//int x = (b = rotate32(b + 0x6C8E9CF7, 7)) ^ 0xC3564E95;
					//return rotate32(y, 11) ^ rotate32(x, 9) ^ (x + y >> 12);
					//int y = (a = rotate32(a + 0x9E3779B9, 2));// ^ 0x632BE5AD;
					//										  //y -= (b = rotate32(b ^ b << 5, 15));// ^ 0xC3564E95;
					//y ^= (b = rotate32(b + 0x6C8E9CF7, 7));
					//return (rotate32(y, 9) ^ rotate32(y, 17) ^ y);
					//return ((y + (y << 4)) ^ (y + (y << 16)) ^ y);

					//return (a = rotate32(a * 0xC6D5, 5)) ^ (b = rotate32(b * 0xA3A9, 9));

					// works very well, passes 32TB
					//return (a = rotate32(a * 0x89A7, 13)) ^ (b = rotate32(b * 0xBCFD, 17));
					//public class Mover32 extends java.util.Random{public int a=1, b=1; public Mover32(){} public Mover32(int seed){for(int i=(seed&0xFFFF);i>0;i--)a=Integer.rotateLeft(a*0x89A7,13);for(int i=seed>>>16;i>0;i--)b=Integer.rotateLeft(b*0xBCFD,17);} protected int next(int bits){return((a=Integer.rotateLeft(a*0x89A7,13))^(b=Integer.rotateLeft(b*0xBCFD,17)))>>> -bits;}}

					uint32_t r = ((a = rotate32(a + 0xAA78EDD7UL, 1)) ^ (b = rotate32(b + 0xC4DE9951UL, 25))) * 0x9E37BUL;
					return r ^ r >> 11 ^ r >> 21;

					// fails at 32GB, various 1/64 bit tests
					//return (a = (a ^ 0x632BE5AD) * 0x9E373) ^ (b = rotate32(b * 0xBCFD, 17));
				}
				std::string mover32::get_name() const { return "mover32"; }
				void mover32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					uint32_t r = a & 0xFFFF;
					a = 1;
					for (; r; r--)
					{
						a += UINT32_C(0xAA78EDD7);
						a = rotate32(a, 1);
					}
					walker->handle(b);
					r = b & 0xFFFF;
					b = 0;
					for (; r; r--)
					{
						b += UINT32_C(0xC4DE9951);
						b = rotate32(b, 25);
						//b *= UINT32_C(0xBCFD);
						//b = rotate32(b, 17);
					}
				}
				Uint64 mover64::raw64() {
					// this is Romu Duo
					uint64_t ap = a;
					//0xD3833E804F4C574Bu == 15241094284759029579u
					a = 0xD3833E804F4C574Bu * b;
					b = rotate64(b,36) + rotate64(b,15) - ap;
					return ap;
					
					//return (a = rotate64(a, 28) * UINT64_C(0x41C64E6B)) + (b = rotate64(b, 37) + UINT64_C(0x9E3779B97F4A7C15)); //0x9E3779B97F4A7C15 //  * UINT64_C(0x9E3779B9)
					
					//const uint64_t aa = a * UINT64_C(0x41C64E6B);
					//const uint64_t bb = b * UINT64_C(0x9E3779B9);
					//return (a = rotate64(aa, 42)) ^ (b = rotate64(bb, 23));
					
					// works, passes 32TB with 2 anomalies
					//return (a = rotate64(a, 26) * UINT64_C(0x41C64E6B)) ^ (b = rotate64(b, 37) + UINT64_C(0x9E3779B97F4A7C15));
				}

				std::string mover64::get_name() const { return "mover64"; }
				void mover64::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					if((a | b) == 0u) a = 1u;
					//walker->handle(a);
					//uint64_t r = a;
					//a = 1;
					//b = 1;
					//for (uint64_t ra = (r & 0xFFFF); ra; ra--)
					//{
					//	a *= UINT64_C(0x41C64E6B);
					//	a = rotate64(a, 42);
					//}
					//for (uint64_t rb = (r >> 16 & 0xFFFF); rb; rb--)
					//{
					//	b *= UINT64_C(0x9E3779B9);
					//	b = rotate64(b, 23);
					//}
				}

				Uint32 multimover32::raw32() {
					// passes 32TB
					//a = rotate32(a * 0x9E37, 17);
					//b = rotate32(b * 0x4E6D, 14);
					//a = rotate32(a * 0xACED, 28);
					//b = rotate32(b * 0xBA55, 19);
					//a *= 0x89A7;
					//b *= 0xBCFD;
					//c += 0xC0EF50EB;
					//const Uint32 p = a;
					//const Uint32 q = b;
					//const Uint32 r = c;
					//const Uint32 s = d;
					//a = rotate32(a, 13) * 0x89A7;
					//b = rotate32(b, 17) * 0xBCFD;
					//c = rotate32(c, 28) * 0xA01B;
					//d = rotate32(d, 16) * 0xC2B9;
					//a += UINT32_C(0x9E3779B9);
					//d = (a += UINT32_C(0x9E3779B9)) - rotate32(d, 9);
					//a = (a << 9) - rotate32(a, 5);
					//b = UINT32_C(0xC3E157C1) - rotate32(b, 19);
					//c = (c << 6) - rotate32(c, 4);

					//a = (a << 6) - rotate32(a, 4);

					//a = (a << 9) + rotate32(a, 8);
					//b = (b << 5) + rotate32(b, 1);
					//c = (c << 27) + rotate32(c, 20);
					//return a ^ b ^ c;
					////passes 32TB with one "unusual" anomaly at 8TB
					//// the subcycle generator (c = rotate32(c, 25) + 0xA5F152BFu) is a poor choice; the period of its subcycle generator
					//// shares a factor with the longer period of (a = rotate32(a, 7) + 0xC0EF50EBu). should try (c = rotate32(c, 30) + 0x9E3779B9u)
					//return (b = rotate32(b, 11) ^ (a = rotate32(a, 7) + 0xC0EF50EBu) ^ (d += (c = rotate32(c, 25) + 0xA5F152BFu)));

					//// passes to 32TB with no anomalies, but the period will be smaller due to the aforementioned shared factor.
					//return (b ^= b << 5 ^ b >> 11 ^ (a = rotate32(a, 7) + 0xC0EF50EBu) + (c = rotate32(c, 25) + 0xA5F152BFu));
					//// passes 32TB, I'm pretty sure? Speed isn't very good, probably doesn't warrant more checks.
					//return (b ^= b << 5 ^ b >> 11 ^ (a = rotate32(a, 7) + 0xC4DE9951u) + (c = rotate32(c, 1) + 0xAA78EDD7u));

					//return (a = rotate32(a, 23) * 0x402AB) ^ (b = rotate32(b, 28) * 0x01621) ^ (c = rotate32(c, 24) * 0x808E9) ^ (d = rotate32(d, 29) * 0x8012D);
					const uint32_t s = (a += 0xC1C64E6DU), t = (b += -((s | -s) >> 31) & 0x9E3779BBU),
					               u = (c += -((s | -s | t | -t) >> 31) & 0xC4DE9951U), v = (d += -((s | -s | t | -t | u | -u) >> 31) & 0x6C8E9CF7U);
				    uint32_t x = (s ^ s >> 17) * (t >> 12 | 1U) ^ v + u;
					x ^= x >> 16;
					x *= 0xAC451U;
					return x ^ x >> 15;

					//const uint32_t s = (stateA += 0xC1C64E6DU);
				    //uint32_t x = (s ^ s >> 17) * ((stateB += -((s | -s) >> 31) & 0x9E3779BBU) >> 12 | 1U);
					//x ^= x >> 16;
					//x *= 0xAC451U;
					//return x ^ x >> 15;



					//return (a = rotate32(a, 13)) + (b = rotate32(b, 17)) ^ (c = rotate32(c, 7));

					//return a ^ b;
					// fails all over
					//a = rotate32(a + 0x9E3779B9, 2);
					//b = rotate32(b + 0x6C8E9CF7, 7);
					//Uint32 y = 0;// (a = rotate32(a + 0x9E3779B9, 2));// +0xC3564E95;
					// +0x6C8E9CF7;
					// +0x632BE5AD;
					//Uint32 y = (b = rotate32(b - (b << 12) ^ 0xC68E9CB7, 21)) ^ (c = rotate32(c + 0xC0EF50EB, 7));

					//Uint32 y = (b = rotate32(b + 0xC68E9CB7 ^ 0xB2B386E5, 24)) ^ (c = rotate32(c - 0x9E3779B9 ^ 0xE541440F, 22));
					//Uint32 y = (a = rotate32(a - 0x9E3779B9 ^ 0xE541440F, 22)) ^ (c = rotate32(c + 0xC0EF50EB, 7));
					//y ^= y << 9;
					//return y ^ y >> 13;
					//return ((y + (y << 4)) ^ (y + (y << 16)) ^ y);
				}
				std::string multimover32::get_name() const { return "multimover32"; }
				void multimover32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(d);
					//if (a == 0) a = 1;
					//if (b == 0) b = 1;
					//if (c == 0) c = 1;
					//if (d == 0) d = 1;

					//walker->handle(d);
					//uint32_t r;
					//r = a & 0xFFFF;
					//a = 1;
					//for (; r; r--)
					//{
					//	//a = rotate32(a, 13) * UINT32_C(0x89A7);
					//	a = (a << 9) - rotate32(a, 5);
					//	//d = (a += UINT32_C(0x9E3779B9)) - rotate32(d, 9);
					//}
					//walker->handle(b);
					//r = b & 0xFFFF;
					//b = 1;
					//for (; r; r--)
					//{
					//	//b = rotate32(b, 17) * UINT32_C(0xBCFD);
					//	b = UINT32_C(0xC3E157C1) - rotate32(b, 19);
					//}
					//walker->handle(c);
					//r = c & 0xFFFF;
					//c = 1;
					//for (; r; r--)
					//{
					//	//c = rotate32(c, 28) * UINT32_C(0xA01B);
					//	c = (c << 6) - rotate32(c, 4);
					//}
					////r = d & 0xFFFF;
					////d = 1;
					////for (; r; r--)
					////{
					////	//d = rotate32(d, 16) * UINT32_C(0xC2B9);
					////	d = (d << 6) - rotate32(d, 4);
					////}
				}

				
				// Cloud values (shift constants).
				uint8_t clouds[8] = {3, 6, 1, 5, 3, 7, 3, 5};
				// Tested values, all good.
				//{1, 2, 6, 5, 4, 1, 6, 1};
				//{1, 2, 3, 4, 5, 1, 2, 3};
				//{1, 1, 1, 1, 1, 1, 1, 1};
				//{3, 1, 3, 1, 3, 1, 3, 1};
				//{3, 7, 5, 1, 4, 3, 5, 1};
				//{7, 3, 6, 1, 3, 1, 3, 1};
				//{1, 3, 5, 2, 1, 5, 7, 3};
				//{3, 6, 1, 5, 3, 7, 3, 5};
				uint8_t masks[8] = {
				    0, 1, 3, 7, 15, 31, 63, 127
				};
				uint32_t cloud::raw32()
				{
					int line = 0;
				    int nbytes = 0;
				    uint8_t bytes[4] = {0, 0, 0, 0};
				    while (nbytes < 4) {
				        if (line == 7) {
				            bytes[nbytes] = result();
				            nbytes++;
						}
					    board[line] ^= board[((line + 1) & 7)];
					    board[((line + 7) & 7)] ^= (board[line] & masks[clouds[line]]) << (8 - clouds[line]);
				   		board[line] >>= clouds[line];

				        line = ((line + 1) & 7);
				    }
				    return (bytes[0] << 24) | (bytes[1] << 16) | (bytes[2] << 8) | (bytes[3]);
				}

				uint8_t cloud::result()
				{
				    uint8_t total = 0;
				    for (int x = 0; x < 8; ++x) {
				        total |= ((board[x] & 1U) << x);
				    }
				    return total;
				}
				std::string cloud::get_name() const { return "cloud"; }
				void cloud::walk_state(StateWalkingObject *walker) {
					walker->handle(board[0]);
					walker->handle(board[1]);
					walker->handle(board[2]);
					walker->handle(board[3]);
					walker->handle(board[4]);
					walker->handle(board[5]);
					walker->handle(board[6]);
					walker->handle(board[7]);
				}

				Uint64 acorn64_10::raw64() {
				    return (j += i += h += g += f += e += d += c += b += a += stream);
//					t = (t ^ t >> 31) * (stream += 0x9E3779B97F4A7C16UL);
//					return t ^ t >> 26 ^ t >> 6;
				}
				std::string acorn64_10::get_name() const { return "acorn64_10"; }
				void acorn64_10::walk_state(StateWalkingObject *walker) {
					walker->handle(stream);
					stream |= 1U;
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(d);
					walker->handle(e);
					walker->handle(f);
					walker->handle(g);
					walker->handle(h);
					walker->handle(i);
					walker->handle(j);
				}

// Chonk8x32
// Passes 64TB of PractRand without anomalies.
// Period is at minimum 2 to the 64, max is much higher (a multiple of 2 to the 64).
				Uint32 chonk8x32::raw32() {
					const Uint32 fa = a;
					const Uint32 fb = b;
					const Uint32 fc = c;
					const Uint32 fd = d;
					const Uint32 fe = e;
					const Uint32 ff = f;
					const Uint32 fg = g;
					const Uint32 fh = h;
					a = fa + 0x9e3779b9;
					b = fa ^ ff;
					c = fb + fd;
					d = rotate32(fc, 25);
					e = fb - fc;
					f = fe ^ fh;
					g = __builtin_ctzll(fa);
					h = fh + fg;
					return ff;
				}
				std::string chonk8x32::get_name() const { return "chonk8x32"; }
				void chonk8x32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
					walker->handle(d);
					walker->handle(e);
					walker->handle(f);
					g = __builtin_ctzll(a);
					walker->handle(h);
				}
				Uint64 spangled_varqual::raw64() {
//					// SpangledRandom
//					// Passes 64TB with no anomalies when using rounds=4
//					// (Really that means 7 rounds total, since there's a hardcocoded round before and two after)
//					// Period is 2 to the 64; has 2 to the 64 different streams, plus a variable amount of "keys" that may be used
//					// Here, they keys are equal to the current round iteration, starting at 1
//					// This is an ARX generator with constant-time skipping through the state and constant-time stream switching
//					// It uses the round function from the Speck cipher
//					// If all streams are appended to one another, the resulting generator would have a period of 2 to the 128, 1D-equidistributed
//					// You can do that by changing stateB's assignment to:
//					// Uint64 b = (stateB += __builtin_ctzll(a));
//					// This eliminates the fast skipping, makes this have one stream, and makes it no longer an ARX generator
//					// It may also slow it down somewhat
//					Uint64 a = (stateA += 0x9E3779B97F4A7C15UL);
//					Uint64 b = (stateB += 0xD1B54A32D192ED03UL);
//					b = (rotate64(b, 56) + a ^ 0xA62B82F58DB8A985UL); a = (rotate64(a, 3) ^ b);
//					for (int i = 1; i <= rounds; i++) {
//						b = (rotate64(b, 56) + a ^ i);
//						a = (rotate64(a, 3) ^ b);
//					}
//					b = (rotate64(b, 56) + a ^ 0xE35E156A2314DCDAL); a = (rotate64(a, 3) ^ b);
//					return (rotate64(a, 3) ^ (rotate64(b, 56) + a ^ 0xBEA225F9EB34556DL));

//// The next variant passes ReMort to 179 PB, somehow.
//// It doesn't get any results marked suspect, but it isn't in the strongest position at the end.
//   6	                   0  - 2.29283e-01	                   0  - 2.29283e-01	                   1  + 2.59071e+00	                   0  - 2.29283e-01	                   0  - 2.29283e-01	
//   7	                4366  - 5.46697e-01	                4482  + 1.01280e+00	                4329  - 1.68021e+00	                4324  - 1.88095e+00	                4424  + 1.78206e-02	
//   8	              360712  - 5.52733e-02	              360870  + 7.79479e-04	              361736  + 2.15956e+00	              359955  - 2.23585e+00	              360331  - 7.55772e-01	
//   9	            22346785  - 4.61707e+00	            22346511  - 4.86946e+00	            22361989  + 1.13803e+00	            22363776  + 2.08722e+00	            22347112  - 4.32465e+00	
//  11	                4533  + 3.14677e+00	                4433  + 7.23294e-02	                4445  + 2.02084e-01	                4354  - 8.46374e-01	                4511  + 2.08173e+00	
//  12	            85008292  - 1.32253e+00	            85017509  - 2.26202e-02	            85022480  + 1.51104e-01	            85017423  - 2.55128e-02	            85018395  - 2.94966e-03	
//  13	          6948700313  + 3.70374e-02	          6948681618  - 1.01254e-03	          6948710573  + 9.95614e-02	          6948709724  + 9.32378e-02	          6948599192  - 1.04169e+00	
//  14	        430509784628  - 4.55994e+00	        430509794206  - 4.49781e+00	        430510189969  - 2.30320e+00	        430510195966  - 2.27554e+00	        430509877592  - 3.97491e+00	
//  16	              361618  + 1.62081e+00	              360551  - 2.53128e-01	              360281  - 9.07421e-01	              361157  + 2.55719e-01	              361127  + 2.07704e-01	
//  17	          6948743354  + 5.02377e-01	          6948704209  + 5.72113e-02	          6948671245  - 2.44167e-02	          6948695568  + 1.83680e-02	          6948741199  + 4.66398e-01	
//  18	        567923660173  + 2.10298e-01	        567922783250  - 4.97098e-01	        567923502019  + 6.18618e-02	        567923723078  + 2.93824e-01	        567922498705  - 1.17209e+00	
//  19	      35186128925378  - 8.93609e-01	      35186129842513  - 6.25200e-01	      35186133564508  - 2.66440e-02	      35186133318250  - 4.19204e-02	      35186130101405  - 5.58085e-01	
//  21	            22361897  + 1.09690e+00	            22346883  - 4.52843e+00	            22347140  - 4.30005e+00	            22362537  + 1.39874e+00	            22362410  + 1.33593e+00	
//  22	        430510168797  - 2.40218e+00	        430509773490  - 4.63272e+00	        430509801636  - 4.44990e+00	        430510207494  - 2.22284e+00	        430510160791  - 2.44015e+00	
//  23	      35186133416420  - 3.54174e-02	      35186129876698  - 6.16119e-01	      35186129128108  - 8.30162e-01	      35186133344861  - 4.01035e-02	      35186134679390  + 6.11102e-04	
//  24	    2179984579837982  + 3.50106e-02	    2179984583788025  + 7.38274e-02	    2179984579654789  + 3.35577e-02	    2179984575016781  + 7.03115e-03	    2179984578568664  + 2.55761e-02	
//  
//                21.082 15 =>   1-9.942362e-02           21.761 15 =>   1-8.382222e-02           18.368 15 =>     8.081723e-01           13.723 15 =>     5.301854e-01           18.406 15 =>     8.110533e-01
//  
//  --- Finished -- ReMort trials: 2251799813685248	2^51.0 (out of 2^51) trials -- spangled128 -- 64 bits: 	ps: 1-9.942e-02  1-8.382e-02*   8.082e-01    5.302e-01    8.111e-01   => p <   3.05979e-02   2^54.32 calls, 2^57.32 bytes	2^36.34 bytes/second	used:  23::20:35:12.35
					// Also passes 64TB with no anomalies when this variant uses rounds=2
					// So really 5 rounds, with 3 of the states hardcoded.
					// This adds a rotation of b to a before all rounds except the first.
					// Here, a left rotation of 41 is the only amount used in the new code.
					// Uint64 a = (stateA += 0x9E3779B97F4A7C15UL);
					// Uint64 b = (stateB += 0xD1B54A32D192ED03UL);
					// b = (rotate64(b, 56) + a ^ 0xA62B82F58DB8A985UL); a = (rotate64(a, 3) ^ b);
					// for (int i = 1; i <= rounds; i++) {
					// 	a += rotate64(b, 41);
					// 	b = (rotate64(b, 56) + a ^ i);
					// 	a = (rotate64(a, 3) ^ b);
					// }
					// a += rotate64(b, 41);
					// b = (rotate64(b, 56) + a ^ 0xE35E156A2314DCDAL); a = (rotate64(a, 3) ^ b);
					// a += rotate64(b, 41);
					// return (rotate64(a, 3) ^ (rotate64(b, 56) + a ^ 0xBEA225F9EB34556DL));

					// This passes 64TB with rounds=4, which really does mean 4 rounds now.
					// It has a third state, which it uses like the key in Speck, and changes per-round and per-result.
					// A variant is tested many lines below this that does not change the third state per-round, but does per-result.
					// Uint64 a = (stateA += 0x9E3779B97F4A7C15UL);
					// Uint64 b = (stateB += 0xD1B54A32D192ED03UL);
					// Uint64 c = (stateC += 0xDE916ABCC965815BUL);
					// for (int i = 1; i < rounds; i++) {
					// 	// c += (rotate64(a, 41)) ^ i;
					// 	b = (rotate64(b, 56) + a ^ (c += 0xBEA225F9EB34556DUL));
					// 	a = (rotate64(a, 3) ^ b);
					// }
					// // c += (rotate64(a, 41));
					// b = (rotate64(b, 56) + a ^ c + 0xF1357AEA2E62A9C5UL);
					// a = (rotate64(a, 3) ^ b);
					// return a;

					// MarshRandom
					// Passes 64TB of PractRand; not terribly fast.
					// Period is 2 to the 64. Has 2 to the 128 streams.
					// Uint64 a = (stateA += 0xDE916ABCC965815BUL);
					// Uint64 b = (stateB += 0xF1357AEA2E62A9C5UL);
					// Uint64 c = (stateC += 0xBEA225F9EB34556DUL);
					// a = (rotate64(a, 3) ^ (rotate64(b, 56) + a ^ c));
					// a ^= a >> 27;
					// a *= 0x3C79AC492BA7B653UL;
					// a ^= a >> 33;
					// a *= 0x1C69B3F74AC4AE35UL;
					// a ^= a >> 27;
					// return a;

//					Uint64 a = (stateA += 0xDE916ABCC965815BUL);
//					Uint64 b = (stateB += 0xF1357AEA2E62A9C5UL);
//					Uint64 c = (stateC += 0xBEA225F9EB34556DUL);
					
					// DraculaRandom
					// Passes 64TB of PractRand with no anomalies.
					// Has a period of 2 to the 192.
					// Needs at least two rounds of the Speck cipher; fails with one.
//					Uint64 a = (stateA = stateA * 0xD1342543DE82EF95UL + 0x9E3779B97F4A7C15UL);
//					Uint64 b = (stateB = stateB * 0x2C6FE96EE78B6955UL + __builtin_ctzll(a));
//					Uint64 c = (stateC = stateC * 0x369DEA0F31A53F85UL + __builtin_ctzll(a&b));
//					b = rotate64(b, 56) + a ^ c;
//					a = (rotate64(a, 3) ^ b);
//					return (rotate64(a, 3) ^ (rotate64(b, 56) + a ^ c));

					
					// Respect192
					// 3 uint64_t for state, period of 2 to the 192.
					// Passes 64TB of PractRand with no anomalies.
					// Runs 4 rounds of the Speck cipher on its three states as a, b, and c, returning a.
					// This one might have issues with its increments at extremely long test lengths, for complicated reasons.
					// At least with uint32_t states, XORing the large constant and the ctz result works better than adding.
					// Uint64 a = (stateA += 0x9E3779B97F4A7C15UL);
					// Uint64 b = (stateB += 0xD1342543DE82EF95UL + __builtin_ctzll(a));
					// Uint64 c = (stateC += 0xA62B82F58DB8A985UL + __builtin_ctzll(a&b));
					// b = rotate64(b, 56) + a ^ c;
					// a = (rotate64(a, 3) ^ b);
					// b = rotate64(b, 56) + a ^ c;
					// a = (rotate64(a, 3) ^ b);
					// b = rotate64(b, 56) + a ^ c;
					// a = (rotate64(a, 3) ^ b);
					// return (rotate64(a, 3) ^ (rotate64(b, 56) + a ^ c));

					// // InsectRandom
					// // Passes 64TB of PractRand with one anomaly at 64TB:
					// // [Low1/64]TMFn(2+1):wl             R= +20.1  p~=   6e-6    unusual
					// Uint64 a = (stateA += 0xBEA225F9EB34556DUL);
					// Uint64 b = (stateB += 0xD1342543DE82EF95UL);// ^ __builtin_ctzll(a));
					// Uint64 c = (stateC += 0xA62B82F58DB8A985UL);// ^ __builtin_ctzll(a&b));
					// b = rotate64(b, 56) + a ^ c;
					// a = (rotate64(a, 3) ^ b);
					// b = rotate64(b, 56) + a ^ c;
					// a = (rotate64(a, 3) ^ b);
					// b = rotate64(b, 56) + a ^ c;
					// a = (rotate64(a, 3) ^ b);
					// return (rotate64(a, 3) ^ (rotate64(b, 56) + a ^ c));

					// RespectfulRandom
					// Passes 64TB of PractRand with no anomalies.
					// Period is 2 to the 64; permits O(1) skip to any point; 2 to the 128 streams.
					Uint64 a = (stateA += 0xBEA225F9EB34556DUL);
					Uint64 b = (stateB += 0xD1342543DE82EF95UL);
					Uint64 c = (stateC += 0xA62B82F58DB8A985UL);
					b = rotate64(b, 56) + a ^ c;
					a = (rotate64(a, 3) ^ b);
					b = rotate64(b, 51) + a ^ c;
					a = (rotate64(a, 5) ^ b);
					b = rotate64(b, 43) + a ^ c;
					a = (rotate64(a, 8) ^ b);
					return (rotate64(a, 13) ^ (rotate64(b, 30) + a ^ c));
				}
				std::string spangled_varqual::get_name() const {
					std::ostringstream str;
					str << "spangled" << rounds;
					return str.str();
				}
				void spangled_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(stateA);
					walker->handle(stateB);
					//stateC = 1UL;//stateA ^ stateB + 0xA62B82F58DB8A985UL;
					walker->handle(stateC);
				}

			}
		}
	}
}
