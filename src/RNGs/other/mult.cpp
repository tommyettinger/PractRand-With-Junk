#include <string>
#include <sstream>
#include <cstdint>
#include "PractRand/config.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include "PractRand/rng_internals.h"
#if defined _MSC_VER && _MSC_VER >= 1800
#include <intrin.h>
#if defined(_M_X64) && !defined(_M_ARM64EC)
#pragma intrinsic(_umul128)
#endif
#endif

#include "PractRand/RNGs/other/mult.h"
//#include "PractRand/test_helpers.h"

namespace PractRand {
	using namespace Internals;
	namespace RNGs {
		namespace Polymorphic {
			namespace NotRecommended {
				Uint32 lcg32of64_varqual::raw32() {
					state = state * 0x5DEECE66D + 11;
					//state = state * 1103515245 + 12345; 0x5DEECE66DL
					return Uint32(state >> outshift);
				}
				std::string lcg32of64_varqual::get_name() const {
					std::ostringstream str;
					str << "lcg(32," << (32 + outshift) << ")";
					return str.str();
				}
				void lcg32of64_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint16 lcg16of64_varqual::raw16() {
					state = state * 1103515245 + 12345;
					return Uint16(state >> outshift);
				}
				std::string lcg16of64_varqual::get_name() const {
					std::ostringstream str;
					str << "lcg(16," << (16 + outshift) << ")";
					return str.str();
				}
				void lcg16of64_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint8 lcg8of64_varqual::raw8() {
					state = state * 1103515245 + 12345;
					return Uint8(state >> outshift);
				}
				std::string lcg8of64_varqual::get_name() const {
					std::ostringstream str;
					str << "lcg(8," << (8 + outshift) << ")";
					return str.str();
				}
				void lcg8of64_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}

				Uint32 lcg32of128_varqual::raw32() {
					//large multiplication is much harder to write in C (than asm)
					//made some compromises here... restricting the multiplier to 32 bits
					//which would hurt quality some, being so small compared to the state
					//but I think the correct comparison is to the output window, not the state
					//which is only 16 bits, so it should all be good
					const Uint32 multiplier = 1103515245;
					const Uint64 adder = 1234567;
					Uint64 a = Uint32(low) * Uint64(multiplier);
					Uint64 b = (low >> 32) * Uint64(multiplier);
					Uint64 old = low;
					low = a + (b << 32);
					b += a >> 32;
					high *= multiplier;
					high += b >> 32;
					high += old;//adds 2**64 to the multiplier
					low += adder;
					if (a < adder) high++;
					if (outshift >= 64) return Uint32(high >> (outshift-64));
					if (outshift > 32) return Uint32( (low >> outshift) | (high << (64-outshift)) );
					return Uint32(low >> outshift);
				}
				std::string lcg32of128_varqual::get_name() const {
					std::ostringstream str;
					str << "lcg(32," << (32 + outshift) << ")";
					return str.str();
				}
				void lcg32of128_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(low);
					walker->handle(high);
				}
				Uint16 lcg16of128_varqual::raw16() {
					const Uint32 multiplier = 1103515245;
					const Uint64 adder = 1234567;
					Uint64 a = Uint32(low) * Uint64(multiplier);
					Uint64 b = (low >> 32) * Uint64(multiplier);
					Uint64 old = low;
					low = a + (b << 32);
					b += a >> 32;
					high *= multiplier;
					high += b >> 32;
					high += old;//adds 2**64 to the multiplier
					low += adder;
					if (a < adder) high++;
					if (outshift >= 64) return Uint16(high >> (outshift-64));
					if (outshift > 48) return Uint16( (low >> outshift) | (high << (64-outshift)) );
					return Uint16(low >> outshift);
				}
				std::string lcg16of128_varqual::get_name() const {
					std::ostringstream str;
					str << "lcg(16," << (16 + outshift) << ")";
					return str.str();
				}
				void lcg16of128_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(low);
					walker->handle(high);
				}
				Uint8 lcg8of128_varqual::raw8() {
					const Uint32 multiplier = 1103515245;
					const Uint64 adder = 1234567;
					Uint64 a = Uint32(low) * Uint64(multiplier);
					Uint64 b = (low >> 32) * Uint64(multiplier);
					Uint64 old = low;
					low = a + (b << 32);
					b += a >> 32;
					high *= multiplier;
					high += b >> 32;
					high += old;//adds 2**64 to the multiplier
					low += adder;
					if (a < adder) high++;
					if (outshift >= 64) return Uint8(high >> (outshift-64));
					if (outshift > 56) return Uint8( (low >> outshift) | (high << (64-outshift)) );
					return Uint8(low >> outshift);
				}
				std::string lcg8of128_varqual::get_name() const {
					std::ostringstream str;
					str << "lcg(8," << (8 + outshift) << ")";
					return str.str();
				}
				void lcg8of128_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(low);
					walker->handle(high);
				}


				//similar to lcg16of32, but with a longer period
				Uint16 lcg16of32_extended::raw16() {
					state = state * 1103515245 + add;
					if (!(state & 0x7fff)) add += 2;
					return Uint16(state >> 16);
				}
				std::string lcg16of32_extended::get_name() const {return "lcg16of32_extended";}
				void lcg16of32_extended::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
					walker->handle(add);
					add |= 1;
				}

				//similar to lcg32, but with a longer period
				Uint32 lcg32_extended::raw32() {
					state = state * 1103515245 + add;
					if (!(state & 0x7fff)) add += 2;
					return state;
				}
				std::string lcg32_extended::get_name() const {return "lcg32_extended";}
				void lcg32_extended::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
					walker->handle(add);
					add |= 1;
				}

				Uint32 clcg32of96_varqual::raw32() {
					lcg1 = lcg1 * 1103515245 + 12345;
					Uint64 tmp = Uint64(lcg2) * 1579544716;
					lcg2 = Uint32(tmp & 0x7FffFFff) + Uint32(tmp >> 33) + 1;
					return lcg2 + Uint32(lcg1 >> outshift);
				}
				std::string clcg32of96_varqual::get_name() const {
					std::ostringstream str;
					str << "clcg(32," << (32 + outshift + 32) << ")";
					return str.str();
				}
				void clcg32of96_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg1); walker->handle(lcg2);
					if (!lcg2) lcg2 = 1;
				}
				Uint16 clcg16of96_varqual::raw16() {
					lcg1 = lcg1 * 1103515245 + 12345;
					Uint64 tmp = Uint64(lcg2) * 1579544716;
					lcg2 = Uint32(tmp & 0x7FffFFff) + Uint32(tmp >> 33) + 1;
					return Uint16(lcg2 >> 12) + Uint16(lcg1 >> outshift);
					//return Uint16((lcg2 + lcg1) >> outshift);
				}
				std::string clcg16of96_varqual::get_name() const {
					std::ostringstream str;
					str << "clcg(16," << (16 + outshift + 32) << ")";
					return str.str();
				}
				void clcg16of96_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg1); walker->handle(lcg2);
					if (!lcg2) lcg2 = 1;
				}
				Uint8 clcg8of96_varqual::raw8() {
					lcg1 = lcg1 * 1103515245 + 12345;
					Uint64 tmp = lcg2 * Uint64(1579544716);
					lcg2 = Uint32(tmp & 0x7FffFFff) + Uint32(tmp >> 33) + 1;
					return Uint8(lcg2) + Uint8(lcg1 >> outshift);
				}
				std::string clcg8of96_varqual::get_name() const {
					std::ostringstream str;
					str << "clcg(8," << (8 + outshift + 32) << ")";
					return str.str();
				}
				void clcg8of96_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg1); walker->handle(lcg2);
					if (!lcg2) lcg2 = 1;
				}

				void pcg32::seed(Uint64 s) { state = 0; raw32(); state += s; raw32(); }
				Uint32 pcg32::raw32() {
					Uint64 oldstate = state;
					state = state * 0x5851f42d4c957f2dULL + inc;
					Uint32 xorshifted = Uint32(((oldstate >> 18u) ^ oldstate) >> 27u);
					Uint32 rot = Uint32(oldstate >> 59u);
					return (xorshifted >> rot) | (xorshifted << ((~rot + 1u) & 31));
				}
				std::string pcg32::get_name() const {
					if (inc == 0xda3e39cb94b95bdbULL) return "pcg32";
					std::ostringstream str;
					str << "pcg32(" << std::hex << inc << ")";
					return str.str();
				}
				void pcg32::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}

				void pcg64::seed(Uint64 s) { state = 0; raw64(); state += s; raw64(); }
				Uint64 pcg64::raw64() {
					//// Passes at least 32TB with no anomalies.
					//// This has probably been tested way more thoroughly than I can do; I'm sure it's fine.
					//// ... Well, I say that, but if you test this on only odd numbers for inc and only even numbers for state,
					//// there is a detectable correlation between streams. I'm not sure how bad it is, but it makes this sometimes
					//// fail Juniper's ICE test.

					//Uint64 oldstate = state;
					//state = state * 0x5851f42d4c957f2dULL + inc;
					//Uint64 word = (oldstate ^ (oldstate >> ((oldstate >> 59u) + 5u))) * 12605985483714917081ULL;
					//return word ^ word >> 43u;

					//// Gets a "suspicious" result at 64TB:
					//// BCFN(2+0,13-0,T)                  R= +12.9  p =  2.0e-6   suspicious

					//Uint64 i = inc ^ rotate64(state, 23);
					//state += 0x5851F42D4C957F2DULL;
					//inc += 0xDA3E39CB94B95BDBULL;
					//Uint64 word = (i ^ i >> (i >> 59u) + 5u) * 12605985483714917081ULL;
					//return (word >> 43u) ^ word;

					//// Passes at least 64TB with no anomalies, but unusually slow...
					//// This takes 3738 seconds to test 1TB, where others took much less than 3000 seconds.

					//Uint64 i = inc ^ rotate64(state, 23);
					//state += 0x5851F42D4C957F2DULL;
					//inc += 0xDA3E39CB94B95BDBULL;
					//i = (i ^ i >> 40) * 12605985483714917081ULL;
					//return (i ^ rotate64(i, 53) ^ rotate64(i, 19));

					//// Passes at least 64TB with no anomalies, but fails Juniper's ICE test badly.

					//Uint64 i = inc ^ state;
					//inc += 0xDA3E39CB94B95BDBULL;
					//state += __lzcnt64(inc);
					//Uint64 r = (i ^ i >> (i >> 59u) + 5u) * 0xF1357AEA2E62A9C5ULL;
					//return (r >> 43u) ^ r;

					//// PactRandom
					//// Passes at least 64TB with no anomalies. It also passes Juniper's ICE test.
					//Uint64 z = (rotate64(state, 13) ^ rotate64(inc, 41) + state) * 0xF1357AEA2E62A9C5ULL;
					//inc += __lzcnt64(state);
					//state += 0xDA3E39CB94B95BDBULL;
					//z = (z ^ rotate64(z, 23) ^ rotate64(z, 47)) * 0xF1357AEA2E62A9C5L;
					//return z ^ z >> 43;

					//i = (i ^ rotate64(i, 34) ^ rotate64(i, 19)) * 12605985483714917081ULL;
					//i = (i ^ i >> 23 ^ i >> 31) * 12605985483714917081ULL;
					//return (i ^ i >> 43);

					// 0xF1357AEA2E62A9C5ULL


//					Uint64 word = (i ^ rotate64(i, 53) ^ rotate64(i, 19)) * 12605985483714917081ULL;
//					Uint64 word = (i ^ i >> 6u ^ i >> 47u) * 12605985483714917081ULL;
//					return word ^ word >> 44u;

					// Modified PCG64, returns 64 bits, 64 bits state.
					// An LCG can be permuted without a variable-distance shift needed at all, and still pass 64TB without anomalies.
					// This one is rather fast, too...
					//Uint64 oldstate = state;
					//state = oldstate * 0x5851f42d4c957f2dULL + inc;
					//Uint64 word = (oldstate ^ rotate64(oldstate, 19) ^ rotate64(oldstate, 41)) * 12605985483714917081ULL;
					//return word ^ word >> 36;

					// Weyl Sequence PCG with two XRXR steps and no variable rotations; 64-bit state, no streams, returns 64 bits.
					// Passes 64TB with no anomalies!
					// Finishes the first TB in less than 45 minutes. Only one multiplication.
					//Uint64 oldstate = state;
					//state = oldstate + 0x5851f42d4c957f2dULL;
					//Uint64 word = (oldstate ^ rotate64(oldstate, 19) ^ rotate64(oldstate, 41)) * 12605985483714917081ULL;
					//return word ^ rotate64(word, 47) ^ rotate64(word, 21);

					// Not quite enough!
					//rng = pcg64(92777c99f3347789), seed = 0x0
					//	length = 4 terabytes(2 ^ 42 bytes), time = 13866 seconds
					//	Test Name                         Raw       Processed     Evaluation
					//	BRank(12) :48K(1)                  R = +6137  p~= 1e-1848    FAIL !!!!!!!!
					//	...and 1051 test result(s) without anomalies
					//Uint64 oldstate = state ^ rotate64(inc, 29);
					//state += 0xD1B54A32D192ED03ULL;
					//inc += 0x8CB92BA72F3D8DD7ULL;
					//Uint64 word = oldstate * 12605985483714917081ULL;
					//return word ^ rotate64(word, 47) ^ rotate64(word, 21);

					//rng = pcg64(c31406d9f3347789), seed = 0x0
					//	length = 2 terabytes(2 ^ 41 bytes), time = 9649 seconds
					//	Test Name                         Raw       Processed     Evaluation
					//	[Low4 / 16]Gap - 16:A                 R = +5.8  p = 3.9e-4   unusual
					//	[Low4 / 16]Gap - 16 : B                 R = +8.8  p = 1.9e-7   very suspicious
					//	...and 1018 test result(s) without anomalies
					//Uint64 oldstate = state ^ rotate64(inc, 37);
					//state += 0xD1B54A32D192ED03ULL;
					//inc += 0x8CB92BA72F3D8DD7ULL;
					//Uint64 word = oldstate * 12605985483714917081ULL;
					//return word ^ rotate64(word, 47) ^ rotate64(word, 13);

					// Passes 64TB!
					// Using a xorshift on state and a rotation on inc lets it pass.
					// 2 to the 64 streams!
					Uint64 oldstate = state ^ state >> 31 ^ rotate64(inc, 37);
					state += 0xD1B54A32D192ED03ULL;
					inc += 0x8CB92BA72F3D8DD7ULL;
					Uint64 word = oldstate * 12605985483714917081ULL;
					return word ^ rotate64(word, 47) ^ rotate64(word, 23);

				}
				std::string pcg64::get_name() const {
					if (inc == 0xda3e39cb94b95bdbULL) return "pcg64";
					std::ostringstream str;
					str << "pcg64(" << std::hex << inc << ")";
					return str.str();
				}
				void pcg64::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
					walker->handle(inc);
				}
				void pcg32_norot::seed(Uint64 s) { state = 0; raw32(); state += s; raw32(); }
				Uint32 pcg32_norot::raw32() {
					Uint64 oldstate = state;
					state = state * 0x5851f42d4c957f2dULL + inc;
					Uint32 xorshifted = Uint32(((oldstate >> 18u) ^ oldstate) >> 27u);
					return xorshifted;
					//Uint32 rot = Uint32(oldstate >> 59u);
					//return (xorshifted >> rot) | (xorshifted << ((~rot + 1u) & 31));
				}
				std::string pcg32_norot::get_name() const {
					if (inc == 0xda3e39cb94b95bdbULL) return "pcg32_norot";
					std::ostringstream str;
					str << "pcg32_norot(" << std::hex << inc << ")";
					return str.str();
				}
				void pcg32_norot::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				void cmrg32of192::seed(Uint64 s) {
					n1m0 = Uint32(s);
					n1m1 = n1m2 = 1;
					n2m0 = Uint32(s >> 32);
					n2m1 = n2m2 = 1;
				}
				Uint32 cmrg32of192::raw32() {
					Uint64 n1 = (Uint64(n1m1) * 1403580 - Uint64(n1m2) * 810728) % ((1ull << 32) - 209);
					Uint64 n2 = (Uint64(n2m0) * 527612 - Uint64(n2m2) * 1370589) % ((1ull << 32) - 22853);
					n1m2 = n1m1;
					n1m1 = n1m0;
					n1m0 = n1;
					n2m2 = n2m1;
					n2m1 = n2m0;
					n2m0 = n2;
					return n1 + n2;
				}
				std::string cmrg32of192::get_name() const {return "cmrg32of192";}
				void cmrg32of192::walk_state(StateWalkingObject *walker) {
					walker->handle(n1m0);
					walker->handle(n1m1);
					walker->handle(n1m2);
					walker->handle(n2m0);
					walker->handle(n2m1);
					walker->handle(n2m2);
				}

				Uint32 xsh_lcg_bad::raw32() {
					Uint64 tmp = x1 ^ (x1 << 11);
					x1 = x2;
					x2 = x3;
					x3 = x0;
					x0 = (x0 >> 19) ^ tmp ^ (tmp >> 8);
					lcg = (lcg * 279470273) % 4294967291;
					return x0 ^ lcg;
				}
				std::string xsh_lcg_bad::get_name() const { return "xsh_lcg_bad"; }
				void xsh_lcg_bad::seed(Uint64 s) {
					x1 = s;
					x0 = x2 = x3 = 0xFFffFFffFFffFFffull;//changed to prevent the bad all-zeroes case
					lcg = 2233445566;
					for (int i = 0; i < 64; i++) raw32();
				}
				void xsh_lcg_bad::walk_state(StateWalkingObject *walker) {
					walker->handle(x0);
					walker->handle(x1);
					walker->handle(x2);
					walker->handle(x3);
					walker->handle(lcg);
				}


				Uint16 mmr16::raw16() {
					Uint16 old = a;
					a = b * 0x69ad;
					b = rotate16(b, 7) ^ c;
					c = rotate16(c, 5) + old;
					return old;
				}
				std::string mmr16::get_name() const { return "mmr16"; }
				void mmr16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}
				Uint32 mmr32::raw32() {
					Uint32 old = a;
					a = b * 0xAC4969AD;
					b = rotate16(b, 13) ^ c;
					c = rotate16(c, 9) + old;
					return old;
				}
				std::string mmr32::get_name() const { return "mmr32"; }
				void mmr32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(c);
				}

				Uint16 garthy16::raw16() {
					if (!counter) scale += 2;
					scale += 2;
					Uint16 temp = value * scale;
					value += ((temp << 7) | (temp >> 9)) ^ counter++;
					return value;
				}
				std::string garthy16::get_name() const { return "garthy16"; }
				void garthy16::walk_state(StateWalkingObject *walker) {
					walker->handle(value); walker->handle(counter); walker->handle(scale);
					scale |= 1;
				}
				Uint32 garthy32::raw32() {
					if (!counter) scale += 2;
					scale += 2;
					Uint32 temp = value * scale;
					value += ((temp << 13) | (temp >> 19)) ^ counter++;
					return value;
				}
				std::string garthy32::get_name() const {return "garthy32";}
				void garthy32::walk_state(StateWalkingObject *walker) {
					walker->handle(value); walker->handle(counter); walker->handle(scale);
					scale |= 1;
				}

				Uint16 binarymult16::raw16() {
					//with Bays-Durham shuffle (size 16) fails @ 32 GB
					Uint16 old = a;
					a = b * (c | 1);
					b = c ^ (old >> 7);
					c ^= old + d++;
					return a;
				}
				std::string binarymult16::get_name() const {return "binarymult16";}
				void binarymult16::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c); walker->handle(d);
				}
				Uint32 binarymult32::raw32() {
					Uint32 old = a;
					a = b * (c | 1);
					b = c ^ (old >> 13);
					c ^= old + d++;
					return a;
				}
				std::string binarymult32::get_name() const {return "binarymult32";}
				void binarymult32::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c); walker->handle(d);
				}

				Uint16 rxmult16::raw16() {
					if (!a) {c++; if (!c) {c=1; d+=2;}}
					a = a * 0x9ad + d;
					b = (((b<<7)|(b>>9)) + a) ^ c;
					Uint16 tmp = b * 5245;
					tmp ^= tmp >> 8;
					return tmp + a;
				}
				std::string rxmult16::get_name() const {return "rxmult16";}
				void rxmult16::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c); walker->handle(d);
					d |= 1;
				}

				Uint64 multish2x64::raw64() {
					Uint64 old = ~a;
					a = (a * 0xa536c4b9) + b;
					b += (old << 21) | (old >> 43);
					return old;
				}
				std::string multish2x64::get_name() const {return "multish2x64";}
				void multish2x64::walk_state(StateWalkingObject *walker ) {
					walker->handle(a); walker->handle(b);
				}
				Uint32 multish3x32::raw32() {
					Uint32 old = a;
					a = (b * 0xa536c4b9) + c++;
					b = ((b << 7) | (b >> 25)) + old;
					return old;
				}
				std::string multish3x32::get_name() const {return "multish3x32";}
				void multish3x32::walk_state(StateWalkingObject *walker ) {
					walker->handle(a); walker->handle(b); walker->handle(c);
				}
				Uint16 multish4x16::raw16() {
					Uint16 old = a;
					if (!c++) d++;
					a = (b^d) * 0x96b9 + c;
					b = ((b << 5) | (b >> 11)) ^ old;
					return old;
				}
				std::string multish4x16::get_name() const {return "multish4x16";}
				void multish4x16::walk_state(StateWalkingObject *walker ) {
					walker->handle(a); walker->handle(b); walker->handle(c); walker->handle(d);
				}

				Uint16 old_mwlac16::raw16() {
					Uint16 oa;
					oa = a;
					a = (b * 0x9785) ^ (a >> 7);
					b = c + (oa >> 2);
					c = d;
					d += ~oa;
					return c;
				}
				std::string old_mwlac16::get_name() const { return "old_mwlac16"; }
				void old_mwlac16::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c); walker->handle(d);
				}
				Uint16 mwlac_varA::raw16() {
					Uint16 oa;
					oa = a * 0x9785;//   1001011110000101
					a = b ^ rotate16(a, 7);
					b += c;
					c = oa;
					return c;
				}
				std::string mwlac_varA::get_name() const { return "mwlac_varA"; }
				void mwlac_varA::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c);
				}
				Uint16 mwlac_varB::raw16() {
					Uint16 oa;
					oa = a * 0x9785;//   1001011110000101
					b = rotate(b, 13);
					a = b ^ rotate16(a, 7);
					b += c;
					c = oa;
					return c;
				}
				std::string mwlac_varB::get_name() const { return "mwlac_varB"; }
				void mwlac_varB::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c);
				}
				Uint16 mwlac_varC::raw16() {
					a *= 0x9785;//   1001011110000101
					b = rotate16(b, 5);
					c = rotate16(c, 13);
					b ^= a;
					a ^= c;
					c += b;
					return b;
				}
				std::string mwlac_varC::get_name() const { return "mwlac_varC"; }
				void mwlac_varC::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c);
				}
				Uint16 mwlac_varD::raw16() {
					a += b;
					b -= c;
					c += a;
					a *= 0x9785;//   1001011110000101
					b = rotate16(b, 7);
					c = rotate16(c, 4);
					return a;
				}
				std::string mwlac_varD::get_name() const { return "mwlac_varD"; }
				void mwlac_varD::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c);
				}
				Uint16 mwlac_varE::raw16() {
					c ^= a;
					a += b;
					b -= c;
					a += c;
					c *= 0x9785;//   1001011110000101
					//shift:	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15
					//						37	39-	41	34-	34	38-	39	38-	35	20
					b = rotate16(b, 6);
					return a;
				}
				std::string mwlac_varE::get_name() const { return "mwlac_varE"; }
				void mwlac_varE::walk_state(StateWalkingObject *walker) {
					walker->handle(a); walker->handle(b); walker->handle(c);
				}

				Uint32 mwc64x::raw32() {
					Uint32 c = state >> 32;
					Uint32 x = Uint32(state);
					state = x * Uint64(4294883355U) + c;
					return x ^ c;
				}
				std::string mwc64x::get_name() const { return "mwc64x"; }
				void mwc64x::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}

				Uint32 xlcg32of64_varqual::raw32() {
					enum {
						X = 0xC74EAD55,//must end in 5 or D
						M = 0x947E3DB3,//must end in 3 or B
					};
					state = (state ^ X) * M;
					return Uint32(state >> outshift);
				}
				std::string xlcg32of64_varqual::get_name() const {
					std::ostringstream str;
					str << "xlcg(32," << (32 + outshift) << ")";
					return str.str();
				}
				void xlcg32of64_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint16 xlcg16of64_varqual::raw16() {
					enum {
						X = 0xC74EAD55,//must end in 5 or D
						M = 0x947E3DB3,//must end in 3 or B
					};
					state = (state ^ X) * M;
					return Uint16(state >> outshift);
				}
				std::string xlcg16of64_varqual::get_name() const {
					std::ostringstream str;
					str << "xlcg(16," << (16 + outshift) << ")";
					return str.str();
				}
				void xlcg16of64_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint8 xlcg8of64_varqual::raw8() {
					enum {
						X = 0xC74EAD55,//must end in 5 or D
						M = 0x947E3DB3,//must end in 3 or B
					};
					state = (state ^ X) * M;
					return Uint8(state >> outshift);
				}
				std::string xlcg8of64_varqual::get_name() const {
					std::ostringstream str;
					str << "xlcg(8," << (8 + outshift) << ")";
					return str.str();
				}
				void xlcg8of64_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint32 xlcg32of128_varqual::raw32() {
					enum {
						X = 0xC74EAD55,//must end in 5 or D
						M = 0x947E3DB3,//must end in 3 or B
					};
					Uint64 a = Uint32(low) * Uint64(M);
					Uint64 b = (low >> 32) * Uint64(M);
					low = a + (b << 32);
					b += a >> 32;
					high = high * M + (b >> 32);
					low ^= X;
					if (outshift >= 64) return Uint16(high >> (outshift - 64));
					if (outshift > 48) return Uint16((low >> outshift) | (high << (64 - outshift)));
					return Uint16(low >> outshift);
				}
				std::string xlcg32of128_varqual::get_name() const {
					std::ostringstream str;
					str << "xlcg(32," << (32 + outshift) << ")";
					return str.str();
				}
				void xlcg32of128_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(low); walker->handle(high);
				}
				Uint16 xlcg16of128_varqual::raw16() {
					enum {
						X = 0xC74EAD55,//must end in 5 or D
						M = 0x947E3DB3,//must end in 3 or B
					};
					Uint64 a = Uint32(low) * Uint64(M);
					Uint64 b = (low >> 32) * Uint64(M);
					low = a + (b << 32);
					b += a >> 32;
					high = high * M + (b >> 32);
					low ^= X;
					if (outshift >= 64) return Uint16(high >> (outshift - 64));
					if (outshift > 48) return Uint16((low >> outshift) | (high << (64 - outshift)));
					return Uint16(low >> outshift);
				}
				std::string xlcg16of128_varqual::get_name() const {
					std::ostringstream str;
					str << "xlcg(16," << (16 + outshift) << ")";
					return str.str();
				}
				void xlcg16of128_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(low); walker->handle(high);
				}
				Uint8 xlcg8of128_varqual::raw8() {
					enum {
						X = 0xC74EAD55,//must end in 5 or D
						M = 0x947E3DB3,//must end in 3 or B
					};
					Uint64 a = Uint32(low) * Uint64(M);
					Uint64 b = (low >> 32) * Uint64(M);
					low = a + (b << 32);
					b += a >> 32;
					high = high * M + (b >> 32);
					low ^= X;
					if (outshift >= 64) return Uint16(high >> (outshift - 64));
					if (outshift > 48) return Uint16((low >> outshift) | (high << (64 - outshift)));
					return Uint16(low >> outshift);
				}
				std::string xlcg8of128_varqual::get_name() const {
					std::ostringstream str;
					str << "xlcg(8," << (8 + outshift) << ")";
					return str.str();
				}
				void xlcg8of128_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(low); walker->handle(high);
				}

				Uint32 cxlcg32of96_varqual::raw32() {
					enum {
						X = 0xC74EAD55,//must end in 5 or D
						M = 0x947E3DB3,//must end in 3 or B
					};
					lcg1 = (lcg1 ^ X) * M;
					Uint64 tmp = lcg2 * Uint64(1579544716);
					lcg2 = Uint32(tmp & 0x7FffFFff) + Uint32(tmp >> 33) + 1;
					return lcg2 + Uint32(lcg1 >> outshift);
				}
				std::string cxlcg32of96_varqual::get_name() const {
					std::ostringstream str;
					str << "cxlcg(32," << (32 + outshift + 32) << ")";
					return str.str();
				}
				void cxlcg32of96_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg1); walker->handle(lcg2);
					if (!lcg2) lcg2 = 1;
				}
				Uint16 cxlcg16of96_varqual::raw16() {
					enum {
						X = 0xC74EAD55,//must end in 5 or D
						M = 0x947E3DB3,//must end in 3 or B
					};
					lcg1 = (lcg1 ^ X) * M;
					Uint64 tmp = lcg2 * Uint64(1579544716);
					lcg2 = Uint32(tmp & 0x7FffFFff) + Uint32(tmp >> 33) + 1;
					return Uint16(lcg2) + Uint16(lcg1 >> outshift);
				}
				std::string cxlcg16of96_varqual::get_name() const {
					std::ostringstream str;
					str << "cxlcg(16," << (16 + outshift + 32) << ")";
					return str.str();
				}
				void cxlcg16of96_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg1); walker->handle(lcg2);
					if (!lcg2) lcg2 = 1;
				}
				Uint8 cxlcg8of96_varqual::raw8() {
					enum {
						X = 0xC74EAD55,//must end in 5 or D
						M = 0x947E3DB3,//must end in 3 or B
					};
					lcg1 = (lcg1 ^ X) * M;
					Uint64 tmp = lcg2 * Uint64(1579544716);
					lcg2 = Uint32(tmp & 0x7FffFFff) + Uint32(tmp >> 33) + 1;
					return Uint8(lcg2) + Uint8(lcg1 >> outshift);
				}
				std::string cxlcg8of96_varqual::get_name() const {
					std::ostringstream str;
					str << "cxlcg(8," << (8 + outshift + 32) << ")";
					return str.str();
				}
				void cxlcg8of96_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(lcg1); walker->handle(lcg2);
					if (!lcg2) lcg2 = 1;
				}

				Uint64 cxm64_varqual::raw64() {
					const Uint64 K = 0x6595a395a1ec531b;
					Uint64 tmp = high >> 32;
					low += K;
					high += K + ((low < K) ? 1 : 0);
					tmp ^= high ^ 0;//(Uint64)this;
					for (int i = 1; i < num_mult; i++) {
						tmp *= K;
						tmp ^= tmp >> 32;
					}
					tmp *= K;
					return tmp + low;
				}
				std::string cxm64_varqual::get_name() const {
					std::ostringstream str;
					str << "cxm" << num_mult << "n64";
					return str.str();
				}
				void cxm64_varqual::walk_state(StateWalkingObject *walker) {
					walker->handle(low);
					walker->handle(high);
					//walker->handle(num_mult);
				}


				Uint32 mo_Cmfr32::raw32() {
					state = ~(2911329625u*state); state = rotate32(state,17);
					return state;
				}
				std::string mo_Cmfr32::get_name() const {return "mo_Cmfr32";}
				void mo_Cmfr32::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint32 mo_Cmr32::raw32() {
					state = 4031235431u * state; state = rotate32(state, 15);
					return state;
				}
				std::string mo_Cmr32::get_name() const { return "mo_Cmr32"; }
				void mo_Cmr32::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint32 mo_Cmr32of64::raw32() {
					state = 38217494031235431ull * state; state = rotate64(state, 37);
					return Uint32(state);
				}
				std::string mo_Cmr32of64::get_name() const { return "mo_Cmr32of64"; }
				void mo_Cmr32of64::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}

				Uint32 murmlac32::raw32() {
					Uint32 tmp = state1;
					for (int i = 0; i < rounds; i++) {
						tmp *= 4031235431u;
						tmp ^= tmp >> 16;
					}
					state1 += state2;
					state2 = tmp;
					return state1;
				}
				std::string murmlac32::get_name() const {
					std::ostringstream str;
					str << "murmlac32(" << rounds << ")";
					return str.str();
				}
				void murmlac32::walk_state(StateWalkingObject *walker) {
					walker->handle(state1); walker->handle(state2);
				}

				Uint64 mulcr64::raw64() {
					Uint64 rv = a * count;
					a = rotate64(a, 24) + b;
					count += 2;
					b = rotate64(b, 37) ^ rv;
					return rv;
				}
				std::string mulcr64::get_name() const { return "mulcr64"; }
				void mulcr64::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(count);
					count |= 1;
				}
				Uint32 mulcr32::raw32() {
					Uint32 rv = a * 2911329625u;
					a = b ^ count++;
					b = rotate32(b, 11) + rv;
					return rv;
				}
				std::string mulcr32::get_name() const { return "mulcr32"; }
				void mulcr32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(count);
				}
				Uint16 mulcr16::raw16() {
					Uint16 rv = a * 2911329625u;
					a = b ^ count++;
					b = rotate16(b, 6) + rv;
					return rv;
				}
				std::string mulcr16::get_name() const { return "mulcr16"; }
				void mulcr16::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					walker->handle(count);
				}






				//0xC5B3CF786C806DF7
				Uint64 mo_Cmres64::raw64() {
					const uint64_t m = state * UINT64_C(0x7CDD0CFB6AD2A499), // 0xC5B3CF786C806DF7 does well up to at least 256GB with rotation by 33, no anomalies
						n = state2 * UINT64_C(0x2C40E78566614E13);
					return (state -= rotate64(m, 35)) ^ (state2 -= rotate64(n, 33));// ^ (weyl += UINT64_C(0x9E3779B97F4A7C15));
				}
				std::string mo_Cmres64::get_name() const { return "mo_Cmres64"; }
				void mo_Cmres64::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
					walker->handle(state2);
					printf("Seed is 0x%016llX, seed2 is 0x%016llX\r\n", state, state2);
				}

				Uint64 tiptoe64::raw64() {
					//WORKS WELL, FAST, NOT EQUIDISTRIBUTED
					//const uint64_t s = state;
					//const uint64_t z = (s ^ s >> 26) * (state += 0x6A5D39EAE116586AULL);
					//return z ^ (z >> 28);

					// passes 2TB no anomalies, does pretty well in gjrand huge, miserably fails gjrand --tera
					//const uint64_t s = (state += 0x6A5D39EAE12657B7ULL); //0x6A5D39EAE12657A9ULL
					//const uint64_t z = (s ^ (s >> 25)) * (s | 1ULL);
					//return z ^ (z >> 19);

					// passes 32TB no failures, passes gjrand huge and tera, but was tested with seed 0 only, which may bias tests. Not equidistributed!
					//const uint64_t s = (state += 0x6C8E9CF570932BD5ULL); //0x6A5D39EAE12657A9ULL
					//const uint64_t z = (s ^ (s >> 25)) * (s | 0xA529ULL);
					//return z ^ (z >> 22);

					//uint64_t z = (state += UINT64_C(0x352E9CF570932BDD)); //UINT64_C(0x6A5D39EAE12657A9)
					//z = ((z >> 25) ^ z) * (z | UINT64_C(1));// (s | 1ULL); //UINT64_C(0x6C8E9CF570932BD5) // UINT64_C(0x2545F4914F6CDD1D)
					//return (z ^ (z >> 22));
					//UINT64_C(0x6C8E9CF570932BD5)

					// left off here, similar to the one that passes 32TB, but not equidistributed
					//state += UINT64_C(0x6C8E9CF570932BD5);
					//uint64_t z = (state ^ state >> 25) * (state | UINT64_C(0xA529));
					//return (z ^ (z >> 22));


					//Uint64 z = (state += 0x9E3779B97F4A7C15);
					//z ^= z >> 27;
					//z *= (Uint64(0x6C8E9CF570932BD5) ^ (z & Uint64(0x8000000000000000)));
					//return z ^ z >> 31;

					// good to 32TB, 2 unusual anomalies
					//Uint64 z = state++;
					//z = ((z ^ rotate64(z, 59) ^ rotate64(z, 23))) * Uint64(0x6C8E9CF570932BD3);
					//z = (z ^ z >> (5 + (z >> 59))) * Uint64(0xAEF17502108EF2D9);
					//return z ^ z >> 25;

					//public static long tiptoe(long state) { return (state = ((state = (state ^ (state << 23 | state >>> 41) ^ (state << 59 | state >>> 5)) * 0x6C8E9CF570932BD5L) ^ state >>> (5 + (state >>> 59))) * 0xAEF17502108EF2D9L) ^ state >>> 25); }
					//Long.rotateLeft(result + 0x6C8E9CF570932BD3L, 17) + 0x9E3779B97F4A7C15
					// seems to work very well; 8TB with one "unusual" anomaly at 1TB
					//Uint64 z = ++state;
					//z = ((z << ((z & 31) + 5u)) ^ z) * Uint64(0x6C8E9CF570932BD3) ^ Uint64(0x9E3779B97F4A7C15);
					//z = (z ^ z >> 24u) * Uint64(0xAEF17502108EF2D9);
					//return (z >> 25u) ^ z;

					//Uint64 z = ++state;
					//z = (z ^ z << 12u) * Uint64(0x6C8E9CF570932BD3);
					//z = (z ^ z >> 25u) * Uint64(0xAEF17502108EF2D9) + Uint64(0x9E3779B97F4A7C15);
					//return z ^ z >> 26u;

					//z ^= z >> 27u;

					// Excellent! Passes 32TB with only 1 anomaly, at 8TB!
					//uint64_t z = state;
					//z = (z ^ rotate64(z, 13) ^ rotate64(z, 31) ^ rotate64(z, 41) ^ rotate64(z, 59)) * UINT64_C(0x6C8E9CF570932BD3);
					//state += UINT64_C(0x9E3779B97F4A7C15);
					//return z ^ z >> 26u;

					// fails at 8TB
					//Uint64 z = (state += Uint64(0x9E3779B97F4A7C15));
					//z = (z ^ rotate64(z, 13) ^ rotate64(z, 41)) * UINT64_C(0x6C8E9CF570932BD3);
					//return z ^ z >> 26u;
					// 0x41C64E6D
					// 0x6F3E235 works pretty well, L'Ecuyer, errata

					//uint64_t z = (state = state * UINT64_C(0x41C64E6D) + UINT64_C(1));//stream);
					//z = (z ^ z >> 32u) * UINT64_C(0xAEF17502108EF2D9);
					//return (z ^ z >> 28u);


					//passes 32TB no anomalies!
					//uint64_t z = state;
					//z = (z ^ z >> 27) * UINT64_C(0xAEF17502108EF2D9);
					//state = state * UINT64_C(1103515245) + UINT64_C(1);
					//return (z ^ z >> 25u);


					//state = (state ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0x6C8E9CF570932BD3);


					//z = (z ^ z << 13u) * Uint64(0x6C8E9CF570932BD3) ^ Uint64(0x6C8E9CF570932BD3) << ((z & 63u) | 8u);

					//public static long tiptoe(long state) { return (state = ((state = (state ^ (state << (state | 8L))) * 0x6C8E9CF570932BD5L) ^ state >>> 24) * 0xAEF17502108EF2D9L) ^ state >>> 25; }

					//uint32_t rot1 = (uint32_t)(state >> 59u);
					//uint32_t high = (uint32_t)(state >> 32u);
					//uint32_t low = (uint32_t)state;
					//uint32_t xored = high ^ low;
					//Uint64 next = rotate32(xored, rot1);
					//return next | (rotate32(high, next & 31u) << 32u);

					// seems to be working well, no anomalies up to 128GB
					//Uint64 z = ++state;
					//z = (z ^ rotate64(z, 58) ^ rotate64(z, 23)) * Uint64(0x6C8E9CF570932BD3);
					//z = (z ^ z >> 26) * Uint64(0xAEF17502108EF2D9);
					//return z ^ z >> 25;


					//return (z ^ z >> (5 + (z >> 59)));

					//return z + (z << (5 + (z & 31)));
					//z += z << (5 + (z & 31));

					//Uint64 z = (state += Uint64(0xBF58476D1CE4E5B9));

					// left off with this at 512GB doing well
					//Uint64 z = ++state ^ Uint64(0x6C8E9CF570932BD5);
					//z *= Uint64(0x94D049BB133111EB);
					//z ^= z >> 22;
					//z = (z ^ z >> 24) * Uint64(0xBF58476D1CE4E5B9);
					//return z ^ z >> 25;

					//passes 32TB but is likely not equidistributed (it could be?)
					//Uint64 z = ++state;
					//z ^= rotate64(z, 58) ^ rotate64(z, 23);
					//Uint64 y = z * Uint64(0x9E3779B97F4A7C15);
					//z = ~(z * Uint64(0x6C8E9CF570932BD6));
					//z *= (y ^= y >> (5 + (y >> 59)));
					//return z ^ z >> 25;

					//Uint64 z = (state += Uint64(0x9E3779B97F4A7C15)), a = z ^ z + 0x94D049BB133111EB, x = 3 * a ^ 2;
					//x *= 2 - a * x;        // 10
					//x *= 2 - a * x;        // 20
					//z = (z ^ z >> 25) * x;
					//x *= 2 - a * x;        // 40
					//x *= 2 - a * x;        // 80
					//z = (z ^ z >> 25) * x;
					//return z ^ z >> 28;

					//Uint64 z = ++state;
					//z = (z ^ rotate64(z, 58) ^ rotate64(z, 23)) * Uint64(0xAEF17502108EF2D9);
					//z = (z ^ z >> 24) * Uint64(0x6C8E9CF570932BD5);
					////z = (z ^ rotate64(z, 42) ^ rotate64(z, 17));
					//return z ^ z >> 25;



					//Uint64 z = (state += Uint64(0x6C8E9CF570932BD5));
					// + Uint64(0x9E3779B97F4A7C15)
					//z = (z ^ z >> 25) * Uint64(0x352E9CF570932BDD);

					//Uint64 z = (state += Uint64(0x6C8E9CF570932BD5));
					////z = (z ^ rotate64(z, 22) ^ rotate64(z, 47)) * Uint64(0x6C8E9CF570932BD5);
					////z *= Uint64(0x6C8E9CF570932BD5);
					//z = (z ^ z >> 26) * Uint64(0x352E9CF570932BDD); //352E9CF5947E3DB3
					//return (z ^ rotate64(z, 22) ^ rotate64(z, 47));

					//public static long tiptoe(long seed) { return ((seed = ((seed *= 0x6C8E9CF570932BD5L) ^ seed >>> 25) * 0x352E9CF570932BDDL) ^ (seed << 22 | seed >>> 42) ^ (seed << 47 | seed >>> 17)); }

					//* Uint64(0x6C8E9CF570932BD5) + Uint64(0x9E3779B97F4A7C15);
					//z = rotate64(z, 37);


					//z -= (rotate64(z, 39) & rotate64(z, 14));
					//uint64_t z = (UINT64_C(0x9E3779B97F4A7C15) ^ (state >> 25)) * (state << 3 ^ UINT64_C(0x6C8E9CF570932BD5));

					//        return (int) ((((state = (x ^ (state += 0x352E9CF570932BDDULL) ^ (state >>> 25)) * (state * x | 1L)) ^ (state >>> 22)) ^
					//((state = (y ^ (state += 0x352E9CF570932BDDL) ^ (state >> > 25)) * (state * y | 1L)) ^ (state >> > 22))) >> > 59);

					//const uint64_t s = (state *= 0x369DEA0F31A53F85ULL);
					//const uint64_t z = (s ^ (s >> 30)) * 0x61C8864680B583EBULL;
					//state += 0x9E3779B97F4A7C15ULL;
					//return z ^ (z >> 26);

					//does well to at least 512 GB, several seeds tried, all good
					//uint64_t z = (state = state * 0x2545F4914F6CDD1DULL + stream),
					//	y = (z ^ z >> 28);
					//return (y ^ (y >> 16)) * 0x5851F42D4C957F2DULL + (0x61C8864680B583EBULL);

					//state = state * 0x369DEA0F31A53F85ULL + stream;
					//return (state ^ state >> 25) + ((state << 19) | (state >> 45)) + ((state << 27) | (state >> 37));

					//state = state * 0x369DEA0F31A53F85ULL + stream;
					//return (rotate64(state, 23) - rotate64(state, 19)) - ((state >> 27) ^ state);

					/*
					uint64_t state = 0ULL; const uint64_t stream = 1ULL; inline uint64_t random() { state = state * 0x369DEA0F31A53F85ULL + stream; return ((state << 23) | (state >> 64 - 23)) - ((state << 19) | (state >> 64 - 19)) - ((state >> 27) ^ state); }
					*/
					//(rotl(state, 27) + rotl(~state, 22)) + ((state >> 26) ^ state);
					//return ((rotate64(state, 37) - rotate64(state, 20))) * 0x2545F4914F6CDD1DULL;
					//return (state ^ state >> 27) + ((state << 27) | (state >> 37)) - ((state << 19) | (state >> 45));
					//return state ^ (rotate64(state, 27) + rotate64(state, 19) + rotate64(-state, 37));
					//return state ^ ((state << 27) | (state >> 37)) + ((state << 19) | (state >> 45)) + ((-state << 37) | (-state >> 27));
					//const uint64_t z = (state ^ state >> 26) + ((state << 27) | (state >> 37));
					//return z ^ z >> 37;

					//uint64_t z = (state = state * 0x369DEA0F31A53F85ULL + 0x1ULL);
					//z = (-z ^ (z >> 25)) * 0x2545F4914F6CDD1DULL;
					//return z ^ z >> 28;
					//0x61C8864680B583EDULL - z

					//return (state ^ (state >> 28)) ^ (((state >> 6) ^ (state >> 33)) * 0x5851F42D4C957F2DL - stream);
					//state = state * 0x2545F4914F6CDD1DULL + stream;
					//const uint64_t z = (state ^ (state >> 28)) * (state | 0x4A529ULL);
					//return z ^ (z >> 28);

					/*
					public ulong state = 1UL; public ulong randomULong() { ulong z = (state *= 0x369DEA0F31A53F85UL); z = (z ^ (z >> 30)) * 0x61C8864680B583EBUL; state += 0x9E3779B97F4A7C15UL;return z ^ (z >> 26); }
					*/
					//z = (j *= 0x369DEA0F31A53F85ULL);
					//z = (z ^ z >> 30) * (z - (j += 0x9E3779B97F4A7C15ULL));
					//z ^= (z >> 26);



					//Uint64 z = state;
					//z ^= (((state += 0xA99635D5B8597AE5ULL) ^ z >> 23) * 0xAD5DE9A61A9C3D95ULL);
					//return z ^ (z >> 29);

					//z ^= (((state += 0x9E3779B97F4A7BB5ULL) ^ z >> 23) * 0x6A5D39EAE116586DULL);
					//z ^= ((state += 0x9E3779B97F4A7C15ULL) ^ (z >> 25)) * 0x6A5D39EAE116586DULL;
					//Uint64 z = (state += 0x9E3779B97F4A7C15ULL);
					//z ^= (z ^ (z >> 25)) * 0x6A5D39EAE116586DULL;
					//return z ^ (z >> 24);

					//z = (z ^ z >> 26) * 0x5851F42D4C957F2DULL;
					//z *= (z << 1) ^ 0x5851F42D4C957F2D; // ^ (z << 11) //0x2545F4914F6CDD1D

					//const uint64_t s = state;
					//const uint64_t z = (s >> 29) * (s | (state += 0x9E3779B97F4A7A55ULL));
					//return z ^ (z >> 29);

					//LinnormRNG normal
					//uint64_t z = (state = state * UINT64_C(0x41C64E6D) + UINT64_C(1));
					//z = (z ^ z >> 27u) * UINT64_C(0xAEF17502108EF2D9);
					//return z ^ z >> 25u;

					//uint64_t z = (state = (state ^ UINT64_C(0x6C8E9CF570932BD5)) * UINT64_C(0x41C64E6B));
					//uint64_t z = (state += UINT64_C(0x9E3779B97F4A7C15));
					//return (stream = (stream >> 1 ^ (-(stream & 1) & UINT64_C(0xD800000000000000))));
					//(stream = (stream >> 1 ^ (-(stream & 8) & UINT64_C(0x1C80000000000000))));
					//((stream = (stream >> 1 ^ (-(stream & 8) & UINT64_C(0x1C80000000000000)))) << 3 ^ UINT64_C(5))

					//QuixoticRNG attempt for speed maybe
					//uint64_t z = (state = (state ^ UINT64_C(0x6C8E9CF570932BD5)) * UINT64_C(0x42042E4DD58B));
					//uint64_t z = (state = (state ^ UINT64_C(0x6C8E9CF570932BD5)) * UINT64_C(0x5DA942042E4DD58B));// 0x5888120408461063 // 0x4120508604029043 //0x482050128A101603
					
					//QuixoticRNG normal
					//uint64_t z = (state = (state ^ UINT64_C(0x6C8E9CF570932BD5)) * UINT64_C(0x41C64E6B));
					//z = (z ^ z >> 27u) + UINT64_C(0xAEF17502108EF2D9);
					//return z ^ z >> 25u;
					
					//uint64_t z = (state = (state ^ UINT64_C(0x6C8E9CF570932BD5)) * UINT64_C(0xC6BC279692B5CC83));
					//z = rotate64(z, 27);
					//z = (z << 26) - z;
					//return z ^ z >> 25u;

					// fails at 512GB
					//					uint64_t z = (state = state * UINT64_C(0x41C64E6B) + UINT64_C(1));
					//					z = (z ^ z >> 27u) + UINT64_C(0xAEF17502108EF2D9);
					//					return z ^ z >> 25u;

					// works really well for some reason???
					//uint64_t z = (state += UINT64_C(0x9E3779B97F4A7C15)) + (stream = (stream >> 1 ^ (-(stream & 1) & UINT64_C(0xD800000000000000))));
					//z = (z ^ z >> 26u) * UINT64_C(0xAEF17502108EF2D9);
					//return (z ^ z >> 29u);

					//// LinnormRNG.determine()
					//uint64_t z = ((state += 0x632BE59BD9B4E019UL) ^ 0x9E3779B97F4A7C15UL) * 0xC6BC279692B5CC83UL;
					//z = (z ^ z >> 27u) * 0xAEF17502108EF2D9UL;
					//return (z ^ z >> 25u);

					//^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xE19B01AA9D42C633);
					//z = (z ^ z >> 28 ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xE19B01AA9D42C633);

					//uint64_t z = ((++state) ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0xD1B54A32D192ED03);
					//z = (z ^ rotate64(z, 23) ^ rotate64(z, 42)) * UINT64_C(0xDB4F0B9175AE2165);
					//return (z ^ z >> 30);

					//uint64_t z = ++state;
					//z = (z ^ rotate64(z, 21) ^ rotate64(z, 45) ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xD1B54A32D192ED03);
					//z = (z ^ rotate64(z, 21) ^ rotate64(z, 45) ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xD1B54A32D192ED03);
					//return (z ^ z >> 27);

					//uint64_t z = ((++state) ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0xD1B54A32D192ED03);
					////z = (z ^ rotate64(z, 21) ^ rotate64(z, 11)) * UINT64_C(0xDB4F0B9175AE2165);
					//z = (z ^ rotate64(z, 38) ^ rotate64(z, 59) ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xC6BC279692B5CC83);
					//return z ^ z >> 27;


					// Quixotic determine attempt?
					// where we left off
					//uint64_t z = ((state += UINT64_C(0x9E3779B97F4A7C15)) ^ UINT64_C(0xAC51F42D4C957F2D)) + UINT64_C(0xC6BC279692B5CC85);// ^ UINT64_C(0x5851F42D4C957F2D)) * UINT64_C(0x41C64E6B); //0xC6BC279692B5CC83
					//z = (z ^ z >> 28) * UINT64_C(0xAEF17502108EF2DB);// *UINT64_C(0xAEF17502108EF2D9); // ^ UINT64_C(0x6C8E9CF570932BD5)
					//z = (z ^ z >> 17) + UINT64_C(0x6C8E9CF570932BD5) ^ UINT64_C(0x9E3779B97F4A7C15);
					//return (z ^ z >> 31u);

					// PaperweightRNG; rather slow
					//uint64_t z = (state -= UINT64_C(0x9E3779B97F4A7C15));
					//z ^= z >> 13;
					//z = (z << 19) - z;
					//z ^= z >> 12;
					//z = (z << 17) - z;
					//z ^= z >> 14;
					//z = (z << 13) - z;
					//return z ^ z >> 24;


					//uint64_t s = (state -= UINT64_C(0x9E3779B97F4A7C15));// ^ UINT64_C(0x6A5D39EAE127586A);
					//s = ((s << 10) - rotate32(s, 6));
					//s -= s << 7;
					//s ^= s >> 15;// ^ UINT64_C(0x6C8E9CF570932BD5);
					//s = (s << 15) - rotate64(s, 11);
					//s -= s << 5;
					//s ^= s >> 11;
					//s = (s << 9) - rotate64(s, 5);
					//s -= s << 13;
					//return s ^ s >> 27;

					// 0xFF51AFD7ED558CCD 0xC4CEB9FE1A85EC53

					//uint64_t x = (state += UINT64_C(0x9E3779B97F4A7C15));
					//x = (x ^ x >> 16) * UINT64_C(0x85EBCA6B);
					//x = (x ^ x >> 13) * UINT64_C(0xC2B2AE35);
					//return x ^ x >> 16;

					//  ^ UINT64_C(0x2545F4914F6CDD1D)) * UINT64_C(0x41C64E6B)
					//const uint64_t z = (state += UINT64_C(0x9E3779B97F4A7C15)), x = (z ^ rotate64(z, 19) ^ rotate64(z, 43)) * UINT64_C(0x9E3779BD);
					//return (x ^ rotate64(x, 29) ^ rotate64(x, 37));

					//const uint64_t z = (state += UINT64_C(0x9E3779B97F4A7C15)), x = (z ^ rotate64(z, 18) ^ rotate64(z, 47)) * UINT64_C(0x41C64E6D);
					//return (x ^ rotate64(x, 25) ^ rotate64(x, 38));

					////WHERE WE LEFT OFF
					//uint64_t y = state, z = y ^ (state += UINT64_C(0xC6BC279692B5CC85));
					//y -= rotate64(z, 29);
					//y *= z;
					//return y ^ y >> 28;

					//x = (x ^ x >> 31) * UINT64_C(0x369DEA0F31A53F85) + UINT64_C(0xAEF17502108EF2D9);
					//x = ((x ^ x >> 34) * UINT64_C(0x41C64E6D) ^ UINT64_C(0x6C8E9CF570932BD5)) * UINT64_C(0x41C64E6B);
					// (z ^ rotate64(z, 42) ^ rotate64(z, 23))
					//z -= (z << 19);// -UINT64_C(0x6C8E9CF570932BD5);
					//z ^= rotate64(z, 39) ^ rotate64(z, 11);// ^ UINT64_C(0xC6BC279692B5CC85);
					//z ^= z >> 11;
					//z -= (z << 13); - UINT64_C(0x6C8E9CF570932BD5);
					//return x ^ x >> 31;

					// seems pretty good, not too fast though
					//uint64_t z = (state -= UINT64_C(0x9E3779B97F4A7C15));
					//z ^= rotate64(z, 42) ^ rotate64(z, 23);// ^ UINT64_C(0x6C8E9CF570932BD5);
					//z -= (z << 19);
					//z ^= rotate64(z, 39) ^ rotate64(z, 11);// ^ UINT64_C(0xC6BC279692B5CC85);
					//z -= (z << 13);
					//return z ^ z >> 31;
					//

					//s = (s ^ s >> 31 ^ UINT64_C(0xC6BC279692B5CC85)) * UINT64_C(0xAEF17502108EF2DB);
					//s = (s ^ s >> 30 ^ UINT64_C(0xAC51F42D4C957F2D)) * UINT64_C(0x41C64E6B);

					//const uint64_t t = ((s << 11) - rotate64(s, 9)) * 127;// UINT64_C(0x6C8E9CF570932BD5);
					//return ((s << 5) - rotate64(s, 3)) ^ (t ^ t >> 30);// ((t << 11) - rotate64(t, 9));


					//s = (s ^ s >> 33 ^ UINT64_C(0xC6BC279692B5CC85)) * UINT64_C(0xAEF17502108EF2DB);
					//return s ^ rotate64(s, 17) ^ rotate64(s, 42);
					//return (s << 7) - rotate64(s, 3);
					//s ^= s >> 30;


					//s = ((s << 14) - rotate64(s, 12)) * UINT64_C(0xC6BC279692B5CC85);
					//s = ((s >> 31) ^ s ^ UINT64_C(0xC6BC279692B5CC85));
					//uint64_t t = s * UINT64_C(3188803096312630803);
					//t = rotate64(t, 33) - s;
					//const uint64_t t = state * 3188803096312630803L;
					//state -= rotate64(t, 33);
					//return (state << 23) - rotate64(state, 19);
					//uint64_t z = ((state += UINT64_C(0x9E3779B97F4A7C15)) ^ UINT64_C(0x5851F42D4C957F2D)) + UINT64_C(0xAEF17502108EF2DB);// ^ UINT64_C(0x5851F42D4C957F2D)) * UINT64_C(0x41C64E6B); //0xC6BC279692B5CC83
					//z = (z ^ z >> 28) * UINT64_C(0xC6BC279692B5CC83);// *UINT64_C(0xAEF17502108EF2D9); // ^ UINT64_C(0x6C8E9CF570932BD5)
					//return (z ^ z >> 27u);


					//uint64_t z = (state += UINT64_C(0x9E3779B97F4A7C15));
					//uint64_t z = ((state += UINT64_C(0x632BE59BD9B4E019)) ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0xC6BC279692B5CC83); 
					//z = (z ^ z >> 26u) * UINT64_C(0xAEF17502108EF2D9);
					//return (z ^ z >> ((z >> 59u)+5u));

					//uint64_t z = (++state ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0xC6BC279692B5CC83);
					//uint64_t z = (state += UINT64_C(0xC6BC279692B5CC83));
					//z ^= -(z >> 63) & UINT64_C(0x632BE59BD9B4E019);
					//z = (z ^ z >> 27u) * UINT64_C(0xAEF17502108EF2D9);
					//z ^= z >> 26u;
					//z = (z ^ UINT64_C(0x6A5D39EAE116586D)) * UINT64_C(0x94D049BB133111EB);// +(-(z >> 63) & UINT64_C(0x9E3779B97F4A7C15));// UINT64_C(0xBC6EF372FE94F82A));
					//z *= (-(z >> 63u) | UINT64_C(0x00FF00FF00FF00FF)) & UINT64_C(0x369DEA0F31A53F85);
					//return (z ^ z >> 30u);// 26u);


					//uint64_t s = (stream = (stream >> 1 ^ (-(stream & 1) & UINT64_C(0xD800000000000000)))),
					//uint64_t z = (state = (state ^ 0xC74EAD55u) * 0x947E3DB3u);
					//^ ~((stream = (stream >> 1 ^ (-(stream & 1) & UINT64_C(0xD800000000000000)))) & UINT16_C(2)); // & UINT64_C(0x6A5D39EAE116586D)
					//z = (z ^ z >> 26u) * UINT64_C(0xAEF17502108EF2D9);
					//return (z ^ z >> 29u);

					//uint64_t z = (state = state * UINT64_C(0x41C64E6D) + UINT64_C(1));
					//z = ((stream = (stream >> 1 ^ (-(stream & 1) & UINT64_C(0xD800000000000000)))) + (z ^ z >> 32u)) * UINT64_C(0xAEF17502108EF2D9);
					//return (z ^ z >> 30u);

					//philox, roughly
					//uint64_t l = (state += UINT64_C(0x9E3779B97F4A7C15)), r;
					//l = _mulx_u64(l, UINT64_C(0xD2B74407B1CE6E93), &r) + stream;
					//r = _mulx_u64(r ^ l, UINT64_C(0xD2B74407B1CE6E93), &l);
					//return l ^ r;

					//Passes 32TB with only one unusual anomaly, around 4GB
					//const uint64_t l = (state += UINT64_C(0x9E3779B97F4A7C15));
					//uint64_t x = (l ^ l >> 28) * UINT64_C(0xCA5A826395121157), y = (x >> 32 ^ l) * UINT64_C(0xD2B74407B1CE6E93);
					//x += rotate64(y, 19);// (y ^ y >> 32);// +UINT64_C(0xD2B74407B1CE6E93);
					//y += rotate64(x, 13);// (x ^ x >> 32);// +UINT64_C(0xCA5A826395121157);
					//return x ^ y;
					//const uint64_t l = (state += UINT64_C(0x9E3779B97F4A7C15));
					//uint64_t x = (l ^ l >> 28), y = (x ^ l >> 17) * UINT64_C(0xD2B74407B1CE6E93);
					//y += rotate64(x, 19);// (y ^ y >> 32);// +UINT64_C(0xD2B74407B1CE6E93);
					//x += rotate64(y, 13);// (x ^ x >> 32);// +UINT64_C(0xCA5A826395121157);
					//return x ^ y;

					//// needs work
					//uint64_t z = ((++state) ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0xD1B54A32D192ED03);
					//z = (z ^ z >> 30u ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xE60E2B722B53AEEB);
					//return (z ^ z >> 28u);

					//uint64_t z = (state = (state ^ UINT64_C(0x6C8E9CF570932BD5)) * UINT64_C(0xD1B54A32D192ED03));
					//z -= z >> 28;
					//return z ^ z >> 26;
                    //uint64_t z = (state += UINT64_C(0x6C8E9CF570932BD5));
					//z ^= z >> 25;
					//z = rotate64(z, 21) * UINT64_C(0xDB4F0B9175AE2165);
					//return rotate64(z, 21) * UINT64_C(0x9E3779B97F4A7C15);
					
					//uint64_t z = (state += UINT64_C(0x9E3779B97F4A7C15));
					//z = (UINT64_C(0xDB4F0B9175AE2165) ^ (z ^ z >> 28)) * UINT64_C(0xD1B54A32D192ED03);
					//z ^= z << 15;
					//return z ^ z >> 21;

					//return s ^ rotate64(s, 11) ^ rotate64(s, 21);
					//s = (state = (z ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0x41C64E6B)) - (z >> 28);// 0xD1B54A32D192ED03 0x41C64E6B   0x6C8E9CF570932BD5
					//               s = (state = (z * UINT64_C(0x369DEA0F31A53F85)) + UINT64_C(0x6C8E9CF570932BD5)) - (z >> 27);
					//const uint64_t h = rotate64(l, 29);// +UINT64_C(0xD2B74407B1CE6E93);
					// philox 0xCA5A826395121157 0xD2B74407B1CE6E93
					// neely  0xC6BC279692B5CC83
					// golden 0xD1B54A32D192ED03 0xDB4F0B9175AE2165
                    // even   0xC6BC2796 0xD1B54A32
					// odd    0xDB4F0B91 0xCA5A8263 0xD2B74407
					//const uint64_t z = (l ^ l >> 25 ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xC6BC279692B5CC83) + UINT64_C(0x632BE59BD9B4E019);// *UINT64_C(0xDB4F0B9175AE2165);// +UINT64_C(0xCA5A826395121157);
					//const uint64_t z = (l ^ l >> 27 ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xC6BC279692B5CC83);// +UINT64_C(0xCA5A826395121157);
					//return (z ^ z >> 26);
					//return (z ^ rotate64(z, 46) ^ rotate64(z, 19));
					// return (state = ((state *= 0xD2B74407B1CE6E93L) ^ (state << 19 | state >>> 45) ^ (state << 37 | state >>> 27)) * 0xCA5A826395121157L) ^ state >>> 28;
					// return ((s = ((s *= 0x9E3779B97F4A7C15L) ^ s >>> 27 ^ 0xDB4F0B9175AE2165L) * 0xC6BC279692B5CC83L) ^ (s << 11 | s >>> 53) ^ (s << 23 | s >>> 41));
				    //uint64_t r = (state += 0x6C8E9CF570932BD5ULL);
					//r = (r | 0xA529ULL) * (r ^ (r >> 25));//(r ^ 0xe7037ed1a0b428dbull); // (r | 0xA529ull)
					//return r ^ (r >> 23);
					
					//uint64_t z = ((state += 0x632BE59BD9B4E019UL) ^ 0x9E3779B97F4A7C15UL) * 0xC6BC279692B5CC83UL;

					//// MizuchiRNG retry, passes 32TB (when PractRand seed is 0xafec9d14, 3 "unusual" anomalies)
					//// (when PractRand seed is 0x339d439a and stream is fixed at 1, no anomalies)
					//uint64_t z = (state = state * 0x369DEA0F31A53F85ULL + stream);
					//z = (z ^ z >> 23u ^ z >> 47u) * 0xAEF17502108EF2D9UL;
					//return (z ^ z >> 25u);

					//uint64_t z = (state += 0xEB44ACCAB455D165ULL);
					//z = (z ^ z >> 23 ^ z >> 47) * (z ^ z * 0x369DEA0F31A53F85ULL + stream);
					//return (z ^ z >> 25);
					
					
					//uint64_t z = (state = state * UINT64_C(0x41C64E6D) + UINT64_C(1));
					//z = (z ^ z >> 27u) * UINT64_C(0xAEF17502108EF2D9);
					//return z ^ z >> 25u;

					//// passes 32TB with no anomalies; slower than SplitMix64 on recent Intel
					//uint64_t z = (state += 0xEB44ACCAB455D165ULL);
					//z = (z ^ z >> 29 ^ z >> 43 ^ z << 7 ^ z << 53) * 0xDB4F0B9175AE2165ULL;
					//return z ^ z >> 26;

					//// passes 32TB, one unusual anomaly at 16TB
					//uint64_t z = (state += 0xEB44ACCAB455D165ULL);
					//z = (z ^ z >> 29 ^ z >> 43 ^ z << 7 ^ z << 53) * 0xDB4F0B9175AE2165ULL;
					//return z ^ z >> 23 ^ z >> 41;

					//uint64_t z = state++ ^ 0xEB44ACCAB455D165ULL;
					//z = (z ^ z >> 29 ^ z >> 43 ^ z << 7 ^ z << 53) * 0xDB4F0B9175AE2165ULL;
					//return z ^ z >> 26;


					//return z ^ z >> 19 ^ z >> 53 ^ z << 24;
					//z = (z ^ z >> 23u ^ z >> 47u) * 0xDB4F0B9175AE2165ULL;//0xD1B54A32D192ED03ULL;
					//return (z ^ z >> 26);

//					s = (s ^ rotate64(s, 41) ^ rotate64(s, 17) ^ 0xD1B54A32D192ED03UL) * 0xAEF17502108EF2D9UL;
//					z = (z ^ z >> 43 ^ z >> 31 ^ z >> 23) * 0xDB4F0B9175AE2165UL;
//					return s ^ s >> 28;

//                    //// VibrantRNG, passes 64TB with no anomalies.
//					//// This is the same as MizuchiRNG with different multipliers.
//					uint64_t z = (state = state * 0xD1342543DE82EF95UL + stream);
//					z = (z ^ z >> 23 ^ z >> 47) * UINT64_C(0xDB4F0B9175AE2165);
//					return z ^ z >> 25;

//					// DiverRNG, verbatim
//					uint64_t z = (state = (state ^ UINT64_C(0x6C8E9CF570932BD5)) * UINT64_C(0xC6BC279692B5CC83));
//					z = rotate64(z, 27) * UINT64_C(0xDB4F0B9175AE2165);
//					return z ^ (z >> 25u);

//					// One anomaly at 64TB,
//					//   [Low4/32]Gap-16:B                 R=  -5.0  p =1-3.9e-4   unusual
//					uint64_t z = (state = -(state ^ (state * state | 5U)));
//					z ^= z >> 31;
//					z *= 0xDB4F0B9175AE2165;
//					z ^= z >> 30;
//					return z;

//					uint64_t z = (state = -(state ^ (state * state | 5U)));
//					uint64_t z = (state = state * 0xD1342543DE82EF95 + 0xDB4F0B9175AE2165);
//					z ^= z >> 30;
//					z *= 0xF1357AEA2E62A9C5;
//					z ^= z >> 31;
//					return z;

//					// Some kind of XQO-based random?
//					// Passes 64TB with no anomalies.
//					uint64_t z = (state += 0x9E3779B97F4A7C15UL);
//					z ^= z * z | 1UL;
//					z ^= z >> 21 ^ z >> 3;
//					z ^= z * z | 1UL;
//					return z ^ z >> 31 ^ z >> 5;

					// uint64_t z = (++state);
					//& 0xFFFFFFFFFFFFFFFEUL;
					// z *= 0xF1357AEA2E62A9C5UL;
					// return z ^ z >> 26;

					// uint64_t z = (state += 0x9E3779B97F4A7C15UL);
					//					z = 0xDE916ABCC965815BUL * (z ^ (z * z | 5UL));
					// z = -(z ^ (z * z | 5UL));
					// z = 0xC6BC279692B5CC83UL * (z ^ (z * z | 5UL));

					// uint64_t z = (state = -(state ^ (state * state | 5U))) + (stream += 0x9E3779B97F4A7C16UL);

					//// Prototype of QuestRandom
					//// The simplest 64-bit Xorshift for state, Weyl sequence for stream.
					//// Passes 64TB with one anomaly at 2TB:
					////   [Low1/8]FPF-14+6/16:cross         R=  +4.9  p =  1.8e-4   unusual
//					state ^= state << 7;
//					uint64_t z = state ^= state >> 9;
//					z *= (stream += 0x9E3779B97F4A7C16UL);
//					z ^= z >> 31;
//					z *= stream;
//					return z ^ z >> 29;

					//// QuestRandom
					//// LFSR for state, odd Weyl sequence for stream.
					//// Passes 64TB with no anomalies. Long enough period?
//					uint64_t z = (state = state >> 1 ^ (-(state & 1UL) & 0x9C82435DE0D49E35UL));
//					// z ^= z >> 31;
//					z *= (stream += 0x9E3779B97F4A7C16UL);
//					z ^= z >> 31;
//					z *= stream;
//					return z ^ z >> 29;
					// z ^= z >> 29;
					// z *= 0xF1357AEA2E62A9C5UL;
					// return z ^ z >> 31;

					//0xD1342543DE82EF95UL is an LCG constant, 0xF1357AEA2E62A9C5UL is an MCG constant.
					//0x2ECBDABC217D106BUL is - 0xD1342543DE82EF95UL

					//uint64_t z = (state = -(state ^ (state * state | 5U)));
					// uint64_t z = (state = state * 0xD1342543DE82EF95 + 0xDB4F0B9175AE2165);
					// z ^= rotate64(z, 25) ^ rotate64(z, 40);
					// return z ^ rotate64(z, 13) ^ rotate64(z, 53);

					//// NarwhalRandom
					//// Passes 64TB with no anomalies.
					//// State has two parts, with state (a counter with a variable step) dependent on stream (an LFSR with many taps).
					//// Period is (2 to the 128) - (2 to the 64). While state can be anything, stream can't be 0.
					// uint64_t z = (state += 0x9E3779B97F4A7C15UL + (stream = stream >> 1 ^ (-(stream & 1UL) & 0x9C82435DE0D49E35UL)));
					// z ^= z >> 30;
					// z *= 0xF1357AEA2E62A9C5UL;
					// return z ^ z >> 29;

					//// Passes 64TB with no anomalies, but near-instantly fails Remortality.
					// uint64_t z = (state += 0xDB4F0B9175AE2165UL);
					// z *= (stream += 0x9E3779B97F4A7C16UL);
					// z ^= z >> 33;
					// z *= stream;
					// return z ^ z >> 30;

					// IceRandom
					// Passes 64TB with no anomalies.
					// Two states, updated with the same period (2 to the 64), which enables constant-time skipping/jumping.
					// ARX; mostly an effort to see if this type of generator could work with so few operation types.

					// Passes over 179PB of ReMort testing:
// 6	                   0  - 2.29283e-01	                   0  - 2.29283e-01	                   0  - 2.29283e-01	                   1  + 2.59071e+00	                   0  - 2.29283e-01	
// 7	                4394  - 1.01122e-01	                4356  - 7.91898e-01	                4492  + 1.33836e+00	                4420  + 5.37216e-03	                4440  + 1.40092e-01	
// 8	              361337  + 6.48559e-01	              359768  - 3.26371e+00	              361689  + 1.93573e+00	              361906  + 3.07141e+00	              360498  - 3.49692e-01	
// 9	            22348946  - 2.86186e+00	            22350553  - 1.82746e+00	            22361959  + 1.12454e+00	            22361813  + 1.06000e+00	            22349745  - 2.31868e+00	
//11	                4544  + 3.76150e+00	                4501  + 1.67010e+00	                4456  + 3.78329e-01	                4437  + 1.08333e-01	                4472  + 7.32531e-01	
//12	            85018021  - 9.00074e-03	            85026661  + 7.09239e-01	            85000873  - 3.82057e+00	            85011585  - 6.28654e-01	            85005707  - 2.04594e+00	
//13	          6948569009  - 1.91190e+00	          6948557891  - 2.29853e+00	          6948667171  - 4.20789e-02	          6948520692  - 3.85079e+00	          6948720506  + 1.88958e-01	
//14	        430511514277  + 2.50723e-01	        430511516798  + 2.54586e-01	        430511162929  - 1.20825e-03	        430511298715  + 2.96490e-02	        430511375185  + 8.33680e-02	
//16	              361058  + 1.16200e-01	              360196  - 1.19702e+00	              360652  - 1.12215e-01	              360816  - 3.84082e-03	              360546  - 2.61573e-01	
//17	          6948698221  + 2.80076e-02	          6948723979  + 2.26916e-01	          6948729194  + 2.90432e-01	          6948765715  + 9.54599e-01	          6948772005  + 1.10774e+00	
//18	        567923733319  + 3.08741e-01	        567922020292  - 2.94967e+00	        567923029789  - 1.42813e-01	        567922574776  - 9.63708e-01	        567922977780  - 1.99737e-01	
//19	      35186130431231  - 4.78100e-01	      35186132119362  - 1.65533e-01	      35186122108108  - 4.38729e+00	      35186122526436  - 4.09683e+00	      35186131118667  - 3.31266e-01	
//21	            22362538  + 1.39924e+00	            22349986  - 2.16605e+00	            22349575  - 2.42946e+00	            22362886  + 1.57878e+00	            22363122  + 1.70670e+00	
//22	        430511112748  - 1.23743e-02	        430511350874  + 6.33445e-02	        430511371311  + 7.99933e-02	        430511051664  - 4.17535e-02	        430511051232  - 4.20230e-02	
//23	      35186121582972  - 4.76599e+00	      35186132291047  - 1.42819e-01	      35186131170349  - 3.21313e-01	      35186122789263  - 3.91943e+00	      35186122187853  - 4.33115e+00	
//24	    2179984587582633  + 1.24598e-01	    2179984576648984  + 1.41158e-02	    2179984587002701  + 1.15983e-01	    2179984595690123  + 2.77337e-01	    2179984587033490  + 1.16433e-01	
//
//              16.778 15 =>     7.319156e-01           17.741 15 =>     7.814574e-01           16.520 15 =>     7.168986e-01           20.590 15 =>     8.874546e-01           13.956 15 =>     5.468045e-01
//
//--- Finished -- ReMort trials: 2251799813685248	2^51.0 (out of 2^51) trials -- Ice -- 64 bits: 	ps:   7.319e-01    7.815e-01    7.169e-01    8.875e-01*   5.468e-01   => p <   1.98124e-01   2^54.32 calls, 2^57.32 bytes	2^37.23 bytes/second	used:  12::21:04:21.12

//					uint64_t z = (state += 0xDB4F0B9175AE2165UL);
//					uint64_t y = (stream += 0x9E3779B97F4A7C15UL);
//					y += z ^ rotate64(z, 11) ^ rotate64(z, 50);
//					z += rotate64(y, 17);
//					y += z ^ rotate64(z, 46) ^ rotate64(z, 21);
//					z += rotate64(y, 47);
//					y += z ^ rotate64(z, 41) ^ rotate64(z, 25);
//					return y;


//// Gets very suspicious at 8TB, but ReMort doesn't detect anything in the same timeframe (12 hours).
//					uint64_t z = (state += 0xDB4F0B9175AE2165UL);
//					uint64_t y = (stream += 0x9E3779B97F4A7C15UL);
//					y += z ^ rotate64(z, 11) ^ rotate64(z, 50);
//					z += y ^ rotate64(y, 46) ^ rotate64(y, 21);
//					y += z ^  rotate64(z, 41) ^ rotate64(z, 25);
//					return y;

//					uint64_t z = (state += 0xDB4F0B9175AE2165UL);
//					uint64_t y = (stream += 0x9E3779B97F4A7C15UL);
//					uint64_t a = y + (z ^ rotate64(z, 11) ^ rotate64(z, 50));
//					uint64_t b = z + (y ^ rotate64(y, 37) ^ rotate64(y, 21));
//					return (a ^ rotate64(a, 44) ^ rotate64(a, 26)) + (b ^ rotate64(b, 53) ^ rotate64(b, 17));

					// uint64_t y = (stream ^ stream >> 7);
					// stream = y ^ y << 9;
					// uint64_t z = (state += 0x9E3779B97F4A7C15UL);
//					uint64_t y = (state += 0xDB4F0B9175AE2165UL);
//					uint64_t z = (stream = stream >> 1 ^ (-(stream & 1UL) & 0x9C82435DE0D49E35UL));
//					uint64_t a = (y + (z ^ rotate64(z, 10) ^ rotate64(z, 50)));
//					uint64_t b = (z + (y ^ rotate64(y, 46) ^ rotate64(y, 21)));
//					return (a ^ rotate64(a, 17) ^ rotate64(a, 41)) + (b ^ rotate64(b, 23) ^ rotate64(b, 47));

					// uint64_t z = (state = state >> 1 ^ (-(state & 1UL) & 0x9C82435DE0D49E35UL));
					// uint64_t y = (stream += 0x9E3779B97F4A7C15UL);
					// uint64_t y = (stream ^ stream >> 7);
					// stream = y ^ y << 9;
					// uint64_t z = (state += 0x9E3779B97F4A7C15UL);

					//uint64_t y = (state += 0xDB4F0B9175AE2165UL);
					//uint64_t z = (stream += 0x9E3779B97F4A7C15UL);
					// z += y ^ rotate64(y, 6)  ^ rotate64(y, 25);
					// y += z ^ rotate64(z, 41) ^ rotate64(z, 59);
					// z += rotate64(y, 47) ^ rotate64(y, 11) ^ rotate64(y, 25);
					// y += rotate64(z, 21) ^ rotate64(z, 50) ^ rotate64(z, 59);
					// z += z ^ rotate64(z, 11) ^ rotate64(z, 50);
					// // y += rotate64(y, 17);
					// z += z ^ rotate64(z, 46) ^ rotate64(z, 21);
					// // y += rotate64(y, 47);
					// z += z ^ rotate64(z, 41) ^ rotate64(z, 25);
					// return y ^ z;

					// WrangleRNG
					// Passes 64TB with no anomalies.
					// Period is exactly 2 to the 64, and there are 2 to the 63 streams of unknown correlation.
					// The increments this uses are meant to be a serious stress-test for the generator; 1 and 3.
					// Even with that pair, this seems strong. It uses only add, rotate, and XOR operations.
					// Instruction pipelining has a lot of possibility here; half of the rotations can be pipelined
					// to compute at the same time as another rotation (I think).
					// Yes, really strong, even with 179PB of Remortality:
// 6	                   0  - 2.29283e-01	                   0  - 2.29283e-01	                   1  + 2.59071e+00	                   0  - 2.29283e-01	                   1  + 2.59071e+00	
// 7	                4468  + 6.33109e-01	                4290  - 3.54632e+00	                4467  + 6.09386e-01	                4429  + 4.35734e-02	                4364  - 5.92113e-01	
// 8	              360561  - 2.36655e-01	              360888  + 3.35052e-03	              361728  + 2.12060e+00	              360170  - 1.29360e+00	              360253  - 9.98396e-01	
// 9	            22353849  - 4.28707e-01	            22353700  - 4.70966e-01	            22354975  - 1.73570e-01	            22356572  - 6.21969e-03	            22354263  - 3.21716e-01	
//11	                4370  - 4.61300e-01	                4544  + 3.76150e+00	                4486  + 1.13758e+00	                4365  - 5.69179e-01	                4345  - 1.11394e+00	
//12	            85030827  + 1.67438e+00	            85004176  - 2.54851e+00	            85004198  - 2.54090e+00	            85024562  + 3.77635e-01	            85009250  - 1.09436e+00	
//13	          6948718452  + 1.68143e-01	          6948705506  + 6.48966e-02	          6948714551  + 1.31954e-01	          6948675932  - 1.00063e-02	          6948583601  - 1.45846e+00	
//14	        430511445289  + 1.56483e-01	        430511484712  + 2.07629e-01	        430511425086  + 1.33071e-01	        430511443462  + 1.54288e-01	        430511602019  + 4.02525e-01	
//16	              360786  - 1.25250e-02	              360613  - 1.59926e-01	              361817  + 2.57405e+00	              361545  + 1.32616e+00	              361337  + 6.48559e-01	
//17	          6948600869  - 1.00103e+00	          6948636831  - 3.23875e-01	          6948632410  - 3.87053e-01	          6948607239  - 8.53954e-01	          6948631380  - 4.02581e-01	
//18	        567923996348  + 8.18430e-01	        567922545572  - 1.04130e+00	        567922579942  - 9.50296e-01	        567922461572  - 1.28120e+00	        567924052629  + 9.59133e-01	
//19	      35186134868091  + 3.19590e-03	      35186136283078  + 8.70694e-02	      35186131697873  - 2.28401e-01	      35186131841686  - 2.05815e-01	      35186134800648  + 2.03965e-03	
//21	            22356002  - 3.97665e-02	            22353724  - 4.64025e-01	            22352577  - 8.53361e-01	            22355248  - 1.28795e-01	            22355475  - 9.66412e-02	
//22	        430511511736  + 2.46860e-01	        430511553918  + 3.14877e-01	        430511558140  + 3.22139e-01	        430511511670  + 2.46760e-01	        430511502906  + 2.33668e-01	
//23	      35186130192222  - 5.35444e-01	      35186136234028  + 8.22578e-02	      35186136189773  + 7.80340e-02	      35186131769909  - 2.16941e-01	      35186130271100  - 5.16160e-01	
//24	    2179984573881378  + 3.54432e-03	    2179984567799668  - 5.00162e-03	    2179984572443224  + 8.25542e-04	    2179984576906887  + 1.54589e-02	    2179984573791677  + 3.31926e-03	
//
//               6.420 15 =>     4.504753e-02           13.082 15 =>     4.808914e-01           12.241 15 =>     4.145231e-01            6.730 15 =>     5.531633e-02            8.844 15 =>     1.579587e-01
//
//--- Finished -- ReMort trials: 2251799813685248	2^51.0 (out of 2^51) trials -- Ice -- 64 bits: 	ps:   4.505e-02*   4.809e-01    4.145e-01    5.532e-02    1.580e-01   => p <   9.14651e-02   2^54.32 calls, 2^57.32 bytes	2^36.61 bytes/second	used:  19::20:20:42.83

//					uint64_t z = (state += 1UL);
//					uint64_t y = (stream += 3UL);
//					z += y ^ rotate64(y, 11) ^ rotate64(y, 50);
//					y += z ^ rotate64(z, 46) ^ rotate64(z, 21);
//					z += y ^ rotate64(y,  5) ^ rotate64(y, 14);
//					y += z ^ rotate64(z, 25) ^ rotate64(z, 41);
//					z += y ^ rotate64(y, 53) ^ rotate64(y,  3);
//					y += z ^ rotate64(z, 31) ^ rotate64(z, 37);
//					return y ^ z;

					//// Testing with some of the worst possible increments (gammas). Still passes!
					//uint64_t z = (state += 1UL);
					//uint64_t y = (stream += 3UL);
					//z += y ^ rotate64(y, 11) ^ rotate64(y, 50);
					//y += z ^ rotate64(z, 46) ^ rotate64(z, 21);
					//z += y ^ rotate64(y,  5) ^ rotate64(y, 14);
					//y += z ^ rotate64(z, 25) ^ rotate64(z, 41);
					//z += y ^ rotate64(y, 53) ^ rotate64(y,  3);
					//y += z ^ rotate64(z, 31) ^ rotate64(z, 37);
					//return y ^ z;

////JacinthRandom, passes 64TB of PractRand with no anomalies.
// 6	                   2  + 1.36750e+01	                   0  - 2.29283e-01	                   0  - 2.29283e-01	                   0  - 2.29283e-01	                   0  - 2.29283e-01	
// 7	                4441  + 1.51585e-01	                4403  - 3.33245e-02	                4392  - 1.21172e-01	                4443  + 1.75929e-01	                4362  - 6.39342e-01	
// 8	              361197  + 3.27498e-01	              360766  - 2.10857e-02	              360703  - 6.25425e-02	              360457  - 4.35072e-01	              360941  + 2.13489e-02	
// 9	            22354389  - 2.92196e-01	            22354860  - 1.94427e-01	            22355636  - 7.66301e-02	            22355831  - 5.54982e-02	            22354712  - 2.23011e-01	
//11	                4483  + 1.04331e+00	                4431  + 5.70454e-02	                4416  + 1.71509e-04	                4483  + 1.04331e+00	                4373  - 4.02009e-01	
//12	            85001224  - 3.67320e+00	            85020868  + 4.57506e-02	            85028807  + 1.15542e+00	            85013859  - 2.98394e-01	            85012281  - 5.14653e-01	
//13	          6948697884  + 2.66708e-02	          6948705170  + 6.28592e-02	          6948744859  + 5.28296e-01	          6948730308  + 3.05015e-01	          6948613275  - 7.25369e-01	
//14	        430512231145  + 2.53856e+00	        430512204267  + 2.40971e+00	        430512382520  + 3.32696e+00	        430512411952  + 3.49260e+00	        430512304550  + 2.90758e+00	
//16	              361413  + 8.68342e-01	              360938  + 1.99144e-02	              360706  - 6.00695e-02	              361171  + 2.79833e-01	              360515  - 3.17023e-01	
//17	          6948547801  - 2.68021e+00	          6948615290  - 6.84779e-01	          6948680478  - 2.06991e-03	          6948760170  + 8.29039e-01	          6948622456  - 5.49893e-01	
//18	        567923149999  - 4.76956e-02	        567922341912  - 1.66587e+00	        567921576800  - 5.31742e+00	        567922499266  - 1.17047e+00	        567923748196  + 3.31068e-01	
//19	      35186131333876  - 2.90820e-01	      35186132074949  - 1.71681e-01	      35186134297814  - 1.56870e-03	      35186133295191  - 4.35274e-02	      35186130631027  - 4.32655e-01	
//21	            22354822  - 2.01579e-01	            22354646  - 2.36389e-01	            22354893  - 1.88321e-01	            22355066  - 1.57904e-01	            22355832  - 5.53986e-02	
//22	        430512607163  + 4.69315e+00	        430512293918  + 2.85258e+00	        430512220802  + 2.48858e+00	        430512382157  + 3.32494e+00	        430512521530  + 4.14471e+00	
//23	      35186132710200  - 9.44037e-02	      35186131954346  - 1.88943e-01	      35186132679832  - 9.75759e-02	      35186133329249  - 4.11646e-02	      35186132196868  - 1.55071e-01	
//24	    2179984573965209  + 3.76133e-03	    2179984575034484  + 7.09488e-03	    2179984572632590  + 1.07506e-03	    2179984571821645  + 2.37759e-04	    2179984574594330  + 5.59564e-03	
//
//              16.933 15 =>     7.413861e-01            8.651 15 =>     1.473665e-01           13.428 15 =>     5.070419e-01           11.653 15 =>     3.656163e-01           11.425 15 =>     3.480435e-01
//
//--- Finished -- ReMort trials: 2251799813685248	2^51.0 (out of 2^51) trials -- Ice -- 64 bits: 	ps:   7.414e-01    1.474e-01*   5.070e-01    3.656e-01    3.480e-01   => p <   6.89827e-01   2^54.32 calls, 2^57.32 bytes	2^37.18 bytes/second	used:  13::09:26:19.48

					//// Using harmonious number-based increments (R2 sequence), we can elide two lines.
					//// We also probably only need to return y, not y ^ z, but the above test used y ^ z .
					//// This passes PractRand to 64TB with no anomalies, returning either y or y ^ z .
//					uint64_t z = (state += 0xC13FA9A902A6328FUL);
//					uint64_t y = (stream += 0x91E10DA5C79E7B1DUL);
//					z += y ^ rotate64(y,  5) ^ rotate64(y, 14);
//					y += z ^ rotate64(z, 25) ^ rotate64(z, 41);
//					z += y ^ rotate64(y, 53) ^ rotate64(y, 3);
//					y += z ^ rotate64(z, 31) ^ rotate64(z, 37);
//					return y;

					//// Passes 64TB with one anomaly,
					//// [Low4/32]FPF-14+6/16:cross        R=  -2.5  p =1-1.6e-4   unusual
					//// Has a period of 2 to the 128.
					//// 0xAEF17502108EF2DAU is a constant that must be even; 0U works as one extreme.
// 6	                   1  + 2.59071e+00	                   0  - 2.29283e-01	                   0  - 2.29283e-01	                   0  - 2.29283e-01	                   0  - 2.29283e-01	
// 7	                4382  - 2.48596e-01	                4421  + 7.80479e-03	                4320  - 2.04970e+00	                4335  - 1.45427e+00	                4517  + 2.35045e+00	
// 8	              359746  - 3.39738e+00	              361851  + 2.75887e+00	              360270  - 9.42643e-01	              361125  + 2.04681e-01	              360921  + 1.27280e-02	
// 9	            22359936  + 4.00175e-01	            22357793  + 3.21724e-02	            22354381  - 2.94028e-01	            22353511  - 5.27427e-01	            22358618  + 1.25208e-01	
//11	                4434  + 8.06509e-02	                4493  + 1.37341e+00	                4461  + 4.76560e-01	                4340  - 1.27844e+00	                4356  - 7.91898e-01	
//12	            85020921  + 4.82426e-02	            85001836  - 3.42319e+00	            85024822  + 4.13086e-01	            85029843  + 1.40959e+00	            84994963  - 6.73706e+00	
//13	          6948668420  - 3.61563e-02	          6948660158  - 8.36724e-02	          6948766788  + 9.79917e-01	          6948756924  + 7.59644e-01	          6948670656  - 2.66748e-02	
//14	        430511173752  - 3.33602e-04	        430511201040  + 5.44024e-04	        430511451749  + 1.64369e-01	        430511456713  + 1.70561e-01	        430511194598  + 1.82417e-04	
//16	              360659  - 1.04543e-01	              359954  - 2.24083e+00	              361024  + 8.08164e-02	              361085  + 1.48864e-01	              359885  - 2.59792e+00	
//17	          6948762738  + 8.86088e-01	          6948636108  - 3.33823e-01	          6948607096  - 8.57127e-01	          6948763977  + 9.14292e-01	          6948631164  - 4.05876e-01	
//18	        567923576709  + 1.20986e-01	        567922828439  - 4.16138e-01	        567924375331  + 1.98123e+00	        567924084128  + 1.04275e+00	        567923497760  + 5.90824e-02	
//19	      35186148795572  + 5.78148e+00	      35186149671177  + 6.51313e+00	      35186137243744  + 2.08874e-01	      35186137378005  + 2.30075e-01	      35186149016572  + 5.96204e+00	
//21	            22353863  - 4.24839e-01	            22359609  + 3.17460e-01	            22358571  + 1.18272e-01	            22353532  - 5.20996e-01	            22354716  - 2.22212e-01	
//22	        430511461917  + 1.77175e-01	        430511222208  + 3.08981e-03	        430511228335  + 4.21514e-03	        430511451803  + 1.64436e-01	        430511619314  + 4.36666e-01	
//23	      35186137980329  + 3.37797e-01	      35186149654933  + 6.49916e+00	      35186148002992  + 5.15678e+00	      35186137383027  + 2.30888e-01	      35186138055867  + 3.52762e-01	
//24	    2179984552801869  - 1.53618e-01	    2179984541361228  - 4.05735e-01	    2179984553541364  - 1.41453e-01	    2179984563942900  - 2.35087e-02	    2179984552561341  - 1.57682e-01	
//
//              12.198 15 =>     4.086815e-01           24.409 15 =>   1-4.087684e-02           13.869 15 =>     5.412265e-01            9.080 15 =>     1.746725e-01           20.238 15 =>     8.768767e-01
//
//--- Finished -- ReMort trials: 2251799813685248	2^51.0 (out of 2^51) trials -- Ice -- 64 bits: 	ps:   4.087e-01  1-4.088e-02*   5.412e-01    1.747e-01    8.769e-01   => p <   7.22545e-01   2^54.32 calls, 2^57.32 bytes	2^37.29 bytes/second	used:  12::09:02:35.56

//					uint64_t z = (state += 0xC13FA9A902A6328FUL);
//					uint64_t y = (stream += 0x91E10DA5C79E7B1DUL + ((z | 0xAEF17502108EF2DAU - z) >> 63));
//					// z += y ^ rotate64(y,  5) ^ rotate64(y, 14); // This line appears to not be needed here.
//					y += z ^ rotate64(z, 25) ^ rotate64(z, 41);
//					z += y ^ rotate64(y, 53) ^ rotate64(y, 3);
//					y += z ^ rotate64(z, 31) ^ rotate64(z, 37);
//					return y;

//					// Passes 64TB with no anomalies.
//					uint64_t z = (state += 0xC13FA9A902A6328FUL);
//					uint64_t y = (stream += 0x91E10DA5C79E7B1DUL + ((z | 0xAEF17502108EF2DAU - z) >> 63));//0xF1357AEA2E62A9C6U
//					y ^= y >> 32;
//					y *= z | 1U;
//					y ^= y >> 33;
//					y *= rotate64(z, 29) | 1U;
//					y ^= y >> 31;
//					return y;
					//state += 0xC13FA9A902A6328FUL;
					//state ^= state << 7;
					//state ^= state >> 9;

					// we know this one is at least fairly good; adding the rotation is stronger than XORing it.
					//state = -(state ^ (state * state | 7ULL));
					//stream = (stream + 1ULL >> 1) * (stream | 1ULL) ^ 7ULL;
					//return state + rotate64(stream, 7);

					// TruxRandom
					// Passes 64TB with no anomalies, seed 0.
					// Period is guaranteed to be at least 2 to the 64, but only because state has that guarantee.
					// The stream is independent, and could have a much lower period.
					// The names are not very accurate here; this generator doesn't have a known amount of streams because:
					// `stream = (stream + 1ULL >> 1) * (stream | 1ULL) ^ 7ULL;` has an unknown period.
					// If stream is confirmed to have the maximum period, which it does for 8-bit variables, then this generator
					// will have a period of 2 to the 64 and 2 to the 64 possible streams.
					//state = state * 0xC6BC279692B5CC83ULL ^ 0x9E3779B97F4A7C15ULL;
					//stream = (stream + 1ULL >> 1) * (stream | 1ULL) ^ 7ULL;
					//return stream + rotate64(state, 20);

					//state = -(state ^ (state * state | 7ULL));
					//stream = (stream + 1ULL >> 1) * (stream | 1ULL) ^ 7ULL;
					//return state + rotate64(stream, 17);

					//stream = _lzcnt_u64(state = state * 0xC6BC279692B5CC83ULL ^ 0x9E3779B97F4A7C15ULL) - (stream ^ (stream * stream | 0x17ULL));
					//stream = _lzcnt_u64(state = state * 0xD1342543DE82EF95ULL - 0xC6BC279692B5CC83ULL) - (stream ^ (stream * stream | 37ULL));
//uint64_t x = state;
//state = -(x ^ (x * x | 5ULL));
//stream = -(_lzcnt_u64(x) - (stream ^ (stream * stream | 23ULL)));
//return stream - rotate64(x, 20);

// CongaRandom
// Passes 64TB with no anomalies, seed 0.
// Period is guaranteed to be 2 to the 128, with 128 bits of state and 64-bit output.
// 1D-equidistributed over its period.
// Uses an XLCG for one state, and an LCG-like second state that updates dependent on the XLCG (using count-leading-zeros).
// The number of rotations this does may be able to be cut down.
//uint64_t x = state;
//state = x * 0xC6BC279692B5CC83ULL ^ 0x9E3779B97F4A7C15ULL;
//x = rotate64(x, 20) + (stream = stream * 0xD1342543DE82EF95ULL + _lzcnt_u64(x));
//return x ^ rotate64(x, 29) ^ rotate64(x, 47);

// SpinyRandom
// Passes 64TB with no anomalies, seed 0.
// Period is guaranteed to be 2 to the 128, with 128 bits of state and 64-bit output.
// 1D-equidistributed over its period.
// Uses an LCG for state, and an LCG-like stream that depends on the first LCG's value (using count-leading-zeros).
// Also uses a construction like the round function for the Speck cipher (without XORing in a key) for output.
// Only requires two rotates, an add, and a XOR for output, none of which depend on mutating state.
// Surprisingly fast in PractRand.
//uint64_t x = state, y = stream;
//state = x * 0x369DEA0F31A53F85ULL + 0xC6BC279692B5CC83ULL;
//stream = y * 0xD1342543DE82EF95ULL + _lzcnt_u64(x);
//return (rotate64(x, 20) + y) ^ rotate64(y, 47);

// Fails at 8TB with [Low1/64]TMFn(2+7):wl, [Low1/64]TMFn(2+8):wl
//uint64_t x = state, y = stream;
//state = x * 0x369DEA0F31A53F85ULL + 0xC6BC279692B5CC83ULL;
//stream = y * 0xD1342543DE82EF95ULL + _lzcnt_u64(x);
//return rotate64(x, 20) + y;

// Fails a lot of tests right away...
//uint64_t x = state, y = stream;
//state = x + 0xD1342543DE82EF95ULL;
//stream = y + (x ^ _lzcnt_u64(x));
//return rotate64(x, 20) + (y * (y + 0xC6BC279692B5CC83ULL) >> 1);

// WeaselRandom
// Passes 64TB of PractRand with no anomalies, seed 0.
// Period is guaranteed to be 2 to the 128, with 128 bits of state and 64-bit output.
// 1D-equidistributed over its period.
// Uses a Weyl sequence and a variable-increment sequence dependent on the other sequence.
// Runs a round of the Speck cipher round function (with key 0),
// multiplies by an LCG constant (see Vigna and Steele), then xor-rotate-xor-rotate and return.
// Only one multiply means this should be pretty fast.
//uint64_t x = state, y = stream, z = (rotate64(x, 3) + y ^ rotate64(y, 56)) * 0xD1342543DE82EF95ULL;
//state = x + 0xC6BC279692B5CC83ULL;
//stream = y + (x ^ _lzcnt_u64(x));
//return z ^ rotate64(z, 29) ^ rotate64(z, 40);
uint64_t z = (state ^ state >> 30) * 0xF1357AEA2E62A9C5UL;
state += 0x9E3779B97F4A7C15UL + (stream = stream >> 1 ^ (-(stream & 1UL) & 0x9C82435DE0D49E35UL));
return z ^ z >> 29;
				}
				std::string tiptoe64::get_name() const { return "tiptoe"; }

				uint64_t mmi(const uint64_t a) {
					uint64_t x = 2 ^ a * 3;
					x *= 2 - a * x;
					x *= 2 - a * x;
					x *= 2 - a * x;
					x *= 2 - a * x;
					return x;
				}

				// uint64_t fixGamma(uint64_t gamma) {
				// 	uint64_t inverse = mmi(gamma |= 1UL), add = 0UL;
				// 	while (abs((long long)__popcnt64(gamma) - 32) > 8
				// 		|| abs((long long)__popcnt64(gamma ^ gamma >> 1) - 32) > 8
				// 		|| abs((long long)__popcnt64(inverse) - 32) > 8
				// 		|| abs((long long)__popcnt64(inverse ^ inverse >> 1) - 32) > 8) {
				// 		inverse = mmi(gamma = gamma * 0xD1342543DE82EF95L + (add += 2L));
				// 	}
				// 	return gamma;
				// }

				void tiptoe64::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
					walker->handle(stream);
					//stream |= 1ULL;
					//stream = (stream ^ UINT64_C(0x369DEA0F31A53F85)) * UINT64_C(0x6A5D39EAE116586D) + (state ^ state >> 17) * UINT64_C(0x9E3779B97F4A7C15);
					//stream = stream << 3 ^ UINT64_C(0x369DEA0F31A53F89);
					printf("Seed is 0x%016llX, Stream is 0x%016llX\r\n", state, stream);
					//printf("Seed is 0x%I64X\r\n", state);
				}


				// UINT64_C(0x369DEA0F31A53F85); // UINT64_C(0x632BE59BD9B4E019); UINT64_C(0x6C8E9CF570932BD5) UINT64_C(0xDB4F0B9175AE2165)
				Uint64 linnormA::raw64() {
					//(state = state * UINT64_C(0x41C64E6D) + UINT64_C(1));
					//// DiverRNG, passes 32TB with one anomaly
					//uint64_t z = (state = (state ^ UINT64_C(0x6C8E9CF570932BD5)) * UINT64_C(0xC6BC279692B5CC83));
					//z = rotate64(z, 27) * UINT64_C(0xDB4F0B9175AE2165);
					//return z ^ (z >> 25u);
					
					//uint64_t z = (rotate64(state, 21) ^ rotate64(state, 41) ^ state++) * UINT64_C(0xC6BC279692B5CC83);
					//z = (z ^ z >> 29) * UINT64_C(0xAEF17502108EF2D9);
					//return z ^ z >> 30;
					//z = rotate64(z, 29) * UINT64_C(0xD1B54A32D192ED03);
					//z = (z ^ z >> 29u) * UINT64_C(0xAEF17502108EF2D9);
					////Early version of DiverRNG.determine(), passes 32TB no anomalies, not very fast
					//uint64_t z = (state++ ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xC6BC279692B5CC83); // golden ratio, neely
					//z = (z ^ z >> 28u ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0xD2B74407B1CE6E93); // phi, used in philox
					//return z ^ rotate64(z, 19) ^ rotate64(z, 37); // primes close to 1/3 and 3/5 of 64

					//uint64_t z = (rotate64(state, 19) ^ rotate64(state, 37) ^ state++) * UINT64_C(0xC6BC279692B5CC83);
					//z = (z ^ z >> 30) * UINT64_C(0xDB4F0B9175AE2165); // phi, used in philox
					//return z ^ z >> 29; // primes close to 1/3 and 3/5 of 64

					// fails at 16TB suddenly wih BRank issue
					//uint64_t z = (rotate64(state, 21) ^ rotate64(state, 33) ^ state++ ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xC6BC279692B5CC83); // generalized golden ratio, neely
					//z = (z ^ z >> 28 ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0xD2B74407B1CE6E93); // phi, used in philox
					//return z ^ z >> 28;
					
					//uint64_t z = (state++ ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xD1B54A32D192ED03);
					////uint64_t z = state++;
					////z = rotate64(z, 3) * UINT64_C(0xC6BC279692B5CC83);
					//z = (z ^ z >> 26u) * UINT64_C(0x9E3779B97F4A7C15) + UINT64_C(0x369DEA0F31A53F85);
					//return (z ^ z >> 25u);

					//uint64_t z = state++;
					//z = (z ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0xC6BC279692B5CC83);
					//z = (z ^ z >> 27u) * UINT64_C(0xAEF17502108EF2D9);
					//return (z ^ z >> 25u) * UINT64_C(0x2545F4914F6CDD1D);        
					
					////passes 32TB no anomalies
					//uint64_t z = (rotate64(state, 5) ^ rotate64(state, 8) ^ state++) * 0x7FFFFFFF;
					//z = (z ^ z >> 30) * 0xBF58476D1CE4E5B9;
					//z = (z ^ z >> 27) * 0x94D049BB133111EB;
					//return z ^ z >> 31;
					
					////passes 32TB no anomalies
					//uint64_t z = state++;
					//z = ((z << 32) ^ rotate64(z, 29)) * UINT64_C(127);//UINT64_C(0x7FFFFFFF);
					//z = (z ^ z >> 30) * UINT64_C(0xBF58476D1CE4E5B9);
					//z = (z ^ z >> 27) * UINT64_C(0x94D049BB133111EB);
					//return z ^ z >> 31;
					
					//// passes 32TB no anomalies
					//uint64_t z = state;
					//z = ((z << (((++state) & 31) + 5)) ^ rotate64(z, 3)) * UINT64_C(0xAEF17502108EF2D9);
					//z = (z ^ (z >> ((z >> 60) + 16))) * UINT64_C(0x369DEA0F31A53F85);
					//return z ^ z >> 27;

					//// passes 32TB no anomalies
					//uint64_t z = ++state;
					//z = ((z << ((z & 31) + 5)) ^ z ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xD1B54A32D192ED03);
					//z = (z ^ (z >> ((z >> 60) + 16))) * UINT64_C(0x369DEA0F31A53F85);
					//return z ^ z >> 27;
					
					
					//// passes 32TB with one anomaly (16GB, [Low16/64]DC6-9x1Bytes-1, unusual)
					//uint64_t z = ++state;
					//z ^= rotate64(z, 21) ^ rotate64(z, 41);
					//z += UINT64_C(0x9E3779B97F4A7C15);
					//z = (z ^ z >> 30) * UINT64_C(0xBF58476D1CE4E5B9);
					//z = (z ^ z >> 27) * UINT64_C(0x94D049BB133111EB);
					//return z ^ z >> 31;

					//// passes 32TB with 4 anomalies, worsening toward the end
					//uint64_t z = ++state;
					//z = ((z << 11) ^ z ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xD1B54A32D192ED03);
					//z ^= z >> 25;
					//z = ((z >> 27) ^ z ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xD1B54A32D192ED03);
					//return ((z >> 29) ^ z);// ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xD1B54A32D192ED03);

					//uint64_t z = ++state;
					//z = ((z << 21) ^ z ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xD1B54A32D192ED03);
					//z = (z ^ (z >> ((z >> 60) | 16))) * UINT64_C(0x369DEA0F31A53F85);
					//return z ^ z >> 27;

					//uint64_t z = ++state;
					//z = (z ^ z << 11 ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xD1B54A32D192ED03);
					//z ^= z >> 25;
					//z = (z ^ z >> 27 ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xD1B54A32D192ED03);
					//return (z ^ z >> 29);// ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xD1B54A32D192ED03);

					// passes 32TB no anomalies, not super-fast but faster than variable-shift unary hashes
					// similar to http://mostlymangling.blogspot.com/2018/07/on-mixing-functions-in-fast-splittable.html 
					//uint64_t z = ++state;
					//z = (z ^ rotate64(z, 12) ^ rotate64(z, 43) ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xE0A28963);
					//z = (z ^ rotate64(z, 52) ^ rotate64(z, 21) ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0x81383173);
					//return z ^ z >> 28;
					//z = (z ^ rotate64(z, 12) ^ rotate64(z, 43) ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0x4823A80B2006E21B); //0xE0A28963
					//z = (z ^ rotate64(z, 52) ^ rotate64(z, 21) ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0x4C8148C08017D353); //0x81383173
					// * UINT64_C(0x9043);//UINT64_C(0x9fb21c651e98df25);

					////passes 2TB without anomalies, then has issues at 4TB:
					//// mildly suspicious : [Low1/64]BCFN(2+3,13-0,T)
					//// unusual           : [Low1/64]BCFN(2+5,13-0,T)
					//uint64_t z = ((++state) ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0xC6BC279692B5CC83);
					//z = (z ^ z >> 27 ^ z >> 12) * UINT64_C(0xAEF17502108EF2D9);
					//return (z ^ z >> 25);

					////where we left off
					//uint64_t z = ((++state) ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0xC6BC279692B5CC83);
					//z = (z ^ z >> 27 ^ z >> 37 ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xD1B54A32D192ED03);
					//return (z ^ z >> 25 ^ z >> 39);


////passes 32TB no anomalies!
//uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
//z = (z ^ rotate64(z, 23) ^ rotate64(z, 41)) * 0x369DEA0F31A53F85ULL;
//return z ^ z >> 34 ^ z >> 26;

//uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
//z ^= z >> 25u ^ z >> 13u;
//z *= 0x2545F4914F6CDD1DULL;
//return z ^ z >> 26u ^ z >> 31u;

//fails at 4TB
//uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
//z = (z ^ z >> 30 ^ z >> 5 ^ 0xDB4F0B9175AE2165ULL) * 0xC6BC279692B5CC83ULL;
//return z ^ z >> 28u;
////return z ^ z >> 31u ^ z >> 23u;

//VelvetRNG, should also work as a simple determine()
//passes 32TB, no anomalies, could be fast in Java, seems fast in MSVC...
//uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
//z = (z ^ rotate64(z, 38) ^ rotate64(z, 23)) * 0x369DEA0F31A53F85ULL;
//return z ^ z >> 39 ^ z >> 26;


//z = ((z << ((++state & 31u) + 5u)) ^ rotate64(z, 4)) * UINT64_C(0xAEF17502108EF2D9);
// 0xE7037ED1A0B428DBULL;//0xAEF17502108EF2D9ULL;
////return (z ^ z >> 47 ^ z >> 23);

////passes 32TB no anomalies, not equidistributed
//uint64_t y = state;
//uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
////z = rotate64(z, 27u);
//y ^= y >> 26u;
//z ^= z >> 28u;
//y *= 0xAEF17502108EF2D9ULL;
//z *= 0x2545F4914F6CDD1DULL;
//return y ^ z ^ y >> 28u ^ z >> 26u;


////passes at least 8TB without anomalies, but can't produce 3/4 of all uint64_t
//uint64_t z = (state *= 0x59071D96D81ECD35ULL);
//z = (z ^ z >> 26u ^ 0xDB4F0B9175AE2165ULL) * 0xC6BC279692B5CC83ULL;
//return z ^ z >> 23u;

//uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
//z = (z ^ z >> 53u ^ z >> 29u ^ z >> 23u ^ z >> 42u) * 0xDB4F0B9175AE2165ULL;
//return z ^ z >> 25u;
//0x818102004182A025ULL
//return z ^ z >> 28u;

//0x369DEA0F31A53F85ULL 0xC6BC279692B5CC83ULL 0x6C8E9CF570932BD5ULL 0x59071D96D81ECD35ULL
					//z = ((z << ((++state & 31u) + 5u)) ^ rotate64(z, 4)) * UINT64_C(0xAEF17502108EF2D9);
					//z = ((z >> 30) ^ rotate64(z, 37)) * UINT64_C(0x369DEA0F31A53F85);
					//z = ((z >> 26) ^ z) * UINT64_C(0x9E3779B97F4A7C15);
					//return z ^ z >> 26;
					//uint64_t z = state++;
					//z = (rotate64(z, 21) ^ rotate64(z, 35) ^ z ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xC6BC279692B5CC83);
					//z = (z ^ z >> 26) * UINT64_C(0xD1B54A32D192ED03);
					// 0xC6BC279692B5CC83  0x6C8E9CF570932BD5
                    //uint64_t z = (++state ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0x4823A80B2006E21B);
					//z = (z ^ z << 7 ^ z >> 23) * UINT64_C(0x4C8148C08017D353);
					//return z ^ z >> 28;

					//return z ^ z >> 26;

					//z = (z ^ z >> 28 ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0xD2B74407B1CE6E93); // phi, used in philox

					//uint64_t z = state++;
					//z = (z << 4) - rotate64(z, 2);
					//z = ((z << 8) - rotate64(z, 6)) * UINT64_C(0xC6BC279692B5CC83);
					//z = (z ^ z >> 27) * UINT64_C(0x9E3779B97F4A7C15);//) * UINT64_C(0xD2B74407B1CE6E93); // phi, used in philox
					//return z ^ z >> 29;
					//z = ((z ^ rotate64(z, 21) ^ rotate64(z, 33) ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xC6BC279692B5CC83));
						// ^ ((rotate64(z, 39) ^ rotate64(z, 11) ^ rotate64(z, 47) ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0xD2B74407B1CE6E93));
					//z = (z ^ z >> 29 ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0xD2B74407B1CE6E93);
					//z = rotate64(z, 27) * UINT64_C(0x9E3779B97F4A7C15);
					//z = rotate64(z, 29) * UINT64_C(0xD2B74407B1CE6E93);

					//z = (rotate64(z, 43) ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0xD2B74407B1CE6E93);
					//z = (rotate64(z, 51) ^ UINT64_C(0x6C8E9CF570932BD5)) * UINT64_C(0xD1B54A32D192ED03);


					//return (z ^ z >> 30u);

					////Fails TMFn after 4 to 16 TB; this is very close to  DiverRNG but uses an LCG, not an XLCG.
					uint64_t z = (state = state * 0xD1342543DE82EF95UL + 1UL);
					z = rotate64(z, 41) * 0xDB4F0B9175AE2165UL;
					return z ^ (z >> 28u);

				}
				std::string linnormA::get_name() const { return "linnormA"; }
				void linnormA::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
					printf("Seed is 0x%016llX\r\n", state);// stream);
				}


				linnormB::linnormB(int rotation, int x)
				{
					R = rotation;
					X = (uint64_t)x;
				}
				//uint64_t reverse(uint64_t n)
				//{
				//	n = _byteswap_uint64(n);
				//	n = ((n & UINT64_C(0xaaaaaaaaaaaaaaaa)) >> 1) | ((n & UINT64_C(0x5555555555555555)) << 1);
				//	n = ((n & UINT64_C(0xcccccccccccccccc)) >> 2) | ((n & UINT64_C(0x3333333333333333)) << 2);
				//	n = ((n & UINT64_C(0xf0f0f0f0f0f0f0f0)) >> 4) | ((n & UINT64_C(0x0f0f0f0f0f0f0f0f)) << 4);
				//	return n;
				//}
				Uint64 linnormB::raw64() {
					//UINT64_C(0x9E3779B97F4A7C15)
					//uint64_t z = (state = state * UINT64_C(0x41C64E6D) + UINT64_C(3));

					//uint64_t z = (state = (state ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0x41C64E6B));
					//z = (z ^ z >> 27u) * UINT64_C(0xAEF17502108EF2D9);
					//return z ^ (z >> 30u);


					//static const uint64_t P1 = 11400714785074694791ULL;  /* 0b1001111000110111011110011011000110000101111010111100101010000111 */
					//static const uint64_t P2 = 14029467366897019727ULL;  /* 0b1100001010110010101011100011110100100111110101001110101101001111 */
					//static const uint64_t P3 = 1609587929392839161ULL;   /* 0b0001011001010110011001111011000110011110001101110111100111111001 */
					//static const uint64_t P4 = 9650029242287828579ULL;   /* 0b1000010111101011110010100111011111000010101100101010111001100011 */
					//static const uint64_t P5 = 2870177450012600261ULL;   /* 0b0010011111010100111010110010111100010110010101100110011111000101 */

					//uint64_t input = (state += 1);
					//uint64_t hash = P5 + 8;
					//uint64_t hash = (state++);
					//hash ^= rotate64(hash, 21) ^ rotate64(hash, 41);
					//hash += P1;
					//hash ^= (hash << ((hash & 31) + 5));
					//hash *= P5;
					//input *= P2;
					//input = rotate64(input, 31);
					//input *= P1;
					//hash ^= input;
					//hash = rotate64(hash, 27) * P1 + P4;

					//hash ^= hash >> 33;
					//hash *= P2;
					//hash ^= hash >> 29;
					//hash *= P3;
					//hash ^= hash >> 32;
					//return hash;

					//uint64_t s0 = ++state * UINT64_C(0x369DEA0F31A53F85);
					//uint64_t s1 = rotate64(s0, 31) + UINT64_C(0xAEF17502108EF2D9);
					//s0 = (s0 ^ rotate64(s0, 24) ^ rotate64(s0, 51));
					//s1 = (s1 ^ s1 >> 29) * UINT64_C(0xDB4F0B9175AE2165);// UINT64_C(0x7FF8A3ED) + UINT64_C(0x9E3779B97F4A7C15);
					//return (s0 ^ s1 ^ s1 >> 29) * UINT64_C(0x41C64E6D);

					//uint64_t z = ((state += UINT64_C(0x9E3779B97F4A7C15)) ^ UINT64_C(0x369DEA0F31A53F85) * UINT64_C(0x61C8864680B583EB));
					//z = (z ^ z >> 30u) * UINT64_C(0x632BE59BD9B4E019);// UINT64_C(0x369DEA0F31A53F85); // UINT64_C(0xAEF17502108EF2D9);
					//return z ^ (z >> 31u);

					//// passes 32TB with one anomaly, unusual at only 2GB: [Low16/64]BCFN(2+2,13-2,T)
					// used command line:
					// RNG_test.exe linnormB -multithreaded -seed 0xa6948409
					//uint64_t z = ++state;
					//z = (z ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0x4823A80B2006E21B);
					//z = (z ^ rotate64(z, 52) ^ rotate64(z, 21) ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0x81383173);
					//return z ^ z >> 28;
					//// passes 32TB with two anomalies, unusual at 8GB and 16GB, both Gap-16:B, earlier anomaly only on low 1/64
					//uint64_t z = ++state;
					//z = (z ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0x4823A80B2006E21B);
					//z = (z ^ rotate64(z, 12) ^ rotate64(z, 43) ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0x81383173);
					//return z ^ z >> 28;

					//uint64_t z = ++state;
					//z = rotate64(z, 1);
					//z = (z ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0x4823A80B2006E21B);
					//z = (z ^ rotate64(z, 52) ^ rotate64(z, 21) ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0x81383173);
					////z = (z ^ rotate64(z, 12) ^ rotate64(z, 43) ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0x81383173);
					//return z ^ z >> 28;

					//uint64_t s = state++;
					//s = reverse(s);
					//s = rotate64(s, 43);
					//s = (s ^ rotate64(s, 41) ^ rotate64(s, 17) ^ 0xD1B54A32D192ED03UL) * 0x2127599BF4325C37UL;
					//s = (s ^ rotate64(s, 42) ^ rotate64(s, 18)) * 0x880355F21E6D1965UL;
					//return (s ^ s >> 28);
					//v *= 0x2127599BF4325C37UL;
					//v ^= v >> 47;
					//return v;
					//v = (v ^ rotate64(v, 41) ^ rotate64(v, 17) ^ 0xD1B54A32D192ED03UL) * 0x2127599BF4325C37UL; //multiplier from Fast-Hash by Zilong Tan, https://github.com/ZilongTan/fast-hash
					//v = (v ^ rotate64(v, 42) ^ rotate64(v, 18)) * 0x880355F21E6D1965UL; // multiplier also from Fast-Hash
					//v = (v ^ rotate64(v, 41) ^ rotate64(v, 18) ^ 0xD1B54A32D192ED03UL) * 0xAEF17502108EF2D9UL;
					//v = (v ^ rotate64(v, 42) ^ rotate64(v, 19)) * 0xDB4F0B9175AE2165UL;

					//v = (v ^ v >> 20 ^ v << 28 ^ 0xD1B54A32D192ED03UL) * 0xAEF17502108EF2D9UL;
					//v = (v ^ v >> 27 ^ v << 21) * 0xDB4F0B9175AE2165UL;
					//v = (v ^ v >> 23 ^ v >> 19 ^ v << 18) * 0xDB4F0B9175AE2165UL;
					//return s ^ s >> 28;

					// 0xA24BAED4963EE407UL;

					// 0x9FB21C651E98DF25UL;

					//uint64_t z = ((++state) ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0x4823A80B2006E21B);
					//z = (z ^ z >> 23 ^ z >> 27) * UINT64_C(0xDB4F0B9175AE2165);// UINT64_C(0xD1B54A32D192ED03);// ^ UINT64_C(0x9E3779B97F4A7C15) UINT64_C(0x81383173)
					//return z ^ z >> 25 ^ z >> 21;

					//uint64_t z = (++state ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0x4823A80B2006E21B);
					//z = (z ^ z >> 12 ^ z >> 43 ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0x81383173);
					//return z ^ z >> 28;

					//uint64_t z = (++state ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0x4823A80B2006E21B);
					//z = (z ^ rotate64(z, 52) ^ rotate64(z, 21) ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0x81383173);
					//return z ^ z >> 28;

					//uint64_t z = (++state ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0x4823A80B2006E21B);
					//z = (z ^ rotate64(z, 52) ^ rotate64(z, 21)) * UINT64_C(0x9E3779B97F4A7C15);
					//return z ^ z >> 28;

					//uint64_t z = ++state;
					////passes 32TB no anomalies when R == 0 and not reversed, seed=0xbd75081b
					//uint64_t v = ++state;
					// comment out next line to disable reversal and rotation
					//v = reverse(rotate64(v, R));
					// comment out next line to disable reversal but not rotation
					//v = rotate64(v, R);

                    //v ^= rotate64(v, 39) ^ rotate64(v, 14);
                    //v *= 0xAEF17502108EF2D9UL;//0xA24BAED4963EE407UL; // second number was used by Evensen
                    //v ^= rotate64(v, 40) ^ rotate64(v, 15);
                    //v *= 0xDB4F0B9175AE2165UL;//0x9FB21C651E98DF25UL; // second number was used by Evensen
                    //return v ^ v >> 28;

					//uint64_t s = state++;
					//s = reverse(s); // this line was commented in and out for different test types
					//s = rotate64(s, R);

					/*
					s = (s ^ (s << 39 | s >> 25) ^ (s << 14 | s >> 50)) * 0xAEF17502108EF2D9UL + 0xD1B54A32D192ED03UL;
					s = (s ^ (s << 40 | s >> 24) ^ (s << 15 | s >> 49)) * 0xDB4F0B9175AE2165UL;
					return s ^ s >> 28;
					*/
					//multipliers from Fast-Hash by Zilong Tan, https://github.com/ZilongTan/fast-hash
					//s = (s ^ (s << 41 | s >> 23) ^ (s << 18 | s >> 46) ^ 0xD1B54A32D192ED03UL) * 0x2127599BF4325C37UL;
					//s = (s ^ s >> 43 ^ s >> 34 ^ s >> 19) * 0x880355F21E6D1965UL;
					//return s ^ s >> 28;

					//s = (s ^ rotate64(s, 41) ^ rotate64(s, 17) ^ 0xD1B54A32D192ED03UL) * 0x2127599BF4325C37UL;
					//s = (s ^ rotate64(s, 40) ^ rotate64(s, 15)) * 0x880355F21E6D1965UL;
					//return (s ^ s >> 28);



					//IT WORKS!
//					uint64_t s = state++;
//					//s = reverse(s); // this line was commented in and out for different test types
//					s = rotate64(s, R);
//					//s = ~s; // this line was commented in and out for different test types
//					s = (s ^ rotate64(s, 41) ^ rotate64(s, 17) ^ 0xD1B54A32D192ED03UL) * 0xAEF17502108EF2D9UL;
//					s = (s ^ s >> 43 ^ s >> 31 ^ s >> 23) * 0xDB4F0B9175AE2165UL;
//					return s ^ s >> 28;


					////PulleyRNG
					////WORKS, 32 TB no anomalies!
					//uint64_t s = state++;
					//s = (s ^ rotate64(s, 41) ^ rotate64(s, 17)) * 0x369DEA0F31A53F85ULL;
					//s = (s ^ s >> 25 ^ s >> 37) * 0xDB4F0B9175AE2165UL;
					//return s ^ s >> 28;


					// 0xBE21F44C6018E14DULL * 0x369DEA0F31A53F85ULL == 1ULL

					////ToppingRNG, passes 32TB with no anomalies.
					////Literally just a different start to SplitMix64 that allows it to pass tests with gamma=1.
					//uint64_t s = state++;
					//s = (s ^ rotate64(s, 23) ^ rotate64(s, 47)) * 0xEB44ACCAB455D165ULL;
					//s = (s ^ s >> 30) * 0xBF58476D1CE4E5B9ULL;
					//s = (s ^ s >> 27) * 0x94D049BB133111EBULL;
					//return (s ^ s >> 31);

//0x2127599BF4325C37ULL;//0x9E3779B97F4A7C15ULL;//0xEB44ACCAB455D165ULL;
//(s ^ s >> 29 ^ s >> 43 ^ s << 7 ^ s << 53)
//					s = (s ^ s << 23 ^ s << 47 ^ s >> 26) * 0x369DEA0F31A53F85ULL;

					uint64_t s = state++;

					//s = reverse_bits64(s); // this line was commented in and out for different test types
					s = rotate64(s, R);
					uint64_t x = s ^ -X; // change X as a parameter to the generator to 0 to disable any bit flips, 1 to flip all, other (32-bit max, assigned to 64-bit) numbers are allowed
					
					x ^= rotate64(x, 39) ^ rotate64(x, 14); // rotl; they correspond to right rotations of 25 and 50
					x *= 0x3C79AC492BA7B653UL;
					x ^= x >> 33;
					x ^= x >> 13;
					x *= 0x1C69B3F74AC4AE35UL;
					x ^= x >> 27;
					return x;

////works well, except for BRank issues on rotation 47, reversed and all flipped
//					x ^= rotate64(x, 39) ^ rotate64(x, 17);
//					x *= 0x3C79AC492BA7B653UL;
//					x ^= x >> 33;
//					x ^= x >> 11;
//					x *= 0x1C69B3F74AC4AE35UL;
//					x ^= x >> 27;
//					return x;



					//s = (s ^ rotate64(s, 23) ^ rotate64(s, 47)) * 0xE7037ED1A0B428DBULL;
					//s = (s ^ rotate64(s, 39) ^ rotate64(s, 14)) * 0xEB44ACCAB455D165ULL;//0xA0761D6478BD642FULL;
					//s = (s ^ s >> 25 ^ s >> 37) * 0xA0761D6478BD642FULL;
					//return s ^ s >> 28;

					////PulleyRNG
					//s = (s ^ rotate64(s, 41) ^ rotate64(s, 17)) * 0x369DEA0F31A53F85ULL;
					//s = (s ^ s >> 25 ^ s >> 37) * 0xDB4F0B9175AE2165UL;
					//return s ^ s >> 28;

//					s = (s ^ rotate64(s, 23) ^ rotate64(s, 51)) * 0xDB4F0B9175AE2165UL;
//					s = (s ^ s >> 28) * 0xEB44ACCAB455D165ULL;
//					return (s ^ s >> 26);
					//s = (s ^ s >> 23 ^ s << 42);// * 0xEB44ACCAB455D165ULL;
//					s = (s ^ rotate64(s, 23) ^ rotate64(s, 47));
//					s = (s ^ s >> 30) * 0xBF58476D1CE4E5B9ULL;
//					s = (s ^ s >> 27) * 0x94D049BB133111EBULL;
//					return (s ^ s >> 26);

//					s = (s ^ rotate64(s, 23) ^ rotate64(s, 47)) * 0xEB44ACCAB455D165ULL;
//					s = (s ^ s >> 30) * 0xBF58476D1CE4E5B9ULL;
//					s = (s ^ s >> 27) * 0x94D049BB133111EBULL;
//					return (s ^ s >> 31);


					//test correlation between gamma 3 and gamma (multiplicative inverse of 3, mod 2 to the 64)
					//uint64_t s = state * 0xAAAAAAAAAAAAAAABULL, t = state++ * 3ULL;
					//s = (s ^ rotate64(s, 41) ^ rotate64(s, 17)) * 0x369DEA0F31A53F85ULL;
					//s = (s ^ s >> 25 ^ s >> 37) * 0xDB4F0B9175AE2165UL;
					//t = (t ^ rotate64(t, 41) ^ rotate64(t, 17)) * 0x369DEA0F31A53F85ULL;
					//t = (t ^ t >> 25 ^ t >> 37) * 0xDB4F0B9175AE2165UL;
					//return s ^ s >> 28 ^ t ^ t >> 28;

					//return s ^ s >> 43 ^ s >> 21;

//					s = (s ^ s >> 43 ^ s >> 31 ^ s >> 23) * 0xDB4F0B9175AE2165UL;
//					return s ^ s >> 28;

    //// score = 1.8856880829987859
    //x ^= x >> 31;
    //x *= UINT64_C(0x1efc5d7c770fc6e9);
    //x ^= x >> 28;
    //x *= UINT64_C(0xda03b996b1182b23);
    //x ^= x >> 30;

    //// score = 1.8849632548546522
    //x ^= x >> 32;
    //x  = ~x * UINT64_C(0xFD3865EC82F016D1);
    //x ^= x >> 28;
    //x *= UINT64_C(0x4E82519D4F727C6D);
    //x ^= x >> 30;
	
	//// score = 1.868895541670708
	//x ^= x >> 33;
    //x  = ~x;
    //x *= UINT64_C(0x906c1c3e34da24ed);
    //x ^= x >> 28;
    //x *= UINT64_C(0x91fcf1c613ff4f85);
    //x ^= x >> 30;
    //return x;
	/* // java
	long x = (state += 0x9E3779B97F4A7C15L);
	x = (x ^ x >>> 33 ^ -1L) * 0x906C1C3E34DA24EDL;
	x = (x ^ x >>> 28) * 0x91FCF1C613FF4F85L;
	return x ^ x >>> 30;
	*/
					//s = (s ^ (s << 39 | s >> 25) ^ (s << 14 | s >> 50) ^ 0xD1B54A32D192ED03UL) * 0xAEF17502108EF2D9UL;
					//s = (s ^ (s << 40 | s >> 24) ^ (s << 15 | s >> 49)) * 0xDB4F0B9175AE2165UL;
					//return s ^ s >> 28;

					//s = (s ^ rotate64(s, 41) ^ rotate64(s, 17) ^ 0xD1B54A32D192ED03UL) * 0x2127599BF4325C37UL;
					//s = (s ^ s >> 34 ^ s >> 19) * 0x880355F21E6D1965UL;
					//return s ^ s >> 43 ^ s >> 28;


					//tested all 64 rotations and all 64 reversed rotations of a counter with -tf 2, up to 1TB for each test.
					//results are in the etc folder.
					//all passed using the same rotations as Pelle Evensen's rrxmrrxmsx_0, using the equivalent of 64-bit right rotations
					//[25,50], then later [24,49], ending with a xorshift by 28, only changing the multipliers after each rotation pair.
					//(the code I tested used the equivalent 64-bit left rotations [39,14] and later [40,15].)
					//the multipliers used are 0xAEF17502108EF2D9UL and 0xDB4F0B9175AE2165UL, in that order.
					//0xAEF17502108EF2D9UL (which replaces 0xA24BAED4963EE407UL) is a multiplier used in PCG-Random.
					//0xDB4F0B9175AE2165UL (which replaces 0x9FB21C651E98DF25UL) is 2 to the 64 divided by the fourth harmonious number,
					//which is the solution to pow(x, 5) = x + 1 (harmonious numbers include the golden ratio and generalize it).
					//some intermediate results were mildly suspicious, and one was suspicious, but the suspicious one cleared up over
					//time, soon having no current anomalies at the 1TB mark. The seed per test was random instead of sequential as in the
					//original test done by Evensen; because the state of the generator always starts as a 32-bit number for each test
					//and the state is only incremented 2^37 times per test, very large internal states weren't tested, but the rotations
					//should have changed the range of bits being operated on to match all possible high and low affected bits.

					//z = (z ^ rotate64(z, 42) ^ rotate64(z, 21) ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xD1B54A32D192ED03);
					////z = (z ^ z << 6 ^ z >> 21 ^ z >> 37 ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xD1B54A32D192ED03);
					//z = (z ^ z << 6 ^ z >> 21 ^ z >> 37) * UINT64_C(0xC6BC279692B5CC83);
					//return (z ^ z << 6 ^ z >> 21 ^ z >> 37);
					
					//z = (z ^ z >> 26 ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0x81383173);
					//return z ^ z >> 28;

					//z = (z ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0x4823A80B2006E21B);
					//z = (z ^ z << 6 ^ z >> 21 ^ z >> 37) * UINT64_C(0xC6BC279692B5CC83);
					//return z ^ z >> 28;

					//uint64_t z = (++state ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xD1B54A32D192ED03);
					//z = (z ^ z << 5 ^ z >> 21 ^ z >> 23 ^ z >> 47) * UINT64_C(0xC6BC279692B5CC83);
					//return z ^ z >> 31;
					//uint64_t z = ++state;
					//z = rotate64(z, R);
					//z = (z ^ rotate64(z, 21) ^ rotate64(z, 42)) * UINT64_C(0x4823A80B2006E21B);
					//uint64_t z = (++state ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xD1B54A32D192ED03);
					//z = (z ^ z << 6 ^ z >> 21 ^ z >> 37) * UINT64_C(0x4C8148C08017D353);
					//return z ^ z >> 28;

				}
				std::string linnormB::get_name() const { 
					std::ostringstream str;
					str << "linnormB(" << R << "," << X << ")";
					return str.str();
				}
				void linnormB::walk_state(StateWalkingObject *walker) {
					//walker->handle(state);
					state = 0UL;
					printf("Seed is 0x%016llX\r\n", state);// stream);
				}
				linnorm32::linnorm32(int rotation, int x)
				{
					R = rotation;
					X = (uint32_t)x;
				}

				Uint32 linnorm32::raw32() {
					//uint32_t z = (stateA = stateA * UINT32_C(0x41C64E6D) + UINT32_C(1));// + (stateC >> 1 ^ (-(stateC & 1) & UINT32_C(0xA3000000)));//stream);
					//uint32_t y = (stateB = (stateB * UINT32_C(163490618)) % UINT32_C(0x7FFFFFFF)) ^ z;
					//z = (z ^ z >> 16u) * UINT32_C(0xAEF1F2D9); //0xAEF17502108EF2D9
					//y = (y ^ y >> 14u) * UINT32_C(0x94D041EB) + UINT32_C(0x6A5D39E1);
					//return (z ^ z >> 15u) + (y ^ y >> 13u);

					//stateB += UINT64_C(0x9E3779B97F4A7C15);
					//stateA += stateA >> 21;

					//stateA = rotate32(stateA, 13) * UINT32_C(0x89A7);
					//stateB = rotate32(stateB, 17) * UINT32_C(0xBCFD);
					//return stateA ^ stateB;
					//return (stateA = rotate32(stateA, 17) * UINT32_C(0xBCFD)) + (stateB = (stateB ^ UINT32_C(0x9E3779BD)) * UINT32_C(0xA952B));
					// ((stateB += UINT32_C(0x9E3779B9)) | UINT32_C(1));
					//stateC = UINT32_C(0xC3E157C1) - rotate32(stateC, 19);
					//stateB += stateB >> 19;
					//stateA += stateA << 8;
					//stateB += stateB << 11;
					//public class SimpleRandom extends Random {public int a=1,b=1; public int next(int bits){return ((a = (a << 13 | a >>> 19) * UINT32_C(0x89A7)) ^ (b = (b << 17 | b >>> 15) * UINT32_C(0xBCFD))) >>> -bits;}}

					//z = (z ^ rotate32(z, 11) ^ rotate32(z, 21)) * (z | UINT32_C(0xFFFE0001)) + (z ^ z >> 14);

					//uint32_t z = ((++stateA) ^ UINT32_C(0xD1B54A35)) * UINT32_C(0x102473);

					//// passes as much as a 32-bit generator can
					//uint32_t z = (stateA += UINT32_C(0x62BD5));
					//z = (z ^ z >> 11 ^ z >> 21) * (z | UINT32_C(0xFFE00001));
					//return z ^ z >> 13 ^ z >> 19;

					uint32_t x = stateA + 1U;
					//x = reverse_bits32(x); // this line was commented in and out for different test types
					x = rotate32(x, R);
					x ^= -X; // change X as a parameter to the generator to 0 to disable any bit flips, 1 to flip all, other (32-bit max, assigned to 64-bit) numbers are allowed
//					x ^= x >> 17;
//    				x *= UINT32_C(0xed5ad4bb);
//    				x ^= x >> 11;
//    				x *= UINT32_C(0xac4c1b51);
//    				x ^= x >> 15;
//    				x *= UINT32_C(0x31848bab);
//    				x ^= x >> 14;
//					x += -(x << 6);
//					x ^= 0xD1B54A35U;
//					x *= 0x94D041EBU;
//					x ^= x >> 17;
//					x += -(x << 9);
//					x ^= 0xD1B54A35U;
//					x *= 0x94D041EBU;
//					x ^= x << 4;
//					x ^= 0xD1B54A35U;
//					x *= 0x94D041EBU;
//					x += -(x << 3);
//					x ^= x << 10;
//					x ^= x >> 15;

					x ^= x >> 17;
					x *= UINT32_C(0xed5ad4bb);
					x ^= x >> 11;
					x *= UINT32_C(0xac4c1b51);
					x ^= x >> 15;
					x *= UINT32_C(0x31848bab);
					x ^= x >> 14;
					stateA = x;

					x += -(x << 6);
					x ^= 0xD1B54A35U;
					x *= 0x94D041EBU;
					x ^= x >> 17;
					x += -(x << 9);
					x ^= 0xD1B54A35U;
					x *= 0x94D041EBU;
					x ^= x << 4;
					x ^= 0xD1B54A35U;
					x *= 0x94D041EBU;
					x += -(x << 3);
					x ^= x << 10;
					x ^= x >> 15;
//					x ^= x >> 17;
////					x *= UINT32_C(0xed5ad4bb);
////					x ^= x >> 11;
////					x *= UINT32_C(0xac4c1b51);
////					x ^= x >> 15;
////					x *= UINT32_C(0x31848bab);
////					x ^= x >> 14;
//
//					x ^= x >> 17;
//					x *= UINT32_C(0xed5ad4bb);
//					x ^= x >> 11;
//					x *= UINT32_C(0xac4c1b51);
//					x ^= x >> 15;
//					x *= UINT32_C(0x31848bab);
//					x ^= x >> 14;
//
//					x += -(x << 6);
//					x ^= 0xD1B54A35U;
//					x *= 0x94D041EBU;
//					x ^= x >> 17;
//					x += -(x << 9);
//					x ^= 0xD1B54A35U;
//					x *= 0x94D041EBU;
//					x ^= x << 4;
//					x ^= 0xD1B54A35U;
//					x *= 0x94D041EBU;
//					x += -(x << 3);
//					x ^= x << 10;
//					x ^= x >> 15;

    				return x;
//					s = (s ^ rotate32(s, 21) ^ rotate32(s, 11) ^ 0xD1B54A35U) * 0x1D2473U;
//					s = (s ^ s >> ((s >> 28) + 4) ^ 0x75AE2165U) * 0x190413U;
//					s = (s ^ rotate32(s, 19) ^ rotate32(s, 13) ^ 0xD1B54A35U) * 0x1D2473U;
//					s = (s ^ s >> ((s >> 28) + 4) ^ 0x9E3779BDU) * 0x1A952BU;
////					return s ^ rotate32(s, 13) ^ rotate32(s, 19);
//					return s ^ s >> 23 ^ s >> 17 ^ s >> 11 ^ s >> 5;




					//uint32_t z = ((++stateA) ^ UINT32_C(0xD1B54A35)) * UINT32_C(0x12D453);
					//uint32_t z = (stateA += UINT32_C(0x62BD5));
					//z = (z ^ z >> 11 ^ z >> 21) * (z | UINT32_C(0xFFE00001));
						//* ((z & UINT32_C(0xFFFF8)) ^ UINT32_C(0xCD7B5));
					//return rotate32(z, 21) ^ rotate32(z, 7) * UINT32_C(0x62BD5);
					//return z ^ z >> 13 ^ z >> 19;

					//z = (z ^ z >> 11 ^ z >> 21) * ((z ^ z >> 15) | UINT32_C(0xFFE00001));
					//return z ^ z >> 15 ^ z >> 20 ^ z >> 13;

					//uint32_t z = ((++stateA) ^ UINT32_C(0xD1B54A35)) * UINT32_C(0x90413);
					//z = (z ^ rotate32(z, 13) ^ rotate32(z, 19) ^ UINT32_C(0x9E3779BD)) * UINT32_C(0x80893);
					//z = (z ^ rotate32(z, 11) ^ rotate32(z, 21) ^ UINT32_C(0x75AE2165)) * UINT32_C(0x90413);
					//return (z ^ z >> 5 ^ z >> 11 ^ z >> 19 ^ z >> 29);
					// *UINT32_C(0x81905);
					//z = (z ^ rotate32(z, 9) ^ rotate32(z, 22) ^ UINT32_C(0x1B54A32D)) * UINT32_C(0xA2203);

					// Java
					//return (z = ((z = (z ^ 0xD1B54A35) * 0x102473) ^ (z << 11 | z >>> 21) ^ (z << 21 | z >>> 11)) * ((z ^ z >>> 15) | 0xFFE00001) + z) ^ z >>> 14;
					//uint32_t z = (stateA += UINT32_C(0x62BD5));
					//z = (z ^ z >> 13) * ((z & UINT32_C(0xFFFF8)) ^ UINT32_C(0xCD7B5));
					//return (((z << 21) | (z >> 11)) ^ (((z << 7) | (z >> 25)) * UINT32_C(0x62BD5)));

					//uint64_t z = ++state;
					//z = ((z << ((z & 31) + 5)) ^ z ^ UINT64_C(0xDB4F0B9175AE2165)) * UINT64_C(0xD1B54A32D192ED03);
					//z = (z ^ (z >> ((z >> 60) + 16))) * UINT64_C(0x369DEA0F31A53F85);
					//return z ^ z >> 27;
					// * UINT32_C(0x1268C3) ^ UINT32_C(0x369DEA0D)
					//z = (z ^ rotate32(z, 11) ^ rotate32(z, 21) ^ UINT32_C(0x369DEA0D)) * UINT32_C(0x102473);

					//z = (z ^ (z >> ((z >> 28) + 4)) ^ UINT32_C(0x369DEA0D)) * UINT32_C(0x102473);
					// * UINT32_C(0x188A23) ^ UINT32_C(0xD1B54A35)
					// * UINT32_C(0x1400CB)
					//z = ((z ^ (z << (z & 15) + 4) ^ UINT32_C(0xDB4F0B95))) * UINT32_C(0x102473);
					//z = ((z ^ (z >> ((z >> 28) + 4)) ^ UINT32_C(0x9E3779BD))) * (z & 0x1FFFF8 | 0xA523u);
					//return z ^ rotate32(z, 11) ^ rotate32(z, 21);

				}
				std::string linnorm32::get_name() const {
					std::ostringstream str;
					str << "linnorm32(" << R << "," << X << ")";
					return str.str();
				}

				void linnorm32::walk_state(StateWalkingObject *walker) {
					walker->handle(stateA);
					//walker->handle(stateB);
					//walker->handle(stateC);
					//printf("Seed is 0x%X, 0x%X, 0x%X\r\n", stateA, stateB, stateC);// stream);
					printf("Seed is 0x%08X\r\n", stateA);// stream);
				}

				Uint16 linnormBounded::raw16() {
					Uint64 z = (state = state * UINT64_C(0x41C64E6D) + UINT64_C(1));
					z = (z ^ z >> 27u) * UINT64_C(0xAEF17502108EF2D9);
					return z ^ (z >> 25u);
					//Uint64 rlow = z & UINT64_C(0xFFFFFFFF);
					//return z;
					//Uint64 low = 0, top = UINT64_C(0x10000);
					//return (rlow >> 16) + (z << 16);

				}
				std::string linnormBounded::get_name() const { return "linnormBounded"; }
				void linnormBounded::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
					printf("Seed is 0x%016llX\r\n", state);
				}

				Uint64 mingler::raw64() {
					//const uint64_t y = -(stateB & 1);
					//const uint64_t y = (stateB * 0x09BEAB3A);

					// excellent. period is 0xFFFFFFFFFFFFFFFF0000000000000000; a little less than xoroshiro but without BRank issues. Not as fast as xoroshiro128++
					//const uint64_t y = -(stateB & 1);
					//const uint64_t z = (stateA = stateA * UINT64_C(0x41C64E6D) + UINT64_C(1) + (stateB = stateB >> 1u ^ (y & UINT64_C(0xD800000000000000))));
					//return (z ^ z >> 28) + (y ^ UINT64_C(0x9E3779B97F4A7C15));

					// high quality (usually?), not terribly fast
					//stateB ^= stateB >> 21;
					//uint64_t z = (stateA = (stateA ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0x41C64E6B) + (stateB ^= stateB << 9));
					//return (z ^ z >> 28) + (stateB ^= stateB >> 29);

					// where we left off
					//uint64_t y = (stateB += (stateB == UINT64_C(0x61C8864680B583EB) ? UINT64_C(0x3C6EF372FE94F82A) : UINT64_C(0x9E3779B97F4A7C15)));
					//y ^= y >> 28;
					//uint64_t z = (stateA = (stateA ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0x41C64E6B) + y);
					//return (z ^ z >> 28) + y;

					//const uint64_t y = stateB * 0x41C64E6B;
					//return (stateB = rotate64(y, 23)) ^ (stateA += UINT64_C(0x9E3779B97F4A7C15));

					//stateA *= stateA;
					//stateA += (stateB += UINT64_C(0x9E3779B97F4A7C15));
					//stateA ^= stateA >> 31;
					////stateA = rotate64(stateA, 32);
					//return stateA ^ stateB;// ^ (stateB ^ stateB >> 32);

					//return (stateA = (stateA << 5) + rotate32(stateA, 1)) * ((stateB += 0x41C64E6D) | 1);
					//return (stateA = (stateA << 5) + rotate64(stateA, 1)) * ((stateB += UINT64_C(0x9E3779B97F4A7C15)) | 1);

//					const uint64_t b = stateB += 0x9E3779B97F4A7C16ULL;
//					const uint64_t a = stateA += ((b > 0x3C6EF372FE94F82CULL) ? 0x6C8E9CF570932BD5ULL : 0ULL); //will be 0x3C6EF372FE94F82CL if signed
//					const uint64_t z = (a ^ a >> 28) * b;
//					return z ^ z >> 28u;

//const uint64_t s = (stateA += 0x6C8E9CF570932BD5ULL);
//if (s == 0ULL) stateB -= 0x9E3779B97F4A7C15ULL;
//const uint64_t z = (s ^ s >> 27) * ((stateB += 0x9E3779B97F4A7C15ULL) | 1ULL);
//return z ^ z >> 25;

//const uint64_t s = (stateA += 0xC6BC279692B5C323ULL);
//const uint64_t z = (s ^ s >> 27) * ((stateB += 0x9E3779B97F4A7C15ULL) | 1ULL);
//if (s == 0ULL) stateB -= 0x9E3779B97F4A7C15ULL;
//return z ^ z >> 27;
					// where we left off
					//stateA *= UINT64_C(0x41C64E6B);
					//stateA = rotate64(stateA, 28);
					//stateB = UINT64_C(0xC6BC279692B5CC8B) - rotate64(stateB, 35);
					//return stateA ^ stateB;
					//end where we left off

					//stateA = stateA * UINT64_C(0x41C64E6D) + UINT64_C(1);
					//stateA += stateA >> 48;
					//stateB += UINT64_C(0x9E3779B97F4A7C15);
					//stateA = UINT64_C(0x9E3779B97F4A7C15) - rotate64(stateA, 42);
					//stateA = rotate64(stateA, 43) - stateA;
					//uint64_t s0 = s[0] * UINT64_C(0x41C64E6B), s1 = s[1];
					//return (s[0] = rotl(s0, 28)) ^ (s[1] = UINT64_C(0xC6BC279692B5CC8B) - rotate64(s1, 35));
					//stateB = UINT64_C(0xC13FA9A902A6328F) - rotate64(stateB, 19);
					//stateB = 3286325185u - rotate32(stateB, 19);
					//stateA += (stateA << 16);
					//stateB += stateB >> 22;
					//stateB += stateB << 5;
					//return z + (stateB ^= stateB >> 29);
					//return (stateA = stateA * UINT64_C(0x41C64E6D) + UINT64_C(1) + (stateB = stateB >> 1u ^ (y & UINT64_C(0xD800000000000000))))
					//	* UINT64_C(0x369DEA0F31A53F85) + (y ^ UINT64_C(0x9E3779B97F4A7C15));
					//return (z ^ z >> 27) + (y ^ UINT64_C(0x9E3779B97F4A7C15));
					// WORKING WELL
					// ((y ^ UINT64_C(0x369DEA0F31A53F85)) | UINT64_C(1))
					//uint64_t y = stateB ^ stateB >> 31;
					//const uint64_t z = (stateA = (stateA * UINT64_C(0x41C64E6D)) + UINT64_C(1)) + (y ^= y << 25);
					//return (z ^ z >> 27u) + (stateB = y ^ y >> 37);

					//YEP
					//uint64_t z = (stateA = (stateA ^ UINT64_C(0x6C8E9CF570932BD5)) * UINT64_C(0x41C64E6B)) + (stateB = stateB >> 1u ^ (y & UINT64_C(0xD800000000000000)));
					//const uint64_t z = (stateA = (stateA ^ UINT64_C(0x6C8E9CF570932BD5)) * UINT64_C(0x41C64E6B)) + (y ^= y << 25);
					//???
					//const uint64_t z = (stateA = (stateA * UINT64_C(0x41C64E6D)) + UINT64_C(1)) + (y ^= y << 25);
					//NOPE
					//uint64_t z = (stateA = (stateA ^ UINT64_C(0x6C8E9CF570932BD5)) * UINT64_C(0x5DA942042E4DD58B)) + (stateB = (y & UINT64_C(0x7FFFFFFF))+(y>>31));

					//YEP?
					//uint64_t z = (stateA = (stateA * UINT64_C(0x41C64E6D)) + UINT64_C(1)) + (stateB = stateB >> 1u ^ (y & UINT64_C(0xD800000000000000)));
					//NOPE
					//uint64_t z = (stateA += UINT64_C(0x632BE59BD9B4E019)) + (stateB = stateB >> 1u ^ (y & UINT64_C(0xD800000000000000)));

					//Uint64 z = (stateA += UINT64_C(0x632BE59BD9B4E019)) ^ y;// (stateA += UINT64_C(0x9E3779B97F4A7C15));// (stateA = stateA * UINT64_C(0x41C64E6D) + UINT64_C(1));
					//z *= y | UINT64_C(1);
					//return (z ^ z >> 33) + y;
					//???
					//return (z ^ z >> 27u) + (stateB = y ^ y >> 37);

					// + UINT64_C(0x6A5D39EAE12657AA)
					//const uint64_t
					//	y = (stateB = stateB >> 1 ^ (-(stateB & 1) & UINT64_C(0xD800000000000000))),
					//	x = ((stateA += UINT64_C(0x9E3779B97F4A7C15)) ^ y) * UINT64_C(0xAEF17502108EF2D9);
					//z = ((x ^ x >> 27) + (y ^ y >> 30)) * UINT64_C(0xAEF17502108EF2D9);
					//return (x ^ x >> 25) + y;
					//(stateA = stateA * UINT64_C(0x41C64E6D))
					//z = (z ^ z >> 27u) * UINT64_C(0xC6BC279692B5CC83);
					//z = (z ^ z >> 27u) * UINT64_C(0xAEF17502108EF2D9);
					//return z ^ (z >> 25u);
					//					stateB ^= stateB << 13;
					//					stateB ^= stateB >> 7;
					//					stateB ^= stateB << 17;
					//					uint64_t z = (stateA += UINT64_C(0x9E3779B97F4A7C15));// stateA * UINT64_C(0x41C64E6D) + UINT64_C(1);
					//					uint64_t z = ((stateA += UINT64_C(0x632BE59BD9B4E019)) ^ UINT64_C(0x9E3779B97F4A7C15)) * UINT64_C(0xC6BC279692B5CC83);
					//					z = ((z ^ z >> 27) + _pext_u64(z + 0x9E3779B97F4A7C15, rotate64(z, 41))) * UINT64_C(0xAEF17502108EF2D9);
					//					return z ^ z >> 25;


					//// TroutRNG, passes 32TB with no anomalies, period of (2 to the 128) minus (2 to the 64)
					//uint64_t s = (stateA += (stateB = (stateB >> 1 ^ (-(stateB & 1ULL) & 0xD800000000000000ULL))) + 0x9E3779B97F4A7C15ULL);
					//s = (s ^ rotate64(s, 41) ^ rotate64(s, 17)) * 0x369DEA0F31A53F85ULL;
					//return s ^ s >> 28;

					//// BellRNG, passes 32TB with one "unusual" anomaly. Same period as Trout.
					//uint64_t s = (stateA += (stateB = (stateB >> 1 ^ (-(stateB & 1ULL) & 0xD800000000000000ULL))) + 0x9E3779B97F4A7C15ULL);
					//s = (s ^ s >> 30) * 0x369DEA0F31A53F85ULL;
					//return s ^ s >> 28;

					//// OrbitRNG from Sarong.
					//// Passes 32TB with no anomalies on two seeds, 3 "unusual" on another, one "mildly suspicious" on another.
					//// Tested heavily due to similarity with ThrustAltRNG; no TMFn issues here, or any others.
					//// Should have a period of exactly 2 to the 128, and allow all states; one-dimensionally equidistributed.
					//const uint64_t s = (stateA += 0xC6BC279692B5C323UL);
					//const uint64_t z = (s ^ s >> 27) * ((stateB += 0x9E3779B97F4A7C15UL) | 1UL);
					//if (!s) stateB -= 0x9E3779B97F4A7C15UL;
					//return z ^ z >> 27;

					//s = (s ^ rotate64(s, 41) ^ rotate64(s, 17)) * 0x369DEA0F31A53F85ULL;
					
					//s = (s ^ s >> 25 ^ s >> 37) * 0xDB4F0B9175AE2165UL;

					//fails TMFn at 2TB
					//const uint64_t s = (stateA += 0x6C8E9CF570932BD5UL);
        			//const uint64_t z = (s ^ s >> 25) * (stateB += 0x9E3779B97F4A7C16UL);
					//return z ^ z >> 22;

					//this one got "very suspicious" at 32TB; switching to shifts of 29 and 27 helps a little
					//const uint64_t s = (stateA += 0xC6BC279692B5C323UL);
        			//const uint64_t z = (s ^ s >> 27) * (stateB += 0x9E3779B97F4A7C16UL);
					//return z ^ z >> 27;

					////works; passes 32TB with one minor anomaly at 16GB: Low8/64]DC6-9x1Bytes-1 , unusual
					//const uint64_t s = (stateA += 0xC6BC279692B5C323UL);
        			//const uint64_t z = (s ^ s >> 31) * (stateB += 0x9E3779B97F4A7C16UL);
					//return z ^ z >> 26;

					//// passes 32TB, one anomaly at 8TB, extremely minor
					//// Low8/64]FPF-14+6/16:cross        R=  -2.5  p =1-1.6e-4   unusual
					//const uint64_t s = (stateA += 0xC6BC279692B5C323UL);
        			//const uint64_t z = (s ^ s >> 31) * (stateB += 0x9E3779B97F4A7C16UL);
					//return z ^ z >> 26 ^ z >> 6;

					//// New OrbitRNG
					//// passes 32TB with one anomaly at 32GB, "[Low8/64]BCFN(2+4,13-1,T)" rated as "unusual"
					//// also unusually fast in Quick-Bench benchmarks... Faster than SplitMix64, Xoroshiro128+, Romu...
					//uint64_t s = (stateA += 0xC6BC279692B5C323u);
  					//uint64_t t = ((s == 0u) ? stateB : (stateB += 0x9E3779B97F4A7C15u));
  					//uint64_t z = (s ^ s >> 31) * ((t ^ t >> 22) | 1u);
					//return z ^ z >> 26;

					//// Passes 32TB, one "unusual" at 32TB, BCFN(2+1,13-0,T)
					//uint64_t s = (stateA += 0xC6BC279692B5C323u);
  					//uint64_t t = ((s == 0u) ? stateB : (stateB += 0x9E3779B97F4A7C15u));
  					//uint64_t z = (s ^ s >> 29) * (t | 1u);
					//return z ^ z >> 28;

					//uint64_t s = (stateA += 0x9E3779B97F4A7C15u);

					//// passes 32TB but just barely, getting "unusual" at 32TB
					//// exact anomaly:
					//// length= 32 terabytes (2^45 bytes), time= 113425 seconds
  					////   Test Name                         Raw       Processed     Evaluation
  					////   BCFN(2+1,13-0,T)                  R= +10.2  p =  5.5e-5   unusual          
  					////   ...and 1133 test result(s) without anomalies
					//uint64_t s = (stateA += 0xC6BC279692B5C323u);
  					//uint64_t t = ((s < 0x18DAF2Du) ? stateB : (stateB += 0x9479D2858AF899E6u));
  					//uint64_t z = (s ^ s >> 31) * t;
					//return z ^ z >> 26;

					////A new variant on Orbit that passes 32TB, and not just by a hair!
					////
					////length= 16 terabytes (2^44 bytes), time= 55924 seconds
					////  Test Name                         Raw       Processed     Evaluation
					////  BCFN(2+1,13-0,T)                  R=  +9.0  p =  2.2e-4   unusual          
					////  ...and 1106 test result(s) without anomalies
					////...
					////length= 32 terabytes (2^45 bytes), time= 112975 seconds
					////  no anomalies in 1134 test result(s)
					//const uint64_t s = (stateA += 0xC6BC279692B5C323u);
  					//const uint64_t z = ((s < 0x6F17146Du) ? stateB : (stateB += 0x9479D2858AF899E6u)) * (s ^ s >> 31);
					//return z ^ z >> 25;

					////has no anomalies before 32TB, but...
					////has trouble at 32TB, getting a "suspicious" BCFN(2+1,13-0,T) as well as two "unusual" TMFn
					//const uint64_t s = (stateA += 0xC6BC279692B5C323u);
                    //const uint64_t t = ((s == 0u) ? stateB : (stateB += 0x9479D2858AF899E6u)) * (s ^ s >> 31);
                	//return t ^ t >> 25;

					//// meant to mimic the behavior of signed comparisons on long values in Java
					//// passes 32TB with no anomalies (!)
					//// other than an offset in the comparison, this is identical to the generator above that also compares against 0x6F17146Du
					//const uint64_t s = (stateA += 0xC6BC279692B5C323u);
  					//const uint64_t z = ((s + 0x8000000000000000u < 0x6F17146Du) ? stateB : (stateB += 0x9479D2858AF899E6u)) * (s ^ s >> 31);
					//return z ^ z >> 25;

					//// passes 32TB with one "unusual" anomaly at 16TB. Can sometimes get a "mildly suspicious" at 32TB.
					//uint64_t s = (stateA += 0xC6BC279692B5C323u);
					//s ^= s >> 31;
  					//s *= ((s < 0x6F17146Du) ? stateB : (stateB += 0x9479D2858AF899E6u));
					//return s ^ s >> 25;

					//0x9E3779B97F4A7C15u weight 38 // 0xB1E131D6149D9795u weight 32 // 0xC15E014828F41863u weight 24 // 0xD1342543DE82EF95u weight 32

					//// has frustrating late issues at 32TB and later, even though it passes 32TB with one "unusual," it gets a "very suspicious" at 64TB.
					//// BCFN(2+1,13-0,T)
					//uint64_t s = (stateA += 0xB1E131D6149D9795u);
					//s ^= s >> 31;
  					//s *= ((s < 0xB0BAFE77u) ? stateB : (stateB += 0x9479D2858AF899E6u));
					//return s ^ s >> 25;

					////GoatRNG (I just like goats...)
					////Passes 32TB with one "unusual" anomaly, [Low1/16]DC6-9x1Bytes-1 , at 512GB. Gets "VERY SUSPICIOUS" BCFN at 64 TB...
					////Also allows all 2 to the 128 states possible and should mitigate the seed correlation in similar versions.
					////This is also very fast in Quick-Bench benchmarks using GCC; faster than SplitMix64 and on-par with Gear.
					//uint64_t s = (stateA += 0xD1342543DE82EF95u);
					//const uint64_t t = ((s ^= s >> 31) < 0xB1E131D6149D9795u ? (stateB += 0xC6BC279692B5C323u) : stateB);
					//s *= ((t ^ t >> 29 ^ t << 11) | 1u);
					//return s ^ s >> 25;

					////GhoulRNG (running out of G-related words...)
					////Passes 64TB with no anomalies. Ran out of memory (had 16GB available) trying to test 128TB.
					////This is a very slight tweak on GoatRNG and should have similar speed, or be a bit faster due to its strict conditional.
					////Also allows all 2 to the 128 states possible and should mitigate the seed correlation in similar versions.
					////Where GoatRNG updates stateB roughly 11/16 times, this only updates it less than 1/4 times.
					////Instead of using one of the Evensen-discovered invertible shift pairs, this just uses a paired right shift on stateA.
					//uint64_t s = (stateA += 0xD1342543DE82EF95u);
					//const uint64_t t = ((s ^= s >> 31 ^ s >> 23 ) >= 0xC6BC279692B5C323u ? (stateB += 0xB1E131D6149D9795u) : stateB);
					//s *= ((t ^ t << 9) | 1u);
					//return s ^ s >> 25;


//					uint64_t s = (stateA += 0xD1342543DE82EF95u);
//					s ^= s >> 31 ^ s >> 23;
//  					const uint64_t t = (stateB += 0xB1E131D6149D9796u);
////  					const uint64_t t = (((s ^= s >> 31) == 0u) ? stateB : (stateB += 0xB1E131D6149D9795u));
//					s *= (t ^ t << 13);
//					return s ^ s >> 25;


//					const uint64_t s = (stateA += 0xD1342543DE82EF95u);
//    				const uint64_t z = (s ^ s >> 31) * (stateB += 0xB1E131D6149D9796u);
//      			return z ^ z >> 25;

////  BCFN(2+1,13-0,T)                  R= +17.9  p =  4.3e-9   very suspicious  
//					const uint64_t s = (stateA += 0xCB9C59B3F9F87D4Du);
//    				const uint64_t z = (s ^ s >> 31) * (stateB += 0xF468D97FAB12104Au);
//      			return z ^ z >> 25;

//// No anomalies through 32TB, but at 64TB it gets 4 TMFn anomalies plus:
////  BCFN(2+1,13-0,T)                  R= +13.8  p =  6.6e-7   suspicious
//					const uint64_t s = (stateA += 0xCB9C59B3F9F87D4DL);
//    				const uint64_t z = (s ^ s >> 31) * (stateB += 0x9738B367F3F0FA9Au);
//					return z ^ z >> 26;

//// This one's very strong, with no anomalies in the first 32TB and just one "unusual" at 64TB:

////rng=mingler, seed=0x0
////length= 64 terabytes (2^46 bytes), time= 234778 seconds
////  Test Name                         Raw       Processed     Evaluation
////  BCFN(2+1,13-0,T)                  R= +10.2  p =  5.7e-5   unusual          
////  ...and 1158 test result(s) without anomalies

//// But, on a different seed, it's a different story.
////rng=mingler, seed=0x1
////length= 64 terabytes (2^46 bytes), time= 235817 seconds
////  Test Name                         Raw       Processed     Evaluation
////  BCFN(2+1,13-0,T)                  R= +17.2  p =  9.5e-9   very suspicious  
////  [Low1/64]TMFn(2+12):wl            R= +22.1  p~=   8e-7    unusual          
////  ...and 1157 test result(s) without anomalies

//					const uint64_t s = (stateA += 0xCB9C59B3F9F87D4Du);
//    				const uint64_t z = (s ^ s >> 31) * (stateB += 0x3463A64C060782B2u);
//					return z ^ z >> 27;


//// Partly tried this, seemed fine to at least 256GB

//	  				uint64_t a = (stateA += 0xCB9C59B3F9F87D4Du);
//					if(a) a = (a ^ a >> 31) * (stateB += 0x3463A64C060782B2u);
//					return a ^ a >> 27;

                    // promising. no anomalies in the first 32GB, but that isn't much.
					//const uint64_t a = ((stateA = stateA * 0xD1342543DE82EF95u + 0xCC62FCEB9202FAADu) >> 31) * (stateB += 0x3463A64C060782B2u);
					//return a ^ a >> 27;

					//const uint64_t a = ((stateA += 0xCC62FCEB9202FAADu)) * (stateB += 0x3463A64C060782B2u);
					//return a ^ a >> 27;


//					const uint64_t s = (stateA += 0xCB9C59B3F9F87D4Du);
//					const uint64_t t = (stateB += 0x9E3779B97F4A7C16u) * (s ^ s >> 31);
////					const uint64_t t = (s ? (stateB += 0x91E10DA5C79E7B1Eu) * (s ^ s >> 31) : 0);
//					return t ^ t >> 27;

//					uint64_t s = (stateA += 0xEA1FBB3AA44F0B9Du);
//					const uint64_t t = ((s ^= s >> 31) ? (stateB += 0x841C6FFBC5B6AF25u) : 0u);
//					s *= ((t ^ t << 9) | 1u);
//					return s ^ s >> 27;

//					const uint64_t s = (stateA += 0xCB9C59B3F9F87D4Du);
//    				uint64_t z = (s ^ s >> 31) * (stateB += 0x3463A64C060782B2u);
//					z ^= z >> 23;
//					return z ^ z >> 27;


//// Current TangleRNG; passes 64TB with no anomalies on seed=1
//// Also quite fast on the JVM; period is 2 to the 64, with 2 to the 63 possible streams.
//// This particular set of large constants and shift amounts does well; it probably isn't unique.
//// The first large constant must be odd to guarantee full period; other than that, a lot of constants were tested.
//// The first is a probable prime according to Java's BigInteger, which uses Lucas-Lehmer. I don't think primality is required.
//// The second large constant is 2 to the 64 divided by the golden ratio, rounded to the nearest integer.
//// The second constant needs to be 2 times any odd integer; equivalently, it must be equal to 2, mod 4.
//// This requirement, along with stateB always starting as an odd number, keeps stateB odd, and ensures stateB will have
////   been all possible odd numbers twice over the period. As a multiplier, stateB must be odd.
//// The double-xorshift is a bijection, as is the earlier single-xorshift, but I don't know how easy it is to recover the state from 2 or more outputs.
//// This generator is not equidistributed over its period, but if you sequentially changed streams each time period was exhausted (every several years),
////   then you'd have a period of 2 to the 127, one-dimensionally equidistributed. You'd also be dead for a few quintillion years.
//					const uint64_t s = (stateA += 0xC6BC279692B5C323UL);
//        			const uint64_t z = (s ^ s >> 31) * (stateB += 0x9E3779B97F4A7C16UL);
//					return z ^ z >> 26 ^ z >> 6;

//// passes 32TB, no anomalies, but fails with BCFN at 64TB
//					const uint64_t s = (stateA += 0xC6BC279692B5C323UL);
//        			const uint64_t z = (s ^ s >> 31) * (stateB += 0x74ED9428DE96EA4AUL);
//					return z ^ z >> 27;

//// mildly suspicious at 32TB, some other anomalies on the way.
//					const uint64_t s = (stateA += 0xC6BC279692B5C323UL);
//        			const uint64_t z = (s ^ s >> 31) * (stateB += 0x61C8864680B583EAUL);
//					return z ^ z >> 27;

//					const uint64_t s = (stateA += 0xC6BC279692B5C323UL);
//        			const uint64_t z = (s ^ s >> 31) * (stateB += 0x9A1AD7561A3AC55AUL);
//					return z ^ z >> 25;

//0x9A1AD7561A3AC55AUL 32
//0x9794A2A79714EA96UL 32

//0x95B534A1ACCD52DAUL 32

//0xde01abbf8f022f55u 34 //problems
//0x61C8864680B583E9u 26
//0x9E3779B97F4A7C16u 38
//0xC13FA9A902A6328Fu 30
//0x91E10DA5C79E7B1Du 34

//0xCC62FCEB9202CAB9u 32
//0x05CB93402A76F7DAu 32

//0xD1342543DE82EF23u 31
//0x74ED9428DE96EA4Au 33
//0xCB9C59B3F9F87D4DUL  0x3463A64C060782B2L

//0xD1342543DE82EF95u

////tmfn at 32, failures at 64
//0x98C5F9D72405F55Au

//0xFA346CBFD5890825u
//0xCC62FCEB9202FAADu
//0x9738B367F3F0FA9Au
//0xCB9C59B3F9F87D4DL 0xF468D97FAB12104AL

//0x9E3779B97F4A7C15u

					//const uint64_t s = (stateA += 0xEB44ACCAB455D165ULL);
        			//const uint64_t z = (s ^ s >> 31 ^ s >> 9) * 0xE7037ED1A0B428DBULL;
					//return z ^ z >> 27;
					// mul:C6BC279692B5C323,xorr:31,xorr:17,mul,xorr:26

					// passes 32TB with no anomalies, but gets VERY SUSPICIOUS Gap test at 64 TB
	  				//uint64_t a = (stateA += 0x9E3779B97F4A7C15UL); //0xC6BC279692B5C323UL //0xCB9C59B3F9F87D4DUL
					//a = (a ^ rotate64(a, 38) ^ rotate64(a, 19)) * 0xCC62FCEB9202FAADUL;
					//return a ^ a >> 43 ^ a >> 27;

					//0xC6BC279692B5C323UL //0xCB9C59B3F9F87D4DUL

					// Passes 64TB with just one early (and minor) anomaly! On seed 0:
					// [Low4/64]DC6-9x1Bytes-1           R=  -5.2  p =1-1.4e-3   unusual
					// On seed 1, no anomalies at all in 64TB of testing.
	  				//uint64_t a = (stateA += 0x9E3779B97F4A7C15UL);
					//a = (a ^ rotate64(a, 47) ^ rotate64(a, 23)) * 0xCC62FCEB9202FAADUL;
					//return a ^ a >> 30 ^ a >> 26;

					//fails at 4TB, even with extra multiply, when using 1 as the increment.
					  // 0xF1DE83E19937733DUL is the modular multiplicative inverse of 0x9E3779B97F4A7C15UL.
					  // 0x82B3C6EF372FE94FUL is that rotated right by 11; when used as the increment, it passes 64TB with two anomalies, both DC6-9.
					//uint64_t a = ++state;
					//a = rotate64(a, 11) * 0x9E3779B97F4A7C15UL;
					//a = (a ^ rotate64(a, 47) ^ rotate64(a, 23)) * 0xCC62FCEB9202FAADUL;
					//return a ^ a >> 30 ^ a >> 26;

					// Passes 64TB with one anomaly (4/16 BCFN, rated "unusual").
					// Should be a good determine() method, not sure.
	  				//uint64_t a = (++stateA ^ 0x9E3779B97F4A7C15UL) * 0xC6BC279692B5C323UL;
					//a = (a ^ rotate64(a, 47) ^ rotate64(a, 23)) * 0xCC62FCEB9202FAADUL;
					//return a ^ a >> 30 ^ a >> 26;
	  				
					//return stateA = (stateA ^ rotate64(stateA, 47) ^ rotate64(stateA, 23) ^ 0x9E3779B97F4A7C15UL) * 0xC6BC279692B5C323UL;
					
					//// passes 64TB with no anomalies!
					//// ShinyRNG
					//uint64_t a = (stateA = stateA * 0xCC62FCEB9202FAADUL + 0x9E3779B97F4A7C15UL);
					//a = (a ^ a >> 31) * 0xC6BC279692B5C323UL;
					//return a ^ a >> 25;

					//const Uint64 b = (stateB = (stateB >> 1ULL ^ (0ULL - (stateB & 1ULL) & 0xD800000000000000ULL)));
					//const Uint64 b = (stateB += 0xCC62FCEB9202FAADUL);

					////GingerRNG; passes 64TB with no anomalies.
					////Period is very good: (2 to the 64) * ((2 to the 64) - 1)
					// 1-dimensionally equidistributed, and without the caveat that xoroshiro has (xoroshiro and xoshiro generate 0 less often).
//				    Uint64 s = (stateA += 0xC6BC279692B5C323UL ^ (stateB = (stateB >> 1ULL ^ (0ULL - (stateB & 1ULL) & 0xD800000000000000ULL))));
//					s = (s ^ s >> 31) * 0xCC62FCEB9202FAADUL;
//					return s ^ s >> 28;

					//// nothing special. slower than above, gets 3 unusual anomalies over 64TB
					//stateB ^= stateB << 7;
					//const Uint64 s = (stateA += 0xC6BC279692B5C323UL ^ (stateB ^= stateB >> 17));
					//return (s ^ s >> 31) * ((stateB ^= stateB >> 13) | 1UL);

					////Problems with this one.
					////It has this anomaly twice; "unusual" at 32TB and "very suspicious" at 64 TB:
					////[Low4/64]FPF-14+6/16:(2,14-0)
					////It also has two other "mildly suspicious" anomalies at 64TB, one FPF, one DC6-9.
					//// I am also not yet certain that it has a good period; it's at minimum 0xFFFFFFFFFFFFFFFF0000000000000000, but should ideally be the square of that.
//				    const uint64_t s = (stateA += 0xCC62FCEB9202FAADUL ^ (stateB = (stateB >> 1UL ^ (0UL - (stateB & 1UL) & 0xD800000000000000UL))));
//				    const uint64_t i = s * 0xC6BC279692B5C323UL ^ ((s < 0xD1342543DE82EF95UL) ? incA : (incA += 0x9E3779B97F4A7C15UL ^ (incB = (incB >> 1UL ^ (0UL - (incB & 1UL) & 0x8000000000000212UL)))));
//					return i ^ i >> 30;

					////passes 64TB with one "unusual" anomaly at 64TB:
					//// BCFN(2+0,13-0,T)
//					uint64_t s = (stateA += 0xCC62FCEB9202FAADUL);
//					s = (s ^ s >> 31 ^ (stateB = s < 0xD1342543DE82EF95UL ? stateB : (stateB >> 1UL ^ (0UL - (stateB & 1UL) & 0xD800000000000000UL)))) * 0xC6BC279692B5C323UL;
//					return s ^ s >> 28;
					
					////again, passes 64TB with one anomaly at 64TB:
					//// [Low8/32]DC6-9x1Bytes-1           R=  -6.2  p =1-1.2e-3   unusual
//					uint64_t s = (stateA += 0xCC62FCEB9202FAADUL);
//					s = (s ^ s >> 31 ^ (stateB = s < 0xBA987654321FEDCBUL ? stateB : (stateB >> 1UL ^ (0UL - (stateB & 1UL) & 0xD800000000000000UL)))) * 0xC6BC279692B5C323UL;
//					return s ^ s >> 28;
					
					////again, passes 64TB with one anomaly, but this time "unusual" at just 2TB, and nothing before or after:
					//// [Low4/32]FPF-14+6/16:all          R=  -5.3  p =1-1.0e-4   unusual
					////I think this is an improvement!
					//// Speed is good, on-par with SplitMix64: https://quick-bench.com/q/8bBSdbMZh5THGlKzXSq3y350Iiw
					////1D equidistributed, period is 0xFFFFFFFFFFFFFFFF0000000000000000 .
//					uint64_t s = (stateA += 0xCC62FCEB9202FAADUL);
//					s = (s ^ s >> 31 ^ (stateB = s < 0xD1342543DE82EF95UL ? stateB : (stateB >> 1UL ^ (0UL - (stateB & 1UL) & 0xD800000000000000UL)))) * 0xC6BC279692B5C323UL;
//					return s ^ s >> 28 ^ s >> 5;

					////
//					uint64_t s = (stateA += 0xCC62FCEB9202FAADUL);
//					s = (s ^ s >> 31 ^ (stateB = s < 0xDB4F0B9175AE2165UL ? stateB : (stateB >> 1UL ^ (0UL - (stateB & 1UL) & 0x8000000000001797UL)))) * 0xC6BC279692B5C323UL;
//					return s ^ s >> 28 ^ s >> 6;

//					//0xDB4F0B9175AE2165
//					uint64_t s = (stateA += 0xCC62FCEB9202FAADUL);
//					s = (s ^ s >> 31) * (s < 0xD1342543DE82EF95UL ? incB : (incB *= 0xF1357AEA2E62A9C5UL));
//					return s ^ s >> 28 ^ s >> 5;
					
					//// GrognardRNG
					//// Passes 64TB with no anomalies!
					//// Grognard is an old French military position for a kind of heavy-weight-carrier.
					//// This is certainly able to handle heavy loads of generation.
					//// Period is 2 to the 128, exactly.
					//// This should be able to step backwards, but not skip. It can leap by updating stateB without stateA.
					//// The zero check for stateB's update allows Java to run this a bit faster without a conditional, I think.
					//uint64_t s = (stateA += 0xD1342543DE82EF95u);
					//const uint64_t t = (s ? (stateB += 0xB1E131D6149D9795u) : stateB);
					//s = (s ^ s >> 31 ^ s >> 21) * ((t ^ t << 9) | 1u);
					//return s ^ s >> 28;

					//// GrogRNG
					//// Passes 64TB with no anomalies!
					//// Period is 2 to the 128, exactly.
					//// Like GrognardRNG, but stripped down a little.
					//// This can step backwards, but not skip, though it can leap in 2-to-the-64 distances by updating stateB without stateA.
					//// The zero check for stateB's update allows Java to run this a bit faster without a conditional, I think.
					//uint64_t s = (stateA += 0xD1342543DE82EF95u);
					//s = (s ^ s >> 31 ^ s >> 21) * ((s ? (stateB += 0xC6BC279692B5C323u) : stateB) | 1u);
					//return s ^ s >> 28;

//*= 0xF1357AEA2E62A9C5UL
					//// TremorRNG
					//// Passes 32TB with no anomalies, but then has a significant anomaly at 64TB:
					//// BCFN(2+1,13-0,T)                  R= +14.0  p =  5.1e-7   suspicious
					//// Period is 2 to the 128, exactly.
					//// GrogRNG should be about the same speed or faster, and higher quality.
					//const uint64_t s = (stateA += 0xC6BC279692B5C323u);
					//const uint64_t t = (stateB += (s >= 0xC6BC279692B5C323u) ? 0u : 0x9E3779B97F4A7C15u);
					//const uint64_t z = (s ^ s >> 31) * (t | 1U);
					//return z ^ z >> 27;

//					return s ^ s >> (s >> 59) + 6;
					//s = (s ^ s >> 31) * 0xCC62FCEB9202FAADUL;
//					return s ^ s >> 28;
//				    Uint64 s = (stateA += 0xC6BC279692B5C323UL);
//					s = (s ^ s >> 31) * ((stateB += 0xCC62FCEB9202FAADUL & 0UL - (s < 0x9E3779B97F4A7C15UL)) | 1UL);
//					s = (s ^ s >> 31) * ((stateB = (stateB >> 1UL ^ (0UL - (stateB & 1UL) & 0xD800000000000000UL))) | 1UL);
//					return s ^ s >> 28;

					//// WaltzRNG
					//// Passes 64TB with no anomalies!
					//// Period is 2 to the 128 minus 2 to the 64.
					//// stateA can be any uint64_t. stateB must not be 0.
					//// Should be 1D equidistributed.
					//// Thanks to Pelle Evensen for suggesting the two-xorshift form for stateB, first found by Richard P. Brent
					//// in https://www.jstatsoft.org/article/view/v011i05 .
//					uint64_t s = (stateA += 0xC6BC279692B5C323u), b = stateB;
//					s = (s ^ s >> 31 ^ s >> 5) * (b | 1u);
//					b ^= b >> 7;
//					stateB = b ^ b << 9;
//					return s ^ s >> 28;
					//// Looked promising (one unusual anomaly at 512GB) until 64TB, where it gets this big problem:
					////BCFN(2+1,13-0,T)                  R= +16.1  p =  3.6e-8   very suspicious
					//uint64_t s = (stateA += 0xC6BC279692B5C323u), b = stateB;
					//s = (s ^ s >> 31) * (b | 1u);
					//b ^= b << 9;
					//stateB = b ^ b >> 7;
					//return s ^ s >> 28;

					//// TingleRNG
					//// Passes 64TB with no anomalies. Period is 2 to the 64, allows skipping.
//					uint64_t z = (stateA += 0xC6BC279692B5C323UL) * (stateB += 0x9E3779B97F4A7C16UL);
//					z = (z ^ z >> 31) * 0xD1342543DE82EF95UL; //0xACBD2BDCA2BFF56DUL;
//					return z ^ z >> 26;
					//// TroubleRNG
					////
//        			stateA = stateA * (stateB += 0x9E3779B97F4A7C16UL) + 0xC6BC279692B5C323UL;
//					return (stateA -= rotate64(stateA, 23));

					// I dunno. Passes 16TB without anomalies. State A must not be 0, state B must be odd.
					// uint64_t r = stateA;
					// r ^= r << 9;
					// stateA = (r ^= r >> 7);
					// r *= (stateB += 0x9E3779B97F4A7C16UL);
					// r ^= r >> 31;
					// r *= (stateB);
					// return (r ^ r >> 31);

//  6	                   1  + 2.59071e+00	                   1  + 2.59071e+00	                   0  - 2.29283e-01	                   0  - 2.29283e-01	                   0  - 2.29283e-01	
//  7	                4506  + 1.87025e+00	                4435  + 8.94254e-02	                4365  - 5.69179e-01	                4303  - 2.84773e+00	                4452  + 3.07898e-01	
//  8	              361488  + 1.11662e+00	              360427  - 5.03448e-01	              360956  + 2.92694e-02	              360942  + 2.18381e-02	              360894  + 4.60659e-03	
//  9	            22354248  - 3.25324e-01	            22355380  - 1.09537e-01	            22346149  - 5.21321e+00	            22346225  - 5.14007e+00	            22354908  - 1.85578e-01	
// 11	                4482  + 1.01280e+00	                4440  + 1.40092e-01	                4352  - 9.02663e-01	                4471  + 7.06996e-01	                4328  - 1.71945e+00	
// 12	            85019477  + 3.97348e-03	            85026225  + 6.31830e-01	            85027391  + 8.48856e-01	            85010292  - 8.70688e-01	            85030286  + 1.52598e+00	
// 13	          6948620577  - 5.83832e-01	          6948691394  + 7.30270e-03	          6948616022  - 6.70322e-01	          6948708691  + 8.58235e-02	          6948669486  - 3.14566e-02	
// 14	        430511726017  + 6.78039e-01	        430511648494  + 4.97420e-01	        430511234586  + 5.54297e-03	        430511158897  - 1.67322e-03	        430511665907  + 5.35559e-01	
// 16	              360861  + 1.67363e-04	              359449  - 5.46443e+00	              360559  - 2.39905e-01	              360103  - 1.55976e+00	              360596  - 1.83362e-01	
// 17	          6948721131  + 1.95533e-01	          6948586479  - 1.37626e+00	          6948675996  - 9.85330e-03	          6948602610  - 9.59669e-01	          6948727236  + 2.65667e-01	
// 18	        567923827094  + 4.62508e-01	        567925092849  + 5.56807e+00	        567923184151  - 2.99550e-02	        567922266472  - 1.93430e+00	        567922953162  - 2.30003e-01	
// 19	      35186142187561  + 1.66532e+00	      35186141057870  + 1.21005e+00	      35186133401504  - 3.63701e-02	      35186134393025  - 5.54878e-04	      35186143061073  + 2.06707e+00	
// 21	            22346137  - 5.22480e+00	            22356364  - 1.50934e-02	            22355343  - 1.14778e-01	            22346907  - 4.50685e+00	            22346557  - 4.82662e+00	
// 22	        430511136379  - 5.65868e-03	        430511752868  + 7.47108e-01	        430511662255  + 5.27443e-01	        430511264288  + 1.43327e-02	        430511119519  - 1.01849e-02	
// 23	      35186132811106  - 8.42397e-02	      35186140957397  + 1.17308e+00	      35186142940938  + 2.00925e+00	      35186134284160  - 1.75633e-03	      35186133636723  - 2.28178e-02	
// 24	    2179984564204183  - 2.18239e-02	    2179984555431176  - 1.12646e-01	    2179984563510681  - 2.64331e-02	    2179984572573862  + 9.94154e-04	    2179984563390121  - 2.72793e-02	
// 
//               13.251 15 =>     4.956004e-01           17.646 15 =>     7.769483e-01           11.233 15 =>     3.327710e-01           18.653 15 =>     8.218331e-01           11.944 15 =>     3.888542e-01
// 
// --- Finished -- ReMort trials: 2251799813685248	2^51.0 (out of 2^51) trials -- Wave -- 64 bits: 	ps:   4.956e-01    7.769e-01    3.328e-01    8.218e-01*   3.889e-01   => p <   6.94893e-01   2^54.32 calls, 2^57.32 bytes	2^37.03 bytes/second	used:  14::19:23:47.41

//					//// WaveRandom
//					//// Passes 64TB with no anomalies.
//					//// Not as clearly correlated between similar starting states as Tangle.
//					uint64_t r = (stateA += 0xC6BC279692B5C323UL);
//					r ^= r >> 31;
//					r *= (stateB += 0x9E3779B97F4A7C16UL);
//					r ^= r >> 29;
//					r *= stateB;
//					return (r ^ r >> 30);

//					//// [Low1/64]Gap-16:A and [Low1/64]Gap-16:B fail at 32TB
//					uint64_t a = (stateA = (stateA ^ 0xF1357AEA2E62A9C5UL) * 0xBB67AE8584CAA73BUL);
//					uint64_t b = (stateB = (stateB ^ 0xD1342543DE82EF95UL) * 0xC6BC279692B5C323UL);
//					return (a ^ a >> 41) + (b ^ b >> 19);

					//// PulseRNG
					//// Passes 64TB with one anomaly,
					////   [Low1/8]BCFN(2+2,13-0,T)          R=  +8.9  p =  2.6e-4   unusual 
					//// Period should be (2 to the 128) - (2 to the 64).
					//// Should be 1D equidistributed, producing each result (2 to the 64) - 1 times.
//					uint64_t a = (stateA += 0xD1342543DE82EF95UL);
//					stateB ^= stateB << 7;
//					uint64_t b = stateB ^= stateB >> 9;
//					b ^= b * b | 1UL;
//					b += a ^ a >> 31;
//					//b ^= b * b | 1UL;
//					return (b ^ b >> 30);

//					//uint64_t a = (stateA = -(stateA ^ (stateA * stateA | 5U)));
//					uint64_t a = (stateA = stateA * 0xD1342543DE82EF95UL + 0xDB4F0B9175AE2165UL);
//					a ^= a >> 32;
//					a *= (stateB += 0x9E3779B97F4A7C15UL) | 1UL;
//					a ^= a >> 31;
//					return a;

//uint64_t z = (state += 0xC13FA9A902A6328FUL);
//uint64_t y = (stream += 0x91E10DA5C79E7B1DUL);
//z += y ^ rotate64(y, 5) ^ rotate64(y, 14);
//y += z ^ rotate64(z, 25) ^ rotate64(z, 41);
//z += y ^ rotate64(y, 53) ^ rotate64(y, 3);
//y += z ^ rotate64(z, 31) ^ rotate64(z, 37);
//return y;

// LeaderRandom
// Passes 64TB of PractRand with no anomalies.
// Period is exactly 2 to the 128, all on one stream.
// This is one-dimensionally equidistributed, returning every uint64_t exactly as often.
// Requires access to the count leading zeroes function, which may vary based on compiler.
// This assumes behavior is the same as Java's Long.numberOfLeadingZeros() method, that is,
// an input of 0 produces an output of 64. The same period length is guaranteed if the
// count leading zeroes function here returns any even constant, but if it doesn't return
// 64 when given 0, the output will eventually be different for the same seed.
// Otherwise, this is a fairly ordinary RNG; it will probably be slower than a SplitMix64
// generator, but has a much longer period.
					// uint64_t r = (stateA += 0x9E3779B97F4A7C15ULL);
					// uint64_t q = (stateB += __lzcnt64(r));
					// r ^= r >> 27 ^ rotate64(q, 21);
					// r *= 0x3C79AC492BA7B653ULL;
					// r ^= r >> 33 ^ rotate64(q, 41);
					// r *= 0x1C69B3F74AC4AE35ULL;
					// r ^= r >> 27;
					// return r;

					//uint64_t b = (stateB += __builtin_ctzll(a));
					//return a ^ rotate64(a, 59) ^ rotate64(a, 8);

					//0xC6BC279692B5C323UL //0xCB9C59B3F9F87D4DUL

//// FleetRandom
//// Passes 64TB of PractRand with no anomalies.
//// Period is exactly 2 to the 128, all on one stream.
//// This is one-dimensionally equidistributed, returning every uint64_t exactly as often.
//// Requires access to the count leading zeroes function, which may vary based on compiler.
//// This assumes behavior is the same as Java's Long.numberOfLeadingZeros() method, that is,
//// an input of 0 produces an output of 64. The same period length is guaranteed if the
//// count leading zeroes function here returns any even constant, but if it doesn't return
//// 64 when given 0, the output will eventually be different for the same seed.

//					uint64_t q = (stateA += 0x9E3779B97F4A7C15ULL);
//					uint64_t r = (stateB += 0xC6BC279692B5C323ULL ^ __lzcnt64(q));
//					r ^= r >> 28;
//					r *= q | 1ULL;
//					r ^= r >> 30 ^ r >> 6;
//					return r;

//// FowlRandom
//// Passes at least 64TB with no anomalies. Passes Juniper's ICE test, typically.
//// Period is exactly 2 to the 128; one stream.
uint64_t x = (stateA += __lzcnt64(stateB));
uint64_t y = (stateB += 0xD1B54A32D192ED03ULL);
x = (x ^ rotate64(y, 37)) * 0x3C79AC492BA7B653ULL;
x = (x ^ x >> 33) * 0xF1357AEA2E62A9C5ULL;
x ^= x >> 27;
return x;

				}
				std::string mingler::get_name() const { return "mingler"; }
				void mingler::walk_state(StateWalkingObject *walker) {
					walker->handle(stateA);
					walker->handle(stateB);
//					if(stateB == 0) stateB = 1U;
//					walker->handle(incA);
//					walker->handle(incB);
					// if(incB == 0) incB = 1U;
//					stateA |= 1U;
//					stateB |= 1U;
					// incB |= 1U;

//					walker->handle(incA);
//					incA |= 1U;
//					walker->handle(incB);
//					incB *= 2U;
//					incB |= 2U;
					//stateB = (stateB & UINT64_C(0x7FFFFFFF)) + (stateB >> 31);
					//if (stateB == 0) stateB = 1;
					//stateA = 0;
					//stateB = 0;
//					printf("stateA is 0x%016LX, stateB is 0x%016LX, incA is 0x%016LX, incB is 0x%016LX\r\n", stateA, stateB, incA, incB);
					//printf("stateA is 0x%016LX, stateB is 0x%016LX\r\n", stateA, stateB);
				}

				Uint32 ta32::raw32() {
//					const Uint32 INCR = 0x6C8E9CF5, ADJ = 0x1;
//					enum { SHIFT_A = 13, SHIFT_B = 11 };
//					stateB ^= stateB >> 6;
//					stateA += INCR;
//					Uint64 z = (stateA ^ stateA >> SHIFT_A) * ((stateB ^= stateB << 1) | ADJ);
//					return (z ^ (z >> SHIFT_B)) + (stateB ^= stateB >> 11);

//					const uint32_t s = (stateA += 0xC1C64E6DU);
//				    uint32_t x = (s ^ s >> 16) * ((stateB += 0x9E3779BBU) | 1U);
//					if (!s) stateB -= 0x9E3779BBU;
//				    x = (x ^ x >> 15) * 0x846CA68BU;
//				    return x ^ x >> 16;

                    ////YOWZA, passes 32TB with no anomalies.
					////Variants came close to failing at 32TB; this doesn't have issues.
					//const uint32_t s = (stateA += 0xC1C64E6DU);
				    //uint32_t x = (s ^ s >> 17) * ((stateB += 0x9E3779BBU) | 1U);
					//if (!s) stateB -= 0x9E3779BBU;
					//x ^= x >> 16;
					//x *= UINT32_C(0xAC4C1B51);
					//x ^= x >> 15;
					//return x;

					////GWT-optimized SilkRNG; has passed 32TB with no anomalies.
					//const uint32_t s = (stateA += 0xC1C64E6DU);
				    //if (s) stateB += 0x9E3779BBU;
					//uint32_t x = (s ^ s >> 17) * (stateB >> 12 | 1U);
					//x ^= x >> 16;
					//x *= 0xAC451U;
					//return x ^ x >> 15;

					////Same as above but with branchless opt
					//const uint32_t s = (stateA += 0xC1C64E6DU);
				    //uint32_t x = (s ^ s >> 17) * ((stateB += -((s | -s) >> 31) & 0x9E3779BBU) >> 12 | 1U);
					//x ^= x >> 16;
					//x *= 0xAC451U;
					//return x ^ x >> 15;


					//const uint32_t s = (stateA = stateA + 0xC1C64E6DU | 0);
				    //uint32_t x = (s ^ s >> 17) * ((stateB = stateB - (-((s | -s) >> 31) & 0x279BU) | 0U) >> 12 | 1U);
					//x ^= x >> 16;
					//x *= 0xAC451U;
					//return x ^ x >> 15;


					//uint32_t x = (s ^ s >> 17) * (stateB | 0xFFF00001U);

					//x *= 0x7feb352dU;
					//x ^= x >> 15;
					//x *= 0x846ca68bU;
				    //uint32_t x = (s ^ s >> 17) * (t | 1U);
					////if (s >= 0xAC4C1B51U) stateB++;
					//return x ^ x >> 16 ^ x >> 11 ^ x >> 21;

//					x *= UINT32_C(0x7feb352d);
//					x ^= x >> 15;
//					x ^= x >> 11;
//					x *= UINT32_C(0x31848bab);
//					x ^= x >> 14;

//					// SpinelRandom, passes 64TB with no anomalies.
//					// Period is 2 to the 64. All states are valid.
//					uint32_t z = (stateA += 0xC13FA9A9U);
//					uint32_t y = (stateB += -((z | -z) >> 31) & 0x91E10DA5U);
//					z += y ^ rotate32(y, 3) ^ rotate32(y, 14);
//					y += z ^ rotate32(z, 28) ^ rotate32(z, 21);
//					z += y ^ rotate32(y, 12) ^ rotate32(y, 23);
//					y += z ^ rotate32(z, 10) ^ rotate32(z, 17);
//					return y;

					//uint32_t z = (stateA += 0xC13FA9A9U);
					////uint32_t y = (stateB += ((z | 0xAEF17502U - z) >> 31) + 0x91E10DA5U); //0xAEF17502 or 0xF1357AEA //0x91E10DA5U
					//uint32_t y = (stateB += ((z | 0xF1357AEAU - z) >> 31) + 0x91E10DA5U); //0xAEF17502 or 0xF1357AEA //0x91E10DA5U
					//z += y ^ rotate32(y, 3) ^ rotate32(y, 14);
					//y += z ^ rotate32(z, 28) ^ rotate32(z, 21);
					//z += y ^ rotate32(y, 12) ^ rotate32(y, 23);
					//y += z ^ rotate32(z, 10) ^ rotate32(z, 17);
					//return y;
					
					// Respite96
					// Passes 64TB of PractRand.
					// Period is 2 to the 96. Has 96 bits of state, with 32-bit output.
					// Does not use multiplication, but does use the count-leading-zeros instruction.
//					uint32_t a = (stateA += 0x91E10DA5U);
//					uint32_t b = (stateB += 0x6C8E9CF5U ^ __lzcnt(a));
//					uint32_t c = (stateC += 0x7FEB352DU ^ __lzcnt(a&b));
//					b = rotate32(b, 24) + a ^ c;
//					a = rotate32(a, 3) ^ b;
//					b = rotate32(b, 24) + a ^ c;
//					a = rotate32(a, 3) ^ b;
//					b = rotate32(b, 24) + a ^ c;
//					a = rotate32(a, 3) ^ b;
//					b = rotate32(b, 24) + a ^ c;
//					a = rotate32(a, 3) ^ b;
//					return a;

 					// Evasive96
					// Passes 64TB of PractRand with no anomalies.
					// Period is 2 to the 96. Has 96 bits of state, with 32-bit output.
					// Does not use multiplication, but does use the count-leading-zeros instruction.
					// Structured so that the step for w can be repeated for additional states, with the last step calling the variable w.
					// That means in an extended form, the w line would change to use z, the new line would use w, and this would use (x&=z) on the extra line.
					// This could also use |= , or mix that with &= .
					//uint32_t x = (stateA += 0x9E3779BB);
					//uint32_t y = (stateB += 0xC13FA9A9 ^ __lzcnt(x));
					//uint32_t w = (stateC += 0x915F77F5 ^ __lzcnt(x&=y));
					//w += x ^ rotate32(x, 3) ^ rotate32(x, 14);
					//x += w ^ rotate32(w, 28) ^ rotate32(w, 21);
					//w += x ^ rotate32(x, 12) ^ rotate32(x, 23);
					//x += w ^ rotate32(w, 10) ^ rotate32(w, 17);
					//return x;

					// Trying to actually extend Evasive96 to Evasive128...
					// Period is 2 to the 128. Has 128 bits of state, with 32-bit output.
					// This gets one anomaly at only 8GB,
					//   [Low4/16]BCFN(2+2,13-1,T)         R=  +9.0  p =  2.5e-4   unusual
					// It needs to not repeat "x &= whatever" operations too much, so this uses "|=".
					// With two "&=" ops in a row, this starts to see serious trouble at 8TB or earlier.
//					uint32_t x = (stateA += 0x9E3779BBU);
//					uint32_t y = (stateB += (0x9E3779BBU * 421U) ^ __lzcnt(x));
//					uint32_t z = (stateC += (0x9E3779BBU * 421U * 412U) ^ __lzcnt(x &= y));
//					uint32_t w = (stateD += (0x9E3779BBU * 421U * 412U * 421U) ^ __lzcnt(x |= z));
//					w += x ^ rotate32(x, 3) ^ rotate32(x, 14);
//					x += w ^ rotate32(w, 28) ^ rotate32(w, 21);
//					w += x ^ rotate32(x, 12) ^ rotate32(x, 23);
//					x += w ^ rotate32(w, 10) ^ rotate32(w, 17);
//					return x;

					// I'm not sure what to call this. It's an experiment, at least.
					// Passes 64TB with no anomalies. Period is 2 to the 128. 128 state bits. 32-bit output.
					//uint32_t n, x, y, z, w;
					//n = (x = (stateA += 0x9E3779BBU));
					//n = (y = (stateB += 0x6C8E9CF5U ^ __lzcnt(x))) + (n ^ rotate32(n, 3) ^ rotate32(n, 14));
					//n = (z = (stateC += 0x7FEB352DU ^ __lzcnt(x &= y))) + (n ^ rotate32(n, 10) ^ rotate32(n, 17));
					//n = (w = (stateD += 0x91E10DA5U ^ __lzcnt(x &= z))) + (n ^ rotate32(n, 21) ^ rotate32(n, 28));
					//return x + (n ^ rotate32(n, 12) ^ rotate32(n, 23));

					// Passes at least 16TB with no anomalies.
					// Without the rotations by 21, this has at least some minor anomalies at 16TB.
					//uint32_t x, y, z, w;
					//x = (stateA += 0x9E3779BBU);
					//y = (stateB += rotate32(x, 21) + __lzcnt(x));
					//z = (stateC += rotate32(y, 21) + __lzcnt(x &= y));
					//w = (stateD += rotate32(z, 21) + __lzcnt(x &= z));
					//w = (w + rotate32(x, 5)) * 0x21f0aaad;
					//w = (w ^ w >> 15) * 0x735a2d97;
					//return w ^ w >> 15;

// Eve32Random
// Passes 64TB with no anomalies; period is 2 to the 128; 32-bit output.
// This took a while to get to pass, and it still uses three multiplications...
//uint32_t x, y, z, w;
//x = (stateA += 0xDB4F0B91);
//y = (stateB += 0x9E3779BD * (rotate32(x, 21) + __lzcnt(x)));
//z = (stateC += 0x9E3779BD * (rotate32(y, 21) + __lzcnt(x &= y)));
//w = (stateD += 0x9E3779BD * (rotate32(z, 21) + __lzcnt(x &= z)));
//return x ^ w ^ rotate32(w, 7) ^ rotate32(w, 24);

// Berry32Random
// Passes 64TB with no anomalies; period is 2 to the 128; 32-bit output.
// Uses a unary hash from skeeto/hash-prospector, prospector32 .
// Because of https://github.com/skeeto/hash-prospector/issues/28 , I tried a different hash that didn't rate as highly.
// Does not use any multiplication before the unary hash.
// Does use the 32-bit count-leading-zeros instruction (which should be an intrinsic).
//uint32_t x, y, z, w;
//x = (stateA += 0xDB4F0B91);
//y = (stateB += (rotate32(x, 21) + __lzcnt(x)));
//z = (stateC += (rotate32(y, 21) + __lzcnt(x &= y)));
//w = (stateD += (rotate32(z, 21) + __lzcnt(x &= z)));
//x += rotate32(w, 21);
//x ^= x >> 15;
//x *= 0x2c1b3c6d;
//x ^= x >> 12;
//x *= 0x297a2d39;
//x ^= x >> 15;
//return x;

// Bear32Random, 1.0
// Passes 64TB with no anomalies; period is 2 to the 128; 32-bit output.
// Uses (most of) a unary hash from skeeto/hash-prospector, prospector32, but elides the first xorshift.
// Does not use any multiplication before the unary hash.
// Does use the 32-bit count-leading-zeros instruction (which should be an intrinsic).
//uint32_t x, y, z, w;
//x = (stateA += 0x9E3779B9);
//y = (stateB += (x + __lzcnt(x)));
//z = (stateC += (y + __lzcnt(x &= y)));
//w = (stateD += (z + __lzcnt(x &= z)));
//w +=  rotate32(x, 21);
//w *= 0x2C1B3C6D;
//w ^= w >> 12;
//w *= 0x297A2D39;
//w ^= w >> 15;
//return w;

// Bear32Random, 2.0
// Meant to improve decorrelation between similar initial states more quickly.
// The main difference with 1.0 is that this uses XOR where 1.0 used addition, and a rotation changed.
// Passes 64TB with no anomalies; period is 2 to the 128; 32-bit output.
// Uses (most of) a unary hash from skeeto/hash-prospector, prospector32, but elides the first xorshift.
// Does not use any multiplication before the unary hash.
// Does use the 32-bit count-leading-zeros instruction (which should be an intrinsic).
//uint32_t x, y, z, w;
//x = (stateA += 0x9E3779B9);
//y = (stateB += (x ^ __lzcnt(x)));
//z = (stateC += (y ^ __lzcnt(x &= y)));
//w = (stateD += (z ^ __lzcnt(x &= z)));
//w += rotate32(x, 13);
//w *= 0x2C1B3C6D;
//w ^= w >> 12;
//w *= 0x297A2D39;
//w ^= w >> 15;
//return w;

// Rawr32Random
// Passes 64TB with no anomalies. Minimum period is 2 to the 64. Maximum period is 2 to the 128, but almost guaranteed to be less than that.
// This one's quite finicky about rotation amounts.
// All states are updated independently of the returned value calculation.
// Uses no multiplication; does use  32-bit count-leading-zeros instruction (which should be an intrinsic).
//uint32_t x = stateA, y = stateB, z = stateC, w = stateD;
//stateA = x + 0x9E3779B9;
//stateB = y + __lzcnt(x);
//stateC = rotate32(w, 3)  - x;
//stateD = rotate32(z, 24) ^ y;
//x ^= rotate32(w, 29) + y;
//y += rotate32(z, 19) ^ x;
//return x - rotate32(y, 21);


//					z += y ^ rotate32(y, 3) ^ rotate32(y, 14);
//					y += z ^ rotate32(z, 28) ^ rotate32(z, 21);
//					z += y ^ rotate32(y, 12) ^ rotate32(y, 23);
//					y += z ^ rotate32(z, 10) ^ rotate32(z, 17);
//					return y;
					//x ^= x >> 17;
					//x *= 0xed5ad4bb;
					//x ^= x >> 11;
					//x *= 0xac4c1b51;
					//x ^= x >> 15;
					//x *= 0x31848bab;
					//x ^= x >> 14;
					//return x;

					//x = (x ^ x >> 15) * 0x2c1b3c6d;
					//x = (x ^ x >> 12) * 0x297a2d39;
					//x ^= x >> 15;
					//return x;

//					x = (x ^ x >> 15) * 0xd168aaad;
//					x = (x ^ x >> 15) * 0xaf723597;
//					x ^= x >> 15;
//					return x;

// Passes 2TB without anomalies, but fails around 8TB.
//uint32_t x = stateA, y = stateB, z = stateC, w = stateD, n = x;
//stateA = x + 0x9E3779B9;
//stateB += x + (__lzcnt(n));
//stateC += y + (__lzcnt(n &= y));
//stateD += z + (__lzcnt(n & z));
//w = (w ^ rotate32(z, 19) + (y ^ rotate32(x, 7) + w) + rotate32(y, 27)) + rotate32(w, 11);
//return w;

//stateA += 0x6C8E9CF5U;
//stateB += 0x7FEB352DU + (__lzcnt(x));
//stateC += 0x91E10DA5U + (__lzcnt(x & y));
//uint32_t x = stateA, y = stateB, z = stateC;
////0xD1B54A32
//stateA = 4 ^ x + 0x9E3779BD;
//stateB = x ^ y + (__lzcnt(x));
//stateC = y ^ z + (__lzcnt(x & y));
//x = (y ^ rotate32(x, 7) + z);
//y = x + rotate32(y, 27);
//x = rotate32(x, 29) ^ (z * y | 1);
//y = rotate32(y,  5) ^ (x * z | 1);
//z = rotate32(z, 14) ^ (y * x | 1);
//x = (y ^ rotate32(x, 19) + z);
//y = x + rotate32(y, 11) ^ z;

//x = (y ^ rotate32(x, 7) + z);
//y = (x + rotate32(y, 27));
//
//x = (y ^ rotate32(x, 19) + z);
//y = (x + rotate32(y, 11));

// Passes 64TB with no anomalies.
//x ^= rotate32(z, 3) + rotate32(y, 24);
//y ^= rotate32(x, 5) + rotate32(z, 15);
//z ^= rotate32(y, 9) + rotate32(x, 29);
//y = (z ^ rotate32(y, 19) + x);
//z = (y + rotate32(z, 11));

//x ^= (z * y | 1);
//y ^= (x * z | 1);
//z ^= (y * x | 1);

//// This block seems to have bad problems starting at 64TB, which is really annoying to test.
//x ^= rotate32(z, 11) + rotate32(y, 19);
//y ^= rotate32(x, 7) + rotate32(z, 23);
//z ^= rotate32(y, 15) + rotate32(x, 21);
//return x^y^z;


//x ^= rotate32(z, 13) + rotate32(y, 23);
//y ^= rotate32(x, 7) + rotate32(z, 19);
//z ^= rotate32(y, 5) + rotate32(x, 29);
//return z;

//z = ((y = (z ^ rotate32(y, 3) + (x ^ 0x6C8E9CF5))) + rotate32(z, 24));
//z = ((y = (z ^ rotate32(y, 3) + (x ^ 0x7FEB352D))) + rotate32(z, 24));
//z = ((y = (z ^ rotate32(y, 3) + (x ^ 0x91E10DA5))) + rotate32(z, 24));

//z = ((y = (z ^ rotate32(y,  7) + x)) + rotate32(z, 24));
//z = ((y = (z ^ rotate32(y,  5) + x)) + rotate32(z, 24));
//z = ((y = (z ^ rotate32(y, 12) + x)) + rotate32(z, 24));

//x ^= rotate32(z, 13) + rotate32(y, 23);
//y ^= rotate32(x, 17) + rotate32(z, 19);
//z ^= rotate32(y,  5) + rotate32(x, 29);


//uint32_t x = stateA, y = stateB, z = stateC;
////0xD1B54A32
//stateA = 4 ^ x + 0x9E3779BD;
//stateB = x ^ y + (__lzcnt(x));
//stateC = y ^ z + (__lzcnt(x & y));
//z = ((y = (z ^ rotate32(y, 3) + (x ^ 0x6C8E9CF5))) + rotate32(z, 24));
//z = ((y = (z ^ rotate32(y, 29) + (x ^ 0x7FEB352D))) + rotate32(z, 8));
//return z;


//uint32_t x, y;
//x = (stateA += 0xDB4F0B91);
//y = (stateB += (rotate32(x, 21) + __lzcnt(x)));
//x ^= rotate32(y, 21);
//x ^= x >> 15;
//x *= 0x2c1b3c6d;
//x ^= x >> 12;
//x *= 0x297a2d39;
//x ^= x >> 15;
//return x;

// 0x39E2D

// Passes at least 4TB with no anomalies.
//uint32_t x = (stateA = stateA + 0x9E3779BD ^ 0xD1B54A32);
//uint32_t t = x & 0xDB4F0B96 - x;
//uint32_t y = (stateB = stateB + rotate32(t, 1) ^ 0xAF723597);
////uint32_t y = (stateB = stateB + 0x91E10DA5 ^ (x & 0xDB4F0B92 - x) >> 31);
//y = rotate32(y, x >> 27);
//y = (y ^ (y >> ((y >> 28u) + 4u))) * 0xB45ED;
//return y ^ y >> 21;

// Passes at least 16TB with no anomalies (ongoing).
//uint32_t x = (stateA = stateA + 0x9E3779BD ^ 0xD1B54A32);
//uint32_t t = x & 0xDB4F0B96 - x;
//uint32_t y = (stateB = stateB + rotate32(t, 1) ^ 0xAF723597);
//y ^= rotate32(x, 17);
//y = (y ^ (y >> ((y >> 28u) + 4u))) * 0xB45ED; //0x39E2D;
//return y ^ y >> 21;

// Gets "mildly suspicious" at 32TB
//uint32_t x = (stateA = stateA + 0x9E3779BD ^ 0xD1B54A32);
//uint32_t t = x & 0xDB4F0B96 - x;
//uint32_t y = (stateB = stateB + rotate32(t, 1) ^ 0xAF723597);
//y += x;
//y = (y ^ (y >> ((y >> 28u) + 4u))) * 0xB45ED;
//return y ^ y >> 22;

/*
rng=ta32, seed=0x0
length= 128 gigabytes (2^37 bytes), time= 412 seconds
  Test Name                         Raw       Processed     Evaluation
  DC6-9x1Bytes-1                    R= +43.2  p =  8.1e-16    FAIL !
  ...and 879 test result(s) without anomalies
*/
//uint32_t x = (stateA = stateA + 0x9E3779BD ^ 0xD1B54A32);
//uint32_t t = x & 0xDB4F0B96 - x;
//uint32_t y = (stateB = stateB + rotate32(t, 1) ^ 0xAF723597);
//y = (y ^ y >> 13 ^ x ^ x >> 19) * 0xB45ED;
//return y ^ y >> 21;

// TaxonRandom
// Completes 64TB with no anomalies; 64 bits of state over two 32-bit words.
// Every pair of states appears once per cycle, which has length 2 to the 64.
// Uses at most 32-bit math. Does use one integer multiply, but follows JS rules for it.
// Uses only operations that are safe on JS Number types, and won't lose precision.
// The state transition is hard to classify, but has a full period.
// Interestingly, (x^y) is guaranteed to take every possible value equally often, but
// it "clumps" into irregularly-sized groups of even and odd values. The other steps
// performed help erase any pattern from the "clumping."
// Uses an operation from PCG-Random to randomly shift right.
// However, this fails juniper's ICE test for correlation of similar initial states.

//uint32_t x = (stateA = stateA + 0x9E3779BD ^ 0xD1B54A32);
//uint32_t t = x & 0xDB4F0B96 - x;
//uint32_t y = (stateB = stateB + rotate32(t, 1) ^ 0xAF723597);
//y = (x ^ y ^ (y >> ((y >> 28u) + 4u))) * 0xB45ED;
//return y ^ y >> 21;

// TaxonRandom, take 2
// This unusual variant on the previous TaxonRandom also passes 64TB with no anomalies.
// It has the same state size (64 bits) and cycle length (2 to the 64).
// This doesn't use anything from PCG-Random anymore; it does a random rotation of x using the
// lower 5 bits of y, and adds that to y, which doesn't interfere with being 1D equidistributed.
// There's a mixed left/right shift pair, which uses data found by Pelle Evensen:
// https://github.com/pellevensen/bijections
// This "take 2" was needed because "take 1" failed juniper's ICE test. This passes it.

//uint32_t x = (stateA = stateA + 0x9E3779BD ^ 0xD1B54A32);
//uint32_t t = x & 0xDB4F0B96 - x;
//uint32_t y = (stateB = stateB + rotate32(t, 1) ^ 0xAF723597);
//y += rotate32(x, y);
//y = (y ^ y >> 22 ^ y << 5) * 0xB45ED;
//return y ^ y >> 21;

// Choo32Random
// Passes 64TB with no anomalies.
// Period is at minimum 2 to the 32. State size is 128 bits over four 32-bit variables.
// This is an adapted version of ChopRandom to be used essentially as a hash of four inputs that
// can be run for longer to get more hashed bits.
//const uint32_t fa = stateA;
//const uint32_t fb = stateB;
//const uint32_t fc = stateC;
//const uint32_t fd = stateD;
//stateA = fb - fc;
//stateB = fa ^ fd;
//stateC = rotate32(fb, fa);
//stateD = fd + 0xADB5B165;
//uint32_t res = stateA + stateB + stateC + stateD;
//res = (res ^ res >> 15) * 0x735a2d97;
//return res ^ res >> 16;

// rng=ta32, seed=0x1
// length = 256 gigabytes(2 ^ 38 bytes), time = 838 seconds
// Test Name                         Raw       Processed     Evaluation
// [Low1 / 16]FPF - 14 + 6 / 16:(0, 14 - 0)     R = +16.7  p = 4.6e-15    FAIL
// [Low1 / 16]FPF - 14 + 6 / 16 : (1, 14 - 0)     R = +8.0  p = 5.1e-7   mildly suspicious
// [Low1 / 16]FPF - 14 + 6 / 16 : all          R = +9.3  p = 3.6e-8   very suspicious
// [Low4 / 16]FPF - 14 + 6 / 16 : (0, 14 - 0)     R = +7.6  p = 1.1e-6   unusual
// [Low4 / 16]FPF - 14 + 6 / 16 : all          R = +6.2  p = 2.7e-5   unusual
// ...and 910 test result(s) without anomalies

//uint32_t x = (stateA = stateA + 0x9E3779BD ^ 0xD1B54A32);
//uint32_t y = (stateB = stateB + __lzcnt(x) ^ x);
//y += rotate32(x, y);
//y = (y ^ rotate32(y, 21) ^ rotate32(y, 9)) * 0x91E10DA5;
//return y ^ y >> 16;

// Passes 32TB with no anomalies, but gets an "unusual" at 64TB:
// rng=ta32, seed=0x0
// length= 64 terabytes (2^46 bytes), time= 179223 seconds
// Test Name                         Raw       Processed     Evaluation
// [Low1/16]FPF-14+6/16:(0,14-0)     R=  +7.5  p =  1.3e-6   unusual
// ...and 1158 test result(s) without anomalies
//uint32_t x = (stateA = stateA + 0x9E3779BD ^ 0xD1B54A32);
//uint32_t y = (stateB = stateB + __lzcnt(x) ^ x * 0x91E10DA5);
//y += rotate32(x, y);
//y = (y ^ rotate32(y, 10) ^ rotate32(y, 23)) * 0xC5F768E7;
//return y ^ y >> 16;


// Taupe32Random
// Passes 64TB with no anomalies.
// State size is just 64 bits in two 32-bit variables.
// Period is exactly 2 to the 64. No streams.
//uint32_t x = (stateA = stateA + 0x9E3779BD ^ 0xD1B54A32);
//uint32_t y = (stateB = stateB + __lzcnt(x) * 0x91E10DA5 ^ x);
//y += rotate32(x, y);
//y = (y ^ y >> 15) * 0xD168AAAD;
//return y ^ rotate32(y, 11) ^ rotate32(y, 23);

// Nu32Random
// Passees 64TB of PractRand, but fails ICE tests badly.
// Uses no multiplication!
// GWT-safe!
// Period is 2 to the 64. No streams.
//uint32_t x = (stateA = stateA + 0x9E3779BD ^ 0xD1B54A32);
//uint32_t t = x & 0xDB4F0B96 - x;
//uint32_t y = (stateB = stateB + rotate32(t, 1) ^ 0xAF723596);
//x += y = rotate32(y, x);
//y +=     rotate32(x, y);
//return (y ^ rotate32(y, 3) ^ rotate32(y, 24));

// Should be more robust, but is not at all!
// rng=ta32, seed=0x0
// length= 4 terabytes (2^42 bytes), time= 10580 seconds
//   Test Name                         Raw       Processed       Evaluation
//   [Low1/8]FPF-14+6/16:(0,14-0)      R= +10.8  p =  1.3e-9     very suspicious
//   [Low1/8]FPF-14+6/16:all           R=  +5.5  p =  1.1e-4     unusual
//   [Low1/16]FPF-14+6/16:(0,14-0)     R= +25.2  p =  5.8e-23    FAIL !!
//   [Low1/16]FPF-14+6/16:(1,14-0)     R= +16.0  p =  1.9e-14    FAIL
//   [Low1/16]FPF-14+6/16:all          R= +14.4  p =  4.9e-13    FAIL
//   [Low1/32]FPF-14+6/16:cross        R=  +7.2  p =  2.0e-6     mildly suspicious
//   [Low1/32]TMFn(2+5):wl             R= +21.4  p~=    1e-6     unusual
//   [Low4/16]FPF-14+6/16:(0,14-0)     R= +17.1  p =  1.8e-15    FAIL
//   [Low4/16]FPF-14+6/16:(1,14-0)     R= +10.8  p =  1.4e-9     very suspicious
//   [Low4/16]FPF-14+6/16:all          R=  +9.7  p =  1.2e-8     very suspicious
//   ...and 1042 test result(s) without anomalies

//uint32_t x = (stateA = stateA + 0x9E3779BD ^ 0xD1B54A32);
//uint32_t t = x & 0xDB4F0B96 - x;
//uint32_t y = (stateB = stateB + rotate32(t, 1) ^ 0xAF723596);
//x += y ^= rotate32(y, x) ^ rotate32(y, 31-x);
//y += x ^= rotate32(x, y) ^ rotate32(x, 31-y);
//return x ^ rotate32(y, 13);

// Gnome32Random
// Uses no multiplication, and still manages to look hash-like immediately.
// Actually uses only ARX operations plus subtraction, though the variable-distance rotations may have a different cost than fixed ones.
// Passes 64TB of PractRand with no anomalies, and passes both ICE tests in juniper.
// Period is 2 to the 64, no streams. 1D-equidistributed.
// The XRR construction is interesting, using two dependent random rotations that are always different.
// Among other things, XORing y with two different rotations of y means the result is never effectively rotated by 0.
// This doesn't mean it can't return y ever; if y is 0 or 0xFFFFFFFF, then any rotation will be the same.
// The state transition is the same as Taxon32Random except 0xAF723597 was changed to 0xAF723596, out of an abundance of caution.

//uint32_t x = (stateA = stateA + 0x9E3779BD ^ 0xD1B54A32);
//uint32_t t = x & 0xDB4F0B96 - x;
//uint32_t y = (stateB = stateB + rotate32(t, 1) ^ 0xAF723596);
//x += y ^= rotate32(y, x) ^ rotate32(y, 31-x);
//y += x ^  rotate32(x, y) ^ rotate32(x, 31-y);
//return (y ^ rotate32(y, 3) ^ rotate32(y, 24));

// This version gets a serious anomaly at 4TB, and most similar generators fail by then or earlier.
// rng=ta32, seed=0x0
// length= 4 terabytes (2^42 bytes), time= 10724 seconds
//   Test Name                         Raw       Processed     Evaluation
//   [Low4/16]FPF-14+6/16:(0,14-0)     R=  +9.4  p =  2.8e-8   suspicious
//   ...and 1051 test result(s) without anomalies
//
//uint32_t x = (stateA = stateA + 0x9E3779BD ^ 0xD1B54A32);
//uint32_t t = x & 0xDB4F0B96 - x;
//uint32_t y = (stateB = stateB + rotate32(t, 1) ^ 0xAF723596);
//x ^= x >> 15 ^ y;
//x *= 0x2c1b3c6d;
//x ^= x >> 12;
//x *= 0x297a2d39;
//x ^= x >> 15;
//return x;

// Fluff32Random, v1
// Passes 64TB with no anomalies.
// Has a period of 2 to the 64, no streams. Two 32-bit states.
// Very similar to FlowRandom, but with 32-bit states and an extended period.
//uint32_t x = (stateA = stateA + 0x9E3779BD ^ 0xD1B54A32);
//uint32_t y = (stateB = stateB + __lzcnt(x) + 0x91E10DA5 ^ x);
//x ^= rotate32(y, 14);
//x *= 0x2c1b3c6d;
//x ^= x >> 12;
//x *= 0x297a2d39;
//x ^= x >> 15;
//return x;

// NOOOOOOO
//rng=ta32, seed=0x0
//length= 64 terabytes (2^46 bytes), time= 197858 seconds
//  Test Name                         Raw       Processed     Evaluation
//  FPF-14+6/16:cross                 R=  -2.5  p =1-1.6e-4   unusual
//  ...and 1158 test result(s) without anomalies
//uint32_t x = (stateA = stateA + 0x9E3779BD ^ 0xD1B54A32);
//uint32_t y = (stateB = stateB + __lzcnt(x) + 0xC5F768E7 ^ x);
//x ^= rotate32(y, 17);
//x *= 0x2c1b3c6d;
//x ^= x >> 12;
//x *= 0x297a2d39;
//x ^= x >> 15;
//return x;

// ADDITIONAL NOOOOOOOO
//rng=ta32, seed=0x0
//length= 64 terabytes (2^46 bytes), time= 177964 seconds
//  Test Name                         Raw       Processed     Evaluation
//  [Low4/16]FPF-14+6/16:(0,14-0)     R=  +8.2  p =  3.5e-7   mildly suspicious
//  ...and 1158 test result(s) without anomalies
// This one only had a problem at the very end... 64TB...
//uint32_t x = (stateA = stateA + 0x9E3779BD ^ 0xD1B54A32);
//uint32_t t = x & 0xDB4F0B96 - x;
//uint32_t y = (stateB = stateB + rotate32(t, 1) ^ 0xAF723596);
//x += (y ^= rotate32(x, 22));
//x ^= rotate32(x, 3) ^ rotate32(x, 24);
//y += (x ^= rotate32(y, 10));
//y ^= rotate32(y, 11) ^ rotate32(y, 26);
//return x ^ y;

//uint32_t x = (stateA = stateA + 0x9E3779BD ^ 0xD1B54A32);
//uint32_t t = x & 0xDB4F0B96 - x;
//uint32_t y = (stateB = stateB + rotate32(t, 1) ^ 0xAF723596);
//x = rotate32(x, 3) ^ (y = rotate32(y, 24) + x ^ 0x2C1B3C6D);
//x = rotate32(x, 3) ^ (y = rotate32(y, 24) + x ^ 0x297A2D39);
//x = rotate32(x, 3) ^ (y = rotate32(y, 24) + x ^ 0x91E10DA5);
//return x;

// Passes at least 4TB with no anomalies
//return stateA = (stateB = rotate32(stateB, 23) + (stateC = stateC + 0xD192ED03U | 0U) | 0U) ^ rotate32(stateC, 19) + stateA ^ rotate32(stateA, 22);

//// Chock32Random
//// Passes 64TB with no anomalies.
//const uint32_t fa = stateA;
//const uint32_t fb = stateB;
//const uint32_t fc = stateC;
//const uint32_t fd = stateD;
//stateA = fb - fc;
//stateB = fa ^ fd;
//stateC = rotate32(fb, 11);
//stateD = fd + 0xADB5B165;
//uint32_t res = stateA + stateB;
//return res ^ rotate32(res, 14) ^ rotate32(res, 23);

// Chip32Random
// Passes 64TB with no anomalies.
const uint32_t fa = stateA;
const uint32_t fb = stateB;
const uint32_t fc = stateC;
const uint32_t fd = stateD;
stateA = fb + fc;
stateB = fd ^ fa;
stateC = rotate32(fb, 11);
stateD = fd + 0x9E3779B9;
return rotate32(fa, 14) ^ rotate32(fb, 23) + fc;

				}
				std::string ta32::get_name() const { return "ta32"; }
				void ta32::walk_state(StateWalkingObject *walker) {
					walker->handle(stateA);
					walker->handle(stateB);
					walker->handle(stateC);
					walker->handle(stateD);
					//stateB |= (stateB == 0);
					printf("Seed is 0x%08X, 0x%08X, 0x%08X, 0x%08X\r\n", stateA, stateB, stateC, stateD);
				}


				Uint64 thrust63::raw64() {
					//WORKS WELL, FAST, NOT EQUIDISTRIBUTED
					// fails at 8TB on binrank
					// 0x3C6EF372FE94F82A , 21 , 28 
					// 0x4F1BBCDCBFA53E0A , 23 , 28 these
					// 0x4F1BBCDCBFA53E0A , 21 , 28 all
					// 0x4F1BBCDCBFA53E0A , 23 , 30 fail binrank
					//uint64_t z = (state += UINT64_C(0xBC6EF372FE94F82A));
					//z *= ((z >> 25) ^ z);
					//return (z >> 28) ^ z;

					//// passes 32TB, but with issues at the end, no failures though
					//uint64_t z = (state += UINT64_C(0x3C6EF372FE94F82A));
					//z *= ((z >> 21) ^ z);
					//return z - (z >> 29);
					state *= state + (state >> 26); state ^= state >> 31; return state;

					//uint64_t z = (state += UINT64_C(0x9E3779B97F4A7C16)); z *= (z >> 25 ^ z); return z - (z >> 22); // 0x6C6EF372FB7CD85A

                    //uint64_t z = (state += UINT64_C(0x3C6EF372FE94F82A));
                    //z *= ((z >> 21) ^ z);
                    //return z - (z >> 25);
                    // 0xD91D39EAE12657AA
                    //0xAD91D39EAE12657A
                    // 0x4F1BBCDCBFA53E0A maybe?
                    // 0x6A5D39EAE127586AULL works, but fails some gjrand tests
                    // 0x6A5D39EAE12657AAULL works on both gjrand and PractRand, with shifts 25 then 22
                    //const uint64_t s = state;
                    //const uint64_t z = (s ^ (s >> 25)) * (state += 0x6A5D39EAE12657AAULL);// 0x6A5D39EAE116586AULL);// 0x6C138EB769B56EBAULL); // 0x6A5D39EAE126579AULL
                    //return z ^ (z >> 22);
                    //z = (s ^ (s >> 25)) * (j += 0x7C139EB769B97F32ULL);
                    //z ^= (z >> 25);

				}
				std::string thrust63::get_name() const { return "thrust63"; }
				void thrust63::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
					//state |= 1;
				}
				void thrust63::seed(Uint64 seed)
				{
//					shift2 = seed & 63u;
//					shift1 = seed >> 8 & 63u;
					autoseed();
					state |= 1;
				}

				Uint64 splitmix64::raw64() {
					//const Uint64 inc = 0x9E3779B97F4A7C15;
					//const Uint64 stage1 = 0xBF58476D1CE4E5B9;
					//const Uint64 stage2 = 0x94D049BB133111EB;
					Uint64 z = (state += 0x222233B8B7C7F28B);
					z = (z ^ z >> 30) * 0xBF58476D1CE4E5B9;
					z = (z ^ z >> 27) * 0x94D049BB133111EB;
					return z ^ z >> 31;
				}
				std::string splitmix64::get_name() const { return "splitmix64"; }
				void splitmix64::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				
				//// Be careful, has late-occuring TMFn issues that can sometimes result in failures.
				//// Not equidistributed at all; cannot produce some outputs and produces others more often.
				//// Period is 2 to the 64; all states are allowed.
				Uint64 thrustAlt64::raw64() {
					Uint64 z = (state += UINT64_C(0x6C8E9CF570932BD5));
					//z = (z ^ z >> 25) * z - (z ^ z >> 30) * (z + 0xEB44ACCAB455D165ULL);
					z = (z ^ z >> 25) * (z | 0xA529L);
					return z ^ (z >> 23);
				}
				std::string thrustAlt64::get_name() const { return "thrustAlt64"; }
				void thrustAlt64::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}
				Uint32 mulberry32::raw32() {
					//uint32_t z = (j += 0x6D2B79F5UL) ^ 0x9E3779B5;
					//z = (z ^ z >> 11) * 0xFFF + 0x9E3779B5;
					//z = (z ^ z >> 13) * 0x1FF - 0x9E3779B5;
					//z = (z ^ z >> 12) * 0x3FF + 0x9E3779B5;
					//z = (z ^ z >> 14) * 0x7FF - 0x9E3779B5;
					//return z ^ z >> 15;
					//uint32_t z = (j += 0x9E3779B9);
					//        z = (z ^ z >> 13) * 0x7FFFF;
					//        z = (z ^ z >> 12) * 0x1FFFF;
					//        z = (z ^ z >> 14) * 0x1FFF;
					//return z ^ z >> 12;

					uint32_t z = (j += 0xB79F5u);
					z = (z ^ (z >> 15)) * (z | 0xFFF0003D);
					z ^= z >> 8;
					z = ((z ^ (z >> 7)) * (j >> 12 | 1u));
					return z ^ (z >> 14);

					//uint32_t z = (j += 0x6D2B79F5u);
					//z = (z ^ (z >> 15)) * (z | 61u);
					//z ^= z >> 8;
					//z += ((z ^ (z >> 7)) * (j | 0xA529u));
					//return z ^ (z >> 14);

				}
				std::string mulberry32::get_name() const { return "mulberry32"; }
				void mulberry32::walk_state(StateWalkingObject *walker) {
					walker->handle(j);
				}

				Uint64 vortex::raw64() {
					/*
					uint64_t z = (state += UINT64_C(0x6C8E9CF570932BD5)) ^ ((stream *= UINT64_C(0x2545F4914F6CDD1D)) >> 30);
					z ^= (z >> 25) * (stream + UINT64_C(0x6A5D39EAE12657BA));
					return z ^ (z >> 28) ^ (state & stream);
					*/
					uint64_t z = (state += UINT64_C(0x6C8E9CF570932BD5));
					z ^= (z >> 26) * (stream *= UINT64_C(0x2545F4914F6CDD1D));
					return z ^ (~z) >> 26; //UINT64_C(0x9E3779B97F4A7BB5) - 
				}
				std::string vortex::get_name() const { return "vortex"; }
				void vortex::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
					walker->handle(stream);
					stream |= 1;
				}
				void vortex::seed(Uint64 seed)
				{
					state = seed & 0xFFFFUL;
					state += UINT64_C(0x6C8E9CF570932BD5);
					state = (state ^ state >> 25) * (state | UINT64_C(0xA529));
					state ^= state >> 22;

					seed >>= 16;

					seed += UINT64_C(0x6C8E9CF570932BD5);
					seed = (seed ^ seed >> 25) * (seed | UINT64_C(0xA529));
					stream = (seed ^ seed >> 22) | UINT16_C(1);
				}

				Uint64 vortex2::raw64() {
					// 0xA529 UINT64_C(0x9E3779B97F4A7BB5)
					state += UINT64_C(0x6C8E9CF570932BD5);
					uint64_t z = (state ^ state >> 25) * ((stream += UINT64_C(0x9E3779B97F4A7BB5)) | UINT64_C(1));
					return (z ^ (z >> 28));
				}
				std::string vortex2::get_name() const { return "vortex2"; }
				void vortex2::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
					walker->handle(stream);
				}
				void vortex2::seed(Uint64 seed)
				{
					state = seed;
					seed += UINT64_C(0x6C8E9CF570932BD5);
					seed = (seed ^ seed >> 25) * (seed | UINT64_C(0xA529));
					stream = seed ^ seed >> 22;
				}
				Uint64 vortex64::raw64() {
					//with just multiply
					//0x2545F4914F6CDD1D
					//with odd increment
					//0x369DEA0F31A53F85

					//0x6C8E9CF570932BD5 0x6C8E9CF970932BD5

					Uint64 z = (state += UINT64_C(0x9E3779B97F4A7C15)); //0x6C8E9CF570932BD5
																		//z = (z ^ z >> 25);// *UINT64_C(0x2545F4914F6CDD1D);
																		//z ^= z >> 13;
																		//z ^= z << 7;

																		//// passes at least 8TB no anomalies
																		//z = (z ^ z >> 13) * UINT64_C(0x7FFFF);// *UINT64_C(0xAEF17502108EF2D9);// 0x2545F4914F6CDD1D 0x9E3779B97F4A7C15
																		//z = (z ^ z >> 12) * UINT64_C(0x1FFFF);
																		//z = (z ^ z >> 14) * UINT64_C(0x1FFF);
																		//return z ^ z >> 24;

																		//z ^= z >> 9;
																		//z -= (z << 13) + UINT64_C(0xAEF17502108EF2D9);
																		//z ^= z >> 7;
																		//z -= (z << 11) + UINT64_C(0x6C8E9CF570932BD5);
					z ^= z >> 9 ^ 0x94D049BB133111EB;
					z -= (z << 13) + UINT64_C(0xAEF17502108EF2D9);
					z ^= z >> 7 ^ 0xBF58476D1CE4E5B9;
					z -= (z << 11) + UINT64_C(0x6C8E9CF570932BD5);

					//z -= z << 17;
					//z -= z << 19;
					//z = ((z << 14) - rotate64(z, 12));
					//z = ((z << 20) - rotate64(z, 17));
					return z ^ z >> 25; //((z >> 22) ^ rotate64(z, 44))
										//z ^= ((z << 19) | (z >> 45)) ^ ((z << 53) | (z >> 11));
										/*
										long z = (state += 0x9E3779B97F4A7C15L);
										z = (z ^ z >> 13) * 0x7FFFF;
										z = (z ^ z >> 12) * 0x1FFFF;
										z = (z ^ z >> 14) * 0x1FFF;
										return z ^ z >> 24;
										*/

										//Uint64 a = (state = state * Uint64(0x369DEA0F31A53F85) + Uint64(0x6C8E9CF970932BD5));

										//Uint64 a = (state += stream);
										////a ^= rotate64(a, a & 63) ^ rotate64(a, a + 39 & 63);
										//a = (a ^ a >> 24) * Uint64(0x2545F4914F6CDD1D);
										//return a ^ a >> 25;


										//Uint64 a = (state += 0x6C8E9CF970932BD5L);
										//Uint64 b = _byteswap_uint64(a);
										//a = (a ^ a >> 25) * Uint64(0x2545F4914F6CDD1D);
										//b ^= b >> 22 ^ a ^ a >> 28;
										//return (b ^ b >> 25);



										//z = (z ^ _byteswap_uint64(z) ^ z >> 25) * Uint64(0x2545F4914F6CDD1D);
										//return (z ^ z >> 22);

										//					Uint64 z = ((state ^ state >> 25)) * UINT64_C(0x59071D96D81ECD35);
										//					return z ^ ((z << 12) | (z >> 52)) ^ ((z << 47) | (z >> 17)) ^ _byteswap_uint64(state += stream);
										//					return z ^ z >> 25 ^ ((z << 12) | (z >> 52)) ^ ((z << 47) | (z >> 17));
										//Uint64 z = (state ^ state >> 25) ^ ((state << 43) | (state >> 21)) ^ ((state << 17) | (state >> 47));
										//return -z ^ ((z << 19) | (z >> 45)) ^ ((z << 53) | (z >> 11));

										//return z ^ z >> 25;

										//passes 32TB with no anomalies (!!!) with stream 12 and seed 0x33f3ca7e
										//state += stream; //0x6C8E9CF570932BD5 0x6C8E9CF970932BD5
										//Uint64 z = (state ^ state >> 25) * UINT64_C(0x2545F4914F6CDD1D);
										//z ^= ((z << 19) | (z >> 45)) ^ ((z << 53) | (z >> 11));
										//return z ^ z >> 25;

										//works at 16TB, probably further (very few anomalies, none worse than unusual)
										//state += stream;
										//Uint64 z = (state ^ state >> 25) * (state | Uint64(0xA529));
										//return z ^ z >> 22;

										// works, passes PractRand at 32TB with stream upper bits 12
										//state += stream;
										//Uint64 z = (state ^ state >> 25) * (state | Uint64(0xA529));
										//return ((z << 43) | (z >> 21)) ^ (((z << 17) | (z >> 47)) * Uint64(0x6C8E9CF570932BD5));

										// also works at 32TB
										//state += stream;
										//Uint64 z = (state ^ state >> 25) * (state | Uint64(0xA529));
										//return z ^ ((z << 42) | (z >> 22)) ^ ((z << 25) | (z >> 39));


										//older
										//Uint64 z = (state ^ state >> 25) * ((state & Uint64(0xFFFFFFF8)) ^ Uint64(0x9E3779B97F4A7C15));

				}
				std::string vortex64::get_name() const {
					std::ostringstream tmp;
					tmp << "vortex64x" << (int)((stream ^ Uint64(0x6C8E9CF570932BD5)) >> 32);
					return tmp.str();
				}
				void vortex64::walk_state(StateWalkingObject *walker) {
					walker->handle(state);
				}


				Uint64 molerat64::raw64() {
					//surprisingly, passes 32TB with one anomaly, seed = 0xde9b705d
					// passes TestU01 in forwards and reverse as well
					//return (a = rotate64(a, 29) * UINT64_C(0xAC564B05)) * UINT64_C(0x818102004182A025);
					return (a += 0x589965CC75374CC3ULL - rotate64(a, 37)) * 0xE7037ED1A0B428DBULL;
					//a = 0x9E3779B97F4A7AF5ULL - rotate64(a, 41);
					//return (a ^ a >> 26) * 0xE7037ED1A0B428DBULL;
					

					//has serious issues at 16TB
					//return (a = rotate64(a, 29) * 0x10005ULL) * 0x818102004182A025ULL;
					//const uint64_t b = a;
					//return rotate64(b, 11) + (a = rotate64(b, 21) * UINT64_C(0x9E3779B9));
					// 0xAEF17502108EF2D9ULL; 0x9E3779B97F4A7AF5ULL 0xE7037ED1A0B428DBULL 0x8EBC6AF09C88C6E3ULL 0xA0761D6478BD642FULL 0x589965CC75374CC3ULL
					//return (z << 29) + rotate64(z, 20);
					//z = (z ^ z >> 28);
					//return z;
				}
				std::string molerat64::get_name() const { return "molerat64"; }
				void molerat64::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					uint64_t r = a;

					//a = UINT64_C(0x9B1B51BEB2EFF7A1); //0x41C64E6B
					//for (uint64_t ra = (r & 0xFFFF); ra; ra--)
					//{
					//	a = rotate64(a, 20) + (a << 27);
					//}
					a &= 0x1fffff; //0x41C64E6B
					++a;
					for (uint64_t ra = (r & 0x7FF); ra; ra--)
					{
						//a = rotate64(a, 29) * 0xAC564B05ULL;
						a += 0x589965CC75374CC3ULL - rotate64(a, 37);
					}

				}

				Uint64 moverCounter64::raw64() {
					//passes 32TB with no anomalies, seed = 0x8b2bfcb
					//return (a = rotate64(a, 21) * UINT64_C(0x9E3779B9)) * (b += UINT64_C(0x9E3779B97F4A7AF6));// *((b += UINT64_C(0x9E3779B97F4A7AF5)) | UINT64_C(1));
					// 0xAEF17502108EF2D9; 0x9E3779B97F4A7AF5 0xC6BC279692B5CC83 0x6C8E9CF570932BD5 0xDB4F0B9175AE2165 0x369DEA0F31A53F85
					//return (b += (a = rotate64(a, 21) * UINT64_C(0x9E3779B9)) ^ UINT64_C(0x9E3779B97F4A7AF5));//* UINT64_C(0x41C64E6D)); //
					//return (a = (rotate64(a, 21) ^ (b += UINT64_C(0x9E3779B97F4A7AF5))) * UINT64_C(0xC6BC279692B5CC83)) * UINT64_C(0x41C64E6B);
					
					//return (a = rotate64(a, 21) * (b += UINT64_C(0x9E3779B97F4A7AF6))) * UINT64_C(0x41C64E6B);

					//passes 32TB one anomaly, TestU01 BigCrush both ways
					//return (a = rotate64(a, 29) * UINT64_C(0xAC564B05) ^ ++b) * UINT64_C(0x818102004182A025);
					////passes at least 8TB and BigCrush both ways
					//return (a = rotate64(a, 35) ^ (b = b * UINT64_C(0x369DEA0F31A53F85) + UINT64_C(0x9E3779B97F4A7AF5))) * UINT64_C(0xDB4F0B9175AE2165);
					
					//return (a = rotate64(a, 29) ^ (b = b * UINT64_C(0xD1B54A32D192ED03) ^ UINT64_C(0xDB4F0B9175AE2165))) * UINT64_C(0x2545F4914F6CDD1D);
					
					//const uint64_t a0 = a * 127ULL;
					//return (a = (rotate64(a0, 35) ^ (b += UINT64_C(0xDB4F0B9175AE2165)))) + a0;

					//const uint64_t a0 = a - (a << 19);
					//return (a += (rotate64(a0, 35) ^ (b += 0xDB4F0B9175AE2165ULL)));

					//return (a = (rotate64(a, 35) ^ (b += 0xDB4F0B9175AE2165ULL))) + (c = 0x9E3779B97F4A7AF5ULL - rotate64(c, 11));
					//return (a = (rotate64(a, 35) ^ (b += 0xDB4F0B9175AE2165ULL))) + (c = (c >> 1 ^ (-(c & 1u) & 0xD800000000000000ULL)));
					//// I think this one does well?
					//return (a = (rotate64(a, 21) + (c = rotate64(c, 35) ^ (b += 0x9E3779B97F4A7AF5ULL))));
					
					//const uint64_t result = ~(a << 1) * (b = rotate64(b, 11) + 0x9E3779B97F4A7AF5ULL);
					//return (result ^ result >> 29) + (a = rotate64(a, 21) + 0xDB4F0B9175AE2165ULL);

					//// passes 32TB with one anomaly
					//return (a = rotate64(a, 29) * UINT64_C(0x9E3779B9) ^ ++b) * UINT64_C(0x8181020042010415);
					//UINT64_C(0x6C8E9CF570932BD5) UINT64_C(0xDB4F0B9175AE2165)
					////passes until 16TB, then mildly suspicious
					//return (a = rotate64(a, 41) ^ (b += UINT64_C(0x6C8E9CF570932BD5))) * UINT64_C(0xDB4F0B9175AE2165);
					//return (a = rotate64(a, 29) ^ (b = b * UINT64_C(0xD1B54A32D192ED03) ^ UINT64_C(0xDB4F0B9175AE2165)))* UINT64_C(0x2545F4914F6CDD1D);
					////passes 32TB with one anomaly
					//const uint64_t a0 = a * UINT64_C(0x2545F4914F6CDD1D);
					//return (a = rotate64(a0, 35) ^ (b += UINT64_C(0xDB4F0B9175AE2165))) - a0;

					//const uint64_t a0 = a * UINT64_C(0xAC564B05);

					//passes 2TB at least with one anomaly at 512GB
					//const uint64_t a0 = a + b + UINT64_C(0xD1B54A32D192ED03);
					//return (a = rotate64(a0, 19) ^ (b += UINT64_C(0xDB4F0B9175AE2165))) * UINT64_C(0xDE4D);
					// passes great up to 8TB, one anomaly at 16TB and 3 mildly suspicious + 1 unusual at 32TB
					//const uint64_t a0 = a + b;
					//return (a = rotate64(a0, 19) ^ (b += UINT64_C(0xDB4F0B9175AE2165))) * UINT64_C(0x8A35);
					
					//return (a = (rotate64(a0, 47) ^ (b += 0x9E3779B97F4A7AF5ULL)));
					
					//b += a;
					//return (a = rotate64(a, 29) * UINT64_C(0xAC564B05)) * (b | 1);
					
					//return (b += (a = rotate64(a, 47) * UINT64_C(0x818102004182A025)) ^ UINT64_C(0x9E3779B97F4A7AF5));
					///return (b += (a = rotate64(a, 29) * UINT64_C(0xAC564B05)));
					
					//final long ab = a + b; return (a = (ab << 19 | ab >>> 45) ^ (b += 0xDB4F0B9175AE2165L)) * 0x8A35L;

					//0xAC564B05UL 0x818102004182A025UL
					//return (a = rotate64(a, 28) ^ (b += 0xD1B54A32D192ED03ULL)) * 0x818102004182A025UL;

					//PartyRandom
					//Passes 64TB with no anomalies.
					//Period is a non-zero multiple of (2 to the 64), minimum guarantee of (2 to the 64).
					//Given all possible subcycles, all results are equally frequent, but any given subcycle may have a different frequency for some results.
					//Uses only ARX operations. The output step is based on 2 rounds of the Speck cipher with key 0,0, and different rotations.
					//uint64_t t = rotate64(a, 37) ^ (b += 0xD1B54A32D192ED03ULL), s = rotate64(t, 47) + b;
					//a = s ^ rotate64(b, 23);
					//return (rotate64(s, 47) + a) ^ rotate64(a, 23);

					// fails BRank immediately.
					return (a = a * 0xD1342543DE82EF95ULL + 0xDE916ABCC965815BULL ^ (b = (b << 1) ^ (0u - (b >> 63) & 0xFEEDBABEDEADBEEFL))) ^ a >> 27;
				}
				std::string moverCounter64::get_name() const { return "moverCounter64"; }
				void moverCounter64::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					b += (b == 0);
					//b |= UINT64_C(1);
					printf("Seed is 0x%016llX, 0x%016llX\r\n", a, b);
					//b = b << 3 | UINT64_C(5);

					//uint64_t r = a;

					//////a = UINT64_C(0x9B1B51BEB2EFF7A1); //0x41C64E6B
					//////for (uint64_t ra = (r & 0xFFFF); ra; ra--)
					//////{
					//////	a = rotate64(a, 20) + (a << 27);
					//////}
					//a = UINT64_C(0x9E3779B9); //0x41C64E6B 0x9E3779B9
					//for (uint64_t ra = (r & 0xFFFF); ra; ra--)
					//{
					//	a = rotate64(a, 21) * UINT64_C(0x9E3779B9);
					//}

				}
				Uint32 moverCounter32::raw32() {
					//passes 32TB with no anomalies, seed = 0x8b2bfcb
					//return (a = rotate64(a, 21) * UINT64_C(0x9E3779B9)) * (b += UINT64_C(0x9E3779B97F4A7AF6));// *((b += UINT64_C(0x9E3779B97F4A7AF5)) | UINT64_C(1));
					// 0xAEF17502108EF2D9; 0x9E3779B97F4A7AF5 0xC6BC279692B5CC83 0x6C8E9CF570932BD5
					//return (b += (a = rotate64(a, 21) * UINT64_C(0x9E3779B9)) ^ UINT64_C(0x9E3779B97F4A7AF5));//* UINT64_C(0x41C64E6D)); //
					//return (a = (rotate64(a, 21) ^ (b += UINT64_C(0x9E3779B97F4A7AF5))) * UINT64_C(0xC6BC279692B5CC83)) * UINT64_C(0x41C64E6B);
					
					//return (a ^= (rotate32(a, 11) ^ (a << 13)) * (b += UINT32_C(0x9E3779BA)) + UINT32_C(0x6C8E9CF5)) * UINT32_C(0x41C64E6B);
					//return ((a = rotate32(a, 17) * UINT32_C(0xBCFD) ^ b) * (b += UINT32_C(0x9E3779BA)));
					//const uint32_t b0 = b;
					//return (a = rotate32(a, 19) + (b0 ^ b0 >> 11)) * (b0 ^ (b += 0x6C8E9CF5u));
					////left off here
					//return ((a = rotate32(a, 17) * UINT32_C(0xBCFD) ^ b) * (b += UINT32_C(0x9E3779BA)));
					//return ((a = rotate32(a, 12) * (b *= UINT32_C(0x12D32D)) + UINT64_C(0x9E3779B9)) * UINT32_C(0xB1AD3));
					////passes 32TB with 2 anomalies, unusual, worsening at the end 
					//const uint64_t a0 = a * UINT32_C(0xA529D), b0 = b;
					//return (a = rotate32(a0, 19) ^ (b += UINT32_C(0x9E3779BD))) + a0 ^ b0;
					//b = rotate32(b, 13) * UINT32_C(0x89A7); // period is a multiple of 16
					//return ((a = rotate32(a, 17) * UINT32_C(0xBCFD))) ^ (b = rotate32(b, 28) * UINT32_C(0xA01B));// * UINT32_C(0xA529D);
					//return ((a = rotate32(a, 17) * UINT32_C(0xBCFD))) ^ (b = rotate32(b, 7) + UINT32_C(0xC0EF50EB));// * UINT32_C(0xA529D);
					//// passes 32TB with one early anomaly, [Low8/32]FPF-14+6/16:(6,14-0), unusual, at 32GB.
					//const uint32_t result = ((a = rotate32(a, 1) + 0x929E7143u) >> 11u | 1u) * (b = rotate32(b, 25) + 0xC4DE9951u);
					//return (result ^ result >> 17) + a;
					// no anomalies up to 8TB; other versions have gone to "suspicious" (on the LSB only) at 32TB, but it isn't especially fast.
					//const uint32_t result = (b = rotate32(b, 25) + 0xC4DE9951u) * (a >> 11 | 1u);
					//return result ^ (result >> 16) + (a = rotate32(a, 1) + 0xAA78EDD7u);
					//// passes 32TB no anomalies
					//const uint32_t b1 = (b += 0x9E3779BDu);
					//return (a = rotate32(a, 21) * (b1 | 0xFFE00001u)) * 0xA5295u ^ b1;
					
					////passes 32TB no anomalies
					//const uint32_t a1 = rotate32(a, 1) + 0xAA78EDD7u;
					//const uint32_t b1 = rotate32(b, 25) + 0xC4DE9951u;
					//const uint32_t r = a1 ^ b1;
					//a = b1 ^ rotate32(b1, 13) ^ rotate32(b1, 19);
					//b ^= a1 + a;
					//return r;
					
					//const uint32_t b1 = (b += 0x9E3779BDu);
					//const uint32_t c = a + 0xC4DE995Au;
					return (a += 0x9E3779BDu + (b = rotate32(b, 17) * 0xBCFDu ^ rotate32(a, 11)));
					


					//(a = rotate32(a, 1) + 0xAA78EDD7u);
					
					//uint32_t r = (a += 0xAA78EDD7u) + (b = rotate32(b, 25) + 0xC4DE9951u);
					//r ^= r >> 7;
					//r ^= r << 1;
					//return r ^ r >> 9;
					
					//return r ^ rotate32(r, 13) ^ rotate32(r, 19);
					//return r;
				}
				std::string moverCounter32::get_name() const { return "moverCounter32"; }
				void moverCounter32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					walker->handle(b);
					//b = b << 3 | UINT64_C(5);

					//uint64_t r = a;
					//////a = UINT64_C(0x9B1B51BEB2EFF7A1); //0x41C64E6B
					//////for (uint64_t ra = (r & 0xFFFF); ra; ra--)
					//////{
					//////	a = rotate64(a, 20) + (a << 27);
					//////}
					//a = UINT32_C(0x89A7); //0x41C64E6B 0x9E3779B9
					//for (uint32_t ra = (r & 0xFFFF); ra; ra--)
					//{
					//	a = rotate32(a, 17) * UINT32_C(0xBCFD);
					//}
					//b = UINT32_C(0x3F10AF16);
					//for (uint32_t rb = (r >> 16); rb; rb--)
					//{
					//	b = rotate32(b, 25) + UINT32_C(0xC0EF50EB);
					//}

				}

				Uint64 xrsr_rev_mul::raw64() {
					const Uint64 s0 = state0;
					Uint64 s1 = state1;
					const Uint64 result = s1 * UINT64_C(0x369DEA0F31A53F85);
					s1 ^= s0;
					state0 = rotate64(s0, 24) ^ s1 ^ (s1 << 16); // a, b
					state1 = rotate64(s1, 37); // c
					return (result ^ result >> 25) + s0;
				}
				std::string xrsr_rev_mul::get_name() const { return "xrsr_rev_mul"; }
				void xrsr_rev_mul::walk_state(StateWalkingObject *walker) {
					walker->handle(state0);
					walker->handle(state1);
				}
				Uint32 twirl32::raw32() {

					//a += 0x62BD5;
					//Uint32 z = (a ^ a >> 13) * ((a & 0xFFFF8) ^ 0xCD7B5);// +0x64B05;
					//return (((z << 7) | (z >> 25)) - a) ^ (z >> 13);

					//const Uint32 INCR = 0x70932BD5, ADJ = 0x63C5, MUL = 0x2C9277B5, PLUS = 0xAC564B05;
					//enum { SHIFT_A = 13, SHIFT_B = 13, ROT = 7 };

					//return rotate32(z, ROT) ^ (z >> SHIFT_B);
					//return (z ^ z >> SHIFT_B);

					//enum {
					//	X = 0xC74EAD55,//0x632BE5AD,//must end in 5 or D
					//	Y = 0x632BE5AD,
					//	M = 0x93DB3,//must end in 3 or B
					//	N = 0xDAD7B,
					//};
					enum {
						X = 0xC74EAD55,//0x632BE5AB,//must end in 5 or D
						M = 0xD3AB3,//must end in 3 or B
					};
					Uint32 z = (a = (a ^ X) * M);
					z ^= z >> 16;
					z ^= z << 7;
					z = (z ^ z >> 11) * Uint32(0xEF52D);
					return z ^ z >> 13;
					/*
					Uint32 z = (a += Uint32(0x6A52D));
					z = (z ^ ((z << 11) | (z >> 21)) ^ ((z << 28) | (z >> 4)) ^ ((z << 5) | (z >> 27)) ^ ((z << 19) | (z >> 13))) * 0xAD7B5;
					z ^= z >> 15;
					z ^= z << 13;
					z ^= z >> 14;
					z = z * Y + N;
					z = (z ^ ((z << 12) | (z >> 20)) ^ ((z << 17) | (z >> 15)) ^ ((z << 22) | (z >> 10)) ^ ((z << 7) | (z >> 25))) * 0x62BD5;
					z ^= z >> 3;
					z ^= z << 23;
					z ^= z >> 25;
					z = (z ^ X) * M;
					return z ^ z >> 16;
					*/
					/* // left off here, fails tests early
					Uint32 z = a += Uint32(0x62BD5);
					z ^= ((z << 11) | (z >> 21)) ^ ((z << 18) | (z >> 14));
					//z = z * Uint32(0xBE5AD) + Uint32(0x6A529);
					z = (z ^ z >> 14) * Uint32(0xAD7B5);
					z = (z ^ z >> 15) * Uint32(0xCDD1D);
					z ^= ((z << 9) | (z >> 23)) ^ ((z << 19) | (z >> 13));
					z = (z ^ z >> 13) * Uint32(0x74C6D);
					return z ^ (z >> 15);
					*/
					//z ^= ((z << 13) | (z >> 19)) ^ ((z << 7) | (z >> 25));
					//z ^= z >> 11;
					//z ^= z << 21;
					//z = (z ^ z >> 13) * 0xCD7B5;
					//return z ^ z >> 13;

					//					z ^= ((z << 21) | (z >> 11)) ^ ((z << 15) | (z >> 17));


					//state += stream; //0x6C8E9CF570932BD5 0x6C8E9CF970932BD5
					//Uint64 z = (state ^ state >> 25) * UINT64_C(0x2545F4914F6CDD1D);
					//z ^= ((z << 19) | (z >> 45)) ^ ((z << 53) | (z >> 11));
					//return z ^ z >> 25;


					//Uint32 state = ++a * 0x62BD5 ^ 0x9E3779B9;
					//Uint32 z = (state ^ state >> 13) * ((state & 0xFFFF8) ^ 0xCD7B5) ^ 0x632BE5AB;
					//return ((z << 21) | (z >> 11)) ^ (((z << 7) | (z >> 25)) * 0x62BD5);

					//return Integer.rotateLeft((state = ((state = state * 0x62BD5 ^ 0x9E3779B9) ^ state >>> 13) * ((state & 0xFFFF8) ^ 0xCD7B5) ^ 0x632BE5AB), 21) ^ (Integer.rotateLeft(state, 7) * 0x62BD5 + 0x932BD);
				}
				std::string twirl32::get_name() const { return "twirl32"; }
				void twirl32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					//if (!a) a = 1;
				}
				Uint32 zig32::raw32() {
					/*
					enum {
					X = 0xC74EAD55,//must end in 5 or D
					Y = 0x70932BD5,
					M = 0x947E3DB3,//0x94DB3,//must end in 3 or B
					N = 0xDB3947E3,//0xC564B,
					};
					a = (a ^ X) * M;
					b = (b ^ Y) * N;
					*/
					enum {
						X = 0xC74EAD55,//0x632BE5AB,//must end in 5 or D
						M = 0xA5CB3,//must end in 3 or B
					};
					// doesn't work well at all
					//a ^= a >> 13;
					//a ^= a << 5;
					//a ^= a >> 19;
					//return (a ^ a >> 10) + (b *= 0xE5CB5);

					//LCG style
					//Uint32 z = (b = b * 0xA5CB5 + 0xC74EAD53);
					//XLCG style
					//Uint32 z = (b = (b ^ 0xC74EAD55) * 0xA5CB3);

					//// Zog, passes 32TB with 1 "unusual"
					//Uint32 z = (b += 0xC74EAD55);
					//a ^= a >> 14;
					//z = (z ^ z >> 10) * 0xA5CB3;
					//a ^= a >> 15;
					//return (z ^ z >> 20) + (a ^= a << 13);


					// passes 32TB with 2 "unusual"
					//Uint32 z = (b = (b ^ 0xC74EAD55) * 0xA5CB3);
					//a ^= a >> 14;
					//z ^= z >> 10;
					//a ^= a >> 15;
					//return (z ^ z >> 20) + (a ^= a << 13);

					// passes 32TB with 3 "unusual"
					// X is 0xC74EAD55
					// M is 0xA5CB3
					/*
					Uint32 z = (b = (b ^ X) * M);
					a ^= a >> 14;
					z ^= z >> 13;
					a ^= a >> 15;
					z = (z ^ z >> 11) + (a ^= a << 13);
					return (z ^ z >> 7);
					*/

					/*
					Uint32 z = (b = (b ^ X) * M);
					a ^= a >> 14;
					z ^= z >> 16;
					a ^= a << 1;
					z ^= z << 7;
					return (z ^ z >> 11) * 0x94DB5 + (a ^= a >> 15);
					*/

					//z = (z ^ z << 7) + (a ^= a >> 15);
					//return z ^ z >> 13; // does not work well with 14
					//z += (z << 16) + (a ^= a >> 15); // unsure if this works
					//z = (z ^ z >> 12) + (a ^= a >> 15); // 64GB pass
					//z = (z ^ z >> 13) + (a ^= a >> 15); // has worked with z >> 12

					// where we left off, passes at least 4 TB
					/*
					Uint32 z = (b += 0x632BE5AB);
					a ^= a >> 14;
					z = (z ^ z >> 15) * Uint32(0xAD7B5);
					a ^= a << 1;
					z = (z ^ z >> 12) * Uint32(0xCDD1D) + (a ^= a >> 15);
					return z ^ z >> 13;
					*/
					/*
					enum {
					X = 0xC74EAD55,//0x632BE5AB,//must end in 5 or D
					M = 0x93DB3,//must end in 3 or B
					};
					Uint32 z = (b = (b ^ X) * M);
					a ^= a >> 14;
					z ^= z >> 15;
					a ^= a << 1;
					z += (a ^= a >> 15) ^ (z ^ z >> 12) * UINT32_C(0xCDD1D);
					return z ^ z >> 13;
					*/

					/* // passes 32TB with only one anomaly, "unusual" on [Low1/32]FPF-14+6/16:all2 at 1TB
					Uint32 z = (b += 0x632BE5AB);
					a ^= a >> 14;
					z = (z ^ z >> 14) * Uint32(0xAD7B5);
					a ^= a << 1; // 13 maybe?
					z = (z ^ z >> 15) * Uint32(0xCDD1D);
					z += (a ^= a >> 15);
					return z ^ z >> 13; // maybe 16 ?

					*/

					//a += 0x62BD5;
					//Uint32 z = (a ^ a >> 13) * ((a & 0xFFFF8) ^ 0xCD7B5);// +0x64B05;
					//return (((z << 7) | (z >> 25)) - a) ^ (z >> 13);

					//const Uint32 INCR = 0x70932BD5, ADJ = 0x63C5, MUL = 0x2C9277B5, PLUS = 0xAC564B05;
					//enum { SHIFT_A = 13, SHIFT_B = 13, ROT = 7 };

					//return rotate32(z, ROT) ^ (z >> SHIFT_B);
					//return (z ^ z >> SHIFT_B);


					// works, passes 32TB of tests
					/*
					a ^= a << 7;
					a ^= a >> 1;
					a ^= a << 9;
					b += 0x632BE5AB;
					Uint32 z = (a + b ^ a - b >> 14) * ((b & 0xFFFF8) ^ 0xCD7B5);
					return z ^ z >> 13 ^ a;
					*/
					//b = b * ((a & 0xFFFF8) ^ 0xCD7B5) + (a & 0x62BD5);
					//return b ^ (b - a) >> 13 ^ a;

					// where we left off
					/*
					a ^= a << 1;
					a ^= a >> 27;
					b += 0x632BE5AB;
					a ^= a << 27;
					Uint32 z = (b + a ^ b - a >> 14) * ((a & 0xFFFFE) ^ 0xCD7B5);
					return z ^ z >> 13;
					*/

					//return (z + b ^ a >> 13) ^ ((z + a ^ b >> 11) ^ (b + a ^ z >> 12));
					//Uint32 z = (b ^ b - a >> 13);
					//return ((z << 21) | (z >> 11)) ^ ((z << 7) | (z >> 25)) * ((a & 0xFFFF8) ^ 0xCD7B5);

					//return (a ^ b >> 10) ^ (b - a >> 13) ^ a * (0xA529 | (b >> 12));
					//Uint32 z = (a ^ a >> 13) * 0x62BD5;
					//return (a ^ z) + ((a << 7) | (a >> 25)) ^ ((z << 22) | (z >> 10));

					//Uint32 state = ++a * 0x62BD5 ^ 0x9E3779B9;
					//Uint32 z = (state ^ state >> 13) * ((state & 0xFFFF8) ^ 0xCD7B5) ^ 0x632BE5AB;
					//return ((z << 21) | (z >> 11)) ^ (((z << 7) | (z >> 25)) * 0x62BD5);

					//return Integer.rotateLeft((state = ((state = state * 0x62BD5 ^ 0x9E3779B9) ^ state >>> 13) * ((state & 0xFFFF8) ^ 0xCD7B5) ^ 0x632BE5AB), 21) ^ (Integer.rotateLeft(state, 7) * 0x62BD5 + 0x932BD);

					//const uint32_t y = (a += 0x9E3779B9);
					//const uint32_t z = ((b += (-(a > 0x52345670) & 0x632BE5AB)) ^ a) * 0xABCFD;
					//return (z ^ rotate32(z, 21) ^ rotate32(z, 9)) * (a >> 11 | 1u) ^ a;
					//const uint32_t s = (a += 0x632BE5ABU);
					//const uint64_t z = (s ^ s >> 13 ^ s >> 20) * (b >> 12 | 0x1U);
					//if (s != 0U)
					//    b += 0x9E3779BDU;
					//return z ^ z >> 21 ^ z >> 11;


					//// passes 32TB with no anomalies
                    //uint32_t s = (a += 0xC1C64E6Du);
                    //uint32_t t = (s == 0u) ? b : (b += 0x9E3779BDu);// (b += -((s | -s) >> 31) & 0x9E3779BBu);
                    //uint32_t x = (s ^ s >> 17) * ~((t ^ t >> 12) & 0x1FFFFEu);
                    //x = (x ^ x >> 16) * 0xAC451u;
                    //return (x ^ x >> 15);
					

					//(a = a * 0x7E57Du + 1u);//(b = b * 0x6EBD5u - 1u)
					//const uint32_t z = (((s + 0x80000000u < 0x65C5u) ? b : (b += 0xD8AB2475u)) | 1u) * (s ^ s >> 15);

					const uint32_t s = (a += 0xC1C64E6Du);
				    uint32_t z = ((s + 0x80000000u < 0x5C5u) ? b : (b += 0xD8AB2475u)) ^ s ^ s >> 15;
					z ^= rotate32(z, (s >> 28) + 16) ^ rotate32(z, (s >> 24 & 15) + 1);
					//const uint32_t s = (a = a * 0xFB85u + 1u);
				    //uint32_t z = ((s == 0u) ? b : (b = b * 0xD09Du + 1u)) ^ s ^ s >> 14;
					//z ^= z >> 15;

					//z ^= z >> (s >> 28) + 4;
					//s ^= s >> (z >> 28) + 6;
					//z ^= z << (s >> 24 & 15) + 4;
					//s ^= s << (z >> 24 & 15) + 1;
					//z ^= z >> (s >> 20 & 15) + 4;
					//s ^= s >> (z >> 20 & 15) + 11;
					//return z;
					//s ^= rotate32(s, (z >> 28) + 4) ^ rotate32(s, (z >> 24 & 15) + 16);
//					s ^= rotate32(s, (z >> 20 & 15) + 8) ^ rotate32(s, (z >> 16 & 15) + 12);
//					z ^= rotate32(z, (s >> 20 & 15) + 8) ^ rotate32(z, (s >> 16 & 15) + 12);
// rotate32(z, (s >> 20 & 15) + 11) ^ rotate32(z, (s >> 16 & 15) + 6)
					//return z ^ rotate32(z, 11) ^ rotate32(z, 23) ^ rotate32(z, 7) ^ rotate32(z, 19);

					//return rotate32(z, (s >> 29)) ^ rotate32(z, (s >> 23 & 7) + 16) ^ rotate32(z, (s >> 17 & 7) + 8);
					return z ^ rotate32(z, (s >> 20 & 15) + 16) ^ rotate32(z, (s >> 16 & 15) + 1);
					
					
					//uint64_t r = rotate32(a, 16) ^ b;
					//r = rotate32(r, a >> 27);
					//r ^= r >> 13;
					//r ^= r << 17;
					//r ^= r >> 15;
					//a = a * 0xFB83u ^ r;
					//b = b * 0xD09Bu ^ 0xC1C64E6Du;
					////a = a * 0xFB85u + r;
					////b = b * 0xD09Du + 1u;
					//return r;


					//uint32_t x, y, z, w;
//x = (stateA += 0xDB4F0B91);
//y = (stateB += (rotate32(x, 21) + __lzcnt(x)));
//z = (stateC += (rotate32(y, 21) + __lzcnt(x &= y)));
//w = (stateD += (rotate32(z, 21) + __lzcnt(x &= z)));
//x += rotate32(w, 21);
//x ^= x >> 15;
//x *= 0x2c1b3c6d;
//x ^= x >> 12;
//x *= 0x297a2d39;
//x ^= x >> 15;
//return x;

				}
				std::string zig32::get_name() const { return "zig32"; }
				void zig32::walk_state(StateWalkingObject *walker) {
					walker->handle(a);
					//if (a == 0) a = 1;
					walker->handle(b);
				}


				Uint64 twinLinear::raw64() {
//					uint64_t r = rotate64(s0, 32) ^ s1;
//					r = rotate64(r, s0 >> 58);
//					//r += 0x2545F4914F6CDD1DULL;
//					s0 = s0 * 0xFF2826ADULL + r;
//					s1 = s1 * 0xFF1CD035ULL + 1ULL;
////					s0 = s0 * 0xB67A49A5466DULL + r;
////					s1 = s1 * 0x87338161EF95ULL + 1ULL;
////					s0 = s0 * 0x2C6FE96EE78B6955ULL + r;
////					s1 = s1 * 0x369DEA0F31A53F85ULL + 1ULL;
//					return r ^ r >> 32;

					// TwinsiesRandom
					// Passes 256TB with no anomalies.
					// Period is 2 to the 64, with 2 to the 64 streams.
					// Stepping backwards is equivalent to a subtraction and a multiplication by an inverse, for each state.
					// Identifying the stream would probably require a variable-length jump through an LCG.
					//uint64_t x = rotate64(s0, 33) + s1;
					//s0 = s0 * 0x369DEA0F31A53F85ULL + 0x2C6FE96EE78B6955ULL;
					//s1 = s1 * 0xD1342543DE82EF95ULL + 0x9E3779B97F4A7C15ULL;
					//return x ^ x >> (x >> 59) + 6 ^ x >> 44;


					// Passes at least 2TB, but a restart was required before it could finish.
//					uint64_t x = (rotate64(s0, 33) ^ s1) * 0xF1357AEA2E62A9C5ULL;
//					s0 += 0x369DEA0F31A53F85ULL;
//					s1 += 0x9E3779B97F4A7C15ULL;
//					return x ^ x >> (x >> 59) + 6 ^ x >> 44;

					// Passes at least 64TB, but fails InitialCorrelationEvaluator test in Juniper.
					//uint64_t x = (rotate64(s0, 33) ^ s1) * 0xF1357AEA2E62A9C5ULL;
					//s0 += 0x369DEA0F31A53F85ULL;
					//s1 += 0x9E3779B97F4A7C15ULL;
					//return x ^ x >> 19 ^ x >> 44;

					// Fails many TMFn tests at 512GB.
					//uint64_t x = rotate64(s0, 33) ^ s1;
					//s0 = s0 * 0x369DEA0F31A53F85ULL + 0x2C6FE96EE78B6955ULL;
					//s1 = s1 * 0xD1342543DE82EF95ULL + 0x9E3779B97F4A7C15ULL;
					//return x ^ x >> 19 ^ x >> 44;

					// Has an anomaly ([Low1/16]Gap-16:A, "unusual") at 64TB...
					//uint64_t x = (rotate64(s0, 33) ^ s1) * 0xF1357AEA2E62A9C5ULL;
					//s0 += 0x369DEA0F31A53F85ULL;
					//s1 += 0x9E3779B97F4A7C15ULL;
					//x ^= x >> (x >> 59) + 9 ^ x >> 47;
					//return x ^ x << (x & 31u) + 6 ^ x << 43;

					// Passes at least 4TB with no anomalies.
					//uint64_t x = (rotate64(s0, 33) ^ s1) * 0xF1357AEA2E62A9C5ULL;
					//s0 += 0x369DEA0F31A53F85ULL;
					//s1 += 0x9E3779B97F4A7C15ULL;
					//x ^= x >> (x >> 59) + 7 ^ x >> 47;
					//return x ^ x << (x & 31u) + 6 ^ x << 41;

					// Passes at least 128GB here, but fails ICE test in Juniper.
					//uint64_t x = (rotate64(s0, 33) ^ s1);
					//s0 += 0x369DEA0F31A53F85ULL;
					//s1 += 0x9E3779B97F4A7C15ULL;
					//x ^= x >> (x >> 59) + 7;
					//x *= 0xF1357AEA2E62A9C5ULL;
					//return x ^ x >> (x >> 59) + 6;

					// Fails fairly quickly, at 64GB.
					//uint64_t x = (rotate64(s0, 21) ^ s1);
					//s0 += 0x369DEA0F31A53F85ULL;
					//s1 += 0x9E3779B97F4A7C15ULL;
					//x ^= x >> 25 ^ x >> 50;
					//x *= 0xF1357AEA2E62A9C5ULL;
					//return x ^ x >> (x >> 59) + 6;

					//fails [Low8/32]Gap-16:A and B early...
					//uint64_t x = s0 ^ s0 >> 31;
					//uint64_t y = s1 ^ s1 >> 29;
					//s0 += 0x369DEA0F31A53F85ULL;
					//s1 += 0x9E3779B97F4A7C15ULL;
					//x ^= rotate64(y, 30);
					//x *= 0xF1357AEA2E62A9C5ULL;
					//return x ^ x >> (x >> 59) + 6;

					// fails ICE test when using rotate64(s0, 33), but does pass 64TB of PractRand
					//uint64_t x = (rotate64(s0, 21) ^ s1);
					//s0 += 0x369DEA0F31A53F85ULL;
					//s1 += 0x9E3779B97F4A7C15ULL;
					//x ^= x >> (x >> 60) + 14;
					//x *= 0xF1357AEA2E62A9C5ULL;
					//return x ^ x >> (x >> 59) + 6 ^ x >> 44;

					// This passes 64TB with no anomalies!
					// Period is exactly 2 to the 128, and all states are valid.
					//uint64_t x = rotate64(s0, 33) + s1;
					//s1 = s1 * 0xD1342543DE82EF95ULL + __lzcnt64(s0);
					//s0 = s0 * 0x369DEA0F31A53F85ULL + 0x2C6FE96EE78B6955ULL;
					//return x ^ x >> 26 ^ x >> 37;

					// Finishes 64 TB with no anomalies. Period is (2 to the 128) - (2 to the 64).
					//uint64_t a = s0, b = s1;
					//uint64_t z = (a + b);
					//s0 = a * 0xD1342543DE82EF95ULL + 1ULL;
					//s1 = (b << 1) ^ (0u - (b >> 63) & 0xFEEDBABEDEADBEEFULL);
					//z = (z ^ rotate64(z, 25) ^ rotate64(z, 50));// *0xBEA225F9EB34556DULL;
					//return (z ^ rotate64(z, 11) ^ rotate64(z, 42)) + b;

// Has a problem only at the very end, at 64TB!
//rng=twinLinear, seed=0x0
//length= 64 terabytes (2^46 bytes), time= 174079 seconds
//  Test Name                         Raw       Processed     Evaluation
//  [Low4/64]FPF-14+6/16:(0,14-0)     R=  +9.6  p =  1.7e-8   suspicious
//  ...and 1158 test result(s) without anomalies
					//uint64_t a = s0, b = s1;
					//uint64_t z = (a + b);
					//s0 = a + 0xD1342543DE82EF95ULL;
					//s1 = (b << 1) ^ (0u - (b >> 63) & 0xFEEDBABEDEADBEEFULL);
					//z = (z ^ rotate64(z, 25) ^ rotate64(z, 50));
					//return (z ^ rotate64(z, 11) ^ rotate64(z, 42)) + b;

// Has trouble earlier than the last one...
// But no problems before 1TB, at least.
//rng=twinLinear, seed=0x0
//length= 2 terabytes (2^41 bytes), time= 4906 seconds
//  Test Name                         Raw       Processed     Evaluation
//  [Low8/32]Gap-16:B                 R=  +7.4  p =  3.3e-6   suspicious
//  ...and 1019 test result(s) without anomalies
					//uint64_t a = s0, b = s1;
					//uint64_t z = (a ^ b);
					//s0 = a + 0xDE916ABCC965815BULL;
					//s1 = b + 0x9E3779B97F4A7C15ULL + __lzcnt64(a);
					//z = (z ^ rotate64(z, 25) ^ rotate64(z, 50)) + a;
					//return (z ^ rotate64(z, 11) ^ rotate64(z, 42)) + b;

// Somehow rotating b does worse... We didn't add a here...
//rng=twinLinear, seed=0x0
//length= 512 gigabytes (2^39 bytes), time= 1272 seconds
//  Test Name                         Raw       Processed     Evaluation
//  [Low1/64]Gap-16:A                 R=  +7.3  p =  4.0e-5   mildly suspicious
//  [Low1/64]Gap-16:B                 R= +15.9  p =  1.5e-13    FAIL
//  ...and 950 test result(s) without anomalies
//uint64_t a = s0, b = s1;
//uint64_t z = (a ^ rotate64(b, 29));
//s0 = a + 0xDE916ABCC965815BULL;
//s1 = b + 0x9E3779B97F4A7C15ULL + __lzcnt64(a);
//z = (z ^ rotate64(z, 25) ^ rotate64(z, 50));
//return (z ^ rotate64(z, 11) ^ rotate64(z, 42)) + b;

//Another failure...
//rng=twinLinear, seed=0x0
//length= 1 terabyte (2^40 bytes), time= 3629 seconds
//  Test Name                         Raw       Processed     Evaluation
//  [Low1/8]BCFN(2+19,13-7,T)         R= +32.6  p =  1.8e-10  very suspicious
//  [Low1/8]BCFN(2+20,13-8,T)         R= +58.7  p =  7.8e-16    FAIL
//  [Low1/8]BCFN(2+21,13-8,T)         R= +68.0  p =  3.4e-18    FAIL
//  ...and 985 test result(s) without anomalies
//uint64_t a = s0, b = s1, z = (b ^ rotate64(b, 11) ^ rotate64(b, 40)) + (a ^ rotate64(a, 25) ^ rotate64(a, 50));
//s0 = a + 0xDE916ABCC965815BULL;
//s1 = b + 0x9E3779B97F4A7C15ULL + __lzcnt64(a);
//return (z ^ rotate64(z, 19) ^ rotate64(z, 47));
 
					//uint64_t a = s0, b = s1, x = rotate64(a, 30) + b;
					//s0 = a * 0xD1342543DE82EF95ULL + 0x9E3779B97F4A7C15ULL;
					//s1 = b * 0xF1357AEA2E62A9C5ULL + __lzcnt64(a);
					//return x ^ x >> 28;

//x = (x ^ x >> 27 ^ y) * 0x3C79AC492BA7B653UL;
//x = (x ^ x >> 33) * 0x1C69B3F74AC4AE35UL;
//x ^= x >> 27 ^ y;

// MorbRandom
// Passes 64TB with no anomalies. Period is 2 to the 128. 1D equidistributed, exactly.
//uint64_t x = s0, y = s1;
//s0 += 0x9E3779B97F4A7C15u;
//s1 += (x + (x >> 1)) >> 63;
//x = (x ^ rotate64(x, 23) ^ rotate64(x, 47) ^ y) * 0xF1357AEA2E62A9C5ULL;
//x = (x ^ rotate64(x, 25) ^ rotate64(x, 50) ^ y);
//return x;

// MorbitalRandom
// Optimization on OrbitalRandom, with a different second multiplier.
// Passes 64TB with no anomalies. Passes ICE and IICE tests.
// Period is 2 to the 128. 1D equidistributed.
uint64_t x = s0;
uint64_t y = s1;
s0 += 0xD1B54A32D192ED03UL;
s1 += 0x8CB92BA72F3D8DD7UL + __lzcnt64(x);
y = (y ^ rotate64(x, 37)) * 0x3C79AC492BA7B653UL;
y = (y ^ y >> 33) * 0xF1357AEA2E62A9C5UL;
y ^= y >> 27;
return y;

				}
				std::string twinLinear::get_name() const { return "twinLinear"; }
				void twinLinear::walk_state(StateWalkingObject *walker) {
					walker->handle(s0);
					walker->handle(s1);
				}

				Uint64 moremur64::raw64() { // named incorrectly, may name later... quarterback64 , due to use of 25 as a rotation?
//				    const uint64_t s = (state ^ 0x9E3779B97F4A7C15u) * 0xC6BC279692B5C323u;
//					//return state += s ^ s >> 41 ^ s >> 23; 0xD1342543DE82EF23
//					return state += (s ^ rotate64(s, 23) ^ rotate64(s, 41));

					//Uint64 x = (state++);
					//x = (x ^ x >> 27) * 0x3C79AC492BA7B653UL;
					//x = (x ^ x >> 33) * 0x1C69B3F74AC4AE35UL;
					//x ^= x >> 27;
					//return x;
					// 0x9E3779B97F4A7C15UL

				    // passes 32TB with two "unusual" anomalies, both BCFN, at 2TB and 4TB.
					// mostly interesting because it is almost as good as the last one here, but
					// not quite, and the only difference is the last uses a shift by 11, while
					// this one uses 13.
                    // Note that this won't pass the full rotated/reversed/flipped battery.
//					Uint64 x = (state++);
//					x ^= x >> 27;
//					x *= 0x3C79AC492BA7B653UL;
//					x ^= x >> 33;
//					x ^= x >> 13;
//					x *= 0x1C69B3F74AC4AE35UL;
//					x ^= x >> 27;
//					return x;

					//x ^= 0x9E3779B97F4A7C15UL;

				// passes at least 32TB with one early "unusual" anomaly
				// passes 32GB for each test of 64 rotations, with or without bit flips, with or without bit-reversal.
				// gets 5 "mildly suspicious" in those 8TB of tests, nothing worse
				// slower than SplitMix64, but should allow any gamma
//					Uint64 x = (state++);
//					x ^= rotate64(x, 39) ^ rotate64(x, 14); // rotl; they correspond to right rotations of 25 and 50
//					x *= 0x3C79AC492BA7B653UL; // from Evensen's improvement on MurmurHash3's mixer, moremur64 
//					x ^= x >> 33;
//					x ^= x >> 13;
//					x *= 0x1C69B3F74AC4AE35UL; // also from moremur64
//					x ^= x >> 27; // one less shift than NASAM, seems to do well
//					return x;
//				// passes at least 32TB with no anomalies
//				// slightly slower than SplitMix64, but should allow any odd gamma
//              // Note that this won't pass the full rotated/reversed/flipped battery.
//					x ^= x >> 27;
//					x *= 0x3C79AC492BA7B653UL;
//					x ^= x >> 33;
//					x ^= x >> 11; // only change between this and Pelle Evensen's moremur64
//					x *= 0x1C69B3F74AC4AE35UL;
//					x ^= x >> 27;
//					return x;

					// wyrand
					// passes at least 16TB, but not even close to equidistributed.
					//uint64_t x = (state += 0x2D358DCCAA6C78A5ull), y = x ^ 0x8BB84B93962EACC9ull;
					//x ^= _umul128(x, y, &y);
					//return x ^ y;
					
					// I'm trying various small transformations on either an XQO sequence or an LCG.
					//uint64_t x = (state = -(state ^ (state * state | 5ull)));
					//uint64_t x = (state = state * 0x3C79AC492BA7B653ULL ^ 0xF1357AEA2E62A9C5ULL);
					//uint64_t x = (state += 0x9E3779B97F4A7C15ULL);
					//x ^= x >> (x >> 59) + 6 ^ x >> 44;
					//x *= 0xD1342543DE82EF95ULL;
					//uint64_t x = (state = state * 0xF1357AEA2E62A9C5ULL + 0x9E3779B97F4A7C15ULL);
					//uint64_t x = (state = state << 1 ^ (0ULL - (state >> 63) & 0xB5E1107E81BC107BULL));
					//uint64_t x = (state = state >> 1 ^ (0ULL - (state & 1ULL) & 0xB5E1107E81BC107BULL));

					// Passes at least 128TB of PractRand.
					// Hardware issues stopped the test earlier than expected.
					// Period is 2 to the 64, state and output sizes are each 64 bits.
					// One-dimensionally equidistributed.
					// May work as a unary hash?
					//uint64_t x = ++state;
					//x ^= x << (x & 31) + 5 ^ x << 41;
					//x *= 0xBEA225F9EB34556DULL;
					//x ^= x >> (x >> 59) + 7 ^ x >> 47;
					//x *= 0xD1342543DE82EF95ULL;
					//return x ^ x >> (x >> 59) + 6 ^ x >> 44;

					//x ^= x << 28 ^ 0xD1342543DE82EF95ULL;
					//x *= 0x3C79AC492BA7B653ULL;
					//x ^= x >> 34 ^ 0xBEA225F9EB34556DULL;
					//x *= 0xB5E1107E81BC107BULL;
					//return x ^ x >> 31;
					//return x ^ x >> (x >> 59) + 6 ^ x >> 44;
					//uint64_t x = (state = state * 0xF1357AEA2E62A9C5ULL + 0x9E3779B97F4A7C15ULL);
					//uint64_t x = (state = state + 0xF1357AEA2E62A9C5ULL ^ 0x9E3779B97F4A7C16ULL);

					//x ^= x << (x & 31) + 5 ^ x << 41;
					//x *= 0xBEA225F9EB34556DULL;
					//x ^= x >> (x >> 59) + 7 ^ x >> 47;
					//x *= 0xD1342543DE82EF95ULL;
					//return x ^ x >> (x >> 59) + 6 ^ x >> 44;

					// memorable constant, but has an "unusual" anomaly at 1TB
					//uint64_t x = (state -= 987654321987654321UL);

					// Not a failure, but close after 8TB and very close after 16TB; failure was likely imminent.
//rng=moremur64, seed=0x0
//length= 8 terabytes (2^43 bytes), time= 22566 seconds
//  Test Name                         Raw       Processed     Evaluation
//  [Low1/8]FPF-14+6/16:cross         R=  +7.1  p =  2.4e-6   mildly suspicious
//  ...and 1080 test result(s) without anomalies
// 
//rng=moremur64, seed=0x0
//length= 16 terabytes (2^44 bytes), time= 42079 seconds
//  Test Name                         Raw       Processed     Evaluation
//  [Low1/8]FPF-14+6/16:cross         R= +11.8  p =  2.8e-10   VERY SUSPICIOUS
//  ...and 1106 test result(s) without anomalies
					//uint64_t x = (state += 5555555555555555555UL);
					//x ^= x * x | 5UL; x = rotate64(x, 35); // round 1, can repeat
					//x ^= x * x | 1UL; x ^= x >> 27; // finisher is different from the round line
					//return x;

					// Sad...
//rng=moremur64, seed=0x0
//length= 4 terabytes (2^42 bytes), time= 20506 seconds
//  Test Name                         Raw       Processed     Evaluation
//  [Low4/32]BCFN(2+0,13-0,T)         R= +11.5  p =  1.1e-5   mildly suspicious
//  [Low8/32]Gap-16:B                 R=  +6.4  p =  2.0e-5   mildly suspicious
//  ...and 1050 test result(s) without anomalies
//uint64_t x = (state += 5555555555555555555UL);
//x ^= x * x | 5UL; x = rotate64(x, 42); // round 1, can repeat
//x ^= x * x | 5UL; x ^= x >> 29; // finisher is different from the round line
//return x;

// HornRandom
// Human-memorable! Passes 64TB with no anomalies! Period is 2 to the 64, 1D-equidistributed.
//uint64_t x = (state += 5555555555555555555UL);
//x ^= x * x | 7UL; x = rotate64(x, 37); // round 1, can repeat
//x ^= x * x | 7UL; x ^= x >> 27; // finisher is different from the round line
//return x;

// WoolRandom
// Also human-memorable, maybe more-so! Passes 64TB with no anomalies (on the fifth seed tried)!
//uint64_t x = (state += 5555555555555555555UL);
//x ^= x * x | 25UL; x = rotate64(x, 39); // round 1, can repeat
//x ^= x * x | 25UL; x ^= x >> 25; // finisher is different from the round line
//return x;

// fails early!
//uint64_t x = (state += 5555555555555555555UL), y = x | 1UL;
//x = x * x + y;
//x = rotate64(x, 35);
//y = x | 1UL;
//x = x * x + y;
//x ^= x >> 26;
//return x;


// Hasher.randomize1()
// Passes at least 128GB without anomalies (interrupted)
//uint64_t x = (state += 0x632BE59BD9B4E019UL);
//x = ((x ^ 0x9E3779B97F4A7C15UL) * 0xC6BC279692B5CC83UL);
//x = (x ^ x >> 27) * 0xAEF17502108EF2D9UL;
//return x ^ x >> 25;


// QomStage1
// Passes 64TB with one anomaly:
//rng=moremur64, seed=0x0
//length= 8 terabytes (2^43 bytes), time= 19660 seconds
//  Test Name                         Raw       Processed     Evaluation
//  [Low1/32]Gap-16:B                 R=  +6.0  p =  4.5e-5   unusual
//  ...and 1080 test result(s) without anomalies
// As a generator with only 64 bits of state, just one "unusual" anomaly is pretty good.
//uint64_t x = (state -= state * state | 1111111111111111111UL);
//x = (x ^ x >> 28) * 5555555555555555555UL;
//return x ^ rotate64(x, 25) ^ rotate64(x, 50);

// Has one "unusual" anomaly at 64TB, nothing before that.
//rng=moremur64, seed=0x0
//length= 64 terabytes (2^46 bytes), time= 168639 seconds
//  Test Name                         Raw       Processed     Evaluation
//  [Low4/32]FPF-14+6/16:cross        R=  -2.6  p =1-9.7e-5   unusual
//  ...and 1158 test result(s) without anomalies
uint64_t x = (state ^ state >> 28) * 5555555555555555555UL;
state -= state * state | 1111111111111111111UL;
return x ^ x >> 28;

					//x ^= x * x | 1UL; x = rotate64(x, 32); // round 1
					//x ^= x * x | 1UL; x = rotate64(x, 32); // round 2, can repeat if desired
					//x ^= x * x | 1UL; x ^= x >> 31; // last part of finisher is different from the round line
					//return x;

				}
				std::string moremur64::get_name() const { return "moremur64"; }
				void moremur64::walk_state(StateWalkingObject *walker) {
					//state = 0UL; // should only be used during the intensive many-short-tests period
					walker->handle(state);
				}

				// Surprisingly strong, considering how bad nr3q1 and nr3q2 are.
				// Passes 64TB with one anomaly, "mildly suspicious" at 128GB:
				// [Low4/64]BCFN(2+3,13-0,T)         R= -10.6  p =1-9.4e-6   mildly suspicious
				// Note that this has slightly less restrictive seeding than Troschuetz.Random, because
				// only u is affected by seeding, and all Uint64 values are valid for u. I'm not 100%
				// certain the way this seeds w is valid for that MWC generator, though.
				Uint64 nr3::raw64() {
					u = u * SeedU1 + SeedU2;
            		v ^= v >> 17;
            		v ^= v << 31;
            		v ^= v >> 8;
            		w = SeedU3 * (w & 0xFFFFFFFFUL) + (w >> 32);
            		Uint64 x = u ^ (u << 21);
            		x ^= x >> 35;
            		x ^= x << 4;
            		return (x + v) ^ w;
				}
				std::string nr3::get_name() const { return "nr3"; }
				void nr3::walk_state(StateWalkingObject *walker) {
					v = SeedV;
					w = SeedW;
					walker->handle(u);
					u ^= v;
					raw64();
					if(u == 0UL) v = SeedV;
					else v = u;
					raw64();
					w = v;
					raw64();
				}

				// Fails BRank immediately, on many tests.
				Uint64 nr3q1::raw64() {
					v ^= v >> 21;
            		v ^= v << 35;
            		v ^= v >> 4;
            		return v * SeedU;
				}
				std::string nr3q1::get_name() const { return "nr3q1"; }
				void nr3q1::walk_state(StateWalkingObject *walker) {
					walker->handle(v);
					v &= 0xFFFFFFFFUL;
					v ^= SeedV;
					v = raw64();
				}

				// Fails BRank immediately, on many tests.
				Uint64 nr3q2::raw64() {
		            v ^= v >> 17;
		            v ^= v << 31;
		            v ^= v >> 8;
		            w = SeedU * (w & 0xFFFFFFFFUL) + (w >> 32);
					return v ^ w;
				}
				std::string nr3q2::get_name() const { return "nr3q2"; }
				void nr3q2::walk_state(StateWalkingObject *walker) {
					walker->handle(v);
            		w = SeedW;
					v &= 0xFFFFFFFFUL;
            		v ^= SeedV;
            		w = raw64();
            		v = raw64();
				}

				Uint64 nova::raw64() {
				// Excellent considering this uses no multiplication.
				// Barely gets one "unusual" anomaly at 1TB, but otherwise passes 64TB:
				//   [Low1/64]Gap-16:B                 R=  +4.9  p =  4.0e-4   unusual
				// Period is (2 to the 128) - (2 to the 64). 1D equidistributed.
		            //v ^= v >> 17;
		            //v ^= v << 31;
		            //v ^= v >> 8;
					//
					//Uint64 x = (w += 0xC6BC279692B5C323UL ^ v);
					//x ^= x >> 41;
					//x += ~(x << 31);
					//return x ^ x >> 23;

				// Let's try something else.
				// Passes 64TB with no anomalies!
				// Also no multiplication, period is (2 to the 128) - (2 to the 64). 1D equidistributed.
				// Seems faster than the above version, probably because v and w (or x) are processed independently.
		            v ^= v >> 4;
		            v ^= v << 35;
		            v ^= v >> 21;

					uint64_t x = (w += 0xC6BC279692B5C323UL);
					x ^= x >> 28;
					x += ~(x << 31);
					x ^= x >> 26;
					return x ^ x >> 11 ^ v;
				}
				std::string nova::get_name() const { return "nova"; }
				void nova::walk_state(StateWalkingObject *walker) {
					walker->handle(v);
            		walker->handle(w);
					if(!v) v = 0x9E3779B97F4A7C15UL;
				}

				Uint64 mars256::raw64() {
					// Passes 64TB with no anomalies.
					// Has notably non-random behavior if it starts with an all-zero state,
					// returning each result twice successively after some initial returns of 0.
					// If on that cycle (which has unknown length), it immediately fails hwd.
					// Otherwise, it does pretty well with hwd; I don't know when it fails.
/*
					const uint64_t fa = stateA;
					const uint64_t fb = stateB;
					const uint64_t fc = stateC;
					const uint64_t fd = stateD;
					stateA = 0xD1342543DE82EF95UL * fd;
					stateB = fa + 0xC6BC279692B5C323UL;
					stateC = rotate64(fb, 47);
					stateD = fa - fc;
					return fd;
*/
//					// Passes 64TB with no anomalies, is a little faster, and does not have the
//					// non-random behavior with the all-zero initial state. It's passed 5PB
//					// in hwd without issues.
//					const uint64_t fa = stateA;
//					const uint64_t fb = stateB;
//					const uint64_t fc = stateC;
//					const uint64_t fd = stateD;
//					stateA = 0xD1342543DE82EF95UL * fd;
//					stateB = fa + 0xC6BC279692B5C323UL;
//					stateC = rotate64(fb, 47) - fd;
//					stateD = fb ^ fc;
//					return fd;

					// QuartzRandom
					// Passes 64TB with no anomalies.
					// Period is at least 2 to the 64, with 2 the 64 guaranteed streams.
					//uint64_t x = (stateA += 0xC13FA9A902A6328FUL);
					//uint64_t y = (stateB += 0x91E10DA5C79E7B1DUL);
					//uint64_t z = stateC + x;
					//uint64_t w = stateD + y;
					//stateC = z + (w ^ rotate64(w, 25) ^ rotate64(w, 43));
					//stateD = w + (z ^ rotate64(z, 47) ^ rotate64(z, 13));
					//return stateC ^ stateD;

					// JasperRandom
					// Passes 64TB with no anomalies.
					// Period is at least 2 to the 64, with 2 to the 64 guaranteed streams.
					// Unlike the above QuartzRandom, this can step backward, and doesn't have shrinking cycles.
					// Very strong results on 179PB of ReMort:
//  6	                   0  - 2.29283e-01	                   0  - 2.29283e-01	                   0  - 2.29283e-01	                   1  + 2.59071e+00	                   0  - 2.29283e-01	
//  7	                4332  - 1.56520e+00	                4318  - 2.13679e+00	                4454  + 3.42208e-01	                4324  - 1.88095e+00	                4507  + 1.91164e+00	
//  8	              361109  + 1.81290e-01	              361378  + 7.63149e-01	              360773  - 1.78373e-02	              361467  + 1.04396e+00	              361461  + 1.02365e+00	
//  9	            22358465  + 1.03355e-01	            22358210  + 7.15877e-02	            22357868  + 3.81142e-02	            22357303  + 5.73588e-03	            22357939  + 4.42027e-02	
// 11	                4415  - 3.81636e-06	                4324  - 1.88095e+00	                4426  + 2.67628e-02	                4475  + 8.11854e-01	                4426  + 2.67628e-02	
// 12	            85018781  - 1.54948e-04	            85008621  - 1.24174e+00	            85015129  - 1.66888e-01	            85014570  - 2.20096e-01	            85016354  - 7.59905e-02	
// 13	          6948709632  + 9.25650e-02	          6948764289  + 9.21463e-01	          6948620865  - 5.78564e-01	          6948794346  + 1.74373e+00	          6948649396  - 1.75030e-01	
// 14	        430511310046  + 3.58944e-02	        430511265640  + 1.48303e-02	        430512102800  + 1.95351e+00	        430511929829  + 1.28609e+00	        430511374270  + 8.25647e-02	
// 16	              362299  + 5.79253e+00	              361730  + 2.13031e+00	              361013  + 7.07403e-02	              361135  + 2.20020e-01	              359879  - 2.63021e+00	
// 17	          6948653054  - 1.40238e-01	          6948707252  + 7.60070e-02	          6948718671  + 1.70305e-01	          6948689892  + 4.54779e-03	          6948708761  + 8.63162e-02	
// 18	        567925318768  + 7.07272e+00	        567922414979  - 1.42499e+00	        567923080149  - 9.67713e-02	        567924251181  + 1.54461e+00	        567924557626  + 2.72072e+00	
// 19	      35186123978936  - 3.16554e+00	      35186126829096  - 1.68664e+00	      35186137739838  + 2.92314e-01	      35186136597463  + 1.21156e-01	      35186124690347  - 2.75316e+00	
// 21	            22356386  - 1.39718e-02	            22357853  + 3.68856e-02	            22358468  + 1.03764e-01	            22357489  + 1.32418e-02	            22358795  + 1.53101e-01	
// 22	        430512067419  + 1.80568e+00	        430511324255  + 4.45691e-02	        430511306192  + 3.37032e-02	        430512034800  + 1.67454e+00	        430512013964  + 1.59336e+00	
// 23	      35186135508102  + 2.70364e-02	      35186126775967  - 1.70998e+00	      35186126254826  - 1.94747e+00	      35186136490617  + 1.08942e-01	      35186136329128  + 9.17112e-02	
// 24	    2179984577673504  + 1.98114e-02	    2179984587147336  + 1.18103e-01	    2179984575399776  + 8.47410e-03	    2179984564436356  - 2.03795e-02	    2179984576898395  + 1.54137e-02	
// 
//               20.016 15 =>     8.706312e-01           14.258 15 =>     5.700576e-01            5.847 15 =>     2.971410e-02           10.700 15 =>     2.881081e-01           13.384 15 =>     5.027062e-01
// 
// --- Finished -- ReMort trials: 2251799813685248	2^51.0 (out of 2^51) trials -- Aqua256 -- 64 bits: 	ps:   8.706e-01    5.701e-01    2.971e-02*   2.881e-01    5.027e-01   => p < 1-7.79087e-02   2^54.32 calls, 2^57.32 bytes	2^37.50 bytes/second	used:  10::16:35:57.89

//					uint64_t x = (stateA += 0xC13FA9A902A6328FUL);
//					uint64_t y = (stateB += 0x91E10DA5C79E7B1DUL);
//					uint64_t z = stateC + x;
//					uint64_t w = stateD + y;
//					stateC = x + (w ^ rotate64(w, 25) ^ rotate64(w, 38));
//					stateD = y + (z ^ rotate64(z, 47) ^ rotate64(z, 19));
//					return stateC + stateD;
					// TernRandom
					// Passes 64TB with no anomalies.
					// Again, period is at least 2 to the 64, with 2 to the 64 guaranteed streams.
					//uint64_t x = (stateA += 0xC13FA9A902A6328FUL);
					//uint64_t y = (stateB += 0x91E10DA5C79E7B1DUL);
					//uint64_t z = stateC ^ x;
					//uint64_t w = stateD ^ y;
					//stateC = x + (rotate64(w, 47) ^ rotate64(w, 13) ^ rotate64(w, 53));
					//stateD = y + (rotate64(z, 19) ^ rotate64(z, 37) ^ rotate64(z, 23));
					//return z;


					// Passes 64TB with no anomalies.
//					uint64_t x = (stateA += 0xC13FA9A902A6328FUL);
//					uint64_t y = (stateB += 0x91E10DA5C79E7B1DUL);
//					uint64_t z = stateC ^ x;
//					uint64_t w = stateD ^ y;
//					stateC = x + (rotate64(w, 47) ^ rotate64(w, 13) ^ rotate64(w, 53));
//					stateD = y + (rotate64(z, 19) ^ rotate64(z, 37) ^ rotate64(z, 23));
//					return z;

//Flops on ReMort around 100 to 170PB.
// 7	                4400  - 5.18470e-02	                4421  + 7.80479e-03	                4457  + 3.97069e-01	                4464  + 5.40934e-01	                4461  + 4.76560e-01	
// 8	              361114  + 1.88447e-01	              361813  + 2.55273e+00	              360692  - 7.20367e-02	              360394  - 5.84423e-01	              360426  - 5.05813e-01	
// 9	            22358321  + 8.47010e-02	            22357601  + 1.92544e-02	            22357218  + 3.33608e-03	            22357509  + 1.42332e-02	            22358927  + 1.75727e-01	
//11	                4502  + 1.70922e+00	                4578  + 6.00814e+00	                4445  + 2.02084e-01	                4334  - 1.49079e+00	                4457  + 3.97069e-01	
//12	            85022540  + 1.56205e-01	            85032187  + 2.07785e+00	            84993202  - 7.76498e+00	            85031598  + 1.89777e+00	            85011701  - 6.08862e-01	
//13	          6948693805  + 1.30826e-02	          6948695706  + 1.88195e-02	          6948688485  + 2.55616e-03	          6948791903  + 1.66719e+00	          6948486784  - 5.61271e+00	
//14	        430510572217  - 8.74323e-01	        430510560593  - 9.07767e-01	        430510157633  - 2.45521e+00	        430510015930  - 3.17865e+00	        430510788585  - 3.66376e-01	
//16	              360441  - 4.70919e-01	              360108  - 1.53904e+00	              360656  - 1.07798e-01	              361243  + 4.21007e-01	              361681  + 1.89885e+00	
//17	          6948643492  - 2.39310e-01	          6948513388  - 4.20235e+00	          6948442937  - 8.38171e+00	          6948500944  - 4.83669e+00	          6948650827  - 1.60961e-01	
//18	        567923307985  - 7.66235e-05	        567921493488  - 5.83949e+00	        567921589565  - 5.23959e+00	        567923486122  + 5.18135e-02	        567922692573  - 6.81245e-01	
//19	      35186112554203  - 1.37286e+01	      35186114499137  - 1.14064e+01	      35186123616266  - 3.38684e+00	      35186121661115  - 4.70865e+00	      35186113137001  - 1.30102e+01	
//21	            22357427  + 1.03960e-02	            22359128  + 2.13175e-01	            22358713  + 1.39831e-01	            22356793  - 1.03203e-03	            22356232  - 2.27323e-02	
//22	        430510172086  - 2.38667e+00	        430510741531  - 4.58335e-01	        430510850931  - 2.60375e-01	        430510305512  - 1.79971e+00	        430510175529  - 2.37048e+00	
//23	      35186121654661  - 4.71337e+00	      35186114291075  - 1.16445e+01	      35186114203340  - 1.17457e+01	      35186121379146  - 4.91720e+00	      35186122477782  - 4.13010e+00	
//24	    2179984607618054  + 6.11676e-01	    2179984614410494  + 8.60396e-01	    2179984605696708  + 5.49001e-01	    2179984599068241  + 3.58776e-01	    2179984606818282  + 5.85176e-01	
//
//              25.239 15 =>   1-3.226669e-02           47.756 15 =>   1-1.442888e-05           40.708 15 =>   1-1.979981e-04           26.469 15 =>   1-2.258465e-02           31.003 15 =>   1-5.546968e-03
//
//--- Finished -- ReMort trials: 2251799813685248	2^51.0 (out of 2^51) trials -- Aqua256 -- 64 bits: 	ps: 1-3.227e-02  1-1.443e-05* 1-1.980e-04  1-2.258e-02  1-5.547e-03   => p <   3.95971e-08  **** SUSPECT **** PRNG: Aqua256	  2^54.32 calls, 2^57.32 bytes	2^37.67 bytes/second	used:   9::12:13:29.09
//					uint64_t x = (stateA += 0xC13FA9A902A6328FUL);
//					uint64_t y = (stateB += 0x91E10DA5C79E7B1DUL);
//					uint64_t z = stateC + x;
//					uint64_t w = stateD + y;
//					stateC = x + rotate64(w, 21);
//					stateD = y + (z ^ rotate64(z, 25) ^ rotate64(z, 50));
//					return w;
					// Fails ReMort after about a minute.
					//uint64_t x = (stateA += 0xC13FA9A902A6328FUL);
					//uint64_t y = (stateB += 0x91E10DA5C79E7B1DUL);
					//uint64_t z = stateC + x;
					//uint64_t w = stateD + y;
					//stateC = x + w;
					//return
					//stateD = y + (z ^ rotate64(z, 44) ^ rotate64(z, 19));

					//// No good.

					//uint64_t x = (stateA += 0xC13FA9A902A6328FUL);
					//uint64_t y = (stateB += 0x91E10DA5C79E7B1DUL);
					//uint64_t z = rotate64(stateC, 23);
					//uint64_t w = rotate64(stateD, 22);
					//stateC = w ^ x;
					//stateD = z ^ y;
					//return z + w;
					
					//Passes at least 32TB with no anomalies, BUT
					//it won't escape if states B, C, and D are all 0.
//					const uint64_t fa = stateA;
//					const uint64_t fb = stateB;
//					const uint64_t fc = stateC;
//					const uint64_t fd = stateD;
//					stateA = fa + 0xC13FA9A902A6328EL;
//					stateB = fa * fc;
//					stateC = fb + fd;
//					stateD = rotate64(fb, 41);
//					return fc;

//// ScruffRandom, passes 64TB of PractRand with no anomalies.
//// Period is at minimum 2 to the 64; all operations have no data dependency on each other.
//// Passes 89PB of ReMort (half the length of normal because of a planned power outage), with a strong final result:
// 7	                2142  - 1.94728e+00	                2221  + 8.17651e-02	                2233  + 2.93058e-01	                2247  + 7.04454e-01	                2197  - 5.05612e-02	
// 8	              181377  + 5.00610e+00	              180722  + 4.83591e-01	              180004  - 9.89892e-01	              180930  + 1.40443e+00	              180728  + 5.03436e-01	
// 9	            11183510  + 2.27016e+00	            11184086  + 2.81899e+00	            11178740  + 6.40369e-03	            11177800  - 4.04517e-02	            11184123  + 2.85627e+00	
//11	                2279  + 2.31158e+00	                2265  + 1.49431e+00	                2240  + 4.76559e-01	                2244  + 6.01349e-01	                2232  + 2.70467e-01	
//12	            42522200  + 3.82542e+00	            42499978  - 2.10962e+00	            42504866  - 4.93860e-01	            42506383  - 2.20975e-01	            42514171  + 5.24772e-01	
//13	          3474332627  - 2.60213e-02	          3474401770  + 1.02359e+00	          3474348730  + 1.25177e-02	          3474399858  + 9.59006e-01	          3474363506  + 1.31452e-01	
//14	        215255398707  - 1.75134e-01	        215255351800  - 2.69976e-01	        215255089196  - 1.17853e+00	        215255036547  - 1.43779e+00	        215255374027  - 2.22486e-01	
//16	              181016  + 1.92530e+00	              180627  + 2.22553e-01	              180317  - 6.65938e-02	              180302  - 8.60668e-02	              180911  + 1.30041e+00	
//17	          3474412641  + 1.43079e+00	          3474421840  + 1.82850e+00	          3474314402  - 2.21375e-01	          3474403315  + 1.07732e+00	          3474361701  + 1.10184e-01	
//18	        283962005391  + 4.26726e-01	        283962199744  + 1.03625e+00	        283962423605  + 2.06802e+00	        283961810830  + 8.30192e-02	        283961756123  + 3.43983e-02	
//19	      17593074156687  + 2.69858e+00	      17593073953524  + 2.54179e+00	      17593068880195  + 1.48036e-01	      17593069404072  + 2.59747e-01	      17593074464838  + 2.94536e+00	
//21	            11177669  - 5.77477e-02	            11184156  + 2.88973e+00	            11184491  + 3.24042e+00	            11178418  - 2.65216e-04	            11177821  - 3.79646e-02	
//22	        215255007254  - 1.59319e+00	        215255329897  - 3.21264e-01	        215255432435  - 1.19573e-01	        215255032292  - 1.45987e+00	        215255066168  - 1.28876e+00	
//23	      17593069274176  + 2.29139e-01	      17593073981337  + 2.56298e+00	      17593073811234  + 2.43477e+00	      17593069401953  + 2.59232e-01	      17593069493214  + 2.81861e-01	
//24	    1089992277004948  - 6.70028e-02	    1089992271968657  - 1.69245e-01	    1089992277309936  - 6.23057e-02	    1089992282125433  - 1.07648e-02	    1089992276720864  - 7.15314e-02	
//
//              23.990 15 =>   1-4.586178e-02           19.854 15 =>     8.651477e-01           11.812 15 =>     3.783828e-01            8.605 15 =>     1.442465e-01           10.630 15 =>     2.852517e-01
//
//--- Finished -- ReMort trials: 1125899906842624	2^50.0 (out of 2^50) trials -- Scruff256 -- 64 bits: 	ps: 1-4.586e-02*   8.651e-01    3.784e-01    1.442e-01    2.853e-01   => p <   7.97358e-01   2^53.32 calls, 2^56.32 bytes	2^37.48 bytes/second	used:   5::10:08:40.45

//					const uint64_t fa = stateA;
//					const uint64_t fb = stateB;
//					const uint64_t fc = stateC;
//					const uint64_t fd = stateD;
//					stateA = fa + 0x9E3779B97F4A7C15UL;
//					stateB = fd * 0xD1342543DE82EF95UL;
//					stateC = fa ^ fb;
//					stateD = rotate64(fc, 21);
//					return fd - fc;

					// SplurgeRandom
					// Passes 64TB with no anomalies.
					// Period is exactly 2 to the 64; there are 2 to the 192 possible streams.
					// It permits O(1) skip to any point in its period.
					// ReMort results are excellent:
//  7	                4443  + 1.75929e-01	                4280  - 4.13579e+00	                4412  - 2.21866e-03	                4444  + 1.88780e-01	                4493  + 1.37341e+00	
//  8	              361632  + 1.68070e+00	              360150  - 1.37045e+00	              360172  - 1.28604e+00	              360104  - 1.55560e+00	              360712  - 5.52733e-02	
//  9	            22355365  - 1.11647e-01	            22357010  + 1.89571e-04	            22347579  - 3.92362e+00	            22347615  - 3.89351e+00	            22356216  - 2.37641e-02	
// 11	                4428  + 3.75169e-02	                4440  + 1.40092e-01	                4339  - 1.31270e+00	                4508  + 1.95348e+00	                4307  - 2.64818e+00	
// 12	            85023227  + 2.20651e-01	            85010250  - 8.79210e-01	            85011373  - 6.65642e-01	            85011783  - 5.95063e-01	            85024042  + 3.11503e-01	
// 13	          6948701511  + 4.27757e-02	          6948659782  - 8.63022e-02	          6948699154  + 3.18792e-02	          6948723387  + 2.20200e-01	          6948766329  + 9.69046e-01	
// 14	        430510196988  - 2.27084e+00	        430510251682  - 2.02656e+00	        430511201160  + 5.52589e-04	        430511176348  - 2.04726e-04	        430510130819  - 2.58495e+00	
// 16	              360741  - 3.49041e-02	              360025  - 1.90095e+00	              361602  + 1.55370e+00	              360978  + 4.31419e-02	              360831  - 1.36929e-03	
// 17	          6948724022  + 2.27407e-01	          6948685554  + 2.37074e-04	          6948614618  - 6.98186e-01	          6948619383  - 6.05926e-01	          6948725284  + 2.42076e-01	
// 18	        567924232857  + 1.48476e+00	        567922936739  - 2.51381e-01	        567923966950  + 7.49370e-01	        567923959309  + 7.31918e-01	        567922604557  - 8.87682e-01	
// 19	      35186138110766  + 3.63841e-01	      35186139446068  + 6.86084e-01	      35186137714471  + 2.87708e-01	      35186137717971  + 2.88341e-01	      35186139730899  + 7.67937e-01	
// 21	            22347010  - 4.41483e+00	            22356956  + 5.51265e-06	            22355480  - 9.59848e-02	            22346693  - 4.70106e+00	            22347041  - 4.38733e+00	
// 22	        430511166900  - 8.24136e-04	        430510225413  - 2.14215e+00	        430510295094  - 1.84256e+00	        430511282982  + 2.19663e-02	        430511164773  - 1.02077e-03	
// 23	      35186137335772  + 2.23296e-01	      35186139464900  + 6.91354e-01	      35186138395295  + 4.24009e-01	      35186137588972  + 2.65459e-01	      35186138900174  + 5.42099e-01	
// 24	    2179984564759586  - 1.84508e-02	    2179984563561999  - 2.60769e-02	    2179984564353549  - 2.08890e-02	    2179984564180771  - 2.19723e-02	    2179984563204771  - 2.86064e-02	
// 
//               11.308 15 =>     3.381111e-01           14.337 15 =>     5.760747e-01           12.895 15 =>     4.652238e-01           15.087 15 =>     6.275381e-01           14.824 15 =>     6.091584e-01
// 
// --- Finished -- ReMort trials: 2251799813685248	2^51.0 (out of 2^51) trials -- splurge256 -- 64 bits: 	ps:   3.381e-01*   5.761e-01    4.652e-01    6.275e-01    6.092e-01   => p <   3.53287e-01   2^54.32 calls, 2^57.32 bytes	2^37.06 bytes/second	used:  14::13:12:32.98

//					uint64_t fa = (stateA += 0xDB4F0B9175AE2165L);
//					uint64_t fb = (stateB += 0xBBE0563303A4615FL);
//					uint64_t fc = (stateC += 0xA0F2EC75A1FE1575L);
//					uint64_t fd = (stateD += 0x89E182857D9ED689L);
//					fb += fa ^ rotate64(fa, 11) ^ rotate64(fa, 46);
//					fc += fb ^ rotate64(fb, 17) ^ rotate64(fb, 34);
//					fd += fc ^ rotate64(fc,  5) ^ rotate64(fc, 58);
//					fa += fd ^ rotate64(fd, 47) ^ rotate64(fd, 38);
//					return fa;

//					uint64_t fa = (stateA += 0xDB4F0B9175AE2165L);
//					uint64_t fb = (stateB += 0xBBE0563303A4615FL);
//					uint64_t fc = (stateC += 0xA0F2EC75A1FE1575L);
//					uint64_t fd = (stateD += 0x89E182857D9ED689L);
//					fb += rotate64(fa, 25);// ^ __rolq(fb, 58);
//					fc += rotate64(fb, 46);// ^ __rolq(fc, 11);
//					fd += rotate64(fc, 37);// ^ __rolq(fd, 21);
//					fa += rotate64(fd, 18);// ^ __rolq(fa, 37);
//					return fa ^ fb ^ fc ^ fd;

//					uint64_t fa = (stateA += 0xDB4F0B9175AE2165L);
//					uint64_t fb = (stateB += 0xBBE0563303A4615FL);
//					fb += fa ^ rotate64(fa, 11) ^ rotate64(fa, 43);
//					fa += fb ^ rotate64(fb, 17) ^ rotate64(fb, 37);
//					return fa ^ fb;

//SparkleRandom
// Period is 2 to the 64.
// Has 2 to the 64 streams.
// Passes 64TB of PractRand without anomalies.
//uint64_t n = (stateA += 0xDB4F0B9175AE2165L);
//uint64_t o = (stateB += 0xBBE0563303A4615FL);
//n = (n ^ rotate64(n, 17) ^ rotate64(n, 53)) * 0xC6BC279692B5C323L;
//o = (o ^ rotate64(o, 11) ^ rotate64(o, 47)) * 0xABC98388FB8FAC03L;
//return (n ^ rotate64(n, 19) ^ rotate64(n, 43)) + (o ^ rotate64(o, 13) ^ rotate64(o, 53));

// ^ (stateC += 0x89E182857D9ED689L)
//uint64_t n = ((stateA += 0xDB4F0B9175AE2165L) ^ 0xF1357AEA2E62A9C5L) * 0xC6BC279692B5C323L;
//uint64_t o = ((stateB += 0xBBE0563303A4615FL) ^ 0x91E10DA5C79E7B1DL) * 0xABC98388FB8FAC03L;
//return (rotate64(o, 19) ^ rotate64(o, 29) ^ rotate64(o, 47)) + (rotate64(n, 17) ^ rotate64(n, 37) ^ rotate64(n, 43));

// SportyRandom
// Period is 2 to the 64.
// Has 2 to the 192 streams, but most likely only 2 to the 189 are decorrelated well.
// (This issue is simply that the lowest bit of states B, C, and D shouldn't be used to help reduce correlation.)
// Nearby stateA and stateB aren't visually correlated, but I don't know about states C and D yet.
// Passes 64TB of PractRand without anomalies.
//  6	                   1  + 2.59071e+00	                   0  - 2.29283e-01	                   0  - 2.29283e-01	                   0  - 2.29283e-01	                   1  + 2.59071e+00	
//  7	                4420  + 5.37216e-03	                4340  - 1.27844e+00	                4474  + 7.84960e-01	                4570  + 5.43241e+00	                4399  - 5.89271e-02	
//  8	              360636  - 1.30769e-01	              360665  - 9.81841e-02	              359723  - 3.53999e+00	              360493  - 3.59605e-01	              359977  - 2.12767e+00	
//  9	            22351780  - 1.19319e+00	            22351832  - 1.16929e+00	            22354463  - 2.75522e-01	            22353597  - 5.01340e-01	            22352469  - 8.96082e-01	
// 11	                4389  - 1.54643e-01	                4455  + 3.60042e-01	                4337  - 1.38258e+00	                4441  + 1.51585e-01	                4289  - 3.60323e+00	
// 12	            85026329  + 6.49889e-01	            85012121  - 5.39852e-01	            85010223  - 8.84710e-01	            85013245  - 3.75579e-01	            85025498  + 5.12702e-01	
// 13	          6948531227  - 3.37076e+00	          6948701909  + 4.47734e-02	          6948670651  - 2.66944e-02	          6948715925  + 1.44201e-01	          6948791117  + 1.64293e+00	
// 14	        430512314021  + 2.95701e+00	        430512157481  + 2.19341e+00	        430511718857  + 6.60187e-01	        430511670457  + 5.45757e-01	        430512056322  + 1.76051e+00	
// 16	              360377  - 6.28493e-01	              362142  + 4.60279e+00	              360844  - 2.36019e-04	              360646  - 1.19006e-01	              359562  - 4.62036e+00	
// 17	          6948872104  + 5.07742e+00	          6948745308  + 5.36156e-01	          6948737183  + 4.02915e-01	          6948629542  - 4.31047e-01	          6948678040  - 5.58656e-03	
// 18	        567922791136  - 4.82451e-01	        567922290145  - 1.84791e+00	        567923860179  + 5.24149e-01	        567923138025  - 5.48882e-02	        567922606621  - 8.82528e-01	
// 19	      35186137233389  + 2.07281e-01	      35186137859411  + 3.14517e-01	      35186135494756  + 2.63015e-02	      35186136324749  + 9.12646e-02	      35186137603665  + 2.68017e-01	
// 21	            22353877  - 4.20988e-01	            22350249  - 2.00542e+00	            22351665  - 1.24692e+00	            22353557  - 5.13391e-01	            22354792  - 2.07317e-01	
// 22	        430511500931  + 2.30767e-01	        430512115457  + 2.00780e+00	        430512125346  + 2.05074e+00	        430511756427  + 7.56515e-01	        430511695847  + 6.04428e-01	
// 23	      35186136761669  + 1.41194e-01	      35186137895169  + 3.21315e-01	      35186136357335  + 9.46139e-02	      35186136230225  + 8.18905e-02	      35186136686953  + 1.31887e-01	
// 24	    2179984565218962  - 1.58747e-02	    2179984563474564  - 2.66852e-02	    2179984566275212  - 1.06859e-02	    2179984566769349  - 8.60984e-03	    2179984565105696  - 1.64919e-02	
// 
//               15.666 15 =>     6.659296e-01           17.347 15 =>     7.622138e-01           11.911 15 =>     3.860725e-01            9.567 15 =>     2.057609e-01           17.339 15 =>     7.611422e-01
// 
// --- Finished -- ReMort trials: 2251799813685248	2^51.0 (out of 2^51) trials -- sporty256 -- 64 bits: 	ps:   6.659e-01    7.622e-01    3.861e-01    2.058e-01*   7.611e-01   => p <   7.32983e-01   2^54.32 calls, 2^57.32 bytes	2^36.60 bytes/second	used:  20::00:24:06.92
// 
// The results above show this passes 179PB of ReMortality testing, with a strong final p-value of 7.32983e-01 .

//					uint64_t x = (stateA += 0xDB4F0B9175AE2165L);
//					x ^= x >> 32;
//					x *= (stateB += 0xBBE0563303A4615FL) | 1L;
//					x ^= x >> 33;
//					x *= (stateC += 0xA0F2EC75A1FE1575L) | 1L;
//					x ^= x >> 32;
//					x *= (stateD += 0x89E182857D9ED689L) | 1L;
//					x ^= x >> 31;
//					return x;

//SprawlRandom
// Passes 64TB of PractRand without anomalies.
// Nearby stateA and stateB initial values are visually correlated, like LaserRandom.
// Has a period of 2 to the 64, with 2 to the 192 streams.
// 190 of those should be decorrelated, but apparently aren't.
//uint64_t x = (stateA += 0xDB4F0B9175AE2165L);
//uint64_t y = (stateB += 0xBBE0563303A4615FL);
//x ^= x >> 35;
//y ^= y >> 29;
//x *= (stateC += 0xA0F2EC75A1FE1575L) | 1L;
//y *= (stateD += 0x89E182857D9ED689L) | 1L;
//x ^= y;
//return x ^ x >> 31 ^ x >> 17;

// SprawlRandom (variant)
// Period is 2 to the 64, still.
// Has 2 to the 192 streams, but nearby initial values for streams are very correlated...
// Still, passes 64TB of PractRand with no anomalies.
//				uint64_t x = (stateA += 0xDB4F0B9175AE2165L);
//				uint64_t y = (stateB += 0xBBE0563303A4615FL);
//				x ^= x >> 41;
//				y ^= y >> 23;
//				x *= (stateC += 0xA0F2EC75A1FE1575L) | 1L;
//				y *= (stateD += 0x89E182857D9ED689L) | 1L;
//				x ^= y;
//				return x ^ x >> 31 ^ x >> 17;

// SpoonRandom
// Passes 64TB of PractRand with no anomalies.
// Period is 2 to the 64, with 2 to the 126 streams.
// Has a smaller state than SportyRandom, and should be faster, but is otherwise similar.
// States B and C must be odd, but stateA has no restrictions.
// Streams do not appear to be correlated at all when stateA and stateB are similar (and stateC is 1).
// NOTE: This is sensitive to the choice of increment for at least states B and C.
//       All of the increments are from the R4 sequence, but stateB had 1 subtracted, and stateC had 1 added.
		uint64_t x = (stateA += 0xDB4F0B9175AE2165L);
		x ^= x >> 32;
		x *= (stateB += 0xBBE0563303A4615EL);
		x ^= x >> 33;
		x *= (stateC += 0xA0F2EC75A1FE1576L);
		x ^= x >> 31;
		return x;

				}
				std::string mars256::get_name() const { return "mars256"; }
				void mars256::walk_state(StateWalkingObject *walker) {
					walker->handle(stateA);
					stateA |= 1L;
					walker->handle(stateB);
					stateB |= 1L;
					walker->handle(stateC);
					stateC |= 1L;
					walker->handle(stateD);
					stateD |= 1L;
				}

				Uint64 marsgwt::raw64() {
//					const uint32_t fa = stateA;
//					const uint32_t fb = stateB;
//					const uint32_t fc = stateC;
//					const uint32_t fd = stateD;
//					stateA = 0xDEF95U * fd;
//					stateB = fa + 0x92B5C323U;
//					stateC = rotate32(fb, 23) - fd;
//					stateD = fb ^ fc;

//					stateA = rotate32(fd, 23) * 0xDEF95U;
//					stateB = fa + 0x92B5C323U;
//					stateC = rotate32(fb, 11) + fd;
//					stateD = fa ^ fb - fc;
//					return fd;

//    //// Passess 64TB with no anomalies.
//	stateA = fb ^ fc ^ fd;
//	stateB = rotate32(fa, 21);
//	stateC = fa + fb;
//	stateD = fc + 0x9E3779B9U;
//	return fc;
//					const uint32_t fa = stateA;
//					const uint32_t fb = stateB;
//					const uint32_t fc = stateC;
//					const uint32_t fd = stateD;
//					stateA = fb ^ fc ^ fd;
//					stateB = rotate32(fa, 21);
//					stateC = fa + fb;
//					stateD = fd + 0x3943D869U;
//					return fc;

					//// ChopRandom
					//// Passes 64TB of PractRand without any anomalies.
					//// Generally quite fast on GWT.
					//const uint32_t fa = stateA;
					//const uint32_t fb = stateB;
					//const uint32_t fc = stateC;
					//const uint32_t fd = stateD;
            	    //stateA = rotate32(fb ^ fc, 26);
            	    //stateB = rotate32(fc ^ fd, 11);
            	    //stateC = fa ^ fb + fc;
            	    //stateD = fd + 0xADB5B165U;
            	    //return fc;

					//// Modified TycheI2; passes 4TB without anomalies
					//stateB = rotate32(stateB,  7) ^ stateC; stateC += 0xADB5B165U;//stateD;
					//stateD = rotate32(stateD,  8) ^ stateA; stateA += stateB;
					//stateB = rotate32(stateD, 12) ^ stateC;// stateC += stateD;
					//stateD = rotate32(stateD, 16) ^ stateA; stateA += stateB;
					//return stateB;

					//// Passes 32TB with no anomalies, but has two anomalies at 64TB:
					////  [Low8/32]DC6-9x1Bytes-1           R= +12.5  p =  3.4e-5   mildly suspicious
					////  [Low8/32]FPF-14+6/16:(19,14-0)    R=  +7.1  p =  3.3e-6   unusual
					const uint32_t fa = stateA;
					const uint32_t fb = stateB;
					const uint32_t fc = stateC;
					const uint32_t fd = stateD;
            	    stateA = rotate32(fb ^ fc, 26);
            	    stateB = rotate32(fc ^ fd, 11);
            	    stateC = fa ^ fb + fc;
            	    stateD = fd + 0xADB5B165U;
            	    return (uint64_t)(fa + fd) << 32 ^ (int32_t)(fb + fc);
				}
				//0x3943D8696D4A3CDDUL;
				std::string marsgwt::get_name() const { return "lizard128"; }
				void marsgwt::walk_state(StateWalkingObject *walker) {
					walker->handle(stateA);
					walker->handle(stateB);
					walker->handle(stateC);
					walker->handle(stateD);
				}
				Uint64 lizard256::raw64() {
//					const uint64_t fa = stateA;
//					const uint64_t fb = stateB;
//					const uint64_t fc = stateC;
//					const uint64_t fd = stateD;
//					stateA = rotate64(fd, 37) + 0xD1342543DE82EF95UL;
//					stateB = rotate64(fc, 47) + fd ^ fa;
//					stateC = rotate64(fa, 12) + 0xC6BC279692B5C323UL;
//					stateD = fa ^ fb + fc;
//					return fd;
					
					////Passes 64TB without anomalies, as well as 500TB of hwd so far.
					////This is almost the same as mars256 above, but attempts to be
					////a little faster by removing `- fd` after the rotation.
					////it also changes `fb ^ fc` to `fb + fc`. It isn't actually any
					////faster, because the pipeline stalls finishing the multiply earlier.
					//const uint64_t fa = stateA;
					//const uint64_t fb = stateB;
					//const uint64_t fc = stateC;
					//const uint64_t fd = stateD;
					//stateA = 0xD1342543DE82EF95UL * fd;
					//stateB = fa + 0xC6BC279692B5C323UL;
					//stateC = rotate64(fb, 47);
					//stateD = fb + fc;
					//return fd;

////hwd issues at 100TB
//					const uint64_t fa = stateA;
//					const uint64_t fb = stateB;
//					const uint64_t fc = stateC;
//					const uint64_t fd = stateD;
//	stateA = fa + 0xC6BC279692B5C323UL;
//	stateB = rotate64(fa, 23) + fd;
//	stateC = rotate64(fd, 35) ^ fb;
//	stateD = rotate64(fb, 44) + fc;
//	return fd;
//	stateA = fa + 0xD1342543DE82EF95UL;
//	stateB = rotate64(fa, 44) + fd;
//	stateC = fb ^ fb << 7;
//	stateD = fc ^ fc >> 9;
//	return fd;

//const uint64_t fa = stateA;
//const uint64_t fb = stateB;
//const uint64_t fc = stateC;
//stateA = fc + 0xC6BC279692B5C323UL;
//stateB = rotate64(fc, 25) ^ fa;
//return stateC ^= rotate64(fa, 44) + fb;

	// const uint64_t fa = stateA;
	// const uint64_t fb = stateB;
	// const uint64_t fc = stateC;
	// const uint64_t fd = stateD;
	// stateA = fc * fd;
	// stateB = fc + rotate64(fa, 39);
	// stateC = fb ^ fd;
	// stateD = fd + 0x9E3779B97F4A7C16UL;
	// return fc;
	
	//// Stranger256 passes a ridiculous amount of extsat (over 2 exabytes) and hwd (at least 7 PB) testing, as well as 64TB of PractRand.
	//// This generator is very solid and uses no multiplication. It isn't astonishingly fast, but is decent.
	//// stateA and stateB form an interleaved 64-bit xorshift RNG. The worst one. That gets fed into the rest, which is chaotic.
	//const uint64_t fa = stateA;
	//const uint64_t fb = stateB;
	//const uint64_t fc = stateC;
	//const uint64_t fd = stateD;
	//stateA = fb ^ fb << 7;
	//stateB = fa ^ fa >> 9;
	//stateC = rotate64(fd, 39) - fb;
	//stateD = fa - fc + 0xC6BC279692B5C323UL;
	//return fc;

////passes at least 8TB with no anomalies, but not especially fast.
//	const uint64_t fa = stateA;
//	const uint64_t fb = stateB;
//	const uint64_t fc = stateC;
//	const uint64_t fd = stateD;
//	stateA = rotate64(fd, 17) ^ fc;
//	stateB = rotate64(fa, 59) + fd;
//	stateC = rotate64(fb, 41) ^ fa;
//	stateD = rotate64(fc, 11) + 0xC6BC279692B5C323UL;
//	return fc;

//	const uint64_t fa = stateA;
//	const uint64_t fb = stateB;
//	const uint64_t fc = stateC;
//	const uint64_t fd = stateD;
//	stateA = fd ^ fc;
//	stateB = rotate64(fa, 41);
//	stateC = fa + fb;
//	stateD = fc + 0x9E3779B97F4A7C15UL;
//	return fc;

	////Passes 64TB with no anomalies.
//	const uint64_t fa = stateA;
//	const uint64_t fb = stateB;
//	const uint64_t fc = stateC;
//	const uint64_t fd = stateD;
//	stateA = fb ^ fc ^ fd;
//	stateB = rotate64(fa, 42);
//	stateC = fa + fb;
//	stateD = fc + 0x9E3779B97F4A7C15UL;
//	return fc;

	// this whole structure seems good and then has problems later.
//	const uint64_t fa = stateA;
//	const uint64_t fb = stateB;
//	const uint64_t fc = stateC;
//	const uint64_t fd = stateD;
//	stateA = fc + fd ^ fb;
//	stateB = rotate64(fa, 36);
//	stateC = fa + fb;
//	stateD = fd + 0xC259492BFC623C6BULL;
//	return fc;
//0x3943D8696D4A3CDDUL is ~0xC6BC279692B5C323ULL + 1ULL;
//0xFC0072FA0B15F4FDULL
//	const uint64_t fa = stateA;
//	const uint64_t fb = stateB;
//	const uint64_t fc = stateC;
//	const uint64_t fd = stateD;
    //// Good rotations: 11-13, 19, 22, 26, 38, 41-42, 46, 50, 53-55
	//// 23 gets two mildly suspicious at 4TB
	//// 39 gets a mildly suspicious at 8TB
	//// 45 gets a mildly suspicious at 512GB, very suspicious at 1TB; both [Low4/16]DC6-9x1Bytes-1
	//// 37 gets a mildly suspicious at 2TB, VERY SUSPICIOUS at 4TB; both [Low8/32]DC6-9x1Bytes-1
//	stateA = rotate64(fc + fb, 37);
//	stateB = fc ^ fd;
//	stateC = fa + fb;
//	stateD = fd + 0x9E3779B97F4A7C15UL;
//	return fc;

    //// 45 gets a mildly suspicious at 8TB, [Low4/16]DC6-9x1Bytes-1
	//// 39 gets an unusual at 256GB, [Low4/32]Gap-16:B
	//// 37 gets an unusual at 2TB, [Low8/64]BCFN(2+1,13-0,T)
//	stateA = rotate64(fc - fb, 37);
//	stateB = fc ^ fd;
//	stateC = fa + fb;
//	stateD = fd + 0x9E3779B97F4A7C15UL;
//	return fa ^ fb;
    //// 37 gets an unusual at 4TB, suspicious at 8TB, both [Low8/32]DC6-9x1Bytes-1
//	stateA = rotate64(fc - fb, 37);
//	stateB = fc ^ fa;
//	stateC = fd + fb;
//	stateD = fd + 0x9E3779B97F4A7C15UL;
//	return fa ^ fc;
    //// 41 gets an unusual at 4TB, BCFN(2+0,13-0,T)
//	stateA = rotate64(fc + fb, 41);
//	stateB = fc ^ fa;
//	stateC = fd + fb;
//	stateD = fd + 0x9E3779B97F4A7C15UL;
//	return fa ^ fc;
			//// Good rotations: 9-13, 15, 18-21, 23, 25-27, 30, 36-37, 42, 44, 47, 49, 51, 53, 55, 57, 59
			//// but this isn't super-fast.
//			stateA = rotate64(fb * 0xF1357AEA2E62A9C5UL, 20);
//			stateB = fc ^ fd;
//			stateC = fa + fb;
//			stateD = fd + 0x9E3779B97F4A7C15UL;
//			return fc;

			//// Good rotations: 2, 4-6, 9-10, 13, 16, 18-21, 23, 27, 30, 33, 36-37, 40, 43-44, 47, 49-50, 53-56, 58, 61
			//// rotation 37 (then 26) passes 64TB with no anomalies.
			//// Calling this one TrimRandom.
			//// Passes 2PB of hwd.
			//// Fails after 2^58.24 bytes of ExtSat.
			//stateA = rotate64(fb + fc, 37);
			//stateB = rotate64(fc ^ fd, 26);
			//stateC = fa + fb;
			//stateD = fd + 0x9E3779B97F4A7C15UL;
			//return fc;

			//// Works pretty well; seed 0 gets a "mildly suspicious" at 4TB, seed 1 is fine at 32TB
			//stateA = rotate64(fb + fc, 7);
			//stateB = rotate64(fc ^ fd, 55);
			//stateC = fa + fb;
			//stateD = fd + 0xC6BC279692B5C323UL;
			//return fc;

			//// Trim256 (original)
			//// Passes 64TB without anomalies.
			//// Remortality test (seed 0) has a phase of suspect results, but recovers after a few hours.
			//// --- Finished -- trials: 2251799813685248	2^51.0 (out of 2^51) trials -- Trim256 -- 64 bits: 	chi2p: 1-2.796e-03 1-5.488e-03 1-1.648e-04 1-4.027e-01 => p < 8.522e-03   2^57.32 calls, 2^60.32 bytes	2^41.15 bytes/second	used:   6::19:15:33.93
			//// A much longer remortality test (with seed 1) passes at least 2 exabytes and is still running.
			//// It is currently at a p-value of 4.019e-01 10 days and 22 hours into testing.
//			stateA = rotate64(fb + fc, 35);
//			stateB = rotate64(fc ^ fd, 46);
//			stateC = fa + fb;
//			stateD = fd + 0x06A0F81D3D2E35EFL;
//			return fc;

		  //// Trim256 (original)
          //// Passes 64TB without anomalies (seed 1)
		  //// On the new boolbin test, seems strong at about 4 days into testing.
		  //// This is slightly faster than the previous Trim256 when run in Java.
		  //// It also has stronger avalanche properties.
//		  stateA = rotate64(fb ^ fc, 57);
//		  stateB = rotate64(fc ^ fd, 11);
//		  stateC = fa + fb;
//		  stateD = fd + 0xADB5B12149E93C39UL;
//		  return fc;

//					const uint64_t fa = stateA;
//					const uint64_t fb = stateB;
//					const uint64_t fc = stateC;
//					const uint64_t fd = stateD;
//					stateA = fd * 0xD1342543DE82EF95UL;
//					stateB = rotate64(fa + fd, 29);
//					stateC = fc + 0xC6BC279692B5C323UL;
//					stateD = fb ^ fc;
//					return fd;
					//// with seed 0, passes 4TB without anomalies, but not especially fast. Not slow either...
//					stateA = 0xD1342543DE82EF95UL * fc;
//					stateB = fa ^ fb ^ fc;
//					stateC = rotate64(fb, 41) + fd;
//					stateD = fd + 0x9E3779B97F4A7C15UL;
//					return fa;

//// SlashRandom, passes 4TB without anomalies (seed 1), fails at 16TB:
//// BCFN(2+0,13-0,T)                  R= +23.9  p =  2.7e-12    FAIL
//// Extremely fast.
//// No period guarantees.
//stateA = rotate64(fc, 44);
//stateB = fa + fc;
//stateC = fb ^ fd;
//stateD = fb + 0xDE916ABCC965815BUL;
//return stateB;

//stateA = rotate64(fc, 12); // try 12, 11, 9
//stateB = fa + fc;
//stateC = fb ^ fd;
//stateD = fb + 0xC6BC279692B5C323UL;//0x9E3779B97F4A7C15UL;//0xDE916ABCC965815BUL;
//return stateB;

//stateA = fb ^ rotate64(fc, 19); // 21, 19
//stateB = fa ^ fd;
//stateC = fb + fd;
//stateD = fd + 0xC6BC279692B5C323UL;//0x9E3779B97F4A7C15UL;//0xDE916ABCC965815BUL;
//return stateA;

// bad rotations: 42 (BRank at 8TB)
//stateA = fb + rotate64(fc, 11); // 39, 36, 28, 25, 23, 22, 20, 19, 13
//stateB = fa ^ fc;
//stateC = fa ^ fd;
//stateD = fd + 0xDE916ABCC965815BUL;
//return stateA;

// Whisker256
// Passes 64TB of PractRand with no anomalies. Passes over 179 PB of Remortality!
// Very fast on the JVM.
// Remortality results after over 10 days of GPU testing, using the WIP name Slash256:
//
// 6	                   0  - 2.29283e-01	                   0  - 2.29283e-01	                   1  + 2.59071e+00	                   0  - 2.29283e-01	                   0  - 2.29283e-01	
// 7	                4293  - 3.37831e+00	                4466  + 5.86116e-01	                4512  + 2.12538e+00	                4453  + 3.24827e-01	                4438  + 1.18467e-01	
// 8	              360824  - 2.36749e-03	              360318  - 7.93868e-01	              360909  + 8.61968e-03	              360076  - 1.67404e+00	              360708  - 5.84486e-02	
// 9	            22356023  - 3.80149e-02	            22356356  - 1.55120e-02	            22358386  + 9.28917e-02	            22359279  + 2.43684e-01	            22356006  - 3.94298e-02	
//11	                4574  + 5.71665e+00	                4490  + 1.26962e+00	                4496  + 1.48127e+00	                4343  - 1.17838e+00	                4371  - 4.41083e-01	
//12	            85007506  - 1.52586e+00	            85015567  - 1.30333e-01	            85021257  + 6.55781e-02	            85002709  - 3.08181e+00	            85011804  - 5.91554e-01	
//13	          6948736276  + 3.89221e-01	          6948738348  + 4.20853e-01	          6948587234  - 1.35509e+00	          6948719548  + 1.79099e-01	          6948705114  + 6.25228e-02	
//14	        430509873970  - 3.99695e+00	        430509863921  - 4.05842e+00	        430511730558  + 6.89485e-01	        430511616945  + 4.31908e-01	        430509899668  - 3.84188e+00	
//16	              361240  + 4.14551e-01	              360544  - 2.64990e-01	              360390  - 5.94648e-01	              360925  + 1.42748e-02	              360218  - 1.11823e+00	
//17	          6948641135  - 2.67773e-01	          6948625204  - 5.02088e-01	          6948656819  - 1.08450e-01	          6948650158  - 1.67465e-01	          6948777970  + 1.26349e+00	
//18	        567924866382  + 4.24016e+00	        567922845548  - 3.87363e-01	        567922710497  - 6.42549e-01	        567924655404  + 3.16558e+00	        567920876906  - 1.04631e+01	
//19	      35186129635851  - 6.81509e-01	      35186131673312  - 2.32376e-01	      35186140882947  + 1.14605e+00	      35186138944166  + 5.53075e-01	      35186133510776  - 2.96832e-02	
//21	            22358005  + 5.02670e-02	            22356118  - 3.05838e-02	            22356265  - 2.06764e-02	            22358551  + 1.15381e-01	            22359230  + 2.33560e-01	
//22	        430511692619  + 5.96803e-01	        430509975720  - 3.40093e+00	        430509938369  - 3.61413e+00	        430511688233  + 5.86519e-01	        430511551341  + 3.10484e-01	
//23	      35186138657076  + 4.83430e-01	      35186131581656  - 2.47512e-01	      35186131867230  - 2.01927e-01	      35186138885530  + 5.38470e-01	      35186142677830  + 1.88547e+00	
//24	    2179984571129474  + 3.53681e-07	    2179984579923680  + 3.57008e-02	    2179984568845378  - 2.33535e-03	    2179984560074928  - 5.57756e-02	    2179984567228868  - 6.88027e-03	
//
//              21.782 15 =>   1-8.293888e-02           12.376 15 =>     4.234616e-01           12.149 15 =>     4.083835e-01           12.310 15 =>     4.182843e-01           20.464 15 =>     8.837629e-01
//
//--- Finished -- ReMort trials: 2251799813685248	2^51.0 (out of 2^51) trials -- Slash256 -- 64 bits: 	ps: 1-8.294e-02*   4.235e-01    4.084e-01    4.183e-01    8.838e-01   => p <   3.88610e-01   2^54.32 calls, 2^57.32 bytes	2^37.55 bytes/second	used:  10::08:30:06.64
//
// Boolbin results after 6 days of GPU testing:
// 1	                   0  - 9.76563e-04	                   0  - 9.76563e-04	                   0  - 9.76563e-04	                   1  - 1.27592e-01	                   0  - 9.76563e-04	                   1  - 1.27592e-01	
// 2	                   0  - 3.07617e-02	                   0  - 3.07617e-02	                   0  - 3.07617e-02	                  38  + 2.65640e-01	                   0  - 3.07617e-02	                  27  - 1.80953e+00	
// 3	                   0  - 6.35742e-01	                   0  - 6.35742e-01	                   1  + 2.08707e-01	                 599  + 2.45839e+00	                   1  + 2.08707e-01	                 584  + 8.74405e-01	
// 4	                  12  + 5.47981e-01	                  16  + 4.10025e+00	                  14  + 1.91153e+00	                6707  + 2.77571e-01	                   7  - 7.49184e-01	                6738  + 8.21918e-01	
// 5	                 103  - 1.52979e+00	                 119  + 6.07804e-02	                 103  - 1.52979e+00	               62076  - 2.36386e-01	                 127  + 9.76597e-01	               62534  + 1.82320e+00	
// 6	                1158  + 1.70884e-01	                1131  - 1.48136e-01	                1167  + 4.61679e-01	              477314  + 5.51827e+00	                1167  + 4.61679e-01	              475040  - 8.98628e-01	
// 7	                9531  + 2.85187e-01	                9499  + 4.21697e-02	                9505  + 7.12780e-02	             3062427  - 3.24775e+00	                9334  - 2.21827e+00	             3064330  - 5.11610e-01	
// 8	               67736  + 5.80921e-01	               67572  + 1.71932e-02	               67175  - 1.95022e+00	            16984615  - 8.59199e-01	               67396  - 2.98237e-01	            16994688  + 2.30118e+00	
// 9	              420304  + 1.10131e-02	              421199  + 2.20692e+00	              419592  - 9.86819e-01	            82222703  + 6.07009e-01	              419395  - 1.68294e+00	            82227525  + 1.71848e+00	
//10	             2312285  + 4.21624e-01	             2310965  - 4.79290e-02	             2310623  - 1.97032e-01	           351718611  + 9.60458e-01	             2311503  + 1.82120e-02	           351706192  + 1.01005e-01	
//11	            11349609  + 9.23948e-01	            11343540  - 7.06446e-01	            11342654  - 1.21779e+00	          1342863830  + 5.25370e-02	            11347070  + 4.30398e-02	          1342809147  - 1.59524e+00	
//12	            50114657  + 4.59588e-02	            50107562  - 6.20741e-01	            50112879  - 1.35300e-03	          4612996819  + 3.42317e-01	            50130065  + 5.71659e+00	          4612942594  - 4.54975e-02	
//13	           200453593  + 5.34859e-03	           200442386  - 5.16135e-01	           200441749  - 5.82806e-01	         14351386621  - 8.73651e-02	           200484710  + 5.15723e+00	         14351520804  + 6.79811e-01	
//14	           730227545  + 7.73169e-02	           730254681  + 1.64418e+00	           730210409  - 1.26791e-01	         40662278069  - 1.74978e-01	           730161066  - 4.76142e+00	         40662314567  - 5.63140e-02	
//15	          2434091299  + 2.47180e-01	          2433962684  - 4.45098e+00	          2434064139  - 2.84466e-03	        105420838577  - 9.68145e-02	          2434060082  - 1.83784e-02	        105421030828  + 7.89405e-02	
//16	          7454310439  - 4.86592e-02	          7454357115  + 1.02418e-01	          7454339267  + 1.28385e-02	        251106424016  - 6.28927e-01	          7454281785  - 3.05221e-01	        251106834260  + 6.56813e-04	
//17	         21047663386  + 9.96758e-01	         21047419882  - 4.62484e-01	         21047377145  - 9.49926e-01	        551449161716  - 2.24540e+00	         21047435619  - 3.26713e-01	        551450086724  - 6.39221e-02	
//18	         54957742024  + 2.01159e+00	         54957346947  - 7.12675e-02	         54957043168  - 2.44228e+00	       1119921866826  + 1.38996e+00	         54957813505  + 2.96949e+00	       1119920796366  + 2.80359e-02	
//19	        133054785685  + 1.67122e-04	        133055315680  + 2.14885e+00	        133055222076  + 1.46237e+00	       2108858490638  + 2.13125e+00	        133054508611  - 5.57508e-01	       2108858509424  + 2.16919e+00	
//20	        299373423415  + 9.22909e-02	        299373819216  + 1.05510e+00	        299372149410  - 4.09918e+00	       3690493304628  - 7.73787e+00	        299372106602  - 4.42211e+00	       3690497113848  - 6.38139e-01	
//21	        627258798131  + 4.73499e-01	        627258079421  - 4.81160e-02	        627257395774  - 1.17191e+00	       6014144044649  - 6.00989e-01	        627257511659  - 8.76523e-01	       6014145991940  + 3.53750e-04	
//22	       1226003424937  - 1.47022e+00	       1226005444762  + 3.74120e-01	       1226004792676  + 5.16664e-04	       9142713816822  - 9.84294e-01	       1226005277753  + 2.12356e-01	       9142718404920  + 2.75905e-01	
//23	       2238791389915  + 2.53479e-03	       2238790104123  - 6.54467e-01	       2238792211428  + 3.59270e-01	      12985309342107  + 1.50901e-01	       2238791914699  + 1.60863e-01	      12985310080750  + 3.52169e-01	
//24	       3824598409276  - 3.05778e+00	       3824605220060  + 3.00659e+00	       3824601957643  + 4.32422e-03	      17253630089069  + 6.43554e-01	       3824602639675  + 1.71816e-01	      17253630854289  + 9.73068e-01	
//25	       6119359680403  - 1.72208e+00	       6119358689170  - 2.93431e+00	       6119363800347  + 1.24747e-01	      21471179862209  - 4.82333e-04	       6119361863273  - 1.84780e-01	      21471186819905  + 2.18916e+00	
//26	       9179044800184  + 1.83468e-02	       9179042627170  - 3.38478e-01	       9179046128893  + 3.29490e-01	      25049706121827  - 5.87500e-01	       9179044544141  + 2.59478e-03	      25049703406899  - 1.71330e+00	
//27	      12918659310311  + 1.39382e+00	      12918650957235  - 1.30738e+00	      12918656509380  + 1.61058e-01	      27420670097521  - 1.41965e-04	      12918653624089  - 1.61147e-01	      27420666950652  - 3.75606e-01	
//28	      17071085613780  + 1.90577e+00	      17071082964084  + 5.46403e-01	      17071081487264  + 1.45738e-01	      28182352566485  - 2.93418e-01	      17071080709393  + 3.74375e-02	      28182353727384  - 1.04331e-01	
//29	      21191687118816  + 1.38517e-01	      21191691544110  + 1.77817e+00	      21191678220998  - 2.43573e+00	      27210545425782  - 7.96763e-01	      21191677361720  - 3.05321e+00	      27210545586429  - 7.42732e-01	
//30	      24723630420601  - 2.63513e-01	      24723639853758  + 1.91493e+00	      24723639513080  + 1.73000e+00	      24691061500791  + 1.87096e+00	      24723643509771  + 4.49054e+00	      24691056897065  + 1.94784e-01	
//31	      27116243007705  + 5.66270e-03	      27116242008692  - 1.35948e-02	      27116237611983  - 9.23383e-01	      21062620770530  + 1.87705e-02	      27116239819473  - 2.88378e-01	      21062620312793  + 1.38887e-03	
//32	      27963622063766  - 3.51180e-01	      27963619305861  - 1.24130e+00	      27963618449052  - 1.62860e+00	      16893977036140  + 1.27467e-02	      27963626758344  + 8.71218e-02	      16893974665143  - 2.15251e-01	
//33	      27116239643870  - 3.25733e-01	      27116242339072  - 2.82508e-03	      27116243068007  + 7.53964e-03	      12741583916432  - 9.63070e-02	      27116235516031  - 1.85894e+00	      12741585076523  + 2.15036e-04	
//34	      24723630362256  - 2.75698e-01	      24723638678962  + 1.31685e+00	      24723640310934  + 2.17786e+00	       9035702611645  + 1.30037e+00	      24723633560230  + 1.39453e-02	       9035702708416  + 1.37483e+00	
//35	      21191687268183  + 1.63722e-01	      21191682109294  - 5.12704e-01	      21191686675907  + 7.61574e-02	       6023799068866  - 2.48554e-02	      21191688398213  + 4.22631e-01	       6023801710376  + 8.43832e-01	
//36	      17071080315186  + 9.61923e-03	      17071078042081  - 2.04378e-01	      17071078427738  - 1.28696e-01	       3774173960883  + 2.91009e+00	      17071078983174  - 5.03147e-02	       3774169494667  - 3.51707e-01	
//37	      12918652842341  - 3.83075e-01	      12918655139445  + 4.07002e-04	      12918658644818  + 9.90913e-01	       2221434478666  + 2.23518e-01	      12918659878035  + 1.79173e+00	       2221432584643  - 6.36801e-01	
//38	       9179050049923  + 3.49022e+00	       9179045022126  + 4.35581e-02	       9179041095589  - 1.18225e+00	       1227635126645  + 3.68491e-01	       9179041204889  - 1.10510e+00	       1227630612649  - 1.20202e+01	
//39	       6119362548585  - 2.33556e-02	       6119360426524  - 1.02144e+00	       6119362202805  - 8.56183e-02	        636552165326  + 1.46878e+00	       6119364851310  + 6.05353e-01	        636550609885  - 5.44097e-01	
//40	       3824602369118  + 7.62649e-02	       3824600698435  - 3.34223e-01	       3824600753464  - 3.02480e-01	        309434838207  + 1.67821e-01	       3824605535247  + 3.59147e+00	        309435649255  + 3.48821e+00	
//41	       2238791293801  - 1.92919e-04	       2238789699201  - 1.16557e+00	       2238792330967  + 4.61426e-01	        140880648934  - 1.58153e-01	       2238788885493  - 2.63556e+00	        140880643001  - 1.70975e-01	
//42	       1226004841734  + 4.49387e-03	       1226005591088  + 5.53248e-01	       1226004505962  - 5.57961e-02	         60004683788  - 1.68764e-01	       1226004858778  + 6.79461e-03	         60004882265  + 1.59550e-01	
//43	        627258600429  + 1.92272e-01	        627259014029  + 9.22969e-01	        627258971687  + 8.23103e-01	         23877717186  - 1.01506e+00	        627257968168  - 1.29474e-01	         23877852151  - 1.79776e-02	
//44	        299372841139  - 5.78214e-01	        299373467429  + 1.47638e-01	        299373199604  - 1.10785e-02	          8863756244  + 1.86021e-05	        299372829834  - 6.10063e-01	          8863725973  - 1.00625e-01	
//45	        133054492347  - 6.26080e-01	        133054044265  - 4.07902e+00	        133054825160  + 1.46767e-02	          3064064477  + 8.19625e-01	        133054951486  + 2.18526e-01	          3064051703  + 4.55031e-01	
//46	         54957512627  + 1.93403e-01	         54957451346  + 3.18164e-02	         54957679740  + 1.32854e+00	           984292096  - 1.70949e+00	         54957625504  + 8.48741e-01	           984379927  + 2.22606e+00	
//47	         21047534554  + 1.21786e-02	         21047555646  + 6.54035e-02	         21047293152  - 2.41365e+00	           293233239  + 2.60364e+00	         21047624322  + 5.31609e-01	           293205416  - 1.27413e-04	
//48	          7454177750  - 3.08858e+00	          7454389308  + 4.80108e-01	          7454419763  + 1.09336e+00	            80759869  - 6.81096e-01	          7454379407  + 3.34340e-01	            80764466  - 9.84530e-02	
//49	          2434137829  + 2.07444e+00	          2434094277  + 3.10844e-01	          2434112689  + 8.66254e-01	            20513197  + 3.69358e-02	          2434049840  - 1.17761e-01	            20514978  + 3.42723e-01	
//50	           730177105  - 2.52342e+00	           730249661  + 1.20228e+00	           730254152  + 1.59436e+00	             4791852  + 6.65191e+00	           730222598  + 9.02312e-03	             4788283  + 8.98260e-01	
//51	           200413063  - 7.78149e+00	           200454132  + 1.23663e-02	           200478418  + 3.33626e+00	             1022231  + 1.12305e-01	           200462344  + 4.77791e-01	             1022055  + 2.59259e-02	
//52	            50121046  + 1.24747e+00	            50108610  - 4.09381e-01	            50113368  + 1.04289e-03	              198594  - 5.79074e-02	            50094608  - 6.85274e+00	              198027  - 2.28804e+00	
//53	            11344517  - 3.03004e-01	            11348268  + 3.17098e-01	            11346213  - 2.20526e-03	               35027  + 3.61823e-02	            11349808  + 1.04101e+00	               35084  + 2.44957e-01	
//54	             2313927  + 2.99075e+00	             2310378  - 3.66069e-01	             2310768  - 1.21457e-01	                5613  + 8.61077e-01	             2311215  - 2.96863e-03	                5627  + 1.24539e+00	
//55	              420242  + 8.65326e-05	              420138  - 2.28397e-02	              421885  + 6.47089e+00	                 760  - 7.33908e-01	              419976  - 1.60825e-01	                 812  + 1.00095e+00	
//56	               67764  + 7.56767e-01	               67545  + 7.41419e-04	               67502  - 1.91080e-02	                 102  + 1.63401e-01	               67313  - 7.49071e-01	                  78  - 4.08104e+00	
//57	                9376  - 1.11936e+00	                9491  + 1.51742e-02	                9474  - 2.64462e-03	                  12  + 1.58543e-01	                9500  + 4.64936e-02	                   6  - 2.06289e+00	
//58	                1202  + 2.93868e+00	                1148  + 1.38597e-02	                1164  + 3.49013e-01	                   2  + 9.87505e-01	                1175  + 8.39043e-01	                   0  - 1.00419e+00	
//59	                 123  + 3.81162e-01	                 101  - 2.02286e+00	                 121  + 1.86589e-01	                   0  - 7.94273e-02	                 115  - 1.54529e-02	                   0  - 7.94273e-02	
//60	                  10  + 9.59078e-03	                  15  + 2.90274e+00	                  11  + 1.75640e-01	                   0  - 5.14807e-03	                   6  - 1.40830e+00	                   0  - 5.14807e-03	
//61	                   3  + 8.79242e+00	                   1  + 2.08707e-01	                   0  - 6.35742e-01	                   0  - 2.62561e-04	                   2  + 2.92760e+00	                   0  - 2.62561e-04	
//
//              50.369 57 =>     3.740456e-01           51.038 57 =>     2.532920e-01           51.300 57 =>     3.440578e-01           57.807 56 =>     6.306025e-01           66.885 57 =>     8.154200e-01           58.033 56 =>     6.655878e-01
//
//--- Finished -- BBin trials: 281474976710656	2^48.0 (out of 2^48) trials -- Slash256 -- 	ps:   3.740e-01    2.533e-01    3.441e-01    6.306e-01    8.154e-01*   6.656e-01   => p <   7.62337e-01   2^52.58 calls, 2^55.58 bytes	2^36.59 bytes/second	used:   6::00:59:24.89
//
//// I am very confident in this generator now.
	//const uint64_t fa = stateA;
	//const uint64_t fb = stateB;
	//const uint64_t fc = stateC;
	//const uint64_t fd = stateD;
	// stateA = fd * 0xF1357AEA2E62A9C5UL;
	// stateB = rotate64(fa, 44);
	// stateC = fb + 0x9E3779B97F4A7C15UL;
	// return (stateD = fa ^ fc);


	// stateA = fd ^ fb + fc;
	// stateB = rotate64(fa, 37);
	// stateC = fc + 0x9E3779B97F4A7C15UL;
	// return (stateD = fa + fb);

	// SquawkRandom
	// Passes 64TB with no anomalies (seed 0).
	// Guaranteed minimum period of 2 to the 64.
	// stateA = fd * 0xF1357AEA2E62A9C5UL;
	// stateB = rotate64(fa, 44);
	// stateC = fc + 0x9E3779B97F4A7C15UL;
	// return (stateD = fa + fb ^ fc);

	//stateA = fd * 0xF1357AEA2E62A9C5UL;
	//stateB = rotate64(fa, 10);
	//stateC = fc + 0x9E3779B97F4A7C15UL;//0xDE916ABCC965815BUL;
	//return (stateD = fc + fb);

// has just one anomaly, at 64TB:
// [Low4/32]Gap-16:B                 R=  +4.9  p =  3.8e-4   unusual
//uint64_t fa = (stateA += 0xDB4F0B9175AE2165L);
//uint64_t fb = (stateB += 0xBBE0563303A4615FL);
//uint64_t fc = (stateC += 0xA0F2EC75A1FE1575L);
//fb += fa ^ rotate64(fa, 7) ^ rotate64(fa, 46);
//fc += fb ^ rotate64(fb, 17) ^ rotate64(fb, 34);
//fa += fc ^ rotate64(fc, 5) ^ rotate64(fc, 41);
//return fa;

// StaunchRandom
// Passes 64TB with no anomalies.
// Period is 2 to the 64; there are 2 to the 128 streams possible.
// How correlated those streams are is unknown.
//uint64_t fa = (stateA += 0xD1B54A32D192ED03L);
//uint64_t fb = (stateB += 0xABC98388FB8FAC03L);
//uint64_t fc = (stateC += 0x8CB92BA72F3D8DD7L);
//fb += fa ^ rotate64(fa, 29) ^ rotate64(fa, 53);
//fc += fb ^ rotate64(fb,  7) ^ rotate64(fb, 19);
//fa += fc ^ rotate64(fc, 37) ^ rotate64(fc, 47);
//return fa;

// SkinkRandom
// Passes 64TB with no anomalies.
// Period is 2 to the 64; there are 2 to the 128 streams possible.
// How correlated those streams are is unknown.
//uint64_t fa = (stateA += 0xD1B54A32D192ED03L);
//uint64_t fb = (stateB += 0xABC98388FB8FAC03L);
//uint64_t fc = (stateC += 0x8CB92BA72F3D8DD7L);
//fb += rotate64(fa, 31) ^ fc;
//fc += rotate64(fb, 19) ^ fa;
//fa += rotate64(fc, 47) ^ fb;
//return fa ^ rotate64(fa, 29) ^ rotate64(fa, 53);

// SquashRandom
// Passes 64TB with no anomalies.
// Period is 2 to the 64, there arre 2 to the 192 streams possible.
// How correlated those streams are is unknown.
// Probably slower than SkinkRandom, above.
				// uint64_t fa = (stateA += 0xDB4F0B9175AE2165L);
				// uint64_t fb = (stateB += 0xBBE0563303A4615FL);
				// uint64_t fc = (stateC += 0xA0F2EC75A1FE1575L);
				// uint64_t fd = (stateD += 0x89E182857D9ED689L);
				// fb += rotate64(fa, 25) ^ fc;
				// fc += rotate64(fb, 46) ^ fd;
				// fd += rotate64(fc, 37) ^ fa;
				// fa += rotate64(fd, 18) ^ fb;
				// return fa ^ rotate64(fa, 13) ^ rotate64(fa, 53);

// just awful.
//				uint64_t fa = (stateA += 0xDB4F0B9175AE2165L);
//				uint64_t fb = (stateB += 0xBBE0563303A4615FL);
//				uint64_t fc = (stateC += 0xA0F2EC75A1FE1575L);
//				uint64_t fd = (stateD += 0x89E182857D9ED689L);
//				uint64_t fx = rotate64(fa, 25) + rotate64(fb, 46);
//				uint64_t fy = rotate64(fc, 18) + rotate64(fd, 37);
//				uint64_t fz = rotate64(fa, 14) + rotate64(fc, 53);
//				uint64_t fw = rotate64(fb, 22) + rotate64(fd, 41);
//				return fx ^ fy ^ fz ^ fw;

// PouchRandom
// Passes 64TB without any anomalies. Passes 179PB of ReMort without suspicion.
// Minimum period is 2 to the 63.
// Multiplies a variable by (an odd) variable.
// 3 other (quicker) operations; 4 math/bit ops total.
// 7	                4411  - 3.86292e-03	                4471  + 7.06996e-01	                4530  + 2.98862e+00	                4357  - 7.65340e-01	                4415  - 3.81636e-06	
// 8	              360979  + 4.38362e-02	              361818  + 2.57940e+00	              361696  + 1.96829e+00	              361293  + 5.35949e-01	              361123  + 2.01679e-01	
// 9	            22352038  - 1.07697e+00	            22351139  - 1.50774e+00	            22365460  + 3.24315e+00	            22366036  + 3.69675e+00	            22351884  - 1.14563e+00	
//11	                4452  + 3.07898e-01	                4497  + 1.51813e+00	                4336  - 1.41820e+00	                4535  + 3.25446e+00	                4357  - 7.65340e-01	
//12	            85016119  - 9.06914e-02	            85014340  - 2.44123e-01	            85001483  - 3.56632e+00	            85012695  - 4.52248e-01	            85015464  - 1.38523e-01	
//13	          6948602685  - 9.57907e-01	          6948658091  - 9.86326e-02	          6948570988  - 1.84681e+00	          6948552119  - 2.51328e+00	          6948659853  - 8.58025e-02	
//14	        430511647206  + 4.94655e-01	        430511593534  + 3.86283e-01	        430512001617  + 1.54621e+00	        430512009075  + 1.57461e+00	        430511593152  + 3.85560e-01	
//16	              361090  + 1.55356e-01	              361483  + 1.09909e+00	              360507  - 3.32197e-01	              361778  + 2.36994e+00	              361050  + 1.07298e-01	
//17	          6948742604  + 4.89704e-01	          6948677018  - 7.56962e-03	          6948549950  - 2.59646e+00	          6948697312  + 2.44766e-02	          6948632635  - 3.83702e-01	
//18	        567923726307  + 2.98487e-01	        567923407380  + 1.51632e-02	        567922566389  - 9.85683e-01	        567922561654  - 9.98198e-01	        567924300246  + 1.71068e+00	
//19	      35186136998801  + 1.72835e-01	      35186137382921  + 2.30871e-01	      35186134335031  - 1.11107e-03	      35186134191133  - 3.31677e-03	      35186136558074  + 1.16578e-01	
//21	            22366120  + 3.76538e+00	            22351442  - 1.35447e+00	            22352579  - 8.52579e-01	            22365349  + 3.15915e+00	            22366255  + 3.87701e+00	
//22	        430511817020  + 9.25689e-01	        430511576997  + 3.55589e-01	        430511716863  + 6.55258e-01	        430511865790  + 1.07424e+00	        430511927640  + 1.27853e+00	
//23	      35186133147045  - 5.45723e-02	      35186137424716  + 2.37692e-01	      35186138352932  + 4.14759e-01	      35186134361950  - 8.29127e-04	      35186132515794  - 1.15617e-01	
//24	    2179984568538371  - 3.01410e-03	    2179984564515401  - 1.98990e-02	    2179984567140887  - 7.19642e-03	    2179984570970172  - 7.93647e-06	    2179984569033306  - 1.96253e-03	
//
//               8.841 15 =>     1.589194e-01           10.362 15 =>     2.642660e-01           22.423 15 =>   1-7.055889e-02           20.423 15 =>     8.828091e-01           10.314 15 =>     2.611097e-01
//
//--- Finished -- ReMort trials: 2251799813685248	2^51.0 (out of 2^51) trials -- Pouch256 -- 64 bits: 	ps:   1.589e-01    2.643e-01  1-7.056e-02*   8.828e-01    2.611e-01   => p <   5.76370e-01   2^54.32 calls, 2^57.32 bytes	2^37.57 bytes/second	used:  10::03:48:02.88

	// const uint64_t a = stateA;
	// const uint64_t b = stateB;
	// const uint64_t c = stateC;
	// const uint64_t d = stateD;
//stateA = c * d;
//stateB = rotate64(a, 47);
//stateC = b - a;
//stateD = d + 0xE35E156A2314DCDAL;// 0x9E3779B97F4A7C16L;
//return c;

// // Alternate PouchRandom
// // Also passes 64TB with no anomalies.
// // Different rotation, that's all.
// stateA = c * d;
// stateB = rotate64(a, 52);
// stateC = b - a;
// stateD = d + 0xE35E156A2314DCDAL;// 0x9E3779B97F4A7C16L;
// return c;

// SnoutRandom
// rot 7, fails 8GB,
//  [Low8/64]DC6-9x1Bytes-1           R= +27.3  p =  2.6e-14    FAIL           
// rot 8, fails 16GB,
//  [Low1/8]DC6-9x1Bytes-1            R= +37.6  p =  3.5e-20    FAIL !!        
//  [Low1/8]BRank(12):4K(1)           R=+152.8  p~=  5.1e-47    FAIL !!! 
// rot 9, fails 32GB,
//  BCFN(2+0,13-0,T)                  R= +27.5  p =  3.1e-14    FAIL           
// rot 10, fails 128,
//  BCFN(2+0,13-0,T)                  R= +23.8  p =  2.8e-12    FAIL           
// rot 11, fails 256GB,
//  BCFN(2+0,13-0,T)                  R= +33.5  p =  1.9e-17    FAIL !         
//  BCFN(2+1,13-0,T)                  R= +27.0  p =  5.9e-14    FAIL
// rot 12, fails 128GB,
//  BCFN(2+0,13-0,T)                  R= +25.4  p =  3.8e-13    FAIL           
// rot 13, fails 128GB,
//  BCFN(2+0,13-0,T)                  R= +39.3  p =  1.4e-20    FAIL !!    
// rot 17, fails 128GB,
//  [Low4/16]BCFN(2+1,13-0,T)         R= +32.2  p =  9.7e-17    FAIL !         
// rot 20, fails 128GB,
//  BRank(12):4K(1)                   R=+540.3  p~=  1.1e-163   FAIL !!!!!     
// rot 23, fails 64GB,
//  BCFN(2+0,13-0,T)                  R= +23.1  p =  7.0e-12    FAIL
// rot 42, fails 256 GB,
//  BCFN(2+0,13-0,T)                  R= +44.5  p =  2.5e-23    FAIL !!        
// rot 53, fails 2GB,
//  BRank(12):4K(1)                   R=+540.3  p~=  1.1e-163   FAIL !!!!!     

//	const uint64_t a = stateA;
//	const uint64_t b = stateB;
//	const uint64_t c = stateC;
//	const uint64_t d = stateD;
//		stateA = a + 0xD1B54A32D192ED03UL;
//		stateB = a ^ d;
//		stateC = rotate64(b, 42);
//		stateD = c + b;
//		return d;

// FloodRandom
// Passes 64TB of PractRand with no anomalies.
// Period is 2 to the 64. There are 2 to the 128 possible streams.
					// uint64_t fa = (stateA += 0xD1B54A32D192ED03L);
					// uint64_t fb = (stateB += 0xABC98388FB8FAC03L);
					// uint64_t fc = (stateC += 0x8CB92BA72F3D8DD7L);
					// uint64_t x = fa ^ rotate64(fb, 19) ^ rotate64(fc, 47);
					// x = (x ^ x >> 27) * 0x3C79AC492BA7B653UL;
					// x = (x ^ x >> 33) * 0x1C69B3F74AC4AE35UL;
					// x ^= x >> 27;
					// return x;

//fb ^= (fc += rotate64(fa, 31));
//fc ^= (fa += rotate64(fb, 19));
//fa ^= (fb += rotate64(fc, 47));
//return fa;
// return (fa ^ rotate64(fa, 29) ^ rotate64(fa, 50)) +
//        (fb ^ rotate64(fb, 19) ^ rotate64(fb, 37)) +
//        (fc ^ rotate64(fc, 13) ^ rotate64(fc, 44));

					// FlowRandom
					// Period is 2 to the 64, has 2 to the 64 streams.
					// Passes 64TB of PractRand with no anomalies.
					// Uses one more numerical operation than DistinctRandom/SplitMix64.
					// Based on Moremur for its constants.
					//uint64_t x = (stateA += 0xD1B54A32D192ED03L);
					//uint64_t y = (stateB += 0x8CB92BA72F3D8DD7L);
					//x = (x ^ rotate64(y, 37)) * 0x3C79AC492BA7B653UL;
					//x = (x ^ x >> 33) * 0x1C69B3F74AC4AE35UL;
					//x ^= x >> 27;
					//return x;

// OrbitalRandom
// Period is 2 to the 128, no streams.
// Passes 64TB of PractRand with no anomalies.
// Uses __lzcnt64 in C++, and would use some other technique to count leading zeros elsewhere.
// In most ways, very similar to FlowRandom, but stateB uses a different increment and adds __lzcnt64(stateA).
// Because stateA and stateB are different in where entropy occurs, this doesn't do (x ^ rotate64(y, 37)), and
// instead does (y ^ rotate64(x, 37)) .
					//uint64_t x = (stateA += 0xD1B54A32D192ED03L);
					//uint64_t y = (stateB += 0x9E3779B97F4A7C15L + __lzcnt64(x));
					//x = (y ^ rotate64(x, 37)) * 0x3C79AC492BA7B653UL;
					//x = (x ^ x >> 33) * 0x1C69B3F74AC4AE35UL;
					//x ^= x >> 27;
					//return x;

// Passes 32TB, but has real trouble with ICE tests.
					//uint64_t x = (stateA += 0xD1B54A32D192ED03UL);
					//uint64_t y = (stateB += 0x9E3779B97F4A7C15UL + __lzcnt64(x));
					//x = (x ^ x >> 32) * 0xBEA225F9EB34556DUL;
					//y = (y ^ y >> 29) * 0xF1357AEA2E62A9C5UL;
					//x ^= x >> 29;
					//y ^= y >> 32;
					//return (y ^ rotate64(x, (int)y & 63));
					
// ThrushRandom
// Unlike ThrashRandom, this passes 64TB with no anomalies.
// It is also faster than SoloRandom, somehow, and has a minimum period of 2 to the 64.
// No multiplication used, either. Incorporating stateD with XOR helps quality, but not period.
//uint64_t fd = stateD;
//stateD ^= (stateA = (stateB = rotate64(stateB, 41) ^ (stateC += 0xBEA225F9EB34556DUL)) + rotate64(stateA, 26));
//return fd;

// ThrooshRandom
// A slightly faster variant on ThrushRandom; also passes 64TB with no anomalies.
// Whoosh!
uint64_t fa = stateA;
uint64_t fb = stateB;
uint64_t fc = stateC;
uint64_t fd = stateD;
stateA = rotate64(fa, 26) + fb;
stateB = rotate64(fb, 41) ^ fc;
stateC = fc + 0xBEA225F9EB34556DUL;
stateD = fd ^ fa;
return fd;
				}
				std::string lizard256::get_name() const { return "lizard256"; }
				void lizard256::walk_state(StateWalkingObject *walker) {
					walker->handle(stateA);
					walker->handle(stateB);
					walker->handle(stateC);
					walker->handle(stateD);
					// stateC |= 1UL;
					// stateD |= 1UL;
				}
				Uint64 plum256::raw64() {
					//const uint64_t fa = stateA;
					//const uint64_t fb = stateB;
					//const uint64_t fc = stateC;
					//const uint64_t fd = stateD;
					//// Good rotations: 23, 25, 29, 46
					//stateA = fd * 0xD1342543DE82EF95UL;
					//stateB = rotate64(fa + fd, rotation);
					//stateC = fc + 0xC6BC279692B5C323UL;
					//stateD = fb ^ fc;
					//return fd;
					//// Good rotations: most; 6, 8-13, 15-18, 20, 26,, 30-31, 33-37, 40, 42-44, 46-48, 51-56
//					stateA = 0xD1342543DE82EF95UL * fc;
//					stateB = fa ^ fc;
//					stateC = rotate64(fb, rotation) + fd;
//					stateD = fd + 0xC6BC279692B5C323UL;
//					return fa;

	const uint64_t fa = stateA;
	const uint64_t fb = stateB;
	const uint64_t fc = stateC;
	const uint64_t fd = stateD;
    //// Good rotations: 23, 41, 53-55
//	stateA = fb ^ fc;
//	stateB = rotate64(fa + fd, rotation);
//	stateC = fa + fb;
//	stateD = fd + 0x9E3779B97F4A7C15UL;
//	return fc;
	//// No good rotations.
//	stateA = fc ^ fd;
//	stateB = rotate64(fa, rotation) + fc;
//	stateC = fa + fb;
//	stateD = fd + 0x9E3779B97F4A7C15UL;
//	return fc;
    //// Good rotations: 11-13, 19, 22-23, 26, 37-39, 41-42, 45-46, 50, 53-55
//	stateA = rotate64(fc + fb, rotation);
//	stateB = fc ^ fd;
//	stateC = fa + fb;
//	stateD = fd + 0x9E3779B97F4A7C15UL;
//	return fc;
			//// Good rotations: 9-13, 15, 18-21, 23, 25-27, 30, 36-37, 42, 44, 47, 49, 51, 53, 55, 57, 59
//			stateA = rotate64(fb * 0xF1357AEA2E62A9C5UL, rotation);
//			stateB = fc ^ fd;
//			stateC = fa + fb;
//			stateD = fd + 0x9E3779B97F4A7C15UL;
//			return fc;

			//// Good rotations: 2, 4-6, 9-10, 13, 16, 18-21, 23, 27, 30, 33, 36-37, 40, 43-44, 47, 49-50, 53-56, 58, 61
			//stateA = rotate64(fb + fc, rotation);
			//stateB = rotate64(fc ^ fd, 63 ^ rotation);
			//stateC = fa + fb;
			//stateD = fd + 0x9E3779B97F4A7C15UL;
			//return fc;

			//stateA = rotate64(fb + fc, rotation);
			//stateB = rotate64(fc ^ fd, 62 - rotation);
			//stateC = fa + fb;
			//stateD = fd + 0xC6BC279692B5C323UL;
			//return fc;
			
//stateA = rotate64(fc, rotation); // try 53, 44, 42, 23, 22, 21, 20 12, 11, 9
//stateB = fa + fc;
//stateC = fb ^ fd;
//stateD = fb + 0xC6BC279692B5C323UL;//0x9E3779B97F4A7C15UL;//0xDE916ABCC965815BUL;
//return stateB;

//	stateA = fd * 0xF1357AEA2E62A9C5UL;
//	stateB = fb + 0xDE916ABCC965815BUL;
//	stateC = fa ^ fb;
//	stateD = rotate64(fc, rotation);
//	return fc;

//stateA = fb + rotate64(fc, rotation); // 21, 19
//stateB = fa ^ fd;
//stateC = fb ^ fd;
//stateD = fd + 0xC6BC279692B5C323UL;//0x9E3779B97F4A7C15UL;//0xDE916ABCC965815BUL;
//return stateA;

//stateA = fb + rotate64(fc, rotation); // 42, 39, 36, 28, 25, 23, 22, 20, 19, 13, 11
//stateB = fa ^ fc;
//stateC = fa ^ fd;
//stateD = fd + 0xDE916ABCC965815BUL;
//return stateA;

//rng=plum256x10, seed=0x0
//length= 1 terabyte (2^40 bytes), time= 4592 seconds
//  no anomalies in 988 test result(s)
//
//rng=plum256x10, seed=0x0
//length= 2 terabytes (2^41 bytes), time= 9253 seconds
//  Test Name                         Raw       Processed     Evaluation
//  BCFN(2+2,13-0,T)                  R=  +9.3  p =  1.6e-4   unusual
//  ...and 1019 test result(s) without anomalies
//
//rng=plum256x10, seed=0x0
//length= 4 terabytes (2^42 bytes), time= 19046 seconds
//  Test Name                         Raw       Processed     Evaluation
//  BCFN(2+2,13-0,T)                  R= +22.2  p =  2.1e-11   VERY SUSPICIOUS
//  ...and 1051 test result(s) without anomalies
	//stateA = rotate64(fc, rotation);
	//stateB = fa + fc;
	//stateC = fb ^ fd;
	//stateD = fd + 0xDE916ABCC965815BUL;
	//return fc + fb;


//rng=plum256x10, seed=0x0
//length= 1 terabyte (2^40 bytes), time= 4671 seconds
//  Test Name                         Raw       Processed     Evaluation
//  BCFN(2+0,13-0,T)                  R= +42.6  p =  2.7e-22    FAIL !!
//  ...and 987 test result(s) without anomalies
//	stateA = rotate64(fc, rotation);
//	stateB = fa + fc;
//	stateC = fb ^ fd;
//	stateD = fd + 0xDE916ABCC965815BUL;
//	return fa + fb + fc;

	// We can't use too many XORs for states...
//rng=plum256x10, seed=0x0
//length= 8 gigabytes (2^33 bytes), time= 18.2 seconds
//  Test Name                         Raw       Processed     Evaluation
//  BCFN(2+0,13-0,T)                  R= +24.3  p =  1.6e-12    FAIL
//  ...and 729 test result(s) without anomalies
	//stateA = rotate64(fc, rotation);
	//stateB = fa ^ fc;
	//stateC = fb ^ fd;
	//stateD = fd + 0xDE916ABCC965815BUL;
	//return fa + fb + fc;
	stateD = fd * 0xDE916ABCC965815BUL + 0x9E3779B97F4A7C15UL;
	stateA = rotate64(fc, rotation);
	stateB = fa + fc;
	stateC = fb ^ fd;
	return fa + fb + fc;

				}

				std::string plum256::get_name() const {
					std::ostringstream tmp;
					tmp << "plum256x" << rotation;
					return tmp.str();
				}
				void plum256::walk_state(StateWalkingObject *walker) {
					walker->handle(stateA);
					walker->handle(stateB);
					walker->handle(stateC);
					walker->handle(stateD);
				}
				
				Uint64 overload320::raw64() {
//					//// Overload320; XORs an unscrambled xoshiro with a Weyl sequence, then runs that through MX3.
//					//// Passes 64TB without anomalies.
//					uint64_t x = stateB ^ (stateE += stream);
//					const uint64_t t = stateB << 17;
//
//					stateC ^= stateA;
//					stateD ^= stateB;
//					stateB ^= stateC;
//					stateA ^= stateD;
//
//					stateC ^= t;
//
//					stateD = rotate64(stateD, 45);
//					x ^= x >> 32;
//					x *= 0xbea225f9eb34556dUL;
//					x ^= x >> 29 ^ stream;
//					x *= 0xbea225f9eb34556dUL;
//					x ^= x >> 32;
//					x *= 0xbea225f9eb34556dUL;
//					x ^= x >> 29;
//
//					return x;

					//const uint64_t x = rotate64(stateB + (stateE += stream), 29);

					//// Passes 64TB with no anomalies.
//					const uint64_t x = stateB * 5UL;
//					const uint64_t t = stateB << 17;
//
//					stateC ^= stateA;
//					stateD ^= stateB;
//					stateB ^= stateC;
//					stateA ^= stateD;
//
//					stateC ^= t;
//
//					stateD = rotate64(stateD, 45);
//					return rotate64(x, 17) * ((stateE += stream) | 1UL);

//					//// Passes 64TB with no anomalies.
//					//// Also passes 64TB with no anomalies when stream is forced to 1 (using seed=1).
//					const uint64_t t = stateB << 17;
//					const uint64_t x = t - stateB;
//
//					stateC ^= stateA;
//					stateD ^= stateB;
//					stateB ^= stateC;
//					stateA ^= stateD;
//
//					stateC ^= t;
//
//					stateD = rotate64(stateD, 45);
//					return rotate64(x, 7) * ((stateE += stream) | 1UL);

    //// Pasar320
	//// Passes 64TB of PractRand with no anomalies.
	//// Fails ReMort rather quickly, but it also runs very quickly (over twice as fast as the generator immediately above).
	//// Various simple changes to this so far have all failed ReMort after different spans of time.
//	const uint64_t fa = stateA;
//	const uint64_t fb = stateB;
//	const uint64_t fc = stateC;
//	const uint64_t fd = stateD;
//	const uint64_t fe = stateE;
//	stateA = rotate64(fc, 41);
//	stateB = fd ^ fc;
//	stateC = fe ^ fb;
//	stateD = fa + fc;
//	stateE = fe + 0xDE916ABCC965815BUL;
//	return stateD;

	const uint64_t fa = stateA;
	const uint64_t fb = stateB;
	const uint64_t fc = stateC;
	const uint64_t fd = stateD;
	const uint64_t fe = stateE;

//	stateA = fd + fc ^ fe;//0xF1357AEA2E62A9C5UL;//0xD1342543DE82EF95UL
//	stateB = fb + 0xDE916ABCC965815BUL;
//	stateC = fa + fe;
//	stateD = rotate64(fe, 42);
//	stateE = fb ^ fd - fa;
//	return fa;

//// Pasar320.
//// Passes 64TB with no anomalies (using seed 3).
//// Passes 179PB of ReMort, 107 PB of BBin.
//  stateA = fe * 0xF1357AEA2E62A9C5UL;
//  stateB = rotate64(fa, 44);
//  stateC = fb + fd;
//  stateD = fd + 0x9E3779B97F4A7C15UL;
//  return stateE = fa ^ fc;

  // Finch320.
  // Passes 64TB with one anomaly at 128GB,
  // [Low4/16]BCFN(2+0,13-0,T)         R=  +8.6  p =  4.1e-4   unusual
//   stateA = rotate64(fe, 50);
//   stateB = fb ^ fa + fd;
//   stateC = rotate64(fb, 25);
//   stateD = fd + 0xF1357AEA2E62A9C5UL;
//   return stateE = fa + fc;

//  // Lantern320.
//  // Passes 64TB with no anomalies.
//  stateA = fa + 0x9E3779B97F4A7C15UL;
//  stateB = rotate64(fe, 41);
//  stateC = fa ^ fb;
//  stateD = rotate64(fc, 17);
//  return stateE = fc + fd;

//  stateA = fa + 0x9E3779B97F4A7C15UL;
//  stateB = rotate64(fe, 49);
//  stateC = fa ^ fb;
//  stateD = rotate64(fc, 25);
//  return stateE = fc + fd;

//  // Ace320, or AceRandom.
//  // Passes 64TB with no anomalies.
//  // Passes over 179PB of ReMort testing without suspicion.

//  7	                4375  - 3.64746e-01	                4371  - 4.41083e-01	                4399  - 5.89271e-02	                4352  - 9.02663e-01	                4500  + 1.63142e+00	
//  8	              359692  - 3.73684e+00	              360668  - 9.50793e-02	              361641  + 1.71977e+00	              360716  - 5.21866e-02	              361894  + 3.00179e+00	
//  9	            22353840  - 4.31204e-01	            22352868  - 7.43442e-01	            22356386  - 1.39718e-02	            22357358  + 7.63311e-03	            22351500  - 1.32607e+00	
// 11	                4392  - 1.21172e-01	                4415  - 3.81636e-06	                4452  + 3.07898e-01	                4435  + 8.94254e-02	                4430  + 5.00829e-02	
// 12	            85004296  - 2.50713e+00	            85027838  + 9.40536e-01	            85013786  - 3.07106e-01	            85014750  - 2.02160e-01	            85009759  - 9.81907e-01	
// 13	          6948758494  + 7.92830e-01	          6948685953  + 4.07383e-04	          6948782218  + 1.38065e+00	          6948724468  + 2.32539e-01	          6948757464  + 7.70979e-01	
// 14	        430511154429  - 2.27668e-03	        430511203405  + 7.25159e-04	        430511083853  - 2.41113e-02	        430511140656  - 4.72048e-03	        430511150850  - 2.82697e-03	
// 16	              361549  + 1.34154e+00	              360430  - 4.96386e-01	              360437  - 4.80102e-01	              361827  + 2.62775e+00	              361452  + 9.93554e-01	
// 17	          6948652253  - 1.47527e-01	          6948576663  - 1.66641e+00	          6948694372  + 1.46848e-02	          6948670385  - 2.77473e-02	          6948545760  - 2.76098e+00	
// 18	        567924722199  + 3.48883e+00	        567922977586  - 1.99967e-01	        567927186798  + 2.64016e+01	        567923911076  + 6.26503e-01	        567926883023  + 2.24216e+01	
// 19	      35186138463734  + 4.39168e-01	      35186140285056  + 9.40398e-01	      35186144862752  + 3.03270e+00	      35186148161071  + 5.27853e+00	      35186136414154  + 1.00598e-01	
// 21	            22356486  - 9.41934e-03	            22353049  - 6.78895e-01	            22353005  - 6.94317e-01	            22356165  - 2.72059e-02	            22356545  - 7.15298e-03	
// 22	        430511220896  + 2.87151e-03	        430511313631  + 3.79946e-02	        430511209946  + 1.36145e-03	        430511192333  + 1.01086e-04	        430511321801  + 4.30039e-02	
// 23	      35186147272591  + 4.61271e+00	      35186140180182  + 9.06421e-01	      35186135873732  + 5.11060e-02	      35186148116716  + 5.24423e+00	      35186145110595  + 3.17997e+00	
// 24	    2179984552996022  - 1.50375e-01	    2179984559999133  - 5.65450e-02	    2179984555537471  - 1.11123e-01	    2179984543308940  - 3.54332e-01	    2179984555051521  - 1.18170e-01	
// 
//               18.149 15 =>     8.000902e-01            7.204 15 =>     7.364054e-02           34.599 15 =>   1-1.683486e-03           15.678 15 =>     6.669651e-01           37.390 15 =>   1-6.457611e-04
// 
// --- Finished -- ReMort trials: 2251799813685248	2^51.0 (out of 2^51) trials -- lantern320 -- 64 bits: 	ps:   8.001e-01    7.364e-02  1-1.683e-03    6.670e-01  1-6.458e-04*  => p <   2.57712e-02   2^54.32 calls, 2^57.32 bytes	2^38.05 bytes/second	used:   7::07:19:07.78

// 1	                   0  - 9.76563e-04	                   0  - 9.76563e-04	                   0  - 9.76563e-04	                   2  + 2.30422e-01	                   0  - 9.76563e-04	                   1  - 1.27592e-01	
// 2	                   0  - 3.07617e-02	                   0  - 3.07617e-02	                   0  - 3.07617e-02	                  37  + 1.19894e-01	                   0  - 3.07617e-02	                  28  - 1.38308e+00	
// 3	                   0  - 6.35742e-01	                   0  - 6.35742e-01	                   0  - 6.35742e-01	                 597  + 2.20092e+00	                   2  + 2.92760e+00	                 579  + 5.24398e-01	
// 4	                   6  - 1.40830e+00	                  12  + 5.47981e-01	                  11  + 1.75640e-01	                6598  - 6.53494e-01	                  12  + 5.47981e-01	                6674  + 1.50315e-02	
// 5	                 110  - 3.45588e-01	                 111  - 2.45179e-01	                 108  - 5.97978e-01	               62655  + 3.36882e+00	                 124  + 5.04234e-01	               62168  - 1.37594e-02	
// 6	                1117  - 6.38081e-01	                1158  + 1.70884e-01	                1133  - 1.06115e-01	              476301  + 7.75028e-01	                1180  + 1.13171e+00	              476194  + 5.25941e-01	
// 7	                9464  - 2.37583e-02	                9615  + 1.95106e+00	                9336  - 2.15750e+00	             3066148  + 1.04371e-01	                9328  - 2.40564e+00	             3069912  + 6.11494e+00	
// 8	               67494  - 2.85661e-02	               67566  + 1.16716e-02	               67678  + 2.90524e-01	            16985546  - 4.91474e-01	               66996  - 4.34839e+00	            16993953  + 1.79195e+00	
// 9	              419802  - 4.48152e-01	              420069  - 6.63410e-02	              420630  + 3.69459e-01	            82222810  + 6.25536e-01	              419528  - 1.19271e+00	            82214683  - 1.11073e-02	
//10	             2309036  - 2.21343e+00	             2312138  + 3.05404e-01	             2312509  + 6.34676e-01	           351720402  + 1.15677e+00	             2310194  - 5.27171e-01	           351685747  - 5.96561e-01	
//11	            11348307  + 3.30272e-01	            11347939  + 2.16638e-01	            11347093  + 4.59196e-02	          1342797508  - 2.49843e+00	            11345184  - 1.24216e-01	          1342840219  - 1.72314e-01	
//12	            50120169  + 9.86077e-01	            50115541  + 1.15094e-01	            50103393  - 1.89555e+00	          4613008095  + 5.64152e-01	            50118381  + 5.48249e-01	          4612991191  + 2.52220e-01	
//13	           200458847  + 1.97339e-01	           200456623  + 8.24525e-02	           200434705  - 1.58997e+00	         14351339495  - 4.74661e-01	           200469336  + 1.40440e+00	         14351599803  + 2.20209e+00	
//14	           730216508  - 1.69981e-02	           730202730  - 4.09916e-01	           730234359  + 2.81132e-01	         40662372965  + 2.73492e-03	           730173836  - 2.92239e+00	         40662312701  - 6.07915e-02	
//15	          2434097745  + 3.94167e-01	          2434039209  - 3.12082e-01	          2434051448  - 9.64538e-02	        105421630845  + 4.53245e+00	          2434074021  + 2.15983e-02	        105420909320  - 8.69910e-03	
//16	          7454308212  - 6.07042e-02	          7454281975  - 3.02794e-01	          7454239072  - 1.09659e+00	        251107108191  + 3.27506e-01	          7454292520  - 1.83297e-01	        251106570684  - 2.50361e-01	
//17	         21047779175  + 3.22740e+00	         21047616672  + 4.57496e-01	         21047372008  - 1.02020e+00	        551450778435  + 4.60562e-01	         21047721293  + 1.95307e+00	        551450245438  - 1.52880e-03	
//18	         54957203340  - 7.73589e-01	         54957273934  - 3.34557e-01	         54957776992  + 2.45696e+00	       1119921027929  + 1.49192e-01	         54957430705  + 8.15840e-03	       1119921089584  + 1.97593e-01	
//19	        133055537422  + 4.30064e+00	        133055355419  + 2.48012e+00	        133053890754  - 5.95607e+00	       2108857219150  + 3.41425e-01	        133054505716  - 5.69423e-01	       2108853411464  - 4.15228e+00	
//20	        299372627461  - 1.32465e+00	        299372901495  - 4.22622e-01	        299373132992  - 5.15281e-02	       3690500015195  + 5.06151e-01	        299374033228  + 2.01163e+00	       3690497355520  - 4.52977e-01	
//21	        627255312308  - 1.37878e+01	        627259345647  + 1.90281e+00	        627257984395  - 1.15149e-01	       6014145174215  - 9.89944e-02	        627257615864  - 6.47470e-01	       6014144939601  - 1.68348e-01	
//22	       1226005548200  + 4.97127e-01	       1226003654221  - 1.01093e+00	       1226005860700  + 9.74767e-01	       9142713531310  - 1.18057e+00	       1226004409358  - 1.04626e-01	       9142717710578  + 8.73988e-02	
//23	       2238790478859  - 3.11970e-01	       2238790408399  - 3.66792e-01	       2238790729152  - 1.53087e-01	      12985305573020  - 4.32291e-01	       2238789976762  - 7.99434e-01	      12985310592582  + 5.40924e-01	
//24	       3824602409938  + 8.82291e-02	       3824598841957  - 2.33297e+00	       3824599885920  - 9.87219e-01	      17253627246033  + 1.38692e-02	       3824602725951  + 2.10335e-01	      17253629241265  + 3.57739e-01	
//25	       6119365463942  + 1.05206e+00	       6119364770296  + 5.55464e-01	       6119365224402  + 8.62791e-01	      21471184546695  + 9.78117e-01	       6119362314096  - 6.13142e-02	      21471184587188  + 9.95479e-01	
//26	       9179045709107  + 1.89621e-01	       9179040096541  - 2.00807e+00	       9179049466655  + 2.80795e+00	      25049710221917  + 2.77920e-03	       9179047587669  + 1.11409e+00	      25049711563958  + 1.02951e-01	
//27	      12918653531776  - 1.82427e-01	      12918658500275  + 9.12466e-01	      12918656706166  + 2.08000e-01	      27420660055297  - 3.72359e+00	      12918650936178  - 1.32081e+00	      27420668668347  - 8.11348e-02	
//28	      17071079030477  - 4.53097e-02	      17071083805144  + 8.88783e-01	      17071077283390  - 4.04125e-01	      28182362525012  + 1.78010e+00	      17071079006973  - 4.77638e-02	      28182360387340  + 8.67752e-01	
//29	      21191683547618  - 1.62883e-01	      21191680232722  - 1.26265e+00	      21191684215020  - 6.68787e-02	      27210552124532  + 1.53321e-01	      21191683292681  - 2.10651e-01	      27210544470835  - 1.15709e+00	
//30	      24723636464464  + 4.93049e-01	      24723632616581  - 5.13964e-03	      24723641787393  + 3.14244e+00	      24691057970862  + 4.32230e-01	      24723638351518  + 1.17005e+00	      24691054935510  + 2.17024e-03	
//31	      27116234883460  - 2.20495e+00	      27116240409962  - 1.79447e-01	      27116237932044  - 8.09036e-01	      21062618526880  - 1.23813e-01	      27116246678663  + 6.08730e-01	      21062616027448  - 8.03677e-01	
//32	      27963626095171  + 2.88166e-02	      27963624890014  - 3.38104e-03	      27963622682912  - 2.26120e-01	      16893972656919  - 9.07339e-01	      27963618190203  - 1.75593e+00	      16893979969119  + 6.83072e-01	
//33	      27116240047263  - 2.43309e-01	      27116241192391  - 7.47240e-02	      27116246945346  + 6.91266e-01	      12741581547721  - 9.48529e-01	      27116244568037  + 1.40544e-01	      12741580890749  - 1.34090e+00	
//34	      24723639737536  + 1.85079e+00	      24723632285785  - 1.91046e-02	      24723635406222  + 2.39460e-01	       9035702873160  + 1.50636e+00	      24723631067343  - 1.46893e-01	       9035702561031  + 1.26225e+00	
//35	      21191689054763  + 6.28408e-01	      21191686566377  + 6.35913e-02	      21191679038054  - 1.91323e+00	       6023797991824  - 3.55797e-01	      21191683236881  - 2.21925e-01	       6023797361644  - 7.28033e-01	
//36	      17071081965684  + 2.47554e-01	      17071082914367  + 5.28759e-01	      17071078996170  - 4.89135e-02	       3774170178288  - 5.81584e-02	      17071081748002  + 1.97903e-01	       3774171527005  + 2.05282e-01	
//37	      12918650579750  - 1.55858e+00	      12918658208456  + 7.63947e-01	      12918655306805  + 4.45390e-03	       2221436368389  + 3.02992e+00	      12918655656978  + 2.69496e-02	       2221432942945  - 3.10917e-01	
//38	       9179047770347  + 1.24501e+00	       9179047156191  + 8.33731e-01	       9179044545145  + 2.62865e-03	       1227633816797  - 3.30800e-01	       9179045857880  + 2.34798e-01	       1227633176459  - 1.32960e+00	
//39	       6119363907311  + 1.57161e-01	       6119361992675  - 1.42544e-01	       6119361630952  - 2.74341e-01	        636552434108  + 2.39884e+00	       6119361385525  - 3.88115e-01	        636550696200  - 3.96200e-01	
//40	       3824596564613  - 7.24630e+00	       3824602779845  + 2.36372e-01	       3824600090046  - 7.90698e-01	        309433252090  - 5.96186e+00	       3824602546605  + 1.34628e-01	        309434829934  + 1.55857e-01	
//41	       2238791846879  + 1.26559e-01	       2238790866426  - 8.97114e-02	       2238790478016  - 3.12599e-01	        140880848060  + 1.76454e-02	       2238792363415  + 4.91358e-01	        140880472409  - 7.53407e-01	
//42	       1226004661753  - 9.12240e-03	       1226003041496  - 2.42994e+00	       1226004705790  - 3.10692e-03	         60004808368  + 9.55822e-03	       1226006125948  + 1.50518e+00	         60004553093  - 8.91793e-01	
//43	        627259991397  + 4.81701e+00	        627259151810  + 1.28750e+00	        627258218420  - 1.92273e-03	         23878050755  + 1.32521e+00	        627258552906  + 1.43250e-01	         23878095586  + 2.07734e+00	
//44	        299372725757  - 9.43388e-01	        299373502732  + 2.01384e-01	        299372990560  - 2.37475e-01	          8863736501  - 4.21850e-02	        299373278207  + 1.47491e-03	          8863846111  + 9.19387e-01	
//45	        133054834577  + 2.15984e-02	        133054745045  - 9.69951e-03	        133054978150  + 2.92212e-01	          3064012951  - 6.51353e-04	        133053939006  - 5.32790e+00	          3063932164  - 2.20521e+00	
//46	         54957431060  + 8.43425e-03	         54957271563  - 3.46359e-01	         54957103244  - 1.70698e+00	           984294300  - 1.53073e+00	         54957610916  + 7.37956e-01	           984357834  + 6.20660e-01	
//47	         21047462939  - 1.46900e-01	         21047789073  + 3.47718e+00	         21047683229  + 1.28857e+00	           293191997  - 6.31960e-01	         21047437540  - 3.11752e-01	           293215016  + 3.01789e-01	
//48	          7454174351  - 3.22850e+00	          7454345264  + 3.34035e-02	          7454442874  + 1.72480e+00	            80757904  - 1.08980e+00	          7454273042  - 4.27366e-01	            80765406  - 4.37552e-02	
//49	          2434039837  - 2.98022e-01	          2434015380  - 1.08500e+00	          2434034474  - 4.28524e-01	            20509833  - 3.03131e-01	          2434022277  - 8.13314e-01	            20513859  + 1.14484e-01	
//50	           730214894  - 3.61398e-02	           730214632  - 3.99202e-02	           730238374  + 4.60767e-01	             4783673  - 1.34428e+00	           730234312  + 2.79290e-01	             4782060  - 3.59755e+00	
//51	           200446420  - 1.87923e-01	           200443091  - 4.47067e-01	           200456713  + 8.61435e-02	             1021730  - 2.57553e-02	           200478731  + 3.41751e+00	             1021788  - 1.06315e-02	
//52	            50119832  + 8.93798e-01	            50100565  - 3.15517e+00	            50108103  - 5.06159e-01	              199129  + 9.20756e-01	            50105595  - 1.13579e+00	              198468  - 2.73846e-01	
//53	            11345567  - 5.69970e-02	            11346470  + 8.60618e-04	            11351901  + 2.69504e+00	               35196  + 1.19611e+00	            11354522  + 5.85525e+00	               34750  - 1.66563e+00	
//54	             2309935  - 8.03581e-01	             2309007  - 2.27055e+00	             2311906  + 1.60025e-01	                5378  - 4.96498e+00	             2310812  - 1.02122e-01	                5547  + 1.72471e-03	
//55	              420468  + 1.28114e-01	              420020  - 1.10992e-01	              419675  - 7.48834e-01	                 790  + 4.61189e-02	              420503  + 1.69679e-01	                 839  + 3.86031e+00	
//56	               67495  - 2.72801e-02	               67434  - 1.59912e-01	               67161  - 2.10358e+00	                  97  - 1.01710e-02	               67370  - 4.17519e-01	                 114  + 2.61282e+00	
//57	                9440  - 1.60516e-01	                9378  - 1.07631e+00	                9727  + 6.48809e+00	                  21  + 9.92158e+00	                9464  - 2.37583e-02	                  10  - 4.55009e-02	
//58	                1122  - 4.23765e-01	                1118  - 5.91721e-01	                1089  - 2.64593e+00	                   0  - 1.00419e+00	                1183  + 1.32829e+00	                   1  - 1.74675e-05	
//59	                 125  + 6.44498e-01	                  86  - 7.91266e+00	                 131  + 1.84709e+00	                   0  - 7.94273e-02	                 138  + 4.03229e+00	                   0  - 7.94273e-02	
//60	                   5  - 2.27370e+00	                  14  + 1.91153e+00	                   8  - 2.96363e-01	                   0  - 5.14807e-03	                  17  + 5.50404e+00	                   0  - 5.14807e-03	
//61	                   1  + 2.08707e-01	                   0  - 6.35742e-01	                   1  + 2.08707e-01	                   0  - 2.62561e-04	                   1  + 2.08707e-01	                   0  - 2.62561e-04	
//
//              64.175 57 =>     7.750012e-01           49.161 57 =>     2.157539e-01           57.579 57 =>     6.134693e-01           67.151 56 =>     8.197123e-01           61.972 57 =>     6.558949e-01           50.300 56 =>     3.481626e-01
//
//--- Finished -- BBin trials: 281474976710656	2^48.0 (out of 2^48) trials -- lantern320 -- 	ps:   7.750e-01    2.158e-01    6.135e-01    8.197e-01*   6.559e-01    3.482e-01   => p <   6.93024e-01   2^52.58 calls, 2^55.58 bytes	2^36.61 bytes/second	used:   5::23:32:48.36

    // stateA = fa + 0x9E3779B97F4A7C15UL;
    // stateB = fa ^ fe;
    // stateC = fb + fd;
    // stateD = rotate64(fc, 52);
    // return stateE = fb - fc;


// LaceRandom (may replace AceRandom)
// Passes 179 PB of ReMort testing without suspicion.
// Ends with a p-value of p <   5.76107e-01
// On rather terrible hardware for the task, finishes in 7::18:15:57.98, where other generators have taken over 20 days.
// Passes 64TB of PractRand with no anomalies.


// 7	                4419  + 3.39252e-03	                4384  - 2.19487e-01	                4519  + 2.44365e+00	                4454  + 3.42208e-01	                4442  + 1.63530e-01	
// 8	              360777  - 1.61030e-02	              360826  - 2.05458e-03	              361963  + 3.41300e+00	              361747  + 2.21372e+00	              359706  - 3.64728e+00	
// 9	            22358092  + 5.88561e-02	            22358078  + 5.74282e-02	            22356245  - 2.19108e-02	            22356526  - 7.84883e-03	            22359152  + 2.17887e-01	
//11	                4413  - 1.02739e-03	                4466  + 5.86116e-01	                4400  - 5.18470e-02	                4461  + 4.76560e-01	                4442  + 1.63530e-01	
//12	            85013579  - 3.32492e-01	            85021531  + 8.16807e-02	            85018834  - 4.48875e-05	            85010833  - 7.64634e-01	            85021625  + 8.76119e-02	
//13	          6948719976  + 1.83471e-01	          6948842230  + 3.59078e+00	          6948628026  - 4.55258e-01	          6948765441  + 9.48186e-01	          6948670136  - 2.87514e-02	
//14	        430511380241  + 8.78773e-02	        430511249982  + 9.58751e-03	        430510938729  - 1.41721e-01	        430510809254  - 3.29234e-01	        430511420459  + 1.27975e-01	
//16	              360066  - 1.71740e+00	              360278  - 9.16960e-01	              359577  - 4.51363e+00	              360266  - 9.55617e-01	              360725  - 4.55659e-02	
//17	          6948642591  - 2.50001e-01	          6948488622  - 5.50872e+00	          6948847034  + 3.81251e+00	          6948687648  + 1.64167e-03	          6948788994  + 1.57829e+00	
//18	        567922137522  - 2.43954e+00	        567923883627  + 5.70170e-01	        567924706853  + 3.41317e+00	        567922066309  - 2.74365e+00	        567923490576  + 5.45391e-02	
//19	      35186132336810  - 1.37047e-01	      35186130744462  - 4.07864e-01	      35186125668319  - 2.23321e+00	      35186128467560  - 1.04548e+00	      35186130860277  - 3.83307e-01	
//21	            22358255  + 7.67711e-02	            22358556  + 1.16100e-01	            22359323  + 2.52958e-01	            22358007  + 5.04568e-02	            22357567  + 1.73105e-02	
//22	        430510928504  - 1.53697e-01	        430511602125  + 4.02730e-01	        430511246275  + 8.51303e-03	        430510886158  - 2.08466e-01	        430510774032  - 3.93719e-01	
//23	      35186128372074  - 1.07866e+00	      35186130413889  - 4.82151e-01	      35186129803730  - 6.35582e-01	      35186128396852  - 1.07000e+00	      35186127069931  - 1.58283e+00	
//24	    2179984580707929  + 4.23303e-02	    2179984577992192  + 2.17794e-02	    2179984583381421  + 6.91708e-02	    2179984585149732  + 9.05268e-02	    2179984582143184  + 5.59243e-02	
//
//               6.579 15 =>     5.025017e-02           12.974 15 =>     4.716522e-01           21.466 15 =>   1-9.052175e-02           11.248 15 =>     3.332448e-01            8.548 15 =>     1.409673e-01
//
//--- Finished -- ReMort trials: 2251799813685248	2^51.0 (out of 2^51) trials -- lantern320 -- 64 bits: 	ps:   5.025e-02*   4.717e-01  1-9.052e-02    3.332e-01    1.410e-01   => p <   5.76107e-01   2^54.32 calls, 2^57.32 bytes	2^37.96 bytes/second	used:   7::18:15:57.98
//    stateA = fa + 0x9E3779B97F4A7C15UL;
//    stateB = fa ^ fe;
//    stateC = fb + fd;
//    stateD = rotate64(fc, 52);
//    stateE = fb + fc;
//    return fb;


// crand64
// Passes 64TB with one anomaly at 8 TB:
// [Low1/32]BCFN(2+7,13-0,T)         R= +11.4  p =  1.3e-5   unusual
// 256 bits of mutable state, 63 bits of unchanging stream.
// Minimum period is 2 to the 64, ARX, about as fast or a little faster than Xoshiro256**.
// Unless you need many streams, consider LaceRandom instead, since it's faster and still ARX.
// from https://github.com/stclib/STC/blob/master/include/stc/crand.h
//
// Passes ReMort (179 PB of it) rather well; it spends a lot of time between 1E-1 and 1E-2, but rarely gets any worse.
//
// 6	                   0  - 2.29283e-01	                   0  - 2.29283e-01	                   1  + 2.59071e+00	                   0  - 2.29283e-01	                   0  - 2.29283e-01	
// 7	                4445  + 2.02084e-01	                4365  - 5.69179e-01	                4433  + 7.23294e-02	                4517  + 2.35045e+00	                4349  - 9.90492e-01	
// 8	              360108  - 1.53904e+00	              360747  - 3.12718e-02	              360361  - 6.71434e-01	              360765  - 2.15719e-02	              360369  - 6.49786e-01	
// 9	            22350788  - 1.69555e+00	            22350229  - 2.01742e+00	            22352952  - 7.13122e-01	            22352465  - 8.97685e-01	            22350631  - 1.78313e+00	
//11	                4434  + 8.06509e-02	                4402  - 3.90457e-02	                4293  - 3.37831e+00	                4551  + 4.18124e+00	                4446  + 2.15842e-01	
//12	            85012453  - 4.88237e-01	            85016188  - 8.62403e-02	            85016839  - 4.97575e-02	            85029017  + 1.20490e+00	            85010333  - 8.62410e-01	
//13	          6948647997  - 1.89355e-01	          6948679614  - 3.12046e-03	          6948764803  + 9.33340e-01	          6948653994  - 1.31920e-01	          6948623502  - 5.31440e-01	
//14	        430511379502  + 8.72108e-02	        430511344182  + 5.83146e-02	        430511890182  + 1.15269e+00	        430511988555  + 1.49710e+00	        430511406211  + 1.12910e-01	
//16	              359625  - 4.18050e+00	              361473  + 1.06447e+00	              360789  - 1.14321e-02	              360862  + 2.13206e-04	              360388  - 5.99794e-01	
//17	          6948531755  - 3.34754e+00	          6948415899  - 1.03650e+01	          6948575376  - 1.70651e+00	          6948642018  - 2.56923e-01	          6948658529  - 9.53598e-02	
//18	        567922234663  - 2.05349e+00	        567923531404  + 8.27786e-02	        567923629325  + 1.74431e-01	        567925252139  + 6.61027e+00	        567924926203  + 4.57337e+00	
//19	      35186133600185  - 2.47167e-02	      35186132417452  - 1.27167e-01	      35186141848070  + 1.52088e+00	      35186140158541  + 8.99487e-01	      35186130797067  - 3.96615e-01	
//21	            22353668  - 4.80301e-01	            22349474  - 2.49651e+00	            22350266  - 1.99525e+00	            22352314  - 9.59220e-01	            22352893  - 7.34353e-01	
//22	        430512125854  + 2.05296e+00	        430511608040  + 4.14253e-01	        430511447844  + 1.59579e-01	        430511998955  + 1.53614e+00	        430512001296  + 1.54500e+00	
//23	      35186143151257  + 2.11102e+00	      35186132170422  - 1.58603e-01	      35186131987698  - 1.84087e-01	      35186140127127  + 8.89470e-01	      35186140483951  + 1.00655e+00	
//24	    2179984563568514  - 2.60318e-02	    2179984575071357  + 7.22855e-03	    2179984565092016  - 1.65673e-02	    2179984556399428  - 9.91553e-02	    2179984566345080  - 1.03787e-02	
//
//              18.559 15 =>     8.173645e-01           17.521 15 =>     7.707619e-01           12.740 15 =>     4.533926e-01           21.536 15 =>   1-8.891471e-02           14.107 15 =>     5.587667e-01
//
//--- Finished -- ReMort trials: 2251799813685248	2^51.0 (out of 2^51) trials -- crand320 -- 64 bits: 	ps:   8.174e-01    7.708e-01    4.534e-01  1-8.891e-02*   5.588e-01   => p <   1.94859e-01   2^54.32 calls, 2^57.32 bytes	2^37.47 bytes/second	used:  10::21:40:48.42

    //const uint64_t result = (stateA ^ (stateD += stateE)) + stateB;
    //stateA = stateB ^ (stateB >> 11);
    //stateB = stateC + (stateC << 3);
    //stateC = ((stateC << 24) | (stateC >> (64 - 24))) + result;
    //return result;
	
	// Tested with the "worst-case" stream, this still passses 64TB fine!
	// Period is 2 to the 64 minimum, but likely 2 to the 319 expected. 2 to the 24 possible streams.
    //stateA = fa + stream;
    //stateB = fa ^ fe;
    //stateC = fb + fd;
    //stateD = rotate64(fc, 44);
    //stateE = fb + fc;
    //return fb;


// BassD in https://quick-bench.com/q/E-wxFyZQMlGJaN-Rflh_IkObTOU
// Passes 64TB with no anomalies!
// Minimim period is 2 to the 64, expected period is a much higher multiple.
stateA = fb ^ fd;
stateB = fe ^ fc;
stateC = fa + fb;
stateD = rotate64(fc, 52);
stateE = fe + 0x9E3779B97F4A7C15UL;
return fa;

				}
				std::string overload320::get_name() const { return "overload320"; }
				void overload320::walk_state(StateWalkingObject *walker) {
					walker->handle(stateA);
					walker->handle(stateB);
					walker->handle(stateC);
					walker->handle(stateD);
					walker->handle(stateE);
					//stateE |= 1UL;
					walker->handle(stream);
					// worst case?
					//stream = 0x9E3779B97F4A7C15UL;
					//stream = (stream & 0x003569CA5369AC00UL) ^ 0x9E3779B97F4A7C15UL;
					printf("stateA: %016llX\n", stateA);
					printf("stateB: %016llX\n", stateB);
					printf("stateC: %016llX\n", stateC);
					printf("stateD: %016llX\n", stateD);
					printf("stateE: %016llX\n", stateE);
					printf("stream: %016llX\n", stream);
					//stream = 1UL;
				}

// Joker40, a shrunken-down AceRandom.
// Running this command,
// ./RNG_test joker40 -te 1 -tlmin 10 -tlmax 50 -multithreaded -seed 0
// This version fails at the 4GB mark.
// The length before failure depends on the length of the subcycle.
				Uint8 joker40::raw8() {
	const uint8_t fa = stateA;
	const uint8_t fb = stateB;
	const uint8_t fc = stateC;
	const uint8_t fd = stateD;
	const uint8_t fe = stateE;
    stateA = fa + 0x97U;
    stateB = fa ^ fe;
    stateC = fb + fd;
    stateD = rotate8(fc, 6);
    stateE = fb - fc;
	return fe;
//    stateA = fa + 0x95U;
//    stateB = fa ^ rotate8(fe, 5);
//    stateC = fb + rotate8(fd, 2);
//    stateD = fe ^ rotate8(fc, 6);
//    return stateE = fc + rotate8(fb, 3);
	

				}
				std::string joker40::get_name() const { return "joker40"; }
				void joker40::walk_state(StateWalkingObject *walker) {
					walker->handle(stateA);
					walker->handle(stateB);
					walker->handle(stateC);
					walker->handle(stateD);
					walker->handle(stateE);
				}

				Uint16 suit80::raw16() {
					const uint16_t fa = stateA;
					const uint16_t fb = stateB;
					const uint16_t fc = stateC;
					const uint16_t fd = stateD;
					const uint16_t fe = stateE;
					// Fails at 1TB (due to DC6-9x1Bytes-1) with the command:
					// ./RNG_test suit80 -tf 2 -tlmin 10 -tlmax 50 -multithreaded -seed 0
						// stateA = fa + 0x9E3DU;
						// stateB = fa ^ fe;
						// stateC = fb + fd;
						// stateD = rotate16(fc, 5);
						// return stateE = fb - fc;
					// Shockingly, passes 64TB with just one anomaly (a very-low-significance DC6-9x1Bytes-1 at 2 KB).
					// This was run with the same command as above.
					stateA = fa + 0x9E3DU;
					stateB = fa ^ fe;
					stateC = fb + fd;
					stateD = rotate16(fc, 11);
					stateE = fb + fc;
					return fb;
				}
				std::string suit80::get_name() const { return "suit80"; }
				void suit80::walk_state(StateWalkingObject* walker) {
					walker->handle(stateA);
					walker->handle(stateB);
					walker->handle(stateC);
					walker->handle(stateD);
					walker->handle(stateE);
				}

				union int_to_float_bits {
					uint32_t integer_bits;
					float converted_float_bits;
				};

				static float intBitsToFloat(uint32_t int_value)
				{
					union int_to_float_bits bits;
					bits.integer_bits = int_value;
					return bits.converted_float_bits;
				}

				union long_to_double_bits {
					uint64_t long_bits;
					double converted_double_bits;
				};

				static double longBitsToDouble(uint64_t long_value)
				{
					union long_to_double_bits bits;
					bits.long_bits = long_value;
					return bits.converted_double_bits;
				}

				Uint32 floatHax32::rand() {
					uint32_t x, y;
					x = (stateA += 0xDB4F0B91);
					y = (stateB += (rotate32(x, 21) + __lzcnt(x)));
					x ^= rotate32(y, 21);
					x ^= x >> 15;
					x *= 0x2c1b3c6d;
					x ^= x >> 12;
					x *= 0x297a2d39;
					x ^= x >> 15;
					return x;
				}

				Uint16 floatHax32::raw16() {
					//uint32_t proto_exp_offset = rand();
					//if (proto_exp_offset == 0) {
					//	return 0;
					//}
					//float f = ldexpf((float)(rand() | 0x80000001), -32 - __lzcnt(proto_exp_offset));
					uint32_t bits = rand();
					float f = intBitsToFloat((126u - __lzcnt(rand()) << 23 & 0u - ((bits | 0u - bits) >> 31)) | (bits & 0x7FFFFF));

					// gets the low 16 bits of the random float f after scaling
					return (uint16_t)(f * 0x800000);
				}
				std::string floatHax32::get_name() const { return "floatHax32"; }
				void floatHax32::walk_state(StateWalkingObject* walker) {
					walker->handle(stateA);
					walker->handle(stateB);
				}

				Uint64 floatHax64::rand() {
					//uint64_t x = (state += 0x9E3779B97F4A7C15UL);
					//x ^= x >> 32;
					//x *= 0xbea225f9eb34556dUL;
					//x ^= x >> 29;
					//x *= 0xbea225f9eb34556dUL;
					//x ^= x >> 32;
					//x *= 0xbea225f9eb34556dUL;
					//x ^= x >> 29;
					//return x;
					uint64_t x = (state += 0xD1342543DE82EF95ULL);
					x = (x ^ x >> 27) * 0x3C79AC492BA7B653ULL;
					x = (x ^ x >> 33) * 0x1C69B3F74AC4AE35ULL;
					x ^= x >> 27;
					return x;
				}

				Uint16 floatHax64::raw16() {
					// uses what Godot does, except that it uses Moremur instead of PCG-Random.
					// Fails at 4TB, with the same issue every inclusive-on-1 float generator has so far.
//rng=floatHax64, seed=0x0
//length= 4 terabytes (2^42 bytes), time= 65724 seconds
//  Test Name                         Raw       Processed     Evaluation
//  FPF-14+6/16:cross                 R= +19.2  p =  1.7e-16    FAIL !
//  ...and 1051 test result(s) without anomalies
					uint32_t proto_exp_offset = (uint32_t)rand();
					if (proto_exp_offset == 0) {
						return 0;
					}
					float f = ldexpf((float)((uint32_t)rand() | 0x80000001), -32 - __lzcnt(proto_exp_offset));

					//uint64_t bits = rand();
					//if ((bits & 0xFFFFFFFF00000000ULL) == 0ULL) {
					//	return 0;
					//}
					//float f = ldexpf((float)((uint32_t)(bits) | 0x80000001), -32 - __lzcnt64(bits));

					//uint64_t bits = rand();
					//float f = intBitsToFloat((126u - (uint32_t)__lzcnt64(bits) << 23 & 0u - (uint32_t)((bits | 0ULL - bits) >> 63)) | ((uint32_t)bits & 0x7FFFFF));

					// better, but not immune? gets "unusual" at 256GB, the low bit of the same troublesome test... [Low1/16]FPF-14+6/16:cross
					//uint64_t x = (state += 0x9E3779B97F4A7C15ULL);
					//x = (x ^ x >> 27) * 0x3C79AC492BA7B653ULL;
					//x = (x ^ x >> 33) * 0x1C69B3F74AC4AE35ULL;
					//float f = (0x1000001ULL * (x >> 37) >> 27) * 5.9604645E-8f;

					// even with MX3, this is "very suspicious" at only 256GB -- the construction is likely at fault.
					//uint64_t x = rand();
					//float f = intBitsToFloat((126u - ((uint32_t)__lzcnt64(x)) << 23) | ((uint32_t)x & 0x7FFFFFu) + 1u) - 2.7105058E-20f;


					// gets the low 16 bits of the random float f after scaling
					//return (uint16_t)(f * 0x800000);
					// high 16 bits
					return (uint16_t)(f * 0x10000);
				}
				std::string floatHax64::get_name() const { return "floatHax64"; }
				void floatHax64::walk_state(StateWalkingObject* walker) {
					walker->handle(state);
				}


				Uint32 floatHaxPCG::rand() {
					uint64_t oldstate = state;
					// Advance internal state
					state = oldstate * 6364136223846793005ULL + 1442695040888963407ULL;
					// Calculate output function (XSH RR), uses old state for max ILP
					uint32_t xorshifted = oldstate >> 27u ^ oldstate >> 45u;
					uint32_t rot = oldstate >> 59u;
					return (xorshifted >> rot) | (xorshifted << (32u - rot & 31u));
				}

				Uint16 floatHaxPCG::raw16() {
					//uint32_t proto_exp_offset = rand();
					//if (proto_exp_offset == 0) {
					//	return 0;
					//}
					//float f = ldexpf((float)((uint32_t)rand() | 0x80000001), -32 - __lzcnt(proto_exp_offset));

					//// fails at 2TB: FPF-14+6/16:(16,14-0)
					//uint32_t bits = rand();
					//float f = intBitsToFloat((126u - (uint32_t)__lzcnt64(state) << 23 & (bits | 0u - bits)) | (bits & 0x7FFFFFu));

					//// passes at least 64TB with no anomalies (!), on seed 0x0.
					//// This is dependent on the structure of PCG-Random for how well it works.
					//// This only has to generate one random uint32_t to function, but it also uses PCG's uint64_t state.
					//// floatHaxPCG::rand() should scramble the output enough that that the state after that output won't be correlated much.
					//// There is detectable correlation after about 2TB, though.
					//// An extra 64-bit multiplication (here, by 0xD1342543DE82EF95ULL) allows it to avoid correlation for our purposes.
					// This works!
					//uint32_t bits = rand();
					//float f = intBitsToFloat((126u - (uint32_t)__lzcnt64(state * 0xD1342543DE82EF95ULL) << 23 & 0u - (bits != 0u)) | (bits & 0x7FFFFFu));
					
					//float f = intBitsToFloat((126u - (uint32_t)__lzcnt64(state * 0xD1342543DE82EF95ULL) << 23 & 0u - ((bits | 0u - bits) >> 31)) | (bits & 0x7FFFFFu));

					// passes at least 512 GB with no anomalies, test interrupted.
					//uint32_t bits = rand();
					//float f = intBitsToFloat((126u - (uint32_t)__lzcnt64(state ^ state << 16u) << 23 & 0u - (bits != 0u)) | (bits & 0x7FFFFFu));

					// passes at least 256GB with no anomalies, but seems slower than above xorshifted state version
					//uint64_t oldstate = state;
					//state = oldstate * 6364136223846793005ULL + 1442695040888963407ULL;
					//uint32_t xorshifted = oldstate >> 27u ^ oldstate >> 45u;
					//uint32_t rot = oldstate >> 59u;
					//uint32_t bits = (xorshifted >> rot) | (xorshifted << (32u - rot & 31u));
					//float f = intBitsToFloat((126u - (uint32_t)__lzcnt64(oldstate * 0xD1342543DE82EF95ULL) << 23 & 0u - (bits != 0u)) | (bits & 0x7FFFFFu));

					// passes 4TB without anomalies. Interrupted because I learned Godot's generator is inclusive on both 0f and 1f . (Whaaaa?!)
					//uint64_t oldstate = state;
					//state = oldstate * 6364136223846793005ULL + 1442695040888963407ULL;
					//uint32_t xorshifted = oldstate >> 27u ^ oldstate >> 45u;
					//uint32_t rot = oldstate >> 59u;
					//uint32_t bits = (xorshifted >> rot) | (xorshifted << (32u - rot & 31u));
					//float f = intBitsToFloat((126u - (uint32_t)__lzcnt64(oldstate ^ oldstate << 37u) << 23 & 0u - (bits != 0u)) | (bits & 0x7FFFFFu));

					//// Passes 64TB with no anomalies, took 197985 seconds.
					//// This generator is, like Godot's randf(), inclusive on 0f and 1f.
					//uint64_t oldstate = state;
					//state = oldstate * 6364136223846793005ULL + 1442695040888963407ULL;
					//uint32_t xorshifted = oldstate >> 27u ^ oldstate >> 45u;
					//uint32_t rot = oldstate >> 59u;
					//uint32_t bits = (xorshifted >> rot) | (xorshifted << (32u - rot & 31u));
					//float f = intBitsToFloat((126u - (uint32_t)__lzcnt64(oldstate ^ oldstate << 37u) << 23) | (bits & 0x7FFFFFu) + 1u) - 2.7105058E-20f;

					//// gets the low 16 bits of the random float f after scaling
					//return (uint16_t)(f * 0x800000);

//					// fails at 512GB for the upper 16 bits!
////rng=floatHaxPCG, seed=0x0
////length= 512 gigabytes (2^39 bytes), time= 2637 seconds
////  Test Name                         Raw       Processed     Evaluation
////  FPF-14+6/16:cross                 R= +15.2  p =  3.9e-13    FAIL
////  ...and 951 test result(s) without anomalies
//					uint64_t oldstate = state;
//					state = oldstate * 6364136223846793005ULL + 1442695040888963407ULL;
//					uint32_t xorshifted = oldstate >> 27u ^ oldstate >> 45u;
//					uint32_t rot = oldstate >> 59u;
//					uint32_t bits = (xorshifted >> rot) | (xorshifted << (32u - rot & 31u));
//					float f = intBitsToFloat((126u - (uint32_t)__lzcnt64(oldstate ^ oldstate << 37u) << 23) | (bits & 0x7FFFFFu) + 1u) - 2.7105058E-20f;
//
					// No matter what variants I try using 64-bit oldstate for clz, this fails at or by 512GB... FPF-14+6/16:cross every time.
					// It's the upper 16 bits of the float that fail, and only those.
					uint64_t oldstate = state;
					state = oldstate * 6364136223846793005ULL + 1442695040888963407ULL;
					uint32_t xorshifted = oldstate >> 27u ^ oldstate >> 45u;
					uint32_t rot = oldstate >> 59u;
					uint32_t bits = (xorshifted >> rot) | (xorshifted << (32u - rot & 31u));
					float f = intBitsToFloat((126u - (uint32_t)__lzcnt64(oldstate ^ (oldstate << rot | oldstate >> 64u - rot) ^ (oldstate << 37u | oldstate >> 27u)) << 23) | (bits & 0x7FFFFFu) + 1u) - 2.7105058E-20f;

					// gets the high 16 bits of the random float f after scaling
					return (uint16_t)(f * 0x10000);
				}
				std::string floatHaxPCG::get_name() const { return "floatHaxPCG"; }
				void floatHaxPCG::walk_state(StateWalkingObject* walker) {
					walker->handle(state);
				}


				Uint64 doubleHaxFlow::rand() {
					uint64_t x = (stateA += 0xD1B54A32D192ED03ULL);
					uint64_t y = (stateB += 0x8CB92BA72F3D8DD7ULL);
					x = (x ^ rotate64(y, 37)) * 0x3C79AC492BA7B653ULL;
					y = (x ^ x >> 33) * 0x1C69B3F74AC4AE35ULL;
					return y ^ y >> 27;
				}

				Uint32 doubleHaxFlow::raw32() {
					//// Passes 64TB of testing with no anomalies.
					//// Tested on the low 32 bits of each returned double, after scaling the double by (1ULL << 52).
					//// Again, passes 64TB of testing with no anomalies.
					//// Tested on the "high 32 bits" of each returned double, after scaling the double by (1ULL << 32).
					uint64_t x = (stateA += 0xD1B54A32D192ED03ULL);
					uint64_t y = (stateB += 0x8CB92BA72F3D8DD7ULL);
					x = (x ^ rotate64(y, 37)) * 0x3C79AC492BA7B653ULL;
					y = (x ^ x >> 33) * 0x1C69B3F74AC4AE35ULL;
					y ^= y >> 27;
					double d = longBitsToDouble((1022ULL - __lzcnt64(x) << 52 & 0ULL - (y != 0ULL)) | (y & 0xFFFFFFFFFFFFFULL));

					//// gets the low 32 bits of the random float f after scaling
					//return (uint32_t)(d * (1ULL << 52));
					//// gets the high 32 bits of the random float f after scaling
					return (uint32_t)(d * (1ULL << 32));
				}
				std::string doubleHaxFlow::get_name() const { return "doubleHaxFlow"; }
				void doubleHaxFlow::walk_state(StateWalkingObject* walker) {
					walker->handle(stateA);
					walker->handle(stateB);
				}

			}
		}
	}
}
