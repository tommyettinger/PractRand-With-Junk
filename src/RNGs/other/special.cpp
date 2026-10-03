#include <string>
#include <sstream>
#include <cstdlib>
#include "PractRand/config.h"
#include "PractRand/rng_basics.h"
#include "PractRand/rng_helpers.h"
#include "PractRand/rng_internals.h"

#include "PractRand/RNGs/other/special.h"

using namespace PractRand::Internals;

namespace PractRand {
	namespace RNGs {
		namespace Polymorphic {
			namespace NotRecommended {
                /*
                Passes 64TB with no anomalies with stream constantly 0, PractRand seed 1.
                With the increment constantly changing, this has one anomaly at 512GB instead, and still passes 64TB.
                The increment shouldn't change by adding an odd number; that would result in an even increment and that halves the period, or worse.
                The failures that aesrand normally gets at 8TB when incrementing only the low 64-bits of the state are quick to remedy.
                If you have a 64-bit stream value, you can split it up into two 32-bit sections, then set only  the upper halves of each 64-bit section.
                This can be done with:

                increment = _mm_set_epi32(streamA, 0xF9F87D4D, streamB, 0x9E3779B9);

                Where streamA and streamB are each 32-bit halves of the 64-bit stream. This has been tested with streamA=0 and streamB=0, which had
                no anomalies through 64TB of testing. Because the increment must have each 64-bit half be an odd number to ensure the generator is
                full-period (it seems to not pass testing if one of the increment halves is even, as well), this only allows changing the upper bits,
                which certainly seems to work well.

                The performance seems... erratic at best. Windows 11 might be slowing down the AES instructions when the command-line window is not
                in the foreground; I have no other explanation for why this slows down by a large factor when the window is minimized.
                */
                Uint64 aesdragontamer::raw64() {
                    if((idx++) == 0){
                        state = _mm_add_epi64(state, increment);
                        __m128i penultimate = _mm_aesenc_si128(state, increment);
                        _mm256_storeu_si256((__m256i *) buf, _mm256_inserti128_si256(
                            _mm256_castsi128_si256(_mm_aesdec_si128(penultimate, increment)),
                            _mm_aesenc_si128(penultimate, increment), 1));
                        ////Tested this too, it works to at least 64TB, constantly changing seed.
                        //increment = _mm_add_epi64(increment, _mm_set_epi64x(0u, 0x9E3779B97F4A7C16u));
                    }
                    return buf[idx &= 3];
                }
                std::string aesdragontamer::get_name() const {
                	return "aesdragontamer";
                }
                void aesdragontamer::walk_state(StateWalkingObject *walker) {
                    idx = 0;
                    //increment = _mm_set_epi32(0xCB9C59B3, 0xF9F87D4D, 0x060782B1, 0x9E3779B9); // works
                    
                    //increment = _mm_set_epi64x(0xCB9C59B3F9F87D4Du, 0x3463A64C060782B1u);
//                  increment = _mm_set_epi8(0x2f, 0x2b, 0x29, 0x25, 0x1f, 0x1d, 0x17, 0x13, 
//                        0x11, 0x0D, 0x0B, 0x07, 0x05, 0x03, 0x02, 0x01);
                    Uint64 seed1, seed2;
                	walker->handle(seed1);
                	walker->handle(seed2);
                    state = _mm_set_epi64x(seed1, seed2);
                }

                Uint64 arsenic64::raw64() {
                	// Passes 128TB with no anomalies.
                	// Not equidistributed. Period is 2 to the 64. 2 to the 64 possible streams.
                	// Uses AES-NI instructions and _mm_add_epi64().
                	// k1 and k2 can be changed to any 256 bits of state, in theory.
                 // state = _mm_add_epi64(state, k0);
                 // auto res = _mm_aesenc_si128(state, k1);
                 // res = _mm_aesenclast_si128(res, k2);
                 // return res[0] ^ res[1];

                	// Fails BRank immediately (2GB).
                	// state = _mm_add_epi64(state, k0);
                    // auto res = _mm_aesenc_si128(state, k1);
                    // return res[0] ^ res[1];

                	// Fails BRank immediately (2GB).
                	// state = _mm_add_epi64(state, k0);
                	// auto res = _mm_aesdec_si128(state, k1);
                	// return res[0] ^ res[1];

                	// Fails almost everything immediately (2GB).
                	// state = _mm_add_epi64(state, k0);
                	// auto res = _mm_aesenclast_si128(state, k2);
                	// return res[0] ^ res[1];

                	// Using k1 for both enc and enclast steps has some anomalies:
// rng=arsenic64, seed=0x0
// length= 32 gigabytes (2^35 bytes), time= 40.9 seconds
//   Test Name                         Raw       Processed     Evaluation
//   [Low4/16]BCFN(2+0,13-0,T)         R=  +8.8  p =  2.9e-4   unusual
//   ...and 804 test result(s) without anomalies
//
// rng=arsenic64, seed=0x0
// length= 512 gigabytes (2^39 bytes), time= 633 seconds
//   Test Name                         Raw       Processed     Evaluation
//   BCFN(2+1,13-0,T)                  R=  -8.3  p =1-2.5e-4   unusual
//   ...and 951 test result(s) without anomalies
                	// Both are very borderline BCFN anomalies with opposed p-values.
                	// state = _mm_add_epi64(state, k0);
                	// auto res = _mm_aesenc_si128(state, k1);
                	// res = _mm_aesenclast_si128(res, k1); // same key as the enc step
                	// return res[0] ^ res[1];

                	// Using all 0 for the key in both enc and enclast has one anomaly:
// rng=arsenic64, seed=0x0
// length= 16 gigabytes (2^34 bytes), time= 21.4 seconds
//   Test Name                         Raw       Processed     Evaluation
//   [Low1/32]BCFN(2+2,13-2,T)         R= +10.6  p =  5.5e-5   unusual
//   ...and 768 test result(s) without anomalies
                	// state = _mm_add_epi64(state, k0);
                	// auto res = _mm_aesenc_si128(state, k1);
                	// res = _mm_aesenclast_si128(res, k1); // same key as the enc step
                	// return res[0] ^ res[1];

                	// Gets through 512GB with no anomalies, but then:
// rng=arsenic64, seed=0x0
// length= 1 terabyte (2^40 bytes), time= 1256 seconds
//   Test Name                         Raw       Processed     Evaluation
//   FPF-14+6/16:all                   R=  +8.2  p =  3.8e-7   suspicious
//   ...and 987 test result(s) without anomalies
//
// rng=arsenic64, seed=0x0
// length= 2 terabytes (2^41 bytes), time= 2431 seconds
//   Test Name                         Raw       Processed     Evaluation
//   FPF-14+6/16:(4,14-0)              R=  +8.0  p =  5.5e-7   unusual
//   FPF-14+6/16:(5,14-0)              R=  +8.2  p =  3.2e-7   mildly suspicious
//   FPF-14+6/16:(9,14-0)              R=  +9.3  p =  3.4e-8   suspicious
//   FPF-14+6/16:(10,14-0)             R=  +8.8  p =  8.7e-8   mildly suspicious
//   FPF-14+6/16:(11,14-0)             R=  +7.6  p =  1.2e-6   unusual
//   FPF-14+6/16:all                   R= +16.5  p =  5.3e-15    FAIL
//   ...and 1014 test result(s) without anomalies
                	// It fails at 2TB.
                	// state = _mm_add_epi64(state, k0);
                	// auto res = _mm_aesdec_si128(state, k3);
                	// return state[0] ^ res[0] ^ state[1] ^ res[1];

                	// Passes 128TB with two "unusual" anomalies:
// rng=arsenic64, seed=0x0
// length= 16 terabytes (2^44 bytes), time= 19856 seconds
//   Test Name                         Raw       Processed     Evaluation
//   [Low4/32]FPF-14+6/16:all          R=  -5.0  p =1-1.9e-4   unusual
//   ...and 1106 test result(s) without anomalies
//
// rng=arsenic64, seed=0x0
// length= 32 terabytes (2^45 bytes), time= 40104 seconds
//   Test Name                         Raw       Processed     Evaluation
//   [Low1/32]Gap-16:A                 R=  +7.2  p =  4.3e-5   unusual
//   ...and 1133 test result(s) without anomalies
                	// state = _mm_add_epi64(state, k0);
                	// auto enc = _mm_aesdec_si128(state, k3);
                	// auto res = enc[0] + enc[1];// + state[0] + state[1];
                	// return res ^ rotate64(res, 25) ^ rotate64(res, 50);

                	// Linear Weyl sequences, decrypt, decrypt, return first half of state.
                	// Has two unusual anomalies, both similar BCFN:
// rng=arsenic64, seed=0x0
// length= 64 gigabytes (2^36 bytes), time= 81.5 seconds
//   Test Name                         Raw       Processed     Evaluation
//   [Low1/32]BCFN(2+2,13-1,T)         R=  -8.0  p =1-3.0e-4   unusual
//   ...and 842 test result(s) without anomalies
// rng=arsenic64, seed=0x0
// length= 2 terabytes (2^41 bytes), time= 2425 seconds
//   Test Name                         Raw       Processed     Evaluation
//   [Low1/32]BCFN(2+0,13-0,T)         R=  -8.1  p =1-3.2e-4   unusual
//   ...and 1019 test result(s) without anomalies
                	// state = _mm_add_epi64(state, k0);
                	// auto d0 = _mm_aesdec_si128(state, k1);
                	// auto d1 = _mm_aesdec_si128(d0, k2);
                	// return d1[0];// ^ d1[1];

                	// One "unusual" anomaly at 128TB:
// rng=arsenic64, seed=0x0
// length= 128 terabytes (2^47 bytes), time= 160716 seconds
//   Test Name                         Raw       Processed     Evaluation
//   [Low1/8]Gap-16:A                  R=  -5.0  p =1-3.9e-4   unusual
//   ...and 1180 test result(s) without anomalies
                	// auto res = _mm_aesenc_si128(_mm_aesenc_si128(state = _mm_add_epi64(state, k0), k1), k2);
                	// return res[0];

                	// Gets an "unusual" anomaly at 512GB:
// rng=arsenic64, seed=0x0
// length= 512 gigabytes (2^39 bytes), time= 629 seconds
//   Test Name                         Raw       Processed     Evaluation
//   [Low4/64]DC6-9x1Bytes-1           R=  +5.8  p =  6.8e-4   unusual
//   ...and 951 test result(s) without anomalies
                	// auto res = _mm_aesenc_si128(_mm_aesenc_si128(state = _mm_add_epi64(state, k0), k1), k3);
                	// return res[0];

                	// Gets an unusual anomaly early, 16GB:
// rng=arsenic64, seed=0x0
// length= 16 gigabytes (2^34 bytes), time= 21.4 seconds
//   Test Name                         Raw       Processed     Evaluation
//   [Low1/64]Gap-16:A                 R=  +6.0  p =  1.7e-4   unusual
//   ...and 768 test result(s) without anomalies
                	// state = _mm_add_epi64(state, k0);
                	// __m128i res = _mm_xor_si128(_mm_aesenc_si128(state, k1), state);
                	// return res[0] + res[1];

                	// Should be unsurprising that this gets the same anomaly at 16GB:
// rng=arsenic64, seed=0x0
// length= 16 gigabytes (2^34 bytes), time= 21.2 seconds
//   Test Name                         Raw       Processed     Evaluation
//   [Low1/64]Gap-16:A                 R=  +6.0  p =  1.7e-4   unusual
//   ...and 768 test result(s) without anomalies
                	// state = _mm_add_epi64(state, k0);
                	// __m128i res = _mm_add_epi64(_mm_aesenc_si128(state, k1), state);
                	// return res[0] ^ res[1];

                	// Fails immediately.
                	// state = _mm_add_epi64(state, k0);
                	// __m128i res = _mm_mul_epi32(_mm_aesenc_si128(state, k1), k0);
                	// return res[0] ^ res[0] >> 29;

                	// Also fails immediately.
                	// state = _mm_add_epi64(state, k0);
                	// __m128i res = _mm_add_epi64(_mm_aesenc_si128(state, k1), state);
                	// return res[0] ^ res[0] >> 29;

                	// Grumble...
// rng=arsenic64, seed=0x0
// length= 32 gigabytes (2^35 bytes), time= 40.8 seconds
//   Test Name                         Raw       Processed     Evaluation
//   [Low8/64]BCFN(2+0,13-0,T)         R=  -8.0  p =1-3.7e-4   unusual
//   ...and 804 test result(s) without anomalies
                	// state = _mm_add_epi64(state, k0);
                	// auto res = _mm_aesenc_si128(_mm_aesenc_si128(state, k2), state);
                	// return res[0];

                	auto res = _mm_aesenc_si128(_mm_aesenc_si128(state = _mm_add_epi64(state, k0), k2), k1);
                	return res[0];


                }
                std::string arsenic64::get_name() const {
                	return "arsenic64";
                }
                void arsenic64::walk_state(StateWalkingObject *walker) {
                    Uint64 seed1, seed2;
                	walker->handle(seed1);
                	walker->handle(seed2);
                    state = _mm_set_epi64x(static_cast<long long>(seed1), static_cast<long long>(seed2));
                }

			}//NotRecommended
		}//Polymorphic
	}//RNGs
}//PractRand
