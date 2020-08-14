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
                */
                Uint64 aesdragontamer::raw64() {
                    if((idx++) == 0){
                        state = _mm_add_epi64(state, increment);
                        __m128i penultimate = _mm_aesenc_si128(state, increment);
                        __m128i penultimate1 = _mm_aesenc_si128(penultimate, increment); 
                        __m128i penultimate2 = _mm_aesdec_si128(penultimate, increment); 
                        __m256i full = _mm256_insertf128_si256(_mm256_castsi128_si256(penultimate2), penultimate1, 1);
                        _mm256_storeu_si256((__m256i *) buf, full);
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
                    increment = _mm_set_epi32(0, 0xF9F87D4D, 0, 0x9E3779B9);
                    //increment = _mm_set_epi64x(0xCB9C59B3F9F87D4Du, 0x3463A64C060782B1u);
//                  increment = _mm_set_epi8(0x2f, 0x2b, 0x29, 0x25, 0x1f, 0x1d, 0x17, 0x13, 
//                        0x11, 0x0D, 0x0B, 0x07, 0x05, 0x03, 0x02, 0x01);
                    Uint64 seed1, seed2;
                	walker->handle(seed1);
                	walker->handle(seed2);
                    state = _mm_set_epi64x(seed1, seed2);
                }
                aesdragontamer::aesdragontamer(){
                }

			}//NotRecommended
		}//Polymorphic
	}//RNGs
}//PractRand
