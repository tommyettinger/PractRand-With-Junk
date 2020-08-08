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
                Passes 64TB with one unusual anomaly at 256GB, when tested without the change to increment (which is also the stream).
                The increment shouldn't change by adding an odd number; that would result in an even increment and that halves the period, or worse.
                I use 0x9E3779B97F4A7C16 here instead of the pi-based number; its last bits are 0b10, which allows 2 to the 63 streams.
                More streams are possible by adding a different number to the otherwise-unchanged first 64 bits of increment.
                The other number should probably also end in 0b10, but shouldn't be added every time the later 64 bits are changed.
                */
                Uint64 aesdragontamer::raw64() {
                    if((idx++) == 0){
                        state = _mm_add_epi64(state, increment);
                        __m128i penultimate = _mm_aesenc_si128(state, increment);
                        __m128i penultimate1 = _mm_aesenc_si128(penultimate, increment); 
                        __m128i penultimate2 = _mm_aesdec_si128(penultimate, increment); 
                        __m256i full = _mm256_insertf128_si256(_mm256_castsi128_si256(penultimate2), penultimate1, 1);
                        _mm256_storeu_si256((__m256i *) buf, full);
                        ////Testing this now
                        increment = _mm_add_epi64(increment, _mm_set_epi64x(0u, 0x9E3779B97F4A7C16u));
                    }
                    return buf[idx &= 3];
                }
                std::string aesdragontamer::get_name() const {
                	return "aesdragontamer";
                }
                void aesdragontamer::walk_state(StateWalkingObject *walker) {
                    idx = 0;
                    increment = _mm_set_epi64x(0xCB9C59B3F9F87D4Du, 0x3463A64C060782B1u);
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
