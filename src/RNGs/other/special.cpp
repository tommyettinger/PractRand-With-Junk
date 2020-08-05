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
				Uint64 aesdragontamer::raw64() {
                      state = _mm_add_epi64(state, increment);
                      __m128i penultimate = _mm_aesenc_si128(state, increment);
                      return _mm_cvtsi128_si64(_mm_aesenc_si128(penultimate, increment));
				}
				std::string aesdragontamer::get_name() const {
					return "aesdragontamer";
				}
				void aesdragontamer::walk_state(StateWalkingObject *walker) {
                    increment = _mm_set_epi64x(0xCB9C59B3F9F87D4Du, 
                    		0x3463A64C060782B1u);
//                    increment = _mm_set_epi8(0x2f, 0x2b, 0x29, 0x25, 0x1f, 0x1d, 0x17, 0x13, 
//                    		  0x11, 0x0D, 0x0B, 0x07, 0x05, 0x03, 0x02, 0x01);
                    Uint64 seed1, seed2;
					walker->handle(seed1);
					walker->handle(seed2);
                    state = _mm_set_epi64x(seed1, seed2);
				}

			}//NotRecommended
		}//Polymorphic
	}//RNGs
}//PractRand
