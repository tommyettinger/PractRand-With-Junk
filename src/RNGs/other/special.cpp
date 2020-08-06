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
                    //printf("%d, %d : ", idx, idx & 3);
                    if((idx++) == 0){
                        state = _mm_add_epi64(state, increment);
                        __m128i penultimate = _mm_aesenc_si128(state, increment);
                        __m128i penultimate1 = _mm_aesenc_si128(penultimate, increment); 
                        __m128i penultimate2 = _mm_aesdec_si128(penultimate, increment); 
                        __m256i full = _mm256_insertf128_si256(_mm256_castsi128_si256(penultimate2), penultimate1, 1);
                        _mm256_storeu_si256((__m256i *) buf, full);
                    }
                    //printf("%X  ", buf[0]);
                    return buf[idx &= 3];
                }
                std::string aesdragontamer::get_name() const {
                	return "aesdragontamer";
                }
                void aesdragontamer::walk_state(StateWalkingObject *walker) {
                    idx = 0;
                    increment = _mm_set_epi64x(0xCB9C59B3F9F87D4Du, 
                    		0x3463A64C060782B1u);
//                  increment = _mm_set_epi8(0x2f, 0x2b, 0x29, 0x25, 0x1f, 0x1d, 0x17, 0x13, 
//                        0x11, 0x0D, 0x0B, 0x07, 0x05, 0x03, 0x02, 0x01);
                    Uint64 seed1, seed2;
                	walker->handle(seed1);
                	walker->handle(seed2);
                    state = _mm_set_epi64x(seed1, seed2);
                }
                aesdragontamer::aesdragontamer(){
                    idx = 0;
                }

			}//NotRecommended
		}//Polymorphic
	}//RNGs
}//PractRand
