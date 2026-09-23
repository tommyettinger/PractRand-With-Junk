
/*
RNGs in the mediocre directory are not intended for real world use
only for research; as such they may get pretty sloppy in some areas

This set is of RNGs that use at least one of the following:
1. complex flow control
2. complex math functions (sqrt, log/exp, sin/cos, etc)
3. anything else not covered by simple.h, mult.h, indirection.h, or fibonacci.h
*/
#include <immintrin.h>

namespace PractRand {
	namespace RNGs {
		namespace Polymorphic {
			namespace NotRecommended {
				//class icg32_16;
				//class eicg32_16;                  
                

				class aesdragontamer : public vRNG64 {
					__m128i state = _mm_set_epi32(0, 0, 0, 0);
					const __m128i increment = _mm_set_epi32(0xCB9C59B3, 0xF9F87D4D, 0x060782B1, 0x9E3779B9);
					Uint64 buf[4];
					int idx;
				public:
					aesdragontamer() : buf(), idx(0)
					{}
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class arsenic64 : public vRNG64 {
					__m128i state = _mm_set_epi32(0, 0, 0, 0);
					const __m128i k0 = _mm_set_epi32(0xC13FA9A9, 0x02A6328F, 0x91E10DA5, 0xC79E7B1D);
					const __m128i k1 = _mm_set_epi32(0xD1B54A32, 0xD192ED03, 0xABC98388, 0xFB8FAC03);
					const __m128i k2 = _mm_set_epi32(0xDB4F0B91, 0x75AE2165, 0xBBE05633, 0x03A4615F);
					// const __m128i k3 = _mm_set_epi32(0xE19B01AA, 0x9D42C633, 0xC6D1D6C8, 0xED0C9631);
					// const __m128i k4 = _mm_set_epi32(0xE60E2B72, 0x2B53AEEB, 0xCEBD76D9, 0xEDB6A8EF);
				public:
					arsenic64()
					{}
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};


			}
		}
	}
}
