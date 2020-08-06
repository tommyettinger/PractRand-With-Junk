
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
					__m128i state;
					__m128i increment;
					Uint64 buf[4];
					int idx;
				public:
					aesdragontamer();
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};


			}
		}
	}
}
