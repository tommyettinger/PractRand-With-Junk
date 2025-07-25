
/*
RNGs in the mediocre directory are not intended for real world use
only for research; as such they may get pretty sloppy in some areas

This set is of RNGs that:
1. use multiplication
2. don't use much indirection, flow control, variable shifts, etc
3. have only a few words of state
4. are likely to have easily detectable bias
*/

namespace PractRand {
	namespace RNGs {
		namespace Polymorphic {
			namespace NotRecommended {
				//similar to the classic LCGs, but with a longer period
				class lcg16of32_extended : public vRNG16 {
					Uint32 state, add;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class lcg32_extended : public vRNG32 {
					Uint32 state, add;
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				//simple classic LCGs
				class lcg32of64_varqual : public vRNG32 {
					Uint64 state;
					int outshift;
				public:
					lcg32of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class lcg16of64_varqual : public vRNG16 {
					Uint64 state;
					int outshift;
				public:
					lcg16of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class lcg8of64_varqual : public vRNG8 {
					Uint64 state;
					int outshift;
				public:
					lcg8of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint8 raw8();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class lcg32of128_varqual : public vRNG32 {
					Uint64 low, high;
					int outshift;
				public:
					lcg32of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class lcg16of128_varqual : public vRNG16 {
					Uint64 low, high;
					int outshift;
				public:
					lcg16of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class lcg8of128_varqual : public vRNG8 {
					Uint64 low, high;
					int outshift;
				public:
					lcg8of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint8 raw8();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				//two LCGs combined
				class clcg8of96_varqual : public vRNG8 {
					Uint64 lcg1;
					Uint32 lcg2;
					int outshift;
				public:
					clcg8of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint8 raw8();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class clcg16of96_varqual : public vRNG16 {
					Uint64 lcg1;
					Uint32 lcg2;
					int outshift;
				public:
					clcg16of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class clcg32of96_varqual : public vRNG32 {
					Uint64 lcg1;
					Uint32 lcg2;
					int outshift;
				public:
					clcg32of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				//LCGs modified by suppressing the carries
				class xlcg32of64_varqual : public vRNG32 {
					Uint64 state;
					int outshift;
				public:
					xlcg32of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xlcg16of64_varqual : public vRNG16 {
					Uint64 state;
					int outshift;
				public:
					xlcg16of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xlcg8of64_varqual : public vRNG8 {
					Uint64 state;
					int outshift;
				public:
					xlcg8of64_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint8 raw8();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xlcg32of128_varqual : public vRNG32 {
					Uint64 low, high;
					int outshift;
				public:
					xlcg32of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xlcg16of128_varqual : public vRNG16 {
					Uint64 low, high;
					int outshift;
				public:
					xlcg16of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xlcg8of128_varqual : public vRNG8 {
					Uint64 low, high;
					int outshift;
				public:
					xlcg8of128_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint8 raw8();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				//modified LCG combined with regular LCG
				class cxlcg8of96_varqual : public vRNG8 {
					Uint64 lcg1;
					Uint32 lcg2;
					int outshift;
				public:
					cxlcg8of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint8 raw8();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class cxlcg16of96_varqual : public vRNG16 {
					Uint64 lcg1;
					Uint32 lcg2;
					int outshift;
				public:
					cxlcg16of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class cxlcg32of96_varqual : public vRNG32 {
					Uint64 lcg1;
					Uint32 lcg2;
					int outshift;
				public:
					cxlcg32of96_varqual(int lcg1_discard_bits) : outshift(lcg1_discard_bits) {}
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class pcg32 : public vRNG32 {
					Uint64 state, inc;
				public:
					pcg32() : state(0x853c49e6748fea9bULL), inc(0xda3e39cb94b95bdbULL) {}
					Uint32 raw32();
					std::string get_name() const;
					void seed(Uint64 s);
					void walk_state(StateWalkingObject *);
				};
				class pcg64 : public vRNG64 {
					Uint64 state, inc;
				public:
					pcg64() : state(0x853c49e6748fea9bULL), inc(0xda3e39cb94b95bdbULL) {}
					Uint64 raw64();
					std::string get_name() const;
					void seed(Uint64 s);
					void walk_state(StateWalkingObject *);
				};
				class pcg32_norot : public vRNG32 {
					Uint64 state, inc;
				public:
					pcg32_norot() : state(0x853c49e6748fea9bULL), inc(0xda3e39cb94b95bdbULL) {}
					Uint32 raw32();
					std::string get_name() const;
					void seed(Uint64 s);
					void walk_state(StateWalkingObject *);
				};
				class cmrg32of192 : public vRNG32 {//I originally encountered this under the name lecuyer3by2b
					//presumably by L'Ecuyer, I adjusted it slightly to output a full 32 bits (instead of ~31.9 bits)
					//it is a Combined Multiple Recursive Generator (the moduli are 2**32-209 and 2**32-22853)
					Uint32 n1m0, n1m1, n1m2, n2m0, n2m1, n2m2;//why is one dimension zero-based and the other not?  no idea, it was that way in the code I based this off of
					Uint32 raw32();
					std::string get_name() const;
					void seed(Uint64 s);
					void walk_state(StateWalkingObject *);
				};
				class xsh_lcg_bad : public vRNG32 {//name was xorwowPlus, I changed it because I wasn't sure it actually qualified as an xorwow
					Uint64 lcg, x0, x1, x2, x3;
					Uint32 raw32();
					std::string get_name() const;
					void seed(Uint64 s);// I also changed the seeding function, because the original permitted the bad all-zeroes state
					void walk_state(StateWalkingObject *);
				};



				//
				class garthy16 : public vRNG16 {
					Uint16 value, scale, counter;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class garthy32 : public vRNG32 {
					Uint32 value, scale, counter;
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				//both sides of the multiply are pseudo-random values in this RNG
				class binarymult16 : public vRNG16 {
					Uint16 a, b, c, d;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class binarymult32 : public vRNG32 {
					Uint32 a, b, c, d;
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class mmr16 : public vRNG16 {
					Uint16 a, b, c;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mmr32 : public vRNG32 {
					Uint32 a, b, c;
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				//uses multiplication, rightshifts, xors, that kind of stuff
				class rxmult16 : public vRNG16 {
					Uint16 a, b, c, d;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				//these are similar to my mwlac algorithm, but lower quality
				class multish2x64 : public vRNG64 {
					Uint64 a, b;
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class multish3x32 : public vRNG32 {
					Uint32 a, b, c;
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class multish4x16 : public vRNG16 {
					Uint16 a, b, c, d;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				//the 16 bit variant of the old version of my mwlac algorithm
				class old_mwlac16 : public vRNG16 {
					Uint16 a, b, c, d;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mwlac_varA : public vRNG16 {
					Uint16 a, b, c, d;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mwlac_varB : public vRNG16 {
					Uint16 a, b, c;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mwlac_varC : public vRNG16 {
					Uint16 a, b, c;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mwlac_varD : public vRNG16 {
					Uint16 a, b, c;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mwlac_varE : public vRNG16 {
					Uint16 a, b, c;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};


				class mwc64x : public vRNG32 {
					Uint64 state;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class cxm64_varqual : public vRNG64 {
					Uint64 low, high;
					int num_mult;
				public:
					cxm64_varqual(int num_mult_) : num_mult(num_mult_) {}
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class mo_Cmfr32 : public vRNG32 {
					Uint32 state;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mo_Cmr32 : public vRNG32 {
					Uint32 state;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mo_Cmr32of64 : public vRNG32 {
					Uint64 state;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class murmlac32 : public vRNG32 {
					Uint32 state1, state2;
					int rounds;
				public:
					Uint32 raw32();
					murmlac32(int rounds_) : rounds(rounds_) {}
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				//multiplication (by a counter), rotate
				class mulcr64 : public vRNG64 {
					Uint64 a, b, count;
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mulcr32 : public vRNG32 {
					Uint32 a, b, count;
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mulcr16 : public vRNG16 {
					Uint32 a, b, count;
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};




				class mo_Cmres64 : public vRNG64 {
					Uint64 state, state2;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};


				class tiptoe64 : public vRNG64 {
					Uint64 state, stream;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class linnormA : public vRNG64 {
					Uint64 state;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class linnormB : public vRNG64 {
					Uint64 state;
					int R;
					uint64_t X;
				public:
				    linnormB(int, int);
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class linnorm32 : public vRNG32 {
					Uint32 stateA;
					int R;
					uint32_t X;

				public:
				    linnorm32(int, int);
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class linnormBounded : public vRNG16 {
					Uint64 state;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mingler : public vRNG64 {
					Uint64 stateA, stateB, incA, incB;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class ta32 : public vRNG32 {
					Uint32 stateA, stateB, stateC, stateD;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class thrust63 : public vRNG64 {
					Uint64 state;
					Uint8 shift1, shift2;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
					void seed(Uint64 seed);
				};
				class vortex : public vRNG64 {
					Uint64 state, stream;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
					void seed(Uint64 seed);
				};
				class vortex2 : public vRNG64 {
					Uint64 state, stream;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
					void seed(Uint64 seed);
				};
				class vortex64 : public vRNG64 {
					Uint64 state, stream;
				public:
					vortex64(int streamPortion) : stream((Uint64)streamPortion << 32 ^ Uint64(0x6C8E9CF570932BD5)) {};
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class splitmix64 : public vRNG64 {
					Uint64 state;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class thrustAlt64 : public vRNG64 {
					Uint64 state;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mulberry32 : public vRNG32 {
					Uint32 j;
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class molerat64 : public vRNG64 {
					Uint64 a;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
					//void seed(Uint64 seed);
				};
				class moverCounter64 : public vRNG64 {
					Uint64 a, b;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
					//void seed(Uint64 seed);
				};
				class moverCounter32 : public vRNG32 {
					Uint32 a, b;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
					//void seed(Uint64 seed);
				};
				class xrsr_rev_mul : public vRNG64 {
					Uint64 state0, state1;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
					//void seed(Uint64 seed);
				};
				class twirl32 : public vRNG32 {
					Uint32 a;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class zig32 : public vRNG32 {
					Uint32 a, b;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class twinLinear : public vRNG64 {
					Uint64 s0, s1;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class moremur64 : public vRNG64 {
					Uint64 state;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				// From https://gitlab.com/pomma89/troschuetz-random/-/blob/main/src/Troschuetz.Random/Generators/NR3Generator.cs
				// Has a long period, and just a few "unusual" anomalies early on.
				class nr3 : public vRNG64 {
					Uint64 u, v, w;
				public:
					const Uint64 SeedU1 = 2862933555777941757UL;
					const Uint64 SeedU2 = 7046029254386353087UL;
					const Uint64 SeedU3 = 4294957665UL;
					const Uint64 SeedV = 4101842887655102017UL;
					const Uint64 SeedW = 1UL;

					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				// From https://gitlab.com/pomma89/troschuetz-random/-/blob/main/src/Troschuetz.Random/Generators/NR3Q1Generator.cs
				// This is just a Marsaglia 64-bit xorshift generator with a multiplication at the end. It fails BRank horribly all-around, but nothing else early on.
				class nr3q1 : public vRNG64 {
					Uint64 v;
				public:
					const Uint64 SeedU = 2862933555777941757UL;
					const Uint64 SeedV = 4101842887655102017UL;

					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				// From https://gitlab.com/pomma89/troschuetz-random/-/blob/main/src/Troschuetz.Random/Generators/NR3Q2Generator.cs
				// Not much improvement over nr3q1, if any. Fails from the start on BRank.
				class nr3q2 : public vRNG64 {
					Uint64 v, w;
				public:
					const Uint64 SeedU = 4294957665UL;
					const Uint64 SeedV = 4101842887655102017UL;
					const Uint64 SeedW = 1UL;

					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class nova : public vRNG64 {
					Uint64 v, w;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mars256 : public vRNG64 {
					Uint64 stateA, stateB, stateC, stateD;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class marsgwt : public vRNG64 {
					Uint32 stateA, stateB, stateC, stateD;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class lizard256 : public vRNG64 {
					Uint64 stateA, stateB, stateC, stateD;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class plum256 : public vRNG64 {
					Uint64 stateA = 0, stateB = 1, stateC = 2, stateD = 3;
					int rotation;
				public:
					plum256(int rot) : rotation(rot) {};
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class overload320 : public vRNG64 {
					Uint64 stateA, stateB, stateC, stateD, stateE, stream;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class joker40 : public vRNG8 {
					Uint8 stateA, stateB, stateC, stateD, stateE;
				public:
					Uint8 raw8();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class suit80 : public vRNG16 {
					Uint16 stateA, stateB, stateC, stateD, stateE;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class floatHax32 : public vRNG16 {
					Uint32 stateA, stateB;
				public:
					Uint32 rand();
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject*);
				};

				class floatHax64 : public vRNG16 {
					Uint64 state;
				public:
					Uint64 rand();
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject*);
				};

				class floatHaxPCG : public vRNG16 {
					Uint64 state;
				public:
					Uint32 rand();
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject*);
				};

				class doubleHaxFlow : public vRNG32 {
					Uint64 stateA, stateB;
				public:
					Uint64 rand();
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject*);
				};

			}
		}
	}
}
