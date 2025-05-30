#include <immintrin.h>

/*
RNGs in the mediocre directory are not intended for real world use
only for research; as such they may get pretty sloppy in some areas

This set is of RNGs that do not make any significant use of:
	multiplication/division, arrays, flow control, complex math functions
*/

namespace PractRand {
	namespace RNGs {
		namespace Polymorphic {
			namespace NotRecommended {
				class xsalta16x3 : public vRNG16 {
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xsaltb16x3 : public vRNG16 {
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xsaltc16x3 : public vRNG16 {
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				//xorshift RNGs, a subset of LFSRs proposed by Marsaglia in 2003
				class xorshift32 : public vRNG32 {
					//constants those Marsaglia described as "one of my favorites" on page 4 of his 2003 paper
					Uint32 a;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xorshift64 : public vRNG64 {
					//constants those Marsaglia used in his sample code on page 4 of his 2003 paper
					Uint64 a;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xorshift128plus : public vRNG64 {
					Uint64 seed0, seed1;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xorshift64of128 : public vRNG64 {
					//the constants are still in need of tuning
					Uint64 high, low;
					void xls(int bits);
					void xrs(int bits);
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xorshift32of128 : public vRNG32 {
					xorshift64of128 impl;
				public:
					Uint32 raw32() {return Uint32(impl.raw64());}
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xorshift16of32 : public vRNG16 {
					xorshift32 impl;
				public:
					Uint16 raw16() {return Uint16(impl.raw32());}
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xorshift32of64 : public vRNG32 {
					xorshift64 impl;
				public:
					Uint32 raw32() {return Uint32(impl.raw64());}
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class xorshift32x4 : public vRNG32 {
					//recommended at the top of Marsaglias 2003 xorshift paper
					Uint32 x,y,z,w;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xorwow32of96 : public vRNG32 {
					xorshift64 impl;
					Uint32 a;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xorwow32x6 : public vRNG32 {
					//recommended at the top of Marsaglias 2003 xorshift paper
					Uint32 x,y,z,w,v,d;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xoroshiro128plus : public vRNG64 {
					// from David Blackman and Sebastiano Vigna (vigna@acm.org), see http://vigna.di.unimi.it/xorshift/
					Uint64 state0, state1;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class xoroshiro128plus_2p64 : public vRNG64 {
					// as xoroshiro128plus, but it skips 2**64-1 outputs between each pair of outputs (testing its recommended parallel sequences)
					Uint64 state0, state1;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};


				//by Ilya O. Levin, see http://www.literatecode.com/2004/10/18/sapparot/
				class sapparot : public vRNG32 {
					Uint32 a, b;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				//variants of sapparot created for testing purposes
				class sap16of48 : public vRNG16 {
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class sap32of96 : public vRNG32 {
					Uint32 a, b, c;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				
				//the low quality variant of the FLEA RNG by Robert Jenkins
				//(he named this variant flea0)
				class flea32x1 : public vRNG32 {
					enum {SIZE = 1};
					Uint32 a[SIZE], b, c, d, i;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				//versions of my SFC algorithm, archived here
				class sfc_v1_16 : public vRNG16 {
					Uint16 a, b, counter;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class sfc_v1_32: public vRNG32 {
					Uint32 a, b, counter;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class sfc_v2_16 : public vRNG16 {
					Uint16 a, b, counter;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class sfc_v2_32 : public vRNG32 {
					Uint32 a, b, counter;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class sfc_v3_16 : public vRNG16 {
					Uint16 a, b, counter;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class sfc_v3_32 : public vRNG32 {
					Uint32 a, b, counter;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class jsf8 : public vRNG8 {
					Uint8 a, b, c, d;
				public:
					Uint8 raw8();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class jsf16 : public vRNG16 {
					Uint16 a, b, c, d;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				//a few simple RNGs just for testing purposes
				class simpleA : public vRNG32 {
					Uint32 a, b, c;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class simpleB : public vRNG16 {
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class simpleC : public vRNG16 {
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class simpleD : public vRNG32 {
					Uint32 a, b, c;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class simpleE : public vRNG32 {
					//seems like a good combination of speed & quality
					//but falls flat when used on 16 bit words (irreversible, statespace issues)
					Uint32 a, b, c;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class simpleF : public vRNG16 {
					Uint16 a, b, c, d;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class simpleG : public vRNG32 {
					Uint32 a, b, c, d;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class trivium_weakenedA : public vRNG32 {
					Uint64 a, b;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class trivium_weakenedB : public vRNG16 {
					Uint32 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				//see http://www.drdobbs.com/tools/229625477
				class mo_Lesr32 : public vRNG32 {
					Uint32 state;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mo_ResrRers32 : public vRNG32 {
					Uint32 a;
					Uint32 b;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mo_Rers32of64 : public vRNG32 {
					Uint64 state;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mo_Resr32of64 : public vRNG32 {
					Uint64 state;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mo_Resdra32of64 : public vRNG32 {
					Uint64 state;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class murmlacish : public vRNG32 {
					Uint32 state1, state2, state3;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class gjishA : public vRNG16 {
					Uint16 a, b, c;
				public:
					Uint16 raw16();//broadly similar to gjrand, but 16 bit and no counter
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class gjishB : public vRNG16 {
					Uint16 a, b, c, counter;
				public:
					Uint16 raw16();//broadly similar to gjrand, but 16 bit
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class gjishC : public vRNG32 {
					Uint32 a, b, c;
				public:
					Uint32 raw32();//broadly similar to gjrand, but 32 bit, no counter, and no 3rd shift
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class gjishD : public vRNG32 {
					Uint32 a, b, c;
				public:
					Uint32 raw32();//broadly similar to gjrand, but 32 bit, no counter
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class ara16 : public vRNG16 {//add, bit rotate, add
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class ara32 : public vRNG32 {//add, bit rotate, add
					Uint32 a, b, c;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class arx16 : public vRNG16 {//add, bit rotate, xor
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class arx32 : public vRNG32 {//add, bit rotate, xor
					Uint32 a, b, c;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class hara16 : public vRNG16 {//heterogenous add, bit rotate, xor
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class harx16 : public vRNG16 {//heterogenous add, bit rotate, xor
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class learx16 : public vRNG16 {//LEAs, bit rotates, xor
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class hlearx16 : public vRNG16 {//heterogenous LEAs, bit rotates, xor
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class alearx16 : public vRNG16 {//add, LEA, bit rotate, xor
					Uint16 a, b, c;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class arac16 : public vRNG16 {//add, bit rotate, add (with counter)
					Uint16 a, b, c, counter;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class arxc16 : public vRNG16 {//add, bit rotate, xor (with counter)
					Uint16 a, b, c, counter;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class rarns16 : public vRNG16 {
					Uint16 xs1, xs2, xs3, history;
					enum { S1 = 3, S2 = 7, S3 = 8 };//3,7,8
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *walker);
					void seed(Uint64 s);
				};






				class xoshiro256starstar : public vRNG64 {
					// from David Blackman and Sebastiano Vigna (vigna@acm.org), see http://vigna.di.unimi.it/xorshift/
					Uint32 state0, state1, state2, state3;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class xoshiro256scramjet : public vRNG64 {
					// from David Blackman and Sebastiano Vigna (vigna@acm.org), see http://vigna.di.unimi.it/xorshift/
					const Uint64 first[4] = { 0x9E3779B97F4A7C15ULL, 0xC13FA9A902A6328FULL, 0x8CB92BA72F3D8DD7ULL, 0xABC98388FB8FAC02ULL };
					const Uint64 second[4] = { 0xDB4F0B9175AE2165ULL, 0xBBE0563303A4615FULL, 0xA0F2EC75A1FE1575ULL, 0x89E182857D9ED688ULL };
					Uint64 state[4];
					int current = 0;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class xoshiro4x32 : public vRNG32 {
					// from David Blackman and Sebastiano Vigna (vigna@acm.org), see http://vigna.di.unimi.it/xorshift/
					Uint32 state0, state1, state2, state3;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class xoshiro5x32 : public vRNG32 {
					// from David Blackman and Sebastiano Vigna (vigna@acm.org), see http://vigna.di.unimi.it/xorshift/
					Uint32 state0, state1, state2, state3, state4;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class xoshiro5x32long : public vRNG64 {
					// from David Blackman and Sebastiano Vigna (vigna@acm.org), see http://vigna.di.unimi.it/xorshift/
					Uint32 state0, state1, state2, state3, state4;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class mover32 : public vRNG32 {
					Uint32 a;
					Uint32 b;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class mover64 : public vRNG64 {
					Uint64 a;
					Uint64 b;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class multimover32 : public vRNG32 {
					Uint32 a;
					Uint32 b;
					Uint32 c;
					Uint32 d;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class thrust64a : public vRNG64 {
					Uint64 state0, state1;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
					//void seed(Uint64 seed);
				};
				class oriole64 : public vRNG64 {
					Uint64 state0, state1, state2;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class xoroshiro32_32 : public vRNG32 {
					Uint32 state0, state1, state2;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class xoroshiro16_16plus : public vRNG16 {
					// from David Blackman and Sebastiano Vigna (vigna@acm.org), see http://vigna.di.unimi.it/xorshift/
					Uint16 state0, state1;
				public:
					Uint16 raw16();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class cloud : public vRNG32 {
					uint8_t board[8];
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				private:
					uint8_t result();
				};

				class acorn64_10 : public vRNG64 {
					Uint64 a, b, c, d, e, f, g, h, i, j, stream;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};
				class chonk8x32 : public vRNG32 {
					Uint32 a, b, c, d, e, f, g, h;
				public:
					Uint32 raw32();
					std::string get_name() const;
					void walk_state(StateWalkingObject *);
				};

				class spangled_varqual : public vRNG64 {
					Uint64 stateA, stateB, stateC;
					int rounds;
				public:
					spangled_varqual(int rounds_) : rounds(rounds_) {}
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject*);
				};

				class chill96 : public vRNG64 {
					Uint32 stateA, stateB, stateC;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject*);
				};

				class choppy128 : public vRNG64 {
					Uint32 a;
					Uint32 b;
					Uint32 c;
					Uint32 d;
				public:
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject*);
				};


				class mace : public vRNG64 {
					Uint64 stateA, stateB, stateC, stateD, stateE, stream;
				public:
					mace(Uint64 stream_) : stream(_pdep_u64(stream_, 0x003569CA5369AC00UL) ^ 0x9E3779B97F4A7C15UL) {}
					Uint64 raw64();
					std::string get_name() const;
					void walk_state(StateWalkingObject*);
				};

			}
		}
	}
}
