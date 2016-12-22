#pragma once

// This file is a compilation of randomc.h/cpp and smft.h/cpp into a single header file done by Arech for NNTL project. Some unnecessary
// for NNTL stuff has been removed(commented out) or simplified.
// Warning 6385 has been disabled for SMFT::BRandom() function.
// Corrected SMFT::MotherBits() to make MSVC run-time checker happy.
// All rights belongs to respective authors of original code.
// 
// For smft.h MEXP is (unjustifiedly) chosen to be #define MEXP 19937 instead of default #define MEXP 11213.

#include <stdint.h>
#include <emmintrin.h>                 // Define SSE2 intrinsics

namespace AFog {

#define AF_INT64_SUPPORTED

	constexpr int GEN_ERROR = 0x80000000;//added by Arech

	//////////////////////////////////////////////////////////////////////////

	/*****************************   randomc.h   **********************************
	* Author:        Agner Fog
	* Date created:  1997
	* Last modified: 2008-11-16
	* Project:       randomc.h
	* Source URL:    www.agner.org/random
	*
	* Description:
	* This header file contains class declarations and other definitions for the
	* randomc class library of uniform random number generators in C++ language.
	*
	* Overview of classes:
	* ====================
	*
	* class CRandomMersenne:
	* Random number generator of type Mersenne twister.
	* Source file mersenne.cpp
	*
	* class CRandomMother:
	* Random number generator of type Mother-of-All (Multiply with carry).
	* Source file mother.cpp
	*
	* class CRandomSFMT:
	* Random number generator of type SIMD-oriented Fast Mersenne Twister.
	* The class definition is not included here because it is not
	* portable to all platforms. See sfmt.h and sfmt.cpp for details.
	*
	* Member functions (methods):
	* ===========================
	*
	* All these classes have identical member functions:
	*
	* Constructor(int seed):
	* The seed can be any integer. The time may be used as seed.
	* Executing a program twice with the same seed will give the same sequence
	* of random numbers. A different seed will give a different sequence.
	*
	* void RandomInit(int seed);
	* Re-initializes the random number generator with a new seed.
	*
	* void RandomInitByArray(int const seeds[], int NumSeeds);
	* In CRandomMersenne and CRandomSFMT only: Use this function if you want
	* to initialize with a seed with more than 32 bits. All bits in the seeds[]
	* array will influence the sequence of random numbers generated. NumSeeds
	* is the number of entries in the seeds[] array.
	*
	* double Random();
	* Gives a floating point random number in the interval 0 <= x < 1.
	* The resolution is 32 bits in CRandomMother and CRandomMersenne, and
	* 52 bits in CRandomSFMT.
	*
	* int IRandom(int min, int max);
	* Gives an integer random number in the interval min <= x <= max.
	* (max-min < MAXINT).
	* The precision is 2^-32 (defined as the difference in frequency between
	* possible output values). The frequencies are exact if max-min+1 is a
	* power of 2.
	*
	* int IRandomX(int min, int max);
	* Same as IRandom, but exact. In CRandomMersenne and CRandomSFMT only.
	* The frequencies of all output values are exactly the same for an
	* infinitely long sequence. (Only relevant for extremely long sequences).
	*
	* uint32_t BRandom();
	* Gives 32 random bits.
	*
	*
	* Example:
	* ========
	* The file EX-RAN.CPP contains an example of how to generate random numbers.
	*
	*
	* Library version:
	* ================
	* Optimized versions of these random number generators are provided as function
	* libraries in randoma.zip. These function libraries are coded in assembly
	* language and support only x86 platforms, including 32-bit and 64-bit
	* Windows, Linux, BSD, Mac OS-X (Intel based). Use randoma.h from randoma.zip
	*
	*
	* Non-uniform random number generators:
	* =====================================
	* Random number generators with various non-uniform distributions are
	* available in stocc.zip (www.agner.org/random).
	*
	*
	* Further documentation:
	* ======================
	* The file ran-instructions.pdf contains further documentation and
	* instructions for these random number generators.
	*
	* Copyright 1997-2008 by Agner Fog.
	* GNU General Public License http://www.gnu.org/licenses/gpl.html
	*******************************************************************************/
	/**************************   mersenne.cpp   **********************************
	* Author:        Agner Fog
	* Date created:  2001
	* Last modified: 2008-11-16
	* Project:       randomc.h
	* Platform:      Any C++
	* Description:
	* Random Number generator of type 'Mersenne Twister'
	*
	* This random number generator is described in the article by
	* M. Matsumoto & T. Nishimura, in:
	* ACM Transactions on Modeling and Computer Simulation,
	* vol. 8, no. 1, 1998, pp. 3-30.
	* Details on the initialization scheme can be found at
	* http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
	*
	* Further documentation:
	* The file ran-instructions.pdf contains further documentation and
	* instructions.
	*
	* Copyright 2001-2008 by Agner Fog.
	* GNU General Public License http://www.gnu.org/licenses/gpl.html
	*******************************************************************************/
	/**************************   mother.cpp   ************************************
	* Author:        Agner Fog
	* Date created:  1999
	* Last modified: 2008-11-16
	* Project:       randomc.h
	* Platform:      This implementation uses 64-bit integers for intermediate calculations.
	*                Works only on compilers that support 64-bit integers.
	* Description:
	* Random Number generator of type 'Mother-Of-All generator'.
	*
	* This is a multiply-with-carry type of random number generator
	* invented by George Marsaglia.  The algorithm is:
	* S = 2111111111*X[n-4] + 1492*X[n-3] + 1776*X[n-2] + 5115*X[n-1] + C
	* X[n] = S modulo 2^32
	* C = floor(S / 2^32)
	*
	* Further documentation:
	* The file ran-instructions.pdf contains further documentation and
	* instructions.
	*
	* Copyright 1999-2008 by Agner Fog.
	* GNU General Public License http://www.gnu.org/licenses/gpl.html
	******************************************************************************/

#ifndef RANDOMC_H
#define RANDOMC_H

	// Define integer types with known size: int32_t, uint32_t, int64_t, uint64_t.
	// If this doesn't work then insert compiler-specific definitions here:
/*
#if defined(__GNUC__) || (defined(_MSC_VER) && _MSC_VER >= 1600)
// Compilers supporting C99 or C++0x have stdint.h defining these integer types
#include <stdint.h>
#define AF_INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#elif defined(_WIN16) || defined(__MSDOS__) || defined(_MSDOS) 
// 16 bit systems use long int for 32 bit integer.
	typedef   signed long int int32_t;
	typedef unsigned long int uint32_t;
#elif defined(_MSC_VER)
// Older Microsoft compilers have their own definition
	typedef   signed __int32  int32_t;
	typedef unsigned __int32 uint32_t;
	typedef   signed __int64  int64_t;
	typedef unsigned __int64 uint64_t;
#define AF_INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#else
// This works with most compilers
	typedef signed int          int32_t;
	typedef unsigned int       uint32_t;
	typedef long long           int64_t;
	typedef unsigned long long uint64_t;
#define AF_INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#endif*/


	/***********************************************************************
	System-specific user interface functions
	***********************************************************************/

	//void EndOfProgram(void);               // System-specific exit code (userintf.cpp)

	//void FatalError(const char *ErrorText);// System-specific error reporting (userintf.cpp)

#if defined(__cplusplus)               // class definitions only in C++
									   /***********************************************************************
									   Define random number generator classes
									   ***********************************************************************/

	class CRandomMersenne {                // Encapsulate random number generator
										   // Choose which version of Mersenne Twister you want:
#if 0 
									   // Define constants for type MT11213A:
#define MERS_N   351
#define MERS_M   175
#define MERS_R   19
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   17
#define MERS_A   0xE4BD75F5
#define MERS_B   0x655E5280
#define MERS_C   0xFFD58000
#else    
									   // or constants for type MT19937:
#define MERS_N   624
#define MERS_M   397
#define MERS_R   31
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   18
#define MERS_A   0x9908B0DF
#define MERS_B   0x9D2C5680
#define MERS_C   0xEFC60000
#endif

	public:
		CRandomMersenne(int seed) {         // Constructor
			RandomInit(seed); LastInterval = 0;
		}
		void RandomInit(int seed)          // Re-seed
		{
			// Initialize and seed
			Init0(seed);

			// Randomize some more
			for (int i = 0; i < 37; i++) BRandom();
		}


		void RandomInitByArray(int const seeds[], int NumSeeds) // Seed by more than 32 bits
		{
			// Seed by more than 32 bits
			int i, j, k;

			// Initialize
			Init0(19650218);

			if (NumSeeds <= 0) return;

			// Randomize mt[] using whole seeds[] array
			i = 1;  j = 0;
			k = (MERS_N > NumSeeds ? MERS_N : NumSeeds);
			for (; k; k--) {
				mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525UL)) + (uint32_t)seeds[j] + j;
				i++; j++;
				if (i >= MERS_N) { mt[0] = mt[MERS_N - 1]; i = 1; }
				if (j >= NumSeeds) j = 0;
			}
			for (k = MERS_N - 1; k; k--) {
				mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941UL)) - i;
				if (++i >= MERS_N) { mt[0] = mt[MERS_N - 1]; i = 1; }
			}
			mt[0] = 0x80000000UL;  // MSB is 1; assuring non-zero initial array

								   // Randomize some more
			mti = 0;
			for (int i = 0; i <= MERS_N; i++) BRandom();
		}

		int IRandom(int min, int max)     // Output random integer
		{
			// Output random integer in the interval min <= x <= max
			// Relative error on frequencies < 2^-32
			if (max <= min) {
				if (max == min) return min; else return 0x80000000;
			}
			// Multiply interval with random and truncate
			int r = int((double)(uint32_t)(max - min + 1) * Random() + min);
			if (r > max) r = max;
			return r;
		}


		int IRandomX(int min, int max)     // Output random integer, exact
		{
			// Output random integer in the interval min <= x <= max
			// Each output value has exactly the same probability.
			// This is obtained by rejecting certain bit values so that the number
			// of possible bit values is divisible by the interval length
			if (max <= min) {
				if (max == min) return min; else return 0x80000000;
			}
#ifdef  AF_INT64_SUPPORTED
			// 64 bit integers available. Use multiply and shift method
			uint32_t interval;                    // Length of interval
			uint64_t longran;                     // Random bits * interval
			uint32_t iran;                        // Longran / 2^32
			uint32_t remainder;                   // Longran % 2^32

			interval = uint32_t(max - min + 1);
			if (interval != LastInterval) {
				// Interval length has changed. Must calculate rejection limit
				// Reject when remainder >= 2^32 / interval * interval
				// RLimit will be 0 if interval is a power of 2. No rejection then
				RLimit = uint32_t(((uint64_t)1 << 32) / interval) * interval - 1;
				LastInterval = interval;
			}
			do { // Rejection loop
				longran = (uint64_t)BRandom() * interval;
				iran = (uint32_t)(longran >> 32);
				remainder = (uint32_t)longran;
			} while (remainder > RLimit);
			// Convert back to signed and return result
			return (int32_t)iran + min;

#else
			// 64 bit integers not available. Use modulo method
			uint32_t interval;                    // Length of interval
			uint32_t bran;                        // Random bits
			uint32_t iran;                        // bran / interval
			uint32_t remainder;                   // bran % interval

			interval = uint32_t(max - min + 1);
			if (interval != LastInterval) {
				// Interval length has changed. Must calculate rejection limit
				// Reject when iran = 2^32 / interval
				// We can't make 2^32 so we use 2^32-1 and correct afterwards
				RLimit = (uint32_t)0xFFFFFFFF / interval;
				if ((uint32_t)0xFFFFFFFF % interval == interval - 1) RLimit++;
			}
			do { // Rejection loop
				bran = BRandom();
				iran = bran / interval;
				remainder = bran % interval;
			} while (iran >= RLimit);
			// Convert back to signed and return result
			return (int32_t)remainder + min;

#endif
		}

		double Random()                    // Output random float
		{
			// Output random float number in the interval 0 <= x < 1
			// Multiply by 2^(-32)
			return (double)BRandom() * (1. / (65536.*65536.));
		}

		uint32_t BRandom()                 // Output random bits
		{
			// Generate 32 random bits
			uint32_t y;

			if (mti >= MERS_N) {
				// Generate MERS_N words at one time
				const uint32_t LOWER_MASK = (1LU << MERS_R) - 1;       // Lower MERS_R bits
				const uint32_t UPPER_MASK = 0xFFFFFFFF << MERS_R;      // Upper (32 - MERS_R) bits
				static const uint32_t mag01[2] = { 0, MERS_A };

				int kk;
				for (kk = 0; kk < MERS_N - MERS_M; kk++) {
					y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
					mt[kk] = mt[kk + MERS_M] ^ (y >> 1) ^ mag01[y & 1];
				}

				for (; kk < MERS_N - 1; kk++) {
					y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
					mt[kk] = mt[kk + (MERS_M - MERS_N)] ^ (y >> 1) ^ mag01[y & 1];
				}

				y = (mt[MERS_N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
				mt[MERS_N - 1] = mt[MERS_M - 1] ^ (y >> 1) ^ mag01[y & 1];
				mti = 0;
			}
			y = mt[mti++];

			// Tempering (May be omitted):
			y ^= y >> MERS_U;
			y ^= (y << MERS_S) & MERS_B;
			y ^= (y << MERS_T) & MERS_C;
			y ^= y >> MERS_L;

			return y;
		}

	private:
		void Init0(int seed)               // Basic initialization procedure
		{
			// Seed generator
			const uint32_t factor = 1812433253UL;
			mt[0] = seed;
			for (mti = 1; mti < MERS_N; mti++) {
				mt[mti] = (factor * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
			}
		}

		uint32_t mt[MERS_N];                // State vector
		int mti;                            // Index into mt
		uint32_t LastInterval;              // Last interval length for IRandomX
		uint32_t RLimit;                    // Rejection limit used by IRandomX
	};


	class CRandomMother {                  // Encapsulate random number generator
	public:
		void RandomInit(int seed)          // Initialization
		{
			int i;
			uint32_t s = seed;
			// make random numbers and put them into the buffer
			for (i = 0; i < 5; i++) {
				s = s * 29943829 - 1;
				x[i] = s;
			}
			// randomize some more
			for (i = 0; i < 19; i++) BRandom();
		}

		int IRandom(int min, int max)      // Get integer random number in desired interval
		{
			// Output random integer in the interval min <= x <= max
			// Relative error on frequencies < 2^-32
			if (max <= min) {
				if (max == min) return min; else return 0x80000000;
			}
			// Assume 64 bit integers supported. Use multiply and shift method
			uint32_t interval;                  // Length of interval
			uint64_t longran;                   // Random bits * interval
			uint32_t iran;                      // Longran / 2^32

			interval = (uint32_t)(max - min + 1);
			longran = (uint64_t)BRandom() * interval;
			iran = (uint32_t)(longran >> 32);
			// Convert back to signed and return result
			return (int32_t)iran + min;
		}

		double Random()                    // Get floating point random number
		{
			// returns a random number between 0 and 1:
			return (double)BRandom() * (1. / (65536.*65536.));
		}

		uint32_t BRandom()                 // Output random bits
		{
			uint64_t sum;
			sum = (uint64_t)2111111111UL * (uint64_t)x[3] +
				(uint64_t)1492 * (uint64_t)(x[2]) +
				(uint64_t)1776 * (uint64_t)(x[1]) +
				(uint64_t)5115 * (uint64_t)(x[0]) +
				(uint64_t)x[4];
			x[3] = x[2];  x[2] = x[1];  x[1] = x[0];
			x[4] = (uint32_t)(sum >> 32);                  // Carry
			x[0] = (uint32_t)sum;                          // Low 32 bits of sum
			return x[0];
		}

		CRandomMother(int seed) {           // Constructor
			RandomInit(seed);
		}
	protected:
		uint32_t x[5];                      // History buffer
	};

#endif // __cplusplus
#endif // RANDOMC_H



	/*****************************    sfmt.h    ***********************************
	* Authors:
	* Mutsuo Saito (Hiroshima University)
	* Makoto Matsumoto (Hiroshima University)
	* Agner Fog (Technical University of Denmark)
	* Date created:  2006
	* Last modified: 2009-02-08
	* Project:       randomc
	* Platform:      This C++ version requires an x86 family microprocessor
	*                with the SSE2 or later instruction set and a compiler
	*                that supports intrinsic functions.
	* Source URL:    www.agner.org/random
	*
	* Description:
	* This header file contains class declarations and other definitions for the
	* "SIMD-oriented Fast Mersenne Twister" (SFMT) random number generator.
	*
	* The SFMT random number generator is a modification of the Mersenne Twister
	* with improved randomness and speed, adapted to the SSE2 instruction set.
	* The SFMT was invented by Mutsuo Saito and Makoto Matsumoto.
	* The present C++ implementation is by Agner Fog.
	*
	* Class description:
	* ==================
	* class CRandomSFMT:
	* Random number generator of type SIMD-oriented Fast Mersenne Twister.
	*
	* Member functions (methods):
	* ===========================
	* Constructor CRandomSFMT(int seed, int IncludeMother = 1):
	* The seed can be any integer.
	* Executing a program twice with the same seed will give the same sequence of
	* random numbers. A different seed will give a different sequence.
	* If IncludeMother is 1 then the output of the SFMT generator is combined
	* with the output of the Mother-Of-All generator. The combined output has an
	* excellent randomness that has passed very stringent tests for randomness.
	* If IncludeMother is 0 then the SFMT generator is used alone.
	*
	* void RandomInit(int seed);
	* Re-initializes the random number generator with a new seed.
	*
	* void RandomInitByArray(int seeds[], int NumSeeds);
	* Use this function if you want to initialize with a seed with more than
	* 32 bits. All bits in the seeds[] array will influence the sequence of
	* random numbers generated. NumSeeds is the number of entries in the seeds[]
	* array.
	*
	* double Random();
	* Gives a floating point random number in the interval 0 <= x < 1.
	* The resolution is 52 bits.
	*
	* int IRandom(int min, int max);
	* Gives an integer random number in the interval min <= x <= max.
	* (max-min < MAXINT).
	* The precision is 2^(-32) (defined as the difference in frequency between
	* possible output values). The frequencies are exact if (max-min+1) is a
	* power of 2.
	*
	* int IRandomX(int min, int max);
	* Same as IRandom, but exact. The frequencies of all output values are
	* exactly the same for an infinitely long sequence.
	*
	* uint32_t BRandom();
	* Gives 32 random bits.
	*
	*
	* Example:
	* ========
	* The file EX-RAN.CPP contains an example of how to generate random numbers.
	*
	* Library version:
	* ================
	* An optimized version of this random number generator is provided as function
	* libraries in randoma.zip. These function libraries are coded in assembly
	* language and support only x86 platforms, including 32-bit and 64-bit
	* Windows, Linux, BSD, Mac OS-X (Intel based). Use randoma.h from randoma.zip
	*
	*
	* Non-uniform random number generators:
	* =====================================
	* Random number generators with various non-uniform distributions are available
	* in stocc.zip (www.agner.org/random).
	*
	*
	* Further documentation:
	* ======================
	* The file ran-instructions.pdf contains further documentation and
	* instructions for these random number generators.
	*
	*
	* Copyright notice
	* ================
	* GNU General Public License http://www.gnu.org/licenses/gpl.html
	* This C++ implementation of SFMT contains parts of the original C code
	* which was published under the following BSD license, which is therefore
	* in effect in addition to the GNU General Public License.
	*
	Copyright (c) 2006, 2007 by Mutsuo Saito, Makoto Matsumoto and Hiroshima University.
	Copyright (c) 2008 by Agner Fog.
	All rights reserved.
	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are met:
	> Redistributions of source code must retain the above copyright notice,
	this list of conditions and the following disclaimer.
	> Redistributions in binary form must reproduce the above copyright notice,
	this list of conditions and the following disclaimer in the documentation
	and/or other materials provided with the distribution.
	> Neither the name of the Hiroshima University nor the names of its
	contributors may be used to endorse or promote products derived from
	this software without specific prior written permission.
	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
	"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
	LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
	A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
	OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
	SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
	LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
	THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
	OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	*******************************************************************************/
	/*****************************   sfmt.cpp   ***********************************
	* Authors:
	* Mutsuo Saito (Hiroshima University)
	* Makoto Matsumoto (Hiroshima University)
	* Agner Fog (Technical University of Denmark)
	* Date created:  2006
	* Last modified: 2009-02-08
	* Project:       randomc
	* Platform:      This C++ version requires an x86 family microprocessor
	*                with the SSE2 or later instruction set and a compiler
	*                that supports intrinsic functions.
	* Source URL:    www.agner.org/random
	* Source URL for original C language implementation:
	*                www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/index.html
	*
	* Description:
	* "SIMD-oriented Fast Mersenne Twister" (SFMT) random number generator.
	* The SFMT random number generator is a modification of the Mersenne Twister
	* with improved randomness and speed, adapted to the SSE2 instruction set.
	* The SFMT was invented by Mutsuo Saito and Makoto Matsumoto.
	* The present C++ implementation is by Agner Fog.
	*
	* Class description and member functions: See sfmt.h
	*
	* Example:
	* ========
	* The file EX-RAN.CPP contains an example of how to generate random numbers.
	*
	* Library version:
	* ================
	* An optimized version of this random number generator is provided as function
	* libraries in randoma.zip. These function libraries are coded in assembly
	* language and support only x86 platforms, including 32-bit and 64-bit
	* Windows, Linux, BSD, Mac OS-X (Intel based). Use randoma.h from randoma.zip
	*
	*
	* Further documentation:
	* ======================
	* See the file ran-instructions.pdf for detailed instructions and documentation
	*
	*
	* Copyright notice
	* ================
	* GNU General Public License http://www.gnu.org/licenses/gpl.html
	* This C++ implementation of SFMT contains parts of the original C code
	* which was published under the following BSD license, which is therefore
	* in effect in addition to the GNU General Public License.
	*
	Copyright (c) 2006, 2007 by Mutsuo Saito, Makoto Matsumoto and Hiroshima University.
	Copyright (c) 2008 by Agner Fog.
	All rights reserved.
	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are met:
	> Redistributions of source code must retain the above copyright notice,
	this list of conditions and the following disclaimer.
	> Redistributions in binary form must reproduce the above copyright notice,
	this list of conditions and the following disclaimer in the documentation
	and/or other materials provided with the distribution.
	> Neither the name of the Hiroshima University nor the names of its
	contributors may be used to endorse or promote products derived from
	this software without specific prior written permission.
	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
	"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
	LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
	A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
	OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
	SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
	LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
	THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
	OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	*******************************************************************************/

#ifndef SFMT_H
#define SFMT_H

//#include <emmintrin.h>                 // Define SSE2 intrinsics
//#include "randomc.h"                   // Define integer types etc

	// Choose one of the possible Mersenne exponents.
	// Higher values give longer cycle length and use more memory:
	//#define MEXP   607
	//#define MEXP  1279
	//#define MEXP  2281
	//#define MEXP  4253
//#define MEXP 11213
	#define MEXP 19937
	//#define MEXP 44497

	// Define constants for the selected Mersenne exponent:
#if MEXP == 44497
#define SFMT_N    348                  // Size of state vector
#define SFMT_M    330                  // Position of intermediate feedback
#define SFMT_SL1    5                  // Left shift of W[N-1], 32-bit words
#define SFMT_SL2	  3                  // Left shift of W[0], *8, 128-bit words
#define SFMT_SR1    9                  // Right shift of W[M], 32-bit words
#define SFMT_SR2	  3                  // Right shift of W[N-2], *8, 128-bit words
#define SFMT_MASK	  0xeffffffb,0xdfbebfff,0xbfbf7bef,0x9ffd7bff // AND mask
#define SFMT_PARITY 1,0,0xa3ac4000,0xecc1327a   // Period certification vector

#elif MEXP == 19937
#define SFMT_N    156                  // Size of state vector
#define SFMT_M    122                  // Position of intermediate feedback
#define SFMT_SL1   18                  // Left shift of W[N-1], 32-bit words
#define SFMT_SL2	  1                  // Left shift of W[0], *8, 128-bit words
#define SFMT_SR1   11                  // Right shift of W[M], 32-bit words
#define SFMT_SR2	  1                  // Right shift of W[N-2], *8, 128-bit words
#define SFMT_MASK	  0xdfffffef,0xddfecb7f,0xbffaffff,0xbffffff6 // AND mask
#define SFMT_PARITY 1,0,0,0x13c9e684   // Period certification vector

#elif MEXP == 11213
#define SFMT_N    88                   // Size of state vector
#define SFMT_M    68                   // Position of intermediate feedback
#define SFMT_SL1	14                   // Left shift of W[N-1], 32-bit words
#define SFMT_SL2	 3                   // Left shift of W[0], *8, 128-bit words
#define SFMT_SR1	 7                   // Right shift of W[M], 32-bit words
#define SFMT_SR2	 3                   // Right shift of W[N-2], *8, 128-bit words
#define SFMT_MASK	 0xeffff7fb,0xffffffef,0xdfdfbfff,0x7fffdbfd // AND mask
#define SFMT_PARITY 1,0,0xe8148000,0xd0c7afa3 // Period certification vector

#elif MEXP == 4253
#define SFMT_N    34                   // Size of state vector
#define SFMT_M    17                   // Position of intermediate feedback
#define SFMT_SL1	20                   // Left shift of W[N-1], 32-bit words
#define SFMT_SL2	 1                   // Left shift of W[0], *8, 128-bit words
#define SFMT_SR1	 7                   // Right shift of W[M], 32-bit words
#define SFMT_SR2	 1                   // Right shift of W[N-2], *8, 128-bit words
#define SFMT_MASK	 0x9f7bffff, 0x9fffff5f, 0x3efffffb, 0xfffff7bb // AND mask
#define SFMT_PARITY 0xa8000001, 0xaf5390a3, 0xb740b3f8, 0x6c11486d // Period certification vector

#elif MEXP == 2281
#define SFMT_N    18                   // Size of state vector
#define SFMT_M    12                   // Position of intermediate feedback
#define SFMT_SL1	19                   // Left shift of W[N-1], 32-bit words
#define SFMT_SL2	 1                   // Left shift of W[0], *8, 128-bit words
#define SFMT_SR1	 5                   // Right shift of W[M], 32-bit words
#define SFMT_SR2	 1                   // Right shift of W[N-2], *8, 128-bit words
#define SFMT_MASK	 0xbff7ffbf, 0xfdfffffe, 0xf7ffef7f, 0xf2f7cbbf // AND mask
#define SFMT_PARITY 0x00000001, 0x00000000, 0x00000000, 0x41dfa600  // Period certification vector

#elif MEXP == 1279
#define SFMT_N    10                   // Size of state vector
#define SFMT_M     7                   // Position of intermediate feedback
#define SFMT_SL1	14                   // Left shift of W[N-1], 32-bit words
#define SFMT_SL2	 3                   // Left shift of W[0], *8, 128-bit words
#define SFMT_SR1	 5                   // Right shift of W[M], 32-bit words
#define SFMT_SR2	 1                   // Right shift of W[N-2], *8, 128-bit words
#define SFMT_MASK	  0xf7fefffd, 0x7fefcfff, 0xaff3ef3f, 0xb5ffff7f  // AND mask
#define SFMT_PARITY 0x00000001, 0x00000000, 0x00000000, 0x20000000  // Period certification vector

#elif MEXP == 607
#define SFMT_N     5                   // Size of state vector
#define SFMT_M     2                   // Position of intermediate feedback
#define SFMT_SL1	15                   // Left shift of W[N-1], 32-bit words
#define SFMT_SL2	 3                   // Left shift of W[0], *8, 128-bit words
#define SFMT_SR1	13                   // Right shift of W[M], 32-bit words
#define SFMT_SR2	 3                   // Right shift of W[N-2], *8, 128-bit words
#define SFMT_MASK	  0xfdff37ff, 0xef7f3f7d, 0xff777b7d, 0x7ff7fb2f  // AND mask
#define SFMT_PARITY 0x00000001, 0x00000000, 0x00000000, 0x5986f054  // Period certification vector
#endif

	// Class for SFMT generator with or without Mother-Of-All generator
	class CRandomSFMT {                              // Encapsulate random number generator
	public:
		CRandomSFMT(int seed, int IncludeMother = 0) {// Constructor
			UseMother = IncludeMother;
			LastInterval = 0;
			RandomInit(seed);
		}
		void RandomInit(int seed)                    // Re-seed
		{
			// Re-seed
			uint32_t i;                         // Loop counter
			uint32_t y = seed;                  // Temporary
			uint32_t statesize = SFMT_N * 4;      // Size of state vector
			if (UseMother) statesize += 5;      // Add states for Mother-Of-All generator

												// Fill state vector with random numbers from seed
			((uint32_t*)state)[0] = y;
			const uint32_t factor = 1812433253U;// Multiplication factor

			for (i = 1; i < statesize; i++) {
				y = factor * (y ^ (y >> 30)) + i;
				((uint32_t*)state)[i] = y;
			}

			// Further initialization and period certification
			Init2();
		}

		void RandomInitByArray(int const seeds[], int NumSeeds) // Seed by more than 32 bits
		{
			// Seed by more than 32 bits
			uint32_t i, j, count, r, lag;

			if (NumSeeds < 0) NumSeeds = 0;

			const uint32_t size = SFMT_N * 4; // number of 32-bit integers in state

											  // Typecast state to uint32_t *
			uint32_t * sta = (uint32_t*)state;

			if (size >= 623) {
				lag = 11;
			} else if (size >= 68) {
				lag = 7;
			} else if (size >= 39) {
				lag = 5;
			} else {
				lag = 3;
			}
			const uint32_t mid = (size - lag) / 2;

			if ((uint32_t)NumSeeds + 1 > size) {
				count = (uint32_t)NumSeeds;
			} else {
				count = size - 1;
			}
#if 0
			// Original code. Argument to func1 is constant!
			for (i = 0; i < size; i++) sta[i] = 0x8B8B8B8B;
			r = func1(sta[0] ^ sta[mid] ^ sta[size - 1]);
			sta[mid] += r;
			r += NumSeeds;
			sta[mid + lag] += r;
			sta[0] = r;
#else
			// 1. loop: Fill state vector with random numbers from NumSeeds
			const uint32_t factor = 1812433253U;// Multiplication factor
			r = (uint32_t)NumSeeds;
			for (i = 0; i < SFMT_N * 4; i++) {
				r = factor * (r ^ (r >> 30)) + i;
				sta[i] = r;
			}
#endif

			// 2. loop: Fill state vector with random numbers from seeds[]
			for (i = 1, j = 0; j < count; j++) {
				r = func1(sta[i] ^ sta[(i + mid) % size] ^ sta[(i + size - 1) % size]);
				sta[(i + mid) % size] += r;
				if (j < (uint32_t)NumSeeds) r += (uint32_t)seeds[j];
				r += i;
				sta[(i + mid + lag) % size] += r;
				sta[i] = r;
				i = (i + 1) % size;
			}

			// 3. loop: Randomize some more
			for (j = 0; j < size; j++) {
				r = func2(sta[i] + sta[(i + mid) % size] + sta[(i + size - 1) % size]);
				sta[(i + mid) % size] ^= r;
				r -= i;
				sta[(i + mid + lag) % size] ^= r;
				sta[i] = r;
				i = (i + 1) % size;
			}
			if (UseMother) {
				// 4. loop: Initialize MotherState
				for (j = 0; j < 5; j++) {
					r = func2(r) + j;
					MotherState[j] = r + sta[2 * j];
				}
			}

			// Further initialization and period certification
			Init2();
		}

		int  IRandom(int min, int max)             // Output random integer
		{
			// Output random integer in the interval min <= x <= max
			// Slightly inaccurate if (max-min+1) is not a power of 2
			if (max <= min) {
				if (max == min) return min; else return 0x80000000;
			}
			// Assume 64 bit integers supported. Use multiply and shift method
			uint32_t interval;                  // Length of interval
			uint64_t longran;                   // Random bits * interval
			uint32_t iran;                      // Longran / 2^32

			interval = (uint32_t)(max - min + 1);
			longran = (uint64_t)BRandom() * interval;
			iran = (uint32_t)(longran >> 32);
			// Convert back to signed and return result
			return (int32_t)iran + min;
		}


		int  IRandomX(int min, int max)             // Output random integer, exact
		{
			// Output random integer in the interval min <= x <= max
			// Each output value has exactly the same probability.
			// This is obtained by rejecting certain bit values so that the number
			// of possible bit values is divisible by the interval length
			if (max <= min) {
				if (max == min) {
					return min;                   // max == min. Only one possible value
				} else {
					return 0x80000000;            // max < min. Error output
				}
			}
			// Assume 64 bit integers supported. Use multiply and shift method
			uint32_t interval;                  // Length of interval
			uint64_t longran;                   // Random bits * interval
			uint32_t iran;                      // Longran / 2^32
			uint32_t remainder;                 // Longran % 2^32

			interval = (uint32_t)(max - min + 1);
			if (interval != LastInterval) {
				// Interval length has changed. Must calculate rejection limit
				// Reject when remainder = 2^32 / interval * interval
				// RLimit will be 0 if interval is a power of 2. No rejection then.
				RLimit = (uint32_t)(((uint64_t)1 << 32) / interval) * interval - 1;
				LastInterval = interval;
			}
			do { // Rejection loop
				longran = (uint64_t)BRandom() * interval;
				iran = static_cast<uint32_t>(longran >> 32);
				remainder = static_cast<uint32_t>(longran & 0xFFFFFFFF);
			} while (remainder > RLimit);
			// Convert back to signed and return result
			return static_cast<int32_t>(iran) + min;
		}

		double Random()                              // Output random floating point number
		{
			// Output random floating point number
			if (ix >= SFMT_N * 4 - 1) {
				// Make sure we have at least two 32-bit numbers
				Generate();
			}
			uint64_t r = *(uint64_t*)((uint32_t*)state + ix);
			ix += 2;
			if (UseMother) {
				// We need 53 bits from Mother-Of-All generator
				// Use the regular 32 bits and the the carry bits rotated
				uint64_t r2 = (uint64_t)MotherBits() << 32;
				r2 |= (MotherState[4] << 16) | (MotherState[4] >> 16);
				r += r2;
			}
			// 53 bits resolution:
			// return (int64_t)(r >> 11) * (1./(67108864.0*134217728.0)); // (r >> 11)*2^(-53)
			// 52 bits resolution for compatibility with assembly version:
			return (int64_t)(r >> 12) * (1. / (67108864.0*67108864.0));  // (r >> 12)*2^(-52)
		}

#pragma warning(disable: 6385)
		uint32_t BRandom()                           // Output random bits
		{
			// Output 32 random bits
			uint32_t y;

			if (ix >= SFMT_N * 4) {
				Generate();
			}
			y = ((uint32_t*)state)[ix++];
			if (UseMother) y += MotherBits();
			return y;
		}
#pragma warning(default: 6385)


	private:
		// Functions used by CRandomSFMT::RandomInitByArray
		static uint32_t func1(uint32_t x) {
			return (x ^ (x >> 27)) * 1664525U;
		}

		static uint32_t func2(uint32_t x) {
			return (x ^ (x >> 27)) * 1566083941U;
		}

		// Subfunction for the sfmt algorithm
		static inline __m128i sfmt_recursion(__m128i const &a, __m128i const &b,
			__m128i const &c, __m128i const &d, __m128i const &mask) {
			__m128i a1, b1, c1, d1, z1, z2;
			b1 = _mm_srli_epi32(b, SFMT_SR1);
			a1 = _mm_slli_si128(a, SFMT_SL2);
			c1 = _mm_srli_si128(c, SFMT_SR2);
			d1 = _mm_slli_epi32(d, SFMT_SL1);
			b1 = _mm_and_si128(b1, mask);
			z1 = _mm_xor_si128(a, a1);
			z2 = _mm_xor_si128(b1, d1);
			z1 = _mm_xor_si128(z1, c1);
			z2 = _mm_xor_si128(z1, z2);
			return z2;
		}

		void Init2()                                 // Various initializations and period certification
		{
			// Various initializations and period certification
			uint32_t i, j, temp;

			// Initialize mask
			static const uint32_t maskinit[4] = { SFMT_MASK };
			mask = _mm_loadu_si128((__m128i*)maskinit);

			// Period certification
			// Define period certification vector
			static const uint32_t parityvec[4] = { SFMT_PARITY };

			// Check if parityvec & state[0] has odd parity
			temp = 0;
			for (i = 0; i < 4; i++) {
				temp ^= parityvec[i] & ((uint32_t*)state)[i];
			}
			for (i = 16; i > 0; i >>= 1) temp ^= temp >> i;
			if (!(temp & 1)) {
				// parity is even. Certification failed
				// Find a nonzero bit in period certification vector
				for (i = 0; i < 4; i++) {
					if (parityvec[i]) {
						for (j = 1; j; j <<= 1) {
							if (parityvec[i] & j) {
								// Flip the corresponding bit in state[0] to change parity
								((uint32_t*)state)[i] ^= j;
								// Done. Exit i and j loops
								i = 5;  break;
							}
						}
					}
				}
			}
			// Generate first random numbers and set ix = 0
			Generate();
		}

		void Generate()                              // Fill state array with new random numbers
		{
			// Fill state array with new random numbers
			int i;
			__m128i r, r1, r2;

			r1 = state[SFMT_N - 2];
			r2 = state[SFMT_N - 1];
			for (i = 0; i < SFMT_N - SFMT_M; i++) {
				r = sfmt_recursion(state[i], state[i + SFMT_M], r1, r2, mask);
				state[i] = r;
				r1 = r2;
				r2 = r;
			}
			for (; i < SFMT_N; i++) {
				r = sfmt_recursion(state[i], state[i + SFMT_M - SFMT_N], r1, r2, mask);
				state[i] = r;
				r1 = r2;
				r2 = r;
			}
			ix = 0;
		}


		uint32_t MotherBits()                        // Get random bits from Mother-Of-All generator
		{
			// Get random bits from Mother-Of-All generator
			uint64_t sum;
			sum =
				(uint64_t)2111111111U * (uint64_t)MotherState[3] +
				(uint64_t)1492 * (uint64_t)MotherState[2] +
				(uint64_t)1776 * (uint64_t)MotherState[1] +
				(uint64_t)5115 * (uint64_t)MotherState[0] +
				(uint64_t)MotherState[4];
			MotherState[3] = MotherState[2];
			MotherState[2] = MotherState[1];
			MotherState[1] = MotherState[0];
			MotherState[4] = (uint32_t)(sum >> 32);       // Carry
			MotherState[0] = (uint32_t) (sum & 0xffffffff);               // Low 32 bits of sum
			return MotherState[0];
		}


		uint32_t ix;                                  // Index into state array
		uint32_t LastInterval;                        // Last interval length for IRandom
		uint32_t RLimit;                              // Rejection limit used by IRandom
		uint32_t UseMother;                           // Combine with Mother-Of-All generator
		__m128i  mask;                                // AND mask
		__m128i  state[SFMT_N];                       // State vector for SFMT generator
		uint32_t MotherState[5];                      // State vector for Mother-Of-All generator
	};

	// Class for SFMT generator without Mother-Of-All generator
	// Derived from CRandomSFMT
	class CRandomSFMT0 : public CRandomSFMT {
	public:
		CRandomSFMT0(int seed) : CRandomSFMT(seed, 0) {}
	};

	// Class for SFMT generator combined with Mother-Of-All generator
	// Derived from CRandomSFMT
	class CRandomSFMT1 : public CRandomSFMT {
	public:
		CRandomSFMT1(int seed) : CRandomSFMT(seed, 1) {}
	};

#endif // SFMT_H


}