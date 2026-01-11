var computeShindo = (function () {
	'use strict';

	function getDefaultExportFromCjs (x) {
		return x && x.__esModule && Object.prototype.hasOwnProperty.call(x, 'default') ? x['default'] : x;
	}

	var compute$1 = {};

	var complex;
	var hasRequiredComplex;

	function requireComplex () {
		if (hasRequiredComplex) return complex;
		hasRequiredComplex = 1;
		//-------------------------------------------------
		// Add two complex numbers
		//-------------------------------------------------
		var complexAdd = function (a, b)
		{
		    if( !Array.isArray(a) || !Array.isArray(b)  )
		        return [0,0];
		    return [a[0] + b[0], a[1] + b[1]];
		};

		//-------------------------------------------------
		// Subtract two complex numbers
		//-------------------------------------------------
		var complexSubtract = function (a, b)
		{
		    if( !Array.isArray(a) || !Array.isArray(b)  )
		        return [0,0];
		    return [a[0] - b[0], a[1] - b[1]];
		};

		//-------------------------------------------------
		// Multiply two complex numbers
		//
		// (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
		//-------------------------------------------------
		var complexMultiply = function (a, b) 
		{
		    if( !Array.isArray(a) || !Array.isArray(b)  )
		        return [0,0];
		    return [(a[0] * b[0] - a[1] * b[1]), 
		            (a[0] * b[1] + a[1] * b[0])];
		};

		//-------------------------------------------------
		// Calculate |a + bi|
		//
		// sqrt(a*a + b*b)
		//-------------------------------------------------
		var complexMagnitude = function (c) 
		{
		    if( !Array.isArray(c)  )
		        return 0;
		    return Math.sqrt(c[0]*c[0] + c[1]*c[1]); 
		};

		//-------------------------------------------------
		// Exports
		//-------------------------------------------------
		complex = {
		    add: complexAdd,
		    subtract: complexSubtract,
		    multiply: complexMultiply,
		    magnitude: complexMagnitude
		};
		return complex;
	}

	/*===========================================================================*\
	 * Fast Fourier Transform Frequency/Magnitude passes
	 *
	 * (c) Vail Systems. Joshua Jung and Ben Bryan. 2015
	 *
	 * This code is not designed to be highly optimized but as an educational
	 * tool to understand the Fast Fourier Transform.
	\*===========================================================================*/

	var fftutil;
	var hasRequiredFftutil;

	function requireFftutil () {
		if (hasRequiredFftutil) return fftutil;
		hasRequiredFftutil = 1;
		//-------------------------------------------------
		// The following code assumes a complex number is
		// an array: [real, imaginary]
		//-------------------------------------------------
		var complex = requireComplex();


		//-------------------------------------------------
		// By Eulers Formula:
		//
		// e^(i*x) = cos(x) + i*sin(x)
		//
		// and in DFT:
		//
		// x = -2*PI*(k/N)
		//-------------------------------------------------
		var mapExponent = {},
		    exponent = function (k, N) {
		      var x = -2 * Math.PI * (k / N);

		      mapExponent[N] = mapExponent[N] || {};
		      mapExponent[N][k] = mapExponent[N][k] || [Math.cos(x), Math.sin(x)];// [Real, Imaginary]

		      return mapExponent[N][k];
		};

		//-------------------------------------------------
		// Calculate FFT Magnitude for complex numbers.
		//-------------------------------------------------
		var fftMag = function (fftBins) {
		    var ret = fftBins.map(complex.magnitude);
		    return ret.slice(0, ret.length / 2);
		};

		//-------------------------------------------------
		// Calculate Frequency Bins
		// 
		// Returns an array of the frequencies (in hertz) of
		// each FFT bin provided, assuming the sampleRate is
		// samples taken per second.
		//-------------------------------------------------
		var fftFreq = function (fftBins, sampleRate) {
		    var stepFreq = sampleRate / (fftBins.length);
		    var ret = fftBins.slice(0, fftBins.length / 2);

		    return ret.map(function (__, ix) {
		        return ix * stepFreq;
		    });
		};

		//-------------------------------------------------
		// Exports
		//-------------------------------------------------
		fftutil = {
		    fftMag: fftMag,
		    fftFreq: fftFreq,
		    exponent: exponent
		};
		return fftutil;
	}

	var twiddle = {};

	/**
	 * Bit twiddling hacks for JavaScript.
	 *
	 * Author: Mikola Lysenko
	 *
	 * Ported from Stanford bit twiddling hack library:
	 *    http://graphics.stanford.edu/~seander/bithacks.html
	 */

	var hasRequiredTwiddle;

	function requireTwiddle () {
		if (hasRequiredTwiddle) return twiddle;
		hasRequiredTwiddle = 1;

		//Number of bits in an integer
		var INT_BITS = 32;

		//Constants
		twiddle.INT_BITS  = INT_BITS;
		twiddle.INT_MAX   =  0x7fffffff;
		twiddle.INT_MIN   = -1<<(INT_BITS-1);

		//Returns -1, 0, +1 depending on sign of x
		twiddle.sign = function(v) {
		  return (v > 0) - (v < 0);
		};

		//Computes absolute value of integer
		twiddle.abs = function(v) {
		  var mask = v >> (INT_BITS-1);
		  return (v ^ mask) - mask;
		};

		//Computes minimum of integers x and y
		twiddle.min = function(x, y) {
		  return y ^ ((x ^ y) & -(x < y));
		};

		//Computes maximum of integers x and y
		twiddle.max = function(x, y) {
		  return x ^ ((x ^ y) & -(x < y));
		};

		//Checks if a number is a power of two
		twiddle.isPow2 = function(v) {
		  return !(v & (v-1)) && (!!v);
		};

		//Computes log base 2 of v
		twiddle.log2 = function(v) {
		  var r, shift;
		  r =     (v > 0xFFFF) << 4; v >>>= r;
		  shift = (v > 0xFF  ) << 3; v >>>= shift; r |= shift;
		  shift = (v > 0xF   ) << 2; v >>>= shift; r |= shift;
		  shift = (v > 0x3   ) << 1; v >>>= shift; r |= shift;
		  return r | (v >> 1);
		};

		//Computes log base 10 of v
		twiddle.log10 = function(v) {
		  return  (v >= 1000000000) ? 9 : (v >= 100000000) ? 8 : (v >= 10000000) ? 7 :
		          (v >= 1000000) ? 6 : (v >= 100000) ? 5 : (v >= 10000) ? 4 :
		          (v >= 1000) ? 3 : (v >= 100) ? 2 : (v >= 10) ? 1 : 0;
		};

		//Counts number of bits
		twiddle.popCount = function(v) {
		  v = v - ((v >>> 1) & 0x55555555);
		  v = (v & 0x33333333) + ((v >>> 2) & 0x33333333);
		  return ((v + (v >>> 4) & 0xF0F0F0F) * 0x1010101) >>> 24;
		};

		//Counts number of trailing zeros
		function countTrailingZeros(v) {
		  var c = 32;
		  v &= -v;
		  if (v) c--;
		  if (v & 0x0000FFFF) c -= 16;
		  if (v & 0x00FF00FF) c -= 8;
		  if (v & 0x0F0F0F0F) c -= 4;
		  if (v & 0x33333333) c -= 2;
		  if (v & 0x55555555) c -= 1;
		  return c;
		}
		twiddle.countTrailingZeros = countTrailingZeros;

		//Rounds to next power of 2
		twiddle.nextPow2 = function(v) {
		  v += v === 0;
		  --v;
		  v |= v >>> 1;
		  v |= v >>> 2;
		  v |= v >>> 4;
		  v |= v >>> 8;
		  v |= v >>> 16;
		  return v + 1;
		};

		//Rounds down to previous power of 2
		twiddle.prevPow2 = function(v) {
		  v |= v >>> 1;
		  v |= v >>> 2;
		  v |= v >>> 4;
		  v |= v >>> 8;
		  v |= v >>> 16;
		  return v - (v>>>1);
		};

		//Computes parity of word
		twiddle.parity = function(v) {
		  v ^= v >>> 16;
		  v ^= v >>> 8;
		  v ^= v >>> 4;
		  v &= 0xf;
		  return (0x6996 >>> v) & 1;
		};

		var REVERSE_TABLE = new Array(256);

		(function(tab) {
		  for(var i=0; i<256; ++i) {
		    var v = i, r = i, s = 7;
		    for (v >>>= 1; v; v >>>= 1) {
		      r <<= 1;
		      r |= v & 1;
		      --s;
		    }
		    tab[i] = (r << s) & 0xff;
		  }
		})(REVERSE_TABLE);

		//Reverse bits in a 32 bit word
		twiddle.reverse = function(v) {
		  return  (REVERSE_TABLE[ v         & 0xff] << 24) |
		          (REVERSE_TABLE[(v >>> 8)  & 0xff] << 16) |
		          (REVERSE_TABLE[(v >>> 16) & 0xff] << 8)  |
		           REVERSE_TABLE[(v >>> 24) & 0xff];
		};

		//Interleave bits of 2 coordinates with 16 bits.  Useful for fast quadtree codes
		twiddle.interleave2 = function(x, y) {
		  x &= 0xFFFF;
		  x = (x | (x << 8)) & 0x00FF00FF;
		  x = (x | (x << 4)) & 0x0F0F0F0F;
		  x = (x | (x << 2)) & 0x33333333;
		  x = (x | (x << 1)) & 0x55555555;

		  y &= 0xFFFF;
		  y = (y | (y << 8)) & 0x00FF00FF;
		  y = (y | (y << 4)) & 0x0F0F0F0F;
		  y = (y | (y << 2)) & 0x33333333;
		  y = (y | (y << 1)) & 0x55555555;

		  return x | (y << 1);
		};

		//Extracts the nth interleaved component
		twiddle.deinterleave2 = function(v, n) {
		  v = (v >>> n) & 0x55555555;
		  v = (v | (v >>> 1))  & 0x33333333;
		  v = (v | (v >>> 2))  & 0x0F0F0F0F;
		  v = (v | (v >>> 4))  & 0x00FF00FF;
		  v = (v | (v >>> 16)) & 0x000FFFF;
		  return (v << 16) >> 16;
		};


		//Interleave bits of 3 coordinates, each with 10 bits.  Useful for fast octree codes
		twiddle.interleave3 = function(x, y, z) {
		  x &= 0x3FF;
		  x  = (x | (x<<16)) & 4278190335;
		  x  = (x | (x<<8))  & 251719695;
		  x  = (x | (x<<4))  & 3272356035;
		  x  = (x | (x<<2))  & 1227133513;

		  y &= 0x3FF;
		  y  = (y | (y<<16)) & 4278190335;
		  y  = (y | (y<<8))  & 251719695;
		  y  = (y | (y<<4))  & 3272356035;
		  y  = (y | (y<<2))  & 1227133513;
		  x |= (y << 1);
		  
		  z &= 0x3FF;
		  z  = (z | (z<<16)) & 4278190335;
		  z  = (z | (z<<8))  & 251719695;
		  z  = (z | (z<<4))  & 3272356035;
		  z  = (z | (z<<2))  & 1227133513;
		  
		  return x | (z << 2);
		};

		//Extracts nth interleaved component of a 3-tuple
		twiddle.deinterleave3 = function(v, n) {
		  v = (v >>> n)       & 1227133513;
		  v = (v | (v>>>2))   & 3272356035;
		  v = (v | (v>>>4))   & 251719695;
		  v = (v | (v>>>8))   & 4278190335;
		  v = (v | (v>>>16))  & 0x3FF;
		  return (v<<22)>>22;
		};

		//Computes next combination in colexicographic order (this is mistakenly called nextPermutation on the bit twiddling hacks page)
		twiddle.nextCombination = function(v) {
		  var t = v | (v - 1);
		  return (t + 1) | (((~t & -~t) - 1) >>> (countTrailingZeros(v) + 1));
		};
		return twiddle;
	}

	/*===========================================================================*\
	 * Fast Fourier Transform (Cooley-Tukey Method)
	 *
	 * (c) Vail Systems. Joshua Jung and Ben Bryan. 2015
	 *
	 * This code is not designed to be highly optimized but as an educational
	 * tool to understand the Fast Fourier Transform.
	\*===========================================================================*/

	var fft;
	var hasRequiredFft;

	function requireFft () {
		if (hasRequiredFft) return fft;
		hasRequiredFft = 1;
		//------------------------------------------------
		// Note: Some of this code is not optimized and is
		// primarily designed as an educational and testing
		// tool.
		// To get high performace would require transforming
		// the recursive calls into a loop and then loop
		// unrolling. All of this is best accomplished
		// in C or assembly.
		//-------------------------------------------------

		//-------------------------------------------------
		// The following code assumes a complex number is
		// an array: [real, imaginary]
		//-------------------------------------------------
		var complex = requireComplex(),
		    fftUtil = requireFftutil(),
		    twiddle = requireTwiddle();

		fft = {
		  //-------------------------------------------------
		  // Calculate FFT for vector where vector.length
		  // is assumed to be a power of 2.
		  //-------------------------------------------------
		  fft: function fft(vector) {
		    var X = [],
		        N = vector.length;

		    // Base case is X = x + 0i since our input is assumed to be real only.
		    if (N == 1) {
		      if (Array.isArray(vector[0])) //If input vector contains complex numbers
		        return [[vector[0][0], vector[0][1]]];      
		      else
		        return [[vector[0], 0]];
		    }

		    // Recurse: all even samples
		    var X_evens = fft(vector.filter(even)),

		        // Recurse: all odd samples
		        X_odds  = fft(vector.filter(odd));

		    // Now, perform N/2 operations!
		    for (var k = 0; k < N / 2; k++) {
		      // t is a complex number!
		      var t = X_evens[k],
		          e = complex.multiply(fftUtil.exponent(k, N), X_odds[k]);

		      X[k] = complex.add(t, e);
		      X[k + (N / 2)] = complex.subtract(t, e);
		    }

		    function even(__, ix) {
		      return ix % 2 == 0;
		    }

		    function odd(__, ix) {
		      return ix % 2 == 1;
		    }

		    return X;
		  },
		  //-------------------------------------------------
		  // Calculate FFT for vector where vector.length
		  // is assumed to be a power of 2.  This is the in-
		  // place implementation, to avoid the memory
		  // footprint used by recursion.
		  //-------------------------------------------------
		  fftInPlace: function(vector) {
		    var N = vector.length;

		    var trailingZeros = twiddle.countTrailingZeros(N); //Once reversed, this will be leading zeros

		    // Reverse bits
		    for (var k = 0; k < N; k++) {
		      var p = twiddle.reverse(k) >>> (twiddle.INT_BITS - trailingZeros);
		      if (p > k) {
		        var complexTemp = [vector[k], 0];
		        vector[k] = vector[p];
		        vector[p] = complexTemp;
		      } else {
		        vector[p] = [vector[p], 0];
		      }
		    }

		    //Do the DIT now in-place
		    for (var len = 2; len <= N; len += len) {
		      for (var i = 0; i < len / 2; i++) {
		        var w = fftUtil.exponent(i, len);
		        for (var j = 0; j < N / len; j++) {
		          var t = complex.multiply(w, vector[j * len + i + len / 2]);
		          vector[j * len + i + len / 2] = complex.subtract(vector[j * len + i], t);
		          vector[j * len + i] = complex.add(vector[j * len + i], t);
		        }
		      }
		    }
		  }
		};
		return fft;
	}

	/*===========================================================================*\
	 * Inverse Fast Fourier Transform (Cooley-Tukey Method)
	 *
	 * (c) Maximilian Bügler. 2016
	 *
	 * Based on and using the code by
	 * (c) Vail Systems. Joshua Jung and Ben Bryan. 2015
	 *
	 * This code is not designed to be highly optimized but as an educational
	 * tool to understand the Fast Fourier Transform.
	\*===========================================================================*/

	var ifft;
	var hasRequiredIfft;

	function requireIfft () {
		if (hasRequiredIfft) return ifft;
		hasRequiredIfft = 1;
		//------------------------------------------------
		// Note: Some of this code is not optimized and is
		// primarily designed as an educational and testing
		// tool.
		// To get high performace would require transforming
		// the recursive calls into a loop and then loop
		// unrolling. All of this is best accomplished
		// in C or assembly.
		//-------------------------------------------------

		//-------------------------------------------------
		// The following code assumes a complex number is
		// an array: [real, imaginary]
		//-------------------------------------------------

		var fft = requireFft().fft;


		ifft = {
		    ifft: function ifft(signal){
		        //Interchange real and imaginary parts
		        var csignal=[];
		        for(var i=0; i<signal.length; i++){
		            csignal[i]=[signal[i][1], signal[i][0]];
		        }
		    
		        //Apply fft
		        var ps=fft(csignal);
		        
		        //Interchange real and imaginary parts and normalize
		        var res=[];
		        for(var j=0; j<ps.length; j++){
		            res[j]=[ps[j][1]/ps.length, ps[j][0]/ps.length];
		        }
		        return res;
		    }
		};
		return ifft;
	}

	/*===========================================================================*\
	 * Discrete Fourier Transform (O(n^2) brute-force method)
	 *
	 * (c) Vail Systems. Joshua Jung and Ben Bryan. 2015
	 *
	 * This code is not designed to be highly optimized but as an educational
	 * tool to understand the Fast Fourier Transform.
	\*===========================================================================*/

	var dft_1;
	var hasRequiredDft;

	function requireDft () {
		if (hasRequiredDft) return dft_1;
		hasRequiredDft = 1;
		//------------------------------------------------
		// Note: this code is not optimized and is
		// primarily designed as an educational and testing
		// tool.
		//------------------------------------------------
		var complex = requireComplex();
		var fftUtil = requireFftutil();

		//-------------------------------------------------
		// Calculate brute-force O(n^2) DFT for vector.
		//-------------------------------------------------
		var dft = function(vector) {
		  var X = [],
		      N = vector.length;

		  for (var k = 0; k < N; k++) {
		    X[k] = [0, 0]; //Initialize to a 0-valued complex number.

		    for (var i = 0; i < N; i++) {
		      var exp = fftUtil.exponent(k * i, N);
		      var term;
		      if (Array.isArray(vector[i]))
		        term = complex.multiply(vector[i], exp);//If input vector contains complex numbers
		      else
		        term = complex.multiply([vector[i], 0], exp);//Complex mult of the signal with the exponential term.  
		      X[k] = complex.add(X[k], term); //Complex summation of X[k] and exponential
		    }
		  }

		  return X;
		};

		dft_1 = dft;
		return dft_1;
	}

	/*===========================================================================*\
	 * Inverse Discrete Fourier Transform (O(n^2) brute-force method)
	 *
	 * (c) Maximilian Bügler. 2016
	 *
	 * Based on and using the code by
	 * (c) Vail Systems. Joshua Jung and Ben Bryan. 2015
	 *
	 * This code is not designed to be highly optimized but as an educational
	 * tool to understand the Fast Fourier Transform.
	\*===========================================================================*/

	var idft_1;
	var hasRequiredIdft;

	function requireIdft () {
		if (hasRequiredIdft) return idft_1;
		hasRequiredIdft = 1;
		//------------------------------------------------
		// Note: Some of this code is not optimized and is
		// primarily designed as an educational and testing
		// tool.
		//-------------------------------------------------

		//-------------------------------------------------
		// The following code assumes a complex number is
		// an array: [real, imaginary]
		//-------------------------------------------------
		var dft = requireDft();

		function idft(signal) {
		    //Interchange real and imaginary parts
		    var csignal = [];
		    for (var i = 0; i < signal.length; i++) {
		        csignal[i] = [signal[i][1], signal[i][0]];
		    }

		    //Apply dft
		    var ps = dft(csignal);

		    //Interchange real and imaginary parts and normalize
		    var res = [];
		    for (var j = 0; j < ps.length; j++) {
		        res[j] = [ps[j][1] / ps.length, ps[j][0] / ps.length];
		    }
		    return res;
		}

		idft_1 = idft;
		return idft_1;
	}

	/*===========================================================================*\
	 * Fast Fourier Transform (Cooley-Tukey Method)
	 *
	 * (c) Vail Systems. Joshua Jung and Ben Bryan. 2015
	 *
	 * This code is not designed to be highly optimized but as an educational
	 * tool to understand the Fast Fourier Transform.
	\*===========================================================================*/

	var fftJs;
	var hasRequiredFftJs;

	function requireFftJs () {
		if (hasRequiredFftJs) return fftJs;
		hasRequiredFftJs = 1;
		fftJs = {
		    fft: requireFft().fft,
		    ifft: requireIfft().ifft,
		    fftInPlace: requireFft().fftInPlace,
		    util: requireFftutil(),
		    dft: requireDft(),
		    idft: requireIdft()
		};
		return fftJs;
	}

	var paddedSignal = {};

	var hasRequiredPaddedSignal;

	function requirePaddedSignal () {
		if (hasRequiredPaddedSignal) return paddedSignal;
		hasRequiredPaddedSignal = 1;
		Object.defineProperty(paddedSignal, "__esModule", { value: true });
		paddedSignal.paddedSignal = paddedSignal$1;
		/**
		 *
		 * Pads the input signal array with zeros to the next power of two length.
		 * @param signals {number[]} The input signal array.
		 * @returns {number[]} The padded signal array.
		 */
		function paddedSignal$1(signals) {
		    const targetLength = Math.pow(2, Math.ceil(Math.log2(signals.length)));
		    return signals.concat(new Array(targetLength - signals.length).fill(0));
		}
		return paddedSignal;
	}

	var filterSignal = {};

	var hasRequiredFilterSignal;

	function requireFilterSignal () {
		if (hasRequiredFilterSignal) return filterSignal;
		hasRequiredFilterSignal = 1;
		Object.defineProperty(filterSignal, "__esModule", { value: true });
		filterSignal.filterSignal = filterSignal$1;
		function filterSignal$1(signal) {
		    return signal.map(v => isNaN(v) ? 0 : v);
		}
		return filterSignal;
	}

	var shindoFilter = {};

	var hasRequiredShindoFilter;

	function requireShindoFilter () {
		if (hasRequiredShindoFilter) return shindoFilter;
		hasRequiredShindoFilter = 1;
		Object.defineProperty(shindoFilter, "__esModule", { value: true });
		shindoFilter.TotalFilterImpl = TotalFilterImpl;
		/**
		 *
		 * Applies a total filter to the given signal in the frequency domain.
		 * @param signal {Phasors} The input signal in frequency domain.
		 * @param sample {number} The sample rate of the original signal.
		 * @returns {Phasors} The filtered signal in frequency domain.
		 */
		function TotalFilterImpl(signal, sample) {
		    const dt = 1 / sample;
		    const Nf = signal.length;
		    const Nfhalf = Nf / 2;
		    for (let i = 0; i <= Nfhalf; i++) {
		        const freq = i / (dt * Nf);
		        const x = freq / 10;
		        const fc = Math.sqrt(1 / freq);
		        const fh = 1 / Math.sqrt(1 +
		            0.694 * Math.pow(x, 2) +
		            0.241 * Math.pow(x, 4) +
		            0.0557 * Math.pow(x, 6) +
		            0.009664 * Math.pow(x, 8) +
		            0.00134 * Math.pow(x, 10) +
		            0.000155 * Math.pow(x, 12));
		        const fl = Math.pow((1 - Math.exp(Math.pow((-(freq / 0.5)), 3))), 0.5);
		        const fa = fc * fh * fl;
		        if (signal[i] === undefined)
		            continue;
		        // @ts-ignore
		        const fs = [signal[i][0] * fa, signal[i][1] * fa];
		        if (isNaN(fs[0]) || isNaN(fs[1])) {
		            continue;
		        }
		        signal[i] = fs;
		    }
		    return signal;
		}
		return shindoFilter;
	}

	var vectorTransform = {};

	var hasRequiredVectorTransform;

	function requireVectorTransform () {
		if (hasRequiredVectorTransform) return vectorTransform;
		hasRequiredVectorTransform = 1;
		Object.defineProperty(vectorTransform, "__esModule", { value: true });
		vectorTransform.vectorTransform = vectorTransform$1;
		/**
		 *
		 * Synthesizes three directional components (NS, EW, UD) into a single vector magnitude array.
		 * @param NS {number[]} The North-South component.
		 * @param EW {number[]} The East-West component.
		 * @param UD {number[]} The Up-Down component.
		 * @returns {number[]} The synthesized vector magnitude array.
		 */
		function vectorTransform$1(NS, EW, UD) {
		    const length = Math.min(NS.length, EW.length, UD.length);
		    const synthesized = new Array(length).fill(0);
		    for (let i = 0; i < length; i++) {
		        synthesized[i] = Math.sqrt(Math.pow((NS[i] || 1), 2) +
		            Math.pow((EW[i] || 1), 2) +
		            Math.pow((UD[i] || 1), 2));
		    }
		    return synthesized;
		}
		return vectorTransform;
	}

	var hasRequiredCompute;

	function requireCompute () {
		if (hasRequiredCompute) return compute$1;
		hasRequiredCompute = 1;
		Object.defineProperty(compute$1, "__esModule", { value: true });
		compute$1.default = computeShindo;
		//@ts-ignore
		const fft_js_1 = requireFftJs();
		const paddedSignal_1 = requirePaddedSignal();
		const filterSignal_1 = requireFilterSignal();
		const shindoFilter_1 = requireShindoFilter();
		const vectorTransform_1 = requireVectorTransform();
		const internalSignalStore = [];
		function computeShindo(xAcc, yAcc, zAcc, SamplingRate = 100) {
		    internalSignalStore.push([xAcc, yAcc, zAcc]);
		    if (internalSignalStore.length < 100)
		        return { gal: -1, intensity: -1 };
		    else
		        internalSignalStore.shift();
		    const xSeries = (0, filterSignal_1.filterSignal)((0, paddedSignal_1.paddedSignal)(internalSignalStore.map(([x, _y, _z]) => x)));
		    const ySeries = (0, filterSignal_1.filterSignal)((0, paddedSignal_1.paddedSignal)(internalSignalStore.map(([_x, y, _z]) => y)));
		    const zSeries = (0, filterSignal_1.filterSignal)((0, paddedSignal_1.paddedSignal)(internalSignalStore.map(([_x, _y, z]) => z)));
		    const xPhasors = (0, fft_js_1.fft)(xSeries);
		    const yPhasors = (0, fft_js_1.fft)(ySeries);
		    const zPhasors = (0, fft_js_1.fft)(zSeries);
		    const xFiltered = (0, shindoFilter_1.TotalFilterImpl)(xPhasors, SamplingRate);
		    const yFiltered = (0, shindoFilter_1.TotalFilterImpl)(yPhasors, SamplingRate);
		    const zFiltered = (0, shindoFilter_1.TotalFilterImpl)(zPhasors, SamplingRate);
		    const xInversed = ((0, fft_js_1.ifft)(xFiltered).map(([x, _y]) => x));
		    const yInversed = ((0, fft_js_1.ifft)(yFiltered).map(([x, _y]) => x));
		    const zInversed = ((0, fft_js_1.ifft)(zFiltered).map(([x, _y]) => x));
		    const totalVector = (0, vectorTransform_1.vectorTransform)(xInversed, yInversed, zInversed);
		    totalVector.sort((a, b) => b - a);
		    const a = Math.floor(0.3 / (1 / SamplingRate)) - 1;
		    const gal = totalVector[a];
		    const Intensity = 2 * Math.log10(gal) + 0.94;
		    const kIntensity = (Math.floor(Math.round(Intensity * 100) / 10) / 10).toFixed(1);
		    return {
		        gal,
		        intensity: parseFloat(kIntensity)
		    };
		}
		return compute$1;
	}

	var computeExports = requireCompute();
	var compute = /*@__PURE__*/getDefaultExportFromCjs(computeExports);

	return compute;

})();
