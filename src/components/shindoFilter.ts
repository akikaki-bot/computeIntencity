//@ts-ignore
import { Phasors } from "fft-js";

/**
 * 
 * Applies a total filter to the given signal in the frequency domain.
 * @param signal {Phasors} The input signal in frequency domain.
 * @param sample {number} The sample rate of the original signal.
 * @returns {Phasors} The filtered signal in frequency domain.
 */
export function TotalFilterImpl( signal: Phasors, sample: number ): Phasors {
    const dt = 1 / sample;

    const Nf = signal.length;
    const Nfhalf = Nf / 2;

    for( let i = 0; i <= Nfhalf; i++ ){
        const freq = i / ( dt * Nf );
        const x = freq / 10;

        const fc = Math.sqrt( 1 / freq );
        const fh = 1 / Math.sqrt(
            1 +
            0.694 * x ** 2 +
            0.241 * x ** 4 +
            0.0557 * x ** 6 +
            0.009664 * x ** 8 +
            0.00134 * x ** 10 +
            0.000155 * x ** 12
        );
        const fl = (
            1 - Math.exp((-( freq / 0.5 ))  ** 3 )
        ) ** 0.5;
        const fa = fc * fh * fl;

        if( signal[i] === undefined ) continue;
        // @ts-ignore
        const fs: [ number, number ] = [ signal[i][0] * fa, signal[i][1] * fa ];
        if( isNaN( fs[0] ) || isNaN( fs[1] ) ){
            continue;
        }
        signal[i] = fs;
    }

    return signal;
}