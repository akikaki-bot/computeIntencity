
/**
 * 
 * Pads the input signal array with zeros to the next power of two length.
 * @param signals {number[]} The input signal array.
 * @returns {number[]} The padded signal array.
 */
export function paddedSignal( signals: number[] ): number[] {
    const targetLength = Math.pow( 2, Math.ceil( Math.log2( signals.length ) ) );
    return signals.concat(new Array(targetLength - signals.length).fill(0));
}