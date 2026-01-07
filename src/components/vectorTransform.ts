
/**
 * 
 * Synthesizes three directional components (NS, EW, UD) into a single vector magnitude array.
 * @param NS {number[]} The North-South component.
 * @param EW {number[]} The East-West component.
 * @param UD {number[]} The Up-Down component.
 * @returns {number[]} The synthesized vector magnitude array.
 */
export function vectorTransform(
    NS: number[],
    EW: number[],
    UD: number[]
) : number[] {
    const length = Math.min( NS.length, EW.length, UD.length );
    const synthesized: number[] = new Array( length ).fill( 0 );

    for( let i = 0; i < length; i++ ){
        synthesized[i] = Math.sqrt(
            ( NS[i] || 1 ) ** 2 +
            ( EW[i] || 1 ) ** 2 +
            ( UD[i] || 1 ) ** 2
        );
    }

    return synthesized;
}