export function filterSignal( signal: number[] ): number[] {
    return signal.map( v => isNaN(v) ? 0 : v );
}