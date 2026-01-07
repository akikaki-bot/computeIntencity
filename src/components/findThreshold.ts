

export function findThreshold(
    signal: number[],
    dt: number,
    targetTime: number = 0.3,
    tol: number = 1e-3
): number {
    let low = 0
    let high = Math.max(...signal)

    while (high - low > tol) {
        const mid = (low + high) / 2;
        const timeAboveThreshold = ( 
            signal.reduce((acc, val) => acc + (val >= mid ? dt : 0), 0) 
        ) * dt;
        if( timeAboveThreshold > targetTime ){
            low = mid;
        } else {
            high = mid;
        }
    }
    return (low + high) / 2;
}