//@ts-ignore
import { fft, ifft } from "fft-js";
import { LoadCSV } from "./components/loadCSV";
import { log, duration } from "./components/console";
import { paddedSignal } from "./components/paddedSignal";
import { filterSignal } from "./components/filterSignal";
import { TotalFilterImpl } from "./components/shindoFilter";
import { vectorTransform } from "./components/vectorTransform";

let SAMPLING_RATE = 100;

const argv = process.argv.slice(2);

if (argv.length < 1) {
    console.error("❌ Please provide the path to the Kyoshin data CSV file as a command-line argument.");
    process.exit(1);
}

const dataPath = argv[0];

if (typeof dataPath !== "string" || dataPath.trim() === "") {
    console.error("❌ Invalid file path provided. Please check the path and try again.");
    process.exit(1);
}

console.log("📁 Kyoshin Data Path:", dataPath);
kyoshinCalculation(dataPath);

export function kyoshinCalculation(path: string) {

    const [startInit, startInitTimeout] = log("🚀 Initializing computation...");
    const [startLoading, startLoadingTimeout] = log("📦 Loading Kyoshin data...");

    const kyoshinData = new LoadCSV(path).load();

    const detectSamplingRate = kyoshinData.filter(row => typeof row.find(cell => /SAMPLING/g.test(cell)) !== "undefined")[0] ?? [];
    console.log("🔍 Detected Sampling Rate Info:", detectSamplingRate[0] ?? "unknown");
    if (detectSamplingRate.length > 1) {
        const stringMatch = detectSamplingRate[0] ?? ""
        const HzMatch = stringMatch.match(/(\d+)\s*Hz/);
        if (HzMatch) {
            const hz = parseInt(HzMatch[1] ?? "0");
            if (!isNaN(hz) && hz > 0) {
                SAMPLING_RATE = hz;
                console.log("⚙️ Using detected sampling rate:", SAMPLING_RATE, "Hz");
            }
        }
    }

    clearInterval(startLoadingTimeout);
    console.log("✨ Loaded Kyoshin data:", kyoshinData.length, "rows. \n✅ Done in", duration(startLoading), "ms");
    const [__initing, __initingTimeout] = log("⌛ Initiating FFT processing...");

    const NSSeries: number[] = filterSignal(paddedSignal(kyoshinData.map(row => parseFloat(row[0] ?? "0"))));
    const EWSeries: number[] = filterSignal(paddedSignal(kyoshinData.map(row => parseFloat(row[1] ?? "0"))));
    const UDSeries: number[] = filterSignal(paddedSignal(kyoshinData.map(row => parseFloat(row[2] ?? "0"))));

    clearInterval(__initingTimeout);
    console.log(`✨ FFT processing initiated.\n✅ Done in ${duration(__initing)} ms, with series length: ${NSSeries.length}`);
    const [__fftStart, __fftTimeout] = log("⌛ Performing FFT...");

    const nsPhasors = fft(NSSeries);
    const ewPhasors = fft(EWSeries);
    const udPhasors = fft(UDSeries);

    clearInterval(__fftTimeout);
    console.log("✨ FFT processing completed.\ndone in", duration(__fftStart), "ms");

    const [__filtering, __filteringTimeout] = log("⌛ Filtering significant frequencies...");

    const nsFiltered = TotalFilterImpl(nsPhasors, SAMPLING_RATE);
    const ewFiltered = TotalFilterImpl(ewPhasors, SAMPLING_RATE);
    const udFiltered = TotalFilterImpl(udPhasors, SAMPLING_RATE);


    clearInterval(__filteringTimeout);
    console.log("✨ Filtering completed.\n✅ Done in", duration(__filtering), "ms");
    const [__inverseFFT, __inverseFFTTimeout] = log("⌛ Performing Inverse FFT...");

    const nsInverse = (ifft(nsFiltered)).map((v: [number, number]) => v[0]);
    const ewInverse = (ifft(ewFiltered)).map((v: [number, number]) => v[0]);
    const udInverse = (ifft(udFiltered)).map((v: [number, number]) => v[0]);

    clearInterval(__inverseFFTTimeout);
    console.log("✨ Inverse FFT completed.\n└✅ Done in", duration(__inverseFFT), "ms");
    const [__vectorTransform, __vectorTransformTimeout] = log("⌛ Performing Vector Transformation...");

    const syntheticGal = vectorTransform(
        nsInverse,
        ewInverse,
        udInverse
    );

    syntheticGal.sort((a, b) => b - a);

    clearInterval(__vectorTransformTimeout);
    console.log("✨ Vector Transformation completed.\n✅ Done in", duration(__vectorTransform), "ms");
    const [__calcuateA, __calcuateATimeout] = log("⌛ Calculating A...");

    const a = Math.floor(0.3 / (1 / SAMPLING_RATE)) - 1;

    const gal = syntheticGal[a] as number;

    clearInterval(__calcuateATimeout);
    console.log(`✨ Calculation of A completed.\n✅ Done in ${duration(__calcuateA)} ms.\n\n🏁 Final Result: ${gal}gal`);

    const Intensity = 2 * Math.log10(gal) + 0.94;

    console.log(`🏁 Estimated Seismic Intensity: ${Math.round(Intensity * 10) / 10}`);

    const kInt = (Math.floor(Math.round(Intensity * 100) / 10) / 10).toFixed(1);

    console.log(`🏁 Estimated Seismic JMA Intensity: ${kInt}`);

    clearInterval(startInitTimeout);
    console.log("🚀 Initialization and computation completed.\n✅ Done in", duration(startInit), "ms");

}
