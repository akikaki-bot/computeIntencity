//@ts-ignore
import { fft, ifft } from "fft-js"
import { paddedSignal } from "./components/paddedSignal"
import { filterSignal } from "./components/filterSignal"
import { TotalFilterImpl } from "./components/shindoFilter"
import { vectorTransform } from "./components/vectorTransform"

const internalSignalStore: [Float, Float, Float][] = [];

interface ComputeResult {
    gal: number,
    intensity: number
}

type Float = number
export default function computeShindo(
    xAcc: Float,
    yAcc: Float,
    zAcc: Float,
    SamplingRate: number = 100
): ComputeResult {
    internalSignalStore.push([xAcc, yAcc, zAcc]);

    if (internalSignalStore.length < 100) return { gal: -1, intensity: -1 };
    else internalSignalStore.shift();

    const xSeries = filterSignal(paddedSignal(internalSignalStore.map(([x, _y, _z]) => x)))
    const ySeries = filterSignal(paddedSignal(internalSignalStore.map(([_x, y, _z]) => y)))
    const zSeries = filterSignal(paddedSignal(internalSignalStore.map(([_x, _y, z]) => z)))

    const xPhasors = fft(xSeries);
    const yPhasors = fft(ySeries);
    const zPhasors = fft(zSeries);

    const xFiltered = TotalFilterImpl(xPhasors, SamplingRate);
    const yFiltered = TotalFilterImpl(yPhasors, SamplingRate);
    const zFiltered = TotalFilterImpl(zPhasors, SamplingRate);

    const xInversed = (ifft(xFiltered).map(([x, _y]) => x));
    const yInversed = (ifft(yFiltered).map(([x, _y]) => x));
    const zInversed = (ifft(zFiltered).map(([x, _y]) => x));

    const totalVector = vectorTransform(xInversed, yInversed, zInversed);
    totalVector.sort((a, b) => b - a);

    const a = Math.floor(0.3 / (1 / SamplingRate)) - 1;
    const gal = totalVector[a] as number;

    const Intensity = 2 * Math.log10(gal) + 0.94;
    const kIntensity = (Math.floor(Math.round(Intensity * 100) / 10) / 10).toFixed(1);

    return {
        gal,
        intensity: parseFloat(kIntensity)
    }
}
