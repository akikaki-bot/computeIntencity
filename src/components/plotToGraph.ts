import path from "path"
const SAVEPATH = path.resolve(__dirname, "../../graphs")
import { ChartJSNodeCanvas } from "chartjs-node-canvas"
import { writeFileSync } from "fs";

export async function plotToGraph(
    title: string,
    data: number[],
    freqTrans: boolean = false
) : Promise<void> {
    const canvas = new ChartJSNodeCanvas({ width: 800, height: 600 });
    const buf = await canvas.renderToBuffer({
        type: "line",
        data: {
            labels: data.map( (_, idx) => idx ),
            datasets: [
                {
                    label: title,
                    data: data,
                    borderColor: 'rgba(75, 192, 192, 1)',
                    backgroundColor: 'rgba(75, 192, 192, 0.2)',
                    fill: false,
                    pointRadius: 0,
                    borderWidth: 1,
                }
            ]
        },
        options: {
            scales: {
                x: {
                    title: {
                        display: true,
                        text: 'Frequency (Hz)'
                    }
                },
                y: {
                    title: {
                        display: true,
                        text: 'Amplitude'
                    }
                }
            }
        }
    })
    writeFileSync( path.resolve(SAVEPATH, `${title}.png`), buf );
}