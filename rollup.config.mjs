import typescript from "@rollup/plugin-typescript";
import nodeResolve from "@rollup/plugin-node-resolve";
import cjs from "@rollup/plugin-commonjs";
export default {
    input: './dist/compute.js',
    output: {
        file: './bundle/main.js',
        format: 'iife',
        name: "computeShindo"
    },
    plugins: [
        nodeResolve({
            moduleDirectories: ['node_modules'],
            extensions: ['.js', '.ts']
        }),
        cjs(),
    ]
}
