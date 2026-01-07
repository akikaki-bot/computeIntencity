import * as fs from "fs"

export class LoadCSV {
    constructor( private path : string ){}

    public load(): string[][]{
        const data = fs.readFileSync(this.path, 'utf-8');

        if( data.length === 0 ){
            throw new Error("CSV file is empty");
        }

        const lines = data.split('\n').map( line => line.trim() ).filter( line => line.length > 0 );
        const result: string[][] = lines.map( line => line.split(',').map( item => item.trim() ) );

        return result;
    }
}