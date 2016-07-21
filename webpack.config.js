module.exports = {
    entry: "./src/entry.ts",
    output: {
        path: "dist",
        filename: "beachballs.js",
        libraryTarget: 'umd',
        umdNamedDefine: true
    },
    module: {
        loaders: [
            {test: /\.ts?$/, exclude: /(node_modules)/, loader: 'ts'}
        ]
    }
    
};