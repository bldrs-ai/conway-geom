import esbuild from 'esbuild'
import {commonConfig} from './common.js'


const config = {
  ...commonConfig,
  outfile: './dist/bundle-node.js',
  platform: 'node',
}

console.log('CONFIG:', config)


esbuild
  .build(config)
  .then((result) => {
    console.log('Build for node succeeded.')
  })
  .catch((err) => {
    console.error(err)
    process.exit(1)
  })
