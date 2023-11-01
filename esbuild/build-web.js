import esbuild from 'esbuild'
import {commonConfig} from './common.js'


const config = {
  ...commonConfig,
  outfile: './dist/bundle-web.js',
  platform: 'browser',
}

console.log('CONFIG:', config)

esbuild
  .build(config)
  .then((result) => {
    console.log('Build for web succeeded.')
  })
  .catch((err) => {
    console.error(err)
    process.exit(1)
  })
