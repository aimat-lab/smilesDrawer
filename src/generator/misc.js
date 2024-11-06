const crypto = require('crypto')
const path = require('path')
const fs = require('fs')
const readline = require('readline')
const yargs = require('yargs/yargs')
const { hideBin } = require('yargs/helpers')

const readSmilesFromCsv = async(file, smilesCol, n = 100) => {
  const stream = fs.createReadStream(file)
  const rl = readline.createInterface({
    input: stream,
    crlfDelay: Infinity
  })

  const result = []
  for await (const line of rl) {
    const smiles = line.split(',')[smilesCol]
    result.push(smiles)
    if (result.length >= 1.25 * n) {
      break
    }
  }

  stream.destroy()
  return result
}

const cliParams = () => {
  const {
    outputDirectory,
    amount, batchSize, size, fonts, fontWeights,
    concurrency,
    outputSvg, outputLabels, outputFlat,
    clean,
    minSmilesLength, maxSmilesLength,
    fromCsvFile: csvFile,
    fromCsvColumn: csvColumn
  } = yargs(hideBin(process.argv)).argv

  const config = {
    csvFile: path.resolve(csvFile),
    csvColumn: csvColumn,
    outputDirectory: path.resolve(outputDirectory),
    size: Number(size) || null,
    fonts: fonts ? fonts.split(',') : ['Roboto'],
    fontWeights: fontWeights ? fontWeights.split(',').map(x => Number(x)) : [200],
    concurrency: Number(concurrency) || 4,
    minSmilesLength: Number(minSmilesLength) || 0,
    maxSmilesLength: Number(maxSmilesLength) || 1000,
    outputSvg: !!outputSvg,
    outputLabels: !!outputLabels,
    outputFlat: !!outputFlat,
    amount: Number(amount) || null,
    batchSize: Number(batchSize) || 100,
    clean: !!clean
  }

  const invalid = Object.entries(config).filter(([key, value]) => value === null)

  if (invalid.length) {
    throw new Error(`invalid configuration values: ${JSON.stringify(invalid, null)}`)
  }

  return config
}

const hash = function(x) {
  return crypto.createHash('sha256').update(x).digest('hex')
}

const wait = ms => new Promise((resolve, reject) => {
  setTimeout(resolve, ms)
})

module.exports = {
  readSmilesFromCsv,
  cliParams,
  hash,
  wait
}
