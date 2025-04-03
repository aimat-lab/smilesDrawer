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

const getElementGraph = (graph, labels) => {
  if (!labels.length) {
    return []
  }
  const bondGraph = []
  const nonTerminals = graph.vertices.filter(v => v.neighbours.length > 1)

  const triplets = new Set()
  const elementIds = new Set(labels.map(l => l.id))

  for (const vertex of nonTerminals) {
    for (const s of vertex.neighbours) {
      for (const t of vertex.neighbours) {
        if (s === t) {
          continue
        }
        triplets.add(`${s}_${vertex.id}_${t}`)
      }
    }
  }

  for (const triplet of triplets) {
    const [a, b, c] = triplet.split('_').map(x => Number(x))
    const edgeA = graph.vertexIdsToEdgeId[`${a}_${b}`]
    const edgeB = graph.vertexIdsToEdgeId[`${b}_${c}`]

    bondGraph.push([edgeA, edgeB])
  }

  const elementGraph = []

  for (const [a, b] of bondGraph) {
    const edgeA = graph.edges[a]
    const edgeB = graph.edges[b]

    const sourceA = graph.vertices[edgeA.sourceId]
    const targetA = graph.vertices[edgeA.targetId]
    const sourceB = graph.vertices[edgeB.sourceId]
    const targetB = graph.vertices[edgeB.targetId]

    const counts = [sourceA.id, targetA.id, sourceB.id, targetB.id].reduce((acc, x) => {
      acc[x] = (acc[x] || 0) + 1
      return acc
    }, {})

    let intersectionId = null
    let intersectionCount = -1
    for (const [id, count] in Object.entries(counts)) {
      if (count > intersectionCount) {
        intersectionId = id
        intersectionCount = count
      }
    }

    const sourceADrawn = elementIds.has(`vertex-id-${sourceA.id}`)
    const targetADrawn = elementIds.has(`vertex-id-${targetA.id}`)
    const sourceBDrawn = elementIds.has(`vertex-id-${sourceB.id}`)
    const targetBDrawn = elementIds.has(`vertex-id-${targetB.id}`)
    const intersectionDrawn = elementIds.has(`vertex-id-${intersectionId}`)

    if (!intersectionDrawn) {
      elementGraph.push([`edge-id-${edgeA.id}`, `edge-id-${edgeB.id}`])
    }

    if (sourceADrawn) {
      elementGraph.push([`vertex-id-${sourceA.id}`, `edge-id-${edgeA.id}`])
    }
    if (targetADrawn) {
      elementGraph.push([`vertex-id-${targetA.id}`, `edge-id-${edgeA.id}`])
    }
    if (sourceBDrawn) {
      elementGraph.push([`vertex-id-${sourceB.id}`, `edge-id-${edgeB.id}`])
    }
    if (targetBDrawn) {
      elementGraph.push([`vertex-id-${targetB.id}`, `edge-id-${edgeB.id}`])
    }
  }

  const sorted = new Set(elementGraph.map(pair => JSON.stringify(pair.slice().sort())))
  const pairs = Array.from(sorted).map(str => JSON.parse(str))

  const filtered = []
  for (const [s, t] of pairs) {
    if (!elementIds.has(s)) {
      continue
    }
    if (!elementIds.has(t)) {
      continue
    }

    filtered.push([s, t])
  }

  return filtered
}

module.exports = {
  readSmilesFromCsv,
  cliParams,
  hash,
  wait,
  getElementGraph
}
