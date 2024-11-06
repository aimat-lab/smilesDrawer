(async() => {
  const { v4: uuid } = require('uuid')
  const path = require('path')
  const treekill = require('tree-kill')
  const fs = require('fs-extra')
  const _ = require('lodash')
  const { fork } = require('child_process')

  const oldCwd = process.cwd()
  const newCwd = path.resolve(path.dirname(__filename))

  if (oldCwd !== newCwd) {
    console.log(`changing working directory from ${oldCwd} to ${newCwd}`)
    process.chdir(newCwd)
  }

  const { readSmilesFromCsv, cliParams, wait } = require('./src/generator/misc')

  const conf = cliParams()

  if (conf.clean) {
    console.log(`deleting ${conf.outputDirectory}`)
    await fs.emptyDir(conf.outputDirectory)
  }

  await fs.ensureDir(conf.outputDirectory)

  console.log('reading smiles file')
  const smilesList = (await readSmilesFromCsv(conf.csvFile, conf.csvColumn, conf.amount))
    .filter(s => s.length >= conf.minSmilesLength && s.length <= conf.maxSmilesLength)
    .slice(0, conf.amount)

  console.log(`found ${smilesList.length} SMILES strings with length between ${conf.minSmilesLength} and ${conf.maxSmilesLength} characters`)

  const label = `generating ${smilesList.length} images with concurrency ${conf.concurrency}`
  console.time(label)

  // aneb: clear state after every n images
  const numberOfBatches = Math.round(conf.amount / conf.batchSize)
  const batches = _.chunk(smilesList, Math.round(conf.amount / numberOfBatches))

  console.log(`processing ${conf.amount} images in ${batches.length} batches (batch size ${conf.batchSize}, concurrency ${conf.concurrency})`)

  const debug = typeof v8debug === 'object'
  const children = {}
  let done = 0
  for (const [index, smilesList] of batches.entries()) {
    const tmpDir = path.join('browser', uuid())
    await fs.ensureDir(tmpDir)

    const browserOptions = {
      userDataDir: tmpDir,
      temporaryDirectory: tmpDir,
      headless: true,
      devtools: false,
      protocolTimeout: 100_000_000
    }

    const message = { conf, smilesList, browserOptions }

    const args = { }

    // aneb: inspector error are IDE-related and do not occur when calling node from command line
    if (debug) {
      const offset = (index % conf.concurrency) + 1
      const port = process.debugPort + offset

      console.log(`adding debug port ${port}`)
      args.execArgv = [`--inspect=${port}`]
    }

    const child = fork('src/worker.js', args)
    children[child.pid] = { }

    child.on('message', function({ browserPid }) {
      children[this.pid] = { browserPid }
    })

    child.on('exit', function(code) {
      done += 1
      const state = code ? 'FAIL' : 'SUCCESS'

      console.log(`${new Date().toUTCString()} - ${state} ${done}/${batches.length} done`)

      treekill(children[this.pid].browserPid, 'SIGKILL')
      treekill(this.pid, 'SIGKILL')
      delete children[this.pid]
    })

    child.send(message)

    while (Object.keys(children).length >= conf.concurrency) {
      await wait(1000)
    }

    await fs.remove(tmpDir)
  }

  // aneb: must also wait for last processes to finish
  while (Object.keys(children).length !== 0) {
    await wait(100)
  }

  console.timeEnd(label)
})()
