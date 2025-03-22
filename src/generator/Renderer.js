const crypto = require('crypto')
const fs = require('fs-extra')

const _ = require('lodash')
const { JSDOM } = require('jsdom')
const { xml2js, js2xml } = require('xml-js')

const Parser = require('../drawer/Parser')
const SvgDrawer = require('../drawer/SvgDrawer')
const SVG = require('./SVG')
const { bondLabels } = require('./types')
const { getPositionInfoFromSvg, resizeImage, drawMasksAroundTextElements } = require('./browser')
const { getElementGraph } = require('./misc')

function Renderer({ outputDirectory, size, fonts, fontWeights, concurrency, outputSvg, outputLabels, outputFlat }) {
  // aneb: find out why this does not work in above scope ...
  const colorMap = require('./colors')

  this.parser = Parser
  this.outputDirectory = outputDirectory
  this.size = size
  this.fonts = fonts
  this.fontWeights = fontWeights
  this.colorMap = colorMap
  this.concurrency = concurrency
  this.outputSvg = outputSvg
  this.outputLabels = outputLabels
  this.outputFlat = outputFlat
  this.waitOptions = { waitUntil: 'domcontentloaded', timeout: 10000 }

  this.svgHelper = new SVG()

  const { document, XMLSerializer } = (new JSDOM('')).window
  this.document = document
  this.XMLSerializer = new XMLSerializer()
}

Renderer.prototype.id = function(x) {
  return crypto.createHash('sha256').update(x).digest('hex')
}

Renderer.prototype.color = function(color, circle = false) {
  const fill = circle ? color : 'none'
  return `fill: ${fill}; stroke: ${color};`
}

Renderer.prototype.randomColorMap = function(keys) {
  const map = {}
  for (const key of keys) {
    // https://stackoverflow.com/questions/5092808/how-do-i-randomly-generate-html-hex-color-codes-using-javascript
    map[key] = '#000000'.replace(/0/g, function() { return (~~(Math.random() * 16)).toString(16) })
  }

  return map
}

Renderer.prototype.makeEdgeAttributesNumeric = function(edge) {
  // aneb: one can only read html attributes as strings, postprocessing is done in one place to avoid handling
  // all types of bonds in browser code which cannot be debugged

  // wedge solid bond is drawn as polygon, all others are drawn from single lines which need to be merged
  if (edge.label === bondLabels.wedgeSolid) {
    edge.points = _.chunk(edge.points.split(/,|\s/).map(p => Number(p)), 2)
    return edge
  }

  for (const pos of ['x1', 'y1', 'x2', 'y2']) {
    edge[pos] = Number(edge[pos])
  }

  return edge
}

Renderer.prototype.positionInfoFromSvgXml = async function(page, xml) {
  // aneb: need to open browser, getBBox is not available via jsdom as it does not render
  await page.setContent(xml, this.waitOptions)

  const dom = await page.evaluate(getPositionInfoFromSvg)
  dom.edges = dom.edges.map(e => this.makeEdgeAttributesNumeric(e))

  return {
    dom,
    xml
  }
}

Renderer.prototype.updateXmlAttributes = function(attributes) {
  const update = ['x', 'y', 'r', 'x1', 'y1', 'x2', 'y2', 'cx', 'cy']

  for (const attr of update) {
    if (!attributes[attr]) {
      continue
    }

    attributes[attr] = Number(attributes[attr]).toFixed(4)
  }

  if (attributes.points) {
    attributes.points = _.chunk(attributes.points.split(/,|\s/).map(n => Number(n).toFixed(4)), 2).join(' ')
  }

  return attributes
}

Renderer.prototype.updateXmlNode = function(node) {
  if (node.attributes) {
    node.attributes = this.updateXmlAttributes(node.attributes)
  }

  if (node.elements && node.elements.length) {
    node.elements = node.elements.map(c => this.updateXmlNode(c))
  }

  return node
}

Renderer.prototype.cleanupLabel = function(label) {
  delete label.style

  if (label.points) {
    label.points = _.chunk(label.points.split(/,|\s/).map(p => Number(p)), 2)
    return label
  }

  if (label.cx && label.cy) {
    label.points = [[label.cx, label.cy]]
    delete label.cx && delete label.cy && delete label.r
    return label
  }

  throw new Error('the label is neither a polygon nor a point!')
}

Renderer.prototype.groupLabels = function(labels) {
  const groups = _.groupBy(labels, 'label-id')
  const result = []
  for (const [id, elementLabels] of Object.entries(groups)) {
    const text = elementLabels[0].text
    const label = elementLabels[0].label
    const xy = elementLabels.map(p => p.xy.toString()).join(' ')
    result.push({ id, label, xy, text })
  }

  return _.sortBy(result, 'id')
}

Renderer.prototype.saveResizedImage = async function(page, smiles, graph, svg, fileName, quality, jsonOnly = false) {
  await page.setContent(svg, this.waitOptions)
  // TODO aneb: images are not resized anymore, clean up code later
  let [updatedSvg, labels, matrix] = await page.evaluate(resizeImage)

  await page.setContent(updatedSvg, this.waitOptions)
  updatedSvg = await page.evaluate(drawMasksAroundTextElements)

  const ops = []

  if (!jsonOnly) {
    const updatedSvgElement = await page.$('svg')
    const capture = updatedSvgElement.screenshot({
      path: `${fileName}.jpg`,
      omitBackground: false,
      quality: quality
    })
    ops.push(capture)
  }

  if (this.outputLabels && labels.length) {
    labels = labels
      .map(l => this.cleanupLabel(l))
      .map(l => ({ ...l, xy: this.svgHelper.transformPoints(l, matrix) }))

    labels = this.groupLabels(labels)

    const elementGraph = getElementGraph(graph, labels)
    const result = { labels, smiles, elementGraph }

    // ops.push(fs.writeFile(`${fileName}-meta.json`, JSON.stringify({ smiles }, null, 2)))
    ops.push(fs.writeFile(`${fileName}.json`, JSON.stringify(result, null, 2)))
  }

  if (this.outputSvg) {
    const updatedSvgXml = js2xml(this.updateXmlNode(xml2js(updatedSvg)), {
      spaces: 2,
      compact: false
    })
    ops.push(fs.writeFile(`${fileName}-after.svg`, updatedSvgXml))
  }

  await Promise.all(ops)
}

Renderer.prototype.smilesToSvgXml = function(smiles) {
  const tree = this.parser.parse(smiles)

  // aneb: need to keep layout relatively constant
  const baseValue = Math.round(this.size * 0.1)

  const options = {
    overlapSensitivity: 1e-1,
    overlapResolutionIterations: 50,
    strokeWidth: _.random(5, 10),
    gradientOffset: _.random(0, 100),
    wedgeBaseWidth: baseValue * _.random(0.2, 0.4),
    dashedWedgeSpacing: baseValue * _.random(0.04, 0.05),
    dashedWedgeWidth: baseValue * _.random(0.6, 0.8),
    bondThickness: baseValue * _.random(0.1, 0.2),
    bondLength: baseValue * _.random(2, 4),
    shortBondLength: _.random(0.7, 0.85),
    bondSpacing: baseValue * _.random(0.2, 0.5),
    font: _.sample(this.fonts),
    fontWeight: _.sample(this.fontWeights),
    fontSizeLarge: baseValue * _.random(0.8, 0.95),
    fontSizeSmall: baseValue * _.random(0.5, 0.65),
    padding: baseValue * _.random(2, 6),
    terminalCarbons: _.sample([true, false]),
    explicitHydrogens: _.sample([true, false])
  }

  const mono = { C: '#000', BACKGROUND: '#fff' }
  const random = this.randomColorMap(Object.keys(this.colorMap))

  const backgroundColor = '#' + _.random(200, 255).toString(16).repeat(3)
  const randomWithWhiteBackGround = { ...random, BACKGROUND: backgroundColor }

  const colormaps = [this.colorMap, mono, randomWithWhiteBackGround]
  const colors = _.sample(colormaps)
  const style = `stroke-width: 0px; background-color: ${colors.BACKGROUND}`
  const svg = this.document.createElementNS('http://www.w3.org/2000/svg', 'svg')
  const drawer = new SvgDrawer({ colors, options })

  drawer.draw(tree, svg)

  // aneb: must set other properties after drawing
  this.svgHelper.update(svg, { style, smiles })

  return [this.XMLSerializer.serializeToString(svg), drawer.preprocessor.graph]
}

Renderer.prototype.getCornersAligned = function({ x, y, width: w, height: h }) {
  const dx = w
  const dy = h

  // aneb: they are drawn in the order that is defined, so it is closed at (x,y) again
  return [
    [x, y],
    [x + dx, y],
    [x + dx, y + dy],
    [x, y + dy]
  ]
}

Renderer.prototype.getCornersOriented = function(edge) {
  // aneb: wedge is already drawn as polygon, all others are just lines, get polygon around lines and then treat both equally
  if (edge.label === bondLabels.wedgeSolid) {
    return [edge.points]
  }

  return this.svgHelper.getEdgePointsOfBoxAroundLine(edge)
}

Renderer.prototype.drawPoints = function({ id, label, points, text }) {
  const color = this.svgHelper.randomColor()

  // aneb: try to avoid overlapping points by using different sizes
  const size = _.random(5, 10)
  return points.map(([x, y]) => {
    return this.svgHelper.createElement('circle', {
      'label-id': `${id}`,
      label: label,
      text: text,
      cx: x,
      cy: y,
      r: size,
      style: this.color(color, true)
    })
  })
}

Renderer.prototype.addLabels = function({ dom, xml }) {
  const svg = new JSDOM(xml).window.document.documentElement.querySelector('svg')

  const nodeCorners = dom.nodes.map(n => ({ ...n, points: this.getCornersAligned(n) }))
  const nodeLabels = nodeCorners.map(n => this.drawPoints(n))
  const edgeLabels = []

  const points = dom.edges.map(e => ({ ...e, points: this.getCornersOriented(e) })).filter(e => !!e.points)
  const hull = Object.values(_.groupBy(points, 'id')).map(e => this.svgHelper.hull(e))
  const hullBox = hull.map(edge => this.drawPoints(edge))
  edgeLabels.push(hullBox)

  this.svgHelper.appendChildren(svg, [...nodeLabels, ...edgeLabels])

  return this.XMLSerializer.serializeToString(svg)
}

Renderer.prototype.imageFromSmilesString = async function(page, smiles) {
  const [svgXmlWithoutLabels, graph] = this.smilesToSvgXml(smiles)
  const { dom, xml } = await this.positionInfoFromSvgXml(page, svgXmlWithoutLabels)

  // aneb: these are only at the original size, the final labels are computed after image has been resized
  const svgXmlWithLabels = this.addLabels({ dom, xml })
  const id = this.id(smiles)

  if (!this.outputFlat) {
    const target = `${this.outputDirectory}/${id}`
    await fs.ensureDir(target)
    await this.saveResizedImage(page, smiles, graph, svgXmlWithoutLabels, `${target}/x`, 100, false)
    await this.saveResizedImage(page, smiles, graph, svgXmlWithLabels, `${target}/y`, 100, false)
    return
  }

  // aneb: debugging only
  await this.saveResizedImage(page, smiles, graph, svgXmlWithoutLabels, `${this.outputDirectory}/${id}-x`, 100, false)
  await this.saveResizedImage(page, smiles, graph, svgXmlWithLabels, `${this.outputDirectory}/${id}-y`, 100, false)
}

module.exports = Renderer
