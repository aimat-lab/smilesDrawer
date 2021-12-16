const crypto = require('crypto')
const fs = require('fs-extra')

const _ = require('lodash')
const { JSDOM } = require('jsdom')
const { xml2js, js2xml } = require('xml-js')

const Parser = require('../drawer/Parser')
const SvgDrawer = require('../drawer/SvgDrawer')
const SVG = require('./SVG')
const { bondLabels, labelTypes } = require('./types')
const { getPositionInfoFromSvg, resizeImage, drawMasksAroundTextElements } = require('./browser')

function Renderer({ outputDirectory, size, fonts, fontWeights, preserveAspectRatio, concurrency, labelType, segment, outputSvg, outputLabels, outputFlat }) {
  // aneb: find out why this does not work in above scope ...
  const colorMap = require('./colors')

  this.parser = Parser
  this.outputDirectory = outputDirectory
  this.size = size
  this.fonts = fonts
  this.fontWeights = fontWeights
  this.preserveAspectRatio = preserveAspectRatio
  this.colorMap = colorMap
  this.concurrency = concurrency
  this.labelType = labelType
  this.segment = segment
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
  const fill = this.segment || circle ? color : 'none'
  return `fill: ${fill}; stroke: ${color};`
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

Renderer.prototype.saveResizedImage = async function(page, smiles, svg, fileName, quality, jsonOnly = false) {
  await page.setContent(svg, this.waitOptions)
  let [updatedSvg, labels, matrix] = await page.evaluate(resizeImage, { size: this.size, preserveAspectRatio: this.preserveAspectRatio })

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
    const cleanLabels = labels
      .map(l => this.cleanupLabel(l))
      .map(l => ({ ...l, xy: this.svgHelper.transformPoints(l, matrix) }))

    const finalLabels = this.groupLabels(cleanLabels)

    // ops.push(fs.writeFile(`${fileName}-meta.json`, JSON.stringify({ smiles }, null, 2)))
    ops.push(fs.writeFile(`${fileName}.json`, JSON.stringify(finalLabels, null, 2)))
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
    overlapSensitivity: 1e-5,
    overlapResolutionIterations: 50,
    strokeWidth: _.sample([5, 6, 7, 8, 9, 10]),
    gradientOffset: _.sample([10, 20, 30, 40, 50, 60, 70, 80, 90, 100]),
    wedgeBaseWidth: baseValue * _.sample([0.2, 0.3, 0.4, 0.5]),
    dashedWedgeSpacing: baseValue * _.sample([0.07, 0.08, 0.09, 0.1]),
    dashedWedgeWidth: baseValue * _.sample([0.6, 0.7, 0.7, 0.8, 0.9]),
    bondThickness: baseValue * _.sample([0.1, 0.15, 0.2]),
    bondLength: baseValue * _.sample([2, 2.5, 3, 3.5, 4]),
    shortBondLength: _.sample([0.7, 0.75, 0.8, 0.85]),
    bondSpacing: baseValue * _.sample([0.2, 0.25, 0.5]),
    font: _.sample(this.fonts),
    fontWeight: _.sample(this.fontWeights),
    fontSizeLarge: baseValue * _.sample([0.8, 0.85, 0.9, 0.95]),
    fontSizeSmall: baseValue * _.sample([0.5, 0.55, 0.6, 0.65]),
    padding: baseValue * _.sample([5, 7.5, 10, 12.5, 15]),
    terminalCarbons: _.sample([true, false]),
    explicitHydrogens: _.sample([true, false])
  }

  const colors = this.colorMap
  const style = `stroke-width: 0px; background-color: ${colors.BACKGROUND}`
  const svg = this.document.createElementNS('http://www.w3.org/2000/svg', 'svg')
  const drawer = new SvgDrawer({ colors, options })

  drawer.draw(tree, svg)

  // aneb: must set other properties after drawing
  this.svgHelper.update(svg, { style, smiles })

  return this.XMLSerializer.serializeToString(svg)
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
  const size = _.floor(_.random(true) * 5 + 2) / 10
  return points.map(([x, y]) => {
    return this.svgHelper.createElement('circle', {
      'label-id': `${id}-label`,
      label: label,
      text: text,
      cx: x,
      cy: y,
      r: size,
      style: this.color(color, true)
    })
  })
}

Renderer.prototype.drawSinglePolygon = function({ id, label, points, text }) {
  const color = this.svgHelper.randomColor()
  return this.svgHelper.createElement('polygon', {
    'label-id': `${id}-label`,
    label: label,
    text: text,
    points: points.join(' '),
    style: this.color(color)
  })
}

Renderer.prototype.drawMultiPolygon = function(edgeElements) {
  const color = this.svgHelper.randomColor()
  const { id, label } = edgeElements[0]
  const points = edgeElements.map(e => this.getCornersOriented(e)).filter(e => !!e)
  return points.map(point => {
    return this.svgHelper.createElement('polygon', {
      'label-id': `${id}-label`,
      label: label,
      points: point.join(' '),
      style: this.color(color)
    })
  })
}

Renderer.prototype.addLabels = function({ dom, xml }) {
  const svg = new JSDOM(xml).window.document.documentElement.querySelector('svg')

  const nodeCorners = dom.nodes.map(n => ({ ...n, points: this.getCornersAligned(n) }))
  const nodeLabels = this.labelType === labelTypes.points && !this.segment
    ? nodeCorners.map(n => this.drawPoints(n))
    : nodeCorners.map(n => this.drawSinglePolygon(n))

  const edgeLabels = []

  if (this.labelType === labelTypes.box) {
    const correctedEdges = dom.edges.map(e => ({ ...e, ...this.svgHelper.correctBoundingBox(e) }))
    const merged = this.svgHelper.mergeBoundingBoxes(correctedEdges)
    const mergedWithPoints = merged.map(n => ({ ...n, points: this.getCornersAligned(n) }))
    edgeLabels.push(...mergedWithPoints.map(e => this.drawSinglePolygon(e)))
  }

  if (this.labelType === labelTypes.oriented) {
    const groupedEdges = _.groupBy(dom.edges, 'id')
    const tightBoxes = Object.values(groupedEdges).map(edge => this.drawMultiPolygon(edge))
    edgeLabels.push(tightBoxes)
  }

  if (this.labelType === labelTypes.points) {
    const points = dom.edges.map(e => ({ ...e, points: this.getCornersOriented(e) })).filter(e => !!e.points)
    const hull = Object.values(_.groupBy(points, 'id')).map(e => this.svgHelper.hull(e))
    const hullBox = this.segment
      ? hull.map(edge => this.drawSinglePolygon(edge))
      : hull.map(edge => this.drawPoints(edge))

    edgeLabels.push(hullBox)
  }

  this.svgHelper.appendChildren(svg, [...nodeLabels, ...edgeLabels])

  return this.XMLSerializer.serializeToString(svg)
}

Renderer.prototype.imageFromSmilesString = async function(page, smiles) {
  const svgXmlWithoutLabels = this.smilesToSvgXml(smiles)
  const { dom, xml } = await this.positionInfoFromSvgXml(page, svgXmlWithoutLabels)

  // aneb: these are only at the original size, the final labels are computed after image has been resized
  const svgXmlWithLabels = this.addLabels({ dom, xml })
  const id = this.id(smiles)
  const quality = _.random(10, 80)

  if (!this.outputFlat) {
    const target = `${this.outputDirectory}/${id}`
    await fs.ensureDir(target)
    await this.saveResizedImage(page, smiles, svgXmlWithoutLabels, `${target}/x`, quality, false)
    await this.saveResizedImage(page, smiles, svgXmlWithLabels, `${target}/y`, 100, true)
    return
  }

  // aneb: debugging only
  await this.saveResizedImage(page, smiles, svgXmlWithoutLabels, `${this.outputDirectory}/${id}-x`, quality, false)
  await this.saveResizedImage(page, smiles, svgXmlWithLabels, `${this.outputDirectory}/${id}-y`, 100, true)
}

module.exports = Renderer
