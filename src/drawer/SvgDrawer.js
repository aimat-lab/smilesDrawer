// we use the drawer to do all the preprocessing. then we take over the drawing
// portion to output to svg
const ArrayHelper = require('./ArrayHelper')
const Atom = require('./Atom')
const Drawer = require('./Drawer')
const Line = require('./Line')
const SvgWrapper = require('./SvgWrapper')
const MathHelper = require('./MathHelper')
const Vector2 = require('./Vector2')
const { bondLabels } = require('../generator/types')

class SvgDrawer {
  constructor({ colors, options }) {
    this.colors = colors
    this.preprocessor = new Drawer(options)
    this.opts = this.preprocessor.opts
  }

  /**
   * Draws the parsed smiles data to an svg element.
   *
   * @param {Object} data The tree returned by the smiles parser.
   * @param {(String|HTMLElement)} target The id of the HTML svg element the structure is drawn to - or the element itself.

   * @returns {Object} The dimensions of the drawing in { width, height }
   */
  draw(data, target) {
    const preprocessor = this.preprocessor

    preprocessor.initDraw(data)

    this.svgWrapper = new SvgWrapper(target, this.preprocessor.opts, this.colors)

    preprocessor.processGraph()

    // Set the canvas to the appropriate size
    this.svgWrapper.determineDimensions(preprocessor.graph.vertices)

    // Do the actual drawing
    this.drawEdges()
    this.drawVertices()

    return this.svgWrapper.constructSvg()
  }

  /**
   * Draws a ring inside a provided ring, indicating aromaticity.
   * @param {Ring} ring A ring.
   */
  drawAromaticityRing(ring) {
    const r = MathHelper.apothemFromSideLength(this.opts.bondLength, ring.getSize())
    this.svgWrapper.drawRing(ring.center.x, ring.center.y, r * 0.8)
  }

  /**
   * Draw the actual edges as bonds.
   */
  drawEdges() {
    const preprocessor = this.preprocessor
    const graph = preprocessor.graph
    const rings = preprocessor.rings
    const drawn = Array(this.preprocessor.graph.edges.length).fill(false)

    graph.traverseBF(0, vertex => {
      const edges = graph.getEdges(vertex.id)
      for (let i = 0; i < edges.length; i++) {
        const edgeId = edges[i]
        if (!drawn[edgeId]) {
          drawn[edgeId] = true
          this.drawEdge(edgeId)
        }
      }
    })

    for (let i = 0; i < rings.length; i++) {
      const ring = rings[i]
      if (preprocessor.isRingAromatic(ring)) {
        this.drawAromaticityRing(ring)
      }
    }
  }

  /**
   * Draw the an edge as a bond.
   * @param {Number} edgeId An edge id.
   */
  drawEdge(edgeId) {
    const preprocessor = this.preprocessor
    const opts = preprocessor.opts
    const svgWrapper = this.svgWrapper
    const edge = preprocessor.graph.edges[edgeId]
    const vertexA = preprocessor.graph.vertices[edge.sourceId]
    const vertexB = preprocessor.graph.vertices[edge.targetId]
    const elementA = vertexA.value.element
    const elementB = vertexB.value.element

    if ((!vertexA.value.isDrawn || !vertexB.value.isDrawn) && preprocessor.opts.atomVisualization === 'default') {
      return
    }

    const edgeIdLabel = 'edge-id'
    const edgeIdValue = `edge-id-${edge.id}`

    const a = vertexA.position
    const b = vertexB.position
    const normals = preprocessor.getEdgeNormals(edge)
    const sides = ArrayHelper.clone(normals)

    sides[0].multiplyScalar(10).add(a)
    sides[1].multiplyScalar(10).add(a)

    // aneb: edge.isPartOfAromaticRing is always false when not drawing ring
    if (edge.bondType === '=' ||
      preprocessor.getRingbondType(vertexA, vertexB) === '=' ||
      (edge.isPartOfAromaticRing && preprocessor.bridgedRing)) {
      // Always draw double bonds inside the ring
      const inRing = preprocessor.areVerticesInSameRing(vertexA, vertexB)
      const s = preprocessor.chooseSide(vertexA, vertexB, sides)

      const label = edge.isPartOfAromaticRing ? bondLabels.aromaticDouble : bondLabels.double

      if (inRing) {
        // Always draw double bonds inside a ring
        // if the bond is shared by two rings, it is drawn in the larger
        // problem: smaller ring is aromatic, bond is still drawn in larger -> fix this
        const lcr = preprocessor.getLargestOrAromaticCommonRing(vertexA, vertexB)
        const center = lcr.center

        normals[0].multiplyScalar(opts.bondSpacing)
        normals[1].multiplyScalar(opts.bondSpacing)

        const line = center.sameSideAs(vertexA.position, vertexB.position, Vector2.add(a, normals[0]))
          ? new Line(Vector2.add(a, normals[0]), Vector2.add(b, normals[0]), elementA, elementB)
          : new Line(Vector2.add(a, normals[1]), Vector2.add(b, normals[1]), elementA, elementB)

        line.shorten(opts.bondLength - opts.shortBondLength * opts.bondLength)

        // aneb: this is the inner double bond line, it should not be drawn when ring is drawn
        // this is still not 100% correct, whole ring is drawn with double bonds
        if (!edge.isPartOfAromaticRing) {
          svgWrapper.drawLine(edgeIdLabel, edgeIdValue, label, line)
        }

        svgWrapper.drawLine(edgeIdLabel, edgeIdValue, label, new Line(a, b, elementA, elementB))
        return
      }

      if ((edge.center ||
        (vertexA.isTerminal() && vertexB.isTerminal())) ||
        ((s.anCount === 0 && s.bnCount > 1) || (s.bnCount === 0 && s.anCount > 1))) {
        this.multiplyNormals(normals, opts.halfBondSpacing)

        const lineA = new Line(Vector2.add(a, normals[0]), Vector2.add(b, normals[0]), elementA, elementB)
        const lineB = new Line(Vector2.add(a, normals[1]), Vector2.add(b, normals[1]), elementA, elementB)

        svgWrapper.drawLine(edgeIdLabel, edgeIdValue, bondLabels.double, lineA)
        svgWrapper.drawLine(edgeIdLabel, edgeIdValue, bondLabels.double, lineB)
        return
      }

      if ((s.sideCount[0] > s.sideCount[1]) || (s.totalSideCount[0] > s.totalSideCount[1])) {
        this.multiplyNormals(normals, opts.bondSpacing)

        const line = new Line(Vector2.add(a, normals[0]), Vector2.add(b, normals[0]), elementA, elementB)

        line.shorten(opts.bondLength - opts.shortBondLength * opts.bondLength)

        svgWrapper.drawLine(edgeIdLabel, edgeIdValue, bondLabels.double, line)
        svgWrapper.drawLine(edgeIdLabel, edgeIdValue, bondLabels.double, new Line(a, b, elementA, elementB))
        return
      }

      if ((s.sideCount[0] < s.sideCount[1]) || (s.totalSideCount[0] <= s.totalSideCount[1])) {
        this.multiplyNormals(normals, opts.bondSpacing)

        const line = new Line(Vector2.add(a, normals[1]), Vector2.add(b, normals[1]), elementA, elementB)

        line.shorten(opts.bondLength - opts.shortBondLength * opts.bondLength)
        svgWrapper.drawLine(edgeIdLabel, edgeIdValue, bondLabels.double, line)
        svgWrapper.drawLine(edgeIdLabel, edgeIdValue, bondLabels.double, new Line(a, b, elementA, elementB))
        return
      }
    }

    if (edge.bondType === '#') {
      normals[0].multiplyScalar(opts.bondSpacing / 1.5)
      normals[1].multiplyScalar(opts.bondSpacing / 1.5)

      const lineA = new Line(Vector2.add(a, normals[0]), Vector2.add(b, normals[0]), elementA, elementB)
      const lineB = new Line(Vector2.add(a, normals[1]), Vector2.add(b, normals[1]), elementA, elementB)

      svgWrapper.drawLine(edgeIdLabel, edgeIdValue, bondLabels.triple, lineA)
      svgWrapper.drawLine(edgeIdLabel, edgeIdValue, bondLabels.triple, lineB)
      svgWrapper.drawLine(edgeIdLabel, edgeIdValue, bondLabels.triple, new Line(a, b, elementA, elementB))
      return
    }

    if (edge.bondType === '.') {
      return
    }

    // aneb: from here on everything is a "single" bond
    const isChiralCenterA = vertexA.value.isStereoCenter
    const isChiralCenterB = vertexB.value.isStereoCenter

    if (edge.wedge === 'up') {
      svgWrapper.drawWedge(edgeIdLabel, edgeIdValue, bondLabels.wedgeSolid, new Line(a, b, elementA, elementB, isChiralCenterA, isChiralCenterB))
      return
    }

    if (edge.wedge === 'down') {
      svgWrapper.drawDashedWedge(edgeIdLabel, edgeIdValue, bondLabels.wedgeDashed, new Line(a, b, elementA, elementB, isChiralCenterA, isChiralCenterB))
      return
    }

    const label = edge.isPartOfAromaticRing ? bondLabels.aromaticSingle : bondLabels.single
    svgWrapper.drawLine(edgeIdLabel, edgeIdValue, label, new Line(a, b, elementA, elementB, isChiralCenterA, isChiralCenterB))
  }

  /**
   * Draws the vertices representing atoms to the canvas.
   *
   */
  drawVertices() {
    const preprocessor = this.preprocessor
    const opts = preprocessor.opts
    const graph = preprocessor.graph
    const svgWrapper = this.svgWrapper

    const vertexIdLabel = 'vertex-id'

    for (let i = 0; i < graph.vertices.length; i++) {
      const vertex = graph.vertices[i]
      const atom = vertex.value
      let charge = 0
      let isotope = 0
      const bondCount = vertex.value.bondCount
      const element = atom.element
      let hydrogens = Atom.maxBonds[element] - bondCount
      const dir = vertex.getTextDirection(graph.vertices)
      const isTerminal = opts.terminalCarbons || element !== 'C' || atom.hasAttachedPseudoElements ? vertex.isTerminal() : false
      const isCarbon = atom.element === 'C'

      const vertexIdValue = `vertex-id-${i}`
      const vertexLabel = 'element-text'

      // This is a HACK to remove all hydrogens from nitrogens in aromatic rings, as this
      // should be the most common state. This has to be fixed by kekulization
      if (atom.element === 'N' && atom.isPartOfAromaticRing) {
        hydrogens = 0
      }

      if (atom.bracket) {
        hydrogens = atom.bracket.hcount
        charge = atom.bracket.charge
        isotope = atom.bracket.isotope
      }

      if (opts.atomVisualization === 'allballs') {
        svgWrapper.drawBall(vertexIdLabel, vertexIdValue, vertexLabel, vertex.position.x, vertex.position.y, element)
      } else if ((atom.isDrawn &&
        (!isCarbon || atom.drawExplicit || isTerminal || atom.hasAttachedPseudoElements)) || graph.vertices.length === 1) {
        if (opts.atomVisualization === 'default') {
          svgWrapper.drawText(vertexIdLabel, vertexIdValue, vertexLabel, vertex.position.x, vertex.position.y,
            element, hydrogens, dir, isTerminal, charge, isotope, atom.getAttachedPseudoElements())
        } else if (opts.atomVisualization === 'balls') {
          svgWrapper.drawBall(vertexIdLabel, vertexIdValue, vertexLabel, vertex.position.x, vertex.position.y, element)
        }
      } else if (vertex.getNeighbourCount() === 2 && vertex.forcePositioned === true) {
        // If there is a carbon which bonds are in a straight line, draw a dot
        const a = graph.vertices[vertex.neighbours[0]].position
        const b = graph.vertices[vertex.neighbours[1]].position
        const angle = Vector2.threePointangle(vertex.position, a, b)

        if (Math.abs(Math.PI - angle) < 0.1) {
          svgWrapper.drawPoint(vertexIdLabel, vertexIdValue, vertexLabel, vertex.position.x, vertex.position.y, element)
        }
      }
    }
  }

  /**
   * Returns the total overlap score of the current molecule.
   * @returns {Number} The overlap score.
   */
  getTotalOverlapScore() {
    return this.preprocessor.getTotalOverlapScore()
  }

  /**
   * Returns the molecular formula of the loaded molecule as a string.
   * @returns {String} The molecular formula.
   */
  getMolecularFormula() {
    return this.preprocessor.getMolecularFormula()
  }

  /**
   * @param {Array} normals list of normals to multiply
   * @param {Number} spacing value to multiply normals by
   */
  multiplyNormals(normals, spacing) {
    normals[0].multiplyScalar(spacing)
    normals[1].multiplyScalar(spacing)
  }
}

module.exports = SvgDrawer
