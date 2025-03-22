/* eslint-disable no-mixed-operators,array-callback-return */
// @ts-check
const MathHelper = require('./MathHelper')
const ArrayHelper = require('./ArrayHelper')
const Vector2 = require('./Vector2')
const Edge = require('./Edge')
const Atom = require('./Atom')
const Ring = require('./Ring')
const RingConnection = require('./RingConnection')
const Graph = require('./Graph')
const SSSR = require('./SSSR')

/**
 * The main class of the application representing the smiles drawer
 *
 * @property {Graph} graph The graph associated with this SmilesDrawer.Drawer instance.
 * @property {Number} ringIdCounter An internal counter to keep track of ring ids.
 * @property {Number} ringConnectionIdCounter An internal counter to keep track of ring connection ids.
 * @property {Number} totalOverlapScore The current internal total overlap score.
 * @property {Object} defaultOptions The default options.
 * @property {Object} opts The merged options.
 */
class Drawer {
  /**
   * The constructor for the class SmilesDrawer.
   *
   * @param {Object} options An object containing custom values for different options. It is merged with the default options.
   */
  constructor(options) {
    this.graph = null
    this.doubleBondConfigCount = 0
    this.doubleBondConfig = null
    this.ringIdCounter = 0
    this.ringConnectionIdCounter = 0
    this.totalOverlapScore = 0

    this.defaultOptions = {
      strokeWidth: 1,
      gradientOffset: 25,
      wedgeBaseWidth: 2,
      dashedWedgeSpacing: 1.25,
      dashedWedgeWidth: 1.5,
      bondThickness: 0.6,
      bondLength: 25,
      shortBondLength: 0.85,
      bondSpacing: 0.18 * 20,
      font: 'Helvetica, Arial, sans-serif',
      fontSizeLarge: 7,
      fontSizeSmall: 4,
      padding: 25,
      atomVisualization: 'default',
      isomeric: true,
      terminalCarbons: false,
      explicitHydrogens: true,
      overlapSensitivity: 0.5,
      overlapResolutionIterations: 10,
      compactDrawing: true,
      experimentalSSSR: false,
      kkThreshold: 0.1,
      kkInnerThreshold: 0.1,
      kkMaxIteration: 2000,
      kkMaxInnerIteration: 50,
      kkMaxEnergy: 1e9
    }

    this.opts = Object.assign(this.defaultOptions, options)
    this.opts.halfBondSpacing = this.opts.bondSpacing / 2.0
    this.opts.bondLengthSq = this.opts.bondLength * this.opts.bondLength
  }

  initDraw(data) {
    this.ringIdCounter = 0
    this.ringConnectionIdCounter = 0

    this.graph = new Graph(data, this.opts.isomeric)
    this.rings = []
    this.ringConnections = []

    this.originalRings = []
    this.originalRingConnections = []

    this.bridgedRing = false

    // Reset those, in case the previous drawn SMILES had a dangling \ or /
    this.doubleBondConfigCount = null
    this.doubleBondConfig = null

    this.initRings()
    this.initHydrogens()
  }

  processGraph() {
    this.position()

    // Restore the ring information (removes bridged rings and replaces them with the original, multiple, rings)
    this.restoreRingInformation()

    // Atoms bonded to the same ring atom
    this.resolvePrimaryOverlaps()

    let overlapScore = this.getOverlapScore()

    this.totalOverlapScore = this.getOverlapScore().total

    for (let o = 0; o < this.opts.overlapResolutionIterations; o++) {
      for (let i = 0; i < this.graph.edges.length; i++) {
        const edge = this.graph.edges[i]
        if (this.isEdgeRotatable(edge)) {
          const subTreeDepthA = this.graph.getTreeDepth(edge.sourceId, edge.targetId)
          const subTreeDepthB = this.graph.getTreeDepth(edge.targetId, edge.sourceId)

          // Only rotate the shorter subtree
          let a = edge.targetId
          let b = edge.sourceId

          if (subTreeDepthA > subTreeDepthB) {
            a = edge.sourceId
            b = edge.targetId
          }

          const subTreeOverlap = this.getSubtreeOverlapScore(b, a, overlapScore.vertexScores)
          if (subTreeOverlap.value > this.opts.overlapSensitivity) {
            const vertexA = this.graph.vertices[a]
            const vertexB = this.graph.vertices[b]
            const neighboursB = vertexB.getNeighbours(a)

            if (neighboursB.length === 1) {
              const neighbour = this.graph.vertices[neighboursB[0]]
              const angle = neighbour.position.getRotateAwayFromAngle(vertexA.position, vertexB.position, MathHelper.toRad(120))

              this.rotateSubtree(neighbour.id, vertexB.id, angle, vertexB.position)
              // If the new overlap is bigger, undo change
              const newTotalOverlapScore = this.getOverlapScore().total

              if (newTotalOverlapScore > this.totalOverlapScore) {
                this.rotateSubtree(neighbour.id, vertexB.id, -angle, vertexB.position)
              } else {
                this.totalOverlapScore = newTotalOverlapScore
              }
            } else if (neighboursB.length === 2) {
              // Switch places / sides
              // If vertex a is in a ring, do nothing
              if (vertexB.value.rings.length !== 0 && vertexA.value.rings.length !== 0) {
                continue
              }

              const neighbourA = this.graph.vertices[neighboursB[0]]
              const neighbourB = this.graph.vertices[neighboursB[1]]

              if (neighbourA.value.rings.length === 1 && neighbourB.value.rings.length === 1) {
                // Both neighbours in same ring. TODO: does this create problems with wedges? (up = down and vice versa?)
                if (neighbourA.value.rings[0] !== neighbourB.value.rings[0]) {
                  continue
                }
                // TODO: Rotate circle
              } else if (neighbourA.value.rings.length !== 0 || neighbourB.value.rings.length !== 0) {
                continue
              } else {
                const angleA = neighbourA.position.getRotateAwayFromAngle(vertexA.position, vertexB.position, MathHelper.toRad(120))
                const angleB = neighbourB.position.getRotateAwayFromAngle(vertexA.position, vertexB.position, MathHelper.toRad(120))

                this.rotateSubtree(neighbourA.id, vertexB.id, angleA, vertexB.position)
                this.rotateSubtree(neighbourB.id, vertexB.id, angleB, vertexB.position)

                const newTotalOverlapScore = this.getOverlapScore().total

                if (newTotalOverlapScore > this.totalOverlapScore) {
                  this.rotateSubtree(neighbourA.id, vertexB.id, -angleA, vertexB.position)
                  this.rotateSubtree(neighbourB.id, vertexB.id, -angleB, vertexB.position)
                } else {
                  this.totalOverlapScore = newTotalOverlapScore
                }
              }
            }

            overlapScore = this.getOverlapScore()
          }
        }
      }
    }

    this.resolveSecondaryOverlaps(overlapScore.scores)

    if (this.opts.isomeric) {
      this.annotateStereochemistry()
    }

    // Initialize pseudo elements or shortcuts
    if (this.opts.compactDrawing && this.opts.atomVisualization === 'default') {
      this.initPseudoElements()
    }

    this.rotateDrawing()
  }

  /**
   * Returns the total overlap score of the current molecule.
   *
   * @returns {Number} The overlap score.
   */
  getTotalOverlapScore() {
    return this.totalOverlapScore
  }

  /**
   * Returns the molecular formula of the loaded molecule as a string.
   *
   * @returns {String} The molecular formula.
   */
  getMolecularFormula() {
    let molecularFormula = ''
    const counts = new Map()

    // Initialize element count
    for (let i = 0; i < this.graph.vertices.length; i++) {
      const atom = this.graph.vertices[i].value

      if (counts.has(atom.element)) {
        counts.set(atom.element, counts.get(atom.element) + 1)
      } else {
        counts.set(atom.element, 1)
      }

      // Hydrogens attached to a chiral center were added as vertices,
      // those in non chiral brackets are added here
      if (atom.bracket && !atom.bracket.chirality) {
        if (counts.has('H')) {
          counts.set('H', counts.get('H') + atom.bracket.hcount)
        } else {
          counts.set('H', atom.bracket.hcount)
        }
      }

      // Add the implicit hydrogens according to valency, exclude
      // bracket atoms as they were handled and always have the number
      // of hydrogens specified explicitly
      if (!atom.bracket) {
        let nHydrogens = Atom.maxBonds[atom.element] - atom.bondCount

        if (atom.isPartOfAromaticRing) {
          nHydrogens--
        }

        if (counts.has('H')) {
          counts.set('H', counts.get('H') + nHydrogens)
        } else {
          counts.set('H', nHydrogens)
        }
      }
    }

    if (counts.has('C')) {
      const count = counts.get('C')
      molecularFormula += 'C' + (count > 1 ? count : '')
      counts.delete('C')
    }

    if (counts.has('H')) {
      const count = counts.get('H')
      molecularFormula += 'H' + (count > 1 ? count : '')
      counts.delete('H')
    }

    const elements = Object.keys(Atom.atomicNumbers).sort()

    elements.map(e => {
      if (counts.has(e)) {
        const count = counts.get(e)
        molecularFormula += e + (count > 1 ? count : '')
      }
    })

    return molecularFormula
  }

  /**
   * Returns the aromatic or largest ring shared by the two vertices.
   *
   * @param {Vertex} vertexA A vertex.
   * @param {Vertex} vertexB A vertex.
   * @returns {(Ring|null)} If an aromatic common ring exists, that ring, else the largest (non-aromatic) ring, else null.
   */
  getLargestOrAromaticCommonRing(vertexA, vertexB) {
    const commonRings = this.getCommonRings(vertexA, vertexB)
    let maxSize = 0
    let largestCommonRing = null

    for (let i = 0; i < commonRings.length; i++) {
      const ring = this.getRing(commonRings[i])
      const size = ring.getSize()

      if (ring.isBenzeneLike(this.graph.vertices)) {
        return ring
      } else if (size > maxSize) {
        maxSize = size
        largestCommonRing = ring
      }
    }

    return largestCommonRing
  }

  /**
   * When drawing a double bond, choose the side to place the double bond. E.g. a double bond should always been drawn inside a ring.
   *
   * @param {Vertex} vertexA A vertex.
   * @param {Vertex} vertexB A vertex.
   * @param {Vector2[]} sides An array containing the two normals of the line spanned by the two provided vertices.
   * @returns {Object} Returns an object containing the following information: {
          totalSideCount: Counts the sides of each vertex in the molecule, is an array [ a, b ],
          totalPosition: Same as position, but based on entire molecule,
          sideCount: Counts the sides of each neighbour, is an array [ a, b ],
          position: which side to position the second bond, is 0 or 1, represents the index in the normal array. This is based on only the neighbours
          anCount: the number of neighbours of vertexA,
          bnCount: the number of neighbours of vertexB
      }
   */
  chooseSide(vertexA, vertexB, sides) {
    // Check which side has more vertices
    // Get all the vertices connected to the both ends
    const an = vertexA.getNeighbours(vertexB.id)
    const bn = vertexB.getNeighbours(vertexA.id)
    const anCount = an.length
    const bnCount = bn.length

    // All vertices connected to the edge vertexA to vertexB
    const tn = ArrayHelper.merge(an, bn)

    // Only considering the connected vertices
    const sideCount = [0, 0]

    for (let i = 0; i < tn.length; i++) {
      const v = this.graph.vertices[tn[i]].position

      if (v.sameSideAs(vertexA.position, vertexB.position, sides[0])) {
        sideCount[0]++
      } else {
        sideCount[1]++
      }
    }

    // Considering all vertices in the graph, this is to resolve ties
    // from the above side counts
    const totalSideCount = [0, 0]

    for (let i = 0; i < this.graph.vertices.length; i++) {
      const v = this.graph.vertices[i].position

      if (v.sameSideAs(vertexA.position, vertexB.position, sides[0])) {
        totalSideCount[0]++
      } else {
        totalSideCount[1]++
      }
    }

    return {
      totalSideCount: totalSideCount,
      totalPosition: totalSideCount[0] > totalSideCount[1] ? 0 : 1,
      sideCount: sideCount,
      position: sideCount[0] > sideCount[1] ? 0 : 1,
      anCount: anCount,
      bnCount: bnCount
    }
  }

  /**
   * Checks whether or not two vertices are in the same ring.
   *
   * @param {Vertex} vertexA A vertex.
   * @param {Vertex} vertexB A vertex.
   * @returns {Boolean} A boolean indicating whether or not the two vertices are in the same ring.
   */
  areVerticesInSameRing(vertexA, vertexB) {
    // This is a little bit lighter (without the array and push) than
    // getCommonRings().length > 0
    for (let i = 0; i < vertexA.value.rings.length; i++) {
      for (let j = 0; j < vertexB.value.rings.length; j++) {
        if (vertexA.value.rings[i] === vertexB.value.rings[j]) {
          return true
        }
      }
    }

    return false
  }

  /**
   * Check whether or not a ring is an implicitly defined aromatic ring (lower case smiles).
   *
   * @param {Ring} ring A ring.
   * @returns {Boolean} A boolean indicating whether or not a ring is implicitly defined as aromatic.
   */
  isRingAromatic(ring) {
    for (let i = 0; i < ring.members.length; i++) {
      const vertex = this.graph.vertices[ring.members[i]]

      if (!vertex.value.isPartOfAromaticRing) {
        return false
      }
    }

    return true
  }

  /**
   * Get the normals of an edge.
   *
   * @param {Edge} edge An edge.
   * @returns {Vector2[]} An array containing two vectors, representing the normals.
   */
  getEdgeNormals(edge) {
    const v1 = this.graph.vertices[edge.sourceId].position
    const v2 = this.graph.vertices[edge.targetId].position

    // Get the normalized normals for the edge
    return Vector2.units(v1, v2)
  }

  /**
   * Returns the type of the ringbond (e.g. '=' for a double bond). The ringbond represents the break in a ring introduced when creating the MST. If the two vertices supplied as arguments are not part of a common ringbond, the method returns null.
   *
   * @param {Vertex} vertexA A vertex.
   * @param {Vertex} vertexB A vertex.
   * @returns {(String|null)} Returns the ringbond type or null, if the two supplied vertices are not connected by a ringbond.
   */
  getRingbondType(vertexA, vertexB) {
    // Checks whether the two vertices are the ones connecting the ring
    // and what the bond type should be.
    if (vertexA.value.getRingbondCount() < 1 || vertexB.value.getRingbondCount() < 1) {
      return null
    }

    for (let i = 0; i < vertexA.value.ringbonds.length; i++) {
      for (let j = 0; j < vertexB.value.ringbonds.length; j++) {
        // if(i != j) continue;
        if (vertexA.value.ringbonds[i].id === vertexB.value.ringbonds[j].id) {
          // If the bonds are equal, it doesn't matter which bond is returned.
          // if they are not equal, return the one that is not the default ("-")
          if (vertexA.value.ringbonds[i].bondType === '-') {
            return vertexB.value.ringbonds[j].bond
          } else {
            return vertexA.value.ringbonds[i].bond
          }
        }
      }
    }

    return null
  }

  /**
   * A helper method to extend the default options with user supplied ones.
   */
  extend() {
    const that = this
    const extended = {}
    let deep = false
    let i = 0
    const length = arguments.length

    if (Object.prototype.toString.call(arguments[0]) === '[object Boolean]') {
      deep = arguments[0]
      i++
    }

    const merge = function(obj) {
      for (const prop in obj) {
        if (Object.prototype.hasOwnProperty.call(obj, prop)) {
          if (deep && Object.prototype.toString.call(obj[prop]) === '[object Object]') {
            extended[prop] = that.extend(true, extended[prop], obj[prop])
          } else {
            extended[prop] = obj[prop]
          }
        }
      }
    }

    for (; i < length; i++) {
      const obj = arguments[i]
      merge(obj)
    }

    return extended
  };

  /**
   * Returns the number of rings this edge is a part of.
   *
   * @param {Number} edgeId The id of an edge.
   * @returns {Number} The number of rings the provided edge is part of.
   */
  edgeRingCount(edgeId) {
    const edge = this.graph.edges[edgeId]
    const a = this.graph.vertices[edge.sourceId]
    const b = this.graph.vertices[edge.targetId]

    return Math.min(a.value.rings.length, b.value.rings.length)
  }

  /**
   * Returns an array containing the bridged rings associated with this  molecule.
   *
   * @returns {Ring[]} An array containing all bridged rings associated with this molecule.
   */
  getBridgedRings() {
    const bridgedRings = []

    for (let i = 0; i < this.rings.length; i++) {
      if (this.rings[i].isBridged) {
        bridgedRings.push(this.rings[i])
      }
    }

    return bridgedRings
  }

  /**
   * Rotates the drawing to make the widest dimension horizontal.
   */
  rotateDrawing() {
    // Rotate the vertices to make the molecule align horizontally
    // Find the longest distance
    let a = 0
    let b = 0
    let maxDist = 0
    for (let i = 0; i < this.graph.vertices.length; i++) {
      const vertexA = this.graph.vertices[i]

      if (!vertexA.value.isDrawn) {
        continue
      }

      for (let j = i + 1; j < this.graph.vertices.length; j++) {
        const vertexB = this.graph.vertices[j]

        if (!vertexB.value.isDrawn) {
          continue
        }

        const dist = vertexA.position.distanceSq(vertexB.position)

        if (dist > maxDist) {
          maxDist = dist
          a = i
          b = j
        }
      }
    }

    let angle = -Vector2.subtract(this.graph.vertices[a].position, this.graph.vertices[b].position).angle()

    if (!isNaN(angle)) {
      // Round to 30 degrees
      const remainder = angle % 0.523599

      // Round either up or down in 30 degree steps
      if (remainder < 0.2617995) {
        angle = angle - remainder
      } else {
        angle += 0.523599 - remainder
      }

      // Finally, rotate everything
      for (let i = 0; i < this.graph.vertices.length; i++) {
        if (i === b) {
          continue
        }

        this.graph.vertices[i].position.rotateAround(angle, this.graph.vertices[b].position)
      }

      for (let i = 0; i < this.rings.length; i++) {
        this.rings[i].center.rotateAround(angle, this.graph.vertices[b].position)
      }
    }
  }

  /**
   * Initializes rings and ringbonds for the current molecule.
   */
  initRings() {
    const openBonds = new Map()

    // Close the open ring bonds (spanning tree -> graph)
    for (let i = this.graph.vertices.length - 1; i >= 0; i--) {
      const vertex = this.graph.vertices[i]

      if (vertex.value.ringbonds.length === 0) {
        continue
      }

      for (let j = 0; j < vertex.value.ringbonds.length; j++) {
        const ringbondId = vertex.value.ringbonds[j].id
        const ringbondBond = vertex.value.ringbonds[j].bond

        // If the other ringbond id has not been discovered,
        // add it to the open bonds map and continue.
        // if the other ringbond id has already been discovered,
        // create a bond between the two atoms.
        if (!openBonds.has(ringbondId)) {
          openBonds.set(ringbondId, [vertex.id, ringbondBond])
        } else {
          const sourceVertexId = vertex.id
          const targetVertexId = openBonds.get(ringbondId)[0]
          const targetRingbondBond = openBonds.get(ringbondId)[1]
          const edge = new Edge(sourceVertexId, targetVertexId, 1)
          edge.setBondType(targetRingbondBond || ringbondBond || '-')
          const edgeId = this.graph.addEdge(edge)
          const targetVertex = this.graph.vertices[targetVertexId]

          vertex.addRingbondChild(targetVertexId, j)
          vertex.value.addNeighbouringElement(targetVertex.value.element)
          targetVertex.addRingbondChild(sourceVertexId, j)
          targetVertex.value.addNeighbouringElement(vertex.value.element)
          vertex.edges.push(edgeId)
          targetVertex.edges.push(edgeId)

          openBonds.delete(ringbondId)
        }
      }
    }

    // Get the rings in the graph (the SSSR)
    const rings = SSSR.getRings(this.graph, this.opts.experimentalSSSR)

    if (rings === null) {
      return
    }

    for (let i = 0; i < rings.length; i++) {
      const ringVertices = [...rings[i]]
      const ringId = this.addRing(new Ring(ringVertices))

      // Add the ring to the atoms
      for (let j = 0; j < ringVertices.length; j++) {
        this.graph.vertices[ringVertices[j]].value.rings.push(ringId)
      }
    }

    // Find connection between rings
    // Check for common vertices and create ring connections. This is a bit
    // ugly, but the ringcount is always fairly low (< 100)
    for (let i = 0; i < this.rings.length - 1; i++) {
      for (let j = i + 1; j < this.rings.length; j++) {
        const a = this.rings[i]
        const b = this.rings[j]
        const ringConnection = new RingConnection(a, b)

        // If there are no vertices in the ring connection, then there
        // is no ring connection
        if (ringConnection.vertices.size > 0) {
          this.addRingConnection(ringConnection)
        }
      }
    }

    // Add neighbours to the rings
    for (let i = 0; i < this.rings.length; i++) {
      const ring = this.rings[i]
      ring.neighbours = RingConnection.getNeighbours(this.ringConnections, ring.id)
    }

    // Anchor the ring to one of it's members, so that the ring center will always
    // be tied to a single vertex when doing repositionings
    for (let i = 0; i < this.rings.length; i++) {
      const ring = this.rings[i]
      this.graph.vertices[ring.members[0]].value.addAnchoredRing(ring.id)
    }

    // Backup the ring information to restore after placing the bridged ring.
    // This is needed in order to identify aromatic rings and stuff like this in
    // rings that are member of the superring.
    this.backupRingInformation()

    // Replace rings contained by a larger bridged ring with a bridged ring
    while (this.rings.length > 0) {
      let id = -1
      for (let i = 0; i < this.rings.length; i++) {
        const ring = this.rings[i]

        if (this.isPartOfBridgedRing(ring.id) && !ring.isBridged) {
          id = ring.id
        }
      }

      if (id === -1) {
        break
      }

      const ring = this.getRing(id)

      const involvedRings = this.getBridgedRingRings(ring.id)

      this.bridgedRing = true
      this.createBridgedRing(involvedRings, ring.members[0])

      // Remove the rings
      for (let i = 0; i < involvedRings.length; i++) {
        this.removeRing(involvedRings[i])
      }
    }
  }

  initHydrogens() {
    // Do not draw hydrogens except when they are connected to a stereocenter connected to two or more rings.
    if (!this.opts.explicitHydrogens) {
      for (let i = 0; i < this.graph.vertices.length; i++) {
        const vertex = this.graph.vertices[i]

        if (vertex.value.element !== 'H') {
          continue
        }

        // Hydrogens should have only one neighbour, so just take the first
        // Also set hasHydrogen true on connected atom
        const neighbour = this.graph.vertices[vertex.neighbours[0]]
        neighbour.value.hasHydrogen = true

        if (!neighbour.value.isStereoCenter || neighbour.value.rings.length < 2 && !neighbour.value.bridgedRing ||
          neighbour.value.bridgedRing && neighbour.value.originalRings.length < 2) {
          vertex.value.isDrawn = false
        }
      }
    }
  }

  /**
   * Returns all rings connected by bridged bonds starting from the ring with the supplied ring id.
   *
   * @param {Number} ringId A ring id.
   * @returns {Number[]} An array containing all ring ids of rings part of a bridged ring system.
   */
  getBridgedRingRings(ringId) {
    const involvedRings = []
    const that = this

    const recurse = function(r) {
      const ring = that.getRing(r)

      involvedRings.push(r)

      for (let i = 0; i < ring.neighbours.length; i++) {
        const n = ring.neighbours[i]

        if (involvedRings.indexOf(n) === -1 &&
          n !== r &&
          RingConnection.isBridge(that.ringConnections, that.graph.vertices, r, n)) {
          recurse(n)
        }
      }
    }

    recurse(ringId)

    return ArrayHelper.unique(involvedRings)
  }

  /**
   * Checks whether or not a ring is part of a bridged ring.
   *
   * @param {Number} ringId A ring id.
   * @returns {Boolean} A boolean indicating whether or not the supplied ring (by id) is part of a bridged ring system.
   */
  isPartOfBridgedRing(ringId) {
    for (let i = 0; i < this.ringConnections.length; i++) {
      if (this.ringConnections[i].containsRing(ringId) &&
        this.ringConnections[i].isBridge(this.graph.vertices)) {
        return true
      }
    }

    return false
  }

  /**
   * Creates a bridged ring.
   *
   * @param {Number[]} ringIds An array of ids of rings involved in the bridged ring.
   * @param {Number} sourceVertexId The vertex id to start the bridged ring discovery from.
   * @returns {Ring} The bridged ring.
   */
  createBridgedRing(ringIds, sourceVertexId) {
    const ringMembers = new Set()
    const vertices = new Set()
    const neighbours = new Set()

    for (let i = 0; i < ringIds.length; i++) {
      const ring = this.getRing(ringIds[i])
      ring.isPartOfBridged = true

      for (let j = 0; j < ring.members.length; j++) {
        vertices.add(ring.members[j])
      }

      for (let j = 0; j < ring.neighbours.length; j++) {
        const id = ring.neighbours[j]

        if (ringIds.indexOf(id) === -1) {
          neighbours.add(ring.neighbours[j])
        }
      }
    }

    // A vertex is part of the bridged ring if it only belongs to
    // one of the rings (or to another ring
    // which is not part of the bridged ring).
    const leftovers = new Set()

    for (const id of vertices) {
      const vertex = this.graph.vertices[id]
      const intersection = ArrayHelper.intersection(ringIds, vertex.value.rings)

      if (vertex.value.rings.length === 1 || intersection.length === 1) {
        ringMembers.add(vertex.id)
      } else {
        leftovers.add(vertex.id)
      }
    }

    // Vertices can also be part of multiple rings and lay on the bridged ring,
    // however, they have to have at least two neighbours that are not part of
    // two rings
    const insideRing = []

    for (const id of leftovers) {
      const vertex = this.graph.vertices[id]
      let onRing = false

      for (let j = 0; j < vertex.edges.length; j++) {
        if (this.edgeRingCount(vertex.edges[j]) === 1) {
          onRing = true
        }
      }

      if (onRing) {
        vertex.value.isBridgeNode = true
        ringMembers.add(vertex.id)
      } else {
        vertex.value.isBridge = true
        ringMembers.add(vertex.id)
      }
    }

    // Create the ring
    const ring = new Ring([...ringMembers])
    this.addRing(ring)

    ring.isBridged = true
    ring.neighbours = [...neighbours]

    for (let i = 0; i < ringIds.length; i++) {
      ring.rings.push(this.getRing(ringIds[i]).clone())
    }

    for (let i = 0; i < ring.members.length; i++) {
      this.graph.vertices[ring.members[i]].value.bridgedRing = ring.id
    }

    // Atoms inside the ring are no longer part of a ring but are now
    // associated with the bridged ring
    for (let i = 0; i < insideRing.length; i++) {
      const vertex = this.graph.vertices[insideRing[i]]
      vertex.value.rings = []
    }

    // Remove former rings from members of the bridged ring and add the bridged ring
    for (const id of ringMembers) {
      const vertex = this.graph.vertices[id]
      vertex.value.rings = ArrayHelper.removeAll(vertex.value.rings, ringIds)
      vertex.value.rings.push(ring.id)
    }

    // Remove all the ring connections no longer used
    for (let i = 0; i < ringIds.length; i++) {
      for (let j = i + 1; j < ringIds.length; j++) {
        this.removeRingConnectionsBetween(ringIds[i], ringIds[j])
      }
    }

    // Update the ring connections and add this ring to the neighbours neighbours
    for (const id of neighbours) {
      const connections = this.getRingConnections(id, ringIds)

      for (let j = 0; j < connections.length; j++) {
        this.getRingConnection(connections[j]).updateOther(ring.id, id)
      }

      this.getRing(id).neighbours.push(ring.id)
    }

    return ring
  }

  /**
   * Returns an array of ring ids shared by both vertices.
   *
   * @param {Vertex} vertexA A vertex.
   * @param {Vertex} vertexB A vertex.
   * @returns {Number[]} An array of ids of rings shared by the two vertices.
   */
  getCommonRings(vertexA, vertexB) {
    const commonRings = []

    for (let i = 0; i < vertexA.value.rings.length; i++) {
      for (let j = 0; j < vertexB.value.rings.length; j++) {
        if (vertexA.value.rings[i] === vertexB.value.rings[j]) {
          commonRings.push(vertexA.value.rings[i])
        }
      }
    }

    return commonRings
  }

  /**
   * Returns the closest vertex (connected as well as unconnected).
   *
   * @param {Vertex} vertex The vertex of which to find the closest other vertex.
   * @returns {Vertex} The closest vertex.
   */
  getClosestVertex(vertex) {
    let minDist = 99999
    let minVertex = null

    for (let i = 0; i < this.graph.vertices.length; i++) {
      const v = this.graph.vertices[i]

      if (v.id === vertex.id) {
        continue
      }

      const distSq = vertex.position.distanceSq(v.position)

      if (distSq < minDist) {
        minDist = distSq
        minVertex = v
      }
    }

    return minVertex
  }

  /**
   * Add a ring to this representation of a molecule.
   *
   * @param {Ring} ring A new ring.
   * @returns {Number} The ring id of the new ring.
   */
  addRing(ring) {
    ring.id = this.ringIdCounter++
    this.rings.push(ring)

    return ring.id
  }

  /**
   * Removes a ring from the array of rings associated with the current molecule.
   *
   * @param {Number} ringId A ring id.
   */
  removeRing(ringId) {
    this.rings = this.rings.filter(function(item) {
      return item.id !== ringId
    })

    // Also remove ring connections involving this ring
    this.ringConnections = this.ringConnections.filter(function(item) {
      return item.firstRingId !== ringId && item.secondRingId !== ringId
    })

    // Remove the ring as neighbour of other rings
    for (let i = 0; i < this.rings.length; i++) {
      const r = this.rings[i]
      r.neighbours = r.neighbours.filter(function(item) {
        return item !== ringId
      })
    }
  }

  /**
   * Gets a ring object from the array of rings associated with the current molecule by its id. The ring id is not equal to the index, since rings can be added and removed when processing bridged rings.
   *
   * @param {Number} ringId A ring id.
   * @returns {Ring} A ring associated with the current molecule.
   */
  getRing(ringId) {
    for (let i = 0; i < this.rings.length; i++) {
      if (this.rings[i].id === ringId) {
        return this.rings[i]
      }
    }
  }

  /**
   * Add a ring connection to this representation of a molecule.
   *
   * @param {RingConnection} ringConnection A new ringConnection.
   * @returns {Number} The ring connection id of the new ring connection.
   */
  addRingConnection(ringConnection) {
    ringConnection.id = this.ringConnectionIdCounter++
    this.ringConnections.push(ringConnection)

    return ringConnection.id
  }

  /**
   * Removes a ring connection from the array of rings connections associated with the current molecule.
   *
   * @param {Number} ringConnectionId A ring connection id.
   */
  removeRingConnection(ringConnectionId) {
    this.ringConnections = this.ringConnections.filter(function(item) {
      return item.id !== ringConnectionId
    })
  }

  /**
   * Removes all ring connections between two vertices.
   *
   * @param {Number} vertexIdA A vertex id.
   * @param {Number} vertexIdB A vertex id.
   */
  removeRingConnectionsBetween(vertexIdA, vertexIdB) {
    const toRemove = []
    for (let i = 0; i < this.ringConnections.length; i++) {
      const ringConnection = this.ringConnections[i]

      if (ringConnection.firstRingId === vertexIdA && ringConnection.secondRingId === vertexIdB ||
        ringConnection.firstRingId === vertexIdB && ringConnection.secondRingId === vertexIdA) {
        toRemove.push(ringConnection.id)
      }
    }

    for (let i = 0; i < toRemove.length; i++) {
      this.removeRingConnection(toRemove[i])
    }
  }

  /**
   * Get a ring connection with a given id.
   *
   * @param {Number} id
   * @returns {RingConnection} The ring connection with the specified id.
   */
  getRingConnection(id) {
    for (let i = 0; i < this.ringConnections.length; i++) {
      if (this.ringConnections[i].id === id) {
        return this.ringConnections[i]
      }
    }
  }

  /**
   * Get the ring connections between a ring and a set of rings.
   *
   * @param {Number} ringId A ring id.
   * @param {Number[]} ringIds An array of ring ids.
   * @returns {Number[]} An array of ring connection ids.
   */
  getRingConnections(ringId, ringIds) {
    const ringConnections = []

    for (let i = 0; i < this.ringConnections.length; i++) {
      const rc = this.ringConnections[i]

      for (let j = 0; j < ringIds.length; j++) {
        const id = ringIds[j]

        if (rc.firstRingId === ringId && rc.secondRingId === id ||
          rc.firstRingId === id && rc.secondRingId === ringId) {
          ringConnections.push(rc.id)
        }
      }
    }

    return ringConnections
  }

  /**
   * Returns the overlap score of the current molecule based on its positioned vertices. The higher the score, the more overlaps occur in the structure drawing.
   *
   * @returns {Object} Returns the total overlap score and the overlap score of each vertex sorted by score (higher to lower). Example: { total: 99, scores: [ { id: 0, score: 22 }, ... ]  }
   */
  getOverlapScore() {
    let total = 0.0
    const overlapScores = new Float32Array(this.graph.vertices.length)

    for (let i = 0; i < this.graph.vertices.length; i++) {
      overlapScores[i] = 0
    }

    for (let i = 0; i < this.graph.vertices.length; i++) {
      let j = this.graph.vertices.length
      while (--j > i) {
        const a = this.graph.vertices[i]
        const b = this.graph.vertices[j]

        if (!a.value.isDrawn || !b.value.isDrawn) {
          continue
        }

        const dist = Vector2.subtract(a.position, b.position).lengthSq()

        if (dist < this.opts.bondLengthSq) {
          const weighted = (this.opts.bondLength - Math.sqrt(dist)) / this.opts.bondLength
          total += weighted
          overlapScores[i] += weighted
          overlapScores[j] += weighted
        }
      }
    }

    const sortable = []

    for (let i = 0; i < this.graph.vertices.length; i++) {
      sortable.push({
        id: i,
        score: overlapScores[i]
      })
    }

    sortable.sort(function(a, b) {
      return b.score - a.score
    })

    return {
      total: total,
      scores: sortable,
      vertexScores: overlapScores
    }
  }

  /**
   * Sets the center for a ring.
   *
   * @param {Ring} ring A ring.
   */
  setRingCenter(ring) {
    const ringSize = ring.getSize()
    const total = new Vector2(0, 0)

    for (let i = 0; i < ringSize; i++) {
      total.add(this.graph.vertices[ring.members[i]].position)
    }

    ring.center = total.divide(ringSize)
  }

  /**
   * Position the vertices according to their bonds and properties.
   */
  position() {
    let startVertex = null

    // Always start drawing at a bridged ring if there is one
    // If not, start with a ring
    // else, start with 0
    for (let i = 0; i < this.graph.vertices.length; i++) {
      if (this.graph.vertices[i].value.bridgedRing !== null) {
        startVertex = this.graph.vertices[i]
        break
      }
    }

    for (let i = 0; i < this.rings.length; i++) {
      if (this.rings[i].isBridged) {
        startVertex = this.graph.vertices[this.rings[i].members[0]]
      }
    }

    if (this.rings.length > 0 && startVertex === null) {
      startVertex = this.graph.vertices[this.rings[0].members[0]]
    }

    if (startVertex === null) {
      startVertex = this.graph.vertices[0]
    }

    this.createNextBond(startVertex, null, 0.0)
  }

  /**
   * Stores the current information associated with rings.
   */
  backupRingInformation() {
    this.originalRings = []
    this.originalRingConnections = []

    for (let i = 0; i < this.rings.length; i++) {
      this.originalRings.push(this.rings[i])
    }

    for (let i = 0; i < this.ringConnections.length; i++) {
      this.originalRingConnections.push(this.ringConnections[i])
    }

    for (let i = 0; i < this.graph.vertices.length; i++) {
      this.graph.vertices[i].value.backupRings()
    }
  }

  /**
   * Restores the most recently backed up information associated with rings.
   */
  restoreRingInformation() {
    // Get the subring centers from the bridged rings
    const bridgedRings = this.getBridgedRings()

    this.rings = []
    this.ringConnections = []

    for (let i = 0; i < bridgedRings.length; i++) {
      const bridgedRing = bridgedRings[i]

      for (let j = 0; j < bridgedRing.rings.length; j++) {
        const ring = bridgedRing.rings[j]
        this.originalRings[ring.id].center = ring.center
      }
    }

    for (let i = 0; i < this.originalRings.length; i++) {
      this.rings.push(this.originalRings[i])
    }

    for (let i = 0; i < this.originalRingConnections.length; i++) {
      this.ringConnections.push(this.originalRingConnections[i])
    }

    for (let i = 0; i < this.graph.vertices.length; i++) {
      this.graph.vertices[i].value.restoreRings()
    }
  }

  /**
   * Creates a new ring, that is, positiones all the vertices inside a ring.
   *
   * @param {Ring} ring The ring to position.
   * @param {(Vector2|null)} [center=null] The center of the ring to be created.
   * @param {(Vertex|null)} [startVertex=null] The first vertex to be positioned inside the ring.
   * @param {(Vertex|null)} [previousVertex=null] The last vertex that was positioned.
   * @param {Boolean} [previousVertex=false] A boolean indicating whether or not this ring was force positioned already - this is needed after force layouting a ring, in order to draw rings connected to it.
   */
  createRing(ring, center = null, startVertex = null, previousVertex = null) {
    if (ring.positioned) {
      return
    }

    center = center || new Vector2(0, 0)

    const orderedNeighbours = ring.getOrderedNeighbours(this.ringConnections)
    const startingAngle = startVertex ? Vector2.subtract(startVertex.position, center).angle() : 0

    const radius = MathHelper.polyCircumradius(this.opts.bondLength, ring.getSize())
    const angle = MathHelper.centralAngle(ring.getSize())

    ring.centralAngle = angle

    let a = startingAngle
    const that = this
    let startVertexId = (startVertex) ? startVertex.id : null

    if (ring.members.indexOf(startVertexId) === -1) {
      if (startVertex) {
        startVertex.positioned = false
      }

      startVertexId = ring.members[0]
    }

    // If the ring is bridged, then draw the vertices inside the ring
    // using a force based approach
    if (ring.isBridged) {
      this.graph.kkLayout(ring.members.slice(), center, startVertex.id, ring, this.opts.bondLength,
        this.opts.kkThreshold, this.opts.kkInnerThreshold, this.opts.kkMaxIteration,
        this.opts.kkMaxInnerIteration, this.opts.kkMaxEnergy)
      ring.positioned = true

      // Update the center of the bridged ring
      this.setRingCenter(ring)
      center = ring.center

      // Setting the centers for the subrings
      for (let i = 0; i < ring.rings.length; i++) {
        this.setRingCenter(ring.rings[i])
      }
    } else {
      ring.eachMember(this.graph.vertices, function(v) {
        const vertex = that.graph.vertices[v]

        if (!vertex.positioned) {
          vertex.setPosition(center.x + Math.cos(a) * radius, center.y + Math.sin(a) * radius)
        }

        a += angle

        if (!ring.isBridged || ring.rings.length < 3) {
          vertex.angle = a
          vertex.positioned = true
        }
      }, startVertexId, (previousVertex) ? previousVertex.id : null)
    }

    ring.positioned = true
    ring.center = center

    // Draw neighbours in decreasing order of connectivity
    for (let i = 0; i < orderedNeighbours.length; i++) {
      const neighbour = this.getRing(orderedNeighbours[i].neighbour)

      if (neighbour.positioned) {
        continue
      }

      const vertices = RingConnection.getVertices(this.ringConnections, ring.id, neighbour.id)

      if (vertices.length === 2) {
        // This ring is a fused ring
        ring.isFused = true
        neighbour.isFused = true

        const vertexA = this.graph.vertices[vertices[0]]
        const vertexB = this.graph.vertices[vertices[1]]

        // Get middle between vertex A and B
        const midpoint = Vector2.midpoint(vertexA.position, vertexB.position)

        // Get the normals to the line between A and B
        const normals = Vector2.normals(vertexA.position, vertexB.position)

        // Normalize the normals
        normals[0].normalize()
        normals[1].normalize()

        // Set length from middle of side to center (the apothem)
        const r = MathHelper.polyCircumradius(this.opts.bondLength, neighbour.getSize())
        const apothem = MathHelper.apothem(r, neighbour.getSize())

        normals[0].multiplyScalar(apothem).add(midpoint)
        normals[1].multiplyScalar(apothem).add(midpoint)

        // Pick the normal which results in a larger distance to the previous center
        // Also check whether it's inside another ring
        let nextCenter = normals[0]
        if (Vector2.subtract(center, normals[1]).lengthSq() > Vector2.subtract(center, normals[0]).lengthSq()) {
          nextCenter = normals[1]
        }

        // Get the vertex (A or B) which is in clock-wise direction of the other
        const posA = Vector2.subtract(vertexA.position, nextCenter)
        const posB = Vector2.subtract(vertexB.position, nextCenter)

        if (posA.clockwise(posB) === -1) {
          if (!neighbour.positioned) {
            this.createRing(neighbour, nextCenter, vertexA, vertexB)
          }
        } else {
          if (!neighbour.positioned) {
            this.createRing(neighbour, nextCenter, vertexB, vertexA)
          }
        }
      } else if (vertices.length === 1) {
        // This ring is a spiro
        ring.isSpiro = true
        neighbour.isSpiro = true

        const vertexA = this.graph.vertices[vertices[0]]

        // Get the vector pointing from the shared vertex to the new centpositioner
        const nextCenter = Vector2.subtract(center, vertexA.position)

        nextCenter.invert()
        nextCenter.normalize()

        // Get the distance from the vertex to the center
        const r = MathHelper.polyCircumradius(this.opts.bondLength, neighbour.getSize())

        nextCenter.multiplyScalar(r)
        nextCenter.add(vertexA.position)

        if (!neighbour.positioned) {
          this.createRing(neighbour, nextCenter, vertexA)
        }
      }
    }

    // Next, draw atoms that are not part of a ring that are directly attached to this ring
    for (let i = 0; i < ring.members.length; i++) {
      const ringMember = this.graph.vertices[ring.members[i]]
      const ringMemberNeighbours = ringMember.neighbours

      // If there are multiple, the ovlerap will be resolved in the appropriate step
      for (let j = 0; j < ringMemberNeighbours.length; j++) {
        const v = this.graph.vertices[ringMemberNeighbours[j]]

        if (v.positioned) {
          continue
        }

        v.value.isConnectedToRing = true
        this.createNextBond(v, ringMember, 0.0)
      }
    }
  }

  /**
   * Rotate an entire subtree by an angle around a center.
   *
   * @param {Number} vertexId A vertex id (the root of the sub-tree).
   * @param {Number} parentVertexId A vertex id in the previous direction of the subtree that is to rotate.
   * @param {Number} angle An angle in randians.
   * @param {Vector2} center The rotational center.
   */
  rotateSubtree(vertexId, parentVertexId, angle, center) {
    const that = this

    this.graph.traverseTree(vertexId, parentVertexId, function(vertex) {
      vertex.position.rotateAround(angle, center)

      for (let i = 0; i < vertex.value.anchoredRings.length; i++) {
        const ring = that.rings[vertex.value.anchoredRings[i]]

        if (ring) {
          ring.center.rotateAround(angle, center)
        }
      }
    })
  }

  /**
   * Gets the overlap score of a subtree.
   *
   * @param {Number} vertexId A vertex id (the root of the sub-tree).
   * @param {Number} parentVertexId A vertex id in the previous direction of the subtree.
   * @param {Number[]} vertexOverlapScores An array containing the vertex overlap scores indexed by vertex id.
   * @returns {Object} An object containing the total overlap score and the center of mass of the subtree weighted by overlap score { value: 0.2, center: new Vector2() }.
   */
  getSubtreeOverlapScore(vertexId, parentVertexId, vertexOverlapScores) {
    const that = this
    let score = 0
    const center = new Vector2(0, 0)
    let count = 0

    this.graph.traverseTree(vertexId, parentVertexId, function(vertex) {
      if (!vertex.value.isDrawn) {
        return
      }

      const s = vertexOverlapScores[vertex.id]
      if (s > that.opts.overlapSensitivity) {
        score += s
        count++
      }

      const position = that.graph.vertices[vertex.id].position.clone()
      position.multiplyScalar(s)
      center.add(position)
    })

    center.divide(score)

    return {
      value: score / count,
      center: center
    }
  }

  /**
   * Returns the current (positioned vertices so far) center of mass.
   *
   * @returns {Vector2} The current center of mass.
   */
  getCurrentCenterOfMass() {
    const total = new Vector2(0, 0)
    let count = 0

    for (let i = 0; i < this.graph.vertices.length; i++) {
      const vertex = this.graph.vertices[i]

      if (vertex.positioned) {
        total.add(vertex.position)
        count++
      }
    }

    return total.divide(count)
  }

  /**
   * Resolve primary (exact) overlaps, such as two vertices that are connected to the same ring vertex.
   */
  resolvePrimaryOverlaps() {
    const overlaps = []
    const done = Array(this.graph.vertices.length)

    // Looking for overlaps created by two bonds coming out of a ring atom, which both point straight
    // away from the ring and are thus perfectly overlapping.
    for (let i = 0; i < this.rings.length; i++) {
      const ring = this.rings[i]

      for (let j = 0; j < ring.members.length; j++) {
        const vertex = this.graph.vertices[ring.members[j]]

        if (done[vertex.id]) {
          continue
        }

        done[vertex.id] = true

        const nonRingNeighbours = this.getNonRingNeighbours(vertex.id)

        if (nonRingNeighbours.length > 1) {
          // Look for rings where there are atoms with two bonds outside the ring (overlaps)
          const rings = []

          for (let k = 0; k < vertex.value.rings.length; k++) {
            rings.push(vertex.value.rings[k])
          }

          overlaps.push({
            common: vertex,
            rings: rings,
            vertices: nonRingNeighbours
          })
        } else if (nonRingNeighbours.length === 1 && vertex.value.rings.length === 2) {
          // Look for bonds coming out of joined rings to adjust the angle, an example is: C1=CC(=CC=C1)[C@]12SCCN1CC1=CC=CC=C21
          // where the angle has to be adjusted to account for fused ring
          const rings = []

          for (let k = 0; k < vertex.value.rings.length; k++) {
            rings.push(vertex.value.rings[k])
          }

          overlaps.push({
            common: vertex,
            rings: rings,
            vertices: nonRingNeighbours
          })
        }
      }
    }

    for (let i = 0; i < overlaps.length; i++) {
      const overlap = overlaps[i]

      if (overlap.vertices.length === 2) {
        const a = overlap.vertices[0]
        const b = overlap.vertices[1]

        if (!a.value.isDrawn || !b.value.isDrawn) {
          continue
        }

        const angle = (2 * Math.PI - this.getRing(overlap.rings[0]).getAngle()) / 6.0

        this.rotateSubtree(a.id, overlap.common.id, angle, overlap.common.position)
        this.rotateSubtree(b.id, overlap.common.id, -angle, overlap.common.position)

        // Decide which way to rotate the vertices depending on the effect it has on the overlap score
        let overlapScore = this.getOverlapScore()
        let subTreeOverlapA = this.getSubtreeOverlapScore(a.id, overlap.common.id, overlapScore.vertexScores)
        let subTreeOverlapB = this.getSubtreeOverlapScore(b.id, overlap.common.id, overlapScore.vertexScores)
        const total = subTreeOverlapA.value + subTreeOverlapB.value

        this.rotateSubtree(a.id, overlap.common.id, -2.0 * angle, overlap.common.position)
        this.rotateSubtree(b.id, overlap.common.id, 2.0 * angle, overlap.common.position)

        overlapScore = this.getOverlapScore()
        subTreeOverlapA = this.getSubtreeOverlapScore(a.id, overlap.common.id, overlapScore.vertexScores)
        subTreeOverlapB = this.getSubtreeOverlapScore(b.id, overlap.common.id, overlapScore.vertexScores)

        if (subTreeOverlapA.value + subTreeOverlapB.value > total) {
          this.rotateSubtree(a.id, overlap.common.id, 2.0 * angle, overlap.common.position)
          this.rotateSubtree(b.id, overlap.common.id, -2.0 * angle, overlap.common.position)
        }
      } else if (overlap.vertices.length === 1) {
        if (overlap.rings.length === 2) {
          // TODO: Implement for more overlap resolution
          // console.log(overlap);
        }
      }
    }
  }

  /**
   * Resolve secondary overlaps. Those overlaps are due to the structure turning back on itself.
   *
   * @param {Object[]} scores An array of objects sorted descending by score.
   * @param {Number} scores[].id A vertex id.
   * @param {Number} scores[].score The overlap score associated with the vertex id.
   */
  resolveSecondaryOverlaps(scores) {
    for (let i = 0; i < scores.length; i++) {
      if (scores[i].score > this.opts.overlapSensitivity) {
        const vertex = this.graph.vertices[scores[i].id]

        if (vertex.isTerminal()) {
          const closest = this.getClosestVertex(vertex)

          if (closest) {
            // If one of the vertices is the first one, the previous vertex is not the central vertex but the dummy
            // so take the next rather than the previous, which is vertex 1
            let closestPosition = null

            if (closest.isTerminal()) {
              closestPosition = closest.id === 0 ? this.graph.vertices[1].position : closest.previousPosition
            } else {
              closestPosition = closest.id === 0 ? this.graph.vertices[1].position : closest.position
            }

            const vertexPreviousPosition = vertex.id === 0 ? this.graph.vertices[1].position : vertex.previousPosition

            vertex.position.rotateAwayFrom(closestPosition, vertexPreviousPosition, MathHelper.toRad(20))
          }
        }
      }
    }
  }

  /**
   * Get the last non-null or 0 angle vertex.
   * @param {Number} vertexId A vertex id.
   * @returns {Vertex} The last vertex with an angle that was not 0 or null.
   */
  getLastVertexWithAngle(vertexId) {
    let angle = 0
    let vertex = null

    while (!angle && vertexId) {
      vertex = this.graph.vertices[vertexId]
      angle = vertex.angle
      vertexId = vertex.parentVertexId
    }

    return vertex
  }

  /**
   * Positiones the next vertex thus creating a bond.
   *
   * @param {Vertex} vertex A vertex.
   * @param {Vertex} [previousVertex=null] The previous vertex which has been positioned.
   * @param {Number} [angle=0.0] The (global) angle of the vertex.
   * @param {Boolean} [originShortest=false] Whether the origin is the shortest subtree in the branch.
   * @param {Boolean} [skipPositioning=false] Whether or not to skip positioning and just check the neighbours.
   */
  createNextBond(vertex, previousVertex = null, angle = 0.0, originShortest = false, skipPositioning = false) {
    if (vertex.positioned && !skipPositioning) {
      return
    }

    // If the double bond config was set on this vertex, do not check later
    let doubleBondConfigSet = false

    // Keeping track of configurations around double bonds
    if (previousVertex) {
      const edge = this.graph.getEdge(vertex.id, previousVertex.id)

      if ((edge.bondType === '/' || edge.bondType === '\\') && ++this.doubleBondConfigCount % 2 === 1) {
        if (this.doubleBondConfig === null) {
          this.doubleBondConfig = edge.bondType
          doubleBondConfigSet = true

          // Switch if the bond is a branch bond and previous vertex is the first
          // TODO: Why is it different with the first vertex?
          if (previousVertex.parentVertexId === null && vertex.value.branchBond) {
            if (this.doubleBondConfig === '/') {
              this.doubleBondConfig = '\\'
            } else if (this.doubleBondConfig === '\\') {
              this.doubleBondConfig = '/'
            }
          }
        }
      }
    }

    // If the current node is the member of one ring, then point straight away
    // from the center of the ring. However, if the current node is a member of
    // two rings, point away from the middle of the centers of the two rings
    if (!skipPositioning) {
      if (!previousVertex) {
        // Add a (dummy) previous position if there is no previous vertex defined
        // Since the first vertex is at (0, 0), create a vector at (bondLength, 0)
        // and rotate it by 90°

        const dummy = new Vector2(this.opts.bondLength, 0)
        dummy.rotate(MathHelper.toRad(-60))

        vertex.previousPosition = dummy
        vertex.setPosition(this.opts.bondLength, 0)
        vertex.angle = MathHelper.toRad(-60)

        // Do not position the vertex if it belongs to a bridged ring that is positioned using a layout algorithm.
        if (vertex.value.bridgedRing === null) {
          vertex.positioned = true
        }
      } else if (previousVertex.value.rings.length > 0) {
        const neighbours = previousVertex.neighbours
        let joinedVertex = null
        let pos = new Vector2(0.0, 0.0)

        if (previousVertex.value.bridgedRing === null && previousVertex.value.rings.length > 1) {
          for (let i = 0; i < neighbours.length; i++) {
            const neighbour = this.graph.vertices[neighbours[i]]
            if (ArrayHelper.containsAll(neighbour.value.rings, previousVertex.value.rings)) {
              joinedVertex = neighbour
              break
            }
          }
        }

        if (joinedVertex === null) {
          for (let i = 0; i < neighbours.length; i++) {
            const v = this.graph.vertices[neighbours[i]]

            if (v.positioned && this.areVerticesInSameRing(v, previousVertex)) {
              pos.add(Vector2.subtract(v.position, previousVertex.position))
            }
          }

          pos.invert().normalize().multiplyScalar(this.opts.bondLength).add(previousVertex.position)
        } else {
          pos = joinedVertex.position.clone().rotateAround(Math.PI, previousVertex.position)
        }

        vertex.previousPosition = previousVertex.position
        vertex.setPositionFromVector(pos)
        vertex.positioned = true
      } else {
        // If the previous vertex was not part of a ring, draw a bond based
        // on the global angle of the previous bond
        const v = new Vector2(this.opts.bondLength, 0)

        v.rotate(angle)
        v.add(previousVertex.position)

        vertex.setPositionFromVector(v)
        vertex.previousPosition = previousVertex.position
        vertex.positioned = true
      }
    }

    // Go to next vertex
    // If two rings are connected by a bond ...
    if (vertex.value.bridgedRing !== null) {
      const nextRing = this.getRing(vertex.value.bridgedRing)

      if (!nextRing.positioned) {
        const nextCenter = Vector2.subtract(vertex.previousPosition, vertex.position)

        nextCenter.invert()
        nextCenter.normalize()

        const r = MathHelper.polyCircumradius(this.opts.bondLength, nextRing.members.length)
        nextCenter.multiplyScalar(r)
        nextCenter.add(vertex.position)

        this.createRing(nextRing, nextCenter, vertex)
      }
    } else if (vertex.value.rings.length > 0) {
      const nextRing = this.getRing(vertex.value.rings[0])

      if (!nextRing.positioned) {
        const nextCenter = Vector2.subtract(vertex.previousPosition, vertex.position)

        nextCenter.invert()
        nextCenter.normalize()

        const r = MathHelper.polyCircumradius(this.opts.bondLength, nextRing.getSize())

        nextCenter.multiplyScalar(r)
        nextCenter.add(vertex.position)

        this.createRing(nextRing, nextCenter, vertex)
      }
    } else {
      // Draw the non-ring vertices connected to this one
      const tmpNeighbours = vertex.getNeighbours()
      let neighbours = []

      // Remove neighbours that are not drawn
      for (let i = 0; i < tmpNeighbours.length; i++) {
        if (this.graph.vertices[tmpNeighbours[i]].value.isDrawn) {
          neighbours.push(tmpNeighbours[i])
        }
      }

      // Remove the previous vertex (which has already been drawn)
      if (previousVertex) {
        neighbours = ArrayHelper.remove(neighbours, previousVertex.id)
      }

      const previousAngle = vertex.getAngle()

      if (neighbours.length === 1) {
        const nextVertex = this.graph.vertices[neighbours[0]]

        // Make a single chain always cis except when there's a tribble (yes, this is a Star Trek reference) bond
        // or if there are successive double bonds. Added a ring check because if there is an aromatic ring the ring bond inside the ring counts as a double bond and leads to =-= being straight.
        if ((vertex.value.bondType === '#' || (previousVertex && previousVertex.value.bondType === '#')) ||
          vertex.value.bondType === '=' && previousVertex && previousVertex.value.rings.length === 0 &&
          previousVertex.value.bondType === '=' && vertex.value.branchBond !== '-') {
          vertex.value.drawExplicit = false

          if (previousVertex) {
            const straightEdge1 = this.graph.getEdge(vertex.id, previousVertex.id)
            straightEdge1.center = true
          }

          const straightEdge2 = this.graph.getEdge(vertex.id, nextVertex.id)
          straightEdge2.center = true

          if (vertex.value.bondType === '#' || previousVertex && previousVertex.value.bondType === '#') {
            nextVertex.angle = 0.0
          }

          nextVertex.drawExplicit = true

          this.createNextBond(nextVertex, vertex, previousAngle + nextVertex.angle)
        } else if (previousVertex && previousVertex.value.rings.length > 0) {
          // If coming out of a ring, always draw away from the center of mass
          const proposedAngleA = MathHelper.toRad(60)
          const proposedAngleB = -proposedAngleA

          const proposedVectorA = new Vector2(this.opts.bondLength, 0)
          const proposedVectorB = new Vector2(this.opts.bondLength, 0)

          proposedVectorA.rotate(proposedAngleA).add(vertex.position)
          proposedVectorB.rotate(proposedAngleB).add(vertex.position)

          // let centerOfMass = this.getCurrentCenterOfMassInNeigbourhood(vertex.position, 100);
          const centerOfMass = this.getCurrentCenterOfMass()
          const distanceA = proposedVectorA.distanceSq(centerOfMass)
          const distanceB = proposedVectorB.distanceSq(centerOfMass)

          nextVertex.angle = distanceA < distanceB ? proposedAngleB : proposedAngleA

          this.createNextBond(nextVertex, vertex, previousAngle + nextVertex.angle)
        } else {
          let a = vertex.angle
          // Take the min and max if the previous angle was in a 4-neighbourhood (90° angles)
          // TODO: If a is null or zero, it should be checked whether or not this one should go cis or trans, that is,
          //       it should go into the oposite direction of the last non-null or 0 previous vertex / angle.
          if (previousVertex && previousVertex.neighbours.length > 3) {
            if (a > 0) {
              a = Math.min(1.0472, a)
            } else if (a < 0) {
              a = Math.max(-1.0472, a)
            } else {
              a = 1.0472
            }
          } else if (!a) {
            const v = this.getLastVertexWithAngle(vertex.id)
            a = v.angle

            if (!a) {
              a = 1.0472
            }
          }

          // Handle configuration around double bonds
          if (previousVertex && !doubleBondConfigSet) {
            const bondType = this.graph.getEdge(vertex.id, nextVertex.id).bondType

            if (bondType === '/') {
              if (this.doubleBondConfig === '/') {
                // Nothing to do since it will be trans per default
              } else if (this.doubleBondConfig === '\\') {
                a = -a
              }
              this.doubleBondConfig = null
            } else if (bondType === '\\') {
              if (this.doubleBondConfig === '/') {
                a = -a
              } else if (this.doubleBondConfig === '\\') {
                // Nothing to do since it will be trans per default
              }
              this.doubleBondConfig = null
            }
          }

          if (originShortest) {
            nextVertex.angle = a
          } else {
            nextVertex.angle = -a
          }

          this.createNextBond(nextVertex, vertex, previousAngle + nextVertex.angle)
        }
      } else if (neighbours.length === 2) {
        // If the previous vertex comes out of a ring, it doesn't have an angle set
        let a = vertex.angle

        if (!a) {
          a = 1.0472
        }

        // Check for the longer subtree - always go with cis for the longer subtree
        const subTreeDepthA = this.graph.getTreeDepth(neighbours[0], vertex.id)
        const subTreeDepthB = this.graph.getTreeDepth(neighbours[1], vertex.id)

        const l = this.graph.vertices[neighbours[0]]
        const r = this.graph.vertices[neighbours[1]]

        l.value.subtreeDepth = subTreeDepthA
        r.value.subtreeDepth = subTreeDepthB

        // Also get the subtree for the previous direction (this is important when
        // the previous vertex is the shortest path)
        const subTreeDepthC = this.graph.getTreeDepth(previousVertex ? previousVertex.id : null, vertex.id)
        if (previousVertex) {
          previousVertex.value.subtreeDepth = subTreeDepthC
        }

        let cis = 0
        let trans = 1

        // Carbons go always cis
        if (r.value.element === 'C' && l.value.element !== 'C' && subTreeDepthB > 1 && subTreeDepthA < 5) {
          cis = 1
          trans = 0
        } else if (r.value.element !== 'C' && l.value.element === 'C' && subTreeDepthA > 1 && subTreeDepthB < 5) {
          cis = 0
          trans = 1
        } else if (subTreeDepthB > subTreeDepthA) {
          cis = 1
          trans = 0
        }

        const cisVertex = this.graph.vertices[neighbours[cis]]
        const transVertex = this.graph.vertices[neighbours[trans]]

        // If the origin tree is the shortest, make them the main chain
        let originShortest = false
        if (subTreeDepthC < subTreeDepthA && subTreeDepthC < subTreeDepthB) {
          originShortest = true
        }

        transVertex.angle = a
        cisVertex.angle = -a

        if (this.doubleBondConfig === '\\') {
          if (transVertex.value.branchBond === '\\') {
            transVertex.angle = -a
            cisVertex.angle = a
          }
        } else if (this.doubleBondConfig === '/') {
          if (transVertex.value.branchBond === '/') {
            transVertex.angle = -a
            cisVertex.angle = a
          }
        }

        this.createNextBond(transVertex, vertex, previousAngle + transVertex.angle, originShortest)
        this.createNextBond(cisVertex, vertex, previousAngle + cisVertex.angle, originShortest)
      } else if (neighbours.length === 3) {
        // The vertex with the longest sub-tree should always go straight
        const d1 = this.graph.getTreeDepth(neighbours[0], vertex.id)
        const d2 = this.graph.getTreeDepth(neighbours[1], vertex.id)
        const d3 = this.graph.getTreeDepth(neighbours[2], vertex.id)

        let s = this.graph.vertices[neighbours[0]]
        let l = this.graph.vertices[neighbours[1]]
        let r = this.graph.vertices[neighbours[2]]

        s.value.subtreeDepth = d1
        l.value.subtreeDepth = d2
        r.value.subtreeDepth = d3

        if (d2 > d1 && d2 > d3) {
          s = this.graph.vertices[neighbours[1]]
          l = this.graph.vertices[neighbours[0]]
          r = this.graph.vertices[neighbours[2]]
        } else if (d3 > d1 && d3 > d2) {
          s = this.graph.vertices[neighbours[2]]
          l = this.graph.vertices[neighbours[0]]
          r = this.graph.vertices[neighbours[1]]
        }

        // Create a cross if more than one subtree is of length > 1
        // or the vertex is connected to a ring
        if (previousVertex &&
          previousVertex.value.rings.length < 1 &&
          s.value.rings.length < 1 &&
          l.value.rings.length < 1 &&
          r.value.rings.length < 1 &&
          this.graph.getTreeDepth(l.id, vertex.id) === 1 &&
          this.graph.getTreeDepth(r.id, vertex.id) === 1 &&
          this.graph.getTreeDepth(s.id, vertex.id) > 1) {
          s.angle = -vertex.angle
          if (vertex.angle >= 0) {
            l.angle = MathHelper.toRad(30)
            r.angle = MathHelper.toRad(90)
          } else {
            l.angle = -MathHelper.toRad(30)
            r.angle = -MathHelper.toRad(90)
          }

          this.createNextBond(s, vertex, previousAngle + s.angle)
          this.createNextBond(l, vertex, previousAngle + l.angle)
          this.createNextBond(r, vertex, previousAngle + r.angle)
        } else {
          s.angle = 0.0
          l.angle = MathHelper.toRad(90)
          r.angle = -MathHelper.toRad(90)

          this.createNextBond(s, vertex, previousAngle + s.angle)
          this.createNextBond(l, vertex, previousAngle + l.angle)
          this.createNextBond(r, vertex, previousAngle + r.angle)
        }
      } else if (neighbours.length === 4) {
        // The vertex with the longest sub-tree should always go to the reflected opposide direction
        const d1 = this.graph.getTreeDepth(neighbours[0], vertex.id)
        const d2 = this.graph.getTreeDepth(neighbours[1], vertex.id)
        const d3 = this.graph.getTreeDepth(neighbours[2], vertex.id)
        const d4 = this.graph.getTreeDepth(neighbours[3], vertex.id)

        let w = this.graph.vertices[neighbours[0]]
        let x = this.graph.vertices[neighbours[1]]
        let y = this.graph.vertices[neighbours[2]]
        let z = this.graph.vertices[neighbours[3]]

        w.value.subtreeDepth = d1
        x.value.subtreeDepth = d2
        y.value.subtreeDepth = d3
        z.value.subtreeDepth = d4

        if (d2 > d1 && d2 > d3 && d2 > d4) {
          w = this.graph.vertices[neighbours[1]]
          x = this.graph.vertices[neighbours[0]]
          y = this.graph.vertices[neighbours[2]]
          z = this.graph.vertices[neighbours[3]]
        } else if (d3 > d1 && d3 > d2 && d3 > d4) {
          w = this.graph.vertices[neighbours[2]]
          x = this.graph.vertices[neighbours[0]]
          y = this.graph.vertices[neighbours[1]]
          z = this.graph.vertices[neighbours[3]]
        } else if (d4 > d1 && d4 > d2 && d4 > d3) {
          w = this.graph.vertices[neighbours[3]]
          x = this.graph.vertices[neighbours[0]]
          y = this.graph.vertices[neighbours[1]]
          z = this.graph.vertices[neighbours[2]]
        }

        w.angle = -MathHelper.toRad(36)
        x.angle = MathHelper.toRad(36)
        y.angle = -MathHelper.toRad(108)
        z.angle = MathHelper.toRad(108)

        this.createNextBond(w, vertex, previousAngle + w.angle)
        this.createNextBond(x, vertex, previousAngle + x.angle)
        this.createNextBond(y, vertex, previousAngle + y.angle)
        this.createNextBond(z, vertex, previousAngle + z.angle)
      }
    }
  }

  /**
   * Check whether or not an edge is rotatable.
   *
   * @param {Edge} edge An edge.
   * @returns {Boolean} A boolean indicating whether or not the edge is rotatable.
   */
  isEdgeRotatable(edge) {
    const vertexA = this.graph.vertices[edge.sourceId]
    const vertexB = this.graph.vertices[edge.targetId]

    // Only single bonds are rotatable
    if (edge.bondType !== '-') {
      return false
    }

    // Do not rotate edges that have a further single bond to each side - do that!
    // If the bond is terminal, it doesn't make sense to rotate it
    // if (vertexA.getNeighbourCount() + vertexB.getNeighbourCount() < 5) {
    //   return false;
    // }

    if (vertexA.isTerminal() || vertexB.isTerminal()) {
      return false
    }

    // Ringbonds are not rotatable
    // noinspection RedundantIfStatementJS
    if (vertexA.value.rings.length > 0 && vertexB.value.rings.length > 0 &&
      this.areVerticesInSameRing(vertexA, vertexB)) {
      return false
    }

    return true
  }

  /**
   * Returns an array of vertices that are neighbouring a vertix but are not members of a ring (including bridges).
   *
   * @param {Number} vertexId A vertex id.
   * @returns {Vertex[]} An array of vertices.
   */
  getNonRingNeighbours(vertexId) {
    const nrneighbours = []
    const vertex = this.graph.vertices[vertexId]
    const neighbours = vertex.neighbours

    for (let i = 0; i < neighbours.length; i++) {
      const neighbour = this.graph.vertices[neighbours[i]]
      const nIntersections = ArrayHelper.intersection(vertex.value.rings, neighbour.value.rings).length

      if (nIntersections === 0 && neighbour.value.isBridge === false) {
        nrneighbours.push(neighbour)
      }
    }

    return nrneighbours
  }

  /**
   * Annotaed stereochemistry information for visualization.
   */
  annotateStereochemistry() {
    const maxDepth = 10

    // For each stereo-center
    for (let i = 0; i < this.graph.vertices.length; i++) {
      const vertex = this.graph.vertices[i]

      if (!vertex.value.isStereoCenter) {
        continue
      }

      const neighbours = vertex.getNeighbours()
      const nNeighbours = neighbours.length
      const priorities = Array(nNeighbours)

      for (let j = 0; j < nNeighbours; j++) {
        const visited = new Uint8Array(this.graph.vertices.length)
        const priority = Array([])
        visited[vertex.id] = 1

        this.visitStereochemistry(neighbours[j], vertex.id, visited, priority, maxDepth, 0)

        // Sort each level according to atomic number
        for (let k = 0; k < priority.length; k++) {
          priority[k].sort(function(a, b) {
            return b - a
          })
        }

        priorities[j] = [j, priority]
      }

      let maxLevels = 0
      let maxEntries = 0
      for (let j = 0; j < priorities.length; j++) {
        if (priorities[j][1].length > maxLevels) {
          maxLevels = priorities[j][1].length
        }

        for (let k = 0; k < priorities[j][1].length; k++) {
          if (priorities[j][1][k].length > maxEntries) {
            maxEntries = priorities[j][1][k].length
          }
        }
      }

      for (let j = 0; j < priorities.length; j++) {
        const diff = maxLevels - priorities[j][1].length
        for (let k = 0; k < diff; k++) {
          priorities[j][1].push([])
        }

        // Break ties by the position in the SMILES string as per specification
        priorities[j][1].push([neighbours[j]])

        // Make all same length. Fill with zeroes.
        for (let k = 0; k < priorities[j][1].length; k++) {
          const diff = maxEntries - priorities[j][1][k].length

          for (let l = 0; l < diff; l++) {
            priorities[j][1][k].push(0)
          }
        }
      }

      priorities.sort(function(a, b) {
        for (let j = 0; j < a[1].length; j++) {
          for (let k = 0; k < a[1][j].length; k++) {
            if (a[1][j][k] > b[1][j][k]) {
              return -1
            } else if (a[1][j][k] < b[1][j][k]) {
              return 1
            }
          }
        }

        return 0
      })

      const order = new Uint8Array(nNeighbours)
      for (let j = 0; j < nNeighbours; j++) {
        order[j] = priorities[j][0]
        vertex.value.priority = j
      }

      // Check the angles between elements 0 and 1, and 0 and 2 to determine whether they are
      // drawn cw or ccw
      // TODO: OC(Cl)=[C@]=C(C)F currently fails here, however this is, IMHO, not a valid SMILES.
      const posA = this.graph.vertices[neighbours[order[0]]].position
      const posB = this.graph.vertices[neighbours[order[1]]].position
      // let posC = this.graph.vertices[neighbours[order[2]]].position;

      const cwA = posA.relativeClockwise(posB, vertex.position)
      // let cwB = posA.relativeClockwise(posC, vertex.position);

      // If the second priority is clockwise from the first, the ligands are drawn clockwise, since
      // The hydrogen can be drawn on either side
      const isCw = cwA === -1

      const rotation = vertex.value.bracket.chirality === '@' ? -1 : 1
      const rs = MathHelper.parityOfPermutation(order) * rotation === 1 ? 'R' : 'S'

      // Flip the hydrogen direction when the drawing doesn't match the chirality.
      let wedgeA = 'down'
      let wedgeB = 'up'
      if (isCw && rs !== 'R' || !isCw && rs !== 'S') {
        vertex.value.hydrogenDirection = 'up'
        wedgeA = 'up'
        wedgeB = 'down'
      }

      if (vertex.value.hasHydrogen) {
        this.graph.getEdge(vertex.id, neighbours[order[order.length - 1]]).wedge = wedgeA
      }

      // Get the shortest subtree to flip up / down. Ignore lowest priority
      // The rules are following:
      // 1. Do not draw wedge between two stereocenters
      // 2. Heteroatoms
      // 3. Draw outside ring
      // 4. Shortest subtree

      const wedgeOrder = new Array(neighbours.length - 1)
      const showHydrogen = vertex.value.rings.length > 1 && vertex.value.hasHydrogen
      const offset = vertex.value.hasHydrogen ? 1 : 0

      for (let j = 0; j < order.length - offset; j++) {
        wedgeOrder[j] = new Uint32Array(2)
        const neighbour = this.graph.vertices[neighbours[order[j]]]
        wedgeOrder[j][0] += neighbour.value.isStereoCenter ? 0 : 100000
        // wedgeOrder[j][0] += neighbour.value.rings.length > 0 ? 0 : 10000;
        // Only add if in same ring, unlike above
        wedgeOrder[j][0] += this.areVerticesInSameRing(neighbour, vertex) ? 0 : 10000
        wedgeOrder[j][0] += neighbour.value.isHeteroAtom() ? 1000 : 0
        wedgeOrder[j][0] -= neighbour.value.subtreeDepth === 0 ? 1000 : 0
        wedgeOrder[j][0] += 1000 - neighbour.value.subtreeDepth
        wedgeOrder[j][1] = neighbours[order[j]]
      }

      wedgeOrder.sort(function(a, b) {
        if (a[0] > b[0]) {
          return -1
        } else if (a[0] < b[0]) {
          return 1
        }
        return 0
      })

      // If all neighbours are in a ring, do not draw wedge, the hydrogen will be drawn.
      if (!showHydrogen) {
        const wedgeId = wedgeOrder[0][1]

        if (vertex.value.hasHydrogen) {
          this.graph.getEdge(vertex.id, wedgeId).wedge = wedgeB
        } else {
          let wedge = wedgeB

          for (let j = order.length - 1; j >= 0; j--) {
            if (wedge === wedgeA) {
              wedge = wedgeB
            } else {
              wedge = wedgeA
            }
            if (neighbours[order[j]] === wedgeId) {
              break
            }
          }

          this.graph.getEdge(vertex.id, wedgeId).wedge = wedge
        }
      }

      vertex.value.chirality = rs
    }
  }

  /**
   *
   *
   * @param {Number} vertexId The id of a vertex.
   * @param {(Number|null)} previousVertexId The id of the parent vertex of the vertex.
   * @param {Uint8Array} visited An array containing the visited flag for all vertices in the graph.
   * @param {Array} priority An array of arrays storing the atomic numbers for each level.
   * @param {Number} maxDepth The maximum depth.
   * @param {Number} depth The current depth.
   * @param parentAtomicNumber
   */
  visitStereochemistry(vertexId, previousVertexId, visited, priority, maxDepth, depth, parentAtomicNumber = 0) {
    visited[vertexId] = 1
    const vertex = this.graph.vertices[vertexId]
    const atomicNumber = vertex.value.getAtomicNumber()

    if (priority.length <= depth) {
      priority.push([])
    }

    for (let i = 0; i < this.graph.getEdge(vertexId, previousVertexId).weight; i++) {
      priority[depth].push(parentAtomicNumber * 1000 + atomicNumber)
    }

    const neighbours = this.graph.vertices[vertexId].neighbours

    for (let i = 0; i < neighbours.length; i++) {
      if (visited[neighbours[i]] !== 1 && depth < maxDepth - 1) {
        this.visitStereochemistry(neighbours[i], vertexId, visited.slice(), priority, maxDepth, depth + 1, atomicNumber)
      }
    }

    // Valences are filled with hydrogens and passed to the next level.
    if (depth < maxDepth - 1) {
      let bonds = 0

      for (let i = 0; i < neighbours.length; i++) {
        bonds += this.graph.getEdge(vertexId, neighbours[i]).weight
      }

      for (let i = 0; i < vertex.value.getMaxBonds() - bonds; i++) {
        if (priority.length <= depth + 1) {
          priority.push([])
        }

        priority[depth + 1].push(atomicNumber * 1000 + 1)
      }
    }
  }

  /**
   * Creates pseudo-elements (such as Et, Me, Ac, Bz, ...) at the position of the carbon sets
   * the involved atoms not to be displayed.
   */
  initPseudoElements() {
    for (let i = 0; i < this.graph.vertices.length; i++) {
      const vertex = this.graph.vertices[i]
      const neighbourIds = vertex.neighbours
      const neighbours = Array(neighbourIds.length)

      for (let j = 0; j < neighbourIds.length; j++) {
        neighbours[j] = this.graph.vertices[neighbourIds[j]]
      }

      // Ignore atoms that have less than 3 neighbours, except if
      // the vertex is connected to a ring and has two neighbours
      if (vertex.getNeighbourCount() < 3 || vertex.value.rings.length > 0) {
        continue
      }

      // TODO: This exceptions should be handled more elegantly (via config file?)

      // Ignore phosphates (especially for triphosphates)
      if (vertex.value.element === 'P') {
        continue
      }

      // Ignore also guanidine
      if (vertex.value.element === 'C' && neighbours.length === 3 &&
        neighbours[0].value.element === 'N' && neighbours[1].value.element === 'N' && neighbours[2].value.element === 'N') {
        continue
      }

      // Continue if there are less than two heteroatoms
      // or if a neighbour has more than 1 neighbour
      let heteroAtomCount = 0
      let ctn = 0

      for (let j = 0; j < neighbours.length; j++) {
        const neighbour = neighbours[j]
        const neighbouringElement = neighbour.value.element
        const neighbourCount = neighbour.getNeighbourCount()

        if (neighbouringElement !== 'C' && neighbouringElement !== 'H' &&
          neighbourCount === 1) {
          heteroAtomCount++
        }

        if (neighbourCount > 1) {
          ctn++
        }
      }

      if (ctn > 1 || heteroAtomCount < 2) {
        continue
      }

      // Get the previous atom (the one which is not terminal)
      let previous = null

      for (let j = 0; j < neighbours.length; j++) {
        const neighbour = neighbours[j]

        if (neighbour.getNeighbourCount() > 1) {
          previous = neighbour
        }
      }

      for (let j = 0; j < neighbours.length; j++) {
        const neighbour = neighbours[j]

        if (neighbour.getNeighbourCount() > 1) {
          continue
        }

        neighbour.value.isDrawn = false

        let hydrogens = Atom.maxBonds[neighbour.value.element] - neighbour.value.bondCount
        let charge = ''

        if (neighbour.value.bracket) {
          hydrogens = neighbour.value.bracket.hcount
          charge = neighbour.value.bracket.charge || 0
        }

        vertex.value.attachPseudoElement(neighbour.value.element, previous ? previous.value.element : null, hydrogens, charge)
      }
    }

    // The second pass
    for (let i = 0; i < this.graph.vertices.length; i++) {
      const vertex = this.graph.vertices[i]
      const atom = vertex.value
      const element = atom.element

      if (element === 'C' || element === 'H' || !atom.isDrawn) {
        continue
      }

      const neighbourIds = vertex.neighbours
      const neighbours = Array(neighbourIds.length)

      for (let j = 0; j < neighbourIds.length; j++) {
        neighbours[j] = this.graph.vertices[neighbourIds[j]]
      }

      for (let j = 0; j < neighbours.length; j++) {
        const neighbour = neighbours[j].value

        if (!neighbour.hasAttachedPseudoElements || neighbour.getAttachedPseudoElementsCount() !== 2) {
          continue
        }

        const pseudoElements = neighbour.getAttachedPseudoElements()

        // eslint-disable-next-line no-prototype-builtins
        if (pseudoElements.hasOwnProperty('0O') && pseudoElements.hasOwnProperty('3C')) {
          neighbour.isDrawn = false
          vertex.value.attachPseudoElement('Ac', '', 0)
        }
      }
    }
  }
}

module.exports = Drawer
