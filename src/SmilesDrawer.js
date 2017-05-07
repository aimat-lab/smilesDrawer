/** The main class of the application representing the smiles drawer */
class SmilesDrawer {
    /**
     * The constructor for the class SmilesDrawer.
     *
     * @param {object} options An object containing custom values for different options. It is merged with the default options.
     */
    constructor(options) {
        this.ringIdCounter = 0;
        this.ringConnectionIdCounter = 0;
        this.canvasWrapper = null;
        this.direction = 1;
        this.totalOverlapScore = 0;

        this.maxBonds = {
            'c': 4,
            'C': 4,
            'n': 3,
            'N': 3,
            'o': 2,
            'O': 2,
            'p': 3,
            'P': 3,
            's': 2,
            'S': 2,
            'b': 3,
            'B': 3,
            'F': 1,
            'I': 1,
            'Cl': 1,
            'Br': 1
        };

        this.defaultOptions = {
            width: 500,
            height: 500,
            bondLength: 16,
            shortBondLength: 9,
            bondSpacing: 4,
            atomVisualization: 'default',
            allowFlips: false,
            isomeric: false,
            debug: false,
            terminalCarbons: false,
            compactDrawing: true,
            fontSizeLarge: 6,
            fontSizeSmall: 4,
            themes: {
                dark: {
                    C: '#fff',
                    O: '#e74c3c',
                    N: '#3498db',
                    F: '#27ae60',
                    CL: '#16a085',
                    BR: '#d35400',
                    I: '#8e44ad',
                    P: '#d35400',
                    S: '#f1c40f',
                    B: '#e67e22',
                    SI: '#e67e22',
                    H: '#252525',
                    BACKGROUND: '#141414'
                },
                light: {
                    C: '#222',
                    O: '#e74c3c',
                    N: '#3498db',
                    F: '#27ae60',
                    CL: '#16a085',
                    BR: '#d35400',
                    I: '#8e44ad',
                    P: '#d35400',
                    S: '#f1c40f',
                    B: '#e67e22',
                    SI: '#e67e22',
                    H: '#d5d5d5',
                    BACKGROUND: '#fff'
                }
            }
        };

        this.opts = this.extend(true, this.defaultOptions, options);

        // Set the default theme.
        this.theme = this.opts.themes.dark;
    }

    /**
     * A helper method to extend the default options with user supplied ones.
     *
     */
    extend() {
        let that = this;
        let extended = {};
        let deep = false;
        let i = 0;
        let length = arguments.length;

        if (Object.prototype.toString.call(arguments[0]) === '[object Boolean]') {
            deep = arguments[0];
            i++;
        }

        let merge = function (obj) {
            for (let prop in obj) {
                if (Object.prototype.hasOwnProperty.call(obj, prop)) {
                    if (deep && Object.prototype.toString.call(obj[prop]) === '[object Object]') {
                        extended[prop] = that.extend(true, extended[prop], obj[prop]);
                    } else {
                        extended[prop] = obj[prop];
                    }
                }
            }
        };

        for ( ; i < length; i++) {
            let obj = arguments[i];
            merge(obj);
        }

        return extended;
    };


    /**
     * Draws the parsed smiles data to a canvas element.
     *
     * @param {object} data The tree returned by the smiles parser.
     * @param {string|HTMLElement} target The id of the HTML canvas element the structure is drawn to - or the element itself.
     * @param {string} themeName='dark' The name of the theme to use. Built-in themes are 'light' and 'dark'.
     * @param {boolean} infoOnly=false Only output info on the molecule without drawing anything to the canvas.
     */
    draw(data, target, themeName = 'light', infoOnly = false) {
        this.data = data;
        this.canvasWrapper = new CanvasWrapper(target, this.opts.themes[themeName], this.opts);
        
        this.ringIdCounter = 0;
        this.ringConnectionIdCounter = 0;

        this.vertices = [];
        this.edges = [];
        this.rings = [];
        this.ringConnections = [];

        this.originalRings = [];
        this.originalRingConnections = [];

        this.bridgedRing = false;
        
        this.initGraph(data);
        this.initRings();

        let t = performance.now();
        this.initPathIncludedDistanceMatrix();
        console.log(performance.now() - t);

        // console.log(this.distanceMatrix);
        
        if (this.opts.isomeric) {
            this.annotateChirality();
        }

        if (!infoOnly) {
            this.position();

            // Restore the ring information (removes bridged rings and replaces them with the original, multiple, rings)
            this.restoreRingInformation();
            
            let overlapScore = this.getOverlapScore();

            this.totalOverlapScore = this.getOverlapScore().total;
            
            for (let i = 0; i < this.edges.length; i++) {
                let edge = this.edges[i];
                
                if (this.isEdgeRotatable(edge)) {
                    let subTreeDepthA = this.getTreeDepth(edge.sourceId, edge.targetId);
                    let subTreeDepthB = this.getTreeDepth(edge.targetId, edge.sourceId);
                    
                    // Only rotate the shorter subtree
                    let a = edge.targetId;
                    let b = edge.sourceId;
                    let depth = subTreeDepthA;

                    if (subTreeDepthA > subTreeDepthB) {
                        a = edge.sourceId;
                        b = edge.targetId;
                        depth = subTreeDepthB;
                    }

                    let subTreeOverlap = this.getSubtreeOverlapScore(b, a, overlapScore.vertexScores);
                    
                    if (subTreeOverlap.value > 1.0) {
                        let vertexA = this.vertices[a];
                        let vertexB = this.vertices[b];
                        let neighbours = vertexB.getNeighbours(a);

                        if (neighbours.length === 1) {
                            let neighbour = this.vertices[neighbours[0]];
                            let angle = neighbour.position.getRotateAwayFromAngle(vertexA.position, vertexB.position, MathHelper.toRad(120));
                            
                            // console.log('Rotate ' + neighbour.id + ' by ' + angle + ' away from ' + vertexA.id + ' around ' + vertexB.id);
                            this.rotateSubtree(neighbour.id, vertexB.id, angle, vertexB.position);
                            
                            // If the new overlap is bigger, undo change
                            let newTotalOverlapScore = this.getOverlapScore().total;

                            if (newTotalOverlapScore > this.totalOverlapScore) {
                                this.rotateSubtree(neighbour.id, vertexB.id, -angle, vertexB.position);
                            } else {
                                this.totalOverlapScore = newTotalOverlapScore;
                            }
                        } else if (neighbours.length == 2) {
                            // Switch places / sides
                            // If vertex a is in a ring, do nothing
                            if (vertexB.value.rings.length + vertexA.value.rings.length > 0) {
                                continue;
                            }

                            let neighbourA = this.vertices[neighbours[0]];
                            let neighbourB = this.vertices[neighbours[1]];

                            let angleA = neighbourA.position.getRotateAwayFromAngle(vertexA.position, vertexB.position, MathHelper.toRad(120));
                            let angleB = neighbourB.position.getRotateAwayFromAngle(vertexA.position, vertexB.position, MathHelper.toRad(120));
                            
                            this.rotateSubtree(neighbourA.id, vertexB.id, angleA, vertexB.position);
                            this.rotateSubtree(neighbourB.id, vertexB.id, angleB, vertexB.position);

                            let newTotalOverlapScore = this.getOverlapScore().total;

                            if (newTotalOverlapScore > this.totalOverlapScore) {
                                this.rotateSubtree(neighbourA.id, vertexB.id, -angleA, vertexB.position);
                                this.rotateSubtree(neighbourB.id, vertexB.id, -angleB, vertexB.position);
                            } else {
                                this.totalOverlapScore = newTotalOverlapScore;
                            }
                        }

                        overlapScore = this.getOverlapScore();
                    }
                }
            }
            
            this.resolveSecondaryOverlaps(overlapScore.scores);
            
            // Set the canvas to the appropriate size
            this.canvasWrapper.scale(this.vertices);

            // Initialize pseudo elements or shortcuts
            if (this.opts.compactDrawing) {
                this.initPseudoElements();
            }

            // Do the actual drawing
            this.drawEdges(this.opts.debug);
            this.drawVertices(this.opts.debug);
            
            this.canvasWrapper.reset();
        }
    }

    /**
     * Initialize the adjacency matrix of the molecular graph.
     * 
     * @returns {array} The adjancency matrix of the molecular graph.
     */
    getAdjacencyMatrix() {
        let length = this.vertices.length;
        let adjacencyMatrix = Array(length);
        
        for (let i = 0; i < length; i++) {
            adjacencyMatrix[i] = new Array(length);
            adjacencyMatrix[i].fill(0);
        }

        for (let i = 0; i < this.edges.length; i++) {
            let edge = this.edges[i];

            adjacencyMatrix[edge.sourceId][edge.targetId] = 1;
            adjacencyMatrix[edge.targetId][edge.sourceId] = 1;
        }

        return adjacencyMatrix;
    }

    /**
     * Get the distance matrix (floyd marshall) of a adjacency matrix.
     * 
     * @param {array} adjacencyMatrix An adjacency matrix.
     * @returns {array} The distance matrix of the graph defined by the adjacency matrix.
     */
    getDistanceMatrix(adjacencyMatrix) {
        let length = adjacencyMatrix.length;
        let distanceMatrix = new Array(length);

        for (let i = 0; i < length; i++) {
            distanceMatrix[i] = new Array(length);

            for (let j = 0; j < length; j++) {
                distanceMatrix[i][j] = (i === j || adjacencyMatrix[i][j] === 1) ? adjacencyMatrix[i][j] : Number.POSITIVE_INFINITY;
            }
        }

        for (let i = 0; i < length; i++) {
            distanceMatrix[i][i] = 0;
        }

        for (let k = 0; k < length; k++) {
            for (let i = 0; i < length; i++) {
                for (let j = 0; j < length; j++) {
                    if (distanceMatrix[i][j] > distanceMatrix[i][k] + distanceMatrix[k][j]) {
                        distanceMatrix[i][j] = distanceMatrix[i][k] + distanceMatrix[k][j];
                    }
                }
            }
        }

        return distanceMatrix;
    }

    /**
     * Returns an edge list constructed form an adjacency matrix.
     * 
     * @param {array} adjacencyMatrix An adjacency matrix.
     * @returns {array} An edge list. E.g. [ [ 0, 1 ], ..., [ 16, 2 ] ]
     */
    getEdgeList(adjacencyMatrix) {
        let length = adjacencyMatrix.length;
        let edgeList = [];

        for (let i = 0; i < length - 1; i++) {
            for (let j = i + 1; j < length; j++) {
                if (adjacencyMatrix[i][j] === 1) {
                    edgeList.push([i,j]);
                } 
            }
        }

        return edgeList;
    }

    /**
     * Returnes the two path-included distance matrices used to find the sssr.
     * 
     * @param {array} adjacencyMatrix An adjacency matrix.
     * @returns {object} The path-included distance matrices. { p1, p2 }
     */
    getPathIncludedDistanceMatrices(adjacencyMatrix) {
        let length = adjacencyMatrix.length;
        let d = Array(length);
        let pe1 = Array(length);
        let pe2 = Array(length);

        for (let i = 0; i < length; i++) {
            d[i] = Array(length);
            pe1[i] = Array(length);
            pe2[i] = Array(length);
            
            for (let j = 0; j < length; j++) {
                d[i][j] = (i === j || adjacencyMatrix[i][j] === 1) ? adjacencyMatrix[i][j] : Number.POSITIVE_INFINITY;

                if (d[i][j] === 1) {
                    pe1[i][j] = [[[i, j]]];
                } else {
                    pe1[i][j] = [];
                }

                pe2[i][j] = [];
            }
        }

        for (let k = 0; k < length; k++) {
            for (let i = 0; i < length; i++) {
                for (let j = 0; j < length; j++) {
                    const previousPathLength = d[i][j];
                    const newPathLength = d[i][k] + d[k][j];

                    if (previousPathLength > newPathLength) {
                        if (previousPathLength === newPathLength + 1) {
                            pe2[i][j] = ArrayHelper.deepCopy(pe1[i][j]);
                        } else {
                            pe2[i][j] = [];
                        }

                        d[i][j] = newPathLength;
                        pe1[i][j] = [ pe1[i][k][0].concat(pe1[k][j][0]) ];
                    } else if (previousPathLength === newPathLength) {
                        if (pe1[i][k].length && pe1[k][j].length) {
                            if (pe1[i][j].length) {
                                pe1[i][j].push(pe1[i][k][0].concat(pe1[k][j][0]));
                            } else {
                                pe1[i][j][0] = pe1[i][k][0].concat(pe1[k][j][0]);
                            }
                        }
                    } else if (previousPathLength === newPathLength - 1) {
                        if (pe2[i][j].length) {
                            pe2[i][j].push(pe1[i][k][0].concat(pe1[k][j][0]));
                        } else {
                            pe2[i][j][0] = pe1[i][k][0].concat(pe1[k][j][0]);
                        }
                    }
                }
            }
        }

        return {
            d: d,
            pe1: pe1, 
            pe2: pe2 
        };
    }

    /**
     * Get the ring candidates from the path-included distance matrices.
     * 
     * @param {array} d The distance matrix.
     * @param {array} pe1 A matrix containing the shortest paths.
     * @param {array} pe2 A matrix containing the shortest paths + one vertex.
     * @returns {array} The ring candidates.
     */
    getRingCandidates(d, pe1, pe2) {
        let length = d.length;
        let candidates = [];
        let c = 0;

        for (let i = 0; i < length; i++) {
            for (let j = 0; j < length; j++) {
                if (d[i][j] === 0 || (pe1[i][j].length === 1 && pe2[i][j] === 0)) {
                    continue;
                } else {
                    // c is the number of vertices in the cycle.
                    if (pe2[i][j].length !== 0) {
                        c = 2 * (d[i][j] + 0.5);
                    } else {
                        c = 2 * d[i][j];
                    }
                    
                    candidates.push([c, pe1[i][j], pe2[i][2]]);
                }
            }
        }

        // Candidates have to be sorted by c
        candidates.sort(function(a, b) {
            return a[0] - b[0];
        });

        return candidates;
    }

    /**
     * Searches the candidates for the smallest set of smallest rings.
     * 
     * @param {array} c The candidates.
     * @param {array} d The distance matrix.
     * @param {array} pe1 A matrix containing the shortest paths.
     * @param {array} pe2 A matrix containing the shortest paths + one vertex.
     * @param {number} nsssr The theoretical number of rings in the graph.
     * @returns {array} The smallest set of smallest rings.
     */
    getSSSR(c, d, pe1, pe2, nsssr) {
        let cSssr = [];

        for (let i = 0; i < c.length; i++) {
            if (c[i][0] % 2 !== 0) {
                for (let j = 0; j < c[i][2].length; j++) {
                    let bonds = c[i][1][0].concat(c[i][2][j]);
                    let atoms = this.bondsToAtoms(bonds);
                    
                    if (bonds.length === atoms.size && !this.pathSetsContain(cSssr, atoms)) {
                        cSssr.push(atoms);
                    }

                    if (cSssr.length === nsssr) {
                        return cSssr;
                    }
                }
            } else {
                for (let j = 0; j < c[i][1].length - 1; j++) {
                    let bonds = c[i][1][j].concat(c[i][1][j + 1]);
                    let atoms = this.bondsToAtoms(bonds);

                    if (bonds.length === atoms.size && !this.pathSetsContain(cSssr, atoms)) {
                        cSssr.push(atoms);
                    }

                    if (cSssr.length === nsssr) {
                        return cSssr;
                    }
                }
            }
        }

        return cSssr;
    }

    /**
     * Return a set of vertex indices contained in an array of bonds.
     * 
     * @param {array} bonds An array of bonds.
     * @returns {set} An array of vertices.
     */
    bondsToAtoms(bonds) {
        let atoms = new Set();
        
        for (let i = 0; i < bonds.length; i++) {
            atoms.add(bonds[i][0]);
            atoms.add(bonds[i][1]);
        }

        return atoms;
    }

    /**
     * Checks whether or not a given path already exists in an array of paths.
     * 
     * @param {array} pathSets An array of sets each representing a path.
     * @param {set} pathSet A set representing a path.
     * @returns {boolean} A boolean indicating whether or not a give path is contained within a set.
     */
    pathSetsContain(pathSets, pathSet) {
        for (let i = 0; i < pathSets.length; i++) {
            if (pathSets[i].size !== pathSet.size) {
                continue;
            }
            
            if (this.areSetsEqual(pathSets[i], pathSet)) {
                return true;
            }
        }

        return false;
    }

    /**
     * Checks whether or not two sets are equal (contain the same elements).
     * 
     * @param {set} setA A set.
     * @param {set} setB A set.
     * @returns {boolean} A boolean indicating whether or not the two sets are equal.
     */
    areSetsEqual(setA, setB) {
        if (setA.size !== setB.size) {
            return false;
        }
        
        for (let element of setA) {
            if (!setB.has(element)) {
                return false;
            }
        }

        return true;
    }

    /**
     * Initializes the path-included distance matrix.
     */
    initPathIncludedDistanceMatrix() {
        let adjacencyMatrix = this.getAdjacencyMatrix();

        // Remove vertices that are not members of a ring
        let removed;

        do {
            removed = 0;

            for (let i = 0; i < adjacencyMatrix.length; i++) {
                let nNeighbours = adjacencyMatrix[i].reduce((a, b) => a + b, 0);
                
                if (nNeighbours === 1) {
                    adjacencyMatrix[i].fill(0);

                    for (let j = 0; j < adjacencyMatrix.length; j++) {
                        adjacencyMatrix[j][i] = 0;
                    }

                    removed++;
                }
            }            
        } while (removed > 0);

        // Update the adjacency matrix (remove rows and columns filled with 0s)
        
        // Keep this as a map of new indices to old indices
        let indices = [];
        let indicesToRemove = [];
        let updatedAdjacencyMatrix = [];

        // Only the rows are filtered here, the columns still have their original values
        for (let i = 0; i < adjacencyMatrix.length; i++) {
            if (adjacencyMatrix[i].indexOf(1) >= 0) {
                indices.push(i);
                updatedAdjacencyMatrix.push(adjacencyMatrix[i]);
            } else {
                indicesToRemove.push(i);
            }
        }

        // Remove the unused values from the adjacency matrix

        for (let i = 0; i < updatedAdjacencyMatrix.length; i++) {
            for (let j = indicesToRemove.length - 1; j >= 0; j--) {
                updatedAdjacencyMatrix[i].splice(indicesToRemove[j], 1);
            }
        }

        adjacencyMatrix = updatedAdjacencyMatrix;

        if (adjacencyMatrix.length === 0) {
            return;
        }

        // Get the edge list and the theoretical number of rings in SSSR
        let edgeList = this.getEdgeList(adjacencyMatrix);
        let nSssr = edgeList.length - adjacencyMatrix.length + 1;

        // Get the distance matrix
        let distanceMatrix = this.getDistanceMatrix(adjacencyMatrix);

        if (nSssr === 0) {
            return;
        }

        let {d, pe1, pe2} = this.getPathIncludedDistanceMatrices(adjacencyMatrix);
        
        this.printMatrix(d);
        console.log(pe1);
        console.log(pe2);
        
        let c = this.getRingCandidates(d, pe1, pe2);
        console.log(c);
        let sssr = this.getSSSR(c, d, pe1, pe2, nSssr);
        console.log(sssr);
        
        //console.log(candidates);
    }

    printMatrix(matrix) {
        let str = '';

        for (let i = 0; i < matrix.length; i++) {
            for (let j = 0; j < matrix.length; j++) {
                str += matrix[i][j] + ' ';
            }

            str += '\n';
        }

        console.log(str);
    }

    /**
     * Returns the number of rings this edge is a part of.
     *
     * @param {number} edgeId The id of an edge.
     * @returns {number} The number of rings the provided edge is part of.
     */
    edgeRingCount(edgeId) {
        let edge = this.edges[edgeId];
        let a = this.vertices[edge.sourceId];
        let b = this.vertices[edge.targetId];

        return Math.min(a.value.rings.length, b.value.rings.length);
    }

    /**
     * Returns an array containing the bridged rings associated with this  molecule.
     *
     * @returns {array} An array containing all bridged rings associated with this molecule.
     */
    getBridgedRings() {
        let bridgedRings = [];

        for (let i = 0; i < this.rings.length; i++) {
            if (this.rings[i].isBridged) {
                bridgedRings.push(this.rings[i]);
            }
        }

        return bridgedRings;
    }

    /**
     * Returns an array containing all fused rings associated with this molecule.
     *
     * @returns {array} An array containing all fused rings associated with this molecule.
     */
    getFusedRings() {
        let fusedRings = [];

        for (let i = 0; i < this.rings.length; i++) {
            if (this.rings[i].isFused) {
                fusedRings.push(this.rings[i]);
            }
        }

        return fusedRings;
    }

    /**
     * Returns an array containing all spiros associated with this molecule.
     *
     * @returns {array} An array containing all spiros associated with this molecule.
     */
    getSpiros() {
        let spiros = [];

        for (let i = 0; i < this.rings.length; i++) {
            if (this.rings[i].isSpiro) {
                spiros.push(this.rings[i]);
            }
        }
        
        return spiros;
    }

    /**
     * Returns a string containing a semicolon and new-line separated list of ring properties: Id; Members Count; Neighbours Count; IsSpiro; IsFused; IsBridged; Ring Count (subrings of bridged rings); Insiders Count (the number of vertices contained within a bridged ring)
     *
     * @returns {string} A string as described in the method description.
     */
    printRingInfo() {
        let result = '';
        for (let i = 0; i < this.rings.length; i++) {
            let ring = this.rings[i];
            result += ring.id + ';';
            result += ring.members.length + ';';
            result += ring.neighbours.length + ';';
            result += ring.isSpiro ? 'true;' : 'false;'
            result += ring.isFused ? 'true;' : 'false;'
            result += ring.isBridged ? 'true;' : 'false;'
            result += ring.rings.length + ';';
            result += ring.insiders.length;
            result += '\n';
        }

        return result;
    }

    /**
     * Returns the total overlap score of the current molecule.
     *
     * @returns {number} The overlap score.
     */
    getTotalOverlapScore() {
        return this.totalOverlapScore;
    }

    /**
     * Returns the ring count of the current molecule.
     *
     * @returns {number} The ring count.
     */
    getRingCount() {
        return this.rings.length;
    }

    /**
     * Checks whether or not the current molecule contains a bridged ring.
     *
     * @returns {boolean} A boolean indicating whether or not the current molecule contains a bridged ring.
     */
    hasBridgedRing() {
        return this.bridgedRing;
    }

    /**
     * Returns the number of heavy atoms (non-hydrogen) in the current molecule.
     *
     * @returns {number} The heavy atom count.
     */
    getHeavyAtomCount() {
        let hac = 0;
        
        for (let i = 0; i < this.vertices.length; i++) {
            if (this.vertices[i].value.element.toLowerCase() !== 'h') {
                hac++;
            }
        }

        return hac;
    }

    /**
     * Initializes the graph (vertices and edges) based on the tree supplied by the smiles parser.
     *
     * @param {object} node The current node in the parse tree.
     * @param {number} parentVertexId=null The id of the previous vertex.
     * @param {boolean} isBranch=false Whether or not the bond leading to this vertex is a branch bond. Branches are represented by parentheses in smiles (e.g. CC(O)C).
     */
    initGraph(node, order = 0, parentVertexId = null, isBranch = false) {
        // Create a new vertex object
        let atom = new Atom(node.atom.element ? node.atom.element : node.atom, node.bond);
        
        atom.branchBond = node.branchBond;
        atom.ringbonds = node.ringbonds;
        atom.bracket = node.atom.element ? node.atom : null;
        atom.setOrder(parentVertexId, order);

        let vertex = new Vertex(atom);
        let parentVertex = this.vertices[parentVertexId];
        
        this.addVertex(vertex);

        // Add the id of this node to the parent as child
        if (parentVertexId !== null) {
            vertex.setParentVertexId(parentVertexId);
            this.vertices[parentVertexId].addChild(vertex.id);

            // In addition create a spanningTreeChildren property, which later will
            // not contain the children added through ringbonds
            this.vertices[parentVertexId].spanningTreeChildren.push(vertex.id);

            // Add edge between this node and its parent
            let edge = new Edge(parentVertexId, vertex.id, 1);
            
            if (isBranch) {
                edge.bondType = vertex.value.branchBond;
            } else {
                edge.bondType = this.vertices[parentVertexId].value.bondType;
            }

            let edgeId = this.addEdge(edge);
            vertex.edges.push(edgeId);
            parentVertex.edges.push(edgeId);
        }
        
        if (atom.bracket && this.opts.isomeric) {
            for (let i = 0; i < atom.bracket.hcount; i++) {
                if (this.opts.isomeric) {
                    this.initGraph({ atom: { element: 'H', bond: '-' }, ringbonds: [] }, i + 1, vertex.id);
                }
            }
        }

        let offset = node.ringbondCount + 1;

        if (atom.bracket) {
            offset += atom.bracket.hcount;
        }
        
        for (let i = 0; i < node.branchCount; i++) {
            this.initGraph(node.branches[i], i + offset, vertex.id, true);
        }

        if (node.hasNext) {
            this.initGraph(node.next, node.branchCount + offset, vertex.id);
        }
    }

    /**
     * Returns the type of the ringbond (e.g. '=' for a double bond). The ringbond represents the break in a ring introduced when creating the MST. If the two vertices supplied as arguments are not part of a common ringbond, the method returns null.
     *
     * @param {Vertex} vertexA A vertex.
     * @param {Vertex} vertexB A vertex.
     * @returns {string|null} Returns the ringbond type or null, if the two supplied vertices are not connected by a ringbond.
     */
    getRingbondType(vertexA, vertexB) {
        // Checks whether the two vertices are the ones connecting the ring
        // and what the bond type should be.
        if (vertexA.value.getRingbondCount() < 1 || vertexB.value.getRingbondCount() < 1) {
            return null;
        }

        for (let i = 0; i < vertexA.value.ringbonds.length; i++) {
            for (let j = 0; j < vertexB.value.ringbonds.length; j++) {
                // if(i != j) continue;
                if (vertexA.value.ringbonds[i].id === vertexB.value.ringbonds[j].id) {
                    // If the bonds are equal, it doesn't matter which bond is returned.
                    // if they are not equal, return the one that is not the default ("-")
                    if (vertexA.value.ringbonds[i].bondType === '-') {
                        return vertexB.value.ringbonds[j].bond;
                    } else {
                        return vertexA.value.ringbonds[i].bond;
                    }
                }
            }
        }

        return null;
    }

    /**
     * Initializes rings and ringbonds for the current molecule.
     *
     */
    initRings() {
        let that = this;
        let openBonds = {};
        let ringId = 0;

        for (let i = this.vertices.length - 1; i >= 0; i--) {
            let vertex = this.vertices[i];
            
            if (vertex.value.ringbonds.length === 0) {
                continue;
            }

            for (let r = 0; r < vertex.value.ringbonds.length; r++) {
                let ringbondId = vertex.value.ringbonds[r].id;
                
                if (openBonds[ringbondId] === undefined) {
                    openBonds[ringbondId] = vertex.id;
                } else {
                    let target = openBonds[ringbondId];
                    let source = vertex.id;
                    let edgeId = that.addEdge(new Edge(source, target, 1));

                    let sourceVertex = that.vertices[source];
                    let targetVertex = that.vertices[target];
                    
                    sourceVertex.addChild(target);
                    targetVertex.addChild(source);

                    sourceVertex.edges.push(edgeId);
                    targetVertex.edges.push(edgeId);

                    let ring = new Ring(ringbondId, source, target);
                    that.addRing(ring);

                    // Annotate the ring (add ring members to ring and rings to vertices)
                    let path = that.getRingVertices(ring.sourceId, ring.targetId);

                    for (let j = 0; j < path.length; j++) {
                        ring.members.push(path[j]);
                        that.vertices[path[j]].value.rings.push(ring.id);
                    }

                    openBonds[ringbondId] = undefined;

                    // Add the order to the new neighbour, this is used for chirality
                    // visualization
                    let targetOffset = targetVertex.value.bracket ? targetVertex.value.bracket.hcount : 0;
                    let sourceOffset = sourceVertex.value.bracket ? sourceVertex.value.bracket.hcount : 0;

                    targetVertex.value.setOrder(source, r + 1 + sourceOffset);
                    sourceVertex.value.setOrder(target, r + 1 + targetOffset);
                }
            }
        }

        // Find connection between rings
        // Check for common vertices and create ring connections. This is a bit
        // ugly, but the ringcount is always fairly low (< 100)
        for (let i = 0; i < this.rings.length - 1; i++) {
            for (let j = i + 1; j < this.rings.length; j++) {
                let a = this.rings[i];
                let b = this.rings[j];

                let ringConnection = new RingConnection(a, b);

                // If there are no vertices in the ring connection, then there
                // is no ring connection
                if (ringConnection.vertices.length > 0) {
                    this.addRingConnection(ringConnection);
                }
            }
        }

        // Add neighbours to the rings
        for (let i = 0; i < this.rings.length; i++) {
            let ring = this.rings[i];
            ring.neighbours = RingConnection.getNeighbours(this.ringConnections, ring.id);
        }

        // Backup the ring information to restore after placing the bridged ring.
        // This is needed in order to identify aromatic rings and stuff like this in
        // rings that are member of the superring.
        this.backupRingInformation();

        // Replace rings contained by a larger bridged ring with a bridged ring
        while (this.rings.length > 0) {
            let id = -1;
            for (let i = 0; i < this.rings.length; i++) {
                let ring = this.rings[i];

                if (this.isPartOfBridgedRing(ring.id)) {
                    id = ring.id;
                }
            }
            
            if (id === -1) {
                break;
            }

            let ring = this.getRing(id);
            let involvedRings = this.getBridgedRingRings(ring.id);

            this.bridgedRing = true;
            this.createBridgedRing(involvedRings, ring.sourceId);

            // Remove the rings
            for (let i = 0; i < involvedRings.length; i++) {
                this.removeRing(involvedRings[i]);
            }
        }
    }

    /**
     * Returns all rings connected by bridged bonds starting from the ring with the supplied ring id.
     *
     * @param {number} ringId A ring id.
     * @returns {array} An array containing all ring ids of rings part of a bridged ring system.
     */
    getBridgedRingRings(ringId) {
        let involvedRings = new Array();
        let that = this;

        let recurse = function (r) {
            let ring = that.getRing(r);
            
            involvedRings.push(r);

            for (let i = 0; i < ring.neighbours.length; i++) {
                let n = ring.neighbours[i];
                
                if (involvedRings.indexOf(n) === -1 &&
                    n !== r &&
                    RingConnection.isBridge(that.ringConnections, that.vertices, r, n)) {
                    recurse(n);
                }
            }
        };

        recurse(ringId);

        return ArrayHelper.unique(involvedRings);
    }

    /**
     * Checks whether or not a ring is part of a bridged ring.
     *
     * @param {number} ringId A ring id.
     * @returns {boolean} A boolean indicating whether or not the supplied ring (by id) is part of a bridged ring system.
     */
    isPartOfBridgedRing(ringId) {
        for (let i = 0; i < this.ringConnections.length; i++) {
            if (this.ringConnections[i].rings.contains(ringId) &&
                this.ringConnections[i].isBridge(this.vertices)) {
                
                return true;
            }
        }

        return false;
    }

    /**
     * Creates a bridged ring.
     *
     * @param {array} ringIds An array of ids of rings involved in the bridged ring.
     * @param {number} sourceVertexId The vertex id to start the bridged ring discovery from.
     * @returns {Ring} The bridged ring.
     */
    createBridgedRing(ringIds, sourceVertexId) {
        let bridgedRing = new Array();
        let vertices = new Array();
        let neighbours = new Array();
        let ringConnections = new Array();

        for (let i = 0; i < ringIds.length; i++) {
            let ring = this.getRing(ringIds[i]);
            
            for (let j = 0; j < ring.members.length; j++) {
                vertices.push(ring.members[j]);
            }

            for (let j = 0; j < ring.neighbours.length; j++) {
                neighbours.push(ring.neighbours[j]);
            }
        }

        // Remove duplicates
        vertices = ArrayHelper.unique(vertices);

        // A vertex is part of the bridged ring if it only belongs to
        // one of the rings (or to another ring
        // which is not part of the bridged ring).
        let leftovers = new Array();
        
        for (let i = 0; i < vertices.length; i++) {
            let vertex = this.vertices[vertices[i]];
            let intersection = ArrayHelper.intersection(ringIds, vertex.value.rings);
            
            if (vertex.value.rings.length == 1 || intersection.length == 1) {
                bridgedRing.push(vertex.id);
            } else {
                leftovers.push(vertex.id);
            }
        }

        // Vertices can also be part of multiple rings and lay on the bridged ring,
        // however, they have to have at least two neighbours that are not part of
        // two rings
        let tmp = new Array();
        let insideRing = new Array();

        for (let i = 0; i < leftovers.length; i++) {
            let vertex = this.vertices[leftovers[i]];
            let onRing = false;

            /*
            if (ArrayHelper.intersection(vertex.getNeighbours(), bridgedRing).length > 1) {
                vertex.value.isBridgeNode = true;
                tmp.push(vertex.id);
            } else {
                vertex.value.isBridge = true;
                insideRing.push(vertex.id);
            }
            */

            for(let j = 0; j < vertex.edges.length; j++) {
                if(this.edgeRingCount(vertex.edges[j]) == 1) {
                    onRing = true;
                }
            }

            if(onRing) {
                vertex.value.isBridgeNode = true;
                tmp.push(vertex.id);
            } else {
                vertex.value.isBridge = true;
                insideRing.push(vertex.id);
            }
        }

        // Merge the two arrays containing members of the bridged ring
        let ringMembers = ArrayHelper.merge(bridgedRing, tmp)

        // The neighbours of the rings in the bridged ring that are not connected by a
        // bridge are now the neighbours of the bridged ring
        neighbours = ArrayHelper.unique(neighbours);
        neighbours = ArrayHelper.removeAll(neighbours, ringIds);

        // The source vertex is the start vertex. The target vertex has to be a member
        // of the birdged ring and a neighbour of the start vertex
        let source = this.vertices[sourceVertexId];
        let sourceNeighbours = source.getNeighbours();
        let target = null;

        for (let i = 0; i < sourceNeighbours.length; i++) {
            let n = sourceNeighbours[i];
            
            if (ringMembers.indexOf(n) !== -1) {
                target = n;
            }
        }
        
        // Create the ring
        let ring = new Ring(-1, sourceVertexId, target);
        
        ring.isBridged = true;
        ring.members = ringMembers;
        ring.neighbours = neighbours;
        ring.insiders = insideRing;
        
        for(let i = 0; i < ringIds.length; i++) {
            ring.rings.push(this.getRing(ringIds[i]).clone());
        }

        this.addRing(ring);

        this.vertices[sourceVertexId].value.anchoredRings.push(ring.id);

        // Atoms inside the ring are no longer part of a ring but are now
        // associated with the bridged ring
        for (let i = 0; i < insideRing.length; i++) {
            let vertex = this.vertices[insideRing[i]];
            
            vertex.value.rings = new Array();
            vertex.value.anchoredRings = new Array();
            vertex.value.bridgedRing = ring.id;
        }

        // Remove former rings from members of the bridged ring and add the bridged ring
        for (let i = 0; i < ringMembers.length; i++) {
            let vertex = this.vertices[ringMembers[i]];
            
            vertex.value.rings = ArrayHelper.removeAll(vertex.value.rings, ringIds);
            vertex.value.rings.push(ring.id);
        }

        // Remove all the ring connections no longer used
        for (let i = 0; i < ringIds.length; i++) {
            for (let j = i + 1; j < ringIds.length; j++) {
                this.removeRingConnectionsBetween(ringIds[i], ringIds[j]);
            }
        }

        // Update the ring connections and add this ring to the neighbours neighbours
        for (let i = 0; i < neighbours.length; i++) {
            let connections = this.getRingConnections(neighbours[i], ringIds);
            
            for (let j = 0; j < connections.length; j++) {
                this.getRingConnection(connections[j]).updateOther(ring.id, neighbours[i]);
            }

            this.getRing(neighbours[i]).neighbours.push(ring.id);
        }
        
        return ring;
    }

    /**
     * Returns an array of vertices that are members of the ring specified by the source and target vertex ids. It is assumed that those two vertices share the ringbond (the break introduced when creating the smiles MST).
     *
     * @param {number} sourceId A vertex id.
     * @param {number} targetId A vertex id.
     * @returns {array} An array of vertex ids.
     */
    getRingVertices(sourceId, targetId) {
        let prev = this.dijkstra(sourceId, targetId);

        // Backtrack from target to source
        let tmp = [];
        let path = [];
        let u = targetId;

        while (u != null) {
            tmp.push(u);
            u = prev[u];
        }

        // Reverse the backtrack path to get forward path
        for (let i = tmp.length - 1; i >= 0; i--) {
            path.push(tmp[i]);
        }

        return path;
    }

    /**
     * Dijkstras algorithm for finding the shortest path between two vertices.
     *
     * @param {number} sourceId The id of the source vertex.
     * @param {number} targetId The id of the target vertex.
     * @returns {array} The path (vertex ids) from the source to the target vertex.
     */
    dijkstra(sourceId, targetId) {
        // First initialize q which contains all the vertices
        // including their neighbours, their id and a visited boolean
        let prev = new Array(this.vertices.length);
        let dist = new Array(this.vertices.length);
        let visited = new Array(this.vertices.length);
        let neighbours = new Array(this.vertices.length);

        // Initialize arrays for the algorithm
        for (let i = 0; i < this.vertices.length; i++) {
            dist[i] = i === sourceId ? 0 : Number.MAX_VALUE;
            prev[i] = null;
            visited[i] = false;
            neighbours[i] = this.vertices[i].getNeighbours();
        }

        // Dijkstras alogrithm
        while (ArrayHelper.count(visited, false) > 0) {
            let u = this.getMinDist(dist, visited);

            // if u is the target, we're done
            if (u == targetId) { 
                return prev;
            }

            visited[u] = true; // this "removes" the node from q

            for (let i = 0; i < neighbours[u].length; i++) {
                let v = neighbours[u][i];
                let tmp = dist[u] + this.getEdgeWeight(u, v);

                // Do not move directly from the source to the target
                // this should never happen, so just continue
                if (u == sourceId && v == targetId || u == targetId && v == sourceId) {
                    continue;
                }

                if (tmp < dist[v]) {
                    dist[v] = tmp;
                    prev[v] = u;
                }
            }
        }
    }

    /**
     * Gets the minimal distance from an array containing distances.
     *
     * @param {array} dist An array of distances.
     * @param {array} visited An array indicated whether or not a vertex has been visited.
     * @returns {number} The id with the minimal distance.
     */
    getMinDist(dist, visited) {
        let min = Number.MAX_VALUE;
        let v = null;

        for (let i = 0; i < dist.length; i++) {
            if (visited[i]) {
                continue;
            }
            else if (dist[i] < min) {
                v = i;
                min = dist[v];
            }
        }

        return v;
    }

    /**
     * Checks whether or not two vertices are in the same ring.
     *
     * @param {Vertex} vertexA A vertex.
     * @param {Vertex} vertexB A vertex.
     * @returns {boolean} A boolean indicating whether or not the two vertices are in the same ring.
     */
    areVerticesInSameRing(vertexA, vertexB) {
        // This is a little bit lighter (without the array and push) than
        // getCommonRings().length > 0
        for (let i = 0; i < vertexA.value.rings.length; i++) {
            for (let j = 0; j < vertexB.value.rings.length; j++) {
                if (vertexA.value.rings[i] == vertexB.value.rings[j]) {
                    return true;
                }
            }
        }

        return false;
    }

    /**
     * Returns an array of ring ids shared by both vertices.
     *
     * @param {Vertex} vertexA A vertex.
     * @param {Vertex} vertexB A vertex.
     * @returns {array} An array of ids of rings shared by the two vertices.
     */
    getCommonRings(vertexA, vertexB) {
        let commonRings = [];

        for (let i = 0; i < vertexA.value.rings.length; i++) {
            for (let j = 0; j < vertexB.value.rings.length; j++) {
                if (vertexA.value.rings[i] == vertexB.value.rings[j]) {
                    commonRings.push(vertexA.value.rings[i]);
                }
            }
        }

        return commonRings;
    }

    /**
     * Returns the smallest ring shared by the two vertices.
     *
     * @param {Vertex} vertexA A vertex.
     * @param {Vertex} vertexB A vertex.
     * @returns {Ring|null} If a smallest common ring exists, that ring, else null.
     */
    getSmallestCommonRing(vertexA, vertexB) {
        let commonRings = this.getCommonRings(vertexA, vertexB);
        let minSize = Number.MAX_VALUE;
        let smallestCommonRing = null;

        for (let i = 0; i < commonRings.length; i++) {
            let size = this.getRing(commonRings[i]).getSize();
            
            if (size < minSize) {
                minSize = size;
                smallestCommonRing = this.getRing(commonRings[i]);
            }
        }

        return smallestCommonRing;
    }

    /**
     * Returns the largest ring shared by the two vertices.
     *
     * @param {Vertex} vertexA A vertex.
     * @param {Vertex} vertexB A vertex.
     * @returns {Ring|null} If a largest common ring exists, that ring, else null.
     */
    getLargestCommonRing(vertexA, vertexB) {
        let commonRings = this.getCommonRings(vertexA, vertexB);
        let maxSize = 0;
        let largestCommonRing = null;

        for (let i = 0; i < commonRings.length; i++) {
            let size = this.getRing(commonRings[i]).getSize();
            
            if (size > maxSize) {
                maxSize = size;
                largestCommonRing = this.getRing(commonRings[i]);
            }
        }

        return largestCommonRing;
    }

    /**
     * Returns the aromatic or largest ring shared by the two vertices.
     *
     * @param {Vertex} vertexA A vertex.
     * @param {Vertex} vertexB A vertex.
     * @returns {Ring|null} If an aromatic common ring exists, that ring, else the largest (non-aromatic) ring, else null.
     */
    getLargestOrAromaticCommonRing(vertexA, vertexB) {
        let commonRings = this.getCommonRings(vertexA, vertexB);
        let maxSize = 0;
        let largestCommonRing = null;

        for (let i = 0; i < commonRings.length; i++) {
            let ring = this.getRing(commonRings[i]);
            let size = ring.getSize();

            if (ring.isBenzeneLike(this.vertices)) {
                return ring;
            } else if (size > maxSize) {
                maxSize = size;
                largestCommonRing = ring;
            }
        }

        return largestCommonRing;
    }

     /**
     * Returns an array of vertices positioned at a specified location.
     *
     * @param {Vector2} position The position to search for vertices.
     * @param {number} radius The radius within to search.
     * @param {number} excludeVertexId A vertex id to be excluded from the search results.
     * @returns {array} An array containing vertex ids in a given location.
     */
    getVerticesAt(position, radius, excludeVertexId) {
        let locals = new Array();

        for (let i = 0; i < this.vertices.length; i++) {
            let vertex = this.vertices[i];
            
            if (vertex.id === excludeVertexId || !vertex.positioned) {
                continue;
            }

            let distance = position.distance(vertex.position);
           
            if (distance <= radius) {
                locals.push(vertex.id);
            }
        }

        return locals;
    }

    /**
     * Returns the closest vertex (connected as well as unconnected).
     *
     * @param {Vertex} vertex The vertex of which to find the closest other vertex.
     * @returns {Vertex} The closest vertex.
     */
    getClosestVertex(vertex) {
        let minDist = 99999;
        let minVertex = null;

        for (let i = 0; i < this.vertices.length; i++) {
            let v = this.vertices[i];

            if (v.id === vertex.id) {
                continue;
            }

            let distSq = vertex.position.distanceSq(v.position);

            if (distSq < minDist) {
                minDist = distSq;
                minVertex = v;
            }
        }

        return minVertex;
    }

    /**
     * Returns the closest vertex (connected as well as unconnected), which is an endpoint.
     *
     * @param {Vertex} vertex The vertex of which to find the closest other vertex.
     * @returns {Vertex} The closest endpoint vertex.
     */
    getClosestEndpointVertex(vertex) {
        let minDist = 99999;
        let minVertex = null;

        for (let i = 0; i < this.vertices.length; i++) {
            let v = this.vertices[i];

            if (v.id === vertex.id || v.getNeighbourCount() > 1) {
                continue;
            }

            let distSq = vertex.position.distanceSq(v.position);

            if (distSq < minDist) {
                minDist = distSq;
                minVertex = v;
            }
        }

        return minVertex;
    }

    /**
     * Returns the rings and vertices contained in a sub-graph.
     *
     * @param {number} vertexId The vertex id to start the sub-graph search from
     * @param {number} previousId The vertex id in the opposite of which the search will be started.
     * @returns {object} An object containing two arrays, one with the vertices in the subgraph and one with the rings in the subgraph.
     */
    getBranch(vertexId, previousId) {
        let vertices = new Array();
        let rings = new Array();
        let that = this;
        
        let recurse = function (v, p) {
            let vertex = that.vertices[v];
            
            for (let i = 0; i < vertex.value.rings.length; i++) {
                rings.push(vertex.value.rings[i]);
            }

            for (let i = 0; i < vertex.children.length; i++) {
                let child = vertex.children[i];
                
                if (child !== p && !ArrayHelper.contains(vertices, { value: child })) {
                    vertices.push(child);
                    recurse(child, v);
                }
            }

            let parentVertexId = vertex.parentVertexId;
            
            if (parentVertexId !== p && parentVertexId !== null && 
                !ArrayHelper.contains(vertices, { value: parentVertexId })) {
                vertices.push(parentVertexId);
                recurse(parentVertexId, v);
            }
        }

        vertices.push(vertexId);
        recurse(vertexId, previousId);

        return {
            vertices: vertices,
            rings: ArrayHelper.unique(rings)
        };
    }

    /**
     * Add a vertex to this representation of a molcule.
     *
     * @param {Vertex} vertex A new vertex.
     * @returns {number} The vertex id of the new vertex.
     */
    addVertex(vertex) {
        vertex.id = this.vertices.length;
        this.vertices.push(vertex);
        
        return vertex.id;
    }

    /**
     * Add an edge to this representation of a molecule.
     *
     * @param {Edge} edge A new edge.
     * @returns {number} The edge id of the new edge.
     */
    addEdge(edge) {
        edge.id = this.edges.length;
        this.edges.push(edge);
        
        return edge.id;
    }

    /**
     * Add a ring to this representation of a molecule.
     *
     * @param {Ring} ring A new ring.
     * @returns {number} The ring id of the new ring.
     */
    addRing(ring) {
        ring.id = this.ringIdCounter++;
        this.rings.push(ring);
        
        return ring.id;
    }

    /**
     * Removes a ring from the array of rings associated with the current molecule.
     *
     * @param {number} ringId A ring id.
     */
    removeRing(ringId) {
        this.rings = this.rings.filter(function (item) {
            return item.id !== ringId;
        });

        // Also remove ring connections involving this ring
        this.ringConnections = this.ringConnections.filter(function (item) {
            return item.rings.first !== ringId && item.rings.second !== ringId;
        });

        // Remove the ring as neighbour of other rings
        for (let i = 0; i < this.rings.length; i++) {
            let r = this.rings[i];
            r.neighbours = r.neighbours.filter(function (item) {
                return item !== ringId;
            });
        }
    }
    
    /**
     * Gets a ring object from the array of rings associated with the current molecule by its id. The ring id is not equal to the index, since rings can be added and removed when processing bridged rings.
     *
     * @param {number} ringId A ring id.
     * @returns {Ring} A ring associated with the current molecule.
     */
    getRing(ringId) {
        for (let i = 0; i < this.rings.length; i++) {
            if (this.rings[i].id == ringId) {
                return this.rings[i];
            }
        }
    }
    
    /**
     * Add a ring connection to this representation of a molecule.
     *
     * @param {RingConnection} ringConnection A new ringConnection.
     * @returns {number} The ring connection id of the new ring connection.
     */
    addRingConnection(ringConnection) {
        ringConnection.id = this.ringConnectionIdCounter++;
        this.ringConnections.push(ringConnection);
        
        return ringConnection.id;
    }
    
    /**
     * Removes a ring connection from the array of rings connections associated with the current molecule.
     *
     * @param {number} ringConnectionId A ring connection id.
     */
    removeRingConnection(ringConnectionId) {
        this.ringConnections = this.ringConnections.filter(function (item) {
            return item.id !== ringConnectionId;
        });
    }

    /**
     * Removes all ring connections between two vertices.
     *
     * @param {number} vertexIdA A vertex id.
     * @param {number} vertexIdB A vertex id.
     */
    removeRingConnectionsBetween(vertexIdA, vertexIdB) {
        let toRemove = new Array();
        for (let i = 0; i < this.ringConnections.length; i++) {
            let ringConnection = this.ringConnections[i];

            if (ringConnection.rings.first === vertexIdA && ringConnection.rings.second === vertexIdB ||
                ringConnection.rings.first === vertexIdB && ringConnection.rings.second === vertexIdA) {
                toRemove.push(ringConnection.id);
            }
        }

        for (let i = 0; i < toRemove.length; i++) {
            this.removeRingConnection(toRemove[i]);
        }
    }


    getRingConnection(id) {
        for (let i = 0; i < this.ringConnections.length; i++) {
            if (this.ringConnections[i].id == id) {
                return this.ringConnections[i];
            }
        }
    }

    /**
     * Get the ring connections associated with a ring, the ring connections between two rings or the ring connections between one ring and multiple other rings.
     *
     * @param {number} ringId A ring id.
     * @param {number|array|null} ringIds=null A ring id, an array of ring ids or null.
     * @returns {array} An array of ring connection ids.
     */
    getRingConnections(ringId, ringIds = null) {
        let ringConnections = new Array();
        
        if (ringIds === null) {
            for (let i = 0; i < this.ringConnections.length; i++) {
                let ringConnection = this.ringConnections[i];
                
                if (ringConnection.rings.first === ringId || ringConnection.rings.second === ringId) {
                    ringConnections.push(ringConnection.id);
                }
            }
        } else if (ringIds.constructor !== Array) {
            for (let i = 0; i < this.ringConnections.length; i++) {
                let ringConnection = this.ringConnections[i];
                
                if (ringConnection.rings.first === ringId && ringConnection.rings.second === ringIds ||
                    ringConnection.rings.first === ringIds && ringConnection.rings.second === ringId) {
                    ringConnections.push(ringConnection.id);
                }
            }
        } else {
            for (let i = 0; i < this.ringConnections.length; i++) {
                for (let j = 0; j < ringIds.length; j++) {
                    let id = ringIds[j];
                    let ringConnection = this.ringConnections[i];
                    
                    if (ringConnection.rings.first === ringId && ringConnection.rings.second === id ||
                        ringConnection.rings.first === id && ringConnection.rings.second === ringId) {
                        ringConnections.push(ringConnection.id);
                    }
                }
            }
        }

        return ringConnections;
    }

    /**
     * Check whether or not the two vertices specified span a bond which is a ring connection (fused rings).
     * 
     * @param {number} vertexIdA A vertex id.
     * @param {number} vertexIdB A vertex id.
     * @returns {boolean} Returns a boolean indicating whether or not the two vertices specify a ringbond.
     */
    isRingConnection(vertexIdA, vertexIdB) {
        for (let i = 0; i < this.ringConnections.length; i++) {
            let ringConnection = this.ringConnections[i];
            
            if (ringConnection.vertices.length !== 2) {
                continue;
            }

            if (ringConnection.vertices[0] === vertexIdA && ringConnection.vertices[1] === vertexIdB ||
                ringConnection.vertices[0] === vertexIdB && ringConnection.vertices[1] === vertexIdA) {
                return true;
            }
        }

        return false;
    }

    /**
     * Returns the overlap score of the current molecule based on its positioned vertices. The higher the score, the more overlaps occur in the structure drawing.
     *
     * @returns {object} Returns the total overlap score and the overlap score of each vertex sorted by score (higher to lower). Example: { total: 99, scores: [ { id: 0, score: 22 }, ... ]  }
     */
    getOverlapScore() {
        let total = 0.0;
        let overlapScores = new Float32Array(this.vertices.length);
        
        for (let i = 0; i < this.vertices.length; i++) {
            overlapScores[i] = 0;
        }

        for (let i = 0; i < this.vertices.length; i++) {
            for (let j = i + 1; j < this.vertices.length; j++) {
                let a = this.vertices[i];
                let b = this.vertices[j];

                let dist = Vector2.subtract(a.position, b.position).length();
                
                if (dist < this.opts.bondLength) {
                    let weighted = (this.opts.bondLength - dist) / this.opts.bondLength;
                    total += weighted;
                    overlapScores[i] += weighted;
                    overlapScores[j] += weighted;
                }
            }
        }

        let sortable = [];

        for (let i = 0; i < this.vertices.length; i++) {
            sortable.push({
                id: i,
                score: overlapScores[i]
            });
        }

        sortable.sort(function (a, b) {
            return b.score - a.score;
        });

        return {
            total: total,
            scores: sortable,
            vertexScores: overlapScores
        };
    }
    
    /**
     * When drawing a double bond, choose the side to place the double bond. E.g. a double bond should always been drawn inside a ring.
     *
     * @param {Vertex} vertexA A vertex.
     * @param {Vertex} vertexB A vertex.
     * @param {array} sides An array containing the two normals of the line spanned by the two provided vertices.
     * @returns {object} Returns an object containing the following information: {
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
        let an = vertexA.getNeighbours(vertexB.id);
        let bn = vertexB.getNeighbours(vertexA.id);
        let anCount = an.length;
        let bnCount = bn.length;

        // All vertices connected to the edge vertexA to vertexB
        let tn = ArrayHelper.merge(an, bn);

        // Only considering the connected vertices
        let sideCount = [0, 0];

        for (let i = 0; i < tn.length; i++) {
            let v = this.vertices[tn[i]].position;
            
            if (v.sameSideAs(vertexA.position, vertexB.position, sides[0])) {
                sideCount[0]++;
            } else {
                sideCount[1]++;
            }
        }

        // Considering all vertices in the graph, this is to resolve ties
        // from the above side counts
        let totalSideCount = [0, 0];

        for (let i = 0; i < this.vertices.length; i++) {
            let v = this.vertices[i].position;
            
            if (v.sameSideAs(vertexA.position, vertexB.position, sides[0])) {
                totalSideCount[0]++;
            } else {
                totalSideCount[1]++;
            }
        }

        return {
            totalSideCount: totalSideCount,
            totalPosition: totalSideCount[0] > totalSideCount[1] ? 0 : 1,
            sideCount: sideCount,
            position: sideCount[0] > sideCount[1] ? 0 : 1,
            anCount: anCount,
            bnCount: bnCount
        };
    }
    
    /**
     * Checks whether or not two vertices are connected.
     *
     * @param {number} vertexIdA A vertex id.
     * @param {number} vertexIdA A vertex id.
     * @returns {boolean} A boolean indicating whether or not two vertices are connected.
     */
    areConnected(vertexIdA, vertexIdB) {
        for(let i = 0; i < this.edges.length; i++) {
            let edge = this.edges[i];
            
            if(edge.sourceId === vertexIdA && edge.targetId === vertexIdB || 
               edge.sourceId === vertexIdB && edge.targetId === vertexIdA) {
                return true;
            }
        }
        return false;
    }
    
    /**
     * Returns the weight of the edge between two given vertices.
     *
     * @param {number} vertexIdA A vertex id.
     * @param {number} vertexIdB A vertex id.
     * @returns {number|null} The weight of the edge or, if no edge can be found, null.
     */
    getEdgeWeight(vertexIdA, vertexIdB) {
        for (let i = 0; i < this.edges.length; i++) {
            let edge = this.edges[i];
            
            if (edge.sourceId == vertexIdA && edge.targetId == vertexIdB || 
                edge.targetId == vertexIdA && edge.sourceId == vertexIdB) {
                return edge.weight;
            }
        }
        
        return null;
    }
    
    /**
     * Returns the edge between two given vertices.
     *
     * @param {number} vertexIdA A vertex id.
     * @param {number} vertexIdB A vertex id.
     * @returns {number|null} The edge or, if no edge can be found, null.
     */
    getEdge(vertexIdA, vertexIdB) {
        for (let i = 0; i < this.edges.length; i++) {
            let edge = this.edges[i];
            
            if (edge.sourceId == vertexIdA && edge.targetId == vertexIdB || 
                edge.targetId == vertexIdA && edge.sourceId == vertexIdB) {
                return edge;
            }
        }
        
        return null;
    }

    /**
     * Applies a force-based layout to a set of provided vertices.
     *
     * @param {array} vertices An array containing vertices to be placed using the force based layout.
     * @param {Vector2} center The center of the layout.
     * @param {number} startVertexId A vertex id. Should be the starting vertex - e.g. the first to be positioned and connected to a previously place vertex.
     * @param {Ring} ring The bridged ring associated with this force-based layout.
     */
    forceLayout(vertices, center, startVertexId, ring) {
        // Constants
        const l = this.opts.bondLength;

        let startVertex = this.vertices[startVertexId];
        let startVertexNeighbours = startVertex.getNeighbours();

        // Add neighbours that are already positioned to the vertices to prevent overlap
        for (let i = 0; i < startVertexNeighbours.length; i++) {
            if (this.vertices[startVertexNeighbours[i]].positioned) {
                vertices.push(startVertexNeighbours[i]);
            }
        }

        // Create adjencency matrix
        let totalLength = vertices.length + ring.rings.length;
        let vToId = new Array(vertices.length);
        let idToV = {};
        let adjMatrix = new Array(totalLength);
        let edges = new Array();
        
        for (let i = 0; i < totalLength; i++) {
            adjMatrix[i] = new Array(totalLength);

            for (let j = 0; j < totalLength; j++) {
                adjMatrix[i][j] = 0;
            }
        }

        for (let i = 0; i < vertices.length; i++) {
            vToId[i] = this.vertices[vertices[i]].id; 
            idToV[vToId[i]] = i;
        }

        for (let i = 0; i < vertices.length - 1; i++) { 
            for (let j = i; j < vertices.length; j++) {
                let edge = this.getEdge(vToId[i], this.vertices[vertices[j]].id);
                
                if (edge !== null)  {
                    adjMatrix[i][j] = l;
                    adjMatrix[j][i] = l;
                    edges.push([i, j]);
                }
            }
        }

        for (let i = 0; i < ring.rings.length; i++) {
            let r = ring.rings[i];
            let index = vertices.length + i;

            for (let j = 0; j < r.members.length; j++) {
                let id = idToV[r.members[j]];
                let radius = MathHelper.polyCircumradius(l, r.getSize());
                
                adjMatrix[id][index] = radius;
                adjMatrix[index][id] = radius;
            }
        }

        for (let i = 0; i < edges.length; i++) {
            for (let j = 0; j < totalLength; j++) {
                adjMatrix[j].push(0);
            }

            adjMatrix.push(new Array());
        
            for (let j = 0; j < totalLength + edges.length; j++) {
                adjMatrix[totalLength + i].push(0);
            }
        }

        // Connect ring centers with edges 
        for (let i = 0; i < ring.rings.length; i++) {
            let r = ring.rings[i];
            let ringIndex = vertices.length + i;
            let ringSize = r.getSize();
            
            for (let j = 0; j < edges.length; j++) {
                let a = edges[j][0];

                // If the vertex and the ring are connected, so must the edge be
                if (adjMatrix[ringIndex][a] !== 0) {
                    let apothem = MathHelper.apothem(adjMatrix[ringIndex][a], ringSize);
                    
                    adjMatrix[ringIndex][totalLength + j] = apothem;
                    adjMatrix[totalLength + j][ringIndex] = apothem;
                }
            }

            // Connecting ring centers, let them have a distance of apothem + apothem
            for (let j = 0; j < ring.rings.length; j++) {
                let r2 = ring.rings[j];

                if (r2.id === r.id) {
                    continue;
                }
                
                // If they do not share a vertex, they are not connected
                let intersection = ArrayHelper.intersection(r.members, r2.members).length;

                if (intersection === 0) {
                    continue;
                }

                let ringIndex2 = vertices.length + j;
                let ringSize2 = r2.getSize();
                let dist = MathHelper.apothemFromSideLength(l, ringSize) + MathHelper.apothemFromSideLength(l, ringSize2);
                
                adjMatrix[ringIndex][ringIndex2] = dist;
            }
        }


        
        totalLength += edges.length;
        
        let edgeOffset = totalLength - edges.length;
        
        let forces = new Array(totalLength);
        let positions = new Array(totalLength);
        let positioned = new Array(totalLength);
        let isRingCenter = new Array(totalLength);
        let ringSize = new Array(totalLength);
        let ringCount = new Array(totalLength);

        for (let i = 0; i < totalLength; i++) {
            isRingCenter[i] = i >= vertices.length && i < edgeOffset;

            ringCount[i] = i < vertices.length ? this.vertices[vToId[i]].value.originalRings.length : 1;

            if (isRingCenter[i]) {
                ringSize[i] = ring.rings[i - vertices.length].members.length;
            } else {
                ringSize[i] = 1;
            }
        }
        
        for (let i = 0; i < totalLength; i++) {
            forces[i] = new Vector2();
            positions[i] = new Vector2(center.x + Math.random() * l, center.y + Math.random() * l);
            positioned[i] = false;

            if (i >= vertices.length) {
                continue;
            }

            let vertex = this.vertices[vToId[i]];
            positions[i] = vertex.position.clone();

            // If the ring size is larger than 2, then put all the non-positioned
            // vertices at the center of the ring instead of 0,0
            if (vertex.position.x === 0 && vertex.position.y === 0) {
                // positions[i] = new Vector2(center.x + Math.random() * l, center.y + Math.random() * l);
            }
            
            if (vertex.positioned && ring.rings.length === 2) {
                positioned[i] = true;
            }
        }
        
        let k = l / 1.4;
        let c = 0.005;
        let maxMove = l / 2.0;
        let maxDist = l * 2.0;
        
        for (let n = 0; n < 600; n++) {
            for (let i = 0; i < totalLength; i++) {
                forces[i].set(0, 0);
            }

            // Set the positions of the edge midpoints
            for (let i = 0; i < edges.length; i++) {
                let index = edgeOffset + i;
                let a = positions[edges[i][0]];
                let b = positions[edges[i][1]];

                positions[index] = Vector2.midpoint(a, b)
            }

            // Repulsive forces
            for (let u = 0; u < totalLength - 1; u++) {
                for (let v = u + 1; v < totalLength; v++) {
                    if (n <= 250 && !(isRingCenter[u] && isRingCenter[v])) {
                        continue;
                    }

                    if (n > 250 && isRingCenter[u] && isRingCenter[v]) {
                        continue;
                    }

                    if (ring.rings.length < 3 && 
                        (isRingCenter[u] || isRingCenter[v])) {
                        continue;
                    }

                    let dx = positions[v].x - positions[u].x;
                    let dy = positions[v].y - positions[u].y;

                    if (dx === 0 || dy === 0) {
                        continue;
                    }

                    let dSq = dx * dx + dy * dy;

                    if (dSq < 0.01) {
                        dx = 0.1 * Math.random() + 0.1;
                        dy = 0.1 * Math.random() + 0.1;

                        dSq = dx * dx + dy * dy;
                    }

                    let d = Math.sqrt(dSq);

                    if (d > adjMatrix[u][v] && n > 200) {
                        continue;
                    }

                    let force = k * k / d;

                    if (n <= 200) {
                        force *= ringSize[u] * ringSize[v];
                    }

                    if (n > 250 && (isRingCenter[u] || isRingCenter[v])) {
                        force *= ringSize[u] * ringSize[v];
                    }

                    let fx = force * dx / d;
                    let fy = force * dy / d;

                    if (!positioned[u]) {
                        forces[u].x -= fx;
                        forces[u].y -= fy;
                    }

                    if (!positioned[v]) {
                        forces[v].x += fx;
                        forces[v].y += fy;
                    }
                }
            }

            // Attractive forces
            for (let u = 0; u < totalLength - 1; u++) {
                for (let v = u + 1; v < totalLength; v++) {
                    if (adjMatrix[u][v] <= 0) {
                        continue;
                    }

                    if (n <= 250 && !(isRingCenter[u] && isRingCenter[v])) {
                        continue;
                    }

                    if (n > 250 && isRingCenter[u] && isRingCenter[v]) {
                        continue;
                    }
                    
                    let dx = positions[v].x - positions[u].x;
                    let dy = positions[v].y - positions[u].y;

                    if (dx === 0 || dy === 0) {
                        continue;
                    }
                    
                    let dSq = dx * dx + dy * dy;
                
                    if (dSq < 0.01) {
                        dx = 0.1 * Math.random() + 0.1;
                        dy = 0.1 * Math.random() + 0.1;

                        dSq = dx * dx + dy * dy;
                    }

                    let d = Math.sqrt(dSq);

                    if (d > maxDist) {
                        d = maxDist;
                        dSq = d * d;
                    }

                    let force = (dSq - k * k) / k;
                    let dOptimal = adjMatrix[u][v];
                    
                    force *= d / dOptimal;

                    let fx = force * dx / d;
                    let fy = force * dy / d;

                    if (!positioned[u]) {
                        forces[u].x += fx;
                        forces[u].y += fy;
                    }

                    if (!positioned[v]) {
                        forces[v].x -= fx;
                        forces[v].y -= fy;
                    }
                }
            }

            // Add the edge forces to the vertices
            for (let i = 0; i < edges.length; i++) {
                let index = edgeOffset + i;
                let force = forces[index];

                let a = edges[i][0];
                let b = edges[i][1];

                forces[a].x += force.x;
                forces[a].y += force.y;

                forces[b].x += force.x;
                forces[b].y += force.y;
            }

            // Move the vertex
            for (let u = 0; u < totalLength; u++) {
                if (positioned[u]) {
                    continue;
                }

                let dx = c * forces[u].x;
                let dy = c * forces[u].y;

                if (dx > maxMove) dx = maxMove;
                if (dx < -maxMove) dx = -maxMove;
                if (dy > maxMove) dy = maxMove;
                if (dy < -maxMove) dy = - maxMove;

                let dSq = dx * dx + dy * dy;
                

                positions[u].x += dx;
                positions[u].y += dy;
            }

            // Place the ring centers in the middle of the members
            if (n > 200 && ring.rings.length > 2) {
                for (let i = 0; i < ring.rings.length; i++) {
                    let r = ring.rings[i];
                    let center = new Vector2();

                    for (let j = 0; j < r.members.length; j++) {
                        let pos = positions[idToV[r.members[j]]];
                        center.x += pos.x;
                        center.y += pos.y;
                    }

                    center.x /= r.members.length;
                    center.y /= r.members.length;

                    positions[vertices.length + i] = center;
                }
            }          
        }

        for (let i = 0; i < totalLength; i++) {
            if (i < vertices.length) { 
                if (!positioned[i]) {
                    this.vertices[vToId[i]].setPositionFromVector(positions[i]);
                    this.vertices[vToId[i]].positioned = true;
                }
            } else if (i < vertices.length + ring.rings.length) {
                let index = i - vertices.length;
                ring.rings[index].center = positions[i];
            }
        }
        
        for (let u = 0; u < vertices.length; u++) {
            let vertex = this.vertices[vertices[u]];
            let parentVertex = this.vertices[vertex.parentVertexId];
            let neighbours = vertex.getNeighbours();
            
            for (let i = 0; i < neighbours.length; i++) {
                let currentVertex = this.vertices[neighbours[i]];
                
                if (currentVertex.positioned) {
                    continue;
                }

                center = this.getSubringCenter(ring, vertex);

                if (currentVertex.value.rings.length === 0) {
                    currentVertex.value.isConnectedToRing = true;
                }

                this.createNextBond(currentVertex, vertex, center);
            }
        }

        // This has to be called in order to position rings connected to this bridged ring
        this.createRing(ring, null, null, null, true);
    }

    /**
     * Gets the center of a ring contained within a bridged ring and containing a given vertex.
     *
     * @param {Ring} ring A bridged ring.
     * @param {Vertex} vertex A vertex.
     * @returns {Vector2} The center of the subring that contains the provided vertex.
     */
    getSubringCenter(ring, vertex) {
        // If there are multiple subrings associated with this ring, always
        // take the smaller one
        let size = Number.MAX_VALUE;
        let center = ring.center;

        for (let i = 0; i < ring.rings.length; i++) {
            let subring = ring.rings[i];
            for (let j = 0; j < subring.members.length; j++) {
                if (subring.members[j] === vertex.id) {
                    if (size > subring.members.length) {
                        center = subring.center;
                        size = subring.members.length;
                    }
                }
            }
        }



        return center;
    }

    /**
     * Draw the actual edges as bonds to the canvas.
     *
     * @param {boolean} debug A boolean indicating whether or not to draw debug helpers.
     */
    drawEdges(debug) {
        let that = this;
        
        for (let i = 0; i < this.edges.length; i++) {
            let edge = this.edges[i];
            let vertexA = this.vertices[edge.sourceId];
            let vertexB = this.vertices[edge.targetId];
            let elementA = vertexA.value.element;
            let elementB = vertexB.value.element;

            if ((!vertexA.value.isDrawn || !vertexB.value.isDrawn) && this.opts.atomVisualization === 'default') {
                continue;
            }

            let a = vertexA.position;
            let b = vertexB.position;
            let normals = this.getEdgeNormals(edge);

            // Create a point on each side of the line
            let sides = ArrayHelper.clone(normals);
            
            ArrayHelper.each(sides, function (v) {
                v.multiply(10);
                v.add(a)
            });

            if (edge.bondType === '=' || this.getRingbondType(vertexA, vertexB) === '=') {
                // Always draw double bonds inside the ring
                let inRing = this.areVerticesInSameRing(vertexA, vertexB);
                let s = this.chooseSide(vertexA, vertexB, sides);
                
                if (inRing) {
                    // Always draw double bonds inside a ring
                    // if the bond is shared by two rings, it is drawn in the larger
                    // problem: smaller ring is aromatic, bond is still drawn in larger -> fix this
                    let lcr = this.getLargestOrAromaticCommonRing(vertexA, vertexB);
                    let center = lcr.center;

                    ArrayHelper.each(normals, function (v) {
                        v.multiply(that.opts.bondSpacing);
                    });

                    // Choose the normal that is on the same side as the center
                    let line = null;
                    
                    if (center.sameSideAs(vertexA.position, vertexB.position, Vector2.add(a, normals[0]))) {
                        line = new Line(Vector2.add(a, normals[0]), Vector2.add(b, normals[0]), elementA, elementB);
                    } else {
                        line = new Line(Vector2.add(a, normals[1]), Vector2.add(b, normals[1]), elementA, elementB);
                    }

                    line.shorten(this.opts.bondLength - this.opts.shortBondLength);

                    // The shortened edge
                    this.canvasWrapper.drawLine(line);

                    // The normal edge
                    this.canvasWrapper.drawLine(new Line(a, b, elementA, elementB));
                } else if (edge.center) {
                    ArrayHelper.each(normals, function (v) {
                        v.multiply(that.opts.bondSpacing / 2.0)
                    });

                    let lineA = new Line(Vector2.add(a, normals[0]), Vector2.add(b, normals[0]), elementA, elementB);
                    let lineB = new Line(Vector2.add(a, normals[1]), Vector2.add(b, normals[1]), elementA, elementB);

                    lineA.shorten(this.opts.bondLength - this.opts.shortBondLength);
                    lineB.shorten(this.opts.bondLength - this.opts.shortBondLength);

                    this.canvasWrapper.drawLine(lineA);
                    this.canvasWrapper.drawLine(lineB);
                } else if (s.anCount == 0 && s.bnCount > 1 || s.bnCount == 0 && s.anCount > 1) {
                    // Both lines are the same length here
                    // Add the spacing to the edges (which are of unit length)
                    ArrayHelper.each(normals, function (v) {
                        v.multiply(that.opts.bondSpacing / 2)
                    });

                    let lineA = new Line(Vector2.add(a, normals[0]), Vector2.add(b, normals[0]), elementA, elementB);
                    let lineB = new Line(Vector2.add(a, normals[1]), Vector2.add(b, normals[1]), elementA, elementB);

                    this.canvasWrapper.drawLine(lineA);
                    this.canvasWrapper.drawLine(lineB);
                } else if (s.sideCount[0] > s.sideCount[1]) {
                    ArrayHelper.each(normals, function (v) {
                        v.multiply(that.opts.bondSpacing)
                    });

                    let line = new Line(Vector2.add(a, normals[0]), Vector2.add(b, normals[0]), elementA, elementB);
                    
                    line.shorten(this.opts.bondLength - this.opts.shortBondLength);
                    this.canvasWrapper.drawLine(line);
                    this.canvasWrapper.drawLine(new Line(a, b, elementA, elementB));
                } else if (s.sideCount[0] < s.sideCount[1]) {
                    ArrayHelper.each(normals, function (v) {
                        v.multiply(that.opts.bondSpacing)
                    });

                    let line = new Line(Vector2.add(a, normals[1]), Vector2.add(b, normals[1]), elementA, elementB);
                    
                    line.shorten(this.opts.bondLength - this.opts.shortBondLength);
                    this.canvasWrapper.drawLine(line);
                    this.canvasWrapper.drawLine(new Line(a, b, elementA, elementB));
                } else if (s.totalSideCount[0] > s.totalSideCount[1]) {
                    ArrayHelper.each(normals, function (v) {
                        v.multiply(that.opts.bondSpacing)
                    });

                    let line = new Line(Vector2.add(a, normals[0]), Vector2.add(b, normals[0]), elementA, elementB);
                    
                    line.shorten(this.opts.bondLength - this.opts.shortBondLength);
                    this.canvasWrapper.drawLine(line);
                    this.canvasWrapper.drawLine(new Line(a, b, elementA, elementB));
                } else if (s.totalSideCount[0] <= s.totalSideCount[1]) {
                    ArrayHelper.each(normals, function (v) {
                        v.multiply(that.opts.bondSpacing)
                    });

                    let line = new Line(Vector2.add(a, normals[1]), Vector2.add(b, normals[1]), elementA, elementB);
                    
                    line.shorten(this.opts.bondLength - this.opts.shortBondLength);
                    this.canvasWrapper.drawLine(line);
                    this.canvasWrapper.drawLine(new Line(a, b, elementA, elementB));
                } else {

                }
            } 
            else if(edge.bondType === '#') {
                ArrayHelper.each(normals, function (v) {
                    v.multiply(that.opts.bondSpacing / 1.5)
                });

                let lineA = new Line(Vector2.add(a, normals[0]), Vector2.add(b, normals[0]), elementA, elementB);
                let lineB = new Line(Vector2.add(a, normals[1]), Vector2.add(b, normals[1]), elementA, elementB);

                lineA.shorten(this.opts.bondLength - this.opts.shortBondLength);
                lineB.shorten(this.opts.bondLength - this.opts.shortBondLength);

                this.canvasWrapper.drawLine(lineA);
                this.canvasWrapper.drawLine(lineB);

                this.canvasWrapper.drawLine(new Line(a, b, elementA, elementB));
            } else {
                let isChiralCenterA = vertexA.value.bracket && vertexA.value.bracket.chirality;
                let isChiralCenterB = vertexB.value.bracket && vertexB.value.bracket.chirality;
                
                if (edge.chiral === 'up') {
                    this.canvasWrapper.drawWedge(new Line(a, b, elementA, elementB, isChiralCenterA, isChiralCenterB));
                } else if (edge.chiral === 'down') {
                    this.canvasWrapper.drawDashedWedge(new Line(a, b, elementA, elementB, isChiralCenterA, isChiralCenterB));
                } else {
                    this.canvasWrapper.drawLine(new Line(a, b, elementA, elementB, isChiralCenterA, isChiralCenterB));
                }
            }

            if (debug) {
                let midpoint = Vector2.midpoint(a, b);
                this.canvasWrapper.drawDebugText(midpoint.x, midpoint.y, 'e: ' + i);
            }
        }

        // Draw ring for benzenes
        for (let i = 0; i < this.rings.length; i++) {
            let ring = this.rings[i];

            if (ring.isAromatic(this.vertices)) {
                this.canvasWrapper.drawAromaticityRing(ring);
            }
        }
    }

    /**
     * Draws the vertices representing atoms to the canvas.
     *
     * @param {boolean} debug A boolean indicating whether or not to draw debug messages to the canvas.
     */
    drawVertices(debug) {
        for (let i = 0; i < this.vertices.length; i++) {
            let vertex = this.vertices[i];
            let atom = vertex.value;

            let charge = 0;
            let isotope = 0;
            let bondCount = this.getBondCount(vertex);
            let element = atom.element.length === 1 ? atom.element.toUpperCase() : atom.element;
            let hydrogens = this.maxBonds[element] - bondCount;
            let dir = vertex.getTextDirection(this.vertices);
            let isTerminal = this.opts.terminalCarbons || element !== 'C' || atom.hasAttachedPseudoElements ? vertex.isTerminal() : false;
            let isCarbon = atom.element.toLowerCase() === 'c';

            if (atom.bracket) {
                hydrogens = atom.bracket.hcount;
                charge = atom.bracket.charge;
                isotope = atom.bracket.isotope;
            }

            if ((!isCarbon || atom.explicit || isTerminal || atom.hasAttachedPseudoElements) && atom.isDrawn) {
                if (this.opts.atomVisualization === 'default') {
                    this.canvasWrapper.drawText(vertex.position.x, vertex.position.y,
                            element, hydrogens, dir, isTerminal, charge, isotope, atom.getAttachedPseudoElements());
                } else if (this.opts.atomVisualization === 'balls') {
                    this.canvasWrapper.drawBall(vertex.position.x, vertex.position.y,
                            element);
                }
            }

            if (debug) {
                let value = 'v: ' + vertex.id + ' ' + ArrayHelper.print(atom.ringbonds);
                this.canvasWrapper.drawDebugText(vertex.position.x, vertex.position.y, value);
            }
        }

        // Draw the ring centers for debug purposes
        if (this.opts.debug) {
            for (let i = 0; i < this.rings.length; i++) {
                let center = this.rings[i].center;
                this.canvasWrapper.drawDebugPoint(center.x, center.y, 
                        'r: ' + this.rings[i].id);
            }
        }   
    }

    /**
     * Position the vertices according to their bonds and properties.
     *
     */
    position() {
        let startVertex = this.vertices[0];

        // If there is a bridged ring, alwas start with the bridged ring
        for (let i = 0; i < this.rings.length; i++) {
            if (this.rings[i].isBridged) {
                for (let j = 0; j < this.rings[i].members.length; j++) {
                    startVertex = this.vertices[this.rings[i].members[j]];
                    
                    if (startVertex.value.originalRings.length === 1) {
                        break;
                    }
                }
            }
        }

        this.createNextBond(startVertex);

        // Atoms bonded to the same ring atom
        this.resolvePrimaryOverlaps();
    }

    /**
     * Reset the positions of rings and vertices. The previous positions will be backed up.
     *
     */
    clearPositions() {
        this.vertexPositionsBackup = [];
        this.ringPositionsBackup = [];

        for (let i = 0; i < this.vertices.length; i++) {
            let vertex = this.vertices[i];
            this.vertexPositionsBackup.push(vertex.position.clone());
            vertex.positioned = false;
            vertex.setPositionFromVector(new Vector2());
        }

        for (let i = 0; i < this.rings.length; i++) {
            var ring = this.rings[i];
            this.ringPositionsBackup.push(ring.center.clone());
            ring.positioned = false;
            ring.center = new Vector2();
        }
    }

    /**
     * Restore the positions backed up during the last clearPositions() call.
     *
     */
    restorePositions() {        
        for (let i = 0; i < this.vertexPositionsBackup.length; i++) {
            this.vertices[i].setPositionFromVector(this.vertexPositionsBackup[i]);
            this.vertices[i].positioned = true;
        }

        for (let i = 0; i < this.ringPositionsBackup.length; i++) {
            this.rings[i].center = this.ringPositionsBackup[i];
            this.rings[i].positioned = true;
        }
    }

    /**
     * Stores the current information associated with rings.
     * 
     */
    backupRingInformation() {
        this.originalRings = [];
        this.originalRingConnections = [];

        for (let i = 0; i < this.rings.length; i++) {
            this.originalRings.push(this.rings[i]);
        }

        for (let i = 0; i < this.ringConnections.length; i++) {
            this.originalRingConnections.push(this.ringConnections[i]);
        }

        for (let i = 0; i < this.vertices.length; i++) {
            this.vertices[i].value.backupRings();
        }
    }

    /**
     * Restores the most recently backed up information associated with rings.
     * 
     */
    restoreRingInformation() {
        // Get the subring centers from the bridged rings
        let bridgedRings = this.getBridgedRings();

        this.rings = [];
        this.ringConnections = [];

        for (let i = 0; i < bridgedRings.length; i++) {
            let bridgedRing = bridgedRings[i];

            for (let j = 0; j < bridgedRing.rings.length; j++) {
                let ring = bridgedRing.rings[j];
                this.originalRings[ring.id].center = ring.center;
            }
        }

        for (let i = 0; i < this.originalRings.length; i++) {
            this.rings.push(this.originalRings[i]);
        }

        for (let i = 0; i < this.originalRingConnections.length; i++) {
            this.ringConnections.push(this.originalRingConnections[i]);
        }

        for (let i = 0; i < this.vertices.length; i++) {
            this.vertices[i].value.restoreRings();
        }
    }

    // TODO: This needs some cleaning up

    /**
     * Creates a new ring, that is, positiones all the vertices inside a ring.
     *
     * @param {Ring} ring The ring to position.
     * @param {Vector2|null} [center=null] The center of the ring to be created.
     * @param {Vertex|null} [startVertex=null] The first vertex to be positioned inside the ring.
     * @param {Vertex|null} [previousVertex=null] The last vertex that was positioned.
     * @param {boolean} [previousVertex=false] A boolean indicating whether or not this ring was force positioned already - this is needed after force layouting a ring, in order to draw rings connected to it.
     */
    createRing(ring, center = null, startVertex = null, previousVertex = null, forcePositioned = false) {
        if (ring.positioned && !forcePositioned) {
            return;
        }

        center = center ? center : new Vector2(0, 0);
        
        let orderedNeighbours = ring.getOrderedNeighbours(this.ringConnections);
        let startingAngle = startVertex ? Vector2.subtract(startVertex.position, center).angle() : 0;

        let radius = MathHelper.polyCircumradius(this.opts.bondLength, ring.getSize());
        let angle = MathHelper.centralAngle(ring.getSize());

        ring.centralAngle = angle;
        
        let a = startingAngle;
        let that = this;

        if (!forcePositioned) {
            ring.eachMember(this.vertices, function (v) {
                let vertex = that.vertices[v];

                if (!vertex.positioned) {
                    vertex.setPosition(center.x + Math.cos(a) * radius, center.y + Math.sin(a) * radius);
                }

                a += angle;
                
                if(!ring.isBridged || ring.rings.length < 3) {
                    vertex.positioned = true;
                }
            }, (startVertex) ? startVertex.id : null, (previousVertex) ? previousVertex.id : null);

            // If the ring is bridged, then draw the vertices inside the ring
            // using a force based approach
            if (ring.isBridged) {
                let allVertices = ArrayHelper.merge(ring.members, ring.insiders);

                this.forceLayout(allVertices, center, startVertex.id, ring);
            }

            // Anchor the ring to one of it's members, so that the ring center will always
            // be tied to a single vertex when doing repositionings
            this.vertices[ring.members[0]].value.addAnchoredRing(ring.id);

            ring.positioned = true;
            ring.center = center;
        }

        // Draw neighbours in decreasing order of connectivity
        for (let i = 0; i < orderedNeighbours.length; i++) {
            let neighbour = this.getRing(orderedNeighbours[i].neighbour);
            
            if (neighbour.positioned) {
                continue;
            }

            let vertices = RingConnection.getVertices(this.ringConnections, ring.id, neighbour.id);
            
            if (vertices.length == 2) {
                // This ring is a fused ring
                ring.isFused = true;
                neighbour.isFused = true;

                let vertexA = this.vertices[vertices[0]];
                let vertexB = this.vertices[vertices[1]];

                // Get middle between vertex A and B
                let midpoint = Vector2.midpoint(vertexA.position, vertexB.position);

                // Get the normals to the line between A and B
                let normals = Vector2.normals(vertexA.position, vertexB.position);

                // Normalize the normals
                ArrayHelper.each(normals, function (v) {
                    v.normalize()
                });

                // Set length from middle of side to center (the apothem)
                let r = MathHelper.polyCircumradius(this.opts.bondLength, neighbour.getSize());
                let apothem = MathHelper.apothem(r, neighbour.getSize());
                
                ArrayHelper.each(normals, function (v) {
                    v.multiply(apothem)
                });

                // Move normals to the middle of the line between a and b
                ArrayHelper.each(normals, function (v) {
                    v.add(midpoint)
                });

                // Check if the center of the next ring lies within another ring and
                // select the normal accordingly
                let nextCenter = normals[0];
                
                if (this.isPointInRing(nextCenter)) {
                    nextCenter = normals[1];
                }

                // Get the vertex (A or B) which is in clock-wise direction of the other
                let posA = Vector2.subtract(vertexA.position, nextCenter);
                let posB = Vector2.subtract(vertexB.position, nextCenter);

                if (posA.clockwise(posB) === -1) {
                    this.createRing(neighbour, nextCenter, vertexA, vertexB);
                } else {
                    this.createRing(neighbour, nextCenter, vertexB, vertexA);
                }
            } else if (vertices.length == 1) {
                // This ring is a spiro
                ring.isSpiro = true;
                neighbour.isSpiro = true;

                let vertexA = this.vertices[vertices[0]];
                
                // Get the vector pointing from the shared vertex to the new center
                let nextCenter = Vector2.subtract(center, vertexA.position);
                nextCenter.invert();
                nextCenter.normalize();

                // Get the distance from the vertex to the center
                let r = MathHelper.polyCircumradius(this.opts.bondLength, neighbour.getSize());
                nextCenter.multiply(r);
                nextCenter.add(vertexA.position);
                this.createRing(neighbour, nextCenter, vertexA);
            }
        }

        // Next, draw atoms that are not part of a ring that are directly attached to this ring
        for (let i = 0; i < ring.members.length; i++) {
            let ringMember = this.vertices[ring.members[i]];
            let ringMemberNeighbours = ringMember.getNeighbours();

            // If there are multiple, the ovlerap will be resolved in the appropriate step
            for (let j = 0; j < ringMemberNeighbours.length; j++) {
                if (ring.thisOrNeighboursContain(this.rings, ringMemberNeighbours[j])) {
                    continue;
                }
                
                let v = this.vertices[ringMemberNeighbours[j]];
                v.value.isConnectedToRing = true;

                this.createNextBond(v, ringMember, ring.center);
            }
        }
    }

    /**
     * Rotate an entire subtree by an angle around a center.
     *
     * @param {number} vertexId A vertex id (the root of the sub-tree).
     * @param {number} parentVertexId A vertex id in the previous direction of the subtree that is to rotate.
     * @param {number} angle An angle in randians.
     * @param {Vector2} center The rotational center.
     */
    rotateSubtree(vertexId, parentVertexId, angle, center) {
        let that = this;
        
        this.traverseTree(vertexId, parentVertexId, function (vertex) {
            vertex.position.rotateAround(angle, center);

            for (let i = 0; i < vertex.value.anchoredRings.length; i++) {
                let ring = that.rings[vertex.value.anchoredRings[i]];

                if (ring) {
                    ring.center.rotateAround(angle, center);
                }
            }
        });
    }

    /**
     * Gets the overlap score of a subtree.
     *
     * @param {number} vertexId A vertex id (the root of the sub-tree).
     * @param {number} parentVertexId A vertex id in the previous direction of the subtree.
     * @param {Array} vertexOverlapScores An array containing the vertex overlap scores indexed by vertex id.
     * @returns {number} The overlap score of the subtree.
     */
    getSubtreeOverlapScore(vertexId, parentVertexId, vertexOverlapScores) {
        let that = this;
        let score = 0;
        let center = new Vector2();
        
        this.traverseTree(vertexId, parentVertexId, function (vertex) {
            let s = vertexOverlapScores[vertex.id];
            score += s;

            let position = that.vertices[vertex.id].position.clone();
            position.multiply(s)
            center.add(position);
        });

        center.divide(score);

        return { value: score, center: center };
    }

    /**
     * Returns the current (positioned vertices so far) center of mass.
     * 
     * @returns {Vector2} The current center of mass.
     */
    getCurrentCenterOfMass() {
        let total = new Vector2();
        let count = 0;
        
        for (let i = 0; i < this.vertices.length; i++) {
            let vertex = this.vertices[i];

            if (vertex.positioned) {
                total.add(vertex.position);
                count++;
            }
        }

        return total.divide(count);
    }

    /**
     * Returns the current (positioned vertices so far) center of mass in the neighbourhood of a given position.
     *
     * @param {Vector2} vec The point at which to look for neighbours.
     * @param {number} [r=currentBondLength*2.0] The radius of vertices to include.
     * @returns {Vector2} The current center of mass.
     */
    getCurrentCenterOfMassInNeigbourhood(vec, r = this.opts.bondLength * 2.0) {
        let total = new Vector2();
        let count = 0;
        let rSq = r * r;
        
        for (let i = 0; i < this.vertices.length; i++) {
            let vertex = this.vertices[i];

            if (vertex.positioned && vec.distanceSq(vertex.position) < rSq) {
                total.add(vertex.position);
                count++;
            }
        }

        return total.divide(count);
    }

    /**
     * Resolve primary (exact) overlaps, such as two vertices that are connected to the same ring vertex.
     *
     */
    resolvePrimaryOverlaps() {
        let overlaps = [];
        let sharedSideChains = []; // side chains attached to an atom shared by two rings
        let done = new Array(this.vertices.length);

        for (let i = 0; i < this.rings.length; i++) {
            let ring = this.rings[i];
            
            for (let j = 0; j < ring.members.length; j++) {
                let vertex = this.vertices[ring.members[j]];

                if (done[vertex.id]) {
                    continue;
                }

                done[vertex.id] = true;

                // Look for rings where there are atoms with two bonds outside the ring (overlaps)
                if (vertex.getNeighbourCount() > 2) {
                    let rings = [];
                    
                    for (var k = 0; k < vertex.value.rings.length; k++) {
                        rings.push(vertex.value.rings[k]);
                    }

                    overlaps.push({
                        common: vertex,
                        rings: rings,
                        vertices: this.getNonRingNeighbours(vertex.id)
                    });
                }
            }
        }

        for (let i = 0; i < sharedSideChains.length; i++) {
            let chain = sharedSideChains[i];
            let angle = -chain.vertex.position.getRotateToAngle(chain.other.position, chain.common.position);
            this.rotateSubtree(chain.vertex.id, chain.common.id, angle + Math.PI, chain.common.position);
        }

        for (let i = 0; i < overlaps.length; i++) {
            let overlap = overlaps[i];

            if (overlap.vertices.length === 1) {
                let a = overlap.vertices[0];
                
                if (a.getNeighbourCount() === 1) {
                    a.flippable = true;
                    a.flipCenter = overlap.common.id;

                    for (let j = 0; j < overlap.rings.length; j++) {
                        a.flipRings.push(overlap.rings[j]);
                    }
                }

                // If the vertex comes out of two rings, it has to be rotated to point straight away (angles between it and both rings the same)
                if (overlap.rings.length === 2) {
                    let neighbours = overlap.common.getNeighbours();
                    let positions = [];

                    for (let j = 0; j < neighbours.length; j++) {
                        let vertex = this.vertices[neighbours[j]];

                        if (!this.isRingConnection(vertex.id, overlap.common.id) && vertex.id !== a.id) {
                            positions.push(vertex.position);
                        }
                    }

                    let midpoint = Vector2.midpoint(positions[0], positions[1]);
                    let angle = a.position.getRotateToAngle(midpoint, overlap.common.position);

                    angle *= a.position.clockwise(midpoint);
                    this.rotateSubtree(a.id, overlap.common.id, angle, overlap.common.position);
                }
            } else if (overlap.vertices.length === 2) {
                let angle = (2 * Math.PI - this.getRing(overlap.rings[0]).getAngle()) / 6.0;
                let a = overlap.vertices[0];
                let b = overlap.vertices[1];
                
                a.backAngle -= angle;
                b.backAngle += angle;
                
                this.rotateSubtree(a.id, overlap.common.id, angle, overlap.common.position);
                this.rotateSubtree(b.id, overlap.common.id, -angle, overlap.common.position);

                if (a.getNeighbourCount() === 1) {
                    a.flippable = true;
                    a.flipCenter = overlap.common.id;
                    a.flipNeighbour = b.id;
                    
                    for (let j = 0; j < overlap.rings.length; j++) {
                        a.flipRings.push(overlap.rings[j]);
                    }
                }
                if (b.getNeighbourCount() === 1) {
                    b.flippable = true;
                    b.flipCenter = overlap.common.id;
                    b.flipNeighbour = a.id;
                    
                    for (let j = 0; j < overlap.rings.length; j++) {
                        b.flipRings.push(overlap.rings[j]);
                    }
                }
            }
        }
    }

    /**
     * Resolve secondary overlaps. Those overlaps are due to the structure turning back on itself.
     *
     * @param {array} scores An array of objects sorted descending by score. An object is in the form of { id: 0, score: 22 }.
     */
    resolveSecondaryOverlaps(scores) {
        for (let i = 0; i < scores.length; i++) {
            if (scores[i].score > this.opts.bondLength / (4.0 * this.opts.bondLength)) {
                let vertex = this.vertices[scores[i].id];

                if (vertex.isTerminal()) {
                    let closest = this.getClosestVertex(vertex);
                    
                    if (closest) {
                        // If one of the vertices is the first one, the previous vertex is not the central vertex but the dummy
                        // so take the next rather than the previous, which is vertex 1
                        let closestPosition = null;
                        
                        if (closest.isTerminal()) {
                            closestPosition = closest.id === 0 ? this.vertices[1].position : closest.previousPosition
                        } else {
                            closestPosition = closest.id === 0 ? this.vertices[1].position : closest.position
                        }

                        let vertexPreviousPosition = vertex.id === 0 ? this.vertices[1].position : vertex.previousPosition;

                        vertex.position.rotateAwayFrom(closestPosition, vertexPreviousPosition, MathHelper.toRad(20));
                    }
                }

                if (vertex.flippable) {
                    let a = vertex.flipRings[0] ? this.rings[vertex.flipRings[0]] : null;
                    let b = vertex.flipRings[1] ? this.rings[vertex.flipRings[1]] : null;
                    let flipCenter = this.vertices[vertex.flipCenter].position;

                    // Make a always the bigger ring than b
                    if (a && b) {
                        let tmp = (a.members.length > b.members.length) ? a : b;
                        
                        b = (a.members.length < b.members.length) ? a : b;
                        a = tmp;
                    }
                    
                    if (this.opts.allowFlips) {
                        if (a && a.allowsFlip()) {
                            vertex.position.rotateTo(a.center, flipCenter);
                            a.setFlipped();
                            
                            if (vertex.flipNeighbour !== null) {
                                // It's better to not straighten the other one, since it will possibly overlap
                                // var flipNeighbour = this.vertices[vertex.flipNeighbour];
                                // flipNeighbour.position.rotate(flipNeighbour.backAngle);
                            }
                        } else if (b && b.allowsFlip()) {
                            vertex.position.rotateTo(b.center, flipCenter);
                            b.setFlipped();
                            
                            if (vertex.flipNeighbour !== null) {
                                // It's better to not straighten the other one, since it will possibly overlap
                                // var flipNeighbour = this.vertices[vertex.flipNeighbour];
                                // flipNeighbour.position.rotate(flipNeighbour.backAngle);
                            }
                        }
                    } else {

                    }

                    // Only do a refresh of the remaining!
                    // recalculate scores (how expensive?)
                    // scores = this.getOverlapScore().scores;
                }
            }
        }
    }

    /**
     * Positiones the next vertex thus creating a bond.
     *
     * @param {Vertex} vertex A vertex.
     * @param {Vertex} previousVertex The previous vertex which has been positioned.
     * @param {ring|number} ringOrAngle Either a ring or a number. If the vertex is connected to a ring, it is positioned based on the ring center and thus the ring is supplied. If the vertex is not in a ring, an angle (in radians) is supplied.
     * @param {number} dir Either 1 or -1 to break ties (if no angle can be elucidated.
     */
    createNextBond(vertex, previousVertex, ringOrAngle, dir) {
        if (vertex.positioned) {
            return;
        }

        // If the current node is the member of one ring, then point straight away
        // from the center of the ring. However, if the current node is a member of
        // two rings, point away from the middle of the centers of the two rings
        if (!previousVertex) {
            // Here, ringOrAngle is always an angle

            // Add a (dummy) previous position if there is no previous vertex defined
            // Since the first vertex is at (0, 0), create a vector at (bondLength, 0)
            // and rotate it by 90°
            let dummy = new Vector2(this.opts.bondLength, 0);
            dummy.rotate(MathHelper.toRad(-120));

            vertex.previousPosition = dummy;
            vertex.setPosition(this.opts.bondLength, 0);
            vertex.angle = MathHelper.toRad(-120);
            vertex.globalAngle = vertex.angle;
            vertex.positioned = true;
        } else if (previousVertex.value.rings.length === 0 && !vertex.value.isBridge && !previousVertex.value.isBridge) {
            // Here, ringOrAngle is always an angle

            // If the previous vertex was not part of a ring, draw a bond based
            // on the global angle of the previous bond
            let v = new Vector2(this.opts.bondLength, 0);
            v.rotate(ringOrAngle);
            v.add(previousVertex.position);

            vertex.globalAngle = ringOrAngle;
            vertex.setPositionFromVector(v);

            vertex.previousPosition = previousVertex.position;
            vertex.positioned = true;
        } else if (previousVertex.value.isBridgeNode && vertex.value.isBridge) {
            // If the previous atom is in a bridged ring and this one is inside the ring
            let pos = Vector2.subtract(ringOrAngle, previousVertex.position);
            pos.normalize();

            // Unlike with the ring, do not multiply with radius but with bond length
            pos.multiply(this.opts.bondLength);
            vertex.position.add(previousVertex.position);
            vertex.position.add(pos);

            vertex.previousPosition = previousVertex.position;
            vertex.positioned = true;
        } else if (vertex.value.isBridge) {
            // The previous atom is in a bridged ring and this one is in it as well
            let v = new Vector2(this.opts.bondLength, 0);
            v.rotate(ringOrAngle);
            v.add(previousVertex.position);

            vertex.globalAngle = ringOrAngle;
            vertex.setPositionFromVector(v);
            vertex.previousPosition = previousVertex.position;
            vertex.positioned = true;
        } else if (previousVertex.value.rings.length === 1 || previousVertex.value.isBridge) {
            // Here, ringOrAngle is always a ring (THIS IS CURRENTLY NOT TRUE - WHY?)
            // Use the same approach as with rings that are connected at one vertex
            // and draw the atom in the opposite direction of the center.
            let pos = Vector2.subtract(ringOrAngle, previousVertex.position);

            pos.invert();
            pos.normalize();
            // Unlike with the ring, do not multiply with radius but with bond length
            pos.multiply(this.opts.bondLength);
            vertex.position.add(previousVertex.position);
            vertex.position.add(pos);
            vertex.previousPosition = previousVertex.position;
            vertex.positioned = true;
        } else if (previousVertex.value.rings.length == 2) {
            // Here, ringOrAngle is always a ring
            let ringA = this.getRing(previousVertex.value.rings[0]);
            let ringB = this.getRing(previousVertex.value.rings[1]);

            // Project the current vertex onto the vector between the two centers to
            // get the direction
            let a = Vector2.subtract(ringB.center, ringA.center);
            let b = Vector2.subtract(previousVertex.position, ringA.center);
            let s = Vector2.scalarProjection(b, a);
            
            a.normalize();
            a.multiply(s);
            a.add(ringA.center);

            let pos = Vector2.subtract(a, previousVertex.position);
            pos.invert();
            pos.normalize();
            pos.multiply(this.opts.bondLength);
            
            vertex.position.add(previousVertex.position);
            vertex.position.add(pos);

            vertex.previousPosition = previousVertex.position;
            vertex.positioned = true;
        }

        // Go to next vertex
        // If two rings are connected by a bond ...
        if (vertex.value.rings.length > 0) {
            let nextRing = this.getRing(vertex.value.rings[0]);
            let nextCenter = Vector2.subtract(vertex.previousPosition, vertex.position);

            nextCenter.invert();
            nextCenter.normalize();

            let r = MathHelper.polyCircumradius(this.opts.bondLength, nextRing.getSize());

            nextCenter.multiply(r);
            nextCenter.add(vertex.position);

            this.createRing(nextRing, nextCenter, vertex);
        } else {
            // Draw the non-ring vertices connected to this one        
            let neighbours = vertex.getNeighbours();
            
            if (previousVertex) {
                neighbours = ArrayHelper.remove(neighbours, previousVertex.id);
            }

            let angle = vertex.getAngle();

            if (neighbours.length === 1) {
                let nextVertex = this.vertices[neighbours[0]];

                // Make a single chain always cis except when there's a tribble bond
                if((vertex.value.bondType === '#' || (previousVertex && previousVertex.value.bondType === '#')) ||
                    vertex.value.bondType === '=' && previousVertex && previousVertex.value.bondType === '=') {
                    vertex.value.explicit = true;
                    
                    if (previousVertex) {
                        let straightEdge1 = this.getEdge(vertex.id, previousVertex.id);
                        straightEdge1.center = true;
                    }
                    
                    let straightEdge2 = this.getEdge(vertex.id, nextVertex.id);    
                    straightEdge2.center = true;

                    nextVertex.globalAngle = angle;
                    nextVertex.angle = 0.0;
                    this.createNextBond(nextVertex, vertex, nextVertex.globalAngle, -dir);
                } else if (previousVertex && previousVertex.value.rings.length > 0) {
                    // If coming out of a ring, always draw away from the center of mass
                    let proposedAngleA = MathHelper.toRad(60);
                    let proposedAngleB = -proposedAngleA;

                    let proposedVectorA = new Vector2(this.opts.bondLength, 0);
                    let proposedVectorB = new Vector2(this.opts.bondLength, 0);
                    
                    proposedVectorA.rotate(proposedAngleA).add(vertex.position);
                    proposedVectorB.rotate(proposedAngleB).add(vertex.position);

                    // let centerOfMass = this.getCurrentCenterOfMassInNeigbourhood(vertex.position, 100);
                    let centerOfMass = this.getCurrentCenterOfMass();
                    let distanceA = proposedVectorA.distance(centerOfMass);
                    let distanceB = proposedVectorB.distance(centerOfMass);

                    nextVertex.angle = distanceA < distanceB ? proposedAngleB : proposedAngleA;
                    
                    if (nextVertex.angle > 0) {
                        dir = -1;
                    } else {
                        dir = 1;
                    }
                    
                    nextVertex.globalAngle = angle + nextVertex.angle;
                    this.createNextBond(nextVertex, vertex, nextVertex.globalAngle, dir);
                } else {
                    if (!dir) {
                        let proposedAngleA = MathHelper.toRad(60);
                        let proposedAngleB = -proposedAngleA;

                        let proposedVectorA = new Vector2(this.opts.bondLength, 0);
                        let proposedVectorB = new Vector2(this.opts.bondLength, 0);
                        
                        proposedVectorA.rotate(proposedAngleA).add(vertex.position);
                        proposedVectorB.rotate(proposedAngleB).add(vertex.position);

                        let centerOfMass = this.getCurrentCenterOfMass();
                        let distanceA = proposedVectorA.distance(centerOfMass);
                        let distanceB = proposedVectorB.distance(centerOfMass);

                        nextVertex.angle = distanceA < distanceB ? proposedAngleB : proposedAngleA;
                        
                        if (nextVertex.angle > 0) {
                            dir = -1;
                        } else {
                            dir = 1;
                        }
                    } else {
                        nextVertex.angle = MathHelper.toRad(60) * dir;
                        dir = -dir;
                    }
                    
                    nextVertex.globalAngle = angle + nextVertex.angle;
                    this.createNextBond(nextVertex, vertex, nextVertex.globalAngle, dir);
                }
            } else if (neighbours.length === 2) {
                // Check for the longer subtree - always go with cis for the longer subtree
                let subTreeDepthA = this.getTreeDepth(neighbours[0], vertex.id);
                let subTreeDepthB = this.getTreeDepth(neighbours[1], vertex.id);
                
                let cis = 0;
                let trans = 1;

                if (subTreeDepthA > subTreeDepthB) {
                    cis = 1;
                    trans = 0;
                }

                if (vertex.position.clockwise(vertex.previousPosition) === 1) {
                    let cisVertex = this.vertices[neighbours[cis]];
                    let transVertex = this.vertices[neighbours[trans]];

                    transVertex.angle = MathHelper.toRad(60);
                    cisVertex.angle = -MathHelper.toRad(60);

                    transVertex.globalAngle = angle + transVertex.angle;
                    cisVertex.globalAngle = angle + cisVertex.angle;

                    this.createNextBond(transVertex, vertex, transVertex.globalAngle, -dir);
                    this.createNextBond(cisVertex, vertex, cisVertex.globalAngle, -dir);
                } else {
                    let cisVertex = this.vertices[neighbours[cis]];
                    let transVertex = this.vertices[neighbours[trans]];

                    transVertex.angle = -MathHelper.toRad(60);
                    cisVertex.angle = MathHelper.toRad(60);

                    transVertex.globalAngle = angle + transVertex.angle;
                    cisVertex.globalAngle = angle + cisVertex.angle;

                    this.createNextBond(cisVertex, vertex, cisVertex.globalAngle, -dir);
                    this.createNextBond(transVertex, vertex, transVertex.globalAngle, -dir);
                }
            } else if (neighbours.length === 3) {
                // The vertex with the longest sub-tree should always go straight
                let d1 = this.getTreeDepth(neighbours[0], vertex.id);
                let d2 = this.getTreeDepth(neighbours[1], vertex.id);
                let d3 = this.getTreeDepth(neighbours[2], vertex.id);
                
                let s = this.vertices[neighbours[0]];
                let l = this.vertices[neighbours[1]];
                let r = this.vertices[neighbours[2]];

                if(d2 > d1 && d2 > d3) {
                    s = this.vertices[neighbours[1]];
                    l = this.vertices[neighbours[0]];
                    r = this.vertices[neighbours[2]];
                }
                else if(d3 > d1 && d3 > d2) {
                    s = this.vertices[neighbours[2]];
                    l = this.vertices[neighbours[0]];
                    r = this.vertices[neighbours[1]];
                }
                
                if (this.getTreeDepth(l.id, vertex.id) === 1 && 
                    this.getTreeDepth(r.id, vertex.id) === 1 &&
                    this.getTreeDepth(s.id, vertex.id) > 1) { 
                    
                    if (!dir) {
                        let proposedAngleA = MathHelper.toRad(60);
                        let proposedAngleB = -proposedAngleA;

                        let proposedVectorA = new Vector2(this.opts.bondLength, 0);
                        let proposedVectorB = new Vector2(this.opts.bondLength, 0);
                        
                        proposedVectorA.rotate(proposedAngleA).add(vertex.position);
                        proposedVectorB.rotate(proposedAngleB).add(vertex.position);

                        // let centerOfMass = this.getCurrentCenterOfMassInNeigbourhood(vertex.position, 100);
                        let centerOfMass = this.getCurrentCenterOfMass();
                        let distanceA = proposedVectorA.distance(centerOfMass);
                        let distanceB = proposedVectorB.distance(centerOfMass);

                        s.angle = distanceA < distanceB ? proposedAngleB : proposedAngleA;
                        
                        if (s.angle > 0) {
                            dir = -1;
                        } else {
                            dir = 1;
                        }
                    } else {
                        s.angle = MathHelper.toRad(60) * dir;
                        dir = -dir;
                    }

                    s.globalAngle = angle + s.angle;
                
                    this.createNextBond(s, vertex, s.globalAngle, -dir);

                    // If it's chiral, the order changes - for anticlockwise, switch the draw order around
                    // to keep the drawing the same
                    if (vertex.value.bracket && vertex.value.bracket.chirality === '@@') {
                        r.angle = MathHelper.toRad(30) * dir;
                        l.angle = MathHelper.toRad(90) * dir;

                        r.globalAngle = angle + r.angle;
                        l.globalAngle = angle + l.angle;

                        this.createNextBond(r, vertex, r.globalAngle);
                        this.createNextBond(l, vertex, l.globalAngle); 
                    } else {
                        l.angle = MathHelper.toRad(30) * dir;
                        r.angle = MathHelper.toRad(90) * dir;

                        l.globalAngle = angle + l.angle;
                        r.globalAngle = angle + r.angle;

                        this.createNextBond(l, vertex, l.globalAngle);
                        this.createNextBond(r, vertex, r.globalAngle);
                    }
                } else {
                    s.angle = 0.0;
                    l.angle = MathHelper.toRad(90);
                    r.angle = -MathHelper.toRad(90);

                    s.globalAngle = angle + s.angle;
                    l.globalAngle = angle + l.angle;
                    r.globalAngle = angle + r.angle;

                    this.createNextBond(s, vertex, s.globalAngle);
                    this.createNextBond(l, vertex, l.globalAngle);
                    this.createNextBond(r, vertex, r.globalAngle);
                }
            } else if (neighbours.length === 4) {
                // The vertex with the longest sub-tree should always go to the reflected opposide direction
                let d1 = this.getTreeDepth(neighbours[0], vertex.id);
                let d2 = this.getTreeDepth(neighbours[1], vertex.id);
                let d3 = this.getTreeDepth(neighbours[2], vertex.id);
                let d4 = this.getTreeDepth(neighbours[3], vertex.id);

                let w = this.vertices[neighbours[0]];
                let x = this.vertices[neighbours[1]];
                let y = this.vertices[neighbours[2]];
                let z = this.vertices[neighbours[3]];

                if(d2 > d1 && d2 > d3 && d2 > d4) {
                    w = this.vertices[neighbours[1]];
                    x = this.vertices[neighbours[0]];
                    y = this.vertices[neighbours[2]];
                    z = this.vertices[neighbours[3]];
                }
                else if(d3 > d1 && d3 > d2 && d3 > d4) {
                    w = this.vertices[neighbours[2]];
                    x = this.vertices[neighbours[0]];
                    y = this.vertices[neighbours[1]];
                    z = this.vertices[neighbours[3]];
                }
                else if(d4 > d1 && d4 > d2 && d4 > d3) {
                    w = this.vertices[neighbours[3]];
                    x = this.vertices[neighbours[0]];
                    y = this.vertices[neighbours[1]];
                    z = this.vertices[neighbours[2]];
                }

                w.angle = -MathHelper.toRad(36);
                x.angle = MathHelper.toRad(36);
                y.angle = -MathHelper.toRad(108);
                z.angle = MathHelper.toRad(108);

                w.globalAngle = angle + w.angle;
                x.globalAngle = angle + x.angle;
                y.globalAngle = angle + y.angle;
                z.globalAngle = angle + z.angle;

                this.createNextBond(w, vertex, w.globalAngle);
                this.createNextBond(x, vertex, x.globalAngle);
                this.createNextBond(y, vertex, y.globalAngle);
                this.createNextBond(z, vertex, z.globalAngle);
            }
        }
    }

    /**
     * Gets the vetex sharing the edge that is the common bond of two rings.
     *
     * @param {Vertex} vertex A vertex.
     * @returns {number|null} The id of a vertex sharing the edge that is the common bond of two rings with the vertex provided or null, if none.
     */
    getCommonRingbondNeighbour(vertex) {
        let neighbours = vertex.getNeighbours();
        
        for (let i = 0; i < neighbours.length; i++) {
            let neighbour = this.vertices[neighbours[i]];
            
            if (ArrayHelper.containsAll(neighbour.value.rings, vertex.value.rings)) {
                return neighbour;
            }
        }

        return null;
    }

    /**
     * Check if a vector is inside any ring.
     *
     * @param {Vector2} vec A vector.
     * @returns {boolean} A boolean indicating whether or not the point (vector) is inside any of the rings associated with the current molecule.
     */
    isPointInRing(vec) {
        for (let i = 0; i < this.rings.length; i++) {
            let ring = this.rings[i];

            if (!ring.positioned) {
                continue;
            }

            let polygon = ring.getPolygon(this.vertices);
            let radius = MathHelper.polyCircumradius(this.opts.bondLength, ring.getSize());
            let radiusSq = radius * radius;

            if (vec.distanceSq(ring.center) < radiusSq) {
                return true;
            }
        }

        return false;
    }

    /**
     * Check whether or not an edge is part of a ring.
     *
     * @param {Edge} edge An edge.
     * @returns {boolean} A boolean indicating whether or not the edge is part of a ring.
     */
    isEdgeInRing(edge) {
        let source = this.vertices[edge.sourceId];
        let target = this.vertices[edge.targetId];

        return this.areVerticesInSameRing(source, target);
    }

    /**
     * Check whether or not an edge is rotatable.
     *
     * @param {Edge} edge An edge.
     * @returns {boolean} A boolean indicating whether or not the edge is rotatable.
     */
    isEdgeRotatable(edge) {
        let vertexA = this.vertices[edge.sourceId];
        let vertexB = this.vertices[edge.targetId];

        // Only single bonds are rotatable
        if (edge.bondType !== '-') {
            return false;
        }

        // C-N bonds are not rotatable because of their high rotational energy barrier
        /* 
        let combination = vertexA.value.element.toLowerCase() + vertexB.value.element.toLowerCase();
        
        if (combination === 'cn' || combination === 'nc') {
            return false;
        }
        */

        // Do not rotate edges that have a further single bond to each side
        // If the bond is terminal, it doesn't make sense to rotate it
        if (vertexA.getNeighbourCount() + vertexB.getNeighbourCount() < 5) {
            return false;
        }

        // Ringbonds are not rotatable
        if (vertexA.value.rings.length > 0 && vertexB.value.rings.length > 0 && 
            this.areVerticesInSameRing(vertexA, vertexB)) {
            return false;
        }

        return true;
    }

    /**
     * Check whether or not a ring is an explicity defined aromatic ring (lower case smiles).
     *
     * @param {Ring} ring A ring.
     * @returns {boolean} A boolean indicating whether or not a ring is explicitly defined as aromatic.
     */
    isRingAromatic(ring) {
        for (let i = 0; i < ring.members.length; i++) {
            if (!this.isVertexInAromaticRing(ring.members[i])) {
                return false;
            }
        }

        return true;
    }

    /**
     * Checks whether or not an edge is part of an explicit aromatic ring (lower case smiles).
     *
     * @param {Edge} edge An edge.
     * @returns {boolean} A boolean indicating whether or not the vertex is part of an explicit aromatic ring.
     */
    isEdgeInAromaticRing(edge) {
        return this.isVertexInAromaticRing(edge.sourceId) &&
            this.isVertexInAromaticRing(edge.targetId);
    }

    /**
     * Checks whether or not a vertex is part of an explicit aromatic ring (lower case smiles).
     *
     * @param {number} vertexId A vertex id.
     * @returns {boolean} A boolean indicating whether or not the vertex is part of an explicit aromatic ring.
     */
    isVertexInAromaticRing(vertexId) {
        var element = this.vertices[vertexId].value.element;
        
        return element == element.toLowerCase();
    }

    /**
     * Get the normals of an edge.
     *
     * @param {Edge} edge An edge.
     * @returns {array} An array containing two vectors, representing the normals.
     */
    getEdgeNormals(edge) {
        let v1 = this.vertices[edge.sourceId].position;
        let v2 = this.vertices[edge.targetId].position;

        // Get the normals for the edge
        let normals = Vector2.normals(v1, v2);

        // Normalize the normals
        ArrayHelper.each(normals, function (v) {
            v.normalize()
        });

        return normals;
    }

    /**
     * Get the depth of a subtree in the direction opposite to the vertex specified as the parent vertex.
     *
     * @param {number} vertexId A vertex id.
     * @param {number} parentVertexId The id of a neighbouring vertex.
     * @returns {number} The depth of the sub-tree.
     */
    getTreeDepth(vertexId, parentVertexId) {
        let neighbours = this.vertices[vertexId].getSpanningTreeNeighbours(parentVertexId);
        let max = 0;

        for (let i = 0; i < neighbours.length; i++) {
            let childId = neighbours[i];
            let d = this.getTreeDepth(childId, vertexId);
            
            if (d > max) {
                max = d;
            }
        }

        return max + 1;
    }

    /**
     * Traverse a sub-tree in the graph.
     *
     * @param {number} vertexId A vertex id.
     * @param {number} parentVertexId A neighbouring vertex.
     * @param {function} callback The callback function that is called with each visited as an argument.
     * @param {number} [maxDepth=null] The maximum depth of the recursion. If null, there is no limit.
     * @param {boolean} [ignoreFirst=false] Whether or not to ignore the starting vertex supplied as vertexId in the callback.
     */
    traverseTree(vertexId, parentVertexId, callback, maxDepth = null, ignoreFirst = false, depth = 1, visited = []) {
        if (maxDepth !== null && depth > maxDepth + 1) {
            return;
        }

        for (let j = 0; j < visited.length; j++) {
            if (visited[j] === vertexId) {
                return;
            }
        }

        visited.push(vertexId);

        let vertex = this.vertices[vertexId];
        let neighbours = vertex.getNeighbours(parentVertexId);

        if (!ignoreFirst || depth > 1) {
            callback(vertex);
        }

        for (let i = 0; i < neighbours.length; i++) {
            this.traverseTree(neighbours[i], vertexId, callback, maxDepth, ignoreFirst, depth + 1, visited);
        }
    }

    /**
     * Gets the number of bonds of a vertex.
     *
     * @param {Vertex} vertex A vertex.
     * @returns {number} The number of bonds the vertex participates in.
     */
    getBondCount(vertex) {
        let count = 0;

        for (let i = 0; i < vertex.edges.length; i++) {
            count += this.edges[vertex.edges[i]].getBondCount();
        }

        return count;
    }

    /**
     * Returns an array of vertices that are neighbouring a vertix but are not members of a ring (including bridges).
     *
     * @param {number} vertexId A vertex id.
     * @returns {array} An array of vertices.
     */
    getNonRingNeighbours(vertexId) {
        let nrneighbours = [];
        let vertex = this.vertices[vertexId];
        let neighbours = vertex.getNeighbours();

        for (let i = 0; i < neighbours.length; i++) {
            let neighbour = this.vertices[neighbours[i]];
            let nIntersections = ArrayHelper.intersection(vertex.value.rings, neighbour.value.rings).length;
            
            if (nIntersections === 0 && neighbour.value.isBridge == false) {
                nrneighbours.push(neighbour);
            }
        }

        return nrneighbours;
    }

    /*
    getChiralOrder(vertexIds, chiralCenterVertexId) {
        let sortedVertexIds = Atom.sortByAtomicNumber(vertexIds, this.vertices);
        
        // Initial check whether there are duplicates, if not, all good to go
        if (!Atom.hasDuplicateAtomicNumbers(sortedVertexIds)) {
            return sortedVertexIds;
        }

        let done = new Array(vertexIds.length);
        let duplicates = Atom.getDuplicateAtomicNumbers(sortedVertexIds);
        
        let maxDepth = 1;
        for (let i = 0; i < duplicates.length; i++) {
            let dupl = duplicates[i];
            
            for (let j = 0; j < dupl.length; j++) {
                let index = dupl[j];
                let vertexId = sortedVertexIds[index].vertexId;
                let total = 0;
                
                console.log(vertexId, chiralCenterVertexId);
                this.traverseTree(vertexId, chiralCenterVertexId, function(vertex) {
                    console.log(vertex);
                    total += vertex.value.getAtomicNumber();
                }, maxDepth, true);

                sortedVertexIds[index].atomicNumber += '.' + total;
            }
        }

        sortedVertexIds = ArrayHelper.sortByAtomicNumberDesc(sortedVertexIds);
        console.log(sortedVertexIds);
        return sortedVertexIds;
    }
    */

    annotateChirality() {
        for (let i = 0; i < this.vertices.length; i++) {
            let vertex = this.vertices[i];
            
            if (vertex.value.bracket && 
                vertex.value.element.toLowerCase() === 'c' &&
                vertex.getNeighbourCount() === 4 ||
                vertex.value.bracket &&
                vertex.value.bracket.hcount > 0 && vertex.getNeighbourCount() === 3) {
                
                let chirality = vertex.value.bracket.chirality;
                
                if (chirality === null) {
                    continue;
                }

                let neighbours = vertex.getNeighbours();
                let orderedNeighbours = new Array(neighbours.length);

                this.vertices[vertex.parentVertexId].value.setOrder(vertex.id, 0);
                
                for (let j = 0; j < neighbours.length; j++) {
                    let neighbourId = neighbours[j];
                    
                    if (neighbourId === vertex.parentVertexId) {
                        orderedNeighbours[0] = neighbourId;
                        continue;
                    }

                    orderedNeighbours[this.vertices[neighbourId].value.getOrder(vertex.id)] = neighbourId;
                }
                
                if (chirality === '@') {
                    let edgeUp = this.getEdge(orderedNeighbours[3], vertex.id);
                   
                    // If the bond already points down here, there is no need to point up
                    // into the other direction
                    if (!(edgeUp.chiral === 'down')) {
                        edgeUp.chiral = 'up';
                    }
                    
                    
                    let edgeDown = this.getEdge(orderedNeighbours[1], vertex.id);
                    
                    // If the bond already points up here, there is no need to point down
                    // into the other direction
                    if (!(edgeDown.chiral === 'up')) {
                        edgeDown.chiral = 'down';
                    }
                } else if (chirality === '@@') {
                    let edgeUp = this.getEdge(orderedNeighbours[1], vertex.id);
                   
                    // If the bond already points down here, there is no need to point up
                    // into the other direction
                    if (!(edgeUp.chiral === 'down')) {
                        edgeUp.chiral = 'up';
                    }
                    
                    
                    let edgeDown = this.getEdge(orderedNeighbours[3], vertex.id);
                    
                    // If the bond already points up here, there is no need to point down
                    // into the other directiononsole.log(vertex, ctn);onsole.log(vertex, ctn);
                    if (!(edgeDown.chiral === 'up')) {
                        edgeDown.chiral = 'down';
                    }
                }
            }
        }
    }

    /**
     * Creates pseudo-elements (such as Et, Me, Ac, Bz, ...) at the position of the carbon sets
     * the involved atoms not to be displayed.
     */
    initPseudoElements() {
        for (let i = 0; i < this.vertices.length; i++) {
            const vertex = this.vertices[i];
            const neighbours = vertex.getNeighbours();

            // Ignore atoms that have less than 3 neighbours, except if
            // the vertex is connected to a ring and has two neighbours

            if (vertex.getNeighbourCount() < 3 &&
                !(vertex.value.isConnectedToRing && vertex.getNeighbourCount() === 2)) {
                continue;
            }

            let ctn = 0;

            for(let j = 0; j < neighbours.length; j++) {
                let neighbour = this.vertices[neighbours[j]];

                if (neighbour.getNeighbourCount() > 1) {
                    ctn++;
                }
            }

            if (ctn > 1) {
                continue;
            }

            // Get the previous atom (the one which is not terminal)
            let previous = null;

            for(let j = 0; j < neighbours.length; j++) {
                let neighbour = this.vertices[neighbours[j]];
                if (neighbour.getNeighbourCount() > 1) {
                    previous = neighbour;
                }
            }


            for(let j = 0; j < neighbours.length; j++) {
                let neighbour = this.vertices[neighbours[j]];
                
                if (neighbour.getNeighbourCount() > 1) {
                    continue;
                }

                neighbour.value.isDrawn = false;
                
                let hydrogens = this.maxBonds[neighbour.value.element] - this.getBondCount(neighbour);
                
                if (neighbour.value.bracket) {
                    hydrogens = neighbour.value.bracket.hcount;
                }

                vertex.value.attachPseudoElement(neighbour.value.element, previous ? previous.value.element : null, hydrogens);
            }
        }
    }

    /**
     * Cleans a SMILES string (removes non-valid characters)
     *
     * @static
     * @param {string} smiles A SMILES string.
     * @returns {string} The clean SMILES string.
     */
    static clean(smiles) {
        return smiles.replace(/[^A-Za-z0-9@\.\+\-\?!\(\)\[\]\{\}/\\=#\$:\*]/g,'');
    }

    /**
     * Applies the smiles drawer draw function to each canvas element that has a smiles string in the data-smiles attribute.
     *
     * @static
     * @param {objects} options SmilesDrawer options.
     * @param {string} [themeName='light'] The theme to apply.
     */
    static apply(options, themeName='light') {
        let smilesDrawer = new SmilesDrawer(options);
        let elements = document.querySelectorAll('canvas[data-smiles]');

        for (let i = 0; i < elements.length; i++) {
            let element = elements[i];
            let data = SmilesDrawer.parse(SmilesDrawer.clean(element.getAttribute('data-smiles')));

            smilesDrawer.draw(data, element, themeName, false);
        }

    }

    /**
     * Parses the entered smiles string.
     * 
     * @static
     * @param {string} smiles A SMILES string.
     * @param {Function} successCallback A callback that is called on success with the parse tree.
     * @param {Function} errorCallback A callback that is called with the error object on error.
     */
    static parse(smiles, successCallback, errorCallback) {
        try {
            if (successCallback) {
                successCallback(SMILESPARSER.parse(smiles));
            }
        } catch (err) {
            if (errorCallback) {
                errorCallback(err);
            }
        }
    }
}