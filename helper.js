/** 
 * Get label info from newMuts
 * 
 * Used in: tree.js for tree nodes 
 * 					main.js for threshold slider labels
 * Parameters: 
 * mutLabels -- labels for mutations provided in JSON
 * newMuts -- 
 * phyloTree -- new tree
 *
 */


//TODO: will always be subset now since we changed infrastructure 

function getLabels(mutLabels, phyloTree, newMuts, subset) {
	var treeLabels = {} 
	var cloneSymbols = {} 
	if(subset)	{ 	//subset
		if(mutLabels.length != 0) {	//labeled
			for(var node in phyloTree) {
				if(newMuts[node].length == 1) {
					var labelArray = mutLabels[newMuts[node]-1]
					var label = ""
					for(var mut in labelArray) {
						label += labelArray[mut] + ", "
					}
					label = label.substr(0, label.length-2)
					treeLabels[node] = label
					cloneSymbols[node] = mutSymbols[newMuts[node]-1]
				} else {
					// TODO: to make PNAS work, had to alter newMuts structure...
					// double for loop to get each individual mut

					var label = ""
					for(var i = 0; i < newMuts[node].length; i++) {
						var cluster = mutLabels[parseInt(newMuts[node][i])-1]
						
						if(selectedMutations.includes(labelToIdx[cluster])) {
							for(var mut = 0; mut < cluster.length; mut++)
								label += cluster[mut] + ", "
						}
					}
					treeLabels[node] = label 
					cloneSymbols[node] = mutSymbols[parseInt(newMuts[node][0])-1] 
				}
			}
		} else {	//not labeled
			for(var node in phyloTree) {
				if(newMuts[node].length == 1) {
					treeLabels[node] = "mut_" + node
					cloneSymbols[node] = mutSymbols[newMuts[node]-1]
				} else {
					var label = ""
					var cloneLabel = ""
					for(var i = 0; i < newMuts[node].length; i++) {
						label += "mut_" + parseInt(newMuts[node][i]) + ",\n"
						cloneLabel += mutSymbols[parseInt(newMuts[node][i])-1] 
					}
					treeLabels[node] = label 
					cloneSymbols[node] = cloneLabel
				}
			}
		}
	} else {	//not subset
		// start at i = 1 to match with order
		for(var i = 1; i < numMutations+1; i++) {
			var label = "";
			var cloneLabel = ""
			if(mutLabels.length != 0) {
				label = mutLabels[i-1]
			} else {
				label = "Mut " + i
			}
			cloneSymbols[i] = mutSymbols[i-1]
			treeLabels[i] = label
		}
	}
	return [treeLabels, cloneSymbols];
}




function getCloneSymbols(old2new, mutSymbols) {
	var cloneSymbols = {}
  for(var node in old2new) {
  	var symbols = ""
  	for(var child in old2new[node]) {
  		symbols += mutSymbols[old2new[node][child]-1]
  	}
    cloneSymbols[old2new[node]] = symbols
  }
  return cloneSymbols
}


/** 
 * Get legend info 
 * Recursively (modified getIntersects function)
 *
 * NOTE: no longer in use (replaced with tree)
**/

function getLegend(tree) {
	var union;
	var root = findRoot(tree)
	var clones = recurse(root, root, {})

	function recurse(node, parent, clones) {
		//console.log("Node: ", node)
		//console.log("Parent: ", parent)
		if (node != root) {
			(function () {
				union = clones[parent].slice()
				union.push(node)
				//get intersection of mutation, and union of current node
				//plus parent's old2new value 
				//clones[node] = union.filter(value => -1 !== selectedMutations.indexOf(value));
				clones[node] = union;
			}());
		} else {	//root case
			var union = []
			union.push(node)
			//old2new[node] = union.filter(value => -1 !== selectedMutations.indexOf(value));
			clones[node] = union;
		}
		for(var index in tree[node]) {
			//console.log("child: ", tree[node][index], " in ", tree[node])
			recurse(tree[node][index], node, clones)

		}
		return clones
	}
	var cloneList = {}
	if(mutLabels.length != numMutations) {
		for(var node in clones) {
      if (clones[node] in cloneList) {
          cloneList[clones[node]].push(node)
      } else {
          cloneList[clones[node]] = [node]
      }
  	}
	} else {
		for(var node in clones) {
			//-1 for zero index
      if (clones[node] in cloneList) {
      		try {
      			cloneList[clones[node]].push(mutLabels[node-1])
      		} catch (e) {
            console.log(e);
            console.log("Error occured while mapping mutation to labels")
        }
      } else {
      	try {
      		cloneList[clones[node]] = [mutLabels[node-1]]
      	} catch (e) {
            console.log(e);
            console.log("Error occured while mapping mutation to labels")
        } 
      }
  	}
	}
	return cloneList
}

/** 
 * Find the root of given phylogeny tree	
 * Returns the VALUE of the root
**/

function findRoot(tree) {
	var nodes = []
	for(var node in tree) {
		nodes.push(parseInt(node))
	}
	for(var node in tree) {
		for(var i = 0; i < tree[node].length; i++) {
			if (nodes.includes(tree[node][i])) { //if child in nodes
				nodes.splice(nodes.indexOf(tree[node][i]), 1)
			}
		}
	}
	
	if(nodes.length > 1) {
		return "Error: More than one node in parent list."
	} else {
		return nodes[0]
	}
}

/**
 * Find root of new phylogeny tree,
 * Where nodes are not necessarily a single digit
 */
function findNewRoot(tree) {
	var nodes = []
	for(var node in tree) {
		nodes.push(node)
	}
	for(var node in tree) {
		for(var i = 0; i < tree[node].length; i++) {
			var curr = JSON.stringify(tree[node][i])
			curr = curr.substring(1, curr.length-1)
			if (nodes.includes(curr)) { //if child in nodes
				nodes.splice(nodes.indexOf(curr), 1)
			}
		}
	}
	if(nodes.length > 1) {
		return "Error: More than one node in parent list."
	} else {
		return nodes[0]
	}
}


/*
 * When users select a subset of mutations! 
 * 
 *
 *
 */

function getIntersects(tree, selectedMutations) {
	var union;
	var root = findRoot(tree)
	var old2new = recurse(root, root, {})

	function recurse(node, parent, old2new) {
		//console.log("Node: ", node)
		//console.log("Parent: ", parent)
		if (node != root) {
			(function () {
				union = old2new[parent].slice()
				union.push(node)
				//get intersection of mutation, and union of current node
				//plus parent's old2new value 
				old2new[node] = union.filter(value => -1 !== selectedMutations.indexOf(value));
			}());
		} else {	//root case
			var union = []
			union.push(node)
			old2new[node] = union.filter(value => -1 !== selectedMutations.indexOf(value));
		}
		for(var index in tree[node]) {
			//console.log("child: ", tree[node][index], " in ", tree[node])
			recurse(tree[node][index], node, old2new)

		}
		return old2new
	}
	return old2new
}

/***
 *
 *  If there exists an edge between two nodes in the old phylogeny tree,
 *  but their old2new values are different, then there should be an edge in
 *  the new tree drawn from the parent's old2new value to the child's old2new value. 
 *
 *  params: oldTree is the original phylogeny tree. 
 *					old2new is the mapping of original mutations to the new mutation groupings
 ***/

function getNewTree(oldTree, old2new, newMuts) {
	newTree = {} 
	// initialize new treenewMuts with new mutation groupings
	for(var node in newMuts) {
		newTree[node] = []
	}
	for(var node in oldTree) {
		for(var child = 0; child < oldTree[node].length; child++) {
			//using JSON.stringify for array comparison 
			if(JSON.stringify(old2new[oldTree[node][child]]) != JSON.stringify(old2new[node])) {
				newTree[old2new[node]].push(old2new[oldTree[node][child]])
			}
		}
	}

	return newTree
}

/* Get new U from old U matrix */
function getNewU(newMuts, matrixU, idx, order) {
	var newU = [] 
	var newIdx = {}
	for(var row = 0; row < matrixU.length; row++) {
	    newU[row] = new Array(Object.keys(newMuts).length).fill(0)
	}
	var counter = 0
  for(var mut in order) {
  	
  	//mut is an index
  	//order[mut] represents the node
  	for(var child in newMuts[order[mut]]) {
  		//for (var i in newMuts[mut][child]) {
  			for(var row = 0; row < matrixU.length; row++) {
	  			newU[row][counter] += matrixU[row][idx[newMuts[order[mut]][child]]]
	  			newU[row][counter] = +(newU[row][counter].toFixed(2))
	  		}
  		//}
  	}
  	newIdx[order[mut]] = counter
  	counter += 1
  }

  return [newU, newIdx]
}


/* 0 new matrix F from new matrix U, recursively*/
function getNewF(tree, matrixU, newIdx) {
	var newMatrixF = []
	for(var row = 0; row < matrixU.length; row++) {
	  newMatrixF[row] = matrixU[row].slice()
	}

	var root = findNewRoot(tree)
	var newMatrixF = recurse(root, root, matrixU, tree, newIdx, newMatrixF) 
	function recurse(node, parent, matrixU, tree, newIdx, newF) {
		for(var child in tree[node]) {
			recurse(tree[node][child], node, matrixU, tree, newIdx, newF)
		}
		if (node != root) {	//not root 
			for(var row = 0; row < matrixU.length; row++) {
				newF[row][newIdx[parent]] += newF[row][newIdx[node]]
				newF[row][newIdx[parent]] = +(newF[row][newIdx[parent]].toFixed(2))
			}
		}
		return newF
	}
	return newMatrixF
}

function prettyPrint(matrix) {
	var pretty = '<table width = "200px">'

	for(var row = 0; row < matrix.length; row++) {
		pretty += "<tr>"
		for(var col = 0; col < matrix[0].length; col++) {
			if (col == 0) {
				pretty += ('<td>' + row + ': </td>')
			}
			pretty += ('<td>' + matrix[row][col] + '</td>')
		}
		pretty += '</tr>'
	}
	pretty += '</table>'
	return pretty 

}














/**
 * Overlapping 
 *
 *
 *
 *
 *
 *
 */

// 1. Determine the "order"//hierarchy of all nodes (from tree, then from index)
/**
  assumption:
  last node in ORDER (bfsorder; all the mutations in the clone)
  is the CURRENT node 
**/
// returns ordering for overlaps 
// parent should be largest, and then children, and for siblings, whichever comes first 

function bfsOrder(tree, root) {
	var order = []
	var queue = []

	// var root = findRoot(tree)
	
	queue.push(root)

	while(queue.length != 0) {
		var curr = queue.shift()
		if (tree[curr].length != 0) { 	// add children to queue
			for(child in tree[curr]) {
				queue.push(tree[curr][child])
			}
		} 
		order.push(curr)
	}
	return order 
}


// 2. Swap orders of matrix U according to the bfsOrder 

function removeOverlap(matrix, order, idx) {
	
	//epsilon represents erode radius (the greater, the smaller)
	var eps = 4; 

	var erodeVal = {}
	for(var mut = 0; mut < order.length; mut++) {
		erodeVal[order[mut]] = 0
	}	

	for(var row = 0; row < matrix.length; row++) {
		// use order[col] for new column index to go in order
		// this replaces the step of reordering the matrix 
		var dict = {}

		for(var col = 0; col < matrix[0].length; col++) {		
			var curFreq = matrix[row][idx[order[col]]]
			if (curFreq != 0 && curFreq in dict) {	// there is dup value; increment
				// set epsilon value (accounting for nesting) 
				if(erodeVal[order[col]] < dict[curFreq]) {
					erodeVal[order[col]] = (dict[curFreq] * eps)	
				}
				dict[curFreq] += 1
			} else {	// no dup value; add to dict 
				dict[curFreq] = 1
			}
		}
	}

	var temp = {} 
	var filters = []
	for(var key in erodeVal) {
		if(!(erodeVal[key] in temp)) {
			temp[erodeVal[key]] = 1
			filters.push(erodeVal[key])
		} 
	}
	return [erodeVal, filters];
}




















function addStyle(str) {
    var node = document.createElement('style');
    node.innerHTML = str;
    document.body.appendChild(node);
}
