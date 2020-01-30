var debugging = false; 

var selectedTreeIdx = 0;

/* Initialize JSON inputs */
var fullMatrixF = data.input;
var phyloTree = phylogeny[selectedTreeIdx]; 

var test = test;

var uThreshold = data.u_thresholds;
var fThreshold = data.f_thresholds;

/* Size of grid */
var n = data.cols;  // columns
//var m = fullMatrixF.length/n;
var m = data.rows;

var numMutations = fullMatrixF[0].length;
var mutLabels = data.mut_labels;

var globalContourIdx = -1;
var globalSliderIdx = -1;

var erode = false;

var labelToIdx = {}
var idxToLabel = {}
for (var cluster = 0; cluster < mutLabels.length; cluster++) {
    labelToIdx[mutLabels[cluster]] = cluster+1
}
for(var label in labelToIdx) {
    var newKey = labelToIdx[label]
    var newVal = ""
    var temp = label.split(",")
    for(var mut in temp) {
        newVal += temp[mut] + ", "
    }
    newVal = newVal.substr(0,newVal.length-2)
    idxToLabel[newKey] = newVal
}

//Color Brewer: Qualitative 1 + 2
var colorRanges = 
    ['#a6cee3', '#e31a1c', '#1f78b4', '#b15928', '#b2df8a', '#33a02c', '#fb9a99', '#ffff99']
//['#1a6699','#7cbde9','#d3e9f8','#60a428', '#93d75b', '#e4f5d6', '#930806','#f73e3b','#fdcfce', '#955b04','#faa21e','#feeacd','#ddcde5','#9768b1','#623d76','#cccc00','#ffff00','#ffffb3', '#2a6f62','#8dd3c7','#c8eae4']
/**
    colorIdx assigns each mutation (mutation introduction) to color. this color is associated
    with its corresponding mutation throughout the application
    NOTE: this updates with every mutation selection
    
    indexToColorIdx converts the old, bfsOrder-based Index (used for selecting mutations and slider contour)
    to the color that should be used, so that all colors are uniform between trees. 

    ONE-INDEXED: 
    MUT-->color
 */ 

var colorIdx = {}
var indexToColorIdx={}

// Object.keys(phyloTree).sort((a,b)=>a-b) sorts by NUMERIC value instead of string
arr=Object.keys(phyloTree).sort((a,b)=>a-b)
for(var i in arr) {
    i=parseInt(i)+1
    colorIdx[i]=colorRanges[i-1]
}
//colorIdx[""]=colorRanges[0]

var availableColors = colorRanges;
// HEX unicode symbols to represent mutations in clones

//var mutSymbols = ["&#x2605;", "&#x2660;", "&#x2662;", "&#x2666;", "&#x2688;"]
var mutSymbols = ["&#x2605;", "&#x2660;", "&#x2662;", "&#x2666;", "&#x2688;", 
                  "&#x213D;", "&#x2601;", "&#x2606;", "&#x262F;", "&#x263D;",
                  "&#x2659;", "&#x265C;", "&#x265E;", "&#x2665;", "&#x266B;",
                  "&#x2690;", "&#x2698;", "&#x2744;", "&#x2716;", "&#x2724;", ""]
/*
 * Global variable to keep track of which 
 * nodes are currently hidden from the user
 *
 * Initialized: in drawCanvas()
 * Used: in tree.js
 */
var activatedMutations = {}

/**
 * Get matrix U 
 * Loop through phylogeny tree and substract 
 * children's frequencies from the parent 
 *
 * Used for: both subset and original mutation drawing
 *
 */
function getMatU(f, tree, idx) {
    var matU = [];
    for(var row = 0; row < f.length; row++) {
        matU[row] = f[row].slice();
    }
    for (var mut in tree) {
        for (var child in tree[mut]) {
            for(var row = 0; row < f.length; row++) {
                matU[row][idx[mut]] -= f[row][idx[tree[mut][child]]]
                // account for javascript floating point arithmetic
                matU[row][idx[mut]] = +(matU[row][idx[mut]].toFixed(2))
            }
        }
    }
    return matU
}

/* Given a phylo tree, idx returns the index of a node in matrix F */
var idx = {}
var counter = 0
for(mut in phyloTree) {
    idx[mut] = counter
    counter +=1 
}

// var fullMatrixU = getMatU(fullMatrixF, phyloTree, idx)

// get all mutations for old2new
// for clone symbols
var selectedMutations = []
if(data.selected_mutations.length == 0) {
    for(var node in phyloTree) {
        selectedMutations.push(parseInt(node))
    }
} else {
    selectedMutations = data.selected_mutations
}
var htmlSelectedMutations = []
for (var cluster = 0; cluster < mutLabels.length; cluster++) {
    var currIdx = labelToIdx[mutLabels[cluster]]
    for (var mut = 0; mut < mutLabels[cluster].length; mut++) {
        if (selectedMutations.includes(currIdx)) {
            htmlSelectedMutations.push(mutLabels[cluster][mut])
        } 
    }
}


getSubsetMutations(selectedMutations)
// var old2new = getIntersects(phyloTree, allMutations)
// var clones = getCloneSymbols(old2new, mutSymbols)


// drawTree(numMutations, phyloTree, colorRanges, order, treeLabels, clones)

function getSubsetMutations(selectedMutations) {
    /* 
     * reset some variables
     * phyloTree --> subsets should use the original tree to create new one
     * matrixF --> reset to match phyloTree
     * activatedMutations --> map resets when new mutations are selected
     */
    var phyloTree = phylogeny[selectedTreeIdx];
    
    //pass by reference even if reading from JSON!
    var temp=data.input;
    var matrixF = [];
    for(var row = 0; row < temp.length; row++) {
        matrixF[row] = temp[row].slice();
    }
    

    numMutations = matrixF[0].length;

    activatedMutations = {};

    /* Given a phylo tree, idx returns the column 
     * index of a node in matrix F */
    idx = {}
    var counter = 0
    for(mut in phyloTree) {
        idx[mut] = counter
        counter +=1 
    }

    var matrixU = getMatU(matrixF, phyloTree, idx)

    /* 
     * start subset process  !!! 
     *
     *
     */

    // step 2. recurse(root)
    var old2new = getIntersects(phyloTree, selectedMutations)
    // 3. get dictionary of new mutations
    // labeledMuts is used for F legend 
    var newMuts = {} 
    var labeledMuts = {}
    for(var node in old2new) {
        if (old2new[node] in newMuts) {
            newMuts[old2new[node]].push(node)
            
            if(mutLabels.length != 0) {
                var l_str = [];
                try {
                    for(var l_node in old2new[node]) {
                        l_str.push(mutLabels[old2new[node][l_node]-1])
                    }
                    labeledMuts[l_str].push(mutLabels[node-1])
                } catch(e) {
                    console.log(e);
                    console.log("Error occured while creating labels for legend")
                }
            }
        } else {
            newMuts[old2new[node]] = [node]
            if(mutLabels.length != 0) {
                var l_str = []
                try {
                    for(var l_node in old2new[node]) {
                        l_str.push(mutLabels[old2new[node][l_node]-1])
                    }
                    labeledMuts[l_str] = [mutLabels[node-1]]
                } catch(e) {
                    console.log(e);
                    console.log("Error occured while creating labels for legend")
                }

            }
        }
    }
    //4. generate new tree and newU 
    //newU = add up values of newMuts to the old keys 
    phyloTree = getNewTree(phyloTree, old2new, newMuts)

    /* root and order needed for drawing colors 
     * and eroding 
     * getting order now allows idx to assign indicies
     * according to their order 
     */

    
    var root = findNewRoot(phyloTree)
    order = bfsOrder(phyloTree, root)



    /* index to color idx
     * basically, converts the BFS order to the actual mutation # 
     */

    for(var o_index in order) {
        curr_node=parseInt(order[o_index][order[o_index].length-1])
        indexToColorIdx[o_index]=curr_node
    } 

    var newInfo = getNewU(newMuts, matrixU, idx, order)
    var newU = newInfo[0]
    idx = newInfo[1]
    
    
    matrixF = getNewF(phyloTree, newU, idx)
    matrixU = getMatU(matrixF, phyloTree, idx)

    numMutations = matrixU[0].length
    
    /***
     * Draw mutation legend
     * for Matrix F 
     ***/


    var legendlabels = (Object.keys(labeledMuts).length == 0) ?  newMuts : labeledMuts
    let flegend = document.getElementById("flegend");
    //flegend.innerHTML = "<h3> Legend </h3>"

    if(Object.keys(labeledMuts).length == 0) {
        console.log("Need to provide mutation labels for frequency legend")
        // var index = 0
        // for (var node in order) {
        //     flegend.innerHTML = flegend.innerHTML + "<p><span style = 'background-color:" + colorRanges[index] + "'></span>" + "Mutation " + order[node] + " </p>"
        //     index += 1
        // }
    } else {
        for (var i = 0; i < selectedMutations.length; i++) {
            var mut = idxToLabel[selectedMutations[i]]
            if(mut.length > 30) {
                mut = mut.substr(0,30)
            }
            //flegend.innerHTML = flegend.innerHTML + "<p><span style = 'background-color:" + colorRanges[i] + "'></span>" + mut + " </p>"
        }
    }



    currMatrixU = matrixU
    currMatrixF = matrixF
    

    
    var erodeInfo = removeOverlap(currMatrixU, order, idx);
    erodeDict = erodeInfo[0]
    filters = erodeInfo[1]


    var labels = getLabels(mutLabels, phyloTree, newMuts, true)
    var treeLabels = labels[0]
    var cloneLabels = labels[1]
    
    //reset thresholds to 0.1 
    for(var i = 0; i < numMutations; i++) {
        uThreshold[i] = 0.1
    }
    //drawThresholds(treeLabels); 



    visualize(currMatrixU, currMatrixF, erode)

    

    var clones = getCloneSymbols(old2new, mutSymbols)

    // sort phyloTree to make sure same trees will be drawn the same 
    // (by drawing the edges in the same order)
    for (var node in phyloTree){
        phyloTree[node].sort()
    }
    //TODO: need to redraw when matrix F is toggled
    drawTree(numMutations, phyloTree, colorRanges, order, treeLabels, clones, cloneLabels)
    addTreeInteractivity() 
    //clear canvases 
    //document.getElementById("contours").innerHTML = '<svg></svg>';
    thresholdRefresh();
}

function visualize(matrixU, matrixF, erode) {
    if(debugging) {
        var tree = JSON.stringify(phyloTree)
        document.getElementById("debug-tree").innerHTML = "Phylogeny Tree: <br><br>" + tree + "<br><br>"

        var uprint = prettyPrint(matrixU)
        document.getElementById("debug-u").innerHTML = "Matrix U: <br><br>" + uprint;

        var fprint = prettyPrint(matrixF)
        document.getElementById("debug-f").innerHTML = "Matrix F: <br><br>" + fprint;

        //clear canvas before drawing
        document.getElementById("clones").innerHTML = cloneBoxHTMl;
        //document.getElementById("frequency").innerHTML = freqBoxHTML;
    }

    drawCanvas();

    function drawCanvas() {
        var temp = [];
        /* empty svg initialization */
        //weird values to fit with border
        var xRange = d3.range(0, n+3);   //rows
        var yRange = d3.range(-1, m+4);   //cols

        var width = 18*n;
        var height = width * (yRange.length / xRange.length);

        
        //get bigger outline for matrix U 

        try {
            drawMatrix(matrixU, "#clones", uThreshold);
            //drawMatrix(matrixF, "#frequency", fThreshold);
        } catch (e) {
            console.log(e);
            console.log("Error occured while drawing matrix")
        }


        /** Function to draw Matrix U and Matrix F 
            matrix: input for matrix u or matrix f 
            div: div ID to draw the SVG visualization in 
            thresholds: thresholds (for each mutation)

        */ 
        function drawMatrix(matrix, div, thresholds) {

            var currmatrix = 'u';

            if(div=="#frequency") {
                currmatrix = 'f'
            }

            var viewbox = "0 0 " +(width) + " "+(height)
            /* Grab SVG from dom */
            var svg = d3.select(div+ " svg")
                .attr("id", "matrixUgrid")
                .attr("width", width)
                .attr("height", height)
                .attr("viewBox", viewbox)

            var group1 = svg.append("g")
                .attr("class", "group1")
            
            group1.append("clipPath")
                .attr("id", "rect-clip")
                .append("rect")
                .attr("x", width/(n+2))
                .attr("y", height/(m+1))
                .attr("width", width-(2*(width/(n-1)))-15)
                .attr("height", height-(2*(height/(m+1))))

            if(erode) {
                for(var index = 0; index < filters.length; index++) {
                    group1.append("filter")
                        .attr("id", "erode-"+currmatrix+filters[index])
                        .append("feMorphology")
                        .attr("operator", "erode")
                        .attr("in", "SourceGraphic")
                        .attr("radius", filters[index] + 4)

                    group1.append("filter")
                        .attr("id", "erode-clip-"+currmatrix+filters[index])
                        .append("feMorphology")
                        .attr("operator", "erode")
                        .attr("in", "SourceGraphic")
                        .attr("radius", filters[index])

                    var outline = group1.append("filter")
                        .attr("id", "outline-"+currmatrix+filters[index])

                    outline.append("feMorphology")
                        .attr("operator", "erode")
                        .attr("in", "SourceGraphic")
                        .attr("result", "OUTER")
                        .attr("radius", filters[index])

                    outline.append("feMorphology")
                        .attr("operator", "erode")
                        .attr("in", "SourceGraphic")
                        .attr("result", "INNER")
                        .attr("radius", filters[index]+1)

                    outline.append("feComposite")
                        .attr("operator", "out")
                        .attr("in", "OUTER")
                        .attr("in2", "INNER")
                }
            }
        
            /* Loop for each mutation (i.e., column in input matrix) */
            for(var index = 0; index < numMutations; index++) {
                
                /** Determine threshold */
                var threshold = [];     //MarchingSquares requires an array

                if (thresholds.length < numMutations) {
                    threshold.push(thresholds[0]);
                } else {
                    threshold.push(thresholds[index]);
                }

                var grid = [];
                temp = [];

                /* Turn column into grid (with n columns) */
                for(var i = 0; i < matrixF.length; i++) {
                    if(i!=0 && i%n == 0) {
                        grid.push(temp);
                        temp = [];
                    }
                    temp.push(matrix[i][index]);
                }
                grid.push(temp);

                // create a border of zeros around the grid!! 
                // for creating a closed path and to visualize entire tumor

                var borders = [] 
                

                for(var row = 0; row < grid.length+2; row++)
                    borders[row] = new Array(grid.length+2).fill(0)


                for(var row = 1; row < borders.length-1; row++) {
                    for(var col = 2; col < borders[0].length-1; col++) {
                        borders[row][col] = grid[row-1][col-1]
                    }
                }

                var curr_node=order[index][order[index].length-1]
    
                var isoLines = [];
                MarchingSquaresJS
                  .isoLines(borders,
                            threshold,
                            {
                                polygons: false,
                                linearRing: false
                            }
                  )
                  .forEach(function(isolines, i) {
                    if(i < colorRanges.length) {
                        isoLines.push({
                          "coords": isolines,
                          "level": i + 1,
                          //"color": colorRanges[index],
                          "color": colorIdx[indexToColorIdx[index]],
                          "val": threshold[i]
                      });
                    }
                    
                  });

                drawLines(div, isoLines, thresholds, xRange, yRange, false);

                // helper function
                /** 
                 *   divId: div to draw the svg in
                 *   lines: isoLines to be drawn (given by grid values)
                 *   intervals: thresholds 
                 *   xs: rows (n)
                 *   ys: columns (m) 
                 */
                function drawLines(divId, lines, intervals, xs, ys, drawingGlobal) {

                    var xScale = d3.scale.linear()
                        .range([0, width])
                        .domain([Math.min.apply(null, xs), Math.max.apply(null, xs)]);

                    var yScale = d3.scale.linear()
                        .range([0, height])
                        .domain([Math.min.apply(null, ys), Math.max.apply(null, ys)]);

                    var identifier = currmatrix + index
                    if(drawingGlobal) {
                        identifier = currmatrix + "n"
                    }
                    

                    var def = group1.append("defs")
                        .attr("opacity", 1.0)

                    var fillcolor = "black"

                    var cloudPathID = "cloud-" + identifier 
                    if(!(cloudPathID in activatedMutations)) {
                        activatedMutations[cloudPathID] = "activated"
                    }
                    def.selectAll("path")
                        .data(lines)
                        .enter()
                        .append("path")
                        .attr("id", cloudPathID)
                        .attr("class", activatedMutations[cloudPathID])
                        .style("stroke-linecap", "round")
                        // .style('opacity', 1.0)
                        .attr("d", function (d) {
                            fillcolor=d.color;
                            var p = "";
                            d.coords.forEach(function (aa) {
                                p += (d3.svg.line()
                                    .x(function (dat) {
                                        return xScale(dat[0]);
                                    })
                                    .y(function (dat) {
                                        return yScale(dat[1]);
                                    })
                                    .interpolate("linear")
                                )(aa) + "";
                            });
                            return p;
                        })

                    // def.append("clipPath").attr("id","clip"+currmatrix+index)
                    //     .append("use").attr("xlink:href", "#cloud"+currmatrix+index);

                    // svg.append("g")
                    //     .append("use").attr("xlink:href","#cloud"+currmatrix+index)
                    //     .attr("clip-path","url(#clip"+currmatrix+index+")")
                    
                    var mask = def.append("mask")
                        .attr("id", "mask-"+identifier)
                        
                    mask.append("rect")
                        .attr("x", 0)
                        .attr("y", 0)
                        .attr("width", width)
                        .attr("height", height)
                        .attr("fill", "white")

                    if(erode) {
                        mask.append("use")
                            .attr("xlink:href", "#cloud-"+identifier)
                            .attr("fill", "black")      // must be black to hide mask in middle
                            .attr("filter", "url(#erode-"+currmatrix+erodeDict[order[index]]+")")
                            .attr("class", "fill-mask")

                        group1.append("g")
                            .attr("class", "cloud")
                            .append("use")
                            // .attr("class","useuse")
                            .attr("id", "color-boundary-"+identifier)
                            .attr("xlink:href","#cloud-"+identifier)
                            .attr("mask", "url(#mask-" +identifier +")")
                            .attr("clip-path", "url(#rect-clip)")
                            .attr("filter","url(#erode-clip-"+currmatrix+erodeDict[order[index]]+")")
                            .attr("fill", fillcolor)
                            .attr("opacity", 0.8)

                        //draw outlines
                        group1.append("g")
                            .attr("class", "outline")
                            .append("use")
                            .attr("xlink:href","#cloud-"+identifier)
                            .attr("filter","url(#outline-"+currmatrix+erodeDict[order[index]]+")")
                            .attr("clip-path","url(#rect-clip)")
                            .attr("fill", "black")    
                    } else {
                        var strokeWidth = 1
                        var cloudId = ""
                        if(drawingGlobal) {
                            strokeWidth = 3
                            cloudId = "remove-cloud-click"
                        }
                        group1.append("g")
                            .append("use")
                            // .attr("class","cloud")
                            .attr("class","cloud #cloud-selector-"+identifier)
                            .attr("id", "color-boundary-"+identifier)
                            .attr("id",cloudId)
                            .attr("xlink:href","#cloud-"+identifier)
                            .attr("mask", "url(#mask-" +identifier +")")
                            .attr("clip-path", "url(#rect-clip)")
                            .attr("fill", fillcolor)
                            .attr("stroke", "#000")
                            .attr("stroke-width", strokeWidth)
                            .attr("opacity", 0.8)

                    }

                    if(currmatrix = 'u') {
                        group1.selectAll("#clones .cloud")
                        .on('mouseenter', function () {
                            var g = this.cloneNode(true)
                            g.class = "temp-cloud"
                            g.id = "remove-cloud"
                            
                            //append to svg to change z-index
                            this.parentElement.parentElement.appendChild(g);

                            d3.selectAll("#remove-cloud")
                                .attr("stroke-width", 3)
                            
                            
                            //remove on mouseout
                            d3.selectAll("#remove-cloud")
                              .on("mouseleave", function() { 
                                while (document.getElementById("remove-cloud")) {
                                    var element = document.getElementById("remove-cloud")
                                    element.parentNode.removeChild(element);
                                }
                               
                            })

                            .on('click', function () {
                                removeGlobal()
                                //first click
                                var g = this.cloneNode(true)
                                //can only have one clicked at once
                                // while (document.getElementById("remove-cloud-click")) {
                                //     var element = document.getElementById("remove-cloud-click")
                                //     element.parentNode.removeChild(element);
                                // }
                                g.class = "temp-cloud"
                                g.id = "remove-cloud-click"
                                this.parentElement.appendChild(g);

                                d3.selectAll("#remove-cloud-click")
                                    .on("click", function() { 
                                        removeGlobal()
                                    })
                                
                            });
                        }) // end of mouse enter clause
                    } //end of matrix u conditional

                    //make associated tree node light up
                    def.selectAll("path")
                        .on('mouseenter', function () {
                            d3.select("#color-boundary-"+identifier)
                                .style('stroke-width', 3)

                            var node = d3.select("#tree-node-"+identifier + " circle")
                        
                            node.style('stroke-width', 3)

                        })
                        .on('mouseleave', function () {
                            
                            d3.select("#color-boundary-"+identifier)
                                // .style('fill', fillcolor)
                                .style('stroke-width', 1)

                            var node = d3.select("#tree-node-"+identifier + " circle")
                        
                            node.style('stroke-width', 1)

                            //append to svg to change z-index
                            //var g = document.getElementById("temp-boundary-"+identifier)
                            //this.parentElement.parentElement.removeChild(g);
                        })
                        .on('click', function () {
                            removeGlobal()
                            var index = d3.select(this).attr("id")
                            index = index.replace(/\D/g,'');
                            if(index != "") {
                                if(index!= globalSliderIdx) {
                                    globalSliderIdx = index
                                    drawNewSlider(index)
                                } else {
                                    globalSliderIdx = -1
                                    drawNewSlider(-1)
                                }
                            } else {    // you clicked on a global boundary, should reset everything
                                globalSliderIdx = -1
                                drawNewSlider(-1)
                                //remove all globally drawn paths at the top 
                                removeGlobal()
                            }
                            
                            
                        });
                    // addStyle("svg .useuse:hover #color-boundary-"+identifier+ " {fill:red;opacity:0}")
                } //end of draw lines 
            } // end of for loop for each mutation

            if(globalSliderIdx != -1) {
                /** Determine threshold */
                var threshold = [];     //MarchingSquares requires an array

                if (thresholds.length < numMutations) {
                    threshold.push(thresholds[0]);
                } else {
                    threshold.push(thresholds[globalSliderIdx]);
                }

                var grid = [];
                temp = [];

                /* Turn column into grid (with n columns) */
                for(var i = 0; i < matrixF.length; i++) {
                    if(i!=0 && i%n == 0) {
                        grid.push(temp);
                        temp = [];
                    }
                    temp.push(matrix[i][globalSliderIdx]);
                }
                grid.push(temp);

                // create a border of zeros around the grid!! 
                // for creating a closed path and to visualize entire tumor

                var borders = [] 
                

                for(var row = 0; row < grid.length+2; row++)
                    borders[row] = new Array(grid.length+2).fill(0)


                for(var row = 1; row < borders.length-1; row++) {
                    for(var col = 2; col < borders[0].length-1; col++) {
                        borders[row][col] = grid[row-1][col-1]
                    }
                }

                var isoLines = [];
                MarchingSquaresJS
                  .isoLines(borders,
                            threshold,
                            {
                                polygons: false,
                                linearRing: false
                            }
                  )
                  .forEach(function(isolines, i) {
                    isoLines.push({
                      "coords": isolines,
                      "level": i + 1,
                      //"color": colorRanges[globalSliderIdx],
                      "color": colorIdx[indexToColorIdx[globalSliderIdx]],
                      "val": threshold[i]});
                  });

                drawLines(div, isoLines, thresholds, xRange, yRange, true);

                // helper function
                /** 
                 *   divId: div to draw the svg in
                 *   lines: isoLines to be drawn (given by grid values)
                 *   intervals: thresholds 
                 *   xs: rows (n)
                 *   ys: columns (m) 
                 */
            }

            /* draw grid points f
               after contour lines so that lines don't hide the grid */

            var gridPoints = svg.append("g");
            var tooltip = d3.select("body").append("div")   
                .attr("class", "tooltip")               
                .style("opacity", 0);

            // +2 for padding and borders 
            var xsection = (width)/(n+2)
            var ysection = (height)/(m+2)
            
            var samples = data.coordinates;
            var coordinates = []
            for(var x = 1; x < (n+2)-1; x++) {
                for(var y = 1; y < (m+2)-1; y++) {
                    coordinates.push({
                        "x": x,
                        "y": y
                    });
                }
            }

            //each grid line is horizontal 
            //for loop lays them down vertically 
            gridPoints.selectAll("circle.vertical")
                .data(coordinates)
                .enter().append("svg:circle")
                .attr("cx", function(d){
                    return d.x*xsection;
                })
                .attr("cy", function(d){
                    return d.y*ysection;
                })
                .attr("r", 1.75)
                .style("stroke", function(d) {
                    var x = d.x - 1 
                    var y = d.y - 1
                    var idx = x + y*(n)
                    if(samples[idx] != "Interpolated") {
                        return "#000"
                    } else {
                        return "#aaa"
                    }
                })
                .style("fill", function(d) {
                    var x = d.x - 1 
                    var y = d.y - 1
                    var idx = x + y*(n)
                    if(samples[idx] != "Interpolated") {
                        return "#000"
                    } else {
                        return "#aaa"
                    }
                })
                .on('mouseenter', function (d) {
                    var x = d.x - 1 
                    var y = d.y - 1
                    var idx = x + y*(n)
                    tooltip.transition()
                        .duration(200)
                        .style("opacity", 0.9)
                    // tooltip.html("Threshold: " + interval.toFixed(2))
                    tooltip.html("Sample: <br>" + samples[idx])
                        .style("left", (d3.event.pageX + 5)+"px")
                        .style("top", (d3.event.pageY - 35)+"px")
                    // console.log("x: " + x)
                    // console.log("y: " + y)
                    // console.log("sample: " + samples[idx])
                })
                .on('mouseleave', function (d) {
                    tooltip.transition()
                        .duration(500)
                        .style("opacity", 0)
                });


            //i=2 to prevent truncated grid points (vertical-wise)
            // for(var i = ysection; i < height-3; i+=ysection) {
            //     gridPoints.selectAll("circle.vertical")
            //         .data(yaxiscoorddata)
            //         .enter().append("svg:circle")
            //         .attr("cx", function(d){return d;})
            //         .attr("cy", i)
            //         .attr("r", 1.75)
            //         .on('mouseover', function (d) {
            //                 debugger
            //             })
            //     //.style("stroke", "rgb(6,120,155)")
            //     //.style("stroke-width", 8);  
            // }
        } // end of draw matrix function
    } // end of draw canvas
}   //end of visualize

/**
 *
 * Draw mutation legend
 * 
 * Draw checkbox mutations 
 *
 **/
// let flegend = document.getElementById("flegend");
// flegend.innerHTML = "<h3> Legend </h3>"

let checkboxes = document.getElementById("checkbox-mutations");
// let selections = document.getElementById("mutation-selection");
// selections.innerHTML = ""

let unselected_checkboxes = document.getElementById("unselected-mutations");
let selected_checkboxes = document.getElementById("selected-mutations");
unselected_checkboxes.innerHTML = ""
selected_checkboxes.innerHTML = ""

var mutToCluster = {}
//if labels are not provided 
if (mutLabels.length == 0) {
    console.log("Need to include labels for Frequency legend")

//if we have labels 
} else {
    // mutation legend
    // for (var i = 0; i < order.length; i++) {
    //     var label = ""
    //     if(order[i] != "") {
    //         for(var j in order[i]) {
    //             try {
    //                 if(!(order[i-1]).includes(order[i][j])) {
    //                     label += idxToLabel[order[i][j]] += ", " 
    //                 }    
    //             } catch(e) {
    //                 console.log("Error creating frequency color legend. ")
    //             }
                
    //         }
    //         label = label.substr(0, label.length-2)
    //     } else {
    //         label = "Unselected ancestral mutations"
    //     }
        
    //     if(label.length > 30) {
    //         label = label.substr(0,30)
    //     }
    // }

    //selection boxes 
    for (var cluster = 0; cluster < mutLabels.length; cluster++) {
        var currIdx = labelToIdx[mutLabels[cluster]]
        for (var mut = 0; mut < mutLabels[cluster].length; mut++) {
            mutToCluster[mutLabels[cluster][mut]] = currIdx
            if (selectedMutations.includes(currIdx)) {
                selected_checkboxes.innerHTML += "<option name = 'mutation' value='" + currIdx + "'>" + mutLabels[cluster][mut] + "</option>"
            } else {
                unselected_checkboxes.innerHTML += "<option name = 'mutation' value='" + currIdx + "'>" + mutLabels[cluster][mut] + "</option>"
            }
        }
    }
    unselected_checkboxes = sortDOM("unselected-mutations", true)
    selected_checkboxes = sortDOM("selected-mutations", true)

    checkboxes.innerHTML += "<br><input class = 'submitMutations' type='submit' value='Submit' onclick='pickMutations()'>"
}



// let checkboxes = document.getElementById("checkbox-mutations");
// let selections = document.getElementById("mutation-selection");
// selections.innerHTML = ""

function drawNewSlider(index) {
    if(index == -1) {
        let uThresholdContainer = document.getElementById("uThresholdContainer")
        uThresholdContainer.innerHTML = "<div class = 'threshold-row'><span class = 'threshold-value-label' id = 'uThresholdLabelAll'>10%</span>" +
                "<span class = 'mutationSliderLabel'>all</span>" +
                "<input type='range' min='1' max='100' value='10' id='uThresholdSliderAll' oninput='updateAllUThresholds(this)'>" +
                
                "</div>" 

    } else { 
        // u 
        let uThresholdContainer = document.getElementById("uThresholdContainer")
        uThresholdContainer.innerHTML = "" 

        var blank = " "

        var sliderIdx = "uThresholdSlider" + index 
        var thresholdValIdx = "uThresholdLabel" + index 

        var thresholdValue = uThreshold[index]*100;
        uThresholdContainer.innerHTML += "<div class = 'threshold-row' id = 'uThresholdRow" + index + "'><span class = 'threshold-value-label' id = '" + thresholdValIdx + "'>0.00</span> " +
            "<input type='range' min='1' max='100' value='"+ thresholdValue + "' id='" + sliderIdx + "' oninput='updateUThreshold(this)'>" +
            "<span class = 'mutationSliderLabel' style = 'background-color:" + colorIdx[indexToColorIdx[index]] +"'>  </span></div>"
            //"<span class = 'mutationSliderLabel' style = 'background-color:" + colorRanges[index] +"'>  </span></div>"

        //document.getElementById(sliderIdx).value = uThreshold[index]*100;
        document.getElementById(thresholdValIdx).innerHTML = uThreshold[index]*100 + "%";
    }
    
}

function sortDOM(divId, ascending) {
    var list = document.getElementById(divId);

    var items = list.childNodes;
    var itemsArr = [];
    for (var i in items) {
        if (items[i].nodeType == 1) { // get rid of the whitespace text nodes
            itemsArr.push(items[i]);
        }
    }

    itemsArr.sort(function(a, b) {
        if(ascending) {
            return a.innerHTML == b.innerHTML
              ? 0
              : (a.innerHTML > b.innerHTML ? 1 : -1);    
        } else {
            return a.innerHTML == b.innerHTML
              ? 0
              : (a.innerHTML < b.innerHTML ? 1 : -1);
        }
        
    });

    for (i = 0; i < itemsArr.length; ++i) {
      list.appendChild(itemsArr[i]);
    }
    return list
}


function removeGlobal() {
    while (document.getElementById("remove-cloud-click")) {
        var element = document.getElementById("remove-cloud-click")
        element.parentNode.removeChild(element);
    }
}

function removeGlobalCloudHover() {
    while (document.getElementById("remove-cloud")) {
        var element = document.getElementById("remove-cloud")
        element.parentNode.removeChild(element);
    }
}