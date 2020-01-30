
/*
 * numMutations: int 
 * phyloTree: dictionary
 * colorRanges: array
 * subset: True or False for phylogeny tree format handling
 */

function drawTree(numMutations, phyloTree, colorRanges, order, treeLabels, clones, cloneLabels) {

  document.getElementById("tree").innerHTML = "<svg></svg>";
  
  // Create the input graph
  var g = new dagreD3.graphlib.Graph()
    .setGraph({})
    .setDefaultEdgeLabel(function() { return {}; });

  var treeIdx = {}

  //var index = 0;
  // use order[index] to make sure tree colors and IDs 
  // are loaded accurately
  for(var index = 0; index < numMutations; index++) {
    /**
      assumption:
      last node in ORDER (bfsorder; all the mutations in the clone)
      is the CURRENT node 
    **/

    g.setNode(index,  { 
      //label: treeLabels[order[index]], 
      labelType: "html",
      labelStyle: "font-size: 10px",
      label: clones[order[index]],
      id: "tree-node-u"+index,
      shape: "circle",
      width: 25,
      height: 25,
      style: "fill: " + colorIdx[indexToColorIdx[index]]});
      //style: "fill: " + colorRanges[index]});
    treeIdx[order[index]] = index 
  //  index += 1
  }

  //empty node 
  g.setNode("head",  { 
      //label: treeLabels[order[index]], 
      label: "      ",
      id: "tree-node-u",
      shape: "circle",
      style: "fill: #fff; stroke-dash-array: 5,5" });

  var rootLabel = ""
  if(idxToLabel[order[0]]) {

    var temp = idxToLabel[order[0]].replace(/\s+/g, '');
    var newLabel = temp.split(",")
    var displayLabels = []
    for(var i in temp) {
        if(htmlSelectedMutations.includes(newLabel[i])) {
          // message += newLabel[i] += ", "
          displayLabels.push(newLabel[i])
        } 
        i+=1
      } 
      if(displayLabels.length > 6) {
        for(var i = 0; i < 6; i++) {
          rootLabel += displayLabels[i] += " &nbsp;&nbsp;   <br>" + "&nbsp;&nbsp;   "
        }
        
        rootLabel+= "[" + (displayLabels.length-6) + " more] " 
      } else {
        for(var i in displayLabels) {
          rootLabel += displayLabels[i] += " &nbsp;&nbsp;   <br>" + "&nbsp;&nbsp;   "
        }
      }
  }


  //(Boolean_expression) ? do.somethingForTrue() : do.somethingForFalse();
  g.setEdge("head", 0, {
    labelType: "html",
    label: clones[order[0]] + rootLabel
  })

  for(var node in phyloTree) {
    for(var child in phyloTree[node]) {
      //debugger
      var label = ""
      var temp = treeLabels[phyloTree[node][child]].replace(/\s+/g, '');
      var newLabel = temp.split(",")
      var displayLabels = []
      var i = 0
      //debugger
      while(i < newLabel.length) {
        if(htmlSelectedMutations.includes(newLabel[i])) {
          // message += newLabel[i] += ", "
          displayLabels.push(newLabel[i])
        } 
        i+=1
      } 
      if(displayLabels.length > 6) {
        for(var i = 0; i < 6; i++) {
          label += newLabel[i] += " &nbsp;&nbsp;   <br>" + "&nbsp;&nbsp;   "
        }
        label+= "[" + (displayLabels.length-6) + " more] " 
      } else {
        for(var i in displayLabels) {
          label += newLabel[i] += " &nbsp;&nbsp;   <br>" + "&nbsp;&nbsp;   "
        }
      }

      g.setEdge(treeIdx[node], treeIdx[phyloTree[node][child]], 
        {
          labelType: "html",
          label: cloneLabels[phyloTree[node][child]] + " " + label,
          arrowhead: "normal"
        });
    }
  }

  try {
    g.nodes().forEach(function(v) {
      var node = g.node(v);
      // Round the corners of the nodes
      node.rx = node.ry = 5;
    });
  } catch(e) {
    console.log(e)
    console.log("Error occured setting nodes and edges in tree.js")
  }
  
  // Create the renderer
  var render = new dagreD3.render();

  var svg = d3.selectAll("#tree")
    .append("svg")
  var svgGroup = svg.append("g");
  var tooltip = d3.select("body").append("div")   
      .attr("class", "tooltip-tree")               
      .style("display", "none");
  render(d3.select("#tree svg g"), g);

  // Center the graph
  var xCenterOffset = (svg.attr("width") - g.graph().width) / 2;
  //svgGroup.attr("transform", "translate(" + xCenterOffset + ", 20)");
  svgGroup.attr("transform", "translate(20, 20)");
  
  var height = g.graph().height + 30
  var width = g.graph().width + 30
  svg.attr("height", height);
  svg.attr("width", width);
  svg.attr("id", "tree-svg")
  var viewbox = "0 0 " + width + " " +  height
  svg.attr("viewbox", viewbox)

  var br = document.getElementById("tree")
  br.innerHTML += " "

  // for(var i = 0; i < numMutations; i++) {
 //      var color=colorIdx[indexToColorIdx[i]]
 //      svg.selectAll(".nodes #tree-node-u"+i)
 //      .on('mouseover', function (d) {

 //        //append to svg to change z-index
 //          // this.parentElement.parentElement.appendChild(g);
 //          // d3.selectAll("#remove-cloud")
 //          //     .attr("stroke-width", 3)

 //          var cloudPath = d3.selectAll("#cloud-u"+d);
 //          var node = d3.select(this).select("circle"); 
 //          var boundary = d3.selectAll("#color-boundary-u"+d)
 //          if(cloudPath.classed('activated')) {
 //              node.style('stroke-width', 3)
 //              boundary.style('stroke-width', 3)
 //          }                   

 //          tooltip.transition()
 //           .duration(200)
 //           .style("display", "inline-block")

 //          var message = ""

 //          if(order[d] == "") {
 //            message = "Unselected ancestral mutations"
 //          } else {
 //            message = "Mutations: "
 //            for(var mut in order[d]) {
 //                if(message.length < 50) {
 //                    message += idxToLabel[order[d][mut]] + "<br>"
 //                } else {
 //                 message = message.substr(0, message.length-2)
 //                 message += "..."
 //                }
 //            }
 //            message = message.substr(0, message.length-2)
 //          }

 //          tooltip.html(message)
 //              .style("left", (d3.event.pageX + 18)+"px")
 //              .style("top", (d3.event.pageY - 28)+"px")
 //      })
 //      .on('mouseout', function (d) {
 //            var cloudPath = d3.selectAll("#cloud-u"+d);
 //            var node = d3.select(this).selectAll("circle"); 
 //            var boundary = d3.selectAll("#color-boundary-u"+d)
 //            if(cloudPath.classed('activated')) {
 //                node.style('stroke-width', 1)
 //                // boundary.style('fill', d3.rgb(colorRanges[d]))
 //                boundary.style('stroke-width', 1)
 //            }

 //            tooltip.transition()
  //          .duration(500)
  //          .style("display", "none")

 //      })
 //      .on('click', function (d) {
 //       alert("clicked!")
            
 //      });
 //    }

  // var symbolLegend = document.getElementById("symbols");
  // symbolLegend.innerHTML = "";
  // for (var index = 0; index < numMutations; index++) {
  //  var symbolLabel = treeLabels[order[index]];
  //  var showMore = ""
  //  if(symbolLabel.length > 30) {
  //    //<a href='javascript:void(0);' onclick='openMore(this); return false;'>more</a>
  //    showMore ="<span class = 'treelabel-showmore'>" + symbolLabel + "</span>"
  //    symbolLabel = symbolLabel.substr(0,30) + "..."
  //  }
 //    symbolLegend.innerHTML += "<p class = 'treelabel'><span>"+ cloneLabels[order[index]] + "</span>" + symbolLabel + showMore + " </p>"

 //  }

}


function addTreeInteractivity() {
  var tooltip = d3.select("body").append("div")   
      .attr("class", "tooltip-tree")               
      .style("display", "none");
// hovering effects 
    for(var i = 0; i < numMutations; i++) {
        var color = colorIdx[indexToColorIdx[i]]
        d3.select("#tree-svg").selectAll(".nodes #tree-node-u"+i)


        .on('mouseover', function () {
            var index = d3.select(this).attr("id")
            index = index.replace(/\D/g,'');

            var cloudPath = d3.selectAll("#cloud-u"+index);
            var node = d3.select(this).select("circle"); 

            var boundary = d3.selectAll("#color-boundary-u"+index)

            if(cloudPath.classed('activated')) {

                // Get associated cloudpath to HIGHLIGHT
                var test= document.getElementsByClassName("#cloud-selector-u"+index)[0]

                var g = test.cloneNode(true)
                g.class = "temp-cloud"
                g.id = "remove-cloud"

                test.parentElement.parentElement.appendChild(g);

                d3.selectAll("#remove-cloud")
                  .attr("stroke-width", 3)

                //debugger

                
                node.style('stroke-width', 3)
                //boundary.style('fill', d3.rgb(colorRanges[d]).darker(1.0))
                boundary.style('stroke-width', 3)
            }                   
                
            tooltip.transition()
              .duration(200)
              .style("display", "inline-block")

              var message = ""

              if(order[index] == "") {
                message = "Unselected ancestral mutations"
              } else {
                message = "Mutations: ";

              for(var mut in order[index]) {
                  var newLabel = idxToLabel[order[index][mut]]
                  newLabel = newLabel.replace(/\s+/g, '');
                  newLabel = newLabel.split(",")
                  newLabel = newLabel.filter(function (el) {
                    return el != "";
                  });
                  var displayLabels = []
                  var i = 0
                  //debugger
                  while(i < newLabel.length) {
                    if(htmlSelectedMutations.includes(newLabel[i])) {
                      // message += newLabel[i] += ", "
                      displayLabels.push(newLabel[i])
                    } 
                    i+=1
                  } 
                  if(displayLabels.length > 6) {
                    for(var i = 0; i < 6; i++) {
                      message += displayLabels[i] += ", "
                    }
                    message+= "[" + (displayLabels.length-6) + " more], " 
                  } else {
                    for(var i in displayLabels) {
                      message += displayLabels[i] += ", "
                    }
                    message=message.substr(0, message.length-1)
                  }
                }
                if(message.substr(message.length-2, message.length) ==", ") {
                  message = message.substr(0, message.length-2)
                }
              }

              tooltip.html(message)
                  .style("left", (d3.event.pageX + 18)+"px")
                  .style("top", (d3.event.pageY - 28)+"px")
          })

          .on('mouseout', function () {

            var d = d3.select(this).attr("id")
            d = d.replace(/\D/g,'');
            var cloudPath = d3.selectAll("#cloud-u"+d);
            var node = d3.select(this).selectAll("circle"); 
            var boundary = d3.selectAll("#color-boundary-u"+d)



            if(cloudPath.classed('activated')) {
                node.style('stroke-width', 1)
                // boundary.style('fill', d3.rgb(colorRanges[d]))
                boundary.style('stroke-width', 1)

                // cloudpath HIGHLIGHTING function: remove it 

                while (document.getElementById("remove-cloud")) {
                    var element = document.getElementById("remove-cloud")
                    element.parentNode.removeChild(element);
                }
            }

            tooltip.transition()
                .duration(500)
                .style("display", "none")

          })
          .on('click', function () {
                var d = d3.select(this).attr("id")
                d = d.replace(/\D/g,'');
                var node = d3.select(this).select("circle")

                            //the actual outlined clone path 
                var cloudPath = d3.selectAll("#cloud-u"+d);

                // drawContours(currMatrixU, index)

                if(cloudPath.classed('activated')) {
                    if(globalSliderIdx = d) {
                        globalSliderIdx = -1
                        drawNewSlider(-1)
                    }
                    cloudPath.style('display','none')
                    node.style('fill', '#fff')
                    node.style("stroke-width", 1)
                    activatedMutations["cloud-u"+d] = "deactivated"
                    cloudPath.classed('activated', false)
                } else {
                    //slider.style.display = 'inline'
                    cloudPath.style('display','inline')
                    node.style('fill', d3.rgb(colorIdx[indexToColorIdx[d]]))

                    activatedMutations["cloud-u"+d] = "activated"
                    cloudPath.classed('activated', true)
                }
                
          });
    }
}
