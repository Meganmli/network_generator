

// this function executes our search via an AJAX call
function runSearch( term ) {
  // hide and clear the previous results, if any
  $('#results').hide();
  $('tbody').empty();
    
  // transforms all the form parameters into a string we can send to the server
  var frmStr = $('#gene_search').serialize();

  // call cgi script to generate the network
  $.ajax({
    url: './generate_network.cgi',
    data: frmStr,
    dataType: 'json',
    success: function(data, textStatus, jqXHR) {
      console.log( frmStr );
      console.log( data );
      processJSON(data);
    },
    error: function(jqXHR, textStatus, errorThrown){
        alert("Failed to generate network for " + frmStr + "!\nPlease enter a valid gene! textStatus: (" + textStatus + ")" +
        jqXHR.responseText + "\n and errorThrown: (" + errorThrown.Message + ")");
    }
  });
  
}

// processes a passed JSON structure representing interaction matches to create network and table
function processJSON( data ) {
    // return error message if present
    if ("message" in data) {
      alert( data.message );
      return
    }

    // keep track of row identifiers
    var next_row_num = 1; 
    var mouse_dict = {};

    // __________ NETWORK RESULTS __________
    // initialize graph network settings
    var cy = cytoscape({
        container: $('#cy'),
          style: [ // the stylesheet for the graph
          {
              selector: 'node',
              style: {
                  'background-color': '#666',
                  'label': 'data(label)',
                  'border-width': 'data(score)',
                  'width':'data(pubs)',
                  'height':'data(pubs)',
                  'line-color': '#666'
              }
          },
          {
              selector: 'edge',
              style: {
                'width': 3,

                // 'label': 'data(label)',
                'line-color': '#666',
                'target-arrow-color': '#666',
                'target-arrow-shape': 'triangle',
                'curve-style': 'bezier', // directed edges
                // "text-background-opacity": 1,
                // "color": "#666",
                // "text-background-color": "#fff"
              }
            },
            // Class selectors
            {
                'selector': '.red',
                'style': {
                    'background-color': '#f5a9c7',
                    'line-color': 'black'
                }
            }
          ],

          // initial viewport state:
          zoom: 1,
          pan: { x: 0, y: 0 },

          // interaction options:
          minZoom: 1,
          maxZoom: 2,
          autolock: false
    });

    // iterate through interactants
    $.each( data.interactants, function(i, item) {
      // collect genes
      mouse_dict[item.mouse_id] = {
        'human_symbol': item.human_symbol, 'human_entrezid':item.human_entrezid, 'pubs':item.pubs, 'score':item.score
      }

      // __________ NETWORK RESULTS __________
      // add nodes (genes)
      cy.add([
        { group: "nodes", data: {id: item.mouse_id, label: item.human_symbol+" ("+item.pubs+" papers)", symbol: item.human_symbol, pubs: item.pubs/5,  score: item.score}, 
                          position: { x: next_row_num*200, y: next_row_num*200 }, 
                          locked: false, grabbable: true }
      ]);
      // color input node red
      if (item.input === "true") {
        cy.nodes(`[id='${item.mouse_id}']`).addClass('red');
      } 
    });

    // iterate through interactions
    $.each( data.interactions, function(i, item) {
      // __________ NETWORK RESULTS __________
      // add edges (interaction)
      cy.add([
        { group: "edges", data: { id: item.geneA_interactant+item.geneB_interactant, label: "confidence score = "+mouse_dict[item.geneB_interactant]['score'], source: item.geneA_interactant, target: item.geneB_interactant} }
      ]);

      // __________ TABLE RESULTS __________
      var this_row_id = 'result_row_' + next_row_num++;
  
      // create a row and append it to the body of the table
      $('<tr/>', { "id" : this_row_id } ).appendTo('tbody');
      // add the row counter column
      $('<td/>', { "text" : next_row_num-1 } ).appendTo('#' + this_row_id);
      // add gene A information
      $('<td/>', { "text" : item.geneA_interactant } ).appendTo('#' + this_row_id); 
      $('<td/>', { "text" : mouse_dict[item.geneA_interactant]['human_entrezid'] } ).appendTo('#' + this_row_id);
      $('<td/>', { "text" : mouse_dict[item.geneA_interactant]['human_symbol'] } ).appendTo('#' + this_row_id);
      // add gene B information
      $('<td/>', { "text" : item.geneB_interactant } ).appendTo('#' + this_row_id);
      $('<td/>', { "text" : mouse_dict[item.geneB_interactant]['human_entrezid'] } ).appendTo('#' + this_row_id);
      $('<td/>', { "text" : mouse_dict[item.geneB_interactant]['human_symbol'] } ).appendTo('#' + this_row_id);

    });

    // set the span that displays the match count
    $('#match_count').text( data.match_count );

    // show the previously hidden result section
    $('#results').show();
    $('#cytoscape_results').show();

    let options = {
      name: 'preset',
      boundingBox: {x1: -100, y1: 0, w:100, h:100}, // constrain layout bounds; { x1, y1, x2, y2 } or { x1, y1, w, h }
      padding: 50,
      fit: true
    }
    cy.layout( options );
    cy.resize(); 

    console.log(cy.width())
    console.log(cy.height())
    console.log(cy.extent())

    cy.nodes().animate({
      position: { x: 700, y: 100 }
    }, {
      duration: 5000,
      complete: function(){
        console.log('Animation complete');
      }
    });
    
    console.log('Animating nodes...');
    setTimeout(function(){
      console.log('Stopping nodes animation');
      cy.nodes().stop();
    }, 2500);
}


// run our javascript once the page is ready
$(document).ready( function() {
    // runSearch when user clicks submit on our search form
    $('#submit').click( function() {
        runSearch();
        return false;  // prevents 'normal' form submission
    });
});

