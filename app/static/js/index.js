var options = {
  shouldSort: true,
  threshold: 0.4,
  maxPatternLength: 32,
  keys: [{
    name: 'user',
    weight: 0.5
  }, {
    name: 'prg',
    weight: 0.3
  }]
};
/*var jobs=[
  {
    "user":"Mark",
    "prg":"Intron",
    "date":"Tue 01 Dec 2020 01:02:58 PM EST"
  },
  {
    "user":"Mark",
    "prg":"Intron",
    "date":"Tue 01 Dec 2020 01:02:42 PM EST"
  },
  {
    "user":"Mary",
    "prg":"Splicing",
    "date":"Tue 01 Dec 2020 01:02:52 PM EST"
  },
  {
    "user":"Mary",
    "prg":"Exon",
    "date":"Tue 01 Dec 2020 01:02:54 PM EST"
  },
  {
    "user":"Jim",
    "prg":"Splicing",
    "date":"Tue 01 Dec 2020 01:02:53 PM EST"
  }
]

var fuse = new Fuse(jobs, options)
*/
//import files from './files.js';

var fuse = new Fuse(files, options)

$('.autocomplete').each(function() {
  var ac = $(this);
  
   ac.on('click', function(e) {
    console.log('boxclick');
    e.stopPropagation();
  })
  .on('focus keyup', search)
  .on('keydown', onKeyDown);
  
  var wrap = $('<div>')
    .addClass('autocomplete-wrapper')
    .insertBefore(ac)
    .append(ac);
  
    var list = $('<div>')
      .addClass('autocomplete-results')
      .on('click', '.autocomplete-result', function(e) {
        console.log('dropclick');
        e.preventDefault();
        e.stopPropagation();
        console.log('notselectindex')
        program=results[$(this).data('index')].prg;
        query_id=results[$(this).data('index')].date;
        console.log(program);
        console.log(query_id);
        //var url='app/allJobs/Exonjob_'+query_id+'/query';
        var url='';
        if(program=='Exon'){url='/exonRes?query_id='+query_id;}
        else if(program=='Intron'){url='/intronRes?query_id='+query_id;}
        else{url='/spliceRes?query_id='+query_id;}
        window.location.href = url;
        selectIndex($(this).data('index'), ac);
    })
    .appendTo(wrap);
});

$(document)
  .on('mouseover', '.autocomplete-result', function(e) {
    var index = parseInt($(this).data('index'), 10);
    if (!isNaN(index)) {
      $(this).attr('data-highlight', index);
    }
  })
  .on('click', clearResults);

function clearResults() {
  results = [];
  numResults = 0;
  $('.autocomplete-results').empty();
  $('.autocomplete-results').attr("style","margin-left:25%;width:0px;height:0px;overflow-y:hidden");
}

function selectIndex(index, autoinput) {
  if (results.length >= index + 1) {
    console.log('selectidex');
    console.log(results[index].prg);
    autoinput.val(results[index].prg);
    clearResults();
  }  
}

var results = [];
var numResults = 0;
var selectedIndex = -1;

function search(e) {
  
  if (e.which === 38 || e.which === 13 || e.which === 40) {
    return;
  }
  var ac = $(e.target);
  var list = ac.next();
  if (ac.val().length > 0) {
    results = _.take(fuse.search(ac.val()), 100);
    numResults = results.length;
    if(numResults > 0){
	    var divs = results.map(function(r, i) {
	    	console.log(r.prg);
	    	if(r.prg=='Exon'){
	    		return '<div class="autocomplete-result" data-index="'+ i +'">'
			     + '<div><b>'+ r.user +'</b> - '+ r.prg +' <img src="app/static/js/exon.png" style="width:20px;height:20px;border-radius: 20%;" > </div>'
			     + '<div class="autocomplete-location">'+ r.date +'</div>'
			     + '</div>';
		}
	    	else if(r.prg=='Intron'){
	    		return '<div class="autocomplete-result" data-index="'+ i +'">'
			     + '<div><b>'+ r.user +'</b> - '+ r.prg +' <img src="app/static/js/intron.png" style="width:20px;height:20px;border-radius: 20%;" > </div>'
			     + '<div class="autocomplete-location">'+ r.date +'</div>'
			     + '</div>';
		}
	    	else if(r.prg=='Splice'){
			return '<div class="autocomplete-result" data-index="'+ i +'">'
			     + '<div><b>'+ r.user +'</b> - '+ r.prg +' <img src="https://free.clipartof.com/450/1702-Free-Clipart-Of-A-Pair-Of-Scissors.jpg" style="width:20px;height:20px;" > </div>'
			     + '<div class="autocomplete-location">'+ r.date +'</div>'
			     + '</div>';
		}
	     });
	    
	    selectedIndex = -1;
	    list.html(divs.join(''))
	      .attr({
		"style":"margin-left:25%;width:150px;height:300px;overflow-x:hidden;overflow-y:auto", 
		"data-highlight":selectedIndex
		});

	    console.log(list);
	}
	else{clearResults();}
	   
  } else {
    numResults = 0;
    clearResults();
    list.empty();
    
  }
}

function onKeyDown(e) {
  console.log('keyDown');
  var ac = $(e.currentTarget);
  var list = ac.next();
  switch(e.which) {
    case 38: // up
      selectedIndex--;
      if (selectedIndex <= -1) {
        selectedIndex = -1;
      }
      list.attr('data-highlight', selectedIndex);
      break;
    case 13: // enter
      selectIndex(selectedIndex, ac);
      break;
    case 9: // enter
      selectIndex(selectedIndex, ac);
      e.stopPropagation();
      return;
    case 40: // down
      selectedIndex++;
      if (selectedIndex >= numResults) {
        selectedIndex = numResults-1;
      }
      list.attr('data-highlight', selectedIndex);
      break;

    default: return; // exit this handler for other keys
  }
  e.stopPropagation();
  e.preventDefault(); // prevent the default action (scroll / move caret)
}
