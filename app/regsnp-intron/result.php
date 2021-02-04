<!DOCTYPE html>
<html lang="en">
<head>
  <title>regSNP-intron</title>
  <meta charset="utf-8">
  <script src="https://code.jquery.com/jquery-3.1.0.min.js" integrity="sha256-cCueBR6CsyA4/9szpPfrX3s49M9vUU5BgtiJj06wt/s=" crossorigin="anonymous">
  </script>
  <!-- Latest compiled and minified CSS -->
  <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous">

  <!-- Optional theme -->
  <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap-theme.min.css" integrity="sha384-rHyoN1iRsVXV4nD0JutlnGaslCJuC7uwjduW9SVrLvRYooPp2bWYgmgJQIXwl/Sp" crossorigin="anonymous">

  <!-- Latest compiled and minified JavaScript -->
  <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" crossorigin="anonymous"></script>
  <link rel="stylesheet" href="resources/css/mystyle.css"/>
  <link rel="icon" href="resources/images/iu_tab.jpg">

  <link rel="stylesheet" href="https://cdn.datatables.net/1.10.12/css/dataTables.bootstrap.min.css">
  <script src="https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min.js"></script>
  <script src="https://cdn.datatables.net/1.10.12/js/dataTables.bootstrap.min.js"></script>

  <script language="javascript" src="//www.biodalliance.org/release-0.13/dalliance-compiled.js"></script>
  <script>
  var query_id = "<?php echo $_GET["query_id"] ?>";
  var data_file = "data/" + query_id + "/output/snp.prediction.json";
  $(document).ready(function() {
	$.fn.dataTable.ext.errMode = 'throw';
    $('#resultTable').DataTable( {
	    "processing": true,
        "scrollX": true,
		"ajax": data_file,
        "columns": [
            { "data": "#chrom" },
            { "data": "pos" },
            { "data": "ref" },
            { "data": "alt" },
            { "data": "disease" },
            { "data": "prob",
              "render": function ( data, type, row ) {
                 return parseFloat(data).toFixed(2);;
              }
            },
            { "data": "tpr",
              "render": function ( data, type, row ) {
                 return parseFloat(data).toFixed(2);;
              }
            },
            { "data": "fpr",
              "render": function ( data, type, row ) {
                 return parseFloat(data).toFixed(2);;
              }
            },
            { "data": "splicing_site" },
            { "data": "name" },
            { "data": "strand" }
        ]
    } );
    
    var table = $('#resultTable').DataTable();

    $('#resultTable tbody').on( 'click', 'tr', function () {
        if ( $(this).hasClass('info') ) {
            $(this).removeClass('info');
        }
        else {
            table.$('tr.info').removeClass('info');
            $(this).addClass('info');
		}
		var record = table.row(this).data();
        console.log(record['#chrom'].substring(3));
        browser.setLocation(record['#chrom'].substring(3), record['pos'] - 50, record['pos'] + 50);
    } );

    var browser =  new Browser({
      chr:          '22',
      viewStart:    30700000,
      viewEnd:      30900000,

      coordSystem: {
        speciesName: 'Human',
        taxon: 9606,
        auth: 'GRCh',
        version: '37',
        ucscName: 'hg19'
      },

      sources:     [{name:                 'Genome',
                     twoBitURI:            '//www.biodalliance.org/datasets/hg19.2bit',
                     tier_type:            'sequence'},
                    {name:                 'Genes',
                     desc:                 'Gene structures from GENCODE 19',
                     bwgURI:               '//www.biodalliance.org/datasets/gencode.bb',
                     stylesheet_uri:       '//www.biodalliance.org/stylesheets/gencode.xml',
                     collapseSuperGroups:  true,
                     trixURI:              '//www.biodalliance.org/datasets/geneIndex.ix'},
                    {name:                 'Repeats',
                     desc:                 'Repeat annotation from Ensembl',
                     bwgURI:               '//www.biodalliance.org/datasets/repeats.bb',
                     stylesheet_uri:       '//www.biodalliance.org/stylesheets/bb-repeats.xml'},
                    {name:                 'Conservation',
                     desc:                 'Conservation', 
                     bwgURI:               '//regsnps-intron.ccbb.iupui.edu/annotation/hg19.100way.phyloP100way.bw',
                     noDownsample:         true}],

    });


  } );
  $(function(){
    $.ajax({
      type: 'HEAD',
      url: data_file,
      success: function() {
        $("#message").addClass("alert-success").html("<strong>The full result with all the features can be downloaded here: <a download href='" + "data/" + query_id + "/output/snp.prediction.txt" + "' target='_blank'>snp.prediction.txt</a>.</strong>");
      },
      error: function() {
        $("#message").addClass("alert-danger").html("<strong>The job may still be running or an error has occurred. Please check the <a href='submission.php?query_id=" + query_id + "'>submission page</a>.</strong>");
      }
    });
  });
  </script>
</head>

<body>
<div class="container-fluid">
  <nav class="navbar navbar-default">
    <div class="container-fluid">
      <!-- Brand and toggle get grouped for better mobile display -->
      <div class="navbar-header">
        <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#myNavbar" aria-expanded="false">
          <span class="sr-only">Toggle navigation</span>
          <span class="icon-bar"></span>
          <span class="icon-bar"></span>
          <span class="icon-bar"></span>
        </button>
        <a class="navbar-brand" href="."><img style="width: 52px;" src="resources/images/iu_tab.jpg">regSNP-intron</a>
      </div>
  
      <!-- Collect the nav links and other content for toggling -->
      <div class="collapse navbar-collapse" id="myNavbar">
        <ul class="nav navbar-nav">
          <li><a href=".">Home </a></li>
		  <li><a href="document.html">Document</a></li>
          <li><a href="about.html">About</a></li>
          <li class="dropdown">
            <a href="#" class="dropdown-toggle" data-toggle="dropdown" role="button" aria-haspopup="true" aria-expanded="false">Tools <span class="caret"></span></a>
            <ul class="dropdown-menu">
              <li><a href="http://watson.compbio.iupui.edu/regSNP-splicing/" target="_blank">regSNP-splicing</a></li>
              <li><a href="http://watson.compbio.iupui.edu/ExonImpact/" target="_blank">ExonImpact</a></li>
              <li role="separator" class="divider"></li>
              <li><a href="http://annovar.openbioinformatics.org/en/latest/" target="_blank">ANNOVAR</a></li>
            </ul>
          </li>
        </ul>
        <ul class="nav navbar-nav navbar-right">
          <li><a href="https://github.com/yliulab/regsnp_intron" target="_blank">GitHub</a></li>
        </ul>
      </div><!-- /.navbar-collapse -->
    </div><!-- /.container-fluid -->
  </nav>
  
  <div class="container-fluid">
    <div class="row content">
      <div class="col-sm-3 sidenav">
        <h4>Resources</h4>
        <ul class="nav nav-pills nav-stacked">
          <li><a href="http://www.hgmd.cf.ac.uk/ac/index.php" target="_blank">HGMD</a></li>
          <li><a href="http://www.1000genomes.org/" target="_blank">1000 Genomes</a></li>
          <li><a href="http://www.ncbi.nlm.nih.gov/clinvar/" target="_blank">ClinVar</a></li>
        </ul><br>
      </div>
  
      <div class="col-sm-9">
		<h1>regSNP-intron</h1>
		<!-- <p id="message"></p> -->
		<div id="message" class="alert"></div>

        </p>
        <table id="resultTable" class="table table-striped table-bordered table-hover" cellspacing="0" width="100%">
          <thead>
              <tr>
                  <th>Chrom</th>
                  <th>Pos</th>
                  <th>Ref</th>
                  <th>Alt</th>
                  <th>Disease</th>
		  <th>Prob</th>
                  <th>TPR</th>
                  <th>FPR</th>
                  <th>Splicing_site</th>
                  <th>Name</th>
                  <th>Strand</th>
              </tr>
          </thead>
        </table>
        <div id="svgHolder"></div>
      </div>
    </div>
  </div>
  
  <footer class="container-fluid">
    <small>
    <p class="text-center">&copy; 2016 CCBB</p>
    <address><p class="text-center"><a href="mailto:yunliu@iupui.edu">Contact Us</a></p></address>
    </small>
  </footer>
</div>
</body>
</html>
