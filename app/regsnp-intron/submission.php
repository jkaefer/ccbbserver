<!DOCTYPE html>

<?php
$query_id = $_GET["query_id"];
$page = "submission.php?query_id=".$query_id;
$sec = "10";
?>

<html lang="en">
<head>
  <title>regSNP-intron</title>
  <meta charset="utf-8">
  <meta http-equiv="refresh" content="<?php echo $sec?>;URL='<?php echo $page?>'">
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
          <li><a href=".">Home <span class="sr-only">(current)</span></a></li>
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
        <h1>Submission</h1>
<?php
echo "<p>Your query has been submitted. The job (id=" . trim($query_id, "query_") . ") is currently running. Once the job is done, the result link will be generated for you.</p><br><p>This page will be auto-refreshed to show the current job status.</p><br>";

ini_set("auto_detect_line_endings", true);
$log = file_get_contents("data/" . $query_id . "/log");
echo nl2br( $log );

if (strpos($log, 'final') !== false) {
  echo "<br><strong><p>You job is finished: ";
  echo "<a href='"."result.php?query_id=" . $query_id ."'>Result Page</a></p></strong>";
}
?>  
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
