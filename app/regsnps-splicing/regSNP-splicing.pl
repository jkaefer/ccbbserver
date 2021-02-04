#!/usr/bin/perl
use warnings;
$|++;
   
use CGI qw(:cgi-lib :standard);  # Use CGI modules that let people read data passed from a form
use CGI::Carp 'fatalsToBrowser';

use File::Basename;
use URI;

use Cwd qw();
my $path = Cwd::cwd();



my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
my $ts = "$year-$mon-$mday-$hour-$min-$sec";
mkdir "../allJobs/Splicejob_$ts",0775 if !-e "../allJobs/Splicejob_$ts"; #changed 0775 to 0777
chmod 0777,"../allJobs/Splicejob_$ts"; #added 02-28-2019 since permissions are borked
unless(-e "../allJobs/Splicejob_$ts") { die "cannot dir"; }

my $job_url = "../allJobs/Splicejob_$ts/result.html";
#my $usr_access_result_url = "http://regsnps-splicing.ccbb.iupui.edu/query_jobs/job_$ts/result.html";
my $usr_access_result_url = "../allJobs/Splicejob_$ts/result.html";


my $wait_notice = "Your data has been successfully uploaded. </br> Please wait upto several minutes to finish prediction (About 2 seconds per SNV). </br></br> Or save the following <a href= \"$usr_access_result_url\" > link </a> and retrieve your result later.";

open(OUT, ">$job_url")||die "$!";
print OUT "
<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n
<html xmlns=\"http://www.w3.org/1999/xhtml\">\n
<head>\n
<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\" />\n

<script  type=\"text/javascript\">
function timeReload(timeoutPeriod){
    setTimeout(\"location.reload(true);\",timeoutPeriod);
}
</script>

<script  type=\"text/javascript\"   src=\"http://regsnps-splicing.ccbb.iupui.edu/source/code_package/js/sorttable.js\"></script>
</head>\n
<body bgcolor=\"white\" onload = \"JavaScript:timeReload(2000);\" >\n
<h3>$wait_notice<h3>\n
<p><h3><a href = \"$usr_access_result_url\">$usr_access_result_url</a></h3><p>
</body>\n
</html>\n
";
close OUT;



my $q = CGI->new;
print $q->redirect("$usr_access_result_url");

#print "Location: $usr_access_result_url\n\n";
print "Content-type: text/html\n\n";
print "<head><META HTTP-EQUIV=\"refresh\" content=\"0;URL=$usr_access_result_url\" /></head>\n";
#print "<h3>Your data has been successfully uploaded.</br> Please wait on this page. </br> Or save the <a href= \"$usr_access_result_url\" > link </a> and retrieve your result later.<h3>";

print "<body bgcolor=\"white\">
<h3>$wait_notice<h3>\n
<p><h3><a href = \"$usr_access_result_url\">$usr_access_result_url</a></h3><p>
</body>";

$old_fh = select(OUTPUT_HANDLE);
$| = 1;
select($old_fh);


$cgi_lib'writefiles = 1;
&ReadParse(%in);                 # This grabs the data passed by the form and puts it in an array

my $real_result = "../allJobs/Splicejob_$ts/result";
open(RE, ">$real_result")||die "$!";
print RE "
<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n
<html xmlns=\"http://www.w3.org/1999/xhtml\">\n
<head>\n
<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\" />\n

<script  type=\"text/javascript\">
function timeReload(timeoutPeriod){
        setTimeout(\"location.reload(true);\",timeoutPeriod);
}
</script>
<script  type=\"text/javascript\"   src=\"http://regsnps-splicing.ccbb.iupui.edu/source/code_package/js/sorttable.js\"></script>
</head>\n


<body bgcolor=\"white\" >\n
<style>
tr.tableHeader {
    background-color: #a0a09b;
    color: #ffffff;
    padding: 2px;
    font-family: Verdana,Arial,Helvetica,sans-serif;
}
tr.tableRowProb {
    background-color: #ff99b3;
    padding: 2px;
}
tr.tableRowPoss {
    background-color: #d0e3f0;
    padding: 2px;
}
tr.tableRowBenign {
    background-color: #e9e9e9;
    padding: 2px;
}
tr.tableHeader td p.bodytext {
    color: #ffffff;
}

th, td {
    padding: 5px;
    font-size: 10pt;
}
th {
    text-align: center;
}
td{
    text-align: center;
}
td#probable{
   background-color: #ff99b3; 
}
td#possible{
    background-color: #d0e3f0;
}
table#splicingSite {
    width: 100%;    
}
table#exonBody {
    width: 100%;
}
tr.tableRowBenign {
    background-color: #e9e9e9;
    padding: 2px;
}
tr.tableHeader td p.bodytext {
    color: #ffffff;
}
tr.tableRowEven {
    background-color: #f5f5f5;
    padding: 2px;
    font-family: Verdana,Arial,Helvetica,sans-serif;
}
tr.tableRowOdd {
    background-color: #e9e9e9;
    padding: 2px;
    font-family: Verdana,Arial,Helvetica,sans-serif;
}

ul, ol, dl \{ 
	padding\: 0\;
	margin\: 0\;
\}
h1, h2, h3, h4, h5, h6, p \{
	margin-top\: 0\;	
	padding-right\: 15px\;
	padding-left\: 15px\;
\}
a img \{ 
	border\: none\;
\}

a\:link \{
	color\: \#42413C\;
	text-decoration\: underline\;
\}
a\:visited \{
	color\: \#6E6C64\;
	text-decoration\: underline\;
\}
a\:hover, a\:active, a\:focus \{ 
	text-decoration\: none\;
\}

\.container \{
	width\: 960px\;
	background-color\: \#FFFFFF\;
	margin\: 0 auto\;
\}


table.sortable thead {
    cursor: default;
}


\#header \{
	height\: 100px\;
    text-align\:center\;
    padding\:5px\;
	background\: black url(\"http://regsnps-splicing.ccbb.iupui.edu/pic/banner.png\") center no-repeat\;	
\}
\#header_txt\{
	margin\: 30px 15px 40px 15px\;
	padding-left\:100px\;
	font-size\: 1\.8em\;
	font-weight\: bold\;
	float\:left\;
	\}
	\#header_logo \{
		padding\:1px\;
		float\:left\;
	\}
	\#sidebar1 \{
		float\: left\;
		width\: 180px\;
		background-color\: \#dfd9bd\;
		padding\: 10px\;
	\}
        \#sidebar2 \{  
            width\: 200px\;  
            height\: 100\%\;  
            background\: \#dfd9bd\;
            float\: right\;  
			padding\:10px\;
        \}  
        \#center \{  
            height\: 100\%\;  
            background\: FFEFD5\;  

            margin-left\:10px\;
			margin-right\:10px\;

			margin\: 0 0 0 210px\;     
			padding\: 10px\;     
        \}  
        \#footer \{  
            background\: grey\;  
            height\: 50px\;  
            width\: 100\%\;  
            position\: absolute\;  
            bottom\: 0\;  
        \}  

\.content ul, \.content ol \{ 
	padding\: 0 15px 15px 40px\;
\}

ul\.nav \{
	list-style\: none\;
	border-top\: 1px solid \#666\; 
	margin-bottom\: 15px\;
\}
ul\.nav li \{
	border-bottom\: 1px solid \#666\;
\}
ul\.nav a, ul\.nav a\:visited \{ 
	padding\: 5px 5px 5px 15px\;
	display\: block\; 
	width\: 160px\; 
	text-decoration\: none\;
	background-color\: \#dfd9bd\;
\}
ul\.nav a\:hover, ul\.nav a\:active, ul\.nav a\:focus \{ 
	background-color\: \#ADB96E\;
	color\: \#FFF\;
\}


</style>


<body bgcolor=\"white\">
  <div id = \"header\" >
	<img src=\"http://regsnps-splicing.ccbb.iupui.edu/pic/iu.png\" align=\"left\" name=\"center_logo\" width=\"100\" height=\"100\" id=\"Insert_logo\" style=\"display:block;\"/>
	<img src=\"http://regsnps-splicing.ccbb.iupui.edu/pic/center_logo.png\" width=\"96\" height=\"95\" align=\"right\" \"display:block;\"/>
	<div id=\"header_txt\" style=\"color:white\">regSNP-splicing</div>
    <!-- end .header --></div>
  <div id=\"sidebar1\">
    <ul class=\"nav\">
      <li><a href=\"http://regsnps-splicing.ccbb.iupui.edu/index.html\">Home</a></li>
      <li><a href=\"#\">FAQ</a></li>
      <li><a href=\"#\">Contact us</a></li>
    </ul>
    <p>News</p>
	<p style=\"font-size:15px\";>Jul,21,2015: Version 1.0 released.</p>
    <p>Related links</p>
    <ul style=\"list-style: none;\">
      <li style=\"font-size:15px\";><a href=\"http://www.hgmd.cf.ac.uk\">HGMD</a></li>
      <li style=\"font-size:15px\";><a href=\"http://www.1000genomes.org\">1000 Genomes Project</a></li>
      <li style=\"font-size:15px\";><a href=\"http://www.ncbi.nlm.nih.gov/clinvar/\">ClinVar</a></li>
    </ul>
	
    <!-- end .sidebar1 --></div>

 
<div id= \"center\">


";

my $att_opt = "../allJobs/Splicejob_$ts/checked_attributes";
open(CAT, ">$att_opt")||die "$!";
my @usr_att = param('attribute');
foreach my $att (@usr_att) {
#   if (param($att)) {
      print CAT "$att\n";
 #  }
}
close CAT;

$data_textbox = $in{"data_textbox"};             # Get the info from the textbox and assign to variable
my $usrdata = $data_textbox;
$usrdata =~ s/\r\n/;/g;
my @showdata = split(/;/, $usrdata);
my $shd = "";
my $shlimit = 9;
if($#showdata > $shlimit ){
    $shd = join('<br>', @showdata[0..$shlimit]);
    #print RE "<p>Your input data:</p><br>$shd<br>......";
}else{
    $shd = join('<br>', @showdata);
    #print RE "<p>Your input data:</p><br>$shd";
}
my $uploaded_fl = '';
my $input_desp = '';

if($data_textbox eq '') ## query is uploaded
{
	my $query = new CGI;
	my $safe_filename_characters = "a-zA-Z0-9_.-";
	my $filename = $query->param("data_file");
	my $tmpfilename = $query->tmpFileName($filename);

	my ( $name, $path, $extension ) = fileparse ( $filename, '..*' );
	$filename = $name . $extension;
	$filename =~ tr/ /_/;
	$filename =~ s/[^$safe_filename_characters]//g;
	#print RE "<p>The uploaded file is \"$filename\"</p>";

	my $upload_filehandle = $query->upload("data_file");
	$uploaded_fl = "../allJobs/Splicejob_$ts/$filename";
	$input_desp = $filename;
	#mkdir "../allJobs/Splicejob_$ts",0775 if !-e "../allJobs/Splicejob_$ts";
	
	open ( UPLOADFILE, ">$uploaded_fl" ) or die "$!";
	binmode UPLOADFILE;
	while ( <$upload_filehandle> ){
		print UPLOADFILE $_;
		#print "<p>$_\n</p>";
	}
	close UPLOADFILE;
	
}
else{ ## query is pasted in text-box
	$uploaded_fl = "../allJobs/Splicejob_$ts/query";
	$input_desp = 'query';
	#mkdir "../allJobs/Splicejob_$ts",0775 if !-e "../allJobs/Splicejob_$ts";

	open(IN, ">$uploaded_fl")||die "$!";
	print IN $data_textbox;
	close IN;
	
}

print RE "<br/>";
#print RE "success 0.     ";

close RE;

system("setsid perl bg_model.pl $uploaded_fl $input_desp $ts $real_result");










