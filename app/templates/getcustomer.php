<?php
$mysqli = new mysqli("localhost","root","12345","db20");
if($mysqli->connect_error) {
  exit('Could not connect');
}

$sql = "SELECT pid FROM prot_data where pid=1";

$stmt = $mysqli->prepare($sql);
#$stmt->bind_param("s", $_GET['q']);
$stmt->execute();
$stmt->store_result();
$stmt->bind_result($pid);
$stmt->fetch();
$stmt->close();

echo "<p>" . $pid . "</p>";



?>
