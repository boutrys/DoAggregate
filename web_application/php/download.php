<?php

// Error reporting:
error_reporting(E_ALL^E_NOTICE);

// Including the connection file:
require('connect.php');


if(!is_bot()){
	$User_name_highlander = $_POST['User_name_highlander'];
	$Password_highlander = $_POST['Password_highlander'];
	$encrypted_password = md5($Password_highlander);
	
	
	$sql = "SELECT COUNT(*) FROM Highlander.users WHERE username = '" . $User_name_highlander . "' AND password = '" . $encrypted_password . "'";
	//$result = $con->query($sql);
	$result = mysqli_query($con, $sql);
	if (mysqli_num_rows($result) > 0) {
		// output data of each row
		while($row = mysqli_fetch_assoc($result)) {
			$res = $row["COUNT(*)"];
		}
	} else {
		$res = 0;
	}
	
	if($res == 1){
		$sql_database = "SELECT analysis, icon, sequencing_target FROM analyses WHERE reference = 'GRCh38' AND variant_caller = 'GATK' ORDER BY ordering";
		$result_database = mysqli_query($con, $sql_database);
		if (mysqli_num_rows($result_database) > 0) {
			$count = 0;
			$database = [];
			while($row = mysqli_fetch_assoc($result_database)) {
				$res = array('analysis' => $row["analysis"], 
							'sequencing_target' => $row["sequencing_target"]
						); 
				$database[$count] = $res;
				$count = $count + 1;
			}
		} else {
			$res = 0;
		}

		echo json_encode(array('error' => false, 'msg' => $database), JSON_FORCE_OBJECT);
	}else{
		echo json_encode(array('error' => true, 'msg' => $res));
	}
}


function is_bot()
{
	/* This function will check whether the visitor is a search engine robot */
	
	$botlist = array("Teoma", "alexa", "froogle", "Gigabot", "inktomi",
	"looksmart", "URL_Spider_SQL", "Firefly", "NationalDirectory",
	"Ask Jeeves", "TECNOSEEK", "InfoSeek", "WebFindBot", "girafabot",
	"crawler", "www.galaxy.com", "Googlebot", "Scooter", "Slurp",
	"msnbot", "appie", "FAST", "WebBug", "Spade", "ZyBorg", "rabaz",
	"Baiduspider", "Feedfetcher-Google", "TechnoratiSnoop", "Rankivabot",
	"Mediapartners-Google", "Sogou web spider", "WebAlta Crawler","TweetmemeBot",
	"Butterfly","Twitturls","Me.dium","Twiceler");

	foreach($botlist as $bot)
	{
		if(strpos($_SERVER['HTTP_USER_AGENT'],$bot)!==false)
		return true;	// Is a bot
	}

	return false;	// Not a bot
}
$con->close();
?>