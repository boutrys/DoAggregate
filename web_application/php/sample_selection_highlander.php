<?php

// Error reporting:
error_reporting(E_ALL^E_NOTICE);

// Including the connection file:
require('connect.php');


if(!is_bot()){
	$User_name_highlander = $_POST['User_name_highlander'];
	$Password_highlander = $_POST['Password_highlander'];
    $analysis = $_POST['analysis'];
	$encrypted_password = md5($Password_highlander);
	
	$options = "";
	$options_selected = "";
	//Retrieve all samples available for the analysis
	$sql_sample_list = "SELECT `value` FROM Highlander." . $analysis . "_possible_values WHERE `field` = 'sample' ORDER BY `value`";
	$result_sample_list = mysqli_query($con, $sql_sample_list);
	if (mysqli_num_rows($result_sample_list) > 0) {
		$count = 0;
		while($row = mysqli_fetch_assoc($result_sample_list)) {
			$options .= "<option>" . $row["value"] . "</option>";
			$options_selected .= "<option style='display: none;'>" . $row["value"] . "</option>";
			$res = array('value' => $row["value"]); 
			if($count == 0){
				$sample_list = array($row["value"]);
			}else{
				array_push($sample_list, $row["value"]);
			}
			$count = $count + 1;
		}
	} else {
		$res = 0;
	}
	//Retrieve all list of samples available for user 
	$sql_list ="SELECT DISTINCT(`key`) FROM Highlander.users_data WHERE `type` = 'VALUES' AND username = '" . $User_name_highlander . "' AND analysis = 'sample' ORDER BY `key`";
	$result_list = mysqli_query($con, $sql_list);
	if (mysqli_num_rows($result_list) > 0) {
		$count = 0;
		$name_list = [];
		$options_list = "";
		while($row = mysqli_fetch_assoc($result_list)) {
			$sql_name ="SELECT `value` FROM Highlander.users_data WHERE `type` = 'VALUES' AND username = '" . $User_name_highlander . "' AND analysis = 'sample' AND `key` = '" . $row["key"] . "' ORDER BY `value`;";
			$tmp_sample_list_res = mysqli_query($con, $sql_name);
			if (mysqli_num_rows($tmp_sample_list_res) > 0) {
				$count_tmp = 0;
				$tmp_sample_list = [];
				while($row_tmp = mysqli_fetch_assoc($tmp_sample_list_res)) {
					$index = array_search(strtoupper($row_tmp["value"]), $sample_list);
					if ($index !== false) {
						$res_tmp = array('value' => $index);
					}else{
						$res_tmp = array('value' => -1);
					}
					$tmp_sample_list[$count_tmp] = $res_tmp;
					$count_tmp = $count_tmp + 1;
				}
			} else {
				$tmp_sample_list = [];
			}
			$new = utf8_encode($row["key"]);
			$options_list .= "<option>" . $new . "</option>";
			$res = array('key' => $new, 'selected' => 0, 'value' => $tmp_sample_list);
			$name_list[$count] = $res;
			$count = $count + 1;
		}
	} else {
		$name_list = [];
	}

	echo json_encode(array('error' => false, 'list_of_sample' => $name_list, 'new_option' => $options, 'new_option_selected' => $options_selected, 'new_option_list' => $options_list), JSON_FORCE_OBJECT);

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