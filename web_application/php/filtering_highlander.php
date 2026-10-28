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

    // Retrieve all annotation for filtering (see script index.js --> variable annotation_for_filtering)
	$annotation_for_filtering = $_POST['annotation_for_filtering'];
    // if more annotation, put rest of code in for loop on annotation_for_filtering + make sure in Form.php the select object respect NameAnotation_full and in annotation_for_filtering in index.js NameAnotation
    // For example variant_type_full and variant_type

	$options = "";
	$options_selected = "";
	//Retrieve all values available for an annotation
    $sql = "SELECT `value` FROM Highlander." . $analysis . "_possible_values WHERE `field` = '" . $annotation_for_filtering . "' ORDER BY `value`";
	$result = mysqli_query($con, $sql);
	if (mysqli_num_rows($result) > 0) {
		$count = 0;
		while($row = mysqli_fetch_assoc($result)) {
            // Define a default value for the annotation
            if($row["value"] == "SNV"){
                $options .= "<option style='display: none;'>" . $row["value"] . "</option>";
                $options_selected .= "<option>" . $row["value"] . "</option>";
            }else{
                $options .= "<option>" . $row["value"] . "</option>";
                $options_selected .= "<option style='display: none;'>" . $row["value"] . "</option>";
            }
			$res = array('value' => $row["value"]); 
			if($count == 0){
				$gene_list = array($row["value"]);
			}else{
				array_push($gene_list, $row["value"]);
			}
			$count = $count + 1;
		}
	} else {
		$res = 0;
	}
		
	//Retrieve Another annotation for filtering
	echo json_encode(array('error' => false, 'new_option' => $options, 'new_option_selected' => $options_selected), JSON_FORCE_OBJECT);

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