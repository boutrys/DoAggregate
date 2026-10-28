<?php

require('configuration.php');

$con = mysqli_connect($db_host, $db_user, $db_pass, $db_database);

if (!$con) {
    echo "Error : cannot connect to database." . PHP_EOL;
    echo "Errno debug : " . mysqli_connect_errno() . PHP_EOL;
    echo "Error debug : " . mysqli_connect_error() . PHP_EOL;
    exit;
}

?>