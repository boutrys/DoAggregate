<?php

/* Database config: fill in the connection details of your Highlander database.
   Never commit real credentials. */

$db_host		= getenv('EXCALIBUR_DB_HOST') ?: 'localhost';
$db_user		= getenv('EXCALIBUR_DB_USER') ?: 'CHANGE_ME';
$db_pass		= getenv('EXCALIBUR_DB_PASS') ?: 'CHANGE_ME';
$db_database	= getenv('EXCALIBUR_DB_NAME') ?: 'Highlander';

/* End config */

?>
