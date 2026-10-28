<?php

$fileData = '';
$message = '';

try {
  if (isset($_FILES['file']['name'][0])) 
{

  $pi = $_POST['pi'];
  $email = $_POST['email'];
  $user = strstr($email, '@', true);
  $user = str_replace(".", "_", $user);
  $folderName = strtoupper($pi . '_' . $user);
  $folder = '/storage/service/excalibur/input/' . date("Y") . '_' . date("m") . '_'. date("d") . '_'. date("h") . '_'. date("i") . '_'. date("s") . '_'. $folderName;    
  $file = $_FILES['file']['tmp_name'];
  $path_parts = pathinfo($file);
  $fileName = $path_parts['filename'];

  $fileData .= $fileName . ',';
    
  $fileData = substr($fileData, 0, -1);
  echo json_encode(array('error' => false, 'msg' => $fileData));
  }
}
catch (Exception $e) {
  echo json_encode(array('error' => true, 'msg' => $e->getMessage()));
}
?>