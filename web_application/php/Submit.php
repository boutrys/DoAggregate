<?php
try {
    // Retrieving data from index.js
    $email = $_POST['email'];
    $pi = $_POST['pi'];
    $Highlander_import = $_POST['Highlander_import'];
    $do_coverage = $_POST['do_coverage'];
    $do_gene_list = $_POST['do_gene_list'];
    $mapqv = $_POST['mapqv'];
    $min_cover = $_POST['min_cover'];
    $remove_individu = $_POST['remove_individu'];
    $remove_position = $_POST['remove_position'];
    $min_nbr_var = $_POST['min_nbr_var'];
    $max_nbr_var = $_POST['max_nbr_var'];
    // Specific to Highlander import
    if($Highlander_import == "on"){
        $maf_filter = $_POST['maf_filter'];
        $HighlanderTable = $_POST['HighlanderTable'];
        $patient_list = $_POST['patient_list'];
        $control_list = $_POST['control_list'];
        $gene_list = $_POST['gene_list'];
        $filter_value = $_POST['filter_value'];
    }


    // Parameters
    $user = strstr($email, '@', true);
    $user = str_replace(".", "_", $user);
    $user = strtoupper($user);
    $folderName = strtoupper($pi . '_' . $user);
    $folder = '/storage/service/excalibur/input/' . date("Y") . '_' . date("m") . '_'. date("d") . '_'. date("H") . '_'. date("i") . '_'. date("s") . '_'. $folderName;
    $oldmask = umask(0);
    mkdir($folder, 0777, true);
    umask($oldmask);
    $excalibur = 'sudo -u boutrys /storage/service/excalibur/data/functions/submit_Excalibur.sh';

    
    // Input files creation for Excalibur
    if($Highlander_import == "off"){
        // Data were uploaded via the website
        $nbr_file = 0;
        if (isset($_FILES['file']['name'][0])) {
            foreach ($_FILES['file']['name'] as $keys => $values) {
                $tmpName = $_FILES['file']['tmp_name'][$keys];
                $file = $_FILES['file']['name'][$keys];     
                if ($keys == 0) {
                    $filename = $folder.'/patient_data.tsv';
                }else if ($keys == 1){
                    $filename = $folder.'/control_data.tsv';
                }else{
                    $filename = $folder.'/'. $file;
                }
                move_uploaded_file($tmpName, $filename);
                $nbr_file++;
            }
        }
        else {
            throw new Exception("missing files");
        }
    }else if($Highlander_import == "on"){
        // Using Highlander
        // patient
        $file_patient = fopen("$folder/list_patient.tsv", "w");
        $patient_list = explode(",", $patient_list);
        foreach ($patient_list as $patient) {
            fwrite($file_patient, $patient . "\n");
        }
        fclose($file_patient);
        // control
        $file_control = fopen("$folder/list_control.tsv", "w");
        $control_list = explode(",", $control_list);
        foreach ($control_list as $control) {
            fwrite($file_control, $control . "\n");
        }
        fclose($file_control);
        // gene list
        if ($do_gene_list == 'on') {
            $file_gene = fopen("$folder/gene_list.tsv", "w");
            $gene_list = explode(",", $gene_list);
            foreach ($gene_list as $gene) {
                fwrite($file_gene, $gene . "\n");
            }
            fclose($file_gene);
            $full_name_gene_file = $folder . '/gene_list.tsv';
        }
        // filters file
        $file_filter = fopen("$folder/filter.txt", "w");
        fwrite($file_filter, $filter_value);
        fclose($file_filter);
    }else{
        throw new Exception("Highlander_import variable should be set to 'on' or 'off'!");
    }
    
    
    //Prepare command line to launch
    $to_launch = $excalibur  . ' -E ' . $email;
    if($Highlander_import == "on"){
        $to_launch = $to_launch . ' -F ' . '"filter.txt"';
        $to_launch = $to_launch . ' -V ' . $HighlanderTable;
        // add 
        if($do_gene_list == 'on'){
            $to_launch = $to_launch . ' -b ' . '"' . $full_name_gene_file . '"';
        }
    }else{
        if($nbr_file == 3){
            $to_launch = $to_launch . ' -b ' . '"' . $filename . '"';
        }
    }
    $to_launch = $to_launch . ' -f ' . $min_nbr_var;
    $to_launch = $to_launch . ' -g ' . $max_nbr_var;
    if($Highlander_import == "on"){
        $to_launch = $to_launch . ' -i ' . $maf_filter;
    }
    $to_launch = $to_launch . ' -m ' . $mapqv;
    $to_launch = $to_launch . ' -n ' . $min_cover;
    if($remove_individu == "off"){
        $to_launch = $to_launch . ' -o';
    }
    if($remove_position == "off"){
        $to_launch = $to_launch . ' -p';
    }
    if($Highlander_import == "off"){
        $to_launch = $to_launch . ' -y';
    }
    if($do_coverage == "off"){
        $to_launch = $to_launch . ' -z';
    }
    $to_launch = $to_launch . ' ' . '"' . $folder . '/' . '"';

    if (empty($folder)) {
        header("Location: ../pages/Form.php");
        throw new Exception("nok");
    }
    else {
        // If you wanna the command line that was sent, see in input file, uncomment next three line
        $file_to_launch = fopen("$folder/to_launch.txt", "w");
        fwrite($file_to_launch, $to_launch);
        fclose($file_to_launch);
        //$output = shell_exec($to_launch);
        //header("Location: ../pages/Live.php");
        //exit();

        ob_implicit_flush(true);
        ob_end_flush();
        $descriptorspec = array(
           0 => array("pipe", "r"),   // stdin is a pipe that the child will read from
           1 => array("file", "$folder/php-out.txt", "a"),   // stdout is a pipe that the child will write to
           2 => array("file", "$folder/php-err.txt", "a")    // stderr is a pipe that the child will write to
        );
        flush();
        
        $process = proc_open($to_launch, $descriptorspec, $pipes, realpath('./'), array());
        stream_set_blocking($pipes[2], 0);
        if (is_resource($process)) {
            #has to use 'false !==' because if the character '0' is read, it's considered false and the php script exits or hangs forever
            while (false !== ($s = fgetc($pipes[1]))) {
                print $s;
                #ob_flush(); #erreur dans log apache = ob_flush(): failed to flush buffer.
                flush();
            }
        }
        
        echo "Errors:\r\n";
        fclose($pipes[1]);
        echo stream_get_contents($pipes[2]);
        echo "\r\n";
        fclose($pipes[2]);
        $exit=proc_close($process);
        echo "*exitcode^$exit*";
        header("Location: ../pages/Live.php");
        exit();
        
    }
}
catch (Exception $e) {
    echo json_encode($e->getMessage());
}

?>