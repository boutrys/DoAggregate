<?php

$keyword = " pipeline";
#$jobs = shell_exec('squeue -h -o "%.100j %.15u %.2t %.10M %.12P %i"| grep ' . $keyword . ' | sort -k6');
$jobs = shell_exec('squeue -h -o "%.100j %.15u %.2t %.10M %.12P %i" | grep ' . $keyword . ' | sort -k6');
$arr = explode("\n", $jobs);
echo json_encode($arr);

?>