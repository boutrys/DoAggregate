$(document).ready(function () {
    $(document).keypress(
        function(event){
        if (event.which == '13') {
            event.preventDefault();
        }
    });

    function Highlander_import() {
        var new_val = 'on';
        $('#Highlander_import').val(new_val).trigger('change');
        var new_val2 = 'off';
        $('#manual_import').val(new_val2).trigger('change');
        $('#clean_data_switch_option').show();
    }
    function manual_import() {
        var new_val = 'off';
        $('#Highlander_import').val(new_val).trigger('change');
        var new_val2 = 'on';
        $('#manual_import').val(new_val2).trigger('change');
        if($('#do_coverage').val() == 'on'){
            var newval_cov = 'off';
            $('#do_coverage').val(newval_cov).trigger('change');
        }
        $('#clean_data_switch_option').hide();
    }
    function do_coverage() {
        if($('#do_coverage').val() === 'on'){
            var new_val = 'off';
            $('#do_coverage').val(new_val).trigger('change');
        }else{
            var new_val = 'on';
            $('#do_coverage').val(new_val).trigger('change');
        }
    }
    function gene_list() {
        $('#geneFile').collapse('hide');
        if($('#gene_list').val() === 'on'){
            var new_val = 'off';
            $('#gene_list').val(new_val).trigger('change');
        }else{
            var new_val = 'on';
            $('#gene_list').val(new_val).trigger('change');
        }
    }
    function remove_individu() {
        if($('#remove_individu').val() === 'on'){
            var new_val = 'off';
            $('#remove_individu').val(new_val).trigger('change');
        }else{
            var new_val = 'on';
            $('#remove_individu').val(new_val).trigger('change');
        }
    }
    function remove_position() {
        if($('#remove_position').val() === 'on'){
            var new_val = 'off';
            $('#remove_position').val(new_val).trigger('change');
        }else{
            var new_val = 'on';
            $('#remove_position').val(new_val).trigger('change');
        }
    }


    $('#Highlander_import').click(Highlander_import)
    $('#manual_import').click(manual_import)
    $('#do_coverage').click(do_coverage)
    $('#gene_list').click(gene_list)
    $('#remove_individu').click(remove_individu)
    $('#remove_position').click(remove_position)


    //Next (choose inbetween loading data via website or using Highlander)
    $("#setOptions").on('click', function (event) {
        // add tippy info 
        const selectElement = document.querySelector('#list_sample_full');
        selectElement.addEventListener('click', function() {
            const selectedOption = selectElement.options[selectElement.selectedIndex].text;
            tippy('#list_sample_full', {
                content: `${selectedOption}`
            });
        });
        const selectElement_gene = document.querySelector('#list_gene_full');
        selectElement_gene.addEventListener('click', function() {
            const selectedOption_gene = selectElement_gene.options[selectElement_gene.selectedIndex].text;
            tippy('#list_gene_full', {
                content: `${selectedOption_gene}`
            });
        });
        
        event.preventDefault()
        event.stopPropagation()
        var email_regex = /^[a-zA-Z0-9._-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,4}$/i;
        if (!email_regex.test($("#email").val())) {
        $('#email').addClass("is-invalid");
        return; 
        }
        setCookie("email",$("#email").val(),180);
        setCookie("group",$("#pi option:selected").val(),180);
        $("#options").addClass("visually-hidden ");
        if($('#Highlander_import').val() == "off"){
        //If import your data files
        $("#submit_via_website_card").removeClass("visually-hidden ");
        $("#submit-part-trigger").removeClass("disabled");
        }else{
        $("#information-part-trigger").addClass("visually-hidden ");
        $("#submit-part-trigger").addClass("visually-hidden ");
        $("#Login_Highlander").removeClass("visually-hidden ");
        }
        $("#information-part-trigger").addClass("disabled");
    });


    //  PART FOR SUBMISSION USING HIGHLANDER 
    var nbr_database=0; //note that it will be equal to true number of database + 1
    var database_id_isOn="0";
    var selected_analysis_name="";
    var val_select_private_list_individu_collapse = false;
    var sample_list_user = [];
    var selected_list_samples_as_patient = [];
    var selected_list_samples_as_control = [];
    var val_select_private_list_gene_collapse = false;
    var gene_list_user = [];

    //Data pretreatment and upload via Highlander
    function fillData_loginHighlander() {
        var formData = new FormData();
        formData.append('User_name_highlander', $('#User_name_highlander').val());
        formData.append('Password_highlander', $('#Password_highlander').val());
        return formData;
    }

    function fillData_sample_selection() {
        var formData = new FormData();
        formData.append('User_name_highlander', $('#User_name_highlander').val());
        formData.append('Password_highlander', $('#Password_highlander').val());
        formData.append('analysis', selected_analysis_name);
        return formData;
    }

    function fillDataFinal_highlander() {
        var formData = new FormData();
        formData.append('email', $('#email').val());
        formData.append('pi', $('#pi option:selected').val());
        formData.append('Highlander_import', $('#Highlander_import').val());
        formData.append('do_coverage', $('#do_coverage').val());
        formData.append('do_gene_list', $('#gene_list').val());
        formData.append('mapqv', $('#mapqv').val());
        formData.append('min_cover', $('#min_cover').val());
        formData.append('remove_individu', $('#remove_individu').val());
        formData.append('remove_position', $('#remove_position').val());
        formData.append('min_nbr_var', $('#min_nbr_var').val());
        formData.append('max_nbr_var', $('#max_nbr_var').val());
        //Specific to Highlander import
        formData.append('User_name_highlander', $('#User_name_highlander').val());
        formData.append('Password_highlander', $('#Password_highlander').val());
        formData.append('analysis', selected_analysis_name);
        formData.append('maf_filter', $('#maf_filter').val());
        formData.append('HighlanderTable', selected_analysis_name);
        //Retrieve user's patient selection
        if(val_select_private_list_individu_collapse == false){
        var selectElement_id = document.getElementById("patient_list_selected");
        }else{
        var selectElement_id = document.getElementById("list_patient_selected");
        }
        var selectElement = selectElement_id.getElementsByTagName("option");
        let selected_patient = [];
        for (let i = 0; i < selectElement_id.options.length; i++) {
        if (selectElement[i].style.display == '') {
            selected_patient.push(selectElement[i].value);
        }
        }

        formData.append('patient_list', selected_patient);
        //Retrieve user's control selection
        if(val_select_private_list_individu_collapse == false){
        var selectElement_id = document.getElementById("control_list_selected");
        }else{
        var selectElement_id = document.getElementById("list_control_selected");
        }
        var selectElement = selectElement_id.getElementsByTagName("option");
        let selected_control = [];
        for (let i = 0; i < selectElement_id.options.length; i++) {
        if (selectElement[i].style.display == '') {
            selected_control.push(selectElement[i].value);
        }
        }
        formData.append('control_list', selected_control);
        //Retrieve user's gene selection
        if($('#gene_list').val() == 'on'){
        if(val_select_private_list_gene_collapse == false){
            var selectElement_id = document.getElementById("gene_list_selected");
        }else{
            var selectElement_id = document.getElementById("list_gene_selected");
        }
        var selectElement = selectElement_id.getElementsByTagName("option");
        let selected_gene = [];
        for (let i = 0; i < selectElement_id.options.length; i++) {
            if (selectElement[i].style.display == '') {
            selected_gene.push(selectElement[i].value);
            }
        }

        formData.append('gene_list', selected_gene);
        }
        var filter_value = "(";
        if($('#Filter_pass').val() == 'PASS'){
        filter_value += "`filters`='PASS'";
        filter_value += " AND ";
        }
        filter_value += "`consensus_prediction` >= "; 
        filter_value += String($('#consensus_prediction_filter').val());
        filter_value += " AND (`gnomad_wes_popmax_af` < ";
        filter_value += String($('#maf_filter').val());
        filter_value += " OR `gnomad_wes_popmax_af` IS NULL))";
        var selectElement_id = document.getElementById("variant_type_selected");
        var selectElement = selectElement_id.getElementsByTagName("option");
        let selected_variant_type = [];
        for (let i = 0; i < selectElement_id.options.length; i++) {
        if (selectElement[i].style.display == '') {
            selected_variant_type.push(selectElement[i].value);
        }
        }
        filter_value += " AND `variant_type` IN (";
        filter_value += `'${selected_variant_type.join("', '")}'`;
        filter_value += ")";
        formData.append('filter_value', filter_value);
        
        return formData;
    }

    function LoginHighlanderFunction(form_data) {
        $.ajax({
        url: "../php/download.php",
        method: "POST",
        data: form_data,
        contentType: false,
        cache: false,
        processData: false,
        success: function (data) {
            var resp = $.parseJSON(data);
            if (resp.error) {
            //login not successfull
            $('#card_drop_loginHighlander').addClass("border-danger");
            $('#uploaded_selection_error_loginHighlander').addClass("is-invalid");
            $('#uploaded_selection_error_loginHighlander').show();

            }else{
            //login successfull
            $('#card_drop_loginHighlander').removeClass("border-danger");
            $('#uploaded_selection_error_loginHighlander').removeClass("is-invalid");
            $('#uploaded_selection_error_loginHighlander').hide();
            event.preventDefault()
            event.stopPropagation()
            $("#Login_Highlander").addClass("visually-hidden ");
            $("#Login_Highlander").addClass("disabled");
            $("#setOptions").addClass("disabled");

            //Page Database selection
            var count_database = 1;
            $("#Database_selection").removeClass("visually-hidden ");
            for (var i = 0, row; row = resp.msg[i]; i++){
                var tmp_name_database = "#DB_";
                tmp_name_database = tmp_name_database.concat(count_database);
                $(tmp_name_database).removeClass("visually-hidden ");
                var tmp_database_label = "#database_label_";
                tmp_database_label = tmp_database_label.concat(count_database);
                $(tmp_database_label).text(row.analysis);
                if(count_database == 1){
                //By default the first analysis is checked and taken as selected analysis
                database_id_isOn = "database_1";
                selected_analysis_name = row.analysis;
                }
                count_database = count_database + 1;
            }
            nbr_database = count_database;
            }
            
        },
        error: function (data) {
            console.log("AJAX ERR Login Highlander");
            console.log(data);
        }
        });
    }

    function sampleSelectionFunction(form_data) {
        $.ajax({
        url: "../php/sample_selection_highlander.php",
        method: "POST",
        data: form_data,
        contentType: false,
        cache: false,
        processData: false,
        xhr: function(){
            //Upload progress
            var xhr = new window.XMLHttpRequest();
            xhr.upload.addEventListener("progress", function(evt){
            if (evt.lengthComputable) {  
                var percentComplete = Math.round((evt.loaded / evt.total) * 100);
                $('#progress_bar_sample_selection').width(percentComplete + '%');
                $('#progress_bar_sample_selection').html(percentComplete + '% completed');  
            }
            }, false); 
            return xhr;
        },
        success: function (data) {
            var resp = $.parseJSON(data);
            if (resp.error) {
            //Not working
            }else{
            //is ok
            $("#sample_list_full").html(resp.new_option);
            $("#patient_list_selected").html(resp.new_option_selected);
            $("#control_list_selected").html(resp.new_option_selected);
            $("#list_patient_selected").html(resp.new_option_selected);
            $("#list_control_selected").html(resp.new_option_selected);
            $("#list_sample_full").html(resp.new_option_list);

            sample_list_user = resp.list_of_sample;
            var tempObj = JSON.parse(JSON.stringify(resp.list_of_sample));
            selected_list_samples_as_patient = tempObj;
            var tempObj2 = JSON.parse(JSON.stringify(resp.list_of_sample));
            selected_list_samples_as_control = tempObj2;
            for (i = 0; i < Object.keys(selected_list_samples_as_patient).length; i++) {
                selected_list_samples_as_patient[i].value = 0;
                selected_list_samples_as_control[i].value = 0;
            }

            $("#Database_selection").addClass("visually-hidden ");
            $("#Database_selection").addClass("disabled");
            event.preventDefault()
            event.stopPropagation()
            $("#Sample_selection").removeClass("visually-hidden ");
            }
            
        },
        error: function (data) {
            console.log("AJAX ERR sample selection Highlander");
            console.log(data);
        }
        });
    }

    function geneSelectionFunction(form_data) {
        $.ajax({
        url: "../php/gene_selection_highlander.php",
        method: "POST",
        data: form_data,
        contentType: false,
        cache: false,
        processData: false,
        xhr: function(){
            //Upload progress
            var xhr = new window.XMLHttpRequest();
            xhr.upload.addEventListener("progress", function(evt){
            if (evt.lengthComputable) {  
                var percentComplete = Math.round((evt.loaded / evt.total) * 100);
                $('#progress_bar_gene_selection').width(percentComplete + '%');
                $('#progress_bar_gene_selection').html(percentComplete + '% completed');  
                var msg = percentComplete + '% completed';
            }
            }, false); 
            return xhr;
        },
        success: function (data) {
            var resp = $.parseJSON(data);
            if (resp.error) {
            //Not working
            }else{
            //is ok
            $("#gene_list_full").html(resp.new_option);
            $("#gene_list_selected").html(resp.new_option_selected);
            $("#list_gene_selected").html(resp.new_option_selected);
            $("#list_gene_full").html(resp.new_option_list);

            gene_list_user = resp.list_of_gene;

            $("#Sample_selection").addClass("visually-hidden ");
            $("#Sample_selection").addClass("disabled");
            event.preventDefault()
            event.stopPropagation()
            $("#gene_selection").removeClass("visually-hidden");
            }
            
        },
        error: function (data) {
            console.log("AJAX ERR gene selection Highlander");
            console.log(data);
        }
        });
    }

    function database(event) {
        var n = event.data.param;
        var new_val = 'on';
        var name_database = '#database_';
        var name_database_to_change = name_database.concat(n);
        database_id_isOn = name_database_to_change;
        var tmp_name_analysis = "#database_label_";
        var tmp_name_analysis_to_select = tmp_name_analysis.concat(n);
        selected_analysis_name = $(tmp_name_analysis_to_select).text();
        $(name_database_to_change).val(new_val).trigger('change');
        var new_val2 = 'off';
        var count_database = 1;
        for (var i = 1; i < nbr_database; i++){
        if(i != n){
            var other_database = name_database.concat(i);
            $(other_database).val(new_val2).trigger('change');
        }
        }
    }
    function select_individu_Collapse() {
        $('#private_list_individu').collapse('hide');
        $('#select_individu_Card').collapse('show');
        val_select_private_list_individu_collapse = false;
    }
    function select_private_list_individu_collapse() {
        $('#select_individu_Card').collapse('hide');
        $('#private_list_individu').collapse('show');
        val_select_private_list_individu_collapse = true;
    }
    function select_gene_Collapse() {
        $('#private_list_gene').collapse('hide');
        $('#select_gene_Card').collapse('show');
        val_select_private_list_gene_collapse = false;
    }
    function select_private_list_gene_collapse() {
        $('#select_gene_Card').collapse('hide');
        $('#private_list_gene').collapse('show');
        val_select_private_list_gene_collapse = true;
    }


    $('#database_1').click({param: "1"}, database)
    $('#database_2').click({param: "2"}, database)
    $('#database_3').click({param: "3"}, database)
    $('#database_4').click({param: "4"}, database)
    $('#database_5').click({param: "5"}, database)
    $('#database_6').click({param: "6"}, database)
    $('#database_7').click({param: "7"}, database)
    $('#database_8').click({param: "8"}, database)
    $('#select_individu_Collapse').click(select_individu_Collapse);
    $('#select_private_list_individu_collapse').click(select_private_list_individu_collapse);
    $('#select_gene_Collapse').click(select_gene_Collapse);
    $('#select_private_list_gene_collapse').click(select_private_list_gene_collapse);

    //Login to Highlander
    $('#setLoginHighlander').on('click', function () {
        var formData_LoginHighlander = fillData_loginHighlander();
        LoginHighlanderFunction(formData_LoginHighlander);
    });

    //Next step 1 : samples selection
    $('#select_dataBase').on('click', function () {
        $('#Loading_data_valid').addClass("is-valid");
        $('#Loading_data_valid').show();
        $('#card_drop_Database_selection').addClass("border-success");
        $(this).prop("disabled", true);
        var formData_sample_selection = fillData_sample_selection();
        sampleSelectionFunction(formData_sample_selection);
    });

    // CHOSE SAMPLE BASED ON ALL SAMPLE IN HIGHLANDER
    //search bar for all sample selection
    $("#search_sample_all").on("keyup", function () {
        var filter = $(this).val().toLowerCase();
        var tmp_full_list = document.getElementById("sample_list_full");
        var li_full = tmp_full_list.getElementsByTagName("option");
        var tmp_selected_patient_list = document.getElementById("patient_list_selected");
        var li_patient = tmp_selected_patient_list.getElementsByTagName("option");
        var tmp_selected_control_list = document.getElementById("control_list_selected");
        var li_control = tmp_selected_control_list.getElementsByTagName("option");
        for (i = 0; i < li_full.length; i++) {
        txtValue = tmp_full_list.options[i].text;
        if(li_patient[i].style.display == "none" && li_control[i].style.display == "none"){
            if (txtValue.toLowerCase().indexOf(filter) > -1) {
            li_full[i].style.display = "";
            } else {
            li_full[i].style.display = "none";
            }
        }
        }
    });
    $('#remove_from_search_sample_all').on('click', function (){
        var new_value = '';
        $("#search_sample_all").val(new_value).trigger('change');
        var tmp_full_list = document.getElementById("sample_list_full");
        var li_full = tmp_full_list.getElementsByTagName("option");
        var tmp_selected_patient_list = document.getElementById("patient_list_selected");
        var li_patient = tmp_selected_patient_list.getElementsByTagName("option");
        var tmp_selected_control_list = document.getElementById("control_list_selected");
        var li_control = tmp_selected_control_list.getElementsByTagName("option");
        for (i = 0; i < li_full.length; i++) {
        if(li_patient[i].style.display == "none" && li_control[i].style.display == "none"){
            li_full[i].style.display = "";
        }  
        }
    });
    $('#add_patient_from_select_sample_all').on('click', function (){
        var tmp_full_list = document.getElementById("sample_list_full");
        var tmp_selected_patient_list = document.getElementById("patient_list_selected");
        var li_full = tmp_full_list.getElementsByTagName("option");
        var li = tmp_selected_patient_list.getElementsByTagName("option");
        for (i = 0; i < tmp_full_list.options.length; i++) {
        if(tmp_full_list.options[i].selected){
            li_full[i].style.display = "none";
            li[i].style.display = "";
            tmp_full_list.options[i].selected = false;
        }
        }   
    });
    $('#remove_patient_from_select_sample_all').on('click', function (){
        var tmp_full_list = document.getElementById("sample_list_full");
        var tmp_selected_patient_list = document.getElementById("patient_list_selected");
        var li_full = tmp_full_list.getElementsByTagName("option");
        var li = tmp_selected_patient_list.getElementsByTagName("option");
        for (i = 0; i < tmp_selected_patient_list.options.length; i++) {
        if(tmp_selected_patient_list.options[i].selected){
            li_full[i].style.display = "";
            li[i].style.display = "none";
        }
        }
    });
    $('#add_control_from_select_sample_all').on('click', function (){
        var tmp_full_list = document.getElementById("sample_list_full");
        var tmp_selected_patient_list = document.getElementById("control_list_selected");
        var li_full = tmp_full_list.getElementsByTagName("option");
        var li = tmp_selected_patient_list.getElementsByTagName("option");
        for (i = 0; i < tmp_full_list.options.length; i++) {
        if(tmp_full_list.options[i].selected){
            li_full[i].style.display = "none";
            li[i].style.display = "";
            tmp_full_list.options[i].selected = false;
        }
        }
    });
    $('#remove_control_from_select_sample_all').on('click', function (){
        var tmp_full_list = document.getElementById("sample_list_full");
        var tmp_selected_patient_list = document.getElementById("control_list_selected");
        var li_full = tmp_full_list.getElementsByTagName("option");
        var li = tmp_selected_patient_list.getElementsByTagName("option");
        for (i = 0; i < tmp_selected_patient_list.options.length; i++) {
        if(tmp_selected_patient_list.options[i].selected){
            li_full[i].style.display = "";
            li[i].style.display = "none";
        }
        } 
    });

    // CHOSE SAMPLE BASED ON USER'S SAMPLE LIST IN HIGHLANDER
    //search bar 
    $("#search_list_sample_all").on("keyup", function () {
        var filter = $(this).val().toLowerCase();
        var tmp_full_list = document.getElementById("list_sample_full");
        var li_full = tmp_full_list.getElementsByTagName("option");
        for (i = 0; i < li_full.length; i++) {
        txtValue = tmp_full_list.options[i].text;
        if(selected_list_samples_as_patient[i].value == 0 && selected_list_samples_as_control[i].value == 0){
            if (txtValue.toLowerCase().indexOf(filter) > -1) {
            li_full[i].style.display = "";
            } else {
            li_full[i].style.display = "none";
            }
        }
        }
    });
    //clean search bar
    $('#remove_from_search_list_sample_all').on('click', function (){
        var new_value = '';
        $("#search_list_sample_all").val(new_value).trigger('change');
        var tmp_full_list = document.getElementById("list_sample_full");
        var li_full = tmp_full_list.getElementsByTagName("option");
        for (i = 0; i < li_full.length; i++) {
        if(selected_list_samples_as_patient[i].value == 0 && selected_list_samples_as_control[i].value == 0){
            li_full[i].style.display = "";
        }
        }
    });
    //display sample name selected for patient
    $('#add_from_select_list_patient').on('click', function (){
        var tmp_full_list = document.getElementById("list_sample_full");
        var tmp_selected_patient_list = document.getElementById("list_patient_selected");
        var li_full = tmp_full_list.getElementsByTagName("option");
        var li = tmp_selected_patient_list.getElementsByTagName("option");
        for (i = 0; i < tmp_full_list.options.length; i++) {
        if(tmp_full_list.options[i].selected){
            li_full[i].style.display = "none";
            selected_list_samples_as_patient[i].value = 1;
            tmp_full_list.options[i].selected = false;
            tmp_sample_list_user = sample_list_user[i];
            // var tmp_to_string = JSON.stringify(tmp_sample_list_user.value);
            for(j = 0; j < Object.keys(tmp_sample_list_user.value).length; j++){
            index = tmp_sample_list_user.value[j].value;
            if(index > -1){
                li[index].style.display = "";
            }
            }
        }
        }
    });
    //remove all samples name selected for patient
    $('#remove_from_select_list_patient').on('click', function (){
        var tmp_full_list = document.getElementById("list_sample_full");
        var tmp_selected_patient_list = document.getElementById("list_patient_selected");
        var li_full = tmp_full_list.getElementsByTagName("option");
        var li = tmp_selected_patient_list.getElementsByTagName("option");
        for (i = 0; i < Object.keys(selected_list_samples_as_patient).length; i++) {
        if(selected_list_samples_as_patient[i].value == 1){
            selected_list_samples_as_patient[i].value = 0;
            li_full[i].style.display = "";
            tmp_sample_list_user = sample_list_user[i];
            // var tmp_to_string = JSON.stringify(tmp_sample_list_user.value);
            for(j = 0; j < Object.keys(tmp_sample_list_user.value).length; j++){
            index = tmp_sample_list_user.value[j].value;
            if(index > -1){
                li[index].style.display = "none";
            }
            }
        }
        } 
    });
    //display sample name selected for control
    $('#add_from_select_list_control').on('click', function (){
        var tmp_full_list = document.getElementById("list_sample_full");
        var tmp_selected_control_list = document.getElementById("list_control_selected");
        var li_full = tmp_full_list.getElementsByTagName("option");
        var li = tmp_selected_control_list.getElementsByTagName("option");
        for (i = 0; i < tmp_full_list.options.length; i++) {
        if(tmp_full_list.options[i].selected){
            li_full[i].style.display = "none";
            selected_list_samples_as_control[i].value = 1;
            tmp_full_list.options[i].selected = false;
            tmp_sample_list_user = sample_list_user[i];
            // var tmp_to_string = JSON.stringify(tmp_sample_list_user.value);
            for(j = 0; j < Object.keys(tmp_sample_list_user.value).length; j++){
            index = tmp_sample_list_user.value[j].value;
            if(index > -1){
                li[index].style.display = "";
            }
            }
        }
        }  
    });
    //remove all samples name selected for control
    $('#remove_from_select_list_control').on('click', function (){
        var tmp_full_list = document.getElementById("list_sample_full");
        var tmp_selected_control_list = document.getElementById("list_control_selected");
        var li_full = tmp_full_list.getElementsByTagName("option");
        var li = tmp_selected_control_list.getElementsByTagName("option");
        for (i = 0; i < Object.keys(selected_list_samples_as_control).length; i++) {
        if(selected_list_samples_as_control[i].value == 1){
            selected_list_samples_as_control[i].value = 0;
            li_full[i].style.display = "";
            tmp_sample_list_user = sample_list_user[i];
            // var tmp_to_string = JSON.stringify(tmp_sample_list_user.value);
            for(j = 0; j < Object.keys(tmp_sample_list_user.value).length; j++){
            index = tmp_sample_list_user.value[j].value;
            if(index > -1){
                li[index].style.display = "none";
            }
            }
        }
        }
    });

    //Next step after sample selection : gene list selection OR filtering
    $('#next_step_after_sample_selection').on('click', function () {
        var sample_selection_is_ok = true;
        //If user select sample based on private list on Highlander, check if sample do not appears in both patient and controls
        if(val_select_private_list_individu_collapse){
        var tmp_selected_patient_list = document.getElementById("list_patient_selected");
        var li_patient = tmp_selected_patient_list.getElementsByTagName("option");
        var tmp_selected_control_list = document.getElementById("list_control_selected");
        var li_control = tmp_selected_control_list.getElementsByTagName("option");
        var count_patient = 0;
        var count_control = 0;
        for(i = 0; i < tmp_selected_patient_list.options.length; i++){
            if(li_patient[i].style.display == ""){
            count_patient = count_patient + 1;
            }
            if(li_control[i].style.display == ""){
            count_control = count_control + 1;
            }
            if(li_patient[i].style.display == "" && li_control[i].style.display == ""){
            sample_selection_is_ok = false;
            $('#sample_selection_via_list_invalid').addClass("is-invalid");
            $('#sample_selection_via_list_invalid').show();
            $('#private_list_individu').addClass("border-danger")
            }
        }
        if (count_patient == 0) {
            $('#list_patient_selected_invalid').addClass("is-invalid");
            $('#list_patient_selected_invalid').show();
            $('#private_list_individu').addClass("border-danger");
        }else{
            $('#list_patient_selected_invalid').hide();
            if(sample_selection_is_ok){
            $('#sample_selection_via_list_invalid').hide();
            }
        }
        if (count_control == 0) {
            $('#list_control_selected_invalid').addClass("is-invalid");
            $('#list_control_selected_invalid').show();
            $('#private_list_individu').addClass("border-danger");
        }else{
            $('#list_control_selected_invalid').hide();
            if(sample_selection_is_ok){
            $('#sample_selection_via_list_invalid').hide();
            }
        }
        }else{
        var tmp_selected_patient_list = document.getElementById("patient_list_selected");
        var li_patient = tmp_selected_patient_list.getElementsByTagName("option");
        var tmp_selected_control_list = document.getElementById("control_list_selected");
        var li_control = tmp_selected_control_list.getElementsByTagName("option");
        var count_patient = 0;
        var count_control = 0;
        for(i = 0; i < tmp_selected_patient_list.options.length; i++){
            if(li_patient[i].style.display == ""){
            count_patient = count_patient + 1;
            }
            if(li_control[i].style.display == ""){
            count_control = count_control + 1;
            }
            if(li_patient[i].style.display == "" && li_control[i].style.display == ""){
            sample_selection_is_ok = false;
            }
        }
        if (count_patient == 0) {
            $('#patient_list_selected_invalid').addClass("is-invalid");
            $('#patient_list_selected_invalid').show();
            $('#select_individu_Card').addClass("border-danger");
        }else{
            $('#patient_list_selected_invalid').hide();
        }
        if (count_control == 0) {
            $('#control_list_selected_invalid').addClass("is-invalid");
            $('#control_list_selected_invalid').show();
            $('#select_individu_Card').addClass("border-danger");
        }else{
            $('#control_list_selected_invalid').hide();
        }
        }
        //lauch next
        if(sample_selection_is_ok == true && count_patient > 0 && count_control > 0){
        $('#private_list_individu').removeClass("border-danger");
        $('#select_individu_Card').removeClass("border-danger");
        $('#Loading_gene_valid').addClass("is-valid");
        $('#Loading_gene_valid').show();
        $('#card_drop_sample_selection').addClass("border-success");
        $(this).prop("disabled", true);
        if($('#gene_list').val() == 'on'){
            var formData_gene_selection = fillData_sample_selection();
            geneSelectionFunction(formData_gene_selection);
        }else{
            var formData_filtering = fillData_filtering();
            filteringFunction(formData_filtering);
            $("#Sample_selection").addClass("visually-hidden ");
            $("#Sample_selection").addClass("disabled");
            event.preventDefault()
            event.stopPropagation()
            $("#filtering_selection").removeClass("visually-hidden");
        }
        
        }
    });

    // CHOSE GENE BASED ON ALL GENE IN HIGHLANDER
    //search bar for all gene selection
    $("#search_gene_all").on("keyup", function () {
        var filter = $(this).val().toLowerCase();
        var tmp_full_list = document.getElementById("gene_list_full");
        var li_full = tmp_full_list.getElementsByTagName("option");
        var tmp_selected_gene_list = document.getElementById("gene_list_selected");
        var li_gene = tmp_selected_gene_list.getElementsByTagName("option");
        for (i = 0; i < li_full.length; i++) {
        txtValue = tmp_full_list.options[i].text;
        if(li_gene[i].style.display == "none"){
            if (txtValue.toLowerCase().indexOf(filter) > -1) {
            li_full[i].style.display = "";
            } else {
            li_full[i].style.display = "none";
            }
        }
        }
    });
    $('#remove_from_search_gene_all').on('click', function (){
        var new_value = '';
        $("#search_gene_all").val(new_value).trigger('change');
        var tmp_full_list = document.getElementById("gene_list_full");
        var li_full = tmp_full_list.getElementsByTagName("option");
        var tmp_selected_gene_list = document.getElementById("gene_list_selected");
        var li_gene = tmp_selected_gene_list.getElementsByTagName("option");
        for (i = 0; i < li_full.length; i++) {
        if(li_gene[i].style.display == "none"){
            li_full[i].style.display = "";
        }  
        }
    });
    $('#add_gene_from_select_gene_all').on('click', function (){
        var tmp_full_list = document.getElementById("gene_list_full");
        var tmp_selected_gene_list = document.getElementById("gene_list_selected");
        var li_full = tmp_full_list.getElementsByTagName("option");
        var li = tmp_selected_gene_list.getElementsByTagName("option");
        for (i = 0; i < tmp_full_list.options.length; i++) {
        if(tmp_full_list.options[i].selected){
            li_full[i].style.display = "none";
            li[i].style.display = "";
            tmp_full_list.options[i].selected = false;
        }
        }   
    });
    $('#remove_gene_from_select_gene_all').on('click', function (){
        var tmp_full_list = document.getElementById("gene_list_full");
        var tmp_selected_gene_list = document.getElementById("gene_list_selected");
        var li_full = tmp_full_list.getElementsByTagName("option");
        var li = tmp_selected_gene_list.getElementsByTagName("option");
        for (i = 0; i < tmp_selected_gene_list.options.length; i++) {
        if(tmp_selected_gene_list.options[i].selected){
            li_full[i].style.display = "";
            li[i].style.display = "none";
        }
        }
    });

    // CHOSE gene BASED ON USER'S gene LIST IN HIGHLANDER
    //search bar 
    $("#search_list_gene_all").on("keyup", function () {
        var filter = $(this).val().toLowerCase();
        var tmp_full_list = document.getElementById("list_gene_full");
        var li_full = tmp_full_list.getElementsByTagName("option");
        for (i = 0; i < li_full.length; i++) {
        txtValue = tmp_full_list.options[i].text;
        if(gene_list_user[i].selected == 0){
            if (txtValue.toLowerCase().indexOf(filter) > -1) {
            li_full[i].style.display = "";
            } else {
            li_full[i].style.display = "none";
            }
        }
        }
    });
    //clean search bar
    $('#remove_from_search_list_gene_all').on('click', function (){
        var new_value = '';
        $("#search_list_gene_all").val(new_value).trigger('change');
        var tmp_full_list = document.getElementById("list_gene_full");
        var li_full = tmp_full_list.getElementsByTagName("option");
        for (i = 0; i < li_full.length; i++) {
        if(gene_list_user[i].selected == 0){
            li_full[i].style.display = "";
        }
        }
    });
    //display gene name selected for gene
    $('#add_from_select_list_gene').on('click', function (){
        var tmp_full_list = document.getElementById("list_gene_full");
        var tmp_selected_gene_list = document.getElementById("list_gene_selected");
        var li_full = tmp_full_list.getElementsByTagName("option");
        var li = tmp_selected_gene_list.getElementsByTagName("option");
        for (i = 0; i < tmp_full_list.options.length; i++) {
        if(tmp_full_list.options[i].selected){
            li_full[i].style.display = "none";
            gene_list_user[i].selected = 1;
            tmp_full_list.options[i].selected = false;
            tmp_gene_list_user = gene_list_user[i];
            for(j = 0; j < Object.keys(tmp_gene_list_user.value).length; j++){
            index = tmp_gene_list_user.value[j].value;
            if(index > -1){
                li[index].style.display = "";
            }
            }
        }
        }  
    });
    //remove all genes name selected for gene
    $('#remove_from_select_list_gene').on('click', function (){
        var tmp_full_list = document.getElementById("list_gene_full");
        var tmp_selected_gene_list = document.getElementById("list_gene_selected");
        var li_full = tmp_full_list.getElementsByTagName("option");
        var li = tmp_selected_gene_list.getElementsByTagName("option");
        for (i = 0; i < tmp_full_list.options.length; i++) {
        if(gene_list_user[i].selected == 1){
            gene_list_user[i].selected = 0;
            li_full[i].style.display = "";
            tmp_gene_list_user = gene_list_user[i];
            
            for(j = 0; j < Object.keys(tmp_gene_list_user.value).length; j++){
            index = tmp_gene_list_user.value[j].value;
            if(index > -1){
                li[index].style.display = "none";
            }
            }
        }
        } 
    });

    //Next step after gene selection : filtering
    $('#next_step_after_gene_selection').on('click', function () {
        var gene_selection_is_ok = true;
        if(val_select_private_list_gene_collapse){
        var tmp_selected_gene_list = document.getElementById("list_gene_selected");
        var li_gene = tmp_selected_gene_list.getElementsByTagName("option");
        var count_gene = 0;
        for(i = 0; i < tmp_selected_gene_list.options.length; i++){
            if(li_gene[i].style.display == ""){
            count_gene = count_gene + 1;
            }
        }
        if (count_gene == 0) {
            $('#list_gene_selected_invalid').addClass("is-invalid");
            $('#list_gene_selected_invalid').show();
            $('#private_list_gene').addClass("border-danger");
        }else{
            $('#list_gene_selected_invalid').hide();
        }
        }else{
        var tmp_selected_gene_list = document.getElementById("gene_list_selected");
        var li_gene = tmp_selected_gene_list.getElementsByTagName("option");
        var count_gene = 0;
        for(i = 0; i < tmp_selected_gene_list.options.length; i++){
            if(li_gene[i].style.display == ""){
            count_gene = count_gene + 1;
            }
        }
        if (count_gene == 0) {
            $('#gene_list_selected_invalid').addClass("is-invalid");
            $('#gene_list_selected_invalid').show();
            $('#select_gene_Card').addClass("border-danger");
        }else{
            $('#gene_list_selected_invalid').hide();
        }
        }
        //lauch next
        if(count_gene > 0){
        $('#private_list_gene').removeClass("border-danger");
        $('#select_gene_Card').removeClass("border-danger");
        $('#selection_gene_valid').addClass("is-valid");
        $('#selection_gene_valid').show();
        $('#card_drop_gene_selection').addClass("border-success");
        $(this).prop("disabled", true);
        var formData_filtering = fillData_filtering();
        filteringFunction(formData_filtering);
        $("#gene_selection").addClass("visually-hidden ");
        $("#gene_selection").addClass("disabled");
        event.preventDefault()
        event.stopPropagation()
        $("#filtering_selection").removeClass("visually-hidden");
        }
    });


    //Load data via the website
    $('#search_file_area').on('click', function () {
        document.getElementById('formFile').onchange = function() {
        fileobj = document.getElementById('formFile').files;
        patients_data = fileobj;
        var formData = fillData(patients_data, 'patient');
        type_of_file='patient';
        uploadFormData(formData, 'patient');
    };
    });

    $('#search_file_area2').on('click', function () {
        document.getElementById('formFile2').onchange = function() {
        fileobj = document.getElementById('formFile2').files;
        controls_data = fileobj;
        var formData = fillData(controls_data, 'control');
        type_of_file='control';
        uploadFormData(formData, 'control');
    };
    });

    $('#search_file_area3').on('click', function () {
        document.getElementById('formFile3').onchange = function() {
        fileobj = document.getElementById('formFile3').files;
        gene_list = fileobj;
        var formData = fillData(gene_list, 'gene_list');
        type_of_file='gene_list';
        uploadFormData(formData, 'gene_list');
    };
    });


    // TO Develop further more
    // Filtering page when submitting data via Highlander
    annotation_for_filtering = ['variant_type']; //Put here all names of annotation to retrieve from Highlander
    count_variant_type = 1;

    function fillData_filtering() {
        var formData = new FormData();
        formData.append('User_name_highlander', $('#User_name_highlander').val());
        formData.append('Password_highlander', $('#Password_highlander').val());
        formData.append('analysis', selected_analysis_name);
        formData.append('annotation_for_filtering', annotation_for_filtering);
        return formData;
    }

    function filteringFunction(form_data) {
        $.ajax({
        url: "../php/filtering_highlander.php",
        method: "POST",
        data: form_data,
        contentType: false,
        cache: false,
        processData: false,
        xhr: function(){
            //Upload progress
            var xhr = new window.XMLHttpRequest();
            xhr.upload.addEventListener("progress", function(evt){
            if (evt.lengthComputable) {  
                var percentComplete = Math.round((evt.loaded / evt.total) * 100);
                // if needed you can add a porgress bar here 
                var msg = percentComplete + '% completed';
            }
            }, false); 
            return xhr;
        },
        success: function (data) {
            var resp = $.parseJSON(data);
            if (resp.error) {
            //Not working
            }else{
            //is ok
            $("#variant_type_full").html(resp.new_option);
            $("#variant_type_selected").html(resp.new_option_selected);
            }
            
        },
        error: function (data) {
            console.log("AJAX ERR sample selection Highlander");
            console.log(data);
        }
        });
    }

    function Filter_pass() {
        var new_val = 'PASS';
        $('#Filter_pass').val(new_val).trigger('change');
        var new_val2 = '';
        $('#Filter_pass_no').val(new_val2).trigger('change');
    }
    function Filter_pass_no() {
        var new_val = '';
        $('#Filter_pass').val(new_val).trigger('change');
        var new_val2 = 'on';
        $('#Filter_pass_no').val(new_val2).trigger('change');
    }
    $('#Filter_pass').click(Filter_pass)
    $('#Filter_pass_no').click(Filter_pass_no)

    $('#add_variant_type').on('click', function (){
        var tmp_full_list = document.getElementById("variant_type_full");
        var tmp_selected_list = document.getElementById("variant_type_selected");
        var li_full = tmp_full_list.getElementsByTagName("option");
        var li = tmp_selected_list.getElementsByTagName("option");
        for (i = 0; i < tmp_full_list.options.length; i++) {
        if(tmp_full_list.options[i].selected){
            li_full[i].style.display = "none";
            li[i].style.display = "";
            tmp_full_list.options[i].selected = false;
        }
        }   
    });
    $('#remove_variant_type').on('click', function (){
        var tmp_full_list = document.getElementById("variant_type_full");
        var tmp_selected_list = document.getElementById("variant_type_selected");
        var li_full = tmp_full_list.getElementsByTagName("option");
        var li = tmp_selected_list.getElementsByTagName("option");
        for (i = 0; i < tmp_selected_list.options.length; i++) {
        if(tmp_selected_list.options[i].selected){
            li_full[i].style.display = "";
            li[i].style.display = "none";
        }
        }
    });

    // Submit via Highlander
    $('#submit_via_highlander').on('click', function () {
        // Check if everything is ok before submitting
        // Allele frequency 
        if ($('#maf_filter').val() && !isNaN($('#maf_filter').val()) && $('#maf_filter').val() > 0 && $('#maf_filter').val() < 1) {
        $('#maf_filter_invalid').hide();
        $('#card_drop_filtering').removeClass("border-danger");
        }else{
        $('#maf_filter_invalid').addClass("is-invalid");
        $('#maf_filter_invalid').show();
        $('#card_drop_filtering').addClass("border-danger");
        event.preventDefault()
        event.stopPropagation()
        return;
        }
        // Consensus prediction 
        if ($('#consensus_prediction_filter').val() && Number.isInteger(Number($('#consensus_prediction_filter').val())) && !isNaN($('#consensus_prediction_filter').val()) && $('#consensus_prediction_filter').val() >= 0 && $('#consensus_prediction_filter').val() < 501) {
        $('#consensus_prediction_filter_invalid').hide();
        $('#card_drop_filtering').removeClass("border-danger");
        }else{
        $('#consensus_prediction_filter_invalid').addClass("is-invalid");
        $('#consensus_prediction_filter_invalid').show();
        $('#card_drop_filtering').addClass("border-danger");
        event.preventDefault()
        event.stopPropagation()
        return;
        }
        // Variant Type 
        var tmp_selected_list = document.getElementById("variant_type_selected");
        var li = tmp_selected_list.getElementsByTagName("option");
        var count_variant_type = 0;
        for(i = 0; i < tmp_selected_list.options.length; i++){
        if(li[i].style.display == ""){
            count_variant_type = count_variant_type + 1;
        }
        }
        if (count_variant_type == 0) {
        $('#variant_type_selected_invalid').addClass("is-invalid");
        $('#variant_type_selected_invalid').show();
        $('#card_drop_filtering').addClass("border-danger");
        event.preventDefault()
        event.stopPropagation()
        return;
        }else{
        $('#variant_type_selected_invalid').hide();
        $('#card_drop_filtering').removeClass("border-danger");
        }
        
        var datas = fillDataFinal_highlander();

        document.location.href = '../pages/Live.php';
        $.ajax({
        url: "../php/Submit.php",
        method: "POST",
        data: datas,
        contentType: false,
        cache: false,
        processData: false,
        success: function (data) {
            
        },
        error: function (data) {
            console.log("err in Submit.php");
            console.log(data);
        }
        });
        event.preventDefault()
        event.stopPropagation()
    });


    //  PART FOR SUBMISSION USING WEBSITE
    var testedFiles = [];
    var patients_data; 
    var controls_data;
    var gene_list;
    var type_of_file='';

    //Data pretreatment and upload via the website
    function fillData(files, type_data) {
        var formData = new FormData();
        if (files.length === 0)
        return false;

        if(type_data === 'patient'){
        $('#noFile').hide();
        }
        if(type_data === 'control'){
        $('#noFile2').hide();
        }
        if(type_data === 'gene_list'){
        $('#noFile3').hide();
        }
        formData.append('email', $('#email').val());
        formData.append('pi', $('#pi option:selected').val());
        formData.append('Highlander_import', $('#Highlander_import').val());
        formData.append('do_coverage', $('#do_coverage').val());
        formData.append('do_gene_list', $('#gene_list').val());
        // formData.append('filtering', $('#filtering').val());
        formData.append('mapqv', $('#mapqv').val());
        formData.append('min_cover', $('#min_cover').val());
        formData.append('remove_individu', $('#remove_individu').val());
        formData.append('remove_position', $('#remove_position').val());
        formData.append('min_nbr_var', $('#min_nbr_var').val());
        formData.append('max_nbr_var', $('#max_nbr_var').val());
        // formData.append('rare_maf_threshold', $('#rare_maf_threshold').val());
        formData.append('file', files[0]);
        testedFiles.push(files[0]);

        return formData;
    }

    function fillDataFinal(file1, file2, file3) {
        var formData = new FormData();
        formData.append('email', $('#email').val());
        formData.append('pi', $('#pi option:selected').val());
        formData.append('Highlander_import', $('#Highlander_import').val());
        formData.append('do_coverage', $('#do_coverage').val());
        formData.append('do_gene_list', $('#gene_list').val());
        // formData.append('filtering', $('#filtering').val());
        formData.append('mapqv', $('#mapqv').val());
        formData.append('min_cover', $('#min_cover').val());
        formData.append('remove_individu', $('#remove_individu').val());
        formData.append('remove_position', $('#remove_position').val());
        formData.append('min_nbr_var', $('#min_nbr_var').val());
        formData.append('max_nbr_var', $('#max_nbr_var').val());
        // formData.append('rare_maf_threshold', $('#rare_maf_threshold').val());
        if (file1.length === 0) {
        return false;
        }else{
        formData.append('file[]', file1[0]);
        }
        if (file2.length === 0)
        return false;
        else{
        formData.append('file[]', file2[0]);
        }
        if(file3 === undefined){

        }else{
        formData.append('file[]', file3[0]);
        }
        return formData;
    }

    function uploadFormData(form_data, type_data) {
        $.ajax({
        url: "../php/Upload.php",
        method: "POST",
        data: form_data,
        contentType: false,
        cache: false,
        processData: false,
        xhr: function(){
            //Upload progress
            var xhr = new window.XMLHttpRequest();
            xhr.upload.addEventListener("progress", function(evt){
            if (evt.lengthComputable) {  
                var percentComplete = Math.round((evt.loaded / evt.total) * 100);
                if(type_data === 'patient'){
                $('#progress_bar_process').width(percentComplete + '%');
                $('#progress_bar_process').html(percentComplete + '% completed');
                }
                if(type_data === 'control'){
                $('#progress_bar_process2').width(percentComplete + '%');
                $('#progress_bar_process2').html(percentComplete + '% completed');
                }
                if(type_data === 'gene_list'){
                $('#progress_bar_process3').width(percentComplete + '%');
                $('#progress_bar_process3').html(percentComplete + '% completed');
                }
                var msg = percentComplete + '% completed';
            }
            }, false); 
            return xhr;
        },
        success: function (data) {
            if(data === ''){
            var resp = {error: "TRUE", msg: "No file selected or File too big: Please contact raphael.helaers@uclouvain.be,"};
            }else{
            var resp = $.parseJSON(data);
            }
            
            if (resp.error) {
            var _f = resp.msg.split(",");
            var img = '<img src="../img/failed.png" alt="twbs" width="32" height="32" class="rounded-circle flex-shrink-0">';
            var msg = '<small class="text-danger opacity-80 text-nowrap">' + _f[0] + '  ' + img + '</small>';
            var msg = _f[0];
            var err = '<div class="alert alert-danger d-flex align-items-center" role="alert"><svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" fill="currentColor" class="bi bi-exclamation-triangle-fill flex-shrink-0 me-2" viewBox="0 0 16 16" role="img" aria-label="Warning:"><path d="M8.982 1.566a1.13 1.13 0 0 0-1.96 0L.165 13.233c-.457.778.091 1.767.98 1.767h13.713c.889 0 1.438-.99.98-1.767L8.982 1.566zM8 5c.535 0 .954.462.9.995l-.35 3.507a.552.552 0 0 1-1.1 0L7.1 5.995A.905.905 0 0 1 8 5zm.002 6a1 1 0 1 1 0 2 1 1 0 0 1 0-2z"/></svg><div>' + msg + '</div></div>'
            
            $('#uploaded_error').text("");
            $('#uploaded_error').append(err);
            if(type_of_file === 'patient'){
                $('#card_drop').addClass("border-danger");
                $('#uploaded_selection_error').text("");
                $('#uploaded_selection_error').append(err);
                $('#search_file_area').addClass("is-invalid");
            }
            if(type_of_file === 'control'){
                $('#card_drop2').addClass("border-danger");
                $('#uploaded_selection_error2').text("");
                $('#uploaded_selection_error2').append(err);
                $('#search_file_area2').addClass("is-invalid");
            }   
            if(type_of_file === 'gene_list'){
                $('#card_drop3').addClass("border-danger");
                $('#uploaded_selection_error3').text("");
                $('#uploaded_selection_error3').append(err);
                $('#search_file_area3').addClass("is-invalid");
            }           
            testedFiles = [];
            } else {
            var _f = testedFiles[0].name;
            if(type_of_file === 'patient'){
                testedFiles = [];
                var new_val = 1;
                $('#uploaded_file').val(new_val).trigger('change');
                $('#uploaded_error').empty();
                $('#uploaded_selection_error').empty();
                $('#card_drop').removeClass("border-danger");
                $('#search_file_area').removeClass("is-invalid");
            }else if(type_of_file === 'control'){
                testedFiles = [];
                var new_val = 1;
                $('#uploaded_file2').val(new_val).trigger('change');
                $('#uploaded_error').empty();
                $('#uploaded_selection_error2').empty();
                $('#card_drop2').removeClass("border-danger");
                $('#search_file_area2').removeClass("is-invalid");
            }else if(type_of_file === 'gene_list'){
                testedFiles = [];
                var new_val = 1;
                $('#uploaded_file3').val(new_val).trigger('change');
                $('#uploaded_error').empty();
                $('#uploaded_selection_error3').empty();
                $('#card_drop3').removeClass("border-danger");
                $('#search_file_area3').removeClass("is-invalid");
            }
            }
        },
        error: function (data) {
            console.log("AJAX ERR in Upload.php");
            console.log(data);
        }
        });
    }


    $('#PatientFile').collapse('show')
    $('#ControlFile').collapse('show')
    var new_val = 0;
    $('#uploaded_file').val(new_val).trigger('change');
    $('#uploaded_file2').val(new_val).trigger('change');
    $('#uploaded_file3').val(new_val).trigger('change');


    // Submit via Website
    $('#submit_via_website').on('click', function () {
        //If data were uploaded via the website
        if ($('#uploaded_file').val() === '0') {
        $('#card_drop').addClass("border-danger");
        $('#search_file_area').addClass("is-invalid");
        $('#noFile').show();
        event.preventDefault()
        event.stopPropagation()
        return;
        } 
        if ($('#uploaded_file2').val() === '0') {
        $('#card_drop2').addClass("border-danger");
        $('#search_file_area2').addClass("is-invalid");
        $('#noFile2').show();
        event.preventDefault()
        event.stopPropagation()
        return;
        }
        if($('#gene_list').val() === 'on'){
        if ($('#uploaded_file3').val() === '0') {
            $('#card_drop3').addClass("border-danger");
            $('#search_file_area3').addClass("is-invalid");
            $('#noFile3').show();
            event.preventDefault()
            event.stopPropagation()
            return;
        }
        var datas = fillDataFinal(patients_data, controls_data, gene_list);
        }else{
        var datas = fillDataFinal(patients_data, controls_data);
        }

        $.ajax({
        url: "../php/Submit.php",
        method: "POST",
        data: datas,
        contentType: false,
        cache: false,
        processData: false,
        success: function (data) {
            document.location.href = '../pages/Live.php';
        },
        error: function (data) {
            console.log("err in Submit.php");
            console.log(data);
        }
        });
        event.preventDefault()
        event.stopPropagation()
    });


});