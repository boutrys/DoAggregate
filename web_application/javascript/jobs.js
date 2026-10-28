$(document).ready(function () {

    $("#search").on("keyup", function () {
        var value = $(this).val().toLowerCase();
        $("#jobsList a").filter(function () {
            $(this).toggle($(this).text().toLowerCase().indexOf(value) > -1)
        });
    });

    var refresh = function (row) {
        if(!row) return false;
        var data = row.split(/\s/).filter(n => n);
        var name = data[0];
        //var user = data[1];
        var state = data[2];
        var time = "Running since " +data[3]+"sec";
        var split = name.split(/_/);
        var prediction = "JOB name : " + name;
        var user=split[8]+" "+split[9];
        var lastid = $("#jobsList a:last").attr("id");
        if (!lastid) lastid = "job_0";
        var split_id = lastid.split("_");
        var nextindex = Number(split_id[1]) + 1;
        var a = document.createElement("a");
        var div = document.createElement("div");
        var p = document.createElement("p");
        var small = document.createElement("small");
        var smallUser = document.createElement("small");
        var span = document.createElement("span");


        a.id = "job_" + nextindex;
        a.classList.add("list-group-item");
        a.classList.add("text-wrap");
        div.classList.add("text-wrap");
        p.classList.add("text-wrap");

        div.classList.add("d-flex");
        div.classList.add("w-100");
        div.classList.add("justify-content-between");

        p.classList.add("mb-1");
        p.textContent = prediction;

        if (state === "PD"){
            small.classList.add("text-muted");
            small.textContent = "Pending";

        }
        else {
            small.textContent = time;
            small.classList.add("text-primary");
//            a.classList.add("list-group-item-primary");
        }

        //smallUser.classList.add("text-muted");
        smallUser.classList.add("fw-bold");
        smallUser.classList.add("text-primary");
        smallUser.textContent = user;

        span.classList.add("position-absolute");
        span.classList.add("top-50");
        span.classList.add("start-0");
        span.classList.add("translate-middle");
        span.classList.add("rounded-pill");
        span.classList.add("badge");
        span.classList.add("bg-secondary");
        span.textContent = Number(nextindex)-1;
        
        if(state === "PD") div.append(span);

        div.append(p);
        div.append(small);

        a.append(div);
        a.append(smallUser);

        if (state === "PD")
            $('#jobsList').append(a);
        else
            $('#jobsList').prepend(a);
        
        if(lastid==="0"){
            $('#'+"job_" + nextindex)
            .hide()
            .css('opacity',0.0)
            .slideDown('slow')
            .animate({opacity: 1.0})
        }
    };
    
    var getJobs = function () {
        $.ajax({
            type: 'POST',
            url: '../php/Jobs.php',
            dataType: 'json',
            cache: false,
            success: function (result) {
                result.forEach(element => {
                        refresh(element);
                });
            },
        });
    };
    getJobs();

});