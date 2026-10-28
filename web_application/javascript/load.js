$(document).ready(function () {
    var link = document.querySelector("link[rel~='icon']");
    if (!link) {
        link = document.createElement('link');
        link.rel = 'icon';
        document.getElementsByTagName('head')[0].appendChild(link);
    }
    link.href = '/excalibur/img/excalibur.png';

    window.addEventListener("pageshow", function (event) {
        var historyTraversal =
            event.persisted || (typeof window.performance != "undefined" && window.performance.navigation.type === 2);
        if (historyTraversal)
            window.location.reload(true);
    });

    var email = getCookie("email");
    if (email)
        $("#email").val(email)

    var group = getCookie("group");
    if (group){
        jQuery("#pi option").filter(function(){
            return $.trim($(this).text()) ==  group
        }).prop('selected', true);
    }

  // if(location.search.substring(1)!==""){
  //   var parameters = location.search.substring(1);
  //   var temp = parameters.split("=");
  //   if(temp[1]==='1'){
  //     $("#err").text("An error has occurred. Please retry.");
  //     $("#alert").removeClass("visually-hidden ");
  // }
  // } else {
  //     $("#alert").addClass("visually-hidden ");
  // }
});