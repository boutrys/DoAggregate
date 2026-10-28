$(document).ready(function () {

  function how_to_use_button() {
    $('#Output_description').collapse('hide');
    $('#FAQ').collapse('hide');
    $('#troubleshooting').collapse('hide');
    $('#how_to_use').removeClass('visually-hidden');
    $('#how_to_use').collapse('show');
  }
  function Output_description_button() {
    $('#how_to_use').collapse('hide');
    $('#FAQ').collapse('hide');
    $('#troubleshooting').collapse('hide');
    $('#Output_description').removeClass('visually-hidden');
    $('#Output_description').collapse('show');
  }
  function FAQ_button() {
    $('#how_to_use').collapse('hide');
    $('#Output_description').collapse('hide');
    $('#troubleshooting').collapse('hide');
    $('#FAQ').removeClass('visually-hidden');
    $('#FAQ').collapse('show');
  }
  function troubleshooting_button() {
    $('#how_to_use').collapse('hide');
    $('#FAQ').collapse('hide');
    $('#Output_description').collapse('hide');
    $('#troubleshooting').removeClass('visually-hidden');
    $('#troubleshooting').collapse('show');
  }

  $('#how_to_use_button').click(how_to_use_button);
  $('#Output_description_button').click(Output_description_button);
  $('#FAQ_button').click(FAQ_button);
  $('#troubleshooting_button').click(troubleshooting_button);

  function importing_via_Highlander_button() {
    $('#importing_your_data').collapse('hide');
    $('#importing_via_Highlander').removeClass('visually-hidden');
    $('#importing_via_Highlander').collapse('show');
  }
  function importing_your_data_button() {
    $('#importing_via_Highlander').collapse('hide');
    $('#importing_your_data').removeClass('visually-hidden');
    $('#importing_your_data').collapse('show');
  }

  $('#importing_via_Highlander_button').click(importing_via_Highlander_button);
  $('#importing_your_data_button').click(importing_your_data_button);

});