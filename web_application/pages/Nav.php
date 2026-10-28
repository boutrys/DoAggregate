<!DOCTYPE html>
<html lang="en">

<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="author" content="Bill Brancart">
  <meta name="author" content="Simon Boutry">
  <title>excalibur</title>
  <link href="/excalibur/css/spacelab.bootstrap.min.css" rel="stylesheet" >
  <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js" integrity="sha384-ka7Sk0Gln4gmtz2MlQnikT1wXgYsOg+OMhuP+IlRH9sENBO0LRn5q+8nbTov4+1p" crossorigin="anonymous"></script>
  <script src="https://unpkg.com/@popperjs/core@2"></script>
  <script src="https://unpkg.com/tippy.js@6"></script>
  <script src="https://cdn.jsdelivr.net/npm/bs-stepper/dist/js/bs-stepper.min.js"></script>
  
  <link rel="stylesheet" href="https://unpkg.com/tippy.js@6/themes/light.css"/>
  <link rel="stylesheet" href="https://unpkg.com/tippy.js@6/animations/scale.css"/>
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bs-stepper/dist/css/bs-stepper.min.css">
  <link rel="stylesheet" href="/excalibur/css/index.css">

</head>

<body>
  <div class="container navbar-expand-lg ">
    <nav class=" navbar navbar-dark bg-dark fixed-top">
      <div class="container-fluid  ">
        <ul class="nav nav-pills">
          <a class="navbar-brand" href="/excalibur/index.php">
            <img src="/excalibur/img/excalibur_white.png" alt="" width="45" height="45" class="d-inline-block align-text-center">
            <b>Excalibur</b>
          </a>
          <button class="navbar-toggler" type="button" data-bs-toggle="collapse" data-bs-target="#navbarSupportedContent" aria-controls="navbarSupportedContent" aria-expanded="false" aria-label="Toggle navigation">
            <span class="navbar-toggler-icon"></span>
          </button>

          <div class="collapse navbar-collapse" id="navbarSupportedContent">
            <ul class="navbar-nav me-auto mb-1 mb-lg-0">
              <li class="nav-item ">
                <a class="nav-link" aria-current="page" href="/excalibur/pages/Form.php">Submit</a>
              </li>
              <li class="nav-item">
                <a class="nav-link" href="/excalibur/pages/Live.php">Job queue</a>
              </li>
            </ul>
          </div>
        </ul>
    </nav>
  </div>
  <nav class="navbar navbar-dark bg-dark ">
    <div class="container-fluid">
      <span class="navbar-brand mb-0 h1"><b>excalibur</b> </span>
    </div>
  </nav>
  </br>
  <script src="/excalibur/javascript/cookie.js"></script>
  <script src="/excalibur/javascript/load.js"></script>
  <script src="/excalibur/javascript/index.js"></script>
  <script src="/excalibur/javascript/form.js"></script>
  <script src="/excalibur/javascript/tooltips.js"></script>

</body>

</html>