# Web application

PHP/JavaScript front end that lets users submit a cohort to the parallel implementation of DoAggregate (see
`../cluster`), follow their jobs in the Slurm queue and download their results. It was deployed at UCLouvain
alongside the Highlander variant database, from which users can select samples, filters and gene lists.

## Structure

- `index.php`: main page with the submission form, live job view and documentation.
- `pages/`: form, navigation and live job pages.
- `php/`: submission (`Submit.php`), upload, job listing (`Jobs.php`), result download and Highlander queries.
- `javascript/`, `css/`, `img/`: front-end code, Bootstrap themes and images.

## Deployment

1. Serve the folder with a PHP-enabled web server (PHP with `mysqli`) on a machine that can submit Slurm jobs.
2. Set the database connection in `php/configuration.php` with the environment variables `EXCALIBUR_DB_HOST`,
   `EXCALIBUR_DB_USER`, `EXCALIBUR_DB_PASS` and `EXCALIBUR_DB_NAME`. Never commit real credentials.
3. In `php/Submit.php`, set the path to `submit_Excalibur.sh` and the system user that launches the jobs.
4. Update the contact e-mail addresses in `index.php`, `pages/Form.php` and `javascript/form.js`.

The Highlander-specific parts (`*_highlander.php`) are only needed if you use a Highlander database. Users can
also upload TSV files directly.
