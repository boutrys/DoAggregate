<?php include('Nav.php'); ?>

<div class="container-fluid container-sm user-select-none fs-60 shadow-sm  w-75">
    <div class="form-floating form-floating-sm w-100" id="options">
        <div class="col-md mb-2">
            <div class="row">
                <?php echo "<br>";?>
            </div>
            <label class="label mb-1 fs-6 fw-bold text ">Starting your submission</label>
        </div>
        <div class="fw-light form-floating  mb-2">
            <input type="email" class="form-control input-sm" id="email" name="email"
                placeholder="name@uclouvain.be" autofocus="on" required>
            <label for="floatingInputGrid">Email</label>
            <div class="invalid-feedback mb-1">Please enter a valid email adress.</div>
        </div>
        <div class="form-floating mb-4">
            <select name="pi" class="fw-light form-select form-select-sm mb-2" id="pi"
                aria-label="Floating label select example" required>
                <option value="NON_SSS">NON-SSS - Select this option if you don't find your research group</option>
                <option value="BCHM_VANSCHAFTINGEN">DDUV - BCHM - Emile Van Schaftingen</option>
                <option value="BCHM_LALOUX">DDUV - BCHM - Géraldine Laloux</option>
                <option value="BCHM_BOMMER">DDUV - BCHM - Guido Bommer</option>
                <option value="BCHM_COLLET">DDUV - BCHM - Jean-François Collet</option>
                <option value="BCHM_DEWACHTER">DDUV - BCHM - Liselot Dewachter</option>
                <option value="BCHM_VEIGA">DDUV - BCHM - Maria Veiga da Cunha </option>
                <option value="CBIO_GATTO">DDUV - CBIO - Laurent Gatto</option>
                <option value="CELL_PIERREUX">DDUV - CELL - Christophe Pierreux</option>
                <option value="CELL_TYTECA">DDUV - CELL - Donatienne Tyteca</option>
                <option value="CELL_HENRIET">DDUV - CELL - Patrick Henriet</option>
                <option value="GEDI_LIMAYE">DDUV - GEDI - Nisha Limaye</option>
                <option value="GEHU_VIKKULA">DDUV - GEHU - Miikka Vikkula</option>
                <option value="GEPI_DECOTTIGNIES">DDUV - GEPI - Anabelle Decottignies</option>
                <option value="GEPI_DESMET">DDUV - GEPI - Charles Desmet</option>
                <option value="IMCA_COULIE">DDUV - IMCA - Pierre Coulie</option>
                <option value="IMCA_LUCAS">DDUV - IMCA - Sophie Lucas</option>
                <option value="INIM_STOCKIS">DDUV - INIM - Julie Stockis</option>
                <option value="LICR_VANDENEYNDE">DDUV - LICR - Benoît Van den Eynde</option>
                <option value="LICR_ZHU">DDUV - LICR - Jingjing Zhu</option>
                <option value="LICR_VIGNERON">DDUV - LICR - Nathalie Vigneron</option>
                <option value="LPAD_LEMAIGRE">DDUV - LPAD - Frédéric Lemaigre</option>
                <option value="LPAD_JACQUEMIN">DDUV - LPAD - Patrick Jacquemin</option>
                <option value="MEMI_VANGASTEL">DDUV - MEMI - Nick Van Gastel</option>
                <option value="MEXP_DEMOULIN">DDUV - MEXP - Jean-Baptiste Demoulin</option>
                <option value="MEXP_RENAUD">DDUV - MEXP - Jean-Christophe Renaud</option>
                <option value="MEXP_COUTELIER">DDUV - MEXP - Jean-Paul Coutelier</option>
                <option value="MEXP_DUMOUTIER">DDUV - MEXP - Laure Dumoutier</option>
                <option value="MEXP_HAVELANGE">DDUV - MEXP - Violaine Havelange</option>
                <option value="PHOS_RIDER">DDUV - PHOS - Mark Rider</option>
                <option value="SCBD_LIEN">DDUV - SCBD - Wen-Hui Lien</option>
                <option value="SIGN_CONSTANTINESCU">DDUV - SIGN - Stefan Constantinescu</option>
                <option value="SIGN_STAERK">DDUV - SIGN - Judith Staerk</option>
                <option value="TILS_VANDERBRUGGEN">DDUV - TILS - Pierre van der Bruggen </option>
                <option value="VIRO_MICHIELS">DDUV - VIRO - Thomas Michiels</option>
                <option value="CEMO_TISSIR">IONS - CEMO - Fadel Tissir</option>
                <option value="CEMO_DELDICQUE">IONS - CEMO - Louise Deldicque</option>
                <option value="CEMO_FRANCAUX">IONS - CEMO - Marc Francaux</option>
                <option value="CEMO_TAJEDDINE">IONS - CEMO - Nicolas Tajeddine</option>
                <option value="CEMO_KIENLEN">IONS - CEMO - Pascal Kienlen-Campard</option>
                <option value="CEMO_GAILLY">IONS - CEMO - Philippe Gailly</option>
                <option value="CEMO_MAGGI">IONS - CEMO - Pietro Maggi</option>
                <option value="CEMO_GUALDANI">IONS - CEMO - Roberta Gualdani</option>
                <option value="CEMO_SLECLERCQ">IONS - CEMO - Sophie Leclercq</option>
                <option value="CEMO_VANPESCH">IONS - CEMO - Vincent Van Pesch</option>
                <option value="CEMO_JOSSIN">IONS - CEMO - Yves Jossin</option>
                <option value="COSY_MOURAUX">IONS - COSY - André Mouraux</option>
                <option value="COSY_BERQUIN">IONS - COSY - Anne Berquin</option>
                <option value="COSY_DEWOLF">IONS - COSY - Arthur Dewolf</option>
                <option value="COSY_DELHAYE">IONS - COSY - Benoit Delhaye</option>
                <option value="COSY_SCHEPENS">IONS - COSY - Bénédicte Schepens</option>
                <option value="COSY_DEJAEGER">IONS - COSY - Dominique De Jaeger</option>
                <option value="COSY_CREVECOEUR">IONS - COSY - Frédéric Crevecoeur</option>
                <option value="COSY_LIBERATI">IONS - COSY - Giulia Liberati</option>
                <option value="COSY_BASTIEN">IONS - COSY - Guillaume Bastien</option>
                <option value="COSY_DUQUE">IONS - COSY - Julie Duqué</option>
                <option value="COSY_MISSAL">IONS - COSY - Marcus Missal</option>
                <option value="COSY_PENTA">IONS - COSY - Massimo Penta</option>
                <option value="COSY_COLLIGNON">IONS - COSY - Olivier Collignon</option>
                <option value="COSY_LAVANDHOMME">IONS - COSY - Patricia Lavand'homme</option>
                <option value="COSY_LEFÈVRE">IONS - COSY - Philippe Lefèvre</option>
                <option value="COSY_HARDWICK">IONS - COSY - Robert Hardwick</option>
                <option value="COSY_NOZARADAN">IONS - COSY - Sylvie Nozaradan</option>
                <option value="COSY_LEGRAIN">IONS - COSY - Valéry Legrain</option>
                <option value="COSY_BLEYENHEUFT">IONS - COSY - Yannick Bleyenheuft</option>
                <option value="NEUR_IVANOIU">IONS - NEUR - Adrian Ivanoiu</option>
                <option value="NEUR_JANKOVSKI">IONS - NEUR - Aleksandr Jankovski</option>
                <option value="NEUR_HANSEEUW">IONS - NEUR - Bernard Hanseeuw</option>
                <option value="NEUR_HUART">IONS - NEUR - Caroline Huart</option>
                <option value="NEUR_RAFTOPOULOS">IONS - NEUR - Christian Raftopoulos</option>
                <option value="NEUR_YUKSEL">IONS - NEUR - Demet Yuksel</option>
                <option value="NEUR_MASQUELIER">IONS - NEUR - Etienne Masquelier</option>
                <option value="NEUR_SERMEUS">IONS - NEUR - Luc Sermeus</option>
                <option value="NEUR_CILIO">IONS - NEUR - Maria-Roberta Cilio</option>
                <option value="NEUR_NASSOGNE">IONS - NEUR - Marie-Cécile Nassogne</option>
                <option value="NEUR_DECAT">IONS - NEUR - Monique Decat</option>
                <option value="NEUR_DEPOTTER">IONS - NEUR - Patrick De Potter</option>
                <option value="NEUR_DETIMARY">IONS - NEUR - Philippe de Timary</option>
                <option value="NEUR_ROMBAUX">IONS - NEUR - Philippe Rombaux</option>
                <option value="NEUR_ELTAHRY">IONS - NEUR - Riëm El Tahry</option>
                <option value="NEUR_FERRAOSANTOS">IONS - NEUR - Susana Ferrao Santos</option>
                <option value="NEUR_DELTOMBE">IONS - NEUR - Thierry Deltombe</option>
                <option value="NEUR_DUPREZ">IONS - NEUR - Thierry Duprez</option>
                <option value="NEUR_VANDERMEEREN">IONS - NEUR - Yves Vandermeeren</option>
                <option value="CARD_GERBER">IREC - CARD - Bernhard Gerber</option>
                <option value="CHEX_BONACCORSI">IREC - CHEX - Eliano Bonaccorsi</option>
                <option value="ECLI_HORSMANS">IREC - ECLI - Yves Horsmans</option>
                <option value="EDIN_BRICHARD">IREC - EDIN - Sonia Brichard</option>
                <option value="EDIN_GILON">IREC - EDIN - Patrick Gilon</option>
                <option value="EDIN_GROOTAERT">IREC - EDIN - Mandy Grootaert</option>
                <option value="EDIN_JONAS">IREC - EDIN - Jean-Christophe Jonas</option>
                <option value="EDIN_LAUSE">IREC - EDIN - Pascale Lause</option>
                <option value="EPID_ROBERT">IREC - EPID - Annie Robert</option>
                <option value="FATH_DESSY">IREC - FATH - Chantal Dessy</option>
                <option value="FATH_CORBET">IREC - FATH - Cyril Corbet</option>
                <option value="FATH_BALLIGAND">IREC - FATH - Jean-Luc Balligand</option>
                <option value="FATH_FERON">IREC - FATH - Olivier Feron</option>
                <option value="FATH_SONVEAUX">IREC - FATH - Pierre Sonveaux</option>
                <option value="GAEN_ILECLERCQ">IREC - GAEN - Isabelle Leclercq</option>
                <option value="GYNE_DOLMANS">IREC - GYNE - Marie-Madeleine Dolmans</option>
                <option value="IMAG_LECOUVET">IREC - IMAG - Frédéric Lecouvet</option>
                <option value="LTAP_DEWULF">IREC - LTAP - Joseph Dewulf</option>
                <option value="LTAP_HAUFROID">IREC - LTAP - Vincent Haufroid</option>
                <option value="LTAP_HUAUX">IREC - LTAP - François Huaux</option>
                <option value="MBLG_KABAMBAMUKADI">IREC - MBLG - Benoit Kabamba Mukadi</option>
                <option value="MIRO_MACHIELS">IREC - MIRO - Jean-Pascal Machiels</option>
                <option value="MIRO_LEE">IREC - MIRO - John Lee</option>
                <option value="MONT_OCAK">IREC - MONT - Sebahat Ocak</option>
                <option value="MORF_BEHETSWYDEMANS">IREC - MORF - Catherine Behets Wydemans</option>
                <option value="MORF_VANBOCKSTAL">IREC - MORF - Mieke Van Bockstal</option>
                <option value="NEFR_DEVUYST">IREC - NEFR - Olivier Devuyst</option>
                <option value="NMSK_STOQUART">IREC - NMSK - Gaëtan Stoquart</option>
                <option value="PEDI_SOKAL">IREC - PEDI - Etienne Sokal</option>
                <option value="PEDI_SCHEERS">IREC - PEDI - Isabelle Scheers</option>
                <option value="PEDI_LYSY">IREC - PEDI - Philippe Lysy</option>
                <option value="PNEU_PILETTE">IREC - PNEU - Charles Pilette</option>
                <option value="REPR_AMORIM">IREC - REPR - Christiani Amorim</option>
                <option value="RUMA_DUREZ">IREC - RUMA - Patrick Durez</option>
                <option value="ADDB_DESRIEUX">LDRI - ADDB - Anne des Rieux</option>
                <option value="ADDB_BELOQUI">LDRI - ADDB - Ana Beloqui-Garcia</option>
                <option value="ADDB_VANBEVER">LDRI - ADDB - Rita Vanbever</option>
                <option value="ADDB_LELOUP">LDRI - ADDB - Gaëtane Leloup</option>
                <option value="ADDB_LEPRINCE">LDRI - ADDB - Julian Leprince</option>
                <option value="BPBL_MUCCIOLI">LDRI - BPBL - Giulio Muccioli</option>
                <option value="BPBL_ALHOUAYEK">LDRI - BPBL - Mireille Alhouayek</option>
                <option value="CMFA_FREDERICK">LDRI - CMFA - Raphaël Frédérick</option>
                <option value="CMFA_LAMBERT">LDRI - CMFA - Didier Lambert</option>
                <option value="GNOS_JLECLERCQ">LDRI - GNOS - Joëlle Leclercq</option>
                <option value="MNUT_BINDELS">LDRI - MNUT - Laure Bindels</option>
                <option value="MNUT_DELZENNE">LDRI - MNUT - Nathalie Delzenne</option>
                <option value="MNUT_CANI">LDRI - MNUT - Patrice Cani</option>
                <option value="MNUT_EVERARD">LDRI - MNUT - Amandine Everard</option>
                <option value="REMA_GALLEZ">LDRI - REMA - Bernard Gallez</option>
                <option value="REMA_JORDAN">LDRI - REMA - Bénédicte Jordan</option>
                <option value="TFAR-CLIP_DALLEUR">LDRI - TFAR-CLIP - Olivia Dalleur</option>
                <option value="TFAR-CLIP_SPINEWINE">LDRI - TFAR-CLIP - Anne Spinewine</option>
                <option value="TFAR-CLIP_HENRARD">LDRI - TFAR-CLIP - Séverine Henrard</option>
                <option value="TFAR-FACM_LORENT">LDRI - TFAR-FACM - Joseph Lorent</option>
                <option value="TFAR-FACM_VANBAMBEKE">LDRI - TFAR-FACM - Françoise Van Bambeke</option>
                <option value="TFAR-PMGK_ELENS">LDRI - TFAR-PMGK - Laure Elens</option>
            </select>
            <label for="pi">Group</label>
            <div class="invalid-feedback">
                Please select a research group.
            </div>
        </div>

        <div class="fw-light form-floating  mb-4">
            <div class="form-check">
                <input class="form-check-input" type="radio" id="Highlander_import" name="Highlander_import" value="on" checked="" data-bs-toggle="collapse" href="#PageLoginHighlander" role="button" aria-expanded="false">
                <label class="form-check-label" for="Highlander_import">Import data via Highlander</label>
            </div>
            <div class="form-check">
                <input class="form-check-input" type="radio" id="manual_import" name="Highlander_import" value="off">
                <label class="form-check-label" for="manual_import">Import your data</label>
            </div>
        </div>

        <div class="fw-light form-floating  mb-4">
            <div class="form-check form-switch">
                <input class="form-check-input" type="checkbox" id="gene_list" name="gene_list" value="off" data-bs-toggle="collapse" href="#geneFile" role="button" aria-expanded="false">
                <label class="form-check-label" for="gene_list">Gene list</label>
            </div>
        </div>

        <div class="fw-light form-floating  mb-4">
            <div class="form-check form-switch" id=clean_data_switch_option>
                <input class="form-check-input" type="checkbox" id="do_coverage" name="do_coverage" checked="" value="on">
                <label class="form-check-label" for="do_coverage">Clean data</label>
            </div>
        </div>

        <div class="fw-light form-floating  mb-4">
            <a id="advancedOption" class="btn btn-primary"
            data-bs-toggle="collapse" href="#OptionCard" role="button" aria-expanded="false"
            aria-controls="multiCollapseExample1">Advanced parameters</a>
        </div>

        <div class="fw-light form-floating  mb-4">
            <div class="collapse multi-collapse" id=OptionCard>
                <label class="form-label fw-bold">Options for data cleaning using coverage matrix</label>

                <div class="row justify-content-md-begin ">
                    <div class="col-md-auto">
                        Mapqv threshold to determine if a position is well covered
                    </div>
                    <div class="col col-lg-2">
                        <input type="text" value="20" class="form-control " id="mapqv" name="mapqv">
                    </div>
                </div>

                <div class="row justify-content-md-begin ">
                    <div class="col-md-auto">
                        Minimum coverage to determine if a position is well covered
                    </div>
                    <div class="col col-lg-2">
                        <input type="text" value="20" class="form-control " id="min_cover" name="min_cover">
                    </div>
                </div>

                <div class="form-switch col-md-auto">
                    <input class="form-check-input" type="checkbox" id="remove_individu" name="remove_individu" checked="" value="on">
                    <label class="form-check-label" for="remove_individu">Remove individuals that are badly covered</label>
                </div>

                <div class="form-switch col-md-auto">
                    <input class="form-check-input" type="checkbox" id="remove_position" name="remove_position" checked="" value="on">
                    <label class="form-check-label" for="remove_position">Remove positions that are badly covered</label>
                </div>
                <div class="fw-light form-floating  mb-4">
                </div>
                <label class="form-label fw-bold">Options for aggregation tests</label>
                <div class="row justify-content-md-begin ">
                    <div class="col-md-auto">
                        Minimum number of unique variants
                    </div>
                    <div class="col col-lg-2">
                        <input type="text" value="2" class="form-control " id="min_nbr_var" name="min_nbr_var">
                    </div>
                </div>

                <div class="row justify-content-md-begin ">
                    <div class="col-md-auto">
                        Maximum number of unique variants
                    </div>
                    <div class="col col-lg-2">
                        <input type="text" value="100000" class="form-control " id="max_nbr_var" name="max_nbr_var">
                    </div>
                </div>

            </div>
        </div>
        
        <div class="d-grid gap-2 d-md-flex justify-content-md-end">
            <button id="setOptions" class="shadow-sm  btn btn-primary w-25">Next</button>
        </div>
    </div>
    
    <?php// This is the part for Highlander Login?>
    <div id="Login_Highlander" class="visually-hidden">
        <div class="row g-2 mb-2 ">
            <div class="col-sm">
            </div>
            <div class="col-sm-1 text-center">
            <span>
                <img src="../img/logo_highlander.png" alt="" width="300" height="187" class="">
            </span>
            </div>
            <div class="col-sm">
            </div>
        </div>
        
        <div class="col-md mb-2">
            <label class="label mb-1 fs-6 fw-bold text ">Please login to Highlander</label>
        </div>

        <div class="card card-body " id="card_drop_loginHighlander">
            <div class="fw-light form-floating  mb-2">
                <input type="text" class="form-control" id="User_name_highlander" name="User_name_highlander" autofocus="on" required>
                <label for="floatingInputGrid">Highlander User's Name</label>
            </div>
            <div class="fw-light form-floating  mb-2">
                <input type="password" class="form-control" id="Password_highlander" name="Password_highlander" autofocus="on" required>
                <label for="floatingInputGrid">Highlander Password</label>
            </div>
            <div id="uploaded_selection_error_loginHighlander" class="invalid-feedback">
                Wrong username and/or password : Please contact raphael.helaers@uclouvain.be if this problem persist.
            </div>
        </div>

        <div class="d-grid gap-2 d-md-flex justify-content-md-end">
                <button id="setLoginHighlander" class="shadow-sm  btn btn-primary w-25">Login</button>
        </div>
            
    </div>  


    <?php// This is the part for Database selection?>
    <div id="Database_selection" class="visually-hidden">
        <div class="row g-2 mb-2 ">
            <div class="col-sm">
            </div>
            <div class="col-sm-1 text-center">
            <span>
                <img src="../img/logo_highlander.png" alt="" width="300" height="187" class="">
            </span>
            </div>
            <div class="col-sm">
            </div>
        </div>
        
        <div class="col-md mb-2">
            <label class="label mb-1 fs-6 fw-bold text ">Please select the database for your analysis</label>
        </div>

        <div class="card card-body " id="card_drop_Database_selection">
            <div id="DB_1" class="form-check visually-hidden">
                <input class="form-check-input" type="radio" id="database_1" name="database_1" value="on" checked="" data-bs-toggle="collapse"  role="button" aria-expanded="false">
                <label class="form-check-label" id="database_label_1" for="database_1" value="Import"> </label>
            </div>
            <div id="DB_2" class="form-check visually-hidden">
                <input class="form-check-input" type="radio" id="database_2" name="database_1" value="off">
                <label class="form-check-label" id="database_label_2" for="database_2" value="Import"></label>
            </div>
            <div id="DB_3" class="form-check visually-hidden">
                <input class="form-check-input" type="radio" id="database_3" name="database_1" value="off">
                <label class="form-check-label" id="database_label_3" for="database_3" value="Import"></label>
            </div>
            <div id="DB_4" class="form-check visually-hidden">
                <input class="form-check-input" type="radio" id="database_4" name="database_1" value="off">
                <label class="form-check-label" id="database_label_4" for="database_4" value="Import"></label>
            </div>
            <div id="DB_5" class="form-check visually-hidden">
                <input class="form-check-input" type="radio" id="database_5" name="database_1" value="off">
                <label class="form-check-label" id="database_label_5" for="database_5" value="Import"></label>
            </div>
            <div id="DB_6" class="form-check visually-hidden">
                <input class="form-check-input" type="radio" id="database_6" name="database_1" value="off">
                <label class="form-check-label" id="database_label_6" for="database_6" value="Import"></label>
            </div>
            <div id="DB_7" class="form-check visually-hidden">
                <input class="form-check-input" type="radio" id="database_7" name="database_1" value="off">
                <label class="form-check-label" id="database_label_7" for="database_7" value="Import"></label>
            </div>
            <div id="DB_8" class="form-check visually-hidden">
                <input class="form-check-input" type="radio" id="database_8" name="database_1" value="off">
                <label class="form-check-label" id="database_label_8" for="database_8" value="Import"></label>
            </div>

            <div id="Loading_data_valid" class="valid-feedback">
                Please wait a few seconds, we are retrieving samples available in Highlander and your sample list(s) from your Highlander user's profile.
            </div>
        </div>

        <div class="card card-body " id="card_retrieve_sample">
            <div class="progress">
                <div class="progress-bar" id="progress_bar_sample_selection" role="progressbar" style="width: 0%;" aria-valuenow="0" aria-valuemin="0" aria-valuemax="100"0%>
                </div>
            </div>
        </div>

            
        <div class="d-grid gap-2 d-md-flex justify-content-md-end">
            <button id="select_dataBase" class="shadow-sm  btn btn-primary w-25">Next Step</button>
        </div>
    </div>  


    <?php// This is the part for Sample selection Using Highlander?>
    <div id="Sample_selection" class="visually-hidden">
        <div class="card card-body " id="card_drop_sample_selection">
            <div class="row g-2 mb-2 ">
                <div class="col-sm">
                </div>
                <div class="col-sm-1 text-center">
                <span>
                    <img src="../img/logo_highlander.png" alt="" width="300" height="187" class="">
                </span>
                </div>
                <div class="col-sm">
                </div>
            </div>
            <div class="row g-2 mb-2 ">
                <div class="col-sm">
                    <a id="select_individu_Collapse" class="btn btn-outline-primary btn-sm w-100"
                                data-bs-toggle="collapse" href="#select_individu_Card" role="button" aria-expanded="false"
                                aria-controls="multiCollapseExample1">Select samples from all available samples
                    </a>
                </div>
                <div class="col-sm-1 text-center">
                    OR
                </div>
                <div class="col-sm">
                    <button id="select_private_list_individu_collapse" class="btn btn-outline-primary btn-sm w-100" type="button"
                                data-bs-toggle="collapse" data-bs-target="#private_list_individu" aria-expanded="false"
                                aria-controls="multiCollapseExample2">Import your Highlander lists of sample
                    </button>
                </div>
            </div>
            
            <?php// Retrieve patient/controls from all samples available in highlander?>
            
            <div class="collapse multi-collapse card card-body" id="select_individu_Card">
                <div class="row g-2 mb-2 ">
                    <div class="col-sm">
                        <div class="row justify-content-md-begin align-items-center">
                        <label class="fw-bold text" for="sample_list_full">Select your Samples</label>
                            <div class="col-md d-flex align-items-center">
                                <input class="form-control mb-2" id="search_sample_all" type="text" placeholder="Search..">
                            </div>
                            <div class="col justify-content-md-begin">
                                <button id="remove_from_search_sample_all" class="shadow-sm  btn btn-primary">
                                    <span>
                                        <img src="../img/delete.png" alt="" width="15" height="15" class="">
                                    </span>
                                </button>
                            </div>
                        </div>
                        <div class="row justify-content-md-begin align-items-center">
                            <style>
                                #sample_list_full {
                                    height: 400px;
                                    overflow-y: scroll;
                                }
                            </style>
                            <select multiple="" class="list-group" id="sample_list_full">
                            </select>
                        </div>
                    </div>
                    <div class="col-sm-1 text-center">

                    </div>

                    <div class="col-sm">
                        <?php// Selected patient from all samples available in highlander?>
                        <div class="row justify-content-md-begin align-items-center">
                            <div class="col-sm-2">
                                <div class="row justify-content-end align-items-center">
                                    <button id="add_patient_from_select_sample_all" class="shadow-sm  btn btn-primary">
                                        <span>
                                            <img src="../img/blue-double-arrow-right.png" alt="" width="20" height="20" class="">
                                        </span>
                                    </button>
                                </div>
                                <div class="row">
                                    <?php echo "<br>";?>
                                </div> 
                                <div class="row justify-content-end align-items-center">
                                    <button id="remove_patient_from_select_sample_all" class="shadow-sm  btn btn-primary">
                                        <span>
                                            <img src="../img/blue-double-arrow-left.png" alt="" width="20" height="20" class="">
                                        </span>
                                    </button>
                                </div>
                            </div>
                            <div class="col-sm-1 text-center">

                            </div>
                            <div class="col-sm justify-content-md-begin " id="card_patient_list_selected">
                                <div class="row justify-content-md-center align-items-center">
                                    <label class="fw-bold text" for="patient_list_selected">Selected as Patients</label>
                                </div>
                                <div class="row justify-content-md-center align-items-center">
                                    <style>
                                        #patient_list_selected {
                                            height: 180px;
                                            overflow-y: scroll;
                                        }
                                    </style>
                                    <select multiple="" class="list-group" id="patient_list_selected">
                                    </select>
                                </div>
                                <div id="patient_list_selected_invalid" class="invalid-feedback">
                                    No patients selected! 
                                </div>
                            </div>
                        </div>
                        <div class="row">
                            <?php echo "<br>";?>
                            <?php echo "<br>";?>
                        </div> 
                        <?php// Selected controls from all samples available in highlander?>  
                        <div class="row justify-content-md-begin align-items-center">
                            <div class="col-sm-2">
                                <div class="row justify-content-end align-items-center">
                                    <button id="add_control_from_select_sample_all" class="shadow-sm  btn btn-primary">
                                        <span>
                                            <img src="../img/blue-double-arrow-right.png" alt="" width="20" height="20" class="">
                                        </span>
                                    </button>
                                </div>
                                <div class="row">
                                    <?php echo "<br>";?>
                                </div> 
                                <div class="row justify-content-end align-items-center">
                                    <button id="remove_control_from_select_sample_all" class="shadow-sm  btn btn-primary">
                                        <span>
                                            <img src="../img/blue-double-arrow-left.png" alt="" width="20" height="20" class="">
                                        </span>
                                    </button>
                                </div>
                            </div>
                            <div class="col-sm-1 text-center">

                            </div>
                            <div class="col-sm justify-content-md-begin " id="card_control_list_selected">
                                <div class="row justify-content-md-center align-items-center">
                                    <label class="fw-bold text" for="control_list_selected">Selected as Controls</label>
                                </div>
                                <div class="row justify-content-md-center align-items-center ">
                                    <style>
                                        #control_list_selected {
                                            height: 180px;
                                            overflow-y: scroll;
                                        }
                                    </style>
                                    <select multiple="" class="list-group" id="control_list_selected">
                                    </select>
                                </div>
                                <div id="control_list_selected_invalid" class="invalid-feedback">
                                    No controls selected! 
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>

            <?php// Retrieve patient/controls from highlander list?>
            
            <div class="collapse multi-collapse card card-body" id="private_list_individu">
                <div class="row g-2 mb-2 ">
                    <div class="col-sm">
                        <div class="row justify-content-md-begin align-items-center">
                        <label class="fw-bold text" for="list_sample_full">Select your list of samples</label>
                            <div class="col-md">
                                <input class="form-control mb-2" id="search_list_sample_all" type="text" placeholder="Search..">
                            </div>
                            <div class="col justify-content-end">
                                <button id="remove_from_search_list_sample_all" class="shadow-sm  btn btn-primary">
                                    <span>
                                        <img src="../img/delete.png" alt="" width="20" height="20" class="">
                                    </span>
                                </button>
                            </div>
                        </div>
                        <div class="row justify-content-md-begin align-items-center">
                            <style>
                                #list_sample_full {
                                    height: 400px;
                                    overflow-y: scroll;
                                }
                            </style>
                            <select multiple="" class="list-group" id="list_sample_full">
                            </select>
                        </div>
                    </div>
                    <div class="col-sm-1 text-center">

                    </div>

                    <div class="col-sm">
                        <?php// Selected patient from highlander list?>
                        <div class="row justify-content-md-begin align-items-center">
                            <div class="col-sm-2">
                                <div class="row justify-content-end align-items-center">
                                    <button id="add_from_select_list_patient" class="shadow-sm  btn btn-primary">
                                        <span>
                                            <img src="../img/blue-double-arrow-right.png" alt="" width="20" height="20" class="">
                                        </span>
                                    </button>
                                </div>
                                <div class="row">
                                    <?php echo "<br>";?>
                                </div> 
                                <div class="row justify-content-end align-items-center">
                                    <button id="remove_from_select_list_patient" class="shadow-sm  btn btn-primary">
                                        <span>
                                            <img src="../img/delete.png" alt="" width="20" height="20" class="">
                                        </span>
                                    </button>
                                </div>
                            </div>
                            <div class="col-sm-1 text-center">

                            </div>
                            <div class="col-sm justify-content-md-begin " id="card_list_patient_selected">
                                <div class="row justify-content-md-center align-items-center">
                                    <label class="fw-bold text" for="list_patient_selected">Selected as Patients</label>
                                </div>
                                <div class="row justify-content-md-center align-items-center">
                                    <style>
                                        #list_patient_selected {
                                            height: 180px;
                                            overflow-y: scroll;
                                        }
                                    </style>
                                    <select multiple="" class="list-group" id="list_patient_selected">
                                    </select>
                                </div>
                                <div id="list_patient_selected_invalid" class="invalid-feedback">
                                    No patients selected! 
                                </div>
                            </div>
                        </div>
                        <div class="row">
                            <?php echo "<br>";?>
                            <?php echo "<br>";?>
                        </div> 
                        <?php// Selected controls from highlander list?>  
                        <div class="row justify-content-md-begin align-items-center">
                            <div class="col-sm-2">
                                <div class="row justify-content-end align-items-center">
                                    <button id="add_from_select_list_control" class="shadow-sm  btn btn-primary">
                                        <span>
                                            <img src="../img/blue-double-arrow-right.png" alt="" width="20" height="20" class="">
                                        </span>
                                    </button>
                                </div>
                                <div class="row">
                                    <?php echo "<br>";?>
                                </div> 
                                <div class="row justify-content-end align-items-center">
                                    <button id="remove_from_select_list_control" class="shadow-sm  btn btn-primary">
                                        <span>
                                            <img src="../img/delete.png" alt="" width="20" height="20" class="">
                                        </span>
                                    </button>
                                </div>
                            </div>
                            <div class="col-sm-1 text-center">

                            </div>
                            <div class="col-sm justify-content-md-begin " id="card_list_control_selected">
                                <div class="row justify-content-md-center align-items-center">
                                    <label class="fw-bold text" for="list_control_selected">Selected as Controls</label>
                                </div>
                                <div class="row justify-content-md-center align-items-center">
                                    <style>
                                        #list_control_selected {
                                            height: 180px;
                                            overflow-y: scroll;
                                        }
                                    </style>
                                    <select multiple="" class="list-group" id="list_control_selected">
                                    </select>
                                </div>
                                <div id="list_control_selected_invalid" class="invalid-feedback">
                                    No controls selected! 
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
                <div id="sample_selection_via_list_invalid" class="invalid-feedback">
                    Sample(s) present both in patient and control, please check your list in Highlander. A sample cannot be present in both patient and control!
                </div>
            </div>

            <div id="Loading_gene_valid" class="valid-feedback">
                Please wait a few seconds, we are retrieving genes available in Highlander and your gene list(s) from your Highlander user's profile.
            </div>
        </div>  

        <div class="card card-body " id="card_retrieve_gene">
            <div class="progress">
                <div class="progress-bar" id="progress_bar_gene_selection" role="progressbar" style="width: 0%;" aria-valuenow="0" aria-valuemin="0" aria-valuemax="100"0%>
                </div>
            </div>
        </div>
        
        <div class="d-grid gap-2 d-md-flex justify-content-md-end">
            <button id="next_step_after_sample_selection" class="shadow-sm  btn btn-primary w-25">Next Step</button>
        </div>
    </div>


    <?php// This is the part for gene list selection Using Highlander?>
    <div id="gene_selection" class="visually-hidden">
        <div class="card card-body " id="card_drop_gene_selection">
            <div class="row g-2 mb-2 ">
                <div class="col-sm">
                </div>
                <div class="col-sm-1 text-center">
                <span>
                    <img src="../img/logo_highlander.png" alt="" width="300" height="187" class="">
                </span>
                </div>
                <div class="col-sm">
                </div>
            </div>
            <div class="row g-2 mb-2 ">
                <div class="col-sm">
                    <a id="select_gene_Collapse" class="btn btn-outline-primary btn-sm w-100"
                                data-bs-toggle="collapse" href="#select_gene_Card" role="button" aria-expanded="false"
                                aria-controls="multiCollapseExample1">Select genes from all available genes
                    </a>
                </div>
                <div class="col-sm-1 text-center">
                    OR
                </div>
                <div class="col-sm">
                    <button id="select_private_list_gene_collapse" class="btn btn-outline-primary btn-sm w-100" type="button"
                                data-bs-toggle="collapse" data-bs-target="#private_list_gene" aria-expanded="false"
                                aria-controls="multiCollapseExample2">Import your Highlander lists of genes
                    </button>
                </div>
            </div>
            
            <?php// Retrieve genes from all genes available in highlander?>
            
            <div class="collapse multi-collapse card card-body" id="select_gene_Card">
                <div class="row g-2 mb-2 ">
                    <div class="col-sm">
                        <div class="row justify-content-md-begin align-items-center">
                        <label class="fw-bold text" for="gene_list_full">Select your genes</label>
                            <div class="col-md">
                                <input class="form-control mb-2" id="search_gene_all" type="text" placeholder="Search..">
                            </div>
                            <div class="col justify-content-end">
                                <button id="remove_from_search_gene_all" class="shadow-sm  btn btn-primary">
                                    <span>
                                        <img src="../img/delete.png" alt="" width="20" height="20" class="">
                                    </span>
                                </button>
                            </div>
                        </div>
                        <div class="row justify-content-md-begin align-items-center">
                            <style>
                                #gene_list_full {
                                    height: 400px;
                                    overflow-y: scroll;
                                }
                            </style>
                            <select multiple="" class="list-group" id="gene_list_full">
                            </select>
                        </div>
                    </div>
                    <div class="col-sm-1 text-center">

                    </div>

                    <div class="col-sm">
                        <?php// Selected genes from all genes available in highlander?>
                        <div class="row justify-content-md-begin align-items-center">
                            <div class="col-sm-2">
                                <div class="row justify-content-end align-items-center">
                                    <button id="add_gene_from_select_gene_all" class="shadow-sm  btn btn-primary">
                                        <span>
                                            <img src="../img/blue-double-arrow-right.png" alt="" width="20" height="20" class="">
                                        </span>
                                    </button>
                                </div>
                                <div class="row">
                                    <?php echo "<br>";?>
                                </div> 
                                <div class="row justify-content-end align-items-center">
                                    <button id="remove_gene_from_select_gene_all" class="shadow-sm  btn btn-primary">
                                        <span>
                                            <img src="../img/blue-double-arrow-left.png" alt="" width="20" height="20" class="">
                                        </span>
                                    </button>
                                </div>
                            </div>
                            <div class="col-sm-1 text-center">

                            </div>
                            <div class="col-sm justify-content-md-begin " id="card_gene_list_selected">
                                <div class="row justify-content-md-center align-items-center">
                                    <label class="fw-bold text" for="gene_list_selected">Selected Genes</label>
                                </div>
                                <div class="row justify-content-md-center align-items-center">
                                    <style>
                                        #gene_list_selected {
                                            height: 400px;
                                            overflow-y: scroll;
                                        }
                                    </style>
                                    <select multiple="" class="list-group" id="gene_list_selected">
                                    </select>
                                </div>
                                <div id="gene_list_selected_invalid" class="invalid-feedback">
                                    No gene selected! 
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>

            <?php// Retrieve genes from highlander list of gene?>
            
            <div class="collapse multi-collapse card card-body" id="private_list_gene">
                <div class="row g-2 mb-2 ">
                    <div class="col-sm">
                        <div class="row justify-content-md-begin align-items-center">
                        <label class="fw-bold text" for="list_gene_full">Select your list of genes</label>
                            <div class="col-md">
                                <input class="form-control mb-2" id="search_list_gene_all" type="text" placeholder="Search..">
                            </div>
                            <div class="col justify-content-end">
                                <button id="remove_from_search_list_gene_all" class="shadow-sm  btn btn-primary">
                                    <span>
                                        <img src="../img/delete.png" alt="" width="20" height="20" class="">
                                    </span>
                                </button>
                            </div>
                        </div>
                        <div class="row justify-content-md-begin align-items-center">
                            <style>
                                #list_gene_full {
                                    height: 400px;
                                    overflow-y: scroll;
                                }
                            </style>
                            <select multiple="" class="list-group" id="list_gene_full">
                            </select>
                        </div>
                    </div>
                    <div class="col-sm-1 text-center">

                    </div>

                    <div class="col-sm">
                        <?php// Selected genes from highlander list?>
                        <div class="row justify-content-md-begin align-items-center">
                            <div class="col-sm-2">
                                <div class="row justify-content-end align-items-center">
                                    <button id="add_from_select_list_gene" class="shadow-sm  btn btn-primary">
                                        <span>
                                            <img src="../img/blue-double-arrow-right.png" alt="" width="20" height="20" class="">
                                        </span>
                                    </button>
                                </div>
                                <div class="row">
                                    <?php echo "<br>";?>
                                </div> 
                                <div class="row justify-content-end align-items-center">
                                    <button id="remove_from_select_list_gene" class="shadow-sm  btn btn-primary">
                                        <span>
                                            <img src="../img/delete.png" alt="" width="20" height="20" class="">
                                        </span>
                                    </button>
                                </div>
                            </div>
                            <div class="col-sm-1 text-center">

                            </div>
                            <div class="col-sm justify-content-md-begin " id="card_list_gene_selected">
                                <div class="row justify-content-md-center align-items-center">
                                    <label class="fw-bold text" for="list_gene_selected">Selected Genes</label>
                                </div>
                                <div class="row justify-content-md-center align-items-center">
                                    <style>
                                        #list_gene_selected {
                                            height: 400px;
                                            overflow-y: scroll;
                                        }
                                    </style>
                                    <select multiple="" class="list-group" id="list_gene_selected" name="options[]">

                                    </select>
                                </div>
                                <div id="list_gene_selected_invalid" class="invalid-feedback">
                                    No Genes selected! 
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
            <div id="selection_gene_valid" class="valid-feedback">
                Gene list selection done, just one last step!
            </div>
        </div>  
        <div class="d-grid gap-2 d-md-flex justify-content-md-end">
            <button id="next_step_after_gene_selection" class="shadow-sm  btn btn-primary w-25">Next Step</button>
        </div>
    </div>


    <?php// This is the part for filtering Using Highlander?>
    <div id="filtering_selection" class="visually-hidden">
        <div class="card card-body " id="card_drop_filtering">
            <div class="row g-2 mb-2 ">
                <div class="col-sm">
                </div>
                <div class="col-sm-1 text-center">
                <span>
                    <img src="../img/logo_highlander.png" alt="" width="300" height="187" class="">
                </span>
                </div>
                <div class="col-sm">
                </div>
            </div>
            <div class="row justify-content-md-begin">
                <label class="label mb-1 fs-6 fw-bold text ">Please select the filtering to apply to your data</label>
            </div>
            <div class="row">
                <?php echo "<br>";?>
            </div> 
            <?php// Filters PASS or None?>
            <div class="row g-2 mb-2 " id="filterPass">
                <div class="col-sm">
                    <label class="fw-bold text" for="Filter_pass">Filters = PASS</label>
                    <?php echo "&nbsp; &nbsp; &nbsp; &nbsp;";?>
                    <input class="form-check-input" type="radio" id="Filter_pass" name="Filter_pass" value="PASS" checked="" data-bs-toggle="collapse" role="button" aria-expanded="false">
                    <label class="form-check-label" for="Filter_pass">Yes</label>
                    <?php echo "   ";?>
                    <input class="form-check-input" type="radio" id="Filter_pass_no" name="Filter_pass" value="">
                    <label class="form-check-label" for="Filter_pass_no">No</label>
                </div>
            </div>
            <div class="row">
                <?php echo "<br>";?>
            </div> 
            <?php// Allele frequency?>
            <div class="row justify-content-md-begin " id="Allele_freq">
                <div class="col-md-auto">
                    <label class="fw-bold text" for="maf_filter">Allele frequency <</label>
                </div>
                <div class="col col-lg-2">
                    <input type="text" value="0.03" class="form-control " id="maf_filter" name="maf_filter">
                </div>
                <div id="maf_filter_invalid" class="invalid-feedback">
                    Must be a positive numeric value below 1! 
                </div>
            </div>
            <div class="row">
                <?php echo "<br>";?>
            </div> 
            <?php// consensus prediction?>
            <div class="row justify-content-md-begin " id="consensus_prediction">
                <div class="col-md-auto">
                    <label class="fw-bold text" for="consensus_prediction_filter"> Consensus prediction >=</label>
                </div>
                <div class="col col-lg-2">
                    <input type="text" value="0" class="form-control " id="consensus_prediction_filter" name="consensus_prediction_filter">
                </div>
                <div id="consensus_prediction_filter_invalid" class="invalid-feedback">
                    Must be a positive numeric value below 501, with no decimal! Use dot "." for decimal separator, not the decimal ","!
                </div>
            </div>
            <div class="row">
                <?php echo "<br>";?>
            </div> 
            <?php// Variant Type?>
            <div class="row justify-content-md-begin align-items-center" id="variant_type">
                <div class="col-sm justify-content-md-begin ">
                    <div class="row justify-content-md-center align-items-center">
                        <label class="fw-bold text" for="variant_type_full">Variant Type</label>
                    </div>
                    <div class="row justify-content-md-center align-items-center">
                        <select multiple="" class="list-group" id="variant_type_full">
                        </select>
                    </div>
                </div>
                <div class="col-sm-1 text-center">

                </div>
                <div class="col-sm-2">
                    <div class="row">
                        <?php echo "<br>";?>
                    </div> 
                    <div class="row justify-content-end align-items-center">
                        <button id="add_variant_type" class="shadow-sm  btn btn-primary">
                            <span>
                                <img src="../img/blue-double-arrow-right.png" alt="" width="25" height="15" class="">
                            </span>
                        </button>
                    </div>
                    <div class="row">
                        <?php echo "<br>";?>
                    </div> 
                    <div class="row justify-content-end align-items-center">
                        <button id="remove_variant_type" class="shadow-sm  btn btn-primary">
                            <span>
                                <img src="../img/blue-double-arrow-left.png" alt="" width="25" height="15" class="">
                            </span>
                        </button>
                    </div>
                </div>
                <div class="col-sm-1 text-center">

                </div>
                <div class="col-sm justify-content-md-begin " id="card_variant_type_selected">
                    <div class="row justify-content-md-center align-items-center">
                        <label class="fw-bold text" for="variant_type_selected">Variant Type Selected </label>
                    </div>
                    <div class="row justify-content-md-center align-items-center">
                        <select multiple="" class="list-group" id="variant_type_selected">
                        </select>
                    </div>
                    <div id="variant_type_selected_invalid" class="invalid-feedback">
                        No variant type selected! 
                    </div>
                </div>
            </div>
            <?php// Add new annotation here?>
            <div class="row g-2 mb-2 ">
                <div class="col-sm">  
                </div>
            </div>
            

        </div>
        <div class="d-grid gap-2 d-md-flex justify-content-md-end">
            <button id="submit_via_highlander" name="submit_via_highlander"
            class="shadow-sm  btn btn-primary mb-2 w-25">Submit</button>
        </div>
    </div>



    <?php// This is the part for upload data directly via the website?>
    <div id="submit_via_website_card" class="visually-hidden">    
        <div class="input-group mb-5">
            <?php// This is the part for patient data upload?>
            <div class="form-floating w-100 ">
                <div id="noFile" class="invalid-feedback mb-2">
                        Patients File Missing.
                </div>
                <div class="collapse multi-collapse text-center" id="PatientFile">
                    <div class="card card-body " id="card_drop">
                        <li class="list-group-item list-group-item-action" id="initSeq">
                            <div class="form-group" id="search_file_area">
                                <label class="form-label fw-bold" for=formFile>Select File for patients</label>
                                <input class="form-control" type="file" id="formFile">
                            </div>
                            <ul name="uploaded_file" class="list-group collapse" id="uploaded_file">

                            </ul>
                            <div id="uploaded_selection_error" class="invalid-feedback">
                                    Upload Patients File.
                            </div>
                        </li>
                    </div>
                </div>
                <div class="progress" id="progress_bar1">
                    <div class="progress-bar" id="progress_bar_process" role="progressbar" style="width: 0%;" aria-valuenow="0" aria-valuemin="0" aria-valuemax="100"0%>
                    </div>
                </div>
            </div>
            <?php// This is the part for control data upload?>
            <div class="form-floating w-100 ">
                <div id="noFile2" class="invalid-feedback mb-2">
                    Controls File Missing.
                </div>
                <div class="collapse multi-collapse text-center" id="ControlFile">
                    <div class="card card-body " id="card_drop2">
                        <li class="list-group-item list-group-item-action" id="initSeq">
                            <div class="form-group" id="search_file_area2">
                                <label class="form-label fw-bold" for=formFile2>Select File for controls</label>
                                <input class="form-control" type="file" id="formFile2">
                            </div>
                            <ul name="uploaded_file2" class="list-group collapse" id="uploaded_file2">

                            </ul>
                            <div id="uploaded_selection_error2" class="invalid-feedback">
                            Upload Controls File.
                            </div>
                        </li>
                    </div>
                </div>
                <div class="progress" id="progress_bar2">
                    <div class="progress-bar" id="progress_bar_process2" role="progressbar" style="width: 0%;" aria-valuenow="0" aria-valuemin="0" aria-valuemax="100"0%>
                    </div>
                </div>
            </div>
            <?php// This is the part for gene list data upload?>
            <div class="form-floating w-100 ">
                <div id="noFile3" class="invalid-feedback mb-2">
                    Gene list File Missing.
                </div>
                <div class="collapse multi-collapse text-center" id="geneFile">
                    <div class="card card-body " id="card_drop3">
                        <li class="list-group-item list-group-item-action" id="initSeq">
                            <div class="form-group" id="search_file_area3">
                                <label class="form-label fw-bold" for=formFile3>Select File for gene list</label>
                                <input class="form-control" type="file" id="formFile3">
                            </div>
                            <ul name="uploaded_file3" class="list-group collapse" id="uploaded_file3">

                            </ul>
                            <div id="uploaded_selection_error3" class="invalid-feedback">
                                Upload Gene list File.
                            </div>
                        </li>
                    </div>
                </div>
                <div class="progress" id="progress_bar3">
                    <div class="progress-bar" id="progress_bar_process3" role="progressbar" style="width: 0%;" aria-valuenow="0" aria-valuemin="0" aria-valuemax="100"0%>
                    </div>
                </div>
            </div>
        </div>

        <div class="d-grid gap-2 d-md-flex justify-content-md-end">
            <button id="submit_via_website" name="submit_via_website"
            class="shadow-sm  btn btn-primary mb-2 w-25">Submit</button>
        </div>
    </div>
</div>

