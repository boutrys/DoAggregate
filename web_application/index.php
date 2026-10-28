<?php include('pages/Nav.php'); ?>
<div class="container-fluid">
    <p class="fw-bold fs-5">About</p>

    <p>This website allows you to perform association testing on a cohort of patients and controls using aggregation test  
    <?php //TODO add ref to Excalibur paper <a href="https://doi.org/10.1038/s41586-021-03819-2" target="_blank" rel="nofollow">AlphaFold v2.2.2</a>?>, 
    running on a dedicated node of the SSS Cluster. 
    <br>
    The different pages (accessible through the black band) can be used to:
    </p>

    <ul>
        <li><a href="pages/Form.php" rel="nofollow"><strong>Submit </strong> cohort</a> to the cluster for association testing. When your analysis will be finished, you'll receive an email inviting you to download the results.</li>
        <li><a href="pages/Live.php" rel="nofollow">Display the <strong>job queue </strong></a> of currently running predictions. You can check it regularly to see when your analysis will run.</li>
        <li> Click on "How to use" to see how to fill <a href="pages/Form.php" rel="nofollow">the form</a> with the requested fields.</li>
    </ul>


    <p class="fw-bold fs-5">Citing this work
    <ul>
        Any publication that discloses findings arising from Excalibur must cite the original publication (still under revision, we will add the article reference later). 
        <?php //<a href="https://github.com/deepmind/alphafold/#citing-this-work" target="_blank" rel="nofollow">cite</a> their original <a href="https://doi.org/10.1038/s41586-021-03819-2" target="_blank" rel="nofollow">publication</a>.
        ?>
    </ul>
    </p>


    <p class="fw-bold fs-5">Technical seminar
    <ul>
        You can download the <a href="data/Genetic_seminar_07_02_23.mp4" rel="nofollow">1h recording</a> and the <a href="data/2022_12_20_genetic_seminar_Excalibur_Web.pptx" rel="nofollow">power point slides</a> of the genetic seminar given by Simon Boutry about Excalibur : a new gene/pathway/GO term prioritization tool for rare-variant. Theory and how to use the webpage. The second part of the video gives you a live demonstration of the use of the website and results exploration.
    </ul>
    </p>


    <div class="card card-body " id="card_drop_sample_selection">
        <div class="row g-2 mb-2 ">
            <div class="col-sm text-center">
                <button id="how_to_use_button" class="btn btn-outline-primary btn-sm fw-bold fs-5" type="button"
                            data-bs-toggle="collapse" data-bs-target="#how_to_use" aria-expanded="false"> How to use
                </button>
            </div>
            <div class="col-sm text-center">
                <button id="Output_description_button" class="btn btn-outline-primary btn-sm fw-bold fs-5" type="button"
                            data-bs-toggle="collapse" data-bs-target="#Output_description" aria-expanded="false">Results description
                </button>
            </div>
            <div class="col-sm text-center">
                <button id="FAQ_button" class="btn btn-outline-primary btn-sm fw-bold fs-5" type="button"
                            data-bs-toggle="collapse" data-bs-target="#FAQ" aria-expanded="false">FAQ
                </button>
            </div>
            <div class="col-sm text-center">
                <button id="troubleshooting_button" class="btn btn-outline-primary btn-sm fw-bold fs-5" type="button"
                            data-bs-toggle="collapse" data-bs-target="#troubleshooting" aria-expanded="false">Troubleshooting
                </button>
            </div>
        </div>
        
        <?php// How to use Excalibur?>
        <div class="card card-body visually-hidden" id="how_to_use">
            In this section we describe all the fields to fill in order to perfom an analysis.
            <br>
            You can look at the video at this link , or download the slides of the presentation here.
            <br>
            <br>
            <p class="fw-bold fs-6"><b>Starting your submission</b>
            <ul>
                <li>Enter your <b>email address (and make sure it is valid!)</b> and select your lab in the list (<a href="mailto:simon.boutry@uclouvain.be" target="_new">contact us</a> if your lab is missing). <br> 
                If the email adress is not valid, you will never receive your results!</li>
                <li>Choose how to import your data</li>
                <ul>
                    <li><a href="#importing_via_Highlander_button" rel="nofollow"><b>Import data via Highlander</b></a> if you wanna build your cohort using data available within your Highlander profile</li>
                    <li><a href="#importing_your_data_button" rel="nofollow"><b>Import your data</b></a> allows you to upload via the website your patient and control files and optionnaly a gene list. Make sure all file respect the criterias described <a href="#importing_via_Highlander" rel="nofollow">here</a></li>
                </ul>
                <li>You can restrict your analysis to a particular <b>gene list</b>. If you use this option, you will be given the opportunity within your submission to give the list of gene that should be analyzed</li>
                <li>You can <b>clean your data</b> using the coverage of each position. This is <b>recommended</b> to improve reliability of results. Note that cleaning your data is time consuming.</li>
                <li>For <b>advanced users</b>, you can modify some <b>advanced parameters</b> for cleaning your data and aggregation tests. By default, we <b>encourage</b> the user <b>not to modify</b> them.</li>
            </ul>
            By clicking the Next button you will be redirected to one of the following
            </p>

            <br>

            <div class="row g-2 mb-2 ">
                <div class="col-sm text-center">
                    <button id="importing_via_Highlander_button" class="btn btn-outline-primary btn-sm fw-bold" type="button"
                                data-bs-toggle="collapse" data-bs-target="#importing_via_Highlander" aria-expanded="false"> Import data via Highlander
                    </button>
                </div>
                <div class="col-sm text-center">
                    OR
                </div>
                <div class="col-sm text-center">
                    <button id="importing_your_data_button" class="btn btn-outline-primary btn-sm fw-bold" type="button"
                                data-bs-toggle="collapse" data-bs-target="#importing_your_data" aria-expanded="false"> Import your data
                    </button>
                </div>
            </div>
            
            <br>
            
            <div class="card card-body visually-hidden" id="importing_via_Highlander">
                <p class="fw-bold fs-6" id="importing_via_Highlander"><b>Import data via Highlander</b> </p>
                
                <ul>
                    Step 1) <b>Connect to your Highlander profile</b>. Make sure you have an Highlander profile, if not, please contact <a href="mailto:raphael.helaers@uclouvain.be" target="_new">Raphaël Helaers</a>. Please enter a valid username and password!
                </ul>
                <ul>
                    Step 2) <b>Choose the database for your analysis</b>. Note that Excalibur has been designed to handle exomes, panels and genomes. So far, we only tested it on exomes. We strongly recommend not to use it on rnaseq!
                </ul>
                <ul>
                    Step 3) <b>Select patients and controls from </b>.
                    <ul>
                        <li>All samples availabl in Highlander OR</li>
                    </ul>
                    <ul>
                        <li>Your sample lists within your Highlander profile.</li>
                    </ul>
                </ul>
                <ul>
                    Optional step) If you want restrict your analysis on a particular gene list, you can <b>Select genes from </b>.
                    <ul>
                        <li>All genes availabl in Highlander OR</li>
                    </ul>
                    <ul>
                        <li>Your gene lists within your Highlander profile.</li>
                    </ul>
                </ul>
                <ul>
                    Final Step) <b>Filtering options</b>. Choose the basic filtering options to apply to both your patients and controls data. For more information about each filtering option available, please read the FAQ section of this page.
                </ul>
                
            </div>

            <div class="card card-body visually-hidden" id="importing_your_data">
                <p class="fw-bold fs-6" id="importing_your_data"><b>Import your data</b></p>
                <ul>
                You must select two <b>TSV files</b>: one for with your patients data, and one with your controls data. Both files must have at least these <b>mandatory columns</b> (see how to construct TSV file with Highlander bellow) :
                    <ul>
                        <li> <b>sample</b> : unique identifier for an individual. </li> 
                        <li> <b>pathology</b> : unique identifier for all individuals with same pathology. </li> 
                        <li> <b>chr</b> : chromosome. </li> 
                        <li> <b>pos</b> : position. </li>
                        <li> <b>reference</b> : the reference base that vary in the samples, or in the population in general. </li>
                        <li> <b>alternative</b> : the alternative base that vary in the samples, or in the population in general. </li>
                        <li> <b>zygosity</b> : zygosity of the variant, encoded as Heterozygous or Homozygous (alternative). </li> 
                        <li> <b>gene_symbol</b> : gene annotation form <a href="https://www.genenames.org/" target="_blank" rel="nofollow">HUGO database</a>. </li> 
                        <li> <b>gnomad_wes_af</b> : alternative allele frequency in the Gnomad Whole Exome Samples. </li> 
                        Example of TSV file for Excalibur
                        <br> <br/>
                        <img src="/excalibur/img/example_tsvFile.jpg" alt="" width="800" height="300" class="d-inline-block align-text-center"> 
                        <br> <br/>
                    </ul>
                </ul>
                <ul>
                    Make sure that the exact same pipeline and parameters have been used to procude the two files (e.g. alignement, variant calling, annotation, ...)
                </ul>
                <ul>
                Optionnaly, you can upload a <b>Gene list</b> as txt or tsv file with no column name, and one gene per line. Note that gene name should be HUGO annotations (see gene_symbol in Highlander or <a href="https://www.genenames.org/" target="_blank" rel="nofollow">HUGO website</a>. Here is an example of such file
                <br> <br/>
                <img src="/excalibur/img/example_geneList.jpg" alt="" width="80" height="243" class="d-inline-block align-text-center"> 
                </ul>
            </div>
            
            <br>

            <p class="fw-bold fs-6"><b>Submit your analysis</b>
            <ul>
                <li>After using the submit button, you will be lead to the job queue.</li>
                <li>When the analysis is finished, a link to download a zipped file of the results will be sent by email.</li>
            </ul>
            </p>
        </div>

        <?php// Output description?>
        <div class="card card-body visually-hidden" id="Output_description">
            <ul> 
                In this section we describe the results <b>folders (in bold)</b> and files that Excalibur will return. Note that if the name of a file is <i> name of file (in italic)</i> this means that this file might not be present in your results, but its perfectly normal. For example, if no genes with significant aggregation was found, then you won't see a file called <i> significant general info </i> in your folder <b> gene </b>. 
            </ul>

            <p><b>Folder "General"</b>
            <ul>
                <li>
                    <b>Folder Gene set analysis</b> 
                    <ul>
                        <li>
                            <b>Folder Go over representation</b>
                            <ul>
                                <li>
                                    <b>Folder BP</b> : Biological Process
                                    <ul>
                                        <li>
                                            Files BP Go over representation (excell + <i> short + simplified + plot</i>)
                                        </li>
                                    </ul>
                                </li>
                                <li>
                                    <b>Folder CC</b> : Cellular component
                                    <ul>
                                        <li>
                                            Files CC Go over representation (excell + <i> short + simplified + plot</i>)
                                        </li>
                                    </ul>
                                </li>
                                <li>
                                    <b>Folder MF</b> : Molecular function
                                    <ul>
                                        <li>
                                            Files MF Go over representation (excell + <i> short + simplified + plot</i>)
                                        </li>
                                    </ul>
                                </li>
                            </ul>
                        </li>
                        <li>
                            <b>Others analysis</b> : results of pathway over representation using consensusPathDB
                            <ul>
                                <li>
                                    Files CPDB over representation (excell + <i> short + simplified + plot</i>)
                                </li>
                            </ul>
                        </li>
                        <li>
                            <i>File gene not annotated goTerm</i>
                        </li>
                        <li>
                            <i>File gene not annotated goTerm</i>
                        </li>
                    </ul>

                </li>
                <li>
                    File Global info analysis : contains global information about each type of genetic region explored (i.e gene, pathway, GO terms), with the following annotations :
                    <ul>
                        <li>
                            Name of experiment
                        </li>
                        <li>
                            Nbr region
                        </li>
                        <li>
                            Nbr region analyzed
                        </li>
                        <li> 
                            Nbr significant region
                        </li>
                        <li>
                            Nbr region not analyzed
                        </li>
                        <li>
                            Proportion Nbr significant region / Nbr region analyzed
                        </li>
                        <li>
                            Proportion Nbr region analyzed / Nbr region
                        </li>
                        <li>
                            Nbr position unique
                        </li>
                        <li>
                            Nbr position unique patient
                        </li>
                        <li>
                            Nbr position unique control
                        </li>
                        <li>
                            Nbr variant unique
                        </li>
                        <li>
                            Nbr variant unique patient
                        </li>
                        <li>
                            Nbr variant unique control
                        </li>
                        <li>
                            Total variant patient
                        </li>
                        <li>
                            Total variant control
                        </li>
                    </ul>
                </li>
            </ul>
            </p>

            <p><b>Folder "Gene"</b>
            <ul>
                <li>
                    <b>Folder <i>Data significant region</i></b> 
                    <ul>
                        <li>
                            Contains all <i>files</i> of significant region raw (variants and annotations) data used.
                        </li>
                    </ul>
                </li>
                <li>
                    File general info : all genetic region ranked by their Excalibur qvalue, and annotated with :
                    <ul>
                        <li>
                            Ranking
                        </li>
                        <li>
                            Excalibur qvalue
                        </li>
                        <li>
                            Best qvalue
                        </li>
                        <li> 
                            Best Stat test
                        </li>
                        <li>
                        Sig Stat test

                        </li>
                        <li>
                        Nbr sig test

                        </li>
                        <li>
                        Nbr test applied

                        </li>
                        <li>
                        Name genetic region

                        </li>
                        <li>
                        Nbr patient with mutation

                        </li>
                        <li>
                        Nbr control with mutation

                        </li>
                        <li>
                        Proportion patients with mutation

                        </li>
                        <li>
                        Proportion controls with mutation

                        </li>
                        <li>
                        Total nbr variant

                        </li>
                        <li>
                        Total nbr variant in patients

                        </li>
                        <li>
                            Total nbr variant in controls
                        </li>
                        <li>
                        Nbr unique variant

                        </li>
                        <li>
                        Nbr unique variant in patient

                        </li>
                        <li>
                        Nbr unique variant in controls

                        </li>
                        <li>
                        Nbr unique variant in both

                        </li>
                        <li>
                        Total nbr heterozygous variant patient

                        </li>
                        <li>
                        Total nbr homozygous variant patient

                        </li>
                        <li>
                        Total nbr heterozygous variant control

                        </li>
                        <li>
                            Total nbr homozygous variant control
                        </li>
                        <li>
                        Proportion of reference homozygous

                        </li>
                        <li>
                        Proportion of heterozygous

                        </li>
                        <li>
                        Proportion of alternatif homozygous

                        </li>
                    </ul>
                </li>
                <li>
                    <i>File Significant general info</i> : same as file general info, but only with genetic region having a Excalibur qvalue < 0.05.
                </li>
                <li>
                    <i>File general info not analyzed</i> : same as file general info, but for all region that were not analyzed.
                </li>
            </ul>
            </p>

            <p><b><i>Folder "Over represented pathway" </i></b>
            <ul>
                If gene set analysis found over represented pathways, then this folder will be created and will have the exact sub folders and files as folder <b> gene </b> except it will be for pathways 
            </ul>
            </p>
            <p><b><i>Folder "Over represented GO" </i></b>
            <ul>
                If gene set analysis found over represented GO terms, then this folder will be created and will have the exact sub folders and files as folder <b> gene </b> except it will be for GO terms 
            </ul>
            </p>
        </div>

        <?php// FAQ?>
        <div class="card card-body visually-hidden" id="FAQ">
            In this section we describe some of the most common questions regarding the use of Excalibur website.
            <br>
            <br>

            <p><b>How to create tsv files necessary for Excalibur using Highlander?</b>
            <ul>
                <li>
                    <b>Select the right database and selection of columns</b> : make sure your selection of columns includes the mandatory columns (see section Import your data). We discourage users to include too many columns in your selection of columns to avoid too big file issues (see Troubleshooting).
                    <br> <br/>
                    <img src="/excalibur/img/example_highlander_select_database.jpg" alt="" class="d-inline-block align-text-center"> 
                    <br> <br/>
                </li>
                <li>
                    <b>Run your filter </b>.
                    <br> <br/>
                    <img src="/excalibur/img/example_highlander_run_filter.jpg" alt="" class="d-inline-block align-text-center"> 
                    <br> <br/>
                </li>
                <li>
                    <b>Download TSV </b>. Once Highlander has retrieve your variants, you can select the "Tools" option, and within it, you can select the "export to tsv" option.
                    <br> <br/>
                    <img src="/excalibur/img/example_highlander_download_tsv.jpg" alt="" class="d-inline-block align-text-center"> 
                    <br> <br/>
                </li>
                <li>
                    <b>Download your patients and controls separately! </b>
                </li>
            </ul>
            </p>

            <p><b>What are these filtering options?</b>
            <ul>
                <li> 
                    <b>Filters = PASS Yes</b> : by default, only variants that pass all filters will be kepts for the analysis. In highlander, Filters is defined as : <br> If filtering was performed on the VCF by the pipeline, this field will contain the name of the filter(s) that the variant DIDN'T pass, or PASS if the variant passed all filters. e.g. if this field show "lowQual,GATKStandard" it means that this variant doesn't pass "lowQual" and "GATKStandard" filters. If this field is empty or equal to '-', it means that no filtering was performed by the pipeline. So this field depends heavily on the variant caller used and the pipeline implementation. The following description is for GATK variant caller only: In a perfect world, the Confidence field would be based on a complete model for all error modes present in the data used to call. Unfortunately, we are still far from this ideal, and we have to use orthogonal approaches to determine which called sites, independent of Confidence, are machine errors and which are real variants. Whatever approach is used to filter the variants, the VCFs produced by the GATK carry both the PASSing filter records (the ones that are good have PASS in their FILTER field) as well as those that fail (the filter field is anything but PASS or a dot). A variant won't pass GATK Standard filter if it met at least one of those conditions: - variant_confidence_by_depth < 2. - fisher_strand_bias > 60. - mapping_quality < 40. - haplotype_score > 13. - rank_sum_test_read_mapping_qual < -12.5. - rank_sum_test_read_pos_bias < -8.
                </li>
                <li>
                    <b>Allele frequency < 0.03 </b> : by default, only variants with an allele frequency below 0.03 will be kept for the analysis. In highlander, gnomad_wes_afgnomad_wes_af is definded as : <br> nomad_wes_afgnomad_wes_af : Alternative allele frequency in the whole gnomAD exome samples (125,748 samples)
                </li>
                <li>
                    <b>Consensus prediction >= 0 </b> : by default, only variants with a consensus prediction equal or above zero will be kept for the analysis. In highlander, consensus_prediction is defined as : <br> Consensus between prediction of all available software.  The base score reflects the number software that predict the variant to be damaging:  Each of the following software having a DAMAGING prediction add +1: Mutation Taster, FATHMM, FATHMM-XF, Polyphen2 (HDIV), Provean, SIFT4G, Mutation Assessor, MCAP, LRT, Lists2, Deogen, ClinPred, BayesDel (with MaxMAF), PrimateAI and MetaSVM. Each of the following scores add +1 when above a certain threshold: CADD phred > 20, VEST > 0.5, REVEL > 0.5, MVP > 0.75 and MutPred > 0.75  So max 20 if all software agree.  If the variant could affect splicing, a +1 or +2 could be added to the base score, if ada and/or rf predictions are AFFECTING_SPLICING.  To this, value is added depending on certain variant impact (annotations from SnpEff): +200 if any prediction indicates that splicing is affected. +300 for SPLICE_SITE_ACCEPTOR, SPLICE_SITE_DONOR, the 2 first or 2 last positions of an exon, the 3rd/4th/5th positions in 3' intron of an exon (before the STOP). +400 for FRAME_SHIFT, STOP_GAINED, STOP_LOST, START_LOST or RARE_AMINO_ACID. +500 for high impact structural variants. So, filtering on consensus_prediction > 0 should yield variants potentially damaging. Choosing a higher value like consensus_prediction > 5 should yield variants probably damaging.
                </li>
                <li>
                    <b>Variant type Selected = SNV</b> : by default, only Single Nucleotide Variant (SNV) will be kept for the analysis. In highlander, variant_type is defined as : <br> Type of variant : SNV (Single Nucleotide Variant), MNV (Multiple Nucleotides Variant), INS (Small Insertion), DEL (Small Deletion), SV (Structural Variant)
                </li>
            </ul>
            </p>

            <p><b>Why some genetic region were not analyzed?</b>
            <ul>
                Explain
            </ul>
            </p>
        </div>

        <?php// Troubleshooting?>
        <div class="card card-body visually-hidden" id="troubleshooting">
            In this section we describe some of the most common questions regarding the use of Excalibur website.
            <br>
            <br>
            <p><b>My file is too big!</b>
            <ul>
                If your input files are too big, you will see a message displayed on the submission form. When this happens, please contact <a href="mailto:simon.boutry@uclouvain.be" target="_new">Simon Boutry</a>. He will check what went wrong and eventually correct problems before you resubmit.
            </ul>
            </p>

            <p><b>My job is not anymore in the queue but I did not receive the email with the zipped results!</b>
            <ul>
                (1) Something went wrong during Excalibur computation, and it crashed.
                <br>
                (2) Everything went well, but the automated email was blocked (it seems to happen rarely).
                <br>
                In both cases, please contact <a href="mailto:simon.boutry@uclouvain.be" target="_new">Simon Boutry</a>.
            </ul>
            </p>
            
            <p><b>Does this tool install anything on my computer?</b>
            <ul>
                No, everything is computed using the cluster of the de Duve institute.
                <br>
                At the end of Excalibur analysis, a zip-archive with the results will be automatically sent to you by email.
            </ul>
            </p>
            
            <p><b>Who should I contact if I have questions or problems?</b>
            <ul>
                You can send a mail to <a href="mailto:simon.boutry@uclouvain.be" target="_new">Simon Boutry</a>.
            </ul>
            </p>
        </div>
    </div>
</div>