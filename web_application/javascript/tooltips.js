$(document).ready(function () {
    tippy.setDefaultProps({
        theme: 'light',
        animation: 'scale',
        followCursor: 'horizontal',
        duration: 150,
        arrow: true
    });

    // Starting your submission
    tippy('#email', {
        content: 'E-mail to which the results download notification will be sent'
    });
    tippy('#pi', {
        content: 'To which group do you belong?'
    });
    tippy('#Highlander_import', {
        content: 'Construct your cohort and analysis using data available wihin your Highlander session. By clicking on the "Next" button, you will be asked to login into your Highlander session. Make sure you have a Highlander account!'
    });
    tippy('#manual_import', {
        content: 'If you already have generated the tsv files necessary for an Excalibur analysis, you can upload them by choosing this option and clicking the "Next" button'
    });
    tippy('#gene_list', {
        content: 'Limit your analysis to  a particular set of genes'
    });
    tippy('#do_coverage', {
        content: 'To enhance reliability of your results. This steps check the BAM files of each individu in your cohort to ensure a good coverage and remove position and/or individu with a poor coverage'
    });
    tippy('#advancedOption', {
        content: 'DONT change these options unless you have very specific needs and hypothesis '
    });

    // PART FOR SUBMISSION USING HIGHLANDER 
    // This is the part for Database selection
    tippy('#card_drop_Database_selection', {
        content: 'Excalibur has been designed to handle exomes, panels and genomes. Note that we only tested it on exomes. We strongly recommend not to use it on rnaseq!'
    });
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
    });
    // This is the part for filtering Using Highlander
    tippy('#filterPass', {
        placement: 'bottom',
        maxWidth: 400,
        content: 'If filtering was performed on the VCF by the pipeline, this field will contain the name of the filter(s) that the variant DIDN\'T pass, or PASS if the variant passed all filters. e.g. if this field show \"lowQual,GATKStandard\" it means that this variant doesn\'t pass \"lowQual\" and \"GATKStandard\" filters. If this field is empty or equal to \'-\', it means that no filtering was performed by the pipeline. So this field depends heavily on the variant caller used and the pipeline implementation. The following description is for GATK variant caller only: In a perfect world, the Confidence field would be based on a complete model for all error modes present in the data used to call. Unfortunately, we are still far from this ideal, and we have to use orthogonal approaches to determine which called sites, independent of Confidence, are machine errors and which are real variants. Whatever approach is used to filter the variants, the VCFs produced by the GATK carry both the PASSing filter records (the ones that are good have PASS in their FILTER field) as well as those that fail (the filter field is anything but PASS or a dot). A variant won\'t pass GATK Standard filter if it met at least one of those conditions: - variant_confidence_by_depth < 2. - fisher_strand_bias > 60. - mapping_quality < 40. - haplotype_score > 13. - rank_sum_test_read_mapping_qual < -12.5. - rank_sum_test_read_pos_bias < -8.'
    });
    tippy('#Allele_freq', {
        content: 'gnomad_wes_afgnomad_wes_af : Alternative allele frequency in the whole gnomAD exome samples (125,748 samples)'
    });
    tippy('#consensus_prediction', {
        placement: 'bottom',
        maxWidth: 600,
        content: 'Consensus between prediction of all available software.  The base score reflects the number software that predict the variant to be damaging:  Each of the following software having a DAMAGING prediction add +1: Mutation Taster, FATHMM, FATHMM-XF, Polyphen2 (HDIV), Provean, SIFT4G, Mutation Assessor, MCAP, LRT, Lists2, Deogen, ClinPred, BayesDel (with MaxMAF), PrimateAI and MetaSVM. Each of the following scores add +1 when above a certain threshold: CADD phred > 20, VEST > 0.5, REVEL > 0.5, MVP > 0.75 and MutPred > 0.75  So max 20 if all software agree.  If the variant could affect splicing, a +1 or +2 could be added to the base score, if ada and/or rf predictions are AFFECTING_SPLICING.  To this, value is added depending on certain variant impact (annotations from SnpEff): +200 if any prediction indicates that splicing is affected. +300 for SPLICE_SITE_ACCEPTOR, SPLICE_SITE_DONOR, the 2 first or 2 last positions of an exon, the 3rd/4th/5th positions in 3\' intron of an exon (before the STOP). +400 for FRAME_SHIFT, STOP_GAINED, STOP_LOST, START_LOST or RARE_AMINO_ACID. +500 for high impact structural variants. So, filtering on consensus_prediction > 0 should yield variants potentially damaging. Choosing a higher value like consensus_prediction > 5 should yield variants probably damaging.'

    });
    tippy('#variant_type', {
        placement: 'bottom',
        maxWidth: 600,
        content: 'Type of variant : SNV (Single Nucleotide Variant), MNV (Multiple Nucleotides Variant), INS (Small Insertion), DEL (Small Deletion), SV (Structural Variant)'
    });


    //  PART FOR SUBMISSION USING WEBSITE
    tippy('#card_drop', {
        content: 'Search Patients File Here'
    });
    tippy('#card_drop2', {
        content: 'Search Controls File Here'
    });
    tippy('#card_drop3', {
        content: 'Search Gene list File Here'
    });


});