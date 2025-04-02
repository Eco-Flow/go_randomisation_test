process GO_SUMMARISE {

    label 'process_single'
    container = 'quay.io/ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'
    publishDir "$params.outdir/output_data/go_summary/" , mode: "${params.publish_dir_mode}", pattern:"*go_term_occurrence_summary.txt"
    publishDir "$params.outdir/output_data/go_summary/" , mode: "${params.publish_dir_mode}", pattern:"*bonferroni_pvalues.txt"
    
    input:
    path(go_tables)

    output:
    path( "*go_term_occurrence_summary.txt" ), emit: go_summary_occurance, optional:true
    path( "*bonferroni_pvalues.txt" ), emit: go_summary_top50, optional:true

    script:
    """

    summarise_go.pl
 
    """
    
}
