#!/usr/bin/env nextflow

log.info """\
 =========================================

 A Nextflow script to randomise the GO enrichment (v1.0.0)

 -----------------------------------------

 Authors:
   - Chris Wyatt <c.wyatt@ucl.ac.uk>

 -----------------------------------------

 Copyright (c) 2024

 =========================================""".stripIndent()


def errorMessage() {
    log.info"""
    =============
    synteny error
    =============
    You failed to provide the input parameter
    Please provide this as follows:
      --input /full/path/to/sample/file
    or use the test profile:
      --profile test
    The pipeline has exited with error status 1.
    """.stripIndent()
    exit 1
}



include { GO_JUNCTIONS_INVER_DIST } from './modules/local/go_junctions_inver_score4_dist.nf'
include { GO_SUMMARISE } from './modules/local/go_random_summarise.nf'


repli_ch = Channel.of(params.replicates)


workflow {

  //gene.scores.txt

  species_inver = Channel.fromPath(params.input)
                      .map { file -> 
                          def species = file.name.split('\\.')[0]  // Extract species name
                          tuple(species, file)  // Return as tuple (file, species)
                      }
 
  input_beds = Channel.fromPath(params.beds).collect()

  go_folder = Channel.fromPath(params.go)

  // Split the params.size string into a list of separate entries (1,3)
  size = Channel.from(params.size.split(','))

  go_and_summary_inver = go_folder.combine(species_inver)

  mergedChannel_inver_dist = go_and_summary_inver.combine(size)

  GO_JUNCTIONS_INVER_DIST ( mergedChannel_inver_dist , input_beds , repli_ch.first() )

  //GO_JUNCTIONS_INVER_DIST.out.go_pvals.view()

  GO_SUMMARISE ( GO_JUNCTIONS_INVER_DIST.out.go_pvals )

}



workflow.onComplete {
    println(workflow.success ? "\nDone! Check results in $params.outdir/ \n" : "Hmmm .. something went wrong\n")
}