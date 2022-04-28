#!/usr/bin/env nextflow

/*
 * AGAT: Another GTF/GFF Analysis Toolkit
 */

process getTranslation {
  /*
  Translate nucleotide FASTA sequences based on GTF features

  Returns
  -------
  Returns 1 file:
      1) Protein FASTA sequence 'translated.fa'
  */

  tag "${gtf}"
  container "quay.io/biocontainers/agat:0.9.0--pl5321hdfd78af_0"
  // container "${params.singularity_dir}/agat.simg"
  memory '20 GB'
  publishDir "${params.outdir}"

  input:
    path gtf
    path fasta

  output:
    path 'translated.fa'

  """
  agat_sp_extract_sequences.pl -g $gtf -f $fasta --protein -o translated.fa
  """
}
