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
  container "${params.singularity_dir}/agat.simg"
  memory '8 GB'  

  input:
    path gtf
    path fasta

  output:
    path 'translated.fa'

  """
  agat_sp_extract_sequences.pl -g $gtf -f $fasta --protein -o translated.fa
  """
}
