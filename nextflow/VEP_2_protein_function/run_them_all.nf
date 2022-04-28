/*
 * From VCF to protein prediction
 */

nextflow.enable.dsl=2

include { run_vep } from '../VEP/workflows/run_vep.nf'
include { getTranslation } from '../protein_function/nf_modules/run_agat.nf'
include { pph2 } from '../protein_function/nf_modules/run_polyphen2.nf'

process create_subs {
  /*
  Generate amino acid substitutions
  */

  // tag "$vcf"
  // container "biocontainers/bcftools:v1.9-1-deb_cv1"
  // memory '4 GB'
  // errorStrategy 'ignore'

  input:
    path vcf

  output:
    path '*.subs'

  """
  export BCFTOOLS_PLUGINS=/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/Cellar/bcftools/1.12/libexec/bcftools/
  bcftools +split-vep *.vcf.gz -d -A tab -s :missense \
           -f '%CHROM:%POS %Feature %Consequence %Protein_position %Amino_acids\n' \
           > protein.subs
  """
}

// fasta=../protein_function/input/fasta/Homo_sapiens.GRCh38.dna_sm.toplevel.fa
// gtf=../protein_function/input/annotation/Homo_sapiens.GRCh38.105.gtf
// vcf=input/homo_sapiens_GRCh38.vcf.gz
// nextflow run run_them_all.nf --fasta $fasta --gtf $gtf --vcf $vcf --chros 21,22 -profile lsf -resume

workflow {
  // Get translated FASTA
  getTranslation(file(params.gtf), file(params.fasta))

  // Get substitutions
  run_vep()
  create_subs( run_vep.out )

  // For each transcript, predict protein function
  // pph2(translated, subs)
}
