process filter_missense {
  container "ensemblorg/ensembl-vep:release_106.1"

  input:
    path vcf

  output:
    path 'filtered.vcf'

  """
  filter_vep -i $vcf -o filtered.vcf --filter "Consequence is missense_variant"
  """
}
