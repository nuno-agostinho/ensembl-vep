#!/bin/sh
singularity_dir=singularity-images

mkdir -p $singularity_dir
singularity pull --force --name $singularity_dir/bcftools.sif \
                 docker://quay.io/biocontainers/bcftools:1.13--h3a49de5_0
singularity pull --force --name $singularity_dir/vep.sif \
                 docker://ensemblorg/ensembl-vep

# required for SIFT/PolyPhen-2 protein prediction workflow
singularity pull --force --name $singularity_dir/agat.sif \
                 docker://quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0

docker_dir="../docker"
docker build -f $docker_dir/polyphen2.Dockerfile -t polyphen2
docker save polyphen2 > $docker_dir/polyphen2.tar
singularity build $singularity_dir/polyphen2.sif \
            docker-archive://$docker_dir/polyphen2.tar

docker build -f $docker_dir/sift.Dockerfile -t sift
docker save sift > $docker_dir/sift.tar
singularity build $singularity_dir/sift.sif \
            docker-archive://$docker_dir/sift.tar
