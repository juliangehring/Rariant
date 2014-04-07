#!/bin/sh
tp53_region="chr17:7571720-7590863"
bam_dir="/ebi/research/huber/users/jgehring/projects/39-error-profiles/data/platinum/align/"
out_dir="/ebi/research/huber/users/jgehring/projects/38-rarient/package/Rariant/inst/platinum/"

## control
samtools view -bh -q 5 ${bam_dir}/NA12877_S1.bam ${tp53_region} > ${out_dir}/control.bam
samtools index ${out_dir}/control.bam

## test
samtools view -bh -q 5 ${bam_dir}/NA12878_S1.bam ${tp53_region} > ${out_dir}/test.bam
samtools index ${out_dir}/test.bam

## test2
samtools view -bh -q 5 ${bam_dir}/NA12881_S1.bam ${tp53_region} > ${out_dir}/test2.bam
samtools index ${out_dir}/test2.bam

## mix
samtools merge ${out_dir}/mix.bam ${out_dir}/test.bam ${out_dir}/test2.bam
samtools index ${out_dir}/mix.bam

## slim
#samtools view -h mix.bam | cut -f 1-12 | samtools view -Sbh - > out.bam
