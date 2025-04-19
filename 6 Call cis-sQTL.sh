


/public/home/liangyan/software/QTLtools cis --vcf filter.chr.dose.vcf.gz \
--bed celltype_chr_phenotype.txt.gz \
--cov celltype_ok_PC.txt \
--permute 1000 --normal --grp-best --out celltype_chr_permutation.txt

cat celltype_*_permutation.txt | gzip -c > celltype_permutations_full.txt.gz

Rscript runFDR_cis.R celltype_permutations_full.txt.gz  0.05 celltype_permutations_full.txt

/public/home/liangyan/software/QTLtools cis --vcf filter.chr.dose.vcf.gz \
                --bed celltype_chr_phenotype.txt.gz \
                --cov celltype_ok_PC.txt --normal --grp-best \
                --mapping celltype_permutations_full.txt.thresholds.txt \
                --out conditions_celltype_chr.txt







