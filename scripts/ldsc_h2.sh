#!/bin/sh

gwas_file=$1
N=$2
OUTDIR=$3

LDSC_DIR=../scripts/ldsc
HM3_SNPLIST=${LDSC_DIR}/eur_w_ld_chr/w_hm3.snplist
REF_LD_CHR=${LDSC_DIR}/eur_w_ld_chr/

TRAIT=$(basename $gwas_file | cut -d. -f1)

python $LDSC_DIR/munge_sumstats.py --sumstats $gwas_file --N $N --out $OUTDIR/$TRAIT --merge-allele $HM3_SNPLIST --ignore BETA,OR,SE,BETA_STD

python $LDSC_DIR/ldsc.py --h2 ${OUTDIR}/$TRAIT.sumstats.gz --ref-ld-chr $REF_LD_CHR --w-ld-chr $REF_LD_CHR --out ${OUTDIR}/$TRAIT
