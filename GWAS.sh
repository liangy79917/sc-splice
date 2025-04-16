




#######从贝服务器上把需要的文件copy，并且计算好MAF
cp /data/home/lingw/workSpace/reference/1000Genomes/hg38/snp-only-as-on/ERZ822766_merged-snp-all.* /data4/liangyan/test/test

##/data4/liangyan/test/test
grep EUR ID.txt |awk '{print $1,$1}' > EUR_ID && wc -l EUR_ID
plink --allow-no-sex --keep-allele-order --bfile ERZ822766_merged-snp-all --keep EUR_ID --make-bed --out EUR_ID

gzip EUR_ID.bed
gzip EUR_ID.bim
gzip EUR_ID.fam

##/home/liang_y/project/dor/gwas
plink --allow-no-sex --keep-allele-order --bfile EUR_ID --freq --out EUR_ID

awk '{if(NR>1) print $1}' EUR_ID.frq > col_1
awk '{if(NR>1) print $2}' EUR_ID.frq|awk -F '-' '{print $2}' > col_2
awk '{if(NR>1) print $3}' EUR_ID.frq > col_3
awk '{if(NR>1) print $4}' EUR_ID.frq > col_4
awk '{if(NR>1) print $5}' EUR_ID.frq > col_5

paste -d '_' col_1 col_2 col_3 col_4 > col_1234
paste -d '_' col_1 col_2 col_4 col_3 > col_1243
paste -d '\t' col_1234 col_1243 col_5 > col_12345

awk -v OFS='\t' '{print $1,$2,$3,1-$3}' col_12345 > col_123456
sed '1i SNP_MAF_1\tSNP_MAF_2\tSMP_1\tSMP_2' col_123456 > col_123456_ok

# grep 1-816376 EUR_ID.frq  ##rs28527770
#    1     1-816376    C    T       0.1203     1006
# grep 1_816376 col_12345_ok

# 1_816376_C_T    1_816376_T_C    0.1203  C ##第一列首个字母是MAF

####总结下MAF的问题
这里的MAF是从旺哥的文件中生成的MAF，我和GWAS数据对应了下，发现仍然有不少对不上，因为我想后面的分析和TWAS都是取model有的SNP，所以这里暂时不处理MAF，只处理hg38位置的问题
这里的MAF是从旺哥的文件中生成的MAF，我和GWAS数据对应了下，发现仍然有不少对不上，因为我想后面的分析和TWAS都是取model有的SNP，所以这里暂时不处理MAF，只处理hg38位置的问题

##但是后来我想maf的核对有利于理解A1和A2，还有利于验证SNP对不对，所以搞下，但是这里注意

> head(data_1)
     SNP_MAF_1   SNP_MAF_2   SMP_1   SMP_2
        <char>      <char>   <num>   <num>
1: 1_16103_G_T 1_16103_T_G 0.03877 0.96123
2: 1_51479_A_T 1_51479_T_A 0.19180 0.80820
3: 1_51898_A_C 1_51898_C_A 0.14210 0.85790
4: 1_51928_A_G 1_51928_G_A 0.14610 0.85390
5: 1_51954_C_G 1_51954_G_C 0.00000 1.00000
6: 1_54490_A_G 1_54490_G_A 0.17500 0.82500

data_1_1 <- data_1[, .(SNP_MAF_1, SMP_1)]
data_1_2 <- data_1[, .(SNP_MAF_2, SMP_2)]

> head(data_2)
     chr    pos effect    ref    beta     SE         P        SNP      SNP_POS
   <int>  <int> <char> <char>   <num>  <num>     <num>     <char>       <char>
1:     1 911428      T      C  0.0310 0.0278 0.2651000  rs4475691 1_911428_T_C
2:     1 911484      C      G  0.0322 0.0282 0.2530000   rs950122 1_911484_C_G
3:     1 978509      A      G -0.0643 0.0213 0.0025360  rs2340596 1_978509_A_G
4:     1 979472      C      G -0.0638 0.0212 0.0026700 rs13303368 1_979472_C_G
5:     1 979560      C      T -0.0669 0.0212 0.0015970 rs13303033 1_979560_C_T
6:     1 983004      T      G -0.0696 0.0210 0.0009299 rs13303118 1_983004_T_G

可能要把data_1数据拆成两份，分别merge，再merge，这里到这里为止，先去搞分群，以后在NSCC上搞

##/home/liang_y/project/dor/gwas/col_12345

##传送到服务器上

https://www.ebi.ac.uk/gwas/efotraits/EFO_0005140
/home/users/nus/bliu3/scratch/ly/utr2/single_clip_150/exon_cov/1/GWAS



##/home/liang_y/project/dor


mkdir IBD

mkdir CD
##/home/liang_y/project/dor/IBD/CD

GCST90446792 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90446001-GCST90447000/GCST90446792





mkdir UD


ulcerative colitis

GCST90446794 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90446001-GCST90447000/GCST90446794

mkdir RA

GCST90132223 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132223

##/home/liang_y/project/dor/RA
zcat GCST90132223.h.tsv.gz|head -n 2|awk '{print $8,$9,$10}'
zcat GCST90132223.h.tsv.gz|awk -v FS='\t' -v OFS='\t' '{if(NR>1) print $1,$2,$3,$4,$5,$6,$8,$9,$10,$11}' > t1 && wc -l t1  
awk -v OFS='_' '{print $1,$2,$3,$4}' t1 > col_1

paste t1 col_1|awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$11}' > t3
sed '1i chr\tpos\teffect\tref\tbeta\tSE\tP\t\SNP\tSNP_POS' t3 > t4  ##SNP_POS第一个SNP是effect


mkdir multiple sclerosis

GCST005531 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005531  ##不行
GCST001198 https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST001001-GCST002000/GCST001198/harmonised/ ##不行
GCST003566 https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003566/harmonised/


mkdir systemic lupus erythematosus

GCST007400 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007400
GCST011096 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST011001-GCST012000/GCST011096/harmonised/

zcat GCST011096.h.tsv.gz|awk -v FS='\t' -v OFS='\t' '{if(NR>1) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > t1 && wc -l t1  
##这里$9和$10不相等，检查后发现$10正确
##这里$7全部是NA

zcat GCST011096.h.tsv.gz|awk -v FS='\t' -v OFS='\t' '{if(NR>1) print $1,$2,$3,$4,$5,$6,$8,$10}' > t1 && wc -l t1  
awk -v OFS='_' '{print $1,$2,$3,$4}' t1 > t2
paste -d '\t' t1 t2 > t3
sed '1i chr\tpos\teffect\tref\tbeta\tSE\tP\t\SNP\tSNP_POS' t3 > t4



# # 读取以制表符分隔的数据文件
library(data.table)
data_1 <- fread("/home/liang_y/project/dor/gwas/col_123456_ok", header = TRUE, sep = "\t")
data_2 <- fread("t4", header = TRUE, sep = "\t")

> head(data_1)
     SNP_MAF_1   SNP_MAF_2   SMP_1   SMP_2
        <char>      <char>   <num>   <num>
1: 1_16103_G_T 1_16103_T_G 0.03877 0.96123
2: 1_51479_A_T 1_51479_T_A 0.19180 0.80820
3: 1_51898_A_C 1_51898_C_A 0.14210 0.85790
4: 1_51928_A_G 1_51928_G_A 0.14610 0.85390
5: 1_51954_C_G 1_51954_G_C 0.00000 1.00000
6: 1_54490_A_G 1_54490_G_A 0.17500 0.82500

data_1_1 <- data_1[, .(SNP_MAF_1, SMP_1)]
data_1_2 <- data_1[, .(SNP_MAF_2, SMP_2)]

> head(data_2)
     chr    pos effect    ref    beta     SE         P        SNP      SNP_POS
   <int>  <int> <char> <char>   <num>  <num>     <num>     <char>       <char>
1:     1 911428      T      C  0.0310 0.0278 0.2651000  rs4475691 1_911428_T_C
2:     1 911484      C      G  0.0322 0.0282 0.2530000   rs950122 1_911484_C_G
3:     1 978509      A      G -0.0643 0.0213 0.0025360  rs2340596 1_978509_A_G
4:     1 979472      C      G -0.0638 0.0212 0.0026700 rs13303368 1_979472_C_G
5:     1 979560      C      T -0.0669 0.0212 0.0015970 rs13303033 1_979560_C_T
6:     1 983004      T      G -0.0696 0.0210 0.0009299 rs13303118 1_983004_T_G

merge 这两个数据框，



merged_data <- merge(data_2, data_1, by = "SNP_MAF", all.x = TRUE)
any(is.na(merged_data))
colSums(is.na(merged_data))
# 筛选出 SMP 列中为 NA 的行
na_smp_rows <- merged_data[is.na(merged_data$SMP), ]




ankylosing spondylitis
GCST005529 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005529

Vitiligo
GCST004785 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004785


Psoriatic arthritis
GCST007043 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007043

Juvenile idiopathic arthritis
GCST005528  http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005528


microscopic colitis
GCST009081 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009081


Autoimmune thyroid diseases (Graves disease or Hashimoto's thyroiditis)
GCST005524 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005524


Sjogren syndrome
GCST012796 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST012001-GCST013000/GCST012796

Neuromyelitis optica
GCST005964 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005964


type 1 diabetes mellitus
GCST90013445 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90013001-GCST90014000/GCST90013445



Hb, hemoglobin
http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002384

Hematocrit
GCST90002383 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002383

platelet count
GCST90002402 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002402

Monocyte count
GCST90002393 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002393

leukocyte count
GCST90002407  http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002407

lymphocyte count
GCST90002388 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002388

neutrophil count
GCST90002398 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002398

Eosinophil counts
GCST90002381 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002381

mean corpuscular volume

GCST90002397 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002397


mean corpuscular hemoglobin concentration

GCST90002322 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002322

Mean corpuscular hemoglobin
GCST90002390 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002390


BMI
GCST90435413 http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90435001-GCST90436000/GCST90435413




