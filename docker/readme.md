# SpecImmune Docker Setup and Testing Guide

Place both `Dockerfile` and `.dockerignore` in the same directory.

After completing all make_db operations, move the hg38 files into the `/db` directory. The directory structure should be:
```
db/
├── CYP
├── hg38/
│   ├── hg38_no_alt.fa
│   ├── hg38_no_alt.fa.amb
│   ├── hg38_no_alt.fa.ann
│   ├── hg38_no_alt.fa.bwt
│   ├── hg38_no_alt.fa.fai
│   ├── hg38_no_alt.fa.pac
│   └── hg38_no_alt.fa.sa
├── HLA
├── HLA_CDS
├── IG_TR
├── KIR
```
Once the above steps are completed, build the Docker image:
```
docker build -t specimmune:v1.0.0 .
```
Docker Test Cases
Display the program help documentation:
```
docker run --rm specimmune:v1.0.1
```
Define directories:
```
dir_db=${specimmune_installed_dir}/db/
dir_test=${specimmune_installed_dir}/test
dir_db: The database directory that includes hg38
dir_test: The test directory containing test data
```
Example paths:
```
dir_db=/data10/zq123/SpecImmune/db/
dir_test=/data10/zq123/SpecImmune/test
```
Run Test Cases:

Test Case 1: HLA
```
docker run --rm -v $dir_db:/SpecImmune/db -v $dir_test:/SpecImmune/test \
specimmune:v1.0.0 --db /SpecImmune/db -r /SpecImmune/test/HLA/test_HLA_lite.fastq.gz \
-j 15 -i HLA -n test_HLA -o /SpecImmune/test/test_20250204 \
--align_method_1 minimap2 -y pacbio
```
Test Case 2: KIR
```
docker run --rm -v $dir_db:/SpecImmune/db -v $dir_test:/SpecImmune/test \
specimmune:v1.0.0 --db /SpecImmune/db/ -r /SpecImmune/test/KIR/KIR_dp50_acc98_1.fastq.gz \
-j 15 -i KIR -n test_new_KIR -o /SpecImmune/test/test_20250204 \
--hete_p 0.2 --align_method_1 minimap2 -y pacbio
```
Test Case 3: IG_TR
```
docker run --rm -v $dir_db:/SpecImmune/db -v $dir_test:/SpecImmune/test \
specimmune:v1.0.0 --db /SpecImmune/db/ -r /SpecImmune/test/IG_TR/vdj.fq.gz \
-j 15 -i IG_TR -n test_TR -o /SpecImmune/test/test_20250204 -y pacbio \
--hg38 /SpecImmune/db/hg38/hg38_no_alt.fa --align_method_1 minimap2
```
Test Case 4: CYP
```
docker run --rm -v $dir_db:/SpecImmune/db -v $dir_test:/SpecImmune/test \
specimmune:v1.0.0 --db /SpecImmune/db/ -r /SpecImmune/test/CYP/HG03579.CYP.fastq.0.1.fq.gz \
-j 15 -i CYP -n test_CYP_nano -o /SpecImmune/test/test_20250204 -y nanopore \
--hg38 /SpecImmune/db/hg38/hg38_no_alt.fa --align_method_1 minimap2
```
Run an Interactive Docker Session for Manual Testing
```
docker run -it --rm -v $dir_db:/SpecImmune/db -v $dir_test:/SpecImmune/test \
specimmune:v1.0.0 bash -c "source /opt/conda/bin/activate SpecImmune && \
python /SpecImmune/scripts/main.py --db /SpecImmune/db/ -r /SpecImmune/test/HLA/test_HLA_lite.fastq.gz \
-j 15 -i HLA -n test_HLA -o /SpecImmune/test/test_20250204 --align_method_1 minimap2 -y pacbio"
```
Now your Docker image is built, and you can run various test cases for SpecImmune!