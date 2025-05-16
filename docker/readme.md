1、dockerfile 还有 .dockerignore 两个文件，放到相同的目录下；
2、再执行完所有的 make_db 的操作后，把 hg38 的文件也放到 /db 目录下，目录结构如下：
db/
├──CYP
├──hg38/
│	├── hg38_no_alt.fa
│	├── hg38_no_alt.fa.amb
│	├── hg38_no_alt.fa.ann
│	├── hg38_no_alt.fa.bwt
│	├── hg38_no_alt.fa.fai
│	├── hg38_no_alt.fa.pac
│	└── hg38_no_alt.fa.sa
├──HLA
├──HLA_CDS
├──IG_TR
├──KIR

3、上述完成之后，执行 docker build -t specimmune:v1.0.3  .
4、docker 测试用例
#输出程序帮助文档；
docker run --rm specimmune:v1.0.3  

dir_db= ${specimmue_installed_dir}/db/
dir_test=${specimmue_installed_dir}/test


#$dir_db 是构建的db 目录， 包括hg38
#$dir_test 是测试的目录，里面有测试数据
dir_db=/data10/zq123/SpecImmune/db/
dir_test=/data10/zq123/SpecImmune/test

#输出程序帮助文档；
docker run --rm specimmune:v1.0.3  

#测试用例1：HLA
docker run --rm -v $dir_db:/SpecImmune/db -v $dir_test:/SpecImmune/test specimmune:v1.0.3  --db /SpecImmune/db -r /SpecImmune/test/HLA/test_HLA_lite.fastq.gz -j 15 -i HLA -n test_HLA -o /SpecImmune/test/test_20250206 --align_method_1 minimap2 -y pacbio

#测试用例2：KIR
docker run --rm -v $dir_db:/SpecImmune/db -v $dir_test:/SpecImmune/test specimmune:v1.0.3  --db /SpecImmune/db/ -r /SpecImmune/test/KIR/KIR_dp50_acc98_1.fastq.gz -j 15 -i KIR -n test_new_KIR  -o /SpecImmune/test/test_20250206  --hete_p 0.2 --align_method_1 minimap2 -y pacbio

#测试用例3：IG_TR
docker run --rm -v $dir_db:/SpecImmune/db -v $dir_test:/SpecImmune/test specimmune:v1.0.3 --db /SpecImmune/db/ -r  /SpecImmune/test/IG_TR/vdj.fq.gz   -j 15 -i IG_TR -n test_TR  -o /SpecImmune/test/test_20250206 -y pacbio  --hg38  /SpecImmune/db/hg38/hg38_no_alt.fa  --align_method_1 minimap2

#测试用例4：CYP
docker run --rm -v $dir_db:/SpecImmune/db -v $dir_test:/SpecImmune/test specimmune:v1.0.3 --db /SpecImmune/db/ -r  /SpecImmune/test/CYP/HG03579.CYP.fastq.0.1.fq.gz -j 15 -i CYP -n test_CYP_nano -o /SpecImmune/test/test_20250206 -y nanopore  --hg38  /SpecImmune/db/hg38/hg38_no_alt.fa  --align_method_1 minimap2




#docker run -it --rm -v $dir_db:/SpecImmune/db -v $dir_test:/SpecImmune/test specimmune:v1.0.3 bash -c "source /opt/conda/bin/activate SpecImmune && python /SpecImmune/scripts/main.py --db /SpecImmune/db/ -r /SpecImmune/test/HLA/test_HLA_lite.fastq.gz -j 15 -i HLA -n test_HLA -o /SpecImmune/test/test_20250204 --align_method_1 minimap2 -y pacbio"