sample=$1
bam=$2
outdir=$3

dir=$(cd `dirname $0`; pwd)
bin=$dir/bin
db=$dir/db

echo 'chromosome	start	end	gene	log2' >$outdir/$sample.antitargetcoverage.cnn
cnvkit.py coverage $bam $db/KIR/gene_extend.bed --output $outdir/$sample.targetcoverage.cnn
cnvkit.py fix $outdir/$sample.targetcoverage.cnn $outdir/$sample.antitargetcoverage.cnn $db/KIR/ref/KIR.ref.cnn --output $outdir/$sample.fix.cnr --no-rmask --no-gc --no-edge
cnvkit.py segment $outdir/$sample.fix.cnr --output $outdir/$sample.fix.cns
cnvkit.py call $outdir/$sample.fix.cns -y -m threshold -t=-2,0.1,2 -o $outdir/$sample.called.cns

