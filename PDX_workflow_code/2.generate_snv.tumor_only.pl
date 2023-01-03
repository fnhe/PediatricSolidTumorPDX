#!usr/bin/perl
use strict;
use warnings;
my $cpu = 1;
my $outdir = "/project/gccri/CPRIT_PDX/hef_folder/9.re_sequencing";

my $insert_size = 300;
my $ref = "/home/hef/Data/hg38/resources_broad_hg38_v0/Homo_sapiens_assembly38.fasta";

my $pwd = `pwd`;
chomp $pwd;
`mkdir -p $pwd/script_sh`;
`mkdir -p $outdir/3.snv $outdir/3.snv/tumor_only`;
open IN,"<$ARGV[0]" or die $!;
while(<IN>){
        chomp;
        my ($id, $bam) = split /\t/;
	open OUT,">script_sh/$id.snv.sh" or die $!;	
	my $fout = "$outdir/3.snv/tumor_only/$id";
	`mkdir -p $fout $fout/strelka $fout/varscan $fout/mutect`;	
	#mutect
	print OUT "/bin/rm -rf $fout/mutect/$id\n";
	if ($id =~ /WES/){
                print OUT "sh /home/hef/Tools/somaticwrapper-master/somatic.Mutect2_tumorOnly/somaticMut.tumor-only.mutect2.sh -c /home/hef/2.project/1.PDX/4.DNAseq_mutation/1.somatic/2.without_G/config.mutect2.tumor_only.ini -p s0 -n $id -b $bam -o $fout/mutect\n";
        }
        elsif ($id =~ /WGS/){
                print OUT "sh /home/hef/Tools/somaticwrapper-master/somatic.Mutect2_tumorOnly/somaticMut.tumor-only.mutect2.sh -c /home/hef/2.project/1.PDX/4.DNAseq_mutation/1.somatic/2.without_G/config.mutect2.tumor_only.wgs.ini -p s0 -n $id -b $bam -o $fout/mutect\n";
        }
	print OUT "cp $fout/mutect/$id/filtered.rem_dbSNP_noCOSMIC.vcf $fout/$id.mutect.vcf\n";
	#strelka
	my $add = "--exome";
        if ($id =~ /WGS/){
                $add = "";
        }
        
	print OUT "source /home/hef/Tools/miniconda3/etc/profile.d/conda.sh;conda activate py2\n";
        print OUT "configureStrelkaGermlineWorkflow.py --bam $bam --referenceFasta $ref $add --runDir $fout/strelka;$fout/strelka/runWorkflow.py -m local -j $cpu >> $fout/strelka_run.log 2>&1\n";
        print OUT "bcftools view -f PASS $fout/strelka/results/variants/variants.vcf.gz > $fout/$id.strelka.vcf\n";
	
	#varscan
        open LIST, ">$fout/varscan/sample.list" or die $!;
        print LIST "TUMOR\n";
		
	print OUT "samtools mpileup -B -f $ref $bam > $fout/varscan/$id.pileup\n";
        print OUT "varscan mpileup2snp $fout/varscan/$id.pileup --output-vcf 1 --vcf-sample-list $fout/varscan/sample.list> $fout/varscan/$id.snp.vcf;varscan mpileup2indel $fout/varscan/$id.pileup --output-vcf 1 --vcf-sample-list $fout/varscan/sample.list > $fout/varscan/$id.indel.vcf\n";
        print OUT "varscan filter $fout/varscan/$id.indel.vcf --output-file $fout/varscan/$id.indel.fil.vcf;varscan filter $fout/varscan/$id.snp.vcf --indel-file $fout/varscan/$id.indel.fil.vcf --output-file $fout/varscan/$id.snp.fil.vcf\n";
        print OUT "bgzip $fout/varscan/$id.indel.fil.vcf; bgzip $fout/varscan/$id.snp.vcf; tabix -p vcf $fout/varscan/$id.snp.vcf.gz; tabix -p vcf $fout/varscan/$id.indel.fil.vcf.gz;bcftools concat -a $fout/varscan/$id.snp.vcf.gz $fout/varscan/$id.indel.fil.vcf.gz  > $fout/varscan/$id.filtered.merge.vcf\n";

        print OUT "perl /home/hef/Tools/fpfilter.pl --vcf-file  $fout/varscan/$id.filtered.merge.vcf --bam-file $bam  --sample TUMOR --reference $ref --output $fout/varscan/$id.varscan.result.vcf\n";
        print OUT "bcftools view -f PASS $fout/varscan/$id.varscan.result.vcf > $fout/$id.varscan.vcf\n";
        print OUT "conda deactivate\n";
}
