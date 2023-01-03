#!usr/bin/perl
use strict;
use warnings;
my $cpu = 1;
my $outdir = "/project/gccri/CPRIT_PDX/hef_folder/9.re_sequencing";

my $insert_size = 300;
my $ref = "/home/hef/Data/hg38/resources_broad_hg38_v0/Homo_sapiens_assembly38.fasta";
my $DB_pon="/home/hef/Data/files_liding/mutect/gatk4_mutect2_4136_pon.vcf.gz";
my $DB_gnomad="/home/hef/Data/files_liding/mutect/af-only-gnomad.hg38.vcf.gz";
my $DB_interval="/home/hef/Data/files_liding/mutect/wgs_calling_regions.hg38.interval_list";
my $region_excluded="/home/hef/Data/hg38/hg38.centromere_telomere.bed";

my $pwd = `pwd`;
chomp $pwd;
`mkdir -p $pwd/script_sh`;
`mkdir -p $outdir/3.snv $outdir/3.snv/tumor_normal`;
open IN,"<$ARGV[0]" or die $!;
while(<IN>){
        chomp;
        my ($tumor_id, $normal_id, $tumor_bam, $normal_bam) = split /\t/;
	open OUT,">script_sh/$tumor_id.snv.sh" or die $!;	
	my $fout = "$outdir/3.snv/tumor_normal/$tumor_id";
	`mkdir -p $fout $fout/pindel $fout/strelka $fout/varscan $fout/mutect`;	
	#strelka
	my $add = "--exome";
        if ($tumor_id =~ /WGS$/){
                $add = "";
        }
	#varscan
	open LIST,">$fout/varscan/varscan_bam.ls" or die $!;
        print LIST "$normal_bam\n$tumor_bam\n";
	#pindel
	open CONFIG,">$fout/pindel/pindel.config" or die $!;
        print CONFIG "$tumor_bam\t$insert_size\t$tumor_id\n$normal_bam\t$insert_size\t$normal_id\n";
	#pindel filter
	open CF, ">$fout/pindel/pindel.filter.config" or die $!;
        print CF "pindel.filter.variants_file = $fout/pindel/all.header\npindel.filter.pindel2vcf = /home/hef/Tools/miniconda3/bin/pindel2vcf\npindel.filter.REF = /home/hef/Data/hg38/resources_broad_hg38_v0/Homo_sapiens_assembly38.fasta\npindel.filter.heterozyg_min_var_allele_freq = 0.2\npindel.filter.homozyg_min_var_allele_freq = 0.8\npindel.filter.mode = somatic\npindel.filter.apply_filter = true\npindel.filter.somatic.min_coverages_t = 10\npindel.filter.somatic.min_coverages_n = 3\npindel.filter.somatic.min_var_allele_freq = 0.10\npindel.filter.somatic.require_balanced_reads = true\npindel.filter.somatic.remove_complex_indels = true\npindel.filter.somatic.max_num_homopolymer_repeat_units = 6\npindel.filter.date=20221225\n";

        print OUT "source /home/hef/Tools/miniconda3/etc/profile.d/conda.sh;conda activate py2\n";
	# strelka
        print OUT "/bin/rm -rf $fout/strelka/manta/runWorkflow.py;configManta.py --normalBam=$normal_bam --tumorBam=$tumor_bam --referenceFasta=$ref $add --runDir=$fout/strelka/manta \n";
        print OUT "$fout/strelka/manta/runWorkflow.py -j $cpu -g 8 >>$fout/strelka_run.log 2>&1\n";
        print OUT "/bin/rm -rf $fout/strelka/runWorkflow.py;configureStrelkaSomaticWorkflow.py --normalBam=$normal_bam --tumorBam=$tumor_bam --referenceFasta=$ref --callMemMb=2048 $add --runDir=$fout/strelka --indelCandidates=$fout/strelka/manta/results/variants/candidateSmallIndels.vcf.gz >$fout/strelka_run.log 2>&1\n";
       
        print OUT "$fout/strelka/runWorkflow.py -m local -j $cpu -g 8 >>$fout/strelka_run.log 2>&1\n";
        print OUT "less $fout/strelka/results/variants/somatic.snvs.vcf.gz|bcftools view -i 'FILTER=\"PASS\"'>$fout/strelka/strelka.snvs.vcf;less $fout/strelka/results/variants/somatic.indels.vcf.gz|bcftools view -i 'FILTER=\"PASS\"'>$fout/strelka/strelka.indels.vcf\n";

	print OUT "picard MergeVcfs I=$fout/strelka/results/variants/somatic.snvs.vcf.gz I=$fout/strelka/results/variants/somatic.indels.vcf.gz O=$fout/strelka/somatic.merge.vcf\n";
        print OUT "bcftools view -f PASS $fout/strelka/somatic.merge.vcf >  $fout/strelka/somatic.merge.filter.vcf\n";
        print OUT "vcf-genotype-annotator $fout/strelka/somatic.merge.filter.vcf TUMOR 0/1 -o $fout/$tumor_id.strelka.vcf\n";

	#varscan
	print OUT "samtools mpileup -q 1 -Q 13 -f $ref -b $fout/varscan/varscan_bam.ls|varscan somatic - $fout/varscan/varscan.out.som --mpileup 1 --output-vcf 1 --output-snp $fout/varscan/varscan_snp --output-indel $fout/varscan/varscan_indel > $fout/varscan_run.log 2>&1\n";
        print OUT "varscan processSomatic $fout/varscan/varscan_indel.vcf\nvarscan processSomatic $fout/varscan/varscan_snp.vcf\n";
        print OUT "varscan somaticFilter $fout/varscan/varscan_snp.Somatic.hc.vcf  --indel-file $fout/varscan/varscan_indel.Somatic.hc.vcf --output-file $fout/varscan/varscan.filtered.snv.vcf\nvarscan somaticFilter $fout/varscan/varscan_indel.Somatic.hc.vcf --output-file $fout/varscan/varscan.filtered.indel.vcf\n";

	print OUT "bgzip $fout/varscan/varscan.filtered.indel.vcf; tabix -p vcf $fout/varscan/varscan.filtered.indel.vcf.gz\n";
        print OUT "bgzip $fout/varscan/varscan.filtered.snv.vcf; tabix -p vcf $fout/varscan/varscan.filtered.snv.vcf.gz\n";
        print OUT "bcftools concat -a $fout/varscan/varscan.filtered.indel.vcf.gz $fout/varscan/varscan.filtered.snv.vcf.gz > $fout/varscan/varscan.filtered.merge.vcf\n";
        print OUT "perl /home/hef/Tools/fpfilter.pl --vcf-file  $fout/varscan/varscan.filtered.merge.vcf --bam-file $tumor_bam  --sample TUMOR --reference $ref --output $fout/varscan/$tumor_id.varscan.result.vcf\n";
        print OUT "bcftools view -f PASS $fout/varscan/$tumor_id.varscan.result.vcf > $fout/$tumor_id.varscan.vcf\n";

	#mutect
        print OUT "gatk Mutect2 -R $ref -L $DB_interval -I $tumor_bam -tumor $tumor_id -I $normal_bam -normal $normal_id --germline-resource $DB_gnomad --panel-of-normals $DB_pon -O $fout/mutect/mutect.raw.vcf > $fout/mutect_run.log 2>&1\n";
        print OUT "gatk FilterMutectCalls -R $ref -V $fout/mutect/mutect.raw.vcf -O $fout/mutect/mutect.filter.vcf\n";
        print OUT "gatk SelectVariants -R $ref -V $fout/mutect/mutect.filter.vcf -O $fout/mutect/mutect.raw.snp.vcf --select-type-to-include SNP --select-type-to-include MNP;gatk SelectVariants -R $ref -V $fout/mutect/mutect.filter.vcf -O $fout/mutect/mutect.raw.indel.vcf --select-type-to-include INDEL\n";
        print OUT "less $fout/mutect/mutect.raw.snp.vcf|bcftools view -i 'FILTER=\"PASS\"'>$fout/mutect/mutect.fil.snp.vcf\nless $fout/mutect/mutect.raw.indel.vcf|bcftools view -i 'FILTER=\"PASS\"'>$fout/mutect/mutect.fil.indel.vcf\n";
	
	print OUT "picard MergeVcfs I=$fout/mutect/mutect.fil.snp.vcf I=$fout/mutect/mutect.fil.indel.vcf O=$fout/$tumor_id.mutect.vcf\n";
	         
	#pindel
        print OUT "pindel -T $cpu -f $ref -i $fout/pindel/pindel.config --exclude $region_excluded -o $fout/pindel/pindel_raw >$fout/pindel_run.log 2>&1\n";
        print OUT "cat $fout/pindel/pindel*_D $fout/pindel/pindel*_SI |grep ChrID > $fout/pindel/all.header\n";
        print OUT "perl /home/hef/Tools/somaticwrapper-master/pindel_filter.v0.5.pl $fout/pindel/pindel.filter.config\n"; #Li idng's pipeline
        print OUT "vt normalize $fout/pindel/all.header.CvgVafStrand_pass.Homopolymer_pass.vcf -r $ref -o $fout/pindel/$tumor_id.indel.filter.output.norm.vcf -n\n";
        print OUT "perl /home/hef/2.project/1.PDX/4.DNAseq_mutation/1.somatic/1.paired/add_pindel_flag.pl $fout/pindel/$tumor_id.indel.filter.output.norm.vcf > $fout/pindel/$tumor_id.indel.filter.output.norm.flag.vcf\n";
	print OUT "conda deactivate\n";

	#For pindel
	print OUT "vt normalize $fout/pindel/all.header.CvgVafStrand_pass.Homopolymer_pass.vcf -r $ref -o $fout/pindel/$tumor_id.indel.filter.output.norm.vcf -n\n";
        print OUT "perl /home/hef/2.project/1.PDX/4.DNAseq_mutation/1.somatic/1.paired/add_pindel_flag.pl $fout/pindel/$tumor_id.indel.filter.output.norm.vcf > $fout/pindel/$tumor_id.indel.filter.output.norm.flag.vcf\n";
        print OUT "bcftools view -f PASS $fout/pindel/$tumor_id.indel.filter.output.norm.flag.vcf > $fout/$tumor_id.pindel.vcf\n";


}
