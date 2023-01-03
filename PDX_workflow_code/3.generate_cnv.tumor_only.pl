#!usr/bin/perl
use strict;
use warnings;
my $cpu = 1;
my $outdir = "/project/gccri/CPRIT_PDX/hef_folder/9.re_sequencing";

my $ref = "/home/hef/Data/hg38/resources_broad_hg38_v0/Homo_sapiens_assembly38.fasta";
my $cnvkit_ini = "/home/hef/Tools/cnvkit_liding/config/config.gencode_grch38.mgi.ini";
my $cnvkit_parameter = "-1.1,-0.25,0.2,0.7";
my $cnvkit_ref_wes = "/project/gccri/CPRIT_PDX/hef_folder/5.CNV/cnvkit/tumor_only/ref/Reference.WES.cnn";
my $cnvkit_ref_wgs = "/project/gccri/CPRIT_PDX/hef_folder/5.CNV/cnvkit/tumor_only/ref/Reference.WGS.cnn";

my $DB_gnomad="/home/hef/Data/files_liding/mutect/af-only-gnomad.hg38.vcf.gz";
my $DB_interval="/home/hef/Data/files_liding/mutect/wgs_calling_regions.hg38.interval_list";
my $region_excluded="/home/hef/Data/hg38/hg38.centromere_telomere.bed";
my $my_pon_WES = "/project/gccri/CPRIT_PDX/hef_folder/3.DNA_mutation/1.somatic/2.without_Germline/sm_withoutG/PON/pon.WES.vcf.gz";
my $my_pon_WGS = "/project/gccri/CPRIT_PDX/hef_folder/3.DNA_mutation/1.somatic/2.without_Germline/sm_withoutG/PON/pon.WGS.vcf.gz";
my $add = "--genotype-germline-sites true --genotype-pon-sites true";

my $mapping_bias_wes = "/project/gccri/CPRIT_PDX/hef_folder/5.CNV/pureCN/WES_mapping_bias.rds";
my $mapping_bias_wgs = "/project/gccri/CPRIT_PDX/hef_folder/5.CNV/pureCN/WGS_mapping_bias.rds";
my $blacklist = "/project/gccri/CPRIT_PDX/hef_folder/5.CNV/pureCN/hg38_simpleRepeats.bed";

my %sex;
open INFO, "</home/hef/2.project/1.PDX/6.CNV/sample_info.txt" or die $!;
while(<INFO>){
        my @t = split /\t/;
        next if (/^ID/);
        my $sex = "?";
        if($t[1] =~/^[Ff]emale/){
                $sex = "F";
        }
        elsif($t[1] =~/^[Mm]ale/){
                $sex = "M";
        }
        $sex{$t[0]} = $sex;
}


my $pwd = `pwd`;
chomp $pwd;
`mkdir -p $pwd/script_sh`;
`mkdir -p $outdir/4.cnv`;
open IN,"<$ARGV[0]" or die $!;
while(<IN>){
        chomp;
        my ($tumor_id, $tumor_bam) = split /\t/;
	open OUT,">script_sh/$tumor_id.cnv.sh" or die $!;	
	my $fout = "$outdir/4.cnv";
	`mkdir -p $fout/cnvkit $fout/cnvkit/$tumor_id  $fout/pureCN $fout/pureCN/$tumor_id`;
	
	#cnvkit
	my $cnvkit_ref = $cnvkit_ref_wgs;
        if($tumor_id =~ /WES/){
                $cnvkit_ref = $cnvkit_ref_wes;
        }
	
	my $tumor_bai = "$tumor_bam.bai";
	print OUT "/bin/rm $fout/cnvkit/$tumor_id/$tumor_id.T.ba*;ln -s $tumor_bam $fout/cnvkit/$tumor_id/$tumor_id.T.bam; ln -s $tumor_bai $fout/cnvkit/$tumor_id/$tumor_id.T.bam.bai\n";
	print OUT "bash /home/hef/Tools/cnvkit_liding/src/cnvkit_wxs.tumorOnly.v2.sh -N $tumor_id -C $cnvkit_ini -p $cnvkit_parameter -B $fout/cnvkit/$tumor_id/$tumor_id.T.bam -R $cnvkit_ref -O $fout/cnvkit/$tumor_id\n";
	print OUT "cnvkit.py export seg $fout/cnvkit/$tumor_id/$tumor_id.T.cns --enumerate-chroms -o $fout/pureCN/$tumor_id/$tumor_id.seg\n";



	my $mapping_bias = $mapping_bias_wes;
	if ($tumor_id =~ /WES/){
		print OUT "gatk Mutect2 -R $ref -L $DB_interval -I $tumor_bam -tumor $tumor_id --germline-resource $DB_gnomad --panel-of-normals $my_pon_WES -XL $region_excluded -O $fout/pureCN/$tumor_id/$tumor_id.raw.vcf $add\n";
	}
	elsif ($tumor_id =~ /WGS/){
		print OUT "gatk Mutect2 -R $ref -L $DB_interval -I $tumor_bam -tumor $tumor_id --germline-resource $DB_gnomad --panel-of-normals $my_pon_WGS -XL $region_excluded -O $fout/pureCN/$tumor_id/$tumor_id.raw.vcf $add\n";
		$mapping_bias = $mapping_bias_wgs;
	}
	my @tt = split /_/, $tumor_id;
	print OUT "Rscript /home/hef/Tools/miniconda3/lib/R/library/PureCN/extdata/PureCN.R --out $fout/pureCN/$tumor_id --sampleid $tumor_id --tumor $fout/cnvkit/$tt[0]_$tt[1]_WGS/$tt[0]_$tt[1]_WGS.T.cnr --seg-file $fout/pureCN/$tt[0]_$tt[1]_WGS/$tt[0]_$tt[1]_WGS.seg --mapping-bias-file $mapping_bias --vcf $fout/pureCN/$tumor_id/$tumor_id.raw.vcf --snp-blacklist $blacklist --genome hg38 --fun-segmentation none --force --post-optimize --seed 123 --sex $sex{$tt[0]}\n";
}
