#!usr/bin/perl
use strict;
use warnings;
my $cpu = 4;
my $outdir = "/project02/gccri/hef/CRIPT/analysis/RNA";

my $star_h = "/home/hef/Data/hg38/star";
my $star_m = "/home/hef/Data/mm10/star";

my $ref_h = "/home/hef/Data/hg38/hg38.fa";
my $ref_m = "/home/hef/Data/mm10/mm10.fa";

my $gtf_h = "/home/hef/Data/hg38/optimist_ref/gencodev22/gencode.v22.annotation.gtf";
my $gtf_m = "/home/hef/Data/mm10/gencode.vM19.annotation.gtf";
my $gff_h = "/home/hef/Data/hg38/optimist_ref/gencodev22/gencode.v22.annotation.gff3";

my $rnaseqc_gtf = "/home/hef/Data/hg38/gencode.v22.genes.gtf";
my $kallisto_h = "/home/hef/Data/hg38/optimist_ref/gencodev22/gencode.v22.all_transcripts.fa.idx";

my $ctat_lib_dir = "/home/hef/Data/hg38/optimist_ref/GRCh38_gencode_v22.star-fusion.v1.10";

my $pwd = `pwd`;
chomp $pwd;
`mkdir -p $outdir/0.disambiguate $outdir/1.trim $outdir/1.trim/trimmed_fq $outdir/1.trim/Final_fq $outdir/2.exp $outdir/2.exp/1.kallisto $outdir/3.fusion $outdir/3.fusion/1.star_fusion $outdir/3.fusion/2.PRADA $outdir/2.exp/2.RSEM $outdir/2.exp/3.htseq $outdir/4.rnaseqc`;
`mkdir -p $pwd/script_sh.rna`;

open OT, ">trimmed_fq.rna.ls" or die $!;
open IN,"<$ARGV[0]" or die $!;
while(<IN>){
        chomp;
        my ($id, $fq1, $fq2) = split /\t/;
        open SH, ">script_sh.rna/$id.rna.sh" or die $!;

	print SH "trim_galore --phred33 --fastqc --length 50 -q 20 --basename $id --gzip -o $outdir/1.trim/trimmed_fq --paired $fq1 $fq2\n";
        $fq1 = "$outdir/1.trim/trimmed_fq/$id*val_1.fq.gz";
        $fq2 = "$outdir/1.trim/trimmed_fq/$id*val_2.fq.gz";

	if ($id =~ /PDX/){
		#map to human
		print SH "STAR --runThreadN $cpu --genomeDir $star_h --sjdbGTFfile $gtf_h --sjdbOverhang 100 --readFilesIn $fq1 $fq2 --outFileNamePrefix $outdir/0.disambiguate/$id\.human. --outSAMtype BAM Unsorted --twopassMode Basic --outSAMattributes All --genomeLoad NoSharedMemory --readFilesCommand zcat --outReadsUnmapped Fastx --outSAMunmapped Within; samtools sort -m 3G -@ $cpu -o $outdir/0.disambiguate/$id\.human.sort.bam -n $outdir/0.disambiguate/$id\.human.Aligned.out.bam; /bin/rm -rf $outdir/0.disambiguate/$id\.human.Aligned.out.bam $outdir/0.disambiguate/$id\.human.*STAR* $outdir/0.disambiguate/$id\.human.Log* $outdir/0.disambiguate/$id\.human.*tab\n";
		#map to mouse
		print SH "STAR --runThreadN $cpu --genomeDir $star_m --sjdbGTFfile $gtf_m --sjdbOverhang 100 --readFilesIn $fq1 $fq2 --outFileNamePrefix $outdir/0.disambiguate/$id\.mouse. --outSAMtype BAM Unsorted --twopassMode Basic --outSAMattributes All --genomeLoad NoSharedMemory --readFilesCommand zcat --outReadsUnmapped Fastx --outSAMunmapped Within; samtools sort -m 3G -@ $cpu -o $outdir/0.disambiguate/$id\.mouse.sort.bam -n $outdir/0.disambiguate/$id\.mouse.Aligned.out.bam; /bin/rm -rf $outdir/0.disambiguate/$id\.mouse.Aligned.out.bam $outdir/0.disambiguate/$id\.mouse.*STAR* $outdir/0.disambiguate/$id\.mouse*out*\n";
		#convert to fq
		print SH "disambiguate -s $id -o $outdir/0.disambiguate -a star $outdir/0.disambiguate/$id\.human.sort.bam $outdir/0.disambiguate/$id\.mouse.sort.bam;samtools merge $outdir/0.disambiguate/$id\.disam.merge.bam $outdir/0.disambiguate/$id\.disambiguatedSpeciesA.bam $outdir/0.disambiguate/$id\.ambiguousSpeciesA.bam $outdir/0.disambiguate/$id\.ambiguousSpeciesB.bam; samtools sort -m 3G -@ $cpu -o $outdir/0.disambiguate/$id\.disam.sortbyname.bam -n $outdir/0.disambiguate/$id\.disam.merge.bam\n";
		print SH "samtools fastq $outdir/0.disambiguate/$id\.disam.sortbyname.bam -1 $outdir/0.disambiguate/$id\.R1.gz -2 $outdir/0.disambiguate/$id\.R2.gz -0 /dev/null -s /dev/null -n -F 0x900; /bin/rm -rf $outdir/0.disambiguate/$id\.human*bam $outdir/0.disambiguate/$id\.mouse* $outdir/0.disambiguate/$id\.*bam $outdir/0.disambiguate/$id\.*log\n";
		#merge with unmapped reads
		print SH "gunzip -c $outdir/0.disambiguate/$id\.R1.gz > $outdir/1.trim/Final_fq/$id\.R1.fq;gunzip -c $outdir/0.disambiguate/$id\.R2.gz > $outdir/1.trim/Final_fq/$id\.R2.fq\n";
		print SH "cat $outdir/1.trim/Final_fq/$id\.R1.fq $outdir/0.disambiguate/$id\.human.Unmapped.out.mate1 > $outdir/1.trim/Final_fq/$id\.R1.temp; mv $outdir/1.trim/Final_fq/$id\.R1.temp $outdir/1.trim/Final_fq/$id\.R1.fq; cat $outdir/1.trim/Final_fq/$id\.R2.fq $outdir/0.disambiguate/$id\.human.Unmapped.out.mate2 > $outdir/1.trim/Final_fq/$id\.R2.temp; mv $outdir/1.trim/Final_fq/$id\.R2.temp $outdir/1.trim/Final_fq/$id\.R2.fq;/bin/rm -rf $outdir/0.disambiguate/$id\.R1.gz $outdir/0.disambiguate/$id\.R2.gz\n";
	}
	else{
		print SH "gunzip -c $fq1 > $outdir/1.trim/Final_fq/$id\.R1.fq; gunzip -c $fq2> $outdir/1.trim/Final_fq/$id\.R2.fq\n";
	}

	$fq1 = "$outdir/1.trim/Final_fq/$id\.R1.fq";
        $fq2 = "$outdir/1.trim/Final_fq/$id\.R2.fq";
	print OT "$id\t$fq1\t$fq2\n";

	#kallisto	
	print SH "kallisto quant -i $kallisto_h -o $outdir/2.exp/1.kallisto/$id -t $cpu --plaintext $fq1 $fq2\n";

	#STAR-Fusion
	print SH "/home/hef/Tools/STAR-Fusion-v1.10.0/STAR-Fusion --genome_lib_dir $ctat_lib_dir --left_fq $fq1 --right_fq $fq2 --FusionInspector validate --examine_coding_effect --CPU $cpu --output_dir $outdir/3.fusion/1.star_fusion/$id;/bin/rm -rf $outdir/3.fusion/1.star_fusion/$id/*preliminary*\n";
	#PRADA
	print SH "source /home/hef/Tools/miniconda3/etc/profile.d/conda.sh;conda activate prada\n";
	print SH "python /home/hef/Tools/PRADA2-master/prada2.py --read1 $fq1 --read2 $fq2 --outdir $outdir/3.fusion/2.PRADA\n";
	print SH "python /home/hef/Tools/PRADA2-master/prada2.py --read1 $fq1 --read2 $fq2 --outdir $outdir/3.fusion/2.PRADA --fusion\n";
	print SH "conda deactivate\n";
	print SH "source /home/hef/Tools/miniconda3/etc/profile.d/conda.sh;conda activate py2\n";
	print SH "python /home/hef/Tools/PRADA2-master/prada2.py --read1 $fq1 --read2 $fq2 --outdir $outdir/3.fusion/2.PRADA --rsem\n";
	print SH "conda deactivate\n";
	print SH "ln -s $outdir/3.fusion/2.PRADA/$id/rsem_results/rsem.genes.results $outdir/2.exp/2.RSEM/$id.genes.results ; ln -s $outdir/3.fusion/2.PRADA/$id/rsem_results/rsem.isoforms.results $outdir/2.exp/2.RSEM/$id.isoforms.results;/bin/rm -rf $outdir/3.fusion/2.PRADA/$id/rsem_results/*bam $outdir/3.fusion/2.PRADA/$id/bam_results\n";
	#HTseq
	print SH "samtools sort -m 2G -@ $cpu -o $outdir/2.exp/3.htseq/$id.sorted.bam $outdir/3.fusion/1.star_fusion/$id/Aligned.out.bam; samtools index $outdir/2.exp/3.htseq/$id\.sorted.bam;/bin/rm -rf $outdir/3.fusion/1.star_fusion/$id/Aligned.out.bam\n";
        print SH "htseq-count -s no -f bam -r pos -n $cpu -t exon -i ID -m union --nonunique all $outdir/2.exp/3.htseq/$id\.sorted.bam $gff_h --additional-attr=gene_id --additional-attr=gene_name --additional-attr=transcript_id --additional-attr=exon_number > $outdir/2.exp/3.htseq/$id\.htseq.txt\n";
	#RNAseQC
	print SH "rnaseqc $rnaseqc_gtf $outdir/2.exp/3.htseq/$id\.sorted.bam --coverage $outdir/4.rnaseqc/$id\n";
	print SH "gzip $outdir/1.trim/Final_fq/$id\.R1.fq; gzip $outdir/1.trim/Final_fq/$id\.R2.fq\n"
}



