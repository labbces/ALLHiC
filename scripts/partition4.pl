#!/usr/bin/perl -w

#This script requires a lot less RAM memory than the original version (b49ddea) and requires bioperl
use Getopt::Std;
getopts "g:d:b:r:";
use Bio::Index::Fasta;


if ((!defined $opt_g)|| (!defined $opt_r)) {
    die "************************************************************************
    Usage: perl $0 -g Allele.gene.table -r draft.asm.fasta
      -h : help and usage.
      -g : Allele.gene.table 
      -b : optional,default prunning.bam
      -r : reference ctg assembly
      -d : optional, default wrk_dir
************************************************************************\n";
}

my $bam    = (defined $opt_b)?$opt_b:"prunning.bam";
my $table  = $opt_g;
my $wrkd   = (defined $opt_d)?$opt_d:"wrk_dir";
my $refSeq = $opt_r;

### Read referece ctg fasta
print STDOUT "START: Indexing genome assembly -- ". localtime ."\n";

my $ctgn;
my $refSeq_inxFile='index.file';

my $refdb = Bio::Index::Fasta->new(-filename => $refSeq_inxFile,
                                   -write_flag => 1);
$refdb->make_index($refSeq);
print STDOUT "END: Indexing genome assembly -- ". localtime ."\n";

print STDOUT "START: Creating data structures to keep Allele.gene.table -- ". localtime ."\n";
### Assign ctgs to pre-defined clusters

my %ctgdb;
open(IN, $table) or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	my $chrn = $data[1];
	next if ($chrn !~ /^Chr/);
	foreach my $i(3..$#data){
		my $ctg = (split/,/,$data[$i])[1];
		$ctgdb{$ctg}->{$chrn}++;
		}
	}
close IN;
my %chrdb; ### pre-defined cluster based on chromosomes of close-releative species
foreach my $ctg (sort(keys(%ctgdb))){
	my $count = 0;
#DMRP	foreach my $chrn (sort {$ctgdb{$ctg}->{$b}<=>$ctgdb{$ctg}->{$a}} sort(keys(%{$ctgdb{$ctg}}))){
#DMRP		$count++;
#DMRP		next if($count>1);
#DMRP#		print "$ctg	$chrn	$ctgdb{$ctg}->{$chrn}\n";
#DMRP		$chrdb{$chrn} .= $ctg.",";
#DMRP		}
	my $maxChr = "";
	my $maxCnt = 0;
	foreach my $chrn (sort(keys(%{$ctgdb{$ctg}}))){
		if ($ctgdb{$ctg}->{$chrn} > $maxCnt){
			$maxCnt = $ctgdb{$ctg}->{$chrn};
			$maxChr = $chrn;
			}
		}
	$chrdb{$maxChr} .= $ctg.",";
	}

print "\tLOG:  There are ". scalar(keys(%ctgdb))." contigs assigned to chromossomes\n";
print "\tLOG:  There are ". scalar(keys(%chrdb))." chromossomes chromossomes with assigned contigs\n";
open TEST, ">test3.chr2ctg.tbl";
foreach my $c(sort(keys(%chrdb))){
 foreach my $ctg(sort(split(/,/,$chrdb{$c}))){
  print TEST "$c\t$ctg\n";
 }
}
print STDOUT "END: Creating data structures to keep Allele.gene.table -- ". localtime ."\n";

system("rm -rf $wrkd");
system("mkdir -p $wrkd");


print STDOUT "START: Separating assembled contigs and BAM hits per reference Chromosome -- ". localtime ."\n";

#Index input BAM file
if (-f $bam.".bai"){ $bam=$bam; print STDOUT "Index file for BAM $bam already present\n";}
else{
	print STDOUT "Index file for BAM $bam, does not exists, sorting and indexing -- ". localtime . "\n";
	my $sortedbam=$bam;
	$sortedbam=~s/.bam/.sorted.bam/;
	system("samtools sort -o $sortedbam $bam");
	system("samtools index $sortedbam");
	$bam=$sortedbam;
	print STDOUT "Finished creating index for BAM in $bam -- ". localtime . "\n";
	}

#Separate contigs and alignments for each reference chrosomome
my %fh_samout;
foreach my $chrn (sort(keys %chrdb)){
	next if ($chrn !~ /^Chr/);
	my @ctgdb  = split(/,/,$chrdb{$chrn});
	print "\tLOG:  There are ". scalar(@ctgdb)." contigs for chromosome $chrn\n";
	if (!-d "$wrkd/$chrn"){
		system("mkdir -p $wrkd/$chrn");
		}
	if (!-f "$wrkd/$chrn/ctg.list"){
#		print "$chrn @ctgdb\n";
		open(my $out, ">$wrkd/$chrn/ctg.list") or die"";
		map {print $out "$_\n";} @ctgdb;
		close $out;
		}
	if (!-f "$wrkd/$chrn/seq.fasta"){
		open(my $faout, ">$wrkd/$chrn/seq.fasta") or die"";
		map {chomp;print $faout ">$_\n".$refdb->get_Seq_by_id($_)->seq()."\n"} @ctgdb;
		close $faout;
		system("samtools faidx $wrkd/$chrn/seq.fasta");
		system("awk '{print \$1,1,\$2}' $wrkd/$chrn/seq.fasta.fai > $wrkd/$chrn/seq.fasta.bed");
		}
	my %tmpdb = (); $tmpdb{'='}++; ### need retain intra-contig links
 	map {$tmpdb{$_}++} @ctgdb;
	open(BAM, "samtools view --regions-file $wrkd/$chrn/seq.fasta.bed $bam|") or die"";
 	while(<BAM>){
         	chomp;
 		my ($c1,$c2) = (split/\s+/,$_)[2,6];
 		next if(!exists($tmpdb{$c1}) or !exists($tmpdb{$c2}));
 		$fh_samout{$chrn}=$chrn."FH";
 		open($fh_samout, ">>$wrkd/$chrn/sample.clean.sam") or die"";
 		print $fh_samout "$_\n";
 		}
 	close BAM;

	close $fh_samout{$chrn};
	system("samtools view -bt $wrkd/$chrn/seq.fasta.fai $wrkd/$chrn/sample.clean.sam > $wrkd/$chrn/sample.clean.bam");
	system("samtools flagstat $wrkd/$chrn/sample.clean.bam >  $wrkd/$chrn/sample.clean.bam.flagstat");
	}

print STDOUT "END: Separating assembled contigs per reference Chromosome -- ". localtime ."\n";
#print STDOUT "START: Separating BAM hits per reference Chromosome -- ". localtime ."\n";
### Read prunning BAM file
#my %saw_ctg;
#my %saw_chrn;
#
#
#foreach my $chrn (keys(%fh_samout)){
#	print $chrn."\n";
#	system("samtools view -bt $wrkd/$chrn/seq.fasta.fai $wrkd/$chrn/sample.clean.sam > $wrkd/$chrn/sample.clean.bam");
#	system("samtools flagstat $wrkd/$chrn/sample.clean.bam >  $wrkd/$chrn/sample.clean.bam.flagstat");
#	system("rm $wrkd/$chrn/sample.clean.sam")
#	}
#print STDOUT "END: Separating BAM hits per reference Chromosome -- ". localtime ."\n";
