#!/usr/bin/env python
import sys
import os
from os.path import exists
import argparse
import multiprocessing
#import pysam
from Bio import SeqIO


def get_opt():
	group = argparse.ArgumentParser()
	group.add_argument('-r', '--ref', help='reference contig level assembly', required=True)
	group.add_argument('-g', '--alleletable', help='Allele.ctg.table', required=True)
	group.add_argument('-b', '--bam', help='bam file, default: prunning.bam', default='prunning.bam')
	group.add_argument('-d', '--workdir', help='work directory, default: wrk_dir', default='wrk_dir')
	group.add_argument('-t', '--thread', help='threads, default: 10', type=int, default=10)
	return group.parse_args()


def read_fasta(in_fa,chrn):
	idx_name=''
	if opts.thread==1:
		idx_name = "index.fasta.biopython"
	else:
		idx_name = "index.fasta.biopython"+"."+chrn

	fa_db = SeqIO.index_db(idx_name, in_fa, "fasta")
	return fa_db

def load_allele(allele_table):
	ctg_on_chr = {}
	chr_contain_ctg = {}
	with open(allele_table, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			chrn = data[0]
			if chrn.startswith('tig') or chrn.startswith('scaffold') or chrn.startswith('utg') or chrn.startswith('ctg'):
				continue
			print(f'{chrn}\n')
			for ctg in data[2:]:
				if ctg not in ctg_on_chr:
					ctg_on_chr[ctg] = {}
				if chrn not in ctg_on_chr[ctg]:
					ctg_on_chr[ctg][chrn] = 0
				ctg_on_chr[ctg][chrn] += 1
	for ctg in ctg_on_chr:
		max_chr = ""
		max_cnt = 0
		for chrn in ctg_on_chr[ctg]:
			if ctg_on_chr[ctg][chrn] > max_cnt:
				max_cnt = ctg_on_chr[ctg][chrn]
				max_chr = chrn
		ctg_on_chr[ctg] = max_chr
		if max_chr not in chr_contain_ctg:
			chr_contain_ctg[max_chr] = {}
		chr_contain_ctg[max_chr][ctg] = 1
	return ctg_on_chr, chr_contain_ctg


def split_files(chrn, allele_table, ref, bam_file, wrkdir):
	wrk_dir = os.path.join(wrkdir, chrn)
	if not os.path.exists(wrk_dir):
		os.mkdir(wrk_dir)
	
	print("\tDealing %s"%chrn)
	ctg_on_chr, chr_contain_ctg = load_allele(allele_table)
	fa_db = read_fasta(ref,chrn)

	sub_bam = os.path.join(wrk_dir, chrn+'.bam')
	sub_fa = os.path.join(wrk_dir, chrn+'.fa')
	sub_bed = os.path.join(wrk_dir, chrn+'.bed')
	sub_list = os.path.join(wrk_dir, chrn+'.list')
	sub_fai = os.path.join(wrk_dir, chrn+'.fa.fai')
	with open(sub_fa, 'w') as fout, open(sub_bed, 'w') as fbed, open(sub_list, 'w') as flist:
		for ctg in chr_contain_ctg[chrn]:
			fout.write(">%s\n%s\n"%(ctg, fa_db[ctg].seq))
			fbed.write(f'{ctg}\t1\t{len(fa_db[ctg].seq)}\n')
			flist.write(f'{ctg}\n')
#DMRP	command1="samtools faidx "+ sub_fa
#DMRP	print(f'{command1}\n')
#DMRP	os.system(command1)
#DMRP	command2="samtools view --regions-file " + sub_bed + " " + bam_file + "| awk 'FNR==NR { a[$NF]; next } $7 in a ||  $7 == \"=\"' " + sub_list + " /dev/stdin | samtools view --threads 4 -o " + sub_bam + " -bt " +sub_fai
#DMRP	print(f'{command2}\n')
#DMRP	os.system(command2)
#DMRP	command3="samtools flagstat " + sub_bam + "> "+ os.path.join(wrk_dir, chrn+'.bam.flagstats')
#DMRP	print(f'{command3}\n')
#DMRP	os.system(command3)
#DMRP	with pysam.AlignmentFile(bam_file, 'rb') as fin:
#DMRP		with pysam.AlignmentFile(sub_bam, 'wb', template=fin) as fout:
#DMRP			for ctg in chr_contain_ctg[chrn]:
#DMRP				for line in fin.fetch(contig=ctg):
#DMRP					if line.next_reference_name and line.next_reference_name in ctg_on_chr and ctg_on_chr[line.next_reference_name]==chrn:
#DMRP						fout.write(line)
#DMRPsystem("samtools view --regions-file $wrkd/$chrn/seq.fasta.bed $bam| awk 'FNR==NR { a[\$NF]; next } \$7 in a ||  \$7 == \"=\"' $wrkd/$chrn/ctg.list /dev/stdin > $wrkd/$chrn/sample.clean.sam");
#DMRP        system("samtools view --threads $CPUS -bt $wrkd/$chrn/seq.fasta.fai $wrkd/$chrn/sample.clean.sam > $wrkd/$chrn/sample.clean.bam");
#DMRP        system("samtools flagstat $wrkd/$chrn/sample.clean.bam >  $wrkd/$chrn/sample.clean.bam.flagstat");
	

def partition_gmap(ref, allele_table, bam, wrkdir, threads):
	if not os.path.exists(wrkdir):
		os.mkdir(wrkdir)
	
	print("Getting groups")
	ctg_on_chr, chrn_db = load_allele(allele_table)
#DMRP	with open(allele_table, 'r') as fin:
#DMRP		for line in fin:
#DMRP			chrn_db[line.strip().split()[0]] = 1

	bai = bam+'.bai'
	if not os.path.exists(bai):
		print("BAI file not found, starting index...")
		ret = os.system('samtools index %s'%bam)
		if ret==0:
			print("Index success")
		else:
			print("Fatal: bam file must be sorted")
			sys.exit(-1)

	print("Splitting files")
	if len(chrn_db) < threads:
		threads = len(chrn_db)
	pool = multiprocessing.Pool(processes=threads)
	for chrn in chrn_db:
		pool.apply_async(split_files, (chrn, allele_table, ref, bam, wrkdir,))
	pool.close()
	pool.join()
	print("Notice: If you got errors of \"Length mismatch\" during allhic extract, it is normal because we split bam with the same header, it will not effect the result")
	print("Finished")


if __name__ == '__main__':
	opts = get_opt()
	ref = opts.ref
	allele_table = opts.alleletable
	bam = opts.bam
	wrkdir = opts.workdir
	threads = opts.thread
	partition_gmap(ref, allele_table, bam, wrkdir, threads)

