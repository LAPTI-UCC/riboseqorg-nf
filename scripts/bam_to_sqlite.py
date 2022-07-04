import sys
import pysam
import operator
import os
import time
import sqlite3
from sqlitedict import SqliteDict

def tran_to_genome(tran, pos, transcriptome_info_dict):
	#print ("tran",list(transcriptome_info_dict))
	traninfo = transcriptome_info_dict[tran]
	chrom = traninfo["chrom"]
	strand = traninfo["strand"]
	exons = sorted(traninfo["exons"])
	#print exons
	if strand == "+":
		exon_start = 0
		for tup in exons:
			exon_start = tup[0]
			exonlen = tup[1] - tup[0]
			if pos > exonlen:
				pos = (pos - exonlen)-1
			else:
				break
		genomic_pos = (exon_start+pos)-1
	elif strand == "-":
		exon_start = 0
		for tup in exons[::-1]:
			exon_start = tup[1]
			exonlen = tup[1] - tup[0]
			if pos > exonlen:
				pos = (pos - exonlen)-1
			else:
				break
		genomic_pos = (exon_start-pos)+1
	return (chrom, genomic_pos)


#  Takes a dictionary with a readname as key and a list of lists as value, each sub list has consists of two elements a transcript and the position the read aligns to in the transcript
#  This function will count the number of genes that the transcripts correspond to and if less than or equal to 3 will add the relevant value to transcript_counts_dict
def processor(process_chunk, master_read_dict, transcriptome_info_dict,master_dict,readseq, unambig_read_length_dict):
	readlen = len(readseq)
	ambiguously_mapped_reads = 0
	#get the read name
	read = list(process_chunk)[0]

	read_list = process_chunk[read] # a list of lists of all transcripts the read aligns to and the positions
	#used to store different genomic poistions
	genomic_positions = []

	#This section is just to get the different genomic positions the read aligns to

	for listname in process_chunk[read]:

		tran = listname[0].replace("-","_").replace("(","").replace(")","")

		pos = int(listname[1])
		genomic_pos = tran_to_genome(tran, pos, transcriptome_info_dict)
		#print ("genomic pos",genomic_pos)
		if genomic_pos not in genomic_positions:
			genomic_positions.append(genomic_pos)

	#If the read maps unambiguously
	if len(genomic_positions) == 1:
		if readlen not in unambig_read_length_dict:
			unambig_read_length_dict[readlen] = 0
		unambig_read_length_dict[readlen] += 1
		#assume this read aligns to a noncoding position, if we find that it does align to a coding region change this to True
		coding=False

		# For each transcript this read alings to
		for listname in process_chunk[read]:
			#get the transcript name
			tran = listname[0].replace("-","_").replace("(","").replace(")","")
			#If we haven't come across this transcript already then add to master_read_dict
			if tran not in master_read_dict:
				master_read_dict[tran] = {"ambig":{}, "unambig":{}, "mismatches":{}, "seq":{}}
			#get the raw unedited positon, and read tags
			pos = int(listname[1])
			read_tags = listname[2]
			#If there is mismatches in this line, then modify the postion and readlen (if mismatches at start or end) and add mismatches to dictionary
			nm_tag = 0
		
			for tag in read_tags:
				if tag[0] == "NM":
					nm_tag = int(tag[1])
			if nm_tag > 0:
				md_tag = ""
				for tag in read_tags:
					if tag[0] == "MD":
						md_tag = tag[1]
				pos_modifier, readlen_modifier,mismatches =  get_mismatch_pos(md_tag,pos,readlen,master_read_dict,tran,readseq)
				# Count the mismatches (we only do this for unambiguous)
				for mismatch in mismatches:
					#Ignore mismatches appearing in the first position (due to non templated addition)
					if mismatch != 0:
						char = mismatches[mismatch]
						mismatch_pos = pos + mismatch
						if mismatch_pos not in master_read_dict[tran]["seq"]:
							master_read_dict[tran]["seq"][mismatch_pos] = {}
						if char not in master_read_dict[tran]["seq"][mismatch_pos]:
							master_read_dict[tran]["seq"][mismatch_pos][char] = 0
						master_read_dict[tran]["seq"][mismatch_pos][char] += 1
				# apply the modifiers
				#pos = pos+pos_modifier
				#readlen = readlen - readlen_modifier


			try:
				cds_start = transcriptome_info_dict[tran]["cds_start"]
				cds_stop = transcriptome_info_dict[tran]["cds_stop"]

				if pos >= cds_start and pos <= cds_stop:
					coding=True
			except:
				pass


			if readlen in master_read_dict[tran]["unambig"]:
				if pos in master_read_dict[tran]["unambig"][readlen]:
					master_read_dict[tran]["unambig"][readlen][pos] += 1
				else:
					master_read_dict[tran]["unambig"][readlen][pos] = 1
			else:
				master_read_dict[tran]["unambig"][readlen] = {pos:1}

		if coding == True:
			master_dict["unambiguous_coding_count"] += 1
		elif coding == False:
			master_dict["unambiguous_non_coding_count"] += 1

	else:
		ambiguously_mapped_reads += 1
		for listname in process_chunk[read]:
			tran = listname[0].replace("-","_").replace("(","").replace(")","")
			if tran not in master_read_dict:
				master_read_dict[tran] = {"ambig":{}, "unambig":{}, "mismatches":{}, "seq":{}}
			pos = int(listname[1])
			read_tags = listname[2]
			nm_tag = 0
			for tag in read_tags:
				if tag[0] == "NM":
					nm_tag = int(tag[1])
			if nm_tag > 0:
				md_tag = ""
				for tag in read_tags:
					if tag[0] == "MD":
						md_tag = tag[1]
					pos_modifier, readlen_modifier,mismatches =  get_mismatch_pos(md_tag,pos,readlen,master_read_dict,tran,readseq)
					# apply the modifiers
					#pos = pos+pos_modifier
					#readlen = readlen - readlen_modifier
				if readlen in master_read_dict[tran]["ambig"]:
					if pos in master_read_dict[tran]["ambig"][readlen]:
						master_read_dict[tran]["ambig"][readlen][pos] += 1
					else:
						master_read_dict[tran]["ambig"][readlen][pos] = 1
				else:
					master_read_dict[tran]["ambig"][readlen] = {pos:1}
	return ambiguously_mapped_reads

def get_mismatch_pos(md_tag,pos,readlen,master_read_dict,tran,readseq):
	nucs = ["A","T","G","C"]
	mismatches = {}
	total_so_far = 0
	prev_char = ""
	for char in md_tag:
		if char in nucs:
			if prev_char != "":
				total_so_far += int(prev_char)
				prev_char = ""
			mismatches[total_so_far+len(mismatches)] = (readseq[total_so_far+len(mismatches)])
		else:
			if char != "^" and char != "N":
				if prev_char == "":
					prev_char = char
				else:
					total_so_far += int(prev_char+char)
					prev_char = ""
	readlen_modifier = 0
	pos_modifier = 0
	five_ok = False
	three_ok = False
	while five_ok == False:
		for i in range(0,readlen):
			if i in mismatches:
				pos_modifier += 1
				readlen_modifier += 1
			else:
				five_ok = True
				break
		five_ok = True


	while three_ok == False:
		for i in range(readlen-1,0,-1):
			if i in mismatches:
				readlen_modifier += 1
			else:
				three_ok = True
				break
		three_ok = True


	return (pos_modifier, readlen_modifier, mismatches)



def process_bam(bam_filepath, transcriptome_info_dict_path,outputfile):
	desc = "NULL"
	start_time = time.time()
	study_dict ={}
	nuc_count_dict = {"mapped":{},"unmapped":{}}
	dinuc_count_dict = {}
	threeprime_nuc_count_dict = {"mapped":{},"unmapped":{}}
	read_length_dict = {}
	unambig_read_length_dict = {}
	unmapped_dict = {}
	master_dict = {"unambiguous_non_coding_count":0,"unambiguous_coding_count":0,"current_dir":os.getcwd()}

	transcriptome_info_dict = {}
	connection = sqlite3.connect(transcriptome_info_dict_path)
	cursor = connection.cursor()
	cursor.execute("SELECT transcript,cds_start,cds_stop,length,strand,chrom,tran_type from transcripts;")
	result = cursor.fetchall()
	for row in result:
		transcriptome_info_dict[str(row[0])] = {"cds_start":row[1],"cds_stop":row[2],"length":row[3],"strand":row[4],"chrom":row[5],"exons":[],"tran_type":row[6]}
	#print list(transcriptome_info_dict)[:10]
	
	cursor.execute("SELECT * from exons;")
	result = cursor.fetchall()
	for row in result:
		transcriptome_info_dict[str(row[0])]["exons"].append((row[1],row[2]))

	#it might be the case that there are no multimappers, so set this to 0 first to avoid an error, it will be overwritten later if there is multimappers
	multimappers = 0
	unmapped_reads = 0
	unambiguous_coding_count = 0
	unambiguous_non_coding_count = 0
	trip_periodicity_reads = 0

	final_offsets = {"fiveprime":{"offsets":{}, "read_scores":{}}, "threeprime":{"offsets":{}, "read_scores":{}}}
	master_read_dict = {}
	prev_seq = ""
	process_chunk = {"read_name":[["placeholder_tran","1","28"]]}
	mapped_reads = 0
	ambiguously_mapped_reads = 0
	master_trip_dict = {"fiveprime":{}, "threeprime":{}}
	master_offset_dict = {"fiveprime":{}, "threeprime":{}}
	master_metagene_stop_dict = {"fiveprime":{}, "threeprime":{}}

	infile = pysam.Samfile(bam_filepath, "rb")
	header = infile.header["HD"]
	unsorted = False
	if "SO" in header:
		if header["SO"] != "queryname":
			unsorted = True
	else:
		unsorted = True
	if unsorted == True:
		print ("ERROR: Bam file appears to be unsorted or not sorted by read name. To sort by read name use the command: samtools sort -n input.bam output.bam")
		print (header,bam_filepath)
		sys.exit()
	total_bam_lines = 0
	all_ref_ids = infile.references

	for read in infile.fetch(until_eof=True):
		total_bam_lines += 1
		if not read.is_unmapped:
			ref = read.reference_id
			tran =  (all_ref_ids[ref]).split(".")[0]
			mapped_reads += 1
			if mapped_reads%1000000 == 0:
				print ("{} reads parsed at {}".format(mapped_reads,(time.time()-start_time)))
			pos = read.reference_start
			readname = read.query_name
			read_tags = read.tags
			if readname == list(process_chunk)[0]:
				process_chunk[readname].append([tran,pos,read_tags])
			#if the current read is different from previous reads send 'process_chunk' to the 'processor' function, then start 'process_chunk' over using current read
			else:
				if list(process_chunk)[0] != "read_name":

					#At this point we work out readseq, we do this for multiple reasons, firstly so we don't count the sequence from a read multiple times, just because
					# it aligns multiple times and secondly we only call read.seq once (read.seq is computationally expensive)
					seq = read.seq
					readlen = len(seq)

					# Note if a read maps ambiguously it will still be counted toward the read length distribution (however it will only be counted once, not each time it maps)
					if readlen not in read_length_dict:
						read_length_dict[readlen] = 0
					read_length_dict[readlen] += 1

					if readlen not in nuc_count_dict["mapped"]:
						nuc_count_dict["mapped"][readlen] = {}
					if readlen not in threeprime_nuc_count_dict["mapped"]:
						threeprime_nuc_count_dict["mapped"][readlen] = {}
					if readlen not in dinuc_count_dict:
						dinuc_count_dict[readlen] = {"AA":0, "TA":0, "GA":0, "CA":0,
									"AT":0, "TT":0, "GT":0, "CT":0,
									"AG":0, "TG":0, "GG":0, "CG":0,
									"AC":0, "TC":0, "GC":0, "CC":0}

					for i in range(0,len(seq)):
						if i not in nuc_count_dict["mapped"][readlen]:
							nuc_count_dict["mapped"][readlen][i] = {"A":0, "T":0, "G":0, "C":0, "N":0}
						nuc_count_dict["mapped"][readlen][i][seq[i]] += 1

					for i in range(0,len(seq)):
						try:
							dinuc_count_dict[readlen][seq[i:i+2]] += 1
						except:
							pass

					for i in range(len(seq),0,-1):
						dist = i-len(seq)
						if dist not in threeprime_nuc_count_dict["mapped"][readlen]:
							threeprime_nuc_count_dict["mapped"][readlen][dist] = {"A":0, "T":0, "G":0, "C":0, "N":0}
						threeprime_nuc_count_dict["mapped"][readlen][dist][seq[dist]] += 1
					ambiguously_mapped_reads += processor(process_chunk, master_read_dict, transcriptome_info_dict,master_dict,prev_seq, unambig_read_length_dict)
				process_chunk = {readname:[[tran, pos, read_tags]]}
				prev_seq = read.seq
		else:
			unmapped_reads += 1

			# Add this unmapped read to unmapped_dict so we can see what the most frequent unmapped read is.
			seq = read.seq
			readlen = len(seq)
			if seq in unmapped_dict:
				unmapped_dict[seq] += 1
			else:
				unmapped_dict[seq] = 1

			# Populate the nuc_count_dict with this unmapped read
			if readlen not in nuc_count_dict["unmapped"]:
				nuc_count_dict["unmapped"][readlen] = {}
			for i in range(0,len(seq)):
				if i not in nuc_count_dict["unmapped"][readlen]:
					nuc_count_dict["unmapped"][readlen][i] = {"A":0, "T":0, "G":0, "C":0, "N":0}
				nuc_count_dict["unmapped"][readlen][i][seq[i]] += 1

			if readlen not in threeprime_nuc_count_dict["unmapped"]:
				threeprime_nuc_count_dict["unmapped"][readlen] = {}

			for i in range(len(seq),0,-1):
				dist = i-len(seq)
				if dist not in threeprime_nuc_count_dict["unmapped"][readlen]:
					threeprime_nuc_count_dict["unmapped"][readlen][dist] = {"A":0, "T":0, "G":0, "C":0, "N":0}
				threeprime_nuc_count_dict["unmapped"][readlen][dist][seq[dist]] += 1

	#add stats about mapped/unmapped reads to file dict which will be used for the final report
	master_dict["total_bam_lines"] = total_bam_lines
	master_dict["mapped_reads"] = mapped_reads
	master_dict["unmapped_reads"] = unmapped_reads
	master_read_dict["unmapped_reads"] = unmapped_reads
	master_dict["ambiguously_mapped_reads"] = ambiguously_mapped_reads
	
	if "read_name" in master_read_dict:
		del master_read_dict["read_name"]
	print ("BAM file processed")
	print ("Creating metagenes, triplet periodicity plots, etc.")
	for tran in master_read_dict:
		try:
			cds_start = transcriptome_info_dict[tran]["cds_start"]
			cds_stop = transcriptome_info_dict[tran]["cds_stop"]
		except:
			continue

		tranlen = transcriptome_info_dict[tran]["length"]
		#Use this to discard transcripts with no 5' leader or 3' trailer
		if cds_start > 1 and cds_stop < tranlen and transcriptome_info_dict[tran]["tran_type"] == 1:
			for primetype in ["fiveprime", "threeprime"]:
				# Create the triplet periodicity and metainfo plots based on both the 5' and 3' ends of reads
				for readlength in master_read_dict[tran]["unambig"]:
					#print "readlength", readlength
					# for each fiveprime postion for this readlength within this transcript
					for raw_pos in master_read_dict[tran]["unambig"][readlength]:
						#print "raw pos", raw_pos
						trip_periodicity_reads += 1
						if primetype == "fiveprime":
							# get the five prime postion minus the cds start postion
							real_pos = raw_pos-cds_start
							rel_stop_pos = raw_pos-cds_stop
						elif primetype == "threeprime":
							real_pos = (raw_pos+readlength)-cds_start
							rel_stop_pos = (raw_pos+readlength)-cds_stop
						#get the readcount at the raw postion
						readcount = master_read_dict[tran]["unambig"][readlength][raw_pos]
						#print "readcount", readcount
						frame = (real_pos%3)
						if real_pos >= cds_start and real_pos <= cds_stop:
							if readlength in master_trip_dict[primetype]:
								master_trip_dict[primetype][readlength][str(frame)] += readcount
							else:
								master_trip_dict[primetype][readlength]= {"0":0.0,"1":0.0,"2":0.0}
								master_trip_dict[primetype][readlength][str(frame)] += readcount
						# now populate offset dict with the 'real_positions' upstream of cds_start, these will be used for metainfo dict
						if real_pos > (-600) and real_pos < (601):
							if readlength in master_offset_dict[primetype]:
								if real_pos in master_offset_dict[primetype][readlength]:
									#print "real pos in offset dict"
									master_offset_dict[primetype][readlength][real_pos] += readcount
								else:
									#print "real pos not in offset dict"
									master_offset_dict[primetype][readlength][real_pos] = readcount
							else:
								#initiliase with zero to avoid missing neighbours below
								#print "initialising with zeros"
								master_offset_dict[primetype][readlength]= {}
								for i in range(-600,601):
									master_offset_dict[primetype][readlength][i] = 0
								master_offset_dict[primetype][readlength][real_pos] += readcount

						# now populate offset dict with the 'real_positions' upstream of cds_start, these will be used for metainfo dict
						if rel_stop_pos > (-600) and rel_stop_pos < (601):
							if readlength in master_metagene_stop_dict[primetype]:
								if rel_stop_pos in master_metagene_stop_dict[primetype][readlength]:
									master_metagene_stop_dict[primetype][readlength][rel_stop_pos] += readcount
								else:
									master_metagene_stop_dict[primetype][readlength][rel_stop_pos] = readcount
							else:
								#initiliase with zero to avoid missing neighbours below
								master_metagene_stop_dict[primetype][readlength] = {}
								for i in range(-600,601):
									master_metagene_stop_dict[primetype][readlength][i] = 0
								master_metagene_stop_dict[primetype][readlength][rel_stop_pos] += readcount
								
	# master trip dict is now made up of readlengths with 3 frames and a count associated with each frame
	# create a 'score' for each readlength by putting the max frame count over the second highest frame count
	for primetype in ["fiveprime", "threeprime"]:
		for subreadlength in master_trip_dict[primetype]:
			maxcount = 0
			secondmaxcount = 0
			for frame in master_trip_dict[primetype][subreadlength]:
				if master_trip_dict[primetype][subreadlength][frame] > maxcount:
					maxcount = master_trip_dict[primetype][subreadlength][frame]
			for frame in master_trip_dict[primetype][subreadlength]:
				if master_trip_dict[primetype][subreadlength][frame] > secondmaxcount and master_trip_dict[primetype][subreadlength][frame] != maxcount:
					secondmaxcount = master_trip_dict[primetype][subreadlength][frame]
			# a perfect score would be 0 meaning there is only a single peak, the worst score would be 1 meaning two highest peaks are the same height
			master_trip_dict[primetype][subreadlength]["score"] = float(secondmaxcount)/float(maxcount)
	#This part is to determine what offsets to give each read length
	print ("Calculating offsets")
	for primetype in ["fiveprime", "threeprime"]:
		for readlen in master_offset_dict[primetype]:
			accepted_len = False
			max_relative_pos = 0
			max_relative_count = 0
			for relative_pos in master_offset_dict[primetype][readlen]:
				# This line is to ensure we don't choose an offset greater than the readlength (in cases of a large peak far up/downstream)
				if abs(relative_pos) < 10 or abs(relative_pos) > (readlen-10):
					continue
				if master_offset_dict[primetype][readlen][relative_pos] > max_relative_count:
					max_relative_pos = relative_pos
					max_relative_count = master_offset_dict[primetype][readlen][relative_pos]
			#print "for readlen {} the max_relative pos is {}".format(readlen, max_relative_pos)
			if primetype == "fiveprime":
				# -3 to get from p-site to a-site, +1 to account for 1 based co-ordinates, resulting in -2 overall
				final_offsets[primetype]["offsets"][readlen] = abs(max_relative_pos-2)
			elif primetype == "threeprime":
				# +3 to get from p-site to a-site, -1 to account for 1 based co-ordinates, resulting in +2 overall
				final_offsets[primetype]["offsets"][readlen] = (max_relative_pos*(-1))+2
			#If there are no reads in CDS regions for a specific length, it may not be present in master_trip_dict
			if readlen in  master_trip_dict[primetype]: 
				final_offsets[primetype]["read_scores"][readlen] = master_trip_dict[primetype][readlen]["score"]
			else:
				 final_offsets[primetype]["read_scores"][readlen] = 0.0
	master_read_dict["offsets"] = final_offsets
	master_read_dict["trip_periodicity"] = master_trip_dict
	master_read_dict["desc"] = "Null"
	master_read_dict["mapped_reads"] = mapped_reads
	master_read_dict["nuc_counts"] = nuc_count_dict
	master_read_dict["dinuc_counts"] = dinuc_count_dict
	master_read_dict["threeprime_nuc_counts"] = threeprime_nuc_count_dict
	master_read_dict["metagene_counts"] = master_offset_dict
	master_read_dict["stop_metagene_counts"] = master_metagene_stop_dict
	master_read_dict["read_lengths"] = read_length_dict
	master_read_dict["unambig_read_lengths"] = unambig_read_length_dict
	master_read_dict["coding_counts"] = master_dict["unambiguous_coding_count"]
	master_read_dict["noncoding_counts"] = master_dict["unambiguous_non_coding_count"]
	master_read_dict["ambiguous_counts"] = master_dict["ambiguously_mapped_reads"]
	master_read_dict["frequent_unmapped_reads"] = (sorted(unmapped_dict.items(), key=operator.itemgetter(1)))[-2000:]
	master_read_dict["cutadapt_removed"] = 0
	master_read_dict["rrna_removed"] = 0
	#If no reads are removed by minus m there won't be an entry in the log file, so initiliase with 0 first and change if there is a line
	master_read_dict["removed_minus_m"] = 0
	master_dict["removed_minus_m"] = 0
	# We work out the total counts for 5', cds 3' for differential translation here, would be better to do thisn in processor but need the offsets
	master_read_dict["unambiguous_all_totals"] = {}
	master_read_dict["unambiguous_fiveprime_totals"] = {}
	master_read_dict["unambiguous_cds_totals"] = {}
	master_read_dict["unambiguous_threeprime_totals"] = {}

	master_read_dict["ambiguous_all_totals"] = {}
	master_read_dict["ambiguous_fiveprime_totals"] = {}
	master_read_dict["ambiguous_cds_totals"] = {}
	master_read_dict["ambiguous_threeprime_totals"] = {}
	print ("calculating transcript counts")
	for tran in master_read_dict:
		if tran in transcriptome_info_dict:
			five_total = 0
			cds_total = 0
			three_total = 0

			ambig_five_total = 0
			ambig_cds_total = 0
			ambig_three_total = 0

			cds_start = transcriptome_info_dict[tran]["cds_start"]
			cds_stop = transcriptome_info_dict[tran]["cds_stop"]
			for readlen in master_read_dict[tran]["unambig"]:
				if readlen in final_offsets["fiveprime"]["offsets"]:
					offset = final_offsets["fiveprime"]["offsets"][readlen]
				else:
					offset = 15
				for pos in master_read_dict[tran]["unambig"][readlen]:
					real_pos = pos+offset
					if real_pos <cds_start:
						five_total += master_read_dict[tran]["unambig"][readlen][pos]
					elif real_pos >=cds_start and real_pos <= cds_stop:
						cds_total += master_read_dict[tran]["unambig"][readlen][pos]
					elif real_pos > cds_stop:
						three_total += master_read_dict[tran]["unambig"][readlen][pos]
			master_read_dict["unambiguous_all_totals"][tran] = five_total+cds_total+three_total
			master_read_dict["unambiguous_fiveprime_totals"][tran] = five_total
			master_read_dict["unambiguous_cds_totals"][tran] = cds_total
			master_read_dict["unambiguous_threeprime_totals"][tran] = three_total

			for readlen in master_read_dict[tran]["ambig"]:
				if readlen in final_offsets["fiveprime"]["offsets"]:
					offset = final_offsets["fiveprime"]["offsets"][readlen]
				else:
					offset = 15
				for pos in master_read_dict[tran]["ambig"][readlen]:
					real_pos = pos+offset
					if real_pos <cds_start:
						ambig_five_total += master_read_dict[tran]["ambig"][readlen][pos]
					elif real_pos >=cds_start and real_pos <= cds_stop:
						ambig_cds_total += master_read_dict[tran]["ambig"][readlen][pos]
					elif real_pos > cds_stop:
						ambig_three_total += master_read_dict[tran]["ambig"][readlen][pos]

			master_read_dict["ambiguous_all_totals"][tran] = five_total+cds_total+three_total+ambig_five_total+ambig_cds_total+ambig_three_total
			master_read_dict["ambiguous_fiveprime_totals"][tran] = five_total+ambig_five_total
			master_read_dict["ambiguous_cds_totals"][tran] = cds_total+ambig_cds_total
			master_read_dict["ambiguous_threeprime_totals"][tran] = three_total+ambig_three_total
	print ("Writing out to sqlite file")
	sqlite_db = SqliteDict(outputfile,autocommit=False)
	for key in master_read_dict:
		sqlite_db[key] = master_read_dict[key]
	sqlite_db["description"] = desc
	sqlite_db.commit()
	sqlite_db.close()


if __name__ == "__main__":
	if len(sys.argv) <= 2:
		print ("Usage: python bam_to_sqlite.py <path_to_bam_file> <path_to_organism.sqlite> <file_description (optional)>")
		sys.exit()
	bam_filepath = sys.argv[1]
	annotation_sqlite_filepath = sys.argv[2]
	#try:
	#	desc = sys.argv[3]
	#except:
	#	desc = bam_filepath.split("/")[-1]
	outputfile = bam_filepath+"v2.sqlite"
	process_bam(bam_filepath,annotation_sqlite_filepath,outputfile)
