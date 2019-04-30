#	CLOMP - Clinically Okay Metagenomic Pipeline
#	Copyright (C) 2019 Ryan C. Shean 
#		This program is free software: you can redistribute it and/or modify
#		it under the terms of the GNU General Public License as published by
#		the Free Software Foundation, either version 3 of the License, or
#		(at your option) any later version.

#		This program is distributed in the hope that it will be useful,
#		but WITHOUT ANY WARRANTY; without even the implied warranty of
#		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#		GNU General Public License for more details.

#		You should have received a copy of the GNU General Public License
#		along with this program. If not, see <http://www.gnu.org/licenses/>.

import ast 
import subprocess
import glob
import argparse 
import os
import operator
from collections import Counter
from ete3 import NCBITaxa
import timeit
ncbi = NCBITaxa()

# takes a path to the ini file and reads in all the variables, all variables
# are globals and assumed to exist and be valid 
# If the variable is assigned to <blank>, we are not sure what would happen.
# It would be good to have a validity checker for each of these variables.  #nicetohave
def parse_ini(ini_path):
	for line in open(ini_path):
		if line[0] and line[0] != '#' and '=' in line:
			line_list = line.split('=')
			var = line_list[0].strip()
			val = line_list[1].strip()
			#Of note, all of the paired-end functionality is untested.  Usually only R1 is used in taxonomical analysis.
			if var == 'PAIRED_END':
				global PAIRED_END
				PAIRED_END = ast.literal_eval(val)
			elif var == 'R1':
				global R1
				R1 = val
			elif var == 'R2':
				global R2
				R2 = val
			elif var == 'SAVE_MIDDLE_FILES':
				global SAVE_MIDDLE_FILES
				SAVE_MIDDLE_FILES = ast.literal_eval(val)
			elif var == 'THREADS': 
				global THREADS
				THREADS = val
			elif var == 'BASE_DELIMITER':
				global BASE_DELIMITER
				BASE_DELIMITER = val
			elif var == 'HOST_FILTER_FIRST':
				global HOST_FILTER_FIRST
				HOST_FILTER_FIRST = ast.literal_eval(val)
			elif var == 'INPUT_SUFFIX':
				global INPUT_SUFFIX
				INPUT_SUFFIX = val
			elif var == 'OVERWRITE':
				global OVERWRITE
				OVERWRITE = ast.literal_eval(val)
			elif var == 'TRIMMOMATIC_JAR_PATH':
				global TRIMMOMATIC_JAR_PATH
				TRIMMOMATIC_JAR_PATH = val
			elif var == 'TRIMMOMATIC_ADAPTER_PATH':
				global TRIMMOMATIC_ADAPTER_PATH
				TRIMMOMATIC_ADAPTER_PATH = val
			elif var == 'SEQUENCER':
				global SEQUENCER
				SEQUENCER = val
			elif var == 'TRIMMOMATIC_OPTIONS':
				global TRIMMOMATIC_OPTIONS
				TRIMMOMATIC_OPTIONS = val
			elif var == 'TRIM_SUFFIX':
				global TRIM_SUFFIX
				TRIM_SUFFIX = val
			elif var == 'BWT_DB_LOCATION':
				global BWT_DB_LOCATION
				BWT_DB_LOCATION = val
			elif var == 'BWT_OPTIONS':
				global BWT_OPTIONS
				BWT_OPTIONS = val
			elif var == 'BWT_SUFFIX':
				global BWT_SUFFIX
				BWT_SUFFIX = val
			elif var == 'BWT_CLEAN':
				global BWT_CLEAN
				BWT_CLEAN = ast.literal_eval(val)
			elif var == 'DB_LIST':
				global DB_LIST
				DB_LIST = ast.literal_eval(val)
			elif var == 'SNAP_ALIGNER_LOCTION':
				global SNAP_ALIGNER_LOCTION
				SNAP_ALIGNER_LOCTION = val			
			elif var == 'SNAP_OPTIONS':
				global SNAP_OPTIONS
				SNAP_OPTIONS = val
			elif var == 'KRAKEN_PATH':
				global KRAKEN_PATH
				KRAKEN_PATH = val
			elif var == 'KRAKEN_DB_PATH':
				global KRAKEN_DB_PATH
				KRAKEN_DB_PATH = val
			elif var == 'EDIT_DISTANCE_OFFSET':
				global EDIT_DISTANCE_OFFSET
				EDIT_DISTANCE_OFFSET = int(val)
			elif var == 'FILTER_LIST':
				global FILTER_LIST
				FILTER_LIST = ast.literal_eval(val)
			elif var == 'H_STRICT':
				global H_STRICT
				H_STRICT = ast.literal_eval(val)
			elif var == 'H_TAXID':
				global H_TAXID
				H_TAXID = val
			elif var == 'LOGIC':
				global LOGIC
				LOGIC = val
			elif var =='WRITE_UNIQUES':
				global WRITE_UNIQUES
				WRITE_UNIQUES = ast.literal_eval(val)
			elif var == 'BLAST_CHECK':
				global BLAST_CHECK
				BLAST_CHECK = ast.literal_eval(val)
			elif var == 'INCLUSION_TAXID':
				global INCLUSION_TAXID
				INCLUSION_TAXID = ast.literal_eval(val)
			elif var == 'EXCLUSION_TAXID':
				global EXCLUSION_TAXID
				EXCLUSION_TAXID = ast.literal_eval(val)
			elif var == 'BUILD_SAMS':
				global BUILD_SAMS
				BUILD_SAMS = ast.literal_eval(val)
				
			
# Takes a list of files to host filter and if output files don't exist, host filters and writes output
# Returns a set containing all output files 
def host_filter(to_host_filter_list):
	# Use a set to prevent duplicates when using paired end reads 
	done_host_filtering_list = set()
	
	# Go through every file in the file list that was passed to host_filter and take the unique sample name.
	for r1 in to_host_filter_list:
		base = r1.split(BASE_DELIMITER)[0]
		
		# Check if output file does not exist or if overwrite is turned on
		if not os.path.isfile(base + BASE_DELIMITER + R1 + BWT_SUFFIX) or OVERWRITE:
			if PAIRED_END:
				#Build name of the R2 file by taking sample name from R1 and putting in R2 where the R1 was.
				r2 = r1.split(R1)[0] + R2 + r1.split(R1)[1]
				
				#Create alignment command for paired-end reads for host filtering with bowtie2.
				align_cmd = 'bowtie2 ' + BWT_OPTIONS + ' --threads ' + THREADS + ' -x ' + \
					BWT_DB_LOCATION + ' -q -1 ' + r1 + + ' -2 ' + r2 + ' -S ' + base + \
					'_mappedSam ' + ' 2>&1 | tee -a ' + base + '.log'
			else:
				#If not paired-end, then create alignment command for single-end reads for host filtering with bowtie2.
				align_cmd = 'bowtie2 ' + BWT_OPTIONS + ' --threads ' + THREADS + ' -x ' + \
				BWT_DB_LOCATION + ' -q -U ' + r1 + ' -S ' + base + '_mappedSam' + ' 2>&1 | tee -a ' + \
				base + '.log'
			
			#Call the bowtie2 command from above.
			subprocess.call(align_cmd, shell=True)
			
			#Convert from sam to bam.
			sort_cmd = 'samtools view -Sb -@ ' + THREADS + ' ' + base + '_mappedSam > ' + base + '_mappedBam' 
			subprocess.call(sort_cmd, shell=True)
			#Unclear if this command is needed.
			sort_cmd2 = 'samtools sort -@ ' + THREADS + ' ' + base + '_mappedBam ' + base +  '_sorted'
			subprocess.call(sort_cmd2, shell=True)
			#Unclear if this command is needed.
			index_cmd = 'samtools index ' + base + '_sorted.bam'
			subprocess.call(index_cmd, shell=True)
			
			
			#Write the output FASTQ file that has been host filtered.
			r1_cmd = 'samtools view -@ ' + THREADS + ' -F 0x40 ' + base + '_mappedBam | ' \
				'awk \'{if($3 == \"*\") print \"@\" $1 \"\\n\" $10 \"\\n\" \"+\" $1 \"\\n\" $11}\' > ' + \
				base + BASE_DELIMITER + R1 + BWT_SUFFIX
			subprocess.call(r1_cmd,shell=True)
			
			#Add the output file to our "done" file list.
			done_host_filtering_list.add(base + R1 + BWT_SUFFIX)
			
			#Perform host filtering on R2 as needed and write output FASTQ that has been host filtered.
			if PAIRED_END:
				r2_cmd = 'samtools view -@ ' + THREADS + ' -f 0x40 ' + base + '_mappedBam | ' \
				'awk \'{if($3 == \"*\") print \"@\" $1 \"\\n\" $10 \"\\n\" \"+\" $1 \"\\n\" $11}\' > ' + \
				base + BASE_DELIMITER + R2 + BWT_SUFFIX
				subprocess.call(r2_cmd, shell=True)
			
			if BWT_CLEAN:
				subprocess.call('rm ' + base + '_mappedSam', shell=True)
				subprocess.call('rm ' + base + '_mappedBam', shell=True)
				subprocess.call('rm ' + base + '_sorted.bam', shell=True)
				subprocess.call('rm ' + base + '_sorted.bam.bai', shell=True)
			if not SAVE_MIDDLE_FILES:
				subprocess.call('rm ' + r1, shell=True)
				if PAIRED_END:
					subprocess.call('rm ' + r2, shell=True)
		else:
			done_host_filtering_list.add(base + BASE_DELIMITER + R1 + BWT_SUFFIX)
			
	return done_host_filtering_list
				
		
		
# takes a string matching the end of file names that we want to trim. Uses trimmomatic to clip
# adapter sequences and short/low quality reads. Writes output to TRIM_SUFFIX returns a set containing
# all the output files 
def trim(to_trim_list):
	done_trim_list = set()
	
	for file_name in to_trim_list:
		base = file_name.split(BASE_DELIMITER)[0]
		if not os.path.isfile(base + BASE_DELIMITER + R1 + TRIM_SUFFIX) or OVERWRITE:
			if PAIRED_END:
				r1 = file_name
				r2 = file_name.split(R1)[0] + R2 + file_name.split(R2)[1]
				trim_cmd = 'java -jar ' + TRIMMOMATIC_JAR_PATH + ' PE -threads ' + THREADS + ' ' + r1 + ' ' + \
					r2 + ' ' + base + BASE_DELIMITER + R1 + TRIM_SUFFIX + ' ' + base + BASE_DELIMITER + R2 + TRIM_SUFFIX + ' ' + SEQUENCER + \
					TRIMMOMATIC_ADAPTER_PATH + TRIMMOMATIC_OPTIONS + ' > ' + base + '.trim.log'
				
			else:
				trim_cmd = 'java -jar ' + TRIMMOMATIC_JAR_PATH + ' SE -threads ' + THREADS + ' ' + file_name + \
					' ' + base + BASE_DELIMITER + R1 + TRIM_SUFFIX +  ' ' + SEQUENCER + TRIMMOMATIC_ADAPTER_PATH + \
					TRIMMOMATIC_OPTIONS + ' > ' + base + '.trim.log'
				
			subprocess.call(trim_cmd,shell=True)
			done_trim_list.add(base + BASE_DELIMITER + R1 + TRIM_SUFFIX)
			if not SAVE_MIDDLE_FILES:
				if PAIRED_END:
					subprocess.call('rm ' + r1, shell=True)
					subprocess.call('rm ' + r2, shell=True)
				else:
					subprocess.call('rm ' + file_name)
		else:
			done_trim_list.add(base + BASE_DELIMITER + R1 + TRIM_SUFFIX)
			
	return done_trim_list

# takes a list of files and aligns all of them to all the snap databases found in DB_LIST
# writes output .sam files to disk and returns a set containing all outputs 
def snap(to_snap_list):
	done_snap_list = set()
	
	# Incremental count for each database that is passed in the initialization file.  This count is appended to the output file.
	count = -1 
	#Loop through the databases from the initialization file.  The goal is to perform each database load into RAM only once.
	for db in DB_LIST:
		#The added variable holds whether a given database has been aligned to for a set of reads.
		added = False 
		count +=  1
		
		# Build SNAP command.  It will be long and comma-separated in order to load each database only once
		snap_cmd = SNAP_ALIGNER_LOCTION + ' ' 
		# if you comma separate commands all linking the same database SNAP will load the database
		# once which DRASTICALLY speeds up the whole process
		for file_name in to_snap_list:
			base = file_name.split(BASE_DELIMITER)[0]
			if not os.path.isfile(base + BASE_DELIMITER + str(count) + '.sam') or OVERWRITE:
				added = True
				print('Aligning ' + base + ' to ' + db)
			
				if PAIRED_END:
					r1 = file_name
					r2 = file_name.split(R1)[0] + R2 + file_name.split(R2)[1] 
					#SNAP will perform the alignment of paired-end reads together, adding to comma-separated SNAP command.
					snap_cmd += ' paired ' + db + ' ' + r1 + ' ' + r2 + ' -o ' + base + BASE_DELIMITER + str(count) + \
						'.sam -t ' + THREADS + ' ' + SNAP_OPTIONS + ' , '
				else:
					snap_cmd += ' single ' + db + ' ' + file_name + ' -o ' + base + BASE_DELIMITER + str(count) + \
						'.sam -t ' + THREADS + ' ' + SNAP_OPTIONS + ' , '
			done_snap_list.add(base + BASE_DELIMITER + str(count) + '.sam')
		# only call snap if we've found at least one file without output
		if added:
			#Execute SNAP command built from above.
			subprocess.call(snap_cmd,shell=True)
		
	return done_snap_list  


# The tie-breaking function takes a list of taxid,edit_distances in the form [[taxid,edit_distance],[taxid,edit_distance],..
# references the global variable LOGIC to control the underlying tiebreaking code, breaks ties and 
# returns a single taxid as a string and a Boolean value (as to whether the read needs to re-BLASTEd against host).
def tie_break(taxid_list):
	score_list = [] 
	actual_taxid_list = []
	for id in taxid_list:
		score_list.append(id[1])
	
	# this can filter out any hits that are sufficiently worse than the edit distance for the 'best'
	# alignment - set to zero to only hold scores that have the best edit distance , and increase to a 
	# number greater than the snap option -d to hold all hits
	best_edit_distance = min(score_list) + EDIT_DISTANCE_OFFSET
	
	# Keep taxids that have an edit distance less than the acceptable edit distance defined above 
	for id in taxid_list:
		if id[1] <= best_edit_distance:
			actual_taxid_list.append(id[0])
	#No longer holding edit distances		
	taxid_list = actual_taxid_list
	lineage_list = []
	
	for id in taxid_list:
		# Not all taxids have valid lineages 
		try:
			#Not every taxid has a valid NCBI lineage, so this statement has to be encased in a try statement.
			lineage = ncbi.get_lineage(id)
			# filters out lineages that contain taxids in the FILTER_LIST variable
			# commonly, this is 'other sequences', 'artificial sequences' or 'environmental samples' 
			if any(x in FILTER_LIST for x in lineage):
				lineage = []
		except:
			lineage = []
		
		if lineage:
			lineage_list.append(lineage)
	
	if not lineage_list:
		return ['*',False]
	
	# controls if use any alignment to the human genome as grounds for classification as human source
	if H_STRICT:
		# check if H_TAXID ever shows up in 
		if any(int(H_TAXID) in sl for sl in lineage_list):
			return [H_TAXID,False]
			

	# count all taxids in all lineages 
	taxid_to_count_map = {}
	for each_lineage in lineage_list:
		for each_taxid in each_lineage:
			if each_taxid in taxid_to_count_map:
				taxid_to_count_map[each_taxid] += 1
			else:
				taxid_to_count_map[each_taxid] = 1
	
	#Set the threshold according to the pre-specified LOGIC in the initialization file
	num_assignments = len(lineage_list)
	if LOGIC == 'strict':
		threshold = num_assignments
	elif LOGIC == '90':
		threshold = num_assignments - ((num_assignments /10) + 1) 
	elif LOGIC == 'oneoff':
		threshold = num_assignments - 1
	else:
		print('invalid logic threshold: defaulting to strict')
		threshold = num_assignments
	
	#Now we will find all the taxids that meet threshold/LOGIC specified above.
	surviving_taxids = []
	for taxid_key in taxid_to_count_map:
		# main filtering - everything that passes this step gets a list intersection and the 
		# most specific taxid left is returned 
		if taxid_to_count_map[taxid_key] >= threshold:
			surviving_taxids.append(taxid_key)
			
	if len(surviving_taxids) == 0:
		return ['*',False]

	d = {}

	for the_value in surviving_taxids:
		d[the_value] = len(ncbi.get_lineage(the_value))
	#Find the remaining taxid with the longest lineage.
	#The longest lineage is defined as the longest list.  #nicetohave: this is not pulling the taxonomic rank at any point.
	assigned_hit = max(d.iteritems(), key=operator.itemgetter(1))[0]
	recheck = False 
	
	#Assign a Boolean value for each read as to whether it needs to be searched against a custom BLAST database
	#Here, we are just assigning whether it needs to be searched or not.  The custom BLAST database would need to be made separately.
	#All reads downstream of INCLUSION_TAXID and but not downstream of EXCLUSION_TAXID will be assigned a True value.
	if BLAST_CHECK:
		assigned_lineage = ncbi.get_lineage(assigned_hit)
		if INCLUSION_TAXID in assigned_lineage and EXCLUSION_TAXID not in assigned_lineage:
			recheck = True
	
	return [assigned_hit, recheck]

#Every read has a taxid assignment or is unassigned at this point.

# wrapper for a kraken script that converts tab seperated taxid\tcount file and writes a 
# Pavian output file for it. Requires a copy of ncbi's taxonomy database and some blank files
# This is a map of assigned taxids to number of occurrences as well as the total number of unassigned reads.
def new_write_kraken(basename, final_counts_map, num_unassigned):
	print('Preparing output for ' + basename)
	# we write a file in the form taxid\tcount 
	l = open(basename + '_temp_kraken.tsv', 'w')
	
	# initialize with the number of unassigned, we'll need to add human host filtering in earlier
	# because some reads will get tie broken to human 
	
	l.write('0\t' + str(num_unassigned))
	# write the rest of the taxids to the file
	for key in final_counts_map.keys():
		l.write('\n' + str(key) + '\t' + str(final_counts_map[key]))
	
	# this close is required so we get an EOF before passing the file to kraken-report 
	l.close()
	
	# kraken-report creates a file that Pavian likes - we name the file base_final_report.tsv
	kraken_report_cmd = KRAKEN_PATH + ' --db ' + KRAKEN_DB_PATH + ' --taxon-counts ' + basename + \
		'_temp_kraken.tsv > ' + basename + '_final_report.tsv'
	subprocess.call(kraken_report_cmd, shell=True)
	
if __name__ == '__main__':
	
	# CLOMP assumes the FASTQ sequence files for a given run are in their own folder along with the initialization file below.
	
	parse_ini('CLOMP.ini')
	
	#Assumes that the files coming off the sequencer will be gzipped, but okay if not.  Currently does not handle other compressions.
	subprocess.call('gunzip *.gz',shell=True)
	
	# We assume that if you are submitting paired end reads every R1 has an R2
	# Core data flow is through lists of files.
	# here we generate the input file list using the provided variables and wildcard expansion (glob is basically like unix ls)
	# Of note, R1 here is a variable that assigns how we identify your R1 reads.
	# We identify the list of the samples to be run off the R1 readfile name.
	if PAIRED_END:
		input_list = glob.glob('*' + R1 + '*' + INPUT_SUFFIX)
	else:
		R1 = ''
		input_list = glob.glob('*' + INPUT_SUFFIX)
	
	#Each of these functions takes a file list as input and returns a file list as output.
	
	if HOST_FILTER_FIRST:
		hf_list = host_filter(input_list)
		trim_list = trim(hf_list)
		to_snap_list = trim_list
	else:
		trim_list = trim(input_list)
		hf_list = host_filter(trim_list)
		to_snap_list = hf_list 
		
	sam_list = snap(to_snap_list)
		
	#file_list = glob.glob('*.sam')
	base_col = set()
	main_start_time = timeit.default_timer()
	
	#Create a list of all of the sample names for which there is a SAM file.
	for item in sam_list:
		base_col.add(item.split(BASE_DELIMITER)[0])
	#For every sample in the folder, go through every SAM file.
	for file_base in base_col:
		base_start_time = timeit.default_timer()
		read_to_taxids_map = {}
		reads_seq_map = {}
		#For every SAM file for a given sample, read in the SAM files.
		for sam_file in glob.glob(file_base + '*.sam'):
			file_start_time = timeit.default_timer()
			print('Reading in ' + sam_file)
			
			#For every line in the SAM file
			line_count = 0
			for line in open(sam_file):
				line_count += 1
				#Skip the first three lines, which are header
				if line_count > 3:
					#For each read, pull the SAM information for that read.
					line_list = line.split('\t')
					current_read = line_list[0]
					snap_assignment_of_current_read = line_list[2]
					sequence_of_current_read = line_list[9]
					
					#If read is unassigned, call it unassigned and throw it out.  Unassigned reads do not have an edit distance, assign it 100.
					if snap_assignment_of_current_read == '*':
						# higher edit distance than anything else makes sure this gets parsed out 
						current_read_taxid = [snap_assignment_of_current_read,100]
					else:
						#Pull the taxid and the edit distance from each line.
						current_read_taxid = [snap_assignment_of_current_read.split('#')[-1],
							int(line_list[12].split(':')[-1])]
					#Create map for each sample.
					#The key in each map is the read ID and the values are lists.
					#For every read, append a list of taxid and edit distance from each SAM file.
					if current_read in read_to_taxids_map:
						read_to_taxids_map[current_read].append(current_read_taxid)
					else: 
						# if this is the first time we've seen this read add the sequence to the list
						# and initialize the read -> taxid_list map  
						read_to_taxids_map[current_read] = [current_read_taxid]
						# also store the read and the sequence, this does need to be in a map 
						reads_seq_map[current_read] = sequence_of_current_read
	
			file_runtime = str(timeit.default_timer() - file_start_time)
			print('Reading in file ' + sam_file + ' took ' + file_runtime)
	
		per_base_runtime = str(timeit.default_timer() - base_start_time)
		print(file_base + ' took ' + per_base_runtime + ' in total to read')
		final_assignment_counts = {}
		
		#Now we've read all the reads for all the SAM files for one sample.  We are still within the sample For loop here.
		#We have all the reads with all the taxids and edit distances that have been assigned.
		print('Breaking ties for ' + file_base)
		tie_break_start = timeit.default_timer()
		# now we're done with the loop and we have a map with reads to list of taxids assigned to them
		g = open(file_base + '_assignments.txt', 'w')
		e = open(file_base + '_unassigned.txt','w')
		unass_count = 0
		taxid_to_read_set = {}
		
		
		for read_key in read_to_taxids_map.keys():
			# Now we need to assign only one taxid to each read, which means we need to run the tie-breaking function. 
			loaded_read = reads_seq_map[read_key]
			
			#Create a results list which is a taxid and a boolean, which answers whether I should re-BLAST this or not.
			r_list = tie_break(read_to_taxids_map[read_key])
			tax_assignment = r_list[0]
			
			#If the read is unassigned, write it to the unassigned file.
			if tax_assignment == '*':
				e.write('>' + read_key + '\n' + loaded_read + '\n')
				unass_count += 1
			# otherwise write it out to the read-by-read assignment file 
			else:
				g.write(read_key + '\t' + str(tax_assignment) + '\t' + loaded_read + '\n')
				# create a mapping of taxid -> unique reads.  Unique reads are defined as reads without the exact same sequence.  This can help in debugging.
				if WRITE_UNIQUES:
					if str(tax_assignment) in taxid_to_read_set:
						taxid_to_read_set[str(tax_assignment)].add(loaded_read,)
					else:
						taxid_to_read_set[str(tax_assignment)] = set([loaded_read])
					
			
					if tax_assignment in final_assignment_counts:
						final_assignment_counts[tax_assignment] += 1
					else:
						final_assignment_counts[tax_assignment] = 1
	
		g.close()
		e.close()
		#For each sample, we make a folder and for every taxid, we create a FASTA file that are named by their taxid.  We lose the read ID in this file.  #nicetohave would be hold the read ID here.
		#Here we will write a FASTA of unique reads
		if WRITE_UNIQUES:
			subprocess.call('mkdir ' + file_base.split('.')[0], shell=True)
			for id in taxid_to_read_set.keys():
				f = open(file_base.split('.')[0] + '/' + str(id) + '_uniques.txt', 'w')
				count = 0
				for item in taxid_to_read_set[id]:
					f.write('>' + str(count) + '\n')
					f.write(item + '\n')
					count += 1
				f.close()
			
		tie_break_time = str(timeit.default_timer() - tie_break_start)
		print('Tie breaking ' + file_base + ' took ' + tie_break_time)
		
		#For each sample, write the Pavian output.
		new_write_kraken(file_base, final_assignment_counts, unass_count)
	
		
