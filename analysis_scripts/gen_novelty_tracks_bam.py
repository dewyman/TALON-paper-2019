'''
Given an input bam file and the associated TALON database, generates 
individual bam files separated on novelty type. 
Input: --c: a csv config file with the following fields
1: bam file
2: database file
3: 'n' for just known vs novel or 'n+' for known and different novelty types
	ie. ISM, NIC, NNC, etc
4: 1 to combine all ISMs (ie exclude ISM-suffix and ISM-prefix categories)
5: public-facing directory (accessible from a web browser) to move 
	output files to. if nonexistent, any invalid path will work
6: url prefix for directory location corresponding to public facing output directory.
	if nonexistent, any character set will do
'''
import subprocess
import argparse
from collections import defaultdict
import os
import numpy as np
import pysam 
import shutil
import sqlite3
import datetime

# get and format output directory
def format_odir(odir):
	cwd = os.getcwd()

	# if first character is not /, use cwd to make this an absolute path
	if odir[0] != '/' and odir[0] != '~':
		odir = cwd+odir
	if odir[-1] != '/':
		odir += '/'
	return odir

# make a dated output directory for the files used for the tracks
def make_dated_folder(odir, bname):
	date = datetime.datetime.now()
	date = date.strftime('%y%m%d')
	odir = odir+date+'_'+bname+'_tracks/'

	if not os.path.isdir(odir):
		print('Making output directory '+odir)
		os.makedirs(odir)

	return(odir)

def write_track(key, fname, colors_dict, sfile, url):
	color = colors_dict[key]
	base = ['track type=bam', 'name="$BASENAME Type $CLASS"', 'bigDataUrl="$URL$FNAME"',\
			'name="novelty $CLASS"', 'visibility=full','color=$COLOR', 'coverageColor=$COLOR', 'bamColorTag=$COLOR', \
			'bamColorMode=off']
	string = ' '.join(base)
	string = string.replace('$CLASS', key)
	string = string.replace('$COLOR', color)
	string = string.replace('$FNAME', fname)
	string = string.replace('$BASENAME', get_basename(sfile))
	string = string.replace('$URL', url)
	return string

# get basename from a file and path string
def get_basename(filepath):
	return os.path.basename(os.path.splitext(filepath)[0])

# convert bam file to sam file so it can be parsed
def bam_to_sam(bamfile, odir, bname):
	samfile = odir+bname+'.sam'
	# print()
	# print('Using bamfile: '+bamfile)
	# print('to make samfile: '+samfile)
	pysam.view('-h', bamfile, '-o',samfile, catch_stdout=False)
	return samfile

# 
def query_db(db, sep_type):
	# begin writing query 
	base_query = """
		SELECT 
			read_name 
		FROM observed as o 
		LEFT JOIN transcript_annotations as ta ON ta.ID = o.transcript_ID 
		WHERE ("""

	# classes to look for and their associated colors
	if sep_type == 'n':
		classes = ['KNOWN', 'NOVEL']
		base_query += "ta.attribute = 'transcript_status' and ta.value = "
		colors = ['0,158,115', '230,159,0']
	elif sep_type == 'n+':
		base_query += "ta.attribute = "
		if not combine_isms:
			classes = ['KNOWN', 'ISM_transcript',\
					   'ISM-prefix_transcript', 'ISM-suffix_transcript',\
					   'NIC_transcript', 'NNC_transcript',\
					   'genomic_transcript','antisense_transcript','intergenic_transcript']
			colors = ['0,158,115','0,114,178','86,180,233','105,139,172',\
					  '213,94,0','230,159,0','240,228,66','0,0,0','204,121,167']
		else:
			classes = ['KNOWN', 'ISM_transcript',\
					   'NIC_transcript', 'NNC_transcript',\
					   'genomic_transcript','antisense_transcript','intergenic_transcript']
			colors = ['0,158,115','0,114,178',\
					  '213,94,0','230,159,0','240,228,66','0,0,0','204,121,167']
	# query DB for results 
	conn = sqlite3.connect(db)
	conn.row_factory = sqlite3.Row 
	cursor = conn.cursor()

	# class-specific dictionaries
	ofiles = defaultdict() # file handles and later file names
	colors_dict = defaultdict() # colors for each class
	read_id_dict = defaultdict(list) # novelty type for each read_id

	# loop through classes and build query db for each one
	for i, c in enumerate(classes):
		colors_dict[c] = colors[i]

		# build the rest of the query
		# weird case for FSM_transcript or KNOWN when doing the novelty types separation
		if c == 'KNOWN' and sep_type == 'n+':
			query = base_query+"'transcript_status' and ta.value = "+"'"+c+"'"+")"
		elif sep_type == 'n':
			query = base_query+"'"+c+"')"
		else: 
			query = base_query+"'"+c+"')"

		# print(query)
		cursor.execute(query)

		# assign a novelty status for each read id 
		read_ids = [k[0] for k in cursor.fetchall()]
		# print('number of reads: '+str(len(read_ids)))
		for read_id in read_ids:
			read_id_dict[read_id].append(c)

		# only open a db file if there is data associated with it
		# from the database
		if len(read_ids):
			ofiles[c] = open(odir+c+'.sam', 'w')

	# close db 
	conn.close()

	return ofiles, colors_dict, read_id_dict

# main
print()
# just parse for config file now
parser = argparse.ArgumentParser(description=\
	'Separates a bam file into novelty classifications and'+\
	' generates a text file with track information to paste into genome browser')
parser.add_argument('--c', help='CSV config file with conversion options.'+\
		' field 1: bam file to convert, '+\
		' field 2: database to use,'+\
		' field 3: n=binary (known v. novel) n+=multiclass novelty conversion'+\
		' field 4: combine ISMs 1 = true, 0 = false'+\
		' field 5: public-facing directory to move files to'+\
		' field 6: url prefix for the data to write to track file')
args = parser.parse_args()

# obtain arguments from config file
cfile = open(args.c, 'r')
for ind, line in enumerate(cfile):
	line = line.replace('\n','').split(',')

	# get arguments 
	bfile = line[0]
	db = line[1]

	# known/novel classifications or known/novelty types
	sep_type = line[2]
	if sep_type != 'n' and sep_type != 'n+':
		print('Incorrect option for separation type given. Please choose n or n+.')
		exit()

	combine_isms = bool(int(line[3]))

	# where to output files to
	temp = os.path.dirname(bfile)
	odir = format_odir(temp)

	# get basename
	bname = get_basename(bfile)

	# make dated output directory
	odir = make_dated_folder(odir, bname)

	# get public directory to try to move to automatically
	pubdir = format_odir(line[4])

	# url to write to trackline
	url = line[5]
	if url[-1] != '/':
		url += '/'

	# print()
	# print('Input arguments')
	# print('bamfile: '+bfile)
	# print('db file: '+db)
	# print('combine_isms: '+str(combine_isms))
	# print('sep type: '+sep_type)
	# print('ouput directory: '+odir)
	# print('basename: '+bname)
	# print('public directory: '+pubdir)
	# print('url: '+url)

	sfile = bam_to_sam(bfile, odir, bname)

	# where to output final tracks to
	if not ind:
		# custom tracks file	
		tfile = odir+bname+'_'+sep_type+'_tracks'
		tfile = open(tfile, 'w')

	# get db info 
	ofiles, colors_dict, read_id_dict = query_db(db, sep_type)

	# # check how many elements are in read id dict (maybe symlink messing up)
	# print()
	# print(len(read_id_dict.items()))
	# print()

	# loop through sam file, get novelty status of each read,
	# and dump to new sam files
	infile = open(sfile, 'r')
	qw = 0
	for line in infile:

		# get the read_id from each non-header line
		if line[0] != '@':
			read_id = line.split('\t')[0]

			# write to corresponding novelty sam file
			try:
				for c in read_id_dict[read_id]:
					ofiles[c].write(line)
			except:
				exceptions+=1

		# write all header lines to all output sam files
		else:
			for _, ofile in ofiles.items():
				ofile.write(line)

	# close and remove samfile
	infile.close()
	os.remove(sfile)

	# close all open files and get sam/bam names
	for key, ofile in ofiles.items():
		fname = ofile.name
		ofile.close()
		ofiles[key] = [fname]
		ofiles[key].append(os.path.splitext(ofiles[key][0])[0]+'.bam')
		ofiles[key].append(os.path.splitext(ofiles[key][0])[0]+'_sorted.bam')

	# convert each sam file to bam
	for key, ofile in ofiles.items():
		pysam.view('-Sb', ofile[0], '-o', ofile[1], catch_stdout=False)
		pysam.sort(ofile[1], '-o', ofile[2], catch_stdout=False)
		pysam.index(ofile[2], catch_stdout=False)

		# write track to custom tracks file
		temp = write_track(key, os.path.basename(ofile[2]), colors_dict, sfile, url)
		tfile.write(temp+'\n')

		# remove sam file and unsorted bam file
		os.remove(ofile[0])
		os.remove(ofile[1])

	# check if public directory exists
	if os.path.isdir(pubdir):
		# copy sorted bam and bai files to public directory
		for _, ofile in ofiles.items():
			# sorted bam
			shutil.copy(ofile[2], pubdir)
			# bam index (bai)
			baifile = ofile[2]+'.bai'
			shutil.copy(baifile, pubdir)

		# copy trackfile too but close first 
		tfile_name = tfile.name
		tfile.close()
		shutil.copy(tfile_name, pubdir)
	else:
		print('Public dir '+pubdir+' does not exist.'+\
		      ' Run move_tracks.py or manually move files.')
		tfile.close()

# close config file
cfile.close()

