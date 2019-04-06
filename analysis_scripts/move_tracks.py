'''
Given an input directory, output directory, and (optionally) a url prefix,
will copy the contents of the input directory to the output directory and
(optionally) update the url in the output *_tracks file from either
gen_novelty_tracks_gtf.py or gen_novelty_tracks_bam.py.
Input:
	--i: input directory containing all files that need to be copied
	--o: output/destination directory to which all files from input will be copied
	(optional) --u: url prefix corresponding to output directory to update
		in the *_tracks file
'''
import shutil
import argparse
import os

# get and format output directory
def format_odir(odir):
	cwd = os.getcwd()

	# if first character is not /, use cwd to make this an absolute path
	if odir[0] != '/' and odir[0] != '~':
		if odir[-1] != '/':
			odir = cwd+'/'+odir
	if odir[-1] != '/':
		odir += '/'
	return odir

# set up input parser
parser = argparse.ArgumentParser(description=\
	'Moves a directory into another directory. Also replaces URL with '+\
	' a new URL if provided to the trackfile')
parser.add_argument('--i', help='Directory with stuff that needs to be moved')
parser.add_argument('--o', help='Public-facing directory to move files to.')
parser.add_argument('--u', help='URL to access the public facing directory OPTIONAL')
args = parser.parse_args()

# grab input arguments
indir = format_odir(args.i)
odir = format_odir(args.o)
if args.u:
	url = args.u
else:
	url = None

# get list of files in the input directory
file_list = os.listdir(indir)

# are we working with bam or gtf files?
files = [f for f in file_list if f[-6:] != 'tracks'][0]
ext = os.path.splitext(files)[1]
if ext == '.gtf':
	ftype = 'gtf'
else:
	ftype = 'bam'

print('Found '+ftype+' files')

# if we were given a url, replace the url field in the track file
if url:

	if url[-1] == '/':
		url = url[:-1]

	track_file = [f for f in file_list if f[-6:] == 'tracks'][0]
	track_file = open(indir+track_file, 'r')
	temp_out = open(track_file.name+'_temp', 'w')

	for line in track_file:
		if ftype == 'bam':
			temp = line.split('bigDataUrl="')[1].split()[0].replace('"', '')
			prefix = '/'.join((temp.split('/')[:-1]))
			line = line.replace(prefix, url)
		elif ftype == 'gtf':
			if line[:4] == 'http': # url-containing line
				prefix = '/'.join(line.split('/')[:-1])
				line = line.replace(prefix, url)

		temp_out.write(line)

	tf = track_file.name
	to = temp_out.name
	track_file.close()
	temp_out.close()
	os.rename(to, tf)

# move directory contents to output directory
if os.path.isdir(odir):
	for file in file_list:
		shutil.copy(indir+file, odir+file)
else: 
	print('Output directory does not exist. Check input and try again.')



