'''
Given an input gtf file from a TALON run, generates 
individual gtf files separated on novelty type. 
Input: --c: a csv config file with the following fields
1: gtf file
2: 'n' for just known vs novel or 'n+' for known and different novelty types
    ie. ISM, NIC, NNC, etc
3: 1 to combine all ISMs (ie exclude ISM-suffix and ISM-prefix categories)
4: public-facing directory (accessible from a web browser) to move 
    output files to. if nonexistent, any invalid path will work
5: url prefix for directory location corresponding to public facing output directory.
    if nonexistent, any character set will do
'''
import argparse
import os
from collections import defaultdict
import shutil
import subprocess
import pprint
import datetime

# get and format output directory
def format_odir(odir):
    cwd = os.getcwd()
    if odir != '':
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
    odir = odir+bname+'_tracks/'

    if not os.path.isdir(odir):
        print('Making output directory '+odir)
        os.makedirs(odir)

    return odir

# get basename from a file and path string
def get_basename(filepath):
    return os.path.basename(os.path.splitext(filepath)[0])

# resets all bools in a dictionary to false
def reset_bool_dict(d):
    for key, item in d.items():
        d[key] = False
    return d 

# closes all filehandles in a dictionary and removes files that are empty
def close_file_dict(d):   
    keys_to_remove = []
    for key, item in d.items():
        temp = item.name
        item.close()
        d[key] = temp

        # check if the file was ever written to
        if not os.stat(temp).st_size:
            print()
            os.remove(temp)
            keys_to_remove.append(key)
    
    # remove entries in dict that were empty so stuff
    # that iterates through dict down the line isn't thrown off
    for k in keys_to_remove:
        del d[k]
    return d

def write_track(key, gtffile, fname, colors_dict, url):
    color = colors_dict[key]
    base = ['track name="$BASE Type $CLASS"',\
            'visibility=full','color=$COLOR']
    string = ' '.join(base)
    string = string.replace('$BASE', get_basename(gtffile))
    string = string.replace('$CLASS', key)
    string = string.replace('$COLOR', color)
    url += fname

    string = string+'\n'+url
    return string

# argument parse config stuff
parser = argparse.ArgumentParser(description=\
    'Separates a gtf file into novelty classifications and'+\
    ' generates a text file with track information to paste into genome browser')
parser.add_argument('--c', help='CSV config file with conversion options.'+\
    ' field 1: gtf file to separate'+\
    ' field 2: n=binary (known v. novel) n+=multiclass novelty conversion'+\
    ' field 3: combine ISMS 1 = true, 0 = false'+\
    ' field 4: public-facing directory to move files to'+\
    ' field 5: url prefix for the data to write to track file')
args = parser.parse_args()

# obtain arguments from config file
cfile = open(args.c, 'r')
for ind, line in enumerate(cfile):
    line = line.replace('\n', '').split(',')

    # arguments from config file
    gtffile = line[0]
    sep_type = line[1]
    if sep_type != 'n' and sep_type != 'n+':
        print('Incorrect option for separation type given. Please choose n or n+.')
        exit()
    combine_isms = bool(int(line[2]))
    pubdir = format_odir(line[3])
    temp = os.path.dirname(gtffile)
    odir = format_odir(temp)
    odir = make_dated_folder(odir, get_basename(gtffile))
    odir = odir + get_basename(gtffile)+'_'
    url = line[4]
    if url[-1] != '/':
        url += '/'

    # where to output final tracks to
    if not ind:
        # custom tracks file    
        tfile = odir+sep_type+'_tracks'
        # print(tfile)
        tfile = open(tfile, 'w')

    # classes to look for and their associated colors
    if sep_type == 'n':
        classes = ['KNOWN', 'NOVEL']
        colors = ['0,158,115', '230,159,0']
    elif sep_type == 'n+':
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

    # class-specific dictionaries
    ofiles = defaultdict() # file handles and later file names
    colors_dict = defaultdict() # colors for each class
    gene_written = defaultdict(bool) # has gene line for that novelty file been written
    transcript_class = defaultdict(bool) # whether the current transcript is of that class
    for i, c in enumerate(classes):
        colors_dict[c] = colors[i]
        gene_written[c] = False 
        transcript_class[c] = False
        ofiles[c] = open(odir+c+'.gtf', 'w')

infile = open(gtffile, 'r')
# print(infile.name)

gene_line = []
trans_line = []
i = 0
for line in infile:
    temp = line.split('\t')
    # if this is EBV stuff ignore it
    if temp[0] == 'chrEBV': 
        continue
    # if we've found a new gene, update the stored gene line
    if temp[2] == 'gene':
        gene_line = line
        gene_written = reset_bool_dict(gene_written)
    elif temp[2] == 'transcript':
        trans_line = line
        transcript_class = reset_bool_dict(transcript_class)
        fields = temp[8]
        temp = fields.split('transcript_status "')[1]
        temp = temp.split('";')[0]
        # check for different types of novelty
        if temp == 'NOVEL':
            if sep_type == 'n+':
                for c in classes:
                    if c+' "TRUE"' in fields:
                        transcript_class[c] = True
                        # print(c)
                        # print(gene_line)
                        if not gene_written[c]:
                            ofiles[c].write(gene_line)
                            gene_written[c] = True
                        # print(trans_line)
                        ofiles[c].write(trans_line)
                        # write gene
                        # write transcript
            elif sep_type == 'n':
                transcript_class['NOVEL'] = True
                # print('Novel')
                if not gene_written['NOVEL']:
                    ofiles['NOVEL'].write(gene_line)
                    gene_written['NOVEL'] = True
                ofiles['NOVEL'].write(trans_line)
                # print(trans_line)
        elif temp == 'KNOWN':
            transcript_class['KNOWN'] = True
            if not gene_written['KNOWN']:
                ofiles['KNOWN'].write(gene_line)
                gene_written['KNOWN'] = True
            ofiles['KNOWN'].write(trans_line)
    else:
        for c, fname in ofiles.items():
            if transcript_class[c]:
                # write exon/CDS/5 or 3' UTR
                ofiles[c].write(line)
    i += 1

# close all gtf files
ofiles = close_file_dict(ofiles)
infile.close()

# pprint.pprint(ofiles)

# generate trackfile
for c, fname in ofiles.items():
    s = write_track(c, gtffile, get_basename(fname)+'.gtf', colors_dict, url)
    tfile.write(s+'\n')
    # print(fname)

# move into odir if it exists
if os.path.isdir(pubdir):
    for _, ofile in ofiles.items():
        shutil.copy(ofile, pubdir)

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



