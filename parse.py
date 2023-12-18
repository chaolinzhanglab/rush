import os
import argparse

def parse(homedir, libSeqCSV, unsorted, top5, bot5):

	nameList = [unsorted, top5, bot5]

	# create a new folder under the home directory for output files
	mainFolder = os.path.join(homedir, 'output/')

	if not os.path.exists(mainFolder):
		os.mkdir(mainFolder)

	unsortedList = list()
	topList = list()
	botList = list()

	# loop through all files under the home dir 
	filenames = os.listdir(homedir)
	for file in filenames:
		if file.endswith('.fastq.gz'):
			# store fastq files containing the same key word in the same group
			if nameList[0] in file:
				unsortedList.append(file)

			elif nameList[1] in file:
				topList.append(file)

			elif nameList[2] in file:
				botList.append(file)

	# concatenate fastq files in each group
	compiled_fastqGroups = [unsortedList, topList, botList]

	c = 0
	for fastqFiles in compiled_fastqGroups:
		cat_output = os.path.join(mainFolder, nameList[c] + '.fastq.gz')
		os.system('cat ' + fastqFiles[0] + ' ' + fastqFiles[1] + ' ' + fastqFiles[2] + ' ' + fastqFiles[3] + ' > ' + cat_output)
		# unzip the combined fastq file
		os.system('gunzip ' + cat_output)
		# check # of rows in each combined fastq file
		os.system('wc -l ' + os.path.join(mainFolder, nameList[c] + '.fastq'))

		c+=1

	# loop through all the combined fastq files in the main folder
	mainFolder_files = os.listdir(mainFolder)

	for fastq_input in mainFolder_files:
		# get the prefix of the fastq file
		groupName = fastq_input.split('.')[0]

		# define input and output file directories
		trimmed_input = os.path.join(mainFolder, fastq_input)
		trimmed_output = os.path.join(mainFolder, groupName + '.trimmed.fastq')
		trimmed_logtxt = os.path.join(mainFolder, groupName + '.txt')
		TruSeq3se = os.path.join(homedir, 'TruSeq3-SE.fa')

		cutadapt_output = os.path.join(mainFolder, groupName + 'ready.fastq')

		libCountCSV = os.path.join(mainFolder, groupName + '_library_count.csv')


		# Step 2 - trim adapters from reads
		os.system('trimmomatic SE -phred33 -trimlog ' + trimmed_logtxt + ' ' + trimmed_input + ' ' + trimmed_output + ' ILLUMINACLIP:' + TruSeq3se +':2:0:10 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:36')

		# Step 3 - remove primers
		os.system('cutadapt -g ggaaaggacgaaacaccgtacccctacc -a ttttttgaattcgctagctaggt --discard-untrimmed -o ' + cutadapt_output + ' ' + trimmed_output)

		# Step 4 - count gRNAs
		os.system('python ' + os.path.join(homedir, 'count_spacers_py3.py') + ' -i ' + libSeqCSV + ' -f ' + cutadapt_output + ' -o ' + libCountCSV + ' -no-g')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='process sequencing library data and find matched gRNAs')
    parser.add_argument('-d', '--homedir', type=str, dest='homedir', help='home directory')
    parser.add_argument('--libcsv', type=str, dest='libSeqCSV', help='library sequences csv', default='library_sequences.csv')
    parser.add_argument('-u', type=str, dest='unsorted', help='group prefix', default='unsorted')
    parser.add_argument('-t', type=str, dest='top5', help='group prefix', default='top5')
    parser.add_argument('-b', type=str, dest='bot5', help='group prefix', default='bot5')
    args = parser.parse_args()

    parse(args.homedir, args.libSeqCSV, args.unsorted, args.top5, args.bot5)



