import os
import numpy as np
import pandas as pd
from scipy import stats
import csv
from utils import *
import warnings

# disable warning
warnings.filterwarnings("ignore")

def analyze_screen(inputDir, oligoLibDir):
	# set output directory 
	outputPath = os.path.join(inputDir, 'screen_results/')

	# create empty lists for storing data
	resultList = list()
	nameList = list()

	# load oligo library info
	oligoLibrary = getOligoLibInfo(oligoLibDir)

	# list the folders under the input directory
	inputDir_folders = os.listdir(inputDir)

	counter = 0
	for folder in inputDir_folders:
		# only process relevant folders
		if 'dual' in folder.lower():
			# get folder name
			nameList.append(folder.split('/')[-2])
		else:
			continue

		# loop through the files in each folder to get raw data
		filenames = os.listdir(folder)
		for f in filenames:
			inputRaw = pd.read_csv(os.path.join(folder, f), header = None)

			if 'top' in f:
			 	topBFPpos = getCounts(inputRaw, oligoLibrary)
			 	topBFPpos.name = 'topBFPpos_counts'

			elif 'bot' in f:
				botBFPpos = getCounts(inputRaw, oligoLibrary)
			 	botBFPpos.name = 'botBFPpos_counts'

			else:
				unsorted = getCounts(inputRaw, oligoLibrary)
				unsorted.name = 'unsorted_counts'

		# calculate RPM
		topBFPpos_RPM = getRPM(topBFPpos)
		topBFPpos_RPM.name = 'topBFPpos_RPM'

		botBFPpos_RPM = getRPM(botBFPpos)
		botBFPpos_RPM.name = 'botBFPpos_RPM'

		unsorted_RPM = getRPM(unsorted)
		unsorted_RPM.name = 'unsorted_RPM'

		allCounts = pd.concat([topBFPpos, botBFPpos, unsorted], axis = 1)
		allRPM = pd.concat([topBFPpos_RPM, botBFPpos_RPM, unsorted_RPM], axis = 1)

		# calculate zscores and relevant values for each screening experiment
		TopVsBot = calz(topBFPpos, botBFPpos, topBFPpos_RPM, botBFPpos_RPM, 'TopVsBot')
		TopVsUnsrt = calz(topBFPpos, unsorted, topBFPpos_RPM, unsorted_RPM,'TopVsUnsrt')
		BotVsUnsrt = calz(botBFPpos, unsorted, botBFPpos_RPM, unsorted_RPM, 'BotVsUnsrt')

		# combine processed data with oligo library info
		screenResults = pd.concat([oligoLibrary, allCounts, allRPM, TopVsBot, TopVsUnsrt, BotVsUnsrt], axis = 1)

		# output file as csv and extract data for later calculataions
		screenResults.to_csv(os.path.join(outputPath, nameList[counter] + '-summary.csv'))

		results_z = screenResults[['TopVsUnsrt_normalized_z', 
									'TopVsUnsrt_filter', 
									'BotVsUnsrt_normalized_z', 
									'BotVsUnsrt_filter']]

		results_z.rename(columns = {'TopVsUnsrt_normalized_z': nameList[counter] + '_TopVsUnsrt_nlzd_z', 'BotVsUnsrt_normalized_z': nameList[counter] + '_BotVsUnsrt_nlzd_z'}, inplace = True)

		# store relevant data from each screen to a list
		resultList.append(results_z)

		counter+=1

	# combine filter values of all screens
	nCombined_TopVsUnsrt = sum(resultList[x]['TopVsUnsrt_filter'] for x in range(len(nameList)))
	nCombined_BotVsUnsrt = sum(resultList[x]['BotVsUnsrt_filter'] for x in range(len(nameList)))
	nCombined = nCombined_TopVsUnsrt + nCombined_BotVsUnsrt
	nCombined.name = 'n_combined'

	valid_filter = nCombined >= len(nameList)
	final_filter= valid_filter.astype('int')
	final_filter.name = 'filter'

	######################################################
	# create an empty list to combine normalized zscores #
	######################################################

	zList = list()

	for s in range(len(nameList)):
		exp = resultList[s]

		if 'dualin' in nameList[s].lower():
		 	zList.append(exp[nameList[s] + '_TopVsUnsrt_nlzd_z']*exp['TopVsUnsrt_filter'] - exp[nameList[s] + '_BotVsUnsrt_nlzd_z']*exp['BotVsUnsrt_filter'])

		else:
		 	zList.append(exp[nameList[s] + '_BotVsUnsrt_nlzd_z']*exp['BotVsUnsrt_filter'] - exp[nameList[s] + '_TopVsUnsrt_nlzd_z']*exp['TopVsUnsrt_filter'])

	final_z = sum(zList)/np.sqrt(nCombined)
	final_z_nan2zeros = final_z.fillna(0)
	final_z_nan2zeros.name = 'combined_z'

	# calculate p-vales and FDR
	pVals = calculatePVal(final_z_nan2zeros)
	pVals.name = 'p_value'

	fdrVals = calculateFDR(pVals, final_filter)
	fdrVals.name = 'FDR'

	eachScreenZ = pd.concat([resultList[x][{nameList[x] +'_TopVsUnsrt_nlzd_z', nameList[x] +'_BotVsUnsrt_nlzd_z'}] for x in range(len(nameList))], axis = 1)
	output_df = pd.concat([oligoLibrary, eachScreenZ, final_z_nan2zeros, pVals, nCombined, final_filter, fdrVals], axis = 1)

	# output final summary as csv file
	output_df.to_csv(os.path.join(outputPath, 'combined-summary.csv'))


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='analyze count csv files to interpret screen results')
	parser.add_argument('-i', type=str, dest='inputDir', help='input directory')
	parser.add_argument('-l', type=str, dest='oligoLibDir', help='oligo library directory')

	args = parser.parse_args()

	analyze_screen(args.inputDir, args.oligoLibDir)






