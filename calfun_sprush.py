# sprush.py>

import pandas as pd
import numpy as np
from scipy import stats


# Get oligo library info
def getOligoLibInfo(oligoDataDir):
	# get targeting gRNA data
	oligo_gRNA = pd.read_excel(oligoDataDir, sheet_name = 'No-Capture') # import data
	oligo_gRNA_filtered = oligo_gRNA[oligo_gRNA['Barcode No.'] == 1] # only get oligos with Barcode 1
	oligo_gRNAinfo = oligo_gRNA_filtered[['Unique ID for Oligo', 'Reverse Complement [Oligo sequence] ']]

	# rename columns
	dfgRNA = oligo_gRNAinfo.rename(columns = {'Unique ID for Oligo' : 'OligoID', 'Reverse Complement [Oligo sequence] ' : 'Reverse_Complement_[Oligo sequence]'})

	# NT gRNA data
	oligoNT = pd.read_excel(oligoDataDir, sheet_name = 'No-Capture-nontargeting') # import data
	oligoNT_filtered = oligoNT[oligoNT['Barcode No.'] == 1] # only get oligos with Barcode 1
	oligoNTinfo = oligoNT_filtered[['Unique ID for Oligo', 'nontargeting sequence 5\'-3\'']]

	# rename columns
	dfNT = oligoNTinfo.rename(columns = {'Unique ID for Oligo' : 'OligoID', 'nontargeting sequence 5\'-3\'' : 'Reverse_Complement_[Oligo sequence]'})

	# get AS data
	oligoAS = pd.read_excel(oligoDataDir, sheet_name = 'sense-ctrl') # import data
	oligoAS_filtered = oligoAS[oligoAS['Barcode No.'] == 1] # only get oligos with Barcode 1
	oligoASinfo = oligoAS_filtered[['Unique ID for Oligo', 'gRNA Minigene Target Sequence ']]

	# rename columns
	dfAS = oligoASinfo.rename(columns = {'Unique ID for Oligo' : 'OligoID', 'gRNA Minigene Target Sequence ': 'Reverse_Complement_[Oligo sequence]'})

	# concatenate relevant oligo info df
	oligoLibrary = pd.concat([dfgRNA, dfNT, dfAS])

	# convert datatype
	oligoLibrary['OligoID'] = oligoLibrary['OligoID'].astype('Int64')              

	# convert sequence string to uppercase
	oligoLibrary['Reverse_Complement_[Oligo sequence]'] = oligoLibrary['Reverse_Complement_[Oligo sequence]'].str.upper()

	# reset index
	oligoLibrary.reset_index(drop = True, inplace = True) 

	return oligoLibrary


# Sort spacer sequneces in raw data according to oligo library info and get counts
def getCounts(rawData, oligoLib):
	Idx = [oligoLib[oligoLib.iloc[:,1] == s].index.tolist()[0] for s in rawData.iloc[:,0]]
	rawData['idx'] = Idx
	rawData.sort_values(by = 'idx', inplace = True)
	rawData.reset_index(drop = True, inplace = True)
	counts = rawData.iloc[:,1]

	return counts

def getRPM(counts):
	# make a copy of the original data to avoid unintended modifications
	counts_replace1 = counts.copy()

	# replace values smaller than 1 with 1 for calculating RPM
	counts_replace1[counts_replace1 < 1] = 1

	# calculate count RPM
	counts_RPM = counts_replace1/sum(counts)*10**6

	return counts_RPM


# Calculate zscores
def calz(xCounts, yCounts, xCounts_RPM, yCounts_RPM, comp):

	pct = sum(xCounts)/(sum(xCounts)+ sum(yCounts))

	sumCounts = xCounts + yCounts
	sumCounts[sumCounts < 1] = 1

	# get max count and max RPM
	Counts = pd.concat([xCounts, yCounts], axis = 1)
	RPM = pd.concat([xCounts_RPM, yCounts_RPM], axis = 1)

	maxCount = Counts.max(axis = 1)
	maxRPM = RPM.max(axis = 1)

	# calculate log2FC
	log2FC = np.log2(xCounts_RPM/yCounts_RPM)

	# calculate p
	P = xCounts/sumCounts
	P = P.astype('float')

	calP = P*(1-P)  
	calP[calP < 10**-6] = 10**-6 

	# calculate zcores and normalized zscores
	z = (P - pct)/np.sqrt(calP/(sumCounts))
	zdif = abs(z - np.median(z))
	nlzd_z = (z - np.median(z))/(np.median(zdif)/0.6745)

	# identify the rows with counts passing threshold and label these rows 
	flt = pd.DataFrame(np.repeat(0, z.shape[0]))
	flt_idx = np.where([[maxRPM[rw] >= 20 and xCounts[rw] >= 5 and yCounts[rw] >= 5] for rw in range(z.shape[0])])[0]
	flt.iloc[flt_idx] = 1 

	# store calculated values as df
	results = pd.DataFrame({comp + '_log2FC': log2FC,
							comp + '_P': P,
							comp + '_z': z,
							comp + '_zdif': zdif,
							comp + '_normalized_z': nlzd_z,
							comp + '_filter': flt.squeeze()})
	return results


def calculatePVal(final_z):
	cdf = pd.DataFrame(stats.norm.cdf(z, loc = 0, scale = 1) for z in final_z)
	cdf_1 = pd.DataFrame(1-cdf)
	concat_cdf = pd.concat([cdf, cdf_1], axis = 1)
	pVals = concat_cdf.min(axis = 1)*2

	return pVals

# calculate FDR
def calculateFDR(pVals, filterData):
	# replace p-values of low count data with NaN to avoid ranking errors
	pVals[filterData[filterData == 0].index.values] = np.nan
	
	# rank p-values
	ranked_pVals = pVals.rank(method = 'first', na_option = 'bottom')

	# only consider p-values of valid data for calculating FDR
	adj_pVals = pVals * len(pVals.dropna())/ranked_pVals
	adj_pVals[adj_pVals > 1] = 1

	return adj_pVals
