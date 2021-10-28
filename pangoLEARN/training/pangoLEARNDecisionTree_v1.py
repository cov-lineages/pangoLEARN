import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.datasets import make_classification
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report, confusion_matrix, make_scorer
from datetime import datetime
import joblib
import sys
import os
from sklearn.model_selection import cross_val_score
from Bio import SeqIO
import pickle

# file with lineage assignments
lineage_file = sys.argv[1]
# file with sequences
sequence_file = sys.argv[2]
# how much of the data will be used for testing, instead of training
testing_percentage = 0.0000000001

relevant_positions = pickle.load(open(sys.argv[5], 'rb'))
relevant_positions.add(0)

# the path to the reference file. 
# This reference sequence must be the same as is used in the pangolearn script!!
referenceFile = sys.argv[3]

# data storage
dataList = []
# dict for lookup efficiency
indiciesToKeep = dict()

referenceId = "Wuhan/WH04/2020"
referenceSeq = ""

idToLineage = dict()
idToSeq = dict()

mustKeepIds = []
mustKeepLineages = []


# function for handling weird sequence characters
def clean(x, loc):
	x = x.upper()
	
	if x == 'T' or x == 'A' or x == 'G' or x == 'C' or x == '-':
		return x

	if x == 'U':
		return 'T'

	# otherwise return value from reference
	return referenceSeq[loc]

def findReferenceSeq():
	with open(referenceFile) as f:
		currentSeq = ""

		for line in f:
			if ">" not in line:
				currentSeq = currentSeq + line.strip()

	f.close()
	return currentSeq


def getDataLine(seqId, seq):
	dataLine = []
	dataLine.append(seqId)

	newSeq = ""

	# for each character in the sequence
	for index in range(len(seq)):
		newSeq = newSeq + clean(seq[index], index)

	dataLine.append(newSeq)
	
	return dataLine


def readInAndFormatData():

	# add the data line for the reference seq
	idToLineage[referenceId] = "A"
	dataList.append(getDataLine(referenceId, referenceSeq))

	# create a dictionary of sequence ids to their assigned lineages
	with open(lineage_file, 'r') as f:
		for line in f:
			line = line.strip()

			split = line.split(",")

			idToLineage[split[0]] = split[1]

	# close the file
	f.close()

	seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(sequence_file, "fasta")}

	print("files read in, now processing")

	for key in seq_dict.keys():
		if key in idToLineage:
			dataList.append(getDataLine(key, seq_dict[key]))
		else:
			print("unable to find the lineage classification for: " + key)


# find columns in the data list which always have the same value
def findColumnsWithoutSNPs():

	# for each index in the length of each sequence
	for index in range(len(dataList[0][1])):
		keep = False

		# loop through all lines
		for line in dataList:

			# if there is a difference somewhere, then we want to keep it
			if dataList[0][1][index] != line[1][index] or index == 0:
				keep = True
				break

		# otherwise, save it
		if keep and index in relevant_positions:
			indiciesToKeep[index] = True


# remove columns from the data list which don't have any SNPs. We do this because
# these columns won't be relevant for a logistic regression which is trying to use
# differences between sequences to assign lineages
def removeOtherIndices(indiciesToKeep):

	# instantiate the final list
	finalList = []

	indicies = list(indiciesToKeep.keys())
	indicies.sort()

	# while the dataList isn't empty
	while len(dataList) > 0:

		# pop the first line
		line = dataList.pop(0)
		seqId = line.pop(0)

		line = line[0]
		# initialize the finalLine
		finalLine = []

		for index in indicies:
			if index == 0:
				# if its the first index, then that's the lineage assignment, so keep it
				finalLine.append(seqId)
			else:
				# otherwise keep everything at the indices in indiciesToKeep
				finalLine.append(line[index])

		# save the finalLine to the finalList
		finalList.append(finalLine)

	# return
	return finalList

def allEqual(list):
		entries = dict()

		for i in list:
			if i not in entries:
				entries[i] = True

		return len(entries) == 1

def removeAmbiguous():
	idsToRemove = set()
	lineMap = dict()
	idMap = dict()

	for line in dataList:
		keyString = ",".join(line[1:])

		if keyString not in lineMap:
			lineMap[keyString] = []
			idMap[keyString] = []
 
		if line[0] in idToLineage:
			lineMap[keyString].append(idToLineage[line[0]])
			idMap[keyString].append(line[0])
		else:
			print("diagnostics")
			print(line[0])
			print(keyString)
			print(line)
	for key in lineMap:
		if not allEqual(lineMap[key]):

			skipRest = False

			# see if any protected lineages are contained in the set, if so keep those ids
			for lineage in lineMap[key]:
				if lineage in mustKeepLineages:
					skipRest = True

					for i in idMap[key]:
						if lineage != idToLineage[i] and i not in mustKeepIds:
							idsToRemove.add(i)

			# none of the lineages are protected, fire at will
			if not skipRest:

				lineageToCounts = dict()

				aLineage = False
				# find most common lineage
				for lineage in lineMap[key]:
					if lineage not in lineageToCounts:
						lineageToCounts[lineage] = 0

					lineageToCounts[lineage] = lineageToCounts[lineage] + 1
					aLineage = lineage

				m = aLineage
				for lineage in lineageToCounts:
					if lineageToCounts[lineage] > lineageToCounts[m]:
						m = lineage


				for i in idMap[key]:
					if m != idToLineage[i]:
						idsToRemove.add(i)

	newList = []

	print("keeping indicies:")

	for line in dataList:
		if line[0] not in idsToRemove:
			print(line[0])
			line[0] = idToLineage[line[0]]
			newList.append(line)

	return newList


print("reading in data " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"), flush=True)

referenceSeq = findReferenceSeq()

readInAndFormatData()

print("processing snps, formatting data " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"), flush=True)

findColumnsWithoutSNPs()

dataList = removeOtherIndices(indiciesToKeep)

print("# sequences before blacklisting")
print(len(dataList))

dataList = removeAmbiguous()

print("# sequences after blacklisting")
print(len(dataList))

# headers are the original genome locations
headers = list(indiciesToKeep.keys())
headers[0] = "lineage"

print("setting up training " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"), flush=True)

pima = pd.DataFrame(dataList, columns=headers)

# nucleotide symbols which can appear
categories = ['A', 'C', 'G', 'T', '-']

# one hot encoding of all headers other than the first which is the lineage
dummyHeaders = headers[1:]

# add extra rows to ensure all of the categories are represented, as otherwise 
# not enough columns will be created when we call get_dummies
for i in categories:
	line = [i] * len(dataList[0])
	pima.loc[len(pima)] = line

# get one-hot encoding
pima = pd.get_dummies(pima, columns=dummyHeaders)

# get rid of the fake data we just added
pima.drop(pima.tail(len(categories)).index, inplace=True)

feature_cols = list(pima)
print(feature_cols)

# remove the last column from the data frame. This is because we are trying to predict these values.
h = feature_cols.pop(0)
X = pima[feature_cols]
y = pima[h]

# separate the data frame into testing/training data sets. 25% of the data will be used for training, 75% for test.
X_train,X_test,y_train,y_test=train_test_split(X,y,test_size=testing_percentage,random_state=0)

print("training " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"), flush=True)

header_out = os.path.join(sys.argv[4],"decisionTreeHeaders_v1.joblib")
joblib.dump(headers, header_out, compress=9)

# instantiate the random forest with 1000 trees
dt = DecisionTreeClassifier()

# fit the model
dt.fit(X,y)

print("testing " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"), flush=True)

# classify the test data
y_pred=dt.predict(X_test)

print(y_pred)

# get the scores from these predictions
y_scores = dt.predict_proba(X_test)

print("generating statistics " + datetime.now().strftime("%m/%d/%Y, %H:%M:%S"), flush=True)

#print the confusion matrix
print("--------------------------------------------")
print("Confusion Matrix")
cnf_matrix = metrics.confusion_matrix(y_test, y_pred)
print(cnf_matrix)

print("--------------------------------------------")
print("Classification report")
print(metrics.classification_report(y_test, y_pred, digits=3))

# save the model files to compressed joblib files
# using joblib instead of pickle because these large files need to be compressed
model_out = os.path.join(sys.argv[4],"decisionTree_v1.joblib")
joblib.dump(dt,  model_out, compress=9)

print("model files created", flush=True)

# this method is used below when running 10-fold cross validation. It ensures
# that the per-lineage statistics are generated for each cross-fold
def classification_report_with_accuracy_score(y_true, y_pred):
	print("--------------------------------------------")
	print("Crossfold Classification Report")
	print(metrics.classification_report(y_true, y_pred, digits=3))
	return accuracy_score(y_true, y_pred)

# optionally, run 10-fold cross validation (comment this out if not needed as it takes a while to run)
# cross_validation_scores = cross_val_score(dt, X=X, y=y, cv=10, scoring=make_scorer(classification_report_with_accuracy_score))
