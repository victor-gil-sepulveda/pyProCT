## Author: Victor Gil Sepulveda
## Started: 17/09/10
## First release: 19/Sep/2010
## Version: 1.0
##
## Description:
## ------------
##
## This file includes some functions to work with plop generated trajectory files. Its main aim is to allow the user
## to select and write the different models of the trajectory using a simple but powerful selection languaje.
##
## Functions:
## ----------
##
##-------------------------------------------
## processDir(dir, file_name_commnon_part)
##-------------------------------------------
## This function will create a list of metric records (which are dictionaries) from the trajectory files in a given
## directory 'dir'.
## You can also specify the common part of the name of the files, or nothing if you want to parse ALL the files inside
## the folder.
##
## Examples:
##
## records = processDir("f","traj")
##
## This sentence will return an array of records (stored in 'records') for all the metrics in all the files in the 'f'
## directory where 'traj' is PART OF its name (this means it will also parse traj.01.pdb, traj.02.pdb... but also
## mytraj.pdb, atrajh.pdb...)
##
## records = processDir("f")
##
## In this case it will parse ALL the files inside the 'f' folder.
##
##
##-------------------------------------------
## filterRecords(expression,records)
##-------------------------------------------
## Given a boolean expression written in 'expression', it will choose and return a subset of the 'records' list parameter
## with all the records that fulfill 'expression'
##
## Examples:
##
## Think that you have this metrics stored inside your .traj files : 'energy', 'totale', 'metrop', 'proc';
## and you want to know which models in processor 1 have energy below -26759. Just use this function like this:
##
## selection = filterRecords("Proc == 1 and Energy<-26759",records)
##
## In your boolean expression write any python-compliant sentence and it will do the trick. It's also case-insensitive.
##
## Do you want to extract models with energy between two values globally declared? you can!:
##
## X = -10000
## Y = -20000
## selection = filterRecords("(Energy<X and Energy>Y)",records)
##
## And in general you can use any complex expression using boolean operators and parentheses.
##
##
##-------------------------------------------
## genSingleTraj(name,records,selection)
##-------------------------------------------
## This function will write a single file which name will be 'name' (where you can also define a local or global path)
## using the records stored in 'selection'. The 'records' parameter will be the list obtained when using processDir.
##
##
##-------------------------------------------
## genMetricsFile(name, metrics, selection)
##-------------------------------------------
## This function will write a single file with the metrics in 'metrics' in colums. 'name' and 'selection' are the file
## name to write and the selected records.


import sys
import os
import numpy

def toNumber(n):
	d = 0
	try:
		d = float(n)
	except:
		try:
			d = int(n)
		except:
			return 0
	return d

def processDir(directory, file_name_commnon_part=""):
	dirList = []
	allList=os.listdir(directory)
	for n in allList:
		if file_name_commnon_part in n:
			dirList.append(n)

	records = []
	for i,fname in enumerate(dirList):
		print "Processing %s ( %d of %d )" %(fname, i+1, len(dirList))
		processFile(os.path.join(sys.argv[1],fname), records, i==0)
	return records

def processFile(filename, records, first):
	line_number = 0
	last_was_remark = False
	file = open(filename)
	for l in file:
		if l[0:7] == "REMARK ":
			if last_was_remark == False:
				## Create a new record
				if first == True:
					records.append({"file":filename, "seqres":True})
				else:
					records.append({"file":filename})
			parts = l.split()
			if len(parts) > 2:
				# Then is a remark of type
				#REMARK  TOTALE            -8690.283
				#REMARK  Steps|              626.000
				#REMARK  L1  Binding Ene|    -81.535

				index = "".join(parts[1:-1])
				if index[-1] == "|":
					index = index[:-1]

				records[-1][index.lower()] = toNumber(parts[-1])
			last_was_remark = True
		else:
			last_was_remark = False
			try:
				records[-1]["body"].append(line_number)
			except:
				try:
					records[-1]["body"].append(line_number)
				except :
					records[-1]["body"] = [line_number]

		line_number = line_number + 1
		first = False
	del records[-1]["body"][1:-1]
	file.close()

## Always use spaces
def filterRecords(expression,records):
	## Format the string
	tags = ["not","and","or",">","<",">=","<=","==","+","-","(",")"]
	expression = expression.lower()
	for t in tags:
		expression = expression.replace(t," "+t+" ")

	## Get all keys from records
	all_keys = set([])
	for r in records:
		all_keys = all_keys | set (r.keys())

	## Identify metrics in expression
	words = expression.split()
	for w in all_keys:
		if w in words:
			expression = expression.replace(w,"r[\""+w+"\"]")

	## Delete SEQRES records
	preselection = []
	for r in records:
		if  not "seqres" in r.keys():
			preselection.append(r)

	selection = []
	for r in preselection:
		if eval(expression):
			selection.append(r)
		print(eval(expression))

	return selection


def copyChunck(origin, to, start,end):
	file = open(origin,"r")
	lines = file.readlines()
	file.close()
	print(origin,start,end)

	for i in range(start,end+1):
		to.write(lines[i])


def genSingleTraj(name, records, selection):
	file = open(name,"w")

	## Write seqres
	seq_record = records[0]
	for r in records:
		if "seqres" in r.keys():
			seq_record = r
			break

	copyChunck(r['file'], file, r['body'][0],r['body'][1])

	## Now write selection
	for r in selection:
		if not 'seqres' in r.keys():
			copyChunck(r['file'], file, r['body'][0],r['body'][1])

	file.close()

def genMetricsFile(name, metrics, selection):
	numpy.savetxt(name, genMetrics(metrics, selection))

def genMetrics(metrics, selection):
	filtered_metrics = []
	for r in selection:
		this_metrics = []
		for m in metrics:
			try:
				this_metrics.append(r[m.lower()])
			except KeyError:
				this_metrics.append(0)
		filtered_metrics.append(this_metrics)
	return numpy.array(filtered_metrics)
