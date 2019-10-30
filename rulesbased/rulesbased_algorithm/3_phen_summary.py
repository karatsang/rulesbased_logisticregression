#Counts number of isolates that antibiotic resistance was correctly, over, and under predicted for each antibiotic tested
#Outputs onto terminal
#python phen_summary.py phenotype_comparison.txt

import sys
import csv
import json
import os
import glob
from collections import Counter

correct=[]
corrects=[]
overpred=[]
underpred=[]
correct2 =[]
with open(sys.argv[-1], 'r') as pc_file:
	pc_file = csv.reader(pc_file, delimiter='\t')
	header = next(pc_file, None)
	for row in pc_file:
		if row[3]:
			correct.append(row[3])
		if row[4]:
			corrects.append(row[4])
		if row[5]:
			overpred.append(row[5])
		if row[6]:
			underpred.append(row[6])

correct=str(correct).replace("'", "")
correct=correct.strip("[]")
correct=correct.replace(" ", "_")
correct=correct.replace(",_", " ")
correct=correct.replace(',', '')
correct=(correct.split())
print("CORRECT R PREDICTION")
print(Counter(correct))


corrects=str(corrects).replace("'", "")
corrects=corrects.strip("[]")
corrects=corrects.replace(" ", "_")
corrects=corrects.replace(",_", " ")
corrects=corrects.replace(',', '')
corrects=(corrects.split())
print("CORRECT S PREDICTION")
print(Counter(corrects))

print(">>>>>>>>")

overpred=str(overpred).replace("'", "")
overpred=overpred.strip("[]")
overpred=overpred.replace(" ", "_")
overpred=overpred.replace(",_", " ")
overpred=overpred.replace(',', '')
overpred=(overpred.split())
print("OVER PREDICTION")
print(Counter(overpred))


print(">>>>>>>>")
underpred=str(underpred).replace("'", "")
underpred=underpred.strip("[]")
underpred=underpred.replace(" ", "_")
underpred=underpred.replace(",_", " ")
underpred=underpred.replace(',', '')
underpred=(underpred.split())
print("UNDER PREDICTION")
print(Counter(underpred))