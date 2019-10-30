# Compares genotype predictions with lab-tested phenotype
# Uses gp.txt files from each isolate and looks for phenotypes in the ORF phenotype file
# Make sure all gp.txt are in the same directory
# Create file called phenotype_comparison.txt with correctly, over and under predicted antibiotic resistances 

# To run:
# python compare_phen.py ast.tsv

import sys
import csv
import json
import os
import glob

predp_csvs = glob.glob('*.gp.txt')

with open("phenotype_comparison.txt", "w") as af:
	writer = csv.writer(af, delimiter='\t', dialect='excel')
	writer.writerow(["Sample ID", "Tested R Phenotype", "Predicted R Phenotype", "Correct R Prediction", "Correct S Prediction", "Overprediction", "Underprediction"])
	ampicillin=0
	amox=0
	amikacin=0
	cefazolin=0
	cefalotin=0
	ciprofloxacin=0
	cefixime=0
	ceftazidime=0
	gentamicin_C=0
	meropenem=0
	nitrofurantoin=0
	pip=0
	tetracycline=0
	tobramycin=0
	tri=0
	cefoxitin=0
	ceftriaxone=0
	ertapenem=0
	for fr in predp_csvs:
		all_predp_dict ={}
		predp_dict={}
		predsus_dict={}
		testp_dict ={}
		testsus_dict={}
		if len(fr.split('_')[0]) == 3:
			sample="C0"+(fr.split("_")[0])
		elif len(fr.split('_')[0]) == 2:
			sample="C00"+(fr.split("_")[0])
		elif len(fr.split('.')[0]) == 5:
			sample=fr.split('.')[0]
		else:
			sample=fr.split("_")[0]
		print(sample)
		with open(fr, 'r') as predp_file:
			predp_file = csv.reader(predp_file, delimiter='\t')
			header = next(predp_file, None)
			for row in predp_file:
				if row[0]:
					all_predp_dict.setdefault(sample,[]).append(row[3])
		with open(sys.argv[-1], 'r') as test_phen_file:
			test_phen_file = csv.reader(test_phen_file, delimiter='\t')
			header = next(test_phen_file, None)
			for row in test_phen_file:
				if sample in row[0]:
					if row[1] == "R" or row[1] == "RA" or row[1] == "I":
						testp_dict.setdefault(sample,[]).append("ampicillin")

					if row[1] == "S":
						testsus_dict.setdefault(sample,[]).append("ampicillin")

					if row[2] == "R" or row[2] == "RA" or row[2] == "I":
						testp_dict.setdefault(sample,[]).append("amoxicillin-clavulanic_acid")

					if row[2] == "S":
						testsus_dict.setdefault(sample,[]).append("amoxicillin-clavulanic_acid")
				
					if row[3] == "R" or row[3] == "RA" or row[3] == "I":
						testp_dict.setdefault(sample,[]).append("amikacin")
	
					if row[3] == "S":
						testsus_dict.setdefault(sample,[]).append("amikacin")							
					if row[4] == "R" or row[4] == "RA" or row[4] == "I":
						testp_dict.setdefault(sample,[]).append("cefazolin")

					if row[4] == "S":
						testsus_dict.setdefault(sample,[]).append("cefazolin")											
					if row[5] == "R" or row[5] == "RA" or row[5] == "I":
						testp_dict.setdefault(sample,[]).append("cefalotin")

					if row[5] == "S":
						testsus_dict.setdefault(sample,[]).append("cefalotin")						
					if row[6] == "R" or row[6] == "RA" or row[6] == "I":
						testp_dict.setdefault(sample,[]).append("ciprofloxacin")

					if row[6] == "S":
						testsus_dict.setdefault(sample,[]).append("ciprofloxacin")							
					if row[7] == "R" or row[7] == "RA" or row[7] == "I":
						testp_dict.setdefault(sample,[]).append("cefixime")

					if row[7] == "S":
						testsus_dict.setdefault(sample,[]).append("cefixime")							
					if row[8] == "R" or row[8] == "RA" or row[8] == "I":
						testp_dict.setdefault(sample,[]).append("ceftazidime")

					if row[8] == "S":
						testsus_dict.setdefault(sample,[]).append("ceftazidime")						
					if row[9] == "R" or row[9] == "RA" or row[9] == "I":
						testp_dict.setdefault(sample,[]).append("gentamicin_C")

					if row[9] == "S":
						testsus_dict.setdefault(sample,[]).append("gentamicin_C")						
					if row[10] == "R" or row[10] == "RA" or row[10] == "I":
						testp_dict.setdefault(sample,[]).append("meropenem")

					if row[10] == "S":
						testsus_dict.setdefault(sample,[]).append("meropenem")												
					if row[11] == "R" or row[11] == "RA" or row[11] == "I":
						testp_dict.setdefault(sample,[]).append("nitrofurantoin")

					if row[11] == "S":
						testsus_dict.setdefault(sample,[]).append("nitrofurantoin")						
					if row[12] == "R" or row[12] == "RA" or row[12] == "I":
						testp_dict.setdefault(sample,[]).append("piperacillin-tazobactam")

					if row[12] == "S":
						testsus_dict.setdefault(sample,[]).append("piperacillin-tazobactam")							
					if row[13] == "R" or row[13] == "RA" or row[13] == "I":
						testp_dict.setdefault(sample,[]).append("tetracycline")	

					if row[13] == "S":
						testsus_dict.setdefault(sample,[]).append("tetracycline")						
					if row[14] == "R" or row[14] == "RA" or row[14] == "I":
						testp_dict.setdefault(sample,[]).append("tobramycin")

					if row[14] == "S":
						testsus_dict.setdefault(sample,[]).append("tobramycin")						
					if row[15] == "R" or row[15] == "RA" or row[15] == "I":
						testp_dict.setdefault(sample,[]).append("trimethoprim-sulfamethoxazole")

					if row[15] == "S":
						testsus_dict.setdefault(sample,[]).append("trimethoprim-sulfamethoxazole")								
					if row[16] == "R" or row[16] == "RA" or row[16] == "I":
						testp_dict.setdefault(sample,[]).append("cefoxitin")

					if row[16] == "S":
						testsus_dict.setdefault(sample,[]).append("cefoxitin")						
					if row[17] == "R" or row[17] == "RA" or row[17] == "I":
						testp_dict.setdefault(sample,[]).append("ceftriaxone")

					if row[17] == "S":
						testsus_dict.setdefault(sample,[]).append("ceftriaxone")							
					if row[18] == "R" or row[18] == "RA" or row[18] == "I":
						testp_dict.setdefault(sample,[]).append("ertapenem")

					if row[18] == "S":
						testsus_dict.setdefault(sample,[]).append("ertapenem")
		
		for k,v in all_predp_dict.items():
			v=str(v).replace("'", "")
			v=v.strip("[]")
			v=v.replace(" ", "_")
			v=v.replace(",_", " ")
			v=v.replace(',', '')
			v=(v.split())
			v=set(v)
			v=list(v)

			if "ampicillin" in set(v):
				predp_dict.setdefault(sample,[]).append("ampicillin")
			if "ampicillin" not in set(v):	
				predsus_dict.setdefault(sample,[]).append("ampicillin")
			if "amoxicillin-clavulanic_acid" in set(v):
				predp_dict.setdefault(sample,[]).append("amoxicillin-clavulanic_acid")
			if "amoxicillin-clavulanic_acid" not in set(v):	
				predsus_dict.setdefault(sample,[]).append("amoxicillin-clavulanic_acid")				
			if "amikacin" in set(v):
				predp_dict.setdefault(sample,[]).append("amikacin")
			if "amikacin" not in set(v):	
				predsus_dict.setdefault(sample,[]).append("amikacin")				
			if "cefazolin" in set(v):
				predp_dict.setdefault(sample,[]).append("cefazolin")
			if "cefazolin" not in set(v):	
				predsus_dict.setdefault(sample,[]).append("cefazolin")				
			if "cefalotin" in set(v):
				predp_dict.setdefault(sample,[]).append("cefalotin")
			if "cefalotin" not in set(v):	
				predsus_dict.setdefault(sample,[]).append("cefalotin")	
			if "ciprofloxacin" in set(v):
				predp_dict.setdefault(sample,[]).append("ciprofloxacin")
			if "ciprofloxacin" not in set(v):	
				predsus_dict.setdefault(sample,[]).append("ciprofloxacin")					
			if "cefixime" in set(v):
				predp_dict.setdefault(sample,[]).append("cefixime")
			if "cefixime" not in set(v):	
				predsus_dict.setdefault(sample,[]).append("cefixime")					
			if "ceftazidime" in set(v):
				predp_dict.setdefault(sample,[]).append("ceftazidime")
			if "ceftazidime" not in set(v):	
				predsus_dict.setdefault(sample,[]).append("ceftazidime")					
			if "gentamicin_C" in set(v):
				predp_dict.setdefault(sample,[]).append("gentamicin_C")
			if "gentamicin_C" not in set(v):
				predsus_dict.setdefault(sample,[]).append("gentamicin_C")				
			if "meropenem" in set(v):
				predp_dict.setdefault(sample,[]).append("meropenem")
			if "meropenem" not in set(v):
				predsus_dict.setdefault(sample,[]).append("meropenem")					
			if "nitrofurantoin" in set(v):
				predp_dict.setdefault(sample,[]).append("nitrofurantoin")
			if "nitrofurantoin" not in set(v):
				predsus_dict.setdefault(sample,[]).append("nitrofurantoin")					
			if "piperacillin-tazobactam" in set(v):
				predp_dict.setdefault(sample,[]).append("piperacillin-tazobactam")
			if "piperacillin-tazobactam" not in set(v):
				predsus_dict.setdefault(sample,[]).append("piperacillin-tazobactam")	
			if "tetracycline" in set(v):
				predp_dict.setdefault(sample,[]).append("tetracycline")
			if "tetracycline" not in set(v):
				predsus_dict.setdefault(sample,[]).append("tetracycline")				
			if "tobramycin" in set(v):
				predp_dict.setdefault(sample,[]).append("tobramycin")
			if "tobramycin" not in set(v):
				predsus_dict.setdefault(sample,[]).append("tobramycin")				
			if "trimethoprim-sulfamethoxazole" in set(v):
				for x in (predp_dict.values()):
					if "trimethoprim-sulfamethoxazole" not in x:
						predp_dict.setdefault(sample,[]).append("trimethoprim-sulfamethoxazole")
			if "trimethoprim-sulfamethoxazole" not in set(v):
				if "trimethoprim" not in set(v):
					for x in (predsus_dict.values()):
						if "trimethoprim-sulfamethoxazole" not in x:
					
							predsus_dict.setdefault(sample,[]).append("trimethoprim-sulfamethoxazole")
			if "trimethoprim-sulfamethoxazole" not in set(v):
				if "sulfamethoxazole" not in set(v):
					for x in (predsus_dict.values()):
						if "trimethoprim-sulfamethoxazole" not in x:
							predsus_dict.setdefault(sample,[]).append("trimethoprim-sulfamethoxazole")

			if "trimethoprim" in set(v): 
				if "sulfamethoxazole" in set(v):
					for x in (predp_dict.values()):
						if "trimethoprim-sulfamethoxazole" not in x:
							predp_dict.setdefault(sample,[]).append("trimethoprim-sulfamethoxazole")
			if "cefoxitin" in set(v):
				predp_dict.setdefault(sample,[]).append("cefoxitin")
			if "cefoxitin" not in set(v):
				predsus_dict.setdefault(sample,[]).append("cefoxitin")					
			if "ceftriaxone" in set(v):
				predp_dict.setdefault(sample,[]).append("ceftriaxone")
			if "ceftriaxone" not in set(v):
				predsus_dict.setdefault(sample,[]).append("ceftriaxone")				
			if "ertapenem" in set(v):
				predp_dict.setdefault(sample,[]).append("ertapenem")
			if "ertapenem" not in set(v):
				predsus_dict.setdefault(sample,[]).append("ertapenem")

		correct=[]
		correctsus=[]
		overpred=[]
		underpred=[]
		matchdict={}

		# True postive
		for k,v in predp_dict.items():
			for k1,v1 in testp_dict.items():
				if k==k1:
					for x1 in set(v):
						if x1 in set(v1):
							correct.append(x1)

		# True negative
		for k2,v2 in predsus_dict.items():
			for k3,v3 in testsus_dict.items():
				if k2==k3:
					for x2 in set(v2):
						if x2 in set(v3):
							correctsus.append(x2)	

		# False positive
		for k4,v4 in predp_dict.items():
			for k5,v5 in testsus_dict.items():
				if k4==k5:
					for x3 in set(v4):
						if x3 in set(v5):
							overpred.append(x3)

		# False negative
		for k6,v6 in predsus_dict.items():
			for k7,v7 in testp_dict.items():
				if k6==k7:
					for x4 in set(v6):
						if x4 in set(v7):
							underpred.append(x4)
		matchdict[(fr.split(".")[0])]=[k, ", ".join(map(str, v1)), ", ".join(map(str, v)), ", ".join(map(str, correct)), ", ".join(map(str, correctsus)), ", ".join(map(str, overpred)), ", ".join(map(str, underpred))]
		for key, value in matchdict.items():
			writer.writerow(value)

