import json
import csv
import sys
import os
import re
import operator

def checkKeyExisted(key, my_dict):
	try:
		nonNone = my_dict[key] is not None
	except KeyError:
		nonNone = False
	return nonNone

filepath = os.path.dirname(os.path.realpath(__file__))

def run():
	for file in os.listdir(filepath):
		if file.endswith("effluxpumpmetamodels.json"):
			tsvfilename = os.path.splitext(file)[0]
			with open(filepath+"/"+tsvfilename+".txt", "w") as af:
				writer = csv.writer(af, delimiter='\t', dialect='excel')
				writer.writerow(["META-MODEL_ID", "META-MODEL_NAME", "CATEGORY", "MODEL_ID", "CUT_OFF", "BEST_HIT_BITSCORE", "BITSCORE", "PERCENT_IDENTITY", "ARO_NAME", "snp", "ARO_CATEGORY", "ARO_ACCESSION", "MODEL_TYPE", "Predicted_DNA","Predicted_Protein","CARD_Protein_Sequence","LABEL"])

				if file.endswith("effluxpumpmetamodels.json"):
					with open(file) as epp_file:
						epp_data = json.load(epp_file)

					for entry in epp_data:
						for category in epp_data[entry]:
							for modelid in epp_data[entry][category]:
								order = {}
								dna = 0
								cgList = []
								bestAROcategory = []
								AROnameList = []
								bestAROcategorydict = {}
								hitID = []
								minARO = []
								minevalue = []	

								for info_modelid in epp_data[entry][category][modelid]:
									if "orf_dna_sequence" in epp_data[entry][category][modelid]:
										dna = 1
									if checkKeyExisted("ARO_category",epp_data[entry][category][modelid]):
										for aroctkey in epp_data[entry][category][modelid]["ARO_category"]:
											cgList.append(str(epp_data[entry][category][modelid]["ARO_category"][aroctkey]["category_aro_name"].encode('ascii','replace')))
									if checkKeyExisted("ARO_category", epp_data[entry][category][modelid]):
										for key in epp_data[entry][category][modelid]["ARO_category"]:
											bestAROcategory.append(str(epp_data[entry][category][modelid]["ARO_category"][key]["category_aro_name"].encode('ascii','replace')))
										bestAROcategorydict[str(minARO)+" "+str(minevalue)] = bestAROcategory
									if "orf_from:" in info_modelid:
										hitID.append(epp_data[entry][category][modelid]["orf_from"])

								match_dict = {}

								if dna == 1:
									if len(epp_data[entry][category][modelid]) != 0:
										if epp_data[entry][category][modelid]["model_type_id"] == 41091:
											if checkKeyExisted("snp", epp_data[entry][category][modelid]):
												temp = epp_data[entry][category][modelid]["snp"]["original"] + str(epp_data[entry][category][modelid]["snp"]["position"]) + epp_data[entry][category][modelid]["snp"]["change"]
											else:
												temp = "n/a"
										elif epp_data[entry][category][modelid]["model_type_id"] == 40293:
											if checkKeyExisted("snp", epp_data[entry][category][modelid]):
												temp = epp_data[entry][category][modelid]["snp"]["original"] + str(epp_data[entry][category][modelid]["snp"]["position"]) + epp_data[entry][category][modelid]["snp"]["change"]
											else:
												temp = "n/a"
										elif epp_data[entry][category][modelid]["model_type_id"] == 40292:
											temp = "n/a"
										match_dict[entry] = [re.sub("\D", "", entry), 
										entry.split(',')[1].lstrip(' '),
										category,
										re.sub("model_id: ", "", modelid),
										epp_data[entry][category][modelid]["type_match"],
										epp_data[entry][category][modelid]["bit_score"],
										epp_data[entry][category][modelid]["pass_bitscore"],
										epp_data[entry][category][modelid]["perc_identity"],
										epp_data[entry][category][modelid]["ARO_name"],
										temp,
										"; ".join(epp_data[entry][category][modelid]["ARO_category"][x]["category_aro_name"] for x in epp_data[entry][category][modelid]["ARO_category"]),
										epp_data[entry][category][modelid]["ARO_accession"],	
										epp_data[entry][category][modelid]["model_type"],							
										epp_data[entry][category][modelid]["orf_dna_sequence"],
										epp_data[entry][category][modelid]["orf_prot_sequence"],
										epp_data[entry][category][modelid]["sequence_from_broadstreet"],
										epp_data[entry][category][modelid]["orf_from"]
										]
									for key, value in match_dict.items():						
										writer.writerow(value)

								else:
									if len(epp_data[entry][category][modelid]) != 0:
										if epp_data[entry][category][modelid]["model_type_id"] == 41091:
											if checkKeyExisted("snp", epp_data[entry][category][modelid]):
												temp = epp_data[entry][category][modelid]["snp"]["original"] + str(epp_data[entry][category][modelid]["snp"]["position"]) + epp_data[entry][category][modelid]["snp"]["change"]
											else:
												temp = "n/a"
										elif epp_data[entry][category][modelid]["model_type_id"] == 40293:
											if checkKeyExisted("snp", epp_data[entry][category][modelid]):
												temp = epp_data[entry][category][modelid]["snp"]["original"] + str(epp_data[entry][category][modelid]["snp"]["position"]) + epp_data[entry][category][modelid]["snp"]["change"]
											else:
												temp = "n/a"
										elif epp_data[entry][category][modelid]["model_type_id"] == 40292:
											temp = "n/a"

										match_dict[entry] = [re.sub("\D", "", entry),
										entry.split(',')[1].lstrip(' '), 
										category,
										re.sub("model_id: ", "", modelid),
										epp_data[entry][category][modelid]["type_match"],
										epp_data[entry][category][modelid]["bit_score"],
										epp_data[entry][category][modelid]["pass_bitscore"],
										epp_data[entry][category][modelid]["ARO_name"],
										temp,
										"; ".join(epp_data[entry][category][modelid]["ARO_category"][x]["category_aro_name"] for x in epp_data[entry][category][modelid]["ARO_category"]),
										epp_data[entry][category][modelid]["ARO_accession"],		
										epp_data[entry][category][modelid]["model_type"],									
										" ",
										epp_data[entry][category][modelid]["orf_prot_sequence"],
										epp_data[entry][category][modelid]["sequence_from_broadstreet"],
										epp_data[entry][category][modelid]["orf_from"]
										]
									
									for key, value in match_dict.items():
										writer.writerow(value)
					
			af.close()

def manual():
	h = {}
	h["META-MODEL_ID"] = "Meta-model found in submitted sequence"
	h["CATEGORY"] = "Efflux Pump Predictor Detection Paradigm for meta-models"
	h["MODEL_ID"] = "Model found in submitted sequences that corresponds to a meta-model"
	h["CUT_OFF"] = "RGI Detection Paradigm for models"
	h["BEST_HIT_BITSCORE"] = "Bit score of detected model"
	h["BIT_SCORE"] = "STRICT detection model - Curated bitscore cutoff for specific model in CARD"
	h["PERCENT_IDENTITY"] = "Quantitative measurement of the similarity between predicted and CARD protein sequence"
	h["ARO_name"] = "ARO term of detected model"
	h["snp"] = "Observed mutation (if applicable)"
	h["ARO_category"] = "ARO Categorization of detected model"
	h["ARO_accession"] = "ARO Accession number of detected model"
	h["Model_type"] = "CARD detection model type"	
	h["Predicted_DNA"] = "ORF predicted DNA sequence"
	h["Predicted_Protein"] = "ORF predicted protein sequence"
	h["CARD_Protein_Sequence"] = "Protein sequence of detected model in CARD"
	h["LABEL"] = "ORF label (internal to RGI)"
	print("\n")
	print("COLUMN","\t\t\t","HELP_MESSAGE")
	for i in h:
		print(i,"\t\t\t",h[i])	
	print("\n")

if __name__ == '__main__':
	run()

