#!/usr/bin/env python

import argparse
import csv
import json
import os
import requests
import subprocess
import sys
import urllib
from Bio import SearchIO

#
#Determine paths.
#
top = os.getenv("KB_TOP")

reference_sequence_deployed = os.path.join(top, "lib", "Burke_HA_Reference_Sequences.fasta")
reference_sequence_dev = os.path.join(top, "modules", "bvbrc_ha_subtype_conversion", "lib", "Burke_HA_Reference_Sequences.fasta")
reference_sequence_local = os.path.join("/home", "ac.mkuscuog", "git", "bvbrc_ha_subtype_conversion", "lib", "Burke_HA_Reference_Sequences.fasta")
if os.path.exists(reference_sequence_deployed):
  REFERENCE_SEQUENCE = reference_sequence_deployed
elif os.path.exists(reference_sequence_dev):
  REFERENCE_SEQUENCE = reference_sequence_dev
else:
  REFERENCE_SEQUENCE = reference_sequence_local

reference_alignment_sequence_deployed = os.path.join(top, "lib", "Burke_HA_Reference_Sequences_alignment.fasta")
reference_alignment_sequence_dev = os.path.join(top, "modules", "bvbrc_ha_subtype_conversion", "lib", "Burke_HA_Reference_Sequences_alignment.fasta")
reference_alignment_sequence_local = os.path.join("/home", "ac.mkuscuog", "git", "bvbrc_ha_subtype_conversion", "lib", "Burke_HA_Reference_Sequences_alignment.fasta")
if os.path.exists(reference_alignment_sequence_deployed):
  REFERENCE_ALIGNMENT_SEQUENCE = reference_alignment_sequence_deployed
elif os.path.exists(reference_alignment_sequence_dev):
  REFERENCE_ALIGNMENT_SEQUENCE = reference_alignment_sequence_dev
else:
  REFERENCE_ALIGNMENT_SEQUENCE = reference_alignment_sequence_local

BLAST_SEQ_ANNOTATION_NAME = "sequence_annotation.tsv"
BLAST_SEQ_ANNOTATION_FILE_HEADER = ["QueryId", "Virus Type", "Segment", "Subtype", "Bitscore", "E-value", "Warning Messages", "Sequence Name"]
BLAST_WARNING_MESSAGE_NO_MATCH = "No Similar HA sequences were identified. Check to make sure that the query sequence provided is from an influenza HA protein."
BLAST_WARNING_MESSAGE_NONE = "None"
BLAST_THRESHOLD = 1e-1 

API_BASE_PATH = "https://www.bv-brc.org/api/"
API_GENOME_FEATURE_SELECT = API_BASE_PATH + "genome_feature/?&select(feature_id)&sort(+feature_id)&in(feature_id,FeatureGroup(%s))" 
API_GENOME_FEATURE_DOWNLOAD = API_BASE_PATH + "genome_feature/?&in(feature_id,(%s))&http_accept=application/protein+fasta" 
API_GENOME_FEATURE_DOWNLOAD_LIST = API_BASE_PATH + "genome_feature/?&in(bvbrc_id,(%s))&http_accept=application/protein+fasta"

HA_REFERENCE_TYPES = { 
  'H1_PR34': 'A/Puerto/Rico/8/34', 
  'H1_1933': 'A/United/Kingdom/1/1933',
  'H1post1995': 'A/NewCaledonia/20/1999',
  'H1N1pdm': 'A/California/04/2009', 
  'H2': 'A/Singapore/1/1957', 
  'H3': 'A/AICHI/2/68', 
  'H4': 'A/swine/Ontario/01911-1/99', 
  'H5mEA-nonGsGD': 'A/mallard/Italy/3401/2005 (LPAI)', 
  'H5': 'A/Vietnam/1203/04 (HPAI)', 
  'H5c221': 'A/chicken/Egypt/0915-NLQP/2009 (HPAI)', 
  'H6': 'A/chicken/Taiwan/0705/99', 
  'H7N3': 'A/Turkey/Italy/220158/02/H7N3', 
  'H7N7': 'A/Netherlands/219/03/H7N7', 
  'H8': 'A/turkey/Ontario/6118/1968', 
  'H9': 'A/Swine/HK/9/98', 
  'H10': 'A/mallard/bavaria/3/2006', 
  'H11': 'A/duck/England/1/1956', 
  'H12': 'A/Duck/Alberta/60/1976', 
  'H13': 'A/gull/Maryland/704/1977', 
  'H14': 'A/mallard/Astrakhan/263/1982', 
  'H15': 'A/duck/Australia/341/1983', 
  'H16': 'A/black-headedgull/Turkmenistan/13/76', 
  'H17': 'A/little-yellow-shoulderedbat/Guatemala/060/2010', 
  'H18': 'A/flat-faced/bat/Peru/033/2010', 
  'B/HONG KONG/8/73': 'B/HONGKONG/8/73', 
  'B/FLORIDA/4/2006': 'B/FLORIDA/4/2006', 
  'B/HUMAN/BRISBANE/60/2008': 'B/HUMAN/BRISBANE/60/2008', 
};

def createFASTAFile(output_dir, job_data):
  input_file = os.path.join(output_dir, "input.fasta")
  if job_data["input_source"] == "fasta_file":
    #Fetch input file from workspace
    try:
      fetch_fasta_cmd = ["p3-cp", "ws:%s" %(job_data["input_fasta_file"]), input_file]
      subprocess.check_call(fetch_fasta_cmd, shell=False)
    except Exception as e:
      print("Error copying fasta file from workspace:\n %s" %(e))
      sys.exit(-1)
  elif job_data["input_source"] == "fasta_data":
    #Copy user data to input file
    try:
      with open(input_file, "w+") as input:
        input.write(job_data["input_fasta_data"])
    except Exception as e:
      print("Error copying fasta data to input file:\n %s" %(e))
      sys.exit(-1)
  elif job_data["input_source"] == "feature_group":
    #Retrieve fasta data from feature group
    try:
      isAuthorized = False
      session = requests.Session();
      if "KB_AUTH_TOKEN" in os.environ:
        print("Reading auth key from environment")
        session.headers.update({ 'Authorization' : os.environ.get('KB_AUTH_TOKEN') })
        isAuthorized = True
      else:
        print("Reading auth key from file")
        tokenFile = os.path.join(os.environ.get('HOME'), ".patric_token")
        if os.path.exists(tokenFile):
          with open(tokenFile) as F:
            tokenString = F.read().rstrip()
            session.headers.update({ 'Authorization' : tokenString })
          isAuthorized = True

      if isAuthorized:      
        genome_select_api = API_GENOME_FEATURE_SELECT %(urllib.quote(job_data["input_feature_group"], safe=""))
        print("Requesting feature ids: %s" %(genome_select_api))
        response = session.get(genome_select_api)

        feature_ids = []
        for data in response.json():
          feature_ids.append(data["feature_id"])

        genome_download_api = API_GENOME_FEATURE_DOWNLOAD %(",".join(feature_ids)) 
        print("Requesting fasta data: %s" %(genome_download_api))
        response = session.get(genome_download_api)

        with open(input_file, "w+") as input:
          input.write(response.content)
      else:
        print("Error authorizing the session for api call")
        sys.exit(-1)
    except Exception as e:
      print("Error retrieving data from feature group:\n %s" %(e))
      sys.exit(-1)
  elif job_data["input_source"] == "feature_list":
    #Retrive fasta data from feature list
    try:
      isAuthorized = False
      session = requests.Session();
      if "KB_AUTH_TOKEN" in os.environ:
        print("Reading auth key from environment")
        session.headers.update({ 'Authorization' : os.environ.get('KB_AUTH_TOKEN') })
        isAuthorized = True
      else:
        print("Reading auth key from file")
        tokenFile = os.path.join(os.environ.get('HOME'), ".patric_token")
        if os.path.exists(tokenFile):
          with open(tokenFile) as F:
            tokenString = F.read().rstrip()
            session.headers.update({ 'Authorization' : tokenString })
          isAuthorized = True
      if isAuthorized:
        feature_ids = job_data["input_feature_list"] 

        genome_download_api = API_GENOME_FEATURE_DOWNLOAD_LIST %(",".join(feature_ids))
        print("Requesting feature list data: %s" %(genome_download_api))
        response = session.get(genome_download_api)

        with open(input_file, "w+") as input:
          input.write(response.content)
      else:
        print("Error authorizing the session for api call")
        sys.exit(-1)
    except Exception as e:
      print("Error retrieving data from feature group:\n %s" %(e))
      sys.exit(-1)

  return input_file

def parseFASTAFile(fasta_data):
  import pdb
  for line in fasta_data:
    pdb.set_trace()
    if line[0] == ">":
      header = line[1:].rstrip()
      break

  data = []
  for line in fasta_data:
    if line[0] == ">":
      yield {"header": header, "data": "".join(data).replace(" ", "").replace("\r", "")}
      data = []
      header = line[1:].rstrip()
      continue
    data.append(line.rstrip())

  yield {"header": header, "data": "".join(data).replace(" ", "").replace("\r", "")}

#Generate gap counts after residues
#Example: 3->2 means there are 2 gaps (-) starting from position 3 in the sequence
def getInsertionPositionForGaps(sequence):
  insertion_position_map = {}
  
  aa_counter = 0
  gap_counter = 0
  for ch in aligned_reference_sequence:
    if ch == '-':
      gap_counter = gap_counter + 1
    else:
      if gap_counter > 0:
        insertion_position_map[aa_counter] = gap_counter
      aa_counter = aa_counter + 1
      gap_counter = 0
    if gap_counter > 0:
      insertion_position_map[aa_counter] = gap_counter

  return insertion_position_map

def getPositionForAlignedCoordinates(sequence):
  aligned_position_map = {}

  counter = 1
  for idx in range(0, len(sequence)):
    if sequence[idx] != "-":
      aligned_position_map[counter] = idx
      counter = counter + 1

  return aligned_position_map

def getInsertionIndices(final_sequence, sequence):
  insertion_indices = []
  aligned_position_map = getPositionForAlignedCoordinates(final_sequence)
  insertion_position_map = getInsertionPositionForGaps(sequence)

  clean_final_seq_len = len(final_sequence.replace("-", ""))
  for position, gap_count in insertion_position_map.items():
    if position == 0:
      for idx in range(gap_count):
        insertion_indices.append(idx)    
    else:
      start = len(final_sequence) if position == clean_final_seq_len else aligned_position_map[position + 1] 
      for idx in range(1, gap_count + 1):
        insertion_indices.append(start - idx)

  return insertion_indices

def getReferenceResidueMap(aligned_query, aligned_reference):
  reference_residue_map = {}

  reference_position = 1
  for idx in range(len(aligned_reference)):
    if aligned_reference[idx] != "-":
      reference_residue_map[reference_position] = aligned_query[idx]
      reference_position = reference_position + 1

  return reference_residue_map

def getFinalAlignedQuery(aligned_query, aligned_reference, final_or_pre_aligned_reference, insertion_indices = None):
  final_aligned_query = []

  reference_residue_map = getReferenceResidueMap(aligned_query, aligned_reference) 

  for _ in range(len(final_or_pre_aligned_reference)):
    final_aligned_query.append("-")

  reference_position = 1 
  for idx in range(len(final_or_pre_aligned_reference)):
    if final_or_pre_aligned_reference[idx] != "-":
      final_aligned_query[idx] = reference_residue_map[reference_position]
      reference_position = reference_position + 1

  if insertion_indices is not None:
    insertion_aa = []
    for idx in range(len(aligned_reference)):
      if aligned_reference[idx] == "-":
        insertion_aa.append(aligned_query[idx])

    for idx in range(len(insertion_aa)):
      index = insertion_indices[idx]
      final_aligned_query[index] = insertion_aa[idx]

  return "".join(final_aligned_query)

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="HA Subtype Numbering Conversion Script")
  parser.add_argument("-j", "--jfile", help="json file for job", required=True)
  parser.add_argument("-o", "--output", help="Output directory. defaults to current directory", required=False, default=".")

  args = parser.parse_args()

  #Load job data
  job_data = None
  try:
    with open(args.jfile, "r") as j:
      job_data = json.load(j)
  except Exception as e:
    print("Error in opening job file:\n %s" %(e))
    sys.exit(-1)

  if not job_data:
    print("job_data is null")
    sys.exit(-1)
  print(job_data)

  #Setup output directory
  output_dir = args.output
  output_dir = os.path.abspath(output_dir)
  if not os.path.exists(output_dir):
    os.mkdir(output_dir)
  os.chdir(output_dir)

  output_file = os.path.join(output_dir, job_data["output_file"] + ".txt")

  selected_types = job_data["types"]

  #Create input file
  input_file = createFASTAFile(output_dir, job_data)

  if os.path.getsize(input_file) == 0:
    print("Input fasta file is empty")
    sys.exit(-1)

  #Parse fasta file and create query fasta for blast
  blast_input_file = os.path.join(output_dir, "blast.fasta")
  index = 1
  sequences = {} 
  with open(blast_input_file, "w") as blast_data:
    with open(input_file) as fasta_data:
      for values in parseFASTAFile(fasta_data):
        query_name = "query" + str(index)
        sequences[query_name] = values

        blast_data.write(">" + query_name)
        blast_data.write("\n")
        blast_data.write(values["data"])
        blast_data.write("\n")

        index = index + 1

  #Parse reference sequences
  reference_sequences = {}
  with open(REFERENCE_ALIGNMENT_SEQUENCE) as fasta_data:
    for values in parseFASTAFile(fasta_data):
      subtype = values["header"].split("|")[2]
      reference_sequences[subtype] = values["data"]

  blast_output_file = os.path.join(output_dir, "blast.out")
  #Run BLAST
  try:
    blast_cmd = ["blastall", "-d", REFERENCE_SEQUENCE, "-i", blast_input_file, "-p", "blastp", "-G", "10", "-E", "1", "-o", blast_output_file]
    subprocess.check_call(blast_cmd, shell=False)
  except Exception as e:
    print("Error running blast for %s" %(input_file))
    sys.exit(-1)
  
  #Parse blast result and create sequence annotation file
  qresults = SearchIO.parse(blast_output_file, "blast-text")
  
  sequence_annotation_file_path = os.path.join(output_dir, BLAST_SEQ_ANNOTATION_NAME)
  sequence_annotation_file = open(sequence_annotation_file_path, "w")
  sequence_annotation_writer = csv.DictWriter(sequence_annotation_file, delimiter="\t", fieldnames=BLAST_SEQ_ANNOTATION_FILE_HEADER)
  sequence_annotation_writer.writeheader()
  
  for qresult in qresults:
    print("Search %s has %i hits" % (qresult.id, len(qresult)))
    
    query_id = qresult.id
    sequence_name = sequences[query_id]["header"] 
    if len(qresult) == 0:
      sequence_annotation_writer.writerow({"QueryId": query_id, 
                                           "Virus Type": "N/A", 
                                           "Segment": "N/A", 
                                           "Subtype": "N/A", 
                                           "Bitscore": "N/A",
                                           "E-value": "N/A",
                                           "Warning Messages": BLAST_WARNING_MESSAGE_NO_MATCH, 
                                           "Sequence Name": sequence_name})
    else:
      bestEvalue = BLAST_THRESHOLD
      bestHSP = None
      for hit in qresult: 
        for hsp in hit:
          if hsp.evalue < bestEvalue:
            bestEvalue = hsp.evalue
            bestHSP = hsp  

      if bestHSP is None:
        sequence_annotation_writer.writerow({"QueryId": query_id,
                                             "Virus Type": "N/A",
                                             "Segment": "N/A",
                                             "Subtype": "N/A",
                                             "Bitscore": "N/A",
                                             "E-value": "N/A",
                                             "Warning Messages": BLAST_WARNING_MESSAGE_NO_MATCH,
                                             "Sequence Name": sequence_name})
      else:
        hitReference = bestHSP.hit_id.strip().split("|")
        #Add hit result to sequence map
        sequences[query_id]["subtype"] = hitReference[2]
        sequence_annotation_writer.writerow({"QueryId": query_id,
                                              "Virus Type": hitReference[0],
                                              "Segment": hitReference[1],
                                              "Subtype": hitReference[2],
                                              "Bitscore": bestHSP.bitscore,
                                              "E-value": bestHSP.evalue,
                                              "Warning Messages": BLAST_WARNING_MESSAGE_NONE,
                                              "Sequence Name": sequence_name})
      #print(vars(qresult))
      #print(qresult.description)

  #Close sequence annotation file
  sequence_annotation_file.close()

  for query_name, value in sequences.items(): 
    if "subtype" not in value:
      print("No subtype found from BLAST: %s" %(query_name))
      break

    blast_subtype = value["subtype"]

    reference_sequence = reference_sequences[blast_subtype]
    sequence = value["data"]

    #Replace '-' with '' in the unaligned reference sequence
    reference_sequence = reference_sequence.replace("-", "")

    muscle_input_file = os.path.join(output_dir, query_name + ".muscle.in")
    muscle_output_file = os.path.join(output_dir, query_name + ".muscle.out")

    with open(muscle_input_file, "w") as muscle_data:
      muscle_data.write(">" + query_name)
      muscle_data.write("\n")
      muscle_data.write(sequence)
      muscle_data.write("\n")
      muscle_data.write(">" + blast_subtype)
      muscle_data.write("\n")
      muscle_data.write(reference_sequence)

    #Run muscle
    try:
      muscle_cmd = ["muscle", "-in", muscle_input_file, "-out", muscle_output_file, "-quiet"]
      subprocess.check_call(muscle_cmd, shell=False)
    except Exception as e:
      print("Error running muscle for %s" %(muscle_input_file))
      sys.exit(-1)

    muscle_sequence = {}
    with open(muscle_output_file) as muscle_data:
      for values in parseFASTAFile(muscle_data):
        muscle_sequence[values["header"]] = values["data"]

    #Combine subtypes from user selection and subtype from blast result
    selected_subtypes = [] 
    for type in list(selected_types):
      if type != blast_subtype:
        selected_subtypes.append(type)

    #Final result file 
    result_file = os.path.join(output_dir, query_name + "_result.fasta")
    result_file_writer = open(result_file, "w")

    #Write query and best hit sequences to the final result file
    muscle_output_sequence = muscle_sequence[query_name].strip()
    aligned_reference_sequence = muscle_sequence[blast_subtype].strip()

    if "-" in aligned_reference_sequence:
      print("Dash exists in the aligned reference sequence")
      #Create final aligned reference sequence
      final_aligned_ref_seq = reference_sequences[blast_subtype]
      insertion_position_map = getInsertionPositionForGaps(aligned_reference_sequence)
      for position, gap_count in insertion_position_map.items():
        gap_str = "".join("-" for _ in range(gap_count)) 
        if position == 0:
          final_aligned_ref_seq = gap_str + final_aligned_ref_seq[0:]
        else:
          aligned_position_map = getPositionForAlignedCoordinates(final_aligned_ref_seq)
          aligned_position = aligned_position_map[position] + 1
          final_aligned_ref_seq = final_aligned_ref_seq[:aligned_position] + gap_str + final_aligned_ref_seq[aligned_position:]

      #insertionIndices = [i for i, ltr in enumerate(aligned_reference_sequence) if ltr == "-"]
      #insertion_indices = getInsertionIndices(final_aligned_ref_seq, reference_sequences[blast_subtype])
      insertion_indices = getInsertionIndices(final_aligned_ref_seq, aligned_reference_sequence)

      final_aligned_query = getFinalAlignedQuery(muscle_output_sequence, aligned_reference_sequence, final_aligned_ref_seq, insertion_indices)
      result_file_writer.write(">" + query_name + "\n")
      result_file_writer.write(final_aligned_query + "\n")
      result_file_writer.write(">" + blast_subtype + "\n")
      result_file_writer.write(final_aligned_ref_seq + "\n")

      #Create aligned reference sequence with insertions for each subtype
      for subtype in selected_subtypes:
        subtype_ref_sequence = reference_sequences[subtype] 
        aligned_ref_seq_with_insertions = []

        counter = 0
        for idx in range(len(final_aligned_ref_seq)):
          if idx in insertion_indices:
            aligned_ref_seq_with_insertions.append("-")
          else:
            aligned_ref_seq_with_insertions.append(subtype_ref_sequence[counter])
            counter = counter + 1

        result_file_writer.write(">" + subtype + "\n")
        result_file_writer.write("".join(aligned_ref_seq_with_insertions) + "\n")
    else: 
      print("Dash does not exist in the aligned reference sequence")
      
      final_aligned_query = getFinalAlignedQuery(muscle_output_sequence, aligned_reference_sequence, reference_sequences[blast_subtype])
      result_file_writer.write(">" + query_name + "\n")
      result_file_writer.write(final_aligned_query + "\n")
      result_file_writer.write(">" + blast_subtype + "\n")
      result_file_writer.write(reference_sequences[blast_subtype] + "\n")

      for subtype in selected_subtypes:
        result_file_writer.write(">" + subtype + "\n")
        result_file_writer.write(reference_sequences[subtype] + "\n")

    #Close result file
    result_file_writer.close()
