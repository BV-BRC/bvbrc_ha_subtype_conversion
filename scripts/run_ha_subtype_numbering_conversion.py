#!/usr/bin/env python

import argparse
import csv
import json
import os
import requests
import subprocess
import sys
import urllib.parse
from Bio import SearchIO


def resolve_reference_path(filename):
    top = os.getenv("KB_TOP") or ""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(top, "lib", filename),
        os.path.join(top, "modules", "bvbrc_ha_subtype_conversion", "lib", filename),
        os.path.join(script_dir, "..", "lib", filename),
        os.path.join("/home", "ac.mkuscuog", "git", "bvbrc_ha_subtype_conversion", "lib", filename),
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    return candidates[-1]


REFERENCE_SEQUENCE = resolve_reference_path("Burke_HA_Reference_Sequences.fasta")
REFERENCE_ALIGNMENT_SEQUENCE = resolve_reference_path("Burke_HA_Reference_Sequences_alignment.fasta")

BLAST_SEQ_ANNOTATION_NAME = "sequence_annotation.tsv"
BLAST_SEQ_ANNOTATION_FILE_HEADER = ["QueryId", "Virus Type", "Segment", "Subtype", "Bitscore", "E-value",
                                    "Warning Messages", "Sequence Name"]
BLAST_WARNING_MESSAGE_NO_MATCH = "No Similar HA sequences were identified. Check to make sure that the query sequence provided is from an influenza HA protein."
BLAST_WARNING_MESSAGE_NONE = "None"
BLAST_THRESHOLD = 1e-1

API_BASE_PATH = "https://www.bv-brc.org/api/"
API_GENOME_FEATURE_SELECT = API_BASE_PATH + "genome_feature/?&select(feature_id)&sort(+feature_id)&in(feature_id,FeatureGroup(%s))"
API_GENOME_FEATURE_SELECT_LIST = API_BASE_PATH + "genome_feature/?&select(feature_id)&sort(+feature_id)&in(patric_id,(%s))"
API_GENOME_FEATURE_DOWNLOAD = API_BASE_PATH + "genome_feature/?&in(feature_id,(%s))&http_accept=application/protein+fasta"


def get_authorized_session():
    session = requests.Session()
    if "KB_AUTH_TOKEN" in os.environ:
        print("Reading auth key from environment")
        session.headers.update({'Authorization': os.environ['KB_AUTH_TOKEN']})
        return session
    print("Reading auth key from file")
    token_file = os.path.join(os.environ.get('HOME', ''), ".patric_token")
    if os.path.exists(token_file):
        with open(token_file) as token_fh:
            session.headers.update({'Authorization': token_fh.read().rstrip()})
        return session
    return None


def fetch_feature_fasta(session, select_api):
    print(f"Requesting feature ids: {select_api}")
    response = session.get(select_api)
    feature_ids = [data["feature_id"] for data in response.json()]

    download_api = API_GENOME_FEATURE_DOWNLOAD % (",".join(feature_ids))
    print(f"Requesting fasta data: {download_api}")
    return session.get(download_api).text


def create_fasta_file(output_dir, job_data):
    input_file = os.path.join(output_dir, "input.fasta")
    source = job_data["input_source"]

    if source == "fasta_file":
        try:
            fetch_fasta_cmd = ["p3-cp", f"ws:{job_data['input_fasta_file']}", input_file]
            subprocess.check_call(fetch_fasta_cmd, shell=False)
        except Exception as e:
            print(f"Error copying fasta file from workspace:\n {e}")
            sys.exit(1)
    elif source == "fasta_data":
        try:
            with open(input_file, "w+") as f:
                f.write(job_data["input_fasta_data"])
        except Exception as e:
            print(f"Error copying fasta data to input file:\n {e}")
            sys.exit(1)
    elif source in ("feature_group", "feature_list"):
        try:
            session = get_authorized_session()
            if session is None:
                print("Error authorizing the session for api call")
                sys.exit(1)

            if source == "feature_group":
                select_api = API_GENOME_FEATURE_SELECT % (
                    urllib.parse.quote(job_data["input_feature_group"], safe=""))
            else:
                select_api = API_GENOME_FEATURE_SELECT_LIST % (','.join(job_data["input_feature_list"]))

            fasta_text = fetch_feature_fasta(session, select_api)
            with open(input_file, "w+") as f:
                f.write(fasta_text)
        except Exception as e:
            print(f"Error retrieving data from feature group:\n {e}")
            sys.exit(1)

    return input_file


def parse_fasta_file(fasta_data):
    header = None
    for line in fasta_data:
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


# Generate gap counts after residues
# Example: 3->2 means there are 2 gaps (-) starting from position 3 in the sequence
def get_insertion_position_for_gaps(sequence):
    insertion_position_map = {}

    aa_counter = 0
    gap_counter = 0
    for ch in sequence:
        if ch == '-':
            gap_counter += 1
        else:
            if gap_counter > 0:
                insertion_position_map[aa_counter] = gap_counter
            aa_counter += 1
            gap_counter = 0
        if gap_counter > 0:
            insertion_position_map[aa_counter] = gap_counter

    return insertion_position_map


def get_position_for_aligned_coordinates(sequence):
    aligned_position_map = {}

    counter = 1
    for idx, ch in enumerate(sequence):
        if ch != "-":
            aligned_position_map[counter] = idx
            counter += 1

    return aligned_position_map


def get_insertion_indices(final_sequence, sequence):
    insertion_indices = []
    aligned_position_map = get_position_for_aligned_coordinates(final_sequence)
    insertion_position_map = get_insertion_position_for_gaps(sequence)

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


def get_reference_residue_map(aligned_query, aligned_reference):
    reference_residue_map = {}

    reference_position = 1
    for idx, ch in enumerate(aligned_reference):
        if ch != "-":
            reference_residue_map[reference_position] = aligned_query[idx]
            reference_position += 1

    return reference_residue_map


def get_final_aligned_query(aligned_query, aligned_reference, final_or_pre_aligned_reference, insertion_indices=None):
    reference_residue_map = get_reference_residue_map(aligned_query, aligned_reference)

    final_aligned_query = ["-"] * len(final_or_pre_aligned_reference)

    reference_position = 1
    for idx, ch in enumerate(final_or_pre_aligned_reference):
        if ch != "-":
            final_aligned_query[idx] = reference_residue_map[reference_position]
            reference_position += 1

    if insertion_indices is not None:
        insertion_aa = [aligned_query[idx] for idx, ch in enumerate(aligned_reference) if ch == "-"]
        for idx, aa in enumerate(insertion_aa):
            final_aligned_query[insertion_indices[idx]] = aa

    return "".join(final_aligned_query)


def main():
    parser = argparse.ArgumentParser(description="HA Subtype Numbering Conversion Script")
    parser.add_argument("-j", "--jfile", help="json file for job", required=True)
    parser.add_argument("-o", "--output", help="Output directory. defaults to current directory", required=False,
                        default=".")

    args = parser.parse_args()

    # Load job data
    job_data = None
    try:
        with open(args.jfile, "r") as j:
            job_data = json.load(j)
    except Exception as e:
        print(f"Error in opening job file:\n {e}")
        sys.exit(1)

    if not job_data:
        print("job_data is null")
        sys.exit(1)
    print(job_data)

    # Setup output directory
    output_dir = os.path.abspath(args.output)
    os.makedirs(output_dir, exist_ok=True)
    os.chdir(output_dir)

    selected_types = job_data["types"]

    # Create input file
    input_file = create_fasta_file(output_dir, job_data)

    if os.path.getsize(input_file) == 0:
        print("Input fasta file is empty")
        sys.exit(1)

    # Parse fasta file and create query fasta for blast
    blast_input_file = os.path.join(output_dir, "blast.fasta")
    sequences = {}
    with open(blast_input_file, "w") as blast_data:
        with open(input_file) as fasta_data:
            for index, values in enumerate(parse_fasta_file(fasta_data), start=1):
                query_name = f"query{index}"
                sequences[query_name] = values

                blast_data.write(f">{query_name}\n")
                blast_data.write(values["data"])
                blast_data.write("\n")

    # Parse reference sequences
    reference_sequences = {}
    with open(REFERENCE_ALIGNMENT_SEQUENCE) as fasta_data:
        for values in parse_fasta_file(fasta_data):
            subtype = values["header"].split("|")[2].split()[0]
            reference_sequences[subtype] = values["data"]

    blast_output_file = os.path.join(output_dir, "blast.out")
    # Run BLAST
    try:
        blast_cmd = ["blastall", "-d", REFERENCE_SEQUENCE, "-i", blast_input_file, "-p", "blastp", "-G", "10", "-E",
                     "1", "-m", "7", "-o", blast_output_file]
        subprocess.check_call(blast_cmd, shell=False)
    except Exception as e:
        print(f"Error running blast for {input_file}: {e}")
        sys.exit(1)

    # Parse blast result and create sequence annotation file
    qresults = SearchIO.parse(blast_output_file, "blast-xml")

    sequence_annotation_file_path = os.path.join(output_dir, BLAST_SEQ_ANNOTATION_NAME)
    with open(sequence_annotation_file_path, "w") as sequence_annotation_file:
        sequence_annotation_writer = csv.DictWriter(sequence_annotation_file, delimiter="\t",
                                                    fieldnames=BLAST_SEQ_ANNOTATION_FILE_HEADER)
        sequence_annotation_writer.writeheader()

        for qresult in qresults:
            print(f"Search {qresult.id} has {len(qresult)} hits")

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
                best_evalue = BLAST_THRESHOLD
                best_hsp = None
                for hit in qresult:
                    for hsp in hit:
                        if hsp.evalue < best_evalue:
                            best_evalue = hsp.evalue
                            best_hsp = hsp

                if best_hsp is None:
                    sequence_annotation_writer.writerow({"QueryId": query_id,
                                                         "Virus Type": "N/A",
                                                         "Segment": "N/A",
                                                         "Subtype": "N/A",
                                                         "Bitscore": "N/A",
                                                         "E-value": "N/A",
                                                         "Warning Messages": BLAST_WARNING_MESSAGE_NO_MATCH,
                                                         "Sequence Name": sequence_name})
                else:
                    hit_reference = best_hsp.hit_id.strip().split("|")
                    # Add hit result to sequence map
                    sequences[query_id]["subtype"] = hit_reference[2]
                    sequence_annotation_writer.writerow({"QueryId": query_id,
                                                         "Virus Type": hit_reference[0],
                                                         "Segment": hit_reference[1],
                                                         "Subtype": hit_reference[2],
                                                         "Bitscore": best_hsp.bitscore,
                                                         "E-value": best_hsp.evalue,
                                                         "Warning Messages": BLAST_WARNING_MESSAGE_NONE,
                                                         "Sequence Name": sequence_name})

    for query_name, value in sequences.items():
        if "subtype" not in value:
            print(f"No subtype found from BLAST: {query_name}")
            continue

        blast_subtype = value["subtype"]

        reference_sequence = reference_sequences[blast_subtype]
        sequence = value["data"]

        # Replace '-' with '' in the unaligned reference sequence
        reference_sequence = reference_sequence.replace("-", "")

        muscle_input_file = os.path.join(output_dir, query_name + ".muscle.in")
        muscle_output_file = os.path.join(output_dir, query_name + ".muscle.out")

        with open(muscle_input_file, "w") as muscle_data:
            muscle_data.write(f">{query_name}\n")
            muscle_data.write(sequence)
            muscle_data.write("\n")
            muscle_data.write(f">{blast_subtype}\n")
            muscle_data.write(reference_sequence)

        # Run muscle
        try:
            muscle_cmd = ["muscle", "-in", muscle_input_file, "-out", muscle_output_file, "-quiet"]
            subprocess.check_call(muscle_cmd, shell=False)
        except Exception as e:
            print(f"Error running muscle for {muscle_input_file}: {e}")
            sys.exit(1)

        muscle_sequence = {}
        with open(muscle_output_file) as muscle_data:
            for values in parse_fasta_file(muscle_data):
                muscle_sequence[values["header"]] = values["data"]

        # Combine subtypes from user selection and subtype from blast result
        selected_subtypes = [t for t in selected_types if t != blast_subtype]

        # Final result file
        result_file = os.path.join(output_dir, query_name + "_result.fasta")
        with open(result_file, "w") as result_file_writer:
            # Write query and best hit sequences to the final result file
            muscle_output_sequence = muscle_sequence[query_name].strip()
            aligned_reference_sequence = muscle_sequence[blast_subtype].strip()

            if "-" in aligned_reference_sequence:
                print("Dash exists in the aligned reference sequence")
                # Create final aligned reference sequence
                final_aligned_ref_seq = reference_sequences[blast_subtype]
                insertion_position_map = get_insertion_position_for_gaps(aligned_reference_sequence)
                for position, gap_count in insertion_position_map.items():
                    gap_str = "-" * gap_count
                    if position == 0:
                        final_aligned_ref_seq = gap_str + final_aligned_ref_seq
                    else:
                        aligned_position_map = get_position_for_aligned_coordinates(final_aligned_ref_seq)
                        aligned_position = aligned_position_map[position] + 1
                        final_aligned_ref_seq = (final_aligned_ref_seq[:aligned_position] + gap_str
                                                 + final_aligned_ref_seq[aligned_position:])

                insertion_indices = get_insertion_indices(final_aligned_ref_seq, aligned_reference_sequence)

                final_aligned_query = get_final_aligned_query(muscle_output_sequence, aligned_reference_sequence,
                                                              final_aligned_ref_seq, insertion_indices)
                result_file_writer.write(f">{query_name}\n")
                result_file_writer.write(final_aligned_query + "\n")
                result_file_writer.write(f">{blast_subtype}\n")
                result_file_writer.write(final_aligned_ref_seq + "\n")

                # Create aligned reference sequence with insertions for each subtype
                for subtype in selected_subtypes:
                    subtype_ref_sequence = reference_sequences[subtype]
                    aligned_ref_seq_with_insertions = []

                    counter = 0
                    for idx in range(len(final_aligned_ref_seq)):
                        if idx in insertion_indices:
                            aligned_ref_seq_with_insertions.append("-")
                        else:
                            aligned_ref_seq_with_insertions.append(subtype_ref_sequence[counter])
                            counter += 1

                    result_file_writer.write(f">{subtype}\n")
                    result_file_writer.write("".join(aligned_ref_seq_with_insertions) + "\n")
            else:
                print("Dash does not exist in the aligned reference sequence")

                final_aligned_query = get_final_aligned_query(muscle_output_sequence, aligned_reference_sequence,
                                                              reference_sequences[blast_subtype])
                result_file_writer.write(f">{query_name}\n")
                result_file_writer.write(final_aligned_query + "\n")
                result_file_writer.write(f">{blast_subtype}\n")
                result_file_writer.write(reference_sequences[blast_subtype] + "\n")

                for subtype in selected_subtypes:
                    result_file_writer.write(f">{subtype}\n")
                    result_file_writer.write(reference_sequences[subtype] + "\n")


if __name__ == "__main__":
    main()
