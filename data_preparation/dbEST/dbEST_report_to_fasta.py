import textwrap
import sys
import os
import gzip

if len(sys.argv) < 2:
    exit("run: python dbEST_report_to_fasta.py <dbEST_report> <filter_keyword:optional>")
else:
    full_path = sys.argv[1]
    file_name = os.path.basename(full_path).replace(".gz","")

FILTER = False

if len(sys.argv) > 2:
    FILTER = " ".join(sys.argv[2:])

def parse(ready_line):

    result = {
        "dbest_id" : "NULL",
        "est_name" : "NULL",
        "genbank_acc": "NULL",
        "tissue_type":"NULL",
        "organism" : "NULL", 
        }

    _idx=0
    splitted = ready_line.strip().split("\n")

    for line in splitted:
        _idx += 1
        if "dbEST Id" in line:
            result["dbest_id"] = splitted[2].replace(' ','').split(":")[-1]
        elif "EST name" in line:
            result["est_name"] = splitted[3].replace(' ','').split(":")[-1]
        elif "GenBank Acc" in line:
            result["genbank_acc"] = splitted[4].replace(' ','').split(":")[-1]
            break    
    
    for i in range(len(splitted)):
        if "SEQUENCE" in splitted[i]:
            _idx = i + 1
            break


    sequence = ""
    dna_set = set("ACGT")
    valid_sequence = True
    idx = _idx

    while True:
        __line =  splitted[idx].replace(' ','')
        if len(set(__line).intersection(dna_set)) >= 3:
            sequence += __line
            idx += 1
        else:
            idx += 1
            break
    
    
    # Assert there is no "ACGT" characters
    for i in set(sequence):
        if i not in ["A","C","G","T"]:
            return False
        

    seq = textwrap.fill(sequence, width=60)

    for line in splitted[idx:]:
        if "Organism" in line:
            result["organism"] = line.split(":")[-1].strip()
        elif "Tissue" in line:
            result["tissue_type"] = line.split(":")[-1].strip()
            break
        else:
            continue

    header = [f"{key}:{value}|" for key, value in result.items()]
    header = "".join(header)
    header = header[:-1]
    
    return [header, seq]
    

fasta_output = f"{file_name}.fa"

if FILTER:
    FILTER_NAME = FILTER.replace(" ","-")
    fasta_output = f"filtered_{FILTER_NAME}_" + fasta_output


with gzip.open(full_path, 'rt', errors='ignore') as report, open(fasta_output, 'w') as fasta, open(fasta_output + ".names", 'w') as names:
    ready_line = ""

    for line in report:
        curr_line = line.strip()
        
        while curr_line.strip() != "||":
            ready_line += curr_line
            curr_line = next(report)


        record = parse(ready_line)
        ready_line = ""
        
        # fasta.write(record[0])
        # fasta.write(record[1] + "\n\n")

        if record:
            if FILTER:
                if FILTER.lower() in record[0].lower():
                    fasta.write(">" + record[0] + "\n")
                    fasta.write(record[1] + "\n\n")
                    names.write(record[0]+"\t"+record[0]+"\n")
            else:
                fasta.write(">" + record[0] + "\n")
                fasta.write(record[1] + "\n\n")
                names.write(record[0][1:]+"\t"+record[0][1:]+"\n")