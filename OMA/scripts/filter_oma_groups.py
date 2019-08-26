import sys

SPECIES = list()
oma_groups_path = str()

if len(sys.argv) > 1:
    oma_groups_path = sys.argv[1]

if len(sys.argv) > 2:
    SPECIES = [x.upper() for x in sys.argv[1:]]

with open(oma_groups_path, 'rt') as oma_groups:
    for i in range(3): next(oma_groups)
    
    for line in oma_groups:
        line = line.strip().split()
        oma_group_id = line[0]
        oma_fingerprint = line[1]
        oma_entries = line[1:]
        valid_entries = list()
        for entry in oma_entries:
            for species in SPECIES:
                if species in entry:
                    valid_entries.append(entry)
                    print(valid_entries)
                    
        
        


    