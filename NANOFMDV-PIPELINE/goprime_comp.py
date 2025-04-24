import sys
from Bio import SeqIO


def get_filename_stub(this_name, this_ext):
    this_stub = this_name

    this_pos = this_name.rfind(this_ext)
    if this_pos > 0:
        this_stub = this_name[:this_pos]

    return this_stub


def check_hit(hit, max_hits):
    # max_hits is a dict with keys 1, 2,3 [not 0]
    for i in range(1, len(max_hits)+1):
        if hit[0] > max_hits[i][0]:
            for j in range(len(max_hits),i,-1):
                max_hits[j] = max_hits[j-1]

            max_hits[i] = hit
            break

    return max_hits


def setup_revcodes():

    rev_codes = {}
    rev_codes["M"] = "K"
    rev_codes["R"] = "Y"
    rev_codes["W"] = "W"
    rev_codes["S"] = "S"
    rev_codes["Y"] = "R"
    rev_codes["K"] = "M"
    rev_codes["V"] = "B"
    rev_codes["H"] = "D"
    rev_codes["D"] = "H"
    rev_codes["B"] = "V"
    rev_codes["N"] = "N"
    rev_codes["-"] = "-"
    rev_codes["A"] = "T"
    rev_codes["C"] = "G"
    rev_codes["G"] = "C"
    rev_codes["T"] = "A"

    return rev_codes


def setup_ambiguities():

    ambis = {}
    ambis["M"] = ["A", "C"]
    ambis["R"] = ["A", "G"]
    ambis["W"] = ["A", "T"]
    ambis["S"] = ["C", "G"]
    ambis["Y"] = ["C", "T"]
    ambis["K"] = ["G", "T"]
    ambis["V"] = ["A", "C", "G"]
    ambis["H"] = ["A", "C", "T"]
    ambis["D"] = ["A", "G", "T"]
    ambis["B"] = ["C", "G", "T"]
    ambis["N"] = ["A", "C", "G", "T"]
    ambis["A"] = ["A"]
    ambis["C"] = ["C"]
    ambis["G"] = ["G"]
    ambis["T"] = ["T"]
    ambis["-"] = ["-"]

    return ambis


def is_ambi(base, ambis):
    base = base.upper()
    isambi = False

    if base not in ["A", "C", "G", "T", "-"]:
        if base in ambis:
            isambi = True
        else:
            print("Error - unrecognised base: " + base)

    return isambi


def checkbase(b1, b2, ambis):
    code1 = ambis[b1]
    code2 = ambis[b2]

    match = False

    for c1 in code1:
        for c2 in code2:
            if c1 == c2:
                match = True
                break

    return match

primer_filename = sys.argv[1]
seq_filename = sys.argv[2]
output_filename = get_filename_stub(seq_filename, ".fa") + "_bind.csv"

print("Primer file = " + primer_filename)
print("Seq file = " + seq_filename)
print("Output file = " + output_filename)

ambis = setup_ambiguities()
primer_count = 0
primers = {}
primer_names = []
primer_seqs = []

for record in SeqIO.parse(sys.argv[1], "fasta"):
    primer_count += 1
    primer_name = record.description
    primer_seq = record.seq.upper()
    primer_revcomp = record.seq.upper().reverse_complement()

    primers[primer_count] = [primer_name, primer_seq, primer_revcomp]

    if primer_name in primer_names:
        print("Warning - primer name is duplicated: " + primer_name)
    if primer_seq in primer_seqs:
        print("Warning - primer seq is duplicated: " + primer_seq)

print("Primer count = " + str(primer_count))

with open(output_filename, "w") as file_output:
    for record in SeqIO.parse(sys.argv[2], "fasta"):
        seq_name = record.description
        seq = record.seq.upper()

        for primer in primers:
            primer_name = primers[primer][0]
            primer_seq = primers[primer][1].upper()
            primer_revcomp = primers[primer][2]

            plen = len(primer_seq)
            slen = len(seq)

            max_hits = {}
            # freq, start_pos, end_pos, F/R, seq match string [dots], primer match string
            max_hits[1] = [0.0,0,0,"?","",""]
            max_hits[2] = [0.0,0,0,"?","",""]
            max_hits[3] = [0.0,0,0,"?","",""]

            window_seq = ""

            # would be good to allow going over end
            # just got to end, create the kmer dynaimcally - as want to skip over gaps - if kmer is plen or reach end then test
            for i in range(0, slen-plen):

                if seq[i] == "-":
                    continue

                window_seq += seq[i]

                if len(window_seq) == plen:
                    start_pos = i+1 - plen + 1
                    end_pos = start_pos + plen -1

                    #this_seq = seq[i:i+plen]
                    this_seq = window_seq
                    match_string = ""
                    primer_string = ""
                    pmatch = 0

                    direction = "F"

                    for j in range(0, plen):
                        if checkbase(this_seq[j], primer_seq[j], ambis):
                            pmatch += 1
                            match_string += "."
                            primer_string += primer_seq[j].lower()
                        else:
                            match_string += this_seq[j]
                            primer_string += primer_seq[j]

                    pfreq = pmatch / plen

                    hit = [pfreq, start_pos, end_pos, direction, match_string, primer_string]

                    if this_seq.upper().count("N") < 3:
                        max_hits = check_hit(hit, max_hits)

                    match_string = ""
                    primer_string = ""
                    pmatch = 0
                    pfreq = 0
                    direction = "R"

                    for j in range(0, plen):
                        if checkbase(this_seq[j], primer_revcomp[j], ambis):
                            pmatch += 1
                            match_string += "."
                            primer_string += primer_seq[j].lower()
                        else:
                            match_string += this_seq[j]
                            primer_string += primer_seq[j]

                    pfreq = pmatch / plen

                    hit = [pfreq, start_pos, end_pos, direction, match_string, primer_string]

                    if this_seq.upper().count("N") < 3:
                        max_hits = check_hit(hit, max_hits)

                    window_seq = window_seq[1:]

            # Now we have a end loop
            # Start of the genome
            # Trime the primer 1 base at a time from the 5' end - see if
            out_str = primer_name + "," + str(plen)

            for hit in max_hits:
                for h in max_hits[hit]:
                    out_str += "," + str(h)

            out_str += "\n"
            file_output.write(out_str)



print("finished goprime_comp.py")
