import sys


def get_filename_stub(this_name, this_ext):
    this_stub = this_name

    this_pos = this_name.rfind(this_ext)
    if this_pos > 0:
        this_stub = this_name[:this_pos]

    return this_stub


goprime_filename = sys.argv[1]
output_filename = get_filename_stub(goprime_filename, ".csv") + ".bed"

print("goprime file = " + goprime_filename)
print("Output file = " + output_filename)

line_count = 0
primers = {}

with open(goprime_filename) as file_handler:
    for line in file_handler:
        line = line.strip()
        line_count += 1
        splits = line.split(",")
        primer_name = splits[0]


        primer_len = int(splits[1])
        primer_freq = float(splits[2])
        primer_start = int(splits[3])
        primer_end = int(splits[4])

        if primer_name not in primers:
            primers[primer_name] = [primer_freq, primer_start-1, primer_end, primer_len]

        else:
            if primer_freq > primers[primer_name][0]:
                primers[primer_name] = [primer_freq, primer_start - 1, primer_end, primer_len]


with open(output_filename, "w") as file_output:
    for p in primers:
        file_output.write("fmdvref\t" + str(primers[p][1]) + "\t" + str(primers[p][2]) + "\t" + p + "\n")




print("finished goprime_bed.py")
