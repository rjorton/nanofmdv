import sys
import os
import subprocess


def print_help():

    print("--------------------")
    print("depth_charge.py")
    print("--------------------")
    print("A simple script to give primary and secondary/supplementary read stats from a BAM file")
    print("And depth/coverage summaries in full mode")
    print("--------------------")
    print("Overall depth_charge.py can be run in multiple ways")
    print("--------------------")
    print("python depth_charge.py my.bam")
    print("This will run quick stats (expanded idxstats with primary and secondary mapped reads)")
    print("Output is to the screen - no progress messages outputted to keep clean")
    print("--------------------")
    print("python depth_charge.py my.bam trans")
    print("Same as above but the output is transposed")
    print("--------------------")
    print("python depth_charge.py my.bam full")
    print("This will run full stats = quick stats plus coverage data e.g. average depth, number of sites > 10 etc")
    print("Output is to my_depth_charge.txt file (samtools depth file my_depth.txt also created)")
    print("--------------------")
    print("python depth_charge.py my.bam full-trans")
    print("Same full stats as above but output is transposed")
    print("--------------------")
    print("python depth_charge.py my.txt full")
    print("If not a bam file, assumes a samtools depth file, will run coverage/depth/breadth stats but not mapping stats")
    print("Outputted is to my_depth_charge.txt file, can also be run with full-trans to trasnspose the output")
    print("--------------------")
    print("REQUIRES samtools in your $PATH: http://www.htslib.org")


def get_filename_stub(this_name, this_ext):
    this_stub = this_name

    this_pos = this_name.rfind(this_ext)
    if this_pos > 0:
        this_stub = this_name[:this_pos]

    return this_stub


def get_mapping_stats(filename, seqs):
    res = {}

    for seq in seqs:
        # store as strings as grabbing the output from command lime
        res[seq] = ["", "", "", ""]

        # samtools view is being used 4 times
        # faster (x5?) to open map and loop through alignments keeping count as go
        run_str = "samtools view -c -F4 " + filename + " " + seq
        mapped = subprocess.run(run_str, shell=True, capture_output=True).stdout.decode("utf-8").strip()

        run_str = "samtools view -c -F4 -F256 -F2048 " + filename + " " + seq
        mapped_true = subprocess.run(run_str, shell=True, capture_output=True).stdout.decode("utf-8").strip()

        run_str = "samtools view -c -F4 -f256 " + filename + " " + seq
        secondary = subprocess.run(run_str, shell=True, capture_output=True).stdout.decode("utf-8").strip()

        run_str = "samtools view -c -F4 -f2048 " + filename + " " + seq
        supplementary = subprocess.run(run_str, shell=True, capture_output=True).stdout.decode("utf-8").strip()

        res[seq][0] = mapped
        res[seq][1] = mapped_true
        res[seq][2] = secondary
        res[seq][3] = supplementary

        print(res[seq])
    return res


def get_unmapped(filename):
    run_str = "samtools view -c -f4 " + filename
    res = subprocess.run(run_str, shell=True, capture_output=True).stdout.decode("utf-8").strip()

    return res


def quick_stats(input_filename, transpose):
    if input_filename[-4:].lower() != ".bam":
        print("Error - was expecting a .bam file for quick stats - exiting")
        sys.exit(1)

    # just to be consistent - a dictionary is used for seqs in full_stats as seqs store the cov data - could be a list a here
    seqs = {}
    ref_lens = {}
    idx_tot = {}

    total_unmapped = ""
    # grab from idxstats
    # total_unmapped = get_unmapped(input_filename)

    run_str = "samtools idxstats " + input_filename
    run_out = subprocess.run(run_str, shell=True, capture_output=True).stdout.decode("utf-8").strip()
    lines = run_out.split("\n")

    for line in lines:
        splits = line.split("\t")

        seq = splits[0]

        if seq == "*":
            total_unmapped = splits[3]
        else:
            if seq in seqs:
                print("Error - seq already in seqs idxstats: " + seq)
            else:
                seqs[seq] = {}
                ref_lens[seq] = splits[1]
                idx_tot[seq] = splits[2]

    mapstats = get_mapping_stats(input_filename, seqs)

    # Columns are the seqs/segs - rows are the read numbers
    if transpose:
        out_str = "RefName"
        for seq in seqs:
            out_str += "\t" + seq
        print(out_str)

        out_str = "RefLength"
        for seq in ref_lens:
            out_str += "\t" + ref_lens[seq]
        print(out_str)

        out_str = "MappedReads"
        for m in mapstats:
            out_str += "\t" + mapstats[m][1]
        print(out_str)

        out_str = "Alignments"
        for m in mapstats:
            out_str += "\t" + mapstats[m][0]
        print(out_str)

        out_str = "Secondary"
        for m in mapstats:
            out_str += "\t" + mapstats[m][2]
        print(out_str)

        out_str = "Supplementary"
        for m in mapstats:
            out_str += "\t" + mapstats[m][3]
        print(out_str)

        out_str = "Unmapped[BAM-Total]"
        for i in range(0, len(mapstats)):
            out_str += "\t" + total_unmapped
        print(out_str)

    else:
        print("RefName\tRefLength\tMappedReads\tAlignments\tSecondary\tSupplementary\tUnmapped[BAM-Total]")

        for m in mapstats:
            print(m + "\t" + ref_lens[m] + "\t" + mapstats[m][1] + "\t" + mapstats[m][0] + "\t" + mapstats[m][2] + "\t" + mapstats[m][3] + "\t" + total_unmapped)


def full_stats(input_filename, transpose):
    bam_file = False
    total_unmapped = 0

    if input_filename[-4:] == ".bam":
        bam_file = True

        print("Input BAM file = " + input_filename)
        output_stub = get_filename_stub(input_filename, ".bam")
        depth_filename = output_stub + "_depth.txt"
        output_filename = output_stub + "_depth_charge.txt"
        print("Output depth file = " + depth_filename)

        os.system("samtools depth -aa -d 0 " + input_filename + " > " + depth_filename)

        total_unmapped = get_unmapped(input_filename)
    else:
        print("Input depth file = " + input_filename)
        depth_filename = input_filename
        output_stub = get_filename_stub(input_filename, ".txt")
        output_filename = output_stub + "_depth_charge.txt"

    print("Output filename = " + output_filename)

    seqs = {}
    line_count = 0

    with open(depth_filename) as file_handler:
        for line in file_handler:

            line_count += 1
            line = line.strip()

            splits = line.split("\t")
            seq = splits[0]
            pos = int(splits[1])
            cov = int(splits[2])

            if seq not in seqs:
                seqs[seq] = {}

            if pos not in seqs[seq]:
                seqs[seq][pos] = cov
            else:
                print("Error - position " + str(pos) + " already has data for genome " + seq)

    lengths = {}
    depths = {}
    breadths = {}
    zeroes = {}
    mapstats = {}
    thresholds = [1, 3, 5, 10, 20, 40, 50, 100, 500, 1000, 5000, 10000, 50000, 100000]

    blank = {}
    for thr in thresholds:
        blank[thr] = float(0)

    print("Number of positions (i.e. length) per chr/segment/genome:")
    for seq in seqs:
        print(seq + " = " + str(len(seqs[seq])))
        lengths[seq] = len(seqs[seq])
        depths[seq] = float(0)
        breadths[seq] = blank.copy()
        zeroes[seq] = float(0)

    for seq in seqs:
        for pos in seqs[seq]:
            depths[seq] += seqs[seq][pos]

            if seqs[seq][pos] == 0:
                zeroes[seq] += 1
            else:
                for thr in thresholds:
                    if seqs[seq][pos] >= thr:
                        breadths[seq][thr] += 1

        depths[seq] = depths[seq] / lengths[seq]
        zeroes[seq] = zeroes[seq] / lengths[seq] * 100

        for thr in breadths[seq]:
            breadths[seq][thr] = breadths[seq][thr] / lengths[seq] * 100

    if bam_file:
        mapstats = get_mapping_stats(input_filename, seqs)

    with open(output_filename, "w") as file_output:
        # Seqs/Segs are rows, stats are columns
        if transpose:
            head_str = "RefName"
            for seq in seqs:
                head_str += "\t" + seq
            file_output.write(head_str + "\n")

            out_str = "RefLength"
            for seq in lengths:
                out_str += "\t" + str(lengths[seq])
            file_output.write(out_str + "\n")

            out_str = "AvDepth"
            for seq in depths:
                out_str += "\t" + str(depths[seq])
            file_output.write(out_str + "\n")

            out_str = "Zero"
            for seq in zeroes:
                out_str += "\t" + str(zeroes[seq])
            file_output.write(out_str + "\n")

            for thr in thresholds:
                out_str = "Cov>=" + str(thr)
                for seq in breadths:
                    out_str += "\t" + str(breadths[seq][thr])
                file_output.write(out_str + "\n")

            if bam_file:
                out_str = "MappedReads"
                for seq in seqs:
                    out_str += "\t" + mapstats[seq][1]
                file_output.write(out_str + "\n")

                out_str = "Alignments"
                for seq in seqs:
                    out_str += "\t" + mapstats[seq][0]
                file_output.write(out_str + "\n")

                out_str = "Secondary"
                for seq in seqs:
                    out_str += "\t" + mapstats[seq][2]
                file_output.write(out_str + "\n")

                out_str = "Supplementary"
                for seq in seqs:
                    out_str += "\t" + mapstats[seq][3]
                file_output.write(out_str + "\n")

                out_str = "Unmapped"
                for i in range(0, len(seqs)):
                    out_str += "\t" + total_unmapped
                file_output.write(out_str + "\n")

        # else not transposing - stats are columns, seqs/segs are rows
        else:
            head_str = "RefName\tRefLength\tAvDepth\tZero"
            for thr in thresholds:
                head_str += "\tCov>=" + str(thr)

            if bam_file:
                head_str += "\tMappedReads\tAlignments\tSecondary\tSupplementary\tUnmapped[BAM-Total]"

            file_output.write(head_str + "\n")

            for seq in seqs:
                out_str = seq + "\t" + str(lengths[seq]) + "\t" + str(depths[seq]) + "\t" + str(zeroes[seq])

                for thr in breadths[seq]:
                    out_str += "\t" + str(breadths[seq][thr])

                if bam_file:
                    out_str += "\t" + mapstats[seq][1] + "\t" + mapstats[seq][0] + "\t" + mapstats[seq][2] + "\t" + mapstats[seq][3] + "\t" + total_unmapped

                file_output.write(out_str + "\n")


# WARNING if BAM index does not exist then will get blanks!
# check it exists - chec
arguments = len(sys.argv)
transpose_arg = False

if arguments == 2:
    quick_stats(sys.argv[1], transpose_arg)
elif arguments == 3:
    if sys.argv[2] == "full" or sys.argv[2] == "full-trans":
        if sys.argv[2] == "full-trans":
            transpose_arg = True

        full_stats(sys.argv[1], transpose_arg)
    elif sys.argv[2] == "trans":
        transpose_arg = True
        quick_stats(sys.argv[1], transpose_arg)
    elif sys.argv[2] == "help":
        print_help()
        sys.exit(1)
    else:
        print("Unrecognised 2nd argument [try -> full OR full-trans OR help]: " + sys.argv[2])
        sys.exit(1)
else:
    print("Unrecognised number of arguments: " + str(arguments))
    sys.exit(1)
