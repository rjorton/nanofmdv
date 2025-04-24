import sys

line_count = 0

with open(sys.argv[1]) as file_handler:
    for line in file_handler:
        line_count += 1

if line_count % 4 == 0:
    read_count=int(line_count/4)
else:
    read_count = line_count / 4

print(read_count)