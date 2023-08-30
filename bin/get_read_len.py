#! /usr/bin/env python
import argparse
import glob
import math

def argParser():
    parser = argparse.ArgumentParser(description="Get read length from fastqc data text file.")
    parser.add_argument('--prefix', help="Prefix (or sample id) for fastq files.")
    parser.add_argument('--workdir', help="Path to fastqc results.")
    return(parser.parse_args())


def main():
    args = argParser()
    files = glob.glob(f"{args.workdir}/{args.prefix}*/fastqc_data.txt")
    # print(f"{args.workdir}/{args.prefix}*/fastqc_data.txt")
    read_length_list = []
    for f in files:
        with open(f) as inF:
            for i, line in enumerate(inF):
                if i == 8 and line.startswith("Sequence length"):
                    read_length_list.append([int(_) for _ in line.strip()[16:].split("-")])
    
    read_length_list2 = [math.ceil(sum(i)/len(i)) for i in read_length_list]
    avg_read_length = max(read_length_list2)
    print(f"{avg_read_length}", end='')

if __name__ == "__main__":
    main()