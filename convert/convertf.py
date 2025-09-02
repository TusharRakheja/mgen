import argparse
import json
import subprocess

parser = argparse.ArgumentParser(
    prog='convertf',
    description='Convert sets from one format to another',
    epilog='Copyright (c) 2025 Tushar Rakheja (The MIT License)'
)

args = None


def get_geno_file_ext_from_format(format):
    if format == 'eigenstrat':
        return 'geno'
    if format == 'plink':
        return 'bed'
    if format == 'tgeno':
        return 'tgeno'


def get_snp_file_ext_from_format(format):
    if format == 'eigenstrat':
        return 'snp'
    if format == 'plink':
        return 'bim'
    if format == 'tgeno':
        return 'snp'


def get_indiv_file_ext_from_format(format):
    if format == 'eigenstrat':
        return 'ind'
    if format == 'plink':
        return 'fam'
    if format == 'tgeno':
        return 'ind'


def get_convertf_format_from_format(format):
    if format == 'eigenstrat':
        return 'PACKEDANCESTRYMAP'
    if format == 'plink':
        return 'PACKEDPED'
    if format == 'tgeno':
        return 'TGENO'


def parse_args():
    global args
    parser.add_argument('-i', '--input-set', dest='input_set', type=str, help='Prefix of your input set (e.g. just "prefix" if prefix.geno, prefix.ind, prefix.snp)', required=True)
    parser.add_argument('-f', '--from', dest='from_format', type=str, help='Current format of set', choices=["eigenstrat", "plink", "tgeno"], required=True)
    parser.add_argument('-t', '--to', dest='to_format', type=str, help='The format to convert the set to', choices=["eigenstrat", "plink", "tgeno"], required=True)
    parser.add_argument('-o', '--output-set', dest='output_set', type=str, help='Prefix of your output set (e.g. just "prefix" if prefix.geno, prefix.ind, prefix.snp)', required=True)
    parser.add_argument('-c', '--convertf', dest='convertf', default="convertf", type=str, help='Command used to invoke convertf on your machine (e.g. "convertf")')
    parser.add_argument('--turn-on-wsl-for-admix-tools', default=False, action='store_true', help="(If running on Windows) Use AdmixTools like convertf via WSL")

    args, _ = parser.parse_known_args()


def fix_family_names_ind_to_ind(input_set, output_set):
    indiv_to_family_name = {}

    with open("{}.ind".format(input_set), 'r') as infile:
        l = infile.readlines()
    
    for line in l:
        line = line.strip()
        if not line:
            continue

        parts = line.split()
        indiv_to_family_name[parts[0]] = parts[2]
    
    with open("{}.ind".format(output_set), 'r') as infile:
        l = infile.readlines()

    with open("{}.ind".format(output_set), 'w') as outfile:
        for line in l:
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            if parts[0] in indiv_to_family_name.keys():
                outfile.write("\t{}\t{}\t{}\n".format(parts[0], parts[1], indiv_to_family_name[parts[0]]))
            else:
                outfile.write("{}\n".format(line))


def fix_family_names_ind_to_fam(input_set, output_set):
    indiv_to_family_name = {}

    with open("{}.ind".format(input_set), 'r') as infile:
        l = infile.readlines()
    
    for line in l:
        line = line.strip()
        if not line:
            continue

        parts = line.split()
        indiv_to_family_name[parts[0]] = parts[2]
    
    with open("{}.fam".format(output_set), 'r') as infile:
        l = infile.readlines()

    with open("{}.fam".format(output_set), 'w') as outfile:
        for line in l:
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            if parts[1] in indiv_to_family_name.keys():
                outfile.write("\t{}\t{}\t{}\t{}\t{}\t{}\n".format(indiv_to_family_name[parts[1]], parts[1], parts[2], parts[3], parts[4], parts[5]))
            else:
                outfile.write("{}\n".format(line))


def fix_family_names_fam_to_ind(input_set, output_set):
    indiv_to_family_name = {}

    print("Called")

    print("Opening file: {}.fam".format(input_set))

    with open("{}.fam".format(input_set), 'r') as infile:
        l = infile.readlines()
    
    for line in l:
        line = line.strip()
        if not line:
            continue

        parts = line.split()
        indiv_to_family_name[parts[1]] = parts[0]

    print(json.dumps(indiv_to_family_name))

    with open("{}.ind".format(output_set), 'r') as infile:
        l = infile.readlines()

    with open("{}.ind".format(output_set), 'w') as outfile:
        for line in l:
            line = line.strip()
            if not line:
                continue

            print("Read line: {}".format(line))
            parts = line.split()
            print("Parts: {}".format(json.dumps(parts)))
            if parts[0] in indiv_to_family_name.keys():
                print("parts[0] i.e. [{}] in indiv_to_family_name.keys()? {}".format(parts[0], parts[0] in indiv_to_family_name.keys()))
                print("Writing: [\t{}\t{}\t{}]".format(parts[0], parts[1], indiv_to_family_name[parts[0]]))
                outfile.write("\t{}\t{}\t{}\n".format(parts[0], parts[1], indiv_to_family_name[parts[0]]))
            else:
                print("parts[0] i.e. [{}] in indiv_to_family_name.keys()? {}".format(parts[0], parts[0] in indiv_to_family_name.keys()))
                print("Writing: [{}]".format(line))
                outfile.write("{}\n".format(line))


def fix_family_names(from_format, to_format, input_set, output_set):
    if from_format == 'eigenstrat' or from_format == 'tgeno':
        if to_format == 'eigenstrat' or to_format == 'tgeno':
            fix_family_names_ind_to_ind(input_set, output_set)
        if to_format == 'plink':
            fix_family_names_ind_to_fam(input_set, output_set)
    if from_format == 'plink':
        if to_format == 'eigenstrat' or to_format == 'tgeno':
            fix_family_names_fam_to_ind(input_set, output_set)
        

def main():
    parse_args()

    # convertf -p <(mgen convertf_p -i <input_set> -o <output_set> -f <from_format> -t <to_format>)

    convertf_p_cmd = [
        "mgen", "convertf_p",
        "-i", args.input_set,
        "-o", args.output_set,
        "-f", args.from_format,
        "-t", args.to_format
    ]

    p1 = subprocess.Popen(convertf_p_cmd, stdout=subprocess.PIPE)

    convertf_cmd = [*(['wsl'] if args.turn_on_wsl_for_admix_tools else []), args.convertf, "-p", "/dev/fd/0"]  # Use /dev/fd/0 to read from stdin

    p2 = subprocess.Popen(convertf_cmd, stdin=p1.stdout)

    p1.stdout.close()

    p2.wait()

    fix_family_names(args.from_format, args.to_format, args.input_set, args.output_set)


if __name__ == '__main__':
    main()