import argparse

parser = argparse.ArgumentParser(
    prog='convertf_p',
    description='Generate pipes for convertf on the fly',
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

    args, _ = parser.parse_known_args()


def main():
    parse_args()

    print(f'\
        genotypename:    {args.input_set}.{get_geno_file_ext_from_format(args.from_format)}\n\
        snpname:         {args.input_set}.{get_snp_file_ext_from_format(args.from_format)}\n\
        indivname:       {args.input_set}.{get_indiv_file_ext_from_format(args.from_format)}\n\
        outputformat:    {get_convertf_format_from_format(args.to_format)}\n\
        genotypeoutname: {args.output_set}.{get_geno_file_ext_from_format(args.to_format)}\n\
        snpoutname:      {args.output_set}.{get_snp_file_ext_from_format(args.to_format)}\n\
        indivoutname:    {args.output_set}.{get_indiv_file_ext_from_format(args.to_format)}\n\
        familynames:     NO\
    ')


if __name__ == '__main__':
    main()