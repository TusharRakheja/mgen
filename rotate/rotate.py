import subprocess
import itertools
import argparse
import hashlib
import yaml
import os

from concurrent.futures import ThreadPoolExecutor
from queue import Queue

CORE_RIGHT = 'core_right'
TARGET = 'target'
TARGETS = 'targets'
SOURCES = 'sources'

RESULTS_QUEUE = None

args = None
IND = None
SNP = None
GENO = None
CONFIG_FILE = None
MODELS_POOL = None
PID = None
SUFFIX = None
RESULTS = {}
CACHE = {}

TARGET_LIST = []

parser = argparse.ArgumentParser(
    prog='rotate',
    description='Run rotating qpAdm models, with or without competition',
    epilog='Copyright (c) 2025 Tushar Rakheja (The MIT License)'
)

def is_valid(conf):
    global CORE_RIGHT, TARGET, TARGETS, TARGET_LIST, SOURCES, IND

    if CORE_RIGHT not in conf.keys():
        return False
    if not isinstance(conf[CORE_RIGHT], list):
        return False
    for source in conf[CORE_RIGHT]:
        if not isinstance(source, str):
            return False

    has_target = TARGET in conf
    has_targets = TARGETS in conf

    if not has_target and not has_targets:
        return False

    if has_target and not isinstance(conf[TARGET], str):
        return False

    if has_targets:
        if not isinstance(conf[TARGETS], list):
            return False
        for target in conf[TARGETS]:
            if not isinstance(target, str):
                return False

    if not isinstance(conf[SOURCES], list):
        return False
    for sources in conf[SOURCES]:
        if not isinstance(sources, list):
            return False
        for source in sources:
            if not isinstance(source, str):
                return False

    with open(IND, 'r') as infile:
        l = infile.readlines()

    samples = set()
    for line in l:
        line = line.strip()
        if not line:
            continue

        parts = line.split()
        samples.add(parts[-1])

    if TARGETS in conf.keys():
        TARGET_LIST = conf[TARGETS]
    
    if TARGET in conf.keys():
        TARGET_LIST.append(conf[TARGET])

    for source in conf[CORE_RIGHT]:
        if source not in samples:
            print("Source {} not found in the {} file".format(source, IND))
            return False

    for target in TARGET_LIST:
        if target not in samples:
            print("Target {} not found in the {} file".format(target, IND))
            return False

    for sources in conf[SOURCES]:
        for source in sources:
            if source != "" and source not in samples:
                print("Source {} not found in the {} file".format(source, IND))
                return False

    return True


def run(model_number, target, model, source_sets, core_sources):
    global PID

    try:
        LEFT = './left_{}_{}_{}'.format(target, PID, model_number)
        RIGHT = './right_{}_{}_{}'.format(target, PID, model_number)
        OUTPUT = './output_{}_{}_{}'.format(target, PID, model_number)
        PARQPADM = './parqpadm_{}_{}_{}'.format(target, PID, model_number)

        generate_parqpadm(PARQPADM, LEFT[2:], RIGHT[2:], target)

        with open(LEFT, 'w') as outfile:
            outfile.write("{}\n".format(target))
            for source in model:
                source = source.strip()
                if not source:
                    continue

                outfile.write("{}\n".format(source))

        with open(RIGHT, 'w') as outfile:
            for source in core_sources:
                outfile.write("{}\n".format(source))
            for sources in source_sets:
                for source in sources:
                    source = source.strip()
                    if not source:
                        continue

                    if source not in model:
                        if not args.no_compete:
                            outfile.write("{}\n".format(source))


        with open(OUTPUT, 'w') as outfile:
            print("{} - Starting qpAdm on model {} for {}".format(model_number, model, target))
            subprocess.call([*(['wsl'] if args.turn_on_wsl_for_admix_tools else []), 'qpAdm', '-p', PARQPADM[2:]], stdout=outfile)

        weights, errors, pvalue = weights_errors_pvalue(OUTPUT)
        return [target, model, weights, errors, pvalue]
    
    except:
        return [target, model, None, None, None]

    finally:
        clean_up_model_files(LEFT, RIGHT, OUTPUT, PARQPADM, model)


def write_results(cur_target):
    global RESULTS_QUEUE, RESULTS

    while not RESULTS_QUEUE[cur_target].empty():
        with open(RESULTS[cur_target], 'a') as outfile:
            target, model, weights, errors, pvalue = RESULTS_QUEUE[cur_target].get().result()
            if weights is None:
                print("\tFailed model {}. Will be attempted again if you re-run the script".format(model))
                continue

            outfile.write(result_row(target, model, weights, errors, pvalue))
            add_model_to_cache(model, target)


def generate_parqpadm(parqpadm, left, right, target):
    global IND, SNP, GENO, SUFFIX
    with open(parqpadm, 'w') as outfile:

        if args.fstats:
            fstats_filename = 'fstats_{}_{}'.format(target, SUFFIX) if args.use_fstats_file is None else args.use_fstats_file.split('/')[-1]
            outfile.write('fstatsname:      {}\n'.format(fstats_filename))
        else:
            outfile.write('indivname:       {}\n'.format(IND))
            outfile.write('snpname:         {}\n'.format(SNP))
            outfile.write('genotypename:    {}\n'.format(GENO))

        outfile.write('popleft:         {}\n'.format(left))
        outfile.write('popright:        {}\n'.format(right))
        outfile.write('details:         YES\n')
        outfile.write('allsnps:         YES\n')
        outfile.write('inbreed:         NO\n')


def generate_poplist(conf, target):
    global CORE_RIGHT, SOURCES, SUFFIX
    
    POPLIST = './poplist_{}_{}'.format(target, SUFFIX)

    if os.path.isfile(POPLIST):
        # already exists
        return

    with open(POPLIST, 'w') as outfile:
        for source in conf[CORE_RIGHT]:
            outfile.write("{}\n".format(source))

        for sources in conf[SOURCES]:
            for source in sources:
                outfile.write("{}\n".format(source))
        
        outfile.write(target)


def generate_parqpfstats(target):
    global IND, SNP, GENO, SUFFIX

    PARQPFSTATS = './parqpfstats_{}_{}'.format(target, SUFFIX)

    if os.path.isfile(PARQPFSTATS):
        # already exists
        return

    with open(PARQPFSTATS, 'w') as outfile:
        outfile.write('indivname:       {}\n'.format(IND.split('/')[-1].strip()))
        outfile.write('snpname:         {}\n'.format(SNP.split('/')[-1].strip()))
        outfile.write('genotypename:    {}\n'.format(GENO.split('/')[-1].strip()))
        outfile.write('poplistname:     poplist_{}_{}\n'.format(target, SUFFIX))
        outfile.write('fstatsoutname:   fstats_{}_{}\n'.format(target, SUFFIX))
        outfile.write('allsnps:         YES\n')
        outfile.write('inbreed:         NO\n')
        outfile.write('scale:           NO\n')


def run_qpfstats(target):
    global SUFFIX

    QP_LOG_OUTPUT = './qpfstats_log_{}_{}'.format(target, SUFFIX)

    if os.path.isfile('./fstats_{}_{}'.format(target, SUFFIX)):
        # already exists
        return

    if args.keep_fstats_log:
        with open(QP_LOG_OUTPUT, 'w') as outfile:
            subprocess.call([*(['wsl'] if args.turn_on_wsl_for_admix_tools else []), 'qpfstats', '-p', 'parqpfstats_{}_{}'.format(target, SUFFIX)], stdout=outfile)
    else:
        subprocess.call([*(['wsl'] if args.turn_on_wsl_for_admix_tools else []), 'qpfstats', '-p', 'parqpfstats_{}_{}'.format(target, SUFFIX)])


def clean_up_model_files(left, right, output, parqpadm, model):
    try:
        os.remove(left)
        os.remove(right)

        if not args.keep_output_files:
            os.remove(output)

        os.remove(parqpadm)
    except OSError:
        print("Error removing files for model {}. Runs are unaffected though.".format(model))


def clean_up_fstats(target):
    global SUFFIX
    try:
        if args.use_fstats_file is None:
            os.remove('./parqpfstats_{}_{}'.format(target, SUFFIX))

        if (args.use_fstats_file is None) and (not args.keep_fstats_file):
            os.remove('./fstats_{}_{}'.format(target, SUFFIX))

        if args.use_fstats_file is None:
            os.remove('./poplist_{}_{}'.format(target, SUFFIX))
    except:
        print("Error removing files for fstats. Runs are unaffected though.")
        pass


def weights_errors_pvalue(output):
    with open(output, 'r') as infile:
        l = infile.readlines()
    
    weights = None
    errors = None
    pvalue = None

    read_pvalue = False

    for line in l:
        line = line.strip()
        if not line:
            continue
        
        if line.startswith('best coefficients:'):
            weights = []
            for weight in line.split()[2:]:
                weights.append(weight)
            continue

        if line.startswith('std. errors:'):
            errors = []
            for error in line.split()[2:]:
                errors.append(error)
            continue

        if 'tail prob' in line:
            read_pvalue = True
            continue

        if read_pvalue:
            pvalue = line.split()[4]
            break

    assert weights is not None
    assert errors is not None
    assert pvalue is not None
    assert read_pvalue is True

    return [weights, errors, pvalue]


def result_row(target, model, weights, errors, pvalue):
    res = ""
    res += "{},".format(target)
    for source in model:
        res += "{},".format(source if source.strip() != "" else "-")

    j = 0
    for i in range(len(model)):
        if model[i].strip() == "":
            res += "0%,"
        else:
            res += "{},".format(str("%.1f" % (float(weights[j])*100)) + '%')
            j += 1

    j = 0
    for i in range(len(model)):
        if model[i].strip() == "":
            res += "0%,"
        else:
            res += "{},".format(str("%.1f" % (float(errors[j])*100)) + '%')
            j += 1

    res += "{},{},".format(pvalue, len(list(filter(lambda source: source.strip() != "", model))))

    ispass = float(pvalue) >= args.pmin

    for i in range(len(weights)):
        ispass = ispass and (float(weights[i]) >= 0) and (float(errors[i]) >= 0)

    res += "{}\n".format(1 if ispass else 0)

    return res


def write_headers(num_sources, target):
    global RESULTS

    if os.path.isfile(RESULTS[target]):
        return

    outfile = open(RESULTS[target], 'w')

    outfile.write('Target,')
    for i in range(num_sources if args.rank is None else args.rank):
        outfile.write('Source {},'.format(i + 1))
    for i in range(num_sources if args.rank is None else args.rank):
        outfile.write('Weight {},'.format(i + 1))
    for i in range(num_sources if args.rank is None else args.rank):
        outfile.write('Error {},'.format(i + 1))
    outfile.write('p-value,Complexity,Pass?\n')

    outfile.close()


def add_model_to_cache(model, target):
    global CACHE

    with open(CACHE[target], 'a') as outfile:
        outfile.write(','.join(model) + '\n')


def is_model_in_cache(model, target):
    global CACHE

    key = ','.join(model)

    try:
        infile = open(CACHE[target], 'r')
    except OSError:
        return False

    l = infile.readlines()

    infile.close()
    
    for line in l:
        if line.strip() == key:
            return True
    
    return False


def parse_args():
    global parser, RESULTS_QUEUE, args, IND, SNP, GENO, CONFIG_FILE, MODELS_POOL, PID, SUFFIX, RESULTS, CACHE

    RESULTS_QUEUE = {}
    parser.add_argument('-s', '--set', dest='set', type=str, help='Path to and prefix of your set (e.g. "/path/to/set" if the files are "/path/to/set.ind", "path/to/set.geno" and "path/to/set.snp")', required=True)
    parser.add_argument(
        '-c', '--config', dest='config', type=str, default='./config.yml', help='Path to the YAML config file (default: "./config.yml")'
    )
    parser.add_argument(
        '-n', '--nthreads', dest='nthreads', type=int, default=1, help='The number of models to run in parallel (default: 1)'
    )
    parser.add_argument('-p', '--pmin', dest='pmin', type=float, default=0.05, help='Minimum viable p-value')
    parser.add_argument('--fstats', default=False, action='store_true', help='Pre-compute fstats for the whole config using qpfstats')
    parser.add_argument('--no-compete', default=False, action='store_true', help='Run the models in the config file without model competition')
    parser.add_argument('--keep-fstats-file', default=False, action='store_true', help='Do not delete the pre-computed fstats file')
    parser.add_argument('--keep-fstats-log', default=False, action='store_true', help="Do not delete the qpfstats log file")
    parser.add_argument('--use-fstats-file', dest='use_fstats_file', type=str, default=None, help='Use pre-computed fstats file')
    parser.add_argument('--turn-on-wsl-for-admix-tools', default=False, action='store_true', help="(If running on Windows) Use AdmixTools like qpAdm and qpfstats via WSL")
    parser.add_argument('-r', '--rank', dest='rank', type=int, default=None, help='The rank of the models you want to run (i.e. how many sources per model)?')
    parser.add_argument('-mr', '--min-rank', dest='min_rank', type=int, default=None, help='The min rank of the models you want to run (i.e. how many sources per model)?')
    parser.add_argument('--keep-output-files', default=False, action='store_true', help="Do not delete the qpAdm output for each model")
    parser.add_argument('--dry-run', default=False, action='store_true', help="Dry run (creates an empty results file but a full cache file)")
    parser.add_argument('--max-rank', default=False, action='store_true', help="The --rank argument should be treated as a max rank")

    args, _ = parser.parse_known_args()

    IND = "{}.ind".format(args.set)
    SNP = "{}.snp".format(args.set)
    GENO = "{}.geno".format(args.set)

    CONFIG_FILE = args.config
    MODELS_POOL = ThreadPoolExecutor(args.nthreads)

    PID = os.getpid()

    with open(CONFIG_FILE, 'rb') as tempconfig:
        SUFFIX = hashlib.md5(tempconfig.read()).hexdigest()


def get_model_list(source_sets):
    global args

    if args.rank is None:
        return list(itertools.product(*source_sets))
    else:
        all_sources = list(set(itertools.chain.from_iterable(source_sets)))
        all_sources.remove('')

        if not args.max_rank:
            return list(itertools.combinations(all_sources, args.rank))
        else:
            all_models = []
            for rank in range(1 if args.min_rank is None else args.min_rank, args.rank + 1):
                models_ = list(itertools.combinations(all_sources, rank))
                models = []
                for model in models_:
                    models.append(list(model))

                for model in models:
                    while len(model) < args.rank:
                        model.append("")

                all_models += models

            return all_models


def main():
    global CORE_RIGHT, TARGET_LIST, SOURCES, CONFIG_FILE, MODELS_POOL, RESULTS, CACHE, SUFFIX

    parse_args()

    with open(CONFIG_FILE, 'r') as infile:
        config = None
        try:
            config = yaml.safe_load(infile)
        except yaml.YAMLError as exc:
            print(exc)

    assert config is not None
    assert is_valid(config)

    source_sets = config[SOURCES]

    models = get_model_list(source_sets)

    for target in TARGET_LIST:
        RESULTS[target] = './results_{}_{}.csv'.format(target, SUFFIX)
        CACHE[target] = './cache_{}_{}.txt'.format(target, SUFFIX)

        RESULTS_QUEUE[target] = Queue()

        write_headers(len(source_sets), target)
        
        print("Will try {} models for {}".format(len(models, target)))

        if not args.dry_run and args.fstats and args.use_fstats_file is None:
            print("Running qpfstats (can take a while) ...")
            generate_poplist(config, target)
            generate_parqpfstats(target)
            run_qpfstats(target)

        model_number = 1
        for model in models:
            if is_model_in_cache(model, target):
                print("{} - Skipping model {} for {} because already in cache.".format(model_number, model, target))
                model_number += 1
                continue

            if len(list(filter(lambda source: source.strip() != "", model))) < 2:
                print("{} - Skipping model {} for {} because we need at least two sources.".format(model_number, model, target))
                model_number += 1
                continue

            if args.dry_run:
                print("{} - Running model {} for {}".format(model_number, model, target))
                add_model_to_cache(model, target)
            else:
                task = MODELS_POOL.submit(run, model_number, target, model, source_sets, config[CORE_RIGHT])
                RESULTS_QUEUE[target].put(task)
            model_number += 1


    for target in TARGET_LIST:
        write_results(target)

        if args.fstats:
            clean_up_fstats(target)


if __name__ == '__main__':
    main()