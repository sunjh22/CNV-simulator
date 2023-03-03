# CNV-simulator
A Python-based whole-genome haplotype-resolution copy number variation simulator.

Download the CNV-simulator

    git clone https://github.com/sunjh22/CNV-simulator.git

Then go to the `CNV-simulator` directory.

You are actually good to run CNV-simulator as follow:

    ./cnv_simulator.py -o /path/to/store/simulated/data -a prefix_of_sample -c desired_coverage reference_genome resource/access-excludes.hg38.analysisSet.bed

You can also specify the number of CNVs to be simulated, minimum, maximum and average CNV length.

The detailed parameters of CNV-simulator

    Usage: cnv_simulator.py [-h] [-v] [-o OUTPUT_DIR] [-a PREFIX] [-l READ_LENGTH]
                            [-i CNV_LIST] [-c COVERAGE] [-n CNV_NUMBER]
                            [-b CNV_MINIMUM_LENGTH] [-B CNV_MAXIMUM_LENGTH]
                            [-e CNV_EXPO_SCALE] [-p AMPLIFY_PROP]
                            genome access

    Simulate whole-genome level haplotype-resolved CNVs and generate NGS reads

    Positional arguments:
    genome                reference genome file in fasta format
    access                genomic accessible region file

    Optional arguments:
    -h, --help            show this help message and exit
    -v, --version         Show program's version number and exit.
    -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                            a name to be used to create the output directory
                            (overrides existing directory with the same name).
                            (default: tmp_simulaiton)
    -a PREFIX, --prefix PREFIX
                            prefix of simulated genome and fastq files (default:
                            sample)
    -l READ_LENGTH, --read_length READ_LENGTH
                            read length (bp) (default: 150)
    -i CNV_LIST, --cnv_list CNV_LIST
                            path to a CNV list file in BED format chr | start |
                            end | cn. If not passed, it is randomly generated
                            using CNV list parameters below (default: None)
    -c COVERAGE, --coverage COVERAGE
                            coverage for reads simulation (default: 1)

    CNV list parameters:
    parameters to be used to generate CNV list if CNV list is not passed

    -n CNV_NUMBER, --cnv_number CNV_NUMBER
                            the number of CNVs to be simulated (default: 200)
    -b CNV_MINIMUM_LENGTH, --cnv_minimum_length CNV_MINIMUM_LENGTH
                            minimum length of each CNV region (default: 5000)
    -B CNV_MAXIMUM_LENGTH, --cnv_maximum_length CNV_MAXIMUM_LENGTH
                            maximum length of each CNV region (default: 2000000)
    -e CNV_EXPO_SCALE, --cnv_expo_scale CNV_EXPO_SCALE
                            scale of a exponential distribution for simulating CNV
                            length (default: 200000)
    -p AMPLIFY_PROP, --amplify_prop AMPLIFY_PROP
                            percentage of amplifications in range [0.0: 1.0].
                            (default: 0.5)

