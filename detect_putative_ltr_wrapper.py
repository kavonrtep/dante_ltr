#!/usr/bin/env python
"""This wrapper is intended to be used on large genomes and large DANTE input to
minimize memory usage, It splits input files to pieces and analyze it on by one by
detect_putative_ltr.R
If input does not exceed memory limit, it will run detect_putative_ltr.R directly
"""

import argparse
import os
import sys
import tempfile
from itertools import cycle
import subprocess


def get_arguments():
    """
    Get arguments from command line
    :return:
    args
    """
    parser = argparse.ArgumentParser(
        description="""detect_putative_ltr_wrapper.py is a wrapper for 
    detect_putative_ltr.R""", formatter_class=argparse.RawTextHelpFormatter
        )
    parser.add_argument(
        '-g', '--gff3', default=None, required=True, help="gff3 file", type=str,
        action='store'
        )
    parser.add_argument(
        '-s', '--reference_sequence', default=None, required=True,
        help="reference sequence as fasta file", type=str, action='store'
        )
    parser.add_argument(
        '-o', '--output', default=None, required=True, help="output file path and prefix",
        type=str, action='store'
        )
    parser.add_argument(
        '-c', '--cpu', default=1, required=False, help="number of CPUs", type=int,
        action='store'
        )
    parser.add_argument(
        '-M', '--max_missing_domains', default=0, required=False, type=int
        )
    parser.add_argument(
        '-L', '--min_relative_length', default=0.6, required=False, type=float,
        help="Minimum relative length of protein domain to be "
             "considered for retrostransposon detection"
        )
    parser.add_argument(
        '-S', '--max_chunk_size', default=100000000, required=False, type=int,
        help='If size of reference sequence is greater than this value, reference is '
             'analyzed in chunks of this size. This is just approximate value - '
             'sequences '
             'which are longer are are not split, default is %(default)s'
        )
    args = parser.parse_args()
    return args


def read_fasta_sequence_size(fasta_file):
    """Read size of sequence into dictionary"""
    fasta_dict = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line[0] == '>':
                header = line.strip().split(' ')[0][1:]  # remove part of name after space
                fasta_dict[header] = 0
            else:
                fasta_dict[header] += len(line.strip())
    return fasta_dict


def make_temp_files(number_of_files):
    """
    Make named temporary files, file will not be deleted upon exit!
    :param number_of_files:
    :return:
    filepaths
    """
    temp_files = []
    for i in range(number_of_files):
        temp_files.append(tempfile.NamedTemporaryFile(delete=False).name)
    return temp_files


def sum_up_stats_files(files):
    """
    Sum up statistics files
    :return:
    """
    new_statistics = {}
    for file in files:
        with open(file, 'r') as fh:
            for line in fh:
                items = line.strip().split('\t')
                if items[0] == 'Classification':
                    header = items
                    continue
                else:
                    counts = [int(item) for item in items[1:]]
                    if items[0] in new_statistics:
                        new_statistics[items[0]] = [sum(x) for x in
                                                    zip(new_statistics[items[0]], counts)]
                    else:
                        new_statistics[items[0]] = counts
    # convert to string, first line is header
    statistics_str = []
    for classification, counts in new_statistics.items():
        statistics_str.append(classification + '\t' + '\t'.join([str(x) for x in counts]))
    sorted_stat_with_header = ['\t'.join(header)] + sorted(statistics_str)
    return sorted_stat_with_header


def main():
    """
    Main function
    """
    args = get_arguments()
    # locate directory of current script
    tool_path = os.path.dirname(os.path.realpath(__file__))
    fasta_seq_size = read_fasta_sequence_size(args.reference_sequence)
    total_size = sum(fasta_seq_size.values())
    number_of_sequences = len(fasta_seq_size)
    if total_size > args.max_chunk_size and number_of_sequences > 1:
        # sort dictionary by values
        seq_id_size_sorted = [i[0] for i in sorted(
            fasta_seq_size.items(), key=lambda x: int(x[1]), reverse=True
            )]
        number_of_temp_files = int(total_size / args.max_chunk_size) + 1
        if number_of_temp_files > number_of_sequences:
            number_of_temp_files = number_of_sequences

        temp_files_fasta = make_temp_files(number_of_temp_files)
        file_handles = [open(temp_file, 'w') for temp_file in temp_files_fasta]
        # make dictionary seq_id_sorted as keys and values as file handles
        seq_id_file_handle_dict = dict(zip(seq_id_size_sorted, cycle(file_handles)))

        # write sequences to temporary files
        with open(args.reference_sequence, 'r') as f:
            for line in f:
                if line[0] == '>':
                    header = line.strip().split(' ')[0][1:]
                    print(header)
                    seq_id_file_handle_dict[header].write(line)
                else:
                    seq_id_file_handle_dict[header].write(line)
        # close file handles
        for file_handle in file_handles:
            file_handle.close()

        # split gff3 file to temporary files -
        # each temporary file will contain gff lines matching fasta
        temp_files_gff = make_temp_files(number_of_temp_files)
        file_handles = [open(temp_file, 'w') for temp_file in temp_files_gff]
        # make dictionary seq_id_sorted as keys and values as file handles
        seq_id_file_handle_dict = dict(zip(seq_id_size_sorted, cycle(file_handles)))
        # write gff lines to chunks
        with open(args.gff3, 'r') as f:
            for line in f:
                if line[0] == '#':
                    continue
                else:
                    header = line.strip().split('\t')[0]
                    seq_id_file_handle_dict[header].write(line)
        # close file handles
        for file_handle in file_handles:
            file_handle.close()

        # run retrotransposon detection on each temporary file
        output_files = make_temp_files(number_of_temp_files)
        for i in range(number_of_temp_files):
            print('Running retrotransposon detection on file ' + str(i))
            subprocess.check_call(
                [f'{tool_path}/detect_putative_ltr.R', '-s', temp_files_fasta[i], '-g',
                 temp_files_gff[i], '-o', output_files[i], '-c', str(args.cpu), '-M',
                 str(args.max_missing_domains), '-L', str(args.min_relative_length)]
                )

        # remove all temporary input files
        for temp_file in temp_files_fasta + temp_files_gff:
            os.remove(temp_file)

        # concatenate output files
        output_file_suffixes = ['_D.fasta', '_DL.fasta', '_DLT.fasta', '_DLTP.fasta',
                                '_DLP.fasta', '.gff3', '_statistics.csv']

        for suffix in output_file_suffixes:
            if suffix == '_statistics.csv':
                # sum up line with same word in first column
                stat_files = [output_file + suffix for output_file in output_files]
                new_statistics = sum_up_stats_files(stat_files)
                with open(args.output + suffix, 'w') as f:
                    f.write("\n".join(new_statistics))
                # remove parsed temporary statistics files
                for file in stat_files:
                    os.remove(file)
            else:
                with open(args.output + suffix, 'w') as f:
                    for i in range(number_of_temp_files):
                        # some file may not exist, so we need to check
                        try:
                            with open(output_files[i] + suffix, 'r') as g:
                                for line in g:
                                    f.write(line)
                            # remove parsed temporary output files
                            os.remove(output_files[i])
                        except FileNotFoundError:
                            pass
    else:
        # no need to split sequences into chunks
        subprocess.check_call(
            [f'{tool_path}/detect_putative_ltr.R', '-s', args.reference_sequence, '-g',
             args.gff3, '-o', args.output, '-c', str(args.cpu), '-M',
             str(args.max_missing_domains), '-L', str(args.min_relative_length)]
            )


if __name__ == '__main__':
    # check version of python must be 3.6 or greater
    if sys.version_info < (3, 6):
        print('Python version must be 3.6 or greater')
        sys.exit(1)
    main()
