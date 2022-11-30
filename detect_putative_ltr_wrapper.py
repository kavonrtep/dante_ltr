#!/usr/bin/env python
"""This wrapper is intended to be used on large genomes and large DANTE input to
minimize memory usage, It splits input files to pieces and analyze it on by one by
detect_putative_ltr.R
If input does not exceed specified max-chunk_size, it will run detect_putative_ltr.R
directly
"""

import argparse
import os
import sys
import tempfile
from itertools import cycle
import subprocess


class Gff3Feature:
    """
    Class for gff3 feature
    """

    def __init__(self, line):
        self.line = line
        self.items = line.strip().split('\t')
        self.header = self.items[0]
        self.source = self.items[1]
        self.type = self.items[2]
        self.start = int(self.items[3])
        self.end = int(self.items[4])
        self.score = self.items[5]
        self.strand = self.items[6]
        self.frame = self.items[7]
        self.attributes = self.items[8]
        self.attributes_dict = {}
        for item in self.attributes.split(';'):
            if item != '':
                key, value = item.split('=')
                self.attributes_dict[key] = value

        self.attributes_str = ';'.join(
            ['{}={}'.format(key, value) for key, value in self.attributes_dict.items()]
            )

    def __str__(self):
        return '\t'.join(
            [self.header, self.source, self.type, str(self.start), str(self.end),
             self.score, self.strand, self.frame, self.attributes_str]
            ) + '\n'

    def __repr__(self):
        return '\t'.join(
            [self.header, self.source, self.type, str(self.start), str(self.end),
             self.score, self.strand, self.frame, self.attributes_str]
            ) + '\n'

    def __eq__(self, other):
        return self.line_recalculated() == other.line_recalculated()

    def __hash__(self):
        return hash(self.line_recalculated())

    def get_line(self):
        """returns original line"""
        return self.line

    def overlap(self, other):
        """
        Check if two features overlap
        :param other:
        :return:
        """
        if self.start <= other.end and self.end >= other.start:
            return True
        else:
            return False

    def line_recalculated(self):
        """
        :return:
        string with recalculated line
        """
        return '\t'.join(
            [self.header, self.source, self.type, str(self.start), str(self.end),
             self.score, self.strand, self.frame, self.attributes_str]
            ) + '\n'

    def __lt__(self, other):
        width = self.end - self.start
        other_width = other.end - other.start
        return width < other_width

    def __gt__(self, other):
        width = self.end - self.start
        other_width = other.end - other.start
        return width > other_width

    def identical_region(self, other):
        """
        Check if two features are identical
        :param other:
        :return:
        """
        if self.start == other.start and self.end == other.end and self.header == \
                other.header:
            return True
        else:
            return False


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
        '-g', '--gff3', default=None, required=True,
        help=("gff3 file with full output from Domain Based Annotation of Transposable "
             "Elements (DANTE)"),

        type=str,
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
        help=('If size of reference sequence is greater than this value, reference is '
             'analyzed in chunks of this size. default is %(default)s '
             'Setting this value too small  will slow down the analysis')
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
        os.remove(temp_files[-1])
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


def read_single_fasta_to_dictionary(fh):
    """
    Read fasta file into dictionary
    :param fh:
    :return:
    fasta_dict
    """
    fasta_dict = {}
    for line in fh:
        if line[0] == '>':
            header = line.strip().split(' ')[0][1:]  # remove part of name after space
            fasta_dict[header] = []
        else:
            fasta_dict[header] += [line.strip()]
    fasta_dict = {k: ''.join(v) for k, v in fasta_dict.items()}
    return fasta_dict


def split_fasta_to_chunks(fasta_file, chunk_size=100000000, overlap=100000):
    """
    Split fasta file to chunks, sequences longe than chuck size are split to overlaping
    peaces. If sequences are shorter, chunck with multiple sequences are created.
    :param fasta_file:

    :param fasta_file:
    :param chunk_size:
    :param overlap:
    :return:
    fasta_file_split
    matching_table
    """
    min_chunk_size = chunk_size * 2
    fasta_dict = read_fasta_sequence_size(fasta_file)
    # calculates ranges for splitting of fasta files and store them in list
    matching_table = []
    fasta_file_split = tempfile.NamedTemporaryFile(delete=False).name
    for header, size in fasta_dict.items():
        if size > min_chunk_size:
            number_of_chunks = int(size / chunk_size)
            adjusted_chunk_size = int(size / number_of_chunks)
            for i in range(number_of_chunks):
                start = i * adjusted_chunk_size
                end = ((i + 1) *
                       adjusted_chunk_size
                       + overlap) if i + 1 < number_of_chunks else size
                new_header = header + '_' + str(i)
                matching_table.append([header, i, start, end, new_header])
        else:
            new_header = header + '_0'
            matching_table.append([header, 0, 0, size, new_header])
    # read sequences from fasta files and split them to chunks according to matching table
    # open output and input files, use with statement to close files
    fasta_dict = read_single_fasta_to_dictionary(open(fasta_file, 'r'))
    with open(fasta_file_split, 'w') as fh_out:
        for header in fasta_dict:
            matching_table_part = [x for x in matching_table if x[0] == header]
            for header2, i, start, end, new_header in matching_table_part:
                fh_out.write('>' + new_header + '\n')
                fh_out.write(fasta_dict[header][start:end] + '\n')
    return fasta_file_split, matching_table


def get_new_header_and_coordinates(header, start, end, matching_table):
    """
    Get new header and coordinates for sequence
    :param header:
    :param start:
    :param end:
    :param matching_table:
    :return:
    new_header
    new_start
    new_end
    """
    matching_table_part = [x for x in matching_table if x[0] == header]
    new_coords = []
    for chunk in matching_table_part:
        if chunk[2] <= start < chunk[3]:
            new_header = chunk[4]
            new_start = start - chunk[2]
            new_end = end - chunk[2]
            new_sequence_length = chunk[3] - chunk[2]
            new_coords.append([new_header, new_start, new_end, new_sequence_length])
    return new_coords


def get_original_header_and_coordinates(new_header, new_start, new_end, matching_table):
    """
    Get original header and coordinates for sequence
    :param new_header:
    :param new_start:
    :param new_end:
    :param matching_table:
    :return:
    original_header
    original_start
    original_end
    """
    matching_table_part = [x for x in matching_table if x[4] == new_header]
    ori_header = matching_table_part[0][0]
    start = matching_table_part[0][2]
    ori_start = new_start + start
    ori_end = new_end + start
    return ori_header, ori_start, ori_end


# recalculate gff3 coordinates, use gff3_feature class
def recalculate_gff3_coordinates(gff3_file, matching_table):
    """
    Recalculate gff3 coordinates, use gff3_feature class
    :param gff3_file:
    :param matching_table:
    :return:
    gff3_file_recalculated
    """
    gff3_file_recalculated = tempfile.NamedTemporaryFile(delete=False).name

    with open(gff3_file, 'r') as fh_in:
        with open(gff3_file_recalculated, 'w') as fh_out:
            for line in fh_in:
                if line[0] == '#':
                    fh_out.write(line)
                else:
                    feature = Gff3Feature(line)
                    new_coords = get_new_header_and_coordinates(
                        feature.header, feature.start, feature.end, matching_table
                        )
                    for new_header, new_start, new_end, sequence_length in new_coords:
                        if new_start >= 1 and new_end <= sequence_length:
                            feature.header = new_header
                            feature.start = new_start
                            feature.end = new_end
                            fh_out.write(str(feature))
    return gff3_file_recalculated


# recalculate gff3 back to original coordinates, use gff3_feature class
def recalculate_gff3_back_to_original_coordinates(gff3_file, matching_table):
    """
    Recalculate gff3 back to original coordinates, use gff3_feature class
    :param gff3_file:
    :param matching_table:
    :return:
    gff3_file_recalculated
    """
    gff3_file_recalculated = tempfile.NamedTemporaryFile(delete=False).name
    with open(gff3_file, 'r') as fh_in:
        with open(gff3_file_recalculated, 'w') as fh_out:
            for line in fh_in:
                if line[0] == '#':
                    fh_out.write(line)
                else:
                    feature = Gff3Feature(line)
                    ori_header, ori_start, ori_end = get_original_header_and_coordinates(
                        feature.header, feature.start, feature.end, matching_table
                        )
                    feature.header = ori_header
                    feature.start = ori_start
                    feature.end = ori_end
                    fh_out.write(str(feature))
    return gff3_file_recalculated


def get_feature_attributes(line):
    """
    Get attributes as dictionary from gff3 list
    :param line:
    :return:
    attributes_dict
    """
    attributes_dict = {}
    for item in line[8].split(';'):
        if item.strip():
            key, value = item.strip().split('=')
            attributes_dict[key] = value
    return attributes_dict


def get_unique_features(gff3_file):
    """
    return list of ID of non-ovelaping features.

    :param gff3_file:
    :return:
    duplicated_features
    """
    good_id = []
    feature_list = []
    with open(gff3_file, 'r') as fh:
        for line in fh:
            if line[0] == '#':
                continue
            feature = Gff3Feature(line)
            if feature.type != 'transposable_element':
                continue
            feature_list.append(feature)  # sort by start position and header
    feature_list.sort(key=lambda x: (x.start, x.header))
    i = 0
    while i < len(feature_list) - 1:
        ch = feature_list[i].header == feature_list[i + 1].header
        if not ch:
            good_id.append(feature_list[i].attributes_dict['ID'])
            i += 1
            continue
        if feature_list[i].identical_region(feature_list[i + 1]):
            # identical position
            good_id.append(feature_list[i].attributes_dict['ID'])
            i += 2
            continue
        if feature_list[i].overlap(feature_list[i + 1]):
            # overlap
            if feature_list[i] > feature_list[i + 1]:
                good_id.append(feature_list[i].attributes_dict['ID'])
                i += 2
                continue
            else:
                good_id.append(feature_list[i + 1].attributes_dict['ID'])
                i += 2
                continue
        else:
            good_id.append(feature_list[i].attributes_dict['ID'])
            i += 1
    if i == len(feature_list) - 1:
        good_id.append(feature_list[i].attributes_dict['ID'])
    return good_id


def filter_gff3_file(gff3_file, good_id, output_file):
    """
    Filter gff3 file by good_id
    :param gff3_file:
    :param good_id:
    :param output_file:
    :return:
    filtered_gff3_file
    """

    with open(gff3_file, 'r') as fh, open(output_file, 'w') as fout:
        for line in fh:
            if line[0] == '#':
                fout.write(line)
            else:
                feature = Gff3Feature(line)
                if ('ID' in feature.attributes_dict and
                        feature.attributes_dict['ID'] in good_id):
                    fout.write(line)
                    continue
                if 'Parent' in feature.attributes_dict:
                    if feature.attributes_dict['Parent'] in good_id:
                        fout.write(line)
                        continue
def verify_seqnames(gff3, fasta):
    """
    Verify that ate least some seqnames in gff3 file are in fasta file
    :param gff3:
    :param fasta:
    :return: True/False
    """
    gff_seqnames = set()
    fasta_seqnames = set()
    with open(gff3, 'r') as fh:
        for line in fh:
            if line[0] == '#':
                continue
            gff_seqnames.add(line.split('\t')[0])
    with open(fasta, 'r') as fh:
        for line in fh:
            if line[0] == '>':
                # seqname in fasta could contain additional description after space
                fasta_seqnames.add(line[1:].split()[0])
    if gff_seqnames.intersection(fasta_seqnames):
        return True
    else:
        return False


def main():
    """
    Main function
    """
    args = get_arguments()
    # locate directory of current script
    tool_path = os.path.dirname(os.path.realpath(__file__))
    # verify that at least some seqnames in gff3 match seqnames in reference fasta file.
    # If not, exit
    if not verify_seqnames(args.gff3, args.reference_sequence):
        # print error message to stderr
        sys.stderr.write(
            'ERROR: No seqnames in gff3 file match seqnames in reference fasta file. '
            'Please check that gff3 file was created from reference fasta file.'
            )
        sys.exit(1)


    # split fasta file to chunks
    fasta_file_split, matching_table = split_fasta_to_chunks(
        args.reference_sequence, chunk_size=args.max_chunk_size
        )
    # recalculate gff3 coordinates
    gff3_file_recalculated = recalculate_gff3_coordinates(args.gff3, matching_table)

    fasta_seq_size = read_fasta_sequence_size(fasta_file_split)
    total_size = sum(fasta_seq_size.values())
    number_of_sequences = len(fasta_seq_size)
    if total_size > args.max_chunk_size and number_of_sequences > 1:
        print('running analysis on chunks of ~ {} Mb'.format(args.max_chunk_size))
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
        with open(fasta_file_split, 'r') as f:
            for line in f:
                if line[0] == '>':
                    header = line.strip().split(' ')[0][1:]
                    seq_id_file_handle_dict[header].write(line)
                else:
                    seq_id_file_handle_dict[header].write(line)
        os.remove(fasta_file_split)
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
        with open(gff3_file_recalculated, 'r') as f:
            for line in f:
                if line[0] == '#':
                    continue
                else:
                    header = line.strip().split('\t')[0]
                    seq_id_file_handle_dict[header].write(line)
        # close file handles
        for file_handle in file_handles:
            file_handle.close()
        os.remove(gff3_file_recalculated)
        # run retrotransposon detection on each temporary file
        output_files = make_temp_files(number_of_temp_files)
        for i in range(number_of_temp_files):
            print('Running retrotransposon detection on file ' + str(i))
            subprocess.check_call(
                [f'{tool_path}/detect_putative_ltr.R', '-s', temp_files_fasta[i], '-g',
                 temp_files_gff[i], '-o', output_files[i], '-c', str(args.cpu), '-M',
                 str(args.max_missing_domains), '-L', str(args.min_relative_length)]
                )

        #remove all temporary input files
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
            elif suffix == '.gff3':
                tmp_gff_unfiltered = tempfile.NamedTemporaryFile(delete=False).name
                with open(tmp_gff_unfiltered, 'w') as f:
                    for i in range(number_of_temp_files):
                        tmp_gff = recalculate_gff3_back_to_original_coordinates(
                            output_files[i] + suffix, matching_table
                            )
                        # remove temporary gff3 file
                        os.remove(output_files[i] + suffix)
                        with open(tmp_gff, 'r') as f_tmp:
                            for line in f_tmp:
                                f.write(line)
                        os.remove(tmp_gff)
                # filter overlapping features
                good_id = get_unique_features(tmp_gff_unfiltered)
                filter_gff3_file(
                    tmp_gff_unfiltered, good_id, args.output + suffix
                    )
                # remove temporary gff3 file
                os.remove(tmp_gff_unfiltered)
            else:
                with open(args.output + suffix, 'w') as f:
                    for i in range(number_of_temp_files):
                        # some file may not exist, so we need to check
                        try:
                            with open(output_files[i] + suffix, 'r') as g:
                                for line in g:
                                    f.write(line)
                            # remove parsed temporary output files
                            os.remove(output_files[i] + suffix)
                        except FileNotFoundError:
                            pass
    else:
        print('running analysis on whole input')
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
