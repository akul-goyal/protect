#!/usr/bin/env python2.7
# Copyright 2016 Arjun Arkal Rao
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from __future__ import print_function
import collections
import csv
from math import ceil
import os
import re
import string

from bd2k.util.expando import Expando

from protect.common import get_files_from_filestore, docker_call, export_results, untargz


def run_fusion_caller(job, star_bam, univ_options, fusion_options):
    """
    This module will run a fusion caller on DNA bams.  This module will be
    implemented in the future.

    This module corresponds to node 10 on the tree
    """
    job.fileStore.logToMaster('Running FUSION on %s' % univ_options['patient'])
    job.fileStore.logToMaster('Fusions are currently unsupported.... Skipping.')
    fusion_file = job.fileStore.getLocalTempFile()
    output_file = job.fileStore.writeGlobalFile(fusion_file)
    return output_file


def star_fusion_disk(rna_fastqs, star_tar):
    return int(4 * ceil(sum([f.size for f in rna_fastqs]) + 524288) +
               2 * ceil(star_tar.size + 524288) +
               5242880)

def run_star_fusion(job, fastqs, junction_file, univ_options, star_fusion_options):
    """
    Runs STAR-Fusion and filters fusion calls using FusionInspector

    :param job job: job
    :param tuple fastqs: RNA-Seq FASTQ Filestore IDs
    :param FsID junction_file: Chimeric junction file
    :param dict univ_options: universal arguments used by almost all tools
    :param dict star_fusion_options: STAR-Fusion specific parameters
    :return: Transgene BEDPE file
    :rtype: Filestore ID
    """
    job.fileStore.logToMaster('Running STAR-Fusion on %s' % univ_options['patient'])
    work_dir = job.fileStore.getLocalTempDir()
    input_files = {'rna_1.fq.gz': fastqs[0],
                   'rna_2.fq.gz': fastqs[1],
                   'STAR.junction': junction_file,
                   'tool_index.tar.gz': star_fusion_options['tool_index']}

    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)
    input_files['tool_index'] = os.path.basename(untargz(input_files['tool_index.tar.gz'], work_dir))

    cores = star_fusion_options['n']
    parameters = ['--chimeric_junction', 'STAR.junction',
                  '--output_dir', 'output',
                  '--genome_lib_dir', input_files['tool_index'],
                  '--CPU', str(cores)]

    docker_call(tool='star-fusion:protect', tool_parameters=parameters, work_dir=work_dir, dockerhub='jpfeil')

    fusion_path = os.path.join(work_dir, 'output/star-fusion.fusion_candidates.final.abridged')

    # Check for fusion prediction
    with open(fusion_path, 'r') as f:
        # Skip header
        f.next()
        try:
            f.next()
        except StopIteration:
            return

    parameters = ['--fusions', fusion_path,
                  '--genome_lib', input_files['tool_index'],
                  '--left_fq', 'rna_1.fq.gz',
                  '--right_fq', 'rna_2.fq.gz',
                  '--out_dir', 'FusionInspector',
                  '--out_prefix', 'FusionInspector',
                  '--prep_for_IGV',
                  '--include_Trinity',
                  '--CPU', cores]

    docker_call(tool='fusion-inspector:protect', tool_parameters=parameters, work_dir=work_dir, dockerhub='jpfeil')

    found_fusion = False
    fusion_path = os.path.join(work_dir, 'FusionInspector/FusionInspector.fusion_predictions.final.abridged.FFPM')
    output_path = os.path.join(work_dir, 'fusion.final')

    # Remove fusions without at least 0.1 fusion fragments per million reads
    if os.path.exists(fusion_path):
        with open(fusion_path, 'r') as f, open(output_path, 'w') as g:
            g.write(f.next())
            for line in f:
                j_ffpm, s_ffpm = line.split()[-2:]
                if sum([float(j_ffpm), float(s_ffpm)]) > 0.1:
                    found_fusion = True
                    g.write(line)

    if found_fusion:
        fusion_bed_path = os.path.join(work_dir, 'FusionInspector/FusionInspector.bed')
        transcript_path = os.path.join(work_dir, 'FusionInspector/FusionInspector.gmap_trinity_GG.fusions.fasta')
        transcript_gff_path = os.path.join(work_dir, 'FusionInspector/FusionInspector.gmap_trinity_GG.fusions.gff3')

        fusion_annotation = job.fileStore.writeGlobalFile(fusion_bed_path)
        filtered_fusions = job.fileStore.writeGlobalFile(output_path)
        transcripts = job.fileStore.writeGlobalFile(transcript_path)
        transcript_annotation = job.fileStore.writeGlobalFile(transcript_gff_path)

        export_results(job, filtered_fusions, output_path, univ_options, subfolder='mutations/fusions')
        export_results(job, transcripts, transcript_path, univ_options, subfolder='mutations/fusions')
        export_results(job, transcript_annotation, transcript_gff_path, univ_options, subfolder='mutations/fusions')
        return job.addChildJobFn(reformat_star_fusion_output,
                                 fusion_annotation,
                                 filtered_fusions,
                                 transcripts,
                                 transcript_annotation,
                                 univ_options).rv()


def parse_star_fusion(infile):
    """
    Parses STAR-Fusion format and returns an Expando object with basic features

    :param str infile: path to STAR-Fusion prediction file
    :return: Fusion prediction attributes
    :rtype: Expando
    """
    reader = csv.reader(infile, delimiter='\t')
    header = reader.next()
    header = {key: index for index, key in enumerate(header)}

    features = ['LeftGene', 'LeftLocalBreakpoint', 'LeftBreakpoint',
                'RightGene', 'RightLocalBreakpoint', 'RightBreakpoint']

    for line in reader:
        yield Expando(dict((feature, line[header[feature]]) for feature in features))


def get_transcripts(transcript_file):
    """
    Parses FusionInspector transcript file and returns dictionary of sequences

    :param str transcript_file: path to transcript FASTA
    :return: de novo assembled transcripts
    :rtype: dict
    """
    with open(transcript_file, 'r') as fa:
        transcripts = {}
        regex = re.compile(r"(?P<ID>TRINITY.*)\s(?P<fusion>.*--.*):(?P<left_start>\d+)-(?P<right_start>\d+)")
        while True:
            # Usually the transcript is on one line
            try:
                info = fa.next()
                seq = fa.next()

                assert info.startswith('>')

                m = regex.search(info)
                if m:
                    transcripts[m.group('ID')] = seq.strip()

            except StopIteration:
                break

            except AssertionError:
                print("WARNING: Malformed fusion transcript file")
    return transcripts


def split_fusion_transcript(annotation_path, transcripts):
    """
    Finds the breakpoint in the fusion transcript and splits the 5' donor from the 3' acceptor

    :param str annotation_path: Path to transcript annotation file
    :param dict transcripts: Dictionary of fusion transcripts
    :return: 5' donor sequences and 3' acceptor sequences
    :rtype: tuple
    """
    annotation = collections.defaultdict(dict)
    forward = 'ACGTN'
    reverse = 'TGCAN'

    trans = string.maketrans(forward, reverse)

    # Pull in assembled transcript annotation
    five_pr_splits = collections.defaultdict(dict)
    three_pr_splits = collections.defaultdict(dict)

    regex = re.compile(r'ID=(?P<ID>.*);Name=(?P<Name>.*);Target=(?P<Target>.*)\s(?P<start>\d+)\s(?P<stop>\d+)')

    # print('SLEEPING')
    # import time
    # time.sleep(1800)

    with open(annotation_path, 'r') as gff:
        for line in gff:
            print(line)
            if line.startswith('#'):
                _, eyd, fusion = line.strip().split()
                fusion, start_stop = fusion.split(':')
                left_break, right_break = start_stop.split('-')

                annotation[fusion]['ID'] = eyd
                annotation[fusion]['left_break'] = left_break
                annotation[fusion]['right_break'] = right_break

            else:
                line = line.strip().split('\t')
                fusion = line[0]
                strand = line[6]
                block_start = line[3]
                block_stop = line[4]
                attr = line[8]

                if strand == '-' and block_start == annotation[fusion]['right_break']:
                    m = regex.search(attr)
                    if m:
                        transcript_id = m.group('Name')
                        transcript_split = int(m.group('stop'))

                        # Take the reverse complement to orient transcripts from 5' to 3'
                        five_seq = transcripts[transcript_id][transcript_split:].translate(trans)[::-1]
                        five_pr_splits[fusion][transcript_id] = five_seq
                        three_seq = transcripts[transcript_id][:transcript_split].translate(trans)[::-1]
                        three_pr_splits[fusion][transcript_id] = three_seq

                if strand == '+' and block_stop == annotation[fusion]['left_break']:
                    m = regex.search(attr)
                    if m:
                        transcript_id = m.group('Name')
                        transcript_split = int(m.group('stop'))
                        five_pr_splits[fusion][transcript_id] = transcripts[transcript_id][:transcript_split]
                        three_pr_splits[fusion][transcript_id] = transcripts[transcript_id][transcript_split:]

    return five_pr_splits, three_pr_splits


def get_gene_ids(fusion_bed):
    """
    Parses FusionInspector bed file to ascertain the ENSEMBL gene ids

    :param str fusion_bed: path to fusion annotation
    :return:
    """
    with open(fusion_bed, 'r') as f:
        gene_to_id = {}
        regex = re.compile(r'(?P<gene>ENSG\d*)')
        for line in f:
            line = line.split('\t')
            transcript, gene_bit, name = line[3].split(';')
            m = regex.search(gene_bit)
            if m:
                gene_to_id[name] = m.group('gene')
    return gene_to_id



def reformat_star_fusion_output(job, fusion_annot, fusion_file, transcript_file, transcript_gff_file, univ_options):
    """
    Writes STAR-Fusion results in Transgene BEDPE format

    :param job job: job
    :param FsID fusion_file: STAR-fusion prediction file
    :param FsID transcript_file: Fusion transcript FASTA file
    :param FsID transcript_gff_file: Fusion transcript GFF file
    :param univ_options: universal arguments used by almost all tools
    :return: Transgene BEDPE file
    :rtype: FsID
    """
    job.fileStore.logToMaster('Reformatting STAR-Fusion output for %s' % univ_options['patient'])

    input_files = {'results.tsv': fusion_file,
                   'fusion.bed': fusion_annot,
                   'transcripts.fa': transcript_file,
                   'transcripts.gff': transcript_gff_file}

    work_dir = job.fileStore.getLocalTempDir()
    input_files = get_files_from_filestore(job, input_files, work_dir, docker=False)

    # Pull in assembled transcript file
    hugo_to_gene_ids = get_gene_ids(input_files['fusion.bed'])
    transcripts = get_transcripts(input_files['transcripts.fa'])
    five_pr_splits, three_pr_splits = split_fusion_transcript(input_files['transcripts.gff'], transcripts)

    # Pull in assembled transcript annotation

    # Header for BEDPE file
    header = ['# chr1', 'start1', 'end1',
              'chr2', 'start2', 'end2',
              'name', 'score',
              'strand1', 'strand2',
              'junctionSeq1', 'junctionSeq2',
              'hugo1', 'hugo2']

    output_path = os.path.join(work_dir, 'fusion_results.bedpe')
    with open(input_files['results.tsv'], 'r') as in_f, open(output_path, 'w') as out_f:
        writer = csv.writer(out_f, delimiter='\t')
        writer.writerow(header)
        for record in parse_star_fusion(in_f):

            left_chr, left_break, left_strand = record.LeftBreakpoint.split(':')

            right_chr, right_break, right_strand = record.RightBreakpoint.split(':')

            fusion = ''.join([record.LeftGene, '--', record.RightGene])
            name = '-'.join([hugo_to_gene_ids[record.LeftGene], hugo_to_gene_ids[record.RightGene]])

            for transcript_id in five_pr_splits[fusion].keys():
                five_prime_seq = five_pr_splits[fusion][transcript_id]
                three_prime_seq = three_pr_splits[fusion][transcript_id]

                writer.writerow([left_chr,
                                 '.',  # Donor start position is not necessary
                                 left_break,
                                 right_chr,
                                 right_break,
                                 '.',  # Acceptor end position is not necessary
                                 name,
                                 '0',
                                 left_strand,
                                 right_strand,
                                 five_prime_seq,
                                 three_prime_seq,
                                 record.LeftGene,
                                 record.RightGene])
    return job.fileStore.writeGlobalFile(output_path)
