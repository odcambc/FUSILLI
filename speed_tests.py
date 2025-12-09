import timeit
import csv

import pyfastx
import numpy as np
import re
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool
from itertools import product

input_fastq = "/Users/bartleby/Desktop/Projects/fusilli/test_data/small1.fq"
fusion_sequence_csv = "/Users/bartleby/Desktop/Projects/fusilli/test.csv"
output_csv = "/Users/bartleby/Desktop/Projects/fusilli/test_counts_time.csv"

fq = pyfastx.Fastq(input_fastq)

def fusion_sequence_to_dict(fusion_sequence_csv):
    fusion_dict = {}
    with open(fusion_sequence_csv, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            fusion_dict[row[0]] = row[1]
    return fusion_dict


def find_match_str(read_sequence, fusion_sequence):
    match_pos = read_sequence.find(fusion_sequence)
    if match_pos != -1:
        return match_pos
    else:
        return None


def find_match_str_parallel(read_sequence, fusion_tuple):
    fusion_name = fusion_tuple[0]
    fusion_sequence = fusion_tuple[1]
    match_pos = read_sequence.find(fusion_sequence)
    if match_pos != -1:
        return fusion_name
    else:
        return None


def find_match_re(read_sequence, fusion_sequence):
    match_pos = re.search(fusion_sequence, read_sequence)
    if match_pos:
        return match_pos
    else:
        return None


def find_match_re_parallel(read_sequence, fusion_tuple):
    fusion_name = fusion_tuple[0]
    fusion_sequence = fusion_tuple[1]
    match_pos = re.search(fusion_sequence, read_sequence)
    if match_pos:
        return fusion_name
    else:
        return None

def find_match_np(read_sequence, fusion_sequence):
    match_pos = np.char.find(read_sequence, fusion_sequence)
    if match_pos != -1:
        return match_pos
    else:
        return None


def find_match_np_parallel(read_sequence, fusion_tuple):
    fusion_name = fusion_tuple[0]
    fusion_sequence = fusion_tuple[1]
    match_pos = np.char.find(read_sequence, fusion_sequence)
    if match_pos != -1:
        return fusion_name
    else:
        return None


def update_counts(fusion_name, match_counts):
    if fusion_name not in match_counts:
        match_counts[fusion_name] = 1
    else:
        match_counts[fusion_name] += 1


def write_match_counts_to_csv(output_csv, match_counts):
    with open(output_csv, "a") as f:
        writer = csv.writer(f)
        for fusion_name, count in match_counts.items():
            writer.writerow([fusion_name, count])


# Single thread str find


def process_fastq_single_str(fq, fusion_dict, match_counts):
    for sequence in fq:
        for fusion_name, fusion_sequence in fusion_dict.items():
            match_pos = find_match_str(sequence.seq, fusion_sequence)
            if match_pos:
                #update_counts(fusion_name, match_counts)
                continue


def main_str_single():
    match_counts = {}
    fq = pyfastx.Fastx(input_fastq)
    fusion_dict = fusion_sequence_to_dict(fusion_sequence_csv)
    process_fastq_single_str(fq, fusion_dict, match_counts)


# Multi thread str find


def process_fastq_parallel_str(fq, fusion_dict, match_counts):
    for name, seq, qual in fq:
        args = product([seq], fusion_dict.items())
        results = p.starmap(find_match_str_parallel, args)
        for result in results:
            if result:
                fusion_name = result
                #update_counts(fusion_name, match_counts)


def main_str_parallel():
    match_counts = {}
    fq = pyfastx.Fastx(input_fastq)
    fusion_dict = fusion_sequence_to_dict(fusion_sequence_csv)
    process_fastq_parallel_str(fq, fusion_dict, match_counts)


# Single thread re search


def process_fastq_single_re(fq, fusion_dict, match_counts):
    for name, seq, qual in fq:
        for fusion_name, fusion_sequence in fusion_dict.items():
            match_pos = find_match_re(seq, fusion_sequence)
            if match_pos:
                #update_counts(fusion_name, match_counts)
                continue


def main_re_single():
    match_counts = {}
    fq = pyfastx.Fastx(input_fastq)
    fusion_dict = fusion_sequence_to_dict(fusion_sequence_csv)
    process_fastq_single_re(fq, fusion_dict, match_counts)


# Multi thread re search

def process_fastq_parallel_re(fq, fusion_dict, match_counts):
    for name, seq, qual in fq:
        args = product([seq], fusion_dict.items())
        results = p.starmap(find_match_re_parallel, args)
        for result in results:
            if result:
                fusion_name = result
                #update_counts(fusion_name, match_counts)


def main_re_parallel():
    match_counts = {}
    fq = pyfastx.Fastx(input_fastq)
    fusion_dict = fusion_sequence_to_dict(fusion_sequence_csv)
    process_fastq_parallel_re(fq, fusion_dict, match_counts)


# Single thread numpy search


def process_fastq_single_np(fq, fusion_dict, match_counts):
    for name, seq, qual in fq:
        for fusion_name, fusion_sequence in fusion_dict.items():
            match_pos = find_match_np(seq, fusion_sequence)
            if match_pos:
                #update_counts(fusion_name, match_counts)
                continue


def main_np_single():
    match_counts = {}
    fq = pyfastx.Fastx(input_fastq)
    fusion_dict = fusion_sequence_to_dict(fusion_sequence_csv)
    process_fastq_single_np(fq, fusion_dict, match_counts)


# Multi thread numpy search


def process_fastq_parallel_np(fq, fusion_dict, match_counts):
    for name, seq, qual in fq:
        args = product([seq], fusion_dict.items())
        results = p.starmap(find_match_np_parallel, args)
        for result in results:
            if result:
                fusion_name = result
                update_counts(fusion_name, match_counts)


def main_np_parallel():
    match_counts = {}
    fq = pyfastx.Fastx(input_fastq)
    fusion_dict = fusion_sequence_to_dict(fusion_sequence_csv)
    process_fastq_parallel_np(fq, fusion_dict, match_counts)


if __name__ == "__main__":
    print("Single thread str search time:")
    print(timeit.timeit(main_str_single, number=1))

    with Pool(processes=32) as p:
        print("Parallel str search time:")
        print(timeit.timeit(main_str_parallel, number=1))

    with ThreadPool(processes=32) as p:
        print("Parallel str search time threadpool:")
        print(timeit.timeit(main_str_parallel, number=1))

    print("Single thread re search time:")
    print(timeit.timeit(main_re_single, number=1))

    with Pool(processes=32) as p:
        print("Parallel re search time:")
        print(timeit.timeit(main_re_parallel, number=1))

    with ThreadPool(processes=32) as p:
        print("Parallel re search time threadpool:")
        print(timeit.timeit(main_re_parallel, number=1))

    print("Single thread np search time:")
    print(timeit.timeit(main_np_single, number=1))

    with Pool(processes=32) as p:
        print("Parallel np search time:")
        print(timeit.timeit(main_np_parallel, number=1))

    with ThreadPool(processes=32) as p:
        print("Parallel np search time threadpool:")
        print(timeit.timeit(main_np_parallel, number=1))
