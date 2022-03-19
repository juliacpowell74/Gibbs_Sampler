#!/usr/bin/env python3

"""Unit tests for the gibbs_sampler program."""


import sys
import pytest
from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq
import os
import random
import re
import subprocess


prg = '../src/gibbs_sampler/gibbs_sampler.py'


def test_exists():
    assert os.path.exists(prg)


def test_usage():
    output = subprocess.Popen(['gibbs_sampler','-h'],stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT)
    outs, errs = output.communicate()
    assert re.match(b'usage',outs,re.IGNORECASE)


def test_bad_width_int():
    width = random.randint(-10, 0)
    output = subprocess.Popen(['gibbs_sampler', 'inputs/test_seqs.fa', f'{width}'],
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT)
    outs, errs = output.communicate()
    assert re.search(f"width {width} must be greater than 0.", str(outs))


def test_bad_num_iter():
    num_iter = random.randint(-10,0)
    output = subprocess.Popen(['gibbs_sampler', 'inputs/test_seqs.fa', '10',
                               f'-n {num_iter}'], stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT)
    outs, errs = output.communicate()
    assert re.search(f'num_iter {num_iter} must be greater than 0.', str(outs))


@pytest.fixture
def get_sampler(monkeypatch):
    monkeypatch.setattr("sys.argv",["gibbs_sampler", "inputs/test_seqs.fa", "4"])
    from gibbs_sampler import gibbs_sampler as gb
    global W
    W = int(sys.argv[2])
    return gb


def test_class_exists(get_sampler):
    sequences = [get_sampler.Sequence(record.id, str(record.seq))
                 for record in SeqIO.parse(sys.argv[1], "fasta")]
    assert sequences[0].name == 'mm9_chr11:95292903-95292913(+)'
    assert len(sequences) == 26


def test_class_init_site(get_sampler):
    sequences = get_sampler.get_sequences(sys.argv[1])
    assert len(sequences[0].instances[0]) == W
    assert(len(sequences[0].instances)) == 52



def test_class_init_error(get_sampler):
    with pytest.raises(SystemExit):
        W = 10
        return get_sampler.Sequence('short_seq', 'AGCT')


@pytest.mark.parametrize('file, exception',
                         [('inputs/test_seqs.fa', None),
                          ('inputs/single_seq.fa',
                           SystemExit('Error: More than 1 sequence must be provided.')),
                          ('inputs/bad_file.embl',
                           TypeError("input file is not in fasta format.")
                         )]
                         )


def test_get_sequence(get_sampler, file, exception):
    sys.argv[1] = file
    try:
        sequences = get_sampler.get_sequences(sys.argv[1])
        assert sequences[0].sequence == "GTCACTGTGTACTCTAGGCTCGTTGGTCCCC" \
                                        "AAGCTTCTGGGTGGCTCTTTCTTATCTCCCGT" \
                                        "CTTACTGTAAGAACAGATGGAGTGCT" \
                                        "AGAACAAGTAGGATTGTGTCTG"
        assert type(sequences[0]) == get_sampler.Sequence
    except (SystemExit, TypeError) as error:
        assert isinstance(error, type(exception))
        assert error.args == exception.args


def test_background_freq(get_sampler):
    sequences = get_sampler.get_sequences(sys.argv[1])
    background = get_sampler.background_freq(sequences)
    assert round(sum(background.values()), 2) == 1.0
    expected = {'A': 0.227997227997228, 'C': 0.28586278586278585,
                'G': 0.21067221067221067, 'T': 0.27546777546777546}
    assert background == expected


def test_calc_pwm(get_sampler):
    motif_instances = get_sampler.Sequence.instances
    pwm = get_sampler.calc_pwm(motif_instances)
    for pos in range(W):
        pos_sum = 0
        for nt in range(4):
            pos_sum += pwm[nt][pos]
        assert pytest.approx(pos_sum, 0.1) == 1


@pytest.fixture
def pwm_test_sets():
    inst = [Seq('GTAA'), Seq('ACGT'),Seq('GCGT'),Seq('ACAT'),Seq('GCGA')]
    inst_small = [Seq('ATAA'), Seq('ACGT'),Seq('GCGT'),Seq('ACAT'),Seq('GCGA')]
    inst_big = [Seq('GTTA'), Seq('ATGT'),Seq('GCAT'),Seq('TCAT'),Seq('CCGC')]

    mot, mot_small, mot_big = [motifs.create(seqs) for seqs in (inst, inst_small, inst_big)]
    pwm, pwm_small, pwm_big = [motif.counts.normalize(pseudocounts=1)
                               for motif in (mot, mot_small, mot_big)]

    return pwm, pwm_small, pwm_big


def test_change_in_pwm(get_sampler, pwm_test_sets):
    pwm, small_change, big_change = pwm_test_sets
    assert get_sampler.change_in_pwm(pwm, pwm) is True
    assert get_sampler.change_in_pwm(pwm, small_change) is True
    assert get_sampler.change_in_pwm(pwm, big_change) is False


def test_change_in_motif(get_sampler):
    motif = "ACGT"
    wrong_motif = "GCGT"
    assert get_sampler.change_in_motif(motif, wrong_motif) is False
    assert get_sampler.change_in_motif(motif, motif) is True




















