#!/usr/bin/env python3

"""Integration Tests

Performs a series of integration and recursion tests to ensure that the
gibbs_sampler works as a whole and that the program produces consistent results.
"""

import sys
from collections import Counter
import pytest


@pytest.fixture()
def get_sampler(monkeypatch):
    monkeypatch.setattr("sys.argv",["gibbs_sampler", "inputs/test_seqs.fa", "10"])
    from gibbs_sampler import gibbs_sampler as gb
    global W
    W = int(sys.argv[2])
    return gb


@pytest.fixture()
def get_vars(get_sampler):
    global sequences
    sequences = get_sampler.get_sequences(sys.argv[1])
    global motif_instances
    motif_instances = sequences[0].instances
    global background
    background = get_sampler.background_freq(sequences)
    return sequences, motif_instances, background


def test_sampler(get_sampler, get_vars):
    motif_counts = []
    pwm_counts = [[0]*10]*4
    s, m, b = get_vars
    for _ in range(100):
        for _ in range(1000):
            new_m, old_m, new_p, old_p = get_sampler.run_sampler(s, m, b)
        motif_counts.append(new_m)
        for i in range(4):
            counts = [x + y for x, y in zip(pwm_counts[i], new_p[i])]
            pwm_counts[i] = counts

    pwm_freq = []
    for j in range(4):
        freq = list(map(lambda x: round(x / 100, 1), [v for v in pwm_counts[j]]))
        pwm_freq.append(freq)
    new_p_array = [[round(new_p[nt][pos], 1) for pos in range(W)] for nt in range(4)]
    assert new_p_array == pwm_freq

    motif_1, motif_2, motif_3 = Counter(motif_counts).most_common(3)
    assert 'TTATCTC' in motif_1 or motif_3 or motif_2
