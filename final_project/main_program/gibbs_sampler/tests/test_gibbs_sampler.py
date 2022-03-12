import argparse
import sys
import pytest
from subprocess import getstatusoutput
import os
from gibbs_sampler import gibbs_sampler


prg = '/Users/juliapowell/BB485/final_project/main_program/old_gibbs_sampler.py'


def test_exists():
    assert os.path.isfile(prg)

def test_output():
    rv, out = getstatusoutput(f'{prg} -h')
    assert rv == 0

sys.path.insert(0,'/Users/juliapowell/BB485/final_project/main_program/')

@pytest.fixture
def set_env(monkeypatch):
    monkeypatch.setattr("sys.argv",["gibbs_sampler","test_seqs.fa","7"])
    import gibbs_sampler as gb
    return gb

def test_get(set_env):
    sequences = set_env.get_sequences
    result = sequences(sys.argv[1])
    assert result[0].sequence == "GTCACTGTGTACTCTAGGCTCGTTGGTCCCCAAGCTTCTGGGT" \
                                 "GGCTCTTTCTTATCTCCCGTCTTACTGTAAGAACAGATGGAGTGCT" \
                                 "AGAACAAGTAGGATTGTGTCTG"
def test_get_error(set_env):
    sequences = set_env.get_sequences(sys.argv[1])
    background = set_env.background_freq(sequences)
    assert round(sum(background.values()), 2) == 1.0














