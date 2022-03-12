#!/usr/bin/env python3

""" Gibbs Sampler

This script runs a basic Gibbs Sampler algorithm for de novo motif
finding. Given a set of sequences and motif width, the algorithm will
run for the desired number of iterations. It will then check for convergence
as defined as no change in motif position or less than a one percent change
in the Position Weight Matrix. If convergence has not been reached, the
sampler will continue to run until it reaches convergence or until it
has repeated the desired number of iterations for second a time to
prevent the sampler from running indefinitely.

Input files must be in fasta format.

This script requires that the third-party package 'Biopython' be
installed within the python working environment.
Installation of the command line tool 'Weblogo3' is also required.

Output is writen to the console and includes the the final motif produced
by the algorithm, its corresponding Position Specific Scoring Matrix,
and the relative entropy of the motif instance. A weblogo is also generated
in eps file format.
"""

import sys
import subprocess
import argparse
import random
import itertools
from Bio import SeqIO
from Bio import motifs

# --------------------------------------------------


def get_args(args=None):
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file',
                        help='Input sequence in fasta file format',
                        type=argparse.FileType('rt'))
    parser.add_argument('width',
                        help='width of motif',
                        type=int)
    parser.add_argument('-t',
                        '--title',
                        help='name of output weblogo file',
                        metavar='title',
                        type=argparse.FileType('wt'),
                        default='WebLogo')
    parser.add_argument('-n',
                        '--num_iter',
                        help='initial number of iterations to run algorithm',
                        metavar='num_iter',
                        type=int,
                        default=1000)

    argv = parser.parse_args(args)

    if argv.width < 1:
        print(parser.error(f'width {argv.width} must be greater than 0.'))
    if argv.num_iter < 1:
        print(parser.error(f'num_iter {argv.num_iter} must be greater than 0.'))

    return argv


# --------------------------------------------------

# initialize global variables received from argparse
args = get_args()
logo_title = args.title
num_iter = args.num_iter
W = args.width
fasta_file = args.file



class Sequence:
    """
    Create a 'Sequence' object.

    Parses a given sequence and initializes a motif instance from
    the sequence.
    Contains additional methods helpful for running the Gibbs algorithm.

    Parameters
    ----------
    name : str
        Sequence identity.
    sequence : str
        The DNA sequence.
    """
    instances = []

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        self.init_site()

    def __str__(self):
        return f'{self.name},{self.sequence}'

    def init_site(self):
        """
        Create a motif instance and appends it to the 'instance' attribute.

        Raises
        ------
        SystemExit
            If the desired motif width is equal to or longer than the sequence.
        AssertionError
            If the motif instance is longer than the desired width.
        """
        if W >= len(self.sequence):
            sys.exit('Error: motif width must be less than sequence length')
        pos = random.randint(0, len(self.sequence) - W)
        motif = self.sequence[pos:pos + W]
        assert len(motif) == W
        self.instances.append(motif)

    def find_site_prob(self, pwm, background):
        """
        Find the probability of the motif at each position.

        Passed as an argument to sample_new_motif method.

        Parameters
        ----------
        pwm : Bio.motifs.matrix.PositionWeightMatrix
            A position weight matrix generated from the motif instances.
        background : dict of {str : float}
            Background nucleotide frequencies.

        Returns
        -------
        motif_score : list of float
            Contains the motif score at each position in the sequence.
        total_score : float
            The sum of each motif score.
        """
        site = 0
        motif_score = []
        total_score = 0
        while site in range(0, len(self.sequence) - W + 1):
            site_score = 1
            index = 0
            for nt in self.sequence[site:site + W]:
                nt_score = pwm[nt][index] / background[nt]
                site_score *= nt_score
                index += 1
            motif_score.append(site_score)
            total_score += site_score
            site += 1
        return motif_score, total_score

    def sample_new_motif(self, find_site_prob, pwm, background):
        """
        Select a new motif instance.

        The selection is random but weighted by the probability of
        the motif occurring at each position.
        Probabilities are calculated using the find_site_prob() method.

        Parameters
        ----------
        find_site_prob : Sequence class method
            Method for determining the motif site probability in a sequence.
        pwm, background : See find_site_prob()
            Arguments passed into the find_site_prob() method.

        Returns
        -------
        new_motif : str
            A new motif instance.

        Raises
        ------
        AssertionError
            If normalized motif scores do not sum to 1.
         """
        motif_score, total_score = find_site_prob(pwm, background)
        weights = []
        for score in motif_score:
            norm_score = score / total_score
            weights.append(norm_score)
        assert round(sum(weights), 2) == 1.0
        pos = random.choices(
            range(0, len(self.sequence) - W + 1),
            weights=weights, k=1)
        new_motif = self.sequence[pos[0]:pos[0] + W]
        return new_motif

# -------------------------------------------------------------------------------------


def get_sequences(file):
    """
    Parse each sequence in a fasta file and create a Sequence object.

    Parameters
    ----------
    file : io.TextIOWrapper
        File handle containing DNA sequences in fasta format.

    Returns
    -------
    sequences : list of Sequence objects
        A list containing all parsed sequences.

    Raises
    ------
    SystemExit
        If one or no sequences were found in the file.
    TypeError
        If input file of sequences is not in fasta format.
    """
    sequences = [Sequence(record.id, str(record.seq))
                 for record in SeqIO.parse(file, "fasta")]
    if sequences == []:
        raise TypeError("input file is not in fasta format.")
    if len(sequences) < 2:
        sys.exit('Error: More than 1 sequence must be provided.')
    return sequences


def background_freq(sequences):
    """
    Find the total nucleotide frequencies from all sequences.

    Parameters
    ----------
    sequences : list of Sequence objects
        A list containing the DNA sequences.
    Returns
    -------
    background : dict of {str : float}
        The background frequency for each nucleotide.
    Raises
    ------
    AssertionError
        If  frequencies do not sum to one.
    """
    background_dict = {'A': 0, "C": 0, "G": 0, "T": 0}
    for seq in sequences:
        background_dict['A'] += seq.sequence.count('A')
        background_dict['C'] += seq.sequence.count('C')
        background_dict['G'] += seq.sequence.count('G')
        background_dict['T'] += seq.sequence.count('T')
    total = sum(background_dict.values())
    background = {k: v / total for (k, v) in background_dict.items()}
    assert round(sum(background.values()), 2) == 1.0
    return background


def calc_pwm(motif_instances):
    """
    Calculate position weight matrix from all motif instances.

    Parameters
    ----------
    motif_instances : list of str
        A list containing randomly selected motif instance from each sequence.

    Returns
    -------
    pwm : Bio.motifs.matrix.PositionWeightMatrix
        A position weight matrix generated from motif instances.
    """
    motif = motifs.create(motif_instances)
    pwm = motif.counts.normalize(pseudocounts=1)
    return pwm


def run_sampler(sequences, motif_instances, background):
    """
    The Gibbs Sampler algorithm.

    Iterates through the sampling, predictive update, and likelihood
    calculation steps. During each iteration except the last, the
    current pwm and motif instance are saved to test for convergence
    after the final iteration of the algorithm.

    Parameters
    ----------
    sequences : list of str
        List of all sequences  contained in fasta file.
    motif_instances : list of str
        List of all motif instances generated from the Sequence class.
    background : dict of {str : float}
        Background nucleotide frequencies.

    Returns
    -------
    new_motif : str
        The final motif generated from the sampler.
    old_motif : str
        The second to last motif generated.
    new_pwm : Bio.motifs.matrix.PositionWeightMatrix
        The final pwm generated from the sampler.
    old_pwm : Bio.motifs.matrix.PositionWeightMatrix
        The second to last pwm generated from the sampler.
    """
    # initialize starting positions in sequences
    # completed when a Sequence object was created for each sequence

    # randomly leaver out one sequence
    leave_out_idx = random.randrange(0, len(sequences))
    leave_out_seq = sequences[leave_out_idx]
    del motif_instances[leave_out_idx]

    # calculate pwm of all but left-out sequence
    pwm = calc_pwm(motif_instances)

    # calculate probability of motif occurring at each site of the
    # left-out sequence and sample a new motif position in the sequence
    new_motif = leave_out_seq.sample_new_motif(
        leave_out_seq.find_site_prob,
        pwm, background)
    motif_instances.insert(leave_out_idx, new_motif)

    for _ in range(num_iter-1):
        # store value of current pwm and motif
        old_pwm = pwm
        old_motif = new_motif

    return new_motif, old_motif, pwm, old_pwm


def change_in_pwm(final_pwm, initial_pwm):
    """
    Calculate the percent change between the ultimate and penultimate pwm.

    A test for convergence. Takes the final two matrices
    generated from the run_sampler function and compares the values
    from each column in the matrices.

    Parameters
    ----------
    final_pwm : Bio.motifs.matrix.PositionWeightMatrix
        The final pwm generated from run_sampler().
    initial_pwm : Bio.motifs.matrix.PositionWeightMatrix
        The second to last pwm generated from run_sampler().

    Returns
    -------
    False : bool
        If the percent change between the final and initial pwm is
        greater than or equal to one percent.
    True : bool
        If the percent change is less than one percent.
    """
    def percent_change(final, initial):
        """Simple function for calculating percent change."""

        return abs((final-initial) / initial)
    pos = 0
    while pos < W:
        pwm_change = 0
        for nt in final_pwm:
            pwm_change = (
                percent_change(final_pwm[nt][pos], initial_pwm[nt][pos])
                if initial_pwm[nt][pos] != 0
                else final_pwm[nt][pos]
            )
        if pwm_change!=0 and round(pwm_change, 2)>=0.01:
            return False
        pos += 1
    return True


def change_in_motif(final_motif, initial_motif):
    """
    Compare the ultimate and penultimate motif instances from run_sampler.

    A test for convergence. Determines if any motif subsequence location
    has changed in the final two iterations of the sampler.

    Parameters
    ----------
    final_motif : str
        The final motif generated from run_sampler().
    initial_motif
        The second to last motif instance generated.

    Returns
    -------
    False : bool
        If the two motif instances do not match.
    True : bool
        If the motif instances do match.
    """
    for i, _ in enumerate(final_motif):
        if final_motif[i] != initial_motif[i]:
            return False
    return True


def web_logo(motif_instances):
    """
    Create a weblogo from the motif instances generated by run_sampler.

    Runs the  program 'WebLogo' by executing a terminal command
    using subprocess.
    The input is the list of motif instances which are writen into stdin
    and the result is written to an output file.

    Parameters
    ----------
    motif_instances : list of str
        A list containing the motif instances generated from run_sampler()

    Raises
    ------
    AssertionError
        If process has not finished.
    """
    with subprocess.Popen(
            ['weblogo'],
            stdin=subprocess.PIPE,
            stdout=logo_title,
            stderr=subprocess.STDOUT) as proc:
        content = [f'>{motif_instances.index(i)}\n{i}\n'
                   for i in motif_instances]
        proc.stdin.write("".join(content).encode('utf-8'))
    assert proc.poll() is not None


def main():
    sequences = get_sequences(fasta_file)
    background = background_freq(sequences)
    motif_instances = Sequence.instances

    # run the Gibbs algorithm
    for _ in range(num_iter):
        final_motif, initial_motif, final_pwm, initial_pwm = run_sampler(
            sequences, motif_instances, background)

    # check for convergence after num_iter repetitions of sampler
    calls = 0
    while (change_in_pwm(final_pwm, initial_pwm) is None or
            change_in_motif(final_motif, initial_motif) is None):
        # repeat sampler until one condition returns True
        itertools.repeat(run_sampler(sequences, motif_instances, background), 1)
        calls += 1
        if calls > num_iter:
            # exit the program if neither condition returns True after
            # running the sampler a second time through the
            # desired number of iterations
            # prevents the loop from running indefinitely
            sys.exit('Error: No motif was found in the provided sequences.')

    pssm = final_pwm.log_odds(background)
    rel_entropy = pssm.mean(background)
    web_logo(motif_instances)

    print(f'The PSSM is:\n{pssm}\n\nThe motif is {final_motif}.\n\n'
          f'The relative entropy is {rel_entropy} bits.')

# --------------------------------------------------


if __name__ == '__main__':
    main()
