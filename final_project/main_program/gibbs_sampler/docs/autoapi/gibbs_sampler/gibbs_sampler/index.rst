:py:mod:`gibbs_sampler.gibbs_sampler`
=====================================

.. py:module:: gibbs_sampler.gibbs_sampler

.. autoapi-nested-parse::

   Gibbs Sampler

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



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   gibbs_sampler.gibbs_sampler.Sequence



Functions
~~~~~~~~~

.. autoapisummary::

   gibbs_sampler.gibbs_sampler.get_args
   gibbs_sampler.gibbs_sampler.get_sequences
   gibbs_sampler.gibbs_sampler.background_freq
   gibbs_sampler.gibbs_sampler.calc_pwm
   gibbs_sampler.gibbs_sampler.run_sampler
   gibbs_sampler.gibbs_sampler.change_in_pwm
   gibbs_sampler.gibbs_sampler.change_in_motif
   gibbs_sampler.gibbs_sampler.web_logo
   gibbs_sampler.gibbs_sampler.main



Attributes
~~~~~~~~~~

.. autoapisummary::

   gibbs_sampler.gibbs_sampler.args
   gibbs_sampler.gibbs_sampler.logo_title
   gibbs_sampler.gibbs_sampler.num_iter
   gibbs_sampler.gibbs_sampler.W
   gibbs_sampler.gibbs_sampler.fasta_file


.. py:function:: get_args(args=None)

   Get command-line arguments


.. py:data:: args
   

   

.. py:data:: logo_title
   

   

.. py:data:: num_iter
   

   

.. py:data:: W
   

   

.. py:data:: fasta_file
   

   

.. py:class:: Sequence(name, sequence)

   Create a 'Sequence' object.

   Parses a given sequence and initializes a motif instance from
   the sequence.
   Contains additional methods helpful for running the Gibbs algorithm.

   :param name: Sequence identity.
   :type name: str
   :param sequence: The DNA sequence.
   :type sequence: str

   .. py:attribute:: instances
      :annotation: = []

      

   .. py:method:: __str__(self)

      Return str(self).


   .. py:method:: init_site(self)

      Create a motif instance and appends it to the 'instance' attribute.

      :raises SystemExit: If the desired motif width is equal to or longer than the sequence.
      :raises AssertionError: If the motif instance is longer than the desired width.


   .. py:method:: find_site_prob(self, pwm, background)

      Find the probability of the motif at each position.

      Passed as an argument to sample_new_motif method.

      :param pwm: A position weight matrix generated from the motif instances.
      :type pwm: Bio.motifs.matrix.PositionWeightMatrix
      :param background: Background nucleotide frequencies.
      :type background: dict of {str : float}

      :returns: * **motif_score** (*list of float*) -- Contains the motif score at each position in the sequence.
                * **total_score** (*float*) -- The sum of each motif score.


   .. py:method:: sample_new_motif(self, find_site_prob, pwm, background)

      Select a new motif instance.

      The selection is random but weighted by the probability of
      the motif occurring at each position.
      Probabilities are calculated using the find_site_prob() method.

      :param find_site_prob: Method for determining the motif site probability in a sequence.
      :type find_site_prob: Sequence class method
      :param pwm: Arguments passed into the find_site_prob() method.
      :type pwm: See find_site_prob()
      :param background: Arguments passed into the find_site_prob() method.
      :type background: See find_site_prob()

      :returns: **new_motif** -- A new motif instance.
      :rtype: str

      :raises AssertionError: If normalized motif scores do not sum to 1.



.. py:function:: get_sequences(file)

   Parse each sequence in a fasta file and create a Sequence object.

   :param file: File handle containing DNA sequences in fasta format.
   :type file: io.TextIOWrapper

   :returns: **sequences** -- A list containing all parsed sequences.
   :rtype: list of Sequence objects

   :raises SystemExit: If one or no sequences were found in the file.
   :raises TypeError: If input file of sequences is not in fasta format.


.. py:function:: background_freq(sequences)

   Find the total nucleotide frequencies from all sequences.

   :param sequences: A list containing the DNA sequences.
   :type sequences: list of Sequence objects

   :returns: **background** -- The background frequency for each nucleotide.
   :rtype: dict of {str : float}

   :raises AssertionError: If  frequencies do not sum to one.


.. py:function:: calc_pwm(motif_instances)

   Calculate position weight matrix from all motif instances.

   :param motif_instances: A list containing randomly selected motif instance from each sequence.
   :type motif_instances: list of str

   :returns: **pwm** -- A position weight matrix generated from motif instances.
   :rtype: Bio.motifs.matrix.PositionWeightMatrix


.. py:function:: run_sampler(sequences, motif_instances, background)

   The Gibbs Sampler algorithm.

   Iterates through the sampling, predictive update, and likelihood
   calculation steps. During each iteration except the last, the
   current pwm and motif instance are saved to test for convergence
   after the final iteration of the algorithm.

   :param sequences: List of all sequences  contained in fasta file.
   :type sequences: list of str
   :param motif_instances: List of all motif instances generated from the Sequence class.
   :type motif_instances: list of str
   :param background: Background nucleotide frequencies.
   :type background: dict of {str : float}

   :returns: * **new_motif** (*str*) -- The final motif generated from the sampler.
             * **old_motif** (*str*) -- The second to last motif generated.
             * **new_pwm** (*Bio.motifs.matrix.PositionWeightMatrix*) -- The final pwm generated from the sampler.
             * **old_pwm** (*Bio.motifs.matrix.PositionWeightMatrix*) -- The second to last pwm generated from the sampler.


.. py:function:: change_in_pwm(final_pwm, initial_pwm)

   Calculate the percent change between the ultimate and penultimate pwm.

   A test for convergence. Takes the final two matrices
   generated from the run_sampler function and compares the values
   from each column in the matrices.

   :param final_pwm: The final pwm generated from run_sampler().
   :type final_pwm: Bio.motifs.matrix.PositionWeightMatrix
   :param initial_pwm: The second to last pwm generated from run_sampler().
   :type initial_pwm: Bio.motifs.matrix.PositionWeightMatrix

   :returns: * **False** (*bool*) -- If the percent change between the final and initial pwm is
               greater than or equal to one percent.
             * **True** (*bool*) -- If the percent change is less than one percent.


.. py:function:: change_in_motif(final_motif, initial_motif)

   Compare the ultimate and penultimate motif instances from run_sampler.

   A test for convergence. Determines if any motif subsequence location
   has changed in the final two iterations of the sampler.

   :param final_motif: The final motif generated from run_sampler().
   :type final_motif: str
   :param initial_motif: The second to last motif instance generated.

   :returns: * **False** (*bool*) -- If the two motif instances do not match.
             * **True** (*bool*) -- If the motif instances do match.


.. py:function:: web_logo(motif_instances)

   Create a weblogo from the motif instances generated by run_sampler.

   Runs the  program 'WebLogo' by executing a terminal command
   using subprocess.
   The input is the list of motif instances which are writen into stdin
   and the result is written to an output file.

   :param motif_instances: A list containing the motif instances generated from run_sampler()
   :type motif_instances: list of str

   :raises AssertionError: If process has not finished.


.. py:function:: main()


