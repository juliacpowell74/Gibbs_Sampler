## Gibbs_Sampler

This program runs the Gibbs Sampler algorithm for
de novo motif discovery. Given a set of sequences, the program
will calculate the most likely motif instance as well as 
the position weight matrix and position specific scoring matrix
(the log2 normalized frequency scores).

The input sequence file should be provided in fasta format. Output is written to the console. 
A eps file containing the weblogo of the motif is also created.
#Getting Started
##Prerequisites
Python 3.8 or later is required.
##Installation
The program can be installed with pip by running the following command:
```bash
$ pip install gibbs_sampler
```

## Usage

To run, simply call gibbs_sampler from the command line along with 
the path to the sequence file and expected width of the motif.

````bash
$ gibbs_sampler <primary sequence file> <width>
````
Optional arguments for the title of the weblogo file and number of 
iterations of the sampler can be specified using the -t and -n flags respectively.

## License

`gibbs_sampler` was created by Monty Python. It is licensed under the terms of the CC0 v1.0 Universal license.

## Credits
Thank you to Professor Hendrix for teaching me how to use python to investigate biologically relevant questions.
