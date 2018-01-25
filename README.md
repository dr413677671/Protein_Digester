The Protein_digester script makes it easy to digest proteins.

The program is able to cope with one or more protein sequences, and digest them into small peptide sequences. The script subsequently also automatically
generates the outcome data in fasta format.

As of Python >= 3.5.3, the script require some Python standard library. For users who still need to support Python 3(Anaconda),several packages are required: re(Regular expression), random and argparse, but also may support older Python versions.



Compatibility
-------------

Protein_digester should work on Python >= 3 (it was tested on 3.5.3).


Configration
------------

The program should be config in Linux command line.
The regular command seems like: ~$ python Protein_digester.py 'FILE_NAME' -o 'OUTPUT_FILENAME' -e 1_or_0 -c CLEAVAGE_NUM -t ENZYME_TYPE_NUM
positional arguments:
  files                 input sequence files(fasta format) with and devided by
                        blank
 
optional arguments:
  -h, --help            show this help message and exit
  -t TYPE, --ENZYME_TYPE TYPE
                        option of the digesting enzyme type: 1 for Trypsin:
                        cuts at Lysine (Lys, K) or Arginine (Arg, R) unless
                        the next amino acid is Proline (Pro, P), 2 for
                        Endoproteinase Lys-C: cuts at Lysine (Lys, K) unless
                        the next amino acid is Proline (Pro, P), 3 for
                        Endoproteinase Arg-C: cuts at Arginine (Arg, R) unless
                        the next amino acid is Proline (Pro, P), 4 for V8
                        proteinase (Glu-C): cuts at Glutamic acid (Glu, E)
                        unless the next amino acid is Proline (Pro, P) and
                        default for enzyme Trypsin
  -e ERRORTRAPPING, --ERROR_TRAPPING ERRORTRAPPING
                        option weather to report the unusual amino acid or the
                        unknown one(unusual amino acid as B,O,U,J,Z; unknown
                        amino acid as X) beside peptide infor mation in the
                        fasta format output: 1 for report 0 for not report
                        default 0
  -o OUTPUT, --OUTPUT OUTPUT
                        set the output filename
  -c CLEAVAGE, --CLEAVAGE CLEAVAGE
                        input the number of miss cleavage of the analysis
                
Bugs
----

If you find a bug, please try to reproduce it with python 3.5.3.

If it happens there also, please file a bug in the python.org issue tracker.
If it does not happen in 3.5.3, file a bug to 413677671@qq.com.


Ran Duan
Email: 413677671@qq.com