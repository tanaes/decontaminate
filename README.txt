README -- decontaminate.py 
version -- 0.3 (mothur compatible)
author -- Jon Sanders <jonsan@gmail.com>

This is a simple script for the removal of putative contaminant sequences.

It depends on an installed and working version of QIIME 1.9.x for functionality.

There are two options for running this script: the freestanding python script 
(decontaminate.py in the top-level directory) or a qiime-integrated version. 

For the former, no installation is necessary. Simply run the python script 
provided. 

For the latter, the scripts in the 'qiime_scripts' directory must either be copied to 
the appropriate directories in your base qiime installation, or the folder decontaminate/qiime_scripts/qiime must be added to your $PYTHONPATH variable. I haven't fully worked this out yet.

Additionally, a wrapper bash script (filter_contaminants.sh) has been provided 
as an example of a complete decontamination procedure. As written, this will 
take a post-split_libraries.py format sequence file, demultiplex into unique 
sequences, cluster just the blanks to form a potential contaminant library, 
identify putative contaminants using the decontaminate.py script and both
reference and abundance-based methods, and filter the starting fasta file down
to just the non-contaminant sequences. This file may be modified to fit your 
needs.
