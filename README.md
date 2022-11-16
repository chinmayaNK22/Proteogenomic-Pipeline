# Proteogenomic-Pipeline

In order to execute the proteogenomic pipeline, it is crucial to generates a sixframe translated protein sequence for a given genome sequence file along with their genome coordinates.

This Proteogenomic-Pipeline repository conists of python packages for generating sixframe translated protein sequences and a gtf file generating package for Genome Search-Specific Peptides (GSSPs).


# How to use six_frame_translation.py
```
>python six_frame_translation.py test_genome.fasta

usage: six_frame_translation.py [-h] -i [-i ...]

Generates a sixframe translated protein sequence database (FASTA) for a given
genome sequence database

positional arguments:
  -i          Genome sequence (FASTA format) to be sixframe translated into
              protein sequence

optional arguments:
  -h, --help  show this help message and exit
```

# How to use generate_gtf.py
```
>python generate_gtf.py PeptideGroups.txt Proteome.fasta sixframe_proteome.fasta

usage: generate_gtf.py [-h] -i [-i ...] -f [-f ...] -sf [-sf ...]

Generates GTF file for peptides identified from Sixframe translated protein
sequence database search in Proteome Discoverer

positional arguments:
  -i          PeptideGroups output from Proteome Discoverer
  -f          Proteome database used in the first step (FASTA format)
  -sf         Sixframe translated proteome database used (FASTA format)

optional arguments:
  -h, --help  show this help message and exit
```

# Citation
Please cite this tool using the DOI [![DOI](https://zenodo.org/badge/388470655.svg)](https://zenodo.org/badge/latestdoi/388470655) 
