====================
Mapping ADDA to PFAM
====================

The section describes how to map |ADDA| domains onto the latest |PFAM|
dataset. Briefly, the sequences in the latest |ADDA| and |PFAM| releases
are compared. Sequences present in |PFAM| but not in |ADDA| are
mapped using |BLAT| onto sequences within |ADDA|.

Preparation
===========

For mapping, a file called :file:`adda.ini` needs
to be present in the working directory. In this file,
the following variables need to be set in the section
``[map]``:

.. glossary::

   filename_adda_sequences
       filename with adda sequences, for example :file:`adda.fasta.gz`

   filename_domains
       filename with adda domains, for example :file:`adda.domains`

   filename_target_sequences
       filename with the target (PFAM) sequences, for example 
       :file:`pfam.fasta.gz`

   chunksize
       size of BLAT jobs. :term:`chunksize` queries are collected and
       run together in one batch.

   min_identity
       minimum percent identity. Only alignments of :term:`min_identity`
       will be used for mapping domains.

The mapping procedure requires the following executables to be in
the :envvar:`PATH`:

* blat: the BLAT search program
* gunzip: gzip utility

Running
=======

The mapping pipeline is run via the command::

    adda.py make map




