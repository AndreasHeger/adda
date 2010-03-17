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
       filename with adda sequences, for example :file:`adda.fasta.gz`. This
       is a |fasta| formatted file. Sequences are labeled by |nid|.

   filename_domains
       filename with adda domains, for example :file:`adda.results`. This is
       a tab-separated file with four columns: sequence identifier (|nid|),
       start, end and family. Sequence coordinates are 0-based half-open. 
       The first line is a header. For example::

	   nid     start   end     family
	   27      0       101     AD092610
	   30      0       147     AD005437
	   42      0       115     AD211211
	   43      0       474     AD004688
	   44      0       150     AD123296
	   55      0       320     AD074056
	   55      320     440     AD054547

   filename_target_sequences
       filename with the target (PFAM) sequences, for example 
       :file:`pfam.fasta.gz`. The first characters up to the first
       whitespace will be used as the sequence identifiers.

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

Results
=======

The map command creates the following files:

.. glossary::

   all.domains
      A file with |ADDA| domains mapped to the sequences in 
      :term:`filename_target_sequences`. The format is identical
      to :term:`filename_domains`.

   direct.domains
      A file with |ADDA| domains in sequences that are part of
      :term:`filename_target_sequences`. The format is identical
      to :term:`filename_domains`.

   mapped.domains
      A file with |ADDA| domains that have been mapped via BLAT
      onto :term:`filename_target_sequences`. The format is identical
      to :term:`filename_domains`.

   indirect.domains
      A file with |ADDA| domains that have been mapped via BLAT
      onto redundant :term:`filename_target_sequences`. The format 
      is identical to :term:`filename_domains`.
   
   mapping.summary
      Summary statistics of the mapping process.

   mapping.coverage
      Table delineating the sequence coverage of each sequence in
      :term:`filename_target_sequences`.



