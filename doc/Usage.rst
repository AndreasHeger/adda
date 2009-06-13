=======
Usage
=======

Overview
========

ADDA is controlled with the script adda.py. The script expects a file adda.ini 
with configuration options in the directory from which it is called. An example 
is in the directory ./test.

Input data
==========

ADDA requires three input files

   * a file with sequences (`SequenceFile`_)

   * a file with a pairwise similarity matrix resulting from all-on-all
     BLAST searches (`PairsdbFile`_)

   * a file with domain assignments from a reference database (`ReferenceFile`_)

   * a configuration file (`AddaIni`_)

A toy example with these files can be found at http://genserv.anat.ox.ac.uk/downloads/contrib/adda.

.. _SequenceFile:

Sequence file
-------------

nrdb.fasta.gz

The sequence file is in fasta format. The format is::

   >Identifier Description
   SEQUENCE

The Identifier is delimited by the first white-space character. The
description is optional and is ignored. The sequence can span multiple
lines.
     
Only sequences appearing in this file will be used by ADDA. These should have been 
filtered to be less than 40% identical (Park et al. 2000).

.. _PairsdbFile:

Pairsdb file
------------

pairsdb.links.gz

A list of pairwise alignments. These have been obtained by
running BLASTP (Altschul et al. 1997) all-on-all and parsed into a tab-separated 
table. The columns are:

query
   identifier of the query sequence
sbjct
   identifier of the sbjct sequence
evalue
   natural log of the E-Value
query_start
   first aligned residue in query
query_end
   last aligned residue + 1 in query
query_ali
   compressed alignment of the query
sbjct_start
   first aligned residue in sbjct
sbjct_end
   last aligned residue + 1 in sbjct
sbjct_ali
   compressed alignment of the sbjct
alignment score 
   (not used)
percent identity 
   (not used)

The alignment coordinates are inclusive/exclusive 0-based coordinates "[)".
The alignment is stored in compressed form as alternating integer numbers
with the prefix "+" and "-". Positive numbers signify character emissions 
and negative numbers insertions. For example, "+3-3+2" with the sequence 
ABCDE will result in "ABC---DE".

ADDA can read alternative formats. See the command line option ``--alignment-format``.

This file can be read in compressed form (gzip, suffix ".gz").

.. _ReferenceFile:

Reference file
==============

references.domains.gz

A list of domains for a subset of sequences. The domains in this file were derived 
from structural domain definitions in SCOP (Andreeva et al. 2008).

This file can be read in compressed form (gzip, suffix ".gz").

.. _AddaIni:

adda.ini
--------

The :file:`adda.ini` should be present in the working directory. It contains
various configuration options for the pipeline and is grouped into sections.

Of interest should be the section ``[files]`` which lists the input and 
output filenames. In particular,

input_graph
   the file with pairwise sequence similarity options (`PairsdbFile`_)

input_fasta
   the file with sequence information (`SequenceFile`_)

input_reference
   the file with reference domain information (`ReferenceFile`_)

Running ADDA
============

ADDA is controlled through a single script :file:`adda.py`.
ADDA proceeds in stage and each stage corresponds to a command to
the script adda.py. To run all stages, run adda.py as

   python adda.py --steps=all

The specific stages are:

   0 ``sequences``: read the sequence file and relabel numerically (:class:`AddaSequences.AddaSequences`)

   1 Pre-processing of graph. These steps are performed in parallel

      1 ``segment`` - compute hierachical decomposition of sequences (:class:`AddaSegment.AddaSegment`)
      
      2 ``profiles`` - build sequence profiles (:class:`AddaProfiles.AddaProfiles`)
      
      3 ``graph`` - format and filtering the alignment graph (:class:`AddaGraph.AddaGraph`)
      
      4 ``fit`` - estimate the error parameters (:class:`AddaFit.AddaFit`)
 
   2 ``index`` index the alignment graph (:class:`AddaIndex.AddaIndex`)

   3 ``optimise`` - decompose sequences into domains (:class:`AddaOptimise.AddaOptimise`)

   4 ``convert`` - convert the alignment graph into a domain graph (:class:`AddaConvert.AddaConvert`)

   5  ``mst`` - build minimum spanning tree of domains (:class:`AddaMst.AddaMst`)

   6 ``align`` - align domains (:class:`AddaAlign.AddaAlign`)

   7 ``cluster`` - cluster domains into families (:class:`AddaCluster.AddaCluster`)

   8 ``families`` - process singletons (:class:`AddaFamilies.AddaFamilies`)

Disk usage
----------

ADDA requires a large amount of temporary diskspace. The largest components are

   * the :file:`files:profile library`. The file size depends on the number of sequences. In the
      test set, 220,000 sequences require 2.6 Gb space.

   * the :file:`files:domain_graph`. The file size depends on the number of domains and 
      the number of pairwise alignments. In the test set, this file is 1.2 Gb.

   * the :file:`files:graph`. The file size depends on the number of pairwise alignments.
      In the test set, this file is 605 Mb.

Parallel runs
-------------

ADDA can use several CPU/cores if available for steps that are embarrassingly parallel.
These steps will create several output files with numeric suffixes that will be later
merged into a single file.

Aborted runs
------------

ADDA will pick up from aborted runs and continue without re-computing previously
computed steps. It will check if a step has run to completion by examining the file
contents and not time stamps. In particular, it will check if a file ends with the
token ''#\\''.

Output
======

TODO

TODO
====

  * Be more economical with disc space. Investigate the use of compressed files.






References
==========

Altschul SF, Madden TL, Schäffer AA, Zhang J, Zhang Z, Miller W, Lipman DJ.
(1997) Gapped BLAST and PSI-BLAST: a new generation of protein database search programs.
Nucleic Acids Res. Sep 1;25(17):3389-402.

Park J, Holm L, Heger A, Chothia C. (2000) RSDB: representative protein 
sequence databases have high information content. Bioinformatics. May;16(5):458-64.

Andreeva A, Howorth D, Chandonia JM, Brenner SE, Hubbard TJ, Chothia C, Murzin AG.
(2008) Data growth and its impact on the SCOP database: new developments.
Nucleic Acids Res. Jan;36(Database issue):D419-25. Epub 2007 Nov 13.
