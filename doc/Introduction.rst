========
Overview
========

Synopsis
========

ADDA is a method to find protein domains in protein sequences.
Briefly, ADDA attempts to split domains into segments that 
correspond as closely as possible to all-on-all pairwise
alignments. A detailed description of the method can be found
in

   Heger A, Holm L. (2003) Exhaustive enumeration of protein domain families.
   J Mol Biol. 2003 May 2;328(3):749-67.
   PMID: `12706730 <http://www.ncbi.nlm.nih.gov/pubmed/12706730>`_ 

Installation
============

Getting ADDA
-------------

Adda is available from http://sourceforge.net/projects/adda. The code
is available from the subversion repository::

   svn co https://adda.svn.sourceforge.net/svnroot/adda/trunk adda

Requirements
------------

ADDA requires the following software:

  * python (2.6) http://www.python.org
  * alignlib (install from subversion) http://sourceforge.net/projects/alignlib
  * matplotlib (0.98 or higher) http://matplotlib.sourceforge.net
  * numpy (1.3.0 or higher) http://numpy.scipy.org
  * cython (0.11.1 or higher) http://www.cython.org

Notes:

  * Python 3 will not work as boost.python support is not there yet
  * Python 2.6 is required for the :mod:`multiprocessing` module. If python 2.5
      is used, install it separately.

Compilation
-----------

Compilation works in the usual pythonic way::

   python setup.py build
   python setup.py install


Citing ADDA
-----------

If you ADDA it in your work, please cite

   Heger A, Holm L. (2003) Exhaustive enumeration of protein domain families.
   J Mol Biol. 2003 May 2;328(3):749-67.
   PMID: `12706730 <http://www.ncbi.nlm.nih.gov/pubmed/12706730>`_ 
