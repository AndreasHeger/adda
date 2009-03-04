import sys, os, re, time

import alignlib
import ProfileLibrary

from AddaModule import AddaModule
import AddaIO
import IndexedFasta

class FastaRecord:
    def __init__(self, title, sequence ):
        self.title = title
        self.pid = re.sub( "\s.*", "", title )
        self.sequence = sequence

def _iterate( infile ):

    h = infile.readline()[:-1]

    if not h: raise StopIteration

    while h[0] != ">":
        h = infile.readline()[:-1]
        if not h: raise StopIteration
        continue

    h = h[1:]
    seq = []

    for line in infile:
        if line[0] == "#": continue
        line = line[:-1] # remove newline
        if not line: continue
        if line[0] == '>':
            yield FastaRecord(h,''.join(seq))
            h = line[1:]
            seq = []
            continue

        seq.append(line)
    yield FastaRecord(h,''.join(seq))

class FastaIterator:

    def __init__(self, f, *args, **kwargs):
        self.mIterator = _iterate(f)
    def __iter__(self):
        return self
    def next(self):
        return self.mIterator.next()

class AddaSequences( AddaModule ):
    """output sequence information."""
    
    mName = "Sequences"
    
    def __init__(self, *args, **kwargs ):

        AddaModule.__init__( self, *args, **kwargs )
                
        self.mFilenameNids = self.mConfig.get( "files", "output_nids", "adda.nids" )    
        self.mFilenameInputFasta = self.mConfig.get( "files", "input_fasta", "adda" )
        self.mFilenameOutputFasta = self.mConfig.get( "files", "output_fasta", "adda" )
        self.mMaxSequenceLength = self.mConfig.get( "segments", "max_sequence_length", 10000 )

    def applyMethod(self ):

        self.mInput = 0
        self.mOutput = 0
        self.mRemoved = 0

        outfile = self.openOutputStream(self.mFilenameNids)
        outfile.write( "nid\thid\tpid\tlength\tsequence\n" )

        nid = 1
        hids = set()
        
        fasta = IndexedFasta.IndexedFasta( self.mFilenameOutputFasta, "w" )

        for seq in FastaIterator( open( self.mFilenameInputFasta, "r" ) ):
            
            self.mInput += 1
            if len( seq.sequence ) > self.mMaxSequenceLength:
                self.mRemoved += 1
                continue

            hid = self.getHID( seq.sequence )
            if hid in hids:
                self.mDuplicates += 1
                continue
            
            hids.add(hid)
            outfile.write( "%s\t%s\t%s\t%i\t%s\n" % (nid, hid, seq.pid, len(seq.sequence), seq.sequence) )
            nid += 1
            fasta.addSequence( nid, seq.sequence )
            self.mOutput += 1

        outfile.close()

    def finish(self):
        
        self.info( "sequences: %i input, %i output, %i removed" %\
                   (self.mInput, self.mOutput, self.mRemoved ) )
        
        AddaModule.finish( self )
        
