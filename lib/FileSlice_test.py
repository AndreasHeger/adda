import unittest, os, glob, re, tempfile

import FileSlice

class TestFileSlice(unittest.TestCase):

    mNumLines = 2
    mNumRecords = 20
    def setUp(self):
        fd, self.mFilename = tempfile.mkstemp()
        os.close(fd)
        outfile = open(self.mFilename, "w")
        l = 0
        for x in range(0,self.mNumRecords): 
            for y in range(0,self.mNumLines): 
                outfile.write( "%i\t%i\t%i\t%i\n" % (x,y,l,y) )
                l += 1
        outfile.close()
        
    def tearDown(self):
        if os.path.exists(self.mFilename):
            os.remove( self.mFilename )

    def checkComplete(self, nchunks ):

        all_data = []
        for chunk in range( 0, nchunks):
            iterator = FileSlice.Iterator( self.mFilename, 
                                           nchunks,
                                           chunk, 
                                           FileSlice.iterator )
            
            data = []
            for d in iterator:
                d = d.split("\t")
                self.assertEqual( len(d), 4 )
                data.append(d)
            all_data += data

        self.assertEqual( len(all_data), self.mNumRecords * self.mNumLines )
        self.assertEqual( [int(x[2]) for x in all_data], range( 0, self.mNumRecords * self.mNumLines ) )

    def testComplete(self):
        self.checkComplete( self.mNumRecords )
        
    def testHalf(self):
        self.checkComplete( self.mNumRecords // 3 )
        
    def testThird(self):
        self.checkComplete( self.mNumRecords // 2 )

    def testSingle(self):
        self.checkComplete( 1 )

    def testDouble(self):
        self.checkComplete( 2 * self.mNumRecords )

    def testPlus1(self):
        self.checkComplete( self.mNumRecords + 1)

if __name__ == '__main__':
    unittest.main()
