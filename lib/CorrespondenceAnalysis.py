import Numeric, LinearAlgebra

##---------------------------------------------------------------------
def GetIndices( matrix ):
    """return order (1st eigenvector) of row and column indicies
    """

    nrows, ncols = matrix.shape

    # calculate row and column sums
    row_sums = Numeric.sum( matrix, 1 )
    col_sums = Numeric.sum( matrix, 0 )    

    a = Numeric.zeros( (nrows, nrows), Numeric.Float)
    for x in range( 0, nrows):
        a[x,x] = 1.0 / float(row_sums[x])

    b = Numeric.zeros( (ncols, ncols), Numeric.Float)    
    for x in range( 0, ncols):
        b[x,x] = 1.0 / float(col_sums[x])
        
    M = Numeric.matrixmultiply( \
        a, Numeric.matrixmultiply( \
        matrix, Numeric.matrixmultiply( \
        b, Numeric.transpose( matrix ))))

    try:
        row_eigenvector = LinearAlgebra.eigenvectors(M)[1][1]
    except:
        return None, None
    
    M = Numeric.matrixmultiply( \
        b, Numeric.matrixmultiply( \
        Numeric.transpose(matrix), Numeric.matrixmultiply( \
        a, matrix )))

    try:
        col_eigenvector = LinearAlgebra.eigenvectors(M)[1][1]    
    except:
        return None, None
    
    return row_eigenvector.astype(Numeric.Float) , col_eigenvector.astype(Numeric.Float)

##---------------------------------------------------------------------
def GetPermutatedMatrix( matrix, map_row_new2old, map_col_new2old,
                         row_headers = None, col_headers = None):
    """return a permuted matrix. Note, that currently this is very
    inefficient, as I do not know how to do this in Numeric.
    """

    nrows, ncols = matrix.shape

    result = Numeric.zeros( (nrows, ncols), matrix.typecode())
    for r in range(0, nrows):
        for c in range(0,ncols):
            result[r,c] = matrix[map_row_new2old[r], map_col_new2old[c]]

    if not row_headers or not col_headers:
        return result

    rows = []
    for x in map_row_new2old:
        rows.append( row_headers[x] )
        
    cols = []
    for x in map_col_new2old:
        cols.append( col_headers[x] )

    return result, rows, cols
        
##---------------------------------------------------------------------
def PermuteRows( matrix ):
    pass
    

if __name__ == "__main__":

    num_rows = 6
    num_cols = 5
    matrix = Numeric.zeros( (num_rows,num_cols), Numeric.Int)

    matrix[0,2] = 1
    matrix[0,3] = 1
    matrix[1,0] = 1
    matrix[1,1] = 1
    matrix[1,4] = 1
    matrix[2,1:5] = 1
    matrix[3,2:4] = 1
    matrix[4,0] = 1
    matrix[4,4] = 1
    matrix[5,0] = 1
    matrix[5,2:5] = 1

    print "matrix=", matrix

    row_indices, col_indices =  GetIndices( matrix )

    map_row_new2old = Numeric.argsort(row_indices)
    map_col_new2old = Numeric.argsort(col_indices)

    print map_row_new2old
    print map_col_new2old

    print GetPermutatedMatrix( matrix, map_row_new2old, map_col_new2old)
    

