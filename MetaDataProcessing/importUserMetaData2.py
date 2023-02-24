#
#
def importUserMetaData2(fn,delimChar='\t', *args, SamplesList=None):

    """Import metadata that can be used for among sample/cell group comparisons"""
    """Input arguments are path/filename for the metadata and the delimiter used in the file"""
    """Optionally, can provide a list of sample names to check the meta file sample against"""

    """Top row of input file must be meta variable names and first column must be the sample identifier variable"""
    
    """ return as pandas DataFrame.  Sample identifiers in column 1 are NOT moved to row index position """  
    
    import pandas as pd
    
    try:
        mt = pd.read_table(fn,sep=delimChar,header=0,comment='#')
    except:
        print ("A problem occurred importing the metadata file:",fn)
        return None

#   File should have >1 row and > 1 column
    if mt.shape[0]==1:
        print( 'Only 1 row read from metadata file',fn,' please check file')
        return None
    elif mt.shape[1]==1:
        print( 'Only 1 column read from metadata file',fn,' please check file and delimiter')
        return None
        
#   If SamplesList was provided, check that contents of column 1 is 1-1 match to contents of Samples

    if SamplesList is not None:
#        print( 'Sample names provided')
        if SamplesList.__class__ != list:
            print( 'Sample names object not a list, sample names check not done')
            return mt

        UserMetaDataSamps = list( mt.iloc[:, 0] )
        UserMetaDataSamps.sort()
        SamplesList.sort()
        if UserMetaDataSamps != SamplesList:
            print('Identifiers in column 1 of metadata file',fn,'are not 1-1 match to identifiers in SamplesList or data matrix')
        
    return mt