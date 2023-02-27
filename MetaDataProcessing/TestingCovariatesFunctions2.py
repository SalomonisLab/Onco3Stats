#
#
def importTestCovFile2( fn, delimChar='\t', *args, metaDF=None ):
    """Import TestCov(ariate) file for user-specified statistical testing"""
    """Input arguments are path/filename and the delimiter used in the file"""
    """Optionally, can provide a metadata dataframe to perform additional checks on the TestCov file contents"""

    """TestCov file must: have 3 columns and NOT have a header row"""
    """General form of a row is metavariable name, role/use of metavariable, and (comma-separated) corresponding values"""
    
    
    """ returns a pandas DataFrame """  
    
    import pandas as pd
    
    allowedMetaVarRoles = ['UID','Restrict','Exclude','NumFilter','Unique','Covariate']
    
    try:
        tc = pd.read_table( fn, sep=delimChar, header=None, comment='#')
    except:
        print ("A problem occurred importing the Test-Covariate file:",fn)
        return None

#   File must have 3 columns
    if tc.shape[1] != 3:
        print( 'The Test-Covariate file',fn,' must have 3 columns')
        return None
    
    tc.columns = ['mVar', 'mRole', 'mVals']
    
#   Check to see specified meta variable roles (e.g., Restrict, Covariate) are allowed
    roles = list( tc['mRole'][tc['mRole'].notna()] ) 
    if not all(r in allowedMetaVarRoles  for r in roles) :
        print( 'Second column in Test-Covariate file has a value not among ',allowedMetaVarRoles )
        return None
    
#   Check that any numeric filtering NumF entries have proper syntax
    NumF_df = tc[tc['mRole']=='NumFilter']
    if NumF_df.shape[0] > 0:
        for k in NumF_df.index:
            kVar = NumF_df['mVar'][k]
            kVals =  NumF_df['mVals'][k]
            cutoff = kVals[1:]
            direction = kVals[0]
            try: 
                cutnum = pd.to_numeric(cutoff)
            except:
                print ("Cutoff value in NumF statement could not be converted to numeric:",cutoff)
                return None
            if not any(r==direction  for r in ['<','>']):
                print ("NumF statement must contain > or < as first character in 3rd column")
                return None
    
#   
# --- If metaDF is provided, check that all meta vars in Test-Covariate file are in UserMetaData file (includes designated UID variable)
    if metaDF is not None:
        if not all( v in list( metaDF.columns )  for v in list( tc['mVar'] ) ):
            print( 'Variable in Test-Covariate file not in MetaData file' )
            return None
    
    return tc

#
# ---------------------------------------------------------------------
#
#
def ApplyTestCovToMeta2(tc=None,mt=None):
#    This function applies any filtering specified in tc (pandas dataframe)
#     ['Restrict', 'Exclude', 'NumFilter' or 'Unique'] to mt (pandas dataframe).
#    The new metadata pandas dataframe (filtered as indicated in tc) is returned

    import pandas as pd
    import numpy as np
    
#   Gather TestCovariate (tc) rows with filtering roles
    Restrict_df = tc[tc['mRole']=='Restrict']
    Exclude_df = tc[tc['mRole']=='Exclude']
    NumF_df = tc[tc['mRole']=='NumFilter']
    Unique_df = tc[tc['mRole']=='Unique']
#    X_df = tc[tc['mRole']=='X']


    fmeta = mt.copy()   
#   Create flag for any meta row to remove
    fmeta['Remov'] = 0
  
    
#   Handle any RESTRICT filtering
    if Restrict_df.shape[0] > 0:
        
        for k in Restrict_df.index:
            kVar = Restrict_df['mVar'][k]
            kVals =  Restrict_df['mVals'][k] 

            kVals = list( kVals.split(',') ) 
            kVals = list(map(str.strip, kVals))

            v = fmeta[kVar]
            remv = np.logical_not( v.isin(kVals) )

            fmeta.loc[remv , 'Remov'] = 1


#   Handle any EXCLUDE filtering
    if Exclude_df.shape[0] > 0:
        
        for k in Exclude_df.index:
            kVar = Exclude_df['mVar'][k]
            kVals =  Exclude_df['mVals'][k] 

            kVals = list( kVals.split(',') ) 
            kVals = list(map(str.strip, kVals))

            v = fmeta[kVar]
            remv =  v.isin(kVals)

            fmeta.loc[remv , 'Remov'] = 2


#   Handle any NumF (numeric) filtering
#     Note that the import function has already checked that direction is either > or <
#      and the (single) value can be coerced to numeric.
#     Example of numeric filtering on an age variable:
#       >19 means keep only samples w/ age > 19
#        coding approach is to remove any records for which age > 19 can not be established

    if NumF_df.shape[0] > 0:
        
        for k in NumF_df.index:

            kVar = NumF_df['mVar'][k]
            kVals =  NumF_df['mVals'][k] 

            cutoff = kVals[1:]
            direction = kVals[0]

#     the following will set any missing or character values to pd.nan            
            numv = pd.to_numeric(fmeta[kVar], errors='coerce')
            cutnum = pd.to_numeric(cutoff)
            
            if direction=='>': 
                kp = numv > cutnum
            else:
                kp = numv < cutnum

            fmeta.loc[np.logical_not(kp) , 'Remov'] = 3


#   Before processing any Unique statements, remove records flagged thus far.
    fmeta = fmeta[fmeta['Remov']==0] 

         
#   Handle any Unique (duplicate removal) filtering
    if Unique_df.shape[0] > 0:
        for k in Unique_df.index:
            kVar = Unique_df['mVar'][k]
            fmeta.drop_duplicates(subset=[kVar], keep='first',inplace=True)
        
    fmeta.drop(columns=['Remov'],inplace=True)
    
    return fmeta

#
# ---------------------------------------------------------------------
#
#
def makeTestGroups2(Cov_df,meta_df):
#
#  Function creates groups for each test specified in Cov_df, using meta_df
#   Cov_df is pandas dataframe of "Covariate" rows that specifies 2 groups based on a metadata variable
#   meta_df is a pandas dataframe (possibly filtered according to other meta variables) 
#     used to ID samples with particular values of covariates.
#     sample identifiers must be the (row) index of meta_df
#
#   function returns a dictionary of dictionaries of the form
#     { comparisonA_name: { 'g1': samples_in_group_1,'g2': samples_in_group_2 },
#       comparisonB_name: { 'g1': samples_in_group_1,'g2': samples_in_group_2 }}
#
    import pandas as pd
    import numpy as np
    
    tc_dict = { }
    
    for k in Cov_df.index:
        kVar = Cov_df['mVar'][k]
        kVals = Cov_df['mVals'][k]

        kVals = list( kVals.split(',') ) 
        kVals = list(map(str.strip, kVals))

        if len( kVals ) != 2:
            print('Covariate statement must have 2 comma-separated values in 3rd column to define groups.')
            return tc_dict
        
        g1_name = kVals[0]
        g2_name = kVals[1]
        comp_name = kVar + '__' + g1_name + '_vs_' + g2_name
    
        g1a = list( meta_df.index[meta_df[kVar]==kVals[0]] )
        g2a = list( meta_df.index[meta_df[kVar]==kVals[1]] )
    
        tc_dict[comp_name] = { 'g1':g1a, 'g2':g2a } 

    return tc_dict

    
