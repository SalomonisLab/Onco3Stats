#
#
#  kw3.py
#
#  January 1, 2023
#
#  python3 programs to support kruskal-wallis test
#  test is basically a 2-group t-test on ranked data



def kruskalWallisChecks( dataf, grpdf, *args, min_group_size=6 ): 
#
#   This function performs some checks before calling KruskalWallisCompute().    
#
#   dataf is a pandas data frame with samples as columns and events (features) as rows.  
#   Event names must be the row names (index), not in a column.
#   Missing values (np.nan='NaN') are allowed
#
#   grpsdf is a pandas data frame with a group membership column ('grp' for now ).
#   Sample names must be the row names (index), not in a column.
#   Aside from ordering, must be 1-1 match between dataf.columns and grpdf.index
#
#   min_group_size is minimum group size for kw test to be computed.
#   min_group_size = 6 is based on statistical rule-of-thumb for Mann-Whitney 2-group test.
#
#  Returns integer to indicate course of action:
#     0 = Problem, do not run KWTest
#     1 = Run KW test
#     2 = Two-sample situation, run MWU test instead

    import pandas as pd
    import numpy as np
    
    # --- Check for 1-1 match between sample names provided in grpdf and the data matrix

    if set(dataf.columns)!= set(grpdf.index) :
        print('Group information is inconsistent with names in data matrix')
        return( 0 )
  
        
    grpCounts = grpdf['grp'].value_counts().to_frame()  # returns a pandas df of counts, indexed by grp, decreasing by count
    grpCounts.columns=['grpRawN']
    nGrps = grpCounts.shape[0]
    minGrpN = min( grpCounts['grpRawN'] )
    # print('nGrps, minGrpN',nGrps,minGrpN)
    
    
    # -- Handle 2 or fewer groups --

    if nGrps < 2 :
        print('Number of sample groups is < 2; Kruskal-Wallis test not conducted')
        return( 0 )
        
    if nGrps == 2 :
        print('Only 2 sample groups found, switching to Mann-Whitney-U test')
        return( 2 )
    
    # -- Don't proceed if already know a group size is < minimum --   

    if minGrpN < min_group_size:
        print('Kruskal-Wallis test not conducted: 1 or more groups has fewer samples than minimum group size: ',minGrpN,' < ',min_group_size)
        return( 0 )
    
    return( 1 )


def KruskalWallisCompute(dataf, grpdf, *args, min_group_size=6 ): 
    #
    #   This function performs the Kruskal-Wallis NP ANOVA test.    
    #
    #   dataf is a pandas data frame with samples as columns and events (features) as rows.  
    #   Event names must be the row names (index), not in a column.
    #   Missing values (np.nan='NaN') are allowed
    #
    #   grpsdf is a pandas data frame with a group membership column.
    #   Sample names must be the row names (index), not in a column.
    #   Aside from ordering, must be 1-1 match between dataf.columns and grpdf.index
    #
    #   min_group_size is minimum group size for kw test to be computed.
    #   6 is based on statistical rule-of-thumb for Mann-Whitney 2-group test.
    #
    #   Returns pandas dataframe with columns containing group Ns (xNaN) & medians, value of KW test statistic, and p-value.
    #   Returned df has row for each event, even if test was not computed.
    
    import pandas as pd
    import numpy as np

    import scipy.stats
    

    # Set up little df with grp info, for group Ns and medians work
    
    grpCounts = grpdf['grp'].value_counts().to_frame()  # returns a pandas df of counts, indexed by grp, decreasing by count
    grpCounts.columns=['grpRawN']
    grpCounts['grpID']=grpCounts.index
    grpCounts = grpCounts.sort_index()
    # print(grpCounts) 
    nGrps = grpCounts.shape[0]

    # Compute group N and median for each event, in blocks by group, right-side joining as go. 
    #. Accumulate in resdf
    gindex = 0
    for gID in grpCounts.index:
        gindex = gindex + 1
        gSamps = grpdf.loc[grpdf['grp'] == gID, 'samp'].tolist()
    #    print(gSamps)
        gdf=dataf[gSamps]
    #    display(gdf)
    #
        meds = np.nanmedian(gdf,axis=1)
    #    print(meds)

        okVals = np.sum(~np.isnan(gdf), axis=1)
    #    print(okVals)

        tempdf=pd.DataFrame(zip(okVals,meds),index=okVals.index)
        tempdf.columns=[('N_' + str(gID)),('Median_' + str(gID))]
        if gindex == 1:
            resdf = tempdf.copy()
        else:
            resdf = resdf.set_index(resdf.index).join(tempdf)
 #       display(resdf)
# display(resdf) 
#    return( resdf )

#    Add columns for the kw statistic & p-value 
    
    resdf['kwStat'], resdf['kwPval'] = [ np.nan, np.nan ]
    

#   ID and extract data for events that meet the min_group_size criterion

    colNums = list(range(0,2*nGrps,2))
    Nsdf = resdf.take(colNums, axis=1)
    # display(Nsdf)
    evMins = Nsdf.min(axis=1)
    evKeeps = evMins.index[evMins >= min_group_size]  
    evKeeps = evKeeps.to_list()
#    print(evKeeps)
    
#   Handle situation of no events meeting min_group_size criterion
    if len(evKeeps) == 0:
        print('Kruskal-Wallis test not conducted: No events meet minimum group size criterion.')
        return( resdf)

#   Subset data to events for which kw test will be conducted.   
    kwdata = dataf.loc[evKeeps,:]
#    display(kwdata)

#   Set up objects for the call to scipy.stats.kruskal()   
    y = np.array(kwdata)
    label, idx = np.unique(grpdf['grp'], return_inverse=True)
    groups = [y[0,idx == i] for i, l in enumerate(label)]

    
#  Perform KW test for each event in kwdata
    for indx in range(0, len(kwdata) ):
        ev = kwdata.index[indx]
        groups = [y[indx,idx == i] for i, l in enumerate(label)] 
        resdf.loc[ev,'kwStat'], resdf.loc[ev,'kwPval'] = scipy.stats.kruskal(*groups)

    return( resdf )
    



    