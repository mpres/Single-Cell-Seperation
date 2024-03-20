
#Mpresley 11/23/23, game plan 9:32, change the DF
'''
Test scoreCellidDF, as a single and save the DF, then do it again (the same process you did last time, scoreCellidArray ).
Do a manual resultDFIntersection manually. 


Mpresley 11/27/23, game plan 8:48 PM, Looking at resultDFIntersection, (let's go line by line and make comments to test out the code.
Once we have the intersection working, we'll integrate it to scoreCellidDf function"

'''


#Mpresley let's iniate this function

def resultDFIntersection(DF1,DF2):
    #Change validation for DataFrame type
    assert isinstance(DF1,pd.DataFrame)
    assert isinstance(DF2,pd.DataFrame)
    
    #Mpresley 11/27/23 will need to make sure the "shape" of the DF1 and DF2 are the same
    assert DF1.shape == DF2.shape
    
    #cell_id have to match
    #MPresley 11/27/23 the column name will have to be the same for the dataframe
#    assert array1[0][4] == array2[0][4]
    DF1 = DF1.reset_index()
    DF1 = DF1.drop('index',axis=1)
    DF2 = DF2.reset_index()
    DF2 = DF2.drop('index',axis=1)

    assert DF1['Cell_id'][0] == DF2['Cell_id'][0]
    #initiate result variable
    #Mpresley 11/27/23 this will need to be a dataframe result? we'll see the most efficient way to do this
    DF_temp = DF1.copy()

    for a in DF1.index:
    #both are positive
        if int(DF1.loc[a][0]) == 1 and int(DF2.loc[a][0]) == 1:
            pass
    #one is negative so the result is negative
        elif int(DF1.loc[a][2]) + int(DF2.loc[a][2]) > 0:
            DF_temp.at[a,"Negative"] = 1
            DF_temp.at[a,"Positive"] = 0
            DF_temp.at[a,"Neutral"] = 0
    #the only result left is neutral
        else:
            DF_temp.at[a,"Neutral"] = 1
            DF_temp.at[a,"Negative"] = 0
            DF_temp.at[a,"Positive"] = 0
    return DF_temp
    
    
    
