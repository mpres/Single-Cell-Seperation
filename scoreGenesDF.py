#Mpresley 12/1/23 8:48 PM
'''
Looks like scoreCellidDF and scoreCellid is working pretty,
the game plan to day is to make a function that sums up a DataFrame
and then to outline the scoreGenes to use scoreCellidDF.
'''
#ResultDF3['Positive'].sum() + ResultDF3['Negative'].sum() + ResultDF3['Neutral'].sum()
#sum off all these is 14025
#len(ResultDF3)
#length of the Datafrme is 14025
#ResultDF3['Positive'].sum()
#sum of ResultDF3 is 1534
#ResultDF2['Positive'].sum()
#sum of ResultDF2 is 3592
#ResultDF['Positive'].sum()
#ResultDF4= resultDFIntersection(ResultDF2,ResultDF3)
#ResultDF4['Positive'].sum()
#after the intersection the positive is only 924 (which makes sense, it should reduce the amount)
#ResultDF5 = resultDFIntersection(ResultDF4,ResultDF3)
#ResultDF5['Positive'].sum()

### ResultDF3['Positive'].sum()
''' Mpresley 12/1/23 9:05 summing columns was pretty easy, now let's take a look at scoreGenes,
 Probably the prep work will be pretty similar, but the other elements will changes quite a bit to

Outline 

ScoreGenes(adata,gene_index,thresh_list)

    Notes on the parameters 
    adata is an annData, Datatype = AnnData
    gene_index, wil be a index of the annData.var (or adata.var). This will help gather the genes that are relevant for the analysis, Datatype = Array
    gene_index may also be an Array of Array ( all genes in a single work as "or" statements, while genes between rows work as "and" statements.
    thresh_list, will be an Array of Arrays (or list of list). structure will be [positive_limit,negative_limit,sign]. exampe [[.8,.2,'+'],[[.1,.6,'-']]]
    if gene_index is list of list, then thresh_list will be list of list of list (Yea, thats kinds of annoying, but necessary)
    
    1. Validation of Datatypes
    2. Validation of param shaps (for gene_index and thresh_list
    3. Validation of Cell_list, get from annData.obs (unique values)

    4. get Cell_list_ids (combination of cell and tissue types)
    5. Make the Result Dataframe
        pd.DataFrame("Positive,Nuetral,Negative,Cell_id"), this data will be summed via Cell_id.
    
    6. Loop through Cell_list_ids
        should get a DF back that will be summed and appended onto the result Df
    7. return the result
    
    Mpresley 2/2/24: (add a parameter to scoreGenes called "summarize". this variable will equal to TRUE, and will return the summaries of each cellId, however if summrize is set to
                     "FALSE" then the full Dataframe will return (so each in this expanded dataframe would represent 1 biological cell). These expanded, unsummarized DataFrames
                      will also have the corresponding indexes which will help me isolate the false positive and false negatives index value and and set thme as a "target".
        
    
'''

def scoreGenesDF(adata,gene_index,thresh_list,summarize=True):
        #validate datatype
        assert isinstance(adata,ad.AnnData)
        assert isinstance(gene_index,np.ndarray)
        assert isinstance(thresh_list,np.ndarray)
        #validate the shape
        assert len(gene_index) == len(thresh_list)
        assert len(gene_index) > 0
        assert thresh_list.shape[0] == gene_index.shape[0]
        #get cell_list_id
        cell_id_list = list(adata.obs['cell_id'].unique())
        #validate cell_list-id
        assert len(cell_id_list) > 0
        assert isinstance(cell_id_list[0],str)
        #Mpresley 12/2/23, create result DataFrame
        if summarize == True:
            ResultDF = pd.DataFrame(columns=['Total_Positive','Total_Neutral','Total_Negative','Cell_id'])
            for c in cell_id_list:
                #MPresley 2/26/24 return parametes to test scoreCellidDF
#                return c,gene_index,thresh_list
                DF_temp = scoreCellidDF(adata,c,gene_index,thresh_list)
                DF_temp_total = pd.DataFrame([[pd.to_numeric(DF_temp['Positive']).sum(),pd.to_numeric(DF_temp['Neutral']).sum(),pd.to_numeric(DF_temp['Negative']).sum(),c]], columns=['Total_Positive','Total_Neutral','Total_Negative','Cell_id'])
                ResultDF = pd.concat([DF_temp_total,ResultDF])
        else:
            ResultDF = pd.DataFrame(columns=['Positive','Neutral','Negative','Cell_idx','Cell_id'])
            for c in cell_id_list:
                DF_temp = scoreCellidDF(adata,c,gene_index,thresh_list)
                ResultDF = pd.concat([DF_temp,ResultDF])
            
        return ResultDF
        #return pd.to_numeric(DF_temp['Positive']).sum(),pd.to_numeric(DF_temp['Neutral']).sum(),pd.to_numeric(DF_temp['Negative']).sum()
        
        


        #DF_2 = pd.DataFrame([[DF_temp['Positive'].sum()],[DF_temp['Neutral'].sum()],[DF_temp['Negative'].sum()],cell_id_list[1]], columns = ['Total_Positive','Total_Neutral','Total_Negative','Cell_id'])
