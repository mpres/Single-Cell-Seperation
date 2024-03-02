#Mpresley 11/23/23, remaking scoreCellid using the scoreCellDF function

def scoreCellidDF(annData,cell_id,gene_index=[],thresh=[[.8,.1,'+']]):
    ''' scoreCellid is a function that is called inside of a python function called scoreGenes, scoreCellid, will take an a cell_id, a gene_index and a 
thresh array np.arrays 
annData is data type for scipi and single cell data see (https://anndata.readthedocs.io/en/latest/)

cell_id will be a single text value, example with a format of tissue+_+celltype

genes_index will be a numpy array of numpy arrays, pointing to the gene index of interest, the columns will work as "OR" statement and the rows as "AND" statement

thresh will be a numpy array of numpy arrays, with corresponding index/rows of the index/rows for gene_index, the thresh will determine the thresholds for whats positive, 
and what is negative and what is neutral.

sign, will determine if  the gene thresh should exclude the cell or add the cell if it hids the gene type. sign will have the same

OUTPUT: will be a list of positive, negative and neutral ## 11/13/23 with the cell_id as well

OUTLINE: 1. validate datatypes of parameter (make sure everything checks out)
         2. create the result array
         3. get the cell_id_index
         4. loop through the cell_id_index
             4.a loop through the cell_id_index usingscoreCell to build up the resultarray,
             4.b add the cell_id to each row of the resultarray 
'''
        #validation of datatypes
    assert isinstance(annData,ad.AnnData)
    assert isinstance(gene_index,np.ndarray)
    assert isinstance(thresh,np.ndarray)
    #Mpresley 2/23/24 set sign variabel
#    sign = thresh[0][2]
#    assert isinstance(sign,np.ndarray) or 
    
#    print('gene_index',gene_idx,'len(gene_index)',len(gene_index))
#    print('thres',thresh,'len(thres)',len(thresh))
    assert len(gene_index) == len(thresh)
    assert len(gene_index) > 0
    
    #find the shpae of gene_index and threshold array
    if isinstance(gene_index[0],np.int64):
        gene_shape = "row"
        print("printing gene_shape",gene_shape)
    elif isinstance(gene_index[0][0],np.int64):
        gene_shape = "table"
        print('gene_shape',gene_shape)
    else:
        print("gene_index array is not valid, make sure the shape and values are correct")
        return -1

    if (isinstance(thresh[0][0],str)):
        thresh_shape = "row"
        print("thresh_shape",thresh_shape)
    elif (thresh[0][0][2] == "+" or thresh[0][0][2] == "-"):
        thresh_shape = "table"
        print("thresh_shape",thresh_shape)
    else:
        print("thresh array is not valid, make sure the shape and values are correct")
        return -1
        
        
        
    #Validate the shape of the threshold index and gene_index
    if gene_shape == "row" and thresh_shape == "row":
        assert len(gene_index) == len(thresh)
    elif gene_shape != thresh_shape :
        print("Error, the shape of gene index and threshold array do not match up")
        return -1
    else:
        assert len(gene_index) == len(thresh)
    
    
    resultDF = pd.DataFrame(columns=['Positive','Neutral','Negative','Cell_idx','Cell_id'])
    
    cellid_index = list(np.where(annData.obs['cell_id'] == cell_id))
#Mpresley 11/18/23 prime this 
    if gene_shape == "row":
#        print('cellid_index[0]',cellid_index[0],'len(cellid_index[0])',len(cellid_index[0]))
        print('processing cell_id',cell_id)
        for c in cellid_index[0]:
            r = scoreCellDF(annData,c,gene_index,thresh,cell_id)
            resultDF = pd.concat([resultDF,r])
        return resultDF
    else:
        for c in cellid_index[0]:
            r = scoreCellDF(annData,c,gene_index[0],thresh[0],cell_id)
            resultDF = pd.concat([resultDF,r])
        
        #now that we have resultArray (the first gene set, we'll grab the second one and add it one to to using resultArrayIntersection(array1,array2)
        #Mpresley 11/23/23, work on getting a resultArrayIntersection version using a DataFrame
        for i in range(len(gene_index)-1):
            resultSingleDF = pd.DataFrame(columns=['Positive','Neutral','Negative','Cell_idx','Cell_id'])
            for c in cellid_index[0]:
                r = scoreCellDF(annData,c,gene_index[i+1],thresh[i+1],cell_id) 
                resultSingleDF = pd.concat([resultSingleDF,r])
#                resultSingleArray.append(r)
            resultDF = resultDFIntersection(resultDF,resultSingleDF)
        return resultDF
