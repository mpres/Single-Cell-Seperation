
def makeGeneRow(adata, gene_list, purity_cutoff = 3, total_hit_cutoff=500, std=2,  i=4 ):
    '''   getGeneVariance and getGeneVarianceScores are dependant for this function
    adata is an annData type
    gene_list = gene list
    i = the number of genes that will be outputted with thresholds

    
    1. Validate data types 
    2. process list (sometimes it will be a pair list)
    3. create result dataframe
    4. loop through gene list, using geneVariance and build dataframe
    5. Sort DataFrame and return dataframe  
    6. i must be 2 are greater
    '''
    
    assert isinstance(adata,ad.AnnData)
    assert isinstance(gene_list,list) or isinstance(gene_list,np.ndarray)
    assert isinstance(i,int)
    assert isinstance(std,int)
    assert 'exclude_target' in list(adata.obs.columns)
    assert i >= 2
    
    #look for a list of pairs, and process it into one list
    if len(gene_list[0]) == 2:
        gene_list0 = [a[0] for a in gene_list ]
        gene_list1 = [a[1] for a in gene_list ]
        gene_list = gene_list0 + gene_list1
        gene_list = list(set(gene_list))
    elif len(gene_list[0]) == 15 and gene_list[1][0:3] == "ENS":
        gene_list = list(set(gene_list))
    else:
        print("Error, Gene_list variable is not a simple list")
        return
    
    gene_list_result = np.array([])
    thresh_gene_list_result = np.array([])
    
    
    #initiate threshhold list and gene list
    df_result,gene_purge_list = getGeneVarianceScores(adata,gene_list,len(gene_list),purity_cutoff=purity_cutoff,total_hit_cutoff=total_hit_cutoff)
#    gene_list = [x for x in gene_list if x not in gene_purge_list]
    print("gene_purge_list ",gene_purge_list)
    print('total_hit_cutoff', total_hit_cutoff)
    print('purity_cutoff', purity_cutoff)
    print('gene_count', len(gene_list))
    time.sleep(5)
    if len(df_result) == 0:
        print('df_result == 0 first iteration',df_result)
        return df_result
    gene_of_interest = df_result.index[0]
    gene_list_result = gene_of_interest
    
    variance_result = getGeneVariance(adata,gene_of_interest)
    bottom_thresh = variance_result[0][1]
    #        top_thresh = variance_result[0][2], Mpresley 2/8/24
    top_thresh = (bottom_thresh/2)
    top_thresh = variance_result[0][2]
    thresh_gene_list_result = np.array([[bottom_thresh,top_thresh,'+']])
    
    pos_index = getGeneVariance(adata,gene_of_interest,pos_idx=True)
    
    adata.obs.loc[list(adata.obs.iloc[pos_index].index),'exclude_target'] = 1
    
    #get first gene
    
    for j in range(i-1):
        print('total_hit_cutoff', total_hit_cutoff)
        print('purity_cutoff', purity_cutoff)
        print("gene_purge_list ",gene_purge_list)
        time.sleep(5)
        df_result,gene_purge_list = getGeneVarianceScores(adata,gene_list,len(gene_list),purity_cutoff=purity_cutoff,total_hit_cutoff=total_hit_cutoff)
#        gene_list = [x for x in gene_list if x not in gene_purge_list]
        print('gene_count', len(gene_list))
        if len(df_result) == 0:
            print("df_result = 0, got to iteration",j+2)
            return gene_list_result, thresh_gene_list_result
        gene_of_interest = df_result.index[0]
        print('gene_of_interest',gene_of_interest)
        print('iteration number - ',j+2)
        time.sleep(3)
        gene_list_result = np.append(gene_list_result,gene_of_interest)

        variance_result = getGeneVariance(adata,gene_of_interest)
        bottom_thresh = variance_result[0][1]
#        top_thresh = variance_result[0][2], Mpresley 2/8/24
        top_thresh = (bottom_thresh/2)
        thresh_array = np.array([[bottom_thresh,top_thresh,'+']])
        thresh_gene_list_result = np.vstack((thresh_gene_list_result,thresh_array))
    
        #Now omit the cells from the dataframe
    
        pos_index = getGeneVariance(adata,gene_of_interest,pos_idx=True)

        adata.obs.loc[list(adata.obs.iloc[pos_index].index),'exclude_target'] = 1
        
 #   Mpresley 2/19/24,        
#    adata.obs['exclude_target'] = 0
    
    gene_row_list = []
    
    print('Exited loop')
    print('gene_list_result',gene_list_result)
    
    for a in gene_list_result:        
        #gene_row_list.append(adata.var.loc[a][10])
        gene_row_list.append(adata.var.index.get_indexer([a])[0])
    gene_list_result = np.array(gene_row_list)
    
    return gene_list_result, thresh_gene_list_result



