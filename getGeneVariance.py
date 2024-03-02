
def getGeneVariance(adata,gene_idx,std=2,pos_idx=False):
        assert isinstance(adata,ad.AnnData)
        assert isinstance(gene_idx,str)
        assert 'target' in adata.obs.columns
        assert gene_idx in adata.var.index
        
        if 'exclude_target' in adata.obs.columns:
            #Mpresley 1/5/24 the exclude target field will not include cells that have been marked with a 1. This is useful when you want to find subsquent target cells and they were previous false negatives
            target_idx = list(adata.obs.index.get_indexer(adata.obs[(adata.obs['target'] == '1') & (adata.obs['exclude_target'] != 1) ].index))
            #Mpresley 2/8/24 added the exclude target to non_target_idx, this should speed of calculation (if it's used for negatives and focus purity_scores
            non_target_idx = list(adata.obs.index.get_indexer(adata.obs[(adata.obs['target'] == '0') & (adata.obs['exclude_target'] != 1) ].index))
            
        else:
            #Mpresley 12/11/23
            target_idx = list(adata.obs.index.get_indexer(list(adata.obs[adata.obs['target'] == '1'].index)))
            non_target_idx = list(adata.obs.index.get_indexer(list(adata.obs[adata.obs['target'] == '0'].index)))
        
#        gene_idx_num = adata.var.index.get_loc('ENSG00000166819')
        gene_idx_num = adata.var.index.get_loc(gene_idx)
        
        #Get Target Matrix, we'll convert this in a a numpy.array so we can get the variance and mean of the target
        target_Matrix = adata.X[target_idx,gene_idx_num]
        non_target_Matrix = adata.X[non_target_idx,gene_idx_num]
        #Mpresley 12/12/23 the target_list will only include the rows that have values above 0
        target_list = [ float(x) for x in target_Matrix.A if x > 0 ]
        non_target_list = [ float(x) for x in non_target_Matrix.A if x > 0 ]
        
        #Mpresley 1/5/24
        if pos_idx == True:
            target_idx_pos = [ int(x) for x in target_idx if adata.X[x,gene_idx_num] > 0 ]
            target_idx_false_pos = [ int(x) for x in non_target_idx if adata.X[x,gene_idx_num] > 0 ]
        
        target_list = np.array(target_list)
        non_target_list = np.array(non_target_list)
        
                
        #Mpresley 12/12/23, get the positive and negative range of target 
        if len(target_list) == 0:
            target_neg_range = 0
            target_pos_range = 0
        else:
            target_neg_range = np.mean(target_list) -  (np.std(target_list) * std)
            target_pos_range = np.mean(target_list) +  (np.std(target_list) * std)
        
        
        #Mpresley 12/12/23, get the positive and negative range of target 
        if len(non_target_list) == 0:
            non_target_neg_range = 0
            non_target_pos_range = 0
        else:
            non_target_neg_range = np.mean(non_target_list) -  (np.std(non_target_list) * std)
            non_target_pos_range = np.mean(non_target_list) +  (np.std(non_target_list) * std)
        
        
        total_target_cells = len(adata.obs[adata.obs['target'] == '1'])
        total_active_target_cells = len(target_list)
        
        target_active_ratio = total_active_target_cells/total_target_cells
        
            
        total_non_target_cells = len(adata.obs[adata.obs['target'] == '0'])
        total_non_active_target_cells = len(non_target_list)
        non_target_active_ratio = total_non_active_target_cells/total_non_target_cells
        
        
        if pos_idx != True:
            return ('target',target_neg_range,target_pos_range,target_active_ratio,total_active_target_cells),('non_target',non_target_neg_range,non_target_pos_range,non_target_active_ratio,total_non_active_target_cells)
        else:
            return target_idx_pos
        
