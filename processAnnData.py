def processAnnData(annDataPath,min_genes=200,min_cells=3,gene_name='feature_name',quanthigh=.98,quantlow=.02,mt=20,highvar=.0125,MaxMean=3,MinMean=0.0125,MinDisp=0.5,MaxVal=10,n_neighbors=10,n_pcs=20,resolution=0.25):
    adata = sc.read_h5ad(annDataPath)
    result = gene_name in list(adata.var.columns)
    assert result == True
    sc.pp.filter_cells(adata, min_genes)#    sc.pp.filter_genes(adata, min_cells) 
    adata.var['mt'] =  adata.var[gene_name].str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, quanthigh)
    lower_lim = np.quantile(adata.obs.n_genes_by_counts.values, quantlow)
    adata = adata[(adata.obs.n_genes_by_counts < upper_lim) & (adata.obs.n_genes_by_counts > lower_lim)]
    adata = adata[adata.obs.pct_counts_mt < mt]
    sc.pp.normalize_total(adata, target_sum=1e4) #normalize every cell to 10,000 umi
    sc.pp.log1p(adata) #change to log counts
    sc.pp.highly_variable_genes(adata, min_mean=MinMean, max_mean=MaxMean, min_disp=MinDisp) #these are default values
    adata = adata[:,adata.var.highly_variable]
    return adata
