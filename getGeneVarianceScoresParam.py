''' mpresley 1/7/24 
  working on step 4 from 1/5/24, (create a function that grabs the most interesting genes)'
  This function will get a list of genes are a list of paired genes and the a unique list of them (maybe 100 genes are so) may be more
  1. It will return a dataframe with the top score (probably returb the percentage hit, plus the amount of hit
  2. It will take a list of genes are an amount of Genes, a total iteration rate (how many of the list do you want
  3. getGeneScores(adata,gene_list,i=100)
  
'''

''' 
Mpresley 3/4/24 rewrite- 1. getGeneVarianceScores/2 to allow for a parameter for GeneScore, 
'''

from IPython.display import clear_output

def getGeneVarianceScoresParam(adata,gene_list,i=100,std=2,purity_cutoff=2,total_hit_cutoff=100,score_cutoff=1,GeneScoreParam=None):
    '''
    getGeneVariance is a dependant for this function
    adata is an annData type
    gene_list = gene list
    i = the amount of genes that will be analyzed from the list
    std = will be used for getGeneVariance
    purity_cutoff is the amount of purity ratio needed to be added to the dataframe result
    same with total_hit_cutoff and score_cutoff (except regarding Total_Target_Hit and Score) 
    GeneScoreParam will determine will weight the geneScore variable with respect to purity_score and the total_target_hit.
        the function goes GeneScore = (Total_Target_hit * GeneScoreParam) * purity_score) unless GeneScoreParam = null or GeneScore = 0,
        in that case Genescore = purity_score
    
    1. Validate data types 
    2. process list (sometimes it will be a pair list)
    3. create result dataframe
    4. loop through gene list, using geneVariance and build dataframe
    5. Sort DataFrame and return dataframe 
    '''

    assert isinstance(adata,ad.AnnData)
    assert isinstance(gene_list,list) or isinstance(gene_list,np.ndarray)
    assert isinstance(i,int)
    assert isinstance(std,int)
    
    print('gene_list count',len(gene_list))
    
    print('GeneScoreParam ',GeneScoreParam)
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
    
    
    # Gene_score is the ratio of (ratio_of_target_hits/ratio_of_non_traget_hits) * target_hits
    # ratio_of_target_hits = target_hits/(total target cells)
    # ratio_of_non_target_hits = target_hits/(total target cells)
    ResultDF = pd.DataFrame(columns=['Total_Target_Hit','Total_Non_Target_hit','Purity_Score','Gene_Score'])
    Gene_purge_list = []

    j = 0
#    print("i",i)
#    print('gene_list count',len(gene_list))
    for g in gene_list:
        if j == i:
            break
        r = getGeneVariance(adata,g,std)
 #       print('record',r)
        clear_output(wait=True)
        print(j)
        Total_Target_hit = r[0][4]
        Total_Non_Target_hit = r[1][4]
        target_ratio = r[0][3] + .001
        non_target_ratio = r[1][3] + .001
        purity_ratio = ( target_ratio / non_target_ratio )


        #Mpresley 3/4/24 add new code here
        if GeneScoreParam == None or GeneScoreParam == 0:
            print("GeneScoreParam has been notified",GeneScoreParam)
            Gene_score = purity_ratio
        else:
            Gene_score = (Total_Target_hit * GeneScoreParam * .01) + purity_ratio
#        Gene_score = (purity_ratio * Total_Target_hit)
        if Gene_score >= score_cutoff and purity_ratio >= purity_cutoff and Total_Target_hit >= total_hit_cutoff:
            ResultDF.loc[g] = [Total_Target_hit,Total_Non_Target_hit,purity_ratio,Gene_score]
        elif(Gene_score < score_cutoff):
            Gene_purge_list.append(g)
            
        j += 1
    
    
    return ResultDF.sort_values('Gene_Score',ascending=False),Gene_purge_list
