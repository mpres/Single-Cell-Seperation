''' 12/18/23 remake scoreCellDF to let "or" to actually work
'''

''' Mpresley 12/18/23 , rewriting of scoreCells to use panda DataFrames
'''

def scoreCellDF(adata,scoreCellIndex,geneIndex,threshHold,Cell_id):
        #Validation,
        #Validate dataType
        assert isinstance(adata,ad.AnnData)
        assert isinstance(geneIndex,np.ndarray)
        assert isinstance(threshHold,np.ndarray)
#        assert isinstance(sign,np.ndarray)
        assert len(geneIndex) == len(threshHold)
        assert len(geneIndex) > 0
        
        
        
        #Look for positive cell 
        for i in range(len(geneIndex)):
            g = geneIndex[i]
            thres_positive = float(threshHold[i][0])
            thres_negative = float(threshHold[i][1])
            thres_sign = threshHold[i][2]
            cellGeneValue = adata.X[scoreCellIndex,g]



            #Look for a positive value (either for a "+" or "-" sign
            if thres_sign == "+":
                if cellGeneValue >= thres_positive:
                    return pd.DataFrame(np.array([[1,0,0,scoreCellIndex,Cell_id]]),columns=['Positive','Neutral','Negative','Cell_idx','Cell_id'])
#                    return [1,0,0,scoreCellIndex]
            elif thres_sign == "-":
                if cellGeneValue <= thres_positive:
                    return pd.DataFrame(np.array([[1,0,0,scoreCellIndex,Cell_id]]),columns=['Positive','Neutral','Negative','Cell_idx','Cell_id'])
                    #return [1,0,0,scoreCellIndex]


       #Look for a neutral value
        for i in range(len(geneIndex)):
            g = geneIndex[i]
            thres_positive = float(threshHold[i][0])
            thres_negative = float(threshHold[i][1])
            thres_sign = threshHold[i][2]
            cellGeneValue = adata.X[scoreCellIndex,g]
            
            if thres_sign == "+":
                if cellGeneValue > thres_negative:
                    return pd.DataFrame(np.array([[0,1,0,scoreCellIndex,Cell_id]]),columns=['Positive','Neutral','Negative','Cell_idx','Cell_id'])
                   # return [0,1,0,scoreCellIndex]
            
            elif thres_sign == "-":
                if cellGeneValue < thres_negative:
                    return pd.DataFrame(np.array([[0,1,0,scoreCellIndex,Cell_id]]),columns=['Positive','Neutral','Negative','Cell_idx','Cell_id'])
           #         return [0,1,0,scoreCellIndex]
    
    
    
            #Look for a negative cell value
        for i in range(len(geneIndex)):
            g = geneIndex[i]
            thres_positive = float(threshHold[i][0])
            thres_negative = float(threshHold[i][1])
            thres_sign = threshHold[i][2]
            cellGeneValue = adata.X[scoreCellIndex,g]
          
            if thres_sign == "+":
                if cellGeneValue <= thres_negative:
                    return pd.DataFrame(np.array([[0,0,1,scoreCellIndex,Cell_id]]),columns=['Positive','Neutral','Negative','Cell_idx','Cell_id'])
                   # return [0,0,1,scoreCellIndex]
            
            elif thres_sign == "-":
                if cellGeneValue > thres_negative:
                    return pd.DataFrame(np.array([[0,0,1,scoreCellIndex,Cell_id]]),columns=['Positive','Neutral','Negative','Cell_idx','Cell_id'])
               #     return [0,0,1,scoreCellIndex]
    
