# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 13:12:31 2019

@author: LuoHan
"""
from scipy.optimize import nnls
import pandas as pd
import numpy as np


##-----------------Loading Function----------------
def Get_SC_Meta_Info(fil):
    """
    """
    return pd.read_table(fil, header = 0)

def Get_SC_Count_Info(fil):
    """
    """
    return pd.read_table(fil, header = 0, index_col = 0)

def Get_bulk_Meta_Info(fil):
    """
    """
    return pd.read_table(fil, header = 0)

def Get_bulk_Count_Info(fil):
    """
    """
    return pd.read_table(fil, header = 0, index_col = 0)


##--------------Class for bulk data preparation------------

class BulkMatrix:
    """
        Calculate relative abudance
    """
    def __init__(self, X, by_col = True):
        
        self.X = X
        self.by_col = by_col
        self.DataMatrixConstruction()
    
    def DataMatrixConstruction(self):
        """
            abudance computation.
        """
        if (self.X < 0).sum()  > 0:
            raise Exception ('Negative entry appears')
        
        if self.by_col:
            self.X = self.X / self.X.sum(axis = 0)
        else:
            self.X = self.X.T / self.X.sum(axis = 1)
            self.X = self.X.T

##---------------Class for single cell data preparation------------

class SingleCellMatrix:
    """
        Single Cell Matrix Data preprocessing.
        Do as what MuSiC do.
    """
    
    def __init__(self, Meta, Count, nonzero = True, markers = None,
                 select_ct = None, ct_cov = False):
        """
            Meta : dataframe(pandas)
                   Phenotype Data.
            
            count : dataframe(pandas)
                    Single Cell Count Data.
            
            nonzero : bool
                      if True, remove all gene with zero expression.
            
            markers : list(array)
                      list or array of gene names. Default as None. 
                      If None, then use all genes provided.
            
            select_ct : list(arrat)
                        list or array of cell type. Default as None.
                        if None, then use all cell types provided.
            
        """
        self.Meta = Meta
        self.Count = Count
        self.nonzero = nonzero
        self.markers = markers
        self.select_ct = select_ct
        self.ct_cov = ct_cov
        
        #------Data preparation------
        self.DataConstruction()
        
    
    def DataConstruction(self):
        """
            Prepare Design matrix and Cross-subject Variance for  Deconvolution
            This function is used for generating cell type specific cross-subject 
            mean and variance for each gene. Cell type specific library size is 
            also calcualted.
        """
        
        if self.select_ct is not None:
            ## select cell types
            mask = self.Meta['cellType'].isin(self.select_ct)
            self.Meta = self.Meta[mask]
            self.Count = self.Count.loc[:,self.Meta['sampleName']]        

        
        if self.nonzero == True:
            ## eliminate non expressed genes
            nz_gene = self.Count.sum(axis = 1) != 0
            self.Count = self.Count[nz_gene]
            
        celltypes= self.select_ct
        samples = set(self.Meta['sampleID'])
        
        ##---------Creating Relative Abudance Matrix----------
        print('Creating Relative Abudance Matrix ...')
        self.theta = self.Theta(celltypes, samples)
        
        ##---------Creating Variance Matrix---------
        print('Creating Variance Matrix ...')
        if self.ct_cov == True:
            self.sigma_ct = self.Sigma(celltypes, samples)
        else:
            self.sigma = self.Sigma(celltypes, samples)
        
        ##---------Creating Library Size Matrix-----
        print('Creating Library Size Matrix ...')
        self.LS = self.LibrarySize(celltypes, samples)
        self.MLS = self.LS.mean(axis=0)
        
        ##---------Creating Design Matrix---------
        print('Creating Design Matrix')
        self.DM = np.zeros(self.theta.shape)
        theta = np.array(self.theta)
        for i in range(len(celltypes)):
            self.DM[:,i] = theta[:,i] * self.MLS[i]
        
        print('markers genes selecting ...')
        if self.markers is None:
            Genes = self.Count.index
            if self.ct_cov != True:
                self.sigma = pd.DataFrame(self.sigma, index = Genes, columns = celltypes)
            else:
                self.sigma_ct = pd.DataFrame(self.sigma_ct, columns = Genes)
            self.DM = pd.DataFrame(self.DM, index = Genes, columns = celltypes)
            self.theta = pd.DataFrame(self.theta, index = Genes, columns = celltypes)
            
        else:
            mask = self.Count.index.isin(self.markers)
            self.Count = self.Count[mask]
            Genes = self.Count.index
            self.DM = self.DM[mask]
            self.theta = self.theta[mask]
            if self.ct_cov != True:
                self.sigma = self.sigma[mask]
                self.sigma = pd.DataFrame(self.sigma, index = Genes, columns = celltypes)
            else:
                self.sigma_ct = self.sigma_ct.T[mask].T
                self.sigma_ct = pd.DataFrame(self.sigma_ct, columns = Genes)
            
            self.DM = pd.DataFrame(self.DM, index = Genes, columns = celltypes)
            self.theta = pd.DataFrame(self.theta, index = Genes, columns = celltypes)
        
        print('Data preparation Done!')
 
           
    def Theta(self, celltypes, samples):
        """
            Prepare theta matrix.
        """
        theta = []
        for ct in celltypes:
            cross_subject = []
            
            for sid in samples:
                
                mask_ct = self.Meta['cellType'] == ct
                mask_sid = self.Meta['sampleID'] == sid
                mask = np.array(mask_ct & mask_sid)
                columns = self.Meta[mask]['sampleName']
                
                sub_ct_sid = self.Count.loc[:,columns]
                sub_ct_sid = sub_ct_sid.sum(axis=1) / sub_ct_sid.sum().sum()
                cross_subject.append(sub_ct_sid.tolist())
            
            cross_subject = np.array(cross_subject)
            cross_subject_mean = cross_subject.mean(axis = 0)
            theta .append(cross_subject_mean.tolist())
        
        
        theta = np.array(theta).T  #reshape (gene * cell type)
        
        return theta
        
    
    def Sigma(self, celltypes, samples):
        """
            Prepare sigma matrix.
        """
        
        if self.ct_cov == True:
            
            sigma_ct = []
            
            #use the covariance across cell types
            n_Genes = len(self.Count.index)
            
            sigma = None
            for ct in celltypes:
                cross_subject = np.array([])
                
                for sid in samples:
                    
                    mask_ct = self.Meta['cellType'] == ct
                    mask_sid = self.Meta['sampleID'] == sid
                    mask = np.array(mask_ct & mask_sid)
                    columns = self.Meta[mask]['sampleName']
                    
                    sub_ct_sid = self.Count.loc[:,columns]
                    sub_ct_sid = np.array(sub_ct_sid.sum(axis=1) / sub_ct_sid.sum().sum())
                    cross_subject = np.hstack((cross_subject,sub_ct_sid))
                
                if sigma is None:
                    sigma = cross_subject
                else:
                    sigma = np.vstack((sigma, cross_subject))
            
            n_subs = len(samples)
            
            for i in range(n_Genes):
                
                index = [i + n_Genes*j for j in range(n_subs)]
                
                sub_M = sigma[:,index]
                sub_cov = np.cov(sub_M)
                sub_cov = sub_cov.reshape((sub_cov.shape[0]*sub_cov.shape[1],))
                
                sigma_ct.append(sub_cov)
            
            sigma_ct = np.array(sigma_ct).T
            
            return sigma_ct
        
        else:
            
            Sigma = []
            for ct in celltypes:
                cross_subject = []
            
                for sid in samples:
                
                    mask_ct = self.Meta['cellType'] == ct
                    mask_sid = self.Meta['sampleID'] == sid
                    mask = np.array(mask_ct & mask_sid)
                    columns = self.Meta[mask]['sampleName']
                
                    sub_ct_sid = self.Count.loc[:, columns]
                    sub_ct_sid = sub_ct_sid.sum(axis=1) / sub_ct_sid.sum().sum()
                    cross_subject.append(sub_ct_sid.tolist())
            
                cross_subject = np.array(cross_subject)
                cross_subject_var = cross_subject.var(axis = 0, ddof=1)
                Sigma.append(cross_subject_var.tolist())
                
            Sigma = np.array(Sigma).T
            
            return Sigma
        
        
    def LibrarySize(self, celltypes, samples):
        """
            prepare Library Matrix.
        """
        
        LS = []
        
        for ct in celltypes:
            cross_subject = []
            
            for sid in samples:
                mask_ct = self.Meta['cellType'] == ct
                mask_sid = self.Meta['sampleID'] == sid
                mask = np.array(mask_ct & mask_sid)
                columns = self.Meta[mask]['sampleName']
                
                sub_ct_sid = self.Count.loc[:,columns]
                
                cross_subject.append(sub_ct_sid.sum().mean())
            
            LS.append(cross_subject)
        
        LS = np.array(LS)
        LS = pd.DataFrame(LS.T, index = samples, columns = celltypes)
        
        return LS
      

##-----------class for NNLS Model-------------
class NNLS:
    """
        simple expands for nnls in scipy.
    """
    def __init__(self, X, y):
        """
            Initial X an y
        """
        self.X = X
        self.y = y
        self.fit()
        self.resid()
    
    def fit(self):
        """
            fitted
        """
        self.coef, self.score = nnls(self.X, self.y)
    
    def resid(self):
        """
            Do as R do ”resid(model)“
        """
        self.r = self.y - np.dot(self.X, self.coef)


#--------------class for Results store-------------
class ModelResults:
    """
        Model results in store
    """
    def __init__(self, initial_p, initial_coef, initial_resid,
                 weight, weight_p, weight_coef, weight_resid,
                 R2, var_p):
        """
        """
        self.initial_p = initial_p
        self.initial_coef = initial_coef
        self.initial_resid = initial_resid
        self.weight = weight
        self.weight_p = weight_p
        self.weight_coef = weight_coef
        self.weight_resid = weight_resid
        self.R2 = R2
        self.var_p = var_p


def weight_cal(Sp, Sigma):
    """
        Calculate weight with cross-subject variance for each cell types
    """
    return (Sp**2 * Sigma).sum(axis=1)


def weight_cal_ct(Sp, Sigma_ct):
    """
        Calculate weight with cross cell type covariance
    """
    
    weights = []
    
    nGenes = Sigma_ct.shape[1]
    n_ct = len(Sp)
    Sp = Sp.values
    Sp_2 = Sp * Sp.reshape((Sp.shape[0],1))
    
    for i in range(nGenes):
        sig = Sigma_ct.iloc[:,i].values.reshape((n_ct, n_ct))
        weights.append((Sp_2 * sig).sum())
    
    return np.array(weights)
        

def music_basic(Y, X, S, sigma, iter_max, nu, eps, centered, normalize):
    """
        weight is estimated with cell type vraiance
    """
    if centered:
        X = X - X.values.mean()
        Y = Y - Y.values.mean()
    
    if normalize:
        X = X / X.values.std()
        S = S * S.values.std()
        Y = Y / Y.values.std()
    else:
        Y = Y * 100

    lm_D = NNLS(X, Y)
    r = lm_D.r
    weight_gene = 1 / (nu + r**2 + weight_cal(Sp=lm_D.coef * S, Sigma= sigma))
    Y_weight = Y.mul(np.sqrt(weight_gene), axis=0)
    X_weight = X.mul(np.sqrt(weight_gene), axis=0)
    lm_D_weight = NNLS(X_weight, Y_weight)
    
    p_weight = lm_D_weight.coef
    p_weight = p_weight / p_weight.sum()
    r = lm_D_weight.r     
    
    print ('    Iteraction NNLS start ...')    
    for i in range(iter_max):
        weight_gene = 1 / (nu + r**2 + weight_cal(Sp=lm_D_weight.coef * S, Sigma= sigma))
        Y_weight = Y.mul(np.sqrt(weight_gene), axis=0)
        X_weight = X.mul(np.sqrt(weight_gene), axis=0)
        lm_D_weight = NNLS(X_weight, Y_weight)
        p_weight_new = lm_D_weight.coef / lm_D_weight.coef.sum()
        r_new = lm_D_weight.r

        if (sum(abs(p_weight_new - p_weight))) < eps:
            print ('    Done, save results.')
            p_weight = p_weight_new
            r = r_new
            R2 = 1 - (Y - np.dot(X.values, lm_D_weight.coef)).var() / Y.var()
            var_p = np.linalg.inv(np.dot(X_weight.T, X_weight)).diagonal() \
                    * (r**2).mean() /(lm_D_weight.coef.sum()**2)
            
            results = ModelResults(initial_p = (lm_D.coef / lm_D.coef.sum()),
                                   initial_coef = lm_D.coef,
                                   initial_resid = lm_D.r,
                                   weight = weight_gene,
                                   weight_p = p_weight,
                                   weight_coef = lm_D_weight.coef,
                                   weight_resid = r_new,
                                   R2 = R2,
                                   var_p = var_p)
            return results
        
        p_weight = p_weight_new
        r = r_new
    print ('    Done, save results.')
    R2 = 1 - (Y - np.dot(X.values, lm_D_weight.coef)).var() / Y.var()
    var_p = np.linalg.inv(np.dot(X_weight.T, X_weight)).diagonal() \
            * (r**2).mean() /(lm_D_weight.coef.sum()**2)
    
    results = ModelResults(initial_p = (lm_D.coef / lm_D.coef.sum()),
                           initial_coef = lm_D.coef,
                           initial_resid = lm_D.r,
                           weight = weight_gene,
                           weight_p = p_weight,
                           weight_coef = lm_D_weight.coef,
                           weight_resid = r_new, 
                           R2 = R2,
                           var_p = var_p)
    return results        


def music_basic_ct(Y, X, S, sigma_ct, iter_max, nu, eps, centered, normalize):
    """
        weight is estimated with cell type covariance
    """
    if centered:
        X = X - X.values.mean()
        Y = Y - Y.values.mean()
    
    if normalize:
        X = X / X.values.std()
        S = S * S.values.std()
        Y = Y / Y.values.std()
    else:
        Y = Y * 100
    
    lm_D = NNLS(X, Y)
    weights  = weight_cal_ct(Sp=lm_D.coef * S, Sigma_ct=sigma_ct)
    weight_gene = 1 / (nu + lm_D.r ** 2 + weights)
    Y_weight = Y.mul(np.sqrt(weight_gene), axis=0)
    X_weight = X.mul(np.sqrt(weight_gene), axis=0)
    
    lm_D_weight = NNLS(X_weight, Y_weight)
    
    p_weight = lm_D_weight.coef
    p_weight = p_weight / p_weight.sum()
    r = lm_D_weight.r
    
    print ('    Iteraction NNLS start ...')
    for i in range(iter_max):
        weight_gene = 1 / (nu + r ** 2 + weight_cal_ct(lm_D_weight.coef * S, Sigma_ct=sigma_ct))
        Y_weight = Y.mul(np.sqrt(weight_gene), axis=0)
        X_weight = X.mul(np.sqrt(weight_gene), axis=0)
        lm_D_weight = NNLS(X_weight, Y_weight)
        p_weight_new = lm_D_weight.coef / lm_D_weight.coef.sum()
        r_new = lm_D_weight.r
        if (sum(abs(p_weight_new - p_weight))) < eps:
            print ('    Done, save results.')
            p_weight = p_weight_new
            r = r_new
            R2 = 1 - (Y - np.dot(X.values, lm_D_weight.coef)).var() / Y.var()
            var_p = np.linalg.inv(np.dot(X_weight.T, X_weight)).diagonal() \
                    * (r**2).mean() /(lm_D_weight.coef.sum()**2)
            
            results = ModelResults(initial_p = (lm_D.coef / lm_D.coef.sum()),
                                   initial_coef = lm_D.coef,
                                   initial_resid = lm_D.r,
                                   weight = weight_gene,
                                   weight_p = p_weight,
                                   weight_coef = lm_D_weight.coef,
                                   weight_resid = r_new,
                                   R2 = R2,
                                   var_p = var_p)
            return results
        
        p_weight = p_weight_new
        r = r_new
    print ('    Done, save results.')
    R2 = 1 - (Y - np.dot(X.values, lm_D_weight.coef)).var() / Y.var()
    var_p = np.linalg.inv(np.dot(X_weight.T, X_weight)).diagonal() \
            * (r**2).mean() /(lm_D_weight.coef.sum()**2)
    
    results = ModelResults(initial_p = (lm_D.coef / lm_D.coef.sum()),
                           initial_coef = lm_D.coef,
                           initial_resid = lm_D.r,
                           weight = weight_gene,
                           weight_p = p_weight,
                           weight_coef = lm_D_weight.coef,
                           weight_resid = r_new, 
                           R2 = R2,
                           var_p = var_p)
    return results
       
                    
def music_prop(bulk_Meta, bulk_Count, sc_Meta, sc_Count, markers = None,
               clusters = 'cellType', samples = 'sampleID', select_ct = None,
               ct_cov = False, iter_max = 1000, nu = 0.0001, eps = 0.01, 
               centered = False, normalize = False):
    """
        MuSiC Deconvolution.
    """
    #-----------Data preprocessing-----------
    Nz_genes = bulk_Count.mean(axis=1) != 0
    bulk_Count = bulk_Count[Nz_genes]
    bulk_gene = bulk_Count.index
    
    if markers is None:
        sc_markers = bulk_gene
    else:
        sc_markers = list(set(markers) & set(bulk_gene))
    
    sc_basis = SingleCellMatrix(sc_Meta, sc_Count, markers=sc_markers, select_ct=select_ct, ct_cov=ct_cov)
    
    cm_gene = sc_basis.DM.index
    if markers is None:
        if len(cm_gene) < 0.2 * min(len(bulk_gene), len(sc_Count.index)):
            raise Exception ("Too few common genes!")
    else:
        if len(cm_gene) < 0.2 * len(set(markers)):
            raise Exception ("Too few common genes!")
    
    print ("Used {} common genes ...".format(len(cm_gene)))
    
    bulk_Count = bulk_Count.loc[cm_gene,:]
    Yjg = BulkMatrix(np.array(bulk_Count)).X
    Yjg = pd.DataFrame(Yjg, index = bulk_Count.index, columns = bulk_Count.columns)
    N_bulks = Yjg.shape[1]
    S = sc_basis.MLS
    
    #------------Deconvolution-----------
    if ct_cov:
        Results = ModelResults(initial_p = [],initial_coef = [],initial_resid = [],
                               weight = [],weight_p = [],weight_coef = [],
                               weight_resid = [],R2 = [],var_p = [])
        
        for i in range(N_bulks):
            
            ## remove zero gene in bulk sample i
            Y_i = Yjg.iloc[:,i]
            mask = Y_i != 0
            Y_i = Y_i[mask]
            D_i = sc_basis.DM[mask]
            sigma_ct_i = sc_basis.sigma_ct.loc[:,mask]
            print ('sample {} has {} available genes'.format(bulk_Count.columns[i],len(Y_i)))
            
            sample_results = music_basic_ct(Y_i, D_i, S, 
                                            sigma_ct = sigma_ct_i,
                                            iter_max = iter_max,
                                            nu = nu, eps = eps,
                                            centered = centered,
                                            normalize = normalize)
            
            Results.initial_p.append(sample_results.initial_p.tolist())
            Results.initial_coef.append(sample_results.initial_coef.tolist())
            Results.initial_resid.append(sample_results.initial_resid.reindex(cm_gene).values.tolist())
            Results.weight.append(sample_results.weight.reindex(cm_gene).values.tolist())
            Results.weight_p.append(sample_results.weight_p.tolist())
            Results.weight_coef.append(sample_results.weight_coef.tolist())
            Results.weight_resid.append(sample_results.weight_resid.reindex(cm_gene).values.tolist())
            Results.R2.append(sample_results.R2)
            Results.var_p.append(sample_results.var_p.tolist())
        
        samples = bulk_Count.columns
        cells = sc_basis.DM.columns
        genes = cm_gene
        
        Results.initial_p = pd.DataFrame(Results.initial_p, index = samples, columns = cells)
        Results.initial_coef = pd.DataFrame(Results.initial_coef, index = samples, columns = cells)
        Results.initial_resid = pd.DataFrame(Results.initial_resid,index = samples, columns = genes)
        Results.weight = pd.DataFrame(np.array(Results.weight).T, index = cm_gene, columns = samples)
        Results.weight_p = pd.DataFrame(Results.weight_p, index = samples, columns = cells)
        Results.weight_coef = pd.DataFrame(Results.weight_coef, index = samples, columns = cells)
        Results.weight_resid = pd.DataFrame(Results.weight_resid, index = samples, columns = genes)
        Results.R2 = pd.Series(Results.R2, index = samples, name = 'R2')
        Results.var_p = pd.DataFrame(Results.var_p, index = samples, columns = cells)
    
    else:
        Results = ModelResults(initial_p = [],initial_coef = [],initial_resid = [],
                               weight = [],weight_p = [],weight_coef = [],
                               weight_resid = [],R2 = [],var_p = [])
        
        for i in range(N_bulks):
            
            ## remove zero gene in bulk sample i
            Y_i = Yjg.iloc[:,i]
            mask = Y_i != 0
            Y_i = Y_i[mask]
            D_i = sc_basis.DM[mask]
            sigma_i = sc_basis.sigma[mask]
            print ('sample {} has {} available genes'.format(bulk_Count.columns[i],len(Y_i)))
            sample_results = music_basic(Y_i, D_i, S, 
                                         sigma = sigma_i,
                                         iter_max = iter_max,
                                         nu = nu, eps = eps,
                                         centered = centered,
                                         normalize = normalize)
            
            Results.initial_p.append(sample_results.initial_p.tolist())
            Results.initial_coef.append(sample_results.initial_coef.tolist())
            Results.initial_resid.append(sample_results.initial_resid.reindex(cm_gene).values.tolist())
            Results.weight.append(sample_results.weight.reindex(cm_gene).values.tolist())
            Results.weight_p.append(sample_results.weight_p.tolist())
            Results.weight_coef.append(sample_results.weight_coef.tolist())
            Results.weight_resid.append(sample_results.weight_resid.reindex(cm_gene).values.tolist())
            Results.R2.append(sample_results.R2)
            Results.var_p.append(sample_results.var_p.tolist())
        
        samples = bulk_Count.columns
        cells = sc_basis.DM.columns
        genes = cm_gene
        
        Results.initial_p = pd.DataFrame(Results.initial_p, index = samples, columns = cells)
        Results.initial_coef = pd.DataFrame(Results.initial_coef, index = samples, columns = cells)
        Results.initial_resid = pd.DataFrame(Results.initial_resid,index = samples, columns = genes)
        Results.weight = pd.DataFrame(np.array(Results.weight).T, index = cm_gene, columns = samples)
        Results.weight_p = pd.DataFrame(Results.weight_p, index = samples, columns = cells)
        Results.weight_coef = pd.DataFrame(Results.weight_coef, index = samples, columns = cells)
        Results.weight_resid = pd.DataFrame(Results.weight_resid, index = samples, columns = genes)
        Results.R2 = pd.Series(Results.R2, index = samples, name = 'R2')
        Results.var_p = pd.DataFrame(Results.var_p, index = samples, columns = cells)            
            
    return Results
        
        
        
      
            
    
    
    
                    
            
            
            
                
                
        
                        
        
        
        
        
        
        
        
        
        
        