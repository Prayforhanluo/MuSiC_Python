##
README

####
MuSiC in python.

Reproduce the deconvolution function of MuSiC.

More information about MuSiC here [xuranw's github](https://github.com/xuranw/MuSiC)

Fix some bug in original MuSiC when calculate the weight by co-variance of cross-subject.

Data formats as txt in data dirct

Use the function like :
	
	 sc_Meta = Get_SC_Meta_Info('sc_meta.txt')
	 sc_Count = Get_SC_Count_Info('sc_count.txt')
	 bulk_Meta = Get_bulk_Meta_Info('bulk_meta.txt')
	 bulk_Count = Get_bulk_Count_Info('bulk_Count.txt')
	 
	 select_ct = ['alpha', 'beta', 'delta', 'gamma', 'acinar', 'ductal']
	 
	 #weight gene by variance of cross-subject
	 Results =  music_prop(bulk_Meta, bulk_Count, sc_Meta, sc_Count,ct_cov=False, select_ct=select_ct)
	 
     #weight gene by co-variance of cross-subject
	 Results =  music_prop(bulk_Meta, bulk_Count, sc_Meta, sc_Count,ct_cov=True, select_ct=select_ct)


results :
	
	Results.initial_p    #NNLS results
	Results.weight_p     #Weighted-NNLS results
	Results.R2			 #R-squared for weighted-NNLS
	...
	...
