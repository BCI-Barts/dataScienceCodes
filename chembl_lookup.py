from chembl_webresource_client.new_client import new_client
import pandas as pd
#-------------------------------------
# Note: this is organism independent
#-------------------------------------
# there is a loose search on human only results
#-------------------------------------
#genes=pd.read_csv("gene_names.txt")
#---------------------------------------
genes=['PDE4D']
results=pd.DataFrame(columns=['gene','gene_name','target_id',
                              'molA_id','molA_type','molA_name',
                              'molB_id','molB_type','molB_name','similarity'])
#-------------------------------------
for gene in genes:
    #--- gene cycle
    target = new_client.target
    res = target.filter(target_synonym__icontains=gene, target_organism__istartswith='Human') 
    print(res)
    for items in res:
        #--- item cycle through results-----
        # 
        #-----------------------------------
        activities = new_client.activity.filter(target_chembl_id__in=items['target_chembl_id']).only(['molecule_chembl_id'])
        #--- extracting target ChEMBL IDs 
        # ---------------------------------------------
        #
        #----------------------------------------------
        for mol in activities:
            #--- Get all molecules at 85% similarity for inferred links
            #
            #
            #----------------------------------------------------------
            molecule     = new_client.molecule
            molAA        = molecule.get(mol['molecule_chembl_id'])
            similarity   = new_client.similarity
           
            molsim = []
            try:
                molsim     = similarity.filter(chembl_id=mol['molecule_chembl_id'], similarity=85).only(['molecule_chembl_id','molecule_type','similarity'])
                molecule = new_client.molecule
            except:
                pass

            for molB in molsim:
                #--- Get molecule info attributes
                #
                #------------------------------------
                molBB     = molecule.get(molB['molecule_chembl_id'])
                
                row=        [gene, 
                            items['pref_name'],
                            items['target_chembl_id'], 
                            molAA['molecule_chembl_id'],
                            molAA['molecule_type'],
                            molAA['pref_name'],
                            molB['molecule_chembl_id'],
                            molB['molecule_type'],
                            molBB['pref_name'],
                            molB['similarity']]

                print('\t'.join([str(x) for x in row]))
                print('\n')
#----------------------------------------------