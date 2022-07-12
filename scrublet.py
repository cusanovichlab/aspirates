conda activate python3.7_base
python3.7

import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os

input_dir="/xdisk/darrenc/hwelfley/redo/91bp/without_testrun/"

thirteen_counts_matrix = scipy.io.mmread(input_dir + 'NEA_13/outs/filtered_feature_bc_matrix/matrix.mtx').T.tocsc()
thirteen_genes = np.array(scr.load_genes(input_dir + 'NEA_13/outs/filtered_feature_bc_matrix/features.tsv', delimiter='\t', column=1))

fifteen_counts_matrix = scipy.io.mmread(input_dir + 'NEA_15/outs/filtered_feature_bc_matrix/matrix.mtx').T.tocsc()
fifteen_genes = np.array(scr.load_genes(input_dir + 'NEA_15/outs/filtered_feature_bc_matrix/features.tsv', delimiter='\t', column=1))

twenty_counts_matrix = scipy.io.mmread(input_dir + 'NEA_20/outs/filtered_feature_bc_matrix/matrix.mtx').T.tocsc()
twenty_genes = np.array(scr.load_genes(input_dir + 'NEA_20/outs/filtered_feature_bc_matrix/features.tsv', delimiter='\t', column=1))

twentytwo_counts_matrix = scipy.io.mmread(input_dir + 'NEA_22/outs/filtered_feature_bc_matrix/matrix.mtx').T.tocsc()
twentytwo_genes = np.array(scr.load_genes(input_dir + 'NEA_22/outs/filtered_feature_bc_matrix/features.tsv', delimiter='\t', column=1))

twentythree_counts_matrix = scipy.io.mmread(input_dir + 'NEA_23/outs/filtered_feature_bc_matrix/matrix.mtx').T.tocsc()
twentythree_genes = np.array(scr.load_genes(input_dir + 'NEA_23/outs/filtered_feature_bc_matrix/features.tsv', delimiter='\t', column=1))


twentysix_counts_matrix = scipy.io.mmread(input_dir + 'NEA_26/outs/filtered_feature_bc_matrix/matrix.mtx').T.tocsc()
twentysix_genes = np.array(scr.load_genes(input_dir + 'NEA_26/outs/filtered_feature_bc_matrix/features.tsv', delimiter='\t', column=1))

twentyeight_counts_matrix = scipy.io.mmread(input_dir + 'NEA_28/outs/filtered_feature_bc_matrix/matrix.mtx').T.tocsc()
twentyeight_genes = np.array(scr.load_genes(input_dir + 'NEA_28/outs/filtered_feature_bc_matrix/features.tsv', delimiter='\t', column=1))

twentynine_counts_matrix = scipy.io.mmread(input_dir + 'NEA_29/outs/filtered_feature_bc_matrix/matrix.mtx').T.tocsc()
twentynine_genes = np.array(scr.load_genes(input_dir + 'NEA_29/outs/filtered_feature_bc_matrix/features.tsv', delimiter='\t', column=1))

thirtyone_counts_matrix = scipy.io.mmread(input_dir + 'NEA_31/outs/filtered_feature_bc_matrix/matrix.mtx').T.tocsc()
thirtyone_genes = np.array(scr.load_genes(input_dir + 'NEA_31/outs/filtered_feature_bc_matrix/features.tsv', delimiter='\t', column=1))

thirtythree_counts_matrix = scipy.io.mmread(input_dir + 'NEA_33/outs/filtered_feature_bc_matrix/matrix.mtx').T.tocsc()
thirtythree_genes = np.array(scr.load_genes(input_dir + 'NEA_33/outs/filtered_feature_bc_matrix/features.tsv', delimiter='\t', column=1))


###

thirteen_scrub = scr.Scrublet(thirteen_counts_matrix, expected_doublet_rate=0.06)
fifteen_scrub = scr.Scrublet(fifteen_counts_matrix, expected_doublet_rate=0.06)
twenty_scrub = scr.Scrublet(twenty_counts_matrix, expected_doublet_rate=0.06)
twentytwo_scrub = scr.Scrublet(twentytwo_counts_matrix, expected_doublet_rate=0.06)
twentythree_scrub = scr.Scrublet(twentythree_counts_matrix, expected_doublet_rate=0.06)
twentysix_scrub = scr.Scrublet(twentysix_counts_matrix, expected_doublet_rate=0.06)
twentyeight_scrub = scr.Scrublet(twentyeight_counts_matrix, expected_doublet_rate=0.06)
twentynine_scrub = scr.Scrublet(twentynine_counts_matrix, expected_doublet_rate=0.06)
thirtyone_scrub = scr.Scrublet(thirtyone_counts_matrix, expected_doublet_rate=0.06)
thirtythree_scrub = scr.Scrublet(thirtythree_counts_matrix, expected_doublet_rate=0.06)


thirteen_doublet_scores, thirteen_predicted_doublets = thirteen_scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
#Detected=0.2%
#Estimated=26.1%

fifteen_doublet_scores, fifteen_predicted_doublets = fifteen_scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)  
#Detected=0.4%
#Estimated=2.6%                                                    

twenty_doublet_scores, twenty_predicted_doublets = twenty_scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
#Detected=0.4%
#Estimated=2.0%

twentytwo_doublet_scores, twentytwo_predicted_doublets = twentytwo_scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
#Detected=1.0%
#Estimated=5.2%

twentythree_doublet_scores, twentythree_predicted_doublets = twentythree_scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
#Detected=0.1%
#Estimated=6.6%


twentysix_doublet_scores, twentysix_predicted_doublets = twentysix_scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
#Detected=0.2%
#Estimated=1.9%


twentyeight_doublet_scores, twentyeight_predicted_doublets = twentyeight_scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
#Detected=0.2%
#Estimated=6.7%

twentynine_doublet_scores, twentynine_predicted_doublets = twentynine_scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
#Detected=0.0%
#Estimated=0.0%

thirtyone_doublet_scores, thirtyone_predicted_doublets = thirtyone_scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

#Detected=0.7%
#Estimated=7.1%

thirtythree_doublet_scores, thirtythree_predicted_doublets = thirtythree_scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
#Detected=0.3%
#Estimated=9.1%

thirteen_scrub.plot_histogram();
plt.savefig("/xdisk/darrenc/hwelfley/91bp/without_testrun/scrubplots/thirteen_scrub.pdf")
plt.show()

fifteen_scrub.plot_histogram();
plt.savefig("/xdisk/darrenc/hwelfley/91bp/without_testrun/scrubplots/fifteen_scrub.pdf")
plt.show()

twenty_scrub.plot_histogram();
plt.savefig("/xdisk/darrenc/hwelfley/91bp/without_testrun/scrubplots/twenty_scrub.pdf")
plt.show()

twentytwo_scrub.plot_histogram();
plt.savefig("/xdisk/darrenc/hwelfley/91bp/without_testrun/scrubplots/twentytwo_scrub.pdf")
plt.show()

twentythree_scrub.plot_histogram(); 
plt.savefig("/xdisk/darrenc/hwelfley/91bp/without_testrun/scrubplots/twentythree_scrub.pdf")                 
plt.show()

twentysix_scrub.plot_histogram();
plt.savefig("/xdisk/darrenc/hwelfley/91bp/without_testrun/scrubplots/twentysix_scrub.pdf")
plt.show()

twentyeight_scrub.plot_histogram();
plt.savefig("/xdisk/darrenc/hwelfley/91bp/without_testrun/scrubplots/twentyeight_scrub.pdf")
plt.show()

twentynine_scrub.plot_histogram();
plt.savefig("/xdisk/darrenc/hwelfley/91bp/without_testrun/scrubplots/twentynine_scrub.pdf")
plt.show()

thirtyone_scrub.plot_histogram();
plt.savefig("/xdisk/darrenc/hwelfley/91bp/without_testrun/scrubplots/thirtyone_scrub.pdf")
plt.show()

thirtythree_scrub.plot_histogram();
plt.savefig("/xdisk/darrenc/hwelfley/91bp/without_testrun/scrubplots/thirtythree_scrub.pdf")
plt.show()




#####open the barcodes for each matrix to make list of barcodes to rm

with open(os.path.join(input_dir, 'NEA_13/outs/filtered_feature_bc_matrix/barcodes.tsv')) as f:
        thirteen_barcodes = [barcode.strip() if not barcode.strip().endswith('-1') else barcode.strip()[:-2] for barcode in f.readlines()]

with open(os.path.join(input_dir, 'NEA_15/outs/filtered_feature_bc_matrix/barcodes.tsv')) as f:
        fifteen_barcodes = [barcode.strip() if not barcode.strip().endswith('-1') else barcode.strip()[:-2] for barcode in f.readlines()]        

with open(os.path.join(input_dir, 'NEA_20/outs/filtered_feature_bc_matrix/barcodes.tsv')) as f:
        twenty_barcodes = [barcode.strip() if not barcode.strip().endswith('-1') else barcode.strip()[:-2] for barcode in f.readlines()]

with open(os.path.join(input_dir, 'NEA_22/outs/filtered_feature_bc_matrix/barcodes.tsv')) as f:
        twentytwo_barcodes = [barcode.strip() if not barcode.strip().endswith('-1') else barcode.strip()[:-2] for barcode in f.readlines()]

with open(os.path.join(input_dir, 'NEA_23/outs/filtered_feature_bc_matrix/barcodes.tsv')) as f:
        twentythree_barcodes = [barcode.strip() if not barcode.strip().endswith('-1') else barcode.strip()[:-2] for barcode in f.readlines()]


with open(os.path.join(input_dir, 'NEA_26/outs/filtered_feature_bc_matrix/barcodes.tsv')) as f:
        twentysix_barcodes = [barcode.strip() if not barcode.strip().endswith('-1') else barcode.strip()[:-2] for barcode in f.readlines()]

with open(os.path.join(input_dir, 'NEA_28/outs/filtered_feature_bc_matrix/barcodes.tsv')) as f:
        twentyeight_barcodes = [barcode.strip() if not barcode.strip().endswith('-1') else barcode.strip()[:-2] for barcode in f.readlines()]

with open(os.path.join(input_dir, 'NEA_29/outs/filtered_feature_bc_matrix/barcodes.tsv')) as f:
        twentynine_barcodes = [barcode.strip() if not barcode.strip().endswith('-1') else barcode.strip()[:-2] for barcode in f.readlines()]

with open(os.path.join(input_dir, 'NEA_31/outs/filtered_feature_bc_matrix/barcodes.tsv')) as f:
        thirtyone_barcodes = [barcode.strip() if not barcode.strip().endswith('-1') else barcode.strip()[:-2] for barcode in f.readlines()]     


with open(os.path.join(input_dir, 'NEA_33/outs/filtered_feature_bc_matrix/barcodes.tsv')) as f:
        thirtythree_barcodes = [barcode.strip() if not barcode.strip().endswith('-1') else barcode.strip()[:-2] for barcode in f.readlines()]            


# Save singlet, double csv files.

thirteen_doublet_scores = thirteen_scrub.doublet_scores_obs_
thirteen_predicted_doublets = thirteen_scrub.predicted_doublets_

fifteen_doublet_scores = fifteen_scrub.doublet_scores_obs_
fifteen_predicted_doublets = fifteen_scrub.predicted_doublets_

twenty_doublet_scores = twenty_scrub.doublet_scores_obs_
twenty_predicted_doublets = twenty_scrub.predicted_doublets_

twentytwo_doublet_scores = twentytwo_scrub.doublet_scores_obs_
twentytwo_predicted_doublets = twentytwo_scrub.predicted_doublets_

twentythree_doublet_scores = twentythree_scrub.doublet_scores_obs_
twentythree_predicted_doublets = twentythree_scrub.predicted_doublets_

twentysix_doublet_scores = twentysix_scrub.doublet_scores_obs_
twentysix_predicted_doublets = twentysix_scrub.predicted_doublets_

twentyeight_doublet_scores = twentyeight_scrub.doublet_scores_obs_
twentyeight_predicted_doublets = twentyeight_scrub.predicted_doublets_

twentynine_doublet_scores = twentynine_scrub.doublet_scores_obs_
twentynine_predicted_doublets = twentynine_scrub.predicted_doublets_

thirtyone_doublet_scores = thirtyone_scrub.doublet_scores_obs_
thirtyone_predicted_doublets = thirtyone_scrub.predicted_doublets_

thirtythree_doublet_scores = thirtythree_scrub.doublet_scores_obs_
thirtythree_predicted_doublets = thirtythree_scrub.predicted_doublets_


#singlet and doublet


#thirteen

thirteen_singlet_barcodes = [thirteen_barcode for i,thirteen_barcode in enumerate(thirteen_barcodes) if thirteen_predicted_doublets[i]==False]
f = open(os.path.join(input_dir, 'NEA_13/outs/thirteen_singlet_barcodes.csv'), "w")
for thirteen_barcode in thirteen_singlet_barcodes:
        f.write('{}\n'.format(thirteen_barcode))

f.close()

thirteen_doublet_barcodes=[thirteen_barcode for i,thirteen_barcode in enumerate(thirteen_barcodes) if thirteen_predicted_doublets[i]==True]
f = open(os.path.join(input_dir, 'NEA_13/outs/thirteen_doublet_barcodes.csv'), "w")
for thirteen_barcode in thirteen_doublet_barcodes:
        f.write('{}\n'.format(thirteen_barcode))


f.close()

#fifteen

fifteen_singlet_barcodes = [fifteen_barcode for i,fifteen_barcode in enumerate(fifteen_barcodes) if fifteen_predicted_doublets[i]==False]
f = open(os.path.join(input_dir, 'NEA_15/outs/fifteen_singlet_barcodes.csv'), "w")
for fifteen_barcode in fifteen_singlet_barcodes:
        f.write('{}\n'.format(fifteen_barcode))


f.close()

fifteen_doublet_barcodes=[fifteen_barcode for i,fifteen_barcode in enumerate(fifteen_barcodes) if fifteen_predicted_doublets[i]==True]
f = open(os.path.join(input_dir, 'NEA_15/outs/fifteen_doublet_barcodes.csv'), "w")
for fifteen_barcode in fifteen_doublet_barcodes:
        f.write('{}\n'.format(fifteen_barcode))


f.close()


#twenty

twenty_singlet_barcodes = [twenty_barcode for i,twenty_barcode in enumerate(twenty_barcodes) if twenty_predicted_doublets[i]==False]
f = open(os.path.join(input_dir, 'NEA_20/outs/twenty_singlet_barcodes.csv'), "w")
for twenty_barcode in twenty_singlet_barcodes:
        f.write('{}\n'.format(twenty_barcode))


f.close()

twenty_doublet_barcodes=[twenty_barcode for i,twenty_barcode in enumerate(twenty_barcodes) if twenty_predicted_doublets[i]==True]
f = open(os.path.join(input_dir, 'NEA_20/outs/twenty_doublet_barcodes.csv'), "w")
for twenty_barcode in twenty_doublet_barcodes:
        f.write('{}\n'.format(twenty_barcode))


f.close()

#22

twentytwo_singlet_barcodes = [twentytwo_barcode for i,twentytwo_barcode in enumerate(twentytwo_barcodes) if twentytwo_predicted_doublets[i]==False]
f = open(os.path.join(input_dir, 'NEA_22/outs/twentytwo_singlet_barcodes.csv'), "w")
for twentytwo_barcode in twentytwo_singlet_barcodes:
        f.write('{}\n'.format(twentytwo_barcode))


f.close()

twentytwo_doublet_barcodes=[twentytwo_barcode for i,twentytwo_barcode in enumerate(twentytwo_barcodes) if twentytwo_predicted_doublets[i]==True]
f = open(os.path.join(input_dir, 'NEA_22/outs/twentytwo_doublet_barcodes.csv'), "w")
for twentytwo_barcode in twentytwo_doublet_barcodes:
        f.write('{}\n'.format(twentytwo_barcode))


f.close()

#23

twentythree_singlet_barcodes = [twentythree_barcode for i,twentythree_barcode in enumerate(twentythree_barcodes) if twentythree_predicted_doublets[i]==False]
f = open(os.path.join(input_dir, 'NEA_23/outs/twentythree_singlet_barcodes.csv'), "w")
for twentythree_barcode in twentythree_singlet_barcodes:
        f.write('{}\n'.format(twentythree_barcode))


f.close()

twentythree_doublet_barcodes=[twentythree_barcode for i,twentythree_barcode in enumerate(twentythree_barcodes) if twentythree_predicted_doublets[i]==True]
f = open(os.path.join(input_dir, 'NEA_23/outs/twentythree_doublet_barcodes.csv'), "w")
for twentythree_barcode in twentythree_doublet_barcodes:
        f.write('{}\n'.format(twentythree_barcode))


f.close()



twentysix_singlet_barcodes = [twentysix_barcode for i,twentysix_barcode in enumerate(twentysix_barcodes) if twentysix_predicted_doublets[i]==False]
f = open(os.path.join(input_dir, 'NEA_26/outs/twentysix_singlet_barcodes.csv'), "w")
for twentysix_barcode in twentysix_singlet_barcodes:
        f.write('{}\n'.format(twentysix_barcode))


f.close()

twentysix_doublet_barcodes=[twentysix_barcode for i,twentysix_barcode in enumerate(twentysix_barcodes) if twentysix_predicted_doublets[i]==True]
f = open(os.path.join(input_dir, 'NEA_26/outs/twentysix_doublet_barcodes.csv'), "w")
for twentysix_barcode in twentysix_doublet_barcodes:
        f.write('{}\n'.format(twentysix_barcode))


f.close()

twentyeight_singlet_barcodes = [twentyeight_barcode for i,twentyeight_barcode in enumerate(twentyeight_barcodes) if twentyeight_predicted_doublets[i]==False]
f = open(os.path.join(input_dir, 'NEA_28/outs/twentyeight_singlet_barcodes.csv'), "w")
for twentyeight_barcode in twentyeight_singlet_barcodes:
        f.write('{}\n'.format(twentyeight_barcode))


f.close()

twentyeight_doublet_barcodes=[twentyeight_barcode for i,twentyeight_barcode in enumerate(twentyeight_barcodes) if twentyeight_predicted_doublets[i]==True]
f = open(os.path.join(input_dir, 'NEA_28/outs/twentyeight_doublet_barcodes.csv'), "w")
for twentyeight_barcode in twentyeight_doublet_barcodes:
        f.write('{}\n'.format(twentyeight_barcode))


f.close()

twentynine_singlet_barcodes = [twentynine_barcode for i,twentynine_barcode in enumerate(twentynine_barcodes) if twentynine_predicted_doublets[i]==False]
f = open(os.path.join(input_dir, 'NEA_29/outs/twentynine_singlet_barcodes.csv'), "w")
for twentynine_barcode in twentynine_singlet_barcodes:
        f.write('{}\n'.format(twentynine_barcode))


f.close()

twentynine_doublet_barcodes=[twentynine_barcode for i,twentynine_barcode in enumerate(twentynine_barcodes) if twentynine_predicted_doublets[i]==True]
f = open(os.path.join(input_dir, 'NEA_29/outs/twentynine_doublet_barcodes.csv'), "w")
for twentynine_barcode in twentynine_doublet_barcodes:
        f.write('{}\n'.format(twentynine_barcode))


f.close()

thirtyone_singlet_barcodes = [thirtyone_barcode for i,thirtyone_barcode in enumerate(thirtyone_barcodes) if thirtyone_predicted_doublets[i]==False]
f = open(os.path.join(input_dir, 'NEA_31/outs/thirtyone_singlet_barcodes.csv'), "w")
for thirtyone_barcode in thirtyone_singlet_barcodes:
        f.write('{}\n'.format(thirtyone_barcode))


f.close()

thirtyone_doublet_barcodes=[thirtyone_barcode for i,thirtyone_barcode in enumerate(thirtyone_barcodes) if thirtyone_predicted_doublets[i]==True]
f = open(os.path.join(input_dir, 'NEA_31/outs/thirtyone_doublet_barcodes.csv'), "w")
for thirtyone_barcode in thirtyone_doublet_barcodes:
        f.write('{}\n'.format(thirtyone_barcode))


f.close()

thirtythree_singlet_barcodes = [thirtythree_barcode for i,thirtythree_barcode in enumerate(thirtythree_barcodes) if thirtythree_predicted_doublets[i]==False]
f = open(os.path.join(input_dir, 'NEA_33/outs/thirtythree_singlet_barcodes.csv'), "w")
for thirtythree_barcode in thirtythree_singlet_barcodes:
        f.write('{}\n'.format(thirtythree_barcode))
f.close()


thirtythree_doublet_barcodes=[thirtythree_barcode for i,thirtythree_barcode in enumerate(thirtythree_barcodes) if thirtythree_predicted_doublets[i]==True]
f = open(os.path.join(input_dir, 'NEA_33/outs/thirtythree_doublet_barcodes.csv'), "w")
for thirtythree_barcode in thirtythree_doublet_barcodes:
        f.write('{}\n'.format(thirtythree_barcode))


f.close()






#DONE!



