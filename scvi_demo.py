# %%

import numpy as np
import scipy as sp
import scipy.sparse
import pandas as pd
import os
# import bbknn

# conda install -c bioconda bbknn 
# conda install scvi-tools -c bioconda -c conda-forge

import anndata as ad

import scanpy as sc
import scvi

import matplotlib.pyplot as plt

import hummingbird as hbdx
from hummingbird.io import print, rule

# %%

%matplotlib inline
%load_ext autoreload
%autoreload 2

# %%

# settings
adata_path = "/data/hbdx/test_classifynder/data/LC__ngs__rpm_log-21.3.0"
# adata_path = "LC__ngs__rpm__log-current"
n_features = 2**14

sc.set_figure_params(figsize=(10,10))
#plt.rcParams['figure.figsize'] = [15,10]
#plt.rcParams['figure.dpi'] = 100 # 200 e.g. is really fine, but slower

# %%
def get_data(adata_path, batch_col="Lab_Multiplexing_pool_ID", cell_col=None, condition_col="Diagnosis_Group"):
    adata = hbdx.io.load(adata_path)
    
#     create new column in obs, pool_ID_simple, which is represented in int instead of the pool ID name. 
    pool_id_mapping = {label:f"pool{i}" for i, label in enumerate(adata.obs.Lab_Multiplexing_pool_ID.unique())}
    adata.obs["Lab_Multiplexing_pool_ID_simple"] = adata.obs.Lab_Multiplexing_pool_ID.map(pool_id_mapping)
    
    
    #adata.obs["gender_diagnosis"] = adata.obs.apply(lambda x: x["Gender"]+"__" +x["Diagnosis_Group"], axis=1)
    
#     not sure what/why is happening here, probably to make the column name smaller
    adata.obs["batch"] = adata.obs[batch_col]

#     setting all unknown cell type to blood
    adata.obs["cell_type"] = adata.obs[cell_col] if not cell_col is None else "blood"
    
    #adata.obs["condition"] = adata.obs[condition_col]

    #adata = adata[adata.obs["Sample_Group"]=="Study"].copy()
    
    return adata

def apply_preprocessing(adata, feature_selection=False, n_top_genes=n_features):
    adata = adata.copy()
    """
    creates new columns with qc_metrics, removes low count genes, scales data
    and optionally selects highly variable genes
    """
    
#     what is the difference between adata.raw.X and adata.X?
    adata.layers["counts"] = adata.raw.X
    #adata.raw = adata.copy()
    
#     scanpy qc metrics, quality control? for what? 
# I guess it creates new columns like dropout
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    
    sc.pp.scale(adata) # Scale data to unit variance and zero mean.

    # Plot dispersions or normalized variance versus means for genes. Documentation vague
    
    # sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes) 
    hvg_selector = hbdx.pipeline.HVGSelector(n_top_genes=n_top_genes)
    

    adata = adata[adata.obs.n_genes_by_counts > 5, :] # filtering out low gene count?
    
#     I suppose sc.pp.highly_variable_genes() creates a boolean column called highly_variable in var which you can use to filter. 
# This is not documented at all
    if feature_selection:
        adata = hvg_selector.fit_transform(adata)
        # adata = adata[:, adata.var.highly_variable]

#     print(adata.var.highly_variable)
    #sc.tl.pca(adata)
    #sc.pp.neighbors(adata)
    #sc.tl.umap(adata)
    
    #adata.raw = adata
    #adata.layers["log_rpm_scaled"] = adata
    #adata = adata.raw

    return adata

    # %%
    adata = get_data(adata_path,
                #cell_col="Gender"
                )
# Filter out
baddies = ["468_0019_S__210122_468_Pool6", 
           "468_0020_S__210122_468_Pool6",
           "468_0021_S__210122_468_Pool6",
           "468_0022_S__210122_468_Pool6",
           "468_0023_S__210122_468_Pool6",
           "468_0024_S__210122_468_Pool6"]

adata = adata[~adata.obs.index.isin(baddies),:].copy()

print(adata.obs.groupby(["Sequencer","Lab_RNA_extr_protocol"]).batch.value_counts())

# %%
# print(adata.var["feature_ID"].head())

# %%
holdout_batches=["GF_Set1", "Konstanz_Pool1", "Thorax_Pool2_D2"] #these do not exist in this dataset?

# %%
# print(adata.obs.head())

# for i in adata.obs.batch.unique():
#     print(i)

# print(adata.obs.batch.unique().categories)
holdout_batches=["201123_454A_P3", "201126_454A_P6"]

if True:
    holdout_batches = np.random.choice(adata.obs.batch.unique().categories, 2, replace=False)

print("Heldout batches are: ", holdout_batches)

adata.obs["heldout"] = adata.obs["batch"].isin(holdout_batches).astype(str)

# print(adata.obs.groupby("batch").count())
# print(adata.obs["holdout"])

# print(adata.obs.head())

# %%

#split data

# only apply preprocessing to batches NOT in holdout. returns top 5000 features/HVG.
ref_adata = apply_preprocessing(adata[~adata.obs.batch.isin(holdout_batches),:], feature_selection=True, n_top_genes=n_features) 
print(ref_adata.shape)

# Held-out batches without any pre-processing done.
new_adata = adata[adata.obs.batch.isin(holdout_batches), :].copy()

print(new_adata.shape)
new_adata.layers["counts"] = new_adata.raw.X
# filter down the held out data on the HVG determined by the held in data.
# Note, this data is scaled, only within itself!
sc.pp.scale(new_adata)
new_adata = new_adata[:, ref_adata.var_names].copy()
print(new_adata.shape)

# %%

n_latent_vae = 2**7
# Filter out samples with less than count of 1, reduces samples from 1679 to 661!
sc.pp.filter_cells(ref_adata, 
                min_counts=1
                    )
print(ref_adata.shape)

scvi.data.setup_anndata(ref_adata, batch_key="batch", layer="counts", categorical_covariate_keys=["Sequencer", "Lab_RNA_extr_protocol", "Sample_Group", "Gender", "Diagnosis_Group"])

vae = scvi.model.SCVI(ref_adata, n_latent=n_latent_vae, n_layers=1)

# %%
print(vae.get_reconstruction_error(ref_adata))
vae.train(max_epochs = 1200, check_val_every_n_epoch=10, early_stopping=True)
print(vae.get_reconstruction_error(ref_adata))

# %%

# print(vae.__dict__)

# print(vae.__dict__.keys())

# %%

#vae.save("scvi_model/")
#vae = scvi.model.SCVI.load("scvi_model/", adata, use_gpu=True)


# scvi.data.transfer_anndata_setup(new_adata, ref_adata, extend_categories=True)
# scvi.data.setup_anndata(new_adata, batch_key="batch", layer="counts", categorical_covariate_keys=["Sequencer", "Lab_RNA_extr_protocol", "Sample_Group", "Gender", "Diagnosis_Group"])

# _ = scvi.model.SCVI(new_adata, n_latent=50)
# _.train()
# print(latent.shape)
# %%

latent = vae.get_latent_representation()
# latent_new = vae.get_latent_representation(new_adata)
# print(latent)

# print(ref_adata.obsm)

# %%
# # OKAY, so if you want to do inference, you need to re-initialize the VAE with the new data. 

scvi.data.transfer_anndata_setup(ref_adata, new_adata, extend_categories=True)
scvi.data.setup_anndata(new_adata, batch_key="batch", layer="counts", categorical_covariate_keys=["Sequencer", "Lab_RNA_extr_protocol", "Sample_Group", "Gender", "Diagnosis_Group"])

vae_new_data = scvi.model.SCVI.load_query_data(new_adata, vae)

vae_new_data.train(max_epochs=0, plan_kwargs=dict(weight_decay=0.0))
# %%
latent_predicted = vae_new_data.get_latent_representation()
latent_control = vae_new_data.get_latent_representation(ref_adata)
# %%
print(latent.shape, latent_predicted.shape, latent_control.shape)

# %%
sc.tl.pca(adata)
sc.pp.neighbors(adata, use_rep="X_pca")
sc.tl.umap(adata)

ref_adata.obsm["X_scVI"] = latent
new_adata.obsm["X_scVI"] = latent_predicted

sc.pp.neighbors(ref_adata, use_rep="X_scVI")
sc.tl.umap(ref_adata)

sc.pp.neighbors(new_adata, use_rep="X_scVI")
sc.tl.umap(new_adata)


# %%

import matplotlib.pyplot as plt

train_elbo = vae.history['elbo_train'][1:]
test_elbo = vae.history['elbo_validation']

ax = train_elbo.plot()
test_elbo.plot(ax = ax)

train_recon = vae.history['reconstruction_loss_train'][1:]
test_recon = vae.history['reconstruction_loss_validation']

ax = train_recon.plot()
test_recon.plot(ax = ax)

# %%
concated = ad.concat([ref_adata, new_adata])
new_adata_mapped = sc.tl.ingest(ref_adata, new_adata, inplace=False, embedding_method='umap')




# %%

colors = ["heldout", "Lab_Multiplexing_pool_ID", 
                                        "Diagnosis_Group", 
                                        "Gender",
                                        "Sequencer",
                                        "Lab_RNA_extr_protocol", 
                                        "Lab_Blocking_Protocol",
                                        "Site_Measurement", 
                                        "TumStageSimp", 
                                        "Diagnosis_Group_sub"]

# Anndata raw, no VAE or scaling or preprocessing or HVG selection

sc.pl.umap(adata, color=colors, ncols=2, use_raw=False, )

# ref_adata VAE latent

sc.pl.umap(ref_adata, color=colors, ncols=2, use_raw=False, )

# "ref_adata + new data inferred from VAE latent"
sc.pl.umap(concated, color=colors, ncols=2, use_raw=False, )


# %%

# %%

# %%

# vae.get_reconstruction_error(new_adata)
# %%
