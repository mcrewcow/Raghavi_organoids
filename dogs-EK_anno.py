#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
get_ipython().system('{sys.executable} -m pip -q install palantir fa2')


# In[2]:


import warnings
warnings.filterwarnings("ignore")
from anndata import AnnData
import numpy as np
import pandas as pd
import scanpy as sc
import scFates as scf
import palantir
import matplotlib.pyplot as plt
sc.settings.verbosity = 3
sc.settings.logfile = sys.stdout
## fix palantir breaking down some plots
import seaborn
seaborn.reset_orig()
get_ipython().run_line_magic('matplotlib', 'inline')

sc.set_figure_params()
scf.set_figure_pubready()


# In[4]:


mg = sc.read_h5ad(r"/mnt/c/Users/Emil/10X/dogs/dogs_integrated.h5ad") #1500


# In[3]:


mgtotal = sc.read_h5ad(r"/mnt/c/Users/Emil/10X/dogs/mgtotal_EK_anno.h5ad")


# In[244]:


mgtotal.obs


# In[7]:


sc.pp.pca(mg)
pca_projections = pd.DataFrame(mg.obsm["X_pca"],index=mg.obs_names)


# In[8]:


dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=4)


# In[9]:


# generate neighbor draph in multiscale diffusion space
mg.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(mg,n_neighbors=30,use_rep="X_palantir")


# In[10]:


# draw ForceAtlas2 embedding using 3 first PCs as initial positions
mg.obsm["X_pca2d"]=mg.obsm["X_pca"][:,:3]
sc.tl.draw_graph(mg, layout = 'fa',init_pos='X_pca2d')


# In[11]:


sc.pp.highly_variable_genes(mg, n_top_genes=1500, flavor='cell_ranger')


# In[12]:


sc.pp.pca(mg)
pca_projections = pd.DataFrame(mg.obsm["X_pca"],index=mg.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=4)
mg.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(mg,n_neighbors=30,use_rep="X_palantir")
mg.obsm["X_pca2d"]=mg.obsm["X_pca"][:,:3]
sc.tl.draw_graph(mg, layout = 'fa',init_pos='X_pca2d')
sc.pl.draw_graph(mg,color=['Cell_ident'], layout = 'fa') #WE USE THAT


# In[246]:


sc.pl.draw_graph(mgtotal,color=['leiden'], layout = 'fa') #WE USE THAT


# In[26]:


sc.pp.pca(mg)
pca_projections = pd.DataFrame(mg.obsm["X_pca"],index=mg.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=10)
mg.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(mg,n_neighbors=30,use_rep="X_palantir")
mg.obsm["X_pca2d"]=mg.obsm["X_pca"][:,:3]
sc.tl.draw_graph(mg, layout = 'fa',init_pos='X_pca2d')
sc.pl.draw_graph(mg,color=['Cell_ident'], layout = 'fa')


# In[19]:


mg.obsp


# In[10]:


sc.pp.pca(mg)
pca_projections = pd.DataFrame(mg.obsm["X_pca"],index=mg.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=8)
mg.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(mg,n_neighbors=50,use_rep="X_palantir")
mg.obsm["X_pca2d"]=mg.obsm["X_pca"][:,:3]
sc.tl.draw_graph(mg, layout = 'fa',init_pos='X_pca2d')
sc.pl.draw_graph(mg,color=['Cell_ident'], layout = 'fa')


# In[9]:


sc.set_figure_params()
sc.pl.draw_graph(mg,color=['Cell_ident'], layout = 'fa')


# In[27]:


sc.tl.umap(mg)


# In[28]:


sc.pl.umap(mg, color = ['Cell_ident'])


# In[5]:


mg_cones = mgtotal[mgtotal.obs['EK_anno'].isin(['Cones'])]


# In[7]:


mg_rods = mgtotal[mgtotal.obs['EK_anno'].isin(['Rods'])]


# In[8]:


sc.pp.pca(mg_cones)
pca_projections = pd.DataFrame(mg_cones.obsm["X_pca"],index=mg_cones.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=10)
mg_cones.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(mg_cones,n_neighbors=30,use_rep="X_palantir")
mg_cones.obsm["X_pca2d"]=mg_cones.obsm["X_pca"][:,:3]
sc.tl.draw_graph(mg_cones, layout = 'fa',init_pos='X_pca2d')
sc.pl.draw_graph(mg_cones,color=['Cell_ident'], layout = 'fa')


# In[9]:


sc.pl.draw_graph(mg_cones,color=['EK_anno'], layout = 'fa')


# In[11]:


sc.pl.draw_graph(mg_cones,color=['treatment'], layout = 'fa')


# In[12]:


sc.tl.umap(mg_cones)


# In[13]:


sc.pl.umap(mg_cones, color = ['Cell_ident'])


# In[14]:


sc.pp.pca(mg_rods)
pca_projections = pd.DataFrame(mg_rods.obsm["X_pca"],index=mg_rods.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=10)
mg_rods.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(mg_rods,n_neighbors=30,use_rep="X_palantir")
mg_rods.obsm["X_pca2d"]=mg_rods.obsm["X_pca"][:,:3]
sc.tl.draw_graph(mg_rods, layout = 'fa',init_pos='X_pca2d')
sc.pl.draw_graph(mg_rods,color=['Cell_ident'], layout = 'fa')


# In[15]:


sc.pl.draw_graph(mg_rods,color=['EK_anno'], layout = 'fa')


# In[16]:


sc.pl.draw_graph(mg_rods,color=['treatment'], layout = 'fa')


# In[17]:


sc.tl.umap(mg_rods)


# In[18]:


sc.pl.umap(mg_rods, color = ['Cell_ident'])


# In[19]:


sc.pl.draw_graph(mgtotal,color=['Cell_ident','leiden','EK_anno'], layout = 'fa')


# In[24]:


sc.pl.draw_graph(mgtotal,color=['EK_anno'], layout = 'fa')


# In[20]:


sc.tl.rank_genes_groups(mgtotal, 'EK_anno', method='logreg')
sc.pl.rank_genes_groups(mgtotal, n_genes=25, sharey=False)


# In[21]:


def cluster_small_multiples(
    adata, clust_key, size=60, frameon=False, legend_loc=None, **kwargs
):
    tmp = adata.copy()

    for i, clust in enumerate(adata.obs[clust_key].cat.categories):
        tmp.obs[clust] = adata.obs[clust_key].isin([clust]).astype("category")
        tmp.uns[clust + "_colors"] = ["#d3d3d3", adata.uns[clust_key + "_colors"][i]]

    sc.pl.draw_graph(
        tmp,
        groups=tmp.obs[clust].cat.categories[1:].values,
        color=adata.obs[clust_key].cat.categories.tolist(),
        size=size,
        frameon=frameon,
        legend_loc=legend_loc,
        **kwargs,
    )


# In[22]:


def cluster_small_multiples_umap(
    adata, clust_key, size=60, frameon=False, legend_loc=None, **kwargs
):
    tmp = adata.copy()

    for i, clust in enumerate(adata.obs[clust_key].cat.categories):
        tmp.obs[clust] = adata.obs[clust_key].isin([clust]).astype("category")
        tmp.uns[clust + "_colors"] = ["#d3d3d3", adata.uns[clust_key + "_colors"][i]]

    sc.pl.umap(
        tmp,
        groups=tmp.obs[clust].cat.categories[1:].values,
        color=adata.obs[clust_key].cat.categories.tolist(),
        size=size,
        frameon=frameon,
        legend_loc=legend_loc,
        **kwargs,
    )


# In[25]:


cluster_small_multiples(mgtotal, clust_key = 'EK_anno', layout = 'fa')


# In[26]:


cluster_small_multiples_umap(mgtotal, clust_key = 'EK_anno')


# In[302]:


sc.tl.umap(mgtotal)


# In[101]:


cluster_small_multiples(mgtotal, clust_key = 'Cell_ident', layout = 'fa')


# In[82]:


sc.tl.rank_genes_groups(mgtotal, 'leiden', method='logreg', use_raw = False)
sc.pl.rank_genes_groups(mgtotal, n_genes=25, sharey=False)


# In[81]:


sc.pl.rank_genes_groups(mgtotal, n_genes=25, sharey=False, use_raw = False)


# In[84]:


mgtotal.obs


# In[27]:


import scvelo as scv


# In[28]:


mgtotal_I = mgtotal[mgtotal.obs['treatment'].isin(['Retina_Injected'])]
mgtotal_UI = mgtotal[mgtotal.obs['treatment'].isin(['Retina_Uninjected'])]


# In[29]:


mgtotal.obs.treatment.value_counts()


# In[30]:


I = scv.read("/mnt/c/Users/Emil/10X/U24CPI_counts_human.loom", cache = True)
UI = scv.read("/mnt/c/Users/Emil/10X/U24CUI_counts_human.loom", cache = True)


# In[ ]:


loom


# In[58]:


mg_rods


# In[51]:


mg_cones_I = scv.utils.merge(mg_cones, I)
mg_cones_UI = scv.utils.merge(mg_cones, UI)


# In[52]:


mg_cones_UI


# In[53]:


mg_cones_velo1 = mg_cones_I.concatenate([mg_cones_UI])


# In[54]:


mg_rods_velo_I = scv.utils.merge(mg_rods, I) #smth wrong with rods upon merging, should be separated


# In[55]:


mg_rods_velo_UI = scv.utils.merge(mg_rods, UI)


# In[56]:


mg_rods_velo1 = mg_rods_velo_I.concatenate([mg_rods_velo_UI])


# In[37]:


mg_rods_velo1


# In[66]:


mgtotal_I = scv.utils.merge(mgtotal_I, I)
mgtotal_UI = scv.utils.merge(mgtotal_UI, UI)


# In[67]:


scv.pl.proportions(mgtotal_I)


# In[68]:


scv.pl.proportions(mgtotal_UI)


# In[69]:


scv.pl.proportions(mg_cones_velo1)


# In[70]:


scv.pl.proportions(mg_rods_velo1)


# In[71]:


scv.pp.filter_and_normalize(mgtotal_I, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(mgtotal_I, n_pcs=30, n_neighbors=30)
scv.pp.filter_and_normalize(mgtotal_UI, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(mgtotal_UI, n_pcs=30, n_neighbors=30)


# In[57]:


scv.pp.filter_and_normalize(mg_rods_velo1, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(mg_rods_velo1, n_pcs=30, n_neighbors=30)



# In[58]:


scv.pp.filter_and_normalize(mg_cones_velo1, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(mg_cones_velo1, n_pcs=30, n_neighbors=30)


# In[74]:


scv.tl.velocity(mgtotal_I)
scv.tl.velocity(mgtotal_UI)


# In[59]:


scv.tl.velocity(mg_cones_velo1)


# In[60]:


scv.tl.velocity(mg_rods_velo1)


# In[42]:


mg_cones_velo1.obs


# In[43]:


mgtotal_I.obs


# In[79]:


mgtotal_UI.obs


# In[80]:


scv.tl.velocity_graph(mgtotal_I)
scv.tl.velocity_graph(mgtotal_UI)



# In[61]:


sc.pp.neighbors(mg_cones_velo1)
sc.pp.neighbors(mg_rods_velo1)


# In[62]:


scv.tl.velocity_graph(mg_cones_velo1)
scv.tl.velocity_graph(mg_rods_velo1)


# In[83]:


sc.pl.draw_graph(mgtotal,color=['Cell_ident'], layout = 'fa', legend_loc = 'on data',use_raw=False)


# In[84]:


scv.pl.velocity_embedding_grid(mgtotal_I, basis='draw_graph_fa', color = 'Cell_ident', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[85]:


scv.pl.velocity_embedding_grid(mgtotal_UI, basis='draw_graph_fa', color = 'Cell_ident', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[86]:


scv.pl.velocity_embedding_stream(mgtotal_I, basis='draw_graph_fa', color = 'Cell_ident', dpi = 300, density = 2, add_margin=0.2, arrow_size = 2, legend_loc = 'right margin')


# In[87]:


scv.pl.velocity_embedding_stream(mgtotal_UI, basis='draw_graph_fa', color = 'Cell_ident', dpi = 300, density = 2, add_margin=0.2, arrow_size = 2, legend_loc = 'right margin')


# In[88]:


sc.pl.draw_graph(mgtotal,color=['OTX2'], layout = 'fa', legend_loc = 'on data',use_raw=False)


# In[89]:


sc.pl.draw_graph(mg_cones,color=['Cell_ident'], layout = 'fa', legend_loc = 'on data',use_raw=False)


# In[63]:


sc.pl.draw_graph(mg_cones_velo1,color=['treatment'], layout = 'fa', legend_loc = 'on data',use_raw=False)


# In[224]:


sc.set_figure_params(figsize = [12,12])
sc.pl.draw_graph(mg_cones_velo1,color=['leiden'], layout = 'fa', legend_loc = 'on data',use_raw=False)


# In[223]:


sc.tl.leiden(mg_cones_velo1, resolution = 1)


# In[226]:


bridge_cones_CPI = mg_cones_velo1[mg_cones_velo1.obs['leiden'].isin(['9','0'])]


# In[228]:


sc.tl.leiden(bridge_cones_CPI, resolution = 1)


# In[229]:


sc.set_figure_params(figsize = [12,12])
sc.pl.draw_graph(bridge_cones_CPI,color=['leiden'], layout = 'fa', legend_loc = 'on data',use_raw=False)


# In[64]:


sc.pl.draw_graph(mg_rods_velo1,color=['treatment'], layout = 'fa', legend_loc = 'on data',use_raw=False)


# In[65]:


scv.pl.velocity_embedding_grid(mg_cones_velo1, basis='draw_graph_fa', color = 'Cell_ident', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[210]:


scv.pl.velocity_embedding_grid(mg_cones_velo1, basis='draw_graph_fa', color = 'Cell_ident', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 20, arrow_size = 6, density = 4)


# In[66]:


scv.pl.velocity_embedding_grid(mg_rods_velo1, basis='draw_graph_fa', color = 'Cell_ident', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[67]:


scv.tl.velocity_confidence(mg_rods_velo1)
keys =  'velocity_confidence','velocity_length'
scv.pl.scatter(mg_rods_velo1, c=keys, cmap='coolwarm', perc=[5, 95], basis='draw_graph_fa' )


# In[68]:


scv.tl.velocity_confidence(mg_cones_velo1)
keys =  'velocity_confidence','velocity_length'
scv.pl.scatter(mg_cones_velo1, c=keys, cmap='coolwarm', perc=[5, 95], basis='draw_graph_fa' )


# In[69]:


scv.tl.velocity_pseudotime(mg_rods_velo1)
scv.pl.scatter(mg_rods_velo1, color='velocity_pseudotime', cmap='gnuplot', basis='draw_graph_fa')


# In[70]:


scv.tl.velocity_pseudotime(mg_cones_velo1)
scv.pl.scatter(mg_cones_velo1, color='velocity_pseudotime', cmap='gnuplot', basis='draw_graph_fa')


# In[71]:


scv.tl.recover_dynamics(mg_cones_velo1)


# In[72]:


scv.tl.recover_dynamics(mg_rods_velo1)


# In[73]:


scv.tl.velocity(mg_cones_velo1, mode='dynamical')
scv.tl.velocity_graph(mg_cones_velo1)
scv.tl.velocity(mg_rods_velo1, mode='dynamical')
scv.tl.velocity_graph(mg_rods_velo1)


# In[74]:


scv.pl.velocity_embedding_grid(mg_rods_velo1, basis='draw_graph_fa', color = 'Cell_ident', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[75]:


scv.pl.velocity_embedding_grid(mg_cones_velo1, basis='draw_graph_fa', color = 'Cell_ident', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[76]:


df = mg_rods_velo1.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]

kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
    pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
    pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
    pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)

scv.get_df(mg_rods_velo1, 'fit*', dropna=True).head()


# In[77]:


df = mg_cones_velo1.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]

kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
    pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
    pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
    pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)

scv.get_df(mg_cones_velo1, 'fit*', dropna=True).head()


# In[78]:


scv.tl.latent_time(mg_cones_velo1)
scv.pl.scatter(mg_cones_velo1, color='latent_time', color_map='gnuplot', size=80, basis = 'draw_graph_fa')


# In[79]:


scv.tl.latent_time(mg_rods_velo1)
scv.pl.scatter(mg_rods_velo1, color='latent_time', color_map='gnuplot', size=80, basis = 'draw_graph_fa')


# In[80]:


mg_cones_velo1.obs


# In[81]:


CPI_Velo_cones_subset = mg_cones_velo1[mg_cones_velo1.obs['orig.ident'].isin(['U24CPI'])]


# In[89]:


CPI_Velo_cones_subset_nond = mg_cones_velo1[mg_cones_velo1.obs['orig.ident'].isin(['U24CPI'])]


# In[82]:


CPI_Velo_rods_subset = mg_rods_velo1[mg_rods_velo1.obs['orig.ident'].isin(['U24CPI'])]


# In[83]:


CUI_Velo_cones_subset = mg_cones_velo1[mg_cones_velo1.obs['orig.ident'].isin(['U24CUI'])]


# In[84]:


CUI_Velo_rods_subset = mg_rods_velo1[mg_rods_velo1.obs['orig.ident'].isin(['U24CUI'])]


# In[86]:


sc.pp.pca(CPI_Velo_cones_subset)
pca_projections = pd.DataFrame(CPI_Velo_cones_subset.obsm["X_pca"],index=CPI_Velo_cones_subset.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=10)
CPI_Velo_cones_subset.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(CPI_Velo_cones_subset,n_neighbors=30,use_rep="X_palantir")
CPI_Velo_cones_subset.obsm["X_pca2d"]=CPI_Velo_cones_subset.obsm["X_pca"][:,:3]
sc.tl.draw_graph(CPI_Velo_cones_subset, layout = 'fa',init_pos='X_pca2d')
sc.pl.draw_graph(CPI_Velo_cones_subset,color=['EK_anno'], layout = 'fa')


# In[90]:


#added later for non-dynamical mode
sc.pp.pca(CPI_Velo_cones_subset_nond)
pca_projections = pd.DataFrame(CPI_Velo_cones_subset_nond.obsm["X_pca"],index=CPI_Velo_cones_subset_nond.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=10)
CPI_Velo_cones_subset_nond.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(CPI_Velo_cones_subset_nond,n_neighbors=30,use_rep="X_palantir")
CPI_Velo_cones_subset_nond.obsm["X_pca2d"]=CPI_Velo_cones_subset_nond.obsm["X_pca"][:,:3]
sc.tl.draw_graph(CPI_Velo_cones_subset_nond, layout = 'fa',init_pos='X_pca2d')
sc.pl.draw_graph(CPI_Velo_cones_subset_nond,color=['EK_anno'], layout = 'fa')


# In[91]:


sc.pp.pca(CUI_Velo_cones_subset)
pca_projections = pd.DataFrame(CUI_Velo_cones_subset.obsm["X_pca"],index=CUI_Velo_cones_subset.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=10)
CUI_Velo_cones_subset.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(CUI_Velo_cones_subset,n_neighbors=30,use_rep="X_palantir")
CUI_Velo_cones_subset.obsm["X_pca2d"]=CUI_Velo_cones_subset.obsm["X_pca"][:,:3]
sc.tl.draw_graph(CUI_Velo_cones_subset, layout = 'fa',init_pos='X_pca2d')
sc.pl.draw_graph(CUI_Velo_cones_subset,color=['EK_anno'], layout = 'fa')


# In[92]:


sc.pp.pca(CPI_Velo_rods_subset)
pca_projections = pd.DataFrame(CPI_Velo_rods_subset.obsm["X_pca"],index=CPI_Velo_rods_subset.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=10)
CPI_Velo_rods_subset.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(CPI_Velo_rods_subset,n_neighbors=30,use_rep="X_palantir")
CPI_Velo_rods_subset.obsm["X_pca2d"]=CPI_Velo_rods_subset.obsm["X_pca"][:,:3]
sc.tl.draw_graph(CPI_Velo_rods_subset, layout = 'fa',init_pos='X_pca2d')
sc.pl.draw_graph(CPI_Velo_rods_subset,color=['EK_anno'], layout = 'fa')


# In[93]:


sc.pp.pca(CUI_Velo_rods_subset)
pca_projections = pd.DataFrame(CUI_Velo_rods_subset.obsm["X_pca"],index=CUI_Velo_rods_subset.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=10)
CUI_Velo_rods_subset.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(CUI_Velo_rods_subset,n_neighbors=30,use_rep="X_palantir")
CUI_Velo_rods_subset.obsm["X_pca2d"]=CUI_Velo_rods_subset.obsm["X_pca"][:,:3]
sc.tl.draw_graph(CUI_Velo_rods_subset, layout = 'fa',init_pos='X_pca2d')
sc.pl.draw_graph(CUI_Velo_rods_subset,color=['EK_anno'], layout = 'fa')


# In[94]:


scv.pl.proportions(CPI_Velo_cones_subset)


# In[95]:


scv.pl.proportions(CUI_Velo_cones_subset)


# In[96]:


scv.pl.proportions(CPI_Velo_rods_subset)


# In[97]:


scv.pl.proportions(CUI_Velo_rods_subset)


# In[98]:


scv.tl.velocity(CPI_Velo_cones_subset)
scv.tl.velocity(CUI_Velo_cones_subset)


# In[99]:


scv.tl.velocity(CPI_Velo_rods_subset)
scv.tl.velocity(CUI_Velo_rods_subset)


# In[100]:


scv.tl.velocity(CPI_Velo_cones_subset_nond)


# In[101]:


scv.tl.velocity_graph(CPI_Velo_cones_subset_nond)


# In[102]:


scv.tl.velocity_graph(CPI_Velo_cones_subset)
scv.tl.velocity_graph(CUI_Velo_cones_subset)


# In[103]:


scv.tl.velocity_graph(CPI_Velo_rods_subset)
scv.tl.velocity_graph(CUI_Velo_rods_subset)


# In[104]:


scv.pl.velocity_embedding_grid(CPI_Velo_cones_subset, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[105]:


scv.pl.velocity_embedding_grid(CUI_Velo_cones_subset, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[106]:


scv.pl.velocity_embedding_grid(CPI_Velo_rods_subset, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[107]:


scv.pl.velocity_embedding_grid(CUI_Velo_rods_subset, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[108]:


scv.tl.recover_dynamics(CUI_Velo_cones_subset)
scv.tl.recover_dynamics(CPI_Velo_cones_subset)


# In[109]:


scv.tl.recover_dynamics(CPI_Velo_rods_subset)
scv.tl.recover_dynamics(CUI_Velo_rods_subset)


# In[110]:


scv.tl.velocity(CUI_Velo_cones_subset, mode='dynamical')
scv.tl.velocity_graph(CUI_Velo_cones_subset)
scv.tl.velocity(CPI_Velo_cones_subset, mode='dynamical')
scv.tl.velocity_graph(CPI_Velo_cones_subset)


# In[111]:


scv.tl.velocity(CPI_Velo_rods_subset, mode='dynamical')
scv.tl.velocity_graph(CPI_Velo_rods_subset)
scv.tl.velocity(CUI_Velo_rods_subset, mode='dynamical')
scv.tl.velocity_graph(CUI_Velo_rods_subset)


# In[112]:


scv.pl.velocity_embedding_grid(CUI_Velo_cones_subset, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[113]:


scv.pl.velocity_embedding_grid(CPI_Velo_cones_subset, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[114]:


scv.pl.velocity_embedding_grid(CPI_Velo_rods_subset, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[115]:


scv.pl.velocity_embedding_grid(CUI_Velo_rods_subset, basis='draw_graph_fa', color = 'EK_anno', dpi = 300, add_margin = 0.2, figsize = (4,4), arrow_length = 8, arrow_size = 2, density = 1)


# In[116]:


scv.tl.latent_time(CPI_Velo_cones_subset)
scv.pl.scatter(CPI_Velo_cones_subset, color='latent_time', color_map='gnuplot', size=80, basis = 'draw_graph_fa')


# In[117]:


scv.tl.latent_time(CUI_Velo_cones_subset)
scv.pl.scatter(CUI_Velo_cones_subset, color='latent_time', color_map='gnuplot', size=80, basis = 'draw_graph_fa')


# In[118]:


scv.tl.latent_time(CPI_Velo_rods_subset)
scv.pl.scatter(CPI_Velo_rods_subset, color='latent_time', color_map='gnuplot', size=80, basis = 'draw_graph_fa')


# In[119]:


scv.tl.latent_time(CUI_Velo_rods_subset)
scv.pl.scatter(CUI_Velo_rods_subset, color='latent_time', color_map='gnuplot', size=80, basis = 'draw_graph_fa')


# In[120]:


scv.tl.velocity_confidence(CPI_Velo_cones_subset)
keys =  'velocity_confidence','velocity_length'
scv.pl.scatter(CPI_Velo_cones_subset, c=keys, cmap='coolwarm', perc=[5, 95], basis = 'draw_graph_fa')


# In[121]:


scv.tl.velocity_confidence(CUI_Velo_cones_subset)
keys =  'velocity_confidence','velocity_length'
scv.pl.scatter(CUI_Velo_cones_subset, c=keys, cmap='coolwarm', perc=[5, 95], basis = 'draw_graph_fa')


# In[122]:


scv.tl.velocity_confidence(CPI_Velo_rods_subset)
keys =  'velocity_confidence','velocity_length'
scv.pl.scatter(CPI_Velo_rods_subset, c=keys, cmap='coolwarm', perc=[5, 95], basis = 'draw_graph_fa')


# In[123]:


scv.tl.velocity_confidence(CUI_Velo_rods_subset)
keys =  'velocity_confidence','velocity_length'
scv.pl.scatter(CUI_Velo_rods_subset, c=keys, cmap='coolwarm', perc=[5, 95], basis = 'draw_graph_fa')


# In[124]:


top_genes = CPI_Velo_cones_subset.var['fit_likelihood'].sort_values(ascending=False).index[:160]
scv.pl.heatmap(CPI_Velo_cones_subset, var_names=top_genes, sortby='latent_time', col_color='Cell_ident', n_convolve=300,yticklabels=True, figsize = (8,32))


# In[126]:


sc.tl.score_genes(CPI_Velo_cones_subset, gene_list = ['LEMD1','MCC','IMPG1','SLC1A7','MAOA','RS1','COL5A1','SLC8A1','PLAAT1','FSTL5','SDK2','IGSF21'], score_name = 'Maturation-CPI', use_raw = False)


# In[125]:


sc.tl.score_genes(CPI_Velo_cones_subset, gene_list = ['LGALS3','GAD2','SPP1','AKAP12','SFRP2','UBE2T','PMEPA1','TPM1','ZFP36L1','IFITM3','PRSS23','TTYH1'], score_name = 'Early genes-CPI', use_raw = False)


# In[129]:


sc.pl.draw_graph(CPI_Velo_cones_subset,color=['Maturation-CPI','Early genes-CPI'], vmax = [0.6,1])


# In[131]:


sc.tl.score_genes(CUI_Velo_cones_subset, gene_list = ['LEMD1','MCC','IMPG1','SLC1A7','MAOA','RS1','COL5A1','SLC8A1','PLAAT1','FSTL5','SDK2','IGSF21'], score_name = 'Maturation-CPI', use_raw = False)


# In[132]:


sc.tl.score_genes(CUI_Velo_cones_subset, gene_list = ['LGALS3','GAD2','SPP1','AKAP12','SFRP2','UBE2T','PMEPA1','TPM1','ZFP36L1','IFITM3','PRSS23','TTYH1'], score_name = 'Early genes-CPI', use_raw = False)


# In[135]:


sc.pl.draw_graph(CUI_Velo_cones_subset,color=['Maturation-CPI','Early genes-CPI'], vmax = [0.6,1])


# In[136]:


top_genes = CUI_Velo_cones_subset.var['fit_likelihood'].sort_values(ascending=False).index[:160]
scv.pl.heatmap(CUI_Velo_cones_subset, var_names=top_genes, sortby='latent_time', col_color='Cell_ident', n_convolve=300,yticklabels=True, figsize = (8,32))


# In[137]:


sc.tl.score_genes(CPI_Velo_cones_subset, gene_list = ['KCNV2','RRAD','MAP2','PLAAT1','KCNK15','SALL1','GREM2','CRABP2','PDE6H','PPP2R2B','XRCC4','PPP4R4'], score_name = 'Maturation-CUI', use_raw = False)


# In[138]:


sc.tl.score_genes(CPI_Velo_cones_subset, gene_list = ['CAV1','AP1S2','ANXA2','SERPINE1','TPM1','RND3','LYPD1','SLC4A7','CALB1','GAD2','AKAP12','DCBLD2'], score_name = 'Early genes-CUI', use_raw = False)


# In[139]:


sc.tl.score_genes(CUI_Velo_cones_subset, gene_list = ['KCNV2','RRAD','MAP2','PLAAT1','KCNK15','SALL1','GREM2','CRABP2','PDE6H','PPP2R2B','XRCC4','PPP4R4'], score_name = 'Maturation-CUI', use_raw = False)


# In[140]:


sc.tl.score_genes(CUI_Velo_cones_subset, gene_list = ['CAV1','AP1S2','ANXA2','SERPINE1','TPM1','RND3','LYPD1','SLC4A7','CALB1','GAD2','AKAP12','DCBLD2'], score_name = 'Early genes-CUI', use_raw = False)


# In[146]:


sc.pl.draw_graph(CPI_Velo_cones_subset,color=['Maturation-CUI','Early genes-CUI'], vmax = [1.2,1.2])


# In[145]:


sc.pl.draw_graph(CUI_Velo_cones_subset,color=['Maturation-CUI','Early genes-CUI'], vmax = [1.2,1.2])


# In[147]:


top_genes = CPI_Velo_rods_subset.var['fit_likelihood'].sort_values(ascending=False).index[:160]
scv.pl.heatmap(CPI_Velo_rods_subset, var_names=top_genes, sortby='latent_time', col_color='Cell_ident', n_convolve=300,yticklabels=True, figsize = (8,32))


# In[243]:





# In[196]:


sc.tl.score_genes(CPI_Velo_rods_subset, gene_list = ['TRH','SFRP2','PRKD1','CLU','GADD45A','RTN4','HES5','APCDD1L-DT','FOS','MYBL1','MDK','TGFB2'], score_name = 'Maturation-CPI', use_raw = False)
sc.tl.score_genes(CPI_Velo_rods_subset, gene_list = ['PALMD','NEUROD4','NFIB','SSX2IP','PKM','CASZ1','RLF','RORB','LRRC8D','UBE4B','CCNL2','GPC3'], score_name = 'Early genes-CPI', use_raw = False)
sc.tl.score_genes(CPI_Velo_rods_subset, gene_list = ['PTPRZ1','TRH','TGFB2','MIAT','GADD45A','APCDD1L-DT','HES5','EFCAB7','MAST2','RP1','GNGT2','MYO9A'], score_name = 'Maturation-CUI', use_raw = False)
sc.tl.score_genes(CPI_Velo_rods_subset, gene_list = ['MYL1','MIR503HG','TTN','SLF2','MYL4','CLU','FABP5','RTN4','SFRP2','ERBB4','OLFM3','USP48'], score_name = 'Early genes-CUI', use_raw = False)

sc.tl.score_genes(CUI_Velo_rods_subset, gene_list = ['TRH','SFRP2','PRKD1','CLU','GADD45A','RTN4','HES5','APCDD1L-DT','FOS','MYBL1','MDK','TGFB2'], score_name = 'Maturation-CPI', use_raw = False)
sc.tl.score_genes(CUI_Velo_rods_subset, gene_list = ['PALMD','NEUROD4','NFIB','SSX2IP','PKM','CASZ1','RLF','RORB','LRRC8D','UBE4B','CCNL2','GPC3'], score_name = 'Early genes-CPI', use_raw = False)
sc.tl.score_genes(CUI_Velo_rods_subset, gene_list = ['PTPRZ1','TRH','TGFB2','MIAT','GADD45A','APCDD1L-DT','HES5','EFCAB7','MAST2','RP1','GNGT2','MYO9A'], score_name = 'Maturation-CUI', use_raw = False)
sc.tl.score_genes(CUI_Velo_rods_subset, gene_list = ['MYL1','MIR503HG','TTN','SLF2','MYL4','CLU','FABP5','RTN4','SFRP2','ERBB4','OLFM3','USP48'], score_name = 'Early genes-CUI', use_raw = False)


# In[198]:


sc.pl.draw_graph(CPI_Velo_rods_subset,color=['Maturation-CPI','Early genes-CPI','Maturation-CUI','Early genes-CUI'], vmax = [1.4,1.4,1.2,0.5])


# In[199]:


sc.pl.draw_graph(CUI_Velo_rods_subset,color=['Maturation-CPI','Early genes-CPI','Maturation-CUI','Early genes-CUI'], vmax = [1.4,1.4,1.2,0.5])


# In[239]:


top_genes = CUI_Velo_rods_subset.var['fit_likelihood'].sort_values(ascending=False).index[:160]
scv.pl.heatmap(CUI_Velo_rods_subset, var_names=top_genes, sortby='latent_time', col_color='Cell_ident', n_convolve=300,yticklabels=True, figsize = (8,32))


# In[131]:


sc.pl.draw_graph(CUI_Velo_cones_subset,color=['ANXA2','LGALS3','PDE6H','CA2'], layout = 'fa', legend_loc = 'on data',use_raw=False)


# In[149]:


sc.tl.leiden(CUI_Velo_cones_subset)


# In[150]:


sc.tl.leiden(CPI_Velo_cones_subset)


# In[151]:


sc.tl.leiden(CPI_Velo_rods_subset)


# In[152]:


sc.tl.leiden(CUI_Velo_rods_subset)


# In[153]:


CUI_Velo_cones_subset.uns['neighbors']['distances'] = CUI_Velo_cones_subset.obsp['distances']
CUI_Velo_cones_subset.uns['neighbors']['connectivities'] = CUI_Velo_cones_subset.obsp['connectivities']

scv.tl.paga(CUI_Velo_cones_subset, groups='leiden')
df = scv.get_df(CUI_Velo_cones_subset, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[154]:


CPI_Velo_cones_subset.uns['neighbors']['distances'] = CPI_Velo_cones_subset.obsp['distances']
CPI_Velo_cones_subset.uns['neighbors']['connectivities'] = CPI_Velo_cones_subset.obsp['connectivities']

scv.tl.paga(CPI_Velo_cones_subset, groups='leiden')
df = scv.get_df(CPI_Velo_cones_subset, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[155]:


CPI_Velo_rods_subset.uns['neighbors']['distances'] = CPI_Velo_rods_subset.obsp['distances']
CPI_Velo_rods_subset.uns['neighbors']['connectivities'] = CPI_Velo_rods_subset.obsp['connectivities']

scv.tl.paga(CPI_Velo_rods_subset, groups='leiden')
df = scv.get_df(CPI_Velo_rods_subset, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[156]:


CUI_Velo_rods_subset.uns['neighbors']['distances'] = CUI_Velo_rods_subset.obsp['distances']
CUI_Velo_rods_subset.uns['neighbors']['connectivities'] = CUI_Velo_rods_subset.obsp['connectivities']

scv.tl.paga(CUI_Velo_rods_subset, groups='leiden')
df = scv.get_df(CUI_Velo_rods_subset, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[157]:


scv.pl.paga(CUI_Velo_cones_subset, basis='draw_graph_fa', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)


# In[158]:


scv.pl.paga(CPI_Velo_cones_subset, basis='draw_graph_fa', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)


# In[159]:


scv.pl.paga(CPI_Velo_rods_subset, basis='draw_graph_fa', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)


# In[160]:


scv.pl.paga(CUI_Velo_rods_subset, basis='draw_graph_fa', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)


# In[163]:


#CUI_Velo_cones_subset
scv.tl.score_genes_cell_cycle(CUI_Velo_cones_subset, use_raw = False)
scv.pl.scatter(CUI_Velo_cones_subset, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95], basis = 'draw_graph_fa')


# In[164]:


scv.tl.score_genes_cell_cycle(CPI_Velo_cones_subset, use_raw = False)
scv.pl.scatter(CPI_Velo_cones_subset, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95], basis = 'draw_graph_fa')


# In[165]:


scv.tl.score_genes_cell_cycle(CUI_Velo_rods_subset, use_raw = False)
scv.pl.scatter(CUI_Velo_rods_subset, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95], basis = 'draw_graph_fa')


# In[166]:


scv.tl.score_genes_cell_cycle(CPI_Velo_rods_subset, use_raw = False)
scv.pl.scatter(CPI_Velo_rods_subset, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95], basis = 'draw_graph_fa')


# In[167]:


import cellrank as cr


# In[168]:


CPI_Velo_cones_subset_nond.obs


# In[169]:


sc.tl.leiden(CPI_Velo_cones_subset_nond)


# In[170]:


sc.pl.draw_graph(CPI_Velo_cones_subset_nond,color=['leiden'], layout = 'fa')


# In[171]:


cr.tl.terminal_states(CPI_Velo_cones_subset_nond, cluster_key="leiden", weight_connectivities=0.2, n_states = 1)


# In[172]:


cr.tl.terminal_states(CPI_Velo_cones_subset, cluster_key="leiden", weight_connectivities=0.2)


# In[173]:


cr.pl.terminal_states(CPI_Velo_cones_subset_nond, basis = 'draw_graph_fa')


# In[174]:


cr.pl.terminal_states(CPI_Velo_cones_subset, basis = 'draw_graph_fa')


# In[175]:


cr.tl.initial_states(CPI_Velo_cones_subset_nond, cluster_key="leiden", n_states = 1)
cr.pl.initial_states(CPI_Velo_cones_subset_nond, discrete=True)


# In[176]:


cr.pl.initial_states(CPI_Velo_cones_subset_nond, discrete=True, basis = 'draw_graph_fa')


# In[177]:


cr.tl.initial_states(CPI_Velo_cones_subset, cluster_key="leiden")
cr.pl.initial_states(CPI_Velo_cones_subset, discrete=True)


# In[178]:


cr.pl.initial_states(CPI_Velo_cones_subset, discrete=True, basis = 'draw_graph_fa')


# In[179]:


cr.tl.lineages(CPI_Velo_cones_subset_nond)
cr.pl.lineages(CPI_Velo_cones_subset_nond, same_plot=False, basis = 'draw_graph_fa')


# In[180]:


cr.tl.lineages(CPI_Velo_cones_subset)
cr.pl.lineages(CPI_Velo_cones_subset, same_plot=False, basis = 'draw_graph_fa')


# In[181]:


scv.tl.recover_latent_time(
    CPI_Velo_cones_subset_nond, root_key="initial_states_probs", end_key="terminal_states_probs"
)


# In[182]:


scv.tl.paga(
    CPI_Velo_cones_subset_nond,
    groups="leiden",
    root_key="initial_states_probs",
    end_key="terminal_states_probs",
    use_time_prior="velocity_pseudotime",
)


# In[184]:


sc.pl.draw_graph(CPI_Velo_cones_subset_nond,color=['Cell_ident','leiden','EK_anno'], layout = 'fa')


# In[186]:


sc.pl.draw_graph(CPI_Velo_cones_subset_nond,color=['leiden'], layout = 'fa')


# In[ ]:





# In[237]:


sc.tl.rank_genes_groups(CPI_Velo_cones_subset_nond, 'leiden', method='logreg', use_raw = False)
sc.pl.rank_genes_groups(CPI_Velo_cones_subset_nond, n_genes=10, sharey=False)


# In[187]:


cr.pl.cluster_fates(
    CPI_Velo_cones_subset_nond,
    mode="paga_pie",
    cluster_key="leiden",
    basis="draw_graph_fa",
    legend_kwargs={"loc": "top right out"},
    legend_loc="top left out",
    node_size_scale=1,
    edge_width_scale=1,
    max_edge_width=4,
    title="directed PAGA",
)


# In[238]:


cr.pl.circular_projection(CPI_Velo_cones_subset_nond, keys="leiden", legend_loc="right")


# In[188]:


scv.tl.paga(
    CPI_Velo_cones_subset_nond,
    groups="leiden",
    root_key="initial_states_probs",
    end_key="terminal_states_probs",
    use_time_prior="latent_time",
)


# In[189]:


CPI_Velo_cones_subset_nond.obs


# In[190]:


cr.pl.cluster_fates(
    CPI_Velo_cones_subset_nond,
    mode="paga_pie",
    cluster_key="leiden",
    basis="draw_graph_fa",
    legend_kwargs={"loc": "top right out"},
    legend_loc="top left out",
    node_size_scale=1,
    edge_width_scale=1,
    max_edge_width=4,
    title="directed PAGA",
)


# In[191]:


cr.tl.lineage_drivers(CPI_Velo_cones_subset_nond)


# In[192]:


cr.pl.lineage_drivers(CPI_Velo_cones_subset_nond, lineage="3", n_genes=10, basis = 'draw_graph_fa')


# In[193]:


# compue DPT, starting from CellRank defined root cell
root_idx = np.where(CPI_Velo_cones_subset_nond.obs["initial_states"] == "13")[0][0]
CPI_Velo_cones_subset_nond.uns["iroot"] = root_idx
sc.tl.dpt(CPI_Velo_cones_subset_nond)

scv.pl.scatter(
    CPI_Velo_cones_subset_nond, basis = 'draw_graph_fa',
    color=["leiden", root_idx, "latent_time", "dpt_pseudotime"],
    fontsize=16,
    cmap="viridis",
    perc=[2, 98],
    colorbar=True,
    rescale_color=[0, 1],
    title=["leiden", "root cell", "latent time", "dpt pseudotime"],
)


# In[195]:


model = cr.ul.models.GAM(CPI_Velo_cones_subset_nond)
cr.pl.gene_trends(
    CPI_Velo_cones_subset_nond,
    model=model, 
    data_key="X",
    genes=["RRAD", "PDE6H"],
    ncols=4,
    time_key="latent_time",
    same_plot=True,use_raw = False,
    hide_cells=True,
    figsize=(15, 4),
    n_test_points=200
)


# In[255]:


model


# In[256]:


plt.subplots(layout="constrained")


# In[260]:


scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
cr.settings.verbosity = 2


# In[270]:


cr.pl.heatmap(
    CPI_Velo_cones_subset_nond,
    model,
    genes=CPI_Velo_cones_subset_nond.varm['terminal_lineage_drivers']["5 corr"].sort_values(ascending=False).index[:100],
    show_absorption_probabilities=True,
    lineages="5",
    n_jobs=1,
    backend="loky",
)


# In[271]:


CPI_Velo_cones_subset_nond.varm


# In[268]:


cr.tl.lineage_drivers(CPI_Velo_cones_subset_nond, cluster_key = 'leiden')


# In[274]:


mg


# In[273]:


mgtotal


# In[279]:


old_to_new = {
'0':'Cones',
'1':'Progenitors',
'2':'Progenitors',
'3':'Cones',
'4':'Inner retina neurons',
'5':'Cones',
'6':'Muller glia',
'7':'Progenitors',
'8':'Cones',
'9':'Progenitors',
'10':'Rods',
'11':'Muller glia',
'12':'Progenitors',
'13':'Glia',
'14':'Inner retina neurons',
'15':'Rods',
'16':'Microglia',
'17':'RPE',
'18':'Cones'}
mgtotal.obs['EK_anno'] = (
mgtotal.obs['leiden']
.map(old_to_new).astype('category')
)


# In[275]:


sc.pl.draw_graph(mgtotal,color=['leiden'], layout = 'fa')


# In[278]:


sc.pl.draw_graph(mgtotal,color=['Cell_ident'], layout = 'fa')


# In[280]:


sc.pl.draw_graph(mgtotal,color=['EK_anno'], layout = 'fa')


# In[284]:


mgtotal.write_h5ad('/mnt/c/Users/Emil/10X/dogs/mgtotal_EK_anno.h5ad') #del mgtotal.raw


# In[283]:


del mgtotal.raw


# In[285]:


mgtotal


# In[196]:


sc.tl.score_genes(CUI_Velo_cones_subset, gene_list = ['ACTL6','BADGRB3','AGRN','ANKS1A','APP','B4GALT5','B4GALT6','BCL11A','BCL2','C1QA','C1QL1','C3','CDKN1C','CLN5','CNTN2','CNTNAP2','CX3CL1','EDNRA','EDNRB','EPHA8','FARP2','FEV','GLDN','KCNB1','KCNIP2','KDM1A','LGI4','LRRK2','MAP3K13','MECP2','MTCH1','MYOC','NR4A2','NRCAM','NTN4','PICK1','RAC3','RB1','RET','RND1','SCARF1','SCLT1','SPTBN4','SRRM4','TBX6','VSX1'], score_name = 'GOBP_NEURON_MATURATION', use_raw = False)
sc.tl.score_genes(CPI_Velo_cones_subset, gene_list = ['ACTL6','BADGRB3','AGRN','ANKS1A','APP','B4GALT5','B4GALT6','BCL11A','BCL2','C1QA','C1QL1','C3','CDKN1C','CLN5','CNTN2','CNTNAP2','CX3CL1','EDNRA','EDNRB','EPHA8','FARP2','FEV','GLDN','KCNB1','KCNIP2','KDM1A','LGI4','LRRK2','MAP3K13','MECP2','MTCH1','MYOC','NR4A2','NRCAM','NTN4','PICK1','RAC3','RB1','RET','RND1','SCARF1','SCLT1','SPTBN4','SRRM4','TBX6','VSX1'], score_name = 'GOBP_NEURON_MATURATION', use_raw = False)
sc.tl.score_genes(CUI_Velo_rods_subset, gene_list = ['ACTL6','BADGRB3','AGRN','ANKS1A','APP','B4GALT5','B4GALT6','BCL11A','BCL2','C1QA','C1QL1','C3','CDKN1C','CLN5','CNTN2','CNTNAP2','CX3CL1','EDNRA','EDNRB','EPHA8','FARP2','FEV','GLDN','KCNB1','KCNIP2','KDM1A','LGI4','LRRK2','MAP3K13','MECP2','MTCH1','MYOC','NR4A2','NRCAM','NTN4','PICK1','RAC3','RB1','RET','RND1','SCARF1','SCLT1','SPTBN4','SRRM4','TBX6','VSX1'], score_name = 'GOBP_NEURON_MATURATION', use_raw = False)
sc.tl.score_genes(CPI_Velo_rods_subset, gene_list = ['ACTL6','BADGRB3','AGRN','ANKS1A','APP','B4GALT5','B4GALT6','BCL11A','BCL2','C1QA','C1QL1','C3','CDKN1C','CLN5','CNTN2','CNTNAP2','CX3CL1','EDNRA','EDNRB','EPHA8','FARP2','FEV','GLDN','KCNB1','KCNIP2','KDM1A','LGI4','LRRK2','MAP3K13','MECP2','MTCH1','MYOC','NR4A2','NRCAM','NTN4','PICK1','RAC3','RB1','RET','RND1','SCARF1','SCLT1','SPTBN4','SRRM4','TBX6','VSX1'], score_name = 'GOBP_NEURON_MATURATION', use_raw = False)


# In[199]:


sc.pl.draw_graph(CPI_Velo_cones_subset,color=['GOBP_NEURON_MATURATION'], vmax = [0.6]) 


# In[200]:


sc.pl.draw_graph(CUI_Velo_cones_subset,color=['GOBP_NEURON_MATURATION'], vmax = [0.6])


# In[201]:


sc.pl.draw_graph(CUI_Velo_rods_subset,color=['GOBP_NEURON_MATURATION'], vmax = [0.6])


# In[202]:


sc.pl.draw_graph(CPI_Velo_rods_subset,color=['GOBP_NEURON_MATURATION'], vmax = [0.6])


# In[206]:


sc.tl.score_genes(CPI_Velo_cones_subset, gene_list = ['ADGRB3','ADGRL1','ANAPC2','ARHGEF15','BCAN','C1QL2','CAMK2B','CDC20','CHRDL1','CLSTN1','CNTNAP1','CX3CR1','DAB2IP','DISC1','IGSF21','NEFL','NEURL1','NEUROD2','NFATC4','NFIA','NLGN1','NRXN1','PALM','PFN1','PTEN','RELN','ROCK1','SEZ6','SEZ6L','SEZ6L2','SHANK1','SYBU','VPS35','YWHAZ'], score_name = 'GOBP_SYNAPSE_MATURATION', use_raw = False)
sc.tl.score_genes(CUI_Velo_cones_subset, gene_list = ['ADGRB3','ADGRL1','ANAPC2','ARHGEF15','BCAN','C1QL2','CAMK2B','CDC20','CHRDL1','CLSTN1','CNTNAP1','CX3CR1','DAB2IP','DISC1','IGSF21','NEFL','NEURL1','NEUROD2','NFATC4','NFIA','NLGN1','NRXN1','PALM','PFN1','PTEN','RELN','ROCK1','SEZ6','SEZ6L','SEZ6L2','SHANK1','SYBU','VPS35','YWHAZ'], score_name = 'GOBP_SYNAPSE_MATURATION', use_raw = False)
sc.tl.score_genes(CPI_Velo_rods_subset, gene_list = ['ADGRB3','ADGRL1','ANAPC2','ARHGEF15','BCAN','C1QL2','CAMK2B','CDC20','CHRDL1','CLSTN1','CNTNAP1','CX3CR1','DAB2IP','DISC1','IGSF21','NEFL','NEURL1','NEUROD2','NFATC4','NFIA','NLGN1','NRXN1','PALM','PFN1','PTEN','RELN','ROCK1','SEZ6','SEZ6L','SEZ6L2','SHANK1','SYBU','VPS35','YWHAZ'], score_name = 'GOBP_SYNAPSE_MATURATION', use_raw = False)
sc.tl.score_genes(CUI_Velo_rods_subset, gene_list = ['ADGRB3','ADGRL1','ANAPC2','ARHGEF15','BCAN','C1QL2','CAMK2B','CDC20','CHRDL1','CLSTN1','CNTNAP1','CX3CR1','DAB2IP','DISC1','IGSF21','NEFL','NEURL1','NEUROD2','NFATC4','NFIA','NLGN1','NRXN1','PALM','PFN1','PTEN','RELN','ROCK1','SEZ6','SEZ6L','SEZ6L2','SHANK1','SYBU','VPS35','YWHAZ'], score_name = 'GOBP_SYNAPSE_MATURATION', use_raw = False)
sc.tl.score_genes(CPI_Velo_cones_subset_nond, gene_list = ['ADGRB3','ADGRL1','ANAPC2','ARHGEF15','BCAN','C1QL2','CAMK2B','CDC20','CHRDL1','CLSTN1','CNTNAP1','CX3CR1','DAB2IP','DISC1','IGSF21','NEFL','NEURL1','NEUROD2','NFATC4','NFIA','NLGN1','NRXN1','PALM','PFN1','PTEN','RELN','ROCK1','SEZ6','SEZ6L','SEZ6L2','SHANK1','SYBU','VPS35','YWHAZ'], score_name = 'GOBP_SYNAPSE_MATURATION', use_raw = False)


# In[204]:


sc.pl.draw_graph(CPI_Velo_cones_subset,color=['GOBP_SYNAPSE_MATURATION']) #, vmax = [0.6]


# In[207]:


sc.pl.draw_graph(CPI_Velo_cones_subset_nond,color=['GOBP_SYNAPSE_MATURATION']) #, vmax = [0.6]


# In[205]:


sc.pl.draw_graph(CUI_Velo_cones_subset,color=['GOBP_SYNAPSE_MATURATION'])


# In[208]:


sc.pl.draw_graph(CUI_Velo_rods_subset,color=['GOBP_SYNAPSE_MATURATION'])


# In[209]:


sc.pl.draw_graph(CPI_Velo_rods_subset,color=['GOBP_SYNAPSE_MATURATION'])


# In[ ]:




