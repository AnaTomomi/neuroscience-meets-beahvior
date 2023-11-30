###############################################################################
# This script visualizes the brain visualization and meta-analysis data in    #
# brain surface plots.                                                        #    
#                                                                             #
# Author: ana.trianahoyos@aalto.fi                                            #
###############################################################################

from neuromaps import transforms
from surfplot import Plot
from neuromaps.datasets import fetch_fslr

import os
import numpy as np
import nibabel as nib

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import seaborn as sns


#Create the colormaps to use
def create_discrete_colormap(cmap_name, i):
    # Load the continuous colormap
    cmap = plt.get_cmap(cmap_name)

    # Create discrete colormap by selecting colors from the continuous colormap
    colors = cmap(np.linspace(0, 1, i))
    discrete_cmap = ListedColormap(colors)

    # Create a normalizer for the colormap with i bins
    norm = BoundaryNorm(np.linspace(0, i, i+1), discrete_cmap.N)

    return discrete_cmap, norm


#Load the data to plot
path = "../results/brain_masks"
images = [file for file in os.listdir(path) if file.endswith(".nii.gz")]
cmap_name = 'YlOrRd'

for image in images:
    img = f'{path}/{image}'
    surf_img = transforms.mni152_to_fslr(img, '32k')
    label = image.split('-')[1].split('-')[0]
    
    nifti_img = nib.load(img)
    nifti_data = nifti_img.get_fdata()
    i = int(np.max(nifti_data))
    
    #Create the colormap
    discrete_cmap, norm = create_discrete_colormap(cmap_name, i)

    #Start plotting the surface data
    surfaces = fetch_fslr()
    lh, rh = surfaces['inflated']

    #Create the plot
    p = Plot(lh, rh)

    # shading of the brain
    lh_sulc, rh_sulc = surfaces['sulc']
    p.add_layer({'left': lh_sulc, 'right': rh_sulc}, cmap='binary_r', cbar=False)
    color_range = (0, i+1)

    # add the masks
    fslr_lh, fslr_rh = surf_img
    p.add_layer({'left': fslr_lh, 'right': fslr_rh}, cmap=discrete_cmap,  
                color_range=color_range, cbar_label='')

    cbar_kws = dict(outer_labels_only=True, pad=.02, n_ticks=2, decimals=0)

    fig = p.build(cbar_kws=cbar_kws)
    
    cbar = fig.axes[1]
    cbar.set_xticks(np.linspace(0, i, i)+0.5)
    cbar.set_xticklabels(np.arange(i)+1, fontsize=14)
    
    fig.axes[1].set_xlabel('number of studies', labelpad=10, fontsize=16)
    fig.show()

    plt.savefig(f'{savepath}/{image[:-7]}.pdf', dpi=300)
    
    del p, fslr_lh, fslr_rh, surf_img, img, label, i, fig