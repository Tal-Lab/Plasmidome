"""
Created on 11/05/2023

Author: Lucy Androsiuk

Description: The grid_image() function generates a 2x2 grid of subplots with 3 heatmaps,
             one for the presence of feature 1, one for the presence of feature 2,
             and one for the correlation between the two.
             The heatmaps are sorted in clustering order and colored based on their cluster assignments.
             The function uses the seaborn library for creating the heatmaps and the gridspec module for arranging the subplots.

Requirements: 3 matrices (dataframes): presence_feature1 (presence matrix for feature 1),
                                       presence_feature2 (presence matrix for feature 2),
                                       corr_feature12 (correlation between presence matrices),
              lists of both Features and locations in the clustering order,
              3 dataframes (or dictionaries) with clusters for feature1 (cluster_feature1_df),
                                                               feature2 (cluster_feature2_df),
                                                               locations (cluster_loc_df).

Based on stackoverflow answer: https://stackoverflow.com/questions/73433022/adding-row-colors-to-a-heatmap
"""

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import matplotlib.gridspec as gridspec
import os

visuals = r'path/to/folder/with/visualizations'

def grid_image():
    '''This function puts 3 heatmaps (clustermaps) into a grid:
    1st (presence matrix for feature 1) in upper right corner (ax2 [0,1]),
    2nd (presence matrix for feature 2) in lower left corner (ax3 [1,0]),
    3rd (correlation between presence matrices) in lower right corner (ax4 [1,1]).
    FYI: Upper left corner may remain free or filled with any other graph.'''
    feature1_order = ['A','B','C','D'] # list of feature1 in clustering order
    feature2_order = ['a','b','c','d'] # list of feature2 in clustering order
    loc_order = [1,2,3,4,5] # list of locations in clustering order
    data1 = presence_feature1 # presence matrix for feature 1
    data1 = data1.reindex(loc_order) # sorting rows (locations) in clustering order
    data1 = data1.reindex(columns = feature1_order) # sorting columns (feature1) in clustering order
    data2 = presence_feature2 # presence matrix for feature 2
    data2 = data2.reindex(feature2_order) # sorting rows (feature2) in clustering order
    data2 = data2.reindex(columns = loc_order) # sorting columns (locations) in clustering order
    data3 = corr_feature12 # correlation between presence matrices data1 and data2
    data3 = data3.reindex(feature2_order)  # sorting rows (feature2) in clustering order
    data3 = data3.reindex(columns = feature1_order)  # sorting columns (feature1) in clustering order

    # Create series with cluster colors and sort indexes in clustering order from cluster dataframes (or dictionaries)
    ### for locations based on provided df - cluster_loc_df
    cluster_loc_df = cluster_loc_df.set_index('Locations')
    loc_cluster_color = dict(
        zip(cluster_loc_df['Locations clusters'].unique(),
            sns.color_palette("colorblind", cluster_loc_df['Locations clusters'].nunique())))
    cluster_loc = cluster_loc_df['Locations clusters'].map(loc_cluster_color)
    cl_loc_col = cluster_loc.reindex(loc_order)

    ### for feature1 based on provided df - cluster_feature1_df
    cluster_feature1_df = cluster_feature1_df.set_index('Feature1')
    feature1_cluster_color = dict(
        zip(cluster_feature1_df['Feature1 clusters'].unique(),
            sns.color_palette("Greys", cluster_feature1_df['Feature1 clusters'].nunique())))
    cluster_feature1 = cluster_feature1_df['Feature1 clusters'].map(feature1_cluster_color)
    cl_feature1_col = cluster_feature1.reindex(feature1_order)

    ### for feature2 based on provided df - cluster_feature2_df
    cluster_feature2_df = cluster_feature2_df.set_index('Feature2')
    feature2_cluster_color = dict(
        zip(cluster_feature2_df['Feature2 clusters'].unique(),
            sns.color_palette("Greys", cluster_feature2_df['Feature2 clusters'].nunique())))
    cluster_feature2 = cluster_feature2_df['Feature2 clusters'].map(feature2_cluster_color)
    cl_feature2_col = cluster_feature2.reindex(feature2_order)

    # Create a 2x2 grid of subplots using gridspec
    fig = plt.figure(figsize = (15, 11))
    gs = gridspec.GridSpec(2, 2, figure = fig)
    sns.set_style("white")

    # Create subplot with distribution plot in upper-left free corner (optional)
    ax1 = fig.add_subplot(gs[0, 0]) #set ax location
    distplot = sns.histplot(data = data3, x = data3.max(), bins = 20, ax = ax1)

    # Create the 1st heatmap (feature1 presence) and place it in the upper right subplot
    ax2 = fig.add_subplot(gs[0, 1])
    heatmap1 = sns.heatmap(data1,
                           cmap = sns.color_palette("Blues", as_cmap = True),
                           linewidths = 0.0,
                           xticklabels = False,
                           yticklabels = False,
                           cbar_kws = {"ticks": [0, 1]},
                           ax = ax2)
    heatmap1.tick_params(axis = 'both', which = 'major', pad = 20, length = 0)  # extra padding to leave room for the row and column colors

    ### adding column (x) colors from the cl_feature1_col
    for i, color in enumerate(cl_feature1_col):
        heatmap1.add_patch(plt.Rectangle(xy = (i, -0.05), width = 1, height = 0.05, color = color, lw = 0,
                                   transform = heatmap1.get_xaxis_transform(), clip_on = False))

    ### adding row (y) colors from the cl_loc_col
    for i, color in enumerate(cl_loc_col):
        heatmap1.add_patch(plt.Rectangle(xy = (-0.05, i), width = 0.05, height = 1, color = color, lw = 0,
                                         transform = heatmap1.get_yaxis_transform(), clip_on = False))

    # Create the second heatmap (feature2 presence) and place it in the lower left subplot
    ax3 = fig.add_subplot(gs[1, 0])
    heatmap2 = sns.heatmap(data2,
                           cmap = sns.color_palette("Blues", as_cmap = True),
                           linewidths = 0.0,
                           xticklabels = False,
                           yticklabels = False,
                           cbar_kws = {"ticks": [0, 1]},
                           ax = ax3)

    heatmap2.tick_params(axis = 'both', which = 'major', pad = 20, length = 0)  # extra padding to leave room for the row and column colors

    ### adding row (y) colors from the cl_feature2_col
    for i, color in enumerate(cl_feature2_col):
        heatmap2.add_patch(plt.Rectangle(xy = ( -0.05, i), width = 0.05, height = 1, color = color, lw = 0,
                                         transform = heatmap2.get_yaxis_transform(), clip_on = False))

    ### adding column (x) colors from the cl_loc_col
    for i, color in enumerate(cl_loc_col):
        heatmap2.add_patch(plt.Rectangle(xy = (i, -0.05), width = 1, height = 0.05, color = color, lw = 0,
                                         transform = heatmap2.get_xaxis_transform(), clip_on = False))

    # Create the third heatmap (correlation) and place it in the lower right subplot
    ax4 =  fig.add_subplot(gs[1, 1])
    heatmap3 = sns.heatmap(data3,
                           cmap = "coolwarm",
                           vmin = -1,
                           vmax = 1,
                           center = 0,
                           xticklabels = False,
                           yticklabels = False,
                           cbar_kws = {'label': 'Pearson correlation', "ticks": [-1, 0, 1]},
                           ax = ax4)

    heatmap3.tick_params(axis = 'both', which = 'major', pad = 20, length = 0)  # extra padding to leave room for the row colors

    ### adding column (x) colors from the cl_feature1_col
    for i, color in enumerate(cl_feature1_col):
        heatmap3.add_patch(plt.Rectangle(xy = (i, -0.05), width = 1, height = 0.05, color = color, lw = 0,
                                         transform = heatmap3.get_xaxis_transform(), clip_on = False))

    ### adding row (y) colors from the cl_feature2_col
    for i, color in enumerate(cl_feature2_col):
        heatmap3.add_patch(plt.Rectangle(xy = (-0.05, i), width = 0.05, height = 1, color = color, lw = 0,
                                         transform = heatmap3.get_yaxis_transform(), clip_on = False))


    # Saving the Grid into SVG and PNG files
    svg_name = 'clusters_grid.svg'
    svg_file = f'{visuals}/{svg_name}'
    png_name = 'clusters_grid.png'
    png_file = f'{visuals}/{png_name}'
    if not os.path.isfile(svg_file) and not os.path.isfile(png_file):
        fig.savefig(svg_file, format = 'svg', dpi = gcf().dpi, bbox_inches = 'tight')
        fig.savefig(png_file, format = 'png', dpi = gcf().dpi, bbox_inches = 'tight')

    sns.set_style("white")
    plt.tight_layout()
    plt.show()

grid_image()