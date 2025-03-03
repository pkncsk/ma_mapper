# ma_mapper.plots
matplotlib wrapper for plotting 
## list of functions
mapper.plots.plot()
mapper.plots.plot_heatmap()
mapper.plots.plot_bar()
mapper.plots.plot_colorbar()
mapper.plots.plot_annotation()
mapper.plots.plot_basecount()
mapper.plots.plot_logos()

### mapper.plots.plot()

`mapper.plots.plot(data=None, show_alignment=False, alignment=None, alignment_col='nulc_white', heatmap=True, heatmap_mode='overlay', vlim=None, heatmap_color=None, heatmap_title=None, heatmap_title_fs=None, heatmap_yhighlight=None, heatmap_yhighlight_col=['grey'], heatmap_yhighlight_alpha=[0.5], heatmap_xhighlight=None, heatmap_xhighlight_col=['grey'], heatmap_xhighlight_alpha=[0.5], heatmap_ylabel=None, heatmap_ylabel_fs=6, heatmap_xlabel=None, heatmap_xlabel_fs=6, aggregated=False, aggregated_data=None, agg_h=10, agg_ylim=[None], agg_ylabel=None, agg_ylabel_fs=6, agg_ylabel_ypos=None, agg_ylabel_xpos=None, agg_yhighlight=None, agg_yhighlight_col=['grey'], agg_yhighlight_alpha=[0.5], agg_xhighlight=None, agg_xhighlight_col=['grey'], agg_xhighlight_alpha=[0.5], agg_xlabel=None, agg_xlabel_fs=None, agg_major_tick=50, agg_yscale=None, agg_ylabel_right=None, agg_ylabel_right_fs=6, agg_ylabel_right_pos=[1.1, 0.5], agg_plottext=None, agg_plottext_fs=6, agg_plottext_pos=[0.99, 0.90], agg_plot_title=None, agg_plot_title_fs=None, agg_colset=None, logos=False, base_count=False, annotation=False, annotation_data=None, anno_title=None, anno_w=2, anno_ylabel=None, anno_ylabel_fs=6, anno_col=['Blues'], anno_cbar=True, anno_cbar_title=None, anno_cbar_label=None, anno_cbar_even_pos=-0.25, show_alignment_colbar=False, colorbar=False, colorbar_width=2, colorbar_steps=None, image_res=600, figsize=[50, 50], show_plot=True, save_to_file=False, output_filepath=None)`

One of the main APIs of ma_mapper package. This function is the wrapper for matplotlib functions to visualize data matrix extracted from genome-wide dataset using modules in ma_mapper, in a spcific layout, mainly a heatmap of data points with possible annotation and colorbar, accompanied by an aggregated plot under the heatmap.  

#### Parameters
mapper.plots.plot() has many parameters inherited from different functions, so they can be categorized by their functions and prefixes:
**Main Parameters: Heatmap area** 
**data : list[numpy.ndarray], np.ndarray, default: None**
The data matrix or the list of data matrix to be plotted.
**show_alignment : bool, default=False**
Show alignment matrix on the heatmap area.
**alignment : np.ndarray, default: None**
The alignment matrix which is used for data extraction (optional)
**alignment_col : str, default: 'nulc_white'**
The color palletes for alignment matrix visualization
'nulc_white': all bases are white with gaps as grey, used for showing where are the gaps on the heatmap
'dna': ACTG - green, yellow, blue, red, give CG yellow and red as they are hotspots for mutation
'dna_jaspar' ACTG- green, blue, red, yellow, same color pallete as JASPAR motif database
**heatmap : bool, default: True**
Show heatmap area
**heatmap_mode : str, default: 'overlay'**
The layout of heatmap area.
'overlay': all matrices are overlaid on top of each other
'spread_horizontal': the matrices are plotted horizontally
'spread_vertical': the matrices are plotted vertically
**vlim : list[[float, float]], default: None**
The minimum and maximum value to be shown on the heatmap. There can be more than one set to control more than one matrix.
**heatmap_color : matplotlib.colors.ListedColormap,str,list[matplotlib.colors.ListedColormap,str], default: None**
The colormap or the list of colormap for heatmap areas, can use colormaps from matplotlib, or custom built-in colormaps such as: custom_cmap.vlag_mpl, custom_cmap.vlag_r_mpl, custom_cmap.RdBu_s_mpl, custom_cmap.RdBu_sr_mpl, and custom_cmap.LimeGreens_mpl which emulate colormap from the package seaborn
**heatmap_title : list[str], default: None**
The list of the titles of the heatmap area
**heatmap_title_fs : int, default: None**
The font size of the titles of the heatmap area
**heatmap_yhighlight : list[[int,int]], default: None**
The list of coordinates of the vertical highlight on the heatmap area 
**heatmap_yhighlight_col : list[str], default: ['grey']**
The list of colors of the vertical highlight on the heatmap area 
**heatmap_yhighlight_alpha : list[float], default: [0.5]**
The list of alpha values of the vertical highlight on the heatmap area 
**heatmap_xhighlight : list[[int,int]], default: None**
The list of coordinates of the horizontal highlight on the heatmap area 
**heatmap_xhighlight_col : list[str], default: ['grey']**
The list of colors of the horizontal highlight on the heatmap area 
**heatmap_xhighlight_alpha : list[float], default: [0.5]**
The alpha value of the horizontal highlight on the heatmap area 
**heatmap_ylabel : list[str], default: None**
The list of y-axis labels of the heatmap area 
**heatmap_ylabel_fs : int, default: 6**
The font size of the y-axis label of the heatmap area
**heatmap_xlabel : list[str], default: None**
The list of x-axis labels of the heatmap area 
**heatmap_xlabel_fs : int, default: 6**
The font size of the x-axis label of the heatmap area

**Additional Parameters: Aggregated dataplot area** 
**aggregated : bool, default: False**
Show aggregated dataplot area
**aggregated_data : list[numpy.ndarray], default: None**
The aggregated data array or the list of aggregated data array to be plotted.
**agg_h : int, default: 10**
The height of the aggregated dataplot area
**agg_ylim : list[[float, float]], default: [None]**
The minimum and maximum value of y-axis of the aggregated dataplot. There can be more than one set to control more than one array.
**agg_ylabel : list[str], default: None**
The list of y-axis labels of the aggregated dataplot area 
**agg_ylabel_fs : int, default: 6**
The font size of the y-axis label of the aggregated dataplot area 
**agg_ylabel_ypos : list[float], default: None**
The list of y-position of y-axis labels of the aggregated dataplot area 
**agg_ylabel_xpos : list[float], default: None**
The list of x-position of y-axis labels of the aggregated dataplot area 
**agg_yhighlight : list[[int,int]], default: None**
The list of coordinates of the vertical highlight on the aggregated dataplot area  
**agg_yhighlight_col : list[str], default: ['grey']**
The list of colors of the vertical highlight on the aggregated dataplot area
**agg_yhighlight_alpha : list[float], default: [0.5]**
The list of alpha values of the vertical highlight on the aggregated dataplot area 
**agg_xhighlight : list[[int,int]], default: None**
The list of coordinates of the horizontal highlight on the aggregated dataplot area
**agg_xhighlight_col : list[str], default: ['grey']**
The list of colors of the horizontal highlight on the aggregated dataplot area 
**agg_xhighlight_alpha : list[float], default: [0.5]**
The alpha value of the horizontal highlight on the aggregated dataplot area  
**agg_xlabel : str, default: None**
The x-axis label of the aggregated dataplot area
**agg_xlabel_fs : int, default: 6**
The font size of the x-axis label of the aggregated dataplot area
**agg_major_tick : int, default: 50**
The major tick of the x-axis label of the aggregated dataplot area
**agg_yscale : list[str], default: None**
The list of the scale of y-axis of the aggregated dataplot area. See matplotlib.scale
**agg_ylabel_right : list[str], default: None**
The list of the y-lables on the right side of the aggregated dataplot area
**agg_ylabel_right_fs : int, default: 6**
The font size of the y-lables on the right side of the aggregated dataplot area
**agg_ylabel_right_pos : list[int,int], default: [1.1,0.5]**
The coordinates of the y-lables on the right side of the aggregated dataplot area
**agg_plottext : list[str], default: None**
The list of text inside the aggregated dataplot area
**agg_plottext_fs : int, default: 6**
The font size of the text inside the aggregated dataplot area
**agg_plottext_pos : list[int,int], default: [0.99,0.90]**
The coordinates of the text of the aggregated dataplot area
**agg_plot_title : list[str], default: None**
The list of the titles of the aggregated dataplot area
**agg_plot_title_fs : int, default: None**
The font size of the titles of the aggregated dataplot area
**agg_colset : list[str], default: None**
The list of the colors for the plot in the aggregated dataplot area
**logos : bool, default: False**
The additional panel for logos representation of the alignment
**base_count : bool**
The additional panel for base count of the alignment

**Additional Parameters: Annotation area**
**annotation : bool, default: False**
Show annotation area
**annotation_data : list[[str,int]], default: None**
The list of data annotations to be plotted.
**anno_title : list[str], default: None**
The list of the title of the data annotations 
**anno_w : int**
The width of the annotation area
**anno_ylabel : str, default: None**
The y-axis label of the annotation area 
**anno_ylabel_fs : int, default: 6**
The font size of the y-axis label of the annotation area 
**anno_col : list[matplotlib.colors.ListedColormap,str], default: ['Blues']**
The colormap or the list of colormap for data annotations, can use colormaps from matplotlib, or custom built-in colormaps such as: custom_cmap.vlag_mpl, custom_cmap.vlag_r_mpl, custom_cmap.RdBu_s_mpl, custom_cmap.RdBu_sr_mpl, and custom_cmap.LimeGreens_mpl which emulate colormap from the package seaborn
**anno_cbar : bool, default: True**
Show the colorbar of the data annotations
**anno_cbar_title : list[str], default: None**
The list of the titles of the annotaion colorbars
**anno_cbar_label : list[[str]], default: None**
The list of the labels of the annotaion colorbars
**anno_cbar_even_pos : float, default: -0.25**
Adjust annotaiton colobar titles when the number of the annotations is even

**Additional Parameters: Colorbar area**
**show_alignment_colbar : bool, default: False**
Show the colorbar of the multiple alignment
**colorbar : bool, default: False**
Show the colorbars of the heatmap area
**colorbar_width : int, default: 2**
The width of the colorbars
**colorbar_steps : list[float], default: None**
The list of the value step in the colorbars

**Main Parameters: Overall plot**
**image_res : int, default: 600**
The image resolution of the overall plot in DPI
**figsize : list[float, float]**
The size of the figure
**show_plot : bool, default: True**
Show the overall plot
**save_to_file : bool, default: False**
Save the output as PNG file. The filepath can be specified, using save_to_file = True would save the list to the same directory as the coordinate table or current working directory.


