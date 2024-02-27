#%%
import sys
import pandas as pd
import numpy as np
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import mapper
from ma_mapper import fetch_data
#%%
input_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/custom_id.fasta.aligned'
aligned_parsed = mapper.parse_alignment(input_filepath, save_to_file= False)
metadata_aligned = mapper.extract_metadata_from_alignment(input_filepath)
#%%
metadata_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11a_coord_with_id.txt'
bam_mapped_min, bam_mapped_max, bam_mapped_forward, bam_mapped_reverse=fetch_data.fetch_bam(metadata_input= metadata_filepath, bam_input='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/KZFP-bam_hg38/znf808.sorted.bam', custom_id = True) 
# %%
metadata_df = pd.read_csv(metadata_filepath, sep='\t')
original_order = metadata_df.iloc[:,4].unique()
#%%
bam_mapped_sorted_min = []
bam_mapped_sorted_max = []
bam_mapped_sorted_forward = []
bam_mapped_sorted_reverse = []
for idx, row in metadata_aligned.iterrows():
    bam_mapped_sorted_min.append(bam_mapped_min[np.where(original_order == row.id)[0][0]])
    bam_mapped_sorted_max.append(bam_mapped_max[np.where(original_order == row.id)[0][0]])
    bam_mapped_sorted_forward.append(bam_mapped_forward[np.where(original_order == row.id)[0][0]])
    bam_mapped_sorted_reverse.append(bam_mapped_reverse[np.where(original_order == row.id)[0][0]])
#%%
filters=mapper.create_filter(aligned_parsed)
# %%
aligned_bam_overlay_min=mapper.map_data(bam_mapped_sorted_min, aligned_parsed, filters = filters)
aligned_bam_overlay_max=mapper.map_data(bam_mapped_sorted_max, aligned_parsed, filters = filters)
aligned_bam_overlay_forward=mapper.map_data(bam_mapped_sorted_forward, aligned_parsed, filters = filters)
aligned_bam_overlay_reverse=mapper.map_data(bam_mapped_sorted_reverse, aligned_parsed, filters = filters)
# %%
#%%
import plotly.graph_objects as go
subfamily = 'MER11A'
fig = go.Figure()
fig.add_trace(go.Heatmap(z= aligned_bam_overlay_min,
                    colorscale='Greys',
                    #colorbar=dict(tick0=0,dtick=1)
))
fig.update_coloraxes(showscale=False)
fig.update_layout(xaxis_showgrid=False, xaxis_zeroline=False,violingap=0,violinmode='overlay',
    title=subfamily + ' bam overlay (min)',
    xaxis_title="Position (bp)",
    yaxis_title="Sequences",
    legend_title="base",
    showlegend=False,
    )
fig.update_layout(showlegend=False,)
fig.show()
#%%
import plotly.graph_objects as go
subfamily = 'MER11A'
fig = go.Figure()

fig.add_trace(go.Heatmap(z= aligned_bam_overlay_max,
                    colorscale='Greens',
                    #colorbar=dict(tick0=0,dtick=1)
))

fig.update_coloraxes(showscale=False)
fig.update_layout(xaxis_showgrid=False, xaxis_zeroline=False,violingap=0,violinmode='overlay',
    title=subfamily + ' bam overlay (max)',
    xaxis_title="Position (bp)",
    yaxis_title="Sequences",
    legend_title="base",
    showlegend=False,
    )
fig.update_layout(showlegend=False,)
fig.show()
#%%
import plotly.graph_objects as go
subfamily = 'MER11A'
fig = go.Figure()
fig.add_trace(go.Heatmap(z= aligned_bam_overlay_reverse,
                    colorscale='Blues', zauto = False,zmax = 5,
                    #colorbar=dict(tick0=0,dtick=1)
))
fig.update_coloraxes(showscale=False)
fig.update_layout(xaxis_showgrid=False, xaxis_zeroline=False,violingap=0,violinmode='overlay',
    title=subfamily + ' bam overlay (reverse)',
    xaxis_title="Position (bp)",
    yaxis_title="Sequences",
    legend_title="base",
    showlegend=True,
    )
#fig.update_layout(showlegend=False,)
fig.show()
#%%
import plotly.graph_objects as go
subfamily = 'MER11A'
fig = go.Figure()
fig.add_trace(go.Heatmap(z= aligned_bam_overlay_forward,
                    colorscale='Reds', zauto= False,zmax = 5,
                    #colorbar=dict(tick0=0,dtick=1)
))
fig.update_coloraxes(showscale=False)
fig.update_layout(xaxis_showgrid=False, xaxis_zeroline=False,violingap=0,violinmode='overlay',
    title=subfamily + ' bam overlay (forward)',
    xaxis_title="Position (bp)",
    yaxis_title="Sequences",
    legend_title="base",
    showlegend=False,
    )
fig.update_layout(showlegend=False,)
fig.show()
# %%
new_order=np.array([871, 870, 868, 869, 442, 441, 443, 444, 258, 539,  29, 161, 352,
        83, 856, 400,   2, 160, 158, 159, 225, 282, 328, 631, 602, 615,
       199, 366, 374, 146, 361, 544, 334, 370, 750, 512, 530, 616, 805,
       206, 685, 864,   5, 356, 566, 582, 472, 497, 766, 224, 117, 190,
       223, 373, 396, 436, 422, 240, 465, 529, 581, 633, 559, 562, 486,
       690, 579, 271, 806, 398,  40, 195, 488,  43, 455, 703,  77, 425,
       713, 692, 651, 767, 487, 701, 686, 758, 857, 813,   8, 634, 461,
       550, 468, 721, 329, 526, 611, 452, 641, 788, 319, 412, 617, 111,
       153, 773, 567, 397, 789,  72, 482,  68, 365, 815, 182, 430, 304,
       760, 586, 311, 475, 511, 661, 449, 744, 479, 426, 677, 780, 165,
       295, 725, 763,  59, 494,  12, 181, 540, 654, 742, 372,  11, 739,
       263, 474, 735, 232, 606, 659, 797, 275, 730, 821, 683,  62, 687,
       342, 362, 584, 473, 330, 333, 408, 357, 536, 575, 495, 255, 504,
       228, 812, 778, 646, 774, 445, 814, 702, 140, 230,  78, 115, 359,
       162, 421,  17, 167, 568, 287, 322, 343,  23, 553, 202,  90, 499,
        34, 227, 264, 144, 337, 109, 561, 545, 154, 546, 715, 518, 663,
       637, 846,  36, 786, 765, 776, 754, 746,  85, 800,  91, 777, 490,
       748, 681, 819, 516, 719, 498, 707, 528, 192, 737, 522, 257, 679,
        14, 280, 268, 571, 624, 247, 297, 621, 811, 409, 569, 600, 741,
       749, 804, 726, 747, 500, 632, 253, 565,   3, 508, 433, 434, 368,
       519, 694, 447, 448, 842, 841, 867, 355, 722, 851, 662, 790, 699,
       697, 698, 863,  22,  20,  84,  21,  87,   7, 385, 267, 830, 832,
       849, 196, 429, 608, 835, 669, 525, 300, 843,   4, 655, 850, 388,
       828, 234, 844, 827, 187, 642, 853, 369, 383, 823, 390, 186, 872,
       845, 503, 249,   1, 250, 769, 345, 437, 248, 866, 847, 838, 817,
       837, 839, 353, 858, 213, 873, 262, 732, 734, 208, 209, 184, 595,
       831,   6, 680, 310, 245, 859, 424, 798, 459, 799,  56, 296, 458,
       840,  28,  37, 829, 852, 520, 796, 141, 801, 848, 626, 417, 580,
        13, 386, 557, 273, 313, 137, 236, 795, 535, 628, 371, 502, 179,
       547, 510, 794, 636, 756, 354, 623, 733, 274, 649, 696, 285, 215,
       438, 622, 124, 785, 513, 555, 142, 534, 462, 432, 768, 779, 277,
       607, 755, 170, 229, 347, 175, 563, 133, 564, 720, 231, 648, 269,
       664,  19, 291, 860, 393, 652,  42, 309,  57, 855, 299, 279, 404,
       836, 542, 194, 281, 450, 770, 833, 278, 653, 834,  65,  24, 470,
       673, 861, 327, 467,  63, 554,  39, 394, 185, 591, 684, 736, 340,
       338, 339, 610, 341, 577, 237, 244, 326, 446,  81, 414, 506, 507,
       824, 476, 787, 114, 704, 480, 604, 665,  54, 301, 308, 205, 127,
       130, 772, 627, 235, 233, 629,  98,  94,  97,  95,  96, 105, 104,
       103, 101, 102,  99, 100, 106, 107, 120, 121, 358, 321, 363, 573,
       531, 198, 578, 176, 191, 298, 336, 344, 251, 691, 156, 548, 259,
       570, 415, 427, 588,  16, 155, 380, 485, 463, 382, 431, 505, 625,
       453, 727,  74,  64,  71, 515, 246, 413, 164, 320,  51, 183, 123,
       477, 377,  70, 218, 168, 226, 163, 464, 599, 454, 491, 521,  82,
       399, 759, 407, 598, 221, 618, 640, 483, 552, 266, 289, 678,  67,
       147,   0, 294,  89, 594,  18,  26,  60, 576, 220,  92, 169, 302,
       303, 112, 589, 745, 537, 288, 807, 419, 351, 420, 242, 460, 389,
       293, 435, 306, 307, 210, 416, 551, 543, 587, 219, 572, 384, 825,
       527, 315, 171, 324, 865, 122, 592,  50,  66, 596,  55, 583, 597,
       689, 489, 619, 200, 423, 492, 305, 410, 166,  73, 391, 469, 471,
       428, 451, 331, 517, 816, 284, 524, 110, 116, 177, 188, 241, 149,
        61, 283, 484, 129, 148, 272, 381, 605, 125, 203, 378, 201, 131,
       126, 252, 180, 348, 376,  30,  38, 134, 173, 286, 207, 335,  88,
       152,  10, 113, 197, 108, 514, 211, 392, 222, 613, 532, 549, 560,
        33, 364,  25, 556, 612, 292, 418, 118, 509, 323, 558,  52, 349,
       395,  75, 138, 456, 481, 403, 533, 312, 751, 387,  48, 174, 346,
       538, 360, 135, 212, 254, 609,  80, 193, 818,  35, 239, 265, 688,
       260, 593,  27, 541,   9,  47, 379, 132,  46, 214, 238,  69, 204,
       350,  79, 128,  15, 136, 119, 157,  45, 217,  44, 261, 151,  31,
        32,  49, 667, 668, 802, 803, 723, 724, 705, 708, 700, 761, 523,
       466, 728, 644, 695, 714, 614, 411, 645, 325, 601, 666,  76, 762,
       603, 375, 496, 682, 718, 658, 574, 660, 753, 643, 738, 731, 820,
       178, 457, 675, 792, 740, 764, 650, 752, 639, 656, 635, 439, 672,
       676, 638, 671, 808, 809, 781, 793, 775, 810, 791, 862, 711, 712,
       270, 290,  58, 314, 402, 318,  41, 757, 716, 717, 332, 501, 189,
       630, 406, 440, 317, 478, 256, 709, 670, 143, 401, 854, 150, 405,
       826, 493, 782, 620,  93, 216, 783,  86, 172, 585, 822, 276, 243,
       771, 743, 729, 674, 710, 590,  53, 693, 706, 784, 139, 145, 657,
       367, 316, 647], dtype="int32")
# %%
aligned_bam_overlay_forward_sorted = []
for i in new_order:
    aligned_bam_overlay_forward_sorted.append(aligned_bam_overlay_forward[i])
#%%
import plotly.graph_objects as go
subfamily = 'MER11A'
fig = go.Figure()
fig.add_trace(go.Heatmap(z= aligned_bam_overlay_forward_sorted,
                    colorscale='Reds', zauto= False,zmax = 5,
                    #colorbar=dict(tick0=0,dtick=1)
))
fig.update_coloraxes(showscale=False)
fig.update_layout(xaxis_showgrid=False, xaxis_zeroline=False,violingap=0,violinmode='overlay',
    title=subfamily + ' bam overlay (forward)',
    xaxis_title="Position (bp)",
    yaxis_title="Sequences",
    legend_title="base",
    showlegend=False,
    )
fig.update_layout(showlegend=False,)
fig.show()
# %%
