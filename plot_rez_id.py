# -*- coding: utf-8 -*-
"""
# ===========================================================================
# ===========================================================================
# !==   Safwan ALJBAAE                                                     ==
# !==   November 2020                                                      ==
# ===========================================================================
"""

import sys
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

#i_folder = sys.argv[1]
i_folder = 1


with open('./ast_list') as f:
    id = f.read().splitlines()
    
res_files = glob.glob("res_arg_*")
n_tp = len(res_files)
print(n_tp)
for i in range(1, n_tp):
    f = "res_arg_" + str("{:02d}".format(i))
    data_res = pd.read_csv(f,
                           skiprows=0,
                           header=None,
                           delim_whitespace=True,
                           index_col=None,
                           names=['t', 'res_ang_cos', 'res_ang_sin', 'res_ang'],
                           low_memory=False,
                           dtype={'t': np.float64,
                                  'res_ang_cos': np.float64,
                                  'res_ang_sin': np.float64,
                                  'res_ang': np.float64
                                  }
                           )
    order = n_tp * (int(i_folder) - 1) + i
    order = int(id[order - 1])
    fig = plt.figure(figsize=(1,1))
    ax = fig.add_subplot(111)
    #plt.subplots_adjust(wspace=1.0, hspace=1.0, left=0.15, right=0.95, bottom=0.2, top=0.9)
    ax.scatter(data_res.t, data_res.res_ang, c='black', marker='.', s=0.5)
    #ax.set_title(str("{:07d}".format(order)), fontsize=5)
    #ax.set_xlabel('$Time$ (MY)', fontsize=5)
    #ax.set_ylabel('$\sigma$ ($^{\circ}$)', fontsize=5)
    ax.axis('off')
    ax.set_xlim(0, 0.11)
    ax.set_ylim(0, 360)
  
    image = 'fig_res_' + str("{:07d}".format(order)) + '.png'

    fig.savefig(image, format='png', dpi=100)
    plt.close(fig)
    
    img = Image.open(image).convert('RGB')
    width, height = img.size
    print(f'The size of the image {image} is: {width} X {height} pixels. ')

