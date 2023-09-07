import sys
print(sys.path)#查看python安装位置
import os, sys
os.getcwd()
# os.listdir(os.getcwd()) 
#
import scipy.io
import numpy as np
import os
import pandas as pd
import scanpy as sc
#读取输入文件（输入文件为10X的标准脚本）
input_dir =  sys.argv[1]
adata = sc.read_10x_mtx(
    input_dir,  
    var_names='gene_symbols', 
    cache=True)
 
#计算双细胞
sc.external.pp.scrublet(adata, expected_doublet_rate = 0.05, threshold = 0.25)

#输出数据
adata.obs[["doublet_score","predicted_doublet"]].to_csv('./doublet.method1.csv')
