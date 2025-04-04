#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import os

# In[ ]:
# フォルダを作成する関数を定義
def create_folder(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def kakunou_true(R, W, pi, v, prm, ftol, gtol, method, dic):

    # 作成するフォルダのパスを指定
    folder_path = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
    fr'short/long_ichiyou/{method}/{dic}/RW_proj={prm.RW_proj}/ftol={ftol}/gtol={gtol}/jac=given/E={prm.E}/K={prm.Col}^2'

    # フォルダを作成
    create_folder(folder_path)

    file_path_f = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
    fr'short/long_ichiyou/{method}/{dic}/RW_proj={prm.RW_proj}/ftol={ftol}/gtol={gtol}/jac=given/E={prm.E}/K={prm.Col}^2/jikken_F.xlsx'
    file_path_h = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
    fr'short/long_ichiyou/{method}/{dic}/RW_proj={prm.RW_proj}/ftol={ftol}/gtol={gtol}/jac=given/E={prm.E}/K={prm.Col}^2/jikken_H.xlsx'
    
    data_f = {'R': R,
              'W': W,
              'pi': pi}
    
    data_h = {'v': v.flatten()}
    
    df_f = pd.DataFrame(data_f)
    df_h = pd.DataFrame(data_h)
    
    df_f.to_excel(file_path_f, index=False)
    df_h.to_excel(file_path_h, index=False)
    
    return file_path_f, file_path_h