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

def kakunou_true(R, W, pi, v, prm, algprm):

    # 作成するフォルダのパスを指定
    folder_path_f = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/excelfile/' \
    fr'short/FISTA/long_shuseki/mesh/E={prm.E}/L={algprm.L}/p_proj={algprm.p_proj}/K={prm.Col}^2/truevalue'
    folder_path_h = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/excelfile/' \
    fr'short/FISTA/long_shuseki/mesh/E={prm.E}/L={algprm.L}/p_proj={algprm.p_proj}/K={prm.Col}^2/truevalue'

    # フォルダを作成
    create_folder(folder_path_f)
    create_folder(folder_path_h)
    
    file_path_f = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/excelfile/' \
    fr'short/FISTA/long_shuseki/mesh/E={prm.E}/L={algprm.L}/p_proj={algprm.p_proj}/K={prm.Col}^2/truevalue/jikken_F.xlsx'
    file_path_h = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/excelfile/' \
    fr'short/FISTA/long_shuseki/mesh/E={prm.E}/L={algprm.L}/p_proj={algprm.p_proj}/K={prm.Col}^2/truevalue/jikken_H.xlsx'
    
    data_f = {'R': R,
              'W': W,
              'pi': pi}
    
    data_h = {'v': v.flatten()}
    
    df_f = pd.DataFrame(data_f)
    df_h = pd.DataFrame(data_h)
    
    df_f.to_excel(file_path_f, index=False)
    df_h.to_excel(file_path_h, index=False)
    
    return file_path_f, file_path_h

# In[ ]:

def kakunou_rel(iteration, err_R, err_W, err_v, err_pi, err_obj, prm, algprm):
    
    # 作成するフォルダのパスを指定
    folder_path = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/excelfile/' \
    fr'short/FISTA/long_shuseki/mesh/E={prm.E}/L={algprm.L}/p_proj={algprm.p_proj}/K={prm.Col}^2/rel_err'
    
    # フォルダを作成
    create_folder(folder_path)
    
    file_path = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/excelfile/' \
    fr'short/FISTA/long_shuseki/mesh/E={prm.E}/L={algprm.L}/p_proj={algprm.p_proj}/K={prm.Col}^2/rel_err/jikken.xlsx'
    
    data = {'iteration': iteration,
            'err_R': err_R,
            'err_W': err_W,
            'err_pi': err_pi,
            'err_v': err_v,
            'err_obj': err_obj}
    
    # リストからDataFrameを作成
    df = pd.DataFrame(data)

    # DataFrameをExcelファイルに書き込み
    df.to_excel(file_path, index=False)
    
    return file_path
