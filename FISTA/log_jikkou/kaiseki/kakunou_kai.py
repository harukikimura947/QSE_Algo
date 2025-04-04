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

def kakunou_true_short(R, W, pi, v, prm):

    # 作成するフォルダのパスを指定
    folder_path = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/excelfile/' \
    fr'short/kaisekikai/long_shuseki/mesh/E={prm.E}/v_proj=None/K={prm.Col}^2'

    # フォルダを作成
    create_folder(folder_path)

    file_path_f = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/excelfile/' \
    fr'short/kaisekikai/long_shuseki/mesh/E={prm.E}/v_proj=None/K={prm.Col}^2/jikken_F.xlsx'
    file_path_h = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/excelfile/' \
    fr'short/kaisekikai/long_shuseki/mesh/E={prm.E}/v_proj=None/K={prm.Col}^2/jikken_H.xlsx'

    data_f = {'R': R,
              'W': W,
              'pi': pi}

    data_h = {'v': v.flatten()}

    df_f = pd.DataFrame(data_f)
    df_h = pd.DataFrame(data_h)

    df_f.to_excel(file_path_f, index=False)
    df_h.to_excel(file_path_h, index=False)

    return file_path_f, file_path_h

def kakunou_true_long(m, n, prm):

    # 作成するフォルダのパスを指定
    folder_path = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/excelfile/' \
    fr'shortlong/kaisekikai/theta={prm.theta_firm}/dist=0.15/E={prm.E}/v_proj=0.1/K={prm.Col}^2/truevalue'

    # フォルダを作成
    create_folder(folder_path)

    file_path_f = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/excelfile/' \
    fr'shortlong/kaisekikai/theta={prm.theta_firm}/dist=0.15/E={prm.E}/v_proj=0.1/K={prm.Col}^2/truevalue/jikken_F.xlsx'
    file_path_h = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/excelfile/' \
    fr'shortlong/kaisekikai/theta={prm.theta_firm}/dist=0.15/E={prm.E}/v_proj=0.1/K={prm.Col}^2/truevalue/jikken_H.xlsx'

    data_f = {'m': m}

    data_h = {'n': n.flatten()}

    df_f = pd.DataFrame(data_f)
    df_h = pd.DataFrame(data_h)

    df_f.to_excel(file_path_f, index=False)
    df_h.to_excel(file_path_h, index=False)

    return file_path_f, file_path_h