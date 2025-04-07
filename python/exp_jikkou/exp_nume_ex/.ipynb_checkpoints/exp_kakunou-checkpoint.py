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

def short_kakunou_true(R, W, pi, v, prm, algprm, long, method, err_short, dic):

    # 作成するフォルダのパスを指定
    folder_path = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
                  fr'short/long_{long}/{method}/err_short={err_short}/{dic}/RW_proj={prm.RW_proj}/alter_T_num={prm.alter_T_num}/' \
                  fr'L={algprm.L}/eta={algprm.eta}/p_proj={algprm.p_proj}/E={prm.E}/K={prm.Col}^2/truevalue'

    # フォルダを作成
    create_folder(folder_path)

    file_path_f = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
                  fr'short/long_{long}/{method}/err_short={err_short}/{dic}/RW_proj={prm.RW_proj}/alter_T_num={prm.alter_T_num}/' \
                  fr'L={algprm.L}/eta={algprm.eta}/p_proj={algprm.p_proj}/E={prm.E}/K={prm.Col}^2/truevalue/jikken_F.xlsx'
    
    file_path_h = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
                  fr'short/long_{long}/{method}/err_short={err_short}/{dic}/RW_proj={prm.RW_proj}/alter_T_num={prm.alter_T_num}/' \
                  fr'L={algprm.L}/eta={algprm.eta}/p_proj={algprm.p_proj}/E={prm.E}/K={prm.Col}^2/truevalue/jikken_H.xlsx'
    
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

def short_kakunou_rel(iteration, err_R, err_W, err_v, err_pi, err_obj, grad_mean, prm, algprm, long, method, err_short, dic):
    
    # 作成するフォルダのパスを指定
    folder_path = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
                  fr'short/long_{long}/{method}/err_short={err_short}/{dic}/RW_proj={prm.RW_proj}/alter_T_num={prm.alter_T_num}/' \
                  fr'L={algprm.L}/eta={algprm.eta}/p_proj={algprm.p_proj}/E={prm.E}/K={prm.Col}^2/rel_err'
    
    # フォルダを作成
    create_folder(folder_path)
    
    file_path = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
                fr'short/long_{long}/{method}/err_short={err_short}/{dic}/RW_proj={prm.RW_proj}/alter_T_num={prm.alter_T_num}/' \
                fr'L={algprm.L}/eta={algprm.eta}/p_proj={algprm.p_proj}/E={prm.E}/K={prm.Col}^2/rel_err/jikken.xlsx'
    
    data = {'iteration': iteration,
            'err_R': err_R,
            'err_W': err_W,
            'err_pi': err_pi,
            'err_v': err_v,
            'err_obj': err_obj,
            'grad_mean': grad_mean}
    
    # リストからDataFrameを作成
    df = pd.DataFrame(data)

    # DataFrameをExcelファイルに書き込み
    df.to_excel(file_path, index=False)
    
    return file_path

def short_kakunou_rel_obj(obj_rel_list, prm, algprm, long, method, err_short, dic):
    
    # 作成するフォルダのパスを指定
    folder_path = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
                  fr'short/long_{long}/{method}/err_short={err_short}/{dic}/RW_proj={prm.RW_proj}/alter_T_num={prm.alter_T_num}/' \
                  fr'L={algprm.L}/eta={algprm.eta}/p_proj={algprm.p_proj}/E={prm.E}/K={prm.Col}^2/rel'
    
    # フォルダを作成
    create_folder(folder_path)
    
    file_path_f = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
                  fr'short/long_{long}/{method}/err_short={err_short}/{dic}/RW_proj={prm.RW_proj}/alter_T_num={prm.alter_T_num}/' \
                  fr'L={algprm.L}/eta={algprm.eta}/p_proj={algprm.p_proj}/E={prm.E}/K={prm.Col}^2/rel/jikken_F.xlsx'
    
    data_f = {'obj_rel': obj_rel_list}
    
    # リストからDataFrameを作成
    df_f = pd.DataFrame(data_f)

    # DataFrameをExcelファイルに書き込み
    df_f.to_excel(file_path_f, index=False)
    
    return file_path_f

def short_kakunou_nume_num(iteration, prm, algprm, long, method, err_short, dic):
    
    # 作成するフォルダのパスを指定
    folder_path = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
                  fr'short/long_{long}/{method}/err_short={err_short}/{dic}/RW_proj={prm.RW_proj}/alter_T_num={prm.alter_T_num}/' \
                  fr'L={algprm.L}/eta={algprm.eta}/p_proj={algprm.p_proj}/E={prm.E}/K={prm.Col}^2/nume_num'
    
    # フォルダを作成
    create_folder(folder_path)
    
    file_path = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
                fr'short/long_{long}/{method}/err_short={err_short}/{dic}/RW_proj={prm.RW_proj}/alter_T_num={prm.alter_T_num}/' \
                fr'L={algprm.L}/eta={algprm.eta}/p_proj={algprm.p_proj}/E={prm.E}/K={prm.Col}^2/nume_num/jikken.xlsx'
    
    data = {'iteration': iteration}
    
    # リストからDataFrameを作成
    df = pd.DataFrame(data)

    # DataFrameをExcelファイルに書き込み
    df.to_excel(file_path, index=False)
    
    return file_path

def long_kakunou_true(m, n, R, W, pi, v, prm, algprm, method, err_short, err_long, dic):
    
    # 作成するフォルダのパスを指定
    folder_path = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
                  fr'shortlong/{method}/err_short={err_short}/err_long={err_long}/{dic}/RW_proj={prm.RW_proj}/alter_T_num={prm.alter_T_num}/' \
                  fr'L={algprm.L}/eta={algprm.eta}/p_proj={algprm.p_proj}/E={prm.E}/theta={prm.theta_firm}/K={prm.Col}^2/truevalue'
    
    # フォルダを作成
    create_folder(folder_path)
    
    file_path_f = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
                  fr'shortlong/{method}/err_short={err_short}/err_long={err_long}/{dic}/RW_proj={prm.RW_proj}/alter_T_num={prm.alter_T_num}/' \
                  fr'L={algprm.L}/eta={algprm.eta}/p_proj={algprm.p_proj}/E={prm.E}/theta={prm.theta_firm}/K={prm.Col}^2/truevalue/jikken_F.xlsx'
    
    file_path_h = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
                  fr'shortlong/{method}/err_short={err_short}/err_long={err_long}/{dic}/RW_proj={prm.RW_proj}/alter_T_num={prm.alter_T_num}/' \
                  fr'L={algprm.L}/eta={algprm.eta}/p_proj={algprm.p_proj}/E={prm.E}/theta={prm.theta_firm}/K={prm.Col}^2/truevalue/jikken_H.xlsx'
    
    data_f = {'m': m,
              'R': R,
              'W': W,
              'pi': pi}
    
    data_h = {'n': n.flatten(),
              'v': v.flatten()}
    
    # リストからDataFrameを作成
    df_f = pd.DataFrame(data_f)
    df_h = pd.DataFrame(data_h)

    # DataFrameをExcelファイルに書き込み
    df_f.to_excel(file_path_f, index=False)
    df_h.to_excel(file_path_h, index=False)
    
    return file_path_f, file_path_h

def long_kakunou_rel(obj_rel_list, prm, algprm, method, err_short, err_long, dic):
    
    # 作成するフォルダのパスを指定
    folder_path = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
                  fr'shortlong/{method}/err_short={err_short}/err_long={err_long}/{dic}/RW_proj={prm.RW_proj}/alter_T_num={prm.alter_T_num}/' \
                  fr'L={algprm.L}/eta={algprm.eta}/p_proj={algprm.p_proj}/E={prm.E}/theta={prm.theta_firm}/K={prm.Col}^2/rel'
    
    # フォルダを作成
    create_folder(folder_path)
    
    file_path_f = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
                  fr'shortlong/{method}/err_short={err_short}/err_long={err_long}/{dic}/RW_proj={prm.RW_proj}/alter_T_num={prm.alter_T_num}/' \
                  fr'L={algprm.L}/eta={algprm.eta}/p_proj={algprm.p_proj}/E={prm.E}/theta={prm.theta_firm}/K={prm.Col}^2/rel/jikken_F.xlsx'
    
    data_f = {'obj_rel': obj_rel_list}
    
    # リストからDataFrameを作成
    df_f = pd.DataFrame(data_f)

    # DataFrameをExcelファイルに書き込み
    df_f.to_excel(file_path_f, index=False)
    
    return file_path_f

def long_kakunou_nume_num(iteration, prm, algprm, method, err_short, err_long, dic):
    
    # 作成するフォルダのパスを指定
    folder_path = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
                  fr'shortlong/{method}/err_short={err_short}/err_long={err_long}/{dic}/RW_proj={prm.RW_proj}/alter_T_num={prm.alter_T_num}/' \
                  fr'L={algprm.L}/eta={algprm.eta}/p_proj={algprm.p_proj}/E={prm.E}/theta={prm.theta_firm}/K={prm.Col}^2/nume_num'
    
    # フォルダを作成
    create_folder(folder_path)
    
    file_path = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
                  fr'shortlong/{method}/err_short={err_short}/err_long={err_long}/{dic}/RW_proj={prm.RW_proj}/alter_T_num={prm.alter_T_num}/' \
                  fr'L={algprm.L}/eta={algprm.eta}/p_proj={algprm.p_proj}/E={prm.E}/theta={prm.theta_firm}/K={prm.Col}^2/nume_num/nume_num.xlsx'
    
    data = {'iteration': iteration}
    
    # リストからDataFrameを作成
    df = pd.DataFrame(data)

    # DataFrameをExcelファイルに書き込み
    df.to_excel(file_path, index=False)
    
    return file_path

def long_kakunou_val(iteration, Z_list, prm, algprm, method, err_short, err_long, dic):
    
    # 作成するフォルダのパスを指定
    folder_path = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
                  fr'shortlong/{method}/err_short={err_short}/err_long={err_long}/{dic}/RW_proj={prm.RW_proj}/alter_T_num={prm.alter_T_num}/' \
                  fr'L={algprm.L}/eta={algprm.eta}/p_proj={algprm.p_proj}/E={prm.E}/theta={prm.theta_firm}/K={prm.Col}^2/val'
    
    # フォルダを作成
    create_folder(folder_path)
    
    file_path = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/exp/excelfile/' \
                  fr'shortlong/{method}/err_short={err_short}/err_long={err_long}/{dic}/RW_proj={prm.RW_proj}/alter_T_num={prm.alter_T_num}/' \
                  fr'L={algprm.L}/eta={algprm.eta}/p_proj={algprm.p_proj}/E={prm.E}/theta={prm.theta_firm}/K={prm.Col}^2/val/Z.xlsx'
    
    data = {'iteration': iteration,
            'Z': Z_list}
    
    # リストからDataFrameを作成
    df = pd.DataFrame(data)

    # DataFrameをExcelファイルに書き込み
    df.to_excel(file_path, index=False)
    
    return file_path
