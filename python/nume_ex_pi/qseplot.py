#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import numpy as np
import pandas as pd
import scipy.optimize as optimize
import scipy.sparse as spsp
import matplotlib.pyplot as plt
import csv
import os

# フォルダを作成する関数を定義
def create_folder(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def rel_plot(rel_path_list, Col_list, prm, algprm):

    from matplotlib.ticker import ScalarFormatter
    from matplotlib.font_manager import FontProperties
    
    #フォルダの作成
    folder_path = fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/graph/' \
    fr'short/FISTA/long_shuseki/mesh/E={prm.E}/L={algprm.L}/p_proj={algprm.p_proj}/' \
    fr'K={Col_list[0]}^2,{Col_list[1]}^2,{Col_list[2]}^2'
    
    create_folder(folder_path)

    # excelファイルの読み込み
    dF_10 = pd.read_excel(fr'{rel_path_list[0]}')
    dF_20 = pd.read_excel(fr'{rel_path_list[1]}')
    dF_30 = pd.read_excel(fr'{rel_path_list[2]}')

    iteration_10 = np.array(dF_10['iteration'].tolist())
    iteration_20 = np.array(dF_20['iteration'].tolist())
    iteration_30 = np.array(dF_30['iteration'].tolist())
    
    err_R_10 = np.array(dF_10['err_R'].tolist())
    err_R_20 = np.array(dF_20['err_R'].tolist())
    err_R_30 = np.array(dF_30['err_R'].tolist())
    
    err_W_10 = np.array(dF_10['err_W'].tolist())
    err_W_20 = np.array(dF_20['err_W'].tolist())
    err_W_30 = np.array(dF_30['err_W'].tolist())
    
    err_pi_10 = np.array(dF_10['err_pi'].tolist())
    err_pi_20 = np.array(dF_20['err_pi'].tolist())
    err_pi_30 = np.array(dF_30['err_pi'].tolist())
    
    err_v_10 = np.array(dF_10['err_v'].tolist())
    err_v_20 = np.array(dF_20['err_v'].tolist())
    err_v_30 = np.array(dF_30['err_v'].tolist())
    
    err_obj_10 = np.array(dF_10['err_obj'].tolist())
    err_obj_20 = np.array(dF_20['err_obj'].tolist())
    err_obj_30 = np.array(dF_30['err_obj'].tolist())

    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Times New Roman'
    plt.rcParams['mathtext.fontset'] = 'cm'  # 数式フォントもTimes New Romanに設定

    font_path = r'C:\Windows\Fonts\meiryo.ttc'  # フォントのパスを指定
    jp_font = FontProperties(fname=font_path)
    
    # Rのグラフを描画
    plt.plot(iteration_10, err_R_10, '-.', label = f'$K = {Col_list[0]}^2$')
    plt.plot(iteration_20, err_R_20, '--', label = f'$K = {Col_list[1]}^2$')
    plt.plot(iteration_30, err_R_30, ':', label = f'$K = {Col_list[2]}^2$')
    plt.yscale('log')

    plt.legend(fontsize=20)

    # plt.title('m')
    plt.xlabel('短期の反復回数', fontproperties=jp_font, fontsize=20)
    plt.ylabel('E$[(R_i^* - R_i) / R_i^*]$', fontsize=20)

    # plt.yticks([1e-1, 1e-2, 1e-3])
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    plt.xlim(left=1, right=800)
    plt.ylim(top=1, bottom=1e-8)

    plt.grid(True)

    plt.subplots_adjust(left=0.17, bottom=0.15)
    
    plt.savefig(fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/graph/' \
    fr'short/FISTA/long_shuseki/mesh/E={prm.E}/L={algprm.L}/p_proj={algprm.p_proj}/' \
    fr'K={Col_list[0]}^2,{Col_list[1]}^2,{Col_list[2]}^2/R.jpg', format='jpg')
    
    plt.show()
    
    # Wのグラフを描画
    plt.plot(iteration_10, err_W_10, '-.', label = f'$K = {Col_list[0]}^2$')
    plt.plot(iteration_20, err_W_20, '--', label = f'$K = {Col_list[1]}^2$')
    plt.plot(iteration_30, err_W_30, ':', label = f'$K = {Col_list[2]}^2$')
    plt.yscale('log')

    plt.legend(fontsize=20)

    # plt.title('m')
    plt.xlabel('短期の反復回数', fontproperties=jp_font, fontsize=20)
    plt.ylabel('E$[(W_i^* - W_i) / W_i^*]$', fontsize=20)

    # plt.yticks([1e-1, 1e-2, 1e-3])
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    plt.xlim(left=1, right=800)
    plt.ylim(top=1, bottom=1e-8)

    plt.grid(True)

    plt.subplots_adjust(left=0.17, bottom=0.15)
    
    plt.savefig(fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/graph/' \
    fr'short/FISTA/long_shuseki/mesh/E={prm.E}/L={algprm.L}/p_proj={algprm.p_proj}/' \
    fr'K={Col_list[0]}^2,{Col_list[1]}^2,{Col_list[2]}^2/W.jpg', format='jpg')

    plt.show()
    
    # piのグラフを描画
    plt.plot(iteration_10, err_pi_10, '-.', label = f'$K = {Col_list[0]}^2$')
    plt.plot(iteration_20, err_pi_20, '--', label = f'$K = {Col_list[1]}^2$')
    plt.plot(iteration_30, err_pi_30, ':', label = f'$K = {Col_list[2]}^2$')
    plt.yscale('log')

    plt.legend(fontsize=20)

    # plt.title('m')
    plt.xlabel('短期の反復回数', fontproperties=jp_font, fontsize=20)
    plt.ylabel('E$[(\pi_i^* - \pi_i) / \pi_i^*]$', fontsize=20)

    # plt.yticks([1e-1, 1e-2, 1e-3])
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    plt.xlim(left=1, right=800)
    plt.ylim(top=1, bottom=1e-8)

    plt.grid(True)

    plt.subplots_adjust(left=0.17, bottom=0.15)
    
    plt.savefig(fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/graph/' \
    fr'short/FISTA/long_shuseki/mesh/E={prm.E}/L={algprm.L}/p_proj={algprm.p_proj}/' \
    fr'K={Col_list[0]}^2,{Col_list[1]}^2,{Col_list[2]}^2/pi.jpg', format='jpg')

    plt.show()
    
    # vのグラフを描画
    plt.plot(iteration_10, err_v_10, '-.', label = f'$K = {Col_list[0]}^2$')
    plt.plot(iteration_20, err_v_20, '--', label = f'$K = {Col_list[1]}^2$')
    plt.plot(iteration_30, err_v_30, ':', label = f'$K = {Col_list[2]}^2$')
    plt.yscale('log')

    plt.legend(fontsize=20)

    # plt.title('m')
    plt.xlabel('短期の反復回数', fontproperties=jp_font, fontsize=20)
    plt.ylabel('E$[(v_{ij}^* - v_{ij}) / v_{ij}^*]$', fontsize=20)

    # plt.yticks([1e-1, 1e-2, 1e-3])
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    plt.xlim(left=1, right=800)
    plt.ylim(top=1, bottom=1e-8)

    plt.grid(True)

    plt.subplots_adjust(left=0.17, bottom=0.15)
    
    plt.savefig(fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/graph/' \
    fr'short/FISTA/long_shuseki/mesh/E={prm.E}/L={algprm.L}/p_proj={algprm.p_proj}/' \
    fr'K={Col_list[0]}^2,{Col_list[1]}^2,{Col_list[2]}^2/v.jpg', format='jpg')

    plt.show()

    # グラフを描画
    plt.plot(iteration_10, err_obj_10, '-.', label = f'$K = {Col_list[0]}^2$')
    plt.plot(iteration_20, err_obj_20, '--', label = f'$K = {Col_list[1]}^2$')
    plt.plot(iteration_30, err_obj_30, ':', label = f'$K = {Col_list[2]}^2$')
    plt.yscale('log')

    plt.legend(fontsize=20)

    # plt.title('m')
    plt.xlabel('短期の反復回数', fontproperties=jp_font, fontsize=20)
    plt.ylabel('E$[(Z_{SD}^* - Z_{SD}) / Z_{SD}^*]$', fontsize=20)

    # plt.yticks([1e-1, 1e-2, 1e-3])
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    plt.xlim(left=1, right=800)
    plt.ylim(top=1, bottom=1e-8)

    plt.grid(True)

    plt.subplots_adjust(left=0.17, bottom=0.15)
    
    plt.savefig(fr'C:/Users/kimura/OneDrive - 国立大学法人東北大学/numerical/short_dual/log/graph/' \
    fr'short/FISTA/long_shuseki/mesh/E={prm.E}/L={algprm.L}/p_proj={algprm.p_proj}/' \
    fr'K={Col_list[0]}^2,{Col_list[1]}^2,{Col_list[2]}^2/Z.jpg', format='jpg')

    plt.show()

