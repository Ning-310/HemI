import base64
import xlrd
from io import BytesIO
import time
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.font_manager as fm

import numpy as np
from numpy import mat, shape
from numpy import *
import random


def loadDataSet(file,filename):

    ticks = time.time()
    if str(file).split('.')[-1] == 'txt':
        data = []
        for i in file:
            curLine = i.decode().strip().split('\t')
            lineArr = []
            for j in curLine:
                lineArr.append(j)
            data.append(lineArr)
        Mat = mat(data)
        m, n = shape(Mat)
        display_data = pd.DataFrame(Mat, index=range(0, m), columns=range(0, n))
    elif filename.split('.')[-1] == 'csv':
        display_data = pd.read_csv(file,index_col=0)
        m = display_data.shape[0]
        n = display_data.shape[1]
        index_name = display_data._stat_axis.values.tolist()
    elif filename.split('.')[-1] == 'xls' or filename.split('.')[-1] == 'xlsx':
        display_data = pd.read_excel(file, index_col=0)
        m = display_data.shape[0]
        n = display_data.shape[1]
        index_name = display_data._stat_axis.values.tolist()
    display_data.to_csv('statics/input/' + str(ticks) + '.csv')
    return display_data,ticks,m,n

def heatmap(attribute,ticks):
    my_font = fm.FontProperties(fname="statics/font/arial.ttf")
    df = pd.read_csv('statics/input/'+ticks+'.csv', index_col=0)

    mat = df.iloc[attribute['row_s']:attribute['row_e']+1,attribute['column_s']:attribute['column_e']+1]
    mat.to_csv("statics/output/heatmap_" + ticks + ".csv")
    if (attribute['color'] == 'custom'):
        cmap = sns.light_palette(attribute['custom_color'], as_cmap=True)
    else:
        cmap = attribute['color']
    fig, ax = plt.subplots(figsize=(attribute['width'], attribute['height']))   # 宽 长

    sns.heatmap(mat,xticklabels=attribute['xstep'], yticklabels=attribute['ystep'],fmt = 'f',mask=False, square=attribute['square'], linewidths=attribute['linewidths'], cmap=cmap, annot=attribute['display'], robust=True)
    plt.xticks(fontproperties = my_font,fontsize=attribute['x_fontsize'], rotation=attribute['x_rotation'])
    plt.yticks(fontproperties = my_font,fontsize=attribute['y_fontsize'], rotation=attribute['y_rotation'])

    ax.set_title(attribute['title'], fontsize=6)
    ax.set_ylabel(attribute['ylabel'], fontsize=6)
    ax.set_xlabel(attribute['xlabel'], fontsize=6)
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 设置中文
    sns.set(font='SimHei')
    plt.rcParams['axes.unicode_minus'] = False
    plt.savefig('statics/output_images/heatmap_'+ ticks +'.png', dpi=attribute['dpi'], bbox_inches='tight', transparent=True)
    plt.savefig('statics/output_images/heatmap_'+ ticks +'.'+attribute['img_format'], dpi=attribute['dpi'], bbox_inches='tight', transparent=True)

    plt.close()
    return mat

def clustermap(attribute,ticks):
    my_font = fm.FontProperties(fname="statics/font/arial.ttf")
    df = pd.read_csv('statics/input/' + ticks + '.csv', index_col=0)
    mat = df.iloc[attribute['row_s']:attribute['row_e'] + 1, attribute['column_s']:attribute['column_e'] + 1]

    if (attribute['color'] == 'custom'):
        cmap = sns.light_palette(attribute['custom_color'], as_cmap=True)
    else:
        cmap = attribute['color']
    g = sns.clustermap(mat, method=attribute['method'], metric=attribute['metric'], figsize=(attribute['width'], attribute['height']), linewidths=attribute['linewidths'], cmap=cmap, annot=attribute['display'],
                   row_cluster=attribute['row_cluster'], col_cluster=attribute['col_cluster'], robust=True)
    g.data2d.to_csv("statics/output/cluster_"+ ticks +".csv")
    display_data_cluster = g.data2d  #html表格数据
    plt.tick_params(labelsize=4)
    plt.xticks(fontproperties=my_font, rotation=attribute['x_rotation'])
    plt.yticks(fontproperties=my_font, rotation=attribute['y_rotation'])

    plt.savefig('statics/output_images/cluster_'+ ticks +'.png', dpi=attribute['dpi'], bbox_inches='tight', transparent=True)
    plt.savefig('statics/output_images/cluster_' + ticks + '.' + attribute['img_format'], dpi=attribute['dpi'], bbox_inches='tight',transparent=True)
    plt.close()
    return display_data_cluster

