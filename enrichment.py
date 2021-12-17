import math
from collections import defaultdict, OrderedDict
import json
from PIL import Image
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from distutils.util import strtobool

from matplotlib import font_manager
from wordcloud import WordCloud, wordcloud
from wordcloud import (WordCloud, get_single_color_func)
import matplotlib.font_manager as fm

class SimpleGroupedColorFunc(object):
    """Create a color function object which assigns EXACT colors
       to certain words based on the color to words mapping

       Parameters
       ----------
       color_to_words : dict(str -> list(str))
         A dictionary that maps a color to the list of words.

       default_color : str
         Color that will be assigned to a word that's not a member
         of any value from color_to_words.
    """

    def __init__(self, color_to_words, default_color):
        self.word_to_color = {word: color
                              for (color, words) in color_to_words.items()
                              for word in words}

        self.default_color = default_color

    def __call__(self, word, **kwargs):
        return self.word_to_color.get(word, self.default_color)


def RGB_to_Hex(a):
    color = '#'
    for i in a:
        num = int(i)
        # 将R、G、B分别转化为16进制拼接转换并大写  hex() 函数用于将10进制整数转换成16进制，以字符串形式表示
        color += str(hex(num))[-2:].replace('x', '0').upper()
    return color
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass

    return False

def gene2uniprot(gene_list, species):
    f = open('statics/uniprot/uniprot_' + species +'.fasta', 'r')
    id_list = []
    for l in f:
        l = l.rstrip()
        if l.startswith('>'):
            for i in gene_list:
                name =' GN=' + i +' '
                if name in l:
                    id_list.append(l.split('|')[1])
    return id_list

def uniprot2gene(uniprot_list, species):
    f = open('statics/uniprot/uniprot_' + species +'.fasta', 'r')
    gene_list = []
    for l in f:
        l = l.rstrip()
        if l.startswith('>'):
            id = l.split('|')[1]  # uniprot id
            for i in uniprot_list:
                if i == id:
                    gene = l.split('GN=')[1]
                    gene_list.append(gene.split(' ')[0])
    return gene_list

def picture(go_set,ad,ticks, attribute):
    my_font = fm.FontProperties(fname="statics/font/arial.ttf")
    eValue = []
    y = []
    mmm = []
    pValue = []
    col_len = len(ad[0])
    if col_len == 9:
        for i in go_set:
            y.append(ad[i][0])
            eValue.append(ad[i][7])
            pValue.append(-math.log(ad[i][8], 10))
            mmm.append(ad[i][1])
    elif col_len == 10:
        for i in go_set:
            y.append(ad[i][1])
            eValue.append(ad[i][8])
            pValue.append(-math.log(ad[i][9], 10))
            mmm.append(ad[i][2])
    elif col_len == 11:
        for i in go_set:
            y.append(ad[i][1])
            eValue.append(ad[i][9])
            pValue.append(-math.log(ad[i][10], 10))
            mmm.append(ad[i][3])
    if attribute['pe_bubblechart'] == 'P1E2':
        x = np.array(eValue)
        color = np.array(pValue)
        plt.xlabel('E-ratio')
    else:
        x = np.array(pValue)
        color = np.array(eValue)
        plt.xlabel('-lg(p-value)')


    mmm = np.array(mmm)
    plt.figure(figsize=(attribute['width_bubblechart'], attribute['height_bubblechart']))
    plt.scatter(x, y, c=color, cmap=attribute['color_bubblechart'], s=mmm * attribute['enlargement_bubblechart'], marker='o')
    plt.colorbar()
    plt.xticks(fontproperties=my_font)
    plt.yticks(fontproperties=my_font)
    ax = plt.gca()  # 获取当前的axes
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['right'].set_color('black')
    # plt.xlim(0, max(x))
    if attribute['pe_bubblechart'] == 'P1E2':

        plt.xlabel('E-ratio')
    else:

        plt.xlabel('-lg($\it{p}$-value)')
    plt.savefig('statics/output_images/bubblechart_' + str(ticks) + '.png', bbox_inches='tight', transparent=True, dpi=attribute['dpi_bubblechart'])
    plt.savefig('statics/output_images/bubblechart_' + str(ticks) + '.jpg', bbox_inches='tight', transparent=True, dpi=attribute['dpi_bubblechart'])
    plt.savefig('statics/output_images/bubblechart_' + str(ticks) + '.pdf', bbox_inches='tight', transparent=True, dpi=attribute['dpi_bubblechart'])
    plt.savefig('statics/output_images/bubblechart_' + str(ticks) + '.tiff', bbox_inches='tight', transparent=True, dpi=attribute['dpi_bubblechart'])
    plt.close()

def Bar(go_set, ad, ticks, attribute_bar):
    eValue = []
    y = []
    mmm = []
    pValue = []
    col_len = len(ad[0])
    if col_len == 9:
        for i in go_set:
            y.append(ad[i][0])
            eValue.append(ad[i][7])
            pValue.append(-math.log(ad[i][8], 10))
            mmm.append(ad[i][1])
    elif col_len == 10:
        for i in go_set:
            y.append(ad[i][1])
            eValue.append(ad[i][8])
            pValue.append(-math.log(ad[i][9], 10))
            mmm.append(ad[i][2])
    elif col_len == 11:
        for i in go_set:
            y.append(ad[i][1])
            eValue.append(ad[i][9])
            pValue.append(-math.log(ad[i][10]))
            mmm.append(ad[i][3])

    if attribute_bar['pe_bar'] == 'P':
        color = np.array(pValue)
    else:
        color = np.array(eValue)
    cm = mpl.cm.get_cmap(attribute_bar['color_bar'])
    colors = cm((color - color.min()) / (color.max() - color.min()))
    mmm = np.array(mmm)
    labels = y
    fig, ax = plt.subplots(figsize=(attribute_bar['width_bar'], attribute_bar['height_bar']), dpi=attribute_bar['dpi_bar'])
    b = ax.barh(np.arange(len(mmm)), mmm, height=attribute_bar['column_width'], tick_label=labels,  color=colors)  # 竖着  横着
    for rect in b:
        w = rect.get_width()
        ax.text(w, rect.get_y() + rect.get_height() / 2, '%d' % int(w), ha='left', va='center')

    # plt.title('Bar', loc='center', fontsize='25',fontweight='bold', color='red')
    my_font = fm.FontProperties(fname="statics/font/arial.ttf")
    sm = plt.cm.ScalarMappable(cmap=attribute_bar['color_bar'],
                               norm=plt.Normalize(vmin=color.min(),
                                                  vmax=color.max()))
    sm._A = []
    plt.colorbar(sm)
    plt.xticks(fontproperties=my_font)
    plt.yticks(fontproperties=my_font)
    ax.set_xlabel('GeneNumber', color='k')
  #  ax.legend()
    axs = plt.gca()  # 获取当前的axes
    axs.spines['bottom'].set_color('black')
    axs.spines['top'].set_color('black')
    axs.spines['left'].set_color('black')
    axs.spines['right'].set_color('black')
    plt.savefig('statics/output_images/bar_' + str(ticks) + '.png', bbox_inches='tight', transparent=True)
    plt.savefig('statics/output_images/bar_' + str(ticks) + '.jpg', bbox_inches='tight', transparent=True)
    plt.savefig('statics/output_images/bar_' + str(ticks) + '.pdf', bbox_inches='tight', transparent=True)
    plt.savefig('statics/output_images/bar_' + str(ticks) + '.tiff', bbox_inches='tight', transparent=True)
    plt.close()

def Rose(go_set, ad, ticks, attribute_rose):
    eValue = []
    name = []
    mmm = []
    pValue = []
    col_len = len(ad[0])
    if col_len == 9:
        for i in go_set:
            name.append(ad[i][0])
            eValue.append(ad[i][7])
            pValue.append(-math.log(ad[i][8], 10))
            mmm.append(ad[i][1])
    elif col_len == 10:
        for i in go_set:
            name.append(ad[i][1])
            eValue.append(ad[i][8])
            pValue.append(-math.log(ad[i][9], 10))
            mmm.append(ad[i][2])
    elif col_len == 11:
        for i in go_set:
            name.append(ad[i][1])
            eValue.append(ad[i][9])
            pValue.append(-math.log(ad[i][10], 10))
            mmm.append(ad[i][3])
    if attribute_rose['pe_rose'] == 'P1E2':
        l = np.array(eValue)
        color = np.array(pValue)
    else:
        color = np.array(eValue)
        l = np.array(pValue)

    cm = mpl.cm.get_cmap(attribute_rose['color_rose'])
    colors = cm((color - color.min()) / (color.max() - color.min()))
    N = len(pValue)
    width = 2 * np.pi / N
    rad = np.cumsum([width] * N) - width / 2
    # 设置 标签颜色、位置
    txt_settings = {
        'span': {0: 0.5, 1: 0.5, 2: 0.5, 3: 0.5},
        'color': {0: 'black', 1: 'black', 2: 'black', 3: 'black'},
        'rot_adj': {0: -90, 1: -90, 2: 90, 3: 90},
        'ha': {0: 'right', 1: 'right', 2: 'left', 3: 'left'}
    }
    # 标签值
    txt_label = [x
                 for x in name]
    plt.figure(figsize=(attribute_rose['width_rose'], attribute_rose['height_rose']), dpi=attribute_rose['dpi_rose'])
    ax = plt.subplot(projection='polar')  # 极坐标图绘制

    ax.set_theta_zero_location('N')  # 设置极坐标的起点（即0度）在正上方向
    ax.grid(False)
    ax.spines['polar'].set_visible(strtobool(attribute_rose['show_circles']))  # 不显示极坐标最外的圆形
    ax.set_yticks([])  # 不显示坐标间隔
    ax.set_thetagrids([])
    # ax.set_ylim(-1,np.ceil(l.max()) + 1)  # 中心为空
    bars = ax.bar(rad, l, width=width, color=colors, alpha=1)

    ax.bar(rad, 1.5, width=width, color='white', alpha=0.2)


    sm = plt.cm.ScalarMappable(cmap=attribute_rose['color_rose'],
                               norm=plt.Normalize(vmin=color.min(),
                                                  vmax=color.max()))
    sm._A = []
    # plt.colorbar(sm,orientation='horizontal', ticklocation='bottom')
    my_font = fm.FontProperties(fname="statics/font/arial.ttf")
    for i in np.arange(N):
        direc = rad[i] // (np.pi / 2)
        ax.text(rad[i],
                l[i] + txt_settings['span'][direc],
                txt_label[i],
                rotation=rad[i] * 180 / np.pi + txt_settings['rot_adj'][direc],
                color=txt_settings['color'][direc],
                ha=txt_settings['ha'][direc], va='center',
                rotation_mode='anchor',  #rotation_mode='default' rotation_mode='anchor'
                alpha=1,
                fontproperties=my_font)


    plt.savefig('statics/output_images/rose_' + str(ticks) + '.png', bbox_inches='tight', transparent=True)
    plt.savefig('statics/output_images/rose_' + str(ticks) + '.jpg', bbox_inches='tight', transparent=True)
    plt.savefig('statics/output_images/rose_' + str(ticks) + '.pdf', bbox_inches='tight', transparent=True)
    plt.savefig('statics/output_images/rose_' + str(ticks) + '.tiff', bbox_inches='tight', transparent=True)
    plt.close()

def Pie(go_set, ad, ticks, attribute_pie):
    my_font = fm.FontProperties(fname="statics/font/arial.ttf")
    eValue = []
    name = []
    mmm = []
    pValue = []
    col_len = len(ad[0])
    if col_len == 9:
        for i in go_set:
            name.append(ad[i][0])
            eValue.append(ad[i][7])
            pValue.append(-math.log(ad[i][8], 10))
            mmm.append(ad[i][1])
    elif col_len == 10:
        for i in go_set:
            name.append(ad[i][1])
            eValue.append(ad[i][8])
            pValue.append(-math.log(ad[i][9], 10))
            mmm.append(ad[i][2])
    elif col_len == 11:
        for i in go_set:
            name.append(ad[i][1])
            eValue.append(ad[i][9])
            pValue.append(-math.log(ad[i][10], 10))
            mmm.append(ad[i][3])

    if attribute_pie['pe_pie'] == 'P1E2':
        color = np.array(pValue)
        w = np.array(eValue)
    else:
        color = np.array(eValue)
        w = np.array(pValue)

    cm = mpl.cm.get_cmap(attribute_pie['color_pie'])
    colors = cm((color - color.min()) / (color.max() - color.min()))
    plt.figure(figsize=(attribute_pie['width_pie'], attribute_pie['height_pie']), dpi=attribute_pie['dpi_pie'])

    sm = plt.cm.ScalarMappable(cmap=attribute_pie['color_pie'],
                               norm=plt.Normalize(vmin=color.min(),
                                                  vmax=color.max()))
    sm._A = []
    plt.colorbar(sm, orientation='horizontal', ticklocation='bottom')

    pie = plt.pie(w,
            labels=name,  # 设置饼图标签
            colors=colors,  # 设置饼图颜色

           # explode=(0, 0.2, 0, 0),  # 第二部分突出显示，值越大，距离中心越远
           # autopct='%.2f%%',  # 格式化输出百分比
            )
    for font in pie[1]:
        font.set_fontproperties(mpl.font_manager.FontProperties(fname='statics/font/arial.ttf'))
    #plt.xlabel(name,fontproperties=my_font)
   # plt.ylabel(name,fontproperties=my_font)

    plt.savefig('statics/output_images/pie_' + str(ticks) + '.png', bbox_inches='tight', transparent=True)
    plt.savefig('statics/output_images/pie_' + str(ticks) + '.jpg', bbox_inches='tight', transparent=True)
    plt.savefig('statics/output_images/pie_' + str(ticks) + '.pdf', bbox_inches='tight', transparent=True)
    plt.savefig('statics/output_images/pie_' + str(ticks) + '.tiff', bbox_inches='tight', transparent=True)

    plt.close()

def WordCloudd(go_set, ad, ticks, attribute_wordcloud):
    eValue = []
    name = []
    mmm = []
    pValue = []
    col_len = len(ad[0])
    frequencies = {}
    if col_len == 9:
        for i in go_set:
            name.append(ad[i][1])
            eValue.append(ad[i][7])
            pValue.append(-math.log(ad[i][8], 10))
            mmm.append(ad[i][1])
    elif col_len == 10:
        for i in go_set:
            name.append(ad[i][1])
            eValue.append(ad[i][8])
            pValue.append(-math.log(ad[i][9], 10))
            mmm.append(ad[i][2])
    elif col_len == 11:
        for i in go_set:
            name.append(ad[i][1])
            eValue.append(ad[i][9])
            pValue.append(-math.log(ad[i][10], 10))
            mmm.append(ad[i][3])
    N = len(np.array(pValue))

    if attribute_wordcloud['pe_wordcloud'] == 'P1E2':
        color = np.array(pValue)
        for i in range(0, N):
            frequencies[name[i]] = eValue[i]
    else:
        color = np.array(eValue)
        for i in range(0, N):
            frequencies[name[i]] = pValue[i]

    cm = mpl.cm.get_cmap(attribute_wordcloud['color_wordcloud'])
    norm = plt.Normalize(color.min(), color.max())
    norm_c = norm(color)
    colors = cm(norm_c)
    colors = colors * 255
    colors = colors[:, 0:3]
    default_color = 'grey'
    color_to_words = {}
    for i in range(0, N):
        R2H = RGB_to_Hex(colors[i])
        z = [name[i]]
        if R2H in color_to_words:
            color_to_words[R2H] = color_to_words[R2H] + [name[i]]
        else:
            color_to_words[R2H] = z

    grouped_color_func = SimpleGroupedColorFunc(color_to_words, default_color)

    wc = WordCloud(font_path='statics/font/arial.ttf', mode='RGBA', min_font_size=1, background_color=None, scale=attribute_wordcloud['scale'], color_func=grouped_color_func).generate_from_frequencies(frequencies)
    plt.imshow(wc, interpolation='bilinear')
    plt.axis('off')
    plt.savefig('statics/output_images/wordcloudImg_' + str(ticks) + '.png', bbox_inches='tight', transparent=True, dpi=attribute_wordcloud['dpi_wordcloud'])
    plt.savefig('statics/output_images/wordcloudImg_' + str(ticks) + '.jpg', bbox_inches='tight', transparent=True, dpi=attribute_wordcloud['dpi_wordcloud'])
    plt.savefig('statics/output_images/wordcloudImg_' + str(ticks) + '.pdf', bbox_inches='tight', transparent=True, dpi=attribute_wordcloud['dpi_wordcloud'])
    plt.savefig('statics/output_images/wordcloudImg_' + str(ticks) + '.tiff', bbox_inches='tight', transparent=True, dpi=attribute_wordcloud['dpi_wordcloud'])
    plt.close()

def do_analysis(gene_set, species, ticks, way):
    go = {}  # 归类
    NN = {}  # 此文件里的全部uniprot id
    goo = defaultdict(list)
    if species == 'human':
        with open('statics/goa/goa_human.json', 'r', encoding='UTF-8') as f:
            go = json.load(f)
        f = open('statics/goa/goa_' + species + '.gaf', 'r')
        for l in f:
            sp = l.rstrip().split('\t')
            NN[sp[1]] = ''
    else:
        f = open('statics/goa/goa_' + species +'.gaf', 'r')
        for l in f:
            sp = l.rstrip().split('\t')
#            if sp[1] not in goo[sp[4]]:
 #               goo[sp[4]].append(sp[1])
            if sp[4] in go:
                go[sp[4]] = go[sp[4]] + [sp[1]]
            else:
                go[sp[4]] = [sp[1]]
            NN[sp[1]] = ''
    """
    test_dict = {
        'version': "1.0",
        'results': goo,
        'explain': {
            'used': True,
            'details': "this is for josn test",
        }
    }
    
    json_str = json.dumps(test_dict)
    with open('goa_human.json', 'w') as json_file:
        json_file.write(json_str)
    """

    goid = {}  # ID：name
    gotp = {}  # ID：namespace
    gid = ''
    name = ''
    tp = ''
    f = open('statics/go.obo', 'r')
    for l in f:
        l = l.rstrip()  # 删除每行后的空格
        if l == '[Term]':
            if gid != '' and name != '' and tp != '':
                goid[gid] = name
                gotp[gid] = tp
            gid = ''
            name = ''
            tp = ''
        if l.startswith('id: GO:'):
            gid = l.split('id: ')[1]
        if l.startswith('name: '):
            name = l.split('name: ')[1][0].upper() + l.split('name: ')[1][1:]
        if l.startswith('namespace: '):
            tp = l.split('namespace: ')[1]

    MM = {}
    gene = {}
    for sp in gene_set:
        if sp in NN:
            MM[sp] = ''
        gene[sp] = ''

    if len(MM) == 0:
        return -1, -1, -1
    ad = []
    for gg in go:  # go
        mm, nn = [], []
        sss = []
        # mm1 = {}
        for gx in MM:
            if gx in go[gg]:
                mm.append(gx)

        go[gg] = set(go[gg])
        m = len(mm)
        M = len(MM)
        n = len(go[gg])
        N = len(NN)
        print('m：{} M：{} n：{} N：{}'.format(m, M, n, N))
        p = stats.fisher_exact([[m, M - m], [n - m, N - n - M + m]])
       # q = -math.log(p[1], 10)
        if gg not in goid.keys() or gg not in gotp.keys():
            continue
        elif mm:
            e = float((m / M) / (n / N))
            if p[1] <= 0.05 and e >= 1:
                if way == 'go':
                    sss.append(gg)
                    sss.append(goid[gg])
                    sss.append(gotp[gg])
                    sss.append(m)
                    sss.append(M)
                    sss.append(m / M)
                    sss.append(n)
                    sss.append(N)
                    sss.append(n / N)
                    sss.append(e)
                    sss.append(float(p[1]))
              #      sss.append(mm)
                    ad.append(sss)
                elif way == 'bp':
                    if gotp[gg] == 'Biological Process':

                        sss.append(gg)
                        sss.append(goid[gg])
                        sss.append(gotp[gg])
                        sss.append(m)
                        sss.append(M)
                        sss.append(m / M)
                        sss.append(n)
                        sss.append(N)
                        sss.append(n / N)
                        sss.append(e)
                        sss.append(float(p[1]))
                 #       sss.append(mm)
                        ad.append(sss)
                elif way == 'mf':
                    if gotp[gg] == 'Molecular Function':
                        sss.append(gg)
                        sss.append(goid[gg])
                        sss.append(gotp[gg])
                        sss.append(m)
                        sss.append(M)
                        sss.append(m / M)
                        sss.append(n)
                        sss.append(N)
                        sss.append(n / N)
                        sss.append(e)
                        sss.append(float(p[1]))
                #        sss.append(mm)
                        ad.append(sss)
                elif way == 'cc':
                    if gotp[gg] == 'Cellular Component':
                        sss.append(gg)
                        sss.append(goid[gg])
                        sss.append(gotp[gg])
                        sss.append(m)
                        sss.append(M)
                        sss.append(m / M)
                        sss.append(n)
                        sss.append(N)
                        sss.append(n / N)
                        sss.append(e)
                        sss.append(float(p[1]))
                 #       sss.append(mm)
                        ad.append(sss)



    dffff = pd.DataFrame(ad, columns=['ID', 'Term', 'Ontology', 'm', 'M', 'm / M', 'n', 'N', 'n / N', 'E-ratio', 'P-value'])
    dffff = dffff.sort_values(by='P-value')
    dffff.to_csv('statics/output/csvout_' + str(ticks) + '.csv')
    return 0, dffff, ad

def do_kegg(gene_set, species,ticks):
    abb_dict = {'yeast': 'sce', 'human': 'hsa', 'arabidopsis': 'ath', 'mouse': 'mmu', 'rat': 'rno', 'fly': 'dme', 'pig': 'ssc', 'cow': 'bta', 'dog': 'cfa', 'chicken':'gga','worm':'cel','zebrafish':'dre'}
    abb = abb_dict[species]
    eid = {}
    f = open('statics/genes_uniprot.list', 'r')
    for l in f:
        if abb +':' not in l:
            continue
        sp = l.rstrip().replace('up:', '').split('\t')
        if sp[0] in eid:
            eid[sp[0]] = eid[sp[0]] + [sp[1]]
        else:
            eid[sp[0]] = [sp[1]]
    go = {}
    gog = {}
    f = open('statics/genes_pathway.list', 'r')
    for l in f:
        sp = l.rstrip().split('\t')
        if 'path:' + abb in l:
            if sp[0] in eid:
                sp[1] = sp[1].replace('path:'+abb, '')
                if sp[1] in go:
                    go[sp[1]] = go[sp[1]] + eid[sp[0]]
                else:
                    go[sp[1]] = eid[sp[0]]
                for k in eid[sp[0]]:
                    gog[k] = ''
    goid = {}
    gotp = {}
    tp = ''
    f = open('statics/pathway.list', 'r')
    for l in f:
        l = l.rstrip()
        if l.startswith('##'):
            tp = l[2:]
        sp = l.split('\t')
        if '#' not in l:
            goid[sp[0]] = sp[1]
            gotp[sp[0]] = tp
    NN = {}
    bgene = {}
    f = open('statics/uniprot/uniprot_' + species +'.fasta', 'r')
    for l in f:
        l = l.rstrip()
        if l.startswith('>'):
            id = l.split('|')[1]
            if id in gog:
                NN[id] = ''
            bgene[id] = ''
    MM = {}
    gene = {}

    for sp in gene_set:
        if sp not in bgene:
            # print(l)
            bgene[sp] = ''
            if sp in gog:
                NN[sp] = ''
        if sp in gog:
            MM[sp] = ''
        gene[sp] = ''
    if len(MM) == 0:
        return -1, -1, -1
    ad = []
    for gg in go:
        mm = []
        nn = []
        sss = []
        mm1 = {}
        for gx in MM:
            if gx in go[gg]:
                mm.append(gx)

        for gx in NN:
            if gx in go[gg]:
                nn.append(gx)
        m = len(mm)
        M = len(MM)
        n = len(nn)
        N = len(NN)
        print('m：{} M：{} n：{} N：{}'.format(m, M, n, N))
        p = stats.fisher_exact([[m, M - m], [n - m, N - n - M + m]])
        #q = -math.log(p[1], 10)
        if n == 0:
            continue
        elif mm:
            e = float((m / M) / (n / N))

            sss.append(abb + gg)
            sss.append(goid[gg])
            sss.append(gotp[gg])
            sss.append(m)
            sss.append(M)
            sss.append(m / M)
            sss.append(n)
            sss.append(N)
            sss.append(n / N)
            sss.append(e)
            sss.append(float(p[1]))
       #     sss.append(mm)
            ad.append(sss)
    dffff = pd.DataFrame(ad, columns=['ID', 'Term', 'Ontology', 'm', 'M', 'm / M', 'n', 'N', 'n / N',
                                          'E-ratio', 'P-value'])
    dffff = dffff.sort_values(by='P-value')
    dffff.to_csv('statics/output/csvout_' + str(ticks) + '.csv')
    return 0, dffff, ad

def do_gesa(gene_set, way,ticks):
    f = open('statics/gsea/' + way + '.gmt', 'r')
    way_list = {}
    NN = set()
    for l in f:
        sp = l.rstrip().split('\t')
        way_list[sp[0]] = [sp[2]]
        NN.add(sp[2])
        for i in range(3, len(sp)):
            way_list[sp[0]] = way_list[sp[0]] + [sp[i]]
            NN.add(sp[i])
    MM = set()
    for l in gene_set:
        if l in NN:
            MM.add(l)
    if len(MM) == 0:
        return -1, -1, -1
    ad = []
    for way in way_list:  # go
        mm, nn = [], []
        sss = []
        for gx in MM:
            if gx in way_list[way]:
                mm.append(gx)
        for gx in NN:
            if gx in way_list[way]:
                nn.append(gx)
        m = mm.__len__()
        M = MM.__len__()
        n = nn.__len__()
        N = NN.__len__()
        if m == 0:
            continue
        p = stats.fisher_exact([[m, M - m], [n - m, N - n - M + m]])
        #q = -math.log(p[1], 10)
        #    print(p)
        print('m：{} M：{} n：{} N：{}'.format(m, M, n, N))
        e = float((m / M) / (n / N))
        if p[1] <= 0.05 and e >= 1:
            sss.append(way)
            sss.append(m)
            sss.append(M)
            sss.append(m / M)
            sss.append(n)
            sss.append(N)
            sss.append(n / N)
            sss.append(e)
            sss.append(float(p[1]))
    #        sss.append(mm)
            ad.append(sss)
    dffff = pd.DataFrame(ad, columns=['ID', 'm', 'M', 'm / M', 'n', 'N', 'n / N',
                                      'E-ratio', 'P-value'])
    dffff = dffff.sort_values(by='P-value')
    dffff.to_csv('statics/output/csvout_' + str(ticks) + '.csv')
    return 0, dffff, ad

def do_do(gene_set, ticks):

    file = open('statics/do_allgenes.txt', 'r', encoding='utf-8')
    DO = {}
    NN = set()
    for i in file:
        curLine = i.strip().split('\t')
        genes = curLine[8].strip().split('/')
        DO[curLine[0].split('DOID:')[1]] = genes
        for g in genes:
            NN.add(g)

    Name = {}  # ID：name
    gid = ''
    name = ''


    f = open('statics/HumanDO.obo', 'r')

    for l in f:
        l = l.rstrip()  # 删除每行后的空格

        if l == '[Term]':
            if gid != '' and name != '':
                Name[gid] = name

            gid = ''
            name = ''
        if l.startswith('id: DOID:'):
            gid = l.split('id: DOID:')[1]

        if l.startswith('name: '):
            name = l.split('name: ')[1]


    MM = {}
    for sp in gene_set:
        if sp in NN:
            MM[sp] = ''

    if len(MM) == 0:
       return -1, -1, -1

    ad = []
    for dd in DO:
        mm, nn = [], []
        sss = []
        # mm1 = {}
        for gx in MM:
            if gx in DO[dd]:
                mm.append(gx)

        DO[dd] = set(DO[dd])
        m = len(mm)
        M = len(MM)
        n = len(DO[dd])
        N = 5063
        print('m：{} M：{} n：{} N：{}'.format(m, M, n, N))
        p = stats.fisher_exact([[m, M - m], [n - m, N - n - M + m]])
        if mm:
            e = float((m / M) / (n / N))
            sss.append(dd)
            sss.append(Name[dd])
            sss.append(m)
            sss.append(M)
            sss.append(m / M)
            sss.append(n)
            sss.append(N)
            sss.append(n / N)
            sss.append(e)
            sss.append(float(p[1]))
            #      sss.append(mm)
            ad.append(sss)
    dffff = pd.DataFrame(ad, columns=['DOID', 'Disease Name', 'm', 'M', 'm / M', 'n', 'N', 'n / N',
                                      'E-ratio', 'P-value'])
    dffff = dffff.sort_values(by='P-value')

    dffff.to_csv('statics/output/csvout_' + str(ticks) + '.csv')


    return 0, dffff, ad

