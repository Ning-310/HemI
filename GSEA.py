import pandas as pd
from scipy.stats import stats

def do_gesa(gene_set, way,ticks):
    f = open('gsea/' + way + '.gmt', 'r')
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
    dffff.to_csv('csvout_' + str(ticks) + '.csv')
    return 0, dffff, ad