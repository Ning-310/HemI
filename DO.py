import pandas as pd
from scipy.stats import stats


def do_do(gene_set, ticks):

    file = open('do_allgenes.txt', 'r', encoding='utf-8')
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


    f = open('HumanDO.obo', 'r')

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

    dffff.to_csv('csvout_' + str(ticks) + '.csv')


    return 0, dffff, ad