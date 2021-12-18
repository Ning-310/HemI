
import pandas as pd
from scipy.stats import stats


def do_go(gene_set, species, ticks, way):
    go = {}  # 归类
    NN = {}  # 此文件里的全部uniprot id


    f = open('goa_' + species + '.gaf', 'r')
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
    f = open('go.obo', 'r')
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

    dffff = pd.DataFrame(ad,
                         columns=['ID', 'Term', 'Ontology', 'm', 'M', 'm / M', 'n', 'N', 'n / N', 'E-ratio', 'P-value'])
    dffff = dffff.sort_values(by='P-value')
    dffff.to_csv('csvout_' + str(ticks) + '.csv')
    return 0, dffff, ad