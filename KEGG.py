import pandas as pd
from scipy.stats import stats

def do_kegg(gene_set, species,ticks):
    abb_dict = {'yeast': 'sce', 'human': 'hsa', 'arabidopsis': 'ath', 'mouse': 'mmu', 'rat': 'rno', 'fly': 'dme', 'pig': 'ssc', 'cow': 'bta', 'dog': 'cfa', 'chicken':'gga','worm':'cel','zebrafish':'dre'}
    abb = abb_dict[species]
    eid = {}
    f = open('genes_uniprot.list', 'r')
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
    f = open('genes_pathway.list', 'r')
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
    f = open('pathway.list', 'r')
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
    f = open('uniprot_' + species +'.fasta', 'r')
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
    dffff.to_csv('csvout_' + str(ticks) + '.csv')
    return 0, dffff, ad