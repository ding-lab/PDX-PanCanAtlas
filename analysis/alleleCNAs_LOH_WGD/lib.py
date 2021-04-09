#!/usr/bin/env python2.7

# Author: Dr Simone Zaccaria
# Old affilliation: Princeton University, NJ (USA)
# New affilliation: UCL Cancer Institute, London (UK)
# Correspondence: s.zaccaria@ucl.ac.uk


import os
import sys
import glob
import re
from os.path import *
from collections import defaultdict
from collections import Counter
from collections import OrderedDict
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
import seaborn as sns
from itertools import cycle
from itertools import combinations
from scipy.stats import chi2_contingency
plt.style.use('ggplot')
sns.set_style("whitegrid")
plt.rcParams["axes.grid"] = True
plt.rcParams["axes.edgecolor"] = "k"
plt.rcParams["axes.linewidth"]  = 1.5



wgd = defaultdict(lambda : None)
nclones = defaultdict(lambda : None)
pur = defaultdict(lambda : defaultdict(lambda : None))
cns = defaultdict(lambda : defaultdict(lambda : dict()))
getgenome = (lambda p : (p[0], int(p[1]), int(p[2])))
getsample = (lambda p : p[3])
getpur = (lambda p : 1.0 - float(p[5]))
states = (lambda v : tuple(map(int, v.split('|'))))
getcns = (lambda p : filter(lambda v : v[1] > 0.0, zip(map(states, p[::2]), map(float, p[1:][::2]))))
form = (lambda p : (getgenome(p), getsample(p), getpur(p), getcns(p[6:])))
for din in glob.glob('hatchet/*'):
    pat = basename(din)
    with open(join(din, 'hatchet.seg.ucn'), 'r') as i:
        for l in (l for l in i if l[0] != '#'):
            g, s, u, c = form(l.strip().split())
            assert pur[pat][s] is None or abs(pur[pat][s] - u) <= 0.01
            pur[pat][s] = u
            cns[pat][g][s] = c
    with open(join(din, 'hatchet.wgd.txt'), 'r') as i:
        for l in i:
            search = re.search("The chosen solution is (.*) with ([1-9]*) clones and is written in ./best.bbc.ucn and ./best.seg.ucn", l)
            if search is not None:
                assert wgd[pat] is None and nclones[pat] is None
                wgd[pat] = 'tetraploid' == search.group(1)
                nclones[pat] = int(search.group(2)) - 1
        assert all(wgd[pat] is not None for pat in pur)
        assert all(nclones[pat] is not None for pat in pur)
wgd = dict(wgd)
nclones = dict(nclones)
pur = {pat : dict(pur[pat]) for pat in pur}
cns = {pat : {g : dict(cns[pat][g]) for g in cns[pat]} for pat in cns}

details = pd.read_csv('sampleinfo.tsv', sep='\t')
ctype = [(r['Formal_CaseID'], r['CancerType']) for i, r in details.iterrows()]
ctype = {p : set(map(lambda y : y[1], filter(lambda x : x[0] == p, ctype))) for p in wgd}
assert all(len(ctype[p]) == 1 for p in wgd)
ctype = {p : list(ctype[p])[0] for p in ctype}
ctypecolos = {}
ctypecolos["BLCA"] = "#FAD2D9"
ctypecolos["BRCA"] = "#ED2891"
ctypecolos["CESC"] = "#F6B667"
ctypecolos["COAD"] = "#00A99D"
ctypecolos["COADREAD"] = "#187251"
ctypecolos["CSCC"] = "#fefbd8"
ctypecolos["GBM"] = "#B2509E"
ctypecolos["GIAD"] = "#ccccff"
ctypecolos["GIST"] = "#00E6E6"
ctypecolos["HNSC"] = "#97D1A9"
ctypecolos["KIRC"] = "#F8AFB3"
ctypecolos["LGG"] = "#754C29"
ctypecolos["LIHC"] = "#CACCDB"
ctypecolos["LUAD"] = "#de3bed"
ctypecolos["LUAS"] = "#C1A72F"
ctypecolos["LUSC"] = "#A084BD"
ctypecolos["MESO"] = "#006298"
ctypecolos["MCC"] = "#C8FFBD"
ctypecolos["OV"] = "#D97D25"
ctypecolos["PRAD"] = "#4597ad"
ctypecolos["PAAD"] = "#F4C91F"
ctypecolos["READ"] = "#009444"
ctypecolos["SARC"] = "#7A439E"
ctypecolos["SCLC"] = "#BE1E2D"
ctypecolos["SKCM"] = "#BBD642"
ctypecolos["STAD"] = "#6E7BA2"
ctypecolos["TGCT"] = "#3953A4"
ctypecolos["THCA"] = "#F9ED32"
ctypecolos["THYM"] = "#CEAC8F"
ctypecolos["UCEC"] = "#DAF1FC"
ctypecolos["UCS"] = "#FBE3C7"
ctypecolos["Other"] = "#cccccc"
orig = [(r['Formal_CaseID'], r['Analysis_ID']) for i, r in details[details['Group']=='Human_Tumor'].iterrows()]
orig = {p : set(map(lambda y : y[1], filter(lambda x : x[0] == p, orig))) for p in wgd}

orderchrs = (lambda x : int(''.join([l for l in x if l.isdigit()])))
grid = defaultdict(lambda : dict())
getgenome = (lambda p : (p[0], int(p[1]), int(p[2])))
getsample = (lambda p : p[3])
getpur = (lambda p : 1.0 - float(p[12]))
states = (lambda v : tuple(map(int, v.split('|'))))
getcns = (lambda p : filter(lambda v : v[1] > 0.0, zip(map(states, p[::2]), map(float, p[1:][::2]))))
form = (lambda p : (getgenome(p), getsample(p), getpur(p), getcns(p[13:])))
for d in glob.glob('hatchet/*'):
    pat = os.path.basename(d)
    with open(os.path.join(d, 'hatchet.bbc.ucn'), 'r') as i:
        for l in (l for l in i if l[0] != '#'):
            g, s, u, C = form(l.strip().split())
            grid[g][(pat, s)] = (u, C)

pos = sorted(grid.keys(), key=(lambda x : (orderchrs(x[0]), x[1], x[2])))
chrs = set(g[0] for g in grid)
mkseg = (lambda L : zip(L[:-1], L[1:]))
segs = {c : mkseg(range(0, max(g[2] for g in grid if g[0]==c), 500000)) for c in chrs}
overlap = (lambda a, b : min(a[1], b[1]) - max(a[0], b[0]))
segs = {c : {s : filter(lambda g : g[0] == c and overlap(g[1:], s) > 0, grid.keys()) for s in segs[c]} for c in chrs}
segs = {c : {s : segs[c][s] for s in segs[c] if len(segs[c][s]) > 0} for c in segs}
cover = (lambda c, s : len(set(s for g in segs[c][s] for s in grid[g])) > 250)
segs = {c : {s : segs[c][s] for s in segs[c] if cover(c, s)} for c in segs}
samples = set(s for g in grid for s in grid[g])
nclones = {s : {len(grid[g][s][1]) for g in grid if s in grid[g]} for s in samples}
assert all(len(nclones[s]) == 1 for s in nclones)
nclones = {s : list(nclones[s])[0] for s in nclones}
mkdef = (lambda s : tuple([((2, 2) if wgd[s[0]] else (1, 1), 1.0) for i in xrange(nclones[s])]))
argmax = (lambda D, s : max(D.keys(), key=(lambda x : D[x])) if len(D) > 0 else mkdef(s))
select = (lambda c, e, s : argmax(Counter([tuple(grid[g][s][1]) for g in segs[c][e] if s in grid[g]]), s))
cngrid = {(c,) + e : {s : select(c, e, s) for s in samples} for c in segs for e in segs[c]}
nclones = {s : {len(cngrid[g][s]) for g in cngrid} for s in samples}
assert all(len(nclones[s]) == 1 for s in nclones)


def draw_grid():
    pos = sorted(cngrid.keys(), key=(lambda x : (orderchrs(x[0]), x[1], x[2])))
    avail = [(t - i, i) for t in xrange(7) for i in reversed(xrange(t+1)) if i <= t - i]
    order = (lambda p : (max(p), min(p)))
    convert = (lambda p : order(p) if sum(p) <= 6 else min(avail, key=(lambda x : abs(p[0] - x[0]) + abs(p[1] - x[1]))))
    df = []
    mapc = {}
    found = set()
    c0 = cngrid[pos[20]]
    sel = {e : {x : c[1] >= max(c[1] for c in c0[e]) for x, c in enumerate(c0[e])} for e in samples}
    for o, b in enumerate(pos):
        for x, j in enumerate((e + (x,), c) for e in sorted(samples, key=(lambda s : (wgd[s[0]], ctype[s[0]], s[0], [x for x in sel[s] if sel[s][x]][0], s[1]))) for x, c in enumerate(cngrid[b][e]) if sel[e][x]):
            e, c = j
            df.extend([{'Cell' : x, 'Genome' : o, 'Value' : convert(c[0])}])
            mapc[x] = e
    df = pd.DataFrame(df)
    found = [v for v in avail if v in set(df['Value'])]
    smap = {v : x for x, v in enumerate(found)}
    df['CN states'] = df.apply(lambda r : smap[r['Value']], axis=1)
    table = pd.pivot_table(df, values='CN states', columns=['Genome'], index=['Cell'], aggfunc='first')
    title = 'Copy-number states'
    palette = {}
    palette.update({(0, 0) : 'darkblue'})
    palette.update({(1, 0) : 'lightblue'})
    palette.update({(1, 1) : 'lightgray', (2, 0) : 'dimgray'})
    palette.update({(2, 1) : 'lightgoldenrodyellow', (3, 0) : 'gold'})
    palette.update({(2, 2) : 'navajowhite', (3, 1) : 'orange', (4, 0) : 'darkorange'})
    palette.update({(3, 2) : 'salmon', (4, 1) : 'red', (5, 0) : 'darkred'})
    palette.update({(3, 3) : 'plum', (4, 2) : 'orchid', (5, 1) : 'purple', (6, 0) : 'indigo'})
    colors = [palette[c] for c in found]
    cmap = LinearSegmentedColormap.from_list('multi-level', colors, len(colors))

    chr_palette = cycle(['#525252', '#969696', '#cccccc'])
    chr_colors = {c : next(chr_palette) for c in sorted(set(b[0] for b in cngrid), key=orderchrs)}
    seen = set()
    seen_add = seen.add
    ordclones = [ctype[mapc[x][0]] for x in table.index if not (ctype[mapc[x][0]] in seen or seen_add(ctype[mapc[x][0]]))]
    cell_palette = cycle(sns.color_palette("tab20", len(set(ordclones))))
    clone_colors = {i : next(cell_palette) for i in ordclones}
    bw = cycle(['#000000', '#ffffff'])
    pdx_color = {x : next(bw) for x in list(OrderedDict.fromkeys([mapc[s][0] for s in table.index]))}

    para = {}
    para['data'] = table
    para['cmap'] = cmap
    para['yticklabels'] = False
    para['row_cluster'] = False
    para['method'] = 'single'
    para['metric'] = 'hamming'
    para['xticklabels'] = False
    para['col_cluster'] = False
    para['figsize'] = (12, 6)
    para['rasterized'] = True
    para['col_colors'] = pd.DataFrame([{'index' : s, 'chromosomes' : chr_colors[pos[x][0]]} for x, s in enumerate(table.columns)]).set_index('index')
    para['row_colors'] = pd.DataFrame([{'index' : x, 'Clone' : ctypecolos[ctype[mapc[x][0]]], 'Tumor' : pdx_color[mapc[x][0]]} for x in table.index]).set_index('index')
    g = sns.clustermap(**para)
    addchr(g, pos)
    legend_TN = [mpatches.Patch(color=ctypecolos[t], label=t) for t in sorted(set(ctype.values()))]
    l2=g.ax_heatmap.legend(loc='center left',bbox_to_anchor=(1.31,0.85),handles=legend_TN,frameon=True)
    l2.set_title(title='tissue type',prop={'size':10})
    plt.savefig("grid.png", bbox_inches = 'tight', dpi=900)
    

def draw_table():
    overlap = (lambda a, b : max(0, min(a[1], b[1]) - max(a[0], b[0])))
    tp53 = ('chr17', 7572926, 7687550)
    covertp53 = {p : {g for g in cns[p] if g[0] == 'chr17' and overlap(tp53[1:], g[1:]) > 0} for p in cns}
    assert all(len(covertp53[p]) > 0 for p in covertp53)
    covertp53 = {p : list(covertp53[p])[0] for p in covertp53}
    frac = (lambda C : sum(sum(c[0]) * c[1] for c in C))
    tp53 = {p : {s : frac(cns[p][covertp53[p]][s]) for s in cns[p][covertp53[p]]} for p in cns}
    
    def value(p, s, f):
        if f == 'WGD':
            return wgd[p]
        elif f == 'TP53 LOH':
            return any(0 in c[0] for c in cns[p][covertp53[p]][s])
        elif f == 'Intra-sample subclonality':
            return len({c[0] for c in cns[p][cns[p].keys()[0]][s]}) > 1
        elif f == 'Inter-sample subclonality':
            return any(any(len({c[0] for c in cns[p][g][s]} - {c[0] for c in cns[p][g][O]}) > 0 for g in cns[p]) for O in pur[p])
        elif f == 'Deletion abudance':
            if wgd[p]:
                return (sum((c[0][0] < 2 or c[0][1] < 2) * (g[2] - g[1]) for g in cns[p] for c in cns[p][g][s] if set(c[0]) != {2}) / float(sum(g[2] - g[1] for g in cns[p] for c in cns[p][g][s] if set(c[0]) != {2}))) >= 0.5
            else:
                return (sum((c[0][0] < 1 or c[0][1] < 1) * (g[2] - g[1]) for g in cns[p] for c in cns[p][g][s] if set(c[0]) != {1}) / float(sum(g[2] - g[1] for g in cns[p] for c in cns[p][g][s] if set(c[0]) != {1}))) >= 0.5
        else:
            assert False
    df = pd.DataFrame([{'Sample' : '{}:{}:{}'.format(ctype[p], p, s), 'Feature' : '{}-{}'.format(x, f), 'Value' : value(p, s, f)} for p in pur for s in pur[p] for x, f in enumerate(['WGD', 'Deletion abudance', 'TP53 LOH', 'Intra-sample subclonality', 'Inter-sample subclonality'])])
    cont = pd.DataFrame([dict([('Sample' , '{}:{}:{}'.format(ctype[p], p, s))] + [(f, value(p, s, f)) for f in ['WGD', 'Deletion abudance', 'TP53 LOH', 'Intra-sample subclonality', 'Inter-sample subclonality']]) for p in pur for s in pur[p]])
    opts = list(combinations(['WGD', 'Deletion abudance', 'TP53 LOH', 'Intra-sample subclonality', 'Inter-sample subclonality'], 2))
    count = (lambda f, v1, v2 : len(cont[(cont[f[0]]==v1)&((cont[f[1]]==v2))]))
    pvalue = (lambda f : (f, chi2_contingency([[count(f, False, False), count(f, False, True)], [count(f, True, False), count(f, True, True)]])[1]))
    map(pvalue, opts)

    table = pd.pivot_table(df, index='Sample', columns='Feature', values=['Value'])
    table = table.reindex(table['Value'].sort_values(by=['0-WGD', '1-Deletion abudance', '2-TP53 LOH', '3-Intra-sample subclonality', '4-Inter-sample subclonality'], ascending=True).index)
    para = {}
    para['data'] = table
    para['cmap'] = 'Greys'
    para['yticklabels'] = False
    para['row_cluster'] = False
    para['xticklabels'] = True
    para['col_cluster'] = False
    para['figsize'] = (4, 12)
    para['rasterized'] = True
    para['row_colors'] = pd.DataFrame([{'index' : x, 'Cancer type' : ctypecolos[ctype[x.split(':')[1]]]} for x in table.index]).set_index('index')
    g = sns.clustermap(**para)
    legend_TN = [mpatches.Patch(color=ctypecolos[t], label=t) for t in sorted(set(ctype.values()))]
    l2=g.ax_heatmap.legend(loc='center left',bbox_to_anchor=(1.31,0.85),handles=legend_TN,frameon=True)
    l2.set_title(title='tissue type',prop={'size':10})
    plt.savefig("table.png", bbox_inches = 'tight', dpi=900)


def draw_nowgd():
    df = pd.read_csv('hatchet/PDMR-981375/hatchet.bbc.ucn', sep='\t')
    df['0.5 - BAF'] = np.abs(0.5 - df['BAF'])
    count = Counter(df['CLUSTER'])
    sel = sorted(count.keys(), key=(lambda x : count[x]), reverse=True)[:7]
    df['Selected'] = df['CLUSTER'].isin(set(sel) - {16}) & df['SAMPLE'].isin(['981375_231_R_T92TV2_v1_2_WES'])
    para = {}
    para['data'] = df[df['Selected']==True]
    para['x'] = '0.5 - BAF'
    para['y'] = 'RD'
    para['hue'] = 'CLUSTER'
    para['fit_reg'] = False
    para['palette'] = ['lightgray', 'khaki', 'dimgray', 'lightblue', 'orange', 'gold']
    para['legend'] = False
    para['scatter_kws'] = {'s' : 7}
    para['size'] = 4
    para['aspect'] = 1.1
    g = sns.lmplot(**para)
    g.despine(top=False, bottom=False, left=False, right=False)
    g.set(xlim=(0.0, 0.52), ylim=(0.35, 2.0))
    plt.savefig("nowgd.png", bbox_inches = 'tight', dpi=900)


def draw_wgd():
    df = pd.read_csv('hatchet/PDMR-714841/hatchet.bbc.ucn', sep='\t')
    df['0.5 - BAF'] = np.abs(0.5 - df['BAF'])
    count = Counter(df['CLUSTER'])
    sel = sorted(count.keys(), key=(lambda x : count[x]), reverse=True)[:14]
    df['Selected'] = df['CLUSTER'].isin(set(sel) - {34, 32, 33, 2}) & df['SAMPLE'].isin(['714841_288_R_KV6MF7NX6M067_v1_2_WES'])
    para = {}
    para['data'] = df[df['Selected']==True]
    para['x'] = '0.5 - BAF'
    para['y'] = 'RD'
    para['hue'] = 'CLUSTER'
    para['fit_reg'] = False
    para['palette'] = ['red', 'khaki', 'salmon', 'red', 'dimgray', 'orange', 'gold', 'orchid', 'orchid', 'navajowhite', 'salmon']
    para['legend'] = False
    para['scatter_kws'] = {'s' : 7}
    para['size'] = 4
    para['aspect'] = 1.1
    g = sns.lmplot(**para)
    g.despine(top=False, bottom=False, left=False, right=False)
    g.set(xlim=(0.0, 0.52), ylim=(0.35, 2.0))
    plt.savefig("wgd.png", bbox_inches = 'tight', dpi=900)


def draw_multipassage():
    PAT = 'PDMR-519858'
    pos = filter(lambda g : g[0] in {'chr2', 'chr3', 'chr6', 'chr7', 'chr11', 'chr12'}, sorted(cngrid.keys(), key=(lambda x : (orderchrs(x[0]), x[1], x[2]))))
    avail = [(t - i, i) for t in xrange(7) for i in reversed(xrange(t+1)) if i <= t - i]
    order = (lambda p : (max(p), min(p)))
    convert = (lambda p : order(p) if sum(p) <= 6 else min(avail, key=(lambda x : abs(p[0] - x[0]) + abs(p[1] - x[1]))))
    df = []
    mapc = {}
    found = set()
    c0 = cngrid[pos[20]]
    sel = {e : {x : c[1] >= .22 for x, c in enumerate(c0[e])} for e in samples}
    for o, b in enumerate(pos):
        for z, j in enumerate((e + (x,), c) for e in sorted(samples, key=(lambda s : (s[0], s[1]))) if e[0] == PAT for x, c in enumerate(cngrid[b][e]) if sel[e][x]):
            e, c = j
            df.extend([{'Cell' : e, 'Genome' : o, 'Value' : convert(c[0])}])
            mapc[z] = e
    df = pd.DataFrame(df)
    found = [v for v in avail if v in set(df['Value'])]
    smap = {v : x for x, v in enumerate(found)}
    df['CN states'] = df.apply(lambda r : smap[r['Value']], axis=1)
    table = pd.pivot_table(df, values='CN states', columns=['Genome'], index=['Cell'], aggfunc='first')
    title = 'Copy-number states'
    palette = {}
    palette.update({(0, 0) : 'darkblue'})
    palette.update({(1, 0) : 'lightblue'})
    palette.update({(1, 1) : 'lightgray', (2, 0) : 'dimgray'})
    palette.update({(2, 1) : 'lightgoldenrodyellow', (3, 0) : 'gold'})
    palette.update({(2, 2) : 'navajowhite', (3, 1) : 'orange', (4, 0) : 'darkorange'})
    palette.update({(3, 2) : 'salmon', (4, 1) : 'red', (5, 0) : 'darkred'})
    palette.update({(3, 3) : 'plum', (4, 2) : 'orchid', (5, 1) : 'purple', (6, 0) : 'indigo'})
    colors = [palette[c] for c in found]
    cmap = LinearSegmentedColormap.from_list('multi-level', colors, len(colors))
    chr_palette = cycle(['#525252', '#969696', '#cccccc'])
    chr_colors = {c : next(chr_palette) for c in sorted(set(b[0] for b in pos), key=orderchrs)}
    para = {}
    para['data'] = table
    para['cmap'] = cmap
    para['yticklabels'] = True
    para['row_cluster'] = False
    para['method'] = 'single'
    para['metric'] = 'hamming'
    para['xticklabels'] = False
    para['col_cluster'] = False
    para['figsize'] = (12, 6)
    para['rasterized'] = True
    para['col_colors'] = pd.DataFrame([{'index' : s, 'chromosomes' : chr_colors[pos[x][0]]} for x, s in enumerate(table.columns)]).set_index('index')
    g = sns.clustermap(**para)
    addchr(g, pos)
    plt.savefig("multipassage.png", bbox_inches = 'tight', dpi=900)


orderchrs = (lambda x : int(''.join([l for l in x if l.isdigit()])))
baf = defaultdict(lambda : dict())
getgenome = (lambda p : (p[0], int(p[1]), int(p[2])))
getsample = (lambda p : p[3])
getpur = (lambda p : 1.0 - float(p[12]))
getbaf = (lambda p : float(p[9]))
form = (lambda p : (getgenome(p), getsample(p), getbaf(p)))
for d in glob.glob('hatchet/*'):
    pat = os.path.basename(d)
    with open(os.path.join(d, 'hatchet.bbc.ucn'), 'r') as i:
        for l in (l for l in i if l[0] != '#'):
            g, s, u = form(l.strip().split())
            if g[0] == 'chr17':
                baf[g][(pat, s)] = u
pos = filter(lambda g : g[0] == 'chr17', sorted(baf.keys(), key=(lambda x : (orderchrs(x[0]), x[1], x[2]))))

def draw_notp53loh():
    plt.figure(figsize=(6, 6))
    df = pd.DataFrame([{'Genome' : x, 'WGD' : wgd[p[0]], 'Cancer type' : ctype[p[0]], 'Patient' : p[0], 'BAF' : baf[g][p]} for x, g in enumerate(pos) for p in baf[g]])
    para = {}
    data = df[df['WGD']==False]
    cmap = sns.cubehelix_palette(as_cmap=True, dark=0, light=1, reverse=False, start=2.8, rot=.1)
    g = sns.kdeplot(data['Genome'], data['BAF'], cmap=cmap, n_levels=100, shade=True)
    THRES = 3000000
    tp53 = ('chr17', 7572926 - THRES, 7687550 + THRES)
    left = min(enumerate(pos), key=(lambda p : abs(p[1][1] - tp53[1])))
    right = min(enumerate(pos), key=(lambda p : abs(p[1][2] - tp53[2])))
    plt.plot((left[0], left[0]), (-0.06, 0.56), '--r', linewidth=2)
    plt.plot((right[0], right[0]), (-0.06, 0.56), '--r', linewidth=2)
    g.set(xlim=(0, len(pos)), ylim=(-0.06, 0.56))
    plt.savefig("plot_notp53loh.png", bbox_inches = 'tight', dpi=900)
    plt.figure(figsize=(6, 2))
    g = sns.kdeplot(df[df['WGD']==False]['BAF'], color='k', shade=True, legend=False, alpha=1)
    g.spines['right'].set_visible(False)
    g.spines['top'].set_visible(False)
    plt.xlim(xmin=-0.02, xmax=0.52)
    plt.grid(b=None)
    plt.savefig("dist_notp53loh.png", bbox_inches = 'tight', dpi=900)

def draw_tp53loh():
    plt.figure(figsize=(6, 6))
    df = pd.DataFrame([{'Genome' : x, 'WGD' : wgd[p[0]], 'Cancer type' : ctype[p[0]], 'Patient' : p[0], 'BAF' : baf[g][p]} for x, g in enumerate(pos) for p in baf[g]])
    data = df[df['WGD']==True]
    cmap = sns.cubehelix_palette(as_cmap=True, dark=0, light=1, reverse=False, start=2.8, rot=.1)
    g = sns.kdeplot(data['Genome'], data['BAF'], cmap=cmap, n_levels=100, shade=True)
    THRES = 3000000
    tp53 = ('chr17', 7572926 - THRES, 7687550 + THRES)
    left = min(enumerate(pos), key=(lambda p : abs(p[1][1] - tp53[1])))
    right = min(enumerate(pos), key=(lambda p : abs(p[1][2] - tp53[2])))
    plt.plot((left[0], left[0]), (-0.06, 0.56), '--r', linewidth=2)
    plt.plot((right[0], right[0]), (-0.06, 0.56), '--r', linewidth=2)
    g.set(xlim=(0, len(pos)), ylim=(-0.06, 0.56))
    plt.savefig("plot_tp53loh.png", bbox_inches = 'tight', dpi=900)
    plt.figure(figsize=(6, 2))
    g = sns.kdeplot(df[df['WGD']==True]['BAF'], color='k', shade=True, legend=False, alpha=1)
    g.spines['right'].set_visible(False)
    g.spines['top'].set_visible(False)
    plt.xlim(xmin=-0.02, xmax=0.52)
    plt.grid(b=None)
    plt.savefig("dist_tp53loh.png", bbox_inches = 'tight', dpi=900)

def addchr(g, pos, color=None):
    corners = []
    prev = 0
    for x, b in enumerate(pos):
        if x != 0 and pos[x-1][0] != pos[x][0]:
            corners.append((prev, x))
            prev = x
    corners.append((prev, x))
    ax = g.ax_heatmap
    ticks = []
    for o in corners:
        ax.set_xticks(np.append(ax.get_xticks(), int(float(o[1] + o[0] + 1) / 2.0)))
        ticks.append(pos[o[0]][0])
    ax.set_xticklabels(ticks, rotation=90, ha='center')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.xaxis.tick_top()
    ax.tick_params(axis='x', which='major', pad=20, length=0)
    ax.tick_params(axis='y', which='major', pad=20, length=0)
