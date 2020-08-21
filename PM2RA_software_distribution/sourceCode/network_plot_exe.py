import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm
from matplotlib.colors import ListedColormap

def plotPMnetwork(project_dict,topnum,networktype,edgecolormap,nodecolormap,networklayout,filter_ss):
    fileplace = project_dict.get('fileplace')
    fdr = project_dict.get('fdr')
    confidenceI = [0.99, 0.95, 0.90][project_dict.get('confidenceI')]
    resfilename = fileplace + '/PM_scores.csv'
    pm_res = pd.read_csv(resfilename, header=0, index_col=0)
    def co_PM_2D_fun(raw_pm_2d,pm1,pm2):
        return np.max([0,np.min([raw_pm_2d-pm1,raw_pm_2d-pm2])])

    if filter_ss == 'ON':
        if fdr == 'ON':
            pm_res.loc[pm_res['qvalue']>1-confidenceI,'raw_pm_2d']=0
            pm_res.loc[pm_res['qvalue1'] > 1 - confidenceI, 'pm1'] = 0
            pm_res.loc[pm_res['qvalue2'] > 1 - confidenceI, 'pm2'] = 0
        if fdr == 'OFF':
            pm_res.loc[pm_res['pvalue'] > 1 - confidenceI, 'raw_pm_2d'] = 0
            pm_res.loc[pm_res['pvalue1'] > 1 - confidenceI, 'pm1'] = 0
            pm_res.loc[pm_res['pvalue2'] > 1 - confidenceI, 'pm2'] = 0

    pm_res['co_PM_2D']=pm_res.apply(lambda x: co_PM_2D_fun(x.raw_pm_2d, x.pm1,x.pm2), axis = 1)
    target = 'raw_pm_2d'
    if networktype == 'Interaction PM scores':
        target = 'co_PM_2D'
    if networktype=='2D PM Scores':
        target = 'raw_pm_2d'
    pm_res=pm_res[pm_res[target]>0]
    pm_res=pm_res.sort_values(by=target,ascending=False)
    pm_res=pm_res.reset_index(drop=True)
    taxnum = np.min([len(pm_res), topnum])
    pm_res = pm_res[0:taxnum]

    G = nx.Graph()
    for i in pm_res.index:
        G.add_weighted_edges_from([(pm_res.loc[i,'taxa1'], pm_res.loc[i,'taxa2'], pm_res.loc[i,target])])

    
    if networklayout=='Spring Layout':
        pos_pm = nx.spring_layout(G)
    if networklayout == 'Fruchterman Reingold Layout':
        pos_pm = nx.fruchterman_reingold_layout(G)
    if networklayout == 'Shell Layout':
        pos_pm = nx.shell_layout(G)
    if networklayout == 'Spectral Layout':
        pos_pm = nx.spectral_layout(G)
    if networklayout == 'Kamada Kawai Layout':
        pos_pm = nx.kamada_kawai_layout(G)
    if networklayout == 'Circular Layout':
        pos_pm = nx.circular_layout(G)

    Blues_big = cm.get_cmap('Blues', 512)
    newBlues = ListedColormap(Blues_big(np.linspace(0.3, 0.85, 256)))

    Greys_big = cm.get_cmap('Greys', 512)
    newGreys = ListedColormap(Greys_big(np.linspace(0.3, 0.85, 256)))

    Purples_big = cm.get_cmap('Purples', 512)
    newPurples = ListedColormap(Purples_big(np.linspace(0.3, 0.85, 256)))

    Greens_big = cm.get_cmap('Greens', 512)
    newGreens = ListedColormap(Greens_big(np.linspace(0.3, 0.85, 256)))

    Oranges_big = cm.get_cmap('Oranges', 512)
    newOranges = ListedColormap(Oranges_big(np.linspace(0.3, 0.85, 256)))

    Reds_big = cm.get_cmap('Reds', 512)
    newReds = ListedColormap(Reds_big(np.linspace(0.3, 0.85, 256)))

    YlOrBr_big = cm.get_cmap('YlOrBr', 512)
    newYlOrBr = ListedColormap(YlOrBr_big(np.linspace(0.3, 0.85, 256)))

    YlOrRd_big = cm.get_cmap('YlOrRd', 512)
    newYlOrRd = ListedColormap(YlOrRd_big(np.linspace(0.3, 0.85, 256)))

    OrRd_big = cm.get_cmap('OrRd', 512)
    newOrRd = ListedColormap(OrRd_big(np.linspace(0.3, 0.85, 256)))

    PuRd_big = cm.get_cmap('PuRd', 512)
    newPuRd = ListedColormap(PuRd_big(np.linspace(0.3, 0.85, 256)))

    BuPu_big = cm.get_cmap('BuPu', 512)
    newBuPu = ListedColormap(BuPu_big(np.linspace(0.3, 0.85, 256)))

    GnBu_big = cm.get_cmap('GnBu', 512)
    newGnBu = ListedColormap(GnBu_big(np.linspace(0.3, 0.85, 256)))

    PuBu_big = cm.get_cmap('PuBu', 512)
    newPuBu = ListedColormap(PuBu_big(np.linspace(0.3, 0.85, 256)))

    YlGnBu_big = cm.get_cmap('YlGnBu', 512)
    newYlGnBu = ListedColormap(YlGnBu_big(np.linspace(0.3, 0.85, 256)))

    RdPu_big = cm.get_cmap('RdPu', 512)
    newRdPu = ListedColormap(RdPu_big(np.linspace(0.3, 0.85, 256)))

    PuBuGn_big = cm.get_cmap('PuBuGn', 512)
    newPuBuGn = ListedColormap(PuBuGn_big(np.linspace(0.3, 0.85, 256)))

    BuGn_big = cm.get_cmap('BuGn', 512)
    newBuGn = ListedColormap(BuGn_big(np.linspace(0.3, 0.85, 256)))

    YlGn_big = cm.get_cmap('YlGn', 512)
    newYlGn = ListedColormap(YlGn_big(np.linspace(0.3, 0.85, 256)))
    colormapdict = {
        'Greys': newGreys,
        'Purples': newPurples,
        'Blues': newBlues,
        'Greens': newGreens,
        'Oranges': newOranges,
        'Reds': newReds,
        'Yellow-orange-brown': newYlOrBr,
        'Yellow-orange-red': newYlOrRd,
        'Orange-red': newOrRd,
        'Purple-red': newPuRd,
        'Red-purple': newRdPu,
        'Blue-purple': newBuPu,
        'Green-blue': newGnBu,
        'Purple-blue': newPuBu,
        'Yellow-green-blue': newYlGnBu,
        'Purple-blue-green': newPuBuGn,
        'Blue-green': newBuGn,
        'Yellow-green': newYlGn
    }
    edgecolormap=colormapdict.get(edgecolormap)
    nodecolormap=colormapdict.get(nodecolormap)

    nodes, nodes_degree = zip(*G.degree)
    nodes = list(nodes)
    nodes_degree = np.array(nodes_degree)
    vmax_n = max(nodes_degree)
    vmin_n = min(nodes_degree)
    nodes_size = (nodes_degree - vmin_n) / (vmax_n - vmin_n) * 700 + 300

    edges, weights = zip(*nx.get_edge_attributes(G, 'weight').items())
    weights = np.array(weights)
    vmax_w = max(weights)
    vmin_w = min(weights)
    edge_width = (weights - vmin_w) / (vmax_w - vmin_w) * 5 + 1
    nlabels = dict(zip(nodes, nodes))


    plt.figure(figsize=(20,20))

    if (vmin_n != vmax_n) and (vmin_w != vmax_w):
        nx.draw(G, pos_pm, edge_color=weights, width=edge_width, edge_cmap=edgecolormap, nodelist=nodes,
                node_size=nodes_size, node_color=nodes_degree, cmap=nodecolormap)
        nx.draw_networkx_labels(G, pos_pm, nlabels, font_size=8)
        edgemap = plt.cm.ScalarMappable(cmap=edgecolormap, norm=plt.Normalize(vmin=vmin_w, vmax=vmax_w))
        edgemap._A = []
        plt.colorbar(edgemap, shrink=0.3, label='PM scores')

        nodemap = plt.cm.ScalarMappable(cmap=nodecolormap, norm=plt.Normalize(vmin=vmin_n, vmax=vmax_n))
        nodemap._A = []
        plt.colorbar(nodemap, shrink=0.3, label='Node degree')
    if (vmin_n != vmax_n) and (vmin_w == vmax_w):
        nx.draw(G, pos_pm, width=edge_width, nodelist=nodes,
                node_size=nodes_size, node_color=nodes_degree, cmap=nodecolormap)
        nx.draw_networkx_labels(G, pos_pm, nlabels, font_size=8)
        nodemap = plt.cm.ScalarMappable(cmap=nodecolormap, norm=plt.Normalize(vmin=vmin_n, vmax=vmax_n))
        nodemap._A = []
        plt.colorbar(nodemap, shrink=0.3, label='Node degree')
    if (vmin_n == vmax_n) and (vmin_w != vmax_w):
        nx.draw(G, pos_pm, edge_color=weights, width=edge_width, edge_cmap=edgecolormap, nodelist=nodes,
                node_size=600, node_color='red')
        nx.draw_networkx_labels(G, pos_pm, nlabels, font_size=8)
        edgemap = plt.cm.ScalarMappable(cmap=edgecolormap, norm=plt.Normalize(vmin=vmin_w, vmax=vmax_w))
        edgemap._A = []
        plt.colorbar(edgemap, shrink=0.3, label='PM scores')
    if (vmin_n == vmax_n) and (vmin_w == vmax_w):
        nx.draw(G, pos_pm)
        nx.draw_networkx_labels(G, pos_pm, nlabels, font_size=8)

    plt.axis('off')
    # plt.show()
    plt.savefig(fileplace + '/PM_network.pdf',dpi=600, bbox_inches='tight')