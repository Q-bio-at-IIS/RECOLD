"""
Library for RECOLD algorithm

Copyright @ R. Yokota
"""

import numpy as np
import random


def data_extract(f,CondSampleNum): 
    import csv  
    dataReader = csv.reader(f)
    
    dat = []
    for row in dataReader:
        dat.append(row)

    seq = []
    for i in xrange(1,len(dat)):
        seq.append(dat[i][0])

    cmat = np.zeros((len(dat)-1,CondSampleNum)) 
    for i in range(len(dat)-1):
        for j in range(CondSampleNum):
            try :
                cmat[i][j] = int(dat[i+1][j+1]) #exclude the first column written about sequences
            except:
                pass
    return cmat, seq


def ry_isomap(D,k,n_components):
    import numpy as np
    from scipy.sparse import lil_matrix
    import scipy.sparse.csgraph as csg
    from sklearn import manifold
    from sklearn.decomposition import PCA
    
    ndata = np.shape(D)[0]
 
    # K-nearest neighbours
    indices = D.argsort()
    neighbours = indices[:,:k+1]

    h = np.zeros((ndata,ndata),dtype=float)
    for i in range(ndata):
        h[i,neighbours[i,:]] = D[i,neighbours[i,:]]
    
    Dsub = (h+h.T)/2.0
    A = lil_matrix(Dsub).tocsr()
    
    # floyd warshall algorithm
    Dinf = csg.floyd_warshall(A,directed=False)

    # mds with smacof algorithm
    MDS = manifold.MDS(n_components,metric=True, dissimilarity="precomputed", max_iter=1000, eps=1e-5)
    Rdspace = MDS.fit_transform(Dinf)
    
    # find better rotation by PCA   
    pca = PCA(n_components=2)
    pca.fit(Rdspace)
    Y= pca.transform(Rdspace)
    
    return Y


def ry_se(Dmod,k,d):
    from sklearn import manifold
    from sklearn.neighbors import NearestNeighbors

    n_components = d    
    n_neighbors = k
    
    neigh = NearestNeighbors(n_neighbors=n_neighbors+1,metric='precomputed')
    neigh.fit(Dmod)
    A = neigh.kneighbors_graph(Dmod,n_neighbors+1, mode='connectivity')
    Dsub = A.toarray()
    se = manifold.SpectralEmbedding(n_components=n_components, affinity='precomputed',n_neighbors=n_neighbors)

    Dsub = (Dsub + Dsub.T)/2.0
    Yse = se.fit_transform(Dsub)
    
    return Yse
    
def calc_dismat(seq,score_flag,go,ge):
    import itertools
    import parasail #sequence alignment package 
    
    comnum = []
    for i in itertools.combinations(range(len(seq)), 2):
        comnum.append(i)

    Dmat = np.zeros((len(seq),len(seq)))
    
    if score_flag == 'sw':    
        for i in range(len(comnum)):
            result = parasail.sw_scan_16(seq[comnum[i][0]], seq[comnum[i][1]], go, ge, parasail.blosum62) #go: open gap  penalty, ge: entend gap penalty
            Dmat[comnum[i][0]][comnum[i][1]] = result.score
    elif score_flag == 'striped_sw':
        for i in range(len(comnum)):
            result = parasail.sw_stats_striped_16(seq[comnum[i][0]], seq[comnum[i][1]], go,ge, parasail.blosum62)
            Dmat[comnum[i][0]][comnum[i][1]] = result.score
    
    Dmat = Dmat + Dmat.T

    ########### storage their own homology scores in the diagonal line ####     
    for i in range(len(seq)):
        result = parasail.sw_stats_striped_16(seq[i], seq[i], go, ge, parasail.blosum62)
        Dmat[i][i] = result.score   


    ########### the matrix of calculated scores btwn all pairs of seqs is transformed into the dissimilarity matrix. ####     
    Dmod = np.zeros((len(seq),len(seq))) #dissimilarity matrix Dmod
    for i in range(len(seq)):
        for j in range(len(seq)):
            Dmod[i,j] = 1 - 2*Dmat[i,j]/(Dmat[i,i]+Dmat[j,j]) #Distance matrix Dmat is transformed into dissimilarity matrix Dmod!!
            
    return Dmod
    
def calc_proj(Dmod,manifold_method,n_components):
    from sklearn import manifold
    if manifold_method == 'tsne':      
        tsne = manifold.TSNE(n_components=n_components, metric='precomputed', random_state=0)
        Ytsne = tsne.fit_transform(Dmod)
    elif manifold_method == 'MDS':
        tsne = manifold.MDS(n_components,metric=True, dissimilarity="precomputed", max_iter=400, eps=1e-5)
        Ytsne = tsne.fit_transform(Dmod)
    elif manifold_method == 'SE':
        n_neighbors = 10
        Ytsne = ry_se(Dmod,n_neighbors,n_components)        
    elif manifold_method == 'isomap':
        kval = 10
        Ytsne = ry_isomap(Dmod,kval,n_components)
            
    return Ytsne

def calc_kde(cmat,Ytsne,CondSampleNum,dataname):
    # load modules
    import matplotlib.pyplot as plt
    from sklearn.neighbors import KernelDensity
    from sklearn.grid_search import GridSearchCV
    
    # set the X-Y range of the projected space 
    xmin, xmax = min(Ytsne[:,1]), max(Ytsne[:,1])
    ymin, ymax = min(Ytsne[:,0]), max(Ytsne[:,0])
    lx,ly = abs(xmax-xmin), abs(ymax-ymin)
    xmin, xmax = xmin-lx/10, xmax+lx/10
    ymin, ymax = ymin-ly/10, ymax+ly/10

    lmin = min(xmin,ymin)
    lmax = max(xmax,ymax)
 
    X, Y = np.mgrid[lmin:lmax:400j, lmin:lmax:400j]
    xy = np.vstack([Y.ravel(), X.ravel()]).T
    params = {'bandwidth': np.logspace(-2, 1, 200)}

    log_Z = np.zeros((CondSampleNum,400*400)) #Column length of logZ is the total number of bins. Here is 400*400 (X*Y).
    Z = np.zeros((CondSampleNum,400*400))

    fig = plt.figure(figsize=(14, 6))

    for i in range(CondSampleNum):
        cind = cmat[:,i].nonzero()
        grid = GridSearchCV(KernelDensity(kernel='exponential'), params) # Here, we set the kernel function of KDE as the exponential func. 
    
        Ybuf = Ytsne[cind[:],:][0]
        grid.fit(Ybuf)
    
        kde = grid.best_estimator_
        Z[i,:] = np.exp(kde.score_samples(xy))
        Zbuf = np.exp(kde.score_samples(xy))
        log_Z[i,:] = kde.score_samples(xy)
        Zbuf = Zbuf*(lmax-lmin)*(lmax-lmin)/160000/sum(Zbuf*(lmax-lmin)*(lmax-lmin)/160000) # Normalization 
        Zbuf = Zbuf.reshape(X.shape)
        ax = fig.add_subplot(2,CondSampleNum/2,i+1,aspect='equal',adjustable='box-forced')
        #ax.imshow(np.rot90(Zbuf), cmap=plt.cm.gist_earth_r, extent=[xmin, xmax, ymin, ymax])
        im = ax.imshow(np.rot90(Zbuf), cmap=plt.cm.gist_earth_r, extent=[lmin, lmax, lmin, lmax],clim=(0.0, 0.0002))
        fig.colorbar(im, ax=ax)
        #ax.plot(Ytsne[cind[:],1][0], Ytsne[cind[:],0][0], 'k.', markersize=2)
        ax.plot(Ybuf[:,1], Ybuf[:,0], 'k.', markersize=2)
        ax.set_title(dataname[i], fontsize=12)
        Z[i,:] = Z[i,:]*(lmax-lmin)*(lmax-lmin)/160000/sum(Z[i,:]*(lmax-lmin)*(lmax-lmin)/160000)
       
    plt.show()  
    return Z, Zbuf, lmin, lmax
    
def calc_JSD(Z,CondSampleNum):
    import itertools
    condcom = []   
    for i in itertools.combinations(range(CondSampleNum), 2):
        condcom.append(i)

    JSD = np.zeros((CondSampleNum,CondSampleNum))
    for i in range(len(condcom)):
        for j in range(len(Z[condcom[i][0],:])):
            M = (Z[condcom[i][0],j]+Z[condcom[i][1],j])/2.0
            JSD[condcom[i][0],condcom[i][1]] += (Z[condcom[i][0],j]*np.log(Z[condcom[i][0],j])+ Z[condcom[i][1],j]*np.log(Z[condcom[i][1],j]))/2.0 - M*np.log(M) 

    JSD = JSD + JSD.T
    return JSD