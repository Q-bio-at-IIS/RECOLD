
import os 
### The directory where this script is
rootdir= os.path.dirname(os.path.abspath(__file__))
### Import modules
import sys; sys.path.append(rootdir+'/mymodule') #need to tidy my modules.

# import the modules for load&save
from scipy.io import savemat
from scipy.io import loadmat
    
"""
Calculation of Dissimilarity Matrix
"""
if True:
    from mymodules import data_extract, calc_dismat
    # Load the data of [Rempala+,2011].
    f= open("Please write the absolute pathname of the supplymentary data of Rempala et al., (2011)",'rU')

    # The number of the condition of samples
    CondSampleNum = 8

    # extract the AA-seq-vector(seq) & counts matrix(cmat) from the csv file
    cmat, seq = data_extract(f,CondSampleNum)

    ### set the parameters for calc_dismat func.###
    # Scores for gap open and extention penalties
    go = 10
    ge = 1
    # select the alignment algorithm (sw,sg,nw) # Caution!! Here, we implemented ONLY the sw and "striped"-sw!!
    score_flag = 'striped_sw'
    ###############################################

    ################################################################
    # calculate the dissimilarity matrix among all observed AA seqs.
    Dmod = calc_dismat(seq,score_flag,go,ge)

    ## save data###############
    savefilename = 'Rempala_Dmod_BLOSUM62_10_1.mat'
    savemat('MatData/'+savefilename, {'Dmod':Dmod})
    savemat('MatData/'+'cmat.mat', {'cmat':cmat})
else :
        print('skip Dissimilarity matrix calculation')


"""
Calculation of Projection by tSNE
"""
if True:
    ############################################################################
    # import the "ry_isomap, "ry_se", & "calc_proj"" functions from mymodules.py
    from mymodules import ry_isomap, ry_se, calc_proj

    # load dissimilaroty matrix (Dmod) of AA seq.
    Dmod= loadmat("MatData/Rempala_Dmod_BLOSUM62_10_1.mat")['Dmod']

    ### set the parameters for calc_proj func.###
    # Select the algorithm of manifold learning. Choice from {tsne,MDS,SE,isomap}.
    manifold_method = "tsne"
    # set the dimension of the projected space
    n_components = 2
    ###############################################

    # calculate the projection
    Y = calc_proj(Dmod,manifold_method,n_components)

    ## save data
    savefilename = 'Rempala_Y_tsne_BLOSUM62_10_1.mat'
    savemat('MatData/'+savefilename, {'Y':Y})
else :
        print('skip tSNE calculation')

"""
Calculate KDE in the projected space
"""
if True:
    # import the "calc_kde" function from mymodules.py
    from mymodules import calc_kde

    ### set the parameters for calc_kde func.###
    # load counts matrix
    cmat = loadmat("MatData/cmat.mat")['cmat']
    # load positions of each AA. seq in the projected space.
    Y = loadmat("MatData/Rempala_Y_tsne_BLOSUM62_10_1.mat")['Y']

    # number of the condition of samples
    CondSampleNum = 8
    # set datanames
    dataname = ['EpTN1','EpTR1','EpTN2','EpTR2','WtTN1','WtTR1','WtTN2','WtTR2']
    #############################################

    # calcualte KDE. Here we used an exponential function as a kernel function of KDE.
    Z, Zbuf, lmin, lmax = calc_kde(cmat,Y,CondSampleNum,dataname)

    ## save data
    savefilename = 'Rempala_Z_tsne_BLOSUM62_10_1.mat'
    savemat('MatData/'+savefilename, {'Z':Z})
else:
    print('skip KDE calculation')


"""
Calculate JSD between repertoire
"""
if True:
    from mymodules import calc_JSD

    ### set the parameters for calc_kde func.###
    # load kde function
    Z = loadmat("MatData/Rempala_Z_tsne_BLOSUM62_10_1.mat")['Z']
    # number of the condition of samples
    CondSampleNum = 8
    # set datanames
    dataname = ['EpTN1','EpTR1','EpTN2','EpTR2','WtTN1','WtTR1','WtTN2','WtTR2']
    #############################################

    # calcualte JSD matrix.
    JSD = calc_JSD(Z,CondSampleNum)

    ## save data
    savefilename = 'Rempala_JSD_tsne_BLOSUM62_10_1.mat'
    savemat('MatData/'+savefilename, {'JSD':JSD})
else:
    print('skip JSD calculation')
