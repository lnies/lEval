from read_write_functions import *
import ROOT as root
# import xit
import numpy as np
from array import array

def root_1D_unbinned_gauss_fit(liste):
    '''Function takes 1D list with unbinned data and returns a list with gaussian fit parameters: x, x_unc., sigma, sigma_unc.'''
    import_list = [float(x) for x in liste]
    # get the x where most of the ions arrive
    hist, hilf_edges = np.histogram(import_list, bins=2000, range=(0,200))
    x_edges = []
    for i in range(len(hilf_edges)-1):
        x_edges.append(np.mean([hilf_edges[i], hilf_edges[i+1]]))
    x_max = x_edges[index_max_value(hist)]
    # calculate mean of limited range of x interval
    mean_x = np.mean(import_list[int(x_max)-10:int(x_max)+10])
    std_x = np.std(import_list[int(x_max)-10:int(x_max)+10])
    # import list to ROOT
    # f = root.TFile( 'test.root', 'recreate' )
    tree = root.TTree( 'tree', 'tree' )

    x = array('d', [ 0. ])
    tree.Branch('x', x, 'x/D')
    for i in range(len(import_list)):
        x[0] = import_list[i]
        tree.Fill()
    # unbinned max. likelihood fit
    x = root.RooRealVar('x', 'x', 0, 200, 'us')
    ds = root.RooDataSet('x-data', 'x-data', root.RooArgSet(x), root.RooFit.Import(tree))
    meant = root.RooRealVar('meant', 'meant', mean_x, mean_x-5, mean_x+5, 'us')
    sigmax = root.RooRealVar('sigmax', 'sigmax', 0.3, 0.01, 5)
    gaussx = root.RooGaussian('gaussx', 'Gaussian distribution', x, meant, sigmax)
    x.setRange('range', mean_x-10, mean_x+10)
    result = gaussx.fitTo(ds, root.RooFit.Range('range')) # , root.RooFit.NumCPU(4))
    # fit results
    x_pos = meant.getValV()
    x_pos_err = meant.getError()
    x_sigma = sigmax.getValV()
    x_sigma_err = sigmax.getError()
    # plot

    c = root.TCanvas('c', 'x', 1000, 700)
    tree.Draw('x')
    tframe = x.frame()
    # ds.plotOn(tframe, root.RooFit.Binning(2000))
    gaussx.plotOn(tframe,root.RooFit.LineColor(root.kRed))
    tframe.Draw()
    c.Update()
    pdffile = 'file_name' + '_hist.pdf'
    # # c.SaveAs(pdffile)
    root.gApplication.Run()
    f.Write()
    f.Close()

    return([x_pos, x_pos_err, x_sigma, x_sigma_err])


def index_max_value(liste):
########## find index of maximum in list
    index_max = 0
    max_value = -99999999999.
    for i in range(len(liste)):
        if liste[i] > max_value:
            max_value = liste[i]
            index_max = i
    return(index_max)



if __name__ == '__main__':
    liste = [0,0,99,100,100,0,0]
    print root_1D_unbinned_gauss_fit(liste)
