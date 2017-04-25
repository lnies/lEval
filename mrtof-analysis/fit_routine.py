# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 13:41:54 2016
Modified on Tuesday April 25 2017

example of command line argument : Path_to_file xlow xhigh binning model name output_name [if set tlow_Corr thigh_corr]
E:/Data_Analysis/Data_Analysis_Cr/63Cr_2016_MRTOF/63Cr_352_1000revs/Cr_run_354
27775120 27775500 0 Gaussian ref_133Cs Gaussian_ref_133Cs [tlow_cor thigh_cor]

@author: mamougeo
@author: datanasov
"""

from ROOT import TCanvas, TH2D, gApplication, gPad, TH1,TArrayD,TF1,TObjArray,TMath, TAxis
from ROOT import RooRealVar, RooRandom, RooArgSet, RooGaussian, RooDataSet, RooDataHist, RooFitResult, RooAbsData, \
    RooChiSquarePdf
from ROOT import RooPlot, RooFit, TArrow, RooArgList, RooMinuit, RooChi2Var, RooGaussModel, RooDecay, RooGenericPdf, \
    RooNLLVar, RooProfileLL, RooMCStudy, gStyle
import sys,os
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm
from time import time
from math import *
import peakutils
from peakutils.plot import plot as pplot
#from scipy.signal import find_peaks_cwt
import glob as gb
import ConfigParser

def create_dict(key_list, val_list):
    """
    The function creates dictionary from given key_list and value_list
    :param key_list:
    :param val_list:
    :return: Dictionary of type (key:val)
    """
    temp_dict = {}
    for arg in range(len(key_list) - 1):
        temp_dict[key_list[arg]] = val_list[arg]
    return temp_dict


def load887file(file887):
    """
    Load run-file configuration
    :param file887:
    :return: Config dictionary
    """
    key_list, val_list = np.loadtxt(file887, delimiter='=', usecols=(0, 1), dtype='object', skiprows=7, unpack=True)
    return create_dict(key_list, val_list)



def load_asc(datafile):
    """
    Load parameters in ascii file
    :param datafile:
    :return: Ascii config dictionnary
    """
    key_list, val_list = np.genfromtxt(datafile, delimiter='=', usecols=(0, 1), dtype='object', skip_header=1,
                                       max_rows=10, unpack=True)
    return create_dict(key_list, val_list), np.loadtxt(datafile, skiprows=13)


def create2d_histogram(name, dict887):
    """
    create an empty 2D histogram to contain the whole MCDWIN spectrum
    :param name:
    :param dict887:
    :return: an empty TH2D histogram
    """
    # lower edge of first bin
    xlow = float(dict887.get('caloff'))-0.5* float(dict887.get('calfact'))
    # upper edge of last bin not included in last bin
    xup = float(dict887.get('caloff')) + (float(dict887.get('range'))-0.5) * float(dict887.get('calfact'))
    # TH2D('HISTO_IDENT', name, nbbinsx, xlow, xup, nbbinsy, ulow, yup)
    return TH2D('HISTO_MCDWIN', name, int(dict887.get('range')), xlow, xup, int(dict887.get('cycles')), -0.5,
                float(dict887.get('cycles'))-0.5)


def fill2d_histo(th2d, data):
    """
    fill a 2D histogram created with create2d_histogram(name, dict887)
    with the mcdwin data.
    :param th2d:
    :param data:
    :return: None
    """
    for i in xrange(len(data[:, 0])):
        counts = data[i, 2]
        # the indeces in MCDWIN start at 0. In ROOT the bin 0 is underflow so need to add +1
        # to the mcdwin index to fill the corresponding ROOT bin.
        th2d.SetBinContent(int(data[i, 0]) + 1, int(data[i, 1]) + 1, counts)


def projection_x(th2d):
    """
    create the projection on the x axis of a th2d
    :param th2d:
    :return: TH1D of the projected data on the x axis
    """
    xproj = th2d.ProjectionX('Projection on ToF axis')
    xproj.SetXTitle("Time Of Flight")
    xproj.SetYTitle("Nb of Counts")
    return xproj

def find_peak(tof_spectrum, nbins,min_range, thres_Val):
    """
	find estimated peaks' positions on the projection of time-of-flight spectrum
	"""

    iflag = 1
    wid_bin = 8
    wid_flag = 1
    xbin_range = np.arange(1, nbins+1)
    ybin_value = []
    ybin_nb = []
    indexes = []
    peaks_x = []  # enhancing searching
    for i in range(1, nbins+1):
        ybin_value.insert(i, (float(tof_spectrum.GetXaxis().GetBinCenter(i))))
        ybin_nb.insert(i, (int(tof_spectrum.GetBinContent(i))))
    ybin_value = np.array(ybin_value)
    ybin_nb = np.array(ybin_nb)
    plt.plot(xbin_range, ybin_nb)
    #plt.show()

    while iflag:
        indexes = peakutils.indexes(ybin_nb, thres=thres_Val, min_dist=min_range)
        print "number of peaks=", len(indexes)
        nb_peak = len(indexes)
        for i_dex in range(len(indexes)):
            print indexes[i_dex], ybin_value[indexes[i_dex]], xbin_range[indexes[i_dex]]
        pplot(xbin_range, ybin_nb, indexes)
        plt.title('First estimate')
        plt.show()

        print "press 'c' to do another search or ENTER to quit search test"
        opt = raw_input()
        if opt == "c":
            min_range = float(raw_input("define a fit range:"))
            thres_Val = float(raw_input("define a threshold:"))
            #indexes = []
            continue
        else:
            iflag = 0

        while True:
            try:
                print 'current range: ', min_range
                print 'current width: ', wid_bin
                peaks_x = peakutils.interpolate(xbin_range, ybin_nb, ind=indexes, width=min_range/wid_bin)
                print 'peak_x', peaks_x
                break

            except:
                print 'crash!!!!!! RuntimeErrors'
                print 'please modify the bin width_: options from 1 to 10'

                wid_bin = int(raw_input())
            iflag = 1
            break

    peaks_x = [int(round(ii)) for ii in peaks_x]
    tlow = ybin_value[peaks_x] - 0.5*float(min_range)
    thigh = ybin_value[peaks_x] + 0.5*float(min_range)
    return tlow, thigh, nb_peak

def projection_y(th2d):
    """
    create the projection on the y axis of a th2d
    :param th2d:
    :return: TH1D of the projected data on the y axis
    """
    yproj = th2d.ProjectionY('Projection on y axis')
    yproj.SetXTitle("Slices")
    yproj.SetYTitle("Nb of Counts")
    return yproj


def get_roi(th2d, tlow, thigh):
    """
    restrict the histogram to the roi defined by tlow and thigh.
    from the tof distribution in the roi computes: mean,rms,skewness and integral
    :param th2d:
    :param tlow:
    :param thigh:
    :return: the Roi mean, rms,integral, skewness and TH2D of the roi
    """
    temp_th2d = TH2D(th2d)
    temp_th2d.GetXaxis().SetRangeUser(tlow, thigh)
    roi_rms = temp_th2d.GetStdDev(1)
    roi_mean = temp_th2d.GetMean(1)
    roi_skewness = temp_th2d.GetSkewness(1)
    roi_counts = temp_th2d.Integral()
    return roi_mean, roi_rms, roi_counts, roi_skewness, temp_th2d


def get_1d_pdf(model_name, tlow, thigh, mean, rms, skewness):
    """
    build the 1d pdf used for the analysis.
    Initialise the parameters accordingly.
    :param model_name:
    :param tlow:
    :param thigh:
    :param mean:
    :param rms:
    :param skewness:
    :return: a list containing the pdf object as well as the roorealvar corresponding, the number of such roorealvar
    """
    # initialisation of parameters SKEWNESS A TO BE POSITVE (TAIL ON THE RIGHT SIDE OF THE PEAK) hence the abs value
    ini_rooreal_dict = {'Gaussian': (mean, rms, skewness),
                        'SEG': (
                            mean - pow(abs(skewness) / 2.0, 1 / 3.0),
                            pow(rms - pow(abs(skewness) / 2.0, 1 / 3.0), 1 / 2.0),
                            pow(abs(skewness) / 2.0, 1 / 3.0)),
                        'EGH': (mean, rms, skewness)
                        }
    in_mean, in_rms, in_tau = ini_rooreal_dict.get(model_name)
    # parameter of models :
    t = RooRealVar("t", "t", tlow, thigh, "ns")
    mean = RooRealVar("mean", "mean", in_mean, tlow, thigh, "ns")
    sigma = RooRealVar("sigma", "sigma", in_rms, 0.1 * rms, 2 * rms, "ns")
    tau = RooRealVar("tau", "tau", in_tau, -100, 100, "ns")

    # RooGaussModel for SEG and DEG:
    gaussm1 = RooGaussModel('gaussm1', "Gaussian distribution", t, mean, sigma)
    # dictionary defining the available pdf (feel free to add your own)
    # PLEASE DO NOT FORGET THE NUMBER OF ROOREALVAR IN THE MODEL
    pdf_dict = {
        'Gaussian': [[RooGaussian('gauss', 'Gaussian Model', t, mean, sigma), t, mean, sigma], 3],
        'SEG': [
            [RooDecay('SEG model', 'SEG model', t, tau, gaussm1, RooDecay.SingleSided), t, mean, sigma, tau, gaussm1],
            4],
        'DEG': [
            [RooDecay('DEG model', 'DEG model', t, tau, gaussm1, RooDecay.DoubleSided), t, mean, sigma, tau, gaussm1],
            4],
        'EGH': [[RooGenericPdf('EGH model', 'EGH model', 'exp(((-(@0-@1)^2)/(2*@2^2+@3*(@0-@1))))',
                               RooArgList(t, mean, sigma, tau)), t, mean, sigma, tau], 4]
    }

    pdf_definition = pdf_dict.get(model_name)
    return pdf_definition[0], pdf_definition[1]


def make_pull(roorealvar, frame, filename, t_stamp):
    """
    extract residuals and pulls of the data with respect to the best fit line
    plot the residuals and the pull as function of the tof
    save figures as eps and pdf
    :param roorealvar:
    :param frame:
    :param filename:
    :return: None
    """
    # extract residual and pull histogram
    hresid = frame.residHist()
    hpull = frame.pullHist()
    # make the plots
    hredframe = roorealvar.frame(RooFit.Title("Residual Distribution"))
    hredframe.addPlotable(hresid, "P")
    hpullframe = roorealvar.frame(RooFit.Title("Pull Distribution"))
    hpullframe.addPlotable(hpull, "P")
    c = TCanvas('Pulls', 'Residuals as a function of the time of flight', 200, 10, 700, 500)
    c.Divide(1, 2)
    c.cd(1)
    gPad.SetLeftMargin(0.15)
    hredframe.GetYaxis().SetTitleOffset(1.6)
    hredframe.Draw()
    c.cd(2)
    gPad.SetLeftMargin(0.15)
    hpullframe.GetYaxis().SetTitleOffset(1.6)
    hpullframe.Draw()
    pdffile = filename + t_stamp + '_resid.pdf'
    #epsfile = filename + t_stamp + '_resid.eps'
    c.SaveAs(pdffile)
    #c.SaveAs(epsfile)


def make_profilell(pdf, roodatahist, filename, t_stamp, nb_params, rooreallist):
    """
    build the nnl
    build the profile likelyhood (see roofit tutorial pdf) for each parameters of the pdf
    draw nnl and profile likelyhood in a canvas containing as many panels as there are parameters
    save the plots to eps and pdf
    :param pdf:
    :param roodatahist:
    :param filename:
    :param nb_params:
    :param rooreallist:
    :return: None
    """
    # reduce the range of the parameters for computation of profillell and later for the mcstudy
    for i in xrange(nb_params):
        roorealvar = rooreallist[i + 1]
        roorealvar.setRange(roorealvar.getValV() + 10 * roorealvar.getAsymErrorLo(),
                            roorealvar.getValV() + 10 * roorealvar.getAsymErrorHi())
    nll = pdf.createNLL(roodatahist)
    c = TCanvas('c1', 'Profile Likelyhood', 200, 10, 700, 500)
    c.Divide(nb_params)
    for i in xrange(nb_params):
        roorealvar = rooreallist[i + 1]
        profile_llmean = nll.createProfile(RooArgSet(roorealvar))
        pllframe = roorealvar.frame()
        pllframe.SetLabelSize(0.02)
        pllframe.SetLabelSize(0.02, "Y")
        nll.plotOn(pllframe, RooFit.ShiftToZero())
        profile_llmean.plotOn(pllframe, RooFit.LineColor(RooFit.kRed))
        pllframe.SetMinimum(0)
        pllframe.SetMaximum(3)
        c.cd(i + 1)
        pllframe.Draw()
    pdffile = filename + t_stamp + '_pll.pdf'
    #epsfile = filename + t_stamp + '_pll.eps'
    c.SaveAs(pdffile)
    #c.SaveAs(epsfile)


def plot_mcdwin(histogramme_mcdwin, filename):
    """
    plot the mcdwin 2d histogram to eps and pdf
    :param histogramme_mcdwin:
    :param filename:
    :return: None
    """
    c = TCanvas('c', 'mcdwin data', 200, 10, 700, 500)
    histogramme_mcdwin.SetXTitle("Time Of Flight")
    histogramme_mcdwin.SetYTitle("Cycles")
    histogramme_mcdwin.Draw("col")
    pdffile = filename + '_mcdwin.pdf'
    #epsfile = filename + '_mcdwin.eps'
    c.SaveAs(pdffile)
    #c.SaveAs(epsfile)


def plot_yproj(histogramme_mcdwin, filename):
    """
    plot the y projection of the mcdwin 2d histogram to eps and pdf
    :param histogramme_mcdwin:
    :param filename:
    :return: None
    """
    c = TCanvas('c', 'Y projection', 200, 10, 700, 500)
    # create Yprojection of MM6 data
    y_projection = projection_y(histogramme_mcdwin)
    y_projection.Draw()
    pdffile = filename + '_yproj.pdf'
    #epsfile = filename + '_yproj.eps'
    c.SaveAs(pdffile)
    #c.SaveAs(epsfile)


def plot_xproj(histogramme_mcdwin, filename):
    """
    plot the x projection of the mcdwin 2d histogram to eps and pdf
    :param histogramme_mcdwin:
    :param filename:
    :return: None
    """
    c = TCanvas('c', 'ToF projection', 200, 10, 700, 500)
    # create Yprojection of MM6 data
    x_projection = projection_x(histogramme_mcdwin)
    x_projection.SetLabelSize(0.02, "X")
    x_projection.SetLabelSize(0.02, "Y")
    x_projection.Draw()
    pdffile = filename + '_xproj.pdf'
    #epsfile = filename + '_xproj.eps'
    c.SaveAs(pdffile)
    #c.SaveAs(epsfile)


def plot_peak(roorealvar, pdf, roodatahist, filename, t_stamp):
    """
    plot the peak as well as the best fit function to eps and pdf
    :param roorealvar:
    :param pdf:
    :param roodatahist:
    :param filename:
    :return: the frame object on the roorealvar t
    """
    tframe = roorealvar.frame(RooFit.Title(filename))
    tframe.SetLabelSize(0.02)
    tframe.SetLabelSize(0.02, "Y")
    roodatahist.plotOn(tframe, RooFit.DrawOption("B"), RooFit.FillColor(RooFit.kGray),RooFit.DataError(RooAbsData.Poisson))
    pdf.plotOn(tframe, RooFit.LineColor(RooFit.kRed))
    pdf.paramOn(tframe, RooFit.Layout(0.6))
    #roodatahist.statOn(tframe)
    c = TCanvas('c','MCDWIN data', 200, 10, 700, 500)
    pdffile = filename + t_stamp + '_peak.pdf'
    #epsfile = filename + t_stamp + '_peak.eps'
    tframe.Draw()
    c.SaveAs(pdffile)
    #c.SaveAs(epsfile)
    return tframe


def rebin_1dhist(th1d, nbin=2):
    """
    rebin the data by combinning nbin adjacent bins together
    :param th1d:
    :param nbbin:
    :return: None
    """
    th1d.Rebin(nbin)
    return


def printtofile(pdf_name, target_filename, rooarglist, t_stamp, cov, nbparam, tlow, thigh, p, frame, filename):
    """
    print p_value, red_chi2, fit results with both hese and minos errors, covariance matrix to file
    :param pdf_name:
    :param rooarglist:
    :param cov:
    :param nbparam:
    :param tlow:
    :param thigh:
    :param p:
    :param frame:
    :param filename:
    :return: None
    """
    file_str = target_filename + t_stamp + ".res"
    temp_list = np.array([], dtype='str')
    # one needs to give the nb of fit parameters as the function chiSquare takes 0 as default
    red_chi2 = frame.chiSquare(nbparam)
    head = 'MCDWIN file' + '\t' + 'Range_low' + '\t' + 'Range_high' + '\t' + 'Model Name' + '\t' + \
           'Parameters with Error' + '\t'+ 'Red_Chi2' + '\t' + 'P-value' + '\t' + 'Covariance-Matrix'
    temp_list = np.append(temp_list, [filename, str(tlow), str(thigh), pdf_name])
    for i in xrange(nbparam):
        temp_list = np.append(temp_list, [rooarglist[i+1].getValV(), rooarglist[i + 1].getError()])
    temp_list = np.append(temp_list, [red_chi2, p])

    l = [cov[i][j] for i in xrange(nbparam) for j in xrange(nbparam)]
    temp_list = np.append(temp_list, [l])

    np.savetxt(file_str, (temp_list,), fmt='%s', header=head, newline='\n', delimiter='\t')



def goodness_of_fit(pdf, data, nbparams):
    """
    print p_value, red_chi2, fit results with both hesse and minos errors, covariance matrix to file
    :param rooarglist:
    :param cov:
    :param nbparam:
    :param tlow:
    :param thigh:
    :param p:
    :param frame:
    :param filename:
    :return: None
    """
    # USE 68.3% POISSON CENTRAL INTERVAL AS UNCERTAINTY (CAN BE ASSYMETRIC AT LOW STATISTIC)
    chi2val = RooChi2Var("chi2", "chi2", pdf, data, RooFit.DataError(RooAbsData.Poisson)).getValV()
    ndofval = data.numEntries() - nbparams
    chi2 = RooRealVar("x", "x", chi2val, 0, 500, "")
    ndof = RooRealVar("ndof", "ndof", ndofval, "")
    chi2_pdf = RooChiSquarePdf("chi2_pdf", "chi2 distribution", chi2, ndof)
    # compute integral from 0 to chi2
    chi2_cdf = chi2_pdf.createCdf(RooArgSet(chi2))
    # proba = chi2_cdf.getValV()
    # # p-value is complementary to one of proba
    # proba = 1 - proba
    # print '--------------------- P VALUE ---------------------------'
    # return proba
    print '--------CHI2----------------'
    print ndofval, chi2val, sqrt(2*ndofval), ndofval-sqrt(2*ndofval),  ndofval+sqrt(2*ndofval)
    if ndofval-sqrt(2*ndofval) <= chi2val <= ndofval+sqrt(2*ndofval):
        return 'TRUE'
    else:
        return 'FALSE'


def mcstudy(rooarglist, ini_val, nbparams, model, nbbins, nbevents, filename):
    """
    Monte carlo validation of the fit
    1500 toy spectrums containing nbevents events are created according the the fitted pdf
    each of these spectrum are fitted to extract the parameters
    the distribution of each parameters is plotted
    the pulls are also plotted and the pull distribution is fitted with a gaussian(mean,sigma)
    if fit is unbiased and errors are ok --> mean = 0 sigma =1
    for each parameters save plots to eps and pdf
    :param rooarglist:
    :param ini_val:
    :param nbparams:
    :param model:
    :param nbbins:
    :param nbevents:
    :param filename:
    :return: None
    """
    x = RooRealVar(rooarglist[0])
    x.setBins(nbbins)
    mgr = RooMCStudy(model, RooArgSet(x), RooFit.Binned(RooFit.kTRUE), RooFit.Silence(),
                     RooFit.FitOptions(RooFit.Save(RooFit.kTRUE)))
    mgr.generateAndFit(1500, int(nbevents), RooFit.kFALSE)
    for i in xrange(nbparams):
        c = TCanvas('c' + str(i), '', 200, 10, 700, 500)
        c.Divide(1, 2)
        c.cd(1)
        pframe = rooarglist[i + 1].frame()
        mgr.plotParamOn(pframe)
        pframe.SetLabelSize(0.02)
        pframe.SetLabelSize(0.02, "Y")
        par_val = TArrow(ini_val[i], 0, ini_val[i], 70, 0.02, "<")
        par_val.SetLineColor(2)
        pframe.addObject(par_val)
        pframe.Draw()
        c.cd(2)
        mpframe = mgr.plotPull(rooarglist[i + 1], RooFit.FrameRange(-5, 5), RooFit.FrameBins(50),
                               RooFit.FitGauss(RooFit.kTRUE))
        mpframe.SetLabelSize(0.02)
        mpframe.SetLabelSize(0.02, "Y")
        mpframe.getObject(2).SetX1(0.7)
        mpframe.getObject(2).SetY1(0.8)
        mpframe.getObject(2).SetX2(0.9)
        mpframe.getObject(2).SetY2(0.9)
        mpframe.Draw()
        pdffile = filename + '_mcs_par-{}.pdf'.format(i + 1)
        #epsfile = filename + '_mcs_par-{}.eps'.format(i + 1)
        c.SaveAs(pdffile)
        #c.SaveAs(epsfile)


def tof_variation(data, dict887, tlow, thigh, filename):
    """
    A WLS model (tof = a*slice+b) is fitted to the data
    A t-test is performed to see if the data is compatible with H0:a = 0
    Not used anymore
    :param data:
    :param dict887:
    :param tlow:
    :param thigh:
    :param filename:
    :return: None
    """
    cp_data = np.copy(data)
    reduced_data = list()
    # split the data array in a list containing an array per slices
    splitted_data = np.split(cp_data, np.where(np.diff(cp_data[:, 1]))[0] + 1)
    for arr in splitted_data:
        # calibrate each tof axis of each slices.
        arr[:, 0] = [float(dict887.get('caloff')) + x * float(dict887.get('calfact')) for x in arr[:, 0]]
        # restrict each slice tof to the values between tlow and thigh
        arr = [arr[i, :] for i in np.where(np.logical_and(arr[:, 0] >= tlow, arr[:, 0] <= thigh))]
        reduced_data.append(arr[0])
    # make a unique numpy array out of the list of numpy array; also get rid of empty arrays
    reduced_data = np.vstack(reduced_data)
    # get rid of empty arrays --> not useful as np.vstack already does this job ;)
    # was reduced_data = np.array([y for y in reduced_data if 0 not in y.shape])
    x = np.array(reduced_data[:, 1])
    x = sm.add_constant(x)
    y = np.array(reduced_data[:, 0])
    wls_model = sm.WLS(y, x, weights=np.array([1 / i for i in y]))
    reg_res = wls_model.fit()
    b, a = reg_res.params
    xfit = np.linspace(0, max(reduced_data[:, 1]))
    yfit = eval('{}*xfit+{}'.format(a, b))
    plt.plot(reduced_data[:, 1], reduced_data[:, 0], 'bo', xfit, yfit, 'r', linewidth=3.0)
    plt.savefig(filename + '_tof_drift.pdf')
    plt.savefig(filename + '_tof_drift.eps')
    print reg_res.summary()

def fitslices_x(th2d,nbofbins,xlow,xhigh,filename,cycles):
    """
    implement TH2D::FitSlicesX() routine of ROOT TH2D --> see root manual for details
    :param fit_res:
    :param filename:
    :param xlow: low edge of the subrange to be fitted (usually range of a high statistics reference)
    :param xhigh: upper edge of the subrange to be fitted.
    :return: TObjectArray containing the histograms containing the fit histograms
    """
    temp_th2d = TH2D(th2d)
    temp_th2d.GetXaxis().SetRangeUser(xlow,xhigh)
    slices = TObjArray()
    if nbofbins != 0:
        temp_th2d.RebinY(nbofbins)
    #this will fit with a gaussian (0 means gaus fit), under and over flow bins are including in the fit, only slices whose x projection contains more than 5 filled bins,with "QNR option,
    #store the results in histogram contained in slices
    temp_th2d.FitSlicesX(0,1,cycles,100,"QNR",slices)
    c = TCanvas('c', 'slices', 200, 10, 700, 500)
    slices[1].SetAxisRange(xlow, xhigh,"Y")
    slices[1].Draw()
    #epsfile='__slices.eps'
    pdffile = '_slices.pdf'
    c.SaveAs(filename+pdffile)
    #c.SaveAs(filename+epsfile)
    return slices

def correct_to_peak(xlow,xup,th2d,ref_bin,nb_of_bins,filename,cycles,xbins):
    """
    The spectrum is corrected to the mean tof of the first slice of a well known reference peak
    defined by xlow and xhigh.
    Makes use of function fitslices_x
    :param xlow: low edge of the subrange defining the peak used for correction
    :param xup: upper edge of the subrange defining the peak used as ref for correction
    :param th2d: 2D root hsitogram representing the MCDWIN data
    :param ref_bin: index of the slice used as reference for correction (usually first or last slice)
    :param nb_of_bins: nb of adjacent to combine in the y axes before correcting
    :param filename: output file name
    :param cyles: number of cycles
    :return: corrected TH2D
    """
    slices = fitslices_x(th2d,nb_of_bins,xlow,xup,filename,cycles)
    cor_list = list()
    cor_th2d = TH2D(th2d)
    if ref_bin >= 1 and ref_bin <= cycles:
        for i in xrange(cycles):
            if slices[1].GetBinContent(i+1) != 0.0:
                correction_factor = int(round(slices[1].GetBinContent(i+1)-slices[1].GetBinContent(ref_bin)))
            else:
                correction_factor = 0
                cor_list.append(correction_factor)
    for i,fact in enumerate(cor_list):
        for j in xrange(xbins):
            cor_th2d.SetBinContent(j+1,i+1,th2d.GetBinContent(j+fact,i+1))
    return cor_th2d



def visualise_correlation_matrix(fit_res, filename):
    """
    Build TH2D representing the correlation_matrix between all your parameters
    Print this matrix to eps and pdf
    :param fit_res:
    :param filename:
    :return: None
    """
    gStyle.SetOptStat(0)
    gStyle.SetPalette(1)
    hcorr = fit_res.correlationHist()
    c = TCanvas('c', 'Correlation Matrix', 200, 10, 700, 500)
    hcorr.Draw("colz")
    pdffile = filename + 'cor.pdf'
    #epsfile = filename + '_cor.eps'
    c.SaveAs(pdffile)
    #c.SaveAs(epsfile)


def import_roofit_dataset(th2d,dict887):
    """
    Import TH2D representing as 2d binned roofit dataset
    Function not used at this stage
    :param th2d:
    :param dict887:
    :return: imported 2d RooDataSet , x  and y roorealvars
    """
    # lower edge of first bin
    xlow = float(dict887.get('caloff'))-0.5*float(dict887.get('calfact'))
    # upper edge of last bin not included in last bin
    xup = float(dict887.get('caloff')) + (float(dict887.get('range'))-0.5) * float(dict887.get('calfact'))
    ylow = -0.5
    yup = float(dict887.get('cycles'))-0.5
    x = RooRealVar("x", "x", xlow, xup, "ns")
    y = RooRealVar("y", "y", ylow, yup, "")
    return RooDataHist('roohist_tof', 'title', RooArgList(x,y), RooFit.Import(th2d)) , x, y

def list_887_file(current_path, default_type="887"):
    pattern = current_path + "*." + default_type
    c_887_file = gb.glob(pattern)
    allfile = []
    for eachfile in c_887_file:
        allfile.append(os.path.basename(eachfile))
    return allfile

def list_asc_file(current_path):
    pattern = current_path + "*.asc"
    c_asc_file = gb.glob(pattern)
    allfile = []
    for eachfile in c_asc_file:
        allfile.append(os.path.basename(eachfile))
    return allfile

def get_result(pdf_name, nbparam, frame, rooarglist, cov):
    red_chi2 = frame.chiSquare(nbparam)
    parVal = np.array([])
    for i in xrange(nbparam):
        parVal = np.append(parVal, [rooarglist[i+1].getValV(), rooarglist[i + 1].getError()])
    l = [cov[i][j] for i in xrange(nbparam) for j in xrange(nbparam)]
    return red_chi2, parVal, l

def fill2d_mpant_histo(th2d, data):
    for i in xrange(len(data)):
        th2d.SetBinContent(int(data[i, 0]) + 1, int(data[i, 1]) + 1, int(data[i,2]))
def counts_sum(th1d,tmin,tmax):
    sum = 0.
    for i in range(th1d.FindBin(tmin),th1d.FindBin(tmax)+1):
        sum += th1d.GetBinContent(i)
    return sum

def mrtof_analysis(argv):
    """
    Main function
    :param argv: commad line argument entered by the end user as cf example in header
    :param filename:
    :return: None
    """
    RooRandom.randomGenerator().SetSeed(int(time()))
    print 'current dir = ', os.getcwd()
    c_pat = os.getcwd() + "\\"  # current path

    input_dir = os.getcwd()
    #print input_dir[0], input_dir[1]
    print input_dir
    common_path = "common"
    if not os.path.exists(common_path):
        os.makedirs(common_path)
    f = open(common_path + "/" + "auto-ranges.ini", "w")
    f_res_all = open(common_path + "/" + "all.res", "w", 0)   # "0" means write immediately to files
    f_res_header = "File Name    Range_low    Range_high   Model Name   count  Parameters with Error        Red_Chi2     P-value     Covariance-Matrix \n"
    f_res_all.write(f_res_header)

    file_type = raw_input("please specify the file, 1 for ascii and 2 mpant \n")
    all_887_file = []
    all_asc_file = []
    all_mpant_file = []
    min_range = []
    thres_Val = []
    min_range = float(raw_input("define a fit range:"))
    thres_Val = float(raw_input("define a threshold:"))
    if file_type == "1":
        all_887_file = list_887_file(c_pat) # find all 887 files
        print "list all 887 files", all_887_file

        all_asc_file = list_asc_file(c_pat) # find all asc files
        print all_asc_file

        for each_887_file in all_887_file:
            fname, fextension = os.path.splitext(each_887_file)
            file887 = each_887_file
            print 'file887 = ', file887, os.path.exists(file887)
            datafile = fname + '.asc'
            print 'datafile = ', datafile, os.path.exists(datafile)
            output_dir = fname
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            output_file = output_dir +"/"
            dict887 = load887file(file887)
            # load MCDWIN data file
            ascdict, mcdwin_data = load_asc(datafile)
            c2 = TCanvas()
            hist_mcdwin = create2d_histogram(output_file, dict887)
            fill2d_histo(hist_mcdwin, mcdwin_data)
            hist_mcdwin.Draw("colz")
            Canvas_name = output_file +  fname + "_slice" + ".png"
            c2.SaveAs(Canvas_name)
            total_projectionx = hist_mcdwin.ProjectionX("Total projection x")
            total_projectionx.SetTitle("Total projection x")
            total_projectionx.Draw()
            Pro_jection_name =  output_file + fname + "_projec" + ".png"
            c2.SaveAs(Pro_jection_name)
            tlow = []
            thigh = []
            tlow.append(float(argv[1]))  # default value 0
            thigh.append(float(argv[2])) # default value 0
            combo_bins = int(argv[3])

            nb_peak = 1
    		# first define the peak range, if tlow and thigh are not predefined in ini file
    		# if nb_peak larger than one, tlow[i] and thigh[i] correspond to the ith peak
            if tlow[0] == 0.0 and thigh[0] == 0.0:
                print "range not defined"
                nbins = total_projectionx.GetNbinsX()
                print "total bin range=", nbins
                tlow, thigh, nb_peak = find_peak(total_projectionx, nbins,min_range,thres_Val)
                for i in range(nb_peak):
                    print tlow[i], thigh[i]
                    f.write("%-60s %15.4f %15.4f %1s %10s %20s %20s\n" % (input_dir+"\\"+fname,tlow[i],thigh[i],argv[3],argv[4],argv[5],argv[6]))

    		# if range for correction are passed then apply correction
            if len(argv) == 9:
                tlow_cor = float(argv[7])
                thigh_cor = float(argv[8])
                hist_mcdwin = correct_to_peak(tlow_cor,thigh_cor,hist_mcdwin,1,0,output_file,int(dict887.get('cycles')), int(dict887.get('range')))
            # loop all the peaks
            c2 = TCanvas("c2", "projectionx", 200,10,700,500)

            f_res_name = output_file + fname + ".res"
            f_res = open(f_res_name, "w")
            f_res.write(f_res_header)
            for i_peak in range(nb_peak):
                mean, rms, counts, skewness, histogramme_roi = get_roi(hist_mcdwin, tlow[i_peak], thigh[i_peak])
                print mean, rms, counts, skewness, histogramme_roi
    			#plot_mcdwin(hist_mcdwin,output_file)
    			#plot_yproj(hist_mcdwin,output_file)
    			#plot_xproj(hist_mcdwin,output_file)
    			# create Xprojection of MCDWIN data
                tof_spectrum = projection_x(histogramme_roi)
                peak_th = str(i_peak+1) + "th peak"
                tof_spectrum.SetTitle(peak_th)
                tof_spectrum.Draw()
                t_mark = str(int(mean))
                projpath = output_file+t_mark+".png"
                #print "time mark = ", t_mark
                c2.SaveAs(projpath)
                if combo_bins != 0:
                    rebin_1dhist(tof_spectrum, combo_bins)
                nbins = tof_spectrum.GetNbinsX()
                print nbins

                # create the RooRealVar and model pdf to be fitted
                pdf_def, nb_of_rooreal = get_1d_pdf(argv[4], tlow[i_peak], thigh[i_peak], mean, rms, skewness)
                analysis_pdf = pdf_def[0]
                print analysis_pdf
                t = pdf_def[1]
                print analysis_pdf
                parameters = RooArgList()
                print "Parameters = ", parameters
                for i in xrange(nb_of_rooreal):
                    parameters.add(pdf_def[i + 1])
                roohist_tof = RooDataHist('roohist_tof', 'title', RooArgList(t), RooFit.Import(tof_spectrum))
                # fit the pdf with maximum likelyhood
                result_mlkh = analysis_pdf.fitTo(roohist_tof, RooFit.Minos(1), RooFit.Save(), RooFit.Verbose(False))
                # Print fit results
                result_mlkh.Print("v")

                nb_of_params = result_mlkh.floatParsFinal().getSize()
                print '---------FLOATING PARAMS---', nb_of_params
                # plot fit Result MLKH FIT to TOF data
                tframe = plot_peak(t, analysis_pdf, roohist_tof, output_file, t_mark)
    			# extract COVARIANCE matrix from fit res also visualise the CORELATION matrix as an histogram
                covariance = result_mlkh.covarianceMatrix()
                visualise_correlation_matrix(result_mlkh, output_file)
    			# Perform a series of test to validate the fit
    			# print the residuals of the fit
                make_pull(t, tframe, output_file, t_mark)
    			# compute p value
                p_val = goodness_of_fit(analysis_pdf, roohist_tof, nb_of_params)
                print p_val
                print '------------------'
    			# disabled p_value as it fails for high statistics peaks. We suspect that the issue is related to the fact that
    			# the chi2 distribution is very flat at high ndof. AS the area is conserved we thus have very very small value
    			# that are hard to handle or prone to numerical instabilities
    			# p_val = -1.0
    			# print res to file
                reduce_chi2, parVal, cov_val = get_result(argv[4], nb_of_params, tframe, parameters, covariance)
                print fname, tlow[i_peak]
                f_res.write("%-20s %-12.2f %-12.2f %-10s %-12.3f %-5.4f %-5.4f %-5.4f %-5.4f %-5s %-5.2f %-5.2f %-5.2f %-5.2f  \n" % (fname, tlow[i_peak], thigh[i_peak], argv[4], parVal[0],parVal[1],parVal[2],parVal[3],reduce_chi2, p_val,cov_val[0],cov_val[1],cov_val[2],cov_val[3]))
                f_res_all.write("%-20s %-12.2f %-12.2f %-10s %-12.3f %-5.4f %-5.4f %-5.4f %-5.4f %-5s %-5.2f %-5.2f %-5.2f %-5.2f  \n" % (fname, tlow[i_peak], thigh[i_peak], argv[4], parVal[0],parVal[1],parVal[2],parVal[3],reduce_chi2, p_val,cov_val[0],cov_val[1],cov_val[2],cov_val[3]))
                #printtofile(argv[4],output_file, parameters, t_mark, covariance, nb_of_params, tlow[i_peak], thigh[i_peak], p_val, tframe, fname)
                print parameters[1].getValV(), parameters[1].getError(), parameters[1].getAsymErrorLo()
                print parameters[2].getValV(), parameters[2].getError(), parameters[2].getAsymErrorLo()
    			# create the profile ll
                make_profilell(analysis_pdf, roohist_tof, output_file, t_mark, nb_of_params, parameters)
                print parameters[1].getValV(), parameters[1].getError(), parameters[1].getAsymErrorLo()
                print parameters[2].getValV(), parameters[2].getError(), parameters[2].getAsymErrorLo()

                #ini_par_val = list()
    			#for i in range(nb_of_params):
    			#    ini_par_val.append(parameters[i + 1].getValV())
    			#    # Monte Carlo study
    			#mcstudy(parameters, ini_par_val, nb_of_params, analysis_pdf, nbins, counts, output_file)

    			# resultD_mlkh = analysis_pdf.fitTo(roohist_tof, RooFit.Minimizer('Minuit2', 'simplex'), RooFit.Minos(True),
    			#                                   RooFit.Save(), RooFit.Verbose(False))
    			# resultD_mlkh.Print("v")
            #f_res.close()
        #f_res_all.close()
        #f.close()

    # process mpant files
    if file_type == "2":
        all_mpant_file = list_887_file(c_pat, "mpa")
        config = ConfigParser.RawConfigParser(allow_no_value=True)

        for each_mpant in all_mpant_file:
            pars = {}
            fname, fextension = os.path.splitext(each_mpant)
            print "fname", fname, "extension", fextension
            output_dir = fname
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            output_file = output_dir +"/"

            begin_line = []
            lookup = 'TDAT'
            nbLine = 0
            with open(each_mpant) as myFile:
                for num, line in enumerate(myFile, 1):      # save histogram + header in list
                    if lookup in line:  # loacate the oistion of data
                        begin_line.append(num)
                nbLine = num  # nbLine should be a list?
            #print each_mpant, "beginline", begin_line, "nb of lines", nbLine

            config.readfp(open('%s' % (each_mpant)))          # read MPANT file to process it
            bin_range = config.getint('MPA4A', 'range')     # get bin range
            nb_cycles = config.getint('MPA4A', 'cycles')
            pars['range'] = int(bin_range)
            pars['cycles'] = int(nb_cycles)

            #print "bin_range", bin_range, "cycle=", nb_cycles

            which_begin_line = 0    # sometimes 2 channels are active. here I decide which histogram to fit
            # offset = 0
            terminate = ''
            time_offset = 0 # time offset (integer) in ns unit
            calfact = 0.    # bin width in ns

            for j in ['CHN1', 'CHN2', 'CHN3', 'CHN4', 'CHN5', 'CHN6']:
                if config.getint(j, 'active') == 1 and config.getint(j, 'TOTALSUM') > 10:
                    time_offset = float(config.getfloat(j, 'caloff'))      # in nanoseconds
                    calfact = float(config.getfloat(j, 'calfact'))
                    pars['caloff'] = time_offset
                    pars['calfact'] = calfact
                elif config.getint(j, 'active') == 1 and time_offset == 0:
                    which_begin_line += 1
                    if which_begin_line == len(begin_line):
                        terminate = 'now'

            #print "timeoffset", time_offset, "calfact", calfact, terminate

            load_file = []
            if terminate == '':
                histogram_data = []
                #histogram_data = [['Time of flight (ns)', 'count', 'cycle']]
                #histogram_data_reduced = [['Ch. / 100 ps', '#']]
                with open(each_mpant, 'rb') as infile:
                    load_file = [[str(h) for h in line.strip().split()] for line in infile]

                #print "last data", nbLine, load_file[nbLine-1]

                maxi = 0
                cc = 0
                nbCycle = 1

                for k in range(begin_line[which_begin_line], nbLine, 1):  # get histogram data from file (not possible with configparser with the non-standard mpant files -.-)
                    cc += 1
                    if cc > bin_range:  # go to the next cycle
                        nbCycle += 1
                        cc = 1

                    help_load_file = [k - begin_line[which_begin_line] + 1 - ((nbCycle-1)*bin_range)]
                    #help_load_file = [(float(k - begin_line[which_begin_line] + 1 - (nbCycle-1.)*bin_range))*calfact + time_offset]
                    #print k, load_file[k], cc, nbCycle

                    help_load_file.extend([nbCycle])
                    help_load_file.extend([float(l) for l in load_file[k]])
                    histogram_data.append(help_load_file)
                #print "ddd", pars.get('caloff'), pars.get('calfact'), pars.get('range'), pars.get('cycles')

                c2 = TCanvas()
                hist_mcdwin_mpant = create2d_histogram(output_file, pars)
                for kk in xrange(len(histogram_data)):
                    hist_mcdwin_mpant.SetBinContent(histogram_data[kk][0]+1, histogram_data[kk][1]+1, histogram_data[kk][2])

                hist_mcdwin_mpant.Draw("colz")
                Canvas_name = output_file +  fname + "_slice" + ".png"
                c2.SaveAs(Canvas_name)
                total_projectionx = hist_mcdwin_mpant.ProjectionX("Total projection x")
                total_projectionx.SetTitle("Total projection x")
                total_projectionx.Draw()
                Pro_jection_name =  output_file + fname + "_projec" + ".png"
                c2.SaveAs(Pro_jection_name)

                nb_peak = 1
                tlow = []
                thigh = []
                tlow.append(float(argv[1]))  # default value 0
                thigh.append(float(argv[2])) # default value 0
                combo_bins = int(argv[3])

                if tlow[0] == 0.0 and thigh[0] == 0.0:
                    print "range not defined"
                    nbins = total_projectionx.GetNbinsX()
                    print "total bin range=", nbins
                    tlow, thigh, nb_peak = find_peak(total_projectionx, nbins,min_range,thres_Val)
                    for i in range(nb_peak):
                        print tlow[i], thigh[i]
                        f.write("%-60s %15.4f %15.4f %1s %10s %20s %20s\n" % (input_dir+"\\"+fname,tlow[i],thigh[i],argv[3],argv[4],argv[5],argv[6]))
                if len(argv) == 9:
                    tlow_cor = float(argv[7])
                    thigh_cor = float(argv[8])
                    hist_mcdwin_mpant = correct_to_peak(tlow_cor,thigh_cor,hist_mcdwin_mpant,1,0,output_file,int(dict887.get('cycles')), int(dict887.get('range')))

                c2 = TCanvas("c2", "projectionx", 200,10,700,500)
                f_res_name = output_file + fname + ".res"
                f_res = open(f_res_name, "w")
                f_res.write(f_res_header)
                for i_peak in range(nb_peak):
                    mean, rms, counts, skewness, histogramme_roi = get_roi(hist_mcdwin_mpant, tlow[i_peak], thigh[i_peak])
                    print mean, rms, counts, skewness, histogramme_roi
                    tof_spectrum = projection_x(histogramme_roi)
                    peak_th = str(i_peak+1) + "th peak"
                    tof_spectrum.SetTitle(peak_th)
                    tof_spectrum.Draw()
                    t_mark = str(int(mean))
                    projpath = output_file+t_mark+".png"
                    c2.SaveAs(projpath)

                    if combo_bins != 0:
                        rebin_1dhist(tof_spectrum, combo_bins)
                    nbins = tof_spectrum.GetNbinsX()
                    print nbins

                    pdf_def, nb_of_rooreal = get_1d_pdf(argv[4], tlow[i_peak], thigh[i_peak], mean, rms, skewness)
                    analysis_pdf = pdf_def[0]
                    print analysis_pdf
                    t = pdf_def[1]
                    print analysis_pdf
                    parameters = RooArgList()
                    print "Parameters = ", parameters
                    for i in xrange(nb_of_rooreal):
                        parameters.add(pdf_def[i + 1])
                    roohist_tof = RooDataHist('roohist_tof', 'title', RooArgList(t), RooFit.Import(tof_spectrum))
                    # fit the pdf with maximum likelyhood
                    result_mlkh = analysis_pdf.fitTo(roohist_tof, RooFit.Minos(1), RooFit.Save(), RooFit.Verbose(False))
                    # Print fit results
                    result_mlkh.Print("v")
                    i_count = counts_sum(tof_spectrum,tlow[i_peak],thigh[i_peak])

                    nb_of_params = result_mlkh.floatParsFinal().getSize()
                    print '---------FLOATING PARAMS---', nb_of_params
                    # plot fit Result MLKH FIT to TOF data
                    tframe = plot_peak(t, analysis_pdf, roohist_tof, output_file, t_mark)
		        	# extract COVARIANCE matrix from fit res also visualise the CORELATION matrix as an histogram
                    covariance = result_mlkh.covarianceMatrix()
                    visualise_correlation_matrix(result_mlkh, output_file)
			        # Perform a series of test to validate the fit
			        # print the residuals of the fit
                    make_pull(t, tframe, output_file, t_mark)
			        # compute p value
                    p_val = goodness_of_fit(analysis_pdf, roohist_tof, nb_of_params)
                    print p_val
                    print '------------------'
			        # disabled p_value as it fails for high statistics peaks. We suspect that the issue is related to the fact that
			        # the chi2 distribution is very flat at high ndof. AS the area is conserved we thus have very very small value
			        # that are hard to handle or prone to numerical instabilities
			        # p_val = -1.0
			        # print res to file
                    reduce_chi2, parVal, cov_val = get_result(argv[4], nb_of_params, tframe, parameters, covariance)
                    print fname, tlow[i_peak]
                    f_res.write("%-10s %-12.2f %-12.2f %-10s %6d %-12.3f %-5.4f %-5.4f %-5.4f %-5.4f %-5s %-5.2f %-5.2f %-5.2f %-5.2f  \n" % (fname, tlow[i_peak], thigh[i_peak], argv[4], i_count, parVal[0],parVal[1],parVal[2],parVal[3],reduce_chi2, p_val,cov_val[0],cov_val[1],cov_val[2],cov_val[3]))
                    f_res_all.write("%-10s %-12.2f %-12.2f %-10s %6d %-12.3f %-5.4f %-5.4f %-5.4f %-5.4f %-5s %-5.2f %-5.2f %-5.2f %-5.2f  \n" % (fname, tlow[i_peak], thigh[i_peak], argv[4], i_count, parVal[0],parVal[1],parVal[2],parVal[3],reduce_chi2, p_val,cov_val[0],cov_val[1],cov_val[2],cov_val[3]))
                    #printtofile(argv[4],output_file, parameters, t_mark, covariance, nb_of_params, tlow[i_peak], thigh[i_peak], p_val, tframe, fname)
                    print parameters[1].getValV(), parameters[1].getError(), parameters[1].getAsymErrorLo()
                    print parameters[2].getValV(), parameters[2].getError(), parameters[2].getAsymErrorLo()
			        # create the profile ll
                    make_profilell(analysis_pdf, roohist_tof, output_file, t_mark, nb_of_params, parameters)
                    print parameters[1].getValV(), parameters[1].getError(), parameters[1].getAsymErrorLo()
                    print parameters[2].getValV(), parameters[2].getError(), parameters[2].getAsymErrorLo()

                    #ini_par_val = list()
			        #for i in range(nb_of_params):
			        #    ini_par_val.append(parameters[i + 1].getValV())
			        #    # Monte Carlo study
			        #mcstudy(parameters, ini_par_val, nb_of_params, analysis_pdf, nbins, counts, output_file)

			        # resultD_mlkh = analysis_pdf.fitTo(roohist_tof, RooFit.Minimizer('Minuit2', 'simplex'), RooFit.Minos(True),
			        #                                   RooFit.Save(), RooFit.Verbose(False))
			        # resultD_mlkh.Print("v")

    f_res.close()
    f.close()
    f_res_all.close()



if __name__ == "__main__":
    mrtof_analysis(sys.argv[1:])
