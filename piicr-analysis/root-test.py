# ---------------------------------------------------------------------------
# Written by Jonas Karthein in 2016/2017. Questions to jonas.karthein@cern.ch
# very helpful: https://accserv.lepp.cornell.edu/svn/packages/root/documentation/users-guide/FittingHistograms.md
# ---------------------------------------------------------------------------
from read_write_functions import *
import time
import ROOT as root
import os


# def spot_fit_root(file_name, spot_positions, all_peak_positions, mean_x, mean_y, cut_range):
#     '''
#     Function to fit position to 2D unbinned csv file.

#     :param cut_range: vector containing [0] = cut_x_min, [1] = cut_x_max, [2] = cut_y_min, [3] = cut_y_max
#     '''
#     multi_process_dict = {}
#     timestamp = (int(time.time()*1000))
#     print timestamp
#     # multi_process_dict[str(timestamp)] = 0      # canvas
#     # multi_process_dict[str(timestamp+1)] = 0    # tree
#     # multi_process_dict[str(timestamp+2)] = 0    # inj
#     # multi_process_dict[str(timestamp+3)] = 0    # x
#     # multi_process_dict[str(timestamp+4)] = 0    # y
#     # multi_process_dict[str(timestamp+5)] = 0    # t
#     # multi_process_dict[str(timestamp+5)] = 0    # ds
#     # multi_process_dict[str(timestamp+6)] = 0    # meanx
#     # multi_process_dict[str(timestamp+7)] = 0    # sigmax
#     # multi_process_dict[str(timestamp+8)] = 0    # meany
#     # multi_process_dict[str(timestamp+9)] = 0    # sigmay
#     # multi_process_dict[str(timestamp+10)] = 0    # gaussx
#     # multi_process_dict[str(timestamp+11)] = 0    # gaussy
#     # multi_process_dict[str(timestamp+12)] = 0    # gaussxy
#     # multi_process_dict[str(timestamp+13)] = 0    # result
#     # multi_process_dict[str(timestamp+14)] = 0    # framex
#     # multi_process_dict[str(timestamp+15)] = 0    # framey
#     # multi_process_dict[str(timestamp+16)] = 0    # nllx
#     # multi_process_dict[str(timestamp+17)] = 0    # profile_llmeanx
#     # multi_process_dict[str(timestamp+18)] = 0    # pllframex
#     # multi_process_dict[str(timestamp+19)] = 0    # nlly
#     # multi_process_dict[str(timestamp+20)] = 0    # profile_llmeany
#     # multi_process_dict[str(timestamp+21)] = 0    # pllframey

#     # import data into TTree
#     multi_process_dict.update(str(timestamp+1)= root.TTree(str(multi_process_dict[str(timestamp+1)]), 'tree'))
#     multi_process_dict[str(timestamp+1)].ReadFile('%s_spot_positions.csv' % (file_name), 'inj/D:x/D:y/D:t/D', ',')

#     # set vector properties
#     inj = root.RooRealVar('inj', 'inj', 0, 0.0, 5000.0)
#     multi_process_dict[str(timestamp+3)] = root.RooRealVar(str(multi_process_dict[str(timestamp+3)]), 'x', 0, cut_range[0], cut_range[1])
#     multi_process_dict[str(timestamp+4)] = root.RooRealVar(str(multi_process_dict[str(timestamp+4)]), 'y', 0, cut_range[2], cut_range[3])
#     t = root.RooRealVar('t', 't', 1000000.0, 1000000., 3000000.0)

#     # set dataset
#     multi_process_dict[str(timestamp+5)] = root.RooDataSet('x-y', 'x-y-data', multi_process_dict[str(timestamp+1)], root.RooArgSet(multi_process_dict[str(timestamp+3)], multi_process_dict[str(timestamp+4)]))

#     # set fitting model
#     multi_process_dict[str(timestamp+6)] = root.RooRealVar(str(multi_process_dict[str(timestamp+6)]), 'meanx', mean_x, mean_x-50.0, mean_x+50.0)
#     multi_process_dict[str(timestamp+7)] = root.RooRealVar(str(multi_process_dict[str(timestamp+7)]), 'sigmax', 50.0, 10.0, 100.0)
#     multi_process_dict[str(timestamp+8)] = root.RooRealVar(str(multi_process_dict[str(timestamp+8)]), 'meany', mean_y, mean_y-50.0, mean_y+50.0)
#     multi_process_dict[str(timestamp+9)] = root.RooRealVar(str(multi_process_dict[str(timestamp+9)]), 'sigmay', 50.0, 10.0, 100.0)

#     multi_process_dict[str(timestamp+10)] = root.RooGaussian(str(multi_process_dict[str(timestamp+10)]), 'Gaussian distribution', multi_process_dict[str(timestamp+3)], multi_process_dict[str(timestamp+6)], multi_process_dict[str(timestamp+7)])
#     multi_process_dict[str(timestamp+11)] = root.RooGaussian(str(multi_process_dict[str(timestamp+11)]), 'Gaussian distribution', multi_process_dict[str(timestamp+4)], multi_process_dict[str(timestamp+8)], multi_process_dict[str(timestamp+9)])
#     multi_process_dict[str(timestamp+12)] = root.RooProdPdf('gxy', 'gx*gy', root.RooArgList(multi_process_dict[str(timestamp+10)],multi_process_dict[str(timestamp+11)]))

#     # actual fit:
#     multi_process_dict[str(timestamp+13)] = multi_process_dict[str(timestamp+12)].fitTo(multi_process_dict[str(timestamp+5)]) # , root.RooFit.NumCPU(4))

#     # plot projections and fits
#     multi_process_dict[str(timestamp)] = root.TCanvas(str(multi_process_dict[str(timestamp)]), 'results', 1000, 700)
#     multi_process_dict[str(timestamp)].Divide(2, 2)

#     multi_process_dict[str(timestamp)].cd(1)
#     multi_process_dict[str(timestamp+14)] = multi_process_dict[str(timestamp+3)].frame()
#     multi_process_dict[str(timestamp+5)].plotOn(multi_process_dict[str(timestamp+14)])
#     multi_process_dict[str(timestamp+12)].plotOn(multi_process_dict[str(timestamp+14)])
#     multi_process_dict[str(timestamp+10)].paramOn(multi_process_dict[str(timestamp+14)])
#     multi_process_dict[str(timestamp+14)].Draw()
#     multi_process_dict[str(timestamp+14)].SetTitle('X-Proj. & unbinned 2D-max.likelihood fit')
#     multi_process_dict[str(timestamp+14)].SetXTitle('channels')

#     multi_process_dict[str(timestamp)].cd(2)
#     multi_process_dict[str(timestamp+15)] = multi_process_dict[str(timestamp+4)].frame()
#     multi_process_dict[str(timestamp+5)].plotOn(multi_process_dict[str(timestamp+15)])
#     multi_process_dict[str(timestamp+12)].plotOn(multi_process_dict[str(timestamp+15)])
#     multi_process_dict[str(timestamp+11)].paramOn(multi_process_dict[str(timestamp+15)])
#     multi_process_dict[str(timestamp+15)].Draw()
#     multi_process_dict[str(timestamp+15)].SetTitle('Y-Proj. & unbinned 2D-max.likelihood fit')
#     multi_process_dict[str(timestamp+15)].SetXTitle('channels')

#     # plot nll's
#     multi_process_dict[str(timestamp)].cd(3)
#     multi_process_dict[str(timestamp+16)] = multi_process_dict[str(timestamp+12)].createNLL(multi_process_dict[str(timestamp+5)]) # , root.RooFit.NumCPU(4))
#     multi_process_dict[str(timestamp+17)] = multi_process_dict[str(timestamp+16)].createProfile(root.RooArgSet(multi_process_dict[str(timestamp+6)]))
#     multi_process_dict[str(timestamp+18)] = multi_process_dict[str(timestamp+6)].frame()
#     multi_process_dict[str(timestamp+16)].plotOn(multi_process_dict[str(timestamp+18)], root.RooFit.ShiftToZero())
#     multi_process_dict[str(timestamp+17)].plotOn(multi_process_dict[str(timestamp+18)], root.RooFit.LineColor(root.kRed))
#     multi_process_dict[str(timestamp+18)].SetTitle('X-NLL analysis')
#     multi_process_dict[str(timestamp+18)].SetXTitle('channels')
#     multi_process_dict[str(timestamp+18)].Draw()


#     multi_process_dict[str(timestamp)].cd(4)
#     multi_process_dict[str(timestamp+19)] = multi_process_dict[str(timestamp+12)].createNLL(multi_process_dict[str(timestamp+5)]) # , root.RooFit.NumCPU(4))
#     multi_process_dict[str(timestamp+20)] = multi_process_dict[str(timestamp+19)].createProfile(root.RooArgSet(multi_process_dict[str(timestamp+8)]))
#     multi_process_dict[str(timestamp+21)] = multi_process_dict[str(timestamp+8)].frame()
#     multi_process_dict[str(timestamp+19)].plotOn(multi_process_dict[str(timestamp+21)], root.RooFit.ShiftToZero())
#     multi_process_dict[str(timestamp+20)].plotOn(multi_process_dict[str(timestamp+21)], root.RooFit.LineColor(root.kRed))
#     multi_process_dict[str(timestamp+21)].SetTitle('Y-NLL analysis')
#     multi_process_dict[str(timestamp+21)].SetXTitle('channels')
#     multi_process_dict[str(timestamp+21)].Draw()
#     multi_process_dict[str(timestamp)].Update()

#     pdffile = file_name + '_x-y-projection-fits.pdf'
#     multi_process_dict[str(timestamp)].SaveAs(pdffile)

#     # fit results
#     x_pos = multi_process_dict[str(timestamp+6)].getValV()
#     x_pos_err = multi_process_dict[str(timestamp+6)].getError()
#     y_pos = multi_process_dict[str(timestamp+8)].getValV()
#     y_pos_err = multi_process_dict[str(timestamp+8)].getError()
#     x_sigma = multi_process_dict[str(timestamp+7)].getValV()
#     x_sigma_err = multi_process_dict[str(timestamp+7)].getError()
#     y_sigma = multi_process_dict[str(timestamp+9)].getValV()
#     y_sigma_err = multi_process_dict[str(timestamp+9)].getError()

#     print 'here'
#     only_spot_positions = [['X position / ch.', 'Y position / ch.']]
#     for i in range(len(spot_positions)):
#         only_spot_positions.append([float(spot_positions[i][1]), float(spot_positions[i][2])])
#     python_plot(only_spot_positions, 'MCP-PS', '%s_mcp' % file_name, 'no', '', 'no', 'gauss', 'Utopia', 'red', 'black', 'full', 0, 8, 6, 'mcp', 60, 0.2, 'on', 'off', 'mcp', [x_sigma, x_sigma_err, y_sigma, y_sigma_err], [x_pos, x_pos_err, y_pos, y_pos_err])

#     #root.gApplication.Run()


#     # output
#     peak_position = [file_name, x_pos, x_pos_err, y_pos, y_pos_err, x_sigma, x_sigma_err, y_sigma, y_sigma_err, multi_process_dict[str(timestamp+5)].numEntries()]
#     all_peak_positions.append(peak_position)

#     return(all_peak_positions)

def spot_fit_root(file_name, spot_positions, all_peak_positions, mean_x, mean_y, cut_range):
    '''
    Function to fit position to 2D unbinned csv file.

    :param cut_range: vector containing [0] = cut_x_min, [1] = cut_x_max, [2] = cut_y_min, [3] = cut_y_max
    '''

    # import data into TTree
    tree = root.TTree('tree', 'tree')
    tree.ReadFile('%s_spot_positions.csv' % (file_name), 'inj/D:x/D:y/D:t/D', ',')

    # set vector properties
    inj = root.RooRealVar('inj', 'inj', 0, 0.0, 5000.0)
    x = root.RooRealVar('x', 'x', 0, cut_range[0], cut_range[1])
    y = root.RooRealVar('y', 'y', 0, cut_range[2], cut_range[3])
    t = root.RooRealVar('t', 't', 1000000.0, 1000000., 3000000.0)

    # set dataset
    ds = root.RooDataSet('x-y', 'x-y-data', tree, root.RooArgSet(x, y))

    # # set fitting model
    # meanx = root.RooRealVar('meanx', 'meanx', mean_x, mean_x-50.0, mean_x+50.0)
    # sigmax = root.RooRealVar('sigmax', 'sigmax', 50.0, 10.0, 100.0)
    # meany = root.RooRealVar('meany', 'meany', mean_y, mean_y-50.0, mean_y+50.0)
    # sigmay = root.RooRealVar('sigmay', 'sigmay', 50.0, 10.0, 100.0)

    # gaussx = root.RooGaussian('gaussx', 'Gaussian distribution', x, meanx, sigmax)
    # gaussy = root.RooGaussian('gaussy', 'Gaussian distribution', y, meany, sigmay)
    # gaussxy = root.RooProdPdf('gxy', 'gx*gy', root.RooArgList(gaussx,gaussy))

    # # actual fit:
    # result = gaussxy.fitTo(ds) # , root.RooFit.NumCPU(4))

    # # plot projections and fits
    # c = root.TCanvas('c', 'results', 1000, 700)
    # c.Divide(2, 2)

    # c.cd(1)
    # framex = x.frame()
    # ds.plotOn(framex)
    # gaussxy.plotOn(framex)
    # gaussx.paramOn(framex)
    # framex.Draw()
    # framex.SetTitle('X-Proj. & unbinned 2D-max.likelihood fit')
    # framex.SetXTitle('channels')

    # c.cd(2)
    # framey = y.frame()
    # ds.plotOn(framey)
    # gaussxy.plotOn(framey)
    # gaussy.paramOn(framey)
    # framey.Draw()
    # framey.SetTitle('Y-Proj. & unbinned 2D-max.likelihood fit')
    # framey.SetXTitle('channels')

    # # plot nll's
    # c.cd(3)
    # nllx = gaussxy.createNLL(ds) # , root.RooFit.NumCPU(4))
    # profile_llmeanx = nllx.createProfile(root.RooArgSet(meanx))
    # pllframex = meanx.frame()
    # nllx.plotOn(pllframex, root.RooFit.ShiftToZero())
    # profile_llmeanx.plotOn(pllframex, root.RooFit.LineColor(root.kRed))
    # pllframex.SetTitle('X-NLL analysis')
    # pllframex.SetXTitle('channels')
    # # pllframex.SetMinimum(0)
    # # pllframex.SetMaximum(3)
    # pllframex.Draw()


    # c.cd(4)
    # nlly = gaussxy.createNLL(ds) # , root.RooFit.NumCPU(4))
    # profile_llmeany = nlly.createProfile(root.RooArgSet(meany))
    # pllframey = meany.frame()
    # nlly.plotOn(pllframey, root.RooFit.ShiftToZero())
    # profile_llmeany.plotOn(pllframey, root.RooFit.LineColor(root.kRed))
    # pllframey.SetTitle('Y-NLL analysis')
    # pllframey.SetXTitle('channels')
    # # pllframey.SetMinimum(0)
    # # pllframey.SetMaximum(3)
    # pllframey.Draw()
    # c.Update()

    # pdffile = file_name + '_x-y-projection-fits.pdf'
    # c.SaveAs(pdffile)

    c2 = root.TCanvas('c2', 'results2', 1000, 500)
    c2.Divide(2, 1)

    c2.cd(1)
    hist = ds.createHistogram(x, y, 60, 60)
    hist.Draw()
    c2.Update()


    # print CUTG

    # c2.cd(2)
    # hist.Draw('CUTG')
    # c2.Update()


    root.gApplication.Run()



    # # fit results
    # x_pos = meanx.getValV()
    # x_pos_err = meanx.getError()
    # y_pos = meany.getValV()
    # y_pos_err = meany.getError()
    # x_sigma = sigmax.getValV()
    # x_sigma_err = sigmax.getError()
    # y_sigma = sigmay.getValV()
    # y_sigma_err = sigmay.getError()

    # os.exit()

    # print 'here'
    # only_spot_positions = [['X position / ch.', 'Y position / ch.']]
    # for i in range(len(spot_positions)):
    #     only_spot_positions.append([float(spot_positions[i][1]), float(spot_positions[i][2])])
    # python_plot(only_spot_positions, 'MCP-PS', '%s_mcp' % file_name, 'no', '', 'no', 'gauss', 'Utopia', 'red', 'black', 'full', 0, 8, 6, 'mcp', 60, 0.2, 'on', 'off', 'mcp', [x_sigma, x_sigma_err, y_sigma, y_sigma_err], [x_pos, x_pos_err, y_pos, y_pos_err])



    # # output
    # peak_position = [file_name, x_pos, x_pos_err, y_pos, y_pos_err, x_sigma, x_sigma_err, y_sigma, y_sigma_err, ds.numEntries()]
    # all_peak_positions.append(peak_position)

    return(all_peak_positions)


if __name__ == '__main__':
    spot_positions = read_csv('G:\\Experiments\\ISOLTRAP\\Software\\PI-ICR\\Python-DAQ\\x_in_preparation\\root-tests', '129gCd_007_p2_spot_positions')
    spot_fit_root('129gCd_007_p2', spot_positions, [], 9., -30., [-750, 750, -750, 750])
