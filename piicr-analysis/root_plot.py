import ROOT as root
from array import array


def root_test(file_name):
    tree = root.TTree('tree', 'tree')
    tree.ReadFile('%s.csv' % (file_name), 'inj/D:x/D:y/D:t/D', ',')

    # set vector properties
    inj = root.RooRealVar('inj', 'inj', 0, 0.0, 5000.0)
    x = root.RooRealVar('x', 'x', 0, fit_range[0], fit_range[1])
    y = root.RooRealVar('y', 'y', 0, fit_range[2], fit_range[3])
    t = root.RooRealVar('t', 't', 1000000.0, 1000000., 3000000.0)

    # set dataset
    ds = root.RooDataSet('x-y', 'x-y-data', tree, root.RooArgSet(x, y))


    # plot
    c = root.TCanvas('c', 'results', 1000, 1000)
    g = ROOT.TGraph(len(x), x,y)
    g.Draw()

    root.gApplication.Run()

def root_test_2(list_x):

    '''Function takes 1D list with unbinned data and returns a list with gaussian fit parameters: x, x_unc., sigma, sigma_unc.'''

    load_tree = 'from_tree'



    import_list_x = [x for x in list_x]

    import os
    os.chdir('/Users/jonaskarthein/cernbox/Python/Doktor/2017-11-06_root_TCutG_tests')

    tree = root.TTree( 'tree', 'tree' )


    if load_tree == 'from_list':
        x = array('d', [ 0. ])
        tree.Branch('x', x, 'x/D')
        y = array('d', [ 0. ])
        tree.Branch('y', y, 'y/D')
        for i in range(len(import_list_x)):
            x[0] = import_list_x[i][0]
            y[0] = import_list_x[i][1]
            tree.Fill()

    elif load_tree == 'from_tree':
        file = root.TFile('85Rb_004.root', 'READ')
        tree = file.Get('p2').Clone()
        tree.Print()


    x_lims = [-750,750]
    y_lims = [-750,750]

    x = root.RooRealVar('x', 'x', x_lims[0], x_lims[1])#, 'us')
    y = root.RooRealVar('y', 'y', y_lims[0], y_lims[1])#, 'us')


    # simulating manual cut
    mycut = root.TCutG('CUTG',5)        # 'CUTG' is the default name root gives a cut
    mycut.SetPoint(0,-200,-200)
    mycut.SetPoint(1,200,-200)
    mycut.SetPoint(2,200,50)
    mycut.SetPoint(3,-200,50)
    mycut.SetPoint(4,-200,-200)
    print('\nSimulated cut:\n')
    mycut.Print()


    # load on the fly created TCutG object (instead of loading cut, take the just created one):
    mycutg = root.TCutG()
    mycutg = root.gROOT.GetListOfSpecials().FindObject('CUTG')
    # mycutg.SetName('mycutg')  # one could rename the cut, but then adjust the loading function
    mycutg.SetVarX('x')                  # important
    mycutg.SetVarY('y')                  # important
    print('\nRead cut from memory:\n')
    mycutg.Print()


    # save cut in file
    f = root.TFile('save_cut.root', 'recreate')
    mycut.Write()
    f.Close()


    # read cut from file
    g = root.TFile('save_cut.root', 'READ')
    cut = root.TCutG()
    cut = g.Get('CUTG').Clone()#'mycutg').Clone()
    print('\nRead cut from file:\n')
    cut.Print()


    # set dataset
    ds = root.RooDataSet('x-y', 'x-y-data', tree, root.RooArgSet(x, y))
    c = root.TCanvas('c', 'results', 1000, 1000)
    htemp = root.TH2F('htemp', 'htemp', 1000, -750, 750, 1000, -750, 750)


    # plot cut data
    tree.Draw('x:y>>htemp', mycut.GetName())
    htemp.SetMarkerStyle(8)
    htemp.SetMarkerSize(0.5)
    htemp.SetTitle('PI-ICR')
    root.gPad.Update()
    c.SaveAs('tree_plot.pdf')


    # create cut tree:
    tree_cut = root.TNtuple()
    tree_cut = tree.CopyTree(mycut.GetName())
    print('\nOriginal tree:\n')
    tree.Print()
    print('\nCut tree:\n')
    tree_cut.Print()


    # show plot
    root.gApplication.Run()


if __name__ == '__main__':
    # fit_range = [-750,750,-750,750]
    # nll_plot = False
    # file_name = '129mCd_002_p2_spot_positions'
    # root_test(file_name)
    list_x = [[1.0,-130.0,-36.5,2402457.0],[3.0,151.0,186.0,2482932.0],[5.0,239.0,220.5,2486431.0],[8.0,186.5,16.5,2485881.0],[9.0,17.5,209.5,2474744.0],[25.0,179.5,222.0,2485512.0],[30.0,196.0,156.5,2469798.0],[31.0,258.5,137.0,2474085.0],[35.0,299.0,193.0,2479816.0],[36.0,267.0,200.0,2478193.0],[40.0,-113.0,219.0,2475943.0],[53.0,0.0,154.5,2496631.0],[58.0,186.0,-38.0,2490658.0],[67.0,188.5,135.5,2504323.0],[69.0,-117.5,141.0,2489333.0],[73.0,216.5,219.0,2473443.0],[84.0,183.5,138.0,2496904.0],[97.0,-94.0,159.5,2466104.0],[99.0,-80.5,39.0,2489934.0],[109.0,203.0,182.0,2490617.0],[110.0,209.0,136.0,2494844.0],[111.0,289.5,208.5,2486257.0],[114.0,170.0,186.0,2489714.0],[118.0,325.0,103.0,2500314.0],[120.0,-95.0,139.5,2496229.0],[121.0,-85.0,166.0,2489965.0],[123.0,-131.5,213.5,2479931.0],[130.0,292.0,232.5,2491278.0],[131.0,211.5,218.5,2475633.0],[132.0,-153.0,175.0,2475701.0],[136.0,164.0,200.5,2490194.0],[139.0,144.5,125.5,2488658.0],[149.0,-153.0,189.0,2472362.0],[150.0,-151.5,179.5,2496810.0],[151.0,264.0,140.0,2471638.0],[156.0,189.5,192.5,2490653.0],[157.0,267.5,189.0,2481933.0],[226.0,218.0,200.5,2483474.0]]
    root_test_2(list_x)
