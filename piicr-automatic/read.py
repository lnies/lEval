import numpy as np
import sys

def main(argv):
    ini = open("generators.ini","r")
    key_list, val_list = np.loadtxt(ini, delimiter='=', usecols=(0, 1), dtype='object', unpack=True)
    devdict = {}
    for arg in range(len(key_list)):
        devdict[key_list[arg]] = val_list[arg]
    add_mag = devdict['gen_mag_excitation']
    add_p1 = devdict['gen_pattern1']
    add_p2 = devdict['gen_pattern2']
    add_ej = devdict['gen_ejection']
    print add_mag
    print add_p1
    print add_p2
    print add_ej
    cfg = open(argv[0],"r")
    pars = np.loadtxt(cfg,skiprows=1, dtype = 'object')
    print pars
    elements = pars[:,0]
    f_mag =  np.around(pars[:,1].astype(np.float),3)
    amp_mag = pars[:,2].astype(np.float)
    n_mag = pars[:,3].astype(np.int)
    f_plus = np.around(pars[:,4].astype(np.float),3)
    amp_plus = pars[:,5].astype(np.float)
    n_plus = pars[:,6].astype(np.int)
    f_c = np.around(pars[:,7].astype(np.float),3)
    amp_c = pars[:,8].astype(np.float)
    n_c = pars[:,9].astype(np.int)
    n_acc = pars[:,10].astype(np.int)
    phase_acc_time = np.around(n_acc/f_c,9)
    rec_dat = np.vstack((elements,f_mag,amp_mag,n_mag,f_plus,amp_plus,n_plus,f_c,amp_c,n_c,n_acc,phase_acc_time)).transpose()
    nloops = argv[1]
    #print elements
    print f_mag
    print amp_mag
    print n_mag
    print f_plus
    print amp_plus
    print n_plus
    print f_c
    print amp_c
    print n_c
    print n_acc
    print nloops
    np.savetxt("params.rec",rec_dat, fmt = '%s', header = 'Element nu- nu-_pulse_ampl nu-_pulse_cycl nu+ nu+_pulse_ampl nu+_pulse_cycles nuc pi_pulse_ampl pi_pulse_cycl phase_acc_cycl phase_acc_time')
    cfg.close()
    ini.close()

if __name__=="__main__":
    main(sys.argv[1:])
