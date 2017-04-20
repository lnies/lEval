# -*- coding: utf-8 -*-
"""
Created on Thurs 16 Feb 2017

@author: mamougeo
"""
###################
# Global variables
###################

###################
from time import sleep
import win32com.client
import pydim as dim
import csdatalib as cs
import numpy as np
import scipy.optimize as scop
import time
import sys
import signal

def prepare_data(filename,guiID,LVconfig,LVUserdata):
    """
    Pack the data to be send to the CS
    :param filename: name of the file to which the data will be recorded
    :param guiID: guiID (cf mm8.ini)
    :param LVconfig: Labview config file (see mm8.ini) --> file config.xml
    :param LVUserdata: Labview userdata file (see mm8.ini) --> file userdata.xml
    :return: data ready to be send to the CS via DIM
    """
    try:
        userdata = open(LVUserdata).read()
        config = open(LVconfig).read()
    except Exception as err:
        print err
    else:
        mydata_gui = cs.pack_val(guiID, cs.CHAR)
        mydata_strt = cs.pack_val(guiID, cs.CHAR) + cs.pack_val(filename, cs.CHAR)
        mydata_conf = cs.pack_val(guiID, cs.CHAR) + cs.pack_val(config, cs.CHAR) + cs.pack_val(userdata, cs.CHAR)
        return userdata, config, mydata_conf, mydata_gui, mydata_strt


def connect(dim_dns_node):
    """
    Attempt to connect to dim dns node. If failure raise exception if not print
    dns node and port
    :param dim_dns_node: name of the dim dns node to connect to
    :return: current dim dns node and dim dns port
    """
    try:
        dim.dic_set_dns_node(dim_dns_node)
        my_node = dim.dic_get_dns_node()
        dns_port = dim.dic_get_dns_port()
        MM6_total_counts = dim.dic_info_service('FPGA-MCS_actTotalCounts', 'I:1', fetch_total_counts, dim.MONITORED)
        MM6_Status = dim.dic_info_service('MM6_Status', 'C', fetch_status, dim.MONITORED)
    except Exception as err:
        print err
    else:
        print 'DNS_NODE = {} with DNS_PORT - {}'.format(my_node, dns_port)
        return my_node, dns_port, MM6_total_counts, MM6_Status


def read_ini_file(filename_ini):
    """
    Load analysis ini file
    :param filename_ini:
    :return: numpy array containing the parameters to be used by the optimizer
    """
    return np.genfromtxt(filename_ini, dtype='object')

def fetch_total_counts(tag,service_data_total_counts):
    global service_total_counts
    service_total_counts = service_data_total_counts

def definition():
    global service_lvstatus
    global service_total_counts
    service_lvstatus = "none"
    service_total_counts = -1

def fetch_status(tag,service_data_status):
    global service_lvstatus
    service_lvstatus = service_data_status.replace('\x00', '')
    print service_lvstatus


def total_counts(flag, init, start, stop, conf, reset, quit, sequencer, guiID, filename, config, userdata, mydata_conf, mydata_gui, mydata_strt):
    print flag
    print service_lvstatus
    while True:
        if (flag == 0) and ((service_lvstatus == "stopped") or (service_lvstatus == "aborted") or (service_lvstatus == "NOTConnected")):
            cs.call_process2(sequencer, conf, "C:{},C:{},C:{}".format(len(guiID), len(config), len(userdata)), mydata_conf)
            flag = 1
        if (flag == 1) and (service_lvstatus == "config"):
            flag = 2
            cs.call_process2(sequencer, init, "C:{}".format(len(guiID)), mydata_gui)
        if (flag == 2) and (service_lvstatus) == "init":
            flag =3
            cs.call_process2(sequencer, reset, "C:{}".format(len(guiID)), mydata_gui)
        if (flag == 3) and service_lvstatus == "data_reset":
            flag =4
            cs.call_process2(sequencer, start, "C:{},C:{}".format(len(guiID), len(filename)), mydata_strt)
        if (flag == 4) and service_lvstatus == "running":
            flag = 5
            cs.call_process2(sequencer, stop, "C:{}".format(len(guiID)), mydata_gui)
        if (flag == 5) and (service_lvstatus == "stopped"):
            print '\n'
            print 'run complete'
            #print "\n"+str(service_total_counts)
            return -1*service_total_counts


def objective_func(x, cs_objects, cs_data):
    """
    Define the objective function
    :param x: 1D array containing the voltages to be set
    :param args: tuple containing all extra parameters needed
    :return: average count rate for 100 shots
    """
    x = np.around(x,2)
    try:
        flag_range = 0
        for i in xrange(len(x)):
            if (x[i] <= float(cs_objects[i,4])) or (x[i] >= float(cs_objects[i,5])):
                flag_range = 1
                raise ValueError
        for i in xrange(len(x)):
                if flag_range == 0:
                    if int(cs_objects[i,2]) != -1:
                        cs.call_process2(cs_objects[i,0], cs_objects[i,1], "I:1,D:1", cs.pack_ch_val([int(cs_objects[i,2])], [x[i]]))
                    else:
                        cs.call_process2(cs_objects[i,0], cs_objects[i,1], "D:1", cs.pack_val([x[i]]))
                else:
                    return

        time.sleep(1)
        flag = 0
        value = total_counts(flag, *cs_data)
        # value = scop.rosen(x)
        return value
    except ValueError:
            print "Value error : value went out of bound"

def disconnect(dns_setup):
    dim.dic_release_service(dns_setup[2])
    dim.dic_release_service(dns_setup[3])

def cal_tof(m_value, nrevs):
    '''
    compute time of flight for a single revolution tg1
	and Acq_delay t_Acq_delay
    '''
    a1 =  3.547895982
    b1 = 0.643097964
    a = 2409.268546
    b = 0.517071131

    tg1 = ((a - a1)*np.sqrt(m_value) + (b - b1))/1000.
    tg1 = nrevs*tg1
    t_Acq_delay = round(a1*np.sqrt(m_value) + b1 + tg1 -10)*1000
    return (tg1, t_Acq_delay)
def config_meas(argv):
    ini = open("general.ini","r")
    key_list, val_list = np.loadtxt(ini, delimiter='=', usecols=(0, 1), dtype='object', unpack=True)
    devdict = {}
    for arg in range(len(key_list)):
        devdict[key_list[arg]] = val_list[arg]
    cfg = open(argv[0],"r")
    pars = np.loadtxt(cfg,skiprows=1, dtype = 'object')
    cfg.close()
    ini.close()
    return (devdict,pars)
def config_mag(add, f_mag, amp_mag, n_mag):
    #configure the function generator of the magnetron excitation
    cs.call_process2(add,'SetChannelFrequency', "I:1;D:1", cs.pack_ch_val([1], [f_mag], cs.DOUBLE))
    cs.call_process2(add,'SetChannelAmplitude', "I:1;D:1", cs.pack_ch_val([1], [amp_mag], cs.DOUBLE))
    cs.call_process2(add,'SetBurstCycles', "I:1;I:1", cs.pack_ch_val([1], [n_mag], cs.INT))
    delay2 = np.round(n_mag/f_mag/4e-9)*4e-9
    cs.call_process2(add,'SetTriggerDelay', "I:1;D:1", cs.pack_ch_val([2], [delay2], cs.DOUBLE))
def config_cycl(add1,add2,addej,f_plus,amp_plus,n_plus,f_c,amp_c,n_c,n_acc):
    #configure the function generator of pattern 1
    cs.call_process2(add1,'SetChannelFrequency', "I:1;D:1", cs.pack_ch_val([1], [f_plus], cs.DOUBLE))
    cs.call_process2(add1,'SetChannelAmplitude', "I:1;D:1", cs.pack_ch_val([1], [amp_plus], cs.DOUBLE))
    cs.call_process2(add1,'SetBurstCycles', "I:1;I:1", cs.pack_ch_val([1], [n_plus], cs.INT))
    cs.call_process2(add1,'SetTriggerDelay', "I:1;D:1", cs.pack_ch_val([1], [0], cs.DOUBLE))
    cs.call_process2(add1,'SetChannelFrequency', "I:1;D:1", cs.pack_ch_val([2], [f_c], cs.DOUBLE))
    cs.call_process2(add1,'SetChannelAmplitude', "I:1;D:1", cs.pack_ch_val([2], [amp_c], cs.DOUBLE))
    cs.call_process2(add1,'SetBurstCycles', "I:1;I:1", cs.pack_ch_val([2], [n_c], cs.DOUBLE))
    delay = np.round(n_plus/f_plus/4e-9)*4e-9
    print 'pi-pulse delay pattern 1 = ',delay
    cs.call_process2(add1,'SetTriggerDelay', "I:1;D:1", cs.pack_ch_val([2], [delay], cs.DOUBLE))
    #configure the function generator of pattern 2
    cs.call_process2(add2,'SetChannelFrequency', "I:1;D:1", cs.pack_ch_val([1], [f_plus], cs.DOUBLE))
    cs.call_process2(add2,'SetChannelAmplitude', "I:1;D:1", cs.pack_ch_val([1], [amp_plus], cs.DOUBLE))
    cs.call_process2(add2,'SetBurstCycles', "I:1;I:1", cs.pack_ch_val([1], [n_plus], cs.INT))
    cs.call_process2(add2,'SetTriggerDelay', "I:1;D:1", cs.pack_ch_val([1], [0], cs.DOUBLE))
    cs.call_process2(add2,'SetChannelFrequency', "I:1;D:1", cs.pack_ch_val([2], [f_c], cs.DOUBLE))
    cs.call_process2(add2,'SetChannelAmplitude', "I:1;D:1", cs.pack_ch_val([2], [amp_c], cs.DOUBLE))
    cs.call_process2(add2,'SetBurstCycles', "I:1;I:1", cs.pack_ch_val([2], [n_c], cs.INT))
    phase_acc = np.round(n_acc/f_c/4e-9)*4e-9
    delay = delay + phase_acc
    print 'pi-pulse delay pattern 2 = ', delay
    cs.call_process2(add2,'SetTriggerDelay', "I:1;D:1", cs.pack_ch_val([2], [delay], cs.DOUBLE))
    #configure the function generator for ejection and MCA trigger
    delay = delay + np.round(n_c/f_c/4e-9)*4e-9
    print 'ejection delay = ', delay
    print 'phase accumulation time = ', phase_acc
    cs.call_process2(addej,'SetTriggerDelay', "I:1;D:1", cs.pack_ch_val([1], [delay], cs.DOUBLE))

def main(argv):
    """
    Main function of the tuner
    :param argv: ini file of the optimization
    :return:
    """
    try:
        ##################################################################
        #Load settings
        ##################################################################
        general, parameters = config_meas(argv)
        add_mag = general['gen_mag_excitation']
        add_p1 = general['gen_pattern1']
        add_p2 = general['gen_pattern2']
        add_ej = general['gen_ejection']
        print add_mag
        print add_p1
        print add_p2
        print add_ej
        elements = parameters[:,0]
        f_mag =  np.around(parameters[:,1].astype(np.float),3)
        amp_mag = parameters[:,2].astype(np.float)
        n_mag = parameters[:,3].astype(np.int)
        f_plus = np.around(parameters[:,4].astype(np.float),3)
        amp_plus = parameters[:,5].astype(np.float)
        n_plus = parameters[:,6].astype(np.int)
        f_c = np.around(parameters[:,7].astype(np.float),3)
        amp_c = parameters[:,8].astype(np.float)
        n_c = parameters[:,9].astype(np.int)
        n_acc = parameters[:,10].astype(np.int)
        phase_acc_time = np.around(n_acc/f_c/4e-9)*4e-9
        nloops = int(argv[1])
        print elements
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
        rec_dat = np.vstack((elements,f_mag,amp_mag,n_mag,f_plus,amp_plus,n_plus,f_c,amp_c,n_c,n_acc,phase_acc_time)).transpose()
        np.savetxt(general['file_path']+'\params.rec',rec_dat, fmt = '%s', header = 'Element nu- nu-_pulse_ampl nu-_pulse_cycl nu+ nu+_pulse_ampl nu+_pulse_cycles nuc pi_pulse_ampl pi_pulse_cycl phase_acc_cycl phase_acc_time' )
#        elements = ["85Rb","87Rb"]
#        revno = [rev*50 for rev in range(1,21)]
#        masses = [38.9631579064506, 40.96127668, 84.91124116, 86.90863195, 132.9049034]
#        nloops = 5
#        add_mag = 'UT_Mag'
#        f_mag = 1084.8
#        amp_mag = 0.4
#        n_mag = 3
#        add_p1 = 'UT_P1'
#        add_p2 = 'UT_P2'
#        add_ej = 'UT_PI-ICR_EJ'
#        f_plus = [1071384.0433, 1046735.7990]
#        n_plus = [233,233]
#        amp_plus = [2.,2.]
#        f_c = [1072468.9946,1047820.8135]
#        amp_c =[2.1,2.1]
#        n_c = [1500,1500]
#        n_acc = [100000,100000]
        #print(n_plus[1]/f_plus[1]+n_acc[1]/f_c[1]+n_c[1]/f_c[1])

        ##################################################################
        #Aliases to commands and parameters of the CS
        ##################################################################
        init = "init"
        start = "strt"
        stop = "stop"
        conf = "conf"
        reset = "resd"
        quit = "quit"
        abort = "abpc"
        sequencer = "MM6"
        dim_dns_node = "pcisoltrap04.cern.ch"
        MMpath = "C:\\ISOLTRAP CS\\Settings\\MM6\\"
        guiID = "pcisoltrap21"
        ###################################################################
        #Setup dim connection
        ###################################################################
        dns_setup = connect(dim_dns_node)
        definition()
        ###################################################################
        #Connect to the TDC Labview VI
        ###################################################################
        LabVIEW = win32com.client.Dispatch("Labview.Application")
        VI = LabVIEW.getvireference('E:\PI-ICR\Read TDC 2017\TDC_DAQ\TDC_DAQ.vi')
        ###################################################################
        #Loop over number of cross-checks
        ###################################################################
        j = 0
        print "Measuring center positions"
        for k in range(len(elements)):
            j+=1
            filename = "D:\\ISOLTRAP CS\\Data\\auto_{}.dat".format(time.strftime("%b%d%H%M%S%Y", time.localtime(time.time())))
            LVconfig = MMpath+"config_"+elements[k]+".xml"
            LVUserData = MMpath+"userdata_"+elements[k]+".xml"
        ####################################################################
            mydata = prepare_data(filename, guiID, LVconfig, LVUserData)
            cs_data = [ init, start, stop, conf, reset, quit, sequencer, guiID, filename]
            for l in xrange(len(mydata)):
                cs_data.append(mydata[l])
    ################################################################
            config_mag(add_mag, f_mag[k], amp_mag[k], n_mag[k])
            config_cycl(add_p1,add_p2,add_ej,f_plus[k],0.001,n_plus[k],f_c[k],0.001,n_c[k],n_acc[k])
            sleep(10)
            tdcfile = general['file_path']+'\c_'+elements[k]+'_'+str(j)+'.bin'
            sleep(1)
            VI.setcontrolvalue('data path',str(tdcfile))
            sleep(1)
            VI.setcontrolvalue('Ctrl. start',str('TRUE'))
            sleep(1)
            VI.setcontrolvalue('Ctrl. start',str('FALSE'))
            sleep(1)
            flag = 0
            value = total_counts(flag, *cs_data)
            VI.setcontrolvalue('Ctrl. stop',str('TRUE'))
            sleep(1)
            VI.setcontrolvalue('Ctrl. stop',str('FALSE'))
            sleep(5)
        print "Measuring patterns"
        for i in range(nloops):
        ###################################################################
        #Loop over the nuclides
        ###################################################################
            for k in range(len(elements)):
                j+=1
                filename = "D:\\ISOLTRAP CS\\Data\\auto_{}.dat".format(time.strftime("%b%d%H%M%S%Y", time.localtime(time.time())))
                LVconfig = MMpath+"config_"+elements[k]+".xml"
                LVUserData = MMpath+"userdata_"+elements[k]+".xml"
        ####################################################################
                mydata = prepare_data(filename, guiID, LVconfig, LVUserData)
                cs_data = [ init, start, stop, conf, reset, quit, sequencer, guiID, filename]
                for l in xrange(len(mydata)):
                    cs_data.append(mydata[l])
    ################################################################
                config_mag(add_mag, f_mag[k], amp_mag[k], n_mag[k])
                config_cycl(add_p1,add_p2,add_ej,f_plus[k],amp_plus[k],n_plus[k],f_c[k],amp_c[k],n_c[k],n_acc[k])
                sleep(10)
                tdcfile = general['file_path']+'\p1p2_'+elements[k]+'_'+str(j)+'.bin'
                sleep(1)
                VI.setcontrolvalue('data path',str(tdcfile))
                sleep(1)
                VI.setcontrolvalue('Ctrl. start',str('TRUE'))
                sleep(1)
                VI.setcontrolvalue('Ctrl. start',str('FALSE'))
                sleep(1)
                flag = 0
                value = total_counts(flag, *cs_data)
                VI.setcontrolvalue('Ctrl. stop',str('TRUE'))
                sleep(1)
                VI.setcontrolvalue('Ctrl. stop',str('FALSE'))
                sleep(5)
    except KeyboardInterrupt:
        print('You pressed Ctrl+C!')
    finally:
        if service_lvstatus != "stopped":
            cs.call_process2(sequencer, abort, "C:{}".format(len(guiID)), mydata[3])
        disconnect(dns_setup)
        cs.call_process2(sequencer, quit, "C:{}".format(len(guiID)), mydata[3])

if __name__=="__main__":
    main(sys.argv[1:])
