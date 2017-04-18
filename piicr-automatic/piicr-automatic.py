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


def total_counts(flag, init, start, stop, conf, reset, quit, sequencer, guiID, filename, config, userdata, mydata_conf, mydata_gui, mydata_strt, tg1, tacq_delay, data_setletter):
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
            cs.call_process2('PPG_ABC', 'SetLetterTime',"C:{},D:{},I:{}".format(len(tg1),1,1), data_setletter)
            cs.call_process2('FAST-MCS', 'SetAqnDelay',"D:1",tacq_delay)
            sleep(1)
            cs.call_process2('FAST-MCS', 'StartScan',"C:1",cs.pack_array(["0"]))
            sleep(5)
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

def main(argv):
    """
    Main function of the tuner
    :param argv: ini file of the optimization
    :return:
    """
    try:
        ##################################################################
        #Settings to adapt
        ##################################################################
        elements = ["39K","41K","85Rb","87Rb","133Cs"]
        revno = [rev*50 for rev in range(1,21)]
        masses = [38.9631579064506, 40.96127668, 84.91124116, 86.90863195, 132.9049034]
        nloops = 1
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
        guiID = "pcisoltrap32"
        ###################################################################
        #Setup dim connection
        ###################################################################
        dns_setup = connect(dim_dns_node)
        definition()
        ###################################################################
        #Loop over number of cross-checks
        ###################################################################


        for i in range(len(revno)):

            for j in range(nloops):

        ###################################################################
        #Loop over the nuclides
        ###################################################################
                for k in range(len(elements)):
                    filename = "D:\\ISOLTRAP CS\\Data\\autotuner_{}.dat".format(time.strftime("%b%d%H%M%S%Y", time.localtime(time.time())))
                    LVconfig = MMpath+"config_"+elements[k]+".xml"
                    LVUserData = MMpath+"userdata_"+elements[k]+".xml"
        ####################################################################
                    mydata = prepare_data(filename, guiID, LVconfig, LVUserData)
                    cs_data = [ init, start, stop, conf, reset, quit, sequencer, guiID, filename]
                    for l in xrange(len(mydata)):
                        cs_data.append(mydata[l])
    ################################################################

    #################################################################
    #Load objects to be scanned
    #################################################################
    #    inifile = argv[0]
    #    cs_objects = read_ini_file(inifile)
    ##################################################################
                    tg1, t_Acq_delay = cal_tof(masses[k],revno[i])
                    print elements[k], revno[i], tg1, t_Acq_delay
                    pack_tg1 = cs.pack_val([tg1],cs.DOUBLE)
                    pack_time_pattern = cs.pack_val([0],cs.INT)
                    letter = 'G1'
                    pack_letter = cs.pack_val(letter,cs.CHAR)
                    cs_data.append(pack_letter)
                    pack_tacqn_delay = cs.pack_val([t_Acq_delay],cs.DOUBLE)
                    cs_data.append(pack_tacqn_delay)
                    set_letter_time_data = pack_letter + pack_tg1 + pack_time_pattern
                    cs_data.append(set_letter_time_data)
                    flag = 0
                    value = total_counts(flag, *cs_data)
                    cs.call_process2('FAST-MCS', 'StopScan',"C:1",cs.pack_array(["0"]))
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
