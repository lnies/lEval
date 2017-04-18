#!/usr/bin/env python
#             author: Rosen Matev     (LHCb)
#  prototype library: Dietrich Bech   (GSI)
# package maintainer: Dinko Atanasov  (ISOLTRAP)
# package maintainer: Maxime Mougeout (ISOLTRAP)
#            version: 1.0
#            contact: dinko.atanasov@cern.ch, maxime.mougeot@cern.ch
#               date: 24 September 2016
#  last modification: 29 November  2016
#
# This library contains basic functions to pack data in accordance with the CS Framework (LabVIEW wrapper).
# For more information please look at the CS Documentation (https://wiki.gsi.de/foswiki/bin/view/CSframework/WebHome)
# and the communication layer of DIM and its wrappers for python/LabVIEW (http://dim.web.cern.ch/dim/).
from struct import *

import pydim
# import time
###################################################################################
################    FORMAT characters have this meaning        ####################
###################################################################################
SHORT  = 'h' # C = short        |Python = int          | Standard size = 2 bytes  #
INT    = 'i' # C = int          |Python = int          | Standard size = 4 bytes  #
CHAR   = 'c' # C = char         |Python = str of len 1 | Standard size = 1 bytes  #
UINT   = 'I' # C = unsigned int |Python = int          | Standard size = 4 bytes  #
LONG   = 'l' # C = long         |Python = int          | Standard size = 4 bytes  #
ULONG  = 'L' # C = unsigned long|Python = int          | Standard size = 4 bytes  #
LLONG  = 'q' # C = long long    |Python = int          | Standard size = 8 bytes  #
ULLONG = 'Q' # C = unsg longlong|Python = int          | Standard size = 8 bytes  #
FLOAT  = 'f' # C = float        |Python = float        | Standard size = 4 bytes  #
DOUBLE = 'd' # C = double       |Python = double       | Standard size = 8 bytes  #
###################################################################################
###################################################################################
def pack_size(x):
    """
    Converts unsigned int (x) into a byte array object of type long in big-endian.
    :param x: integer number
    :return: bytearray(x)
    """
    return bytearray(pack('!'+LONG, x))
##################################################################################


def pack_header(x):
    """
    Converts an object into byte array of type INT in big-endian. This is used to convert
    headers. See pack_array() if one wants to convert size of object.
    :param x: object of type array
    :return: bytearray(x)
    """
    return bytearray(pack('!'+len(x)*INT, *x))
##################################################################################


def pack_array(x, FORMAT=CHAR):
    """
    Converts array of type 'format' to byte array object.
    The size of the object is appended to the byte array as well.
    :param x: object of type array
    :param FORMAT: Default=CHAR, Definition of FORMAT is given at the start of this script.
    :return: bytearray(x)
    """
    return bytearray(pack('!'+len(x)*FORMAT, *x))
##################################################################################


def pack_ch_val(ch_idx, ch_val, FORMAT=DOUBLE):
    """
    Prepare channel value in suitable format to be send to CSframework Object.
    :param ch_idx: Channel Index, FORMAT = INT
    :param ch_val: Channel Value,
    :param FORMAT: Format used to convert Channel value
    :return: bytearray(ch_idx+ch_val)
    """
    byte_ch = pack_array(ch_idx, INT)
    byte_val = pack_array(ch_val, FORMAT)
    return pack_size(len(byte_ch)) + byte_ch + pack_size(len(byte_val)) + byte_val
##################################################################################


def pack_val(val, FORMAT=DOUBLE):
    """
    Prepare value in suitable format to be send to CSframework Object.
    :param val: Value,
    :param FORMAT: Format used to convert Channel value
    :return: bytearray(val)
    """
    byte_val = pack_array(val, FORMAT)
    return pack_size(len(byte_val))+byte_val
##################################################################################


def call_process2(objectName, commandSelector, typeDescriptor, data, accessID=-1):
    """
    This function prepares data to be send to the CS Framework. It is based on
    the original library written for C by D. Beck. The FORMAT of the data expected
    by the CS LabVIEW wrapper is a byte array configured in big-endian.
    Prototype function: CSCallProcess2("char", "char", "char", accessID, byteArray, size);

    Example taken from the the C source code to set value of TestSystem\Super object with method SetEvtPresetValue:
    CSCallProcess2("TestSystem\\Super", "SetEvtPresetValue", "", -1, byteArray, size);

    :param objectName: string describing the CS Object name
    :param commandSelector: string describing the CS Selector method
    :param typeDescriptor: string describing the FORMAT of the data (examples: "D:1"; "I:1"; "C:23434")
    :param data: byte array with the values to be send
    :param accessID: DIM accessID. See documentation of DIM for more information. Default = -1
    :return: None
    """
    header = pack_header([
     0,                     # total size of my data, added later
     1,                     # version number
     0,                     # call type: simple call = 0
     0,                     # error number: no error
     0,                     # expire time: never expire
     len(objectName),       # length of object name
     len(commandSelector),  # length of selector
     0,                     # length of async object: unused (simple call)
     0,                     # length of async selector: unused (simple call)
     len(typeDescriptor),   # length of descriptor
     len(data),             # size of net data
    ])
    # Append data to the header
    new_data = header+pack_array(objectName) + pack_array(commandSelector) + pack_array(typeDescriptor) + data
    # add the total size of my data to the header
    new_data[:4] = pack_size(len(new_data))
    # call the send function
    send_command(objectName, commandSelector, accessID, new_data)
##################################################################################


def send_command(objectName, commandSelector, accessID, data):
    """
    This function packs the rest of the information expected by CS LabVIEW wrapper for checksum.
    Actual send is performed by using the pydim and the dic_cmnd_service(). Expected type of the send data is tuple.
    For usage and more information see documentation of DIM and/or PyDIM.

    :param objectName: string describing the CS Object name
    :param commandSelector: string describing the CS Selector method
    :param accessID: string describing the FORMAT of the data (examples: "D:1"; "I:1"; "C:23434")
    :param data: byte array with the values to be send
    :return: None
    """
    header = pack_header([
     0,           # total size of my data, added later
     1,           # version number of data format
     accessID,    # accessID
     len(data),   # size of net data
    ])
    # Append the net data with the rest of the information
    new_data = header + data
    # add the total size of the net data to the header
    #print "before adding", len(map(int, new_data))
    new_data[:4] = pack_size(len(new_data))
    # build the object command
    command = 'CAE_{}_{}'.format(objectName, commandSelector)
    #print "after adding", len(new_data)
    # finally send the command using PyDIM
    pydim.dic_cmnd_service(command, (len(new_data)+1, 1, accessID, len(data), str(data),), "I:1;I:1;I:1;I:1;C" )
##################################################################################


if __name__ == '__main__':
    pydim.dic_set_dns_node("pcisoltrap28.cern.ch")
    my_node = pydim.dic_get_dns_node()
    dns_port = pydim.dic_get_dns_port()
    print 'DNS_NODE = {} with DNS_PORT - {}'.format(my_node, dns_port)

    # Example1: SetChannelVoltage at SimPowerSupply
    myObject = 'SimPowerSupply'
    myMethod = 'SetChannelVoltage'
    presetChannel = 1
    presetValue = 12
    call_process2(myObject, myMethod, "I:1;D:1", pack_ch_val([presetChannel], [presetValue], DOUBLE))
    # Example1: SetFreqeuncy at SimAFG
    myObject2 = 'SimAFG'
    myMethod2 = 'SetFrequency'
    presetValue2 = 999.0
    call_process2(myObject2, myMethod2, "D:1", pack_val([presetValue2], DOUBLE))
