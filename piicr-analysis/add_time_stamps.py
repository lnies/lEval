import subprocess
import time
import os


def add_time_stamps(folder_name):
    '''Adds timestamps for raw piicr data from .tdc to .txt'''
    exe_str = 'G:\\Experiments\\ISOLTRAP\\Software\\PI-ICR\\Python-DAQ\\remote_piicr_add_timestamps.exe'
    folder_name = 'E:\\Jonas\\test'     # folder where raw data .txt and .tdc is placed
    process = subprocess.Popen(exe_str+' '+folder_name,  stderr=subprocess.PIPE)
    time.sleep(2)   # just a waiting time until all files are manipulated
    process.kill()

if __name__ == '__main__':
    add_time_stamps('bla')
