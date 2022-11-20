# -*- coding: utf-8 -*-
"""
Created on Mon 11 January 2022
Modified by Lukas.Nies@cern.ch on 04/02/2021
@author: Lukas Nies
@contact: Lukas.Nies@cern.ch
@license: MIT
"""

import mmap
import os
import sys
import numpy as np
import pandas as pd
from io import StringIO
from configparser import ConfigParser, RawConfigParser
from chardet import detect
import re
import datetime

class CustomParser(ConfigParser):
	"""
	Original from https://gitlab.cern.ch/datanaso/dEval/-/blob/master/etc/utils.py
	ConfigParser with a custom dictionary conversion method.
	"""
	def as_dict(self):
		d = dict(self._sections)
		for k in d:
			d[k] = dict(self._defaults, **d[k])
			d[k].pop('__name__', None)
		return d

### MR-ToF MS RELATED FILE HANDLING

class ProcessorBase():
	'''
	Base-class for MR-ToF MS data proessing
	'''
	def __init__(self):
		self.files = []
		self.data = {}
		self.pars = {}
		self.df_dict = {}
		self.step = 0

	def add_all(self, to_csv=False):
		'''
		Returns sum of all handled files while incrementing sweep numbers from one file to the next to 
		get a rolling sweep number. 
		Parameters:
			- to_csv: File name to store combined .csv in
		Return:
			- flattened and sweep-adjusted dataframe
		'''
		if len(self.files) == 0:
			print(f"(ProcessorBase.add_all): Data not processed yet or empty.")
			return
		# Adjust the sweep numbers
		for i in np.arange(0, len(self.df_dict)):
			if i == 0: 
				continue
			key = list(self.df_dict.keys())[i]
			# print(key)
			key_m1 = list(self.df_dict.keys())[i-1]
			# print(key_m1, self.df_dict[key_m1].iloc[-1]['sweep'])
			self.df_dict[key]['sweep'] += self.df_dict[key_m1].iloc[-1]['sweep'] + 1
		#
		df = pd.concat(self.df_dict, ignore_index=True)
		# Save to file if file name is passed
		if to_csv != False:
			df.to_csv(to_csv, index=False)
		#
		return df

	def parse_header(self, key, txt):
		parser = CustomParser(strict=False)
		parser.read_file(StringIO(txt))
		tmp = parser.as_dict()
		self.pars[key] = dict(**tmp['CHN1'])
		self.pars[key].update(**tmp['MPA4A'])

class MCS6Lst(ProcessorBase):
	'''
	Process each list file contained in array of lst_files
	Param:
		lst_files: array of list files
	'''
	"""
	Created on Mon 23 April 2018
	Modified and adapted to Python 3 on Wed 17 July 2019
	Modified by Lukas.Nies@cern.ch on 21/10/2021
	@author: Maxime Mougeot
	@author: Jonas Karthein
	@contact: maxime.mougeot@cern.ch
	@contact: jonas.karthein@cern.ch
	@license: MIT license
	"""
	def __init__(self):
		"""
		Initialize the conversion dataframe and some other varaibles
		"""
		ProcessorBase.__init__(self) # self.map_files = file_array
		self.files = []
		#-------------------------Create a dataframe containing the conversion table from the MCS6A manual--------------------------------------#
		Time_Patch_Value = ['0', '5', '1', '1a', '2a', '22', '32','2','5b','Db','f3','43','c3','3']
		conversion_dict = {'Data_Length' : pd.Series([2,4,4,6,6,6,6,6,8,8,8,8,8,8], index=Time_Patch_Value),
		'Data_Lost_Bit' : pd.Series([np.nan ,np.nan ,np.nan ,np.nan ,np.nan ,np.nan ,47,np.nan ,63,np.nan ,47,63,np.nan ,63], index=Time_Patch_Value),
		'Tag_Bits' : pd.Series([(np.nan,np.nan) ,(np.nan,np.nan) ,(np.nan,np.nan) ,(np.nan,np.nan) ,(40,47),(40,47),(np.nan,np.nan), (np.nan,np.nan), (48,62),(48,63),(48,63),(48,62),(48,63),(58,62)], index=Time_Patch_Value),
		'Sweep_Counter': pd.Series([(np.nan,np.nan),(24,31),(np.nan,np.nan),(32,47),(32,39),(np.nan,np.nan),(40,46),(np.nan,np.nan),(32,47),(32,47),(40,46),(np.nan,np.nan),(np.nan,np.nan),(np.nan,np.nan)], index=Time_Patch_Value),
		'Time_Bits': pd.Series([12,20,28,28,28,36,36,44,28,28,36,44,44,54], index=Time_Patch_Value),
		'Max_Sweep_Length': pd.Series([0.0000004096,0.000105,0.027,0.027,0.027,6.872,6.872,1759.2,0.027,0.027,6.872,1759.2,1759.2,1801440], index=Time_Patch_Value)}
		self.conversion_df = pd.DataFrame(conversion_dict)
		self.df_dict = {}

	def convert_bytes(self,bytearray,nb_bits,data_lost_bit,tag_bits,sweep_counter,time_bits,verbose=0):
		'''
		Perform the actual conversion of a single event from binary to integer numbers
		See pages  5-19 and 5-20 of FastCom MCS6A for more information on the bits meaning
		Param:
			bytearray : an array of nb_bits/8 bytes encapsulating the data of a single MCS6A stop event
			nb_bits : total number of bits on which the data is encoded
			data_lost_bit : Data lost bit. Indicates if the fifo was fullself. 1 bit index
			tag_bits : Tag bits. Tag info of a single stop event (see manual). Tuple containing the bit indexes.
			sweep_counter: Sweep number of a single stop event. Tuple containing the bit indexes.
			time_bits : Number of bits encoding the time of flight a single stop eventself.
			The tof seems to be given in an unbinned format with 100ps resolution (to be confirmed).
		Return:
			Decoded tof, sweep, channel, edge, tag, fifo
		'''
		bit_word = ''

		for bytes in reversed(bytearray):
			bit_word += '{0:08b}'.format(ord(chr(bytes)))

		#convert data lost bit always first index in the reversed array (last in the manual)
		if np.isnan(data_lost_bit):
			fifo = np.nan
			index_data_lost_bit = -1
		else:
			index_data_lost_bit = nb_bits-1-int(data_lost_bit)
			fifo = int(bit_word[index_data_lost_bit],2)

		#convert tag bit
		if np.isnan(tag_bits[0]):
			tag = np.nan
			index_high_tag = -1
			index_low_tag = -1
		else:
			index_high_tag = nb_bits-1-int(tag_bits[1])
			index_low_tag = nb_bits-1-int(tag_bits[0])
			tag = int(bit_word[index_high_tag:index_low_tag+1],2)

		#convert sweep number
		if np.isnan(sweep_counter[0]):
			sweep = np.nan
			index_high_sweep = -1
			index_low_sweep = -1
		else:
			index_high_sweep = nb_bits-1-int(sweep_counter[1])
			index_low_sweep = nb_bits-1-int(sweep_counter[0])
			sweep = int(bit_word[index_high_sweep:index_low_sweep+1],2)

		if bit_word != "000000000000000000000000000000000000000000000000" and verbose>1:
			print(f"bit_word: {bit_word}")
			print(f"index_data_lost_bit: {fifo}")
			print(f"index_high_tag: {index_high_tag}")
			print(f"index_low_tag: {index_low_tag}")
			print(f"tag: {tag}")
			print(f"index_high_sweep: {index_high_sweep}")
			print(f"index_low_sweep: {index_low_sweep}")
			print(f"sweep: {sweep}")

		#convert time of flight
		index_high_tof = max(index_low_sweep,index_low_tag,index_data_lost_bit)+1
		index_low_tof = index_high_tof+time_bits
		tof = int(bit_word[index_high_tof:index_low_tof],2)

		#these are always there no matter the format
		channel = int(bit_word[index_low_tof+1:],2)
		edge = int(bit_word[index_low_tof],2)

		# if tof != 0:
		#     print(tof, sweep-1, channel, edge, tag, fifo)

		return tof, sweep-1, channel, edge, tag, fifo

	def decode_binary(self,binary,time_patch, verbose = 0):
		'''
		Read the binary part of the file by chunks and decode each chunk according to the format
		given in time_patch
		The length of a chunk is given by the time_patch
		Param : binary part of
		Return: nunpy array containing the converted data : tof, sweep, channel, edge, tag, fifo
		'''

		#-----------extract data from the dataframe--------------------------#
		data_length = int(self.conversion_df.loc[time_patch.decode('ascii'),'Data_Length'])
		nb_bits = 8*data_length #convert nb of bytes into nb of bits
		data_lost_bit = self.conversion_df.loc[time_patch.decode('ascii'),'Data_Lost_Bit']
		tag_bits = self.conversion_df.loc[time_patch.decode('ascii'),'Tag_Bits']
		sweep_counter = self.conversion_df.loc[time_patch.decode('ascii'),'Sweep_Counter']
		time_bits = int(self.conversion_df.loc[time_patch.decode('ascii'),'Time_Bits'])
		max_sweep_length = self.conversion_df.loc[time_patch.decode('ascii'),'Max_Sweep_Length']

		steps = len(binary[binary.tell():])/data_length
		first_it = True

		if verbose>1:
			print(f"Data length: {data_length}\nN bits: {nb_bits}\nData lost bit: {data_lost_bit}\n\
					tag bits: {tag_bits}\nsweep_counter: {sweep_counter}\ntime_bits: {time_bits}\n\
					max sweep length: {max_sweep_length}\nsteps: {steps}\n")
		# !
		# Introduce sweep_counter_overflow: in some cases, MCS6 seems to allocated only a small amount of bits for storing the sweep number.
		# In time_patch=32 this is for example only 7 bits -> can count to 128 and then resets to 0. With the overflow counter we just
		# count how many overflows happen and just at the necessary amount of sweeps to the overall sweep number
		sweep_counter_overflow = 0
		old_sweep = 0 # for detecting when overflow happens
		# loop through all bytewords
		for i in range(int(steps)):
			if verbose>0:
				if (i%(int(steps/10))==0):
					print(f"Step {i} of {steps}.")
			byteword = binary.read(data_length)
			tof, sweep, channel, edge, tag, fifo = self.convert_bytes(byteword,nb_bits,
				data_lost_bit, tag_bits, sweep_counter, time_bits, verbose=verbose)
			# Check whether overflow happened (for example old_sweep = 127, new sweep is 0)
			# Only do for non-zero events:
			if tof != 0: 
				if verbose>1: print(f"old_sweep: {old_sweep}")
				if old_sweep > sweep:
					sweep_counter_overflow += 1
				if verbose>1: print(f"sweep_counter_overflow: {sweep_counter_overflow}")
				old_sweep = sweep 
				# Add overflow to the sweep number (in case sweep has 7bit int -> 2**7=128)
				sweep += sweep_counter_overflow*(2**(sweep_counter[1]-sweep_counter[0]+1))
				if verbose>1: print(f"sweep: {sweep}")
			#
			if channel != 0 :#means for real data
				if first_it:
					converted_data = np.array([tof, sweep, channel, edge, tag, fifo])
					first_it = False
				else :
					converted_data = np.vstack((converted_data, np.array([tof, sweep, channel, edge, tag, fifo])))
		binary.close()
		return converted_data

	def get_time_patch_and_binary(self, listfile, verbose=False):
		'''
		Memory map the list file and isolate the time_patch and the binary part of the data
		Param:
			listfile : input list file
		Return:
			mapped_file : memore map of the input listfile
			time_patch : string code indicating the format in which the data are written (see manual)
		'''

		#
		mapped_file = mmap.mmap(listfile.fileno(), 0, access=mmap.ACCESS_READ)
		search_dict = {'section' : '[DATA]' , 'list_file_type' : 'time_patch', 'file_date' : 'cmline0', 'report-file': 'REPORT-FILE'}
		#-----------------set file index to time patch code -----------------#
		pos_type_from = mapped_file.find(search_dict['list_file_type'].encode('ascii'))+len(search_dict['list_file_type'])+1
		mapped_file.seek(pos_type_from)
		time_patch = mapped_file.readline().strip('\r\n'.encode('ascii'))
		self.pars['time_patch'] = time_patch.decode()
		#-----------------set file index to time patch code -----------------#
		mapped_file = mmap.mmap(listfile.fileno(), 0, access=mmap.ACCESS_READ)
		pos_type_from = mapped_file.find(search_dict['file_date'].encode('ascii'))+len(search_dict['file_date'])+1
		mapped_file.seek(pos_type_from)
		file_date = mapped_file.readline().strip('\r\n'.encode('ascii'))
		self.pars['cmline0'] = file_date.decode()
		#-----------------set file index to time patch code -----------------#
		mapped_file = mmap.mmap(listfile.fileno(), 0, access=mmap.ACCESS_READ)
		pos_type_from = mapped_file.find(search_dict['report-file'].encode('ascii'))+len(search_dict['report-file'])+1
		mapped_file.seek(pos_type_from)
		report_file = mapped_file.readline().strip('\r\n'.encode('ascii'))
		self.pars['report-file'] = report_file.decode()
		#-----------set file index to beginning of DATA-----------------------------------#
		pos_data_from = mapped_file.find(search_dict['section'].encode('ascii'))
		mapped_file.seek(pos_data_from)
		#---readline and there no matter what the file index should point to the beginning of the binary data
		mapped_file.readline()

		if verbose>1:
			print(f"pos_type_from: {pos_type_from}\npos_data_from: {pos_data_from}\ntime_patch: {time_patch}")

		return mapped_file, time_patch

	def process(self,file_array,to_csv = False, verbose=0):
		"""
		Perform the processing of the files 
		Parameters:
			- file_array: Array of file-paths
			- to_csv: if true, saves files under it's file name with .csv extension
			- verbose: verbosity
		"""
		full_info = False   # for regular application the channel, edge, tag and fifo info are constant so they don't have to be saved. In that case keep full_info = False
		self.files = file_array
		for filename in self.files:
			with open(filename,'rb') as listfile:

				binary, time_patch = self.get_time_patch_and_binary(listfile, verbose=verbose)

				if full_info:
					converted_data = self.decode_binary(binary,time_patch,verbose=verbose)    # np.array with tof, sweep, channel, edge, tag, fifo
					header_res ='tof,sweep,channel,edge,tag,fifo'
					if to_csv:
						np.savetxt('{}/{}.csv'.format(os.path.split(filename)[0],os.path.splitext(os.path.basename(filename))[0]),converted_data,
					fmt = '%i,%i,%i,%i,%f,%f', header = header_res)
				else:
					converted_data = pd.DataFrame(self.decode_binary(binary,time_patch,verbose)[:, [0,1]], columns=['tof', 'sweep'])     # saves only tof and sweep info
					converted_data.tof = converted_data.tof/10  # 100ps -> ns
					converted_data.sweep = converted_data.sweep.astype('int64')  # sweep is only int
					if to_csv:
						converted_data.to_csv('{}/{}.csv'.format(os.path.split(filename)[0],os.path.splitext(os.path.basename(filename))[0]), index=False)
				print('File {} loaded successfully!'.format(os.path.splitext(os.path.basename(filename))[0]))
				self.df_dict[os.path.splitext(os.path.basename(filename))[0]] = converted_data
		if full_info == False:
			return(pd.concat(self.df_dict, axis=1))    # convert dict of dataframes into one dataframe with two column name levels

class MPANTMpa(ProcessorBase):
	"""
	Original from https://gitlab.cern.ch/datanaso/dEval/-/blob/master/etc/mpant.py
	Class handling the MPANT (mpa) data files!
	Data format is 'asc'
	"""
	def __init__(self):
		ProcessorBase.__init__(self) # self.map_files = file_array
		self.files = []

	def read(self, mpa):
		"""
		Read mpa data file
		:return:
		"""
		with open(mpa, 'r') as f:
			fs = f.read()
			name = os.path.basename(mpa).split('.')[0]
			raw_header, raw_data = fs.split('[DATA]\n')
			if bool(raw_data) and len(raw_data.split(' ')) >= 9:
				self.parse_header(name, raw_header)
				self.df_dict[name] = pd.read_csv(StringIO(raw_data), delimiter=' ', usecols=(0, 1, 2), header=None, names=['tof', 'sweep', 'counts'])
				self.data[name] = self.df_dict[name].to_numpy()

	def process(self, files, to_csv = False):
		'''
		Processes MPA file.
		Parameters:
		- files: array of files to be processed
		- to_csv: Defaults to False; if true, saves the processed files
		Return:
		- pars: dictionary of measurement parameters
		- data: 2D array of processed data
		- data_dict: array of panda dataframes with data
		'''
		self.files = files
		# for i, f in enumerate(self.files):
		for f in self.files:
			self.read(f)
			# Convert bins to tof in ns
			name = os.path.basename(f).split('.')[0]
			# Check of file name is in dict 
			if name not in self.df_dict.keys():
				print(f"(MPANTMpa.process): WARNING: file {name} not processed. Empty?")
				continue
			#
			self.df_dict[name].tof = (self.df_dict[name].tof - 0.5) * float(self.pars[name]['calfact']) + float(self.pars[name]['caloff'])
			if to_csv:
				self.df_dict[name].to_csv('{}/{}.csv'.format(os.path.split(f)[0],os.path.splitext(os.path.basename(f))[0]), index=False)
		#
		return self.pars, self.data, pd.concat(self.df_dict, axis=1)

class MCDWIN887(ProcessorBase):
	"""
	Original from: https://gitlab.cern.ch/datanaso/dEval/-/blob/master/etc/mcdwin.py
	Class handling the MCDWIN (.887) data files!
	Pass the .887 file and the data format is 'csv, asc' is handled automatically.
	"""
	def __init__(self):
		ProcessorBase.__init__(self) # self.map_files = file_array

	def get_signals(self):
		return self.next_887, self.status_887

	def read(self, p887):
		"""
		Read mpa data file
		:return:
		"""
		if self.parse_header(p887) == -1:
			return -1
		folder = os.path.dirname(p887) + os.path.sep
		name = os.path.basename(p887).split('.')[0]

		if self.pars[name]['sweepmode'] == '80' and self.pars[name]['fmt'] == 'csv':
			self.df_dict[name] = self.read_csv_80(folder, name)
		elif self.pars[name]['sweepmode'] == '84' and self.pars[name]['fmt'] == 'csv':
			self.df_dict[name] = self.read_csv_84(folder, name)
		elif self.pars[name]['sweepmode'] == '84' and self.pars[name]['fmt'] == 'asc':
			self.df_dict[name] = self.read_asc_84(folder, name)

		if self.df_dict[name]['counts'].sum() != 0:
			self.data[name] = self.df_dict[name][self.df_dict[name]['counts'] > 0].to_numpy()
		else:
			self.data[name] = self.df_dict[name][0:1].to_numpy()

	def read_csv_80(self, folder, name):
		self.pars[name]['cycles'] = '1'
		fname = folder + name + f'.{self.pars[name]["fmt"]}'
		df = pd.read_csv(fname, sep='\t', header=None, dtype=float, prefix='Y')

		chan = [y for y in range(int(self.pars[name]['range']))]
		cycles = [np.full((int(self.pars[name]['range'])), i) for i in range(1)]
		slices = np.concatenate(cycles, axis=0)
		return pd.DataFrame({'tof': np.array(chan), 'sweep': slices, 'counts': df['Y1'].to_numpy()})

	def read_csv_84(self, folder, name):
		fname = folder + name + f'.{self.pars[name]["fmt"]}'
		df = pd.read_csv(fname, sep='\t', header=None, dtype=float, prefix='Y')

		chan = [y for x in range(int(self.pars[name]['cycles'])) for y in range(int(self.pars[name]['range']))]
		cycles = [np.full((int(self.pars[name]['range'])), i) for i in range(int(self.pars[name]['cycles']))]
		slices = np.concatenate(cycles, axis=0)

		return pd.DataFrame({'tof': np.array(chan), 'sweep': slices, 'counts': df['Y1'].to_numpy()})

	def read_asc_84(self, folder, name):
		fname = folder + name + f'.{self.pars[name]["fmt"]}'
		with open(fname, 'r') as f:
			fs = f.read()
			_, raw_data = fs.split('[DATA]\n')
			if len(raw_data.split(' ')) >= 3:  # If there is data (more than three lines), load it to df
				df = pd.read_csv(StringIO(raw_data), delimiter=' ', usecols=(0, 1, 2), header=None, names=['tof', 'sweep', 'counts'])
			else:  # Else build an empty dummy df
				df = pd.DataFrame({'tof': np.zeros(5),
								   'sweep': np.zeros(5),
								   'counts': np.zeros(5)})
		return df

	def parse_header(self, key):
		"""
		Parse header if file exists and is .887 file
		"""
		try:
			with open(key, 'r') as f:
				fs = '[root]\n' + f.read()
				key = os.path.basename(key).split('.')[0]
				parser = CustomParser(strict=False)
				parser.read_file(StringIO(fs))
				tmp = parser.as_dict()

				self.pars[key] = dict(**tmp['root'])
				# correct for inconsistent 'fmt' keyword in FASTCom configuration file (.887)
				if self.pars[key]['fmt'] == '3':
					self.pars[key]['fmt_idx'] = '3'
					self.pars[key]['fmt'] = 'asc'
		except(Exception, FileNotFoundError) as err:
			return -1

	def process(self, files, to_csv=True):
		"""
		
		"""
		#
		self.files = files
		#
		for i, f in enumerate(files):
			if self.read(f) == -1:
				print(f"(MCDWIN887.process): WARNING: file {f} not processed. Empty?")
				continue
			# Convert bins to tof in ns
			name = os.path.basename(f).split('.')[0]
			# Check of file name is in dict 
			if name not in self.df_dict.keys():
				print(f"(MCDWIN887.process): WARNING: file {name} not processed. Empty?")
				continue
			#
			self.df_dict[name].tof = self.df_dict[name].tof * float(self.pars[name]['calfact']) + float(self.pars[name]['caloff'])
			if to_csv:
				self.df_dict[name].to_csv('{}/{}.csv'.format(os.path.split(f)[0],os.path.splitext(os.path.basename(f))[0]), index=False)

		#
		return self.pars, self.data, pd.concat(self.df_dict, axis=1)

### ToF-ICR RELATED FILE HANDLING

class MM8():
	'''
	Processes .dat files passed as a list 
	'''
	"""
	Created on Mon 04 April 2022
	Adopted from D. Atanasov https://gitlab.cern.ch/datanaso/dEval
	@author: Lukas Nies
	@author: Dinko Atanasov
	@contact: Lukas.Nies@cern.ch
	@license: MIT license
	"""
	def __init__(self):
		#
		self.files = []
		self.header_dict = {}

	def __parse_header(self, mm_file):
		"""
		@ Author: D. Atanasov (https://gitlab.cern.ch/datanaso/dEval)
		# Modified by L. Nies
		Loads information from MM tmp files to the pars dictionary.
		"""
		try:
			with open(mm_file, 'r+b') as f:
				raw_data = f.read()
				enc = detect(raw_data)['encoding']
				if enc is None:
					header_buffer = raw_data.decode('cp437', errors='ignore')
				else:
					header_buffer = raw_data.decode(enc, errors='ignore')
				cut_from = 8
				cut_to = header_buffer.find('*----------------here the binary part starts --------------*\n')
				header = header_buffer[cut_from: cut_to]
				header = header.replace('Ejection1', 'Ejection_one')
				header = header.replace('Ejection2', 'Ejection_two')
				header = header.replace('cooling1', 'cooling_one')
				header = header.replace('cooling2', 'cooling_two')
				header = re.sub(r':\s+:', ': ' + 12 * ' ' + 'NOT_USED :', header)
				header = re.sub(r'\d+\s:', '', header)
				header = re.sub(r'\s\s:\s0.\d{3}', '', header)
				header = re.sub(r'\s\s:\s\d+.\d{3}', '', header)
				header = header.replace(' ', '')
				header = header.replace(',', '\n')
				header = re.sub(' : ', '=', header)
				header = re.sub('=\n', '=-\n', header)
				# header = header.replace('BROKEN', '')
				parser = RawConfigParser(allow_no_value=True, strict=False)
				parser.read_file(StringIO(header))
				pars = {}
				for section in parser.sections():
					for arg in parser.options(section):
						if (section == 'Clean') or (section == 'MCA') or (section == 'Excit') or ('SCAN' in section):
							pars[section.lower() + '_' + arg] = parser.get(section, arg)
						elif (section == 'DATAFILENAME'):
							pars['filename'] = parser.get(section, arg)
						else:
							pars[arg] = parser.get(section, arg)
				return pars
		except(Exception, FileNotFoundError) as err:
			return -1

	def __month_mapping(self, s):
		'''
		Takes string of first three letters of month name and mapps it to two-digit number of month in gregorian calendar
		'''
		mapping ={
			'Jan': '01',
			'Feb': '02',
			'Mar': '03',
			'Apr': '04',
			'Mai': '05',
			'Jun': '06',
			'Jul': '07',
			'Aug': '08',
			'Sep': '09',
			'Oct': '10',
			'Nov': '11',
			'Dec': '12',
		}
		return mapping[s]

	def __add_dt_fields(self):
		"""
		Add field 'dt' to all entries in 'self.header_dict' as datetime converted time of measurement start
		"""
		for mm_file in self.header_dict:
			file = self.header_dict[mm_file]['filename'].split("_")[1].split(".")[0].split("\\")[-1]
			#
			date_string = file[3:]
			# 
			month = self.__month_mapping(file[0:3])
			year = date_string[-4:]
			day = date_string[:2]
			hour = date_string[2:4]
			minute = date_string[4:6]
			second = date_string[6:8]
			dt = datetime.datetime.strptime(f"{month}/{day}/{year} {hour}:{minute}:{second}", '%m/%d/%Y %H:%M:%S')
			self.header_dict[mm_file]['datetime'] = dt

	def process(self, files, to_csv=True):
		"""
		
		"""
		#
		self.files = files
		#
		for i, f in enumerate(self.files):
			# Parse headers 
			header = self.__parse_header(f)
			if header == -1:
				print(f"(MM8.process): WARNING: file {f} not found or empty?")
				continue
			#
			name = os.path.basename(f).split('.')[0]
			#
			self.header_dict[f] = header
		# Add field for time of measurement start as datetime formate
		self.__add_dt_fields()
		#
		return self.header_dict
