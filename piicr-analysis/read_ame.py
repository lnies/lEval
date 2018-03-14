import pandas as pd
import numpy as np
import re

class AME():
    def __init__(self):
        self.ame_cols = ['A+1 INDICATOR', 'N-Z', 'N', 'Z', 'A', 'ELEMENT', 'O', 'MASS EXCESS', 'MASS EXCESS UNC.', 'BINDING ENERGY/A', 'BINDING ENERGY UNC./A', 'DECAY', 'BETA-DECAY ENERGY', 'BETA-DECAY ENERGY UNC.', 'ATOMIC MASS', 'ATOMIC MASS UNC.']
        self.ame_format = ['a1', 'i3', 'i5', 'i5', 'i5', '1x', 'a3', 'a4', '1x', 'f13.5', 'f11.5', 'f11.3', 'f9.3', '1x', 'a2', 'f11.3', 'f9.3', '1x', 'i3', '1x', 'f12.5', 'f11.5']
        self.dtype = {}
        self.widths = []
        self.atomic_mass_correction_flag = False
        self.ame_cols_correction_flag = True
        self.full_ame_df = {}

    def correct_atomic_mass(self):
        '''Corrects for the whitespace within the atomic mass. Must be deleted by hand in the file!'''
        character_sum = 0
        for i in [-4, -3, -2]:
            character_sum += eval(re.findall('[0-9.]+', self.ame_format[i])[0].replace('.', '+1+'))
        del self.ame_format[-2] # atomic mass w/o int(M/u)
        del self.ame_format[-2] # white space
        del self.ame_format[-2] # just int(M/u)
        self.ame_format.insert(-2, 'f16.5')


    def get_widths(self):
        '''Calculates widths of the columns.'''
        n = 0
        x = 0
        for i in self.ame_format:
            if i == '1x':
                self.widths[-1] += 1
                x += 1
            else:
                # re.search() function gets the numbers from a string (being the first int)
                # replace(...) exchanges the digit points by +1+ to add the number of digits before and after +1 for the point
                # eval function calculates an equation from a string and sums up all digits
                self.widths.append(int(re.search(r'\d+', i).group()))
                if ''.join(re.findall('[a-zA-Z]+', i)) == 'a':
                    self.dtype[self.ame_cols[n-x]] = str
                elif ''.join(re.findall('[a-zA-Z]+', i)) == 'i':
                    self.dtype[self.ame_cols[n-x]] = np.int32
                elif ''.join(re.findall('[a-zA-Z]+', i)) == 'f':
                    self.dtype[self.ame_cols[n-x]] = np.float64
            n += 1


    def get_loading_info(self):
        '''Gets format type and header of AME.'''
        if self.atomic_mass_correction_flag != False and self.ame_cols_correction_flag != True:
            fp = open('AME16.txt')
            for i, line in enumerate(fp):
                if i == 23:
                    self.ame_format = line[16:]
                elif i == 37:
                    self.ame_cols = line
                elif i > 37:
                    break
            fp.close()
            self.ame_cols = re.split('  +', self.ame_cols)
            self.ame_format = re.split(',', self.ame_format)
            self.ame_format[-1] = self.ame_format[-1].replace('\n', '')
            self.ame_cols[0] = self.ame_cols[0].replace('1', '')
            self.ame_cols[4] = self.ame_cols[4].replace('EL', 'ELEMENT')
            self.ame_cols[-1] = self.ame_cols[-1].replace('\n', '')
            self.ame_cols.insert(0, 'A+1 INDICATOR')
            self.ame_cols.insert(8, 'MASS EXCESS UNC.')
            self.ame_cols.insert(10, 'BINDING ENERGY UNC./A')
            self.ame_cols.insert(11, 'DECAY')
            self.ame_cols.insert(13, 'BETA-DECAY ENERGY UNC.')
            self.ame_cols.insert(15, 'ATOMIC MASS UNC.')

            self.atomic_mass_correction_flag = False
            self.ame_cols_correction_flag = True


    def ame_to_df(self):
        '''Imports AME to df'''
        self.full_ame_df = pd.read_fwf('AME16.txt', skiprows=39, widths=self.widths, names=self.ame_cols)


    def get_full_df(self):
        '''Returns full AME as dataframe.'''
        # self.get_loading_info()
        self.correct_atomic_mass()
        self.get_widths()
        self.ame_to_df()
        return(self.full_ame_df)


    def get_mass_excess(self, isotope=''):
        '''Returns mass excess of specific mass.'''
        if isotope != '':
            pass
        else:
            print 'Please enter isotope!'


    def get_number_symbol(self, item):
        """
        Get the Atomic number and the Element symbol from a specific string ()
        Examples: (40Ca19F), (H1H1:O16), 85Rb, 136Cd etc.
        :param item:
        :return: Atomic number, Element Symbol
        """
        return re.findall('[0-9][0-9][0-9]|[0-9][0-9]|[0-9]', item), re.findall('[A-Z][a-z]|[A-Z]', item)

if __name__ == '__main__':
    ame = AME()
    print ame.get_full_df()
    print ame.get_mass_excess('85Rb')

