#!/usr/bin/env python
#####################################################
# How to run this script
#
# From Linux Terminal:$./calculon.py Results_File Symbol_Ref1 Symbol_Ref2 Atomic_Number Output_File_Name
#
# 'Symbol_Ref1' and 'Symbol_Ref2' can take the following format --> '40Ca+19F' or '85Rb' (without the quotes)
#
# v1.0 Created  by D. Atanasov @ 16.09.2014
# v2.0 Modified by D. Atanasov @ 14.02.2017
# v3.0 Modified to have OOP by D. Atanaso @ 21.03.2017
# v4.0 Modified to import full AME table by J. Karthein @ 11.04.2017
#####################################################
import os
import sys
import re
import numpy as np
import pandas as pd
import platform


#   ######################################################
#   #              Function's Definition                 #
#   ######################################################
class Ame():
    def __init__(self):
        self.scale = 931494.061                          # [amu]->[keV]  #
        self.me = 548.57990946                           # Electron mass #
        self.me_unc = 0.00000022                         # uncertainty   #
        self.load_table()

    def get_number_symbol(self, item):
        """
        Get the Atomic number and the Element symbol from a specific string ()
        Examples: (40Ca19F), (H1H1:O16), 85Rb, 136Cd etc.
        :param item:
        :return: Atomic number, Element Symbol
        """
        return re.findall('[0-9][0-9][0-9]|[0-9][0-9]|[0-9]', item), re.findall('[A-Z][a-z]|[A-Z]', item)

    def load_table(self, ame='AME16.txt'):
        '''
        The function reads the AME table and returns a DataFrame of str's.

        Input name:     'AME16.txt'
        Further info:   If a new table will be uploaded, please change the file + name.
                        Before starting please remove all '*' by ' ' and all '#' by '.'
        Questions to:   jonas.karthein@cern.ch
        '''
        if platform.system() == 'Windows':
            os.chdir('G:\\Experiments\\ISOLTRAP\\Software\\PI-ICR\\Python-DAQ')
        elif platform.system() == 'Darwin': # MacOS
            # os.chdir('/Volumes/dfs/Software/PI-ICR/Python-DAQ/')
            os.chdir('/Volumes/Python-DAQ/')
        ame_import = np.genfromtxt(ame, skip_header=39, dtype=['a1', 'int', 'int', 'int', 'int', 'a4', 'a4', 'float', 'float', 'float', 'float', 'a3', 'float', 'float', 'int', 'float', 'float'], delimiter=[1,3,5,5,5,4,4,16,12,12,6,3,11,8,4,13,11])
        ame_table = [list(elem) for elem in ame_import.tolist()]    # convert list of tuples to list of lists

        for i in range(len(ame_table)):
            ame_table[i][15] = ame_table[i][14] * 1E6 + ame_table[i][15]    # calculate full atomic mass (int*1000000+float)
            for j in range(len(ame_table[i])):
                if type(ame_table[i][j]) == str:
                    ame_table[i][j] = ame_table[i][j].replace(" ", "")      # delete all spaces --> makes searching easier
                ame_table[i][j] = str(ame_table[i][j])                      # convert all entries to int (needed for Dinkos class)
        self.df = pd.DataFrame(ame_table, columns=['cc', 'NZ', 'N', 'Z', 'A', 'EL', 'o', 'mass excess / keV', 'mass excess unc / keV', 'binding energy / keV', 'binding energy unc / keV', 'B', 'beta decay energy / keV', 'beta decay energy unc / keV', 'atomic mass (int) / Dalton', 'atomic mass / micro Dalton', 'atomic mass unc'])


    def get_ame_mass(self, el_expr):
        """
        Checks AME Table for existing of the element with the given atomic number.
        Calculates the mass. If a list is provided to the function (such as a molecule) it calculates
        the summed mass of the constituents.
        :param atomic_number:
        :param symbols:
        :param ame:
        :return: el_expr, atom_mass, np.sqrt(unc)
        """
        ame_mass = 0.0
        ame_unc = 0.0
        self.idx = []
        atomic_number, symbols = self.get_number_symbol(el_expr)
        for i in range(len(symbols)):
            self.idx.append(self.df[(self.df['EL'] == symbols[i]) & (self.df['A'] == atomic_number[i])].index.tolist())
        if len(self.idx[0]) == 0:
            print "Element or Atomic number not found in database"
        else:
            for i in range(len(self.idx)):
                ame_mass += float(self.df.get_value(self.idx[i][0], 'atomic mass / micro Dalton'))
                ame_unc += float(self.df.get_value(self.idx[i][0], 'atomic mass unc'))**2
        return el_expr, ame_mass, np.sqrt(ame_unc)

    def get_ion_mass(self, mass_atom, mass_atom_unc, charge=1):
        """
        Calculate the ion mass
        :param mass_atom:
        :param mass_atom_unc:
        :param charge:
        :return: massIon, massIonUnc
        """
        return (mass_atom - charge*self.me), np.sqrt(mass_atom_unc**2 + charge*self.me_unc**2)

    def get_isobars(self, atomic_number):
        """
        Checks AME Table for existing of the element with the given atomic number.
        Calculates the mass. If a list is provided to the function (such as a molecule) it calculates
        the summed mass of the constituents.
        :param atomic_number:
        :param symbols:
        :param ame:
        :return: el_expr, atom_mass, np.sqrt(unc)
        """
        idx = []
        print atomic_number, type(atomic_number)
        print self.df[self.df['A'] == str(atomic_number)]
        idx = self.df[self.df['A'] == str(atomic_number)].index.tolist()
        print len(idx)
        isobars = np.empty((len(idx), 3,), dtype=object)
        if len(idx) == 0:
            print "No isobars found! Please check the atomic Number"
        else:
            for i in range(len(idx)):
                print idx[i], type(idx[i])
                isobars[i][0] = self.df.get_value(idx[i], 'EL')
                isobars[i][1] = float(self.df.get_value(idx[i], 'atomic mass / micro Dalton'))
                isobars[i][2] = float(self.df.get_value(idx[i], 'atomic mass unc'))
        return isobars


if __name__ == "__main__":
    main(sys.argv)
