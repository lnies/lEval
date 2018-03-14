import pandas as pd
import numpy as np


def s2n(df):
    print '\nRandom DataFrame:'
    print df
    mind_these_numbers = [2, 4]                         # test list of numbers to not consider

    for i in range(1, len(df)-1):

        if i not in mind_these_numbers:                 # to filter by indices
            df.loc[i, 's2n'] = (df.loc[i+1, 'a']        # adding the predecessor, ...
                                + df.loc[i, 'a']        # ...the current entry and ...
                                + df.loc[i-1, 'a'])     # ...the successor into a new column

        elif df.loc[i, 'b'] not in mind_these_numbers:  # to filter by entries in list 'b'
            df.loc[i, 's3n'] = (df.loc[i+1, 'a']        # adding the predecessor, ...
                                + df.loc[i, 'a']        # ...the current entry and ...
                                + df.loc[i-1, 'a'])     # ...the successor into a new column

    print '\nModified DataFrame:'
    print df


if __name__ == '__main__':
    df = pd.DataFrame(np.random.randint(low=0, high=10, size=(10, 5)),columns=['a','b','c','d','e']) # random dataframe
    s2n(df)
