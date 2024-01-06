import os
import pandas as pd

temperature = ['BT', 'RT']
hold_dur = [10, 20, 40]

path = '../output/'

for temp in temperature:
    for hold in hold_dur:

        mypath = path + f'{temp}_{hold}'
        files = os.listdir(mypath)

        df = {}
        for f in files:
            if f[:4] == "prop":
                cell = f[5:8]
                data = pd.read_csv(mypath + f'/{f}')
                gleak = data['gleak (nS)']
                Eleak = data['Eleak (mV)']
                I_fixed = - gleak * Eleak
                df[cell] = I_fixed
        
        df = pd.DataFrame(df)
        df.to_csv(f'I_fixed_{temp}_{hold}.csv', index = None)
        