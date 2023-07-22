import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import gridspec
from matplotlib import cm
matplotlib.rc('font', family='arial', size = 8)

import helpers

pathtodatadir = 'raw_data/'
directories = os.listdir(pathtodatadir)
out_f = os.listdir('output/')

# if protocol file not in resources, then create
# Note: Protocol same for all six chips (checked)
if 'protocol.csv' not in os.listdir('resources/'):
    name = 'resources/protocol.csv'
    file = helpers.LoadJson(directories[0], pathtodatadir + directories[0] + '/')
    file.voltage_protocol(name)   

if 'protocol_fig.csv' not in os.listdir('resources/'):
    name = 'resources/protocol_fig.csv'
    file = helpers.LoadJson(directories[0], pathtodatadir + directories[0] + '/')
    file.voltage_protocol(name, fig = True)    

if 'BT_10_sweep_time.csv' not in os.listdir('resources/'):
    for dir in directories:
        file = helpers.LoadJson(directories[0], pathtodatadir + directories[0] + '/')
        sweep_time = file.sweep_time()
        sweep_time = pd.DataFrame(sweep_time)
        sweep_time.to_csv(f'resources/{dir[16:21]}_sweep_time.csv')

# if chip dir not in output, create
for dir in directories:
    if dir[16:21] not in out_f:
        os.mkdir(f'output/{dir[16:21]}')

class ReadData():
    def __init__(self, directoryname):
        """
        Reads data recorded from an individual chip plate (temperature X thold)
        :directoryname: directory name for each chip output
        """
        self.pathtoallcells = pathtodatadir + directoryname + '/'
        self.jsonfile = helpers.LoadJson(directoryname, self.pathtoallcells)

        self.temp = directoryname[16:18] # temperature 
        self.hold_dur = directoryname[19:21] # holding duration
        self.sweep_drug = self.jsonfile.sweep_drug_add() # first sweep drug

        # array of sweep time
        self.sweep_time = pd.DataFrame(self.jsonfile.sweep_time())
        self.n_sweeps = len(self.jsonfile.sweep_time()) # n_sweeps

        # protocol details
        self.dt = self.jsonfile.time_difference_sample() # ms
        self.protocol = pd.read_csv('resources/protocol.csv', delimiter=',')
        self.ind_ramp_start = self.protocol[self.protocol.iloc[:,0] == 110].\
            index.tolist()[0] #Knowldege of protocol needed here
        self.ind_ramp_end = self.protocol[self.protocol.iloc[:,0] == 510].\
            index.tolist()[0] 
        self.ind_step_start = self.protocol[self.protocol.iloc[:,0] == 860].\
            index.tolist()[0]
        self.ind_step_end = self.protocol[self.protocol.iloc[:,0] == 1010].\
            index.tolist()[0]

        # Placeholders to change per cell
        self.rseal, self.cap, self.rseries = None, None, None
        self.all_sweeps = None
        self.all_sweeps_drug_sub = None
        self.drug_sweep_curr = None
        self.gleak, self.Eleak = None, None

    def qc_prop(self, well):
        """
        Loads the Rseal (GOhms), Cap (pF), and Rseries (MOhms) for the cell in {well} 
        :param well: Well ID, e.g. C01
        """
        rseal, cap, rseries =  self.jsonfile.machine_prop(well)
        #this function updates the machine prop per well
        self.rseal, self.cap, self.rseries = rseal, cap, rseries 

    def qc_check(self, cell, f):
        """
        Returns cells selected after quality control
        QC1.0: NaN in rseal, cap, rseries for sweeps upto n_drug + 2
        QC1.1: 0.1 GOhms < Rseal < 8 GOhms for sweeps upto n_drug + 2
        QC1.2: Cap > 1pF for sweeps upto n_drug + 2
        QC1.3: Rseries < 40 MOhms for sweeps upto n_drug + 2
        QC2.1: Checks abs(leak current) < abs(ical) for sweeps upto n_drug - 1
        QC2.2: Checks abs(drug) < abs(Ical) for all sweeps upto n_drug - 1
        QC3: Checks that net current comes in, not go out (n_drug/2)
        QC4: Checks that SNR > 50 for all sweeps (n_drug/2)

        :param outputpath: path to a directory where QC results to be saved
                        path should end with /
        return: array containing well ID of selected cells
        """
        
        f.write('\n')
        qc1, f = self.qc_range(cell, f) # qc1 check for NaN, range 
        qc2_1, qc5 = self.qc_leak(cell)
        f.write(f'{int(qc2_1)}|')
        qc2_2 = self.qc_drug()
        f.write(f'{int(qc2_2)}|')
        qc3_1 = self.qc_ical()
        f.write(f'{int(qc3_1)}|')
        qc4 = self.qc_noise()
        f.write(f'{int(qc4)}|')  
        f.write(f'{int(qc5)}|')

        cell_pass = qc3_1 and qc2_2 and qc2_1 and qc1 and qc4 and qc5
            
        return cell_pass, f

    def qc_range(self, well, f):
        """
        Checks only selected sweeps [n_drug + 2] for qc1.0-1.3
        """
        self.qc_prop(well) # rseal, cap, and rseries are set here
        rseal = self.rseal[:self.sweep_drug + 1] 
        cap = self.cap[:self.sweep_drug + 1]
        rseries = self.rseries[:self.sweep_drug + 1]

        if np.NaN in rseal or np.NaN in cap or np.NaN in rseries:
            f.write(f'{well} |{0}|n/a|n/a|n/a|')
            return False, f
        
        if max(rseal) > 8 or min(rseal) < 0.1:
            a = False
        else:
            a = True
        if min(cap) < 1:
            b = False
        else:
            b = True
        if max(rseries) > 40:
            c = False
        else:
            c = True
        
        f.write(f'{well} |{1}|{int(a)}|{int(b)}|{int(c)}|')
        return a and b and c, f

    def _read_sweeps(self, cell):
        """
        Loads the raw current for all sweeps (n_drug +2)
        Converts current from A to pA
        Capacitative filtering
        """

        for file in self.jsonfile.cell_fnames:
            if file[-8:-5] == cell:
                filename = file
                break

        cell_data = pd.read_csv(self.pathtoallcells + filename, sep =';', \
            usecols=range(1,self.n_sweeps + 2))

        # Convert A to pA 
        all_sweeps = cell_data.iloc[2:, 1:].astype(float) * pow(10, 12)

        #Capacitative filter
        all_sweeps = helpers.cap_filter(self.ind_step_start, self.ind_step_end,\
                self.dt, all_sweeps)
        self.all_sweeps = all_sweeps.iloc[:, :self.sweep_drug + 1]

    def _voltage_array_leak(self):
        """
        Loads volatge protocol array
        """
        return self.protocol.iloc[self.ind_ramp_start: self.ind_ramp_end, 1]

    def _leak_parameters(self, cell):
        """
        Computes leak parameters for all sweeps n_drug + 2
        """
        self._read_sweeps(cell)

        #Compute leak parameters for each sweep
        leak_data = self.all_sweeps.iloc[self.ind_ramp_start: self.ind_ramp_end, :]
        leak_parameters = helpers.linear_leak_fit(leak_data, self._voltage_array_leak(), \
            len(leak_data.iloc[0]))
        
        self.gleak, self.Eleak = leak_parameters[0], leak_parameters[1]

        return leak_parameters

    def _leak_current_step(self, cell):
        """
        Leak current subtracted traced and drug subtracted traces loaded here
        returns peak value of leak current across all sweeps n_drug + 2
        """

        leak_parameters = self._leak_parameters(cell)
        leak_current_step = [] # calculate the leak current at peak (0 mV)
        for i in range(self.sweep_drug):
            leak_current_step.append(leak_parameters[0][i] * (0 - leak_parameters[1][i]))

        # Post-processing
        leak_current = np.zeros((self.sweep_drug + 1, self.jsonfile.no_of_samples()))

        for i in range(self.sweep_drug + 1):
            leak_current[i] = self._leak_model().simulate((leak_parameters[0][i], \
                leak_parameters[1][i]), 0)

        leak_current = pd.DataFrame(leak_current.transpose())
        subtracted_trace = self.all_sweeps - leak_current.values
        #Capacitative filter
        subtracted_trace = helpers.cap_filter(self.ind_step_start, self.ind_step_end,\
             self.dt, subtracted_trace)

        #Drug Subtraction
        drug_subtracted = np.zeros((self.sweep_drug - 1, self.jsonfile.no_of_samples()))
        #Second sweep after adding drug
        sweep_after_drug = subtracted_trace.iloc[:, -1] 
        self.drug_sweep_curr = sweep_after_drug 
        
        for i in range(self.sweep_drug - 1):
            drug_subtracted[i] = subtracted_trace.iloc[:, i] - sweep_after_drug

        drug_subtracted = pd.DataFrame(drug_subtracted.transpose())
        #Capacitative filter
        drug_subtracted = helpers.cap_filter(self.ind_step_start, self.ind_step_end,\
             self.dt, drug_subtracted)

        self.all_sweeps_drug_sub = drug_subtracted

        # check variation of gleak
        g_leak_arr = pd.DataFrame(leak_parameters[0][:self.sweep_drug -1])
        n_std = g_leak_arr.std()/g_leak_arr.mean()
               
        if n_std[0] > 2:
            qc_gleak_var = False
        else:
            qc_gleak_var = True

        return leak_current_step, qc_gleak_var

    def qc_leak(self, well):
        """
        Checks that magnitude of leak current < extracted ICaL for all sweep before n_drug
        """
        leak_current_step, qc_gleak_var = self._leak_current_step(well) #This will store current
        all_sweeps = self.all_sweeps_drug_sub.iloc[self.ind_step_start: self.ind_step_end, :self.sweep_drug - 1]
        all_sweeps = all_sweeps.min(axis=0)


        for i in range(self.sweep_drug - 1):
            if leak_current_step[i] < all_sweeps.iloc[i]:
                return False, qc_gleak_var
        return True, qc_gleak_var

    def qc_drug(self):
        """
        Note: always has to be called after self.all_sweeps has been initialised
        Checks that abs(drug) < abs(Ical) for all sweeps n_drug - 1
        """
        all_sweeps = self.all_sweeps_drug_sub.iloc[self.ind_step_start: self.ind_step_end, :self.sweep_drug + 1]
        all_sweeps = all_sweeps.min(axis=0)
        min_drug = self.drug_sweep_curr.min(axis=0)
        
        for i in range(self.sweep_drug - 1):
            if min_drug < all_sweeps.iloc[i]:
                return False
        return True

    def qc_ical(self):
        """
        always has to be called after self.all_sweeps has been initialised
        Checks for the first half sweeps that the net current brought in is -ve
        """

        all_sweeps = self.all_sweeps_drug_sub.iloc[self.ind_step_start: self.ind_step_end, :self.sweep_drug -1]

        for i in range(round(self.sweep_drug/2)):
            if all_sweeps.iloc[:,i].sum() > 0:
                return False
        return True

    def qc_noise(self):
        """
        Checks that SNR > 50 for first half sweeps
        """
        i_start = self.ind_step_start
        i_end = self.ind_step_start + int(60/self.dt)
        signal = self.all_sweeps_drug_sub.iloc[i_start: i_end,\
             :int((self.sweep_drug - 1)/2)].min()
        noise = self.all_sweeps_drug_sub.iloc[self.ind_ramp_end: self.ind_step_start, \
            :int((self.sweep_drug - 1)/2)].std()
        snr = abs(signal/ noise)

        if min(snr) < 50:    
            return False
        else:
            return True

    def _leak_model(self):
        leakmodel = helpers.LinearLeakModel(self.protocol.iloc[:, 1])
        return leakmodel
    
    def _cal_tarray(self):
        ## Calcualte the time array
        i_min = self.all_sweeps_drug_sub.iloc[self.ind_step_start: self.ind_step_end, :].idxmin(axis = 0)
        t_stamp = [] # in seconds
        for i in range(len(i_min)):
            if np.isnan(i_min[i]):
                t_stamp.append(np.nan)
            else:
                # add time from sweep to the sweep start time
                t = self.protocol.iloc[int(i_min[i]), 0]/1000 + \
                    self.sweep_time.iloc[i].to_numpy()[0] 

                t_stamp.append(t) # in seconds  
        
        return pd.DataFrame(t_stamp)  
    
    def extract_data(self, df):
        """
        :param file: file which stores all 
        """

        # Access each cell individually
        # File of cells that passed and shape

        outpath = f'output/{dir[16:21]}/'
        cell_array = self.jsonfile.all_available_cells()

        f = open(outpath + '_cell_qc.txt', 'w')
        f.write('well | qc1.0 | qc1.1| qc1.2| qc1.3| qc2.1 | qc2.2 | qc3 | qc4| qc5')

        # f3: selected cells

        for cell in cell_array:
            cell_pass, f = self.qc_check(cell, f)

            if cell_pass == True:
                # save drug subtracted current per sweep
                self.all_sweeps_drug_sub.to_csv(outpath + f'{cell}.csv', index=False)
                
                # save gleak, Eleak, Rseries, Cap, Rseal, per sweep
                leak_parameters = self.gleak, self.Eleak
                art = pd.DataFrame(leak_parameters, index = np.array(['gleak (nS)', 'Eleak (mV)']))
                art = art.transpose()
                art['cap (pF)'] = pd.DataFrame(self.cap.transpose()[:self.sweep_drug + 1])
                art['rseries (MOhms)'] = pd.DataFrame(self.rseries.transpose()[:self.sweep_drug + 1])
                art['rseal (GOhms)'] = pd.DataFrame(self.rseal.transpose()[:self.sweep_drug + 1])
                art.to_csv(outpath + f'prop_{cell}.csv', index = True)

                # Calculate rundown array 
                ical_peak = self.all_sweeps_drug_sub.iloc[self.ind_step_start: self.ind_step_end, :].min(axis = 0)
                rundown = -1*ical_peak/ical_peak.iloc[0]
                
                # Calcualte the time array
                t_arr = self._cal_tarray()

                rrate_val = helpers.cal_r_rate(rundown, t_arr) # rundown rate
                
                # Determine shape of rundown
                shape = helpers.rundwon_shape(rundown*10, t_arr) # to make numerical orders similar

                # plot
                cmap = matplotlib.cm.get_cmap('viridis')(np.linspace(0, 1, len(self.all_sweeps_drug_sub.iloc[0])))
                fig = plt.figure(figsize=(6.6, 3))
                gs = gridspec.GridSpec(ncols =2, nrows = 1)
                ax1 = fig.add_subplot(gs[0, 0])
                ax2 = fig.add_subplot(gs[0, 1])
                ax3 = ax2.twinx()
                time = self.protocol.iloc[:,0].iloc[self.ind_step_start: self.ind_step_end] 
                time_x = time - time.iloc[0]
                for i in range(len(self.all_sweeps_drug_sub.iloc[0])):
                    ax1.plot(time_x, self.all_sweeps_drug_sub.iloc[:, i].iloc[self.ind_step_start: self.ind_step_end], color = cmap[i])
                    ax2.scatter(t_arr.iloc[i], ical_peak.iloc[i], color = cmap[i])
                    ax3.scatter(t_arr.iloc[i], -ical_peak.iloc[i]/ical_peak.iloc[0], color = cmap[i])
                
                ax1.set_xlabel('Time from the beginning of the step to 0 mV\n at each sweep (ms)')
                ax2.set_xlabel('Time (s)')
                ax1.set_ylabel('Current (pA)')
                ax2.set_ylabel('Peak Current (pA)')
                ax3.set_ylabel('Rundown')
                plt.suptitle(f'{cell} rundown shape: {shape}')
                plt.tight_layout()
                plt.savefig(outpath + f'{cell}.png')
                plt.close()


                # save select cell, shape, r_rate
                df['Cell ID'].append(cell)
                df['Run rate'].append(rrate_val)
                df['thold'].append(dir[19:21])
                df['shape'].append(shape)

                if dir[16:18] == 'BT':
                    df['Temperature'].append(310)
                elif dir[16:18] == 'RT':
                    df['Temperature'].append(298)
                else:
                    raise ValueError('Temperature not found')

                inaca_stat = helpers.inaca_status(cell)
                df['INaCa'].append(inaca_stat)
        f.close()

        return df

df = {
    'Cell ID': [],
    'Run rate': [],
    'Temperature': [],
    'INaCa': [],
    'thold': [],
    'shape': [],
}

for dir in directories:
    print(f'Extracting temperature: {dir[16:18]} and thold: {dir[19:21]} s')
    readdata = ReadData(dir)

    df = readdata.extract_data(df)

df = pd.DataFrame(df)
df.to_csv('output/r_rate_database.csv', index=False)

print('Current succesfully extracted')
