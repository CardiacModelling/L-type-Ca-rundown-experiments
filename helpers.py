import json
import numpy as np
import pints
from scipy import stats
import pandas as pd
import os

class LoadJson():
    def __init__(self, directoryname, pathtodatadirectory):
        """
        Loads json file that contains details of the experimental conditions at 
        which the voltage-clamp was run

        :param directoryname: One of six depending on temperature and thold
        :param pathtodatadirectory: path to the directory
        :return: python dictionary
        """
        json_name = directoryname[:16] + directoryname[19:] + '.json'
        jsonpath = pathtodatadirectory + json_name
        
        with open(jsonpath) as json_file:
            self.json_data = json.load(json_file)

        # file names of all cells
        cell_fnames = os.listdir(pathtodatadirectory)
        self.cell_fnames = [f for f in cell_fnames if f.endswith('_.csv')]

    def all_available_cells(self):
        """
        Returns the array of all well ID recorded
        """
        cells = []
        for filename in self.cell_fnames:
            cells.append(filename[-8:-5])

        return cells #returns all availabe cells in one directory        

    def time_add_drug(self):
        """
        Returns time in seconds at which dug was added
        """
        time = self.json_data['CompoundAddition']['CompStateProt'][1]['Timestamp_s'] 
        return time #s

    def no_of_samples(self):
        """
        Returns the number of samples recorded per sweep
        """
        return self.json_data['TraceHeader']['MeasurementLayout']['NofSamples']

    def sweep_time(self):
        """
        Returns an array with time in seconds at the beginning of each sweep
        """
        return self.json_data['TraceHeader']['TimeScaling']['SweepTime'][0] #s
    
    def sweep_drug_add(self):
        """
        Returns the sweep # at which drug has been added
        """
        sweeptimearray = self.sweep_time()

        for i in range(len(sweeptimearray)):
            if sweeptimearray[i] > self.time_add_drug():
                return i + 1    

    def time_array(self):
        """
        Returns an array of the time at which current was recorded in a sweep (s)
        """
        return self.json_data['TraceHeader']['TimeScaling']['TR_Time']

    def time_difference_sample(self):
        """
        Returns delta time difference (in ms)
        """
        dt = self.json_data['TraceHeader']['TimeScaling']['TR_Time'][1] - \
            self.json_data['TraceHeader']['TimeScaling']['TR_Time'][0]
        return dt * 1000 #Convert to ms

    def time_measured(self):
        """
        Returns the total time for each sweep (in ms)
        """
        t = self.json_data['TraceHeader']['TimeScaling']['TR_Time'][-1]
        return t * 1000 #Convert to ms

    def time_difference_sweep(self):
        """
        Returns the time between two sweeps (s)
        NOTE: small differences in the time between different sweeps
        """
        dt = self.json_data['TraceHeader']['TimeScaling']['SweepTime'][0][1] - \
            self.json_data['TraceHeader']['TimeScaling']['SweepTime'][0][0]
        return dt  # seconds

    def voltage_protocol(self, protocolfname, fig = False):
        """
        Saves the time (ms) and voltage (mV) at protocolfname
        """
        dt = self.time_difference_sample() #ms
        f = open(protocolfname, 'w')
        f.write('Time(ms), Voltage(mV)\n')
        protocol_array = self.json_data['ExperimentConditions']['VoltageProtocol']
        
        for i in range(len(protocol_array)):
            start = protocol_array[i]['SegmentStart_ms'] #ms
            duration = protocol_array[i]['Duration ms'] #ms
            voltage_start = protocol_array[i]['VoltageStart'] #mV
            voltage_end =  protocol_array[i]['VoltageEnd'] #mV
            n_time_samples = round(duration/dt)
            time_array = np.linspace(start, start + duration, n_time_samples, endpoint=False)

            if voltage_start == voltage_end:
                for i in range(n_time_samples):
                    f.write('%s, %s\n' %(round(time_array[i], 1), voltage_start))
                
            else:
                if fig == True: # change if making the file for figure purposes only
                    voltage_end = -90
                voltage_range = np.linspace(voltage_start, voltage_end, n_time_samples)
                for i in range(n_time_samples):
                    f.write('%s, %s\n' %(round(time_array[i], 1), voltage_range[i]))   
        f.close()

    def machine_prop(self, well):
        """
        :param well: String
        :return: three arrays of the Rseal, Cap, and Rseries recorded per sweep
        """
        row = ord(well[0]) - 65
        column = int(well[1:]) - 1
        RSeal = self.json_data['QCData']['RSeal']
        Capacitance = self.json_data['QCData']['Capacitance']
        RSeries = self.json_data['QCData']['Rseries']
        rseal_array = []
        cap_array = []
        rseries_array = []

        for i in range(len(RSeal)):
            rseal_array.append(RSeal[i][column][row])
            cap_array.append(Capacitance[i][column][row])
            rseries_array.append(RSeries[i][column][row])

        rseal_array = [x*pow(10, -9) if x is not None else None for x in rseal_array] #GOhms
        cap_array = [x*pow(10, 12) if x is not None else None for x in cap_array] #pF
        rseries_array = [x*pow(10, -6) if x is not None else None for x in rseries_array] #MOhms

        rseal_array = np.array(rseal_array, dtype=np.float64)
        cap_array = np.array(cap_array, dtype=np.float64)
        rseries_array= np.array(rseries_array, dtype=np.float64)

        return rseal_array, cap_array, rseries_array

class LinearLeakModel(pints.ForwardModel):
    def __init__(self, voltage_array):
        self.voltage_array = voltage_array
    
    def n_parameters(self):
        return 2

    def simulate(self, parameters, _):
        # times not needed here, 
        g, E = parameters[0], parameters[1]
        return g * (self.voltage_array - E)

def linear_leak_fit(allsweeps, voltage_array, n_sweeps):
    Eleak = []
    gleak = []
    for i in range(n_sweeps):
        slope, intercept , _, _, _ = stats.linregress(voltage_array, allsweeps.iloc[:,i])
        Eleak.append(-1* intercept / slope)
        gleak.append(slope)
    
    return np.array(gleak), np.array(Eleak)

def cap_filter(ind_step_start, ind_step_end, dt, df_current):
    """df_current here should be the subtracted current"""
    n_point_ignore = int(1/dt) #Ignore 1 ms after new voltage step
    df_current.iloc[ind_step_start : ind_step_start + n_point_ignore, :] = np.NaN
    df_current.iloc[ind_step_end : ind_step_end + n_point_ignore, :] = np.NaN

    return df_current

def cal_r_rate(rundown, t_arr):
    rundown = rundown.dropna()
    t_arr = t_arr.dropna()

    # Calculate r_rate of rundown
    rdiff = rundown.iloc[1:] - rundown.iloc[:-1].values
    tdiff = t_arr.iloc[1:] - t_arr.iloc[:-1].values 

    r_rate_arr = rdiff.div(tdiff[0]) 
    r_rate_median = 60 * r_rate_arr.median() # per min

    return r_rate_median # per min


def rundwon_shape(r_data, t_data):
    """
    return answer: 1 (linear), 2(saturating), 3(other)
    """
    rundown = r_data.dropna().to_numpy()
    t_arr = t_data.dropna().values.flatten()

    a = stats.linregress(t_arr, rundown)
    b = stats.linregress(np.log(t_arr), rundown)
    
    if b.rvalue > 0.85:
        answer = 'Saturating'
    elif a.rvalue > 0.95:
        answer = 'Linear'
    else:
        answer = 'Other'

    return answer

def inaca_status(cell):
    if 0 < int(cell[1:]) < 13:
        return 'Off'
    elif 12 < int(cell[1:]) < 25:
        return 'On'
    else:
        raise ValueError('INaCa not determined')

def leak_proportion_calcium(temp, voltage):
    """
    default voltage is the protocol
    concentrations in mM, radius in pm
    Ionic radii: https://doi.org/10.1039/TF9646002075
    Dielectric constant: https://doi.org/10.1063/1.4940432
    ensure that the voltage array does not have volt = 0 mV, change to 0.00001 mV
    """
    K_i, K_o, zk, ak = 110, 3.5, 1, 138
    Na_i, Na_o, zna, ana = 9.1, 78.75, 1, 102
    Ca_i, Ca_o, zca, aca = 0, 2.15, 2, 100
    Cl_i, Cl_o, zcl, acl = 15, 89.05, -1, 181
    F_i, F_o, zf, af = 100, 0, -1, 133
    Mg_i, Mg_o, zmg, amg = 0, 1, 2, 72
    Cs_i, Cs_o, zcs, acs = 0, 0.5, 1, 167
    # ATP_i, ATP_o, zatp, aatp = 4, 0 , -1, 700
    # GTP_i, GTP_o, zgtp, agtp = 0.1, 0 , -1, 700
    Io = 0.5 * (4 * Ca_o + K_o + Na_o + Cl_o + F_o + 4 * Mg_o + Cs_o) * 0.001 # adjusted to Molar
    Ii = 0.5 * (4 * Ca_i + K_i + Na_i + Cl_i + F_i + 4 * Mg_i + Cs_i) * 0.001 # adjusted to Molar

    R = 8.314 #J/mol/K
    F = 96.5 #C/mmol
    if temp == 'BT':
        T = 310
    else:
        T = 298
    frt = F/(R*T)

    def eqn_davies(I, zx):
        A = 1.826 * pow(10, 6) / pow(74 * T, 1.5)
        y = - pow(zx, 2) * A * ((pow(I, 0.5)/ (1 + pow(I, 0.5))) - 0.3 * I)
        return pow(10, y)

    gamma_mono_i = eqn_davies(Ii, 1)
    gamma_di_i = eqn_davies(Ii, 2)
    gamma_mono_o = eqn_davies(Io, 1)
    gamma_di_o = eqn_davies(Io, 2)

    for i in range(len(voltage)):
        if voltage.iloc[i] == 0:
            voltage.iloc[i] = 0.00001

    def gradient_mono(X_i, X_o, zx, ax):
        expon = np.exp(- voltage * zx * frt)
        return pow(zx, 2) * voltage * (X_i * gamma_mono_i - X_o *gamma_mono_o * expon)/ (ax * (1 - expon))

    def gradient_di(X_i, X_o, zx, ax):
        expon = np.exp(- voltage * zx * frt)            
        return pow(zx, 2) * voltage * (X_i*gamma_di_i - gamma_di_o * X_o * expon)/ (ax * (1 - expon))

    grad_K = gradient_mono(K_i, K_o, zk, ak)
    grad_Na = gradient_mono(Na_i, Na_o, zna, ana)
    grad_Ca = gradient_di(Ca_i, Ca_o, zca, aca)
    grad_Cl = gradient_mono(Cl_i, Cl_o, zcl, acl)
    grad_F = gradient_mono(F_i, F_o, zf, af)
    grad_Mg = gradient_di(Mg_i, Mg_o, zmg, amg)
    grad_Cs = gradient_mono(Cs_i, Cs_o, zcs, acs)
    # grad_ATP = gradient(ATP_i, ATP_o, zatp, aatp)
    # grad_GTP = gradient(GTP_i, GTP_o, zgtp, agtp)

    total = grad_K + grad_Na + grad_Ca + grad_Cl + grad_F + grad_Mg + grad_Cs #+ grad_ATP + grad_GTP

    return grad_K, grad_Na, grad_Ca, grad_Cl, grad_F, grad_Mg, grad_Cs, total