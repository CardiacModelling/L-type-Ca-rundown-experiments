import numpy as np

def selected_cells(temp, hold_dur):
    """
    Returns an array of the selected cells for a given temeprature and holding duration
    """
    path = 'output/'

    pathf = path + f'{temp}_{hold_dur}/_cell_qc.txt'
    f = open(pathf, 'r')
    lines = f.readlines()
    total = len(lines) - 1
    ind_qc1_0 = 0
    ind_qc1_1 = 0
    ind_qc1_2 = 0
    ind_qc1_3 = 0
    ind_qc2_1 = 0
    ind_qc2_2 = 0
    ind_qc3_1 = 0
    ind_qc4 = 0
    ind_qc5 = 0
    cells = []
    
    for i in range(1, len(lines)):
        qc1_0 = int(lines[i][5:6])
        qc1_1 = int(lines[i][7:8])
        qc1_2 = int(lines[i][9:10])
        qc1_3 = int(lines[i][11:12])
        qc2_1 = int(lines[i][13:14])
        qc2_2 = int(lines[i][15:16])
        qc3_1 = int(lines[i][17:18])
        qc4 = int(lines[i][19:20])
        qc5 = int(lines[i][21:22])
        ind_qc1_0 += qc1_0
        ind_qc1_1 += qc1_1
        ind_qc1_2 += qc1_2
        ind_qc1_3 += qc1_3
        ind_qc2_1 += qc2_1
        ind_qc2_2 += qc2_2
        ind_qc3_1 += qc3_1
        ind_qc4 += qc4
        ind_qc5 +=qc5
        selected = qc1_0 * qc1_1 * qc1_2 * qc1_3 * qc2_1 * qc2_2 * qc3_1 * qc4 * qc5
        if selected == 1:
            cells.append(lines[i][:3])
    return cells



def select_cell_per_qc():
    """
    Outputs the number of cells selected from each chip
    """
    path = 'output/'
    temperature = ['BT', 'RT']
    holding_duration = [10, 20, 40]
    print('temp', 'hold |', 'QC0', 'QC1.1', 'QC1.2', 'QC1.3', 'QC2.1', 'QC2.2', 'QC3', 'QC4', 'QC5', '| SELECTED', )
    print('-----------------------------------------------------------------')
    for temp in temperature:
        for hold_dur in holding_duration:
            pathf = path + f'{temp}_{hold_dur}/_cell_qc.txt'
            f = open(pathf, 'r')
            lines = f.readlines()
            total = len(lines) - 1
            ind_qc1_0 = 0
            ind_qc1_1 = 0
            ind_qc1_2 = 0
            ind_qc1_3 = 0
            ind_qc2_1 = 0
            ind_qc2_2 = 0
            ind_qc3_1 = 0
            ind_qc4 = 0
            ind_qc5 = 0
            selected = 0
            
            for i in range(1, len(lines)):
                qc1_0 = int(lines[i][5:6])
                qc1_1 = int(lines[i][7:8])
                qc1_2 = int(lines[i][9:10])
                qc1_3 = int(lines[i][11:12])
                qc2_1 = int(lines[i][13:14])
                qc2_2 = int(lines[i][15:16])
                qc3_1 = int(lines[i][17:18])
                qc4 = int(lines[i][19:20])
                qc5 = int(lines[i][21:22])
                ind_qc1_0 += qc1_0
                ind_qc1_1 += qc1_1
                ind_qc1_2 += qc1_2
                ind_qc1_3 += qc1_3
                ind_qc2_1 += qc2_1
                ind_qc2_2 += qc2_2
                ind_qc3_1 += qc3_1
                ind_qc4 += qc4
                ind_qc5 +=qc5
                selected += qc1_0 * qc1_1 * qc1_2 * qc1_3 * qc2_1 * qc2_2 * qc3_1 * qc4 * qc5
                
            print(temp, hold_dur, "|", ind_qc1_0, ind_qc1_1, ind_qc1_2, ind_qc1_3, ind_qc2_1,\
                 ind_qc2_2, ind_qc4, ind_qc5, ind_qc3_1, "|", selected)
            



def simple_beeswarm(y, nbins=None):
        """
        Returns x coordinates for the points in ``y``, so that plotting ``x`` and
        ``y`` results in a bee swarm plot.
        """
        y = np.asarray(y)
        if nbins is None:
            nbins = len(y) // 6

        # Get upper bounds of bins
        x = np.zeros(len(y))
        ylo = np.min(y)
        yhi = np.max(y)
        dy = (yhi - ylo) / nbins
        ybins = np.linspace(ylo + dy, yhi - dy, nbins - 1)

        # Divide indices into bins
        i = np.arange(len(y))
        ibs = [0] * nbins
        ybs = [0] * nbins
        nmax = 0
        for j, ybin in enumerate(ybins):
            f = y <= ybin
            ibs[j], ybs[j] = i[f], y[f]
            nmax = max(nmax, len(ibs[j]))
            f = ~f
            i, y = i[f], y[f]
        ibs[-1], ybs[-1] = i, y
        nmax = max(nmax, len(ibs[-1]))

        # Assign x indices
        dx = 1 / (nmax // 2)
        for i, y in zip(ibs, ybs):
            if len(i) > 1:
                j = len(i) % 2
                i = i[np.argsort(y)]
                a = i[j::2]
                b = i[j+1::2]
                x[a] = (0.5 + j / 3 + np.arange(len(b))) * dx
                x[b] = (0.5 + j / 3 + np.arange(len(b))) * -dx

        return x


