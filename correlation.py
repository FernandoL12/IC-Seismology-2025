#! /usr/bin/env python3

###############################
### The CC correlation code ###
###############################

## Libraries
import argparse
import numpy as np
import seaborn as sns
from obspy.clients import fdsn
import matplotlib.pyplot as plt
from obspy.core import UTCDateTime, AttribDict, Stream

# For correlation
from scipy.signal import correlate, correlation_lags
from obspy.signal.cross_correlation import xcorr_pick_correction

import warnings
warnings.filterwarnings("ignore")
#_______________________________________________________________________

###############
## Functions ##
###############

## 1) Input ____________________________________________________________
def cmdline():
    parser = argparse.ArgumentParser(
                    prog='CC',
                    description='Cross-correlation code',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # General Parameters
    g0 = parser.add_argument_group('General')
    g0.add_argument('-v', '--verbose', action='store_true' )
    g0.add_argument('-s', '--station', type = str,
                        help = 'Indicate station to use.')
    g0.add_argument('events', nargs='+', help='Id of one or more events to process.')

    
    # Correlation parameters
    g1 = parser.add_argument_group('Correlation Parameters')
    g1.add_argument('-c', '--correct', action = 'store_true',
                        help = 'Perform a correction shift of given seconds to search for best alignment.')
    g1.add_argument('-cs', '--correction-shift', default = 2.0, type = float,
                        help = 'Amount of shift allowed while searching for the max correlation in a given window.')
    g1.add_argument('-w', '--window', default = 'P/1/2', type = str,
                        help = 'Window size to perform the correlation. Format is Phase/pre[/pos]. Pre and Pos are positive numbers, Phase is one of P or S. [] are optional.')

    
    # Filters
    g2 = parser.add_argument_group('Data Processing Parameters')
    g2.add_argument('-lp', '--low-pass', type = float, 
                        default = 8.5, help = 'Low pass filter value [Hz].')
    g2.add_argument('-hp', '--high-pass', type = float, 
                        default = 1.2, help = 'High pass filter value [Hz].')

    
    # FDSN server sources
    g3 = parser.add_argument_group('Data Fetching parameters')
    g3.add_argument('-F', '--fdsn', type = str, 
                        default = 'http://seisvl.sismo.iag.usp.br/', help = 'FDSN server to fetch data. If no -E option is given this is the default server to fetch events.')
    g3.add_argument('-E', '--event-fdsn', type = str,
                        default = None, help = 'FDSN serve to fetch event. Defaults to same as -F')

    # Process
    args = parser.parse_args()

    return args


## 2) Data Extraction __________________________________________________
#	2.1) Get phase time pick and waveform_id.id (net, sta, loc, chan)
def evpicks(evid, phases = ['P'], station = None, client = 'http://10.110.0.135:18003'):
    '''
    str, str (default=SeisVL) list['str'] --> str, float

    Connect to SeisVL and get the event corresponded to ID given.
    After that, gets the first event preferred origin and it's arrivals and
    yield station code (net.sta.loc.chan) and the time pick of phases  list
    given.
    '''
    
    cl_e = fdsn.Client(client) if type(client) == str else client

    #Confere se o ID existe na SeisVL
    try:
        evp = cl_e.get_events(eventid = evid, includearrivals=True)
    except:
        print(f'No event with {evid} ID detected. Make sure it is in Vale event list (val2025....)!') # None
        return [(None, None)]

    #Pega o 1o evento e sua origem preferida
    E = evp[0]
    O = E.preferred_origin()

    # Faço um laço para as chegadas, pois delas eu tiro a fase que eu quero + ID
    # Com o ID de A eu comparo com o ID do Pick, se bater eu fico com o Pick (que carrega
    # informação do tempo e da waveform ID, que é o código da estação)
    all_data = []
    for A in O.arrivals:
        if A.time_weight <= 0.0: continue
        if A.phase not in phases: continue
        P = [ P for P in E.picks if P.resource_id == A.pick_id ][0]
        
        if station is not None and station == P.waveform_id.id:
            return [ (P.waveform_id.id, P.time) ] 

        all_data.append((P.waveform_id.id, P.time))

    if station is not None:
        raise Exception("Station {} not found.".format(station))

    return all_data


#	2.2) Get trace information and UTC and Relative time of each sample
def evtrace(sta, tp, t0, t1, fmin = 2.0, fmax = 10.0, margin = 2.0, client = None):
    '''
    str, float, float, float, float (opt), float (opt), float (opt), str (opt) --> numpy.ndarray, numpy.ndarray

    Connect to SeisVL and get the waveform of a given period of time, returning the 
    waverform and the relative time.
    '''

    ## Clients
    cl_e = fdsn.Client(client) if type(client) == str else client
    data = cl_e.get_waveforms(*sta.split("."),
                              starttime = t0 - margin,
                              endtime = t1 + margin,
                              attach_response = True)
    
    ## Errors
    if len(data) == 0:
        raise Exception("No data for station {} at time t = {}-{}".format(sta, t0, t1))
    elif len(data) > 1:
        raise Exception("Data with gaps for station {} at time t = {}-{}".format(sta, t0, t1))
    
    ## Process
    trace = data[0]
    trace.detrend()
    trace.filter("bandpass", freqmin=fmin, freqmax=fmax)
    trace.trim(t0,t1)
    
    return trace


## 3) Processing _______________________________________________________
#	3.1) Pick correction function
def Ppick_cc(trace1, trace2):
    """
    Trace, Trace, float --> list, list, float

    Utilizes cross-correlation between two given traces to 
    correct P phase pick.
    """
    
    dt = trace1.stats.delta
    corr = correlate(trace1.data, trace2.data, mode='valid')
    lags = correlation_lags(len(trace1.data), len(trace2.data), mode='valid')
    OFFSET = lags[corr.argmax()] * dt
    
    return corr, lags, OFFSET
    
    
#	3.2) Trim data with the same amount of samples needed
def npts_cut(tr, t0, length = 2, npts = None):
    '''
    Trace, UTCDateTime, int, int, float --> Trace
    
    tr     = trace
    t0     = UTCDateTime (para o cut)
    length = tempo em segundos de dados
    npts   = tempo em número de amostras desejado
    '''

    # Acha o tempo da amostra mais próxima do tempo dado = t0
    t0real = tr.times("utcdatetime")[(np.abs(tr.times("utcdatetime")-t0).argmin())]-tr.stats.delta/2.0

    # Antes de cortar faz uma copia para não destruir o traco original
    trc = tr.copy()

    # Corta exatamente npts amostras
    if npts is None:
        trc.trim(t0real, t0real + length, nearest_sample=False)    
    else:
        trc.trim(t0real, t0real + npts * tr.stats.delta, nearest_sample=False)
   
    return trc    


# 3.3) Builds post correlation correction matrix
def corr_matrix(ev_id, station, phase, fmin, fmax,
                start1, start2, end1, end2,  maxshift, 
                correction, cl_e, cl_d):
    """
    Creates a corrected P pick correlation matrix.

    ev_id: list of events IDs.

    station: station code (example: VL.SLBO..HHZ).

    t0: initial time to trim data in seconds [s].

    phase: list with the phase you would like to make a graph.

    fmin: minimum frequency [Hz] to cut (high-pass).
    fmax: maximum frequency [Hz] to cut (low-pass).

    start1: beginning of first trace (must be the smaller one). 
            Time in seconds [s] before P arrival.
    start2: beginning of second trace (must be the bigger one).
            Time in seconds [s] before P arrival.

    end1: end of first trace (must be the smaller one).
          Time in seconds [s] after P arrival.
    end2: end of second trace (must be the bigger one).
          Time in seconds [s] after P arrival.

    maxshift: maximum shift allowed in seconds [s].

    n_samples: number of samples to cut traces.

    correction: condition if the finction should make a 
    cross-correlation pick correction.
    """
    
    size = len(ev_id)
    
    # Criando duas listas para arnazenar os lag e as formas de ondas pós correlação
    results = []
    
    # Laço que passa por todos os eventos em ev_id e faz a leitura e correção do pick de P para as réplicas que foram
    # registradas pela variável "station"
    for i in range(size):
        # Passa a lista ev_id e tira os dados da "station"
        for s1, t1 in evpicks(ev_id[i], phases=phase, client = cl_e):
            if s1 == station:
                data1 = evtrace(station, t1, t1-start1, t1+end1, fmin, fmax, client = cl_d)
                dt = data1.stats.delta
    
                # Passa a lista ev_id e tira as formas de ondas dos demais eventos que aconteceram em "station"
                for j in range(i + 1, size):
                    for s2, t2 in evpicks(ev_id[j], phases=phase, client = cl_e):
                        if s2 == station:
                            # Aqui vai ser feito uma análise com um chute inicial
                            data2 = evtrace(station, t2, t2-start2, t2+end2, fmin, fmax, client = cl_d)

                            OFFSET_CORR = 0.0
                            OFFSET= 0.0
                            lags = None
                            corr = None

                            if correction:
                                corr, lags, OFFSET = Ppick_cc(data1, data2)
                                FACTOR1 = 1/np.max(data1.data)
                                FACTOR2 = 1/np.max(data2.data)
                                OFFSET_CORR = (t1 - data1.times('utcdatetime')[0]) - (t2 - data2.times('utcdatetime')[0] + OFFSET)

                            data1 = npts_cut(data1, t0 = t1 - start1, length = (start1 + end1))
                            data2 = npts_cut(data2, t0 = t2 - start1 + OFFSET_CORR, npts = data1.stats.npts)

                            print(f"Append i= {i} j={j} evA={ev_id[i]} evB={ev_id[j]} OFFSET={OFFSET:+5.2f} OFFSET_COR={OFFSET_CORR:+5.2f} {'!' if OFFSET_CORR > maxshift else ''}")

                            results.append(AttribDict({
                                'i': i,
                                'j': j,
                                'data1': data1,
                                'data2': data2,
                                'lags': lags,
                                'corr': corr,
                                'OFFSET': OFFSET,
                                'OFFSET_CORR': OFFSET_CORR,
                                's1':s1,
                                's2':s2,
                                't1':t1,
                                't2':t2,
                                'eid1': ev_id[i],
                                'eid2': ev_id[j],
                                'M': np.abs(np.corrcoef(data1.data, data2.data)[0][1])
                            }))

    return results


def assembly_matrix(results):
    size = max(max([ r.i for r in results ]), max([ r.j for r in results ])) + 1
    Mcorr = np.ones([size, size])

    for r in results:
        # ~ print("i=",r.i, "j=", r.j)
        Mcorr[r.i][r.j] = r.M
   
    for i in range(size):
        for j in range(size):
            if j>i: Mcorr[j][i] = -1

    return Mcorr


## 4) Visualização _____________________________________________________
# 	4.1) Plot a heatmap using a given matrix
def plot_matrix(corr_M, ev_id, figsize=(8,6), cmap="Accent_r"):
    """
    matrix, list, tuple, string --> heatmap
    
    Plot a heatmap using a given matrix
    """
    
    size = len(ev_id)
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize)

    ### Inicial
    Max = np.max(np.abs(corr_M))
    Min = -Max
    
    sns.heatmap(
        corr_M,
        cmap=cmap,
        vmin=Min,
        vmax=Max,
        annot=True,
        ax=ax,
        xticklabels=ev_id,
        yticklabels=ev_id,
        cbar_kws={'label': 'Correlation value'}
        )
    
    ax.figure.axes[-1].yaxis.label.set_size(12)
    ax.set_title(f'Correlation matrix\nStation code: {station} | Nº of events: {size}', fontsize=16)
    ax.tick_params(axis="x", rotation=20, labelsize=12)
    ax.tick_params(axis="y", labelsize=12)
    
    plt.show()
    # ~ plt.savefig("Correlation matrix", *, transparent=None, dpi='figure', 
            # ~ format=None, metadata=None, bbox_inches=None, pad_inches=0.1,
            # ~ facecolor='auto', edgecolor='auto', backend=None,
           # ~ )
           

#	4.2) Print da matriz
def print_matrix(matrix, ev_id):
    """
    matrix, list --> string
    
    Print the matrix values in string format
    """
    
    size = len(ev_id)

    print('Matriz de correlação:')
    for i in range(size):
        print('|', end='   ')
        for j in range(size):
            if matrix[i,j] < 0:
                print(f'{matrix[i,j]:.2f}', end='   ')
            else:
                print(f'{matrix[i,j]:.2f}', end='    ')
        print('|', end ='')        
        print()
    print()


#	4.3) Plot graphs
def plot_graph(info_dict, ncols=5, figsize=(30,10)):
    """
    dict, int, tuple --> graph
    
    Plota os gráficos das formas de ondas sobrepostas pós correlação
    """
    # O número de gráficos deve ser igual ao número de correções que foram feitas
    total_de_graficos = len(info_dict['corr'])
    
    # Gostaria que tivesse, no máximo, 5 colunas. Se tiiver menos de 5 gráficos, uma coluna
    # para cada gráfico
    if total_de_graficos < 5:
        ncols = total_de_graficos
    
    # Ajustando o número de linhas de acordo com o npumero de colunas escolhidas (no caso 5)
    if total_de_graficos % ncols == 0:
        nrows = total_de_graficos//ncols
        deno = ncols
    else:
        nrows = total_de_graficos//ncols + 1
        deno = ncols + 1
    
    # Criando os gráficos
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, squeeze=False)
    
    
    # Criando um contador que passará pelo indice do info_dict
    cont = 0
    for row in range(nrows):
        for col in range(ncols):
            title = info_dict['title'][cont]
            image = ax[row][col]
            # Plot da correlações
            #ax[row][col].plot(info_dict['lags'][cont], info_dict['corr'][cont], color='red')
            
            # Plotando o 1o evento
            Y1 = info_dict['data1'][cont]
            X1 = Y1.times()
            Y1 = Y1.data / np.max(Y1.data)
            
            image.plot(X1, Y1, color='darkorange', label=f'{title[:11]}')
            
            # Plotando o 2o evento
            Y2 = info_dict['data2'][cont]
            X2 = Y2.times() + info_dict['offset'][cont]
            Y2 = Y2.data / np.max(Y2.data)
    
            image.plot(X2, Y2, 'b--', label=f'{title[18:]}')
    
            # Configurando o gráfico
            image.set_title(title, fontsize=14)
            image.set_xlim(0.0,2.5)
            image.grid(alpha=0.25)
            image.legend(fontsize=12)
            cont += 1
            
            if cont == total_de_graficos:
                break
    
    
    plt.tight_layout()
       
######################################################################################################
## Código Principal
######################################################################################################
if __name__ == '__main__':
    # 1) Call arguments
    args = cmdline()

    # Client
    client_data  = fdsn.Client(args.fdsn)
    client_event = fdsn.Client(args.event_fdsn) if args.event_fdsn is not None else client_data
    
    # Events ID
    events = args.events
    
    # Station
    station = args.station
    
    # Phase (P or S), cut pre and after phase pick
    phase, pre, pos = args.window.split("/")
    pre = float(pre)
    pos = float(pos)

    # Band-pass
    fmin = args.high_pass
    fmax = args.low_pass

    # Correlation parameters
    # Correct pick
    correct = args.correct
    
    # Correction shift 
    corr_shift =args.correction_shift

    # 2) Processing & 3) Correlation 

    if phase == "P":
        results=corr_matrix(events, station, phase=phase, fmin=fmin,
                          fmax=fmax,start1=pre, start2=pre+corr_shift,
                          end1=pos, end2=pos + corr_shift,
                          maxshift = corr_shift,
                          correction = correct,
                          cl_e = client_event, cl_d = client_data)

    elif phase == "S":
        results=corr_matrix(ev_id, station, t0=0, phase=phase, fmin=fmin,
                          fmax=fmax, start1=pre, start2=pre, 
                          end1=pos, end2=pos+2*shiftmargin,
                          maxshift=pre + 4*shiftmargin, n_samples=200,
                          correction=correct)
                          

   # Gerar os resultados
    Mcorr = assembly_matrix(results)
    plot_matrix(Mcorr, events)
        
