#! /usr/bin/env python3

#
# The CC correlation code.
#

from obspy.clients import fdsn
import argparse


def cmdline():
    parser = argparse.ArgumentParser(
                    prog='CC',
                    description='Cross-correlation code',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # General Parameters
    g0 = parser.add_argument_group('General')

    g0.add_argument('-v', '--verbose', action='store_true' )
    g0.add_argument('-a', '--action' , choices=[ 'matrix' ])
    g0.add_argument('events', nargs='+', help='Id of one or more events to process.')

    
    # Correlation parameters
    g1 = parser.add_argument_group('Correlation Parameters')
    g1.add_argument('-c', '--correct', action = 'store_true',
                        help = 'Perform a correction shift of given seconds to search for best alignment.')
    g1.add_argument('-cs', '--correction-shift', default = 2.0, type = float,
                        help = 'Amount of shift allowed while searching for the max correlation in a given window.')
    g1.add_argument('-w', '--window', default = 'P/-1/2', type = str,
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
    g3.add_argument('-E', '--event-fdns', type = str,
                        default = None, help = 'FDSN serve to fetch event. Defaults to same as -F')

    # Process
    args = parser.parse_args()

    return args

######################################################################################################
###  FUNCTIONS  ###
# 1) Get phase time pick and waveform_id.id (net, sta, loc, chan)
def evpicks(evid, phases = ['P']):
    '''
    str, list['str'] --> str, float

    Connect to SeisVL and get the event corresponded to ID given.
    After that, gets the first event preferred origin and it's arrivals and
    yield station code (net.sta.loc.chan) and the time pick of phases  list
    given.
    '''
    
    #Connecta na SeisVL
    cl_e = fdsn.Client("http://10.110.0.135:18003") #SeisVL

    #Confere se o ID existe na SeisVL
    try:
        evp = cl_e.get_events(eventid = evid, includearrivals=True)
    except:
        return print(f'No event with {evid} ID detected. Make sure it is in Vale event list (val2025....)!') # None

    #Pega o 1o evento e sua origem preferida
    E = evp[0]
    O = E.preferred_origin()

    # Faço um laço para as chegadas, pois delas eu tiro a fase que eu quero + ID
    # Com o ID de A eu comparo com o ID do Pick, se bater eu fico com o Pick (que carrega
    # informação do tempo e da waveform ID, que é o código da estação)
    for A in O.arrivals:
        if A.time_weight <= 0.0: continue
        if A.phase not in phases: continue
        P = [ P for P in E.picks if P.resource_id == A.pick_id ][0]

        # Yield é tipo um return, mas que continua o loop e não para a função
        yield (P.waveform_id.id, P.time) #waveform_id.id é o código da estação


# 2) Get trace information and UTC and Relative time of each sample
def evtrace(sta, tp, t0, t1, fmin=2.0, fmax=10.0, margin = 2.0):
    '''
    str, float, float, float, float (opt), float (opt), float (opt) --> numpy.ndarray, numpy.ndarray

    Connect to SeisVL and get the waveform of a given period of time, returning the 
    waverform and the relative time.
    '''
    
    ## Clients
    cl_e = fdsn.Client("http://10.110.0.135:18003")
    data = cl_e.get_waveforms(*sta.split("."), starttime=t0 - margin, endtime=t1 + margin, attach_response=True)
    
    ## Errors
    ##
    if len(data) == 0:
        raise Exception("No data")
    elif len(data) > 1:
        raise Exception("Data with gap")
    
    ## Process
    ##
    trace = data[0]
    trace.detrend()
    trace.filter("bandpass", freqmin=fmin, freqmax=fmax)
    trace.trim(t0,t1)
    
    ## Return
    ##
    #OBS: A data tem que ser um obspy.Trace para rodar a correlação depois
    data   = trace #.data
    
    # Estou retornando tempo UTC e Relativo
    # UTC pois precisa para a função de correlação
    # E relativo pelo fato de ser mais conveniente na hora de plotar a waveform 
    times_utc  = trace.times("utcdatetime")
    times_r = times_utc - tp 
    
    return (data, times_utc, times_r)

# 3) Print da matriz
def plot_matrix(ev_id, matrix):
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


    print(f'\n\n')
######################################################################################################

def main:
    args = cmdline()
















if __name__ == '__main__':
    args = cmdline()
    main()