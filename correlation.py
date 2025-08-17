#! /usr/bin/env python3

#############################
## The CC correlation code.##
#############################

from obspy.clients import fdsn
import argparse

######################################################################################################
## Functions
######################################################################################################

## Utils
def cmdline():
    parser = argparse.ArgumentParser(
                    prog='CC',
                    description='Cross-correlation code',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # General Parameters
    g0 = parser.add_argument_group('General')
    g0.add_argument('-a', '--action' , choices=[ 'matrix' ])
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
    g3.add_argument('-E', '--event-fdns', type = str,
                        default = None, help = 'FDSN serve to fetch event. Defaults to same as -F')

    # Process
    args = parser.parse_args()

    return args


## Processamento
def evpicks(evid, phases = ['P'], station = None, client = None):
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
        return print(f'No event with {evid} ID detected. Make sure it is in Vale event list (val2025....)!') # None

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
            return (P.waveform_id.id, P.time)

        all_data.append((P.waveform_id.id, P.time))

    if station is not None:
        raise Exception("Station {} not found.".format(station))

    return all_data


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
    ##
    if len(data) == 0:
        raise Exception("No data for station {} at time t = {}-{}".format(sta, t0, t1))
    elif len(data) > 1:
        raise Exception("Data with gaps for station {} at time t = {}-{}".format(sta, t0, t1))
    
    ## Process
    ##
    trace = data[0]
    trace.detrend()
    trace.filter("bandpass", freqmin=fmin, freqmax=fmax)
    
    trace.trim(t0,t1)
    
    return trace


def correlate(data, args):
    corr_results = {
        'a': {
            'b' : (corr_sem_deslocar, cor_deslocado, lag_deslocamento),
            'c' : (), 
        },
        'b' : {
                'a': (),
                'c': ()
        }
    }

    for i1,d1 in enumerate(data):
        for i2,d2 in enumerate(data[i1+1:]):
            pass

    return corr_results


def same_length(small, big, t0, offset=0, n_samples=None):
    """
    Trace, Trace, float, list, float(opt) --> Trace, Trace

    Arguments:

    small: obspy.core.trace.Trace (must contain .data, .time and
           .stats.delta)

    big: obspy.core.trace.Trace (must contain .data, .time and
           .stats.delta)

    t0: float (initial time to trim data)

    n_samples: float (number of samples new traces must have)

    n_sec: float (number of seconds that new trace should have)

    Gets two traces and the initial time. Cut both traces with the same
    amount of samples, moving, if necessary, the interval's limits to
    the nearest sample.

    By default, if n_samples is not given, is used the number of samples
    is equal 2s of 'small' seismometer acquisition samples.
    """

    if len(small.data) > len(big.data):
        raise Exception(f"The trace 'small' is bigger than 'big'. "
                        f"Please, switch this arguments position.")
    if len(small.data) == len(big.data):
        return small, big

    dt = small.stats.delta
    if n_samples is None:
        n_samples = int(2 * 1/dt)

    # Trimming SMALL
    # finding initial time to trim (small_ti)
    times_small = small.times()
    find=False
    cont=0
    for i, t in enumerate(times_small[:-1]):
        if len(small) == n_samples:
            break
        if find:
            cont += 1
        if t == t0:
            small_ti = t0
            cont = 1
            find=True
        elif t + dt > t0 and not find:
            if abs(t) - abs(t0) <= abs(times_small[i+1]) - abs(t0):
                small_ti = t
                cont = 1
            else:
                small_ti = times_small[i+1]
                cont = 0
            find=True
        # finding final time to trim (small_tf)
        if cont == n_samples:
            small_tf = t
            break
    if find:
        if cont < n_samples:
            raise Exception(f"{n_samples} samples exceeds number of samples"
                        f" after t0(total={cont}) for trace 'small'")
    if find:
        tr = small.stats.starttime
        small.trim(tr+small_ti, tr+small_tf)

    # Trimming BIG
    # finding initial time to trim (big_ti)
    times_big = big.times()+offset
    start = times_big[0]
    find=False
    cont=0
    for i, t in enumerate(times_big[:-1]):
        if find:
            cont+=1
        if t == t0:
            big_ti = np.abs(t0-start)
            cont = 1
            find=True
        elif t + dt > t0 and not find:
            if abs(t) - abs(t0) <= abs(times_big[i+1]) - abs(t0):
                big_ti = np.abs(t-start)
                cont = 1
            else:
                big_ti = times_big[i+1]
                cont = 0
            find=True
        
        # finding final time to trim (big_tf)
        if cont == n_samples:
            big_tf = np.abs(t-start)
            break    
    if cont < n_samples:
        raise Exception(f"{n_samples} samples exceeds number of samples"
                        f" after t0(total={cont}) for trace 'small'")
    tr = big.stats.starttime
    big.trim(tr+big_ti, tr+big_tf)

    return small, big


# Pick correction function
def Ppick_cc(trace1, trace2):
    """
    Trace, Trace --> list, list, float

    Utilizes cross-correlation between two given traces to 
    correct P phase pick.
    """
    corr = correlate(trace1.data, trace2.data, mode='valid')
    lags = correlation_lags(len(trace1.data), len(trace2.data), mode='valid')
    
    if abs(lags[corr.argmax()] * dt) <= maxshift:
        OFFSET = lags[corr.argmax()] * dt
    else:
        OFFSET = 0

    return corr, lags, OFFSET


## Visualização
def plot_matrix(data, corr_results, args):
    return
    '''
    list, matrix --> print
    
    Receives a list o IDs and an already existing matrix with correlation values.
    It is recommended that the values are normalized by maximum lag.
    '''
    
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
## Código Principal
######################################################################################################
if __name__ == '__main__':
    # Argumentos
    args = cmdline()

    phase, pre, pos = arg.window.split("/")
    pre = float(pre)
    pos = float(pos)

    client_data  = fdsn.Client(args.fdsn)
    client_event = fdsn.Client(args.event_fdsn) if args.event_fdsn is not None else client_data

    # Baixar os dados para os eventos na janela desejada
    data = {}
    for i, evid in enumerate(args.events):
        try:
            shift_margin = args.correction_shift + 1.0 if args.correct else 0.0
            stid, t = evpicks(evid, phases=[ phase ], station = args.station, client = client_event)
            data = evtrace(stid, t, t - (pre + shift_margin) , t + (pos + shift_margin),
                           fmin = args.high_pass, fmax = args.low_pass, client = client_data)
        except E:
            print("Error,", E)
            continue

        data[evid] = {
                'stid'  : stid,
                't_ref' : t,
                'trace' : trace,
            }

    # Processar
    corr_results = correlate(data, args)

    # Gerar os resultados
    if args.actions == 'matrix':
        plot_matrix(data, corr_results, args)
    elif args.actions == 'correct':
        pass
        # ~ plot_wf(data, corr_results, args)
    else:
        raise Exception("Bad mode")

    sys.exit(0)

    # ~ # Fazendo a correção do pick da onda P    
    # ~ # Recriando o gráfico de correlação pós ajuste
    # ~ Mcorr = np.ones([size, size])
    
    # ~ # Criando duas listas para arnazenar os lag e as formas de ondas pós correlação
    # ~ info_dict = {'data1':[], 'data2':[],'lags':[],'corr':[], 'offset' : [], 'title': []}
    
    # ~ # Definindo os parâmetros para evpicks e evtrace
    # ~ # 1) Filtro aplicado
    # ~ fmin = args.low_pass
    # ~ fmax = args.high_pass
    
    # ~ # 2) tamanho do 1o e do 2o traço, em tempo relativo ao pick da onda P
    # ~ start1 = 0.2
    # ~ start2 = 4
    # ~ diff_start = abs(start2 - start1)
    
    # ~ end1 = 2
    # ~ end2 = 2
    
    # ~ maxshift = args.correction_shift # [s]
    
    # ~ # Laço que passa por todos os eventos em ev_id e faz a leitura e correção do pick de P para as réplicas que foram
    # ~ # registradas pela variável "station"
    
    # ~ for i in range(len(ev_id)):
        # ~ # Passa a lista ev_id e tira os dados da "station"
        # ~ for s1, t1 in evpicks(ev_id[i], phases=['P']):
            # ~ if s1 == station:
                # ~ data1, _, times1 = evtrace(station, t1, t1-start1, t1+end1, fmin, fmax)
                # ~ dt = data1.stats.delta
    
                # ~ # Passa a lista ev_id e tira as formas de ondas dos demais eventos que aconteceram em "station"
                # ~ for j in range(i + 1, len(ev_id)):
                    # ~ for s2, t2 in evpicks(ev_id[j], phases=['P']):
                        # ~ if s2 == station:
                            # ~ # Aqui vai ser feito uma análise com um chute inicial
                            # ~ data2, _, times2 = evtrace(station, t2, t2-start2, t2+end2, fmin, fmax)
                            # ~ corr = correlate(data1.data, data2.data, mode='valid')
                            # ~ lags = correlation_lags(len(data1.data), len(data2.data), mode='valid')
                            
                            # ~ if abs(lags[corr.argmax()] * dt) <= maxshift:
                                # ~ OFFSET = lags[corr.argmax()] * dt
                            # ~ else:
                                # ~ OFFSET = 0
                            # ~ FACTOR1 = 1/np.max(data1.data)
                            # ~ FACTOR2 = 1/np.max(data2.data)
                            
                            # ~ info_dict['data1'].append(data1)
                            # ~ info_dict['data2'].append(data2)
                            # ~ info_dict['lags'].append(lags)
                            # ~ info_dict['corr'].append(corr)
                            # ~ info_dict['offset'].append(OFFSET)
                            # ~ info_dict['title'].append(f'{ev_id[i]}  -X-  {ev_id[j]}')
                            
    
                            # ~ # Constrói matriz
                            # ~ t1    = data1.times() 
                            
                            # ~ t2    = data2.times() + OFFSET
                            
                            # ~ index1   = bisect.bisect_left(t2, 0) #next(i for i, t in enumerate(t2) if t >= 0)
                            # ~ index2   = bisect.bisect_left(t2, np.max(t1)) #next(i for i, t in enumerate(t2) if t >= np.max(t1))
                    
                            # ~ if len(data1) != len(data2[index1:index2]):
                                # ~ # print(i, j)
                                # ~ # print(t1)
                                # ~ # print(len(t1))
                                # ~ # print()
                                # ~ # print(len(t2[index1:index2]))
                                # ~ # print(t2[index1:index2])
                                
                                # ~ if len(data2[index1:index2]) < len(data1):
                                    # ~ index2+=1
                                # ~ elif len(data2[index1:index2]) > len(data1):
                                    # ~ index2-=1
                    
                            # ~ try:
                                # ~ Mcorr[i,j] = np.corrcoef(data1, data2[index1:index2])[0][1]
                            # ~ except:
                                # ~ Mcorr[i,j] = 0
    
    # ~ for i in range(size):
        # ~ for j in range(size):
            # ~ if j>i:
                # ~ Mcorr[j][i] = -1

    # ~ # Plota o gráfico (está em uma célula a parte, pois se for mudar a paleta de cores não demora tanto)
    # ~ Max = np.max(np.abs(Mcorr))
    # ~ Min = -Max
    # ~ plt.figure(figsize=(8, 8))
    # ~ plt.imshow(Mcorr, cmap='Accent_r', vmin = Min, vmax = Max)
    # ~ plt.colorbar(label='Correlation value')
    # ~ plt.xticks(range(len(ev_id)), ev_id, rotation=20)
    # ~ plt.yticks(range(len(ev_id)), ev_id)
    # ~ plt.title(f'Correlation matrix (corrected)\nStation code: {station}\nNumber of events: {len(ev_id)}')
    # ~ plt.show()

    # ~ # Efetivamente faz o plot
    # ~ plot_matrix(ev_id, Mcorr)


    # ~ # Costruindo o plot das formas de ondas corrigidas sobrepostas
    # ~ # O número de gráficos deve ser igual ao número de correções que foram feitas
    # ~ total_de_graficos = len(info_dict['corr'])
    
    # ~ # Gostaria que tivesse, no máximo, 5 colunas. Se tiiver menos de 5 gráficos, uma coluna
    # ~ # para cada gráfico
    # ~ if total_de_graficos < 5:
        # ~ ncols = total_de_graficos
    # ~ else:
        # ~ ncols = 5
    
    # ~ # Ajustando o número de linhas de acordo com o npumero de colunas escolhidas (no caso 5)
    # ~ if total_de_graficos % ncols ==0:
        # ~ nrows = total_de_graficos//ncols
        # ~ deno = ncols
    # ~ else:
        # ~ nrows = total_de_graficos//ncols + 1
        # ~ deno = ncols + 1
    
    # ~ # Criando os gráficos
    # ~ fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=[15,8], squeeze=False)
    
    
    # ~ # Criando um contador que passará pelo indice do info_dict
    # ~ cont = 0
    # ~ for row in range(nrows):
        # ~ for col in range(ncols):
    
            # ~ image = ax[row][col]
            # ~ # Plot da correlações
            # ~ #ax[row][col].plot(info_dict['lags'][cont], info_dict['corr'][cont], color='red')
            
            # ~ # Plotando o 1o evento
            # ~ Y1 = info_dict['data1'][cont]
            # ~ X1 = Y1.times()
            # ~ Y1 = Y1.data / np.max(Y1.data)
            
            # ~ image.plot(X1, Y1, color='darkorange')
            
            # ~ # Plotando o 2o evento
            # ~ Y2 = info_dict['data2'][cont]
            # ~ X2 = Y2.times() + info_dict['offset'][cont]
            # ~ Y2 = Y2.data / np.max(Y2.data)
    
            # ~ image.plot(X2, Y2, 'b--')
    
            # ~ # Configurando o gráfico
            # ~ image.set_title(info_dict['title'][cont], fontsize=10)
            # ~ image.set_xlim(0.0,2.5)
            # ~ image.grid(alpha=0.25)
            
            # ~ cont += 1
            
            # ~ if cont == total_de_graficos:
                # ~ break
    
    # ~ plt.tight_layout()
