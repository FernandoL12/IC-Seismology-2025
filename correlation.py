#! /usr/bin/env python3

#
# The CC correlation code.
#

from obspy.clients import fdsn
import argparse


def main:
    args = cmdline()
    # Lista de eventos
    ev_id = args.events ###########
    # Número de estações/dimensão da matriz nxn
    size = len(ev_id)
    # Estação que será extraída a forma de onda
    station = 'VL.SLBO..HHZ' ###########


    # Fazendo a correção do pick da onda P
    
    # Recriando o gráfico de correlação pós ajuste
    Mcorr = np.ones([size, size])
    
    # Criando duas listas para arnazenar os lag e as formas de ondas pós correlação
    info_dict = {'data1':[], 'data2':[],'lags':[],'corr':[], 'offset' : [], 'title': []}
    
    # Definindo os parâmetros para evpicks e evtrace
    # 1) Filtro aplicado
    fmin = args.low_pass
    fmax = args.high_pass
    
    # 2) tamanho do 1o e do 2o traço, em tempo relativo ao pick da onda P
    start1 = 0.2
    start2 = 4
    diff_start = abs(start2 - start1)
    
    end1 = 2
    end2 = 2
    
    maxshift = args.correction_shift # [s]
    
    # Laço que passa por todos os eventos em ev_id e faz a leitura e correção do pick de P para as réplicas que foram
    # registradas pela variável "station"
    
    for i in range(len(ev_id)):
        # Passa a lista ev_id e tira os dados da "station"
        for s1, t1 in evpicks(ev_id[i], phases=['P']):
            if s1 == station:
                data1, _, times1 = evtrace(station, t1, t1-start1, t1+end1, fmin, fmax)
                dt = data1.stats.delta
    
                # Passa a lista ev_id e tira as formas de ondas dos demais eventos que aconteceram em "station"
                for j in range(i + 1, len(ev_id)):
                    for s2, t2 in evpicks(ev_id[j], phases=['P']):
                        if s2 == station:
                            # Aqui vai ser feito uma análise com um chute inicial
                            data2, _, times2 = evtrace(station, t2, t2-start2, t2+end2, fmin, fmax)
                            corr = correlate(data1.data, data2.data, mode='valid')
                            lags = correlation_lags(len(data1.data), len(data2.data), mode='valid')
                            
                            if abs(lags[corr.argmax()] * dt) <= maxshift:
                                OFFSET = lags[corr.argmax()] * dt
                            else:
                                OFFSET = 0
                            FACTOR1 = 1/np.max(data1.data)
                            FACTOR2 = 1/np.max(data2.data)
                            
                            info_dict['data1'].append(data1)
                            info_dict['data2'].append(data2)
                            info_dict['lags'].append(lags)
                            info_dict['corr'].append(corr)
                            info_dict['offset'].append(OFFSET)
                            info_dict['title'].append(f'{ev_id[i]}  -X-  {ev_id[j]}')
                            
    
                            # Constrói matriz
                            t1    = data1.times() 
                            
                            t2    = data2.times() + OFFSET
                            
                            index1   = bisect.bisect_left(t2, 0) #next(i for i, t in enumerate(t2) if t >= 0)
                            index2   = bisect.bisect_left(t2, np.max(t1)) #next(i for i, t in enumerate(t2) if t >= np.max(t1))
                    
                            if len(data1) != len(data2[index1:index2]):
                                # print(i, j)
                                # print(t1)
                                # print(len(t1))
                                # print()
                                # print(len(t2[index1:index2]))
                                # print(t2[index1:index2])
                                
                                if len(data2[index1:index2]) < len(data1):
                                    index2+=1
                                elif len(data2[index1:index2]) > len(data1):
                                    index2-=1
                    
                            try:
                                Mcorr[i,j] = np.corrcoef(data1, data2[index1:index2])[0][1]
                            except:
                                Mcorr[i,j] = 0
    
    for i in range(size):
        for j in range(size):
            if j>i:
                Mcorr[j][i] = -1

    # Plota o gráfico (está em uma célula a parte, pois se for mudar a paleta de cores não demora tanto)
    Max = np.max(np.abs(Mcorr))
    Min = -Max
    plt.figure(figsize=(8, 8))
    plt.imshow(Mcorr, cmap='Accent_r', vmin = Min, vmax = Max)
    plt.colorbar(label='Correlation value')
    plt.xticks(range(len(ev_id)), ev_id, rotation=20)
    plt.yticks(range(len(ev_id)), ev_id)
    plt.title(f'Correlation matrix (corrected)\nStation code: {station}\nNumber of events: {len(ev_id)}')
    plt.show()

    # Efetivamente faz o plot
    plot_matrix(ev_id, Mcorr)


    # Costruindo o plot das formas de ondas corrigidas sobrepostas
    # O número de gráficos deve ser igual ao número de correções que foram feitas
    total_de_graficos = len(info_dict['corr'])
    
    # Gostaria que tivesse, no máximo, 5 colunas. Se tiiver menos de 5 gráficos, uma coluna
    # para cada gráfico
    if total_de_graficos < 5:
        ncols = total_de_graficos
    else:
        ncols = 5
    
    # Ajustando o número de linhas de acordo com o npumero de colunas escolhidas (no caso 5)
    if total_de_graficos % ncols ==0:
        nrows = total_de_graficos//ncols
        deno = ncols
    else:
        nrows = total_de_graficos//ncols + 1
        deno = ncols + 1
    
    # Criando os gráficos
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=[15,8], squeeze=False)
    
    
    # Criando um contador que passará pelo indice do info_dict
    cont = 0
    for row in range(nrows):
        for col in range(ncols):
    
            image = ax[row][col]
            # Plot da correlações
            #ax[row][col].plot(info_dict['lags'][cont], info_dict['corr'][cont], color='red')
            
            # Plotando o 1o evento
            Y1 = info_dict['data1'][cont]
            X1 = Y1.times()
            Y1 = Y1.data / np.max(Y1.data)
            
            image.plot(X1, Y1, color='darkorange')
            
            # Plotando o 2o evento
            Y2 = info_dict['data2'][cont]
            X2 = Y2.times() + info_dict['offset'][cont]
            Y2 = Y2.data / np.max(Y2.data)
    
            image.plot(X2, Y2, 'b--')
    
            # Configurando o gráfico
            image.set_title(info_dict['title'][cont], fontsize=10)
            image.set_xlim(0.0,2.5)
            image.grid(alpha=0.25)
            
            cont += 1
            
            if cont == total_de_graficos:
                break
    
    plt.tight_layout()


######################################################################################################
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
    args = cmdline()
    
    #Connecta na SeisVL
    cl_e = fdsn.Client(args.fdsn) #SeisVL

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
    args = cmdline()
    ## Clients
    cl_e = fdsn.Client(args.fdsn)
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

if __name__ == '__main__':
    args = cmdline()
    main()