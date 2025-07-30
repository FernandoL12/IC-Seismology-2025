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

    g0.add_argument('-v', '--verbose',
                        action='store_true')

    g0.add_argument('-a', '--action', choices=[ 'matrix' ])

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


if __name__ == '__main__':
    args = cmdline()
