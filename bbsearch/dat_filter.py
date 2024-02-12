# -*- coding: utf-8 -*-
import os
import numpy as np
import sys
import shutil
import time
from scipy import signal
from argparse import ArgumentParser

def get_dt_from_inf(inf_file):
    """
    Get the time resolution from the inf file
    """
    dt = -1
    mstr = "time series bin (sec)"
    with open(inf_file, 'r') as fin:
        for line in fin:
            if mstr in line:
                dt = float(line.split('=')[-1])
            else: pass
    return dt


def check_files(datfile):
    """
    Make sure we have both a dat file and an inf file
    """
    retval = 0
    
    basename = datfile.rsplit('.dat', 1)[0]
    inffile = "%s.inf" %basename
    
    if not os.path.exists(datfile):
        print("Missing dat file: %s" %datfile)  
        retval += 1
    
    if not os.path.exists(inffile):
        print("Missing inf file: %s" %inffile)  
        retval += 1

    return datfile, inffile, retval

    
def butter_bandstop(nyq, cutoff_freq_start, 
                    cutoff_freq_stop, order=3):
    """
    Create a butterworth bandstop filter
    """
    cutoff_freq_start = cutoff_freq_start / nyq
    cutoff_freq_stop = cutoff_freq_stop / nyq

    sos = signal.butter(order,
             [cutoff_freq_start, cutoff_freq_stop],
             btype="bandstop", analog=False, output='sos')

    return sos


def parse_zap(zlist):
    """
    Parse the list of f0,nh,df strings given as input
    into arrays of f0, nh, and df
    """
    f0 = []
    nh = []
    ww = []
    for zz in zlist:
        f0_ii, nh_ii, ww_ii = zz.split(',')
        f0.append( float(f0_ii) )
        nh.append( int(nh_ii) )
        ww.append( float(ww_ii) )
    f0 = np.array(f0)
    nh = np.array(nh)
    ww = np.array(ww)

    return f0, nh, ww


def get_freqs(f0, nharms, fnyq):
    """
    get array of freqs to zap
    """
    if f0 > fnyq:
        print("Fundamental above Nyquist")
        return np.array([])
    else: pass

    if nharms == -1:
        nh = int(fnyq / f0)
    else:
        nh = nharms

    freqs = f0 * (1 + np.arange(nh+1))

    print("Filter freqs (Hz):")
    for ff in freqs:
        print("  %0.3f" %ff)
    print("")

    return freqs


def get_filter_freq_ranges(freqs, nharms, widths, fnyq):
    """
    get freq starts and stops for all freq,nh,w 
    """
    fcenters = []
    fstarts  = []
    fstops   = []

    for ii, f0 in enumerate(freqs):
        fc_ii = get_freqs(f0, nharms[ii], fnyq)
        flo_ii = fc_ii - widths[ii]
        fhi_ii = fc_ii + widths[ii] 
    
        fcenters += fc_ii.tolist()
        fstarts += flo_ii.tolist()
        fstops += fhi_ii.tolist()

    fcenters = np.array(fcenters)
    fstarts = np.array(fstarts)
    fstops = np.array(fstops)

    return fcenters, fstarts, fstops


def filter_harms(datfile, dt, freqs, nharms, widths, outbase=None):
    """
    Filter out freq f0 and nharms harmonics.
    from the fb file.

    nharms = 0 : just filter fundamental (f0)
            -1 : filter all harmonics up to Nyquist
             N : filter f0, 2f0, ..., (N+1)f0

    Typical attenuation is ~120-150 dB around the filtered
    frequencies.
    """
    tstart = time.time()
    
    basename = datfile.rsplit('.dat', 1)[0]
    inffile = "%s.inf" %basename

    fnyq = 1 / (2.0 * dt)

    fcs, flos, fhis = get_filter_freq_ranges(freqs, nharms, widths, fnyq)

    if outbase is not None:
        out_dat = "%s_filter.dat" %outbase
        out_inf = "%s_filter.inf" %outbase
    else:
        out_dat = "%s_filter.dat" %basename
        out_inf = "%s_filter.inf" %basename

    # read data
    print("Reading data from: %s" %datfile)
    dat = np.fromfile(datfile, dtype='f4')

    for ii in np.arange(0, len(fcs), 1):
        print("Filtering frequency: %.1f Hz" %(fcs[ii]))
        f_start = flos[ii]
        f_stop  = fhis[ii]
        sos = butter_bandstop(fnyq, f_start, f_stop, order=3)
        dat = signal.sosfiltfilt(sos, dat)

    # Write dat file
    dat = dat.astype('f4')
    dat.tofile(out_dat)

    # Copy inf file with new name
    shutil.copyfile(inffile, out_inf)

    tstop = time.time()
    dt = (tstop - tstart) 
    print("Took %.1f min" %(dt/60.))

    return


def parse_input():
    """
    Use argparse to parse input
    """
    prog_desc = "Filter out frequencies from *.dat file"
    parser = ArgumentParser(description=prog_desc)
    parser.add_argument('datfile', help='Time series *.dat file')
    parser.add_argument('-z', '--zap',
           help='Filter frequency - comma separated list giving the '+\
                'center frequency, the number of harmonics beyond ' +\
                'fundamental, and width in Hz (e.g., \'60.0,5,1.0\' to '+\
                '60Hz signal and 5 harmonics with width 1.0Hz.  To zap '+\
                'multiple frequencies, repeat this argument',
           action='append', required=True, default=[])
    parser.add_argument('-o', '--outbase', default=None,
                        help='Basename for output file ' +\
                             '(def: same as input but with \"_filter\" '+\
                             'appended)',
                        required=False)

    args = parser.parse_args()

    return args


def main():
    """
    Run filter
    """
    tstart = time.time()

    # Parse input
    args = parse_input()

    print("\n\n===== PARAMETERS =====")
    
    datfile = args.datfile
    print("  Input data file: %s" %datfile)
    zaplist = args.zap
    print("  zap: %s" %( "; ".join(zaplist)))
    outbase = args.outbase
    if outbase is not None:
        print("  Output file base name: %s" %outbase)
    else:
        dfn = datfile.rsplit('/', 1)[-1]
        basename = dfn.rsplit('.dat', 1)[0]
        print("  Output file base name: %s" %basename)
    print("======================\n\n")
   
    datfile, inffile, retval = check_files(datfile)

    freqs, nharms, widths = parse_zap(zaplist)
    
    if not retval:
        dt = get_dt_from_inf(inffile)
        print("\nSample Time: %.2f us\n" %(dt * 1e6))
        filter_harms(datfile, dt, freqs, nharms, widths, 
                     outbase=outbase)

    return 


debug = 0

if __name__ == "__main__":
    if debug:
        pass
    else:
        main()
