#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import glob
import os
import copy
import sys
import time
import your 
from subprocess import call 
from argparse import ArgumentParser
import write_filterbank as wfil 
import snippet_plots_sp as sp_plt

#scriptdir = "/src/bb_proc/bbsearch"
cur_dir = os.path.realpath(__file__)
scriptdir  = cur_dir.rsplit('/', 1)[0]

############################
## Filterbank Parameters ##
############################

def get_chan_info(data_file):
    yr = your.Your(data_file)
    foff = yr.your_header.foff
    fch1 = yr.your_header.fch1
    dt   = yr.your_header.tsamp
    nchans = yr.your_header.nchans

    return nchans, fch1, foff, dt


##############################
## Making and Reading Cands ##
##############################

class SP_CAND:
    def __init__(self, line, num):
        cols = line.split()
        self.dm = float(cols[0])
        self.snr = float(cols[1])
        self.time = float(cols[2])
        self.samp = int(cols[3])
        self.wbins = int(cols[4])
        
        self.cnum = num

    def __str__(self):
        str = "SP_CAND(t=%.2f s, DM=%.1f pc/cc, SNR=%.2f)" %(\
                                  self.time, self.dm, self.snr)
        return str

    def __repr__(self):
        str = "SP_CAND(t=%.2f s, DM=%.1f pc/cc, SNR=%.2f)" %(\
                                 self.time, self.dm, self.snr)
        return str

def cands_from_spfile(spfile):
    """
    Read in the candidates from a *singlepulse (sp) file
    and return an array of SP_CAND class objects
    """
    ii = 1
    candlist = []
    with open(spfile, 'r') as fin:
        for line in fin:
            if line[0] in ["\n", "#", 'i']:
                continue
            else: pass

            SP = SP_CAND(line, ii)
            candlist.append(SP)
            ii += 1

    return np.array(candlist)


##########################
## PROCESSING PIPELINE  ##
##########################

def fix_file(filfile, tel, dsn):
    """
    Check if we need to update the obs or 
    machine key in the header.  if we do, 
    update them.
    """
    if tel is None and dsn==False:
        print("Don't need to fix header")
        return
    else: pass

    fix_cmd = "python %s/fix_sigproc_header.py " %scriptdir +\
              "--filterbank %s " %filfile

    if tel is not None:
        fix_cmd += "--telescope %s " %tel

    if dsn:
        fix_cmd += "--dsn"

    print(fix_cmd)
    call(fix_cmd, shell=True)

    return


def del_chans_to_string(nums):
    """
    Take list of channel numbers to remove and convert
    them to a string that can be input to PRESTO
    """
    # Get list of ranges
    nums = sorted(set(nums))
    gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s+1 < e]
    edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
    ranges = list(zip(edges, edges))

    # Shorten string using ":" when necessary
    out_str = ""
    for i in np.arange(0, len(ranges), 1):
        if (ranges[i][0] == ranges[i][1]):
            out_str = out_str + str(ranges[i][0]) + ","
        else:
            out_str = out_str + str(ranges[i][0]) +\
                      ":" + str(ranges[i][1]) + ","

    # Remove trailing comma if nec
    if out_str[-1] == ',':
        out_str = out_str.rstrip(',')
    else: pass

    return out_str


def get_zap_chans(edgezap, nchansub, nsub, zchans="", flip=True):
    """
    Generate a list of channels to zap from two options:

      1. Edge zapping of subbands 
      2. Input list of channels via a comma sep string

    Returns a string that can be passed to PRESTO
    """
    zap_list = []

    # First we'll get the edge channels to zap
    if edgezap > 0:
        zlo = np.arange(0, edgezap)
        zhi = np.arange(nchansub-edgezap, nchansub)
        for ii in range(nsub):
            zlo_ii = zlo + ii * nchansub
            zhi_ii = zhi + ii * nchansub
            zap_list += zlo_ii.tolist()
            zap_list += zhi_ii.tolist()

    # Now the manually entered channels
    if zchans != "":
        zz_list = [ int(zz) for zz in zchans.split(',') ]
        zap_list += zz_list 

    # Get unique chans, flip for PRESTO, and make string
    if len(zap_list):
        zap_list = np.unique(zap_list)

        if flip:
            Nchans = nsub * nchansub
            zap_list = Nchans - 1 - zap_list[::-1]

        # Get formatted string
        zap_chan_str = del_chans_to_string(zap_list)
        

    else:
        zap_chan_str = ""

    return zap_chan_str


def dedisperse(filfile, dm, zapstr, zdm=False, outdir='.'):
    """
    Run prepdata to dedisperse the filterbank file 
    at some DM, zapping channels if desired

    Include zerodm option if zdm=True

    Will do nothing if file already exists
    """
    filbase = (filfile.split("/")[-1]).split(".fil")[0] 
    outbase =  "%s/%s_DM%.3f" %(outdir, filbase, dm)
    datfile = "%s.dat" %outbase

    if zapstr != "":
        zap = "-ignorechan %s " %zapstr
    else:
        zap = ""

    if zdm:
        zdm_str = "-zerodm "
    else:
        zdm_str = ""

    dm_cmd = "prepdata -filterbank " +\
             "-dm %.3f " %dm +\
             "-nobary -noclip " +\
             "%s" %zdm_str +\
             zap +\
             "-o %s " %outbase +\
             "%s " %filfile

    if os.path.exists(datfile):
        print("  Found dat file: %s" %datfile)
        print("  Skipping dedispersion")
        pass
    else:
        print(dm_cmd)
        call(dm_cmd, shell=True)

    return datfile


def filter_dat(datfile, zaplist):
    """
    Run filter on dat file 
    """
    if len(zaplist) == 0:
        print("Skipping filter of dat file")
        return datfile
    else: pass
    
    zstr = ""
    for zz in zaplist:
        zstr += "-z %s " %zz
    
    filter_cmd = "python -u %s/dat_filter.py " %scriptdir +\
                 "%s %s" %(zstr, datfile)
    
    print(filter_cmd)
    call(filter_cmd, shell=True)

    basename = datfile.rsplit('.dat', -1)[0]
    outdat = "%s_filter.dat" %basename  

    return outdat


def sp_search(datfile, snr, bb=False, maxwidth=1.0, dtrendlen=8):
    """
    Run single pulse search

      snr:  SNR threshold
    
      maxwidth: do not consider templates wider than 
                this maxwidth (in sec).

      detrendlen: detrending length (thousands of bins)
                  Must be power of two.  So a detrendlen
                  of 8 means that detrending will be done 
                  using a size of 8000 bins.

      bb: Do NOT ignore bad blocks if True
    """
    datbase = datfile.split(".dat")[0]
    spfile = "%s.singlepulse" %datbase

    if bb:
        b_str = ""
    else:
        b_str = "-b "

    sp_cmd = "python -u %s/single_pulse_search_w16ms.py " %scriptdir +\
             "-t %.2f " %snr +\
             "-m %.4f " %maxwidth +\
             "-d %d " %dtrendlen +\
             "%s" %b_str +\
             "%s" %datfile
   
    if os.path.exists(spfile):
        print("  Found singlepulse file: %s" %spfile)
        print("  Skipping single pulse search")

    else:
        print(sp_cmd)
        call(sp_cmd, shell=True)

    return spfile


def extract_snippets(filfile, outbase, splist, nspec):
    """
    Extract data around candidates and write to 
    small filterbank file.  Candidates input as 
    list of SPCAND class objects.

    Output data will have nspec time samples
    """
    # Get number of zeros to bad so cand nums
    # are all the same length
    #nz = int( np.log10(len(splist)) + 0.5 )
    nz = int( np.ceil(np.log10(len(splist)+1)) )
      
    for ii, cc in enumerate(splist):
        print("Cand %d / %d" %(ii+1, len(splist)))
        num_str = str(ii).zfill(nz)  
        outfile = "%s_cand%s.fil" %(outbase, num_str)
        nstart = cc.samp
         
        ext_cmd = "extract %s " %filfile +\
                  "%d %d " %(nstart, nspec) +\
                  "> %s" %(outfile)

        print(ext_cmd)
        call(ext_cmd, shell=True)

    return 


def your_extract_snippets(filfile, outbase, splist, 
                          nspec, nmax=-1):
    """
    Extract data around candidates and write to 
    small filterbank file.  Candidates input as 
    list of SPCAND class objects.

    Output data will have nspec time samples
    """
    # Get number of zeros to bad so cand nums
    # are all the same length
    nz = int( np.ceil(np.log10(len(splist)+1)) )

    # Get hdr info
    yr = your.Your(filfile)
    hdr_bytes = yr.hdrbytes
    nchans = yr.nchans
    dt = yr.tsamp 
    bps = yr.nbits // 8  # bytes per sample

    # Sort sp cands by time
    tt = np.array([ ss.time for ss in splist ])
    xx = np.argsort(tt)
    spsort = splist[xx]

    with open(filfile, 'rb') as fin:
        # Skip header bytes 
        #fin.seek(hdr_bytes)

        # Now loop over cands and extract data
        for ii, cc in enumerate(spsort):
            if nmax > 0 and ii >= nmax:
                break
            else:
                pass
            print("Cand %d / %d" %(ii+1, len(spsort)))
            # Give index from orig file
            #num_str = str(xx[ii]).zfill(nz)  
            num_str = str(cc.cnum).zfill(nz)  
            outfile = "%s_cand%s" %(outbase, num_str)

            # start position in bytes
            nstart = (cc.samp - nspec//2) * nchans * bps + hdr_bytes
            # If cand is at very start of file, make sure we 
            # we don't back up into header
            if nstart < hdr_bytes:
                nstart = hdr_bytes
              
            #print(nstart)
            fin.seek(nstart)
            #print(fin.tell())

            # get start time 
            tstart = cc.samp * dt
            
            # num to read in smaples b/c fromfile
            nread  = nspec * nchans  
            
            # Get the data
            dat = np.fromfile(fin, count=nread, dtype='f4')
            dat = np.reshape(dat, (-1, nchans))
            #print(dat.shape)
    
            # Write to file
            if len(dat):
                wfil.write_snippet(dat, yr, outfile, tstart)
            else:
                print("No data found!")

    return 


def parse_input():
    """
    Use argparse to parse input
    """
    prog_desc = "Single pulse search and cand extraction"
    parser = ArgumentParser(description=prog_desc)

    parser.add_argument('infile', help='Input file name')
    parser.add_argument('outdir', help='Output data directory')
    parser.add_argument('-dm', '--dm', required=True, type=float,
           help='Dispersion Measure (pc/cc)')
    parser.add_argument('-snr', '--snrmin', required=False, 
           help='SNR threshold for single pulse search (def: 6)',
           type=float, default=6.0)
    parser.add_argument('-w', '--width', default=0.1,
           help='Output size for candidate snippets in sec (def: 0.1)',
           required=False, type=float)
    parser.add_argument('-zap', '--zap', required=False, default="", 
           help='Comma separated list of frequency channels ')
    parser.add_argument('-ezap', '--edgezap', required=False, 
           help='Subband edge channels to zap (def: 2)',
           type=int, default=2)
    parser.add_argument('-nsub', '--nsub', required=False, 
           help='Number of subbands in data (def: 1)',
           type=int, default=1)
    parser.add_argument('-ncs', '--nchansub', required=True, 
           help='Number of channels per subband',
           type=int, default=1024)
    parser.add_argument('-f', '--filter',
           help='Filter frequency - comma separated list giving the '+\
                'center frequency, the number of harmonics beyond ' +\
                'fundamental, and width in Hz (e.g., \'60.0,5,1.0\' to '+\
                '60Hz signal and 5 harmonics with width 1.0Hz).  To zap '+\
                'multiple frequencies, repeat this argument',
           action='append', required=False, default=[])
    parser.add_argument('-mw', '--maxwidth', required=False, 
           help='Max boxcar width in ms for SP search (def: 10)',
           type=float, default=10.0)
    parser.add_argument('-mc', '--maxcands', required=False, 
           help='Max cands for plotting. ' +\
                'Do not plot if exceeds this number. (def = -1, no lim)',
           type=int, default=-1)
    parser.add_argument('--zerodm', action='store_true', 
           help='Apply the zero DM filter when dedispersing')
    parser.add_argument('--badblocks', action='store_true', 
           help='Ignore bad blocks in single pulse search')
    parser.add_argument('-tel', '--tel', default=None, required=False, 
           help='DSN Telescope Name GS/RO/CN (def: None)')

    args = parser.parse_args()

    return args


def main():
    """
    Run processing
    """
    tstart = time.time()

    # Parse input
    args = parse_input()

    print("\n\n===== PARAMETERS =====")
    filfile = args.infile
    print("  Filterbank File: %s" %filfile)
    outdir = args.outdir
    print("  Output Directory: %s" %outdir)
    dm = args.dm
    print("  Dispersion Measure: %.2f pc/cc" %dm)
    zdm = args.zerodm
    print("  Zero DM during de-dispersion: %r" %zdm)
    snr = args.snrmin
    print("  Candidate SNR Threshold: %.1f" %snr)
    blocks = args.badblocks
    print("  Ignore bad blocks in SP search: %r" %blocks)
    width = args.width
    print("  Candidate Snippet Size: %.3f sec" %width)
    zap = args.zap
    print("  List of Channels to zap: %s" %zap)
    ezap = args.edgezap
    print("  Edge Channels to zap: %d" %ezap)
    nsub = args.nsub
    print("  Number of subbands: %d" %nsub)
    nchansub = args.nchansub
    print("  Number of channels per subband: %d" %nchansub)
    filter_list = args.filter
    if len(filter_list):
        fzap_str = ';'.join(filter_list)
    else:
        fzap_str = "No filtering"
    print("  Filtering (f0, nh, W): %s" %fzap_str)
    mw = args.maxwidth
    print("  Max single pulse template width: %.1fms" %mw)
    tel = args.tel
    print("  Telescope: %s" %tel)
    max_cands = args.maxcands
    print("  Max cands for plotting: %d" %max_cands)
    print("===================\n\n")

    # Check that fil file and output dir exist
    if not os.path.exists(filfile):
        print("Filterbank file not found!")
        print("   %s" %filfile)
        return 
    # Check if output directory exists
    if not os.path.exists(outdir):
        print("Output directory does not exist!")
        print("   %s" %outdir)
        return 

    ### Run fix file ###
    if tel is not None:
        dsn = True
        fix_file(filfile, tel, dsn)
    else: 
        pass

    ### Get filterbank info ###
    nchans, fch1, foff, dt = get_chan_info(filfile)

    ### Get zap channel string ###
    zstr = get_zap_chans(ezap, nchansub, nsub, zchans=zap)

    ### Dedisperse ###
    print("\n\n===== DEDISPERSION =====")
    datfile = dedisperse(filfile, dm, zstr, zdm=zdm, outdir=outdir)

    ### Frequency Filter ###
    print("\n\n===== FILTERING =====")
    datfile = filter_dat(datfile, filter_list)

    ### Single Pulse Search ###
    print("\n\n===== SP SEARCH =====")
    mw_sec = mw * 1e-3
    mw_bins = min( int(mw_sec/dt), 8000 ) 
    spfile = sp_search(datfile, snr, bb=blocks, 
                       maxwidth=mw_sec, dtrendlen=32)

    # Read cands from SP file
    splist = cands_from_spfile(spfile)
    ncands = len(splist)
    print("  Number of candidates: %d" %ncands) 

    nplot = 10

    if (max_cands > 0) and (ncands > max_cands):
        print("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print(" Found %d single pulse candidates" %ncands)
        print(" which exceeds user defined limit of %d. " %max_cands)
        print("")
        print(" Will NOT extract all cands and make plots!")
        print(" Will plot ~10 for diagnostic purposes")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    else:
        nplot = -1

    ### Extract snippets ###
    print("\n\n===== EXTRACT CAND DATA =====")
    
    # Get output base
    filbase = (filfile.split("/")[-1]).split(".fil")[0] 
    outbase =  "%s/%s_DM%.3f" %(outdir, filbase, dm)

    # Set spec + extract data
    nspec = int( width / dt + 0.5 )
    your_extract_snippets(filfile, outbase, splist, nspec, nmax=nplot)

    ### Make plots from candidates ###
    sbase = "%s_DM%.3f" %(filbase, dm)
    # Set dec factor to get about 100 channels
    f_dec = max( int(nchans/100), 1 )
    nbins=100
    wmax=-1
    rmax=-1

    # Make summary plots
    sp_plt.make_summary_plots(spfile, outdir, sbase, 
                              outdir=outdir, rmax=rmax)

    # Make cand plots
    pzstr = get_zap_chans(ezap, nchansub, nsub, zchans=zap, flip=False)
   
    if nplot > -1:
        plt_splist = splist[: nplot] 
    else:
        plt_splist = spfile
    sp_plt.make_snippet_plots(plt_splist, outdir, sbase, outdir=outdir, 
                              snr_min=snr, wmax=wmax, t_dec=-1, 
                              f_dec=f_dec, outbins=nbins, rmax=rmax, 
                              zstr=pzstr)

    return


debug = 0

if __name__ == "__main__":
    if debug:
        pass
    else:
        main()
