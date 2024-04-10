import numpy as np
import os
import sys
import glob
import time
from subprocess import call
from argparse import ArgumentParser
import bp_rfi as bp_rfi

#srcdir = '/src/bb_proc'
cur_dir = os.path.realpath(__file__)
srcdir  = cur_dir.rsplit('/', 1)[0]

def convert_cs2fil(csdir, bname, outdir, nchan, dm, 
                   nthread, memlim):
    """
    Run bb2fil_chunk.py to convert cs fil to fil

    This will make a DADA file from the cs file, 
    then run digifil to make a channelized filterbank 
    file with intrachannel dispersive delays removed.
    """
    tstart = time.time()

    script_path = "%s/bb2fil/bb2fil_chunk.py" %srcdir
    cmd = "python -u " +\
          "%s " %script_path +\
          "--basename %s " %bname +\
          "--cs_dir %s " %csdir +\
          "--dada_dir %s " %outdir +\
          "--fil_dir %s " %outdir +\
          "--nchan %d " %nchan +\
          "--dm %.4f " %dm +\
          "--nthread %d " %nthread +\
          "--mem_lim %.1f " %memlim

    print(cmd)
    call(cmd, shell=True)

    tstop = time.time()
    tdur = tstop - tstart

    return tdur


def check_rfi_tavg(nchan, tfac):
    """
    digifil will complain if nchan * tfac > 2^19
    
    Need to track down why that is and if you can 
    change it.  In the meantime, we have this. 
    """
    if nchan * tfac >= 2**19:
        tfac_max = max( int( 2**19 / nchan ), 1 )
        if tfac_max > 10:
            # Make it somewhat round 
            tfac_out = (tfac_max // 10) * 10
        else:
            print("\n\n")
            print("  Decimation only by factor of %d" %(tfac_max))
            print("  Probably not worth it!!")
            tfac_out = tfac_max
            if tfac_max < 5:
                print("\n")
                print("  Since tfac < 5, it's definitely not worth it")
                print("  Setting tfac to -1 (ie, skipping averaging)")
                tfac_out = -1
        print("\n\n")
        print("  Reducing tfact from %d to %d b/c of digifil" %(tfac, tfac_out))
        print("  limitation that nchan * tfac < 2**19")
        print("\n\n")
    else: tfac_out = tfac

    return tfac_out


def make_rfi_fil(filfile, outdir, tfac=512, nthread=1):
    """
    Use digifil to make a highly decimated version 
    of filfile that can be used for RFI flagging 
    """
    tstart = time.time()

    # Get filfile base
    filfn = filfile.split('/')[-1]
    filbase = filfn.split('.fil')[0]

    # Make output filename
    outfn = "%s_avg%d.fil" %(filbase, tfac)
    
    # Full output path
    outfile = "%s/%s" %(outdir, outfn)

    # Set up command
    dec_cmd = "digifil -I0 -b-32 -t %d " %tfac +\
              "-threads %d " %nthread +\
              "-o %s %s " %(outfile, filfile)

    print(dec_cmd)

    # If file does not exist, run command
    if os.path.exists(outfile):
        print("  Decimated file already exists!")
        print("  %s" %outfile)
    else:
        call(dec_cmd, shell=True)

    tstop = time.time()
    tdur = tstop - tstart

    return tdur, outfile


def run_search(filfile, outdir, nchan, nsub, dm, snr, mw, max_cands, 
               width=0.1, tel='RO', edgezap=2, zap_str="", 
               fzaps=[], avoid_badblocks=False, apply_zerodm=False):
    """
    Run bbsearch.py to search for cands

    This will de-disperse, run single pulse search, 
    and make snippets 
    """
    tstart = time.time()

    script_path = "%s/bbsearch/bbsearch.py" %srcdir 

    if avoid_badblocks:
        bb_str = "--badblocks "
    else:
        bb_str = ""

    if apply_zerodm:
        zdm_str = "--zerodm "
    else:
        zdm_str = ""

    if len(fzaps):
        fz_str = ""
        for fz in fzaps:
            fz_str += "-f %s " %fz
    else:
        fz_str = ""
    
    cmd = "python -u %s " %script_path +\
          "%s" %bb_str +\
          "%s" %zdm_str +\
          "-dm %.4f " %dm +\
          "-snr %.1f " %snr +\
          "-w %.4f " %width +\
          "-mw %.3f " %mw +\
          "-ezap %d " %edgezap +\
          "-zap %s " %zap_str +\
          "-nsub %d " %nsub +\
          "-ncs %d " %(int(nchan/nsub)) +\
          "%s " %fz_str +\
          "-tel %s " %tel +\
          "-mc %d " %max_cands +\
          "%s %s " %(filfile, outdir)
   
    # Run command 
    print(cmd)   
    call(cmd, shell=True)

    tstop = time.time()
    tdur = tstop - tstart 

    return tdur


def parse_input():
    """
    Use argparse to parse input
    """
    prog_desc = "Pipeline to process and search baseband data"
    parser = ArgumentParser(description=prog_desc)

    parser.add_argument('csdir', 
                        help='Directory containing *.cs files')
    parser.add_argument('basename', 
                        help='Base name of cs file: {basename}*.cs')
    parser.add_argument('outdir', help='Output data directory')
    parser.add_argument('-dm', '--dm', 
                        help='DM for intra-channel coherent dedispersion',
                        required=True, type=float)
    parser.add_argument('-nc', '--nchan', 
                       help='Total number of output channels in filterbank',
                       required=True, type=int)
    parser.add_argument('-m', '--memlim', default=16.0,  
                        help='Max memory to use during DADA conversion in GB (def: 16)',
                        required=False, type=float)
    parser.add_argument('-nt', '--nthread', default=1, 
                        help='Number of threads for digifil processing (def: 1)',
                        required=False, type=int)
    parser.add_argument('-rt', '--rfidec', default=512, 
                        help='Number of samples to decimate filfile ' +\
                             'for RFI diagnostics (def: 512, to skip this: -1)',
                        required=False, type=int)
    parser.add_argument('-snr', '--snrmin', default=6.0,  
                        help='Minimum SNR for single pulse search (def: 6.0)',
                        required=False, type=float)
    parser.add_argument('-ezap', '--edgezap', required=False,
                        help='Subband edge channels to zap (def: 2)',
                        type=int, default=2)
    parser.add_argument('-nsub', '--nsub', required=False,
                        help='Number of subbands in data (def: 1)',
                        type=int, default=1)
    parser.add_argument('-w', '--width', default=0.1,  
                        help='Output size for candidate snippets in sec (def: 0.1)',
                        required=False, type=float)
    parser.add_argument('-mw', '--maxwidth', required=False,
                        help='Max boxcar width in ms for SP search (def: 10)',
                        type=float, default=10.0)
    parser.add_argument('-mc', '--maxcands', required=False,
                        help='Max cands for plotting. ' +\
                        'Do not plot if exceeds this number. (def = -1, no lim)',
                        type=int, default=-1)
    parser.add_argument('-f', '--filter',
           help='Filter frequency - comma separated list giving the '+\
                'center frequency, the number of harmonics beyond ' +\
                'fundamental, and width in Hz (e.g., \'60.0,5,1.0\' to '+\
                '60Hz signal and 5 harmonics with width 1.0Hz).  To zap '+\
                'multiple frequencies, repeat this argument',
           action='append', required=False, default=[])
    parser.add_argument('--zerodm', action='store_true',
                        help='Apply the zero DM filter when dedispersing')
    parser.add_argument('--badblocks', action='store_true',
                        help='Ignore bad blocks in single pulse search')
    parser.add_argument('-tel', '--tel', default='RO',  
                        help='DSN Telescope Name GS/RO/CN (def: RO)',
                        required=False)

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
    csdir = args.csdir
    print("  Baseband file directory: %s" %csdir)
    bname = args.basename
    print("  Input basename: %s" %bname)
    outdir = args.outdir
    print("  Output Directory: %s" %outdir)
    dm = args.dm
    print("  Dispersion Measure: %.2f pc/cc" %dm)
    nchan = args.nchan
    print("  Output Frequency Channels: %d" %nchan)
    memlim = args.memlim
    print("  Memory Limit: %.1f GB" %memlim)
    nthread = args.nthread
    print("  Max Threads: %d" %nthread)
    rfi_tdec = args.rfidec
    print("  Input RFI decimation factor: %d" %rfi_tdec)
    snr = args.snrmin
    print("  Candidate SNR Threshold: %.1f" %snr)
    width = args.width
    print("  Candidate Snippet Size: %.3f sec" %width)
    mw = args.maxwidth
    print("  Max single pulse template width: %.1fms" %mw)
    max_cands = args.maxcands
    print("  Max cands for plotting: %d" %max_cands)
    ezap = args.edgezap
    print("  Edge Channels to zap: %d" %ezap)
    nsub = args.nsub
    print("  Number of subbands: %d" %nsub)
    filter_list = args.filter
    if len(filter_list):
        fzap_str = ';'.join(filter_list)
    else:
        fzap_str = "No filtering"
    print("  Filtering (f0, nh, W): %s" %fzap_str)
    zdm = args.zerodm
    print("  Zero DM during de-dispersion: %r" %zdm)
    blocks = args.badblocks
    print("  Ignore bad blocks in SP search: %r" %blocks)
    tel = args.tel
    print("  Telescope: %s" %tel)
    print("===================") 
    
    # Make sure output directory exists, 
    # if it doesn't, print warning and exit
    if not os.path.exists(outdir):
        print("outdir not found: %s" %outdir)
        print("Make sure outdir exists, then try again")
        sys.exit(0)
    else: pass
        
    # If filterbank does not already exists,
    # then run the baseband to filterbank 
    # conversion
    print("\n\n=== BASEBAND TO FILTERBANK ===")
    filfile = "%s/%s.fil" %(outdir, bname)
    if not os.path.exists(filfile):
        tfil = convert_cs2fil(csdir, bname, outdir, nchan, dm, 
                              nthread, memlim)
    else:
        print("  filfile exists: %s" %filfile)
        print("  Skipping filterbank conversion...")
        tfil = 0

    # If desired, make decimated filterbank for 
    # finding bad channels.  Note that make_rfi_fil 
    # will check to see if the file already exits
    # Also make a plot showing RFI 
    print("\n\n=== FINDING BAD CHANNELS ===")
    rfi_tdec = check_rfi_tavg(nchan, rfi_tdec)
    if rfi_tdec > 0:
        dec_dur, rfi_fil = make_rfi_fil(filfile, outdir, tfac=rfi_tdec, 
                                        nthread=1) 
        # Find bad channels from bandpass       
        nchan_win = bp_rfi.get_win_num(nchan, nsub, wfrac=0.2) 
        zap_str = bp_rfi.bp_bad_chans(rfi_fil, outdir, mode='avg', 
                                      diff_thresh=3, nchan_win=nchan_win)
        # Make plots showing RFI stats calculated either with 
        # 1 min of data of time/4, whichever is smallest 
        rtime = bp_rfi.get_time_chunk(rfi_fil, 60, nt_min=4)
        outbase_rfi = "%s/bp" %outdir
        bp_rfi.rfi_plot(rfi_fil, rtime, outbase_rfi, bpass=True)
        bp_rfi.rfi_plot(rfi_fil, rtime, outbase_rfi, bpass=False)
    else:
        dec_dur = 0
        zap_str = ""

    # Make sure filterbank file exists, then 
    # run search
    print("\n\n=== SINGLE PULSE SEARCH ===")
    if not os.path.exists(filfile):
        print("  filfile not found: %s" %filfile)
        return 
    else: 
        pass
     
    tsearch = run_search(filfile, outdir, nchan, nsub, dm, snr, mw, max_cands, 
                 width=width, tel=tel, edgezap=ezap, zap_str=zap_str, 
                 fzaps=filter_list, avoid_badblocks=blocks, apply_zerodm=zdm)

    tstop = time.time()
    total_time = tstop - tstart
 
    # Now print summary of times
    print("##################################################")
    print("##              TIME SUMMARY                    ##")
    print("##################################################")
    print("")
    print("Filterbank:         %.1f minutes" %(tfil/60.))
    print("Find Bad Chans:     %.1f minutes" %(dec_dur/60.))
    print("Searching:          %.1f minutes" %(tsearch/60.))
    print("")
    print("Total Time:         %.1f minutes" %(total_time/60.))

    return


debug = 0

if __name__ == "__main__":
    if debug:
        pass
    else:
        main()
