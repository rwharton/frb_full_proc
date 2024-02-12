import your 
import numpy as np
import os
from argparse import ArgumentParser
from subprocess import call
import glob
#import sp_spec as splt

import matplotlib.pyplot as plt 
from matplotlib.gridspec import GridSpec

#os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

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

###########################
##  De-disperse and Avg  ##
###########################

# Lorimer & Kramer (2006)
kdm = 4148.808 # MHz^2 / (pc cm^-3)

def deltaT(freqs, f0, dt, kdm=kdm):
    """
    Return array of single unit DM delays in samples for MHz freqs
    """
    return (kdm / dt) * (freqs**-2.0 - f0**-2)


def dmdt(DM, freqs, f0, dt, kdm=kdm):
    """
    Return DM delays (in samples) for DM=DM
    """
    return DM * deltaT(freqs, f0, dt, kdm=kdm)


def dedisperse_dspec(dspec, dm, freqs, f0, dt, kdm=kdm, reverse=False):
    nchans = dspec.shape[0]
    nt = dspec.shape[1]
    dsamps = dmdt(dm, freqs, f0, dt, kdm=kdm)
    #dsamps -= np.min(dsamps)
    if reverse:
        sgn = -1.0
    else:
        sgn = +1.0

    dout = np.zeros( dspec.shape )
    for ii, dd in enumerate(dspec):
        ak = np.fft.rfft(dd)
        bfreq = np.arange(len(ak)) / (1.0 * len(dd))
        shift = np.exp(sgn * 1.0j * 2 * np.pi * bfreq * dsamps[ii])
        dd_shift = np.fft.irfft( ak * shift )
        dout[ii] = dd_shift[:]

    return dout


def dspec_avg_chan(dspec, freqs, avg_chan=1):
    Nchan = dspec.shape[1]
    n = int(Nchan / avg_chan)

    freq_out = np.zeros(n)
    dd_out = np.zeros( (dspec.shape[0], n) )

    for ii in range(n):
        sl = slice(ii * avg_chan, (ii+1) * avg_chan)
        freq_out[ii] = np.mean(freqs[sl])
        dd_out[:, ii] = np.mean(dspec[:, sl], axis=1)

    return freq_out, dd_out


def dspec_avg_time(dspec, tt, avg_tsamp=1):
    Nt = dspec.shape[0]
    n = int(Nt / avg_tsamp)

    tt_out = np.zeros(n)
    dd_out = np.zeros( (n, dspec.shape[1]) )

    for ii in range(n):
        sl = slice(ii * avg_tsamp, (ii+1) * avg_tsamp)
        tt_out[ii] = np.mean(tt[sl])
        dd_out[ii, :] = np.mean(dspec[sl,:], axis=0)

    return tt_out, dd_out


############################
##  Get data / Make Plot  ##
############################

def read_bpass(bp_file):
    """
    Get bandpass from file 
    """
    dat = np.loadtxt(bp_file)
    freqs = dat[:, 0]
    bp = dat[:, 1]
    return freqs, bp


def chans_from_zapstr(zstr):
    """
    Convert comma separated list of channels 
    to array and treat colons as ranges min:max:sep
    """
    chans = []
    for zz in zstr.split(','):
        if ":" in zz:
            zrange = zz.split(':')
            zlo = int(zrange[0])
            zhi = int(zrange[1])
            zlist = np.arange(zlo, zhi+1, 1).tolist()
        else:
            zlist = [int(zz)]
        chans += zlist
    zchans = np.unique( chans )
    return zchans


def get_snippet_data2(filfile, dm, favg=1, tavg=1, 
                      bpass=True, bpfile=None, 
                      zapstr=""):
    """
    Use your to read data and get metadata
    """
    # Get info
    nchans, fch1, foff, dt = get_chan_info(filfile)
    yr = your.Your(filfile)
    nsamps = yr.your_header.nspectra
    freqs = np.arange(nchans) * foff + fch1
    tt = np.arange(nsamps) * dt
    tt -= np.mean(tt)

    # Get data
    dat = yr.get_data(0, yr.your_header.nspectra)

    # Dedisperse
    dout = dedisperse_dspec(dat.T, dm, freqs, freqs[0], dt)
    dout = dout.T

    # Bandpass
    if bpass:
        if bpfile is None:
            Nthird = int( nsamps / 3 )
            bp = np.mean(dout[:Nthird], axis=0)
        else:
            fbp, bp = read_bpass(bpfile)
        dout = dout/bp - 1
    else: pass

    # Zap chans -- replace with median
    if len(zapstr):
        zchans = chans_from_zapstr(zapstr)
        gchans = np.setdiff1d(np.arange(dout.shape[1]), zchans)
        davg = np.mean(dout[:, gchans])
        #print(davg)
        dout[:, zchans] = davg

    # Average in time (if desired)
    if tavg > 1:
        tt_out, dout = dspec_avg_time(dout, tt, avg_tsamp=tavg)
    else:
        tt_out = tt

    # Average in freq (if desired)
    if favg > 1:
        ff_out, dout = dspec_avg_chan(dout, freqs, avg_chan=favg)
    else:
        ff_out = freqs

    return tt_out, ff_out, dout


def get_snippet_data(filfile, dm, favg=1, tavg=1, bpass=True):
    """
    Use your to read data and get metadata
    """
    # Get info
    nchans, fch1, foff, dt = get_chan_info(filfile)
    yr = your.Your(filfile)
    nsamps = yr.your_header.nspectra
    freqs = np.arange(nchans) * foff + fch1
    tt = np.arange(nsamps) * dt
    tt -= np.mean(tt)

    # Get data
    dat = yr.get_data(0, yr.your_header.nspectra)

    # Dedisperse
    dout = dedisperse_dspec(dat.T, dm, freqs, freqs[0], dt)
    dout = dout.T

    if bpass:
        Nthird = int( nsamps / 3 )
        bp = np.mean(dout[:Nthird], axis=0)
        dout = dout/bp - 1
    else: pass

    # Average in time (if desired)
    if tavg > 1:
        tt_out, dout = dspec_avg_time(dout, tt, avg_tsamp=tavg)
    else:
        tt_out = tt

    # Average in freq (if desired)
    if favg > 1:
        ff_out, dout = dspec_avg_chan(dout, freqs, avg_chan=favg)
    else:
        ff_out = freqs

    return tt_out, ff_out, dout


def make_plot(filfile, dm, favg=1, tavg=1, spec_sig=5,
              outbins=128, outfile=None, cnum=None, 
              ctime=None, cdm=None, cwidth=None, 
              csnr=None, bpass=True, bpfile=None, zstr=""):
    """
    make 3 panel plot
    """
    if outfile is not None:
        plt.ioff()
    else: pass

    #tt, freqs, dat = get_snippet_data(filfile, dm,
    #                         favg=favg, tavg=tavg, bpass=True)
    tt, freqs, dat = get_snippet_data2(filfile, dm,
                             favg=favg, tavg=tavg, bpass=bpass,
                             bpfile=bpfile, zapstr=zstr)
    tt *= 1e3
    #tt0, _, dat0 = get_snippet_data(filfile, 0,
    #                         favg=favg, tavg=tavg, bpass=True)
    #tt0 *= 1e3

    # Set tlim if desired
    if outbins > 0:
        Nt = len(tt)
        xt_lo = max( 0, Nt//2 - outbins//2)
        xt_hi = min( Nt-1, Nt//2 + outbins//2)
        tlim = (tt[xt_lo], tt[xt_hi])
    else:
        tlim = (tt[0], tt[-1])

    fig = plt.figure(constrained_layout=True)
    gs = GridSpec(3, 3, figure=fig)

    ax_t  = fig.add_subplot(gs[0, 0:2])
    ax_ds = fig.add_subplot(gs[1:, 0:2])
    ax_f  = fig.add_subplot(gs[1:, 2])
    ax_txt = fig.add_subplot(gs[0, 2])

    # Dynamic Spectrum
    ext = [tt[0], tt[-1], freqs[0], freqs[-1]]
    d_sig = np.std(dat)
    d_med = np.median(dat)
    vmax = d_med + 4 * d_sig
    vmin = d_med - 3 * d_sig
    ax_ds.imshow(dat.T, aspect='auto', interpolation='nearest',
                 origin='lower', extent=ext, vmin=vmin, vmax=vmax)
    ax_ds.set_xlabel("Time (ms)", fontsize=14)
    ax_ds.set_ylabel("Freq (MHz)", fontsize=14)
    ax_ds.set_xlim(tlim)

    # Time series
    ts = np.mean(dat, axis=1)
    Nthird = len(ts) // 3
    avg = np.mean(ts[:Nthird])
    sig = np.std(ts[:Nthird])
    ts = (ts-avg)/sig

    #ts0 = np.mean(dat0, axis=1)
    #ts0 = (ts0-avg)/sig
    #ax_t.plot(tt0, ts0, c='0.5', zorder=-1)

    ax_t.plot(tt, ts)
    ax_t.tick_params(axis='x', labelbottom=False)
    tylim = ax_t.get_ylim()
    #ax_t.set_xlim(tt[0], tt[-1])
    ax_t.set_xlim(tlim)

    # Spectrum
    xpk = np.argmax(ts)
    #xx_below = np.where( ts <= 0.1*np.max(ts) )[0]
    #xx_lo = np.max(xx_below[xx_below <= xpk ])
    #xx_hi = np.min(xx_below[xx_below >= xpk ])
    xx_lo = int( len(ts) // 2 ) - 2
    xx_hi = int( len(ts) // 2 ) + 2

    off_spec = np.mean(dat[:Nthird], axis=0)
    off_sig = np.std(off_spec)
    off_sig = off_sig * np.sqrt(Nthird/(xx_hi-xx_lo))

    spec = np.mean(dat[xx_lo:xx_hi], axis=0) / off_sig
    ax_f.plot(spec, freqs)
    ax_f.set_ylim(freqs[0], freqs[-1])
    ax_f.tick_params(axis='y', labelleft=False)

    # Shade region used for spec
    tlo = tt[xx_lo]
    thi = tt[xx_hi]
    ax_t.fill_betweenx([-10, 100], tlo, thi, color='r', alpha=0.1)
    #print(tlo, thi)
    ax_t.set_ylim(tylim)

    # Add cand info subplot
    ax_txt.axis('off')
    outstr = ""
    if cnum is not None:
        outstr += "Cand: %d\n" %cnum
    if ctime is not None:
        outstr += "Time: %.3f s\n" %ctime
    if cdm is not None:
        outstr += "DM: %.2f\n" %cdm
    if cwidth is not None:
        outstr += "Width: %d bins\n" %cwidth
    if csnr is not None:
        outstr += "SNR: %.1f\n" %csnr
    ax_txt.text(0.00, 0.9, outstr, fontsize=12, ha='left', va='top',
                transform=ax_txt.transAxes)

    if outfile is not None:
        plt.savefig(outfile, dpi=100, bbox_inches='tight')
        plt.ion()
    else:
        plt.show()

    return

#########################
##  Remove Duplicates  ##
#########################

def get_chan_info(data_file):
    yr = your.Your(data_file)
    foff = yr.your_header.foff
    fch1 = yr.your_header.fch1
    dt   = yr.your_header.tsamp
    nchans = yr.your_header.nchans
    
    return nchans, fch1, foff, dt


def t_dm_shifts(dDMs, nchan, fch1, df):
    """
    calc time offset from incorrect dm
    """
    fchans = np.arange(nchan) * df + fch1
    fhi = np.max(fchans)
    
    avg_f2 = np.mean( fchans**-2.0 - fhi**-2.0) 

    dts = 4.15e3 * avg_f2 * dDMs

    return dts


def sift_dm_dupes(clist, nchan, fch1, df, dt):
    """
    sift through candidate list and remove 
    candidates that may be duplicates of 
    higher SNR burts at wrong DMs
    """ 
    snrs = np.array([ cc.snr for cc in clist ])
    xx = np.argsort(snrs)[::-1]

    c_snrs = clist[xx]
    c_out  = []

    while(len(c_snrs)):
        #print(c_snrs[0].snr)
        # add current candidate to output
        c_out.append(c_snrs[0])

        # get dms and times
        dms = np.array([ cc.dm for cc in c_snrs ])
        tts = np.array([ cc.time for cc in c_snrs ])
        wws = np.array([ cc.wbins for cc in c_snrs ]) * dt

        # get dm and time of current candidate 
        dm0 = c_snrs[0].dm
        tt0 = c_snrs[0].time
        ww0 = c_snrs[0].wbins * dt

        # Calc offsets
        dts = t_dm_shifts(dms-dm0, nchan, fch1, df)

        # what cands are within dts?
        cond_xx = (np.abs(tt0-tts) <= np.abs(dts) + (wws + ww0)) & \
                  (np.sign(tt0-tts) * np.sign(dts) >= 0) 

        # cands outside are then just the negation
        yy = np.where( cond_xx == False )[0]

        # Now update the c_snrs array to just keep the 
        # ones that are NOT duplicates 
        c_snrs = c_snrs[yy]

        # print the length as a check
        #print(len(c_snrs))

    # Convert output to array 
    c_out = np.array(c_out)

    return c_out
          


################################
## Cross Matching Cand Lists  ##
################################

def cross_check_clists(clist1, clist2, nchans, fch1, df, dt):
    """
    Check for cross-matches for clist1 in clist2
    """
    # Time Samples
    ts1 = np.array([ cc.samp for cc in clist1 ]) * dt 
    ts2 = np.array([ cc.samp for cc in clist2 ]) * dt 

    # Widths 
    w1 = np.array([ cc.wbins for cc in clist1 ]) * dt 
    w2 = np.array([ cc.wbins for cc in clist2 ]) * dt

    # dms 
    dm1 = np.array([ cc.dm for cc in clist1 ]) 
    dm2 = np.array([ cc.dm for cc in clist2 ]) 

    # match lists 
    midx1 = []
    midx2 = []

    for ii in range(len(clist1)):
        tdiff = ts1[ii] - ts2
        jj = np.argmin( np.abs( tdiff ) )
        td_ij = tdiff[jj]
        dm_dt = t_dm_shifts(dm2[jj]-dm1[ii], nchans, fch1, df)
        cond1 = (np.abs(td_ij) <= np.abs(dm_dt) + (w2[jj] + w1[ii])) & \
                (np.sign(td_ij) * np.sign(dm_dt) >= 0)
        cond2 = np.abs(td_ij) <= (w2[jj] + w1[ii])

        cond = cond1 | cond2

        print(ii, td_ij * 1e3, w2[jj] + w1[ii], dm_dt * 1e3, cond)
        if cond:
            midx1.append(ii)
            midx2.append(jj)
        else:
            pass

    midx1 = np.array(midx1)
    midx2 = np.array(midx2)

    return midx1, midx2
             
        
def cross_check_clists_orig(clist1, clist2):
    """
    Check for cross-matches for clist1 in clist2
    """
    # Samples 
    ts1 = np.array([ cc.samp for cc in clist1 ]) 
    ts2 = np.array([ cc.samp for cc in clist2 ]) 

    # Widths 
    w1 = np.array([ cc.wbins for cc in clist1 ])  
    w2 = np.array([ cc.wbins for cc in clist2 ])  

    # match lists 
    midx1 = []
    midx2 = []

    for ii in range(len(clist1)):
        tdiff = ts1[ii] - ts2
        jj = np.argmin( np.abs( tdiff ) )
        td_ij = np.abs(tdiff[jj])
        if td_ij <= (w1[ii] + w2[jj])//2:
            midx1.append(ii)
            midx2.append(jj)
        else:
            pass

    midx1 = np.array(midx1)
    midx2 = np.array(midx2)

    return midx1, midx2


######################
##  Filter by Rate  ##
######################

def filter_by_rate(clist, rate_max):
    """
    Filter out time periods that have counts per 
    minute above rate_max
    """
    # get times
    tt = np.array([ cc.time for cc in clist ])
    
    # bin to histogram
    bb = np.arange(0, np.max(tt), 60)
    nn, tt_e = np.histogram(tt, bins=bb)

    # find where rate exceeded (if any)
    xxr = np.where( nn >= rate_max )[0]

    if len(xxr) > 0:
        yy_bad = np.array([])
        for xxi in xxr:
            tlo = tt_e[xxi]
            thi = tt_e[xxi+1]
            
            yy = np.where( (tt >= tlo) & (tt <= thi) )[0]
            yy_bad = np.hstack( (yy_bad, yy) )

        yy_good = np.setdiff1d( np.arange(len(tt)), yy_bad )
        yy_good = np.unique(yy_good)

        clist_out = clist[yy_good]

    else:
        clist_out = clist[:]

    return clist_out


########################
##  Candidate Params  ##
########################

def cand_params(candlist):
    """
    return arrays of cand pars
    """
    tts = np.array([ cc.time for cc in candlist ])
    dms = np.array([ cc.dm for cc in candlist ])
    wws = np.array([ cc.wbins for cc in candlist ])
    snrs = np.array([ cc.snr for cc in candlist ])

    return tts, snrs, dms, wws
         

def get_csift(clist, data_file):
    """
    Run extra dm sifting to clist
    """
    yr = your.Your(data_file)
    foff = yr.your_header.foff
    fch1 = yr.your_header.fch1  
    dt   = yr.your_header.tsamp
    nchans = yr.your_header.nchans
    csift = sift_dm_dupes(clist, nchans, fch1, foff, dt)
    return csift
    

def get_fil_dict(fildir, basename):
    """
    Get dictionary of fil paths
    """
    # Get fil paths
    filfiles = glob.glob('%s/%s*.fil' %(fildir, basename))
    
    # get nums 
    nums = [ int(fn.rsplit("_cand")[-1].rsplit(".fil")[0]) for fn in filfiles ]
    
    # make dictionary
    fn_dict = { nums[ii] : filfiles[ii] for ii in range(len(filfiles)) }
    
    return fn_dict

    
def make_snippet_plots(sp, fildir, basename, outdir='.', 
                       snr_min=7, wmax=-1, t_dec=-1, f_dec=1, 
                       outbins=128, rmax=-1,
                       bpass=True, bpfile=None, zstr=""):
    """
    Make plots from snippet filterbanks

    t_dec == - 1 means use width//2

    Assumes cand in middle of file
    """
    # Assuming string is file name  
    if type(sp) == str:
        clist_all = cands_from_spfile(sp)

    # If list or array assume it is clist
    elif type(sp) in [np.ndarray, list]:
        clist_all = [ cc for cc in sp ] 

    # Otherwise we dont know
    else: 
        print("Unknown sp type: must be file name or array/list")
        return

    # Make new array with snr_min cut
    clist = np.array([ cc for cc in clist_all if cc.snr > snr_min ])

    # If rmax > 0, make rate cut off
    if rmax > 0:
        clist = filter_by_rate(clist, rmax)
    
    # Make new array with wmax cut
    if wmax <= 0:
        pass
    else:
        clist = np.array([ cc for cc in clist if cc.wbins <= wmax ])

    # Print number of cands
    print("Processing %d cands" %len(clist))

    # Get fil file dict
    fn_dict = get_fil_dict(fildir, basename)

    # now convert to your format and save
    Nc = len(clist)
    for ii, cc in enumerate(clist):
        print("Cand %d/%d" %(ii+1, Nc))
    
        cnum = cc.cnum
        #print(cnum)
        infile = fn_dict.get(cnum)
        
        if infile is None:
            print("File not found: %s" %infile)
            continue
        else: 
            print("Found file: %s" %infile)

        if t_dec == -1:
            t_dec_fac = max(1, cc.wbins // 2)
        else:
            t_dec_fac = t_dec

        outfile = "%s/cand%04d.png" %(outdir, cnum)

        make_plot(infile, cc.dm, favg=f_dec, tavg=t_dec_fac, 
                   cnum=cnum, ctime=cc.time, cdm=cc.dm, 
                   cwidth=cc.wbins, csnr=cc.snr, outbins=outbins, 
                   outfile=outfile, bpass=bpass, bpfile=bpfile, 
                   zstr=zstr)
        
    return 


def write_cands(clist, outfile):
    """
    Write out candidates
    """
    with open(outfile, 'w') as fout:
        for cc in clist:
            ostr = "%.1f   %.1f   %.3f    %d   %d\n" %(\
                    cc.dm, cc.snr, cc.time,  cc.samp, cc.wbins)
            fout.write(ostr)

    return


#######################
###  SUMMARY PLOTS  ###
#######################

def time_summary(csift, clist, title=None, outfile=None, rmax=-1):
    """
    Summary plots along time axis
    """
    if outfile is not None:
        plt.ioff()
    else: pass

    fig = plt.figure(figsize=(14, 10))
    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2, sharex=ax1)
    ax3 = fig.add_subplot(3, 1, 3, sharex=ax1)

    tt0, snr0, dm0, ww0 = cand_params(clist)
    tt1, snr1, dm1, ww1 = cand_params(csift)

    # Top Panel is DM vs time with SNR size
    ax1.scatter(tt0, dm0, s=25*(snr0/10)**2.0, fc='none', 
                ec='k', lw=1)
    ax1.scatter(tt1, dm1, s=25*(snr1/10)**2.0, ec='r', 
                marker='s', lw=1, fc='none')

    ax1.set_yscale('log')
    ax1.set_ylabel('DM (pc/cc)', fontsize=14)

    # Middle panel is count rate per minute
    bb = np.arange(0, np.max(tt0), 60)
    ax2.hist(tt0, bins=bb, color='k')
    ax2.hist(tt1, bins=bb, color='r')
    if rmax > -1:
        ax2.axhline(y=rmax, c='k', ls='--')

    ax2.set_yscale('log')
    ax2.set_ylabel('Hits/min', fontsize=14)

    # Bottom plot is snr
    ax3.plot(tt0, snr0, marker='.', ls='', c='k')
    ax3.plot(tt1, snr1, marker='o', mfc='none', mec='r', 
             ls='')
    ax3.set_ylim(5)

    ax3.set_ylabel('SNR', fontsize=14)
    ax3.set_xlabel("Time (s)", fontsize=14)
    
    if title is not None:
        fig.suptitle(title, fontsize=16, y=0.925)

    if outfile is not None:
        plt.savefig(outfile, dpi=100, bbox_inches='tight')
        plt.close()
        plt.ion()
    else:
        plt.show()
    
    return


def snr_summary(csift, clist, title=None, outfile=None):
    """
    Summary plots along time axis
    """
    if outfile is not None:
        plt.ioff()
    else: pass
    
    fig = plt.figure(figsize=(12, 8))
    ax1 = fig.add_subplot(2, 3, 1)
    ax2 = fig.add_subplot(2, 3, 2)
    ax3 = fig.add_subplot(2, 3, 3)
    ax4 = fig.add_subplot(2, 3, 4)
    ax5 = fig.add_subplot(2, 3, 5)
    ax6 = fig.add_subplot(2, 3, 6)

    tt0, snr0, dm0, ww0 = cand_params(clist)
    tt1, snr1, dm1, ww1 = cand_params(csift)

    fs = 10

    # Top left plot is SNR Hist
    ax1.hist(snr1, log=True)
    ax1.set_xlabel("SNR", fontsize=fs)

    # Top middle plot is DM hist
    dbins = 10**np.linspace(2, 6, 20)
    ax2.hist(dm1, bins=dbins)
    ax2.set_xscale('log') 
    ax2.set_xlabel("DM (pc/cc)", fontsize=fs)

    # Top right plot is width (bins)
    wbins = 2**np.arange(10) / np.sqrt(2)
    ax3.hist(ww1, bins=wbins)
    ax3.set_xscale('log', base=2)
    ax3.set_xlabel("Width (bins)", fontsize=fs)
    

    # Bottom Left plot DM vs Width
    ax4.plot(dm1, ww1, marker='o', ls='', mfc='none')
    ax4.set_xscale('log', base=10)
    ax4.set_yscale('log', base=2)
    ax4.set_xlabel('DM (pc/cc)', fontsize=fs)
    ax4.set_ylabel('Width (bins)', fontsize=fs)

    # Bottom Middle plot DM vs SNR
    ax5.plot(dm1, snr1, marker='o', ls='', mfc='none')
    ax5.set_xscale('log')
    ax5.set_xlabel('DM (pc/cc)', fontsize=fs)
    ax5.set_ylabel('SNR', fontsize=fs)

    # Bottom Right plot width vs SNR
    ax6.plot(ww1, snr1, marker='o', ls='', mfc='none')
    ax6.set_xscale('log', base=2)
    ax6.set_xlabel('Width (bins)', fontsize=fs)
    ax6.set_ylabel('SNR', fontsize=fs)

    if title is not None:
        fig.suptitle(title, fontsize=16, y=0.925)
    
    if outfile is not None:
        plt.savefig(outfile, dpi=100, bbox_inches='tight')
        plt.close()
        plt.ion()
    else:
        plt.show()
    
    return


def make_summary_plots(spfile, fildir, basename, outdir='.', rmax=-1):
    """
    make and save the two summary plots

    basename will be the title of each plot 
    and the base of the outfile name
    """
    # get one fil file for params for sifting
    fil_files = glob.glob("%s/%s*fil" %(fildir, basename))
    if len(fil_files) == 0:
        print("No snippet *.fil files found!")
        return 
    else: pass

    data_file = fil_files[0]

    clist = cands_from_spfile(spfile)
    csift = get_csift(clist, data_file)
    
    # Time summary
    #time_file = "%s_time.png" %basename
    time_file = "%s/cand_time.png" %outdir
    time_summary(csift, clist, title=basename, 
                 outfile=time_file, rmax=rmax)

    # SNR summary 
    #snr_file = "%s_snr.png" %basename
    snr_file = "%s/cand_snr.png" %outdir
    snr_summary(csift, clist, title=basename, 
                outfile=snr_file)

    return


def parse_input():
    """
    Use 'your' to make pulse candidates + plots from snippet filterbanks
    """
    prog_desc = "Make plots of snippet filterbanks"
    parser = ArgumentParser(description=prog_desc)

    parser.add_argument('infile', help='Input *.singlepulse file')
    parser.add_argument('fildir', help='Filterbank file directory')
    parser.add_argument('-b', '--basename',
                        help='Basename for filfiles.  Will search '+\
                             'on basename*fil (default = all in fildir)',
                        required=False, default="")
    parser.add_argument('-o', '--outdir',
                        help='Output directory for plots (default = .)', 
                        required=False, default=".")
    parser.add_argument('-s', '--snr',
                        help='SNR threshold for plotting (default = all)',
                        required=False, type=float, default=-1)
    parser.add_argument('-w', '--wmax',
                        help='Maximum width to consider (default = -1 ie no max)',
                        required=False, type=int, default=-1)
    parser.add_argument('-t', '--tdec',
                        help='Time decimation factor (default = -1, auto)',
                        required=False, type=int, default=-1)
    parser.add_argument('-f', '--fdec',
                        help='Frequency decimation factor (default = 1)',
                        required=False, type=int, default=1)
    parser.add_argument('-n', '--nbins',
                        help='Number of output bins (default = 128)',
                        required=False, type=int, default=128)
    parser.add_argument('-r', '--rmax',
                        help='Max number of cands per minute '+\
                             '(default = -1, no limit)',
                        required=False, type=float, default=-1)
    parser.add_argument('-z', '--zapchans',
                        help='Comma separated list of channels to zap',
                        required=False, default="")
    parser.add_argument('--no_cands', action='store_true',
                        help='Do NOT make candidate plots (default false)',
                        default=False)
    parser.add_argument('--no_sum', action='store_true',
                        help='Do NOT make summary plots (Default false)',
                        default=False)

    args = parser.parse_args()
    
    infile = args.infile
    fildir = args.fildir
    outdir = args.outdir
    bname = args.basename
    tdec = args.tdec
    fdec = args.fdec
    snr = args.snr
    wmax = args.wmax
    rmax = args.rmax
    nbins = args.nbins
    zstr = args.zapchans
    skip_cands = args.no_cands
    skip_sum = args.no_sum
    
    return infile, fildir, outdir, bname, tdec, fdec, snr, \
           wmax, nbins, zstr, rmax, skip_cands, skip_sum

debug = 0

if __name__ == "__main__":
    if debug:
        pass

    else:
        infile, fildir, outdir, bname, tdec, fdec, snr, wmax, nbins, zstr, rmax, skip_cands, skip_sum = parse_input()

        if not skip_sum:
            make_summary_plots(infile, fildir, bname, 
                               outdir=outdir, rmax=rmax)

        if not skip_cands:
            make_snippet_plots(infile, fildir, bname, outdir=outdir, 
                               snr_min=snr, wmax=wmax, t_dec=tdec, f_dec=fdec, 
                               outbins=nbins, rmax=rmax, zstr=zstr)
