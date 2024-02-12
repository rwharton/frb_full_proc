import numpy as np
import subprocess
import time
from tabulate import tabulate
import os
import sys
import shutil
import glob
import parse_arg as pa
import click
import multiprocessing as mp
import re
from datetime import datetime, timedelta
from astropy.time import Time, TimeDelta

# set top of src directory (where this file is)
cur_dir = os.path.realpath(__file__)
srcdir  = cur_dir.rsplit('/', 1)[0]
#srcdir  = cur_dir


def get_pars(parfile):
    """
    Read in the parameter file, parse it, and 
    return a class called 'par'

    This is a way to simulate what we did before 
    of simply importing the parameter file as a 
    module.
    """
    pdict = pa.get_par_dict(parfile)
    par = pa.Parameters(pdict)
    return par

def vrad_2_cs(par, scan_params):
    """
    Conversion from vrad to cs files
    """
    # Hold dict of params
    print("\n\n---------------------------------------------")
    print("Converting from vrad to cs")
    print("---------------------------------------------")

    scanno = scan_params['scanno']
    data_amount = scan_params['dur']
    skip = scan_params['skip']
    dm = scan_params['dm']
    
    params = {"inf": par.inf_file, "scanno": scanno,  
              "vrad_dir": par.vrad_dir, "nproc": par.nproc, 
              "vrad_base": par.vrad_base, "out": par.outdir, 
              "freq": par.freq_band, "source": par.source,
              "tele": par.telescope, "amount": data_amount, 
              "skip": skip, 
              "dm": dm, "src" : srcdir}


    
    # vrad2cs will make cs_dir and vdr_dir if they do not 
    # already exist.

    # NOTE: 'inf_file' needs to be in the output directory
    
    cmd = "python3 -u %(src)s/vrad2cs/vrad2cs.py" \
           " --inf_file %(out)s/%(inf)s" \
           " --scanno %(scanno)i" \
           " --vrad_dir %(vrad_dir)s" \
           " --vrad_base %(vrad_base)s --out_dir %(out)s" \
           " --freq_band %(freq)s --nproc %(nproc)d" \
           " --dm %(dm)f --source %(source)s --telescope %(tele)s" \
           " --skip %(skip)i" \
           " --data_amount %(amount)i" % params
        
    print(cmd)
    subprocess.run(cmd, shell=True)

    return

def read_scan_table(scan_file,scanno):
    """
    Given scan_table, find MJD time and return it 
    in YY/Day_HH:MM:SS format
    """
    # Define the output format
    output_format = "%y/%j_%H:%M:%S"
    # Open the scan table file
    with open(scan_file, 'r') as file:
        lines = file.readlines()
        nlines = len(lines)
        print("wam>>> nlines in scan table", nlines)

    line = lines[scanno-1]
    #line = file.readline().rstrip()

    # Split the line into parts
    parts = line.split()

    src_name = parts[0]

    # Extract the date and time parts
    date = parts[3] # could be different
    time = parts[4] + parts[5] + parts[6] # could be different
    dur_min = int(parts[1])
    dur_sec = int(parts[2])
    dur = dur_min*60 + dur_sec
    
    date_time_str = date + time.replace(' ', '')
    print(date_time_str)
    date_time_obj = datetime.strptime(date_time_str, '%y%m%d%H%M%S')
    output = date_time_obj.strftime(output_format)
    print(Time(date_time_obj).yday)
    return Time(date_time_obj), output, dur, src_name


def cs_2_fil(par, scan_params):
    """
    Conversion from cs to fil files
    """
    print("\n\n---------------------------------------------")
    print("Converting from cs to fil")
    print("---------------------------------------------")

    scanno = scan_params['scanno']
    data_amount = scan_params['dur']
    dm = scan_params['dm']
    tele = par.telescope

    tel_dict = {"robledo"   : "RO", 
                "canberra"  : "CN", 
                "goldstone" : "GS"}
    
    tel_id = tel_dict.get(tele)

    scan_file = par.outdir + "/" + par.inf_file
    print("cs_2_fil: scan_file", scan_file)

    # Extract inf from scan file
    mjd_start, start_time, dur_sec, src_name = read_scan_table(scan_file, scanno)

    # modify outdir for each scan
    scanno_formatted = "{:02d}".format(scanno)
    scan_outdir = par.outdir + "/s" + scanno_formatted + "-" + src_name
    print("cs_2_fil: scan_outdir", scan_outdir)

    cs_dir = "%s/cs" % (scan_outdir)
    if not os.path.exists(cs_dir):
        print("cs_2_fil: error: cs dir does not exist")
        sys.exit(0)
        
    if par.freq_band == "s":
        cs_infiles = glob.glob("%s/slcp-*.cs" % (cs_dir)) \
                          + glob.glob("%s/srcp-*.cs" % (cs_dir))
    else:
        cs_infiles = glob.glob("%s/xlcp-*.cs" % (cs_dir)) \
                          + glob.glob("%s/xrcp-*.cs" % (cs_dir))

    nsub = len(cs_infiles)/2 # per polarization

    band_labels = False
    for cs_file in cs_infiles:
        print("cs file: ", cs_file)
        basename = os.path.basename(cs_file)
        [band_name,batch,band_no,tmp1] = re.split('[-.]', basename)
        band_id = band_name + "-" + batch
        print("band_id: ", band_id)
        if not band_labels:
            if "lcp" in band_id:
                band1_id = par.freq_band + "lcp-" + batch
                band2_id = par.freq_band + "rcp-" + batch
                band_labels = True
            else:
              band1_id = par.freq_band + "rcp-" + batch
              band2_id = par.freq_band + "lcp-" + batch
              band_labels = True  
        

    band_id_list = [band1_id, band2_id]
    print("band_id_list = ", band_id_list)
    
    params = {"inf": par.inf_file, "scanno": scanno,  
              "vrad_dir": par.vrad_dir, "nproc": par.nproc, 
              "vrad_base": par.vrad_base, "out": par.outdir, 
              "freq": par.freq_band, "source": par.source,
              "tele": par.telescope, "amount": data_amount, 
              "dm": dm, "src" : srcdir}

    ezap = 2
    npol = 2
    with mp.Pool(npol) as pool:
        # Create a list of arguments for each call to vrad_2_vdr
        args_list = [(scan_params, nsub, ezap, scan_outdir, 
                      band_id_name, data_amount, srcdir, 
                      tel_id) for band_id_name in band_id_list]
        print("frb1: args_list ==>")
        print(args_list)
        # Use the pool to map the vrad_2_vdr function to the args_list
        pool.starmap(cs_2_fil_run, args_list)

    return


def cs_2_fil_run(scan_params, nsub, ezap, scan_outdir, band_id_name, 
                      data_amount, srcdir, tel_id):
    print("cs_2_fil_run===>")

    dm = scan_params['dm']
    nchan = scan_params['nchan']
    snr = scan_params['snr']
    w = scan_params['w']
    mw = scan_params['mw']
    nt = scan_params['nt']
    mc = scan_params['mc']

    snr_formatted = "{:.1f}".format(snr)
    dm_formatted = "{:.3f}".format(dm)
    w_formatted = "{:.1f}".format(w)
    mw_formatted = "{:.1f}".format(mw)

    odir = "out-"+ band_id_name
    odir = scan_outdir + "/" + odir
    cs_dir = "%s/cs" % (scan_outdir)
    log_file = band_id_name+".log"
    log_file = scan_outdir + "/" + log_file

    if not os.path.exists(odir):
        os.makedirs(odir, exist_ok=True)

    if tel_id is None:
        tel_str = ""
    else:
        tel_str = "-tel %s" %tel_id

    if mc > 0:
        cmd = "python3 -u %s/bb_proc.py" \
            " -dm %s -nc %d -nt %d -snr %s -ezap %d -w %s -mw %s " \
            " --badblocks %s -nsub %d -mc %d %s %s %s > %s" \
            % (srcdir, dm_formatted, nchan, nt, snr_formatted, ezap, \
               w_formatted, mw_formatted, tel_str, nsub, mc, cs_dir, \
               band_id_name, odir, log_file)
    else:
        cmd = "python3 -u %s/bb_proc.py" \
            " -dm %s -nc %d -nt %d -snr %s -ezap %d -w %s -mw %s " \
            " --badblocks %s -nsub %d %s %s %s > %s" \
            % (srcdir, dm_formatted, nchan, nt, snr_formatted, ezap, \
               w_formatted, mw_formatted, tel_str, nsub, cs_dir, \
               band_id_name, odir, log_file)

    print("cmd>>>", cmd)
    print("")
    subprocess.run(cmd, shell=True)

    return

@click.command()
@click.option("--parfile", type=str,
              help="parfile [required]", required=True)
#@click.option("--data_amount", type=int,
#              help="Amount of data to write (in seconds)", required=False)
@click.option("--scanno", type=int,
              help="Scan number (from scan.table)", required=True)
@click.option("--dm", type=float,
              help="DM", required=True)
@click.option("--nchan", type=int,
              help="Number of channels", required=True)
@click.option("--snr", type=float,
              help="SNR threshold", required=True)
@click.option("--w", type=float,
              help="width", required=True)
@click.option("--nt", type=int,
              help="Number of threads", required=True)
@click.option("--mw", type=float,
              help="Maximum width", required=True)
@click.option("--mc", type=int,
              help="Maximum cands", required=True)
@click.option("--skip", type=int,
              help="Number of seconds to skip", required=True)
@click.option("--dur", type=int,
              help="Duration (seconds)", required=True)

def run_pipeline(parfile, scanno, dur, skip, dm, nchan, snr, w, mw, nt, mc):
    """
    Run the pipeline
    """

    #print("run_pipeline: cur_dir = ", cur_dir)
    #print("run_pipeline: srcdir  = ", srcdir)

    #print("main: run_pipeline: dur", dur)
    if dur is None:
        dur = 0
    #print("main: run_pipeline: dur", dur)
        
    # Get pars 
    par = get_pars(parfile)
    
    print('vrad', par.vrad_to_cs)
    print('cs', par.cs_to_fil)
    print('plotfil', par.plot_fil)
    print('plot_scint', par.plot_scint)
    print('combine', par.combine_pol)

    print('vrad_dir', par.vrad_dir)
    print('vrad_base', par.vrad_base)

    print('cur_dir ',cur_dir)
    print('srcdir ',srcdir)

    print('scanno ', scanno)
    print('dm ', dm)
    print('dur ', dur)
    print('skip ', skip)
    print('nchan ', nchan)
    print('w ', w)
    print('mw ', mw)
    print('snr ', snr)
    print('nt ', nt)
    print('mc ', mc)

    print('outdir', par.outdir)

    scan_params = {"scanno": scanno, "dur": dur, "skip": skip, "dm": dm, 
                   "nchan": nchan, "w": w, "mw": mw, "snr": snr, "nt": nt, "mc": mc}


    print("scan_params >>>", scan_params.values())
    print("scanno >>>", scan_params['scanno'])
    
    times = []

    if par.vrad_to_cs:
        st = time.time()
        vrad_2_cs(par, scan_params)
        times.append(["vrad->cs", round(time.time()-st, 2)])
    

    if par.cs_to_fil:
        st = time.time()
        cs_2_fil(par, scan_params)
        times.append(["cs->fil", round(time.time()-st, 2)])
    
       
    total = sum(item[1] for item in times)
    times.append(["total", total])
    print(tabulate(times, headers =['Step','Time (s)'], 
                   tablefmt = 'fancy_grid'))
 
        
    return


debug = 0

if __name__ == "__main__":
    if not debug: 
        run_pipeline()
