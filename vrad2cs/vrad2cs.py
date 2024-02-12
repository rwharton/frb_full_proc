import numpy as np
import os 
import sys
import glob
import time
import subprocess
import click
from datetime import datetime, timedelta
from astropy.time import Time, TimeDelta
#import sigproc as fb
import multiprocessing as mp
import re

cur_dir = os.path.realpath(__file__)
srcdir  = cur_dir.rsplit('/', 1)[0]

def extract_start_time(inf_file):
    """
    Given inf_file/scan_table, find MJD time and return it 
    in YY/Day_HH:MM:SS format
    """
    if "scan.table" in inf_file:
        # Define the output format
        output_format = "%y/%j_%H:%M:%S"
        # Open the scan table file
        with open(inf_file, 'r') as file:

            line = file.readline().rstrip()

            # Split the line into parts
            parts = line.split()

            # Extract the date and time parts
            date = parts[3] # could be different
            time = parts[4] + parts[5] + parts[6] # could be different

            date_time_str = date + time.replace(' ', '')
            print(date_time_str)
            date_time_obj = datetime.strptime(date_time_str, '%y%m%d%H%M%S')
            output = date_time_obj.strftime(output_format)
            print(Time(date_time_obj).yday)
            return Time(date_time_obj), output
    else:
        epoch_key = "Epoch of observation (MJD)"
        mjd = None

        with open(inf_file, 'r') as file:
            lines = file.readlines()
            for line in lines:
                if epoch_key in line:
                    _, value = line.split('=', 1)
                    mjd = value.strip()
                    break

        t = Time(mjd, scale='utc', format='mjd')
        print(t.yday)

        # convert to required format
        input_format = "%Y:%j:%H:%M:%S.%f"
        dt = datetime.strptime(t.yday, input_format)

        # Format the datetime object using strftime
        output_format = "%y/%j_%H:%M:%S"
        formatted_time = dt.strftime(output_format)
        return t, formatted_time


def extract_start_time_from_scan_table(inf_file,scanno):
    """
    Given scan_table, find MJD time and return it 
    in YY/Day_HH:MM:SS format
    """
    # Define the output format
    output_format = "%y/%j_%H:%M:%S"
    # Open the scan table file
    with open(inf_file, 'r') as file:
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


def vrad_2_vdr(pulse_times, vrad_file, out_dir, 
               data_amount, srcdir=srcdir):
    """
    Converts raw telescope data (*.vrad) given list of times 
    of pulses to vdr format
    
    Executes script in format:

    rdef_io -V -f TIME -n 60 INPUT_FILE OUTPUT_FILE

    rdef_loc: location of rdef_io script
    """
    # Create directory to put vdr files in 
    vdr_dir = "%s/vdr" % (out_dir)
    if not os.path.exists(vdr_dir):
        os.makedirs(vdr_dir, exist_ok=True)

    # replace with vdr suffix to create output path + place 
    # in vdr folder
    basename = os.path.basename(vrad_file)
    basename = os.path.splitext(basename)[0]
    #print("wam: vrad_2_vdr: basename: ", basename)

    # parse basename
    [yy,doy,ch,dss,tmp1,band] = re.split('[-_]', basename)
    band = band.lower()
    ch = int(ch)
    ch = "%02d" % ch
    outfile = band + "-0001-" + str(ch) + ".vdr"
    
    out_file = "%s/%s" % (vdr_dir, outfile)

    # wait for process to finish before executing another
    rdef_loc = "%s/rdef_io" %srcdir
    #rdef_loc = "/home/majid/software/olr_tools/rdef_io"

    t = pulse_times
    
    cmd = "%s -V -f %s -n %d %s %s" % (rdef_loc, t, data_amount, 
                                       vrad_file, out_file)
        
    # if file exists already, we don't need to run
    if os.path.isfile(out_file):
        #os.remove(out_file)
        print("File: %s already exists" %out_file)
    else:
        print(cmd)
        print("")
        subprocess.run(cmd, shell=True)

    return


def vdr_2_cs(vdr_file, out_dir, freq_band, source_name, freq, bw, 
             telescope, srcdir=srcdir, band_num=1):
    """
    Converts *.vdr files to *.cs files
    Executes script in format:

    vsrdump_char_100 {vdr_file} -o {output_file} 
    -name {source_name} -comp -freq {freq} -bw {bw} -{telescope} -vsr
    """
    # Create directory to put cs files in 
    cs_dir = "%s/cs" % (out_dir)
    if not os.path.exists(cs_dir):
        os.makedirs(cs_dir, exist_ok=True)

    # replace with cs suffix to create output path + place in cs folder
    basename = os.path.basename(vdr_file)
    basename = os.path.splitext(basename)[0]

    out_file = "%s/%s.cs" % (cs_dir, basename)

    if freq_band == "s":
        vsrdump_loc = "%s/vsrdump_char_100" %srcdir
    else: 
        vsrdump_loc = "%s/vsrdump_char" %srcdir

    cmd = "%s %s -o %s -name %s -comp -freq %.2f -bw %.2f -%s -vsr" \
        % (vsrdump_loc, vdr_file, out_file, source_name, 
           freq, bw, telescope)

    # if file exists already, don't need to run
    if os.path.isfile(out_file):
        #os.remove(out_file)
        print("File: %s already exists" %out_file)
    else:
        print("cmd>>>: ", cmd)
        print("")
        subprocess.run(cmd, shell=True)

    return


def cleanup_vdr(vdr_dir, cs_dir, freq_band):
    """
    Check that cs files exist, then remove vdr
    """
    # Get cs files
    if freq_band == "s":
        glob_str = "%s/s*cs" %(cs_dir)
    else:
        glob_str = "%s/x*cs" %(cs_dir)
    
    cs_files = glob.glob(glob_str)

    print("\n\n")
    print("########################")
    print(" Cleaning up vdr files ")
    print("########################")
    print("\n\n")
    
    # For each cs file, look for corresponding vdr
    for cs_file in cs_files:
        cs_fn = cs_file.split('/')[-1]
        cs_bn = cs_fn.rsplit('.', 1)[0]
        vdr_file = "%s/%s.vdr" %(vdr_dir, cs_bn)
        # if it exists, remove it
        if os.path.isfile(vdr_file):
            print("%s exists!" %cs_file)
            print("Removing %s" %vdr_file)
            os.remove(vdr_file)
        else:
            pass
    
    return
    

@click.command()
@click.option("--inf_file", type=str, 
              help="Full path of .inf file. ***You can also input a scan table", required=True)
#@click.option("--sp_file", type=str, 
#              help="Full path of single pulse file", required=True)
@click.option("--vrad_dir", type=str,
              help="Directory containing vrad files", required=True)
@click.option("--vrad_base", type=str,
              help="Basename of vrad files {vrad_base}*.vrad", required=True)
@click.option("--out_dir", type=str,
              help="Output directory of *.cs and intermediate *.vdr files", required=True)
@click.option("--freq_band", type=str,
              help="Frequency band to process (s/x)", required=True)
@click.option("--dm", type=float,
              help="Dispersion Measure", required = True)
@click.option("--source", type=str,
              help="Name of the pulse source", required=True)
@click.option("--telescope", type=str,
              help="Name of telescope", required=True)
@click.option("--data_amount", type=int,
              help="Amount of data to write (in seconds)", required=True)
@click.option("--skip", type=int,
              help="Amount of data to skip at the start (in seconds)", required=True)
#@click.option("--data_amount", type=float,
#              help="Amount of data to write (in seconds)", required=True)
@click.option("--nproc", type=int,
              help="Number of processes to use in multiprocessing", 
              required=True)
@click.option("--scanno", type=int,
              help="Scan number (from scan.table)", required=True)

def vrad_2_cs(inf_file, scanno, vrad_dir, vrad_base, out_dir, 
              freq_band, dm, source, telescope, skip, data_amount,
              nproc, srcdir=srcdir):
    # dictionary containing freq for every freq band
    freq_dct = {"s": 2250, "x1": 8224, "x2": 8256, "x3": 8288, "x4": 8320}
    bw_dct = {"s": 100, "x": 32}

    # Extract inf from scan file
    mjd_start, start_time, dur_sec, src_name = extract_start_time_from_scan_table(inf_file, scanno)
    
    # Take each vrad file and convert to vdr
    print("\n\nvrad -> vdr...")
    
    if freq_band == "s":
        vrad_infiles = glob.glob("%s/%s*SLCP.vrad" % (vrad_dir, vrad_base)) +\
                       glob.glob("%s/%s*SRCP.vrad" % (vrad_dir, vrad_base))
    else:
        vrad_infiles = glob.glob("%s/%s*XLCP.vrad" % (vrad_dir, vrad_base)) +\
                       glob.glob("%s/%s*XRCP.vrad" % (vrad_dir, vrad_base))


    print("vrad2cs:vrad_2_cs: nproc      = ", nproc)

    if len(vrad_infiles) == 0:
        print("No *vrad files found!")
        print("vrad_dir: %s" %vrad_dir)
        print("vrad_base: %s" %vrad_base)
        sys.exit(0)

    new_mjd_time = mjd_start + timedelta(seconds=skip)
    output_format = "%y/%j_%H:%M:%S"
    converted_times = new_mjd_time.strftime(output_format)
    print("vrad2cs: vrad_2_cs: converted_times: ", converted_times)

    nproc = len(vrad_infiles)
    if data_amount == 0:
        data_amount  = dur_sec

    # modify outdir for each scan
    scanno_formatted = "{:02d}".format(scanno)
    scan_outdir = out_dir + "/s" + scanno_formatted + "-" + src_name
    print("vrad2cs:vrad_2_cs: outdir      = ", out_dir)
    print("vrad2cs:vrad_2_cs: scan_outdir = ", scan_outdir)

    with mp.Pool(nproc) as pool:
        # Create a list of arguments for each call to vrad_2_vdr
        args_list = [(converted_times, vrad_file, scan_outdir, 
                      data_amount, srcdir) for vrad_file in vrad_infiles]
        print("type of args_list = ", type(args_list))
        print("len  of args_list = ", len(args_list))
        print(args_list)
        # Use the pool to map the vrad_2_vdr function to the args_list
        pool.starmap(vrad_2_vdr, args_list)

    
    # Take vdr files and convert to cs
    print("\n\nvdr -> cs...")
    # We know all vdr files will be contained in vdr directory
    vdr_infiles = glob.glob("%s/vdr/*" % (scan_outdir))
    
    if len(vdr_infiles) == 0:
        print("No *vdr files found!")
        print("vdr_dir: %s/vdr" %scan_outdir)
        sys.exit(0)

    # if frequency band is s they will have same frequency
    if freq_band == "s":
        s_freq = freq_dct[freq_band]
        s_bw   = bw_dct[freq_band]
        with mp.Pool(nproc) as pool:
            # Create a list of arguments for each call to vdr_2_cs
            args_list = [(vdr_file, scan_outdir, freq_band, source, \
                  s_freq, s_bw, telescope, srcdir) for vdr_file in vdr_infiles]
            # Use the pool to map the vdr_2_cs function to the args_list
            pool.starmap(vdr_2_cs, args_list)

    # we have four different subbands for x band
    elif freq_band == "x":
        band_info = []
        for vdr_file in vdr_infiles:
            band_num = int(re.search(r'-00(\d+)', vdr_file).group(1))
            subband = "x" + str(band_num)
            x_freq_sb = freq_dct[subband]
            x_bw_sb = bw_dct[freq_band]
            # associate each file with freq and bw for its 
            # value (1, 2, 3, 4)
            band_info.append((vdr_file, x_freq_sb, x_bw_sb, band_num))

        with mp.Pool(nproc) as pool:
            # Create a list of arguments for each call to vdr_2_cs
            args_list = [(vdr_file, scan_outdir, freq_band, source, \
                          x_freq_sb, x_bw_sb, telescope, srcdir, band_num) 
                         for vdr_file, x_freq_sb, x_bw_sb, band_num in band_info]

            # Use the pool to map the vdr_2_cs function to the args_list
            pool.starmap(vdr_2_cs, args_list)

    else:
        print("freq_band must be s or x")

    # cleanup
    vdr_dir = "%s/vdr" %(scan_outdir)
    cs_dir = "%s/cs" %(scan_outdir)
    cleanup_vdr(vdr_dir, cs_dir, freq_band)

    return 

debug = 0

if __name__ == "__main__":
    if not debug: 
        vrad_2_cs()
