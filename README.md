# frb_full_proc
Convert, process, and search DSN baseband data.

This pipeline will take telescope baseband data in 
`.vrad` form, convert to `.vdr`, then convert to 
`.cs` files in SIGPROC format.  From the 
`.cs` baseband files, it will then convert them to 
a filterbank with arbitrary number of frequency channels, 
then de-disperse and search for single pulses. Small 
snippets of data around each candidate will be written 
to filterbank format and plots will be made for each.


## Usage 
Running with the `--help` option will show the following usage:

    Usage: frb_full_proc.py [OPTIONS]
    
      Run the pipeline
    
    Options:
      --parfile TEXT    parfile [required]  [required]
      --scanno INTEGER  Scan number (from scan.table)  [required]
      --dm FLOAT        DM  [required]
      --nchan INTEGER   Number of channels  [required]
      --snr FLOAT       SNR threshold  [required]
      --w FLOAT         width  [required]
      --nt INTEGER      Number of threads  [required]
      --mw FLOAT        Maximum width  [required]
      --mc INTEGER      Maximum cands  [required]
      --skip INTEGER    Number of seconds to skip  [required]
      --dur INTEGER     Duration (seconds)  [required]
      --help            Show this message and exit.
    

Here's an example of a recent processing run:

    python /src/frb_proc_full/frb_proc_full.py --parfile params_22m295_sband.txt 
           --scanno 1 --dur 10 --dm 219.46 --skip 0 --nchan 400 --nt 16 --snr 6.0 
           --w 0.1 --mw 0.1 --mc 100

The parameter file (`--parfile`) gives a text file with various 
observation specific parameters and options (see below for an 
example).  The scan number (`--scanno`) gives the scan to be 
processed.  The targets for each scan are given by the scan table, 
which is specified in the parameter file.  The duration (`--dur`) 
gives the number of seconds of data to process (to process it all, 
just give 0).  The skip parameter (`--skip`) gives the number of 
seconds to skip after the start of the scan.    

In this case, we are removing intrachannel delays by coherently 
de-dispersing at DM=219.46 (`--dm`) and producing a filterbank with 
a total of 400 channels (`--nchan`), which is how we determine the 
output time resolution.  We are using 16 threads (`--nt`) for the 
filterbank conversion with `digifil`.

For the search, we are searching candidates up to a maximum width 
of 0.1 seconds (`--mw`) and a min snr (`--snr`) of 6.  For each 
candidate we will extract 0.1 seconds (`--w`) of data and write 
a snippet filterbank file that will be used to make candidate 
waterfall plots.  To avoid being swamped by candidates, you can 
set a max number (`--mc`).  If the number of candidates exceeds 
this value, then no plots will be made (except for a few for diagnostic 
purposes).

## Parameter File

Here is an example parameter file. Most of this should be 
self-explanatory.  Right now telescope options are just the 
DSN stations ("robledo", "goldstone", "canberra").

    #####################
    ###   OBS INFO    ###
    #####################
    
    freq_band = "s" # set to x or s
    source = "FRB220912A"
    telescope = "robledo"
    
    ############################
    ###  OUTPUT DIRECTORIES  ###
    ############################
    
    # Output directory
    outdir = "."
    
    ##########################
    ###  INPUT VRAD FILES  ###
    ##########################
    
    # Directory containing the vrad files
    vrad_dir = "/raw/22m295"
    
    # vrad base
    vrad_base = "22-295"
    
    # Names of info file
    #  -->  Place these in output directory
    inf_file = "scan.table.22m295"
    
    ######################
    ###  STEPS TO RUN  ###
    ######################
    
    # Use (0 or False) and (1 or True) to
    # set which steps to run
    
    vrad_to_cs  = 1
    cs_to_fil   = 1
    
    # max number of processes to run at
    # once with parallel processing
    nproc = 8
