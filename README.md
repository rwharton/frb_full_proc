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

    python /src/bb_proc_full/bb_proc_full.py --parfile params_22m295_sband.txt 
           --scanno 1 --dur 10 --dm 219.46 --skip 0 --nchan 400 --nt 16 --snr 6.0 
           --w 0.1 --mw 0.1 --mc 100

The 

In this case, we are removing intrachannel delays by coherently 
de-dispersing at DM=219.46 (`--dm`) and producing a filterbank with 
a total of 128 channels (`--nchan`), which is how we determine the 
output time resolution.  We are using 16 threads (`--nt`) for the 
filterbank conversion with `digifil`.
