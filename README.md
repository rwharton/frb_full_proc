# bb_proc
End to end baseband processing for DSN data.  

This pipeline will take `.cs` baseband files as input, 
convert them to a filterbank with arbitrary number of 
frequency channels, then de-disperse and search for 
single pulses. Small snippets of data around each 
candidate will be written to filterbank format and 
plots will be made for each.


## Usage 
Running with the `-h` option will show the following usage:

    usage: bb_proc.py [-h] -dm DM -nc NCHAN [-m MEMLIM] [-nt NTHREAD] 
                      [-rt RFIDEC] [-snr SNRMIN] [-ezap EDGEZAP] [-nsub NSUB]
                      [-w WIDTH] [-mw MAXWIDTH] [--zerodm] [--badblocks] [-tel TEL]
                      csdir basename outdir
    
    Pipeline to process and search baseband data
    
    positional arguments:
      csdir                 Directory containing *.cs files
      basename              Base name of cs file: {basename}*.cs
      outdir                Output data directory
    
    optional arguments:
      -h, --help            show this help message and exit
      -dm DM, --dm DM       DM for intra-channel coherent dedispersion
      -nc NCHAN, --nchan NCHAN
                            Total number of output channels in filterbank
      -m MEMLIM, --memlim MEMLIM
                            Max memory to use during DADA conversion in GB (def: 16)
      -nt NTHREAD, --nthread NTHREAD
                            Number of threads for digifil processing (def: 1)
      -rt RFIDEC, --rfidec RFIDEC
                            Number of samples to decimate filfile for RFI 
                            diagnostics (def: 512, to skip this: -1)
      -snr SNRMIN, --snrmin SNRMIN
                            Minimum SNR for single pulse search (def: 6.0)
      -ezap EDGEZAP, --edgezap EDGEZAP
                            Subband edge channels to zap (def: 2)
      -nsub NSUB, --nsub NSUB
                            Number of subbands in data (def: 1)
      -w WIDTH, --width WIDTH
                            Output size for candidate snippets in sec (def: 0.1)
      -mw MAXWIDTH, --maxwidth MAXWIDTH
                            Max boxcar width in ms for SP search (def: 10)
      -mc MAXCANDS, --maxcands MAXCANDS
                        Max cands for plotting. Do not plot if exceeds this
                        number. (def = -1, no lim)
      -f FILTER, --filter FILTER
                            Filter frequency - comma separated list giving the
                            center frequency, the number of harmonics beyond
                            fundamental, and width in Hz (e.g., '60.0,5,1.0' to
                            60Hz signal and 5 harmonics with width 1.0Hz). To zap
                            multiple frequencies, repeat this argument
      --zerodm              Apply the zero DM filter when dedispersing
      --badblocks           Ignore bad blocks in single pulse search
      -tel TEL, --tel TEL   DSN Telescope Name GS/RO/CN (def: RO)

For better or worse, it's a lot of options.

Here's an example of a recent processing run:

    python /src/bb_proc/bb_proc.py -dm 219.46 -nc 128 -nt 16 -m 64 
           -rt 500 -snr 8 -ezap 2 -w 0.1 -nsub 4 --badblocks -tel RO 
           -f 60,10,0.5 -mw 20 -mc 9000 /shared/frb20220912A/22m316/xband xlcp-0001 .

In this case, we are removing intrachannel delays by coherently 
de-dispersing at DM=219.46 (`-dm`) and producing a filterbank with 
a total of 128 channels (`-nc`), which is how we determine the output 
time resolution.  We are using 16 threads (`-nt`) for the filterbank 
conversion with `digifil` and are allowing for up to 64 GB of memory 
(`-m 64`) to be used in the conversion. 

The baseband data files are 
found in `/shared/frb20220912A/22m316/xband` and have a basename of 
`xlcp-0001` (meaning the sub-bands are called, e.g., `xlcp-0001-01.cs`). 
It is important to give the basename as the full name before the `-XX.cs` 
at the end of the files.  There are 4 channels in the baseband data (`-nsub`), 
which is important for flagging sub-band edge channels.  In this case, 
we will zap two channels (`-ez`) at the edges of the sub-bands. 

For the searching, we will first de-disperse to DM=219.46.  If you want 
to remove certain periodicities in the de-dispersed time series, you can 
use the filter option.  In this case we use `-f 60,10,0.5` which sets the 
fundamental at 60 Hz, selects 10 harmonics, and sets a zapping width of 
0.5 Hz.  The filter will be run on the dat file before searching.
The search is done with a modified version of PRESTO's `single_pulse_search.py` 
with much larger boxcar widths.  Here, we are setting the minimum SNR limit 
to be 8 (`-snr`), the maximum boxcar template width to be 20ms (`-mw`) and 
are ignoring candidates 
found in bad blocks (`--badblocks`).  For the candidates found we will create 
snippet filterbank files of 0.1 seconds around each candidate (`-w`).  Using 
these, we will make plots.  To avoid bursts of RFI, we will exclude from 
plots any candidates that fall during a minute with more than 500 bursts 
per minute (`-rt`).

We can also optionally set an upper limit to the candidates (`-mc`) beyond 
which candidate filterbanks and plots will **not** be made.  This is a minor 
precaution against observations with lots of RFI or strong 60 Hz signal where 
the code would take forever to make useless plots.  If the number of candidates 
found by the single pulse search exceeds the value given by `-mc` (in our 
example here, 9000), then plots of the first 10 candidates will be made for 
diagnostic purposes and then no other candidates will be plotted.  There will 
still be a summary plot showing the number of candidates over time.

The `-tel` option allows you to set the telescope name in the output filterbank 
header.  Right now, you can just give one of the three DSN dish names.

All of the files produced in the processing will be place in the output 
directory.  In the case of this example, that is just the current working 
directory.



