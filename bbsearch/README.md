# dsn_baseband

baseband_search.py - A FRB search pipeline for channelized baseband data output from digifil / baseband repo.


### RUN NOTES:

Each filterbank file must have a .fil extension. They get moved to sub-directories, searched, and candidates extracted.
Don't forget to update scriptdir at the top of the script! This should coincide with the location of the code/repo.
The python3 branch can be run in python3.

Update the script location path at the top of baseband_search.py:

```
scriptdir = "/data3/pearlman/baseband_pipeline/dev"
```

### Run Example:

```
python baseband_search.py --obsroot sband --nsub 7 --nchansub 32 --edgezap 2 --zapchans 0,1,2 --dm 87.77 --spthresh 6.0 --width 0.100 --telescope RO --dsn
```

```
(base) [pearlman@cenx3:dev]$ python baseband_search.py --help
Usage: baseband_search.py [OPTIONS]

Options:
  --obsroot TEXT      Root name of candidate directory (for each filterbank).
  --nsub INTEGER      Number of subbands in baseband data.
  --nchansub INTEGER  Number of channels per subband.
  --edgezap INTEGER   Number of channels to zap at top and bottom of each
                      subband.

  --zapchans TEXT     Comma separated list of channels to zap. Index 0
                      corresponds to the highest frequency channel.

  --dm FLOAT          Dispersion measure (for incoherent, inter-channel
                      dedispersion and searching).

  --spthresh FLOAT    S/N threshold to use for single pulse search.
  --width FLOAT       Length of data to extract around each candidate burst
                      (units: seconds).

  --telescope TEXT    Telescope (Madrid=RO, Goldstone=GS, Canberra=CN).
  --dsn               Set a dummy machine_id (999) for the DSN.
  --help              Show this message and exit.
```

### SOFTWARE DEPENDENCIES:

python: numpy, h5py, click, logging [pip install these]

sigproc (https://github.com/aaronpearlman/sigproc_pearlman_dsn)

presto (https://github.com/aaronpearlman/presto_dsn)
