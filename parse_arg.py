def get_par_dict(infile):
    """
    Read an input file with parameter names 
    and values and build a dictionary from them

    Skip comments and blank lines, split on "="

    Note that all values will be strings and 
    we'll deal with types later
    """
    par_dict = {}
    with open(infile, 'r') as fin:
        for line in fin:
            if line[0] in ["\"", "#", " ", "\n"]:
                continue
            else: pass

            cols = line.split('=')

            if len(cols) == 1:
                continue
            else: pass

            key = cols[0].strip()
            val = cols[1].strip()

            # remove comments 
            if '#' in val:
                val = val.split('#')[0]
                val = val.strip()

            # remove quotes around strings
            if "\"" in val:
                val = val.strip("\"")
            if "\'" in val:
                val = val.strip("\'")

            par_dict[key] = val

    return par_dict


class Parameters():
    """
    Class object to store the parameters from 
    the par file. Takes as input the dict of 
    parameters
    """
    def __init__(self, pd):
        # Output directory
        self.outdir =  self.pval(pd, 'outdir', 's')

        # Multiprocessing
        self.nproc = self.pval(pd, 'nproc', 'd', 1)
        
        # vrad2cs conversion
        self.inf_file = self.pval(pd, 'inf_file', 's')
        self.sp_file  = self.pval(pd, 'sp_file', 's')
        self.sp_dt_us = self.pval(pd, 'sp_dt_us', 'f')
        self.sp_filter  = self.pval(pd, 'sp_filter', 'b')
        self.sp_snr_min = self.pval(pd, 'sp_snr_min', 'f')
        self.sp_match_fac = self.pval(pd, 'sp_match_fac', 'f')
        self.vrad_dir  = self.pval(pd, 'vrad_dir', 's')
        self.vrad_base = self.pval(pd, 'vrad_base', 's')
        
        self.freq_band = self.pval(pd, 'freq_band', 's')
        self.source = self.pval(pd, 'source', 's')
        self.telescope = self.pval(pd, 'telescope', 's')
        self.data_amount = self.pval(pd, 'data_amount', 'f')

        # inputs for cs2fil conversion
        self.cs_base = self.pval(pd, 'cs_base', 's')
        self.cs_dir  = self.pval(pd, 'cs_dir', 's')
        self.dada_dir  = self.pval(pd, 'dada_dir', 's')
        self.fil_dir   = self.pval(pd, 'fil_dir', 's')
        self.dm = self.pval(pd, 'dm', 'f')
        self.nchans = self.pval(pd, 'nchans', 'dl')
        
        # inputs for combine_pol
        self.lcp_base   = self.pval(pd, 'lcp_base', 's')
        self.rcp_base   = self.pval(pd, 'rcp_base', 's')
        self.comb_base  = self.pval(pd, 'comb_base', 's')
        self.ez = self.pval(pd, 'ez', 'f')
        self.dthresh = self.pval(pd, 'dthresh', 'f')
        self.nwin = self.pval(pd, 'nwin', 'd')
        
        # inputs for plotfil
        self.offset_file = self.pval(pd, 'offset_file', 's')
        self.fil_base = self.pval(pd, 'fil_base', 's')
        self.png_dir  = self.pval(pd, 'png_dir', 's')
        self.npy_dir  = self.pval(pd, 'npy_dir', 's')
        self.tavg = self.pval(pd, 'tavg', 'd')
        self.tdur = self.pval(pd, 'tdur', 'f')

        # scint plot inputs -- update
        self.npy_base = self.pval(pd, 'npy_base', 's')

        # conversion steps
        self.vrad_to_cs = self.pval(pd, 'vrad_to_cs', 'b') 
        self.cs_to_fil = self.pval(pd, 'cs_to_fil', 'b') 
        self.plot_fil = self.pval(pd, 'plot_fil', 'b') 
        self.plot_scint = self.pval(pd, 'plot_scint', 'b') 
        self.combine_pol = self.pval(pd, 'combine_pol', 'b') 
        
        # dm optimization?
        self.dm_opt = self.pval(pd, 'dm_opt', 'b')
        
        self.dm_lo = self.pval(pd, 'dm_lo', 'f')
        self.dm_hi = self.pval(pd, 'dm_hi', 'f')
        self.dm_step = self.pval(pd, 'dm_step', 'f')
        self.dm_dir  = self.pval(pd, 'dm_dir', 's')

        # Cleanup
        self.cleanup = self.pval(pd, 'cleanup', 'b')
        self.vdr_dir  = self.pval(pd, 'vdr_dir', 's')

    def pval(self, pdict, key, dtype, default=None):
        val = pdict.get(key)
        if val is None:
            val = default
        else:
            if dtype == 'd':
                val = int(val)
            elif dtype == 'f':
                val = float(val)
            elif dtype == 'b':
                if "True" in val:
                    val = True
                elif "False" in val:
                    val = False
                else:
                    val = bool( int(val) )
            elif dtype == 's':
                val = str(val)
            elif dtype == 'dl':
                vl = val.strip('[]').split(',')
                val = [ int(vv) for vv in vl ]
            elif dtype == 'fl':
                vl = val.strip('[]').split(',')
                val = [ float(vv) for vv in vl ]
            elif dtype == 'sl':
                vl = val.strip('[]').split(',')
                val = [ vv.strip() for vv in vl ]
            else:
                print("Unrecognized format")
                val = None 
        return val
