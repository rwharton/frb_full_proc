import your 
from your.formats.filwriter import make_sigproc_object


def write_snippet(data, in_yr, outbase, toffset):
    """
    Write snippet data filterbank from input your file 

    toffset in seconds from start of file
    """
    hdr_in = in_yr.your_header

    outfile = "%s.fil" %outbase

    tstart = hdr_in.tstart + (toffset / (24 * 3600.))

    # Setup output filterbank file
    sig_obj = make_sigproc_object(
            rawdatafile   = outfile,
            source_name   = hdr_in.source_name,
            nchans        = hdr_in.nchans, 
            foff          = hdr_in.foff, 
            fch1          = hdr_in.fch1, 
            tsamp         = hdr_in.tsamp, 
            tstart        = tstart, 
            src_raj       = in_yr.src_raj, 
            src_dej       = in_yr.src_dej, 
            machine_id    = in_yr.machine_id,
            nbeams        = in_yr.nbeams,
            ibeam         = in_yr.ibeam,
            nbits         = in_yr.nbits,
            nifs          = in_yr.nifs,
            barycentric   = in_yr.barycentric,
            pulsarcentric = in_yr.pulsarcentric ,
            telescope_id  = in_yr.telescope_id,
            data_type     = in_yr.data_type,
            az_start      = in_yr.az_start,
            za_start      = in_yr.za_start)

    # Write header
    sig_obj.write_header(outfile)

    # write data 
    sig_obj.append_spectra(data, outfile)

    return 
