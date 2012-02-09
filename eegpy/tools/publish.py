"""Tools for plot publishing"""
import os
import tempfile
import ftplib

def uploadfig(pylab_module, upload_configuration,filename=None,dpi=150):
    """Upload plot to some server

    Uploads a matplotlib-figure to a server instead of showing it in a X window or writing 
    it to a file. Intended for situations where plots are created on a, e.g., SSH session.

    Support for various protocols are planned. Atm only FTP upload will work.
    
    :param pylab_module: The Pylab-module to use
    :param upload_configuration: At least containing .protocol to determine appropriate method
    :type upload_configuration: object
    :param dpi: Resolution of plot 
    :type dpi: integer

    .. code-block:: python
    
        from collections import namedtuple
        import matplotlib
        matplotlib.use("Agg")
        import pylab as p

        
        UploadConfiguration = namedtuple( "UploadConfiguration", 
                ["protocol", "server","username","pwd","directory"]) 
        uc = UploadConfiguration("ftp", "www.thorstenkranz.de","foo","bar","/")

        xs = p.arange(0,10,0.1)
        p.plot(xs,p.cos(xs))
        uploadfig(p,uc)
    """
    p = pylab_module
    if not upload_configuration.protocol.upper() in ["FTP"]:
        raise ValueError("Unknown protocol")
    if filename==None:
        fh,fn = tempfile.mkstemp(suffix=".png")
    else:
        fn_without_dir = os.path.split(filename)[1]
        if fn_without_dir!=filename:
            raise ValueError("filename must not contain any path, just filename.")
        fn = os.path.join(tempfile.mkdtemp(),fn_without_dir)
    #fh.close()
    p.savefig(fn,dpi=dpi)
    try:
        if upload_configuration.protocol.upper() == "FTP":
            _uploadfig_ftp(fn,upload_configuration)
    finally:
        os.unlink(fn)


def _uploadfig_ftp(fn, upload_configuration):
    uc = upload_configuration
    ftp = ftplib.FTP(uc.server)
    ftp.login(uc.username,uc.pwd)
    ftp.cwd(uc.directory)
    ftp.storbinary("STOR " + os.path.split(fn)[1], open(fn,"rb"),1024)

__all__ = [uploadfig]

