import os
import tempfile
from ftplib import FTP
from collections import namedtuple

def uploadfig(pylab_module, upload_configuration,dpi=150):
    """Upload plot to some server
    
    :param pylab_module: The Pylab-module to use
    :param upload_configuration: At least containing .protocol to determine appropriate method
    :type upload_configuration: object
    :param dpi: Resolution of plot 
    :type dpi: integer
    """
    p = pylab_module
    if not upload_configuration.protocol.upper() in ["FTP"]:
        raise ValueError("Unknown protocol")
    fh,fn = tempfile.mkstemp(suffix=".png")
    #fh.close()
    p.savefig(fn,dpi=dpi)
    if upload_configuration.protocol.upper() == "FTP":
        _uploadfig_ftp(fn,upload_configuration)

def _uploadfig_ftp(fn, upload_configuration):
    uc = upload_configuration
    ftp = FTP(uc.server)
    ftp.login(uc.username,uc.pwd)
    ftp.cwd(uc.directory)
    ftp.storbinary("STOR " + os.path.split(fn)[1], open(fn,"rb"),1024)

__all__ = [uploadfig]

if __name__=="__main__":
    import pylab as p
    
    UploadConfiguration = namedtuple( "UploadConfiguration", ["protocol", "server","username","pwd","directory"]) 
    uc = UploadConfiguration("ftp", "www.thorstenkranz.de","f006149e","easypeasy","/")

    xs = p.arange(0,10,0.1)
    p.plot(xs,p.cos(xs))
    uploadfig(p,uc)
