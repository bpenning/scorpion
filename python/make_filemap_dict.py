import commands
import os

def filemap_from_single_diroctory(dir,rootfile,experiment):
    #UGLY but this works
    (err, xsec) = commands.getstatusoutput('grep "TOTAL" '+dir+'/pythia-output.log | grep -o "[0-9]\{1,100\}\.[0-9]\{1,100\}[E\-]\{0,2\}[0-9]\{0,3\}"')
    #to be robust against giving the absolute path to the rootfile
    rootfile=os.path.basename(rootfile)
    rootfile_path=os.path.join(dir,rootfile)
    internal_name,_=os.path.splitext(rootfile)
    return {
            'xsec':float(xsec),
            'rootfiles':[rootfile_path],
            'internal_name':internal_name,
            'experiment':experiment,
            }
