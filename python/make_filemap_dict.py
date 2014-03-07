import commands
import os

def filemap_from_pythia6_single_diroctory(dir,rootfile,experiment):
    #UGLY but this works
    (err, xsec) = commands.getstatusoutput('grep "TOTAL" '+dir+'/pythia-output.log | grep -o "[0-9]\{1,100\}\.[0-9]\{1,100\}[E\-]\{0,2\}[0-9]\{0,3\}"')
    #very important: convert from milibarn to barn
    xsec=float(xsec)*1e-3
    #to be robust against giving the absolute path to the rootfile
    rootfile=os.path.basename(rootfile)
    rootfile_path=os.path.join(dir,rootfile)
    internal_name,_=os.path.splitext(rootfile)
    return {
            'xsec':xsec,
            'rootfiles':[rootfile_path],
            'internal_name':internal_name,
            'experiment':experiment,
            }

def filemap_from_pythia8_single_diroctory(dir,rootfile,experiment):
    #get xsection
    pythialog=os.path.join(dir,'pythia-output.log')
    with open(pythialog,'r') as f:
        relevant_block=False
        for line in f:
            if 'PYTHIA Event and Cross Section Statistics' in line:
                relevant_block=True
            if 'End PYTHIA Event and Cross Section Statistics' in line:
                relevant_block=False
            if relevant_block and ('sum' in line):
                xsec_mbarn=line.split('|')[-2].split()[0]
                xsec=0.001*float(xsec_mbarn)
    #to be robust against giving the absolute path to the rootfile
    rootfile=os.path.basename(rootfile)
    rootfile_path=os.path.join(dir,rootfile)
    internal_name,_=os.path.splitext(rootfile)
    return {
            'xsec':xsec,
            'rootfiles':[rootfile_path],
            'internal_name':internal_name,
            'experiment':experiment,
            }
