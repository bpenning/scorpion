import os
import ROOT as r
def print_and_save_CLs(directory):
    f = r.TFile(os.path.join(directory,"outputfile.root"),"READ")
    fout = open(os.path.join(directory,'exclusion.txt'),'w')
    for direc in f.GetListOfKeys():
	print direc.GetName()
	tree = f.Get('{0}/LimitTree'.format(direc.GetName())) 
        for entry in tree:
            excl = 100*(1.0 - entry.cls[0])
            print("{0} : {1}".format(directory,excl))
            fout.write(str(excl))
    fout.close()

