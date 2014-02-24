import os
import ROOT as r
def print_and_save_CLs(directory,lim_tree_dir):
    f = r.TFile(os.path.join(directory,"outputfile.root"),"READ")
    fout = open(os.path.join(directory,'exclusion.txt'),'w')
    print(lim_tree_dir)
    tree = f.Get('{0}/LimitTreec'.format(lim_tree_dir))
    for entry in tree:
        excl = 100*(1.0 - entry.cls[0])
        print("{0} : {1}".format(directory,excl))
        fout.write(str(excl))
    fout.close()

