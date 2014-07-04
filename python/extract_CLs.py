import os
import ROOT as r
def print_and_save_CLs(directory,lim_tree_dir,xsec=1.0):
    f = r.TFile(os.path.join(directory,"outputfile.root"),"READ")
    fout = open(os.path.join(directory,'exclusion.txt'),'w')
    print(lim_tree_dir)
    #tree = f.Get('{0}/LimitTree'.format(lim_tree_dir))
    tree = f.Get('{0}/AnalysisTree'.format(lim_tree_dir))
    total_yeild=0
    for entry in tree:
      for sig in entry.sigeff:
        #excl = 100*(1.0 - entry.cls[0])
        #print("{0} : {1}".format(directory,excl))
        yeild=(sig*xsec*19.5E15*0.3*0.3)/36 #Coupling squared!
        print sig,yeild
        fout.write(str(sig)+" "+str(yeild)+"\n")
    fout.close()

