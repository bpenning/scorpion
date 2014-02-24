import optparse
from python.make_filemap_dict import filemap_from_single_diroctory 
from python.run_limit_code import runlim
from python.extract_CLs import print_and_save_CLs 

def parse_args():
    parser=optparse.OptionParser()
    parser.add_option('--pythia-delphes-dir')
    parser.add_option('--jaf-output-dir')
    parser.add_option('--CM-energy',default=7,type=int)
    parser.add_option('--experiment',default='CMS7')
    parser.add_option('--rootfile',default='delphes-output.root')
    options, args=parser.parse_args()
    return options

if __name__=="__main__":
    args=parse_args()
    pythia_delphes_dir=args.pythia_delphes_dir
    experiment=args.experiment
    rootfile=args.rootfile
    jaf_output_dir=args.jaf_output_dir
    com=args.CM_energy
    filemap_dict=filemap_from_single_diroctory(pythia_delphes_dir,rootfile,experiment)
    runlim(jaf_output_dir,filemap_dict,com,ss5b=1,os5b=1,lp5b=1,alphat7bb=1)
    print_and_save_CLs(jaf_output_dir,filemap_dict['internal_name'])


