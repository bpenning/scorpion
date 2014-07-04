#! /usr/bin/env python
import optparse
from python.make_filemap_dict import filemap_from_pythia6_single_diroctory , filemap_from_pythia8_single_diroctory
from python.run_limit_code import runlim
from python.extract_CLs import print_and_save_CLs 

pythia_versions={
        '6':{'filemap_function':filemap_from_pythia6_single_diroctory},
        '8':{'filemap_function':filemap_from_pythia8_single_diroctory},
        }

def parse_args():
    parser=optparse.OptionParser()
    parser.add_option('--pythia-delphes-dir',help='output dir of pythia-delphes')
    parser.add_option('--jaf-output-dir')
    parser.add_option('--with-cross-section',type=float,help='give xsec (in barns)')
    parser.add_option('--CM-energy',default=7,type=int)
    parser.add_option('--experiment',default='CMS7')
    parser.add_option('--rootfile',default='delphes-output.root')
    parser.add_option('--pythia-version',choices=pythia_versions.keys(),default='8')
    parser.add_option('--gen-info',action="store_true")
    options, args=parser.parse_args()
    return options

if __name__=="__main__":
    args=parse_args()
    pythia_delphes_dir=args.pythia_delphes_dir
    experiment=args.experiment
    rootfile=args.rootfile
    jaf_output_dir=args.jaf_output_dir
    gen_info = args.gen_info
    com=args.CM_energy
    filemap_from_single_diroctory=pythia_versions[args.pythia_version]['filemap_function']
    filemap_dict=filemap_from_single_diroctory(pythia_delphes_dir,rootfile,experiment)
    if args.with_cross_section:
        filemap_dict['xsec']=args.with_cross_section
    runlim(jaf_output_dir,filemap_dict,com,gen_info,ss5b=1,os5b=1,lp5b=1,alphat7bb=1)
    #print_and_save_CLs(jaf_output_dir,filemap_dict['internal_name'])
    print_and_save_CLs(jaf_output_dir,"delphes-output_CMS8_MonoJet_analysis8",filemap_dict['xsec'])
    


