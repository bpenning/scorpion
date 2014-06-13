#! /usr/bin/env python
import optparse
from python.make_filemap_dict import filemap_from_pythia6_single_diroctory 
from python.make_filemap_dict import filemap_from_pythia8_single_diroctory
from python.run_limit_code import runlim
from python.extract_CLs import print_and_save_CLs 
from python.stop_xsections import get_stop_x_section_from_slha_file
from python.stop_xsections import get_sbot_x_section_from_slha_file
from python.gluino_xsections import get_gluino_x_section_from_slha_file

pythia_versions={
        '6':{'filemap_function':filemap_from_pythia6_single_diroctory},
        '8':{'filemap_function':filemap_from_pythia8_single_diroctory},
        }

analyses={
        'all-7tev':['ss5b','os5b','lp5b','alphat7bb'],        
        'lp20b-only':['lp20b'],        
        'ge3lp20b-only':['ge3lp20b'],
        'mt220b-only':['mt220b'],
        'ss820b-only':['ss820b'],
        }

def parse_args():
    parser=optparse.OptionParser()
    parser.add_option('--pythia-delphes-dir',
            help='output dir of pythia-delphes')
    parser.add_option('--jaf-output-dir')
    parser.add_option('--with-cross-section',type=float,
            help='give xsec (in barns)')
    parser.add_option('--with-cms-stop-cross-section',
            help='provide slha file from which to extract mstop')
    parser.add_option('--with-cms-sbot-cross-section',
            help='provide slha file from which to extract msbot')
    parser.add_option('--with-cms-gluino-cross-section',
            help='provide slha file from which to extract mgluino')
    parser.add_option('--CM-energy',default=7,type=int)
    parser.add_option('--experiment',default='CMS7')
    parser.add_option('--rootfile',default='delphes-output.root')
    parser.add_option('--pythia-version',choices=pythia_versions.keys(),
            default='8')
    parser.add_option('--analyses',choices=analyses.keys(),default='all-7tev',
            help='specify which analeses to run')
    options, args=parser.parse_args()
    return options

if __name__=="__main__":
    args=parse_args()
    pythia_delphes_dir=args.pythia_delphes_dir
    experiment=args.experiment
    rootfile=args.rootfile
    jaf_output_dir=args.jaf_output_dir
    com=args.CM_energy
    filemap_from_single_diroctory=pythia_versions[args.pythia_version]['filemap_function']
    filemap_dict=filemap_from_single_diroctory(pythia_delphes_dir,rootfile,experiment)
    if args.with_cross_section:
        filemap_dict['xsec']=args.with_cross_section
    elif args.with_cms_stop_cross_section:
        slhafile=args.with_cms_stop_cross_section
        filemap_dict['xsec']=get_stop_x_section_from_slha_file(slhafile)
    elif args.with_cms_sbot_cross_section:
        slhafile=args.with_cms_sbot_cross_section
        filemap_dict['xsec']=get_sbot_x_section_from_slha_file(slhafile)
    elif args.with_cms_gluino_cross_section:
        slhafile=args.with_cms_gluino_cross_section
        filemap_dict['xsec']=get_gluino_x_section_from_slha_file(slhafile)
    analyses_kwargs={}
    for name in analyses[args.analyses]:
        analyses_kwargs[name]=True
    runlim(jaf_output_dir,filemap_dict,com,**analyses_kwargs)
#    print_and_save_CLs(jaf_output_dir,filemap_dict['internal_name'])


