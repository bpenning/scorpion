#! /usr/bin/env python
import optparse
from python.make_filemap_dict import filepair_dict_from_pythia6_single_directory
from python.make_filemap_dict import filepair_dict_from_pythia8_single_directory
from python.run_limit_code import runlim
from python.stop_xsections import get_stop_x_section_from_slha_file
from python.stop_xsections import get_sbot_x_section_from_slha_file
from python.gluino_xsections import get_gluino_x_section_from_slha_file
from python.extract_CLs import print_and_save_CLs
analyses_dict = {
        'all-7tev': ['ss5b', 'os5b', 'lp5b', 'alphat7bb'],
        'os5b-only': ['os5b'],
        'lp20b-only': ['lp20b'],
        'monojet20b-only': ['monojet20b'],
        'dmbsr1-only': ['dmbsr1'],
        'ge3lp20b-only': ['ge3lp20b'],
        'mt220b-only': ['mt220b'],
        'ss820b-only': ['ss820b'],
        'alphat-only': ['alphat7bb'],
        'alphat8-only': ['alphat20b'],
        'alphat13-only': ['alphat13T'],
        'alphat8valid-only': ['alphat20bvalid'],
        'all-8tev': ['monojet20b', 'mt220b', 'lp20b', 'ss820b', 'os5b',
            'ge3lp20b'],
        'all-dm-8tev': ['monojet20b', 'dmbsr1', 'alphat20b'],
        '8tev-only': ['monojet20b', 'mt220b', 'lp20b', 'ss820b',
            'ge3lp20b'],
        'hinv-8tev': ['hinv20b'],
#        'all-8tev': ['monojet20b', 'lp20b', 'ss820b', 'os5b',],
        }

def vararg_callback(option, opt_str, value, parser):
    assert value is None
    value = []

    def floatable(str):
        try:
            float(str)
            return True
        except ValueError:
            return False

    for arg in parser.rargs:
        # stop on --foo like options
        if arg[:2] == "--" and len(arg) > 2:
            break
        # stop on -a, but not on -3 or -3.0
        if arg[:1] == "-" and len(arg) > 1 and not floatable(arg):
            break
        value.append(arg)

    del parser.rargs[:len(value)]
    setattr(parser.values, option.dest, value)

def parse_args():
    """
    Parse options from command line
    """
    parser = optparse.OptionParser()
    parser.add_option('--pythia-delphes-dirs', action='callback',
            callback=vararg_callback, help='output dir of pythia-delphes',
            dest='pythia_delphes_dirs')
    parser.add_option('--experiments', action='callback', dest='experiments',
            callback=vararg_callback)
    parser.add_option('--analyses', choices=analyses_dict.keys(),
            default='all-7tev', help='specify which analeses to run')
    parser.add_option('--gen-info', action='store_true')
    parser.add_option('--use-event-weights', action='store_true')
    parser.add_option('--expected-limits', action='store_true')
    parser.add_option('--jaf-output-dir')
    parser.add_option('--with-cross-section', type=float,
            help='give xsec (in barns)')
    parser.add_option('--with-cms-stop-cross-section',
            help='provide slha file from which to extract mstop')
    parser.add_option('--with-cms-sbot-cross-section',
            help='provide slha file from which to extract msbot')
    parser.add_option('--with-cms-gluino-cross-section',
            help='provide slha file from which to extract mgluino')
    parser.add_option('--xsec-factor', type=float, default=1.0)
    parser.add_option('--rootfile', default='delphes-output.root')
    parser.add_option('--pythia-version', default='6')
    parser.add_option('--delphes-version',choices=('2','3'), default='3')
    parser.add_option('--skip-limits', action='store_true')
    options, args = parser.parse_args()
    return options

def main(pythia_delphes_dirs, jaf_output_dir, with_cross_section,
        with_cms_stop_cross_section, with_cms_sbot_cross_section,
        with_cms_gluino_cross_section, xsec_factor, experiments, 
        rootfile, analyses, gen_info, pythia_version, 
        use_event_weights, delphes_version, expected_limits, skip_limits):

    """
    Main program
    """
    delphes_int = -1
    if (delphes_version == '2'):
	delphes_int = 0
    elif (delphes_version == '3'):
	delphes_int = 1
    else:
	print "Invalid delphes option"
	exit()

    if not len(pythia_delphes_dirs) == len(experiments):
        print('ERROR: the number of directories does not equal the number of '
                'experiments')
    if pythia_version == '6':
        filepair_dict_from_single_directory = filepair_dict_from_pythia6_single_directory
    elif pythia_version == '8':
        filepair_dict_from_single_directory = filepair_dict_from_pythia8_single_directory
    filemap_dict = {}
    for i in range(len(experiments)):
        pythia_delphes_dir = pythia_delphes_dirs[i]
        experiment = experiments[i]
        filepair_dict = filepair_dict_from_single_directory(pythia_delphes_dir,
                rootfile,with_cross_section)
        if with_cms_stop_cross_section:
            slhafile = with_cms_stop_cross_section
            filepair_dict['xsec'] = get_stop_x_section_from_slha_file(slhafile)
        elif with_cms_sbot_cross_section:
            slhafile = with_cms_sbot_cross_section
            filepair_dict['xsec'] = get_sbot_x_section_from_slha_file(slhafile)
        elif with_cms_gluino_cross_section:
            slhafile = with_cms_gluino_cross_section
            filepair_dict['xsec'] = get_gluino_x_section_from_slha_file(slhafile)
        analyses_kwargs = {}
        for name in analyses_dict[analyses]:
            analyses_kwargs[name] = True
        # if with_cross_section:
        #     filepair_dict['xsec'] = with_cross_section
        filepair_dict['xsec'] = xsec_factor*filepair_dict['xsec']
        filemap_dict[experiment] = filepair_dict

    runlim(jaf_output_dir, filemap_dict, gen_info,delphes_int, use_event_weights, expected_limits, skip_limits, **analyses_kwargs)
    #print_and_save_CLs(jaf_output_dir)
    #jaf_CMS8_alphaT20b_analysis20
    

if __name__ == "__main__":
    main(**vars(parse_args()))
