from __future__ import print_function
import os, sys
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import base.process
from base.process import process


def collect_args():
    """Returns a tb-specific argument parser.
    """
    parser = base.process.collect_args()
    parser.set_defaults(
        json="prepared/tb.json"
    )
    return parser


config = {
    "dir": "tb",
    "in": "prepared/tb.json",
    "newick_tree_options": {"nthreads": 4},
    "clock_filter": {
        "n_iqd": 4,
    },
    "geo_inference": ['country', 'region'], # what traits to perform this on
    "geo_inference_options": {
        "root_state": {
            "region": "southeast_asia",
            "country": "vietnam",
        },
    },
    "auspice": { ## settings for auspice JSON export
        "color_options": {
            "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
            "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete"},
        },
        "controls": {'authors':['authors']},
        "defaults": {'mapTriplicate': True},
    },
    "timetree_options": {
        "Tc": 'opt',
        "confidence":True
    }
}

if __name__=="__main__":
    parser = collect_args()
    params = parser.parse_args()

    if params.clean:
        config["clean"] = True

    if params.json:
        config["in"] = params.json

    config["newick_tree_options"]["raxml"] = not params.no_raxml

    runner = process(config)
#    print("aligning")
    runner.align(fill_gaps=True)
#    print("building tree")
#    runner.build_tree()
#    print("timetree")
#    runner.timetree_setup_filter_run()
#    print("geo_inference")
#    runner.run_geo_inference()
#    print("saving nexus")
#    runner.save_as_nexus()
#    print("auspice_export")
#    runner.auspice_export()
