# script to run analysis
from argparse import ArgumentParser
from cluster_helper import cluster_view
from spp_idr import idr


def main(control, experimental, spp_path, view=None):
    view = unblock_view(view)
    peak_caller = idr.SPPPeakCaller(spp_path)
    peak_caller.set_map_function(view.map)
    peaks = idr.call_peaks(control, experimental, peak_caller)

def unblock_view(view):
    if view:
        view.block = False
    return view

def block_view(view):
    if view:
        view.block = True
    return view

if __name__ == "__main__":
    parser = ArgumentParser(description='run IDR analysis on a set of chip-seq calls')
    parser.add_argument('--controls', nargs="+", help="Control BAM files")
    parser.add_argument('--experimental', nargs="+", help="Experimental BAM files")
    parser.add_argument('--lsf-queue', help="LSF queue name")
    parser.add_argument('--num-jobs', default=1, help="number of parallel jobs to run")
    parser.add_argument('spp_path', help="Path to spp installation.")
    args = parser.parse_args()

    if args.lsf_queue:
        with cluster_view("LSF", args.lsf_queue, args.num_jobs) as view:
            main(args.control, args.experimental, args.spp_path, view)
