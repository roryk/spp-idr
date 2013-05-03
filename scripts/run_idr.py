# script to run analysis
from argparse import ArgumentParser
from cluster_helper import cluster_view
from spp_idr import idr
import os


def main(control, experimental, tool_path, view=None):
    if not view:
        mapper = map
    spp_path = os.path.join(tool_path, "phantompeakqualtools", "run_spp.R")
    idr_runner_path = os.path.join(tool_path, "idrCode", "batch-consistency-analysis.r")
    idr_plotter_path = os.path.join(tool_path, "idrCode", "batch-consistency-plot.r")
    peak_caller = idr.SPPPeakCaller(spp_path)
    peak_caller.set_map_function(view.map)
    peaks = idr.run_analysis(control, experimental, spp_path, idr_runner_path,
                             idr_plotter_path, mapper)


if __name__ == "__main__":
    parser = ArgumentParser(description='run IDR analysis on a set of chip-seq calls')
    parser.add_argument('--controls', nargs="+", help="Control BAM files")
    parser.add_argument('--experimental', nargs="+", help="Experimental BAM files")
    parser.add_argument('--lsf-queue', help="LSF queue name")
    parser.add_argument('--sge-queue', help="SGE queue name")
    parser.add_argument('--torque-queue', help="Torque queue name")
    parser.add_argument('--num-jobs', default=1, help="number of parallel jobs to run")
    parser.add_argument('tool_path', help="Path to spp and idr installation.")
    args = parser.parse_args()

    if args.lsf_queue:
        with cluster_view("LSF", args.lsf_queue, args.num_jobs) as view:
            main(args.control, args.experimental, args.spp_path, view)
    elif args.sge_queue:
        with cluster_view("LSF", args.lsf_queue, args.num_jobs) as view:
            main(args.control, args.experimental, args.spp_path, view)
    elif args.torque_queue:
        with cluster_view("LSF", args.lsf_queue, args.num_jobs) as view:
            main(args.control, args.experimental, args.spp_path, view)
    else:
        main(args.control, args.experimental, args.spp_path, None)
