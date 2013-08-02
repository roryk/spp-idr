# script to run analysis
from argparse import ArgumentParser
from cluster_helper.cluster import cluster_view
from spp_idr import idr
import os

def main(control, experimental, tool_path, mapper, caller, cores_per_job):
    spp_path = os.path.join(tool_path, "phantompeakqualtools", "run_spp.R")
    idr_runner_path = os.path.join(tool_path, "idrCode", "batch-consistency-analysis.r")
    idr_plotter_path = os.path.join(tool_path, "idrCode", "batch-consistency-plot.r")
    peaks = idr.run_analysis(control, experimental, spp_path, idr_runner_path,
                             idr_plotter_path, mapper, caller, cores_per_job)


if __name__ == "__main__":
    parser = ArgumentParser(description='run IDR analysis on a set of chip-seq calls')
    parser.add_argument('--control', nargs="+", help="Control BAM files")
    parser.add_argument('--experimental', nargs="+", help="Experimental BAM files")
    parser.add_argument('--lsf-queue', help="LSF queue name")
    parser.add_argument('--sge-queue', help="SGE queue name")
    parser.add_argument('--torque-queue', help="Torque queue name")
    parser.add_argument('--num-jobs', default=1, help="number of parallel jobs to run",
                        type=int)
    parser.add_argument('--caller', default="spp", help="Peak caller to run "
                        "(spp or clipper)")
    parser.add_argument('--cores-per-job', default=1, type=int,
                        help="Number of cores to run for each job.")
    parser.add_argument('tool_path', help="Path to spp and idr installation.")
    args = parser.parse_args()

    args.control = map(os.path.abspath, args.control)
    args.experimental = map(os.path.abspath, args.experimental)

    if args.lsf_queue:
        with cluster_view("LSF", args.lsf_queue, args.num_jobs,
                          cores_per_job=args.cores_per_job) as view:
            main(args.control, args.experimental, args.tool_path, view.map,
                 args.caller, args.cores_per_job)
    elif args.sge_queue:
        with cluster_view("LSF", args.lsf_queue, args.num_jobs,
                          cores_per_job=args.cores_per_job) as view:
            main(args.control, args.experimental, args.tool_path, view.map,
                 args.caller, args.cores_per_job)
    elif args.torque_queue:
        with cluster_view("LSF", args.lsf_queue, args.num_jobs,
                          cores_per_job=args.cores_per_job) as view:
            main(args.control, args.experimental, args.tool_path, view.map,
                 args.caller, args.cores_per_job)
    else:
        main(args.control, args.experimental, args.tool_path, map,
             args.caller, args.cores_per_job)
