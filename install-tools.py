from argparse import ArgumentParser
from spp_idr import idr
import subprocess
import os

PHANTOMPEAKQUALTOOLS_URL = "http://phantompeakqualtools.googlecode.com/files/ccQualityControl.v.1.1.tar.gz"
IDR_URL = "https://sites.google.com/site/anshulkundaje/projects/idr/idrCode.tar.gz"
FIXED_SPP_URL = "https://dl.dropboxusercontent.com/u/2822886/spp/spp_1.10.1.fixed.tar.gz"

if __name__ == "__main__":
    parser = ArgumentParser(description='Download and install SPP and IDR scripts for use '
                            'with spp-idr')
    parser.add_argument('--tool_path', default="tools",
                        help="path to install tools to")
    args = parser.parse_args()
    idr.safe_makedir(args.tool_path)
    print "Installing phantompeakqualtools."
    idr.download_to_dir(PHANTOMPEAKQUALTOOLS_URL, args.tool_path)
    print "Installing IDR."
    idr.download_to_dir(IDR_URL, args.tool_path)
    print "Installing fixed SPP."
    spp_file = idr.download_to_dir(FIXED_SPP_URL, os.path.join(args.tool_path,
                                                            "phantompeakqualtools"))
    cl = ["R", "CMD", "INSTALL", spp_file]
    cur_dir = os.getcwd()
    os.chdir(os.path.join(args.tool_path, "phantompeakqualtools"))
    subprocess.check_call(cl)
    os.chdir(cur_dir)
    print "Installation of tools complete."
