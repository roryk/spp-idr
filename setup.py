#!/usr/bin/env python
from setuptools import setup, find_packages

setup(name = "spp_idr",
      version = "0.0.1",
      author = "Rory Kirchner",
      author_email = "rory.kirchner@gmail.com",
      description = "Call peaks with SPP and do IDR analysis.",
      license = "MIT",
      url = "https://github.com/roryk/idr-wrapper/",
      namespace_packages = ["spp_idr"],
      packages = find_packages(),
      scripts = ['scripts/run_idr.py'],
      install_requires = ["pysam >= 0.4.1",
                          "nose >= 1.3.0",
                          "ipython-cluster-helper"])
