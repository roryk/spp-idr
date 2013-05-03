# spp-idr: a wrapper around SPP and IDR to perform CHIP-seq analysis
This is a simplified implementation of the workflow sketched out here:
https://sites.google.com/site/anshulkundaje/projects/idr

Anshul explains the purpose of IDR well:
```
Reproducibility is essential to reliable scientiﬁc discovery in high-throughput
experiments. The IDR (Irreproducible Discovery Rate) framework is a uniﬁed
approach to measure the reproducibility of ﬁndings identiﬁed from replicate
experiments and provide highly stable thresholds based on
reproducibility. Unlike the usual scalar measures of reproducibility, the IDR
approach creates a curve, which quantitatively assesses when the ﬁndings are no
longer consistent across replicates. In layman's terms, the IDR method compares
a pair of ranked lists of identifications (such as ChIP-seq peaks). These ranked
lists should not be pre-thresholded i.e. they should provide identifications
across the entire spectrum of high confidence/enrichment (signal) and low
confidence/enrichment (noise). The IDR method then fits the bivariate rank
distributions over the replicates in order to separate signal from noise based
on a defined confidence of rank consistency and reproducibility of
identifications i.e the IDR threshold.
```

# Installation
```
git clone git@github.com:roryk/spp-idr.git
cd spp-idr
python setup.py install
python install-tools.py --tool_path=tools
```
# Usage
## Run locally
'''
run_idr.py --controls BAM1, BAM2, BAM3 --experimental BAM1, BAM2, BAM3 tool_path
'''

## Run on Platform LSF on queue 'default'
'''
run_idr.py --lsf-queue default --controls BAM1, BAM2, BAM3 --experimental BAM1, BAM2, BAM3 tool_path
'''

## Run on Platform LSF on queue 'default'
'''
run_idr.py --lsf-queue default --controls BAM1, BAM2, BAM3 --experimental BAM1, BAM2, BAM3 tool_path
'''

## Run on Sun Grid Engine (SGE) on queue 'default'
'''
run_idr.py --sge-queue default --controls BAM1, BAM2, BAM3 --experimental BAM1, BAM2, BAM3 tool_path
'''

## Run on Torque on queue 'default'
'''
run_idr.py --torque-queue torque --controls BAM1, BAM2, BAM3 --experimental BAM1, BAM2, BAM3 tool_path
'''