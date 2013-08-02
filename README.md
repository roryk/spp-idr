### spp-idr: a wrapper around peak callers and IDR for performing peak analyses
This repo probably should get renamed, as it has expanded past using just SPP.

Two peaks callers can be used:

- SPP for performing peak calls on punctate chromatin marks for CHIP-seq experiments.
- clipper for performing peak calls on RNA binding experiments such as CLIP-seq or
IMPACT-seq, for reads mapped to hg19 only

#### Example situation
You have BAM files from a CHIP-seq experiment made with punctate chromatin marks such as
H3k9ac, H3k27ac, H3k4me3, etc. and want to call a set of peaks reproducible across your
replicates.

This is a simplified implementation of the workflow sketched out here:
https://sites.google.com/site/anshulkundaje/projects/idr

Anshul explains the purpose of IDR well:

>Reproducibility is essential to reliable scientiﬁc discovery in high-throughput experiments. The IDR (Irreproducible Discovery Rate) framework is a uniﬁed approach to measure the reproducibility of ﬁndings identiﬁed from replicate experiments and provide highly stable thresholds based on reproducibility. Unlike the usual scalar measures of reproducibility, the IDR approach creates a curve, which quantitatively assesses when the ﬁndings are no longer consistent across replicates. In layman's terms, the IDR method compares a pair of ranked lists of identifications (such as ChIP-seq peaks). These ranked lists should not be pre-thresholded i.e. they should provide identifications across the entire spectrum of high confidence/enrichment (signal) and low confidence/enrichment (noise). The IDR method then fits the bivariate rank distributions over the replicates in order to separate signal from noise based on a defined confidence of rank consistency and reproducibility of identifications i.e the IDR threshold.

An intuitive overview of IDR can be found here:
http://www.personal.psu.edu/users/q/u/qul12/IDR101.pdf

#### Installation
```
git clone git@github.com:roryk/spp-idr.git
cd spp-idr
python setup.py install
python install-tools.py --tool_path=tools
```
#### Usage
##### Run locally
```
run_idr.py --caller spp --controls BAM1 BAM2 BAM3 --experimental BAM1 BAM2 BAM3 tool_path
```

##### Run on Platform LSF on queue 'default'
You should set num-jobs to be twice the amount of BAM files you have.
```
run_idr.py --num-jobs 8 --lsf-queue default --caller spp --control BAM1 BAM2 BAM3 --experimental BAM1 BAM BAM3 tool_path
```

##### Run on Sun Grid Engine (SGE) on queue 'default'
You should set num-jobs to be twice the amount of BAM files you have.
```
run_idr.py --num-jobs 8 --sge-queue default --caller spp --control BAM1 BAM2 BAM3 --experimental BAM1 BAM2 BAM3 tool_path
```

##### Run on Torque on queue 'default'
You should set num-jobs to be twice the amount of BAM files you have.
```
run_idr.py --caller spp --num-jobs 8 --torque-queue torque --control BAM1 BAM2 BAM3 --experimental BAM1 BAM2 BAM3 tool_path
```

#### Output
This generates quite a lot of output files. The conservative set of peak calls
and the optimum set of peak calls can be found in the peaks directory.

##### SPP specific important information
There are also several diagnostic plots in the "peaks" directory created by SPP
which are important to look at to determine if the correct fragment length was
calculated. You should look at those and make sure they look correct, otherwise
you might end up analyzing a lot of noisy garbage. The documentation for SPP has
more information about these plots.
