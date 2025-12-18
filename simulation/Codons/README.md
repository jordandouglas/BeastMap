
# Testing stochastic mapping (nucleotide) on codon sequences using simulation studies 


## Dependencies

- BEAST 2.7 - assumed to be installed in the `~/beast/bin/` directory
- BeastMap package installed
- CodonSubstModel package installed 
- GammaSpikeModel package installed for tree priors
- R v4, available from command line as `Rscript` 


## Instructions



1. Simulate 100 birth-death trees and then simulate 100 datasets under these trees. The substitution histories of these trees will be recorded.  


```bash prepare.sh```





2. Then to run MCMC on a single replicate:

```
cd templates
cd rep1
~/beast/bin/beast -df var.seq.json ../../run.xml
```

3. Or, to run MCMC on all replicates one at a time:

```
bash run.sh
```



4. Or, you can run MCMC on all replicates in parallel across multiple threads (default 4). **Warning**: these jobs will run in the background, so do not run this script unless you are prepared to monitor/kill the jobs


```
bash runParallel.sh 
```


## XML files

- `truth.xml` - simulating the tree and other parameters. By default, there are 30 taxa.

- `seqsim.xml` - simulates sequences down the tree, which is specified in json format and parsed into the xml file through the `-df` option of beast.

- `run.xml` - performs MCMC on a dataset specified by the `-df` option.


