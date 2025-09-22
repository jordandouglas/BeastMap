# BeastMap
A BEAST 2 package for counting the number of synonymous, non-synonymous, and indel mutations on each branch. The method first performs ancestral sequence reconstruction on the internal nodes, and then uses stochastic mapping to sample a mutation pathway along each branch. This all happens during MCMC, and the package is compatible with a wide range of existing BEAST 2 site, clock, and tree models.

Warning: package is currently in pre-pre-release. It has passed the simulation studies, but the code is still quite volatile.

## Install

Package is currently not released.

To build:

Download the repository and use ant:

```
mkdir ~/.beast/2.7/beastmap/
cd BeastMap/
ant package
cp build/dist/beastmap.package.v*.zip ~/.beast/2.7/beastmap/tmp.zip
cd ~/.beast/2.7/beastmap/
unzip -o tmp.zip
```


Or you can download the zip file directly from the releases section. Please make sure to install the `CodonSubstModels` package too, as that is a dependency.


## Available counters

Each counter requires a ```BranchMutationSampler```, which will stochastically sample the mutations at the time of logging. This ensures that the various mutation summarisers below will be in harmony.

### SubstitutionSum
Counts the total number of substitutions per branch.

```<sampler spec="beastmap.logger.mut.SubstitutionSum" sampler="@mutationsampler" />```

### SynonymousSubstSum
Counts the total number of synonymous substitutions per branch. Requirements: **nucleotide** or **codon** data. Options: code (default: universal); readingFrame (default: 1).

```<sampler spec="beastmap.logger.mut.SynonymousSubstSum" sampler="@mutationsampler" code="universal" readingFrame="1"/>```

### NonSynonymousSubstSum
Counts the total number of non-synonymous substitutions per branch. Requirements: **nucleotide** or **codon** data. Options: code (default: universal); readingFrame (default: 1).

```<sampler spec="beastmap.logger.mut.NonSynonymousSubstSum" sampler="@mutationsampler"  code="universal" readingFrame="1"/>```

### FromToSubstSum
Counts the total number of substitutions per branch from one of the characters in state X to one of the characters in state Y (user defined). Options: from (e.g. 'A' or 'AG'); from (e.g. 'C' or 'ACG').

``` <sampler spec="beastmap.logger.mut.FromToSubstSum" sampler="@mutationsampler" from="A" to="CG" />```

### NucleotideTransitionCounter
Counts the total number of transitions per branch (purine to purine or pyrimidine to pyrimidine). Requirements: **nucleotide** data.

``` <sampler spec="beastmap.logger.mut.NucleotideTransitionCounter" sampler="@mutationsampler" />```

### NucleotideTransversionCounter
Counts the total number of substitutions per branch (purine to pyrimidine or vice versa). Requirements: **nucleotide** data.

``` <sampler spec="beastmap.logger.mut.NucleotideTransitionCounter" sampler="@mutationsampler" />```


### AminoAcidClassChanges
Counts the total number of amino acid substitutions per branch such that the amino acid chnages to a different functional class. Requirements: **amino acid** data. These classes are based on the BLOSUM62 matrix non-polar: {IVLM}, amide/amine: {DENQ}, basic: {HKR}, aromatic: {FWY}, small/polar: {AST}, cysteine: {C}, glycine: {G}, proline: {P}.

``` <sampler spec="beastmap.logger.mut.AminoAcidClassChanges" sampler="@mutationsampler" />```


### AminoAcidClassRemains
Counts the total number of amino acid substitutions per branch such that the amino acid remains in the same functional class. Requirements: **amino acid** data.

``` <sampler spec="beastmap.logger.mut.AminoAcidClassRemains" sampler="@mutationsampler" />```


### SubstitutionSummer

Takes one of the other per-branch counters and adds all the numbers together across the whole tree

``` <log spec="beastmap.logger.SubstitutionSummer" counter="@ID_OF_COUNTER"  />```


## Setting up counters in XML

At every log, the ancestral sequence of each internal node will be stochatsically sampled, and so will the mutations along each branch. These mutations are summarised by the following loggers.

Append the following loggers to the bottom of the XML file to count the number of substitutions along each branch in the tree logger. In this example, the trees will be logged with length `SubstitutionSum` but you can leave in default units by removing the `lengths` input.
```
  <logger id="treelog" spec="Logger" fileName="substitution.trees" logEvery="10000" mode="tree">
      <log id="SampledSubstTreeLogger" spec="beastmap.logger.SampledSubstTreeLogger" lengths="@SubstitutionSum" tree="@tree">
          <sampler id="SubstitutionSum" spec="beastmap.logger.mut.SubstitutionSum" sampler="@mutationsampler" />
          <sampler id="SynonymousSubstSum" spec="beastmap.logger.mut.SynonymousSubstSum" sampler="@mutationsampler" code="universal" readingFrame="1"/>
          <sampler id="NonSynonymousSubstSum" spec="beastmap.logger.mut.NonSynonymousSubstSum" sampler="@mutationsampler"  code="universal" readingFrame="1"/>
          <sampler id="FromToSubstSum" spec="beastmap.logger.mut.FromToSubstSum" sampler="@mutationsampler" from="A" to="G" />
          <sampler id="NucleotideTransitionCounter" spec="beastmap.logger.mut.NucleotideTransitionCounter" sampler="@mutationsampler" />
          <sampler id="NucleotideTransversionCounter" spec="beastmap.logger.mut.NucleotideTransversionCounter" sampler="@mutationsampler" />
          <sampler id="SubstitutionSumFiltered" spec="beastmap.logger.mut.SubstitutionSum" sampler="@mutationsampler" filter="2,1-99\3" />

          <!-- This will log the ancestral sequences onto the tree -- it will make the tree files quite large so turn it off if you dont want it -->
          <sampler id="AncestralSequenceLogger" spec="beastmap.logger.AncestralSequenceLogger" sampler="@mutationsampler"/>

          <sampler id="NucleotideTransversionCounter" spec="beastmap.logger.mut.NucleotideTransversionCounter">
            <sampler id="mutationsampler" spec="beastmap.evolution.BranchMutationSampler" tag="seq" useAmbiguities="true" substModelIsNodeDependent="false" burnin="50000" >
               <tree idref="tree" />
               <siteModel idref="siteModelID" />
               <branchRateModel idref="clockModelID" />
               <data spec="beastmap.evolution.PatternlessAlignment" data="@data" />
            </sampler>
          </sampler>
      </log>
  </logger>
```

You can also place these loggers in the trace file as well as the tree file to estimate their ESS:

```
<logger id="tracelog" spec="Logger" fileName="$(filebase).log" logEvery="10000" model="@posterior" sanitiseHeaders="true" sort="smart">
   <sampler idref="NucleotideTransitionCounter" />
</logger>
```



## Examples

See the xml file in examples.
