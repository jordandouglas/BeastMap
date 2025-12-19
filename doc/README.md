
# Configuring BeastMap XML files

While BEAUti offers many options for setting up an analysis, there are still some features that can only be done by editing XML files. This page is intended for advanced users.



## Subsitution counters

Each counter shares a common ```BranchMutationSampler``` for any given tree likelihood, which will stochastically sample the mutations at the time of logging. This ensures that the various mutation summarisers below will be in harmony.



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



## Changing the tree lengths

By default, all branch lengths are in the default time units. However using the `lengths` parameter of a `SampledSubstTreeLogger` class, you can set the branch length to any of the other terms being mapped (e.g., number of substitutions, number of synonymous substitutions, etc.)
In this example, the trees will be logged with length `SubstitutionSum` but you can leave in default units by removing the `lengths` input.

```
 <logger id="treelog" spec="Logger" fileName="substitution.trees" logEvery="10000" mode="tree">
      <log id="SampledSubstTreeLogger" spec="beastmap.logger.SampledSubstTreeLogger" lengths="@SubstitutionSum" tree="@tree">
      ...
      </log>
</logger>
```

## Setting up counters in XML

At every log, the ancestral sequence of each internal node will be stochatsically sampled, and so will the mutations along each branch. These mutations are summarised by the following loggers.

Append the following loggers to the bottom of the XML file to count the number of substitutions along each branch in the tree logger. 
```
  <logger id="treelog" spec="Logger" fileName="substitution.trees" logEvery="10000" mode="tree">
      <log id="SampledSubstTreeLogger" spec="beastmap.logger.SampledSubstTreeLogger" tree="@tree">
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


## Indels and missing data

A `BranchMutationSampler` has an option called `indel`, which accepts a second `BranchMutationSampler`. These two stochastic mappers must share the same taxa and the latter must encode binary data. If this input is provided, then the missing sites (i.e., the deletions / gaps) will be masked out after ancestral sequence reconstruction, and the substitution counts will be smaller. If you do not use an indel mode, the gaps will be imputed with nucleotides/amino acids/etc, and the reported number of substitutions will be greater.


```
  <sampler id="BeastMap.StochasticMapper.data" spec="beastmap.evolution.BranchMutationSampler" burnin="100000" likelihood="@treeLikelihood.data" indel="@BeastMap.StochasticMapper.data_simpleindel">
```


## Filters

Any of the counters can be filtered to a certain range of sites. In the exampkle below, we will only count the substitutions in sites 1-10. 

```<sampler spec="beastmap.logger.mut.SubstitutionSum" sampler="@mutationsampler" filter="1-10" />```




## Weaving partitions together (and codon partition models)

You can join multiple partitions together using the `StochasticMapJoiner` class. In the example below, the three codon partitions are joined together using the `join="intersperse"` option. THis requires the three datasets are all the same length, and will then produce a single stochastic mapper from an alignment like so: `codonpos1_site1, codonpos2_site1, codonpos3_site1, codonpos1_site2, codonpos2_site2, codonpos3_site2, ...`. Alternatively the `join="concat"`option will simply concatenante them end-to-end.

```
 <log id="gluedPartitionSampler" spec="beastmap.evolution.StochasticMapJoiner" join="intersperse" >

      <plate var="c" range="1,2,3">
          <sampler id="mutationsampler.$(c)" spec="beastmap.evolution.BranchMutationSampler" >
               <tree idref="tree" />
               <siteModel idref="siteModel.$(c)" />
               <branchRateModel idref="clockModel" />
               <data spec="beastmap.evolution.PatternlessAlignment" data="@data.$(c)" />
          </sampler>
      </plate>

  </log>
```

This `StochasticMapJoiner` can then be treated the same as any `BranchMutationSampler`, for example:


```
  <sampler id="SynonymousSubstSum" spec="beastmap.logger.mut.SynonymousSubstSum" sampler="@mutationsampler" code="universal" readingFrame="1"/>
```

See the [simulation/Codons](/simulation/Codons) folder for some examples.
