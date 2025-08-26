# MutationTree
A BEAST 2 package for counting the number of synonymous, non-synonymous, and indel mutations on each branch



## Install

Package is currently not relased.

To build:

Download the repository and use ant:

```
mkdir ~/.beast/2.7/mutationtree/
cd MutationTree/
ant package
cp build/dist/mutationtree.package.v*.zip ~/.beast/2.7/mutationtree/tmp.zip
cd ~/.beast/2.7/mutationtree/
unzip -o tmp.zip
```


## Available counters

Each counter requires a ```BranchMutationSampler```, which will stochastically sample the mutations at the time of logging. This ensures that the various mutation summarisers below will be in harmony.

### SubstitutionSum
Counts the total number of substitutions per branch.

```<sampler spec="mutationtree.logger.mut.SubstitutionSum" sampler="@mutationsampler" />```

### SynonymousSubstSum
Counts the total number of synonymous substitutions per branch. Requirements: **nucleotide** data. Options: code (default: universal); readingFrame (default: 1).

```<sampler spec="mutationtree.logger.mut.SynonymousSubstSum" sampler="@mutationsampler" code="universal" readingFrame="1"/>```

### NonSynonymousSubstSum
Counts the total number of non-synonymous substitutions per branch. Requirements: **nucleotide** data. Options: code (default: universal); readingFrame (default: 1).

```<sampler spec="mutationtree.logger.mut.NonSynonymousSubstSum" sampler="@mutationsampler"  code="universal" readingFrame="1"/>```

### FromToSubstSum
Counts the total number of substitutions per branch from one of the characters in state X to one of the characters in state Y (user defined). Options: from (e.g. 'A' or 'AG'); from (e.g. 'C' or 'ACG').

``` <sampler spec="mutationtree.logger.mut.FromToSubstSum" sampler="@mutationsampler" from="A" to="CG" />```

### NucleotideTransitionCounter
Counts the total number of transitions per branch (purine to purine or pyrimidine to pyrimidine). Requirements: **nucleotide** data.

``` <sampler spec="mutationtree.logger.mut.NucleotideTransitionCounter" sampler="@mutationsampler" />```

### NucleotideTransversionCounter
Counts the total number of substitutions per branch (purine to pyrimidine or vice versa). Requirements: **nucleotide** data.

``` <sampler spec="mutationtree.logger.mut.NucleotideTransitionCounter" sampler="@mutationsampler" />```





## Setting up counters in XML

At every log, the ancestral sequence of each internal node will be stochatsically sampled, and so will the mutations along each branch. These mutations are summarised by the following loggers.

Append the following loggers to the bottom of the XML file to count the number of substitutions along each branch in the tree logger. In this example, the trees will be logged with length `SubstitutionSum` but you can leave in default units by removing the `lengths` input.
```
  <logger id="treelog" spec="Logger" fileName="substitution.trees" logEvery="10000" mode="tree">
      <log id="SampledSubstTreeLogger" spec="mutationtree.logger.SampledSubstTreeLogger" lengths="@SubstitutionSum" tree="@tree">
          <sampler id="SubstitutionSum" spec="mutationtree.logger.mut.SubstitutionSum" sampler="@mutationsampler" />
          <sampler id="SynonymousSubstSum" spec="mutationtree.logger.mut.SynonymousSubstSum" sampler="@mutationsampler" code="universal" readingFrame="1"/>
          <sampler id="NonSynonymousSubstSum" spec="mutationtree.logger.mut.NonSynonymousSubstSum" sampler="@mutationsampler"  code="universal" readingFrame="1"/>
          <sampler id="FromToSubstSum" spec="mutationtree.logger.mut.FromToSubstSum" sampler="@mutationsampler" from="A" to="A" />
          <sampler id="NucleotideTransitionCounter" spec="mutationtree.logger.mut.NucleotideTransitionCounter" sampler="@mutationsampler" />

          <sampler id="NucleotideTransversionCounter" spec="mutationtree.logger.mut.NucleotideTransversionCounter">
            <sampler id="mutationsampler" spec="mutationtree.evolution.BranchMutationSampler" tag="seq" useAmbiguities="true" >
               <tree idref="tree" />
               <siteModel idref="siteModelID" />
               <branchRateModel idref="clockModelID" />
               <data spec="mutationtree.evolution.PatternlessAlignment" data="@data" />
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

See the mutationtree.xml file in examples.
