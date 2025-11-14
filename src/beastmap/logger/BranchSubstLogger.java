package beastmap.logger;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.Nucleotide;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;
import beastmap.evolution.StochasticMapper;
import beastmap.util.Mutation;
import beastmap.util.MutationUtils;
import codonmodels.evolution.datatype.Codon;
import codonmodels.evolution.datatype.GeneticCode;
import starbeast3.evolution.speciation.GeneTreeForSpeciesTreeDistribution;


@Description("Counts the number and type of substitutions on each branch")
public abstract class BranchSubstLogger extends CalculationNode implements Loggable, Function, StochasticMapProperty {

	
	final public Input<StochasticMapper> samplerInput = new Input<>("sampler", "mutation sampler to log");
	final public Input<StochasticMapper> truthInput = new Input<>("truth", "mutation sampler to log", Validate.XOR, samplerInput);
	final public Input<Boolean> includeRootInput = new Input<>("includeRoot", "count the root too?", false);
	
	
	final public Input<GeneTreeForSpeciesTreeDistribution> geneInput = new Input<>("gene", "map branches to a species tree?");
	//final public Input<Tree> speciesTreeInput = new Input<>("speciesTree", "map branches to a species tree?");
	
	
	final public Input<String> filterInput = new Input<>("filter", "specifies which of the sites in the input alignment we will restrict to" +
	            "First site is 1." +
	            "Filter specs are comma separated, either a singleton, a range [from]-[to] or iteration [from]:[to]:[step]; " +
	            "1-100 defines a range, " +
	            "1-100\3 or 1:100:3 defines every third in range 1-100, " +
	            "1::3,2::3 removes every third site. " +
	            "Default for range [1]-[last site], default for iterator [1]:[last site]:[1]", Validate.OPTIONAL);
	

	 protected List<Integer> filter;
	 GeneTreeForSpeciesTreeDistribution geneTreePrior;
	 
	 
    
    @Override
    public void initAndValidate() {
    	
    	DataType dt = this.getDataType();
    	if (!this.canHandleDataType(dt)) {
    		throw new IllegalArgumentException(this.getID() + " cannot support datatype " + dt.getClass());
    	}
    	
    	this.filter = null;
    	if (filterInput.get() != null && !filterInput.get().isEmpty()) {
    		this.filter = MutationUtils.parseFilterSpec(getSiteAndPatternCount(), filterInput.get());
    	}
    	
    	
    	this.geneTreePrior = geneInput.get();
    	if (this.geneTreePrior != null) {
    		
    		
    		Tree geneTreeOther = samplerInput.get() == null ? (Tree) truthInput.get().getTree() : (Tree) samplerInput.get().getTree();
    		if (this.geneTreePrior.getGeneTree() != geneTreeOther) {
    			throw new IllegalArgumentException("mismatching gene trees " + this.geneTreePrior.getGeneTree().getID() + " != " + geneTreeOther);
    		}
    	}
    	
    }
    
    

    protected double getMutationSummary(List<Mutation> mutations, Node node){
    	
    	//Log.warning("getMutationSummary " + mutations.size());
    	
    	
    	if (this.filter == null) {
    		return getFilteredMutationSummary(mutations, node);
    	}else {
    		
    		// Filter the mutations and then give the remainder to the child class
    		List<Mutation> filteredMutations = new ArrayList<Mutation>();
        	for (Mutation mutation : mutations) {
        		if (this.filter.contains(mutation.getSiteNr())){
        			filteredMutations.add(mutation);
        		}
        	}
        	
        	
        	return getFilteredMutationSummary(filteredMutations, node);
    	}
    	
    	
    }
    
    
    
    @Override
    public Object getPropertyOfNode(Node node) {
    	return getArrayValue(node.getNr());
    }
    
    
    /**
     * For example, the number of mutations on the branch
     * @param nodeNr
     * @return
     */
    public abstract double getFilteredMutationSummary(List<Mutation> mutations, Node node);
    
    
    
    /**
     * The name to log onto the tree / log file
     * @return
     */
	//protected abstract String getName();
	
	
	/**
	 * Is this data type supported?
	 * @param dataType
	 * @return
	 */
	protected abstract boolean canHandleDataType(DataType dataType);
    
	
    /**
     * Please call this before calling log or getArrayValue, as the values are sampled once per state at most
     * @param sampleNr
     */
	@Override
    public void sampleMutations(long sampleNr) {
    	StochasticMapper sampler = samplerInput.get();
    	if (sampler == null) return;
    	sampler.sampleMutations(sampleNr);
    }


    @Override
    public double getArrayValue(int dim) {
    	
    	
    	Node speciesNode = this.getTree().getNode(dim);
    	List<Node> nodes = getAllNodesOnBranch(speciesNode);
    	double total = 0;
    	
    	for (Node node : nodes) {
    		
    		List<Mutation> mutations;
    		if (samplerInput.get() != null) {
        		
        		// Make sure to call sampleMutations before calling this method so that it depends on the current stochastic sample
        		mutations = getSpeciesTreeMutationSummary(dim, node);
        		
        	}else {
        		mutations = truthInput.get().getMutationsOnBranch(dim); // TODO implement on MSC for simulation
        	}
        	
    		
    		Collections.sort(mutations);
    		
    		double x = getMutationSummary(mutations, node);
    		if (this.geneTreePrior != null) {
    			//Log.warning("snode=" + dim + " gnode=" + node.getNr() + " nmut=" + mutations.size() + " x=" + x);
    		}
    		
    		total += x;
    	}
    	
    	return total;
    	
    	
    	
    }
    
    
    // Get all nodes within this species branch (or just return the branch itself if there aren't any)
    private List<Node> getAllNodesOnBranch(Node speciesNode){
    	
    	List<Node> nodes = new ArrayList<>();
    	if (this.geneTreePrior == null) {
    		nodes.add(speciesNode);
    	}else {
    		
    		Tree geneTree = (Tree) this.geneTreePrior.getGeneTree();
			for (Node geneNode : geneTree.getRoot().getAllChildNodesAndSelf()) {
				
				if (geneNode.isRoot()) continue;
				
				// Consider only the gene branches that reside within this species branch
				if (!this.geneTreePrior.mapGeneBranchToSpeciesNodes(geneNode.getNr()).contains(speciesNode)) {
					continue;
				}
				nodes.add(geneNode);
			}
    		
    	}
    	
    	return nodes;
    	
    }
    
    
    // If we are counting a species tree, then consider all gene tree branch segments within this species branch
    protected List<Mutation> getSpeciesTreeMutationSummary(int speciesNodeNr, Node geneNode){
    	

    	if (this.geneTreePrior == null) {
    		return samplerInput.get().getMutationsOnBranch(speciesNodeNr);
    	}else {
    		
    		List<Mutation> mutationsSpeciesBranch = new ArrayList<Mutation>();
    		
    		if (geneNode.isRoot()) return mutationsSpeciesBranch;
    		
    		Tree speciesTree = this.getTree();
    		Node speciesNode = speciesTree.getNode(speciesNodeNr);
    		//Tree geneTree = (Tree) this.geneTreePrior.getGeneTree();

    		// Start and end times of interval
    		double start = speciesNode.isRoot() ? Double.POSITIVE_INFINITY : speciesNode.getParent().getHeight();
    		double end = speciesNode.getHeight();
    		
			
			// Consider only the gene branches that reside within this species branch
			if (!this.geneTreePrior.mapGeneBranchToSpeciesNodes(geneNode.getNr()).contains(speciesNode)) {
				return mutationsSpeciesBranch;
			}
			
			
			List<Mutation> mutationsGeneBranch = samplerInput.get().getMutationsOnBranch(geneNode.getNr());
			
			// Renumber the mutations on the branch
			for (Mutation mut : mutationsGeneBranch) {
				
				// Count all mutations that occur within this species branch
				double mutHeight = geneNode.getParent().getHeight() - mut.getTime();
				if (mutHeight < 0) mutHeight = 0;
				if (mutHeight < start && mutHeight >= end) {
					mutationsSpeciesBranch.add(mut);
				}
			}
				
        	
        	return mutationsSpeciesBranch;
    		
    	}
    	
    	
    }





	@Override
	public void init(PrintStream out) {
		Tree tree = getTree();
		for (int i = 0; i < getDimension(); i ++) {
			String id = tree.getNode(i).getID();
			if (id == null) id = "node" + i;
			out.print(this.getID() + "." + id + "\t");
		}
		
	}



	@Override
	public void log(long sampleNr, PrintStream out) {
		
		this.sampleMutations(sampleNr);
		
		for (int nodeNr = 0; nodeNr < getDimension(); nodeNr ++) {
			double nmut = getArrayValue(nodeNr);
			out.print(nmut + "\t");
		}
		
		
	}
	
	protected Tree getTree() {
		
		if (geneInput.get() != null) {
			return geneInput.get().speciesTreeInput.get();
		}
		
		// Locus tree
		if (samplerInput.get() == null) {
			return (Tree) truthInput.get().getTree();
		}else {
			return (Tree) samplerInput.get().getTree();
		}
	}



	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}
	
    @Override
    public int getDimension() {
    	
    	Tree tree = getTree();
    	if (includeRootInput.get()) return tree.getNodeCount();
        return tree.getNodeCount()-1;
    }


    
    protected int[] getSequenceForNode(Node node) {
    	
    	if (samplerInput.get() == null) {
			return truthInput.get().getStatesForNode(getTree(), node);
		}else {
			return samplerInput.get().getStatesForNode(getTree(), node);
		}
    	
    	
    }
    
    
    protected DataType getDataType(){
    	return samplerInput.get() != null ? samplerInput.get().getDataTypeOfMapper() : truthInput.get().getDataTypeOfMapper();
    }
    
    
    
    /**
     * We are assuming that the number of patterns is the number of sites, which will be true if a PatternlessAlignment is used
     * @return
     */
    protected int getSiteAndPatternCount(){
    	int siteCount = samplerInput.get() != null ? samplerInput.get().getPatternCount() : truthInput.get().getPatternCount();
    	return siteCount;
    }
    

    /**
     * Return number of synonymous and non-synonymous substitutions in a array of length 2
     * Assuming that the data type is nucleotide
     * @param mutations
     * @param code
     * @param codon
     * @param openReadingFrame
     * @return
     */
	protected int[] getSynonymousAndNonSynonymousSubstitutionCount(List<Mutation> mutations, GeneticCode code, Codon codon, int openReadingFrame) {
		
		StochasticMapper mapper = samplerInput.get() != null ? samplerInput.get() : truthInput.get();
		
		if (this.getDataType() instanceof Nucleotide) {
			return MutationUtils.getSNCountForNucleotides(mutations, code, codon, openReadingFrame, mapper);
			//return getSNCountForNucleotides(mutations, code, codon, openReadingFrame);
		}
		
		if (this.getDataType() instanceof Codon) {
			return MutationUtils.getSNCountForCodons(mutations, code, codon, mapper);
			//return getSNCountForCodons(mutations, code, codon);
		}
		
		return null;
		
	}
	
	/**
	 *  Return number of synonymous and non-synonymous substitutions in a array of length 2, but using the 
	 */
	protected int[] getUnconditionalSynonymousAndNonSynonymousSubstitutionCount(List<Mutation> mutations, GeneticCode code, Codon codon, int openReadingFrame) {
		if (this.getDataType() instanceof Nucleotide) {
			return MutationUtils.getSNCountForNucleotides(mutations, code, codon, openReadingFrame, samplerInput.get().getUnconditionalData());
			//return getSNCountForNucleotides(mutations, code, codon, openReadingFrame);
		}
		
		if (this.getDataType() instanceof Codon) {
			return MutationUtils.getSNCountForCodons(mutations, code, codon, samplerInput.get().getUnconditionalData());
			//return getSNCountForCodons(mutations, code, codon);
		}
		
		return null;
	}

	
	
	
	
	

	






    
}



