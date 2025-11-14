package beastmap.evolution;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastmap.util.Mutation;
import starbeast3.evolution.speciation.GeneTreeForSpeciesTreeDistribution;

public class BranchMutationSamplerSpeciesTree  extends BEASTObject implements StochasticMapper {

	final public Input<Tree> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);
	final public Input<List<StochasticMapper>> samplerInput = new Input<>("sampler", "mutation sampler to join", new ArrayList<>());
	final public Input<List<GeneTreeForSpeciesTreeDistribution>> geneInput = new Input<>("gene", "gene tree to species tree mapping (one per sampler)", new ArrayList<>());
	
	
	long lastSample = -1;
	int nsites;
	int ngenes;
	List<StochasticMapper> samplers;
	List<GeneTreeForSpeciesTreeDistribution> genes;
	List<List<Mutation>> mutationsAlongEachBranch = new ArrayList<>();
	//int[][] sampledSequences;
	
	@Override
	public void initAndValidate() {
		
		this.samplers = samplerInput.get();
		this.genes = geneInput.get();
		if (this.samplers.isEmpty() || this.samplers.size() != this.genes.size()) {
			throw new IllegalArgumentException("Please provide one stochastic mapper per gene-to-species mapper");
		}
		ngenes = this.genes.size();
		
		
		// Ensure all samplers have the same datatype
		DataType dt = samplerInput.get().get(0).getDataTypeOfMapper();
		for (StochasticMapper mapper : samplerInput.get()) {
			if (mapper.getDataTypeOfMapper().getClass() != dt.getClass()) {
				throw new IllegalArgumentException("Please ensure that all sampler datasets have the same datatype " + dt.getClass().getCanonicalName() + " != " + mapper.getDataTypeOfMapper().getClass().getCanonicalName());
			}
		}
		
		// Ensure the gene and stochastic mappers have the same trees
		for (int g = 0; g < ngenes; g++) {
			StochasticMapper sampler = this.samplers.get(g);
			GeneTreeForSpeciesTreeDistribution genePrior = this.genes.get(g);
			Tree geneTree = sampler.getTree();
			Tree geneTree2 = (Tree) genePrior.getGeneTree();
			if (geneTree != geneTree2) {
				throw new IllegalArgumentException("Tree number " + g + " does not match across the stochastic mapper and gene mapper: " + geneTree.getID() + " != " + geneTree2.getID());
			}
		}
		
		
		// Length
		nsites = 0;
		for (StochasticMapper mapper : samplerInput.get()) {
			nsites += mapper.getPatternCount();
		}
		
		this.lastSample = -1;
		
	}
	
	
	@Override
	public void sampleMutations(long sample) {
		
		
		// Only do this once per logged state
    	if (sample == this.lastSample) return;
		
    	
    	// Sample the mutations
		for (StochasticMapper mapper : samplerInput.get()) {
			mapper.sampleMutations(sample);
		}

		
		
		Tree speciesTree = this.getTree();
    	int nnodes = speciesTree.getNodeCount();
		
		// Join them
		this.mutationsAlongEachBranch.clear();
		//this.sampledSequences = new int[speciesTree.getNodeCount()][nsites];
		

    	
    	for (int speciesNodeNr = 0; speciesNodeNr < nnodes; speciesNodeNr ++) {
    		
    		Node speciesNode = speciesTree.getNode(speciesNodeNr);
    		List<Mutation> mutationsBranchTotal = new ArrayList<>();
    		
    		
    		// Start and end times of interval
    		double start = speciesNode.isRoot() ? Double.POSITIVE_INFINITY : speciesNode.getParent().getHeight();
    		double end = speciesNode.getHeight();
    		
    		
    		int siteNr = 0;
    		for (int g = 0; g < ngenes; g++) {
    			StochasticMapper sampler = this.samplers.get(g);
    			GeneTreeForSpeciesTreeDistribution genePrior = this.genes.get(g);
    			Tree geneTree = sampler.getTree();
    			//Node[] geneNodesInSpecies = genePrior.mapSpeciesNodeToGeneTreeNodes(speciesNode);
    			//List<Node> geneBranchesInSpecies = new ArrayList<>();
    			
    			
    			//Tree geneTree = sampler.getTree();
    			
    			for (Node geneNode : geneTree.getRoot().getAllChildNodesAndSelf()) {
    				
    				
    				if (geneNode.isRoot()) continue;
    				
    				// Consider only the gene branches that reside within this species branch
    				if (!genePrior.mapGeneBranchToSpeciesNodes(geneNode.getNr()).contains(speciesNode)) {
    					continue;
    				}
    				
//    				if (speciesNode.isRoot()) {
//    					Log.warning("root " + geneNode.getNr() + " for " + g + " geneNode root " + geneNode.isRoot());
//    				}
    			
	    				
    				
    				
        			List<Mutation> mutationsBranch = sampler.getMutationsOnBranch(geneNode.getNr());
        			
        			// Renumber the mutations on the branch
        			for (Mutation mut : mutationsBranch) {
        				
        				
        				// Count all mutations that occur within this species branch
        				double mutHeight = geneNode.getParent().getHeight() - mut.getTime();
        				if (mutHeight < 0) mutHeight = 0;
        				
//        				if (speciesNode.isRoot()) {
//        					Log.warning("root " + start + " - " + end + " ? " + mutHeight);
//        				}
        				if (mutHeight < start && mutHeight >= end) {
        					Mutation copy = new Mutation(mut);
	        				copy.setSiteNr(copy.getSiteNr() + siteNr);
	        				mutationsBranchTotal.add(copy);
        				}
        				
        				
        			}
	    			
	    			
	    			// Copy node sequence
	    			//System.arraycopy(sampler.getStatesForNode(geneTree, geneNode), 0, this.sampledSequences[speciesNodeNr], siteNr, sampler.getPatternCount());
	    			
    			
    			}
    			
    			// Limitation: will be challenging to ensure we are not counting across boundaries
    			siteNr += sampler.getPatternCount();
    			
    			
    			
    		}
    		
    		
    		
    		//Log.warning("BranchMutationSampler :" +  speciesNodeNr + " there are " + mutationsBranchTotal.size());
        	this.mutationsAlongEachBranch.add(mutationsBranchTotal);
    		

    		
    	}
    
		
		this.lastSample = sample;
		
	}
	
	

	
	@Override
	public Tree getTree() {
		return treeInput.get();
	}

	
	
	@Override
	public int[] getStatesForNode(Tree tree, Node node) {
		//return this.sampledSequences[node.getNr()];
		
		// Does not make sense. There are multiple sequences not just one.
		return null;
	}

	@Override
	public List<Mutation> getMutationsOnBranch(int nodeNr){
		return this.mutationsAlongEachBranch.get(nodeNr);
	}
	
	@Override
	public List<Mutation> getMutationsOnBranchAtSite(int nodeNr, int siteNr){
		List<Mutation> muts = this.mutationsAlongEachBranch.get(nodeNr);
		List<Mutation> siteMutations = new ArrayList<>();
		for (int i = 0; i < muts.size(); i ++) {
			Mutation m = muts.get(i);
			if (m.getSiteNr() == i) {
				siteMutations.add(m);
			}
		}
		
		return siteMutations;
	}
	

	
	@Override
	public DataType getDataTypeOfMapper() {
		return samplerInput.get().get(0).getDataTypeOfMapper(); // All data types are the same
	}

	@Override
	public int getPatternCount() {
		return nsites;
	}

	@Override
	public StochasticMapper getUnconditionalData() {
		return null;
	}


	
	
	
	

}
