package beastmap.evolution;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastmap.indel.SimpleIndelCodingAlignment;
import beastmap.util.Mutation;

@Description("Glues together two or more branch samplers into a single object. Useful for joining partitions")
public class StochasticMapJoiner extends BEASTObject implements StochasticMapper {

	// concat = generate a concatenated alignment D1-D2-D3
	// intersperse = glue together 1 site at a time D1(site1)-D2(site1)-D3(site1)-D1(site2)-...
			// useful for codon partitions
	String[] JOIN_OPTIONS = new String[] { "concat", "intersperse" };
	
	final public Input<List<StochasticMapper>> samplerInput = new Input<>("sampler", "mutation sampler to join", new ArrayList<>());
	final public Input<String> joinInput = new Input<>("join", "how to join the alignments together", JOIN_OPTIONS[0], JOIN_OPTIONS);
	public Input<BranchMutationSampler> indelInput = new Input<BranchMutationSampler>("indel", "another mapper can be used to remove all sites that are predicted to have gaps");
	
	
	String join;
	int nsites;
	long lastSample;
	List<List<Mutation>> mutationsAlongEachBranch = new ArrayList<>();
	int[][] sampledSequences;
	
	int gapChar;
	SimpleIndelCodingAlignment indelData;
	
	@Override
	public void initAndValidate() {
		this.join = joinInput.get();
		if (samplerInput.get().isEmpty()) {
			Log.warning("Please provide at least one sampler");
		}
		if (this.join.equals("intersperse")) {
			
			// When using a patternless alignment, pattern count = site count
			int length = samplerInput.get().get(0).getPatternCount();
			for (StochasticMapper mapper : samplerInput.get()) {
				if (mapper.getPatternCount() != length) {
					throw new IllegalArgumentException("Please ensure that all sampler datasets have the same length when using join='intersperse' " + length + " != " + mapper.getPatternCount());
				}
			}
			
		}
		
		// Ensure all samplers have the same tree
		Tree tree = samplerInput.get().get(0).getTree();
		for (StochasticMapper mapper : samplerInput.get()) {
			if (mapper.getTree() != tree) {
				throw new IllegalArgumentException("Please ensure that all sampler datasets have the same tree " + tree.getID() + " != " + mapper.getTree().getID());
			}
		}
		
		// Ensure all samplers have the same datatype
		DataType dt = samplerInput.get().get(0).getDataType();
		for (StochasticMapper mapper : samplerInput.get()) {
			if (mapper.getDataType().getClass() != dt.getClass()) {
				throw new IllegalArgumentException("Please ensure that all sampler datasets have the same datatype " + dt.getClass().getCanonicalName() + " != " + mapper.getDataType().getClass().getCanonicalName());
			}
		}
		
		// Length
		nsites = 0;
		for (StochasticMapper mapper : samplerInput.get()) {
			nsites += mapper.getPatternCount();
		}
		
		
		if (indelInput.get() != null) {
			
			PatternlessAlignment data = (PatternlessAlignment)indelInput.get().getData();
			if (! (data.alignmentInput.get() instanceof SimpleIndelCodingAlignment)) {
				throw new IllegalArgumentException("Please ensure that the indel data is of type " + SimpleIndelCodingAlignment.class.getName() + ". Currently it is " + data.alignmentInput.get());
			}
			indelData = (SimpleIndelCodingAlignment)data.alignmentInput.get();
		}
		
		
		
    	gapChar = getDataType().stringToEncoding(""+DataType.GAP_CHAR).get(0);
		
		
		
		lastSample=-1;
		
		
	}

	@Override
	public void sampleMutations(long sample) {
		
		BranchMutationSampler indels = indelInput.get();
		if (indels != null) indelInput.get().sampleMutations(sample);
		
	
		// Only do this once per logged state
    	if (sample == this.lastSample) return;
		
		for (StochasticMapper mapper : samplerInput.get()) {
			mapper.sampleMutations(sample);
		}

		
		Tree tree = this.getTree();
    	int nnodes = tree.getNodeCount();
		
		// Join them
		this.mutationsAlongEachBranch.clear();
		this.sampledSequences = new int[tree.getNodeCount()][nsites];
		

    	
    	for (int nodeNr = 0; nodeNr < nnodes; nodeNr ++) {
    		
    		Node node = tree.getNode(nodeNr);
    		List<Mutation> mutationsBranchTotal = new ArrayList<>();
    		
    		if (this.join.equals("concat")) {
    			
    			
    			int siteNr=0;
        		for (int samplerNr = 0; samplerNr < samplerInput.get().size(); samplerNr++) {
        			StochasticMapper sampler = samplerInput.get().get(samplerNr);
        			
        			if (!node.isRoot()) {
	        			List<Mutation> mutationsBranch = sampler.getMutationsOnBranch(nodeNr);
	        			
	        			// Renumber the mutations on the branch
	        			for (Mutation mut : mutationsBranch) {
	        				Mutation copy = new Mutation(mut);
	        				copy.setSiteNr(copy.getSiteNr() + siteNr);
	        				mutationsBranchTotal.add(copy);
	        			}
        			}
        			
        			// Copy node sequence
        			System.arraycopy(sampler.getStatesForNode(tree, node), 0, this.sampledSequences[nodeNr], siteNr, sampler.getPatternCount());

        			// Increment counter for the next alignment
        			siteNr += sampler.getPatternCount();
        			
        		}
    			
    			
    			
    		}
    		
    		else if (this.join.equals("intersperse")) {
    			
    			int newSiteNr = 0;
    			for (int siteNr = 0; siteNr < samplerInput.get().get(0).getPatternCount(); siteNr ++) {
    				for (int samplerNr = 0; samplerNr < samplerInput.get().size(); samplerNr++) {
    					StochasticMapper sampler = samplerInput.get().get(samplerNr);
    					
    					if (!node.isRoot()) {
    						
	    					List<Mutation> mutationsBranch = sampler.getMutationsOnBranchAtSite(nodeNr, siteNr);
	    					
	    					
	    					// Renumber the mutations on the branch
	            			for (Mutation mut : mutationsBranch) {
	            				
	            				//Log.warning("mutation on branch " + siteNr + " = " + mut.getSiteNr() );
	            				
	            				Mutation copy = new Mutation(mut);
	            				copy.setSiteNr(newSiteNr);
	            				mutationsBranchTotal.add(copy);
	            				
	            				//Log.warning("renumbered " + mut.getSiteNr() + "->" + copy.getSiteNr() + " for sampler " + samplerNr + " on " + siteNr + "/" + nsites);
	            			}
    					}
    					
    					this.sampledSequences[nodeNr][newSiteNr] = sampler.getStatesForNode(tree, node)[siteNr];
    					
    					newSiteNr ++;
    				}
    				
    			}
    			
    			
    			
    			
    		}
    		
    		maskWithGaps(indels, this.sampledSequences[nodeNr], node, nsites);

    		
    		//Log.warning("BranchMutationSampler :" +  nodeNr + " there are " + mutationsBranch.size());
        	this.mutationsAlongEachBranch.add(mutationsBranchTotal);
    		

    		
    	}
		
		
		this.lastSample = sample;
		
	}
	
	
	private void maskWithGaps(BranchMutationSampler indels, int[] siteStates, Node node, int nsites) {
		
		
		// Mask any sites that correspond to deletions
    	if (indels != null) {
    		
    		int[] gapsChildCompact = indels.getStatesForNode(getTree(), node);
    		boolean[] gapsChild = indelData.expandIndelCoding(gapsChildCompact);
    		
    		if (gapsChild.length != siteStates.length) {
    			throw new IllegalArgumentException("Please ensure that the original dataset in 'indel' is the same as this one. The lengths are not matching: " + gapsChild.length  + "!=" + siteStates.length);
    		}
    		
    		for (int siteNr = 0; siteNr < nsites; siteNr ++) {
    			if (!gapsChild[siteNr]) {
    				siteStates[siteNr] = gapChar;
    			}
    		}
    		
    	}
		
		
	}
	
	
	@Override
	public int[] getStatesForNode(Tree tree, Node node) {
		return this.sampledSequences[node.getNr()];
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
	public Tree getTree() {
		return samplerInput.get().get(0).getTree(); // All trees are the same
	}



	@Override
	public DataType getDataType() {
		return samplerInput.get().get(0).getDataType(); // All data types are the same
	}

	@Override
	public int getPatternCount() {
		return nsites;
	}
	
	

	
}



