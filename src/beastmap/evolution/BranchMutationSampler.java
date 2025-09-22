package beastmap.evolution;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import beastmap.indel.SimpleIndelCodingAlignment;
import beastmap.util.Mutation;
import beastmap.util.MutationUtils;

@Description("Stochastically samples a mutation trajectory along each branch for a given state")
public class BranchMutationSampler extends AncestralSequenceTreeLikelihood implements StochasticMapper{
	
	
	
	public Input<Boolean> nodeDependentInput = new Input<Boolean>("substModelIsNodeDependent", "set to false if the substituion model does not vary among branches (and gain some extra efficiency). but this will return the wrong answer for epoch subst models", true);
	public Input<BranchMutationSampler> indelInput = new Input<BranchMutationSampler>("indel", "another mapper can be used to remove all sites that are predicted to have gaps");
	public Input<Integer> burninInput = new Input<Integer>("burnin", "do not print any numbers until we have reached this state number (in case the tree branches are too long at the initial state)", 0);
	
	
	long lastSample;
	List<List<Mutation>> mutationsAlongEachBranch;
	SiteModel.Base siteModel;
	SubstitutionModel substitutionModel; 
	
	
	double[][] qmatrixBase;

	final int MAX_N_LOOPS = 1000000;
	final int MAX_MUT_BRANCH = 10000; // Any more than this and we risk running out of memory
	int gapChar;
	
	
	int n1 = 0;
	int n2 = 0;
	
	SimulatedAlignmentWithMutations unconditionalData; // Data simulated unconditional on the observed data
	SimpleIndelCodingAlignment indelData;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		this.lastSample = -1;
		
		this.siteModel = (SiteModel.Base) siteModelInput.get();
		this.substitutionModel = siteModel.getSubstitutionModel();
		
		qmatrixBase = null;
		
		
		if (! (dataInput.get() instanceof PatternlessAlignment)) {
			throw new IllegalArgumentException("Please ensure that the data is of type " + PatternlessAlignment.class.getName() + ". Currently it is " + getData().getClass());
		}
		
		
		if (indelInput.get() != null) {
			
			PatternlessAlignment data = (PatternlessAlignment)indelInput.get().getData();
			if (! (data.alignmentInput.get() instanceof SimpleIndelCodingAlignment)) {
				throw new IllegalArgumentException("Please ensure that the indel data is of type " + SimpleIndelCodingAlignment.class.getName() + ". Currently it is " + data.alignmentInput.get());
			}
			indelData = (SimpleIndelCodingAlignment)data.alignmentInput.get();
			gapChar = indelData.getGapChar();
		}
		
		
		// Initialise empty arrays
		this.mutationsAlongEachBranch = new ArrayList<>();
		this.mutationsAlongEachBranch.clear();
		Tree tree = (Tree) treeInput.get();
		int nbranches = tree.getNodeCount()-1;
		for (int nodeNr = 0; nodeNr < nbranches; nodeNr ++) {
			this.mutationsAlongEachBranch.add(new ArrayList<>());
		}
		
		
//		Log.warning("datatype " + getDataTypeOfMapper().getClass());		
//		// Build a mock alignment
//		List<Sequence> seqs = new ArrayList<>();
//		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
//			Sequence s = new Sequence(tree.getNode(i).getID(), getDataTypeOfMapper().encodingToString(new int[] { 0 })); // Dummy sequence
//			seqs.add(s);
//		}
//		Alignment mockData = new Alignment(seqs, getDataTypeOfMapper().getTypeDescription());
//		
//		unconditionalData = new SimulatedAlignmentWithMutations();
//	    int nsites = dataInput.get().getPatternCount();
//	    unconditionalData.initByName("data", mockData, "tree", tree, "siteModel", this.siteModel, "branchRateModel", this.branchRateModel, "sequencelength", nsites);
//		
	}

	
	@Override
	public StochasticMapper getUnconditionalData() {
		return unconditionalData;
	}
	




	@Override
	 public void sampleMutations(long sample) {
		
		
		//Log.warning(this.getID() + " sampling regular mutations ");
		
		// Only do this once per logged state
    	if (sample == this.lastSample) return;
    	
    	if (sample < burninInput.get()) {
    		return;
    	}
		
		BranchMutationSampler indels = indelInput.get();
		if (indels != null) indelInput.get().sampleMutations(sample);
	    	
    	
    	
    	Tree tree = (Tree) treeInput.get();
    	int nbranches = tree.getNodeCount()-1;
    	
    	
    	// Reset list of per-branch mutations
    	this.mutationsAlongEachBranch.clear();
    	
    	
    	// We can just do this once at the start and not repeat on every branch
    	if (!nodeDependentInput.get()) {
			qmatrixBase = MutationUtils.getTransitionRates(substitutionModel, null, false);
		}
    	
    	
    	// Set up rate matrix
    	//this.generalSubstitutionModel.setupRateMatrix();
    	//double[][] qmatrix = this.generalSubstitutionModel.getRateMatrix();
    	
    	
    	// Ancestral sequence reconstruction
    	super.calculateLogP();
    	super.redrawAncestralStates();
    	
    	
    
    	
    	// If the PatternlessAlignment is being used, then each site has its own pattern (which is best here)
    	int nsites = dataInput.get().getPatternCount();
    	
    	// Sample trajectory of mutations on each branch
    	
    	int maxOnInternal = 0; 
    	for (int nodeNr = 0; nodeNr < nbranches; nodeNr ++) {
    		
    		
    		List<Mutation> mutationsBranch = new ArrayList<Mutation>();
    		Node node = tree.getNode(nodeNr);
    		if (node.isRoot()) continue;
    		int[] siteStates = super.getStatesForNode(tree, node);
    		int[] parentStates = super.getStatesForNode(tree, node.getParent());
    		
    		
    		// Mask any sites that correspond to deletions
    		maskWithGaps(indels, siteStates, node, nsites);
    		
    		
//        	if (indelInput.get() != null) {
//        		
//        		int[] gapsChildCompact = indels.getStatesForNode(tree, node);
//        		//int[] gapsParentCompact = indels.getStatesForNode(tree, node.getParent());
//        		boolean[] gapsChild = indelData.expandIndelCoding(gapsChildCompact);
//        		//boolean[] gapsParent = indelData.expandIndelCoding(gapsParentCompact);
//        		
//        		
//        		if (gapsChild.length != siteStates.length) {
//        			throw new IllegalArgumentException("Please ensure that the original dataset in 'indel' is the same as this one. The lengths are not matching: " + gapsChild.length  + "!=" + siteStates.length);
//        		}
//        		
//        		for (int siteNr = 0; siteNr < nsites; siteNr ++) {
//        			if (!gapsChild[siteNr]) {
//        				siteStates[siteNr] = gapChar;
//        			}
////        			if (!gapsParent[siteNr]) {
////        				parentStates[siteNr] = gapChar;
////        			}
//        		}
//        		
//        	}
    		
    		
    		double branchRate = this.branchRateModel.getRateForBranch(node);
    		
    		for (int siteNr = 0; siteNr < nsites; siteNr ++) {
    			
    			// What category was sampled for this site?
    			int siteCategory = ((RateCategorySampledLikelihoodCore)this.likelihoodCore).getSampledSiteCategory(siteNr);
    			double siteRate = siteModel.getRateForCategory(siteCategory, node);
    			
    			List<Mutation> mutationsSite = sampleMutationsConditionalEndPoints(parentStates[siteNr], siteStates[siteNr], node.getLength(), siteNr, this.substitutionModel, branchRate*siteRate, node);
    			mutationsBranch.addAll(mutationsSite);
    		}
    		
    		
    		//Log.warning("BranchMutationSampler :" +  nodeNr + " there are " + mutationsBranch.size());
        	this.mutationsAlongEachBranch.add(mutationsBranch);
    		
        	
        	if (nodeNr >+ tree.getLeafNodeCount()) {
        		maxOnInternal = Math.max(maxOnInternal, mutationsBranch.size());
        	}
    		
    	}
    
    	
    	if (maxOnInternal == 0) {
    		// Sometimes we have a state where all mutations appear on the leaves and none on the internal nodes. A bug or feature in the ASR
    		//Log.warning("Warning: there are no mutations on the internal branches for some reason");
    	}
    	
    	
    	// Sample unconditional
    	if (unconditionalData != null) {
    		unconditionalData.simulate();
    	}
    	
    	
    	// Update sample number
		sample = this.lastSample;
		
		
		
//		for (int nodeNr = 0; nodeNr < nbranches; nodeNr ++) {
//			n1 += this.getMutationsOnBranch(nodeNr).size();
//			n2 += unconditionalData.getMutationsOnBranch(nodeNr).size();
//		}
//		
//		if (n1 > n2) {
//			Log.warning("UNCOND " + n1 + " / " + n2);
//		}else {
//			Log.warning("COND " + n1 + " / " + n2);
//		}
//		
		
    	
    	
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
	public List<Mutation> getMutationsOnBranch(int nodeNr){
		return this.mutationsAlongEachBranch.get(nodeNr);
	}
	
	@Override
	public List<Mutation> getMutationsOnBranchAtSite(int nodeNr, int siteNr){
		List<Mutation> muts = this.mutationsAlongEachBranch.get(nodeNr);
		List<Mutation> siteMutations = new ArrayList<>();
		for (int i = 0; i < muts.size(); i ++) {
			Mutation m = muts.get(i);
			if (m.getSiteNr() == siteNr) {
				siteMutations.add(m);
			}
		}
		
		return siteMutations;
	}

	
	/**
	 * Using the modified rejection sampler described below. There are more sophisticated methods but they are not always better
	 * https://pmc.ncbi.nlm.nih.gov/articles/PMC2818752/
	 */
	private List<Mutation> sampleMutationsConditionalEndPoints(int parent, int child, double time, int siteNr,  SubstitutionModel qmatrix, double clockRate, Node node){
		 
		 List<Mutation> arr = null;
		 int nstates = qmatrix.getStateCount();
		 
		 // What do to if gap? or ambiguous?
		 if (parent < 0 || child < 0 || parent >= nstates || child >= nstates) {
			 
			 // TODO
			 return new ArrayList<>();
		 }
		 
		 
		 // Allow each node to have its own q matrix. But if the matrix changes mid branch, we will have problems
		 double[][] qmatrixNode = qmatrixBase != null ? qmatrixBase : MutationUtils.getTransitionRates(qmatrix, node, this.nodeDependentInput.get());
		 
		 
		 int loopNr = 0;
		 while (arr == null) {
			 
			 
			 // Simple rejection sampling
			 if (parent == child) {
				 arr = rejectionSampleCase1(time, parent, clockRate, qmatrixNode, node, siteNr);
			 }
			 
			 // Condition the rejection sampling on knowing there was at least one mutation
			 else {
				 arr = rejectionSampleCase2(time, parent, child, clockRate, qmatrixNode, node, siteNr);
			 }
			 
			 loopNr++;
			 if (loopNr > MAX_N_LOOPS) {
				 arr = new ArrayList<>();
				 Log.warning(this.getID() + ": Taking too long to get mutations on branch. Setting to parsimonious estimate. If this warning happens to often, it is because the rejection sampler is too slow. This is because there are too many mutations along the branches. Maybe one of the clock/substitution rates is too low?");
				 if (parent != child) {
					 Mutation mut = new Mutation(parent, child, -1, siteNr, parent, child, node);
					 arr.add(mut);
					 break;
				 }
			 }
		 }
		 


		 
		 return arr;
		 
	 }
	
	
	
	/**
	 * Rejection sampler where the start and end state are the same (nothing special just a regular rejection sampler)
	 * Algorithm 2 Case 1 from https://pmc.ncbi.nlm.nih.gov/articles/PMC2818752/
	 * @return
	 */
	private List<Mutation> rejectionSampleCase1(double time, int state, double clockRate, double[][] qmatrix, Node node, int siteNr){
		
		
		List<Mutation> arr = new ArrayList<>();
		
		int from = state;
		double t = 0;
		while (true) {
			 
			 
			double lambda = 0;
			double[] outRates = qmatrix[from];
			for (int i = 0; i < outRates.length; i ++) lambda += outRates[i]*clockRate;
			
			//When does the next mutation occur?
			double dt = Randomizer.nextExponential(lambda);
		 
			if (t + dt > time) break;
		 
		 
			// What was the mutation?
			outRates[from] = 0;
			int nextState = Randomizer.randomChoicePDF(outRates);
			//Log.warning("mutated from " + from + " to " + nextState + " with rates (" + outRates[0] + "," + outRates[1] + "," + outRates[2] + "," + outRates[3] + ")" );
		 
			// Make mutation
			Mutation mut = new Mutation(from, nextState, t+dt, siteNr, state, state, node);
			arr.add(mut);
		 
		 
			// Increment time, update parental state
			from = nextState;
			t = t + dt;
			
			if (arr.size() > MAX_MUT_BRANCH) {
				 Log.warning(this.getID() + ": rs1 -- too many mutations on each branch! Perhaps branches are too long or clock rate too slow. Using parsimonious estimate.");
				 return new ArrayList<Mutation>();
			}
			 
		}
		
		
		// Success
		if (from == state) {
			return arr;
		}
		
		// Failure
		else {
			return null;
		}
		 
		
		 
	}
	
	/**
	 * Rejection sampler where the start and end state are different
	 * Algorithm 2 Case 2 from https://pmc.ncbi.nlm.nih.gov/articles/PMC2818752/
	 * @return
	 */
	private List<Mutation> rejectionSampleCase2(double time, int parent, int child, double clockRate, double[][] qmatrix, Node node, int siteNr){
		
		
		List<Mutation> arr = new ArrayList<>();
		
		int from = parent;
		double t = 0;
		boolean first = true;
		while (true) {
			
			
			double lambda = 0;
			double[] outRates = qmatrix[from];
			for (int i = 0; i < outRates.length; i ++) lambda += outRates[i]*clockRate;
			
			//When does the next mutation occur?
			double dt;
			if (first) {
				
				// Sample time conditional on there being at least one change
				dt = sampleTimeConditionalOnAtLeastOneChange(lambda, time);
				first = false;
				
			}else {
				dt = Randomizer.nextExponential(lambda);
			}
			 
			
		 
			if (t + dt > time) break;
		 
		 
			// What was the mutation?
			outRates[from] = 0;
			int nextState = Randomizer.randomChoicePDF(outRates);
			//Log.warning("mutated from " + from + " to " + nextState + " with rates (" + outRates[0] + "," + outRates[1] + "," + outRates[2] + "," + outRates[3] + ")" );
		 
			// Make mutation
			Mutation mut = new Mutation(from, nextState, t+dt, siteNr, parent, child, node);
			arr.add(mut);
		 
			// Increment time, update parental state
			from = nextState;
			t = t + dt;
			
			if (arr.size() > MAX_MUT_BRANCH) {
				 Log.warning(this.getID() + ": rs2 -- too many mutations on each branch! Perhaps branches are too long or clock rate too slow. Using parsimonious estimate.");
				 arr = new ArrayList<>();
				 mut = new Mutation(parent, child, -1, siteNr, parent, child, node);
				 arr.add(mut);
				 return arr;
			}
			 
		}
		
		
		// Success
		if (from == child) {
			return arr;
		}
		
		// Failure
		else {
			return null;
		}
		
	}
	
	
	/**
	 * Sample the waiting time to the first transition in time T, conditional on there being at least one transition during that period
	 * Equation 2.1 from https://pmc.ncbi.nlm.nih.gov/articles/PMC2818752/
	 * @param lambda - the total outgoing rate of initial state a (positive number)
	 * @param t - the total length of time
	 * @return
	 */
	public static double sampleTimeConditionalOnAtLeastOneChange(double lambda, double t) {
		double u = Randomizer.nextFloat();
		return -Math.log(1 - u*(1-Math.exp(-t*lambda))) / lambda;
	}


	@Override
	public Tree getTree() {
		return (Tree) treeInput.get();
	}


	@Override
	public int[] getStatesForNode(Tree tree, Node node) {
		return super.getStatesForNode(tree, node);
	}


	@Override
	public DataType getDataTypeOfMapper() {
		return dataInput.get().getDataType();
	}
	
	
	
	@Override
	public int getPatternCount() {
		return dataInput.get().getPatternCount();
	}


	public Alignment getData() {
		return dataInput.get();
	}
	
	
	
	
	
	
}




