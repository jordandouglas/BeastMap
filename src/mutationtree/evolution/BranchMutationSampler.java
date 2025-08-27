package mutationtree.evolution;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Log;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import mutationtree.util.Mutation;
import mutationtree.util.MutationUtils;

@Description("Stochastically samples a mutation trajectory along each branch for a given state")
public class BranchMutationSampler extends AncestralSequenceTreeLikelihood {
	
	
	
	long lastSample;
	List<List<Mutation>> mutationsAlongEachBranch;
	SiteModel siteModel;
	SubstitutionModel substitutionModel; 
	
	

	final int MAX_N_LOOPS = 100000000;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		this.lastSample = -1;
		this.mutationsAlongEachBranch = new ArrayList<>();
		this.siteModel = (SiteModel) siteModelInput.get();
		this.substitutionModel = siteModel.getSubstitutionModel();
		
	}

	
	
	 public void sampleMutations(long sample) {
	    	
    	// Only do this once per logged state
    	if (sample == this.lastSample) return;
    	
    	Tree tree = (Tree) treeInput.get();
    	int nbranches = tree.getNodeCount()-1;
    	
    	
    	// Reset list of per-branch mutations
    	this.mutationsAlongEachBranch.clear();
    	
    	
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
    	
    	
    	// Update sample number
		sample = this.lastSample;
    	
    	
    }
	 
	
	public List<Mutation> getMutationsOnBranch(int nodeNr){
		return this.mutationsAlongEachBranch.get(nodeNr);
	}
	
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
		 
		 
		 int loopNr = 0;
		 while (arr == null) {
			 
			 
			 // Simple rejection sampling
			 if (parent == child) {
				 arr = rejectionSampleCase1(time, parent, clockRate, qmatrix, node, siteNr);
			 }
			 
			 // Condition the rejection sampling on knowing there was at least one mutation
			 else {
				 arr = rejectionSampleCase2(time, parent, child, clockRate, qmatrix, node, siteNr);
			 }
			 
			 loopNr++;
			 if (loopNr > MAX_N_LOOPS) {
				 arr = new ArrayList<>();
				 Log.warning("Taking too long to get mutations on branch. Setting to parsimonious estimate. If this warning happens to often, it is because the rejection sampler is too slow.");
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
	private List<Mutation> rejectionSampleCase1(double time, int state, double clockRate, SubstitutionModel qmatrix, Node node, int siteNr){
		
		
		List<Mutation> arr = new ArrayList<>();
		
		int from = state;
		double t = 0;
		while (true) {
			 
			 
			double lambda = 0;
			double[] outRates = MutationUtils.getTransitionRates(qmatrix, from, node);
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
	private List<Mutation> rejectionSampleCase2(double time, int parent, int child, double clockRate, SubstitutionModel qmatrix, Node node, int siteNr){
		
		
		List<Mutation> arr = new ArrayList<>();
		
		int from = parent;
		double t = 0;
		boolean first = true;
		while (true) {
			
			
			double lambda = 0;
			double[] outRates = MutationUtils.getTransitionRates(qmatrix, from, node);
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
	
	
	
	
	
	
	
	
	
}




