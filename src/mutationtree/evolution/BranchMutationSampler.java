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

@Description("Stochastically samples a mutation trajectory along each branch for a given state")
public class BranchMutationSampler extends AncestralSequenceTreeLikelihood {
	
	
	
	long lastSample;
	List<List<Mutation>> mutationsAlongEachBranch;
	SiteModel siteModel;
	SubstitutionModel substitutionModel; 
	
	
	final double EPSILON = 1e-10; // The small amount of time considered as 1 dt for computing the instantaneous transition rates from a P matrix
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
		 
		 List<Mutation> arr = new ArrayList<>();
		 int nstates = qmatrix.getStateCount();
		 
		 // What do to if gap? or ambiguous?
		 if (parent < 0 || child < 0 || parent >= nstates || child >= nstates) {
			 
			 // TODO
			 return arr;
		 }
		 
		
		 
		 int loopNr = 0;
		 int finalState = -1;
		 while (finalState != child) {
			 
			 
			 arr.clear();
			 
			 int from = parent;
			 double t = 0;
			 while (true) {
				 
				 
				 //When does the next mutation occur?
				 double dt = Randomizer.nextExponential(clockRate);
				 
				// Log.warning("time " + dt + " out of " + time);
				 
				 if (t + dt > time) break;
				 
				 
				 // What was the mutation?
				 double[] outRates = getTransitionProbabilities(qmatrix, from, node);
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
			 
			 
			 // Hopefully the final state is the child or we will need to repeat that loop again
			 finalState = from;
			 
			 if (finalState != child) {
				//Log.warning("want " + child + " but getting " + finalState);
			 }
			 
			 loopNr++;
			 if (loopNr > MAX_N_LOOPS) {
				 Log.warning("Taking too long to get mutations on branch. Setting to parsimonious estimate. If this error happens to often, the rejection sampler is not working properly.");
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
	 * Many substitution model implementations do not return the Q matrix, and some of these are node/time dependent
	 * The most general approach that works around the setup of beast2 is to callculate P(Qt) for a very small t, which will be a close approximation of Q
	 * @param substModel
	 * @param fromState
	 * @param node
	 * @return
	 */
	private double[] getTransitionProbabilities(SubstitutionModel substModel, int fromState, Node node) {
		
		
		
		int nstates = substModel.getStateCount();
		double[] matrix = new double[nstates*nstates];
		
		// From time 0 to time epsilon, assuming a clock rate of 1 change per unit of time. If epsilon is too small, the probabilities may underflow
		substModel.getTransitionProbabilities(node, 0, EPSILON, 1, matrix);
		
		
		double outSum = 0;
		double[] outProbs = new double[nstates];
		for (int toState = 0; toState < nstates; toState ++) {
			if (toState == fromState) {
				outProbs[toState] = 0;
			}else {
				int index2d = fromState*nstates + toState;
				outProbs[toState] = matrix[index2d];
			}
			outSum += outProbs[toState];
		}
		
		
		if (outSum < 1e-100) {
			Log.warning("Zero sum. Epsilon is too small.");
		}
		for (int toState = 0; toState < nstates; toState ++) {
			outProbs[toState] = outProbs[toState] / outSum;
		}
		
		return outProbs;
		
		
	}
	 
	
	
	
	
	
}




