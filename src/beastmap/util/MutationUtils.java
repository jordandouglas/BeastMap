package beastmap.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BetaDistributionImpl;
import org.apache.commons.math3.special.Beta;

import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.Nucleotide;
import beast.base.evolution.substitutionmodel.ComplexSubstitutionModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import beastmap.evolution.StochasticMapper;
import codonmodels.evolution.datatype.Codon;
import codonmodels.evolution.datatype.GeneticCode;
import starbeast3.evolution.speciation.GeneTreeForSpeciesTreeDistribution;

public class MutationUtils {
	
	
	final static double EPSILON = 1e-6; // The small amount of time considered as 1 dt for computing the instantaneous transition rates from a P matrix
	
	
	
	/**
	 * Create a single instance of uniform frequencies so that we don't need to do this every time and see the verbose printing
	 */
	private static Frequencies frequenciesVector = null;
	private static void prepareFrequencies(int nstates) {
		
		
		if (frequenciesVector != null && frequenciesVector.getFreqs().length == nstates) return;
		
		// Build a CTMC
		List<Double> uniformFreqs = new ArrayList<>();
		for (int i = 0; i < nstates; i ++) uniformFreqs.add(1.0 / nstates);
		RealParameter f = new RealParameter();
		f.initByName("value", uniformFreqs);
		frequenciesVector = new Frequencies();
		frequenciesVector.initByName("frequencies", f);
		
	}
	
	/**
	 * Many substitution model implementations do not return the Q matrix, and some of these are node/time dependent
	 * The most general approach that works around the setup of beast2 is to calculate P(Qt)/t for a very small t, whose off-diagonal elements will be a close approximation of Q
	 * Note that the diagonal element will be 0 not the negative row sum
	 * @param substModel
	 * @param fromState
	 * @param node
	 * @return
	 */
	public static double[] getTransitionRates(SubstitutionModel substModel, int fromState, Node node, int nstates) {
		
		
		
		double[] matrix = new double[nstates*nstates];
		
		// From time 0 to time epsilon, assuming a clock rate of 1 change per unit of time. If epsilon is too small, the probabilities may underflow
		substModel.getTransitionProbabilities(node, MutationUtils.EPSILON, 0, 1, matrix);
		
		
		double outSum = 0;
		double[] outRates = new double[nstates];
		for (int toState = 0; toState < nstates; toState ++) {
			if (toState == fromState) {
				outRates[toState] = 0;
			}else {
				int index2d = fromState*nstates + toState;
				outRates[toState] = matrix[index2d] / MutationUtils.EPSILON;
			}
			outSum += outRates[toState];
			
			
//			int index2d = fromState*nstates + toState;
//			outRates[toState] = matrix[index2d] / EPSILON;
		}
		
		
		if (outSum < 1e-100) {
			Log.warning("Zero sum. Epsilon is too small. " + outSum);
			//Log.warning(outRates[0] + " " + outRates[1] + " " + outRates[2] + " " + outRates[3]) ;
		}
	
		
		//Log.warning(outRates[0] + " " + outRates[1] + " " + outRates[2] + " " + outRates[3]) ;
		
		return outRates;
		
		
	}
	
	
	public static double[][] getTransitionRates(SubstitutionModel substModel, Node node, boolean nodeDependent, int nstates) {
		
		
		
		if (substModel instanceof GeneralSubstitutionModel && !nodeDependent) {
			GeneralSubstitutionModel gensub = (GeneralSubstitutionModel)substModel;
			gensub.setupRelativeRates();
			gensub.setupRateMatrix();
			return gensub.getRateMatrix();
		}
		
		
		double[] matrix = new double[nstates*nstates];
		
		// From time 0 to time epsilon, assuming a clock rate of 1 change per unit of time. If epsilon is too small, the probabilities may underflow
		substModel.getTransitionProbabilities(node, MutationUtils.EPSILON, 0, 1, matrix);
		
		
		double[][] outRates = new double[nstates][nstates];
		for (int fromState = 0; fromState < nstates; fromState ++) {
			
			double outSum = 0;
			for (int toState = 0; toState < nstates; toState ++) {
				if (toState == fromState) {
					outRates[fromState][fromState] = 0;
				}else {
					int index2d = fromState*nstates + toState;
					outRates[fromState][toState] = matrix[index2d] / MutationUtils.EPSILON;
				}
				outSum += outRates[fromState][toState];
			}
			if (outSum < 1e-100) {
				Log.warning("Zero sum. Epsilon is too small. " + outSum);
			}
		}
		
		//Log.warning(outRates[0] + " " + outRates[1] + " " + outRates[2] + " " + outRates[3]) ;
		
		return outRates;
		
		
	}
	
	
	/**
	 * Returns the sequence at any point in the tree, not just on a node
	 * @param mapper
	 * @param tree
	 * @param node
	 * @param height
	 * @return
	 */
	public static int[] getStatesForBranchAtTime(StochasticMapper mapper, Tree tree, Node node, double height, boolean inclusive) {
		
		
		if (height < 0 || node.isRoot() || height <= node.getHeight()) {
			return mapper.getStatesForNode(tree, node);
		}
		
		if (height > node.getParent().getHeight()) {
			throw new IllegalArgumentException("Dev error 3443: invalid height. " + node.getParent().getHeight() + " !> " + height + " !> " + node.getHeight());
			//return null;
		}
		
		
		// Get the sequence on the bottom of this branch
		int[] statesOnNode = mapper.getStatesForNode(tree, node);
		statesOnNode = Arrays.copyOf(statesOnNode, statesOnNode.length);
		
		
		// Go backwards in time until we hit the right height
		List<Mutation> mutationsOnbranch = mapper.getMutationsOnBranch(node.getNr());
		Collections.sort(mutationsOnbranch);
		for (int i = mutationsOnbranch.size()-1; i >= 0; i --) {
			
			Mutation mutation = mutationsOnbranch.get(i);
			double heightOfMutation = node.getParent().getHeight() - mutation.getTime();
			
			//Log.warning("heightOfMutation=" + heightOfMutation);
			
			if (inclusive) {
				if (heightOfMutation >= height) break;
			}else {
				if (heightOfMutation > height) break;
			}
			
			
			
			
			int siteNr = mutation.getSiteNr();
			int from = mutation.getFrom();
			int to = mutation.getTo();
			
			if (statesOnNode[siteNr] != to) {
				throw new IllegalArgumentException("Dev error 274278: " + statesOnNode[siteNr] + " != " + to);
			}
			
			// Mutate backwards
			statesOnNode[siteNr] = from;
			
		}
		
		
		//Log.warning("broken at=" + height);
		
		return statesOnNode;
		
		
	}
	
	

	/**
	 * Return a list of sites that correspond to a string filter
	 * Adapted from FilteredAlignment
	 */
    public static List<Integer> parseFilterSpec(int siteCount, String filterString) {
    	
    	
    	List<Integer> filter = new ArrayList<>();
    	
        // parse filter specification
        String[] filters = filterString.split(",");
        int[] from = new int[filters.length];
        int[] to = new int[filters.length];
        int[] step = new int[filters.length];
        for (int i = 0; i < filters.length; i++) {
            filterString = " " + filters[i] + " ";
            if (filterString.matches(".*-.*")) {
                // range, e.g. 1-100/3
                if (filterString.indexOf('\\') >= 0) {
                	String str2 = filterString.substring(filterString.indexOf('\\') + 1); 
                	step[i] = parseInt(str2, 1);
                	filterString = filterString.substring(0, filterString.indexOf('\\'));
                } else {
                	step[i] = 1;
                }
                String[] strs = filterString.split("-");
                from[i] = parseInt(strs[0], 1) - 1;
                to[i] = parseInt(strs[1], siteCount) - 1;
            } else if (filterString.matches(".*:.*:.+")) {
                // iterator, e.g. 1:100:3
                String[] strs = filterString.split(":");
                from[i] = parseInt(strs[0], 1) - 1;
                to[i] = parseInt(strs[1], siteCount) - 1;
                step[i] = parseInt(strs[2], 1);
            } else if (filterString.trim().matches("[0-9]*")) {
                from[i] = parseInt(filterString.trim(), 1) - 1;
                to[i] = from[i];
            	step[i] = 1;
            } else {
                throw new IllegalArgumentException("Don't know how to parse filter " + filterString);
            }
        }
        
        
        // Calculate filter
        boolean[] isUsed = new boolean[siteCount];
        for (int i = 0; i < to.length; i++) {
            for (int k = from[i]; k <= to[i]; k += step[i]) {
                isUsed[k] = true;
            }
        }
        // count
        int k = 0;
        for (int i = 0; i < isUsed.length; i++) {
            if (isUsed[i]) {
                k++;
            }
        }
        
        
        // set up index set
        k = 0;
        for (int i = 0; i < isUsed.length; i++) {
            if (isUsed[i]) {
            	filter.add(i);
            }
        }
        
        
        return filter;
        
    }

    
    /**
	 * Adapted from FilteredAlignment
	 */
    public static int parseInt(String str, int defaultValue) {
        str = str.replaceAll("\\s+", "");
        try {
            return Integer.parseInt(str);
        } catch (Exception e) {
            return defaultValue;
        }
    }



	    /**
     * Simulate directly using Gillespie's algorithm and return the child state
     * @param node
     * @param clockRate
     * @param substModel
     * @param mutations
     * @return
     */
    public static int simulateMutationsDownBranch(int parentState, Node node, double clockRate, SubstitutionModel qmatrix, List<Mutation> arr, int siteNr, int nstates) {
    	

    	
		int from = parentState;
		double t = 0;
		double time = node.getLength();
		while (true) {
			 

			// Outgoing rate
			double lambda = 0;
			double[] outRates = MutationUtils.getTransitionRates(qmatrix, from, node, nstates);
			for (int r = 0; r < outRates.length; r++) lambda += outRates[r]*clockRate;
			
			//When does the next mutation occur?
			double dt = Randomizer.nextExponential(lambda);
		 
			if (t + dt > time) break;
		 
		 
			// What was the mutation?
			int nextState = Randomizer.randomChoicePDF(outRates);
			//Log.warning("mutated from " + from + " to " + nextState + " with rates (" + outRates[0] + "," + outRates[1] + "," + outRates[2] + "," + outRates[3] + ")" );
		 
			// Make mutation
			Mutation mut = new Mutation(from, nextState, t+dt, siteNr, parentState, -1, node);
			arr.add(mut);
		 
		 
			// Increment time, update parental state
			from = nextState;
			t = t + dt;
			 
		}
    		
    		 
    	return from;
    }

    
    

    public static int[] getSNCountForCodons(List<Mutation> mutations, GeneticCode code, Codon codon, StochasticMapper mapper) {

		
		
		if (mutations.isEmpty()) return new int[] {0, 0};
		
		
		
		// Get child and parent sequences
		Node child = mutations.get(0).getNode();
		//Node parent = child.getParent();
		
		if ( mutations.get(0).getTime() < 0) {
			Log.warning("Unexpected dev error 0435: mutaiton time is -1");
		}
		
		
		// Multispecies coalescent? -- start the count at the point where this gene branch extends into the species node parent
		double maxHeight = child.getParent().getHeight() - mutations.get(0).getTime();
		double minHeight = child.getParent().getHeight() - mutations.get(mutations.size()-1).getTime();
		
		//Log.warning("maxHeight=" + maxHeight + " minHeight=" + minHeight + " t=" +  mutations.get(0).getTime() + " p=" + child.getParent().getHeight());
		
		// Reconstruct sequence at the desired height, in case this is a gene tree within a species tree
		int[] childSequence = MutationUtils.getStatesForBranchAtTime(mapper, mapper.getTree(), child, minHeight, true);
		int[] parentSequence = MutationUtils.getStatesForBranchAtTime(mapper, mapper.getTree(), child, maxHeight, false);
		
		int nNonSyn = 0;
		int nSyn = 0;
		int ncodons = mapper.getSequenceLength(); // Pattern count should equal site count if using PatternlessAlignment
		for (int codonNr = 0; codonNr < ncodons; codonNr ++) {
			
			//Log.warning("Position " + codonNr + " of " + ncodons);
			
			// Find all mutations in this position
			List<Mutation> codonMutations = new ArrayList<>();
			for (Mutation mut : mutations) {
				if (mut.getSiteNr() == codonNr) {
					codonMutations.add(mut);
				}
			}
			
			// Sort by time along branch (forward in time)
			Collections.sort(codonMutations);
			if (codonMutations.isEmpty()) continue;
			
			
			
			// Parent and child states
			int parentCodonNr = parentSequence[codonNr];
			int childCodonNr = childSequence[codonNr];
			
			
			//Log.warning("parentCodonNr " + parentCodonNr + ", childCodonNr " + childCodonNr);
			
			// I don't know what to do about gaps/ambig/stop codons just yet...
			if (parentCodonNr == -1 || childCodonNr == -1) continue;
			
		
			int parentAA = code.getAminoAcidState(parentCodonNr);
			
			
			
			for (Mutation mut : codonMutations) {
				int from = mut.getFrom();
				int to = mut.getTo();
				if (parentCodonNr != from) {
					Log.warning("Unexpected dev error 24244: " + parentCodonNr + "!=" + from);
				}
				parentCodonNr = to;
				
				
				
				// Translate codon. Stop codons, gaps, and ambigs will get -1
				int nextAA = to < 0 ? -1 : code.getAminoAcidState(to);
				
				//Log.warning(from +  " to " + to + " aa " + nextAA);
				
				
				// For debugging
//				DataType aminoacid = new Aminoacid();
//				String fromChar = this.getDataType().getCharacter(from);
//				String toChar = this.getDataType().getCharacter(to);
//				String fromAAChar = aminoacid.getCharacter(parentAA);
//				String toAAChar = aminoacid.getCharacter(nextAA);
				
				
				// Synonymous mutation
				if (nextAA == parentAA) {
					//Log.warning("s: " + fromAAChar + "("+fromChar+")" + "->" + toAAChar + "("+toChar+")" + (" (from " + from + " to " + to + ")"));
					nSyn ++;
				}
				
				// Non-synonymous mutation
				else{
					nNonSyn++;
					//Log.warning("ns: " + fromAAChar + "("+fromChar+")" + "->" + toAAChar + "("+toChar+")" + (" (from " + from + " to " + to + ")"));
				}
				
				parentAA = nextAA;
				
			}
						
						
			
		}
		
		
		return new int[] { nSyn, nNonSyn };
		
	}
		
	
	

	
	public static int[] getSNCountForNucleotides(List<Mutation> mutations, GeneticCode code, Codon codon, int openReadingFrame, StochasticMapper mapper) {

			
		
		if (mutations.isEmpty()) return new int[] {0, 0, 0};
		
		
		// Get child and parent sequences
		Node child = mutations.get(0).getNode();
		
		
		
		// Multispecies coalescent? -- start the count at the point where this gene branch extends into the species node parent
		double maxHeight = child.getParent().getHeight() - mutations.get(0).getTime();
		double minHeight = child.getParent().getHeight() - mutations.get(mutations.size()-1).getTime();
		
		// Reconstruct sequence at the desired height, in case this is a gene tree within a species tree
		int[] childSequence = MutationUtils.getStatesForBranchAtTime(mapper, mapper.getTree(), child, minHeight, true);
		int[] parentSequence = MutationUtils.getStatesForBranchAtTime(mapper, mapper.getTree(), child, maxHeight, false);
		
		
//		Log.warning("ch: " + Arrays.toString(childSequence));
//		Log.warning("pa: " + Arrays.toString(parentSequence));
		
		
		int nNonSyn = 0;
		int nSyn = 0;
		int nStops = 0;
		int nsites = mapper.getSequenceLength(); // Pattern count should equal site count if using PatternlessAlignment
		for (int codonNr = 0; codonNr < nsites/3; codonNr ++) {
			int codon1 = codonNr*3 + openReadingFrame-1;
			int codon2 = codon1+1;
			int codon3 = codon1+2;
			
			if (codon2 >= nsites || codon3 >= nsites) break;
			
			// Find all mutations in this codon
			List<Mutation> codonMutations = new ArrayList<>();
			for (Mutation mut : mutations) {
				if (mut.getSiteNr() == codon1 || mut.getSiteNr() == codon2 || mut.getSiteNr() == codon3) {
					codonMutations.add(mut);
				}
			}
			
			// Sort by time along branch (forward in time)
			Collections.sort(codonMutations);
			if (codonMutations.isEmpty()) continue;
			
			
			
			
			// Parent state
			int[] parentCodon = new int[] { parentSequence[codon1], parentSequence[codon2], parentSequence[codon3] };
			int[] childCodon = new int[] { childSequence[codon1], childSequence[codon2], childSequence[codon3] };
			
			
			
			int parentCodonNr = getCodonNr(parentCodon, codon);
			int childCodonNr = getCodonNr(childCodon, codon);
			
			

			int parentAA = parentCodonNr < 0 ? -1 : code.getAminoAcidState(parentCodonNr);
			
			//Log.warning("Parent " + parentCodon[0] + parentCodon[1] + parentCodon[2] + "->" + parentCodonNr);
			for (Mutation mut : codonMutations) {
				int codonPos = mut.getSiteNr() == codon1 ? 0 : mut.getSiteNr() == codon2 ? 1 : 2;
				int from = mut.getFrom();
				int to = mut.getTo();
				if (parentCodon[codonPos] != from) {
					Log.warning("Unexpected dev error 2424: " + parentCodon[codonPos] + "!=" + from);
				}
				parentCodon[codonPos] = to;
				
				
				// Translate codon. Stop codons, gaps, and ambigs will get -1
				int nextCodonNr = getCodonNr(parentCodon, codon);
				int nextAA = nextCodonNr < 0 ? -1 : code.getAminoAcidState(nextCodonNr);
				if (nextCodonNr < 0 && isStopCodon(parentCodon, codon)) {
					nStops ++;
				}
				
				
				//Log.warning("Mutating from " + from + " to " + to + " at time " + mut.time);
				
				if (nextAA == parentAA) {
					
					// Synonymous mutation
					//Log.warning("s: " + parentAA + "->" + nextAA);
					nSyn ++;
				}else{
					
					// Non-synonymous mutation
					nNonSyn++;
					//Log.warning("ns: " + parentAA + "->" + nextAA);
				}
				
				
				parentAA = nextAA;
				
				//Log.warning("Inter  " + parentCodon[0] + parentCodon[1] + parentCodon[2]);
			}
			//Log.warning("Child  " + childCodon[0] + childCodon[1] + childCodon[2] + "->" + childCodonNr);
			//System.out.println();
			
			if (parentCodon[0] != childCodon[0]) {
				Log.warning("Unexpected dev error 24254.0: " + parentCodon[0] + "!=" + childCodon[0]);
			}
			if (parentCodon[1] != childCodon[1]) {
				Log.warning("Unexpected dev error 24254.1: " + parentCodon[1] + "!=" + childCodon[1]);
			}
			if (parentCodon[2] != childCodon[2]) {
				Log.warning("Unexpected dev error 24254.2: " + parentCodon[2] + "!=" + childCodon[2]);
			}
			
			
		}
		
//		System.out.println("number of stops " + nStops);
//		System.out.println("AAA=" + getCodonNr(new int[] {0,0,0}, codon) + " " + isStopCodon(new int[] {0,0,0}, codon));
//		System.out.println("TGA=" + getCodonNr(new int[] {3,2,0}, codon) + " " + isStopCodon(new int[] {3,2,0}, codon));
//		System.out.println("A--=" + getCodonNr(new int[] {0,-1,-1}, codon) + " " + isStopCodon(new int[] {0,-1,-1}, codon));
		
		return new int[] { nSyn, nNonSyn, nStops };
	}
	
	
	

	
	
	/**
	 * dNdS for a given site
	 * @param codon
	 * @param code
	 * @param mapper
	 * @param codonSiteNr
	 * @return
	 */
	public static double[] getdNdS(Codon codon, GeneticCode code, StochasticMapper mapper, int codonSiteNr, double vN, double vS, double[] frequencies) {
		
		
		boolean wholeSequence = codonSiteNr == -1;

		final int a = 0;
		final int c = 1;
		final int g = 2;
		final int t = 3;
		
		
		// Per branch
		int nS = 0;
		int nN = 0;
		int nTransitions = 0;
		int nTransversions = 0;
		for (int branchNr = 0; branchNr < mapper.getTree().getNodeCount(); branchNr++) {
			List<Mutation> mutations = mapper.getMutationsOnBranch(branchNr);
			
			
			
			// Count number of transitions and transversions from the whole sequence
			for (Mutation mut : mutations) {
				if (mut.getFrom() == a && mut.getTo() == g) nTransitions ++;
				else if (mut.getFrom() == g && mut.getTo() == a) nTransitions ++;
				else if (mut.getFrom() == c && mut.getTo() == t) nTransitions ++;
				else if (mut.getFrom() == t && mut.getTo() == c) nTransitions ++;
				else nTransversions ++;
			}
						
			
			
			// Filter to just one site ?
			if (!wholeSequence) {
				
				
				int pos1 = codonSiteNr*3;
				int pos2 = pos1 + 1;
				int pos3 = pos1 + 2;
				
				List<Mutation> filteredMutations = new ArrayList<Mutation>();
	        	for (Mutation mutation : mutations) {
	        		if (mutation.getSiteNr() == pos1 || mutation.getSiteNr() == pos2 || mutation.getSiteNr() == pos3){
	        			filteredMutations.add(mutation);
	        		}
	        	}
				
	        	mutations = filteredMutations;
				
			}
			
			//Log.warning(branchNr + " " + mutations.size());
			
			// Count nN and nS for the target region only
			int[] counts = getSNCountForNucleotides(mutations, code, codon, 1, mapper);
			nS += counts[0];
			nN += counts[1];
			
			

			
		}
		
		
		
		//if (nN == 0) return new double[] { 0, 0 };
		//if (nS == 0) return new double[] { Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY };
		
		
//		// Estimate the codon frequencies directly from the root sequence
//		int[] rootSequence = mapper.getStatesForNode(mapper.getTree(), mapper.getTree().getRoot());
//		int[] rootCodons = new int[(int) Math.floor(rootSequence.length / 3)];
//		for (int i = 0; i < rootCodons.length; i ++) {
//			int ntPos = i*3;
//			int[] triplet = new int[] { rootSequence[ntPos], rootSequence[ntPos+1], rootSequence[ntPos+2] };
//			int codonNr = MutationUtils.getCodonNr(triplet, codon);
//			rootCodons[i] = codonNr;
//			
//			//System.out.print(i + " " + triplet + " " + codonNr);
//		}
//		
		
//		int sum = 0;
//		double[] frequencies = new double[codon.stateCount];
//		for (int codonNr = 0; codonNr < frequencies.length; codonNr ++) {
//			
//			// How many times do we see this state?
//			double freq = 0;
//			for (int i = 0; i < rootCodons.length; i ++) {
//				if (rootCodons[i] == codonNr) {
//					freq += 1.0;
//					sum += 1;
//				}
//			
//			}
//			
//			frequencies[codonNr] = freq;
//			
//		}
//
//		// Normalise to sum to 1
//		for (int codonNr = 0; codonNr < frequencies.length; codonNr ++) {
//			frequencies[codonNr] = frequencies[codonNr] / sum;
//		}
//		
		
		
		
		
		// Estimate kappa based on frequencies
		double transitionFlux = 0;
		double transversionFlux = 0;
		
		for (int codonFrom = 0; codonFrom < codon.getStateCount(); codonFrom++) {
			for (int codonTo = 0; codonTo < codon.getStateCount(); codonTo++) {
			
				if (codonFrom == codonTo) continue;
				
				int pos = MutationUtils.getMutationPosition(codonFrom, codonTo, codon);
				if (pos == -1) continue;
				
				double rate = frequencies[codonTo];
				if (MutationUtils.isTransition(codonFrom, codonTo, codon)) {
					//Log.warning(codonFrom + " " + codonTo + " is a transition");
					transitionFlux += rate * frequencies[codonFrom];
				}else {
					//Log.warning(codonFrom + " " + codonTo + " is a transversion");
					transversionFlux += rate * frequencies[codonFrom];
				}
			
			}
		}
		
		
		
		// Kappa
		double vT = 4;
		double vV = 1;
		double eTeV = rBetaPrime(vT, vV, nTransitions, nTransversions);
		double phiKappa = transitionFlux / transversionFlux;
		double kappa =  eTeV / phiKappa;
		//Log.warning("We have " + nTransitions  + " ts and " + nTransitions  + " tv and  kappa " + kappa + " transversionFlux=" + transversionFlux + " transitionFlux=" + transitionFlux);
		
		
		
		// Assuming an HKY-like matrix with codon frequencies, we will calculate omega directly from the ratio nN/nS
		double eNeS = rBetaPrime(vN, vS, nN, nS);

		
		//double eNeS = 1.0 * nN / nS;
		double sFlux = 0;
		double nFlux = 0;
		int nSynChanges = 0;
		int nNonsynChanges = 0;
		for (int codonFrom = 0; codonFrom < codon.getStateCount(); codonFrom++) {
			for (int codonTo = 0; codonTo < codon.getStateCount(); codonTo++) {
			
				if (codonFrom == codonTo) continue;
				
				int pos = MutationUtils.getMutationPosition(codonFrom, codonTo, codon);
				if (pos == -1) continue;
				
				double kappaRate = MutationUtils.isTransition(codonFrom, codonTo, codon) ? kappa : 1;
				double rate = frequencies[codonTo] * kappaRate;
				if (MutationUtils.areSynonymous(codonFrom, codonTo, codon)) {
					sFlux += rate * frequencies[codonFrom];
					nSynChanges++;
				}else {
					nFlux += rate * frequencies[codonFrom];
					nNonsynChanges++;
				}
			
			}
		}
		
		double phi = nFlux / sFlux;
		//Log.warning("We have " + nS + " and " + nN + " with " + nTransitions  + " ts and " + nTransversions  + " tv and  kappa " + kappa + " and a gradient of " + phi);
		
		
		double omega_JC = eNeS * nSynChanges / nNonsynChanges;
		double omega_HKY = eNeS / phi;
		
		return new double[] {omega_JC, omega_HKY};
		
		
	}
	

	
	/**
	 * Return dNdS estimated by both phi Jukes Cantor and phi HKY, across the whole sequence
	 * @param codon
	 * @param code
	 * @param mapper
	 * @return
	 */
	public static double[] getdNdS(Codon codon, GeneticCode code, StochasticMapper mapper, double vN, double vS) {
		double[] freqs = getCodonFreqs(mapper, codon, 0.1, -1);
		return getdNdS(codon, code, mapper, -1, vN, vS, freqs);
	}
	
	
	
	public static double[] getCodonFreqs(StochasticMapper mapper, Codon codon, double pseudocount, int codonNr) {

		
		
		
		
		List<Mutation> mutations = new ArrayList<>();
		for (int branchNr = 0; branchNr < mapper.getTree().getNodeCount(); branchNr++) {
			
			List<Mutation> mutationsBranch = mapper.getMutationsOnBranch(branchNr);
			
			// Filter to just one site ?
			if (codonNr >= 0) {
				
				
				int pos1 = codonNr*3;
				int pos2 = pos1 + 1;
				int pos3 = pos1 + 2;
	        	for (Mutation mutation : mutationsBranch) {
	        		if (mutation.getSiteNr() == pos1 || mutation.getSiteNr() == pos2 || mutation.getSiteNr() == pos3){
	        			mutations.add(mutation);
	        		}
	        	}
				
				
			}else {
				mutations.addAll(mutationsBranch);
			}
			
			
			
		}
		
		int nstates = codon.stateCount;
		
		double[][] nchanges = new double[nstates][nstates];
		
		
		// Estimate series of time-non-reversible rates based on how many times we see each change
		for (int i = 0; i < nstates; i ++) {
			for (int j = 0; j < nstates; j ++) {
				nchanges[i][j] = pseudocount; 
			}
		}
		
		for (Mutation mut : mutations) {
			int from = mut.getFrom();
			int to = mut.getTo();
			nchanges[from][to] ++;			
		}
		
		
		// As vector
		List<Double> rateVector = new ArrayList<>();
		for (int i = 0; i < nstates; i ++) {
			for (int j = 0; j < nstates; j ++) {
				if (i == j) continue;
				rateVector.add(nchanges[i][j]);
			}
		}
		

		

		RealParameter r = new RealParameter();
		r.initByName("value", rateVector);
		
		
		prepareFrequencies(nstates);
		ComplexSubstitutionModel model = new ComplexSubstitutionModel();
		model.initByName("rates", r, "frequencies", frequenciesVector);
		
		
		// Find the equilibrium distribution of this Markov chain
		double[] freqs = getFrequencies(model);
		
//		System.out.print("found freqs: ");
//		for (int i = 0; i < nstates; i ++) {
//			System.out.print(freqs[i] + "\t");
//		}
//		System.out.println();
		
		return freqs;
		
		
	}
	
	
	
	/*
	 * Find the equilibrium distribution
	 * From the abyss package
	 */
	private static double[] getFrequencies(GeneralSubstitutionModel model) {
		
		final double DEFAULT_BRANCH_LENGTH = 100000;
		
		
        int nrOfStates = model.getStateCount();
        double t = DEFAULT_BRANCH_LENGTH;
        boolean equilibrium = false;
        double[] p;
        double[] f = new double[nrOfStates];

        while (!equilibrium) {
            p = new double[nrOfStates * nrOfStates];
            model.getTransitionProbabilities(null, 1.0, 0., t, p);
            boolean reached = true;
            for (int i = 0; i < nrOfStates; i++) {
                f[i] = p[i];
                for (int j = 1; j < nrOfStates; j++) {
                    if (p[i] - p[j * nrOfStates + i] > 1e-6) {
                        reached = false;
                        break;
                    }
                }
                if (!reached) break;
            }
            if (reached) equilibrium = true;
            t *= 10;
            
            //Log.warning("t=" + t);
            
        }
        return f;
    }
	
	
	
	/**
	 * Return a random number from a beta prime distribution, based on counts and pseudocounts
	 * @param shape1
	 * @param shape2
	 * @param count1
	 * @param count2
	 * @return
	 */
	public static double rBetaPrime(double shape1, double shape2, int count1, int count2) {
		
		org.apache.commons.math.distribution.BetaDistribution beta = new BetaDistributionImpl(1, 1);
		double x = 0;
		beta.setAlpha(shape1 + count1);
		beta.setBeta(shape2 + count2);
		try {
			double eT = beta.inverseCumulativeProbability(Randomizer.nextDouble());
			x = eT / (1-eT);
		} catch (MathException e) {
			e.printStackTrace();
			x = 1;
		}
		
		return x;
		
	}
	
	
	// Returns the single point where the two codons differ, or -1 if they differ by 0, 2 or 3 positions
	public static int getMutationPosition(int codon1, int codon2, Codon codon) {
		
		String c1 = codon.encodingToString(new int[] { codon1 } );
		String c2 = codon.encodingToString(new int[] { codon2 } );
		String[] bits1 = c1.split("");
		String[] bits2 = c2.split("");
		
		// Where do they differ?
		int differ1 = bits1[0].equals(bits2[0]) ? 0 : 1;
		int differ2 = bits1[1].equals(bits2[1]) ? 0 : 1;
		int differ3 = bits1[2].equals(bits2[2]) ? 0 : 1;
		
		
		if (differ1 + differ2 + differ3 != 1) return -1;
		if (differ1 == 1) return 0;
		if (differ2 == 1) return 1;
		if (differ3 == 1) return 2;
		
		return -2;
		
	}


	public static boolean isTransition(int codon1, int codon2, Codon codon) {
		
		
		
		Nucleotide nt = new Nucleotide();
		
		final int a = 0;
		final int c = 1;
		final int g = 2;
		final int t = 3;
		
		int pos = getMutationPosition(codon1, codon2, codon);
		if (pos < 0) return false;
		
		
		String c1 = codon.encodingToString(new int[] { codon1 } );
		String c2 = codon.encodingToString(new int[] { codon2 } );
		String[] bits1 = c1.split("");
		String[] bits2 = c2.split("");
		
		int nt1 = nt.stringToEncoding(bits1[pos]).get(0);
		int nt2 = nt.stringToEncoding(bits2[pos]).get(0);		
		
		
		//System.out.println(c1 + " " + c2 + " " + pos + " " + bits1[pos] + " " + nt1 + " " + nt2);
		
		if (nt1 == a && nt2 == g) return true;
		if (nt1 == g && nt2 == a) return true;
		if (nt1 == c && nt2 == t) return true;
		if (nt1 == t && nt2 == c) return true;
		return false;
		
	}


	public static boolean areSynonymous(int codon1, int codon2, Codon codon) {
		String aa1 = codon.stateToAminoAcid(new int[] {codon1});
		String aa2 = codon.stateToAminoAcid(new int[] {codon2});
		return aa1.equals(aa2);
	}

	public static boolean isStopCodon(int[] threeNucleotides, Codon codon) {
		if (threeNucleotides[0] < 0 || threeNucleotides[1] < 0 || threeNucleotides[2] < 0) return false; // Gaps or ambig
		Nucleotide nt = new Nucleotide();
		String triplet = nt.encodingToString(threeNucleotides);
		return codon.isTripletStopCodon(triplet);
	}
    
	public static int getCodonNr(int[] threeNucleotides, Codon codon) {
		
		try {
			return codon.getCodonState(threeNucleotides[0], threeNucleotides[1], threeNucleotides[2]);
		} catch (Exception e) {
			return -1;
		}
	}

	
	 

}
