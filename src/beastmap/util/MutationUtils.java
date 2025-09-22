package beastmap.util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import beastmap.evolution.StochasticMapper;
import codonmodels.evolution.datatype.Codon;
import codonmodels.evolution.datatype.GeneticCode;

public class MutationUtils {
	
	
	final static double EPSILON = 1e-6; // The small amount of time considered as 1 dt for computing the instantaneous transition rates from a P matrix
	
	
	/**
	 * Many substitution model implementations do not return the Q matrix, and some of these are node/time dependent
	 * The most general approach that works around the setup of beast2 is to calculate P(Qt)/t for a very small t, whose off-diagonal elements will be a close approximation of Q
	 * Note that the diagonal element will be 0 not the negative row sum
	 * @param substModel
	 * @param fromState
	 * @param node
	 * @return
	 */
	public static double[] getTransitionRates(SubstitutionModel substModel, int fromState, Node node) {
		
		
		
		int nstates = substModel.getStateCount();
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
	
	
	public static double[][] getTransitionRates(SubstitutionModel substModel, Node node, boolean nodeDependent) {
		
		
		
		if (substModel instanceof GeneralSubstitutionModel && !nodeDependent) {
			GeneralSubstitutionModel gensub = (GeneralSubstitutionModel)substModel;
			gensub.setupRelativeRates();
			gensub.setupRateMatrix();
			return gensub.getRateMatrix();
		}
		
		
		int nstates = substModel.getStateCount();
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
    public static int simulateMutationsDownBranch(int parentState, Node node, double clockRate, SubstitutionModel qmatrix, List<Mutation> arr, int siteNr) {
    	

    	
		int from = parentState;
		double t = 0;
		double time = node.getLength();
		while (true) {
			 

			// Outgoing rate
			double lambda = 0;
			double[] outRates = MutationUtils.getTransitionRates(qmatrix, from, node);
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

    
    

	
	// Number of synonymous and non-synonymous, when using codon data
	public static int[] getSNCountForCodons(List<Mutation> mutations, GeneticCode code, Codon codon, StochasticMapper mapper) {
		
		
		if (mutations.isEmpty()) return new int[] {0, 0};
		
		
		
		// Get child and parent sequences
		Node child = mutations.get(0).getNode();
		Node parent = child.getParent();
		int[] childSequence = mapper.getStatesForNode(mapper.getTree(), child);
		int[] parentSequence = mapper.getStatesForNode(mapper.getTree(), parent);
		
		int nNonSyn = 0;
		int nSyn = 0;
		int ncodons = mapper.getPatternCount(); // Pattern count should equal site count if using PatternlessAlignment
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
		
	
	

	
	// Number of synonymous and non-synonymous, when using nucleotide data
	public static int[] getSNCountForNucleotides(List<Mutation> mutations, GeneticCode code, Codon codon, int openReadingFrame, StochasticMapper mapper) {
			
		if (mutations.isEmpty()) return new int[] {0, 0};
		
		// Get child and parent sequences
		Node child = mutations.get(0).getNode();
		Node parent = child.getParent();
		int[] childSequence = mapper.getStatesForNode(mapper.getTree(), child);
		int[] parentSequence = mapper.getStatesForNode(mapper.getTree(), parent);
		
		
		int nNonSyn = 0;
		int nSyn = 0;
		int nsites = mapper.getPatternCount(); // Pattern count should equal site count if using PatternlessAlignment
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
			
			
			
			// I don't know what to do about gaps/ambig/stop codons just yet...
			//if (parentCodonNr == -1 || childCodonNr == -1) continue;
			
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
		
		
		return new int[] { nSyn, nNonSyn };
	}
	
	
    
	private static int getCodonNr(int[] threeNucleotides, Codon codon) {
		try {
			return codon.getCodonState(threeNucleotides[0], threeNucleotides[1], threeNucleotides[2]);
		} catch (Exception e) {
			return -1;
		}
	}



	 

}
