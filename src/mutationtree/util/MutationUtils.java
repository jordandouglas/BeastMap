package mutationtree.util;

import beast.base.core.Log;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;

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
	 

}
