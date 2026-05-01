package beastmap.operator;

import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beastmap.evolution.AncestralSequenceTreeLikelihood;

public class LeafOperatorTest extends Operator {

	
	public Input<IntegerParameter> parameterInput = new Input<IntegerParameter>("x", "parameter to operate on", Input.Validate.REQUIRED);
	public Input<AncestralSequenceTreeLikelihood> likelihoodInput = new Input<AncestralSequenceTreeLikelihood>("likelihood", "ancestral sequence tree likelihood for doing stochastic mapping", Input.Validate.REQUIRED);

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public double proposal() {
		
		
		// Ancestral reconstruction
		AncestralSequenceTreeLikelihood likelihood = likelihoodInput.get();
		likelihood.calculateLogP();
		
		
		double[][][] logProbs = likelihood.redrawAncestralStates();
		
		// Print the probs
		for (int leafNr = 0; leafNr < logProbs.length; leafNr ++) {
			
			System.out.println("leaf " + leafNr);
			double[][] leafProbs = logProbs[leafNr];
			for (int patternNr = 0; patternNr < leafProbs.length; patternNr ++) {
				
				
				System.out.print("\tsite " + patternNr + " has log-probs:\t");
				double[] patternProbs = leafProbs[patternNr];
				for (int stateNr = 0; stateNr < patternProbs.length; stateNr ++) {
					System.out.print(patternProbs[stateNr] + "\t");
				}
				System.out.println();
				
			}
			
			
		}
		
		
		return 0;
	}
	
}
