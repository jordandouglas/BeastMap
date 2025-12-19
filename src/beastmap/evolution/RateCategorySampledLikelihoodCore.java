package beastmap.evolution;

import beast.base.core.Log;
import beast.base.evolution.likelihood.BeerLikelihoodCore;
import beast.base.util.Randomizer;

/**
 * Likelihood core which samples a gamma category rather than integrating over them
 * @author jdou557
 *
 */
public class RateCategorySampledLikelihoodCore extends BeerLikelihoodCore {

	
	protected double[] rootFreqs = null;
	protected double[] siteCategoryPosteriors = null;
	protected int[] patternSiteCategories;
	
	public RateCategorySampledLikelihoodCore(int nrOfStates) {
		super(nrOfStates);
	}
	
	@Override
	public void initialize(int nodeCount, int patternCount, int matrixCount, boolean integrateCategories, boolean useAmbiguities) {
		super.initialize(nodeCount, patternCount, matrixCount, integrateCategories, useAmbiguities);
		siteCategoryPosteriors = new double[matrixCount];
		patternSiteCategories = new int[patternCount];
	}
	
	
	/**
	 * This must be called before 'integratePartials'
	 * @param freqs
	 */
	public void setRootFrequencies(double[] rootFreqs) {
		this.rootFreqs = rootFreqs;
	}
	
	
	@Override
	public void integratePartials(int nodeIndex, double[] proportions, double[] outPartials) {
        calculateAndSamplePartials(partials[currentPartialsIndex[nodeIndex]][nodeIndex], proportions, outPartials);
    }
	
	
	public int getSampledSiteCategory(int patternNum) {
		return patternSiteCategories[patternNum];
	}
	
	/**
     * Samples a category for each pattern (make sure that each pattern = 1 site by using the right alignment class)
     * The probability of selecting a pattern is equal to the site partial * root freq * proportion
     */
	protected void calculateAndSamplePartials(double[] inPartials, double[] proportions, double[] outPartials) {

		
		//Log.warning("nrOfPatterns=" + nrOfPatterns + ",  nrOfMatrices=" + nrOfMatrices + " nrOfStates=" + nrOfStates);
		
        // Sample partial for each pattern from the categories
        for (int k = 0; k < nrOfPatterns; k++) {
        	
        	
        	double sum = 0;
        	for (int l = 0; l < nrOfMatrices; l++) {
        		
        		// What is the site partial for this category?
        		siteCategoryPosteriors[l] = 0;
        		for (int i = 0; i < nrOfStates; i++) {
        			int index = l*this.nrOfStates*this.nrOfPatterns + k*this.nrOfStates + i;
        			siteCategoryPosteriors[l] += this.rootFreqs[i] * inPartials[index];
        		}
        		siteCategoryPosteriors[l] = siteCategoryPosteriors[l] * proportions[l];
        		
        		sum += siteCategoryPosteriors[l];
        		
        	}
        	
        	
        	
        	// Sample a category
        	int siteCategory;
        	try {
        		siteCategory = sum <= 0 ? 0 : Randomizer.randomChoicePDF(siteCategoryPosteriors);
        	}catch(Error e) {
        		
        		// Use MAP state
        		double max = siteCategoryPosteriors[0];
                int choice = 0;
                for (int i = 1; i < siteCategoryPosteriors.length; i++) {
                    if ((siteCategoryPosteriors[i] - max)/(siteCategoryPosteriors[i] + max) > 1e-10) {
                        max = siteCategoryPosteriors[i];
                        choice = i;
                    }
                }
                siteCategory = choice;
        		
        	}
        	
        	
        	
        	patternSiteCategories[k] = siteCategory;
        	
        	
        	
        	//Log.warning("Sampled cat " + siteCategory + " for site " + k);
        	
        	
        	// Partial
        	for (int i = 0; i < nrOfStates; i++) {
    			int index = siteCategory*this.nrOfStates*this.nrOfPatterns + k*this.nrOfStates + i;
    			outPartials[k*this.nrOfStates + i] = inPartials[index];
    		}
        	
        	
        }
        
        
        
    }


}
