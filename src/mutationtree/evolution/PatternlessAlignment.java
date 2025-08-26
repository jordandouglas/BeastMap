package mutationtree.evolution;


import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;

@Description("An alignment for which each site is its own pattern (useful for ancestral stae reconstruction with site categories)")
public class PatternlessAlignment extends Alignment {
	
	
	final public Input<Alignment> dataInput = new Input<>("data", "The alignment.", Input.Validate.REQUIRED);
	
	
	
	@Override
    public void initAndValidate() {
		for (String input : dataInput.get().getInputs().keySet()) {
			this.setInputValue(input, dataInput.get().getInputValue(input));
		}
		super.initAndValidate();
	}
	
	@Override
	protected void calcPatterns(boolean log) {
		
		
		int taxonCount = counts.size();
        int siteCount = counts.get(0).size();
        
        
        Log.warning("initialising patternless alignment with " +  taxonCount + " taxa and " + siteCount + " sites");;

        // convert data to transposed int array
        int[][] data = new int[siteCount][taxonCount];
        for (int i = 0; i < taxonCount; i++) {
            List<Integer> sites = counts.get(i);
            for (int j = 0; j < siteCount; j++) {
                data[j][i] = sites.get(j);
            }
        }
		
        
        // Each pattern is the site 
        sitePatterns = new int[siteCount][taxonCount];
        for (int i = 0; i < siteCount; i++) {
            sitePatterns[i] = data[i];
        }
        
        
        // All weights should be equal
        patternWeight = new int[siteCount];
        for (int i = 0; i < siteCount; i++) {
        	patternWeight[i] = 1;
        }
        
        
        // All site indices point to themselves
        patternIndex = new int[siteCount];
        for (int i = 0; i < siteCount; i++) {
        	patternIndex[i] = i;
        }
        
        
        maxStateCount = 0;
        for (int m_nStateCount1 : stateCounts) {
            maxStateCount = Math.max(maxStateCount, m_nStateCount1);
        }
        
        
        
	}
	 
	 
	 

}
