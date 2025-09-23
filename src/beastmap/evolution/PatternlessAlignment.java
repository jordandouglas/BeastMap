package beastmap.evolution;


import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.FilteredAlignment;
import beast.base.evolution.alignment.Alignment.SiteComparator;
import beast.base.evolution.datatype.DataType;

@Description("An alignment for which each site is its own pattern (useful for ancestral stae reconstruction with site categories)")
public class PatternlessAlignment extends FilteredAlignment {
	
	
	//final public Input<Alignment> dataInput = new Input<>("data", "The alignment.", Input.Validate.REQUIRED);
	
	int[] filter;
    int[] from;
    int[] to;
    int[] step;
    
    
    
    
    
    
    public PatternlessAlignment() {
    	filterInput.setRule(Validate.OPTIONAL);
    	sequenceInput.setRule(Validate.OPTIONAL);
        siteWeightsInput.setRule(Validate.FORBIDDEN);
    }
    

	
	@Override
    public void initAndValidate() {
		
		
		String filterStr = "1-" + alignmentInput.get().getSiteCount();
//		if (alignmentInput.get() instanceof FilteredAlignment) {
//			FilteredAlignment aln = (FilteredAlignment) alignmentInput.get();
//			filterStr = aln.filterInput.get();
//		}
		
		if (alignmentInput.get().userDataTypeInput.get() != null) {
			this.setInputValue(userDataTypeInput.getName(), alignmentInput.get().userDataTypeInput.get());
		}
		
		
		for (String input : alignmentInput.get().getInputs().keySet()) {
			Object value = alignmentInput.get().getInputValue(input);
			if (input.equals(constantSiteWeightsInput.getName())) continue;
			if (input.equals(filterInput.getName())) continue;
			if (input.equals(alignmentInput.getName())) continue;
			if (input.equals(sequenceInput.getName())) continue;
			if (input.equals(taxonSetInput.getName())) continue;
			this.setInputValue(input, value);
		}
		
		
		Log.warning("filter " + filterStr + " " + alignmentInput.get().getTaxonCount());
		this.setInputValue("filter", filterStr);
		
		parseFilterSpec();
		calcFilter();
		super.initAndValidate();
		
	}
		
	
	
    private void parseFilterSpec() {
        // parse filter specification
        String filterString = filterInput.get();
        String[] filters = filterString.split(",");
        from = new int[filters.length];
        to = new int[filters.length];
        step = new int[filters.length];
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
                to[i] = parseInt(strs[1], alignmentInput.get().getSiteCount()) - 1;
            } else if (filterString.matches(".*:.*:.+")) {
                // iterator, e.g. 1:100:3
                String[] strs = filterString.split(":");
                from[i] = parseInt(strs[0], 1) - 1;
                to[i] = parseInt(strs[1], alignmentInput.get().getSiteCount()) - 1;
                step[i] = parseInt(strs[2], 1);
            } else if (filterString.trim().matches("[0-9]*")) {
                from[i] = parseInt(filterString.trim(), 1) - 1;
                to[i] = from[i];
            	step[i] = 1;
            } else {
                throw new IllegalArgumentException("Don't know how to parse filter " + filterString);
            }
        }
    }
    
    private int parseInt(String str, int defaultValue) {
        str = str.replaceAll("\\s+", "");
        try {
            return Integer.parseInt(str);
        } catch (Exception e) {
            return defaultValue;
        }
    }


	
	
	private void calcFilter() {
		
        boolean[] isUsed = new boolean[alignmentInput.get().getSiteCount()];
        
       // Log.warning("isUsed " + isUsed.length);
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
        filter = new int[k];
        k = 0;
        for (int i = 0; i < isUsed.length; i++) {
            if (isUsed[i]) {
                filter[k++] = i;
            }
        }
    }
	

	
	@Override
	protected void calcPatterns() {
		
		//super.calcPatterns();
		
		int taxonCount = counts.size();
        int siteCount = filter.length;
        
        
        Log.warning("initialising patternless alignment with " +  taxonCount + " taxa and " + siteCount + " sites ");
        
        //Log.warning(dataInput.get().getID() + " " + dataInput.get().getClass());
        
        
        

        // convert data to transposed int array
        int[][] data = new int[siteCount][taxonCount];
        for (int i = 0; i < taxonCount; i++) {
            List<Integer> sites = counts.get(i);
            for (int j = 0; j < siteCount; j++) {
                data[j][i] = sites.get(filter[j]);
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
