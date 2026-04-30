package beastmap.evolution;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.inference.parameter.RealParameter;
import codonmodels.GeneralCodonFrequencies;
import codonmodels.evolution.alignment.CodonAlignment;


@Description("Convert a vector of 20 aa frequencies into 61 codons")
public class AAToCodonFrequencies extends GeneralCodonFrequencies {
	
	
	final public Input<RealParameter> aaFreqsInput = new Input<>("aaFreqs","20 aa frequencies", Input.Validate.REQUIRED);
	final public Input<RealParameter> codonFreqsInput = new Input<>("codonFreqs","61 baseline codon frequencies", Input.Validate.REQUIRED);
	
	
	@Override
    public void initAndValidate() {
		
        Alignment data = null;
        if (dataInput.get() != null && dataInput.get() instanceof CodonAlignment) {
        	data = dataInput.get();
        } else {
        	data = getData(this);
        }
		
		codonDataType = CodonAlignment.toCodonAlignment(data).getDataType();
        codonStateCount = codonDataType.getStateCount();
        needsUpdate = true;
		
		if (codonFreqsInput.get() == null || codonFreqsInput.get().getDimension() != codonStateCount) {
        	throw new IllegalArgumentException("Must provide " + codonStateCount + " codon frequencies");
        }
		
		if (aaFreqsInput.get().getDimension() != 20) {
        	throw new IllegalArgumentException("Must provide " + 20 + " aa frequencies");
        }

        
     	
   }
	
	
    private Alignment getData(BEASTInterface o2) {
		if (o2 instanceof GenericTreeLikelihood) {
			GenericTreeLikelihood tl = (GenericTreeLikelihood) o2;
			return tl.dataInput.get();
		}
		
    	for (BEASTInterface o : o2.getOutputs()) {
    		Alignment data = getData(o);
    		if (data != null) {
    			return data;
    		}
    	}
    	return null;
	}
	
	
	@Override
    protected void update() {
        // 60/61 codon frequencies (separated by white space) are in fixed order AAA AAC AAG AAT ... TTA TTC TTG TTT
        RealParameter codonFreqs = codonFreqsInput.get();
        RealParameter aaFreqs = aaFreqsInput.get();
        
        
        Aminoacid aminoAcids = new Aminoacid();
        
        freqs = new double[codonFreqs.getDimension()];
       
       
        // Multiply each coding class by its amino acid frequency
        for (int aa = 0; aa < 20; aa++) {
        	
        	String aaStr = aminoAcids.encodingToString(new int[] {aa});
        	double aaFreq = aaFreqs.getArrayValue(aa);
        	
        	
        	// What codons have this amino acid?
        	double doubleClassFreq = 0;
        	for (int i = 0; i < freqs.length; i++) {
        		String translated = codonDataType.stateToAminoAcid(new int[] { i });
        		boolean codonMapsToAA = translated.equals(aaStr);
        		if (codonMapsToAA) {
        			//Log.warning(aa + " " + aaStr + " is coded by " + codonDataType.encodingToString(new int[] {i}) +  " freq " + codonFreqs.getArrayValue(i));
        			doubleClassFreq += codonFreqs.getArrayValue(i);
        		}
    		}
        	
        	
        	for (int i = 0; i < freqs.length; i++) {
        		
        		
        		String translated = codonDataType.stateToAminoAcid(new int[] { i });
        		boolean codonMapsToAA = translated.equals(aaStr);
        		if (codonMapsToAA) {
        			
        			//Log.warning(aa + " " + aaStr + " is coded by " + codonDataType.encodingToString(new int[] {i}) +  " rate " + aaFreq);
//        			Log.warning("codon freqs " + codonFreqs.getArrayValue(i) + " " + doubleClassFreq);
//        			
        			freqs[i] = codonFreqs.getArrayValue(i) * aaFreq/doubleClassFreq;
        		}
                
            }
        	
        	
        }
        
        
        // Normalise so they sum to 1
        double sum = 0;
        for (int i = 0; i < freqs.length; i++) {
            sum += freqs[i];
        }
        
        
        for (int i = 0; i < freqs.length; i++) {
            freqs[i] = freqs[i] / sum;
           // Log.warning(codonDataType.encodingToString(new int[] {i}) + " has freq " + freqs[i]);
        }
        

    	needsUpdate = true;
        
        
	}

}
