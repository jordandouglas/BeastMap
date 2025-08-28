package mutationtree.logger.mut;

import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.Nucleotide;
import mutationtree.codons.Codon;
import mutationtree.codons.GeneticCode;
import mutationtree.logger.BranchSubstLogger;
import mutationtree.util.Mutation;


@Description("Number of non-synonymous mutations")
public class NonSynonymousSubstSum extends BranchSubstLogger {

	
	final public Input<Integer> readingFrameInput = new Input<>("readingFrame", "reading frame in position 1 2 or 3?", 1, new Integer[] { 1,2,3 });
	final public Input<String> geneticCodeInput = new Input<>("code", "name of genetic code", "universal", GeneticCode.GENETIC_CODE_NAMES);
	
	
	int frame;
	GeneticCode code;
	Codon codon;
	
	@Override
    public void initAndValidate() {
		
		super.initAndValidate();
		
		
		// Parse inputs
		this.code = GeneticCode.findByName(geneticCodeInput.get());
		if (this.code == null) {
			throw new IllegalArgumentException("Cannot find code " + geneticCodeInput.get());
		}
		this.codon = new Codon(code);
		this.frame = readingFrameInput.get();
		
		
		// If there is a filter, make sure none of the codons have been split into bits
//		if (this.filter != null) {
//			
//			for (int i = 0; i < filter.length; i ++){
//				
//				
//				
//			}
//			
//		}
		
		
		
    }
	
	
	
	@Override
	public double getFilteredMutationSummary(List<Mutation> mutations) {
		int[] counts = super.getSynonymousAndNonSynonymousSubstitutionCount(mutations, code, codon, frame);
		return counts[1];
	}
	

	@Override
	protected String getName() {
		return "nonsyn";
	}



	@Override
	protected boolean canHandleDataType(DataType dataType) {
		return dataType instanceof Nucleotide;
	}

}
