package beastmap.logger.mut;

import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.Nucleotide;
import beast.base.evolution.tree.Node;
import beastmap.logger.BranchSubstLogger;
import beastmap.util.Mutation;
import codonmodels.evolution.datatype.Codon;
import codonmodels.evolution.datatype.GeneticCode;


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
		
		
    }
	
	
	
	@Override
	public double getFilteredMutationSummary(List<Mutation> mutations, Node node) {
		int[] counts = super.getSynonymousAndNonSynonymousSubstitutionCount(mutations, code, codon, frame);
		return counts[1];
	}
	

	@Override
	public String getName() {
		return "nonsyn" + "." + this.getID();
	}



	@Override
	protected boolean canHandleDataType(DataType dataType) {
		return dataType instanceof Nucleotide || dataType instanceof Codon;
	}

}
