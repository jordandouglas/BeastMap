package beastmap.logger.mut;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.Nucleotide;
import beast.base.evolution.tree.Node;
import beastmap.logger.BranchSubstLogger;
import beastmap.util.Mutation;
import beastmap.util.MutationUtils;
import codonmodels.evolution.datatype.Codon;
import codonmodels.evolution.datatype.GeneticCode;


@Description("Estimate of dN/dS for a branch")
public class Estimate_dNdS extends BranchSubstLogger {

	final public Input<Integer> readingFrameInput = new Input<>("readingFrame", "reading frame in position 1 2 or 3?", 1, new Integer[] { 1,2,3 });
	final public Input<String> geneticCodeInput = new Input<>("code", "name of genetic code", "universal", GeneticCode.GENETIC_CODE_NAMES);
	
	
	final int a = 0;
	final int c = 1;
	final int g = 2;
	final int t = 3;
	
	int frame;
	GeneticCode code;
	Codon codon;
	
	
	final double pSpN = 1.0 * 138/438; // Number of possible synonymous vs non-synonymous changes, stop codons included
	
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
		
		
		// Count nN and nS
		int[] counts = super.getSynonymousAndNonSynonymousSubstitutionCount(mutations, code, codon, frame);
		int nS = counts[0];
		int nN = counts[1];
		if (nS == 0) return Double.POSITIVE_INFINITY; 
		
		
		// Count number of transitions
		int nTransitions = 0;
		for (Mutation mut : mutations) {
			if (mut.getFrom() == a && mut.getTo() == g) nTransitions ++;
			if (mut.getFrom() == g && mut.getTo() == a) nTransitions ++;
			if (mut.getFrom() == c && mut.getTo() == t) nTransitions ++;
			if (mut.getFrom() == t && mut.getTo() == c) nTransitions ++;
		}
		
		
		// Estimate the codon frequencies directly from the root sequence
		int[] rootSequence = super.getSequenceForNode(super.getTree().getRoot());
		int[] rootCodons = new int[(int) Math.floor(rootSequence.length / 3)];
		for (int i = 0; i < rootCodons.length; i ++) {
			int ntPos = i*3;
			int[] triplet = new int[] { rootSequence[ntPos], rootSequence[ntPos+1], rootSequence[ntPos+2] };
			int codonNr = MutationUtils.getCodonNr(triplet, codon);
			rootCodons[i] = codonNr;
			
			//System.out.print(i + " " + triplet + " " + codonNr);
		}
		
		
		int sum = 0;
		double[] frequencies = new double[codon.stateCount];
		for (int codonNr = 0; codonNr < frequencies.length; codonNr ++) {
			
			// How many times do we see this state?
			double freq = 0;
			for (int i = 0; i < rootCodons.length; i ++) {
				if (rootCodons[i] == codonNr) {
					freq += 1.0;
					sum += 1;
				}
			
			}
			
			frequencies[codonNr] = freq;
			
		}

		// Normalise to sum to 1
		for (int codonNr = 0; codonNr < frequencies.length; codonNr ++) {
			frequencies[codonNr] = frequencies[codonNr] / sum;
		}
		
		
		
		// Log.warning("We have " + nS + " and " + nN + " with " + nTransitions + " transitions " + frequencies[0] + frequencies[1]);
		
		return nN/nS*pSpN;
		
	}

	@Override
	public String getName() {
		return "syn" + "." + this.getID();
	}
	
	@Override
	public boolean canHandleDataType(DataType dataType) {
		return dataType instanceof Nucleotide;
	}
	
	
	@Override
	public boolean isSummable() {
		return false;
	}

}
