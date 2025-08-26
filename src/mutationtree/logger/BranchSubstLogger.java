package mutationtree.logger;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;
import mutationtree.codons.Codon;
import mutationtree.codons.GeneticCode;
import mutationtree.evolution.BranchMutationSampler;
import mutationtree.util.Mutation;


@Description("Counts the number and type of substitutions on each branch")
public abstract class BranchSubstLogger extends CalculationNode implements Loggable, Function {

	
	final public Input<BranchMutationSampler> samplerInput = new Input<>("sampler", "mutation sampler to log", Validate.REQUIRED);
	


    
    @Override
    public void initAndValidate() {
    	
    	DataType dt = samplerInput.get().dataInput.get().getDataType();
    	if (!this.canHandleDataType(dt)) {
    		throw new IllegalArgumentException(this.getID() + " cannot support datatype " + dt.getClass());
    	}
    	
    }
    
    
    
    /**
     * For example, the number of mutations on the branch
     * @param nodeNr
     * @return
     */
    public abstract double getMutationSummary(List<Mutation> mutations);
    
    
    
    /**
     * The name to log onto the tree
     * @return
     */
	protected abstract String getName();
	
	
	/**
	 * Is this data type supported?
	 * @param dataType
	 * @return
	 */
	protected abstract boolean canHandleDataType(DataType dataType);
    
	
    /**
     * Please call this before calling log or getArrayValue, as the values are sampled once per state at most
     * @param sampleNr
     */
    public void sampleMutations(long sampleNr) {
    	BranchMutationSampler sampler = samplerInput.get();
    	sampler.sampleMutations(sampleNr);
    }


    @Override
    public double getArrayValue(int dim) {
    	
    	// Make sure to call sampleMutations before calling this method so that it depends on the current stochastic sample
    	BranchMutationSampler sampler = samplerInput.get();
    	List<Mutation> mutations = sampler.getMutationsOnBranch(dim);
    	return getMutationSummary(mutations);
    }





	@Override
	public void init(PrintStream out) {
		Tree tree = (Tree) samplerInput.get().treeInput.get();
		for (int i = 0; i < getDimension(); i ++) {
			String id = tree.getNode(i).getID();
			if (id == null) id = "" + i;
			out.print(this.getID() + "." + id + "\t");
		}
		
	}



	@Override
	public void log(long sampleNr, PrintStream out) {
		
		this.sampleMutations(sampleNr);
		
		Tree tree = (Tree) samplerInput.get().treeInput.get();
		for (int nodeNr = 0; nodeNr < getDimension(); nodeNr ++) {
			double nmut = getArrayValue(nodeNr);
			out.print(nmut + "\t");
		}
		
		
	}



	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}
	
    @Override
    public int getDimension() {
    	Tree tree = (Tree) samplerInput.get().treeInput.get();
        return tree.getNodeCount()-1;
    }


    

    /**
     * Return number of synonymous and non-synonymous substitutions in a array of length 2
     * Assuming that the data type is nucleotide
     * @param mutations
     * @param code
     * @param codon
     * @param openReadingFrame
     * @return
     */
	protected int[] getSynonymousAndNonSynonymousSubstitutionCount(List<Mutation> mutations, GeneticCode code, Codon codon, int openReadingFrame) {
		
		
		if (mutations.isEmpty()) return new int[] {0, 0};
		
		// Get child and parent sequences
		BranchMutationSampler sampler = samplerInput.get();
		Node child = mutations.get(0).getNode();
		Node parent = child.getParent();
		int[] childSequence = sampler.getStatesForNode(sampler.treeInput.get(), child);
		int[] parentSequence = sampler.getStatesForNode(sampler.treeInput.get(), parent);
		
		
		int nNonSyn = 0;
		int nSyn = 0;
		int nsites = sampler.dataInput.get().getPatternCount(); // Pattern count should equal site count if using PatternlessAlignment
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
			if (parentCodonNr == -1 || childCodonNr == -1) continue;
			
			int parentAA = code.getAminoAcidState(parentCodonNr);
			
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
	
	
	
	protected int getCodonNr(int[] threeNucleotides, Codon codon) {
		try {
			return codon.getCodonState(threeNucleotides[0], threeNucleotides[1], threeNucleotides[2]);
		} catch (Exception e) {
			return -1;
		}
	}

	





    
}



