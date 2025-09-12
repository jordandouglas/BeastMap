package beastmap.logger;

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
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.Nucleotide;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;
import beastmap.evolution.BranchMutationSampler;
import beastmap.evolution.RecordedMutationSimulator;
import beastmap.evolution.SimulatedAlignmentWithMutations;
import beastmap.evolution.StochasticMapper;
import beastmap.util.Mutation;
import beastmap.util.MutationUtils;
import codonmodels.evolution.datatype.Codon;
import codonmodels.evolution.datatype.GeneticCode;


@Description("Counts the number and type of substitutions on each branch")
public abstract class BranchSubstLogger extends CalculationNode implements Loggable, Function, StochasticMapProperty {

	
	final public Input<StochasticMapper> samplerInput = new Input<>("sampler", "mutation sampler to log");
	final public Input<RecordedMutationSimulator> truthInput = new Input<>("truth", "mutation sampler to log", Validate.XOR, samplerInput);
	
	final public Input<String> filterInput = new Input<>("filter", "specifies which of the sites in the input alignment we will restrict to" +
	            "First site is 1." +
	            "Filter specs are comma separated, either a singleton, a range [from]-[to] or iteration [from]:[to]:[step]; " +
	            "1-100 defines a range, " +
	            "1-100\3 or 1:100:3 defines every third in range 1-100, " +
	            "1::3,2::3 removes every third site. " +
	            "Default for range [1]-[last site], default for iterator [1]:[last site]:[1]", Validate.OPTIONAL);
	

	 protected List<Integer> filter;
	 
	 
    
    @Override
    public void initAndValidate() {
    	
    	DataType dt = this.getDataType();
    	if (!this.canHandleDataType(dt)) {
    		throw new IllegalArgumentException(this.getID() + " cannot support datatype " + dt.getClass());
    	}
    	
    	this.filter = null;
    	if (filterInput.get() != null && !filterInput.get().isEmpty()) {
    		this.filter = MutationUtils.parseFilterSpec(getSiteAndPatternCount(), filterInput.get());
    	}
    	
    }
    
    

    public double getMutationSummary(List<Mutation> mutations, Node node){
    	
    	if (this.filter == null) {
    		return getFilteredMutationSummary(mutations, node);
    	}else {
    		
    		// Filter the mutations and then give the remainder to the child class
    		List<Mutation> filteredMutations = new ArrayList<Mutation>();
        	for (Mutation mutation : mutations) {
        		if (this.filter.contains(mutation.getSiteNr())){
        			filteredMutations.add(mutation);
        		}
        	}
        	return getFilteredMutationSummary(filteredMutations, node);
    	}
    	
    	
    	
    	
    }
    
    
    @Override
    public Object getPropertyOfNode(Node node) {
    	return getArrayValue(node.getNr());
    }
    
    
    /**
     * For example, the number of mutations on the branch
     * @param nodeNr
     * @return
     */
    public abstract double getFilteredMutationSummary(List<Mutation> mutations, Node node);
    
    
    
    /**
     * The name to log onto the tree / log file
     * @return
     */
	//protected abstract String getName();
	
	
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
	@Override
    public void sampleMutations(long sampleNr) {
    	StochasticMapper sampler = samplerInput.get();
    	if (sampler == null) return;
    	sampler.sampleMutations(sampleNr);
    }


    @Override
    public double getArrayValue(int dim) {
    	
    	
    	List<Mutation> mutations;
    	if (samplerInput.get() != null) {
    		
    		// Make sure to call sampleMutations before calling this method so that it depends on the current stochastic sample
    		mutations = samplerInput.get().getMutationsOnBranch(dim);
    	}else {
    		mutations = truthInput.get().getMutationsOnBranch(dim);
    	}
    	
    	Node node = this.getTree().getNode(dim);
    	return getMutationSummary(mutations, node);
    }





	@Override
	public void init(PrintStream out) {
		Tree tree = getTree();
		for (int i = 0; i < getDimension(); i ++) {
			String id = tree.getNode(i).getID();
			if (id == null) id = "node" + i;
			out.print(this.getID() + "." + id + "\t");
		}
		
	}



	@Override
	public void log(long sampleNr, PrintStream out) {
		
		this.sampleMutations(sampleNr);
		
		for (int nodeNr = 0; nodeNr < getDimension(); nodeNr ++) {
			double nmut = getArrayValue(nodeNr);
			out.print(nmut + "\t");
		}
		
		
	}
	
	protected Tree getTree() {
		if (samplerInput.get() == null) {
			return (Tree) truthInput.get().getTree();
		}else {
			return (Tree) samplerInput.get().getTree();
		}
	}



	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}
	
    @Override
    public int getDimension() {
    	Tree tree = getTree();
        return tree.getNodeCount()-1;
    }


    
    protected int[] getSequenceForNode(Node node) {
    	
    	if (samplerInput.get() == null) {
			return truthInput.get().getSequenceForNode(node);
		}else {
			return samplerInput.get().getStatesForNode(getTree(), node);
		}
    	
    	
    }
    
    
    protected DataType getDataType(){
    	return samplerInput.get() != null ? samplerInput.get().getDataType() : truthInput.get().getDataTypeOfSimulator();
    }
    
    
    
    /**
     * We are assuming that the number of patterns is the number of sites, which will be true if a PatternlessAlignment is used
     * @return
     */
    protected int getSiteAndPatternCount(){
    	int siteCount = samplerInput.get() != null ? samplerInput.get().getPatternCount() : truthInput.get().getData().getSiteCount();
    	return siteCount;
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
		
		if (this.getDataType() instanceof Nucleotide) {
			return getSNCountForNucleotides(mutations, code, codon, openReadingFrame);
		}
		
		if (this.getDataType() instanceof Codon) {
			return getSNCountForCodons(mutations, code, codon);
		}
		
		return null;
		
	}
	
	
	// Number of synonymous and non-synonymous, when using codon data
	private int[] getSNCountForCodons(List<Mutation> mutations, GeneticCode code, Codon codon) {
		
		
		if (mutations.isEmpty()) return new int[] {0, 0};
		
		
		
		// Get child and parent sequences
		Node child = mutations.get(0).getNode();
		Node parent = child.getParent();
		int[] childSequence = getSequenceForNode(child);
		int[] parentSequence = getSequenceForNode(parent);
		
		int nNonSyn = 0;
		int nSyn = 0;
		int ncodons = this.getSiteAndPatternCount(); // Pattern count should equal site count if using PatternlessAlignment
		for (int codonNr = 0; codonNr < ncodons; codonNr ++) {
			
			//Log.warning("Position " + codonNr + " of " + ncodons);
			
			// Find all mutations in this position
			List<Mutation> codonMutations = new ArrayList<>();
			for (Mutation mut : mutations) {
				if (mut.getSiteNr() == codonNr) {
					codonMutations.add(mut);
				}
			}
			
			// Sort by time along branch (forward in time)
			Collections.sort(codonMutations);
			if (codonMutations.isEmpty()) continue;
			
			
			
			// Parent and child states
			int parentCodonNr = parentSequence[codonNr];
			int childCodonNr = childSequence[codonNr];
			
			
			//Log.warning("parentCodonNr " + parentCodonNr + ", childCodonNr " + childCodonNr);
			
			// I don't know what to do about gaps/ambig/stop codons just yet...
			if (parentCodonNr == -1 || childCodonNr == -1) continue;
			
		
			int parentAA = code.getAminoAcidState(parentCodonNr);
			
			
			
			for (Mutation mut : codonMutations) {
				int from = mut.getFrom();
				int to = mut.getTo();
				if (parentCodonNr != from) {
					Log.warning("Unexpected dev error 24244: " + parentCodonNr + "!=" + from);
				}
				parentCodonNr = to;
				
				
				
				// Translate codon. Stop codons, gaps, and ambigs will get -1
				int nextAA = to < 0 ? -1 : code.getAminoAcidState(to);
				
				//Log.warning(from +  " to " + to + " aa " + nextAA);
				
				
				// For debugging
//				DataType aminoacid = new Aminoacid();
//				String fromChar = this.getDataType().getCharacter(from);
//				String toChar = this.getDataType().getCharacter(to);
//				String fromAAChar = aminoacid.getCharacter(parentAA);
//				String toAAChar = aminoacid.getCharacter(nextAA);
				
				
				// Synonymous mutation
				if (nextAA == parentAA) {
					//Log.warning("s: " + fromAAChar + "("+fromChar+")" + "->" + toAAChar + "("+toChar+")" + (" (from " + from + " to " + to + ")"));
					nSyn ++;
				}
				
				// Non-synonymous mutation
				else{
					nNonSyn++;
					//Log.warning("ns: " + fromAAChar + "("+fromChar+")" + "->" + toAAChar + "("+toChar+")" + (" (from " + from + " to " + to + ")"));
				}
				
				parentAA = nextAA;
				
			}
						
						
			
		}
		
		return new int[] { nSyn, nNonSyn };
		
	}
		
	
	// Number of synonymous and non-synonymous, when using nucleotide data
	private int[] getSNCountForNucleotides(List<Mutation> mutations, GeneticCode code, Codon codon, int openReadingFrame) {
			
		if (mutations.isEmpty()) return new int[] {0, 0};
		
		// Get child and parent sequences
		Node child = mutations.get(0).getNode();
		Node parent = child.getParent();
		int[] childSequence = getSequenceForNode(child);
		int[] parentSequence = getSequenceForNode(parent);
		
		
		int nNonSyn = 0;
		int nSyn = 0;
		int nsites = this.getSiteAndPatternCount(); // Pattern count should equal site count if using PatternlessAlignment
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



