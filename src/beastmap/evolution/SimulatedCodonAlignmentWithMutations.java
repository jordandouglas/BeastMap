package beastmap.evolution;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.parser.XMLProducer;
import beast.base.util.Randomizer;
import beastmap.util.Mutation;
import beastmap.util.MutationUtils;
import codonmodels.evolution.alignment.CodonAlignment;
import codonmodels.evolution.datatype.Codon;
import codonmodels.evolution.datatype.GeneticCode;



/**
 * Adapated from the SimulatedCodonAlignment class in the codonsubstmodels package
 */
@Description("An alignment containing sequences randomly generated using a"
        + "given site model down a given tree.")
public class SimulatedCodonAlignmentWithMutations extends CodonAlignment implements StochasticMapper {
    final public Input<Tree> m_treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);
    
    final public Input<List<SiteModel.Base>> m_pSiteModelInput = new Input<>("siteModel", "site model for leafs in the beast.tree", new ArrayList<>());
    final public Input<IntegerParameter> siteModelNumberInput = new Input<>("siteModelNumber", "the site model number of each site will be saved to this list", Validate.OPTIONAL);
    
    
    final public Input<IntegerParameter> siteModelChangesInput = new Input<>("changes", "a list of internal node numbers where the site model changes, same length as the siteModel list", Validate.OPTIONAL);
    
    final public Input<BranchRateModel.Base> m_pBranchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");
    final public Input<Integer> m_sequenceLengthInput = new Input<>("sequencelength", "nr of samples to generate (default 1000).", 1000);
    final public Input<String> m_outputFileNameInput = new Input<>(
            "outputFileName",
            "If provided, simulated alignment is additionally written to this file.");    
    
    final public Input<Boolean> oneSiteModelPerSiteInput = new Input<>("oneSiteModelPerSite", "if there are multiple site models, do we sample one per site?", false);
    

    /**
     * nr of samples to generate *
     */
    protected int m_sequenceLength;
    /**
     * tree used for generating samples *
     */
    protected Tree m_tree;
    /**
     * site model used for generating samples *
     */
    //protected SiteModel.Base m_siteModel;
    /**
     * branch rate model used for generating samples *
     */
    protected BranchRateModel m_branchRateModel;
    /**
     * nr of categories in site model *
     */
    int m_categoryCount;
    /**
     * nr of states in site model *
     */
    int m_stateCount;

    /**
     * name of output file *
     */
    String m_outputFileName;    

    /**
     * an array used to transfer transition probabilities
     */
    protected double[][] m_probabilities;
    
    long lastSampleNr;
    boolean oneSiteModelPerSite;
    
    int[][] sequencesAll;
    List<List<Mutation>> mutationsAlongEachBranch;
    
    public SimulatedCodonAlignmentWithMutations() {
        
        // Override the sequence input requirement.
        sequenceInput.setRule(Validate.OPTIONAL);
    }
    
    

    @Override
    public void initAndValidate() {
    	
    	mutationsAlongEachBranch = new ArrayList<>();
        m_tree = m_treeInput.get();
        
        m_branchRateModel = m_pBranchRateModelInput.get();
        m_sequenceLength = m_sequenceLengthInput.get();
        
        if (m_pSiteModelInput.get().isEmpty()) {
        	throw new IllegalArgumentException("Please provide a site model");
        }
        
        oneSiteModelPerSite = oneSiteModelPerSiteInput.get();
        if (!oneSiteModelPerSite) {
        	
	        if (m_pSiteModelInput.get().size() > 1 && (siteModelChangesInput.get() == null || m_pSiteModelInput.get().size() != siteModelChangesInput.get().getDimension())) {
	        	throw new IllegalArgumentException("The number of site models must be at the same as the number of changes");
	        }
        
	        if (m_pSiteModelInput.get().size() > 1) {
	        	int rootNr = m_tree.getRoot().getNr();
	        	
	        	boolean rootIncluded = false;
	        	for (int i = 0; i < siteModelChangesInput.get().getDimension(); i ++) {
	        		int nr = siteModelChangesInput.get().getNativeValue(i);
	        		if (nr == rootNr) {
	        			rootIncluded = true;
	        			break;
	        		}
	        	}
	        	if (!rootIncluded) {
	        		throw new IllegalArgumentException("One of the changes must be the root number " + rootNr);
	        	}
	        	
	    	}
        
    	}else {
    		if (siteModelNumberInput.get() != null) {
    			siteModelNumberInput.get().setDimension(m_sequenceLengthInput.get());
    		}
    	}

        GeneticCode geneticCode = GeneticCode.findByName(geneticCodeInput.get());
        setGeneticCode(geneticCode);

        m_stateCount = geneticCode.getStateCount();
        m_categoryCount = this.getSiteModel(m_tree.getRoot().getNr(), 0).getCategoryCount();
        m_probabilities = new double[m_categoryCount][m_stateCount * m_stateCount];
        m_outputFileName = m_outputFileNameInput.get();
        
        
        sequencesAll = new int[m_tree.getNodeCount()][];
        sequenceInput.get().clear();
        
        simulate();        
        
        // Write simulated alignment to disk if requested:
        if (m_outputFileName != null) {
            PrintStream pstream;
			try {
				pstream = new PrintStream(m_outputFileName);
	            pstream.println(new XMLProducer().toRawXML(this));
	            pstream.close();
			} catch (FileNotFoundException e) {
				throw new IllegalArgumentException(e.getMessage());
			}
        }
        
        alignment = alignmentInput.get();
        alignment.sequenceInput.get().clear();
        alignment.sequenceInput.get().addAll(sequenceInput.get());
        
        super.initAndValidate();
    }


    
    /**
     * Convert integer representation of sequence into a Sequence
     *
     * @param seq  integer representation of the sequence
     * @param node used to determine taxon for sequence
     * @return Sequence
     */
    Sequence intArray2Sequence(int[] seq, Node node) {
        DataType dataType = m_dataType;
        
        
        String seqString = dataType.encodingToString(seq);
        
        // Find taxon with same name if tree is labelled
        int taxonNum = node.getNr();
        if (node.getID() != null && !node.getID().isEmpty()) {
	        for (int i = 0; i < alignmentInput.get().getTaxaNames().size(); i ++) {
	        	if (alignmentInput.get().getTaxaNames().get(i).equals(node.getID())) {
	        		taxonNum=i;
	        		break;
	        	}
	        }
        }
        
        String taxon = alignmentInput.get().getTaxaNames().get(taxonNum);
        
        
        return new Sequence(taxon, seqString);
    }

    /**
     * perform the actual sequence generation
     *
     * @return alignment containing randomly generated sequences for the nodes in the
     *         leaves of the tree
     */
    public void simulate() {
        Node root = m_tree.getRoot();
        
        int[] siteModelCategory = new int[m_sequenceLength];
        int[] category = new int[m_sequenceLength];
        int[] seq = new int[m_sequenceLength];
        for (int siteNr = 0; siteNr < m_sequenceLength; siteNr ++) {
        	
        	if (this.oneSiteModelPerSite) {
        		siteModelCategory[siteNr] = Randomizer.nextInt(this.m_pSiteModelInput.get().size());
        		//siteModelCategory[siteNr] = 1; // test
        	}else {
        		siteModelCategory[siteNr] = 0;
        	}
        	
        	if (siteModelNumberInput.get() != null) {
        		siteModelNumberInput.get().setValue(siteNr, siteModelCategory[siteNr]);
        	}
        	
			SiteModel.Base rootSiteModel = this.getSiteModel(root.getNr(), siteModelCategory[siteNr]);
			double[] categoryProbs = rootSiteModel.getCategoryProportions(root);
			category[siteNr] = Randomizer.randomChoicePDF(categoryProbs);
			double[] frequencies = rootSiteModel.getSubstitutionModel().getFrequencies();
			seq[siteNr] = Randomizer.randomChoicePDF(frequencies);
        }
        
        this.sequencesAll[root.getNr()] = seq;
        
        this.mutationsAlongEachBranch = new ArrayList<>();
        for (int i = 0; i < m_tree.getNodeCount(); i++) {
        	this.mutationsAlongEachBranch.add(null);
        }
        traverse(root, seq, category, siteModelCategory);

    } 

    /**
     * recursively walk through the tree top down, and add sequence to alignment whenever
     * a leave node is reached.
     *
     * @param node           reference to the current node, for which we visit all children
     * @param parentSequence randomly generated sequence of the parent node
     * @param category       array of categories for each of the sites
     * @param alignment
     */
    void traverse(Node node, int[] parentSequence, int[] category, int[] siteModelCategory) {
    	
    	
    	
    	
        for (int childIndex = 0; childIndex < 2; childIndex++) {
        	
        	
        	List<Mutation> mutationsBranch = new ArrayList<>();
            Node child = (childIndex == 0 ? node.getLeft() : node.getRight());
            int[] seq = new int[m_sequenceLength];
            for (int i = 0; i < m_sequenceLength; i++) {
            	SiteModel.Base nodeSiteModel = this.getSiteModel(node.getNr(), siteModelCategory[i]);
                double clockRate = (m_branchRateModel == null ? 1.0 : m_branchRateModel.getRateForBranch(child));
                double siteRate = nodeSiteModel.getRateForCategory(category[i], child);
                seq[i] = MutationUtils.simulateMutationsDownBranch(parentSequence[i], child, clockRate*siteRate, nodeSiteModel.getSubstitutionModel(), mutationsBranch, i, this.getDataTypeOfMapper().getStateCount());
            }
            
            this.mutationsAlongEachBranch.set(child.getNr(), mutationsBranch);
            this.sequencesAll[child.getNr()] = seq;

            if (child.isLeaf()) {
                sequenceInput.setValue(intArray2Sequence(seq, child), this);
            } else {
                traverse(child, seq, category, siteModelCategory);
            }
        }
    }

    
    @Override
    public List<Mutation> getMutationsOnBranch(int nodeNr){
		return this.mutationsAlongEachBranch.get(nodeNr);
	}

    
	public int[] getSequenceForNode(Node node) {
		return sequencesAll[node.getNr()];
	}

	@Override
	public Tree getTree() {
		return m_tree;
	}


	@Override
	public void sampleMutations(long sampleNr) {
		if (lastSampleNr == sampleNr) return;
		simulate();
		lastSampleNr = sampleNr;
	}


	@Override
	public List<Mutation> getMutationsOnBranchAtSite(int nodeNr, int siteNr){
		List<Mutation> muts = this.mutationsAlongEachBranch.get(nodeNr);
		List<Mutation> siteMutations = new ArrayList<>();
		for (int i = 0; i < muts.size(); i ++) {
			Mutation m = muts.get(i);
			if (m.getSiteNr() == siteNr) {
				siteMutations.add(m);
			}
		}
		
		return siteMutations;
	}



	@Override
	public int[] getStatesForNode(Tree tree, Node node) {
		return getSequenceForNode(node);
	}


	@Override
	public StochasticMapper getUnconditionalData() {
		return this;
	}

	@Override
	public Codon getDataType() {
		if (m_dataType == null) {
			initDataType();
		}
		return super.getDataType();
	}

	@Override
	public DataType getDataTypeOfMapper() {
		if (m_dataType == null) {
			initDataType();
		}
		return super.getDataType();
	}
	
//	@Override
//	public int getPatternCount() {
//		return m_sequenceLength;
//	}
	
	@Override
	public int getPatternCount() {
        return super.getPatternCount();
    }
	
	@Override
	public int getSequenceLength() {
        return this.m_sequenceLength;
    }
	

	
    private SiteModel.Base getSiteModel(int nodeNr, int siteModelNr){
    	
    	if (m_pSiteModelInput.get().size() == 1) {
    		return m_pSiteModelInput.get().get(0);
    	}
    	
    	Tree tree = m_treeInput.get();
    	
    	if (oneSiteModelPerSite) {
    		return m_pSiteModelInput.get().get(siteModelNr);
    	}
    	
    	
    	
    	List<Integer> changes = new ArrayList<>();
    	for (int i = 0; i < siteModelChangesInput.get().getDimension(); i++) {
    		changes.add(siteModelChangesInput.get().getNativeValue(i));
    	}
    	
    	// Find the most recent site model above this one
    	Node node = tree.getNode(nodeNr);
    	
    	while (!changes.contains(node.getNr())) {
    		node = node.getParent();
    		if (node == null) return null;
    	}
    	
		int index = changes.indexOf(node.getNr());
		
		//System.out.println(nodeNr + " has site model " + index + " and " + m_pSiteModelInput.get().get(index).getSubstitutionModel().getFrequencies()[0]);
		return m_pSiteModelInput.get().get(index);
    	
    }

    
    


//    @Override
//    public int getSiteCount() {
//    	return m_sequenceLength;
//    }
//
//	@Override
//	public DataType getDataTypeOfSimulator() {
//		return super.getDataType();
//	}
//    	
//
//	@Override
//    public List<Mutation> getMutationsOnBranch(int nodeNr){
//		return this.mutationsTree.get(nodeNr);
//	}
//
//    @Override
//	public int[] getSequenceForNode(Node node) {
//		return sequencesAll[node.getNr()];
//	}
//
//	@Override
//	public Tree getTree() {
//		return m_tree;
//	}
//
//	@Override
//	public Alignment getData() {
//		return this;
//	}
//	
//	
//	@Override
//	public Codon getDataTypeOfMapper() {
//		if (m_dataType == null) {
//			initDataType();
//		}
//		return super.getDataType();
//	}
	    
    
} // class SequenceAlignment

