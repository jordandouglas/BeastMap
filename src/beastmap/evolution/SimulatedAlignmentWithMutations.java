package beastmap.evolution;

	
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.parser.XMLProducer;
import beast.base.util.Randomizer;
import beastmap.util.Mutation;
import beastmap.util.MutationUtils;


/**
 * Adapted from SimulatedAlignment in beast.app
 */

@Description("An alignment containing sequences randomly generated using a"
        + "given site model down a given tree.")
public class SimulatedAlignmentWithMutations extends Alignment implements StochasticMapper {
    final public Input<Alignment> m_data = new Input<>("data", "alignment data which specifies datatype and taxa of the beast.tree", Validate.REQUIRED);
    final public Input<Tree> m_treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);
    final public Input<SiteModel.Base> m_pSiteModelInput = new Input<>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);
    final public Input<BranchRateModel.Base> m_pBranchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");
    final public Input<Integer> m_sequenceLengthInput = new Input<>("sequencelength", "nr of samples to generate (default 1000).", 1000);
    final public Input<String> m_outputFileNameInput = new Input<>(
            "outputFileName",
            "If provided, simulated alignment is additionally written to this file.");
    final public Input<Long> localSeedInput = new Input<>(
            "seed",
            "Optional local random seed for simulating this alignment. If not set, global seed is used.",
            Input.Validate.OPTIONAL
    );
    
    
    
    int[][] sequencesAll;
    List<List<Mutation>> mutationsAlongEachBranch;

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
    protected SiteModel.Base m_siteModel;
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
    long lastSampleNr = -1;
    
    public SimulatedAlignmentWithMutations() {
        
        // Override the sequence input requirement.
        sequenceInput.setRule(Validate.OPTIONAL);
    }
    

    @Override
    public void initAndValidate() {
    	
    	mutationsAlongEachBranch = new ArrayList<>();
        m_tree = m_treeInput.get();
        m_siteModel = m_pSiteModelInput.get();
        m_branchRateModel = m_pBranchRateModelInput.get();
        m_sequenceLength = m_sequenceLengthInput.get();
        m_stateCount = m_data.get().getMaxStateCount();
        m_categoryCount = m_siteModel.getCategoryCount();
        m_probabilities = new double[m_categoryCount][m_stateCount * m_stateCount];
        m_outputFileName = m_outputFileNameInput.get();
        
        sequencesAll = new int[m_tree.getNodeCount()][];
        sequenceInput.get().clear();

        Long customSeed = localSeedInput.get();
        long originalSeed = Randomizer.getSeed();
        long seedToUse = customSeed != null ? customSeed : originalSeed;

        if (customSeed != null) {
            Log.info.println();
            Log.info.println("Random number seed for alignment simulation: " + customSeed);
            Log.info.println();
        }

        try {
            Randomizer.setSeed(seedToUse);
            simulate();
        } finally {
            Randomizer.setSeed(originalSeed);
        }

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
        DataType dataType = m_data.get().getDataType();
        String seqString = dataType.encodingToString(seq);
//	    	StringBuilder seq = new StringBuilder();
//	    	String map = m_data.get().getMap();
//	    	if (map != null) {
//	    		for (int i  = 0; i < m_sequenceLength; i++) {
//	    			seq.append(map.charAt(seq[i]));
//	    		}
//	    	} else {
//	    		for (int i  = 0; i < m_sequenceLength-1; i++) {
//	    			seq.append(seq[i] + ",");
//	    		}
//				seq.append(seq[m_sequenceLength-1] + "");
//	    	}
        //String taxon = m_data.get().getTaxaNames().get(node.getNr());
        
        
        
        // Find taxon with same name if tree is labelled
        int taxonNum = node.getNr();
        if (node.getID() != null && !node.getID().isEmpty()) {
	        for (int i = 0; i < m_data.get().getTaxaNames().size(); i ++) {
	        	if (m_data.get().getTaxaNames().get(i).equals(node.getID())) {
	        		taxonNum=i;
	        		break;
	        	}
	        }
        }
        
        String taxon = m_data.get().getTaxaNames().get(taxonNum);
        
        
        return new Sequence(taxon, seqString);
    } // intArray2Sequence

    /**
     * perform the actual sequence generation
     *
     * @return alignment containing randomly generated sequences for the nodes in the
     *         leaves of the tree
     */
    public void simulate() {
        Node root = m_tree.getRoot();

        //Log.warning("simulate m_sequenceLength=" + m_sequenceLength);
        
        double[] categoryProbs = m_siteModel.getCategoryProportions(root);
        int[] category = new int[m_sequenceLength];
        for (int i = 0; i < m_sequenceLength; i++) {
            category[i] = Randomizer.randomChoicePDF(categoryProbs);
        }

        
       
        
        double[] frequencies = m_siteModel.getSubstitutionModel().getFrequencies();
        int[] seq = new int[m_sequenceLength];
        for (int i = 0; i < m_sequenceLength; i++) {
            seq[i] = Randomizer.randomChoicePDF(frequencies);
        }
        this.sequencesAll[root.getNr()] = seq;
        
        
        


        //alignment.setDataType(m_siteModel.getFrequencyModel().getDataType());
        this.mutationsAlongEachBranch = new ArrayList<>();
        for (int i = 0; i < m_tree.getNodeCount(); i++) {
        	this.mutationsAlongEachBranch.add(null);
        }
        traverse(root, seq, category);

    } // simulate

    /**
     * recursively walk through the tree top down, and add sequence to alignment whenever
     * a leave node is reached.
     *
     * @param node           reference to the current node, for which we visit all children
     * @param parentSequence randomly generated sequence of the parent node
     * @param category       array of categories for each of the sites
     * @param alignment
     */
    void traverse(Node node, int[] parentSequence, int[] category) {
        for (int childIndex = 0; childIndex < 2; childIndex++) {
        	
        	
        	List<Mutation> mutationsBranch = new ArrayList<>();
            Node child = (childIndex == 0 ? node.getLeft() : node.getRight());
            int[] seq = new int[m_sequenceLength];
            for (int i = 0; i < m_sequenceLength; i++) {
                double clockRate = (m_branchRateModel == null ? 1.0 : m_branchRateModel.getRateForBranch(child));
                double siteRate = m_siteModel.getRateForCategory(category[i], child);
                seq[i] = MutationUtils.simulateMutationsDownBranch(parentSequence[i], child, clockRate*siteRate, m_siteModel.getSubstitutionModel(), mutationsBranch, i);
               // Log.warning("site " + i);
                
            }
            
            this.mutationsAlongEachBranch.set(child.getNr(), mutationsBranch);
            this.sequencesAll[child.getNr()] = seq;

            if (child.isLeaf()) {
                sequenceInput.setValue(intArray2Sequence(seq, child), this);
            } else {
                traverse(child, seq, category);
            }
        }
    } // traverse
    
    
    


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
	public DataType getDataTypeOfMapper() {
		return super.getDataType();
	}
	
	@Override
	public int getPatternCount() {
		return m_sequenceLength;
	}

//	@Override
//	public Alignment getData() {
//		return this;
//	}
//    
//	@Override
//	public DataType getDataTypeOfSimulator() {
//		return super.getDataType();
//	}
//    
    
    
    
    
    
    
    
    
    


} // class SequenceAlignment


