package beastmap.evolution;


import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.core.ProgramStatus;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.UserDataType;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.likelihood.LikelihoodCore;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.util.Randomizer;
import beastfx.app.beauti.Beauti;



/**
 * Adapted from AncestralStateTreeLikelihood in beast classic, which came from BEAST X
 */
@Description("Ancestral State Tree Likelihood")
public class AncestralSequenceTreeLikelihood extends TreeLikelihood  {
    public static final String STATES_KEY = "states";

    public Input<String> tagInput = new Input<String>("tag","label used to report trait", "seq");
    public Input<Boolean> useMAPInput = new Input<Boolean>("useMAP","whether to use maximum aposteriori assignments or sample", false);
    public Input<Boolean> returnMLInput = new Input<Boolean>("returnML", "report integrate likelihood of tip data", true);
    public Input<Boolean> sampleTipsInput = new Input<Boolean>("sampleTips", "if tips have missing data/ambigous values sample them for logging (default true)", false);
    
    
    public Input<GenericTreeLikelihood> likelihoodInput = new Input<GenericTreeLikelihood>("likelihood", "tree likelihood that we will use to get the site, clock, and tree");
    
    
    public AncestralSequenceTreeLikelihood() {
		treeInput.setRule(Validate.OPTIONAL);
		siteModelInput.setRule(Validate.OPTIONAL);
		branchRateModelInput.setRule(Validate.OPTIONAL);
	}
    

	int[][] storedTipStates;

	/** parameters for each of the leafs **/
	IntegerParameter[] parameters;

	/** and node number associated with parameter **/
	int[] leafNr;

	int traitDimension;

    /**
     * Constructor.
     * Now also takes a DataType so that ancestral states are printed using data codes
     *
     * @param patternList     -
     * @param treeModel       -
     * @param siteModel       -
     * @param branchRateModel -
     * @param useAmbiguities  -
     * @param storePartials   -
     * @param dataType        - need to provide the data-type, so that corrent data characters can be returned
     * @param tag             - string label for reconstruction characters in tree log
     * @param forceRescaling  -
     * @param useMAP          - perform maximum aposteriori reconstruction
     * @param returnML        - report integrate likelihood of tip data
     */
    int patternCount;
    int stateCount;

    int[][] tipStates; // used to store tip states when using beagle
    
    @Override
    public void initAndValidate() {
    	
    	// Beauti?
    	boolean inBEAUti = ProgramStatus.name.equals("BEAUti");
    	if (inBEAUti || dataInput.get().getSiteCount() == 0 || dataInput.get().getTaxonCount() == 0 || dataInput.get().getCounts().isEmpty()) {
    		return;
    	}
    	
//    	Log.warning("sites " + dataInput.get().getSiteCount());
//    	Log.warning("taxa " + dataInput.get().getTaxonCount());
    	
    	
    	
    	
    	if (likelihoodInput.get() != null) {
    		treeInput.setValue(likelihoodInput.get().treeInput.get(), this);
    		siteModelInput.setValue(likelihoodInput.get().siteModelInput.get(), this);
    		branchRateModelInput.setValue(likelihoodInput.get().branchRateModelInput.get(), this);
    	}
    	
    	// Do not use beagle as the likelihood core is not implemented in beagle
    	implementationInput.setValue(AncestralSequenceTreeLikelihood.class.getName(), this);
    	
    	super.initAndValidate();

    	
        this.tag = tagInput.get();
        TreeInterface treeModel = treeInput.get();
        patternCount = dataInput.get().getPatternCount();
        dataType = dataInput.get().getDataType();
        stateCount = dataType.getStateCount();

        reconstructedStates = new int[treeModel.getNodeCount()][patternCount];
        storedReconstructedStates = new int[treeModel.getNodeCount()][patternCount];

        this.useMAP = useMAPInput.get();
        this.returnMarginalLogLikelihood = returnMLInput.get();
      
        
        if (beagle != null) {
            if (!(siteModelInput.get() instanceof SiteModel.Base)) {
            	throw new IllegalArgumentException ("siteModel input should be of type SiteModel.Base");
            }
            m_siteModel = (SiteModel.Base) siteModelInput.get();
        	substitutionModel = (SubstitutionModel.Base) m_siteModel.substModelInput.get();
            int nStateCount = dataInput.get().getMaxStateCount();
            probabilities = new double[(nStateCount + 1) * (nStateCount + 1)];
        }

        int tipCount = treeModel.getLeafNodeCount();
        tipStates = new int[tipCount][];

        Alignment data = dataInput.get();
        for (Node node : treeInput.get().getExternalNodes()) {
            String taxon = node.getID();
            int taxonIndex = data.getTaxonIndex(taxon);
            if (taxonIndex == -1) {
            	if (taxon.startsWith("'") || taxon.startsWith("\"")) {
                    taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
                }
                if (taxonIndex == -1) {
                	throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
                }
            }
            tipStates[node.getNr()] = new int[patternCount];
            if (!m_useAmbiguities.get()) {
            	likelihoodCore.getNodeStates(node.getNr(), tipStates[node.getNr()]);
            } else {
            	int [] states = tipStates[node.getNr()];
	            for (int i = 0; i < patternCount; i++) {
	                int code = data.getPattern(taxonIndex, i);
	                int[] statesForCode = data.getDataType().getStatesForCode(code);
	                if (statesForCode.length==1)
	                    states[i] = statesForCode[0];
	                else
	                    states[i] = code; // Causes ambiguous states to be ignored.
	            }

            }
    	}
        
        
        traitDimension = tipStates[0].length;
        


		storedTipStates = new int[tipStates.length][traitDimension];
		for (int i = 0; i < tipStates.length; i++) {
			System.arraycopy(tipStates[i], 0, storedTipStates[i], 0, traitDimension);
		}


    }
    
    
    @Override
    protected LikelihoodCore createLikelihoodCore(int stateCount) {
    	return new RateCategorySampledLikelihoodCore(stateCount);
    }


    @Override
    public void store() {
        super.store();

        for (int i = 0; i < reconstructedStates.length; i++) {
            System.arraycopy(reconstructedStates[i], 0, storedReconstructedStates[i], 0, reconstructedStates[i].length);
        }

        storedAreStatesRedrawn = areStatesRedrawn;
        storedJointLogLikelihood = jointLogLikelihood;
        
        
        // deal with ambiguous tips
        if (leafNr != null) {
			for (int i = 0; i < leafNr.length; i++) {
				int k = leafNr[i];
				System.arraycopy(tipStates[k], 0, storedTipStates[k], 0, traitDimension);
			}
        }
    }

    @Override
    public void restore() {

        super.restore();

        int[][] temp = reconstructedStates;
        reconstructedStates = storedReconstructedStates;
        storedReconstructedStates = temp;

        areStatesRedrawn = storedAreStatesRedrawn;
        jointLogLikelihood = storedJointLogLikelihood;
        
        // deal with ambiguous tips
        if (leafNr != null) {
			for (int i = 0; i < leafNr.length; i++) {
				int k = leafNr[i];
				int[] tmp = tipStates[k];
				tipStates[k] = storedTipStates[k];
				storedTipStates[k] = tmp;
				// Does not handle ambiguities or missing taxa
				likelihoodCore.setNodeStates(k, tipStates[k]);
			}
        }
    }
    
    @Override
    protected boolean requiresRecalculation() {
    	likelihoodKnown = false;

    	boolean isDirty = super.requiresRecalculation();
    	if (!m_useAmbiguities.get()) {
    		return isDirty;
    	}
    	
    	
    	int hasDirt = Tree.IS_CLEAN;
		
		// check whether any of the leaf trait parameters changed
		for (int i = 0; i < leafNr.length; i++) {
			if (parameters[i].somethingIsDirty()) {
				int k = leafNr[i];
				for (int j = 0; j < traitDimension; j++) {
					tipStates[k][j] = parameters[i].getValue(j);
				}
				likelihoodCore.setNodeStates(k, tipStates[k]);
				isDirty = true;
				// mark leaf's parent node as dirty
				Node leaf = treeInput.get().getNode(k);
				// leaf.makeDirty(Tree.IS_DIRTY);
				leaf.getParent().makeDirty(Tree.IS_DIRTY);
	            hasDirt = Tree.IS_DIRTY;
			}
		}
		isDirty |= super.requiresRecalculation();
		this.hasDirt |= hasDirt;

		return isDirty;
    	
    	
    }


    public DataType getDataTypeOfMapper() {
        return dataType;
    }

    public int[] getStatesForNode(TreeInterface tree, Node node) {
        if (tree != treeInput.get()) {
            throw new RuntimeException("Can only reconstruct states on treeModel given to constructor");
        }

        if (!likelihoodKnown) {
        	try {
        		 calculateLogP();
        	} catch (Exception e) {
				throw new RuntimeException(e.getMessage());
			}
        }

        if (!areStatesRedrawn) {
            redrawAncestralStates();
        }
        return reconstructedStates[node.getNr()];
    }


    public void redrawAncestralStates() {
    	
        jointLogLikelihood = 0;
        TreeInterface tree = treeInput.get();
        traverseSample(tree, tree.getRoot(), null);
        
        areStatesRedrawn = true;
        

    }


    
    @Override
    public double calculateLogP() {
        areStatesRedrawn = false;
        
        // Set the root frequencies prior to sampling gamma categories
        double[] rootFrequencies = substitutionModel.getFrequencies();
        if (rootFrequenciesInput.get() != null) {
            rootFrequencies = rootFrequenciesInput.get().getFreqs();
        }
        ((RateCategorySampledLikelihoodCore) this.likelihoodCore).setRootFrequencies(rootFrequencies);
        
        
        
        double marginalLogLikelihood = super.calculateLogP();
        likelihoodKnown = true;

        if (returnMarginalLogLikelihood) {
            return logP;
        }
        
        // redraw states and return joint density of drawn states
        redrawAncestralStates();
        logP = jointLogLikelihood;
        return logP;
    }


    
    // More stable
    private int drawChoiceLog(double[] logProbs) {
    	
    	// Max prob
    	double max = logProbs[0];
        for (int i = 1; i < logProbs.length; i++) {
            if (logProbs[i] > max) {
                max = logProbs[i];
            }
        }
        
        
        // Return a u.a.r character
        if (Double.isInfinite(max) || Double.isNaN(max)) {
        	int rand = Randomizer.nextInt(logProbs.length);
    		return rand;
        }
		
    	
        // Convert to relative probabilities in a stable fashion
    	double[] probs = new double[logProbs.length];
    	for (int i = 0; i < logProbs.length; i++) {
    		probs[i] = Math.exp(logProbs[i] - max);
    		if (Double.isNaN(probs[i])) probs[i] = 0;
    		//Log.warning(logProbs[i] + " max = " + max + " " + probs[i]);
    	}
    	
    	
    	
    	return drawChoice(probs);
        	
    }

    private int drawChoice(double[] probs) {
    	
        if (useMAP) {
            double max = probs[0];
            int choice = 0;
            for (int i = 1; i < probs.length; i++) {
                if ((probs[i] - max)/(probs[i] + max) > 1e-10) {
                    max = probs[i];
                    choice = i;
                }
            }
            return choice;
            
        } else {
        	
        	if (Randomizer.getTotal(probs) <= 1e-100) {
        		Log.warning("Warning: ancestral state reconstruction probabilities sum to 0, suggesting numerical issues. Selecting a state uniformly at random.");
        		
        		// Return a u.a.r character
        		int rand = Randomizer.nextInt(probs.length);
        		return rand;
        		
        	}
        	
        	try {
        		return Randomizer.randomChoicePDF(probs);
        	}catch(Exception e) {
        		
        		// Use MAP state
        		Log.warning("Warning: ancestral state encountered error, suggesting numerical issues. Selecting a state uniformly at random.");
        		double max = probs[0];
                int choice = 0;
                for (int i = 1; i < probs.length; i++) {
                    if ((probs[i] - max)/(probs[i] + max) > 1e-10) {
                        max = probs[i];
                        choice = i;
                    }
                }
                return choice;
        		
        	}
        }
    }

    public void getStates(int tipNum, int[] states)  {
        // Saved locally to reduce BEAGLE library access
        System.arraycopy(tipStates[tipNum], 0, states, 0, states.length);
    }


    
    /**
     * Traverse (pre-order) the tree sampling the internal node states.
     *
     * @param tree        - TreeModel on which to perform sampling
     * @param node        - current node
     * @param parentState - character state of the parent node to 'node'
     */
    public void traverseSample(TreeInterface tree, Node node, int[] parentState) {

        int nodeNum = node.getNr();

        Node parent = node.getParent();

        // This function assumes that all partial likelihoods have already been calculated
        // If the node is internal, then sample its state given the state of its parent (pre-order traversal).

        double[] conditionalProbabilities = new double[stateCount];
        int[] state = new int[patternCount];

        
        
        
        if (!node.isLeaf()) {

            if (parent == null) {

                double[] rootPartials = m_fRootPartials;

                double[] rootFrequencies = substitutionModel.getFrequencies();
                if (rootFrequenciesInput.get() != null) {
                    rootFrequencies = rootFrequenciesInput.get().getFreqs();
                }

                // This is the root node
                for (int j = 0; j < patternCount; j++) {
            		System.arraycopy(rootPartials, j * stateCount, conditionalProbabilities, 0, stateCount);

                    for (int i = 0; i < stateCount; i++) {
                        //conditionalProbabilities[i] *= rootFrequencies[i];
                        conditionalProbabilities[i] = Math.log(conditionalProbabilities[i]) + Math.log(rootFrequencies[i]);
                    }
                    try {
                        //state[j] = drawChoice(conditionalProbabilities);
                        state[j] = drawChoiceLog(conditionalProbabilities);
                        
                    } catch (Error e) {
                        System.err.println(e.toString());
                        state[j] = 0;
                    }
                    reconstructedStates[nodeNum][j] = state[j];

                    //System.out.println("Pr(j) = " + rootFrequencies[state[j]]);
                    jointLogLikelihood += Math.log(rootFrequencies[state[j]]);
                }

            } else {

                // This is an internal node, but not the root
                double[] partialLikelihood = new double[stateCount * patternCount];



                likelihoodCore.getNodePartials(node.getNr(), partialLikelihood);


                for (int j = 0; j < patternCount; j++) {

                    int parentIndex = parentState[j] * stateCount;
                    int childIndex = j * stateCount;

                    for (int i = 0; i < stateCount; i++) {
                        //conditionalProbabilities[i] = partialLikelihood[childIndex + i] * probabilities[parentIndex + i];
                        conditionalProbabilities[i] = Math.log(partialLikelihood[childIndex + i]) + Math.log(probabilities[parentIndex + i]);
                    }

                    
                    // Sampled ancestor?
                    if (node.getLength() <= 0) {
                    	 state[j] = parentState[j];
                    }else {
                    	 state[j] = drawChoiceLog(conditionalProbabilities);
                    }
                    
                   
                    reconstructedStates[nodeNum][j] = state[j];

                    double contrib = probabilities[parentIndex + state[j]];
                    //System.out.println("Pr(" + parentState[j] + ", " + state[j] +  ") = " + contrib);
                    jointLogLikelihood += Math.log(contrib);
                }
            }

            // Traverse down the two child nodes
            Node child1 = node.getChild(0);
            traverseSample(tree, child1, state);

            Node child2 = node.getChild(1);
            traverseSample(tree, child2, state);
        } else {

            // This is an external leaf
        	getStates(nodeNum, reconstructedStates[nodeNum]);


        	if (sampleTipsInput.get()) {
	            // Check for ambiguity codes and sample them
	            for (int j = 0; j < patternCount; j++) {
	
	                final int thisState = reconstructedStates[nodeNum][j];
	                final int parentIndex = parentState[j] * stateCount;
	                likelihoodCore.getNodeMatrix(nodeNum, 0, probabilities);
	                if (dataType.isAmbiguousCode(thisState)) {
		                    
	                    boolean [] stateSet = dataType.getStateSet(thisState);
	                    for (int i = 0; i < stateCount; i++) {
	                        //conditionalProbabilities[i] =  stateSet[i] ? probabilities[parentIndex + i] : 0;
	                        conditionalProbabilities[i] =  stateSet[i] ? Math.log(probabilities[parentIndex + i]) : 0;
	                    }
	                    
	                    
	                    // Sampled ancestor?
	                    if (node.getLength() <= 0) {
	                    	reconstructedStates[nodeNum][j] = parentState[j];
	                    }else {
	                    	reconstructedStates[nodeNum][j] = drawChoiceLog(conditionalProbabilities);
	                    }
	                    
	                    
	                    //reconstructedStates[nodeNum][j] = drawChoice(conditionalProbabilities);
	                    
	                }
	
	                double contrib = probabilities[parentIndex + reconstructedStates[nodeNum][j]];
	                //System.out.println("Pr(" + parentState[j] + ", " + reconstructedStates[nodeNum][j] +  ") = " + contrib);
	                jointLogLikelihood += Math.log(contrib);
	            }
        	}
        	
        }
    }
    
    
    @Override
    public void log(final long sample, final PrintStream out) {
    	hasDirt = Tree.IS_FILTHY;
    	calculateLogP();
        out.print(getCurrentLogP() + "\t");
    }


    protected DataType dataType;
    private int[][] reconstructedStates;
    private int[][] storedReconstructedStates;

    private String tag;
    private boolean areStatesRedrawn = false;
    private boolean storedAreStatesRedrawn = false;

    private boolean useMAP = false;
    private boolean returnMarginalLogLikelihood = true;

    private double jointLogLikelihood;
    private double storedJointLogLikelihood;

    boolean likelihoodKnown = false;
}
