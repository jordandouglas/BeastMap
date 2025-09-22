package beastmap.logger;

import java.io.PrintStream;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.parameter.RealParameter;
import beastmap.evolution.StochasticMapper;
import beastmap.util.Mutation;
import beastmap.util.MutationUtils;



@Description("Logs the tree such that each branch has a difference reconsturcted sequence. Will no longer be a strictly binary tree.")
public class TypedTreeLogger extends BEASTObject implements Loggable {

	final public Input<Tree> treeInput = new Input<>("tree", "tree to be logged", Validate.REQUIRED);
	final public Input<StochasticMapper> samplerInput = new Input<>("sampler", "mutation sampler to log", Validate.REQUIRED);
	final public Input<Integer> decimalPlacesInput = new Input<>("dp", "the number of decimal places to use writing branch lengths, rates and real-valued metadata, use -1 for full precision (default = full precision)", -1);
	final public Input<String> typeLabelInput = new Input<>("typeLabel", "the name of the branch metadata variable that stores the state", "state");
	
	final public Input<String> filterInput = new Input<>("filter", "specifies which of the sites in the input alignment we will restrict to" +
            "First site is 1." +
            "Filter specs are comma separated, either a singleton, a range [from]-[to] or iteration [from]:[to]:[step]; " +
            "1-100 defines a range, " +
            "1-100\3 or 1:100:3 defines every third in range 1-100, " +
            "1::3,2::3 removes every third site. " +
            "Default for range [1]-[last site], default for iterator [1]:[last site]:[1]", Validate.OPTIONAL);


	protected List<Integer> filter;
	
	
	String typeLabel;
	StochasticMapper sampler;
	
	private DecimalFormat df;

    @Override
    public void initAndValidate() {
    	
    	this.typeLabel = typeLabelInput.get();
    	this.sampler = samplerInput.get();
    	this.filter = null;
    	if (filterInput.get() != null && !filterInput.get().isEmpty()) {
    		this.filter = MutationUtils.parseFilterSpec(sampler.getPatternCount(), filterInput.get());
    	}
    	
    	
        int dp = decimalPlacesInput.get();
        if (dp < 0) {
            df = null;
        } else {
            // just new DecimalFormat("#.######") (with dp time '#' after the decimal)
            df = new DecimalFormat("#."+new String(new char[dp]).replace('\0', '#'));
            df.setRoundingMode(RoundingMode.HALF_UP);
        }

    }


	@Override
	public void init(PrintStream out) {
		Tree tree = this.getCurrentStatedTree(0);
		tree.init(out);
	}

	@Override
	public void log(long sample, PrintStream out) {
		

		Tree tree = this.getCurrentStatedTree(sample);
		
		out.print("tree STATE_" + sample + " = ");
		out.println(this.toNewick(tree.getRoot(), new ArrayList<>(), null));
        out.print(";");
		
		
	}

	@Override
	public void close(PrintStream out) {
		Tree mtTree = this.getCurrentStatedTree(0);
		mtTree.close(out);
	}
	
	
	String toNewick(Node node, List<Function> metadataList, BranchRateModel.Base branchRateModel) {
        StringBuffer buf = new StringBuffer();
        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(toNewick(node.getLeft(), metadataList, branchRateModel));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toNewick(node.getRight(), metadataList, branchRateModel));
            }
            buf.append(")");
        } else {
            buf.append(node.getNr() + 1);
        }
		StringBuffer buf2 = new StringBuffer();
		buf2.append("[&");
		if (metadataList.size() > 0) {
			boolean needsComma = false;
			for (Function metadata : metadataList) {
				if (metadata instanceof Parameter<?>) {
					Parameter<?> p = (Parameter<?>) metadata;
					int dim = p.getMinorDimension1();
					if (p.getMinorDimension2() > node.getNr()) {
						if (needsComma) {
							buf2.append(",");
						}
						buf2.append(((BEASTObject) metadata).getID());
						buf2.append('=');
						if (dim > 1) {
							buf2.append('{');
							for (int i = 0; i < dim; i++) {
								if (metadata instanceof RealParameter) {
									RealParameter rp = (RealParameter) metadata;
									appendDouble(buf2, rp.getMatrixValue(node.getNr(), i));
								} else {
									buf2.append(p.getMatrixValue(node.getNr(), i));
								}
								if (i < dim - 1) {
									buf2.append(',');
								}
							}
							buf2.append('}');
						} else {
							if (metadata instanceof RealParameter) {
								RealParameter rp = (RealParameter) metadata;
								appendDouble(buf2, rp.getArrayValue(node.getNr()));
							} else {
								buf2.append(metadata.getArrayValue(node.getNr()));
							}
						}
						needsComma = true;
					} else {
					
					}
				} else {
					if (metadata.getDimension() > node.getNr()) {
						if (needsComma) {
							buf2.append(",");
						}
						buf2.append(((BEASTObject) metadata).getID());
						buf2.append('=');
						buf2.append(metadata.getArrayValue(node.getNr()));
						needsComma = true;
					}
				}
			}
			if (buf2.length() > 2 && branchRateModel != null) {
				buf2.append(",");
			}
		}
		
//		if (branchRateModel != null) {
//			buf2.append("rate=");
//			appendDouble(buf2, branchRateModel.getRateForBranch(node));
//		}
		
		buf2.append(typeLabel + "='" +  node.getMetaData(typeLabel) + "'");
		
		buf2.append(']');
		if (buf2.length() > 3) {
			buf.append(buf2.toString());
		}
        buf.append(":");
        appendDouble(buf, node.getLength());
        return buf.toString();
    }

	
	
	
	public Tree getCurrentStatedTree(long sampleNr) {
		
		
		sampler.sampleMutations(sampleNr);
		
		// Convert the tree into a tree where each node has a difference sequence
		Tree originalTree = treeInput.get();
		Node newRoot = originalTree.getRoot().copy(); // Clone
		List<Node> newNodes = newRoot.getAllChildNodesAndSelf();
		
		for (int i = 0; i < newNodes.size(); i ++) {
			
			Node node = newNodes.get(i);
			int nodeNr = node.getNr();
			
			if (!node.isRoot()) {
				
				
				double startAge = node.getParent().getHeight();
				int[] seq = sampler.getStatesForNode(sampler.getTree(), node.getParent());
				int[] seqCopy = new int[seq.length];
				for (int j = 0; j < seq.length; j ++) seqCopy[j] = seq[j];
				
				// Sort by time along branch (forward in time)
				List<Mutation> mutations = sampler.getMutationsOnBranch(nodeNr);
				mutations = this.filterMutations(mutations);
				Collections.sort(mutations);
				for (Mutation mutation : mutations) {
					
					
					// Insert a new branch between this node and its parent
					Node older = new Node();
					Node p = node.getParent();
					p.removeChild(node);
					p.addChild(older);
					older.addChild(node);
					
					double height = startAge - mutation.getTime();
					
					
					
					String newState = sampler.getDataTypeOfMapper().encodingToString(filterSequence(seqCopy));
					
					// Set age and state
					older.setHeight(height);
					older.setMetaData(typeLabel, newState); 
	
					if (older.getHeight() <= node.getHeight()) {
						Log.warning("negative branch length " + older.getHeight() +  "<=" + node.getHeight());
					}
					
					// Mutation happens AFTER this child
					int siteNr = mutation.getSiteNr();
					seqCopy[siteNr] = mutation.getTo();
					
				}
				
			}
			
			// height unchanged
			int[] childSeq = sampler.getStatesForNode(sampler.getTree(), node);
			//int[] childSeq = sampler.getStatesForNode(sampler.getTree(), node);
			String state = sampler.getDataTypeOfMapper().encodingToString(filterSequence(childSeq)); 
			node.setMetaData(typeLabel, state);
			
//			if (node.isLeaf()) {
//				Log.warning("Setting " + node.getID() + " to " + state);
//			}
			
			
		}
		

		
		//Log.warning("new tree " + newTree.getNodeCount());
		for (Node node : newRoot.getAllChildNodesAndSelf()) {
			String type = (String) node.getMetaData(typeLabel);
			node.metaDataString = "&" + typeLabel + "='" + type + "'";
			
			
			
		}
		
		List<String> labels = new ArrayList<>();
		for (int i = 0; i < newRoot.getLeafNodeCount(); i ++) {
			labels.add(null);
		}
		for (Node leaf : newRoot.getAllLeafNodes()) {
			labels.set(leaf.getNr(), leaf.getID());
		}
		
		String newick = newRoot.toString(labels);
		//int[] dummy = new int[1];
        //String newick = newRoot.toSortedNewick(dummy, true);
		Tree t2 = new Tree(newick);
		
		
//		for (Node leaf : t2.getRoot().getAllLeafNodes()) {
//			String id = labels.get(leaf.getNr());
//			leaf.setID(id);
//		}
//		
		
		//Log.warning(newick);
		
		// This seems to duplicate each leaf name from 'x' into 'xx'. We will correct this now
		for (Node leaf : t2.getRoot().getAllLeafNodes()) {
			String leafName = leaf.getID();
			int len = leafName.length();
			String first = leafName.substring(0, len/2);
			String last = leafName.substring(len/2, len);
			if (first.equals(last)) {
				//Log.warning("Taking first");
				leaf.setID(first);
			}
			
		}

		return t2;
		

	
	
	}
	
	private int[] filterSequence(int[] fullSeq) {
		
		if (this.filter == null) return fullSeq;
		
		int[] filtered = new int[this.filter.size()];
		int k = 0;
		for (int j = 0; j < fullSeq.length; j ++) {
			if (this.filter.contains(j)) {
				filtered[k] = fullSeq[j];
				k++;
			}
		}
		
		return filtered;
	}
	
	
    /**
     * Appends a double to the given StringBuffer, formatting it using
     * the private DecimalFormat instance, if the input 'dp' has been
     * given a non-negative integer, otherwise just uses default
     * formatting.
     * @param buf
     * @param d
     */
    private void appendDouble(StringBuffer buf, double d) {
        if (df == null) {
            buf.append(d);
        } else {
            buf.append(df.format(d));
        }
    }
	

    private List<Mutation> filterMutations(List<Mutation> mutations){
    	
    	
    	if (this.filter == null) {
    		return mutations;
    	}
    	
    		
		// Filter the mutations and then give the remainder to the child class
		List<Mutation> filteredMutations = new ArrayList<Mutation>();
    	for (Mutation mutation : mutations) {
    		if (this.filter.contains(mutation.getSiteNr())){
    			filteredMutations.add(mutation);
    		}
    	}
    	
    	return filteredMutations;
	
    	
    }
    

	
}
