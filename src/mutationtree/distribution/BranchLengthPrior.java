package mutationtree.distribution;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;


@Description("A tree prior for a non-time tree")
public class BranchLengthPrior extends SpeciesTreeDistribution {

	final public Input<RealParameter> rateInput = new Input<>("rate", "exponential distribution prior rate of branch lengths", Input.Validate.REQUIRED);
	final public Input<RealParameter> midpointInput = new Input<>("midpoint", "the difference between the root's two furtherest tips is an exponential with this mean", Input.Validate.OPTIONAL);
	
	
	@Override
	public double calculateTreeLogLikelihood(TreeInterface tree) {
		
		double rate = rateInput.get().getArrayValue();
		
		double logPTree = 0;
		for (Node node : treeInput.get().getNodesAsArray()) {
			if (node.isRoot()) continue;
			double length = node.getLength();
			
//			// So that this process is generative, the two branches coming out of the root will count as a single branch
			// Edit: since we are using an exponential prior, this step is redundant and gives the same answer
//			if (midpointInput.get() != null && node.getParent().isRoot()) {
//				
//				
//				// Skip one branch
//				if (node == node.getParent().getRight()){
//					continue;
//				}else {
//					Node left = node;
//					Node right = node.getParent().getRight();
//					length = left.getHeight() + right.getHeight();
//				}
//				
//			}
			
			if (length <= 0) {
				return Double.NEGATIVE_INFINITY;
			}
			
			logPTree += Math.log(rate) - rate*length; // Exponential distribution density in log space
			
			if (node.getHeight() < 0) {
				Log.warning("height=" + node.getHeight());
			}
			
		}
		
		if (midpointInput.get() != null) {
			Node root = treeInput.get().getRoot();
			
			// Find the furtherest tip on either side
			Node left = root.getLeft();
			Node right = root.getRight();
			
			double leftHeight = Double.POSITIVE_INFINITY;
			double rightHeight = Double.POSITIVE_INFINITY;
			
			for (Node node : left.getAllLeafNodes()) {
				leftHeight = Math.min(leftHeight, node.getHeight());
			}
			for (Node node : right.getAllLeafNodes()) {
				rightHeight = Math.min(rightHeight, node.getHeight());
			}
			
			
			// Difference between heights
			double rootToLeft = root.getHeight() - leftHeight;
			double rootToRight = root.getHeight() - rightHeight;
			double diff = Math.abs(rootToLeft - rootToRight);
			double midpointRate = 1 / midpointInput.get().getArrayValue();
			logPTree += Math.log(midpointRate) - midpointRate*diff; // Exponential distribution density in log space
			
			
		}
		
		return logPTree;
	}
	
    @Override
    public List<String> getConditions() {
        List<String> conditions = new ArrayList<>();
        conditions.add(rateInput.get().getID());
        return conditions;
    }

    @Override
    public List<String> getArguments() {
        List<String> arguments = new ArrayList<>();
        arguments.add(treeInput.get().getID());
        return arguments;
    }
    
    @Override
    protected boolean requiresRecalculation() {
    	if (InputUtil.isDirty(rateInput)) return true;
        return treeInput.get().somethingIsDirty();
    }

}
