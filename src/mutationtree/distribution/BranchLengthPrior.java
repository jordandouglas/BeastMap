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
	
	
	@Override
	public double calculateTreeLogLikelihood(TreeInterface tree) {
		
		double rate = rateInput.get().getArrayValue();
		
		double logPTree = 0;
		for (Node node : treeInput.get().getNodesAsArray()) {
			if (node.isRoot()) continue;
			double length = node.getLength();
			if (length <= 0) {
				return Double.NEGATIVE_INFINITY;
			}
			logPTree += Math.log(rate) - rate*length; // Exponential distribution density in log sapce
			
			if (node.getHeight() < 0) {
				Log.warning("height=" + node.getHeight());
			}
			
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
