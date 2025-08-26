package mutationtree.operator;

import beast.base.core.Description;
import beast.base.evolution.operator.TipDatesRandomWalker;
import beast.base.evolution.tree.Node;


@Description("Moves a random tip, and then adjusts the height of each node so that all heights are non-negative")
public class UnconstrainedTipRandomWalker extends TipDatesRandomWalker {

	
	@Override
    public double proposal() {
		double logHR = super.proposal();
		
		
		// Get youngest tip
		double minHeight = Double.POSITIVE_INFINITY;
		for (Node node : treeInput.get().getNodesAsArray()) {
			minHeight = Math.min(minHeight, node.getHeight());
		}
		
		for (Node node : treeInput.get().getNodesAsArray()) {
			node.setHeight(node.getHeight() - minHeight);
		}
		
		
		return logHR;
		 
	 }
}
