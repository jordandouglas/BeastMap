package mutationtree.evolution;

import beast.base.evolution.tree.Tree;

public class TimelessTree extends Tree {
	
	
	@Override
    public int scale(final double scale) {
		
		// As we are parameterising by branch length, not node height, the number of terms scaled is actually the number of branches not the number of internal nodes
        root.scale(scale);
        return this.getNodeCount() - 1;
    }
	

}
