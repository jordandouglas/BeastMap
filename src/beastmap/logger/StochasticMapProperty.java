package beastmap.logger;

import beast.base.evolution.tree.Node;

public interface StochasticMapProperty {
	
	
	public String getName();
	public void sampleMutations(long sampleNr);
	public Object getPropertyOfNode(Node node);
	

}
