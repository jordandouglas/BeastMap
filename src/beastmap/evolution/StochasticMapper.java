package beastmap.evolution;

import java.util.List;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastmap.util.Mutation;

public interface StochasticMapper {
	
	public void sampleMutations(long sampleNr);

	public List<Mutation> getMutationsOnBranch(int dim);

	public List<Mutation> getMutationsOnBranchAtSite(int nodeNr, int siteNr);
	
	public Tree getTree();

	public int[] getStatesForNode(Tree tree, Node node);

	public DataType getDataType();
	
	public int getPatternCount();
	

}
