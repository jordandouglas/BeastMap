package mutationtree.evolution;

import java.util.List;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import mutationtree.util.Mutation;

public interface RecordedMutationSimulator {

	public List<Mutation> getMutationsOnBranch(int dim);

	public Tree getTree();

	public int[] getSequenceForNode(Node node);

	public DataType getDataTypeOfSimulator();

	public Alignment getData();
	
	
	

}
