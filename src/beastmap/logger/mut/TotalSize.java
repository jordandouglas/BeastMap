package beastmap.logger.mut;

import java.util.List;

import beast.base.core.Description;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.tree.Node;
import beastmap.logger.BranchSubstLogger;
import beastmap.util.Mutation;


@Description("Records the length of this node. The length is the number of non-gap sites.")
public class TotalSize extends BranchSubstLogger {

	int gapChar;
	
	@Override
    public void initAndValidate() {
		gapChar = samplerInput.get().getDataType().stringToEncoding(""+DataType.GAP_CHAR).get(0);
		super.initAndValidate();
	}
	
	@Override
	public String getName() {
		return "size" + "." + this.getID();
	}

	@Override
	public double getFilteredMutationSummary(List<Mutation> mutations, Node node) {
		
		int ungappedSize = 0;
		int[] sequence = samplerInput.get().getStatesForNode(getTree(), node);
		for (int i = 0; i < sequence.length; i ++) {
			int c = sequence[i];
			if (c != gapChar) ungappedSize++;
		}
		return ungappedSize;
		
	}

	@Override
	protected boolean canHandleDataType(DataType dataType) {
		return true;
	}

}
