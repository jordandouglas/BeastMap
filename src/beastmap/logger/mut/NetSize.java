package beastmap.logger.mut;

import java.util.List;

import beast.base.core.Description;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.tree.Node;
import beastmap.logger.BranchSubstLogger;
import beastmap.util.Mutation;


@Description("Records the change in length from parent to child. The length is the number of non-gap sites.")
public class NetSize extends BranchSubstLogger {

	
	int gapChar;
    public void initAndValidate() {
		gapChar = samplerInput.get().getDataType().stringToEncoding(""+DataType.GAP_CHAR).get(0);
		super.initAndValidate();
	}
	
	@Override
	public String getName() {
		return "netsize" + "." + this.getID();
	}

	@Override
	public double getFilteredMutationSummary(List<Mutation> mutations, Node node) {
		
		
		// Child size
		int ungappedSizeChild = 0;
		int[] sequenceChild = samplerInput.get().getStatesForNode(getTree(), node);
		for (int i = 0; i < sequenceChild.length; i ++) {
			int c = sequenceChild[i];
			if (c != gapChar) ungappedSizeChild++;
		}
		
		// Parent size
		int ungappedSizeParent = 0;
		if (node.getParent() != null) {
			int[] sequenceParent = samplerInput.get().getStatesForNode(getTree(), node.getParent());
			for (int i = 0; i < sequenceParent.length; i ++) {
				int c = sequenceParent[i];
				if (c != gapChar) ungappedSizeParent++;
			}
			
		}
				
		
		
		return ungappedSizeParent - ungappedSizeChild;
		
	}

	@Override
	protected boolean canHandleDataType(DataType dataType) {
		return true;
	}

}
