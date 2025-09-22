package beastmap.logger.mut;

import java.util.List;

import beast.base.core.Description;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.tree.Node;
import beastmap.logger.BranchSubstLogger;
import beastmap.util.Mutation;


@Description("Number of mutations")
public class SubstitutionSum extends BranchSubstLogger {

	
	@Override
    public void initAndValidate() {
		super.initAndValidate();
	}
	
	@Override
	public double getFilteredMutationSummary(List<Mutation> mutations, List<Mutation> mutationsUnconditional, Node node) {
		return mutations.size();
	}

	@Override
	public String getName() {
		return "substitutions" + "." + this.getID();
	}

	@Override
	protected boolean canHandleDataType(DataType dataType) {
		return true;
	}

}
