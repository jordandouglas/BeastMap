package mutationtree.logger.mut;

import java.util.List;

import beast.base.core.Description;
import beast.base.evolution.datatype.DataType;
import mutationtree.logger.BranchSubstLogger;
import mutationtree.util.Mutation;


@Description("Number of mutations")
public class SubstitutionSum extends BranchSubstLogger {

	
	@Override
    public void initAndValidate() {
		super.initAndValidate();
	}
	
	@Override
	public double getMutationSummary(List<Mutation> mutations) {
		return mutations.size();
	}

	@Override
	protected String getName() {
		return "substitutions";
	}

	@Override
	protected boolean canHandleDataType(DataType dataType) {
		return true;
	}

}
