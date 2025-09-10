package beastmap.logger.mut;

import java.util.List;

import beast.base.core.Description;
import beast.base.core.Log;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.Nucleotide;
import beastmap.logger.BranchSubstLogger;
import beastmap.util.Mutation;

@Description("Count the number of transversions")
public class NucleotideTransversionCounter extends BranchSubstLogger {
	
	
	final int a = 0;
	final int c = 1;
	final int g = 2;
	final int t = 3;

	
	@Override
    public void initAndValidate() {
		super.initAndValidate();
	}
	
	@Override
	public double getFilteredMutationSummary(List<Mutation> mutations) {
		
		
		int count = 0;
		for (Mutation mut : mutations) {
			if (mut.getFrom() == mut.getTo()) {
				Log.warning("Dev error: 24294 from and to are the same " + mut.getFrom());
			}
			if (mut.getFrom() == a && mut.getTo() != g) count ++;
			if (mut.getFrom() == g && mut.getTo() != a) count ++;
			if (mut.getFrom() == c && mut.getTo() != t) count ++;
			if (mut.getFrom() == t && mut.getTo() != c) count ++;
		}
		return count;
		
	}

	@Override
	public String getName() {
		return "transversion" + "." + this.getID();
	}

	@Override
	protected boolean canHandleDataType(DataType dataType) {
		return dataType instanceof Nucleotide;
	}

}
