package mutationtree.logger.mut;

import java.util.Collections;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.datatype.DataType;
import mutationtree.logger.BranchSubstLogger;
import mutationtree.util.Mutation;

@Description("Count the number of mutations from X to Y where X and Y are sets")
public class FromToSubstSum extends BranchSubstLogger {
	
	
	final public Input<String> fromInput = new Input<>("from", "'from' character(s), e.g., 'AG' for nucleotides. Order has no meaning and duplicates are ignored.", Input.Validate.REQUIRED);
	final public Input<String> toInput = new Input<>("to", "'to' character(s), e.g., 'T' for nucleotides. Order has no meaning and duplicates are ignored.", Input.Validate.REQUIRED);

	
	List<Integer> from;
	List<Integer> to;
	
	@Override
    public void initAndValidate() {
		
		// Check that we can handle this data type
		super.initAndValidate();
		
		DataType dt = samplerInput.get().dataInput.get().getDataType();
		this.from = dt.stringToEncoding(fromInput.get());
		this.to = dt.stringToEncoding(toInput.get());
		
		
	}
	
	@Override
	public double getMutationSummary(List<Mutation> mutations) {
		
		int count = 0;
		for (Mutation mut : mutations) {
			if (this.from.contains(mut.getFrom()) && this.to.contains(mut.getTo())) {
				count ++;
			}
		}
		return count;
		
	}

	@Override
	protected String getName() {
		return "from" + fromInput.get() + "to" + toInput.get();
	}

	@Override
	protected boolean canHandleDataType(DataType dataType) {
		
		try {
			dataType.stringToEncoding(fromInput.get());
			dataType.stringToEncoding(toInput.get());
		}catch (Exception e) {
			e.printStackTrace();
			return false;
		}
		
		return true;
	}

}
