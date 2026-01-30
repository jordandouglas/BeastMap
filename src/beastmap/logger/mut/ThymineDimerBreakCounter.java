package beastmap.logger.mut;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.Nucleotide;
import beast.base.evolution.tree.Node;
import beastmap.evolution.StochasticMapper;
import beastmap.logger.BranchSubstLogger;
import beastmap.util.Mutation;
import beastmap.util.MutationUtils;

@Description("Count the number of TT or AA dimers that get broken")
public class ThymineDimerBreakCounter extends BranchSubstLogger {
	
	
	final int a = 0;
	final int c = 1;
	final int g = 2;
	final int t = 3;
	final int gap = 4;

	
	@Override
    public void initAndValidate() {
		super.initAndValidate();
	}
	
	@Override
	public double getFilteredMutationSummary(List<Mutation> mutations, Node child) {
		
		
		// Multispecies coalescent? -- start the count at the point where this gene branch extends into the species node parent
		double maxHeight = child.getParent().getHeight() - mutations.get(0).getTime();
		double minHeight = child.getParent().getHeight() - mutations.get(mutations.size()-1).getTime();
		
		StochasticMapper mapper = samplerInput.get() != null ? samplerInput.get() : truthInput.get();
		int[] childSequence = MutationUtils.getStatesForBranchAtTime(mapper, mapper.getTree(), child, minHeight, true);
		int[] parentSequence = MutationUtils.getStatesForBranchAtTime(mapper, mapper.getTree(), child, maxHeight, false);
		int[] nextSequence = Arrays.copyOf(parentSequence, parentSequence.length);
		
		int count = 0;
		
		// Going forward in time
		for (Mutation mut : mutations) {
			
			//double height = child.getParent().getHeight() - mut.getTime();
			int siteNr = mut.getSiteNr();
			int from = mut.getFrom();
			int to = mut.getTo();
			nextSequence[siteNr] = to;
			
			if (parentSequence[siteNr] != from) {
				Log.warning("Dev error 98765 " + parentSequence[siteNr] + " != " + from);
			}
			
			
			int nextSite = siteNr + 1;
			int prevSite = siteNr - 1;
			
			
			if (prevSite >= 0) {
				
				
				// ?
			}
			
			if (nextSite < parentSequence.length) {
				
				// ?
			}
			
			// TODO -- not sure how to handle indels
			
			// Did we do a TT -> xT/tY ?
			
			
			
			
		}
		

		
		return count;
		
	}

	@Override
	public String getName() {
		return "transition" + "." + this.getID();
	}

	@Override
	public boolean canHandleDataType(DataType dataType) {
		return dataType instanceof Nucleotide;
	}
	
	
	@Override
	public boolean isSummable() {
		return true;
	}

}
