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
	
	
	final public Input<Boolean> incodonInput = new Input<>("incodon", "require that we only consider two contiguous bases if they are in the same codon, starting from site 1", false);
	
	
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
		
		
		final int adenine = 0;
		final int thymine = 3;
		
		if (mutations.isEmpty()) return 0;
		
		
		// Multispecies coalescent? -- start the count at the point where this gene branch extends into the species node parent
		double maxHeight = child.getParent().getHeight() - mutations.get(0).getTime();
		double minHeight = child.getParent().getHeight() - mutations.get(mutations.size()-1).getTime();
		
		StochasticMapper mapper = samplerInput.get() != null ? samplerInput.get() : truthInput.get();
		int[] childSequence = MutationUtils.getStatesForBranchAtTime(mapper, mapper.getTree(), child, minHeight, true);
		int[] parentSequence = MutationUtils.getStatesForBranchAtTime(mapper, mapper.getTree(), child, maxHeight, false);
		
		int count = 0;
		
		
		// Iterate across all sites, and pay attention to this site and the one immediately downstream from it
		for (int siteNr = 0; siteNr < childSequence.length-1; siteNr++) {
			
			int nextSiteNr = siteNr+1;

			// Discount any AA or TT that cross from one codon into the next, just in case there are introns/gaps
			if (incodonInput.get()) {
				if (siteNr % 3 == 2) {
					continue;
				}
			}
			
			boolean dimerHasBroken = false;
			int stateOfBreak = -1;
			
			int state1 = parentSequence[siteNr];
			int state2 = parentSequence[nextSiteNr];
			
			// Going forward in time
			for (Mutation mut : mutations) {
				
				
				// We only care about mutations in these two sites here
				if (mut.getSiteNr() != siteNr && mut.getSiteNr() != nextSiteNr) {
					continue;
				}
				
				
				// One of these just changed
				int newState1 = state1;
				int newState2 = state2;
				if (mut.getSiteNr() == siteNr) newState1 = mut.getTo();
				if (mut.getSiteNr() == nextSiteNr) newState2 = mut.getTo();
				
				
				// Sanity check
				if (mut.getSiteNr() == siteNr && state1 != mut.getFrom()) {
					Log.warning("Dev error 98764.1 " + state1 + " != " + mut.getFrom());
				}
				if (mut.getSiteNr() == nextSiteNr && state2 != mut.getFrom()) {
					Log.warning("Dev error 98764.2 " + state2 + " != " + mut.getFrom());
				}
					
				
				// Count all occurrences of TT->XY or AA->XY where X differs from the initial state, and Y can be anything
				if (!dimerHasBroken) {
					
					
					// Second position has broken. Wait until the first position breaks
					if (state1 == state2 && (state2 == thymine || state2 == adenine) && newState2 != state2) {
						stateOfBreak = state1;
						dimerHasBroken = true;
					}
					
					
					// The first position has already broken, increment count
					else if (state1 == state2 && (state2 == thymine || state2 == adenine) && newState1 != state1) {
						count ++;
					}
					
					
				}
				
				else {
					
					
					
					// Break has repaired back to TT or AA
					if (newState1 == newState2 && newState1 == stateOfBreak) {
						stateOfBreak = -1;
						dimerHasBroken = false;
					}
					
					// Break has completed
					else if (newState1 != stateOfBreak) {
						count ++;
						stateOfBreak = -1;
						dimerHasBroken = false;
					}
					
				}
				
				
				if (state1 == thymine && state2 == thymine && newState1 != thymine) {
					count ++;
				}
				

				// This is a more stringent test that requires both parts mutate not just one of them
//				if (!dimerHasBroken) {
//					
//					// Did it break just now?
//					if (state1 == thymine && state2 == thymine && (newState1 != thymine || newState2 != thymine)) {
//						dimerHasBroken = true;
//					}
//					
//				}
//				
//				else {
//					
//					
//					
//					// Did the break repair back to TT?
//					if (newState1 == thymine && newState2 == thymine) {
//						
//						
//						// Go back to the start. No thymine breaks
//						dimerHasBroken = false;
//						
//					}
//					
//					
//					
//					// Has the double mutation finished? Increment the count
//					else if (newState1 != thymine && newState2 != thymine) {
//						
//						count ++;
//						dimerHasBroken = false;
//						
//					}
//					
//				}
//				
				
				state1 = newState1;
				state2 = newState2;
				
				
			}
			
			
			// Sanity check now that we are at the bottom of the branch
			if (childSequence[siteNr] != state1) {
				Log.warning("Dev error 98765.1 " + childSequence[siteNr] + " != " + state1);
			}
			if (childSequence[nextSiteNr] != state2) {
				Log.warning("Dev error 98765.2 " + childSequence[nextSiteNr] + " != " + state2);
			}
			
			
			
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
