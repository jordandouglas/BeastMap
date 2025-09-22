package beastmap.logger.mut;

import java.util.List;

import beast.base.core.Description;
import beast.base.evolution.tree.Node;
import beastmap.util.Mutation;

@Description("Counts the number of times an amino acid mutates into a different member of the same functional class")
public class AminoAcidClassRemains extends AminoAcidClassChanges {
	
	
	
	
	@Override
	public double getFilteredMutationSummary(List<Mutation> mutations, List<Mutation> mutationsUnconditional, Node node) {
		
		int count = 0;
		for (Mutation mut : mutations) {
			
			int from = mut.getFrom();
			int to = mut.getTo();
			
			if (!aaClasses.containsKey(from) || !aaClasses.containsKey(to)) continue;
			
			int fromClass = aaClasses.get(from);
			int toClass = aaClasses.get(to);
			
			
			// Count the similarities
			if (fromClass == toClass) count ++;
			
		}
		
		return count;
		
	}

	@Override
	public String getName() {
		return "AAClassRemains" + "." + this.getID();
	}

	
	

}
