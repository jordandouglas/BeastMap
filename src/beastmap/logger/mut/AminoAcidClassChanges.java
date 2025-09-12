package beastmap.logger.mut;

import java.util.HashMap;
import java.util.List;

import beast.base.core.Description;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.tree.Node;
import beastmap.logger.BranchSubstLogger;
import beastmap.util.Mutation;


@Description("Counts the number of changes between amino functional groups, largely based on BLOSUM62 classes")
public class AminoAcidClassChanges extends BranchSubstLogger {

	
	
	HashMap<Integer, Integer> aaClasses = new HashMap<>();
	
	final int nonpolar = 0;
	final int amide = 1;
	final int basic = 2;
	final int aromatic = 3;
	final int small = 4;
	final int cys = 5;
	final int pro = 6;
	final int gly = 7;
	
	
	 @Override
	 public void initAndValidate() {
		 super.initAndValidate();
		 
		 
		 // ACDEFGHIKLMNPQRSTVWY
		 
		 // Large non-polar
		 aaClasses.put(7, nonpolar); // I
		 aaClasses.put(9, nonpolar); // L
		 aaClasses.put(10, nonpolar); // M
		 aaClasses.put(17, nonpolar); // V
		 
		 // Amide/amine 
		 aaClasses.put(2, amide); // D
		 aaClasses.put(3, amide); // E
		 aaClasses.put(11, amide); // N
		 aaClasses.put(13, amide); // Q
		 
		 // Basic 
		 aaClasses.put(6, basic); // H
		 aaClasses.put(8, basic); // K
		 aaClasses.put(14, basic); // R
		 
		 // Aromatic 
		 aaClasses.put(4, aromatic); // F
		 aaClasses.put(18, aromatic); // W
		 aaClasses.put(19, aromatic); // Y
		 
		 // Small/polar 
		 aaClasses.put(0, small); // A
		 aaClasses.put(15, small); // S
		 aaClasses.put(16, small); // T
		 
		 // Cysteine
		 aaClasses.put(1, cys); // C
		 
		 // Glycine 
		 aaClasses.put(5, gly); // G
		 
		 // Proline 
		 aaClasses.put(12, pro); // P
		 
		 
	 }
	
	@Override
	public double getFilteredMutationSummary(List<Mutation> mutations, Node node) {
		
		int count = 0;
		for (Mutation mut : mutations) {
			
			int from = mut.getFrom();
			int to = mut.getTo();
			
			if (!aaClasses.containsKey(from) || !aaClasses.containsKey(to)) continue;
			
			int fromClass = aaClasses.get(from);
			int toClass = aaClasses.get(to);
			
			
			// Count the differences
			if (fromClass != toClass) count ++;
			
		}
		
		return count;
		
	}

	@Override
	public String getName() {
		return "AAClassChanges" + "." + this.getID();
	}

	@Override
	protected boolean canHandleDataType(DataType dataType) {
		return dataType instanceof Aminoacid;
	}

}
