package beastmap.indel;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.TreeSet;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.datatype.Binary;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.TwoStateCovarion;


@Description("Represents an alignment using its gaps https://www.jstor.org/stable/pdf/2585224.pdf")
public class SimpleIndelCodingAlignment extends Alignment {
	
    final public Input<Alignment> dataInput = new Input<>("data", "the original alignment", Validate.REQUIRED);
    

    
    SimpleIndelCoding indels;
    int gapChar;
	
	@Override
    public void initAndValidate() {
		
		//dataTypeInput.set("binary");
		
		
		if (userDataTypeInput.get() != null) {
			DataType dt = userDataTypeInput.get();
			if (! (dt instanceof Binary) && !(dt instanceof TwoStateCovarion)) {
				throw new IllegalArgumentException("datatype " + userDataTypeInput.get().getID() + " must be Binary or TwoStateCovarion");
			}
		}else {
			dataTypeInput.setValue("binary", this);
		}
		
		stateCountInput.setValue(2, this);
		
		gapChar = dataInput.get().getDataType().stringToEncoding(""+DataType.GAP_CHAR).get(0);
		List<Sequence> sequencesOrig = dataInput.get().sequenceInput.get();
		indels = new SimpleIndelCoding(sequencesOrig, dataInput.get().getDataType(), this.gapChar);
		
		//int[][] indelCodingMatrix = prepareIndelMatrix(sequences);
		int[][] indelCodingMatrix = indels.getSICMatrix();
		
		for (int seqNr = 0; seqNr < sequencesOrig.size(); seqNr++) {
			Sequence sequence = sequencesOrig.get(seqNr);
			String seq = "";
			for (int pos = 0; pos < indelCodingMatrix[seqNr].length; pos++) {
				int val = indelCodingMatrix[seqNr][pos];
				String c = val == 1 ? "1" : val == 0?  "0" : "-";
				seq += c;
			}
			//System.out.println("indels " + sequence.getTaxon() + ":\n\t" + seq);
			System.out.println(">" + sequence.getTaxon() + ":\n" + seq.replaceAll("0", "A").replaceAll("1", "G"));
			
			Sequence indelSeq = new Sequence(sequence.getTaxon(), seq);
			sequenceInput.setValue(indelSeq, this);
		}
		
    	if (types.size() == 0) {
    		findDataTypes();    		
    	}
    	
    	
    	
    	
    	// Test the back-conversion
//    	List<Integer> x = sequenceInput.get().get(0).getSequence(new Binary());
//    	System.out.println(">recon" + sequenceInput.get().get(0).getTaxon());
//    	expandIndelCoding(x);
    	
    	
        super.initAndValidate();
        
        
	}
	
	
	
	public boolean[] expandIndelCoding(List<Integer> simpleIndelCoding ){
		int[] arr = new int[simpleIndelCoding.size()];
		for (int i = 0; i < arr.length; i ++) arr[i] = simpleIndelCoding.get(i);
		return expandIndelCoding(arr);
	}
	
	/**
	 * Expand a compact binary representation into an alignment of the original size filled with 0's (for gaps) and 1's (for not gaps)
	 * @return
	 */
	public boolean[] expandIndelCoding(int[] simpleIndelCoding ){
		
		boolean[] expanded = indels.expandIndelSequence(simpleIndelCoding);
		
//		for (int pos = 0; pos < expanded.length; pos++) {
//			System.out.print(expanded[pos] ? "N" : "-");
//		}
//		System.out.println();
		
		return expanded;
		
	}
	
	public class SimpleIndelCoding {

	    public class Interval {
	        public final int start;
	        public final int end;
	        public Interval(int s, int e) { start = s; end = e; }
	    }

	    private final List<Interval> blocks;
	    private final int[][] sicMatrix;
	    private final int alnLength;

	    /** Constructor: compute superset-aware SIC from sequences */
	    public SimpleIndelCoding(List<Sequence> sequences, DataType dtOriginal, int gapChar) {
	        if (sequences == null || sequences.isEmpty())
	            throw new IllegalArgumentException("Sequences cannot be empty");

	        this.alnLength = sequences.get(0).getSequence(dtOriginal).size();
	        

	        // 1. collect gap runs for each sequence
	        List<List<Interval>> gapsPerSeq = new ArrayList<>();
	        for (Sequence seq : sequences) {
	            List<Integer> s = seq.getSequence(dtOriginal);
	            List<Interval> runs = new ArrayList<>();
	            boolean inGap = false;
	            int start = -1;
	            for (int i = 0; i < alnLength; i++) {
	                //boolean isGap = dtOriginal.isAmbiguousCode(s.get(i));
	                boolean isGap = s.get(i) == gapChar; // Gaps only not ambigs
	                if (isGap && !inGap) {
	                    inGap = true;
	                    start = i;
	                } else if (!isGap && inGap) {
	                    inGap = false;
	                    runs.add(new Interval(start, i - 1));
	                }
	            }
	            if (inGap) runs.add(new Interval(start, alnLength - 1));
	            gapsPerSeq.add(runs);
	        }

	        // 2. collect unique blocks
	        TreeSet<Interval> uniqueBlocks = new TreeSet<>(Comparator
	                .comparingInt((Interval iv) -> iv.start)
	                .thenComparingInt(iv -> iv.end));
	        for (List<Interval> seqRuns : gapsPerSeq) {
	            uniqueBlocks.addAll(seqRuns);
	        }
	        this.blocks = new ArrayList<>(uniqueBlocks);

	        // 3. build SIC matrix (superset-aware)
	        int nSeq = sequences.size();
	        sicMatrix = new int[nSeq][blocks.size()];

	        for (int b = 0; b < blocks.size(); b++) {
	            Interval block = blocks.get(b);
	            for (int s = 0; s < nSeq; s++) {
	                List<Interval> runs = gapsPerSeq.get(s);
	                boolean exactMatch = false;
	                boolean strictSuperset = false;

	                for (Interval r : runs) {
	                    if (r.start == block.start && r.end == block.end) {
	                        exactMatch = true;
	                        break;
	                    } else if (r.start <= block.start && r.end >= block.end
	                               && !(r.start == block.start && r.end == block.end)) {
	                        strictSuperset = true;
	                    }
	                }

	                if (exactMatch) sicMatrix[s][b] = 1;
	                else if (strictSuperset) sicMatrix[s][b] = -1;
	                else sicMatrix[s][b] = 0;
	            }
	        }
	    }

	    public List<Interval> getBlocks() { return blocks; }
	    public int[][] getSICMatrix() { return sicMatrix; }
	    public int getAlignmentLength() { return alnLength; }

	    /** Expand a short SIC sequence into full-length alignment */
	    public boolean[] expandIndelSequence(int[] indelStates) {
	        if (indelStates.length != blocks.size())
	            throw new IllegalArgumentException("Length mismatch with blocks");
	        boolean[] full = new boolean[alnLength];
	        Arrays.fill(full, true); // residues

	        for (int b = 0; b < blocks.size(); b++) {
	            Interval block = blocks.get(b);
	            if (indelStates[b] == 1 || indelStates[b] == -1) {
	                for (int i = block.start; i <= block.end; i++)
	                    full[i] = false; // gap
	            }
	        }
	        return full;
	    }
	}

	public int getGapChar() {
		return gapChar;
	}
	
}




