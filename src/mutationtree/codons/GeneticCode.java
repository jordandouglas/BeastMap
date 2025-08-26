

package mutationtree.codons;


import beast.base.core.Description;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.Nucleotide;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A set of standard genetic codes.
 * Note: the codeState here is same the states in {@link Codon Codon} data type.
 * See mapCodeToStateSet.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 *
 * Merged and modified from BEAST 1 GeneticCode, AminoAcids, Nucleotides.
 */
@Description("A set of standard genetic codes.")
public class GeneticCode {

//    public static final String GENETIC_CODE = "geneticCode";
    
    /**
     * Constants used to refer to the built in code tables
     */
    public static final int UNIVERSAL_ID = 0;
    public static final int VERTEBRATE_MT_ID = 1;
    public static final int YEAST_ID = 2;
    public static final int MOLD_PROTOZOAN_MT_ID = 3;
    public static final int MYCOPLASMA_ID = 4;
    public static final int INVERTEBRATE_MT_ID = 5;
    public static final int CILIATE_ID = 6;
    public static final int ECHINODERM_MT_ID = 7;
    public static final int EUPLOTID_NUC_ID = 8;
    public static final int BACTERIAL_ID = 9;
    public static final int ALT_YEAST_ID = 10;
    public static final int ASCIDIAN_MT_ID = 11;
    public static final int FLATWORM_MT_ID = 12;
    public static final int BLEPHARISMA_NUC_ID = 13;
    public static final int NO_STOPS_ID = 14;

    /**
     * Standard genetic code tables from GENBANK
     * Nucleotides go A, C, G, T - Note: this is not the order used by the Genbank web site
     * With the first codon position most significant (i.e. AAA, AAC, AAG, AAT, ACA, etc.).
     */
    public static final String[] GENETIC_CODE_TABLES = {
        // Universal
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
        // Vertebrate Mitochondrial
        "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
        // Yeast
        "KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
        // Mold Protozoan Mitochondrial
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
        // Mycoplasma
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
        // Invertebrate Mitochondrial
        "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
        // Ciliate
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF",
        // Echinoderm Mitochondrial
        "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
        // Euplotid Nuclear
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF",
        // Bacterial
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
        // Alternative Yeast
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
        // Ascidian Mitochondrial
        "KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
        // Flatworm Mitochondrial
        "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF",
        // Blepharisma Nuclear
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF",
        // No stops
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYYQYSSSSWCWCLFLF"
    };

    /**
     * Names of the standard genetic code tables from GENBANK
     */
    public static final String[] GENETIC_CODE_NAMES = {
        "universal", "vertebrateMitochondrial", "yeast", "moldProtozoanMitochondrial",
        "mycoplasma", "invertebrateMitochondrial", "ciliate", "echinodermMitochondrial",
        "euplotidNuclear", "bacterial", "alternativeYeast", "ascidianMitochondrial",
        "flatwormMitochondrial", "blepharismaNuclear", "noStops"
    };

    /**
     * Descriptions of the standard genetic code tables from GENBANK
     */
    public static final String[] GENETIC_CODE_DESCRIPTIONS = {
        "Universal", "Vertebrate Mitochondrial", "Yeast", "Mold Protozoan Mitochondrial",
        "Mycoplasma", "Invertebrate Mitochondrial", "Ciliate", "Echinoderm Mitochondrial",
        "Euplotid Nuclear", "Bacterial", "Alternative Yeast", "Ascidian Mitochondrial",
        "Flatworm Mitochondrial", "Blepharisma Nuclear", "Test case with no stop codons"
    };

    public static final GeneticCode UNIVERSAL = new GeneticCode(UNIVERSAL_ID);
    public static final GeneticCode VERTEBRATE_MT = new GeneticCode(VERTEBRATE_MT_ID);
    public static final GeneticCode YEAST = new GeneticCode(YEAST_ID);
    public static final GeneticCode MOLD_PROTOZOAN_MT = new GeneticCode(MOLD_PROTOZOAN_MT_ID);
    public static final GeneticCode MYCOPLASMA = new GeneticCode(MYCOPLASMA_ID);
    public static final GeneticCode INVERTEBRATE_MT = new GeneticCode(INVERTEBRATE_MT_ID);
    public static final GeneticCode CILIATE = new GeneticCode(CILIATE_ID);
    public static final GeneticCode ECHINODERM_MT = new GeneticCode(ECHINODERM_MT_ID);
    public static final GeneticCode EUPLOTID_NUC = new GeneticCode(EUPLOTID_NUC_ID);
    public static final GeneticCode BACTERIAL = new GeneticCode(BACTERIAL_ID);
    public static final GeneticCode ALT_YEAST = new GeneticCode(ALT_YEAST_ID);
    public static final GeneticCode ASCIDIAN_MT = new GeneticCode(ASCIDIAN_MT_ID);
    public static final GeneticCode FLATWORM_MT = new GeneticCode(FLATWORM_MT_ID);
    public static final GeneticCode BLEPHARISMA_NUC = new GeneticCode(BLEPHARISMA_NUC_ID);
    public static final GeneticCode NO_STOPS = new GeneticCode(NO_STOPS_ID);

//    public static final GeneticCode[] GENETIC_CODES = {
//        UNIVERSAL, VERTEBRATE_MT, YEAST, MOLD_PROTOZOAN_MT, MYCOPLASMA, INVERTEBRATE_MT,
//        CILIATE, ECHINODERM_MT, EUPLOTID_NUC, BACTERIAL, ALT_YEAST, ASCIDIAN_MT,
//        FLATWORM_MT, BLEPHARISMA_NUC, NO_STOPS
//    };


    public GeneticCode(int geneticCodeId) {
        this.geneticCodeId = geneticCodeId;
        codeTable = GENETIC_CODE_TABLES[geneticCodeId];

        assert getCodeTableNoStopCodons().length() == getStateCount();
    }

    // rm stop codons from states
    public int getStateCount() {
        return getCodeTableLength() - getStopCodonCount();
    }

    /**
     * Return the genetic code table currently used
     * including stop codons
     */
    public String getCodeTable() {
        return codeTable;
    }

    /**
     * Return the genetic code table currently used
     * excluding stop codons
     */
    public String getCodeTableNoStopCodons() {
        return codeTable.replaceAll("\\*", "");
    }

    /**
     * 64, which is also the state count of Codon.
     * Return the length of the code string in the currently used genetic code table,
     * also the number of defined {@link Codon#CODON_TRIPLETS triplets}.
     * @see DataType#getStateCount()
     * @see DataType#stringToEncoding(String)
     */
    public int getCodeTableLength() {
        return codeTable.length();
    }

    /**
     * Returns the name of the genetic code
     */
    public String getName() {
        return GENETIC_CODE_NAMES[geneticCodeId];
    }
    
    /**
     * Returns the description of the genetic code
     */
    public String getDescription() {
        return GENETIC_CODE_DESCRIPTIONS[geneticCodeId];
    }

    public static GeneticCode findByDescription(String desc) {
        int idx = Arrays.asList(GENETIC_CODE_DESCRIPTIONS).indexOf(desc);
        if (idx >= 0) return findByName(GENETIC_CODE_NAMES[idx]);
        throw new RuntimeException("Unknown genetics code");
    }

    public static GeneticCode findByName(String codeStr) {
        GeneticCode geneticCode = null;
        if (codeStr.equals(GeneticCode.UNIVERSAL.getName())) {
            geneticCode = GeneticCode.UNIVERSAL;
        } else if (codeStr.equals(GeneticCode.VERTEBRATE_MT.getName())) {
            geneticCode = GeneticCode.VERTEBRATE_MT;
        } else if (codeStr.equals(GeneticCode.YEAST.getName())) {
            geneticCode = GeneticCode.YEAST;
        } else if (codeStr.equals(GeneticCode.MOLD_PROTOZOAN_MT.getName())) {
            geneticCode = GeneticCode.MOLD_PROTOZOAN_MT;
        } else if (codeStr.equals(GeneticCode.INVERTEBRATE_MT.getName())) {
            geneticCode = GeneticCode.INVERTEBRATE_MT;
        } else if (codeStr.equals(GeneticCode.CILIATE.getName())) {
            geneticCode = GeneticCode.CILIATE;
        } else if (codeStr.equals(GeneticCode.ECHINODERM_MT.getName())) {
            geneticCode = GeneticCode.ECHINODERM_MT;
        } else if (codeStr.equals(GeneticCode.EUPLOTID_NUC.getName())) {
            geneticCode = GeneticCode.EUPLOTID_NUC;
        } else if (codeStr.equals(GeneticCode.BACTERIAL.getName())) {
            geneticCode = GeneticCode.BACTERIAL;
        } else if (codeStr.equals(GeneticCode.ALT_YEAST.getName())) {
            geneticCode = GeneticCode.ALT_YEAST;
        } else if (codeStr.equals(GeneticCode.ASCIDIAN_MT.getName())) {
            geneticCode = GeneticCode.ASCIDIAN_MT;
        } else if (codeStr.equals(GeneticCode.FLATWORM_MT.getName())) {
            geneticCode = GeneticCode.FLATWORM_MT;
        } else if (codeStr.equals(GeneticCode.BLEPHARISMA_NUC.getName())) {
            geneticCode = GeneticCode.BLEPHARISMA_NUC;
        } else if (codeStr.equals(GeneticCode.NO_STOPS.getName())) {
            geneticCode = GeneticCode.NO_STOPS;
        } else {
            throw new RuntimeException("Unknown genetics code");
        }
        return geneticCode;
    }

    protected Nucleotide nucleotide = new Nucleotide(); // char <=> state

    /**
     * Get a Nucleotide state from {@link Nucleotide#codeMap codeMap}.
     * <code>codeMap = "ACGTURYMWSKBDHVNX" + GAP_CHAR + MISSING_CHAR</code>
     * @see DataType#stringToEncoding(String)
     * @param c a Nucleotide character in codeMap
     * @return Nucleotide state, also the index in {@link Nucleotide#codeMap codeMap}.
     */
    public int getNucleotideState(char c) {
        return nucleotide.stringToEncoding(String.valueOf(c)).get(0);
    }

    /**
     * Infer char from {@link Nucleotide#codeMap codeMap}.
     * @see DataType#encodingToString(int[])
     * @param state Nucleotide state, also the index in {@link Nucleotide#codeMap codeMap}.
     * @return a Nucleotide character
     */
    public char getNucleotideChar(int state) {
        return nucleotide.encodingToString(new int[]{state}).charAt(0);
    }

    /**
     * This table maps amino acid characters into state codes (0-25).
     * Amino Acids go ACDEFGHIKLMNPQRSTVWYBZX*?-,
     * Other letters; j, o, and u are mapped to ?
     * *, ? and - are mapped to themselves
     * All other chars are mapped to -
     */
    public static final int[] AMINOACID_STATES = {
            25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,    // 0-15
            25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,    // 16-31
    //                                     *        -
            25,25,25,25,25,25,25,25,25,25,23,25,25,25,25,25,    // 32-47
    //                                                    ?
            25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,24,    // 48-63
    //          A  B  C  D  E  F  G  H  I  j  K  L  M  N  o
            25, 0,20, 1, 2, 3, 4, 5, 6, 7,24, 8, 9,10,11,24,    // 64-79
    //       P  Q  R  S  T  u  V  W  X  Y  Z
            12,13,14,15,16,24,17,18,22,19,21,25,25,25,25,25,    // 80-95
    //          A  B  C  D  E  F  G  H  I  j  K  L  M  N  o
            25, 0,20, 1, 2, 3, 4, 5, 6, 7,24, 8, 9,10,11,24,    // 96-111
    //       P  Q  R  S  T  u  V  W  X  Y  Z
            12,13,14,15,16,24,17,18,22,19,21,25,25,25,25,25     // 112-127
    };

    /**
     * A table to translate state numbers (0-25) into one letter codes.
     */
    public static final char[] AMINOACID_CHARS= {
            'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
            'S','T','V','W','Y','B','Z','X',
            Codon.STOP_CHARACTER,DataType.MISSING_CHAR,DataType.GAP_CHAR
    };

    /**
     * Get character corresponding to a given state mapping to
     * {@link #AMINOACID_STATES AMINOACID_STATES}.
     * @param aaState the state numbers in {@link #AMINOACID_STATES AMINOACID_STATES},
     *              which is different to the state in {@link Codon Codon}.
     * @return
     */
    public char getAminoAcidChar(int aaState) {
        return AMINOACID_CHARS[aaState];
    }


    /**
     * @param codeState the character index from the string of
     *                  the currently used genetic code table.
     *                  Note: after excluding stop codons.
     * @return an Amino Acid character in codeTableNoStopCodons.
     *         if codeState >= getStateCount(), return -.
     */
    public char getAminoAcidNoStopCodons(int codeState) {
        if (codeState >= getStateCount())
            return DataType.GAP_CHAR;
        return getCodeTableNoStopCodons().charAt(codeState);
    }

    /**
     * Return the Amino Acid state ({@link #AMINOACID_STATES AMINOACID_STATES})
     * corresponding to a given codon (triplet) state.
     * Note: the state is the canonical state (generated combinatorially),
     * which is different to BEAST 2 Aminoacid states, and which can be
     * from {@link Codon#getCodonState(int,int,int)} getCodonState}.
     * @param codeState the state in {@link Codon Codon} triplets
     * @see #AMINOACID_STATES
     * @return Amino Acid state
     */
    public int getAminoAcidState(int codeState) {
//        if (codeState == Codon.UNKNOWN_STATE)
//            return AMINOACID_STATES[Aminoacid.MISSING_CHAR];
//        else if (codeState == Codon.GAP_STATE)
//            return AMINOACID_STATES[Aminoacid.GAP_CHAR];

        char aa = getAminoAcidNoStopCodons(codeState);
        return AMINOACID_STATES[aa];
    }

    //++++++++ Note: codeState != index +++++++++

    /**
     * Get amino acid corresponding to a given index which includes stop codons.
     * @see #GENETIC_CODE_TABLES
     * @param index the Amino Acid character index in the currently used genetic code table
     *              including stop codons, also the same index to the {@link Codon Codon} triplets.
     * @return an Amino Acid character in codeTable.
     *         if codeState >= codeTable.length, return -.
     */
    public char getAminoAcidFromCodeTable(int index) {
        if (index >= getCodeTableLength()) // Ambiguous handled in Codon TODO >= 64
            return DataType.GAP_CHAR;
        return codeTable.charAt(index);
    }


//    /**
//     * Note that the state is the canonical state (generated combinatorially)
//     * @return whether the codeState is a stop codon
//     */
//    public boolean isStopCodon(int codeState) {
//        return (getAminoAcidState(codeState) == Codon.STOP_STATE);
//    }

    /**
     * @param index the Amino Acid character index in the currently used genetic code table
     *        including stop codons, also the same index to the {@link Codon Codon} triplets.
     */
    public boolean isStopCodonIndex(int index) {
        if (index >= getCodeTableLength())
            // TODO >= 64,
            // which must expand the standard genetic code table to handle ambiguous nucleotides
            throw new UnsupportedOperationException("Ambiguous nucleotides are not supported at the moment !");
//            return false;
        return (getAminoAcidFromCodeTable(index) == Codon.STOP_CHARACTER);
    }

    /**
     * Return the codon states of all stop codon in
     * the currently used genetic code table.
     */
    public int[] getStopCodonIndices() {
        List<Integer> indices = new ArrayList<>();
        for (int i = 0; i < getCodeTableLength(); i++) {
            if (getAminoAcidFromCodeTable(i) == Codon.STOP_CHARACTER) {
                indices.add(i);
            }
        }
        // >= Java 1.8
        return indices.stream().mapToInt(i -> i).toArray();
    }

    /**
     * Return the number of stop codon "*"
     * in the currently used genetic code table.
     */
    public int getStopCodonCount() {
        int count = 0;
        for (int i = 0; i < getCodeTableLength(); i++) {
            if (getAminoAcidFromCodeTable(i) == Codon.STOP_CHARACTER)
                count++;
        }
        return count;
    }

    protected int geneticCodeId;
    protected String codeTable;

    @Override
    public String toString() {
        return getName();
    }

	public int getGeneticCodeID() {
		return geneticCodeId;
	}
}
