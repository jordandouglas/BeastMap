package mutationtree.logger;

import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import mutationtree.evolution.AncestralSequenceTreeLikelihood;

@Description("Reconstructs sequences at internal nodes and logs them in NEXUS format")
public class AncestralSequenceLogger extends AncestralSequenceTreeLikelihood implements Function, Loggable {
	final public Input<Boolean> logIndividualInput = new Input<>("logIndividualSites", "if true, tree log gets one entry for every site, "
			+ "if false complete sequence is logged", false);
	
	
	int [] siteStates;
	
	@Override
	public void init(PrintStream out) {
		((Tree)treeInput.get()).init(out);
		
		dataType =  dataInput.get().getDataType();
		siteStates = new int[dataInput.get().getSiteCount()];	
	}
	
	@Override
	public void log(long sample, PrintStream out) {
		calculateLogP();
		redrawAncestralStates();
        out.print("tree STATE_" + sample + " = ");
        TreeInterface tree = treeInput.get();
        tree.getRoot().sort();
        out.print(toNewick(tree.getRoot()));
        out.print(";");
	}

	
    String toNewick(Node node) {
        StringBuffer buf = new StringBuffer();
        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(toNewick(node.getLeft()));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toNewick(node.getRight()));
            }
            buf.append(")");
        } else {
            buf.append(node.getNr() + 1);
        }

	    int [] patternstates = getStatesForNode(treeInput.get(), node);
	    for (int i = 0; i < siteStates.length; i++) {
	    	siteStates[i] = patternstates[dataInput.get().getPatternIndex(i)];
	    }
	    String seq = dataType.encodingToString(siteStates);
	    if (logIndividualInput.get()) {
	    	buf.append("[&");
	    	for (int k = 0; k < seq.length(); k++) {
	    		buf.append((k > 0 ? "," : "") + tagInput.get()
	    		+ (k < 10 ? "0":"")
	    		+ (k < 100 ? "0":"")
	    		+ k + "=\"" + seq.charAt(k) + "\"");
	    	}
	    	buf.append("]");
	    	
	    } else {
	    	buf.append("[&" + tagInput.get() + "=\"" + seq + "\"]");
	    }

	    buf.append(':');
        buf.append(node.getLength());
        return buf.toString();
    }
    
	@Override
	public void close(PrintStream out) {
		((Tree)treeInput.get()).close(out);
	}	

}
