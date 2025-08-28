package mutationtree.operator;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;


@Description("Scales one branch lengths in a tree (in contrast with Scale which scales a height not a length)")
public class BranchLengthScaler extends TreeOperator {
	
	//public final Input<RealParameter> rateInput = new Input<>("rate", "if specified, this parameter is scaled in the opposite direction as the tree");
	public final Input<KernelDistribution> kernelDistributionInput = new Input<>("kernelDistribution", "provides sample distribution for proposals", 
	    		KernelDistribution.newDefaultKernelDistribution());
	public final Input<Double> scaleFactorInput = new Input<>("scaleFactor", "scaling factor: range from 0 to 1. Close to zero is very large jumps, close to 1.0 is very small jumps.", 0.75);
	final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);
	
	
	protected KernelDistribution kernelDistribution;
	private double scaleFactor;

    @Override
	public void initAndValidate() {
    	scaleFactor = scaleFactorInput.get();
    	kernelDistribution = kernelDistributionInput.get();
	}
    
	protected double getScaler(int i, double value) {
    	return kernelDistribution.getScaler(i, value, getCoercableParameterValue());
	}

	@Override
	public double proposal() {
		
		
		double logHR = 0;
		
		// Scale branch lengths
		final Tree tree = (Tree)InputUtil.get(treeInput, this);
		final double scale = getScaler(0, Double.NaN);
		int nbranches = tree.getNodeCount()-1;
		
		// Sample a branch to move
		int branchNr = Randomizer.nextInt(nbranches);
		Node x = tree.getNode(branchNr);
		double oldLength = x.getLength();
		double newLength = oldLength * scale;
		
		// Adjust height of clade
		double dt = oldLength - newLength;
        List<Node> childNodes = new ArrayList<>();
        x.getAllChildNodesAndSelf(childNodes);
        for (Node node : childNodes) {
        	node.setHeight(node.getHeight() + dt);
        }
        
        logHR = Math.log(scale);
       
		
		
		// Translate the tree so that youngest tip is at height 0
		double minHeight = Double.POSITIVE_INFINITY;
		for (Node node : treeInput.get().getNodesAsArray()) {
			minHeight = Math.min(minHeight, node.getHeight());
		}
		for (Node node : treeInput.get().getNodesAsArray()) {
			node.setHeight(node.getHeight() - minHeight);
		}
				
       
		return logHR;
		
	}
	
	@Override
	public void optimize(double logAlpha) {
	    // must be overridden by operator implementation to have an effect
		if (optimiseInput.get()) {
	        double delta = calcDelta(logAlpha);
	        double scaleFactor = getCoercableParameterValue();
	        delta += Math.log(scaleFactor);
	        scaleFactor = Math.exp(delta);
	        setCoercableParameterValue(scaleFactor);
		}
	}
	  

	    
	@Override
	public double getTargetAcceptanceProbability() {
		return 0.3;
	}
	    


	@Override
	public String getPerformanceSuggestion() {
	    double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
	    double targetProb = getTargetAcceptanceProbability();
	
	    double ratio = prob / targetProb;
	    if (ratio > 2.0) ratio = 2.0;
	    if (ratio < 0.5) ratio = 0.5;
	
	    // new scale factor
	    double newWindowSize = getCoercableParameterValue() * ratio;
	
	    DecimalFormat formatter = new DecimalFormat("#.###");
	    if (prob < 0.10 || prob > 0.40) {
	        return "Try setting scale factor to about " + formatter.format(newWindowSize);
	    } else return "";
	}
	
	
	
    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        scaleFactor = Math.max(value, 0.0);
    }



}
