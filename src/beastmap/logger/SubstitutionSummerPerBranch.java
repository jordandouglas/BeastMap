package beastmap.logger;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;

@Description("Sums multiple stochastic mappers together, on a per-branch basis")
public class SubstitutionSummerPerBranch extends CalculationNode implements StochasticMapProperty, Loggable, Function {

	
	final public Input<Boolean> weightInput = new Input<>("weight", "weight sums by alignment length?", false);
	final public Input<List<BranchSubstLogger>> counterInput = new Input<>("counter", "counter to sum", new ArrayList<>());
	
	Tree tree;
	
	
	@Override
	public void initAndValidate() {
		
		if (counterInput.get().isEmpty()) {
			throw new IllegalArgumentException("Please provide at least one counter");
		}
		
		Tree tree = counterInput.get().get(0).getTree();
		for (BranchSubstLogger logger : counterInput.get()) {
			Tree tree2 = logger.getTree();
			if (tree != tree2) {
				throw new IllegalArgumentException("Please ensure all trees are the same " + tree.getID() + " != " + tree2.getID());
			}
		}
		
		this.tree = tree;
		
	}

	@Override
	public void init(PrintStream out) {
		
		for (int i = 0; i < this.getDimension(); i ++) {
			out.print(getName() + "." + i + "\t");
		}
	}

	@Override
	public void log(long sample, PrintStream out) {
		
		sampleMutations(sample);
		for (int i = 0; i < this.getDimension(); i ++) {
			out.print(getArrayValue(i) + "\t");
		}
		
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public String getName() {
		// TODO Auto-generated method stub
		return this.getID();
	}

	@Override
	public void sampleMutations(long sampleNr) {
		for (BranchSubstLogger logger : counterInput.get()) {
			logger.sampleMutations(sampleNr);
		}
	}

	@Override
	public Object getPropertyOfNode(Node node) {
		return getArrayValue(node.getNr());
	}

	@Override
	public int getDimension() {
		return this.tree.getNodeCount();
	}

	@Override
	public double getArrayValue(int dim) {
		
		//double wsum = 0;
		double total = 0;
		for (BranchSubstLogger logger : counterInput.get()) {
			
			//double w = 1;
//			if (weightInput.get()) {
//				w = logger.getSiteAndPatternCount();
//			}
			total += logger.getArrayValue(dim);
			//wsum += w;
		}
		
		
		//if (weightInput.get()) total = total / wsum;
		return total;
		
	}

	
}

