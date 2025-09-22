package beastmap.logger;

import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;


@Description("Sums an entire vector of BranchSubstLogger terms")
public class SubstitutionSummer extends CalculationNode implements  StochasticMapProperty, Loggable {

	
	final public Input<BranchSubstLogger> counterInput = new Input<>("counter", "counter to log", Input.Validate.REQUIRED);
	
	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}


	@Override
	public void init(PrintStream out) {
		out.print(getName() + "\t");
	}

	@Override
	public void log(long sample, PrintStream out) {
		
		counterInput.get().sampleMutations(sample);
		
		double sum = 0 ;//(Double) getPropertyOfNode(null);
		for (int i = 0; i < counterInput.get().getDimension(); i ++) {
			sum += counterInput.get().getArrayValue(i);
		}
		
		out.print(sum + "\t");
	}

	@Override
	public void close(PrintStream out) {
		
	}


	@Override
	public String getName() {
		return "sum." + counterInput.get().getID();
	}


	@Override
	public void sampleMutations(long sampleNr) {
		counterInput.get().sampleMutations(sampleNr);
	}


	@Override
	public Object getPropertyOfNode(Node node) {
		double sum = 0;
		for (int i = 0; i < counterInput.get().getDimension(); i ++) {
			sum += counterInput.get().getArrayValue(i);
		}
		return sum;
	}



}
