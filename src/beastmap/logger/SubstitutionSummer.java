package beastmap.logger;

import java.io.PrintStream;
import beast.base.core.Function;
import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;


@Description("Sums an entire vector of BranchSubstLogger terms")
public class SubstitutionSummer extends CalculationNode implements StochasticMapProperty, Loggable {

	
	final public Input<Function> counterInput = new Input<>("counter", "counter to log", Input.Validate.REQUIRED);
	
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
		
		sampleMutations(sample);
		
		
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
		
		String id = "";
		Function f = counterInput.get();
		if (f instanceof BEASTInterface) {
			id = ((BEASTInterface)f).getID();
		}
		return "sum." + id;
	}


	@Override
	public void sampleMutations(long sampleNr) {
		
		Function f = counterInput.get();
		if (f instanceof StochasticMapProperty) {
			((StochasticMapProperty)f).sampleMutations(sampleNr);
		}
		
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
