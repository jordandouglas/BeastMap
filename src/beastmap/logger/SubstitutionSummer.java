package beastmap.logger;

import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;


@Description("Sums an entire vector of BranchSubstLogger terms")
public class SubstitutionSummer extends CalculationNode implements Loggable {

	
	final public Input<BranchSubstLogger> counterInput = new Input<>("counter", "counter to log", Input.Validate.REQUIRED);
	
	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}


	@Override
	public void init(PrintStream out) {
		out.print("sum." + counterInput.get().getID() + "\t");
	}

	@Override
	public void log(long sample, PrintStream out) {
		
		counterInput.get().sampleMutations(sample);
		
		double sum = 0;
		for (int i = 0; i < counterInput.get().getDimension(); i ++) {
			sum += counterInput.get().getArrayValue(i);
		}
		
		out.print(sum + "\t");
	}

	@Override
	public void close(PrintStream out) {
		
	}

}
