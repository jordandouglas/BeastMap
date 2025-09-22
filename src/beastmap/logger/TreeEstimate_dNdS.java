package beastmap.logger;

import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;
import codonmodels.evolution.datatype.Codon;
import codonmodels.evolution.datatype.GeneticCode;


@Description("Estimate of dN/dS for the whole tree")
public class TreeEstimate_dNdS extends CalculationNode implements Loggable {
	
	final public Input<StochasticMapProperty> nNInput = new Input<>("nN", "number of non-synonymous substituions", Input.Validate.REQUIRED);
	final public Input<StochasticMapProperty> nSInput = new Input<>("nS", "number of synonymous substituions", Input.Validate.REQUIRED);
	
	
	final double pSpN = 1.0 * 138/438; // Number of possible synonymous vs non-synonymous changes
	int frame;
	GeneticCode code;
	Codon codon;
	
	@Override
    public void initAndValidate() {
		
		
	
		
		
    }
	
	




	@Override
	public void init(PrintStream out) {
		out.print(this.getID() + "\t");
	}


	@Override
	public void log(long sample, PrintStream out) {
		nNInput.get().sampleMutations(sample);
		nSInput.get().sampleMutations(sample);
		
		double nN = (Double)nNInput.get().getPropertyOfNode(null);
		double nS = (Double)nSInput.get().getPropertyOfNode(null);
		
		
		double val = 0;
		if (nN == 0) val = 0;
		else if (nS == 0) val = Double.POSITIVE_INFINITY; 
		else val = nN/nS*pSpN;
		
		out.print(val + "\t");
		
	}


	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}



	
}


