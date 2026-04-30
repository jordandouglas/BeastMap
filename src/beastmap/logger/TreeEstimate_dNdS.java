package beastmap.logger;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;
import beastmap.evolution.StochasticMapper;
import beastmap.util.Mutation;
import beastmap.util.MutationUtils;
import codonmodels.evolution.datatype.Codon;
import codonmodels.evolution.datatype.GeneticCode;


@Description("Estimate of dN/dS for the whole tree")
public class TreeEstimate_dNdS extends CalculationNode implements Loggable, StochasticMapProperty {
	
	final public Input<StochasticMapper> mapperInput = new Input<>("mapper", "stochastic mapper", Input.Validate.REQUIRED);
	final public Input<String> geneticCodeInput = new Input<>("code", "name of genetic code", "universal", GeneticCode.GENETIC_CODE_NAMES);
	final public Input<Boolean> perSiteInput = new Input<>("eachSite", "log one value for each site?", false);
	final public Input<String> filterInput = new Input<>("filter", "specifies which of the sites in the input alignment we will restrict to" +
            "First site is 1." +
            "Filter specs are comma separated, either a singleton, a range [from]-[to] or iteration [from]:[to]:[step]; " +
            "1-100 defines a range, " +
            "1-100\3 or 1:100:3 defines every third in range 1-100, " +
            "1::3,2::3 removes every third site. " +
            "Default for range [1]-[last site], default for iterator [1]:[last site]:[1]", Validate.OPTIONAL);
	
	
	final public Input<Double> vnInput = new Input<>("vN", "Shape parameter for non-synonymous prior", 1.0);
	final public Input<Double> vsInput = new Input<>("vS", "Shape parameter for synonymous prior", 2.0);
	
	
	
	protected List<Integer> filter;
	final double pSpN = 1.0 * 138/438; // Number of possible synonymous vs non-synonymous changes
	int frame;
	GeneticCode code;
	Codon codon;
	
	@Override
    public void initAndValidate() {
		
		// Parse inputs
		this.code = GeneticCode.findByName(geneticCodeInput.get());
		if (this.code == null) {
			throw new IllegalArgumentException("Cannot find code " + geneticCodeInput.get());
		}
		this.codon = new Codon(code);
		
		
		this.filter = new ArrayList<>();
		
		
		// Log all sites?
		int nsites = mapperInput.get().getPatternCount();
		int ncodons = nsites / 3;
		if (perSiteInput.get()) {
			for (int i = 0; i < ncodons; i ++) {
				this.filter.add(i);
			}
		} 
		
		
		// Log just some sites?
		else if (filterInput.get() != null && !filterInput.get().isEmpty()) {
			
    		this.filter = MutationUtils.parseFilterSpec(ncodons, filterInput.get());
    	}
		
		
		
		
    }
	
	




	@Override
	public void init(PrintStream out) {
		
		
		String id = (this.getID() == null ? "" : "."+this.getID());
		if (!this.filter.isEmpty()) {
			for (int i : filter) {
				int s = i + 1;
				out.print("dNdS.JC" + id + "." + s + "\t" + "dNdS.HKY" + id + "." + s + "\t");
			}
				
				
		}else {
			out.print("dNdS.JC" + id + "\t" + "dNdS.HKY" + id + "\t");
		}
		
		
	}


	@Override
	public void log(long sample, PrintStream out) {
		
		
		StochasticMapper mapper = mapperInput.get();
		mapper.sampleMutations(sample);
		
		if (!this.filter.isEmpty()) {
			
			
			
			double[] freqs = MutationUtils.getCodonFreqs(mapper, codon, 0.1, -1);
			for (int codonNr : filter) {
				double[] vals = MutationUtils.getdNdS(codon, code, mapper, codonNr, vnInput.get(), vsInput.get(), freqs);
				out.print(vals[0] + "\t" + vals[1] + "\t");
			}
				
				
		}else {
			
			double[] vals = MutationUtils.getdNdS(codon, code, mapper, vnInput.get(), vsInput.get());
			out.print(vals[0] + "\t" + vals[1] + "\t");
			
		}
		
		
		//System.out.println("XXX " + sample);
		
		
		
		
	}


	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}



	@Override
	public String getName() {
		return "dNdS";
	}



	@Override
	public void sampleMutations(long sampleNr) {
		StochasticMapper mapper = mapperInput.get();
		mapper.sampleMutations(sampleNr);
	}



	@Override
	public Object getPropertyOfNode(Node node) {
		// TODO Auto-generated method stub
		return null;
	}



	
}


