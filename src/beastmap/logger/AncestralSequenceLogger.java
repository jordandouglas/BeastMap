package beastmap.logger;


import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;
import beastmap.evolution.StochasticMapper;

@Description("Logs the sequence at each internal node")
public class AncestralSequenceLogger extends CalculationNode implements StochasticMapProperty, Function {
	
	final public Input<StochasticMapper> samplerInput = new Input<>("sampler", "mutation sampler to log", Input.Validate.OPTIONAL);
	final public Input<StochasticMapper> truthInput = new Input<>("truth", "mutation sampler to log", Validate.XOR, samplerInput);
	
	int [] siteStates;
	
	StochasticMapper sampler;
	StochasticMapper truth;
	
	@Override
	public void initAndValidate() {
		this.sampler = samplerInput.get();
		this.truth = truthInput.get();
	}
	
	
	@Override
	public void sampleMutations(long sampleNr) {
		if (this.sampler != null) {
			this.sampler.sampleMutations(sampleNr);
		}
		
	}


	@Override
	public Object getPropertyOfNode(Node node) {
		
		
		StochasticMapper mapper = this.sampler != null ? this.sampler : this.truth;
		
		DataType dt = mapper.getDataTypeOfMapper();
		int[] seqInt = null;
		try {
			seqInt = mapper.getStatesForNode(mapper.getTree(), node);
			
		}catch (Exception e) {
			seqInt = new int[mapper.getPatternCount()];
			for (int i = 0; i < mapper.getPatternCount();  i++) {
				seqInt[i] = 0;
			}
		}
		
		
		
		return dt.encodingToString(seqInt);
	}


	@Override
	public String getName() {
		return "sequence" + "." + this.getID();
	}


	@Override
	public int getDimension() {
		if (this.sampler != null) {
			return sampler.getTree().getNodeCount();
		}else {
			return truth.getTree().getNodeCount();
		}
	}


	@Override
	public double getArrayValue(int dim) {
		return 0; // Not implemented -- use getPropertyOfNode so we can return a String
	}	
	






}
