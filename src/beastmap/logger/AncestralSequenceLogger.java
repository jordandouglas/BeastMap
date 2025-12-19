package beastmap.logger;


import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;
import beastmap.evolution.StochasticMapper;

@Description("Logs the sequence at each internal node")
public class AncestralSequenceLogger extends CalculationNode implements StochasticMapProperty, Function {
	
	final public Input<StochasticMapper> samplerInput = new Input<>("sampler", "mutation sampler to log", Input.Validate.REQUIRED);
	
	int [] siteStates;
	
	StochasticMapper sampler;
	
	@Override
	public void initAndValidate() {
		this.sampler = samplerInput.get();
	}
	
	
	@Override
	public void sampleMutations(long sampleNr) {
		this.sampler.sampleMutations(sampleNr);
	}


	@Override
	public Object getPropertyOfNode(Node node) {
		
		DataType dt = sampler.getDataTypeOfMapper();
		int[] seqInt = null;
		try {
			seqInt = sampler.getStatesForNode(sampler.getTree(), node);
			
		}catch (Exception e) {
			seqInt = new int[sampler.getPatternCount()];
			for (int i = 0; i < sampler.getPatternCount();  i++) {
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
		return sampler.getTree().getNodeCount();
	}


	@Override
	public double getArrayValue(int dim) {
		return 0; // Not implemented -- use getPropertyOfNode so we can return a String
	}	
	






}
