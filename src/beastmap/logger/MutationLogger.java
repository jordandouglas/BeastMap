package beastmap.logger;


import java.util.List;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;
import beastmap.evolution.StochasticMapper;
import beastmap.util.Mutation;

@Description("Logs the time and nature of mutations along each branch")
public class MutationLogger extends CalculationNode implements StochasticMapProperty, Function {
	
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
		
		if (node.isRoot()) return "";
		
		List<Mutation> mutations = sampler.getMutationsOnBranch(node.getNr());
		String result = "";
		
		for (int i = 0; i < mutations.size(); i++) {
			Mutation mut = mutations.get(i);
			result += mut.getString(sampler.getDataTypeOfMapper());
			if (i < mutations.size()-1) {
				result += ",";
			}
		}
		
		
		result += "";
		return result;
		
	}


	@Override
	public String getName() {
		return "mut" + "." + this.getID();
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


