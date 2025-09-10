package beastmap.util;

import beast.base.evolution.tree.Node;

/**
 * A mutation along a branch
 */
public class Mutation implements Comparable<Mutation> {
	
	int from;
	int to;
	double time;
	int siteNr;
	int branchParent;
	int branchChild;
	Node node;
	
	
	public Mutation(Mutation other) {
		this.from = other.from;
		this.to = other.to;
		this.time = other.time;
		this.siteNr = other.siteNr;
		this.branchParent = other.branchParent;
		this.branchChild = other.branchChild;
		this.node = other.node;
	}
	
	public Mutation(int from, int to, double time, int siteNr, int branchParent, int branchChild, Node node) {
		this.from = from;
		this.to = to;
		this.time = time;
		this.siteNr = siteNr;
		this.branchParent = branchParent;
		this.branchChild = branchChild;
		this.node = node;
	}

	public int getSiteNr() {
		return siteNr;
	}
	
	public Node getNode() {
		return node;
	}
	
	 @Override
     public int compareTo(Mutation other) {
         return Double.compare(this.time, other.time);
     }

	 public int getFrom() {
		return from;
	 }
	
	 public int getTo() {
		return to;
	 }
	 
	 public void setSiteNr(int siteNr) {
		 this.siteNr = siteNr;
	 }

}


