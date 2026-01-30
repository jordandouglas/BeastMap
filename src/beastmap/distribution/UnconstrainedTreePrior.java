package beastmap.distribution;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BetaDistributionImpl;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import beastmap.evolution.TimelessTree;


@Description("A tree prior for a non-time tree")
@Citation(value="Douglas, J., & Bromham, L. (2025). Reconstructing substitution histories on phylogenies, with accuracy, precision, and coverage. bioRxiv, 2025-12.", DOI="10.64898/2025.12.21.695861")
public class UnconstrainedTreePrior extends SpeciesTreeDistribution {

	final public Input<RealParameter> meanInput = new Input<>("mean", "mean branch length under a gamma distribution prior", Input.Validate.REQUIRED);
	final public Input<RealParameter> shapeInput = new Input<>("shape", "gamma distribution prior shape of branch lengths. Set to 1 for exponential distribution.", Input.Validate.REQUIRED);
	
	
	final public Input<RealParameter> midpointInput = new Input<>("midpoint", "the relative position of the root between the two furtherest tips is a "
			+ "beta(alpha, alpha). When alpha is large (e.g. 50), this approaches the midpoint. When alpha is small (e.g., 1), all root placements "
			+ "are equally likely along this path.", Input.Validate.OPTIONAL);
	
	
	org.apache.commons.math.distribution.GammaDistribution gammaDist = new GammaDistributionImpl(1, 1);
	org.apache.commons.math.distribution.BetaDistribution betaDist = new BetaDistributionImpl(1, 1);
	
	public void initAndValidate() {
		
		if (! (treeInput.get() instanceof TimelessTree)) {
			Log.warning("Warning: assuming that " + treeInput.get().getID() + " is an unconstrained substitution tree.");
		}
		
		super.initAndValidate();
	}
	
	@Override
	public double calculateTreeLogLikelihood(TreeInterface tree) {
		
		double mean = meanInput.get().getArrayValue();
		double shape = shapeInput.get().getArrayValue();
		double scale = mean/shape;
		this.gammaDist = new GammaDistributionImpl(shape, scale);
		
		double logPTree = 0;
		for (Node node : treeInput.get().getNodesAsArray()) {
			if (node.isRoot()) continue;
			double length = node.getLength();
			
			// The two branches coming out of the root will count as a single branch
			if (midpointInput.get() != null && node.getParent().isRoot()) {
				
				
				// Skip one branch
				if (node == node.getParent().getRight()){
					continue;
				}else {
					Node left = node;
					Node right = node.getParent().getRight();
					length = left.getLength() + right.getLength();
				}
				
			}
			
			if (length <= 0) {
				return Double.NEGATIVE_INFINITY;
			}
			
			logPTree += this.gammaDist.logDensity(length);
			
			if (node.getHeight() < 0) {
				Log.warning("height=" + node.getHeight());
			}
			
		}
		
		if (midpointInput.get() != null) {
			
			
			Node root = treeInput.get().getRoot();
			
			// Find the furtherest tip on either side
			Node left = root.getLeft();
			Node right = root.getRight();
			
			double leftHeight = getYoungestLeafHeight(left);
			double rightHeight = getYoungestLeafHeight(right);
			
			
			// Difference between heights
			double rootToLeft = root.getHeight() - leftHeight;
			double rootToRight = root.getHeight() - rightHeight;
			double alpha = midpointInput.get().getArrayValue(); 
			this.betaDist = new BetaDistributionImpl(alpha, alpha);
			double proportion = rootToLeft / (rootToLeft + rootToRight); 
			
			//Log.warning("diff" + diff + " midpointRate " + midpointRate + " rootToLeft " + rootToLeft + " rootToRight " + rootToRight);
			
			logPTree += this.betaDist.logDensity(proportion);
			
			
		}
		
		return logPTree;
	}
	
	
	private double getYoungestLeafHeight(Node node) {
		
		if (node.isLeaf()) return node.getHeight();
		
		double youngest = Double.POSITIVE_INFINITY;
		for (Node leaf : node.getAllLeafNodes()) {
			youngest = Math.min(youngest, leaf.getHeight());
		}

		return youngest;
	}
	
    @Override
    public List<String> getConditions() {
        List<String> conditions = new ArrayList<>();
        conditions.add(meanInput.get().getID());
        conditions.add(shapeInput.get().getID());
        if (midpointInput.get() != null) conditions.add(midpointInput.get().getID());
        return conditions;
    }

    @Override
    public List<String> getArguments() {
        List<String> arguments = new ArrayList<>();
        arguments.add(treeInput.get().getID());
        return arguments;
    }
    
    @Override
    protected boolean requiresRecalculation() {
    	if (InputUtil.isDirty(meanInput)) return true;
    	if (InputUtil.isDirty(midpointInput)) return true;
    	if (InputUtil.isDirty(shapeInput)) return true;
    	
        return treeInput.get().somethingIsDirty();
    }
    
    
    @Override
    public void sample(State state, Random random) {
    	
    	 if (sampledFlag) return;

         sampledFlag = true;

         // Cause conditional parameters to be sampled
         sampleConditions(state, random);
         
         
         
         double mean = meanInput.get().getArrayValue();
         double shape = shapeInput.get().getArrayValue();
         double scale = mean/shape;
         this.gammaDist = new GammaDistributionImpl(shape, scale);
         Tree tree = (Tree) treeInput.get();
         
         

         // Simulate tree conditional on new parameters
         List<Node> activeLineages = new ArrayList<>();
         for (Node oldLeaf : tree.getExternalNodes()) {
             Node newLeaf = new Node(oldLeaf.getID());
             newLeaf.setNr(oldLeaf.getNr());
             newLeaf.setHeight(0.0);
             activeLineages.add(newLeaf);
         }

         int nextNr = activeLineages.size();

         while (activeLineages.size() > 1) {
			int k = activeLineages.size();
			
			// Sample 2 branch lengths
			double branchLength1=0, branchLength2=0;
			try {
				branchLength1 = this.gammaDist.inverseCumulativeProbability(Randomizer.nextFloat());
				branchLength2 = this.gammaDist.inverseCumulativeProbability(Randomizer.nextFloat());
			} catch (MathException e) {
				e.printStackTrace();
			}
			
			
			// To make the root, we will sample just a single branch length and then split it in half
			if (activeLineages.size() == 2) {
				double sum = branchLength1;
				branchLength1 = sum/2;
				branchLength2 = sum/2;
			}
             
             // Sample 2 nodes
             Node node1 = activeLineages.get(random.nextInt(k));
             Node node2;
             do {
                 node2 = activeLineages.get(random.nextInt(k));
             } while (node2.equals(node1));

             // Join them
             Node newParent = new Node();
             newParent.setNr(nextNr++);
             newParent.addChild(node1);
             newParent.addChild(node2);
             
             // Set parent height to 0 and children heights to negative lengths
             newParent.setHeight(0);
             
             // Push all the left nodes down
             List<Node> childNodes = new ArrayList<>();
             node1.getAllChildNodesAndSelf(childNodes);
             for (Node node : childNodes) {
            	 node.setHeight(node.getHeight() - branchLength1);
             }
             
             // Push all the right nodes down
             childNodes = new ArrayList<>();
             node2.getAllChildNodesAndSelf(childNodes);
             for (Node node : childNodes) {
            	 //Log.warning("setting " + node.isLeaf() + " height to " + node.getHeight() + "-" + branchLength2 + " " + node.toNewick());
            	 node.setHeight(node.getHeight() - branchLength2);
            	 
             }
             
             
             // Update list
             activeLineages.remove(node1);
             activeLineages.remove(node2);
             activeLineages.add(newParent);
             
         }
         
         
        


		tree.assignFromWithoutID(new Tree(activeLineages.get(0)));
		
		// Reset tree height so that the youngest leaf has height 0
		setNodeHeights(tree);
		
		
		if (midpointInput.get() != null) {
			
			double alpha = midpointInput.get().getValue();
			this.betaDist = new BetaDistributionImpl(alpha, alpha);
			double rootProportion;
			try {
				rootProportion = betaDist.inverseCumulativeProbability(Randomizer.nextDouble());
				//Log.warning("alpha = " + alpha + " prop = " + rootProportion);
			} catch (MathException e) {
				e.printStackTrace();
				rootProportion = 0.5;
			}
			
			// Convert to an unrooted tree
			RootnessNode unrootedNode = convert(tree.getRoot());
			MidpointResult result = findMidpoint(unrootedNode, rootProportion);
			Tree tree2 = convert(result.nodeA, result.nodeB, result.distFromAToRoot);
			 
			//Log.warning("before rooting " + tree.getRoot().toNewick());
			//Log.warning(" after rooting" + tree2.getRoot().toNewick());
			 
			//tree.assignFromWithoutID(new Tree(tree2.getRoot()));
			tree.assignFromWithoutID(tree2);
		
			 
		}
		
	
         
        
         
    }
    
  
    
    
    
    // Reset tree height so that the youngest leaf has height 0
    private static void setNodeHeights(Tree tree) {
    	setNodeHeights(tree.getRoot());
    }
    
    
    // Reset tree height so that the youngest leaf has height 0
    private static void setNodeHeights(Node root) {
    	if (root.isLeaf()) {
    		root.setHeight(0);
    		return;
    	}
    	
    	double minHeight = Double.POSITIVE_INFINITY;
		for (Node node : root.getAllChildNodesAndSelf()) {
			minHeight = Math.min(minHeight, node.getHeight());
		}
		for (Node node : root.getAllChildNodesAndSelf()) {
			node.setHeight(node.getHeight() - minHeight);
			//Log.warning("Setting node height to " + node.getHeight());
		}
    }
    
    
    
    /**
     * Convert a BEAST2 rooted Node (with heights) into an unrooted RootnessNode structure.
     * @param beastRoot The BEAST2 root node
     * @return The corresponding RootnessNode (unrooted representation)
     */
    private RootnessNode convert(Node beastRoot) {
        Map<Node, RootnessNode> mapping = new HashMap<>();
        
        RootnessNode left = convertRecursive(beastRoot.getLeft(), mapping);
        RootnessNode right = convertRecursive(beastRoot.getRight(), mapping);
        left.addNeighbour(right, beastRoot.getLeft().getLength() + beastRoot.getRight().getLength());
        return left;
    }

    private RootnessNode convertRecursive(Node beastNode, Map<Node, RootnessNode> mapping) {
    	
        // If we already built this node, return it
        if (mapping.containsKey(beastNode)) {
            return mapping.get(beastNode);
        }

        // Make new RootnessNode (label only if leaf)
        String label = beastNode.isLeaf() ? beastNode.getID() : null;
        RootnessNode thisNode = new RootnessNode(label);
        mapping.put(beastNode, thisNode);

        // Process children
        for (int i = 0; i < beastNode.getChildCount(); i++) {
            Node child = beastNode.getChild(i);
            RootnessNode childNode = convertRecursive(child, mapping);

            // branch length = parent.height - child.height
            double length = beastNode.getHeight() - child.getHeight();

            thisNode.addNeighbour(childNode, length);
        }

        return thisNode;
    }
    
    
    public static Tree convert(RootnessNode nodeA, RootnessNode nodeB, double rootHeightFromALeaf) {
       
    
        
        double connectingLength = nodeA.getBranchLengthTo(nodeB);
        double heightOfA = nodeA.lengthToFurthestLeaf(nodeB);
        double heightOfB = nodeB.lengthToFurthestLeaf(nodeA);
        double spanDistance = heightOfA + heightOfB + connectingLength;

        
        double lenToA = rootHeightFromALeaf - heightOfA;
        double lenToB = connectingLength - lenToA;
        
        
    	// Create new root
        Node root = new Node();
        root.setHeight(rootHeightFromALeaf);
        
        //Log.warning("rootHeightFromALeaf " + rootHeightFromALeaf +  ", connectingLength " + connectingLength + ", heightOfA " + heightOfA + ", heightOfB " + heightOfB + " , lenToA" + lenToA + " , lenToB" + lenToB);

        // Create two subtrees, one descending from nodeA and one from nodeB
        Node left  = buildSubtree(nodeA, rootHeightFromALeaf, lenToA, nodeB, new HashMap<>(), new HashSet<>());
        Node right = buildSubtree(nodeB, rootHeightFromALeaf, lenToB, nodeA,  new HashMap<>(), new HashSet<>());
        
        root.addChild(left);
        root.addChild(right);

        
        
        
        // Renumber nodes
        int nr = 0;
        for (Node leaf : root.getAllChildNodesAndSelf()) {
        	if (!leaf.isLeaf()) continue;
        	leaf.setNr(nr);
        	nr++;
        }
        for (Node internal : root.getAllChildNodesAndSelf()) {
        	if (internal.isRoot() || internal.isLeaf()) continue;
        	internal.setNr(nr);
        	nr++;
        }
        root.setNr(nr);
        
        Tree tree = new Tree(root);
        
        return tree;
    }

    /**
     * Recursively build a subtree below a given neighbour.
     *
     * @param current The current LengthedNode being converted
     * @param parentHeight The height of the BEAST parent node
     * @param cameFrom The neighbour we came from (so we don't backtrack)
     * @param mapping Cache from LengthedNode -> Node
     * @param visited Set of visited LengthedNodes
     */
    private static Node buildSubtree(RootnessNode current, double parentHeight, double branchLen, RootnessNode cameFrom, Map<RootnessNode, Node> mapping, Set<RootnessNode> visited) {
        if (mapping.containsKey(current)) {
            return mapping.get(current);
        }
        visited.add(current);

        Node thisNode = new Node();
        if (current.getLabel() != null) {
            thisNode.setID(current.getLabel());
        }

        // Height = parentHeight - branchLength
        double height = parentHeight - branchLen;
        thisNode.setHeight(height);
        mapping.put(current, thisNode);

        // Recurse into neighbours (except the one we came from)
        for (RootnessNode neigh : current.getNeighbours()) {
            if (neigh == cameFrom) continue;
            Node child = buildSubtree(neigh, height, current.getBranchLengthTo(neigh), current, mapping, visited);
            thisNode.addChild(child);
        }

        return thisNode;
    }
    
    
    // Tree nodes where they store their neighbours and the length between them, but there is no direction
    public class RootnessNode {

        private final List<RootnessNode> neighbours = new ArrayList<>();
        private final List<Double> branchLengths = new ArrayList<>();

        private final String label; // optional, for leaf names

        public RootnessNode(String label) {
            this.label = label;
        }

        
        // Length to furthest leaf on this side of the other node
        public double lengthToFurthestLeaf(RootnessNode otherNode) {
        	Set<RootnessNode> visited = new HashSet<>();
			return lengthToFurthestLeaf(otherNode, visited);
		}
        
        private double lengthToFurthestLeaf(RootnessNode cameFrom, Set<RootnessNode> visited) {
        	
        	visited.add(this);
        	
        	double maxDistFromLeaves = 0;
        	for (RootnessNode neigh : this.getNeighbours()) {
                if (neigh == cameFrom) continue;
                double d1 = neigh.lengthToFurthestLeaf(this, visited);
                double d2 = neigh.getBranchLengthTo(this);
                double d = d1 + d2;
                //Log.warning("dist = " + d1 + " + " + d2 + "=" + d);
                if (d > maxDistFromLeaves) {
                	maxDistFromLeaves = d;
                }
            }
        	
			return maxDistFromLeaves;
		}

		public String getLabel() {
            return label;
        }

        /** Add a connection between RootnessNode node and another, with a branch length */
        public void addNeighbour(RootnessNode other, double length) {
            if (neighbours.size() >= 3) {
                throw new IllegalStateException("A RootnessNode can have at most 3 neighbours.");
            }
            neighbours.add(other);
            branchLengths.add(length);

            // Add reciprocal link if not already there
            if (!other.neighbours.contains(this)) {
                other.addNeighbour(this, length);
            }
        }

        /** Number of neighbours (1 = leaf, 3 = internal node) */
        public int degree() {
            return neighbours.size();
        }
        

        /** Return an unmodifiable view of neighbours */
        public List<RootnessNode> getNeighbours() {
            return Collections.unmodifiableList(neighbours);
        }

        /** Return branch length to a specific neighbour */
        public double getBranchLengthTo(RootnessNode neighbour) {
            int idx = neighbours.indexOf(neighbour);
            if (idx < 0) {
                throw new IllegalArgumentException("Not a neighbour of this node");
            }
            return branchLengths.get(idx);
        }

        /** Set branch length to a specific neighbour */
        public void setBranchLengthTo(RootnessNode neighbour, double length) {
            int idx = neighbours.indexOf(neighbour);
            if (idx < 0) {
                throw new IllegalArgumentException("Not a neighbour of this node");
            }
            branchLengths.set(idx, length);

            // also update reciprocal link
            int idxOther = neighbour.neighbours.indexOf(this);
            if (idxOther >= 0) {
                neighbour.branchLengths.set(idxOther, length);
            }
        }

        @Override
        public String toString() {
            return label != null ? label : super.toString();
        }
        
    }
    
    
    public static MidpointResult findMidpoint(RootnessNode anyNode, double midpointProportion) {
    	
        // Step 1: find one farthest leaf
    	RootnessNode leaf1 = farthestLeaf(anyNode).node;

        // Step 2: from leaf1, find farthest leaf (this is leaf2)
        FarthestResult res = farthestLeaf(leaf1);
        RootnessNode leaf2 = res.node;
        double diameter = res.distance;

        // Step 3: get path leaf1 -> leaf2
        List<RootnessNode> path = getPath(leaf1, leaf2);

        // Step 4: walk along the path
        double halfway = diameter * midpointProportion;
        double distSoFar = 0.0;
        
       // Log.warning("midpointProportion " + midpointProportion);

        for (int i = 0; i < path.size() - 1; i++) {
        	RootnessNode a = path.get(i);
        	RootnessNode b = path.get(i + 1);
            double len = a.getBranchLengthTo(b);

            if (distSoFar + len >= halfway) {
                double distFromA = halfway - distSoFar;
                return new MidpointResult(a, b, distFromA, halfway);
            }
            distSoFar += len;
        }

        throw new IllegalStateException("Midpoint not found (should not happen)");
    }

    /** Helper: get farthest leaf from start using BFS/DFS */
    private static FarthestResult farthestLeaf(RootnessNode start) {
        Set<RootnessNode> visited = new HashSet<>();
        FarthestResult result = new FarthestResult(start, 0.0);

        dfsFarthest(start, null, 0.0, visited, result);
        return result;
    }

    private static void dfsFarthest(RootnessNode node, RootnessNode parent,
                                    double dist, Set<RootnessNode> visited,
                                    FarthestResult result) {
        visited.add(node);
        if (node.degree() == 1 && dist > result.distance) {
            result.node = node;
            result.distance = dist;
        }
        for (RootnessNode neigh : node.getNeighbours()) {
            if (neigh == parent) continue;
            dfsFarthest(neigh, node, dist + node.getBranchLengthTo(neigh), visited, result);
        }
    }

    /** Helper: get path between two nodes */
    private static List<RootnessNode> getPath(RootnessNode start, RootnessNode target) {
        Map<RootnessNode, RootnessNode> parentMap = new HashMap<>();
        Deque<RootnessNode> stack = new ArrayDeque<>();
        Set<RootnessNode> visited = new HashSet<>();

        stack.push(start);
        visited.add(start);

        // DFS until we reach target
        while (!stack.isEmpty()) {
        	RootnessNode node = stack.pop();
            if (node == target) break;
            for (RootnessNode neigh : node.getNeighbours()) {
                if (!visited.contains(neigh)) {
                    visited.add(neigh);
                    parentMap.put(neigh, node);
                    stack.push(neigh);
                }
            }
        }

        // Reconstruct path
        List<RootnessNode> path = new ArrayList<>();
        for (RootnessNode n = target; n != null; n = parentMap.get(n)) {
            path.add(n);
            if (n == start) break;
        }
        Collections.reverse(path);
        return path;
    }

    /** Result container */
    public static class MidpointResult {
        public final RootnessNode nodeA, nodeB;
        public final double distFromA; // distance from nodeA along branch
        public final double distFromAToRoot;

        public MidpointResult(RootnessNode a, RootnessNode b, double distFromA, double distFromAToRoot) {
            this.nodeA = a;
            this.nodeB = b;
            this.distFromA = distFromA;
            this.distFromAToRoot = distFromAToRoot;
        }
    }

    private static class FarthestResult {
    	RootnessNode node;
        double distance;
        FarthestResult(RootnessNode n, double d) { node = n; distance = d; }
    }
    
    

}



