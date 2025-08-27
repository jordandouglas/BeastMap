package mutationtree.distribution;

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


@Description("A tree prior for a non-time tree")
public class UnconstrainedTreePrior extends SpeciesTreeDistribution {

	final public Input<RealParameter> rateInput = new Input<>("rate", "exponential distribution prior rate of branch lengths", Input.Validate.REQUIRED);
	final public Input<RealParameter> midpointInput = new Input<>("midpoint", "the difference between the root's two furtherest tips is an exponential with this mean", Input.Validate.OPTIONAL);
	
	
	@Override
	public double calculateTreeLogLikelihood(TreeInterface tree) {
		
		double rate = rateInput.get().getArrayValue();
		
		double logPTree = 0;
		for (Node node : treeInput.get().getNodesAsArray()) {
			if (node.isRoot()) continue;
			double length = node.getLength();
			
//			// So that this process is generative, the two branches coming out of the root will count as a single branch
			// Edit: since we are using an exponential prior, this step is redundant and gives the same answer
//			if (midpointInput.get() != null && node.getParent().isRoot()) {
//				
//				
//				// Skip one branch
//				if (node == node.getParent().getRight()){
//					continue;
//				}else {
//					Node left = node;
//					Node right = node.getParent().getRight();
//					length = left.getHeight() + right.getHeight();
//				}
//				
//			}
			
			if (length <= 0) {
				return Double.NEGATIVE_INFINITY;
			}
			
			logPTree += Math.log(rate) - rate*length; // Exponential distribution density in log space
			
			if (node.getHeight() < 0) {
				Log.warning("height=" + node.getHeight());
			}
			
		}
		
		if (midpointInput.get() != null) {
			Node root = treeInput.get().getRoot();
			
			// Find the furtherest tip on either side
			Node left = root.getLeft();
			Node right = root.getRight();
			
			double leftHeight = Double.POSITIVE_INFINITY;
			double rightHeight = Double.POSITIVE_INFINITY;
			
			for (Node node : left.getAllLeafNodes()) {
				leftHeight = Math.min(leftHeight, node.getHeight());
			}
			for (Node node : right.getAllLeafNodes()) {
				rightHeight = Math.min(rightHeight, node.getHeight());
			}
			
			
			// Difference between heights
			double rootToLeft = root.getHeight() - leftHeight;
			double rootToRight = root.getHeight() - rightHeight;
			double diff = Math.abs(rootToLeft - rootToRight);
			double midpointRate = 1 / midpointInput.get().getArrayValue();
			logPTree += Math.log(midpointRate) - midpointRate*diff; // Exponential distribution density in log space
			
			
		}
		
		return logPTree;
	}
	
    @Override
    public List<String> getConditions() {
        List<String> conditions = new ArrayList<>();
        conditions.add(rateInput.get().getID());
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
    	if (InputUtil.isDirty(rateInput)) return true;
    	if (InputUtil.isDirty(midpointInput)) return true;
        return treeInput.get().somethingIsDirty();
    }
    
    
    @Override
    public void sample(State state, Random random) {
    	
    	 if (sampledFlag) return;

         sampledFlag = true;

         // Cause conditional parameters to be sampled
         sampleConditions(state, random);
         
         
         
         double lambda = rateInput.get().getValue();
         double midpointMean = midpointInput.get().getValue();
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
             double branchLength1 = Randomizer.nextExponential(lambda);
             double branchLength2 = Randomizer.nextExponential(lambda);

             
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
//		
		
		//Log.warning("before " + tree.getRoot().toNewick());
		
//		
//		// Midpoint rooting: find the two leaves that are the furtherest apart
//		int ntaxa = tree.getLeafNodeCount();
//		double[][] dists = new double[ntaxa][ntaxa];
//		for (int i = 0; i < ntaxa; i ++) {
//			
//			Node leaf1 = tree.getNode(i);
//			for (int j = i+1; j < ntaxa; j ++) {
//				
//				// The root is at height 0 and all leaves are negative
//				Node leaf2 = tree.getNode(j);
//				double distance = -(leaf1.getHeight() + leaf2.getHeight());
//				dists[i][j] = distance;
//				
//			}
//		}
//		
//		// Max distance pair
//		double maxD = 0;
//		int younger = 0;
//		for (int i = 0; i < ntaxa; i ++) {
//			for (int j = i+1; j < ntaxa; j ++) {
//				if (dists[i][j] > maxD) {
//					maxD = dists[i][j];
//					
//					Node leaf1 = tree.getNode(i);
//					Node leaf2 = tree.getNode(j);
//					if (leaf1.getHeight() < leaf2.getHeight()) {
//						younger = i;
//					}else {
//						younger = j;
//					}
//				}
//			}
//		}
//		
//		
//		// Put the root in the middle of these two leaves, plus a Laplace distributed amount
//		double w = Randomizer.nextExponential(1/midpointMean);
//		if (Randomizer.nextBoolean()) w = -w;
//		double distAboveYounger = maxD/2 + w;
//		
//		
//		
		
		
		
		// NOT WORKING
		
		 // Convert to an unrooted tree
		// RootnessNode unrootedNode = convert(tree.getRoot());
		//// MidpointResult result = findMidpoint(unrootedNode);
		// Tree tree2 = convert(result.nodeA, result.nodeB, maxD/2);
		 
//		 Log.warning("before " + tree.getRoot().toNewick());
//		 Log.warning("after " + tree2.getRoot().toNewick());
//		 
		 
		// tree.assignFromWithoutID(new Tree(tree2.getRoot()));
		
		
		// Find the node whose parent will become the root
		//double t = 0;
		//Node node = tree.getNode(younger);
		
		// Go up the tree from younger towards the root until we are half way to the other node.
		// Because 'younger' is the younger of the two nodes, we will stop before we reach the tentative root
		// TODO ACCOUNT FOR RANDOMNESS, IN SOME CASE WE WANT THE OLDER ONE
//		while (true) {
//			
//			double length = node.getLength();
//			if (t + length >= distAboveYounger) {
//				
//				// We will put the root above this node
//				
//				
//				// Free the current root node and keep track of which of its children will become the parent and child later 
//				Node root = tree.getRoot();
//				Node oldRootElder, oldRootYounger;
//				boolean goingUpLeft = root.getLeft().getAllChildNodesAndSelf().contains(node);
//				if (goingUpLeft) {
//					oldRootElder = root.getLeft();
//					oldRootYounger = root.getRight();
//				}else {
//					oldRootYounger = root.getLeft();
//					oldRootElder = root.getRight();
//				}
//				root.removeChild(oldRootYounger);
//				root.removeChild(oldRootElder);
//				
//				
//				// Move the root above our current node
//				Node parent = node.getParent();
//				Node grandParent = parent.getParent(); // Will be null if parent is the root
//				if (grandParent != null) {
//					grandParent.removeChild(parent); // Redundant if grandparent is the root
//				}
//				Node sister = parent.getLeft() == node ? parent.getRight() : parent.getLeft();
//				
//				root.setHeight(distAboveYounger);
//				root.addChild(node);
//				root.addChild(parent);
//				
//				// Adjust the whole clade under parent by a fixed distance
//				double dParentClade =  2 * (t + length - distAboveYounger);
//				for (Node n : parent.getAllChildNodesAndSelf()) {
//					n.setHeight(n.getHeight() - dParentClade);
//				}
//				
//				
//				// Now we need to reconnect the nodes around the old root
//				if (grandParent != null) {
//					
//					
//					
//					
//				}else {
//					
//					
//					
//				}
//				
//				
//				break;
//			}
//			
//			t = t + length;
//			
//			
//		}
//		
		
		
		
	
         
        
         
    }
    
  
    
    
    
    // Reset tree height so that the youngest leaf has height 0
    private void setNodeHeights(Tree tree) {
    	double minHeight = Double.POSITIVE_INFINITY;
		for (Node node : tree.getNodesAsArray()) {
			minHeight = Math.min(minHeight, node.getHeight());
		}
		for (Node node : tree.getNodesAsArray()) {
			node.setHeight(node.getHeight() - minHeight);
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
    
    
    public static Tree convert(RootnessNode nodeA, RootnessNode nodeB, double rootHeight) {
        // Create new root
        Node root = new Node();
        root.setHeight(rootHeight);

        // Create two subtrees, one descending from nodeA and one from nodeB
        Node left = buildSubtree(nodeA, rootHeight, nodeB, new HashMap<>(), new HashSet<>());
        Node right = buildSubtree(nodeB, rootHeight, nodeA, new HashMap<>(), new HashSet<>());

        root.addChild(left);
        root.addChild(right);

        return new Tree(root);
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
    private static Node buildSubtree(RootnessNode current, double parentHeight, RootnessNode cameFrom, Map<RootnessNode, Node> mapping, Set<RootnessNode> visited) {
        if (mapping.containsKey(current)) {
            return mapping.get(current);
        }
        visited.add(current);

        Node thisNode = new Node();
        if (current.getLabel() != null) {
            thisNode.setID(current.getLabel());
        }

        // Height = parentHeight - branchLength
        double length = current.getBranchLengthTo(cameFrom);
        double height = parentHeight - length;
        thisNode.setHeight(height);
        mapping.put(current, thisNode);

        // Recurse into neighbours (except the one we came from)
        for (RootnessNode neigh : current.getNeighbours()) {
            if (neigh == cameFrom) continue;
            Node child = buildSubtree(neigh, height, current, mapping, visited);
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
    
    
    public static MidpointResult findMidpoint(RootnessNode anyNode) {
        // Step 1: find one farthest leaf
    	RootnessNode leaf1 = farthestLeaf(anyNode).node;

        // Step 2: from leaf1, find farthest leaf (this is leaf2)
        FarthestResult res = farthestLeaf(leaf1);
        RootnessNode leaf2 = res.node;
        double diameter = res.distance;

        // Step 3: get path leaf1 -> leaf2
        List<RootnessNode> path = getPath(leaf1, leaf2);

        // Step 4: walk halfway along the path
        double halfway = diameter / 2.0;
        double distSoFar = 0.0;

        for (int i = 0; i < path.size() - 1; i++) {
        	RootnessNode a = path.get(i);
        	RootnessNode b = path.get(i + 1);
            double len = a.getBranchLengthTo(b);

            if (distSoFar + len >= halfway) {
                double distFromA = halfway - distSoFar;
                return new MidpointResult(a, b, distFromA);
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

        public MidpointResult(RootnessNode a, RootnessNode b, double distFromA) {
            this.nodeA = a;
            this.nodeB = b;
            this.distFromA = distFromA;
        }
    }

    private static class FarthestResult {
    	RootnessNode node;
        double distance;
        FarthestResult(RootnessNode n, double d) { node = n; distance = d; }
    }
    
    

}



