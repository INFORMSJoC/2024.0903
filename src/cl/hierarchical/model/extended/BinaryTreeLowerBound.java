package cl.hierarchical.model.extended;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Keep track of the best lowerbound when branching in branch-and-price
 * @author rickw
 *
 */
public class BinaryTreeLowerBound {
	private Node root;
	private Map<Branch, Node> nodeMap;

	public BinaryTreeLowerBound(Branch rootBranch, double value) {
		this.root = new Node(null, value);
		this.nodeMap = new LinkedHashMap<>();
		nodeMap.put(rootBranch, root);
	}

	/**
	 * Go through the tree and update the best found LP values
	 * @param start
	 * @return
	 */
	public void updateBestLP(Branch start) {
		Node curNode = nodeMap.get(start);

		while(true) {
			Node parentNode = curNode.getParent();
			if(parentNode==null || !parentNode.hasTwoChilds()) {
				break;
			}
	
			Node left = parentNode.getLeft();
			Node right = parentNode.getRight();
			
			double bestLP = Math.min(left.getValue(), right.getValue());
			parentNode.setValue(bestLP);
			curNode = parentNode;
		}
	}
	
	public double getBestLP() {
		return root.getValue();
	}

	/**
	 * Forbid is left
	 * @param curParent
	 * @param left
	 * @param value
	 */
	public void addLeft(Branch curParent, Branch left, double value) {
		Node parentNode = nodeMap.get(curParent);
		Node leftNode = new Node(parentNode, value);
		parentNode.addLeft(leftNode);
		nodeMap.put(left, leftNode);
	}

	/**
	 * Enforce is right
	 * @param curParent
	 * @param right
	 * @param value
	 */
	public void addRight(Branch curParent, Branch right, double value) {
		Node parentNode = nodeMap.get(curParent);
		Node rightNode = new Node(parentNode, value);
		parentNode.addRight(rightNode);
		nodeMap.put(right, rightNode);
	}

	private class Node {
		private double value;
		private Node left;
		private Node right;
		private Node parent;

		private Node(Node parent, double value) {
			this.parent = parent;
			this.value = value;
		}

		public void setValue(double value) {
			this.value = value;
		}

		private void addLeft(Node left) {
			this.left = left;
		}

		private void addRight(Node right) {
			this.right = right;
		}

		public boolean hasTwoChilds() {
			if(left!=null && right!=null) {
				return true;
			}
			return false;
		}

		public double getValue() {
			return value;
		}

		public Node getLeft() {
			return left;
		}

		public Node getRight() {
			return right;
		}

		public Node getParent() {
			return parent;
		}
	}
}