package cl.hierarchical.model.extended;

public class Edge_MIP_Colgen {
	private Node_MIP_Colgen origin, destination;

	public Edge_MIP_Colgen(Node_MIP_Colgen origin, Node_MIP_Colgen destination) {
		this.origin = origin;
		this.destination = destination;
	}

	public Node_MIP_Colgen getOrigin() {
		return origin;
	}

	public Node_MIP_Colgen getDestination() {
		return destination;
	}
}
