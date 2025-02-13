package cl.hierarchical.model.extended;

public class Edge {
	private int origin, destination;

	public Edge(int origin, int destination) {
		this.origin = origin;
		this.destination = destination;
	}

	public int getOrigin() {
		return origin;
	}

	public int getDestination() {
		return destination;
	}
}
