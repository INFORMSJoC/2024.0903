package cl.hierarchical.model.pdcg;

import org.apache.commons.math3.linear.RealMatrix;

public class PrimalDualPoint {
	private RealMatrix pointX, pointPi, pointS;

	public PrimalDualPoint(RealMatrix pointX, RealMatrix pointPi, RealMatrix pointS) {
		this.pointX = pointX;
		this.pointPi = pointPi;
		this.pointS = pointS;
	}

	public RealMatrix getPointX() {
		return pointX;
	}

	public RealMatrix getPointPi() {
		return pointPi;
	}

	public RealMatrix getPointS() {
		return pointS;
	}
	
	
}
