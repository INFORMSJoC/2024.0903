package cl.hierarchical.model.extended;

import java.util.List;

import cl.data.Instance;

public class BranchInformationExtra {
	private boolean isBranch;
	private Instance reducedInstance;
	private DualInformation reducedDuals;
	private List<List<Integer>> branchMapping;
	private List<Integer> branchMappingWeights;
	
	public boolean isBranch() {
		return isBranch;
	}
	public void setBranch(boolean isBranch) {
		this.isBranch = isBranch;
	}
	public Instance getReducedInstance() {
		return reducedInstance;
	}
	public void setReducedInstance(Instance reducedInstance) {
		this.reducedInstance = reducedInstance;
	}
	public DualInformation getReducedDuals() {
		return reducedDuals;
	}
	public void setReducedDuals(DualInformation reducedDuals) {
		this.reducedDuals = reducedDuals;
	}
	public List<List<Integer>> getBranchMapping() {
		return branchMapping;
	}
	public void setBranchMapping(List<List<Integer>> branchMapping) {
		this.branchMapping = branchMapping;
	}
	public List<Integer> getBranchMappingWeights() {
		return branchMappingWeights;
	}
	public void setBranchMappingWeights(List<Integer> branchMappingWeights) {
		this.branchMappingWeights = branchMappingWeights;
	}
	
	
}
