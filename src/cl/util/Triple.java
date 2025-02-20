package cl.util;

import java.util.Objects;

public class Triple<E1, E2, E3> {

	public final E1 first;
	public final E2 second;
	public final E3 third;


	public Triple(E1 first, E2 second, E3 third) {
		super();
		this.first = first;
		this.second = second;
		this.third = third;
	}
	
	public static <E1,E2,E3> Triple<E1,E2,E3> of(E1 first, E2 second, E3 third) {
		return new Triple<>(first,second,third);
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) return true;
		if (o == null || getClass() != o.getClass()) return false;
		Triple<?, ?, ?> triple = (Triple<?, ?, ?>) o;
		return first.equals(triple.first) &&
				second.equals(triple.second) &&
				Objects.equals(third, triple.third);
	}

	@Override
	public int hashCode() {
		return Objects.hash(first, second, third);
	}

	@Override
	public String toString() {
		return "Triple [first=" + first + ", second=" + second + ", third=" + third + ']';
	}
}
