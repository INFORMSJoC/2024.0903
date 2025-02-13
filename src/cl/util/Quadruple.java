package cl.util;

import java.util.Objects;

public class Quadruple<E1, E2, E3, E4> {

	public final E1 first;
	public final E2 second;
	public final E3 third;
	public final E4 fourth;


	public Quadruple(E1 first, E2 second, E3 third, E4 fourth) {
		super();
		this.first = first;
		this.second = second;
		this.third = third;
		this.fourth = fourth;
	}
	
	public static <E1,E2,E3,E4> Quadruple<E1,E2,E3,E4> of(E1 first, E2 second, E3 third, E4 fourth) {
		return new Quadruple<>(first,second,third, fourth);
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) return true;
		if (o == null || getClass() != o.getClass()) return false;
		Quadruple<?, ?, ?, ?> triple = (Quadruple<?, ?, ?, ?>) o;
		return first.equals(triple.first) &&
				second.equals(triple.second) &&
				third.equals(triple.third) &&
				fourth.equals(triple.fourth);
	}

	@Override
	public int hashCode() {
		return Objects.hash(first, second, third, fourth);
	}

	@Override
	public String toString() {
		return "Quadruple [first=" + first + ", second=" + second + ", third=" + third + ", fourth=" + fourth + "]";
	}
}
