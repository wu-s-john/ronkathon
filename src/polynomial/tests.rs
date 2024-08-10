use std::collections::{BTreeMap, HashMap};

use crate::polynomial::multivariate_polynomial::MultivariatePolynomial;

use super::*;

#[fixture]
fn poly() -> Polynomial<Monomial, PlutoBaseField, 4> {
  // Coefficients of the polynomial 1 + 2x + 3x^2 + 4x^3
  let a = PlutoBaseField::new(1);
  let b = PlutoBaseField::new(2);
  let c = PlutoBaseField::new(3);
  let d = PlutoBaseField::new(4);
  Polynomial::<Monomial, PlutoBaseField, 4>::new([a, b, c, d])
}

#[rstest]
fn evaluation(poly: Polynomial<Monomial, PlutoBaseField, 4>) {
  // Should get: 1 + 2*(2) + 3*(2)^2 + 4*(2)^3 = 49
  let y = poly.evaluate(PlutoBaseField::new(2));
  let r = PlutoBaseField::new(49);
  assert_eq!(y, r);
}

#[test]
fn evaluation_with_zero() {
  // Coefficients of the polynomial 1 + 3x^2
  let a = PlutoBaseField::new(1);
  let b = PlutoBaseField::new(0);
  let c = PlutoBaseField::new(3);
  let polynomial = Polynomial::<Monomial, PlutoBaseField, 3>::new([a, b, c]);
  let y = polynomial.evaluate(PlutoBaseField::new(0));

  // Should get: 1 + 3(0)^2 = 1
  let r = PlutoBaseField::new(1);
  assert_eq!(y, r);
}

#[rstest]
fn lagrange_evaluation(poly: Polynomial<Monomial, PlutoBaseField, 4>) {
  // Convert to Lagrange basis using roots of unity
  let lagrange = poly.dft();
  println!("{}", lagrange);

  // Should get: 1 + 2*(2) + 3*(2)^2 + 4*(2)^3= 49
  let r = lagrange.evaluate(PlutoBaseField::new(2));
  assert_eq!(r, PlutoBaseField::new(49));
}

#[test]
#[should_panic]
fn no_roots_of_unity() {
  // Coefficients of the polynomial 1 + 2x
  let a = PlutoBaseField::new(1);
  let b = PlutoBaseField::new(2);
  let c = PlutoBaseField::new(3);
  let polynomial = Polynomial::<Monomial, PlutoBaseField, 3>::new([a, b, c]);
  polynomial.dft();
}

#[rstest]
fn check_coefficients(poly: Polynomial<Monomial, PlutoBaseField, 4>) {
  assert_eq!(poly.coefficients, [
    PlutoBaseField::new(1),
    PlutoBaseField::new(2),
    PlutoBaseField::new(3),
    PlutoBaseField::new(4)
  ]);

  assert_eq!(poly.dft().coefficients, [
    PlutoBaseField::new(10),
    PlutoBaseField::new(79),
    PlutoBaseField::new(99),
    PlutoBaseField::new(18)
  ]);
}

#[rstest]
fn degree(poly: Polynomial<Monomial, PlutoBaseField, 4>) {
  assert_eq!(poly.degree(), 3);
}

#[rstest]
fn leading_coefficient(poly: Polynomial<Monomial, PlutoBaseField, 4>) {
  assert_eq!(poly.leading_coefficient(), PlutoBaseField::new(4));
}

#[rstest]
fn pow_mult(poly: Polynomial<Monomial, PlutoBaseField, 4>) {
  assert_eq!(
    poly.pow_mult::<2>(PlutoBaseField::new(5)),
    Polynomial::<Monomial, PlutoBaseField, 6>::new([
      PlutoBaseField::new(0),
      PlutoBaseField::new(0),
      PlutoBaseField::new(5),
      PlutoBaseField::new(10),
      PlutoBaseField::new(15),
      PlutoBaseField::new(20)
    ])
  );
}

#[test]
fn trim_zeros() {
  let mut coefficients = vec![
    PlutoBaseField::new(1),
    PlutoBaseField::new(2),
    PlutoBaseField::new(3),
    PlutoBaseField::new(4),
    PlutoBaseField::ZERO,
  ];
  Polynomial::<Monomial, PlutoBaseField, 5>::trim_zeros(coefficients.as_mut());
  assert_eq!(coefficients, [
    PlutoBaseField::new(1),
    PlutoBaseField::new(2),
    PlutoBaseField::new(3),
    PlutoBaseField::new(4)
  ]);
}

#[rstest]
fn dft(poly: Polynomial<Monomial, PlutoBaseField, 4>) {
  assert_eq!(poly.dft().coefficients, [
    PlutoBaseField::new(10),
    PlutoBaseField::new(79),
    PlutoBaseField::new(99),
    PlutoBaseField::new(18)
  ]);
  // let poly =
  //   Polynomial::<Monomial, PlutoBaseField>::new(vec![PlutoBaseField::ZERO,
  // PlutoBaseField::ZERO]); assert_eq!(poly.coefficients, [PlutoBaseField::ZERO]);
}


#[test]
fn test_multivariate_polynomial_creation() {
    let mut poly = MultivariatePolynomial::<PlutoBaseField>::new();
    poly.insert_term(BTreeMap::from([(0, 2), (1, 1)]), PlutoBaseField::new(3)); // 3x_0^2*x_1
    poly.insert_term(BTreeMap::from([(0, 1)]), PlutoBaseField::new(2)); // 2x_0
    poly.insert_term(BTreeMap::new(), PlutoBaseField::new(1)); // 1

    assert_eq!(poly.degree(), 3);
    assert_eq!(poly.variables(), vec![0, 1]);
}

#[test]
fn test_multivariate_polynomial_addition() {
    let mut poly1 = MultivariatePolynomial::<PlutoBaseField>::new();
    poly1.insert_term(BTreeMap::from([(0, 2)]), PlutoBaseField::new(1)); // x_0^2
    poly1.insert_term(BTreeMap::from([(1, 1)]), PlutoBaseField::new(2)); // 2x_1

    let mut poly2 = MultivariatePolynomial::<PlutoBaseField>::new();
    poly2.insert_term(BTreeMap::from([(0, 2)]), PlutoBaseField::new(3)); // 3x_0^2
    poly2.insert_term(BTreeMap::from([(2, 1)]), PlutoBaseField::new(4)); // 4x_2

    let result = poly1 + poly2;

    assert_eq!(result.coefficient(&BTreeMap::from([(0, 2)])), Some(&PlutoBaseField::new(4))); // 4x_0^2
    assert_eq!(result.coefficient(&BTreeMap::from([(1, 1)])), Some(&PlutoBaseField::new(2))); // 2x_1
    assert_eq!(result.coefficient(&BTreeMap::from([(2, 1)])), Some(&PlutoBaseField::new(4))); // 4x_2
}

#[test]
fn test_multivariate_polynomial_multiplication() {
    let mut poly1 = MultivariatePolynomial::<PlutoBaseField>::new();
    poly1.insert_term(BTreeMap::from([(0, 1)]), PlutoBaseField::new(2)); // 2x_0
    poly1.insert_term(BTreeMap::new(), PlutoBaseField::new(1)); // 1

    let mut poly2 = MultivariatePolynomial::<PlutoBaseField>::new();
    poly2.insert_term(BTreeMap::from([(1, 1)]), PlutoBaseField::new(3)); // 3x_1
    poly2.insert_term(BTreeMap::new(), PlutoBaseField::new(1)); // 1

    let result = poly1 * poly2;

    assert_eq!(result.coefficient(&BTreeMap::from([(0, 1), (1, 1)])), Some(&PlutoBaseField::new(6))); // 6x_0*x_1
    assert_eq!(result.coefficient(&BTreeMap::from([(0, 1)])), Some(&PlutoBaseField::new(2))); // 2x_0
    assert_eq!(result.coefficient(&BTreeMap::from([(1, 1)])), Some(&PlutoBaseField::new(3))); // 3x_1
    assert_eq!(result.coefficient(&BTreeMap::new()), Some(&PlutoBaseField::new(1))); // 1
}

#[test]
fn test_multivariate_polynomial_evaluation() {
    let mut poly = MultivariatePolynomial::<PlutoBaseField>::new();
    poly.insert_term(BTreeMap::from([(0, 2), (1, 1)]), PlutoBaseField::new(3)); // 3x_0^2*x_1
    poly.insert_term(BTreeMap::from([(0, 1)]), PlutoBaseField::new(2)); // 2x_0
    poly.insert_term(BTreeMap::new(), PlutoBaseField::new(1)); // 1

    let points = BTreeMap::from([
        (0, PlutoBaseField::new(2)),
        (1, PlutoBaseField::new(3)),
    ]);

    let result = poly.evaluate(&points);
    // 3*(2^2)*(3) + 2*(2) + 1 = 3*4*3 + 4 + 1 = 36 + 4 + 1 = 41
    assert_eq!(result, PlutoBaseField::new(41));
}

#[test]
fn test_multivariate_polynomial_display() {
    let mut poly = MultivariatePolynomial::<PlutoBaseField>::new();
    poly.insert_term(BTreeMap::from([(0, 2), (1, 1)]), PlutoBaseField::new(3)); // 3x_0^2*x_1
    poly.insert_term(BTreeMap::from([(0, 1)]), PlutoBaseField::new(2)); // 2x_0
    poly.insert_term(BTreeMap::new(), PlutoBaseField::new(1)); // 1

    assert_eq!(poly.to_string(), "3*x_0^2*x_1 + 2*x_0 + 1");
}