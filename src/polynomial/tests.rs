use std::collections::{BTreeMap, HashMap};

use self::sumcheck::Random;

use super::*;
use crate::polynomial::{
  multivariate_polynomial::{MultivariatePolynomial, MultivariateTerm, MultivariateVariable}, sumcheck::verify_sumcheck_first_round,
};

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

  assert_eq!(
    poly.variables().into_iter().collect::<std::collections::HashSet<_>>(),
    vec![0, 1].into_iter().collect::<std::collections::HashSet<_>>()
  );
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

  println!("Addition Result polynomial: {}", result);

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

  println!("Multiplication Result polynomial: {}", result);

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

  let points = BTreeMap::from([(0, PlutoBaseField::new(2)), (1, PlutoBaseField::new(3))]);

  println!("{}", poly);

  let result = poly.evaluate(&points);
  // 3*(2^2)*(3) + 2*(2) + 1 = 3*4*3 + 4 + 1 = 36 + 4 + 1 = 41
  assert_eq!(result, PlutoBaseField::new(41));
}

#[test]
fn test_apply_variables_single_variable() {
  let mut poly = MultivariatePolynomial::<PlutoBaseField>::new();
  poly.insert_term(BTreeMap::from([(0, 2)]), PlutoBaseField::new(3)); // 3x_0^2
  poly.insert_term(BTreeMap::from([(0, 1)]), PlutoBaseField::new(2)); // 2x_0
  poly.insert_term(BTreeMap::new(), PlutoBaseField::new(1)); // 1

  let variables = vec![(0, PlutoBaseField::new(2))];
  let result = poly.apply_variables(&variables);

  assert_eq!(result.coefficient(&BTreeMap::new()), Some(&PlutoBaseField::new(17))); // 3*2^2 + 2*2 +
                                                                                    // 1 = 17
}

#[test]
fn test_apply_variables_multiple_variables() {
  let mut poly = MultivariatePolynomial::<PlutoBaseField>::new();
  poly.insert_term(BTreeMap::from([(0, 2), (1, 1)]), PlutoBaseField::new(3)); // 3x_0^2*x_1
  poly.insert_term(BTreeMap::from([(0, 1)]), PlutoBaseField::new(2)); // 2x_0
  poly.insert_term(BTreeMap::from([(1, 1)]), PlutoBaseField::new(4)); // 4x_1
  poly.insert_term(BTreeMap::new(), PlutoBaseField::new(1)); // 1

  println!("Apply Multiple Variables Polynomial: {}", poly);

  let variables = vec![(0, PlutoBaseField::new(2)), (1, PlutoBaseField::new(3))];
  let result = poly.apply_variables(&variables);

  println!("Reduced Multiple Variable Polynomial: {}", result);

  assert_eq!(result.coefficient(&BTreeMap::new()), Some(&PlutoBaseField::new(53))); // 3*2^2*3 + 2*2
                                                                                    // + 4*3 + 1 =
                                                                                    // 53
}

#[test]
fn test_apply_variables_partial_application() {
  let mut poly = MultivariatePolynomial::<PlutoBaseField>::new();
  poly.insert_term(BTreeMap::from([(0, 2), (1, 1), (2, 1)]), PlutoBaseField::new(3)); // 3x_0^2*x_1*x_2
  poly.insert_term(BTreeMap::from([(0, 1), (2, 1)]), PlutoBaseField::new(2)); // 2x_0*x_2
  poly.insert_term(BTreeMap::from([(1, 1)]), PlutoBaseField::new(4)); // 4x_1

  let variables = vec![(0, PlutoBaseField::new(2)), (1, PlutoBaseField::new(3))];
  let result = poly.apply_variables(&variables);

  assert_eq!(result.coefficient(&BTreeMap::from([(2, 1)])), Some(&PlutoBaseField::new(40))); // 3*2^2*3*x_2 + 2*2*x_2 = 36x_2 + 4x_2 = 40x_2
  assert_eq!(result.coefficient(&BTreeMap::new()), Some(&PlutoBaseField::new(12))); // 4*3 = 12
}

#[test]
fn test_apply_variables_no_effect() {
  let mut poly = MultivariatePolynomial::<PlutoBaseField>::new();
  poly.insert_term(BTreeMap::from([(0, 2), (1, 1)]), PlutoBaseField::new(3)); // 3x_0^2*x_1
  poly.insert_term(BTreeMap::from([(2, 1)]), PlutoBaseField::new(2)); // 2x_2

  let variables = vec![(3, PlutoBaseField::new(5))];
  let result = poly.apply_variables(&variables);

  assert_eq!(result.coefficient(&BTreeMap::from([(0, 2), (1, 1)])), Some(&PlutoBaseField::new(3)));
  assert_eq!(result.coefficient(&BTreeMap::from([(2, 1)])), Some(&PlutoBaseField::new(2)));
}

#[test]
fn test_apply_variables_empty() {
  let mut poly = MultivariatePolynomial::<PlutoBaseField>::new();
  poly.insert_term(BTreeMap::from([(0, 2), (1, 1)]), PlutoBaseField::new(3)); // 3x_0^2*x_1

  let variables = vec![];
  let result = poly.apply_variables(&variables);

  assert_eq!(result.coefficient(&BTreeMap::from([(0, 2), (1, 1)])), Some(&PlutoBaseField::new(3)));
}

impl Random for PlutoBaseField {
  fn random<R: Rng + ?Sized>(rng: &mut R) -> Self {
      let value = rng.gen_range(0..PlutoPrime::Base as usize);
      PlutoBaseField::new(value)
  }
}

#[test]
fn test_first_round_sumcheck() {
  // Create the polynomial f(x1, x2, x3) = x1 * (x2 + x3) - (x2 * x3)
  let poly = MultivariatePolynomial::<PlutoBaseField>::from_terms(vec![
    MultivariateTerm::new(
      vec![MultivariateVariable::new(0, 1), MultivariateVariable::new(1, 1)],
      PlutoBaseField::ONE,
    ), // x1 * x2
    MultivariateTerm::new(
      vec![MultivariateVariable::new(0, 1), MultivariateVariable::new(2, 1)],
      PlutoBaseField::ONE,
    ), // x1 * x3
    MultivariateTerm::new(
      vec![MultivariateVariable::new(1, 1), MultivariateVariable::new(2, 1)],
      -PlutoBaseField::ONE,
    ), // -x2 * x3
  ]);

  // Compute the first round of sumcheck
  let (claimed_sum, first_univariate_poly) = poly.prove_first_sumcheck_round();

  // Print the claimed sum and the first univariate polynomial
  println!("Claimed sum: {:?}", claimed_sum);
  println!("First univariate polynomial: {}", first_univariate_poly);

  // Verify the first round
  let (valid, challenge) = verify_sumcheck_first_round(claimed_sum, &first_univariate_poly);

  assert!(valid, "First round verification failed");
  println!("First round verification passed");
  println!("Challenge for next round: {:?}", challenge);

  // Additional checks
  // 1. Check that the claimed sum is correct
  let expected_sum = poly.generate_sumcheck_claim();
  assert_eq!(claimed_sum, expected_sum, "Claimed sum is incorrect");

  // 2. Check that the first univariate polynomial is correct
  let x1 = first_univariate_poly.variables()[0];
  let g_0 = first_univariate_poly.evaluate(&[(x1, PlutoBaseField::ZERO)].into_iter().collect());
  let g_1 = first_univariate_poly.evaluate(&[(x1, PlutoBaseField::ONE)].into_iter().collect());
  assert_eq!(g_0 + g_1, claimed_sum, "g(0) + g(1) should equal the claimed sum");

  // 3. Check that the degree of the first univariate polynomial is at most 1
  assert!(
    first_univariate_poly.degree() <= 1,
    "Degree of first univariate polynomial should be at most 1"
  );
}
