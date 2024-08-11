use std::collections::BTreeMap;

use rand::thread_rng;

use self::boolean_array::get_all_possible_boolean_values;
use super::{multivariate_polynomial::MultivariatePolynomial, *};

pub struct SumcheckProof<F: FiniteField> {
  // Define the structure of your proof here
  pub claimed_sum: F,

  // Vector of univariate polynomials, one for each round
  pub round_polynomials: Vec<MultivariatePolynomial<F>>,

  // Vector of evaluations of the round polynomials at the challenge points
  pub round_evaluations: Vec<F>,

  // The final evaluation point (all challenges combined)
  pub final_point: Vec<F>,

  // The final evaluation of the original multivariate polynomial
  pub final_evaluation: F,
}

pub trait Random {
  fn random<R: Rng + ?Sized>(rng: &mut R) -> Self;
}

impl<F: FiniteField + Random> MultivariatePolynomial<F> {
  pub fn generate_sumcheck_claim(&self) -> F {
    // Compute and return the sum of the polynomial over all boolean inputs
    let variables = self.variables();
    let num_variables = variables.len();
    let mut sum = F::ZERO;

    // Iterate over all possible boolean assignments
    for i in 0..(1 << num_variables) {
      let mut assignment = BTreeMap::new();
      for (j, &var) in variables.iter().enumerate() {
        assignment.insert(var, if (i & (1 << j)) != 0 { F::ONE } else { F::ZERO });
      }
      sum += self.evaluate(&assignment);
    }

    sum
  }

  pub fn sumcheck_prover(&self) -> SumcheckProof<F> {
    let variables = self.variables();
    let num_variables = variables.len();
    let mut round_polynomials = Vec::new();
    let mut round_evaluations = Vec::new();
    let mut partial_assignment = Vec::new();

    // Compute the claimed sum (H in the protocol description)
    let claimed_sum = self.generate_sumcheck_claim();

    let mut rng = thread_rng();

    // Iterate through all rounds
    for round in 0..num_variables {
      // Compute the univariate polynomial for this round
      let univariate_poly = self.compute_univariate_polynomial(round, partial_assignment.clone());
      round_polynomials.push(univariate_poly.clone());

      // In a real implementation, we would receive a challenge from the verifier here
      // For now, we'll generate a random challenge

      let challenge = F::random(rng);
      round_evaluations
        .push(univariate_poly.evaluate(&BTreeMap::from([(variables[round], challenge)])));

      // Add the challenge to our partial assignment for the next round
      partial_assignment.push(challenge);
    }

    // Compute the final evaluation
    let final_point = partial_assignment.clone();
    let final_evaluation =
      self.evaluate(&final_point.iter().cloned().zip(variables).map(|(v, k)| (k, v)).collect());

    SumcheckProof {
      claimed_sum,
      round_polynomials,
      round_evaluations,
      final_point,
      final_evaluation,
    }
  }

  pub fn sumcheck_verifier(&self, proof: &SumcheckProof<F>) -> bool {
    let variables = self.variables();
    let num_variables = variables.len();
    let mut partial_sum = proof.claimed_sum;
    let mut partial_assignment = Vec::new();

    // Verify each round
    for round in 0..num_variables {
      let univariate_poly = &proof.round_polynomials[round];

      // Check if the polynomial is of degree at most 1
      if univariate_poly.degree() > 1 {
        return false;
      }

      // Check if g_i(0) + g_i(1) equals the partial sum from the previous round
      let sum_at_endpoints = univariate_poly
        .evaluate(&BTreeMap::from([(variables[round], F::ZERO)]))
        + univariate_poly.evaluate(&BTreeMap::from([(variables[round], F::ONE)]));
      if sum_at_endpoints != partial_sum {
        return false;
      }

      // Update partial sum for the next round
      partial_sum = proof.round_evaluations[round];

      // Add the challenge to our partial assignment for the next round
      partial_assignment.push(proof.final_point[round]);
    }

    // Perform the final round check
    let final_evaluation =
      self.evaluate(&proof.final_point.iter().cloned().zip(variables).collect());
    if final_evaluation != proof.final_evaluation {
      return false;
    }

    true
  }

  fn compute_univariate_polynomial(
    &self,
    round: usize,
    partial_assignment: Vec<F>,
  ) -> MultivariatePolynomial<F> {
    let variables = self.variables();
    let num_variables = variables.len();

    // Create a polynomial to store the result
    // First create a partial evaluation
    let partial_poly = self.apply_variables(
      &partial_assignment.iter().enumerate().map(|(i, &v)| (i, v)).collect::<Vec<_>>(),
    );

    let result_polynomial = get_all_possible_boolean_values(num_variables - round - 1)
      .map(|bool_values| {
        let further_assignments: Vec<F> = bool_values.iter().map(|&b| if b { F::ONE } else { F::ZERO }).collect();
        partial_poly.clone().apply_variables(
            &(0..num_variables).skip(round).zip(further_assignments).collect::<Vec<_>>()
        )
      })
      .fold(MultivariatePolynomial::new(), |acc, poly| acc + poly);

    // Assert that the resulting polynomial has only one variable
    assert_eq!(
      result_polynomial.variables().len(),
      1,
      "The univariate polynomial should have only one variable"
    );

    result_polynomial
  }

  fn test_polynomial_degree(&self) -> bool { todo!("Implement polynomial degree test") }
}
