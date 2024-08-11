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

impl<F: FiniteField + Display> MultivariatePolynomial<F> {
  pub fn generate_sumcheck_claim(&self) -> F {
    // Compute and return the sum of the polynomial over all boolean inputs
    let variables = self.variables();
    let num_variables = variables.len();
    let mut sum = F::ZERO;

    // Use the get_all_possible_boolean_values function to iterate over all boolean assignments
    for boolean_assignment in get_all_possible_boolean_values(num_variables) {
      let assignment: BTreeMap<_, _> = variables
        .iter()
        .zip(boolean_assignment.iter())
        .map(|(&var, &b)| (var, if b { F::ONE } else { F::ZERO }))
        .collect();
      sum += self.evaluate(&assignment);
    }

    sum
  }

  pub fn prove_first_sumcheck_round(&self) -> (F, MultivariatePolynomial<F>) {
    let variables = self.variables();
    let num_variables = variables.len();

    let sum = get_all_possible_boolean_values(num_variables)
      .map(|bool_values| {
        let mut assignment = BTreeMap::new();
        for (i, &b) in bool_values.iter().enumerate() {
          assignment.insert(variables[i].clone(), if b { F::ONE } else { F::ZERO });
        }
        self.evaluate(&assignment)
      })
      .sum();

    // Compute the univariate polynomial g_1(X1)
    let univariate_poly = self.compute_univariate_polynomial(0, vec![]);

    (sum, univariate_poly)
  }

  pub fn prove_sumcheck_round_i(
    &self,
    i: usize,
    partial_assignment: Vec<F>, 
  ) -> MultivariatePolynomial<F> {
    return self.compute_univariate_polynomial(i, partial_assignment);
  }

  pub fn prove_sumcheck_last_round(
    &self,
    i: usize,
    partial_assignment: Vec<F>,
  ) -> MultivariatePolynomial<F> {
    return self.compute_univariate_polynomial(i, partial_assignment);
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
        let further_assignments: Vec<F> =
          bool_values.iter().map(|&b| if b { F::ONE } else { F::ZERO }).collect();
        let further_variables = ((round + 1)..num_variables).zip(further_assignments).collect::<Vec<_>>();
        let poly = partial_poly.clone().apply_variables(
          &further_variables,
        );
        poly
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

pub fn verify_sumcheck_first_round<F: FiniteField + Random>(
    claimed_sum: F,
    univariate_poly: &MultivariatePolynomial<F>
) -> (bool, F) {
    // Step 1: Verify that the polynomial is univariate (has only one variable)
    if univariate_poly.variables().len() != 1 {
        return (false, F::ZERO);
    }

    // Step 2: Verify that g(0) + g(1) = claimed_sum
    let var = 0;
    let sum_at_endpoints = univariate_poly.evaluate(&[(var, F::ZERO)].into_iter().collect())
        + univariate_poly.evaluate(&[(var, F::ONE)].into_iter().collect());


    if sum_at_endpoints != claimed_sum {
        return (false, F::ZERO);
    }

    // Step 3: Generate a random challenge
    let mut rng = rand::thread_rng();
    let challenge: F = F::random(&mut rng);

    // Step 4: Evaluate the polynomial at the challenge point
    let eval_at_challenge = univariate_poly.evaluate(&[(var, challenge)].into_iter().collect());

    // Return true (verification passed) and the evaluation at the challenge point
    (true, eval_at_challenge)
}

// Verify the i-th round of the sumcheck protocol
pub fn verify_sumcheck_round_i<F: FiniteField + Random>(
    round: usize,
    previous_random_value: F,
    previous_univariate_poly: &MultivariatePolynomial<F>,
    current_univariate_poly: &MultivariatePolynomial<F>,
) -> (bool, F) {
    // Step 1: Verify that the current polynomial is univariate
    if current_univariate_poly.variables().len() != 1 {
        return (false, F::ZERO);
    }

    // Step 2: Verify that g_i(r_{i-1}) = g_{i-1}(0) + g_{i-1}(1)
    let prev_var = round - 1;
    let sum_at_endpoints = previous_univariate_poly.evaluate(&[(prev_var, F::ZERO)].into_iter().collect())
        + previous_univariate_poly.evaluate(&[(prev_var, F::ONE)].into_iter().collect());

    let eval_at_previous_challenge = current_univariate_poly.evaluate(&[(prev_var, previous_random_value)].into_iter().collect());

    if eval_at_previous_challenge != sum_at_endpoints {
        return (false, F::ZERO);
    }

    // Step 3: Generate a new random challenge
    let mut rng = rand::thread_rng();
    let new_challenge: F = F::random(&mut rng);

    // Step 4: Evaluate the current polynomial at the new challenge point
    let eval_at_new_challenge = current_univariate_poly.evaluate(&[(round, new_challenge)].into_iter().collect());

    // Return true (verification passed) and the evaluation at the new challenge point
    (true, eval_at_new_challenge)
}

// You first verify the last rounds by ensuring that the last values fulfill the equality of g_{v−1}(r_{v−1}) = gv(0) +gv(1).
// Then, it randomly gets a random value. Then plugs it into the current univariate polynomial
pub fn verify_sumcheck_last_round<F: FiniteField + Random>(
    round: usize,
    previous_random_value: F,
    previous_univariate_poly: &MultivariatePolynomial<F>,
    current_univariate_poly: &MultivariatePolynomial<F>,
) -> (bool, F) {
    // Step 1: Verify that the current polynomial is univariate
    if current_univariate_poly.variables().len() != 1 {
        return (false, F::ZERO);
    }

    // Step 2: Verify that g_{v-1}(r_{v-1}) = g_v(0) + g_v(1)
    let prev_var = round - 1;
    let sum_at_endpoints = current_univariate_poly.evaluate(&[(prev_var, F::ZERO)].into_iter().collect())
        + current_univariate_poly.evaluate(&[(prev_var, F::ONE)].into_iter().collect());

    let eval_at_previous_challenge = previous_univariate_poly.evaluate(&[(prev_var, previous_random_value)].into_iter().collect());

    if eval_at_previous_challenge != sum_at_endpoints {
        return (false, F::ZERO);
    }

    // Step 3: Generate a new random challenge
    let mut rng = thread_rng();
    let new_challenge: F = F::random(&mut rng);

    // Step 4: Evaluate the current polynomial at the new challenge point
    let eval_at_new_challenge = current_univariate_poly.evaluate(&[(round, new_challenge)].into_iter().collect());

    // Return true (verification passed) and the evaluation at the new challenge point
    (true, eval_at_new_challenge)
}

pub fn simulate_sumcheck_protocol<F: FiniteField + Random + Display>(
    polynomial: &MultivariatePolynomial<F>
) -> bool {
    let num_variables = polynomial.variables().len();
    let mut challenges = Vec::new();

    // First round
    let (claimed_sum, first_univariate_poly) = polynomial.prove_first_sumcheck_round();
    let (valid, challenge) = verify_sumcheck_first_round(claimed_sum, &first_univariate_poly);
    if !valid {
        return false;
    }
    challenges.push(challenge);

    let mut previous_univariate_poly = first_univariate_poly;

    // Intermediate rounds
    for i in 1..num_variables - 1 {
        let univariate_poly = polynomial.prove_sumcheck_round_i(i, challenges.clone());
        let (valid, challenge) = verify_sumcheck_round_i(
            i,
            *challenges.last().unwrap(),
            &previous_univariate_poly,
            &univariate_poly,
        );
        if !valid {
            return false;
        }
        challenges.push(challenge);
        previous_univariate_poly = univariate_poly;
    }

    // Last round
    let last_univariate_poly = polynomial.prove_sumcheck_last_round(num_variables - 1, challenges.clone());
    let (valid, final_challenge) = verify_sumcheck_last_round(
        num_variables - 1,
        *challenges.last().unwrap(),
        &previous_univariate_poly,
        &last_univariate_poly,
    );
    if !valid {
        return false;
    }
    challenges.push(final_challenge);

    // Final check
    
    let challenges_with_indices: BTreeMap<usize, F> = challenges.iter().cloned().enumerate().collect();
    let final_evaluation = polynomial.evaluate(&challenges_with_indices);
    let expected_evaluation = last_univariate_poly.evaluate(&[(num_variables - 1, final_challenge)].into_iter().collect());

    final_evaluation == expected_evaluation
}