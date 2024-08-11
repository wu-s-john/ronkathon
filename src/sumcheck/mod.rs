use std::{fmt::Display, hash::{Hash, Hasher}};

use rand::{Rng, SeedableRng};
use crate::{algebra::field::FiniteField, polynomial::{multivariate_polynomial::MultivariatePolynomial, to_bytes::ToBytes}, random::{Random, RandomOracle}};

pub mod boolean_array;
#[cfg(test)] mod tests;

use self::boolean_array::get_all_possible_boolean_values;

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

impl<F: FiniteField + Random> RandomOracle for F {
    fn random_oracle<R: Rng + ?Sized>(rng: &mut R, input: &[u8]) -> Self {
        // This is a simplified example. In a real implementation,
        // you'd want to use a cryptographic hash function here.
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        input.hash(&mut hasher);
        let hash = hasher.finish();

        // Use the hash to seed a new RNG
        let mut seeded_rng = rand::rngs::StdRng::seed_from_u64(hash);

        // Generate a random field element using the seeded RNG
        Self::random(&mut seeded_rng)
    }
}

impl<F: FiniteField + Display> MultivariatePolynomial<F> {
  pub fn generate_sumcheck_claim(&self) -> F {
    // Compute and return the sum of the polynomial over all boolean inputs
    let variables = self.variables();
    let num_variables = variables.len();
    let mut sum = F::ZERO;

    // Use the get_all_possible_boolean_values function to iterate over all boolean assignments
    for boolean_assignment in get_all_possible_boolean_values(num_variables) {
      let assignment: Vec<(usize, F)> = variables
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
        let assignment: Vec<(usize, F)> = variables.iter().enumerate()
          .map(|(i, &var)| (var, if bool_values[i] { F::ONE } else { F::ZERO }))
          .collect();
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
    assert!(
      result_polynomial.variables().len() <= 1,
      "The univariate polynomial should have at most one variable"
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
    let sum_at_endpoints = univariate_poly.evaluate(&[(var, F::ZERO)]) + univariate_poly.evaluate(&[(var, F::ONE)]);


    if sum_at_endpoints != claimed_sum {
        return (false, F::ZERO);
    }

    // Step 3: Generate a random challenge
    let mut rng = rand::thread_rng();
    let challenge: F = F::random(&mut rng);

    // Return true (verification passed) and the evaluation at the challenge point
    (true, challenge)
}

// Verify the i-th round of the sumcheck protocol
pub fn verify_sumcheck_univariate_poly_sum<F: FiniteField + Random>(
    round: usize,
    challenge: F,
    previous_univariate_poly: &MultivariatePolynomial<F>,
    current_univariate_poly: &MultivariatePolynomial<F>,
) -> (bool, F) {
    // Step 1: Verify that the current polynomial is univariate
    if current_univariate_poly.variables().len() > 1 {
        return (false, F::ZERO);
    }

    // Step 2: Verify that g_i(r_{i-1}) = g_{i-1}(0) + g_{i-1}(1)
    let prev_var = round - 1;
    let sum_at_endpoints = previous_univariate_poly.evaluate(&[(prev_var, challenge)]);

    let eval_at_previous_challenge = current_univariate_poly.evaluate(&[(round, F::ZERO)])
        + current_univariate_poly.evaluate(&[(round, F::ONE)]);
    if eval_at_previous_challenge != sum_at_endpoints {
        return (false, F::ZERO);
    }

    // Step 3: Generate a new random challenge
    let mut rng = rand::thread_rng();
    let new_challenge: F = F::random(&mut rng);

    // Return true (verification passed) and the evaluation at the new challenge point
    (true, new_challenge)
}

// You would just need to verify that the last unit polynomial applied that you created with the last challenge is equal to applying 
// all the challenges applied
// Then, it randomly gets a random value. Then plugs it into the current univariate polynomial
pub fn verify_sumcheck_last_round<F: FiniteField + Random>(
    challenges: Vec<F>,
    univariate_poly: &MultivariatePolynomial<F>,
    poly: &MultivariatePolynomial<F>,
) -> bool {
    // Step 1: Apply all challenges to the original polynomial
    let mut challenges_with_indices = Vec::new();
    for (i, challenge) in challenges.iter().enumerate() {
        challenges_with_indices.push((i, *challenge));
    }
    let poly_evaluation = poly.evaluate(&challenges_with_indices);

    // Step 2: Generate a random challenge for the last variable
    let mut rng = rand::thread_rng();
    let last_challenge: F = F::random(&mut rng);

    // Step 3: Evaluate the univariate polynomial at the last challenge
    let last_var = challenges.len();
    let univariate_evaluation = univariate_poly.evaluate(&[(last_var, last_challenge)]);

    // Step 4: Compare the evaluations
    poly_evaluation == univariate_evaluation
}

impl<F: FiniteField + Display> ToBytes for F {
    fn to_bytes(&self) -> Vec<u8> {
        // Implement this based on how your field elements are represented
        // This is just an example:
        self.to_string().into_bytes()
    }
}

impl<F: FiniteField + Display> ToBytes for MultivariatePolynomial<F> {
    fn to_bytes(&self) -> Vec<u8> {
        // Implement this based on how your polynomials are represented
        // This is just an example:
        self.to_string().into_bytes()
    }
}


pub fn non_interactive_sumcheck_prove<F: FiniteField + Random + RandomOracle + Display + ToBytes>(
    polynomial: &MultivariatePolynomial<F>
) -> SumcheckProof<F> {
    let num_variables = polynomial.variables().len();
    let mut challenges = Vec::new();
    let mut round_polynomials = Vec::new();
    let mut round_evaluations = Vec::new();

    // First round
    let (claimed_sum, first_univariate_poly) = polynomial.prove_first_sumcheck_round();
    round_polynomials.push(first_univariate_poly.clone());

    // Generate challenge using random oracle
    let mut rng = rand::thread_rng();
    let challenge: F = F::random_oracle(&mut rng, &claimed_sum.to_bytes());
    challenges.push(challenge);
    round_evaluations.push(first_univariate_poly.evaluate(&[(0, challenge)]));

    let mut previous_univariate_poly = first_univariate_poly;

    // Intermediate rounds
    for i in 1..num_variables {
        let univariate_poly = polynomial.prove_sumcheck_round_i(i, challenges.clone());
        round_polynomials.push(univariate_poly.clone());

        // Generate challenge using random oracle
        let challenge: F = F::random_oracle(&mut rng, &previous_univariate_poly.to_bytes());
        challenges.push(challenge);
        round_evaluations.push(univariate_poly.evaluate(&[(i, challenge)]));

        previous_univariate_poly = univariate_poly;
    }

    // Final evaluation
    let final_point = challenges.clone();
    let final_evaluation = polynomial.evaluate(&final_point.iter().cloned().enumerate().collect::<Vec<_>>());

    SumcheckProof {
        claimed_sum,
        round_polynomials,
        round_evaluations,
        final_point,
        final_evaluation,
    }
}

pub fn non_interactive_sumcheck_verify<F: FiniteField + Random + RandomOracle + Display>(
    proof: &SumcheckProof<F>,
    polynomial: &MultivariatePolynomial<F>
) -> bool {
    let num_variables = polynomial.variables().len();

    // Verify first round
    let (valid, challenge) = verify_sumcheck_first_round(proof.claimed_sum, &proof.round_polynomials[0]);
    if !valid {
        return false;
    }

    // Verify that the challenge matches what the prover used
    let mut rng = rand::thread_rng();
    let expected_challenge: F = F::random_oracle(&mut rng, &proof.claimed_sum.to_bytes());
    if challenge != expected_challenge {
        return false;
    }

    // Verify intermediate rounds
    for i in 1..num_variables {
        let (valid, challenge) = verify_sumcheck_univariate_poly_sum(
            i,
            proof.final_point[i-1],
            &proof.round_polynomials[i-1],
            &proof.round_polynomials[i],
        );
        if !valid {
            return false;
        }

        // Verify that the challenge matches what the prover used
        let expected_challenge: F = F::random_oracle(&mut rng, &proof.round_polynomials[i-1].to_bytes());
        if challenge != expected_challenge {
            return false;
        }
    }

    // Verify last round
    verify_sumcheck_last_round(
        proof.final_point.clone(),
        &proof.round_polynomials.last().unwrap(),
        polynomial,
    )
}