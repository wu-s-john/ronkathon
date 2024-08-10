use std::collections::{HashMap, BTreeMap};
use std::hash::Hash;
use std::ops::{Add, Mul};

use crate::algebra::field::FiniteField;

use super::{Monomial, Polynomial};

use std::hash::Hasher;

type Hi = BTreeMap<usize, usize>;
#[derive(PartialEq, Eq)]
struct SortedBTreeMap<K: Ord, V>(BTreeMap<K, V>);

impl<K: Ord + Hash, V: Hash> Hash for SortedBTreeMap<K, V> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        for (k, v) in self.0.iter() {
            k.hash(state);
            v.hash(state);
        }
    }
}

pub struct MultivariatePolynomial<F: FiniteField> {
    terms: HashMap<BTreeMap<usize, usize>, F>,
}

impl<F: FiniteField> MultivariatePolynomial<F> {
    pub fn new() -> Self {
        Self { terms: HashMap::new() }
    }

    pub fn insert_term(&mut self, exponents: BTreeMap<usize, usize>, coefficient: F) {
        if coefficient != F::ZERO {
            let entry = self.terms.entry(exponents.clone()).or_insert(F::ZERO);
            *entry += coefficient;
            if *entry == F::ZERO {
                self.terms.remove(&exponents);
            }
        }
    }

    pub fn evaluate(&self, points: &HashMap<usize, F>) -> F {
        self.terms.iter().map(|(exponents, coeff)| {
            let term_value = exponents.iter().map(|(&var, &exp)| {
                points.get(&var).unwrap_or(&F::ONE).pow(exp)
            }).product::<F>();
            *coeff * term_value
        }).sum()
    }

    pub fn degree(&self) -> usize {
        self.terms.keys()
            .map(|exponents| exponents.values().sum::<usize>())
            .max()
            .unwrap_or(0)
    }

    pub fn variables(&self) -> Vec<usize> {
        self.terms.keys()
            .flat_map(|exponents| exponents.keys().cloned())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect()
    }
}

impl<F: FiniteField> Add for MultivariatePolynomial<F> {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        for (exponents, coeff) in rhs.terms {
            self.insert_term(exponents, coeff);
        }
        self
    }
}

impl<F: FiniteField> Mul for MultivariatePolynomial<F> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut result = MultivariatePolynomial::new();
        for (exp1, coeff1) in &self.terms {
            for (exp2, coeff2) in &rhs.terms {
                let mut new_exp = exp1.clone();
                for (&var, &exp) in exp2 {
                    *new_exp.entry(var).or_insert(0) += exp;
                }
                result.insert_term(new_exp, *coeff1 * *coeff2);
            }
        }
        result
    }
}

// Implement From for univariate polynomials
impl<F: FiniteField, const D: usize> From<Polynomial<Monomial, F, D>> for MultivariatePolynomial<F> {
    fn from(poly: Polynomial<Monomial, F, D>) -> Self {
        let mut result = MultivariatePolynomial::new();
        for (i, &coeff) in poly.coefficients.iter().enumerate() {
            if coeff != F::ZERO {
                let mut exponents = BTreeMap::new();
                exponents.insert(0, i);
                result.insert_term(exponents, coeff);
            }
        }
        result
    }
}

// Extend Polynomial to support conversion to multivariate
impl<F: FiniteField, const D: usize> Polynomial<Monomial, F, D> {
    pub fn to_multivariate(self, variable_index: usize) -> MultivariatePolynomial<F> {
        let mut result = MultivariatePolynomial::new();
        for (i, &coeff) in self.coefficients.iter().enumerate() {
            if coeff != F::ZERO {
                let mut exponents = BTreeMap::new();
                exponents.insert(variable_index, i);
                result.insert_term(exponents, coeff);
            }
        }
        result
    }
}