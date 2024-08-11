use std::collections::{HashMap, BTreeMap};
use std::hash::Hash;
use std::ops::{Add, Mul};
use itertools::Itertools;
use nalgebra as na;

use na::DMatrix;
use crate::algebra::field::FiniteField;
use super::*;
use super::{Monomial, Polynomial};

use std::ops::Sub;

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MultivariateVariable {
    pub index: usize,
    pub exponent: usize,
}

impl MultivariateVariable {
    pub fn new(index: usize, exponent: usize) -> Self {
        MultivariateVariable { index, exponent }
    }
}

#[derive(PartialEq, Eq)]
pub struct MultivariateTerm<F: FiniteField> {
    pub variables: Vec<MultivariateVariable>,
    pub coefficient: F,
}

impl<F: FiniteField> MultivariateTerm<F> {
    pub fn new(variables: Vec<MultivariateVariable>, coefficient: F) -> Self {
        MultivariateTerm { variables, coefficient }
    }
}

#[derive(PartialEq, Eq, Debug)]
pub struct MultivariatePolynomial<F: FiniteField> {
    terms: HashMap<BTreeMap<usize, usize>, F>,
}


impl<F: FiniteField + Clone> Clone for MultivariatePolynomial<F> {
    fn clone(&self) -> Self {
        MultivariatePolynomial {
            terms: self.terms.clone(),
        }
    }
}

impl<F: FiniteField> MultivariatePolynomial<F> {
    pub fn new() -> Self {
        Self { terms: HashMap::new() }
    }

    pub fn from_terms(terms: Vec<MultivariateTerm<F>>) -> Self {
        let mut poly = MultivariatePolynomial::new();
        for term in terms {
            let mut btree_map = BTreeMap::new();
            for var in term.variables {
                btree_map.insert(var.index, var.exponent);
            }
            poly.insert_term(btree_map, term.coefficient);
        }
        poly
    }

    fn insert_term(&mut self, exponents: BTreeMap<usize, usize>, coefficient: F) {
        if coefficient != F::ZERO {
            let entry = self.terms.entry(exponents.clone()).or_insert(F::ZERO);
            *entry += coefficient;
            if *entry == F::ZERO {
                self.terms.remove(&exponents);
            }
        }
    }

    pub fn coefficient(&self, exponents: &BTreeMap<usize, usize>) -> Option<&F> {
        self.terms.get(exponents)
    }

    pub fn evaluate(&self, points: &[(usize, F)]) -> F {
        self.terms.iter().map(|(exponents, coeff)| {
            let term_value = exponents.iter().map(|(&var, &exp)| {
                points.iter()
                    .find(|&&(v, _)| v == var)
                    .map(|&(_, value)| value.pow(exp))
                    .unwrap_or(F::ONE)
            }).product::<F>();
            *coeff * term_value
        }).sum()
    }

    pub fn apply_variables(&self, variables: &[(usize, F)]) -> Self {
        let mut result = MultivariatePolynomial::new();
        
        for (exponents, coeff) in &self.terms {
            let mut new_exponents = exponents.clone();
            let mut new_coeff = *coeff;
            
            for &(var, value) in variables {
                if let Some(exp) = new_exponents.get(&var) {
                    new_coeff *= value.pow(*exp);
                    new_exponents.remove(&var);
                }
            }
            
            if !new_exponents.is_empty() {
                result.insert_term(new_exponents, new_coeff);
            } else {
                result.insert_term(BTreeMap::new(), new_coeff);
            }
        }
        
        result
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

impl<F: FiniteField> MultivariatePolynomial<F> {
    pub fn mul_matrix(self, matrix: &DMatrix<F>) -> MultivariatePolynomial<F> {

        let mut poly = MultivariatePolynomial::new();
        for i in 0..matrix.nrows() {
            
            for j in 0..matrix.ncols() {
                if let Some(coefficient) = matrix.get((i, j)) {
                    for (exp, coef) in &self.terms {
                        let mut new_exp = exp.clone();
                        new_exp.insert(j, *new_exp.get(&j).unwrap_or(&0) + 1);
                        poly.insert_term(new_exp, coefficient.clone() * coef.clone());
                    }
                }
            }
        }

        return poly;
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

impl<F: FiniteField> Sub for MultivariatePolynomial<F> {
    type Output = Self;

    fn sub(mut self, rhs: Self) -> Self::Output {
        for (exponents, coeff) in rhs.terms {
            // Negate the coefficient and insert
            self.insert_term(exponents, -coeff);
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

impl<F: FiniteField + Display> Display for MultivariatePolynomial<F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let mut first = true;
        for (exponents, coeff) in self.terms.iter().sorted_by(|(a_exp, _), (b_exp, _)| {
            a_exp.iter()
                .zip(b_exp.iter())
                .find(|((a_var, a_pow), (b_var, b_pow))| {
                    a_var.cmp(b_var).then_with(|| b_pow.cmp(a_pow)).is_ne()
                })
                .map_or(std::cmp::Ordering::Equal, |(_, _)| std::cmp::Ordering::Less)
        }) {
            if !first {
                write!(f, " + ")?;
            }
            first = false;

            if *coeff != F::ONE || exponents.is_empty() {
                write!(f, "{}", coeff)?;
            }

            let mut first_var = true;
            for (&var, &exp) in exponents {
                if exp > 0 {
                    if !first_var || *coeff != F::ONE {
                        write!(f, "*")?;
                    }
                    write!(f, "x_{}", var)?;
                    if exp > 1 {
                        write!(f, "^{}", exp)?;
                    }
                    first_var = false;
                }
            }
        }

        if first {
            write!(f, "0")?;
        }

        Ok(())
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
