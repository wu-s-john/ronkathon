use rand::Rng;

pub trait Random {
  fn random<R: Rng + ?Sized>(rng: &mut R) -> Self;
}

pub trait RandomOracle: Random {
  fn random_oracle<R: Rng + ?Sized>(rng: &mut R, input: &[u8]) -> Self;
}
