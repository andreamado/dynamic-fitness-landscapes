use rand::prelude::*;
use rand_distr::Normal;

use super::linear_algebra::{SquareMatrix, Vector};

#[allow(dead_code)]
pub struct MultivariateNormal<const S: usize> {
    mean: Vector<S>,
    covariance_matrix: SquareMatrix<S>,
    l_matrix: SquareMatrix<S>
}

impl<const S: usize> MultivariateNormal<S> {
    pub fn generate(&self) -> Vector<S> {
        match self.l_matrix {
            SquareMatrix::Null => Vector::new(),
            SquareMatrix::NonNull(l_matrix) => {
                let mut rng = thread_rng();
                let normal = Normal::new(0., 1.).unwrap();

                let mut temp = [0.; S];
                for i in 0..S {
                    temp[i] = normal.sample(&mut rng);
                }

                let mut res = [0.; S];
                for i in 0..S {
                    res[i] = (0..=i).fold(self.mean[i], |acc, j| {
                        acc + l_matrix[i][j] * temp[j]
                    });
                }
                Vector::from(res)
            }
        }
    }

    pub fn new(mean: Vector<S>, covariance_matrix: SquareMatrix<S>) -> std::result::Result<Self, &'static str> {
        let l_matrix = match covariance_matrix {
            SquareMatrix::Null => SquareMatrix::Null,
            SquareMatrix::NonNull(covariance_matrix) => {
                // Cholesky–Banachiewicz algorithm
                // https://en.wikipedia.org/wiki/Cholesky_decomposition
                let mut l_matrix = [[0.; S]; S];
                for i in 0..S {
                    for j in 0..=i {
                        let sum = (0..j).fold(0., |acc, k| {
                            acc + l_matrix[i][k] * l_matrix[j][k]
                        });

                        l_matrix[i][j] = if i == j {
                            (covariance_matrix[i][i] - sum).sqrt()
                        } else {
                            (covariance_matrix[i][j] - sum) / l_matrix[j][j]
                        };

                        if l_matrix[i][j].is_nan() {
                            return Err("covariance_matrix is not positive-definite")
                        }
                    }
                }
                SquareMatrix::NonNull(l_matrix)
            }
        };
        Ok(MultivariateNormal {
            mean,
            covariance_matrix,
            l_matrix
        })
    }
}

