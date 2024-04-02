use super::math::linear_algebra::{SquareMatrix, Vector};

#[derive(Copy, Clone, Debug)]
pub enum FitnessModel<const S: usize> {
    HoC{
        cb: SquareMatrix<S>
    },
    Additive{
        mu: Vector<S>,
        ca: SquareMatrix<S>
    },
    RoughMountFuji {
        mu: Vector<S>,
        ca: SquareMatrix<S>,
        cb: SquareMatrix<S>
    }
}

impl<const S: usize> FitnessModel<S> {
    pub fn new_hoc(params: Vec<f64>) -> Self {
        let cb_diagonal    = params[0];
        let cb_offdiagonal = params[1];
        let mut cb = [[0_f64; S]; S];
        for i in 0..S {
            for j in 0..S {
                if i == j {
                    cb[i][i] = cb_diagonal;
                } else {
                    cb[i][j] = cb_offdiagonal;
                }
            }
        }
        FitnessModel::HoC {
            cb: SquareMatrix::from(cb)
        }
    }

    pub fn new_additive(params: Vec<f64>) -> Self {
        let mu = Vector::from([params[0]; S]);

        let ca_diagonal    = params[1];
        let ca_offdiagonal = params[2];

        let mut ca = [[0_f64; S]; S];
        for i in 0..S {
            for j in 0..S {
                if i == j {
                    ca[i][i] = ca_diagonal;
                } else {
                    ca[i][j] = ca_offdiagonal;
                }
            }
        }
        FitnessModel::Additive {
            mu, ca: SquareMatrix::from(ca)
        }
    }

    pub fn new_rmf(params: Vec<f64>) -> Self {
        let mu = Vector::from([params[0]; S]);

        let ca_diagonal    = params[1];
        let ca_offdiagonal = params[2];
        let cb_diagonal    = params[3];
        let cb_offdiagonal = params[4];

        let ca = if ca_diagonal > 0. {
            let mut ca = [[0_f64; S]; S];
            for i in 0..S {
                for j in 0..S {
                    if i == j {
                        ca[i][i] = ca_diagonal;
                    } else {
                        ca[i][j] = ca_offdiagonal;
                    }
                }
            }
            SquareMatrix::from(ca)
        } else {
            SquareMatrix::<S>::Null
        };

        let cb = if cb_diagonal > 0. {
            let mut cb = [[0_f64; S]; S];
            for i in 0..S {
                for j in 0..S {
                    if i == j {
                        cb[i][i] = cb_diagonal;
                    } else {
                        cb[i][j] = cb_offdiagonal;
                    }
                }
            }
            SquareMatrix::from(cb)
        } else {
            SquareMatrix::<S>::Null
        };

        FitnessModel::RoughMountFuji {
            mu, ca, cb
        }
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut r = Vec::<u8>::new();
        match self {
            Self::HoC {cb} => {
                r.extend(cb.to_bytes());
                r.push(0);
            },
            Self::Additive {mu, ca} => {
                r.extend(ca.to_bytes());
                r.extend(mu.to_bytes());
                r.push(1);
            }
            Self::RoughMountFuji {mu, ca, cb} => {
                r.extend(ca.to_bytes());
                r.extend(cb.to_bytes());
                r.extend(mu.to_bytes());
                r.push(2);
            }
        }
        r
    }
    pub fn from_bytes(vec: &Vec<u8>) -> Self {
        // const mat_size: usize = S*S*8+1;
        match vec.last() {
            Some(&0) => Self::HoC {
                cb: SquareMatrix::<S>::from_bytes(&vec[0..(S*S*8+1)]).unwrap()
            },
            Some(&1) => {
                Self::Additive {
                    ca: SquareMatrix::<S>::from_bytes(&vec[0..(S*S*8+1)]).unwrap(),
                    mu: Vector::<S>::from_bytes(&vec[(S*S*8+1)..(vec.len()-1)]).unwrap()
                }
            },
            Some(&2) => {
                Self::RoughMountFuji {
                    ca: SquareMatrix::<S>::from_bytes(&vec[0..(S*S*8+1)]).unwrap(),
                    cb: SquareMatrix::<S>::from_bytes(&vec[(S*S*8+1)..2*(S*S*8+1)]).unwrap(),
                    mu: Vector::<S>::from_bytes(&vec[2*(S*S*8+1)..(vec.len()-1)]).unwrap()
                }
            }
            Some(&_) => panic!("Model type not recognized"),
            None     => panic!("Could not load fitness model: empty vector")
        }
    }

    fn t(a: f64) -> f64 {
        (a * 100_000.).trunc() / 100_000.
    }

    pub fn get_name(&self) -> String {
        match self {
            Self::HoC {cb} => {
                format!(
                    "HoC_S{}_cd{:.5}_co{:.5}",
                    S, Self::t(cb[(0, 0)]), Self::t(cb[(0, 1)])
                )
            },
            Self::Additive {mu, ca} => {
                format!(
                    "additive_S{}_mu{:.5}_cd{:.5}_co{:.5}",
                    S, Self::t(mu[0]), Self::t(ca[(0, 0)]), Self::t(ca[(0, 1)])
                )
            }
            Self::RoughMountFuji {mu, ca, cb} => {
                format!(
                    "RMF_S{}_mu{:.5}_cad{:.5}_cao{:.5}_cbd{:.5}_cbo{:.5}",
                    S, Self::t(mu[0]), Self::t(ca[(0, 0)]), Self::t(ca[(0, 1)]), Self::t(cb[(0, 0)]), Self::t(cb[(0, 1)])
                )
            }
        }
    }
}

impl<const S: usize> ::std::convert::From<&str> for FitnessModel<S> {
    fn from(s: &str) -> Self {
        let s: Vec<&str> = match s.split_once("--model ") {
            None       => panic!("it was not possible to find option --model"),
            Some((_, s)) => s.split_whitespace().collect()
        };

        match &s[..] {
            ["HoC", cb_diagonal, cb_offdiagonal, ..] => {
                let parameters = vec![cb_diagonal.parse::<f64>().unwrap(), cb_offdiagonal.parse::<f64>().unwrap()];
                FitnessModel::new_hoc(parameters)
            },
            ["Additive", mu, ca_diagonal, ca_offdiagonal, ..] => {
                let parameters = vec![
                    mu.parse::<f64>().unwrap(),
                    ca_diagonal.parse::<f64>().unwrap(),
                    ca_offdiagonal.parse::<f64>().unwrap()
                ];
                FitnessModel::new_additive(parameters)
            },
            ["RoughMountFuji", mu, ca_diagonal, ca_offdiagonal, cb_diagonal, cb_offdiagonal, ..] => {
                let parameters = vec![
                    mu.parse::<f64>().unwrap(),
                    ca_diagonal.parse::<f64>().unwrap(),
                    ca_offdiagonal.parse::<f64>().unwrap(),
                    cb_diagonal.parse::<f64>().unwrap(),
                    cb_offdiagonal.parse::<f64>().unwrap()
                ];
                FitnessModel::new_rmf(parameters)
            },
            [model, ..] => panic!("Did not recognize the model: {}", model),
            [..]        => panic!("No model found")
        }
    }
}
