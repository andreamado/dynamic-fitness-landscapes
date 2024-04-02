use std::ops::Index;
use std::collections::HashMap;

use super::{
    genotype::{Genotype, possible_sequences},
    fitness_model::FitnessModel,
    math::{
        multivariate_normal::MultivariateNormal,
        linear_algebra::Vector
    }
};

use serde::{Serialize, Deserialize};

#[derive(Clone)]
pub struct MultidimensionalRoughMountFuji<const L: usize, const S: usize> {
    phenotype: HashMap<Genotype<L>, Vector<S>>,
    fitness_model: FitnessModel<S>
}

impl<const L: usize, const S: usize> MultidimensionalRoughMountFuji<L, S> {
    pub fn new(fitness_model: FitnessModel<S>) -> Self {
        let mut phenotype = HashMap::new();

        match fitness_model {
            FitnessModel::HoC { cb } => {
                let mvn_b = MultivariateNormal::new(Vector::new(), cb).unwrap();
                for seq in possible_sequences::<L>() {
                    let g = Genotype::<L>::from_sequence(&seq);
                    phenotype.insert(g, mvn_b.generate());
                }
            },
            FitnessModel::Additive { mu, ca } => {
                let mvn_a = MultivariateNormal::new(mu, ca).unwrap();
                let additive_component: Vec::<Vector<S>> = (0..L).map(|_| mvn_a.generate()).collect();
                for seq in possible_sequences::<L>() {
                    let g = Genotype::<L>::from_sequence(&seq);

                    let mut p = Vector::new();
                    for i in 0..L {
                        let gi = g[i] as f64;
                        let ai = additive_component[i];
                        for r in 0..S {
                            p[r] += ai[r] * gi;
                        }
                    }
                    phenotype.insert(g, p);
                }
            },
            FitnessModel::RoughMountFuji { mu, ca, cb } => {
                let mvn_a = MultivariateNormal::new(mu,            ca).unwrap();
                let mvn_b = MultivariateNormal::new(Vector::new(), cb).unwrap();

                let additive_component: Vec::<Vector<S>> = (0..L).map(|_| mvn_a.generate()).collect();
                for seq in possible_sequences::<L>() {
                    let g = Genotype::<L>::from_sequence(&seq);

                    let mut p = mvn_b.generate();
                    for i in 0..L {
                        let gi = g[i] as f64;
                        let ai = additive_component[i];
                        for r in 0..S {
                            p[r] += ai[r] * gi;
                        }
                    }
                    phenotype.insert(g, p);
                }
            }
        }

        MultidimensionalRoughMountFuji {
            phenotype, fitness_model
        }
    }

    #[inline]
    pub fn get_multiplicative(&self, g: Genotype<L>) -> Vector<S> {
        let mut phenotype = self[g];
        for p in phenotype.iter_mut() {
            *p = p.exp();
        }
        phenotype
    }
    pub fn to_vec(&self) -> VecRMF {
        VecRMF {
            v: (self.phenotype.iter().map(|(g, p)| (g.to_vec(), p.to_vec())).collect(), self.fitness_model.to_bytes())
        }
    }
    pub fn from_vec(vec: &(Vec<(Vec<u8>, Vec<f64>)>, Vec<u8>)) -> Self {
        let mut phenotype = HashMap::<Genotype<L>, Vector<S>>::with_capacity(vec.0.len());
        for (g, p) in &vec.0 {
            phenotype.insert(Genotype::<L>::from_sequence(&g[..]), Vector::<S>::from_vec(&p));
        }
        Self {
            phenotype,
            fitness_model: FitnessModel::<S>::from_bytes(&vec.1)
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct VecRMF {
    pub v: (Vec<(Vec<u8>, Vec<f64>)>, Vec<u8>)
}

impl<const L: usize, const S: usize> Index<Genotype<L>> for MultidimensionalRoughMountFuji<L, S> {
    type Output = Vector<S>;
    #[inline]
    fn index(&self, g: Genotype<L>) -> &Self::Output {
        match self.phenotype.get(&g) {
            Some(p) => &p,
            None    => panic!("phenotype value not found for the required genotype!")
        }
    }
}
