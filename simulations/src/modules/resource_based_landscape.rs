use super::{
    multidimensional_rough_mount_fuji::{MultidimensionalRoughMountFuji, VecRMF},
    population::FixedSizePopulation,
    genotype::{Genotype, possible_sequences},
    fitness_landscape::{FitnessLandscape, FitnessType},
    fitness_model::FitnessModel,
    math::linear_algebra::Vector
};

use std::{
    collections::HashMap,
    fs::File,
    error::Error
};

#[derive(Clone)]
pub struct ResourceBasedFitnessLandscape<const L: usize, const S: usize> {
    phenotypic_landscape: MultidimensionalRoughMountFuji<L, S>,
    null_model: bool
}

impl<const L: usize, const S: usize> ResourceBasedFitnessLandscape<L, S> {
    pub fn new(fitness_model: FitnessModel<S>) -> Self {
        ResourceBasedFitnessLandscape {
            phenotypic_landscape: MultidimensionalRoughMountFuji::<L,S>::new(fitness_model),
            null_model: false
        }
    }

    pub fn save(&self, model_name: &str, l: usize) -> Result<(), Box<dyn Error>> {
        let filename = format!(
            "landscapes/L{}_{}_{}.dat",
            L, model_name, l
        );
        let file = File::create(filename).unwrap();
        serde_cbor::to_writer(file, &self.to_vec())?;
        Ok(())
    }

    pub fn get_occupied_fitness_landscape(&self, population: &FixedSizePopulation<L>, resources: &Vector<S>) -> HashMap<Genotype<L>,f64> {
        let mut fitness_landscape = HashMap::<Genotype<L>, f64>::with_capacity(population.n_genotypes());
        if self.null_model {
            // Warning: this fitness is not correctly normalized!
            let (genotypes, n) = population.to_vector();

            for (i, &g) in genotypes.iter().enumerate() {
                let fitness: f64 = self.phenotypic_landscape
                                  .get_multiplicative(g)
                                  .iter().sum();
                fitness_landscape.insert(g, n[i] as f64 * fitness);
            }
        } else {
            let (genotypes, n) = population.to_vector();
            let n_genotypes = population.n_genotypes();
            let pop_size = population.size() as f64;

            let mut alpha = vec![[0_f64; S]; n_genotypes];
            let mut n_alpha = [0_f64; S];

            for (i, &g) in genotypes.iter().enumerate() {
                let phenotype = self.phenotypic_landscape.get_multiplicative(g);

                for j in 0..S {
                    alpha[i][j] = phenotype[j];
                    n_alpha[j] += (n[i] as f64) * phenotype[j]
                }
            }

            let resources: Vec<f64> = (0..S).map(|j| resources[j] / n_alpha[j]).collect();
            let fitness = (0..n_genotypes).map(|i| {
                (0..S).fold(0_f64, |acc, j| {
                    acc + alpha[i][j] * resources[j] * (n[i] as f64) / pop_size
                })
            });
            let fitness_iter = fitness.into_iter();
            let total: f64 = fitness_iter.clone().sum();

            for (i, f) in fitness_iter.enumerate() {
                fitness_landscape.insert(genotypes[i], f/total);
            }
        }
        fitness_landscape
    }

    pub fn mean_phenotypic_distance(&self, population: &FixedSizePopulation<L>) -> f64 {
        let mut mean_distance = 0f64;
        for (&g1, &n1) in population.iter() {
            let p1 = self.phenotypic_landscape.get_multiplicative(g1);

            for (&g2, &n2) in population.iter() {
                let p2 = self.phenotypic_landscape.get_multiplicative(g2);

                let mut dist = 0.;
                for r in 0..S {
                    let dif = p1[r] - p2[r];
                    dist += dif * dif;
                }

                mean_distance += dist.sqrt() * (n1 * n2) as f64;
            }
        }
        mean_distance / ((population.size() * (population.size() - 1)) as f64)
    }

    pub fn as_null_model(&mut self) {
        self.null_model = true;
    }

    pub fn get_full_fitness_landscape(&self, population: &FixedSizePopulation<L>, resources: &Vector<S>) -> FitnessLandscape<L> {
        let mut fitness_landscape = FitnessLandscape::<L>::new(FitnessType::Multiplicative);

        if self.null_model {
            for &g in &possible_sequences::<L>() {
                let g = Genotype::from_sequence(&g);
                let fitness: f64 = self.phenotypic_landscape
                                  .get_multiplicative(g)
                                  .iter().sum();
                fitness_landscape.add_genotype(g, fitness / S as f64);
            }

            let mean_fitness = population.iter().map(|(g, &n)| {
                                   (n as f64) * fitness_landscape.get(g).unwrap()
                               }).sum::<f64>() / population.size() as f64;

            fitness_landscape.normalize(mean_fitness);
        } else {
            let sum_r: Vec<f64> = (0..S).map(|r| {
                population.iter().map(|(&g, &n)| {
                    let ar = self.phenotypic_landscape.get_multiplicative(g)[r];
                    (n as f64) * ar
                }).sum()
            }).collect();

            let mean_fitness = resources.iter().sum::<f64>() / population.size() as f64;
            for &g in &possible_sequences::<L>() {
                let g = Genotype::from_sequence(&g);

                let a = self.phenotypic_landscape.get_multiplicative(g);
                let fitness: f64 = (0..S).map(|r| {
                    a[r] * resources[r] / sum_r[r]
                }).sum();

                fitness_landscape.add_genotype(g, fitness / mean_fitness);
            }
        }
        fitness_landscape
    }

    pub fn to_vec(&self) -> VecRMF {
        self.phenotypic_landscape.to_vec()
    }

    pub fn from_vec(vec: &VecRMF) -> Self {
        Self {
            phenotypic_landscape: MultidimensionalRoughMountFuji::<L, S>::from_vec(&vec.v),
            null_model: false
        }
    }

    pub fn load(filename: &str) -> Self {
        let loaded = if let Ok(file) = File::open(filename) {
            serde_cbor::from_reader(file).unwrap()
        } else {
            panic!("Could not open {}", filename)
        };

        Self::from_vec(&loaded)
    }
}
