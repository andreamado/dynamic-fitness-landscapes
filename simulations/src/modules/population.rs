use rand_distr::{Binomial, Bernoulli, Distribution, WeightedAliasIndex};
use rand::prelude::IteratorRandom;

use std::{
    collections::HashMap,
    ops::{Index, Deref},
    fmt,
    fs::File,
    io::{BufWriter, Write},
    error::Error
};

use super::{
    genotype::Genotype,
    resource_based_landscape::ResourceBasedFitnessLandscape,
    math::linear_algebra::Vector
};

#[allow(dead_code)]
pub enum InitialPopulation<const L: usize> {
    NeutralSFS,
    Binomial(f64), // The value is the probability of derived allele
    SingleGenotype(Genotype<L>)
}

#[derive(Clone)]
pub struct FixedSizePopulation<const L: usize> {
    population: HashMap<Genotype<L>, usize>,
    pop_size:   usize,
    binomial_coefficients: [f64; L]
}

impl<const L: usize> FixedSizePopulation<L> {
    /// Creates a new, empty population
    pub fn new(size: usize) -> Self {
        Self {
            population: HashMap::new(),
            pop_size:   size,
            binomial_coefficients: (1..=L).map(|n| {
                binomial(L, n)
            }).collect::<Vec<f64>>().try_into().unwrap()
        }
    }

    /// Adds one individual with the specified genotype
    pub fn add_individual(&mut self, genotype: Genotype<L>) {
        match self.population.get_mut(&genotype) {
            // if the new genotype is already present in the population add an individual
            Some(g) => { *g += 1 },
            // otherwise add the new genotype with one individual
            None    => { self.population.insert(genotype, 1); }
        }
    }

    pub fn add_genotype(&mut self, genotype: Genotype<L>, n: usize) {
        match self.population.get_mut(&genotype) {
            // if the new genotype is already present in the population add n individuals
            Some(g) => { *g += n },
            // otherwise add the new genotype with n individuals
            None    => { self.population.insert(genotype, n); }
        }
    }

    pub fn initialize(&mut self, initial_population: InitialPopulation<L>) {
        self.population.clear();
        match initial_population {
            InitialPopulation::SingleGenotype(genotype) => {
                self.population.insert(genotype, self.pop_size);
            },
            InitialPopulation::NeutralSFS => {
                unimplemented!();
                // Check haploid_recombination2 for a reference implementation
            },
            // Generates an initial population where each individual has a probability equal to
            // minor_allele_probability of carrying the minor allele form for each allele
            InitialPopulation::Binomial(minor_allele_probability) => {
                let mut rng = rand::thread_rng();
                let allele_type = Bernoulli::new(minor_allele_probability).unwrap();

                for _ in 0..self.pop_size {
                    let genotype: Vec<u8> = (0..L)
                        .map(|_| allele_type.sample(&mut rng) as u8)
                        .collect();
                    self.add_individual(Genotype::<L>::from_sequence(&genotype));
                }
            },
        }
    }

    /// Cleans population, retaining only the genotypes with one or more individuals
    fn clean_population(&mut self) {
        self.population.retain(|_, &mut n| n > 0)
    }

    pub fn mutation(&mut self, mutation_rate_per_locus: f64) {
        let mut rng = rand::thread_rng();

        // The probability of a genotype acquiring one or more mutations is one minus the probability
        // of not acquiring any mutation.
        let genotype_mutation_probability = 1. - (1. - mutation_rate_per_locus).powi(L as i32);

        // Distribution that checks the number of mutations
        let m = mutation_rate_per_locus;
        let weights: Vec<f64> = (1..=L).map(|n| {
            self.binomial_coefficients[n-1] * (1. - m).powi(L as i32 - n as i32) * m.powi(n as i32)
        }).collect();
        let number_of_mutations = WeightedAliasIndex::new(weights).unwrap();

        // For each genotype present in the population
        for (genotype, n) in self.population.clone().drain() {
            // Count how many individuals will carry mutations
            let bin = Binomial::new(n as u64, genotype_mutation_probability).unwrap();
            let individuals_with_mutations = bin.sample(&mut rng) as usize;

            // Remove the mutated individuals from the population
            match self.population.get_mut(&genotype) {
                Some(g) => *g -= individuals_with_mutations,
                None    => unreachable!("Genotype not found!")
            }

            // Generate new genotypes for the mutated individuals
            for _ in 0..individuals_with_mutations {
                let mut new_genotype = genotype.clone();

                // How many mutations?
                let n_mutations = number_of_mutations.sample(&mut rng) + 1;

                // which mutations?
                for i in (0..L).choose_multiple(&mut rng, n_mutations) {
                    new_genotype.mutate(i);
                }
                self.add_individual(new_genotype)
            }
        }

        // clean the genotypes that have no individuals
        self.clean_population();
    }

    pub fn to_vector(&self) -> (Vec<Genotype<L>>, Vec<usize>) {
        let mut genotypes = Vec::<Genotype<L>>::with_capacity(self.population.len());
        let mut ns = Vec::<usize>::with_capacity(self.population.len());

        for (&g, &n) in &self.population {
            genotypes.push(g);
            ns.push(n);
        }
        (genotypes, ns)
    }

    pub fn to_vec(&self) -> Vec<(Vec<u8>, usize)> {
        self.population.iter().map(|(&g, &n)| (g.to_vec(), n)).collect()
    }

    pub fn save(&self, filename: &str) -> Result<(), Box<dyn Error>> {
        let file = File::create(filename).unwrap();
        let mut file = BufWriter::new(file);

        write!(file, "#")?;
        for i in 0..L {
            write!(file, "g{}\t", i)?;
        }
        write!(file, "n\n")?;

        for (&g, &n) in self.population.iter() {
            write!(file, "{}\t{}\n", g, n)?;
        }
        Ok(())
    }

    pub fn from_vec(vec: &Vec<(Vec<u8>, usize)>) -> Self {
        let pop_size = vec.iter().fold(0, |acc, (_, n)| acc + n);
        let mut population = Self::new(pop_size);

        for (g, n) in vec {
            population.add_genotype(Genotype::from_sequence(&g[..]), *n);
        }
        population
    }


    pub fn wright_fisher<const S: usize>(&mut self, landscape: &ResourceBasedFitnessLandscape<L,S>, resources: &Vector<S>) {
        let mut rng = rand::thread_rng();

        // Get the fitnesses of the genotypes
        let fitness_landscape = landscape.get_occupied_fitness_landscape(&self, resources);
        let n_genotypes = fitness_landscape.len();

        let mut genotypes = Vec::<Genotype<L>>::with_capacity(n_genotypes);
        let mut fitnesses = Vec::<f64>::with_capacity(n_genotypes);

        for (g, f) in fitness_landscape {
            genotypes.push(g);
            fitnesses.push(f);
        }

        // Create the new population
        let mut new_population = vec![0_usize; n_genotypes];
        let new_indices = rand::distributions::WeightedIndex::new(&fitnesses).unwrap();
        for _ in 0..self.pop_size {
            new_population[new_indices.sample(&mut rng)] += 1;
        }

        self.population.clear();
        for i in 0..n_genotypes {
            if new_population[i] > 0 {
                self.population.insert(genotypes[i], new_population[i]);
            }
        }
    }

    /// Returns the number of genotypes *currently* present in the population
    #[inline]
    pub fn n_genotypes(&self) -> usize {
        self.population.len()
    }

    pub fn distribution(&self) -> HashMap<Genotype<L>, f64> {
        self.population.iter()
                       .map(|(g, n)| ((*g).clone(), (*n as f64)/(self.pop_size as f64)))
                       .collect()
    }

    /// Returns the population size
    #[inline]
    pub fn size(&self) -> usize {
        self.pop_size
    }

    /// Returns the absolute Shannon entropy of the population
    pub fn shannon_entropy(&self) -> f64 {
        let mut entropy = 0_f64;
        let size = self.pop_size as f64;
        for &n in self.population.values() {
            if n > 0 {
                let f: f64 = (n as f64) / size;
                entropy -= f * f.ln();
            }
        }
        entropy
    }

    pub fn haplotype_diversity(&self) -> f64 {
        let mut h = 1_f64;
        let size = self.pop_size as f64;
        for &n in self.population.values() {
            let f: f64 = (n as f64) / size;
            h -= f * f;
        }
        // No correction for sample size since we sample the full population
        h
    }

    pub fn nucleotide_diversity(&self) -> f64 {
        let mut pi = 0_f64;
        let size = self.pop_size as f64;
        for (&gi, &ni) in &self.population {
            let xi = ni as f64 / size;
            for (gj, &nj) in &self.population {
                let xj = nj as f64 / size;
                let kij = gi.n_differences(gj) as f64;
                pi += xi * xj * kij;
            }
        }
        // No correction for sample size since we sample the full population
        pi
    }

}

impl<const L: usize> Index<Genotype<L>> for FixedSizePopulation<L> {
    type Output = usize;
    fn index(&self, g: Genotype<L>) -> &Self::Output {
        match self.population.get(&g) {
            Some(i) => &i,
            None    => &0
        }
    }
}
impl<const L: usize> IntoIterator for FixedSizePopulation<L> {
    type Item = (Genotype<L>, usize);
    type IntoIter = std::collections::hash_map::IntoIter<Genotype<L>, usize>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.population.into_iter()
    }
}

impl<const L: usize> Deref for FixedSizePopulation<L> {
    type Target = HashMap<Genotype<L>, usize>;

    fn deref(&self) -> &Self::Target {
        &self.population
    }
}

impl<const L: usize> fmt::Display for FixedSizePopulation<L> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (&g, &n) in &self.population {
            write!(f, "{}\t{}\n", g, n)?
        };
        write!(f, "")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn general() {
        const L: usize = 5;
        let size = 100;
        let mut population = FixedSizePopulation::<L>::new(size);
        population.initialize(InitialPopulation::SingleGenotype(Genotype::<L>::new()));

        population.mutation(0.);
        assert_eq!(population[Genotype::<L>::new()], size);
        assert_eq!(population[Genotype::<L>::from_sequence(&[0, 1, 0, 1, 0])], 0);

        population.mutation(1.);
        assert_eq!(population[Genotype::<L>::new()], 0);
    }
}

/// Computes the binomial coefficient
fn binomial(n: usize, k: usize) -> f64 {
      (n-k..=n).product::<usize>() as f64
    / (1..=k  ).product::<usize>() as f64
}
