
use std::{
    collections::HashMap,
    fmt,
    fs::File,
    io::{BufWriter, Write},
    cmp::Ordering,
    error::Error
};

use super::genotype::Genotype;

pub type VecLandscape = Vec<(Vec<u8>, f64)>;

#[derive(Eq, PartialEq, Copy, Clone)]
pub enum FitnessType {
    Multiplicative,
    Additive
}

pub struct FitnessLandscape<const L: usize> {
    pub landscape: HashMap<Genotype<L>,f64>,
    tp: FitnessType
}

impl<const L: usize> FitnessLandscape<L> {
    #[inline]
    pub fn new(tp: FitnessType) -> Self {
        FitnessLandscape {
            landscape: HashMap::<Genotype<L>,f64>::new(),
            tp
        }
    }

    /// Adds a genotype fitness pair
    #[inline]
    pub fn add_genotype(&mut self, g: Genotype<L>, f: f64) {
        self.landscape.insert(g, f);
    }

    /// Returns the fitness of genotype g using fitness type tp
    #[inline]
    pub fn get_fitness(&self, g: &Genotype<L>, tp: FitnessType) -> Option<f64> {
        match self.landscape.get(g) {
            None => None,
            Some(&fitness) => {
                if self.tp == tp {
                    Some(fitness)
                } else {
                    match self.tp {
                        FitnessType::Multiplicative => Some(fitness.ln()),
                        FitnessType::Additive       => Some(fitness.exp())
                    }
                }
            }
        }
    }

    /// Returns the fitness of genotype g in the current form
    #[inline]
    pub fn get(&self, g: &Genotype<L>) -> Option<&f64> {
        self.landscape.get(g)
    }

    pub fn normalize(&mut self, norm: f64) {
        self.landscape = self.landscape.iter().map(|(&g, &f)| {
            (g, f / norm)
        }).collect();
    }

    /// Returns the gamma statistics of epistasis
    pub fn gamma(&self) -> f64 {
        let (mut cov, mut var) = (0., 0.);
        for g in self.landscape.keys() {
            for j in 0..L {
                let sj = match self.get_fitness_effect(&g, j, FitnessType::Additive) {
                    Some(s) => s,
                    None    => continue
                };
                for i in 0..L {
                    if i == j { continue };
                    let sij = match self.get_fitness_effect(&g.cmutate(i), j, FitnessType::Additive) {
                        Some(s) => s,
                        None    => continue
                    };
                    cov += sj * sij;
                    var += sj * sj;
                }
            }
        }
        cov / var
    }

    /// Returns the maximum fitness in the landscape
    pub fn max(&self) -> Option<(&Genotype<L>, &f64)> {
        self.landscape.iter().reduce(|(g_a, f_a), (g_b, f_b)| {
            if f_a > f_b {
                (g_a, f_a)
            } else {
                (g_b, f_b)
            }
        })
    }
    /// Returns the minimum fitness in the landscape
    pub fn min(&self) -> Option<(&Genotype<L>, &f64)> {
        self.landscape.iter().reduce(|(g_a, f_a), (g_b, f_b)| {
            if *f_a < *f_b {
                (g_a, f_a)
            } else {
                (g_b, f_b)
            }
        })
    }
    // Returns a tuple with the mean and the variance of the fitness landscape
    pub fn mean_var(&self) -> (f64, f64) {
        let size = self.landscape.len() as f64;
        let mean = self.landscape.iter().fold(0_f64, |acc, (_, &f)| acc + f)   / size;
        let var  = self.landscape.iter().fold(0_f64, |acc, (_, &f)| acc + f*f) / size - mean * mean;
        (mean, var)
    }

    // Returns a vector with the additive fitness effects
    pub fn fitness_effects(&self, tp: FitnessType) -> Vec<f64> {
        let mut fitness_effects = Vec::<f64>::with_capacity(self.landscape.len() * L);
        for g in self.landscape.keys() {
            for i in 0..L {
                match self.get_fitness_effect(&g, i, tp) {
                    Some(s) => fitness_effects.push(s),
                    None    => continue
                }
            }
        }
        fitness_effects
    }

    pub fn get_fitness_effect(&self, g: &Genotype<L>, i: usize, tp: FitnessType) -> Option<f64> {
        let f  = self.get_fitness(&g, tp)?;
        let fi = self.get_fitness(&g.cmutate(i), tp)?;
        Some(match tp {
            FitnessType::Additive       => fi - f,
            FitnessType::Multiplicative => fi / f
        })
    }

    pub fn rank_order(&self) -> Vec<Genotype<L>> {
        let mut vec = self.to_vec_raw();
        vec.sort_by(|(_, f1), (_, f2)|
            if f1 > f2 { Ordering::Greater } else { Ordering::Less }
        );
        vec.iter().map(|(g, _)| *g).collect()
    }

    // Requires all fitnesses to be distinct
    pub fn spearman_rho(&self, other: &FitnessLandscape<L>) -> f64 {
        let r1 = self.rank_order();
        let r2 = other.rank_order();
        let n = r1.len();

        let mut d2 = 0;
        for i in 0..n {
            let j = r2.iter().position(|&g| g == r1[i]).unwrap();
            d2 += (i - j) * (i - j);
        }
        1. - 6. / ((n * (n*n - 1)) as f64) * (d2 as f64)
    }

    /// Returns a vector listing all local maxima genotypes in the landscape
    pub fn maxima(&self) -> Vec<Genotype<L>> {
        self.landscape.iter().filter_map(|(&g, &f)| {
            // Checks if fitness is larger than all neighbors and returns the genotype if so
            if (0..L).fold(true, |acc, i| {
                match self.get(&g.cmutate(i)) {
                    Some(&fi) => acc && (f > fi),
                    // If neighboring genotype doesn't exist, ignore it
                    None      => acc
                }
            }) { Some(g) } else { None }
      }).collect()
    }
    /// Returns a vector listing all local minima genotypes in the landscape
    pub fn minima(&self) -> Vec<Genotype<L>> {
        self.landscape.iter().filter_map(|(&g, &f)| {
            // Checks if fitness is larger than all neighbors and returns the genotype if so
            if (0..L).fold(true, |acc, i| {
                match self.get(&g.cmutate(i)) {
                    Some(&fi) => acc && (f < fi),
                    // If neighboring genotype doesn't exist, ignore it
                    None      => acc
                }
            }) { Some(g) } else { None }
      }).collect()
    }
    pub fn strains_selected(&self) -> Vec<Genotype<L>> {
        self.landscape.iter().filter_map(|(&g, &f)| {
            if f > 1. { Some(g) } else { None }
        }).collect()
    }

    pub fn to_vec(&self) -> VecLandscape {
        let mut v = Vec::with_capacity(self.landscape.len());
        for &g in self.landscape.keys() {
            let f = self.get_fitness(&g, FitnessType::Additive).unwrap_or(f64::NAN);
            let g = g.to_vec();
            v.push((g, f));
        }
        v
    }
    fn to_vec_raw(&self) -> Vec<(Genotype<L>, f64)> {
        let mut v = Vec::with_capacity(self.landscape.len());
        for &g in self.landscape.keys() {
            let f = self.get(&g).unwrap_or(&f64::NAN);
            v.push((g, *f));
        }
        v
    }
    pub fn from_vec(&self, v: VecLandscape) -> Self {
        let mut landscape = Self::new(FitnessType::Additive);
        for (g, f) in &v {
            landscape.add_genotype(Genotype::<L>::from_sequence(g), *f);
        }
        landscape
    }

    pub fn save(&self, filename: &str) -> Result<(), Box<dyn Error>> {
        let file = File::create(filename).unwrap();
        write!(BufWriter::new(file), "{}", self)?;
        Ok(())
    }
}

impl<const L: usize> fmt::Display for FitnessLandscape<L> {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (&g, &fitness) in &self.landscape {
            write!(f, "{}\t{}\n", g, fitness)?
        }
        write!(f, "")
    }
}
