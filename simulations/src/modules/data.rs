use std::{
    collections::HashMap,
    error::Error,
    fs::File,
    io::{BufWriter, Write}
};

use serde::{Serialize, Deserialize};
use rand::Rng;

use super::{
    genotype::Genotype,
    math::linear_algebra::Vector,
    population::FixedSizePopulation,
    resource_based_landscape::ResourceBasedFitnessLandscape,
    fitness_landscape::VecLandscape,
    parameters::Parameters
};

const BUFFER_SIZE: usize = 5000;
const MAX_GENERATIONS: usize = 500;
const THRESHOLD: f64 = 0.1;

pub struct Data<'a, const S: usize> {
    summary: BufWriter<File>,
    parameters: &'a Parameters<S>,
    buffer: Vec<DataPoint>,
    pos: usize,
    past_top_genotypes: Vec<[i64; MAX_TOPGENOTYPES]>
}

impl<'a, const S: usize> Data<'a, S> {
    pub fn from_parameters(parameters: &'a Parameters<S>, l: usize) -> Self {
        let unique_id = rand::thread_rng().gen_range(0..10000);

        let folder_name = if parameters.folder_name.len() > 0 {
            parameters.folder_name.clone()
        } else {
            "data/".to_string()
        };

        let filename = if parameters.null_model { format!(
            "{}L{}_{}_m{:e}_r[{}]_null_{}.dat",
            folder_name,
            l, parameters.model.get_name(), parameters.mutation_rate_per_locus,
            &parameters.resources.iter().fold(String::new(), |acc, r| format!("{},{:.3}", acc, r))[1..],
            unique_id
        ) } else { format!(
            "{}L{}_{}_m{:e}_r[{}]_full_{}.dat",
            folder_name,
            l, parameters.model.get_name(), parameters.mutation_rate_per_locus,
            &parameters.resources.iter().fold(String::new(), |acc, r| format!("{},{:.3}", acc, r))[1..],
            unique_id
        ) };

        let file = File::create(filename).unwrap();
        let mut summary = BufWriter::new(file);

        summary.write(b"#n_pop\tlandscape_idx\treplicate\tt\tentropy\thaplotype_diversity\tnucleotide_diversity\tstrains\tn_maxima\tn_minima\tmaximum\tminimum\tgamma\tmean\tvar\tfitness_wildtype\tmean_phenotypic_distance").unwrap();
        for i in 0..MAX_TOPGENOTYPES {
            summary.write(format!("\ttg{}\tn{}", i, i).as_bytes()).unwrap();
        }
        summary.write(b"\n").unwrap();

        Self {
            summary,
            parameters,
            buffer: vec![DataPoint::empty(); BUFFER_SIZE],
            pos: 0,
            past_top_genotypes: vec![[-1; MAX_TOPGENOTYPES]; BUFFER_SIZE]
        }
    }

    pub fn save_landscape<const L: usize>(&self, landscape: &ResourceBasedFitnessLandscape<L,S>, l: usize) -> Result<(), Box<dyn Error>> {
        landscape.save(&self.parameters.model.get_name()[..], l)?;
        Ok(())
    }

    pub fn save_datapoint<const L: usize>(&mut self,
        l: usize,
        r: usize,
        population: &FixedSizePopulation<L>,
        landscape:  &ResourceBasedFitnessLandscape<L,S>,
        resources:  &Vector<S>,
        t: usize,
        write_to_file: bool
    ) -> Result<(), Box<dyn Error>> {
            self.buffer[self.pos] = DataPoint::new(&population, &landscape, &resources, l, r, t);
            self.past_top_genotypes[self.pos] = self.top_genotypes();
            self.pos = (self.pos + 1) % BUFFER_SIZE;

            if write_to_file {
                self.buffer[self.pos].save(&mut self.summary)?;
            }
            Ok(())
        }

    pub fn write_to_file(&mut self) -> Result<(), Box<dyn Error>> {
        let beg = self.pos - 1 - MAX_GENERATIONS + BUFFER_SIZE;
        for i in 0..MAX_GENERATIONS {
            self.buffer[(beg + i) % BUFFER_SIZE].save(&mut self.summary)?;
        }
        Ok(())
    }

    pub fn top_genotypes(&self) -> [i64; MAX_TOPGENOTYPES] {
        let mut tg = HashMap::<i64, usize>::with_capacity(MAX_TOPGENOTYPES);
        for i in 0..MAX_GENERATIONS {
            for &g in &self.buffer[(self.pos - i + BUFFER_SIZE) % BUFFER_SIZE].top_genotypes {
                if g == -1 { break }
                let count = tg.entry(g).or_insert(0);
                *count += 1;
            }
        }
        tg.retain(|_, n| (*n as f64) / (MAX_GENERATIONS as f64) > THRESHOLD);

        let mut arr: Vec<i64> = tg.into_keys().collect();
        arr.sort_unstable();
        arr.resize(MAX_TOPGENOTYPES, -1);
        arr.try_into().unwrap()
    }

    pub fn flush(&mut self) -> Result<(), Box<dyn Error>> {
        self.summary.flush()?;
        Ok(())
    }

    pub fn stable_state(&self) -> bool {
        let tg1 = self.past_top_genotypes[(self.pos - 1 + BUFFER_SIZE) % BUFFER_SIZE];
        for i in 1..MAX_GENERATIONS {
            let tg2 = self.past_top_genotypes[(self.pos - i + BUFFER_SIZE) % BUFFER_SIZE];
            if tg1 != tg2 {
                return false
            }
        }
        true
    }
}

impl<'a, const S: usize> Drop for Data<'a, S> {
    fn drop(&mut self) {
        self.flush().unwrap();
    }
}

const MAX_TOPGENOTYPES: usize = 10;

#[derive(Serialize, Deserialize, Clone)]
pub struct DataPoint {
    size: usize,
    l: usize,
    r: usize,
    t: usize,
    entropy:  f64,
    nucleotide_diversity: f64,
    haplotype_diversity: f64,
    strains:  usize,
    n_maxima: usize,
    n_minima: usize,
    maximum_minimum: [f64; 2],
    top_genotypes:   [i64;   MAX_TOPGENOTYPES],
    n_top_genotypes: [usize; MAX_TOPGENOTYPES],
    gamma: f64,
    mean:  f64,
    var:   f64,
    fitness_wildtype: f64,
    mean_phenotypic_distance: f64,
    landscape: Option<VecLandscape>
}

impl DataPoint {
    pub fn new<const S: usize, const L: usize>(
        population: &FixedSizePopulation<L>,
        landscape:  &ResourceBasedFitnessLandscape<L,S>,
        resources:  &Vector<S>,
        l: usize,
        r: usize,
        t: usize
    ) -> Self {
            let fitness_landscape = landscape.get_full_fitness_landscape(&population, &resources);
            let (_, &max) = fitness_landscape.max().unwrap_or((&Genotype::new(), &f64::NAN));
            let (_, &min) = fitness_landscape.min().unwrap_or((&Genotype::new(), &f64::NAN));
            let (mean, var) = fitness_landscape.mean_var();

            let fitness_wildtype = *fitness_landscape.get(&Genotype::new()).unwrap_or(&f64::NAN);

            let mean_phenotypic_distance = landscape.mean_phenotypic_distance(population);

            let mut top_genotypes   = [-1; MAX_TOPGENOTYPES];
            let mut n_top_genotypes = [ 0; MAX_TOPGENOTYPES];
            let mut k = 0;

            for (g, n) in population.iter() {
                if let Some(f) = fitness_landscape.get(g) {
                    if *f > 1. {
                        top_genotypes[k] = g.index() as i64;
                        n_top_genotypes[k] = *n;
                        k += 1;
                    }
                }
                if k >= MAX_TOPGENOTYPES { break; }
            }

            Self {
                size: population.size(),
                l, r, t,
                entropy: population.shannon_entropy(),
                haplotype_diversity: population.haplotype_diversity(),
                nucleotide_diversity: population.nucleotide_diversity(),
                strains: population.n_genotypes(),
                n_maxima: fitness_landscape.maxima().len(),
                n_minima: fitness_landscape.minima().len(),
                maximum_minimum: [max, min],
                top_genotypes,
                n_top_genotypes,
                gamma: fitness_landscape.gamma(),
                mean:  mean,
                var:   var,
                fitness_wildtype,
                mean_phenotypic_distance,
                landscape: None
            }
        }

    pub fn empty() -> Self {
        Self {
            size: 0,
            l: 0,
            r: 0,
            t: 0,
            entropy: f64::NAN,
            haplotype_diversity: f64::NAN,
            nucleotide_diversity: f64::NAN,
            strains: 0,
            n_maxima: 0,
            n_minima: 0,
            maximum_minimum: [f64::NAN, f64::NAN],
            top_genotypes: [-1; MAX_TOPGENOTYPES],
            n_top_genotypes: [0; MAX_TOPGENOTYPES],
            gamma: f64::NAN,
            mean:  f64::NAN,
            var:   f64::NAN,
            fitness_wildtype: f64::NAN,
            mean_phenotypic_distance: f64::NAN,
            landscape: None
        }
    }

    pub fn get(&self, property: &str) -> f64 {
        match property {
            "entropy"  => self.entropy,
            "strains"  => self.strains as f64,
            "n_maxima" => self.n_maxima as f64,
            "n_minima" => self.n_minima as f64,
            "maximum"  => self.maximum_minimum[0],
            "minimum"  => self.maximum_minimum[1],
            "gamma"    => self.gamma,
            "mean"     => self.mean,
            "var"      => self.var,
            "fitness_wildtype" => self.fitness_wildtype,
            "phenotypic_distance" => self.mean_phenotypic_distance,
            _ => panic!("unrecognized property: {}", property)
        }
    }

    pub fn save(&self, file: &mut BufWriter<File>) -> Result<(), Box<dyn Error>> {
        let mut top_genotypes = Vec::<String>::with_capacity(MAX_TOPGENOTYPES);
        for i in 0..MAX_TOPGENOTYPES {
            top_genotypes.push( self.top_genotypes[i].to_string()   );
            top_genotypes.push( self.n_top_genotypes[i].to_string() );
        }

        file.write(
            format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                self.size, self.l, self.r, self.t,
                self.entropy, self.haplotype_diversity,
                self.nucleotide_diversity, self.strains,
                self.n_maxima, self.n_minima,
                self.maximum_minimum[0], self.maximum_minimum[1],
                self.gamma,
                self.mean, self.var,
                self.fitness_wildtype,
                self.mean_phenotypic_distance,
                top_genotypes.join("\t")
            ).as_bytes()
        )?;
        Ok(())
    }
}
