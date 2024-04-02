//! ecoevo_landscapes simulates a population evolving on the fitness landscape
//! and records statistical information about the population and fitness 
//! landscape
//! 
//! For information on the parameters, run `ecoevo_landscape --help`

pub mod modules;
use modules::{
    population::{
        FixedSizePopulation, 
        InitialPopulation
    },
    resource_based_landscape::ResourceBasedFitnessLandscape,
    genotype::Genotype,
    data::Data,
    parameters::Parameters
};

use std::time::Instant;

fn main() {
    const L: usize = 10;
    const S: usize = 2;

    let t_max = 100_000;
    let t_min = 15_000;
    let params = Parameters::<S>::from_command_line();

    let mut data = Data::from_parameters(&params, L);

    let mut output = String::new();
    output.push_str(&format!("#{}\t{} model\n", params.model.get_name(), if params.null_model {"null"} else {"full"}));
    output.push_str(&format!("#landscape_id\tpop_size\treplicate\ttime(s)\n"));

    for l in params.landscapes[0]..params.landscapes[1] {
        let landscape = {
            let landscape_filename = format!(
                "landscapes/L{}_{}_{}.dat",
                L, params.model.get_name(), l
            );
            let mut landscape = ResourceBasedFitnessLandscape::<L, S>::load(&landscape_filename[..]);
            if params.null_model { landscape.as_null_model(); }
            landscape
        };

        for &pop_size in &params.pop_size {
            let mut population = FixedSizePopulation::<L>::new(pop_size);
            for r in 0..params.replicates {
                let start = Instant::now();

                population.initialize(InitialPopulation::SingleGenotype(Genotype::<L>::random()));
                for t in 0..t_max {
                    population.mutation(params.mutation_rate_per_locus);
                    population.wright_fisher(&landscape, &params.resources);

                    if t > t_min - 501 {
                        let _ = data.save_datapoint(l, r, &population, &landscape, &params.resources, t, false);
                    }
                    if t > t_min && data.stable_state() {
                        break
                    }
                }
                let _ = data.write_to_file();
                output.push_str(&format!("{}\t{}\t{}\t{:.3}\n", l, pop_size, r, start.elapsed().as_secs_f32()));
            }
        }
        data.flush().unwrap();
    }
    println!("{}\n", output);
}
