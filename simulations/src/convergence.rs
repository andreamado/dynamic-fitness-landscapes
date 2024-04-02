//! convergence records detailed information about fitness landscape evolution
//! over time, including plots of the fitness landscape.
//! 
//! For information on the parameters, run `convergence --help`
//! 
//! The user can convert a time series of plots to a video using ffmpeg with a 
//! command like `ffmpeg -i example/%06d.svg -vf format=yuv420p output.mp4`


pub mod modules;
use modules::{
    population::{FixedSizePopulation, InitialPopulation},
    resource_based_landscape::ResourceBasedFitnessLandscape,
    genotype::Genotype,
    data::Data,
    parameters::Parameters,
    plot_landscape::FitnessLandscapePlot
};

use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    const L: usize = 10;
    const S: usize = 2;

    let t_max = 100_000;
    let params = Parameters::<S>::from_command_line_convergence();

    let mut data = Data::from_parameters(&params, L);
    let l = params.landscapes[0];

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
        population.initialize(InitialPopulation::SingleGenotype(Genotype::<L>::random()));

        for t in 0..t_max {
            population.mutation(params.mutation_rate_per_locus);
            population.wright_fisher(&landscape, &params.resources);

            data.save_datapoint(l, 0, &population, &landscape, &params.resources, t, true).unwrap();

            let fitness_landscape = landscape.get_full_fitness_landscape(&population, &params.resources);
            let filename = format!("{}landscape_data_{:06}.dat", params.folder_name, t);
            fitness_landscape.save(&filename)?;

            let colors = population.distribution();
            let filename = format!("{}{:06}.svg", params.folder_name, t);
            let res = FitnessLandscapePlot::new(&fitness_landscape.landscape, None, Some(&colors)).plot(&filename);
            if res.is_err() {
              println!("Could not save file {}. Skipping...", filename);
            }
        }
    }
    data.flush()?;

    Ok(())
}
