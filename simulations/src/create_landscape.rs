//! create_landscapes creates landscapes landscapes that will be used by
//! convergence and ecoevo_landscapes programs
//! 
//! For information on the parameters, run `create_landscape --help`


pub mod modules;
use modules::{
    resource_based_landscape::ResourceBasedFitnessLandscape,
    parameters::Parameters
};

use std::path::Path;

fn main() {
    const L: usize = 10;
    const S: usize = 2;

    let params = Parameters::<S>::from_command_line_landscape();

    for l in 0..params.landscapes[0] {
        let landscape_filename = format!(
            "landscapes/L{}_{}_{}.dat",
            L, params.model.get_name(), l
        );

        if Path::new(&landscape_filename[..]).exists() {
            println!("{} already exists. Skipping...", landscape_filename)
        } else {
            let res = ResourceBasedFitnessLandscape::<L, S>::new(params.model).save(&params.model.get_name()[..], l);
            if res.is_err() {
                println!("Could not save file {}", landscape_filename);
            }
        }
    }
}
