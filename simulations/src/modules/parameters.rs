use clap::{Arg, App, AppSettings, ArgGroup, values_t, value_t};

use super::{
    fitness_model::FitnessModel,
    math::linear_algebra::Vector
};


pub struct Parameters<const S: usize> {
    pub pop_size: Vec<usize>,
    pub mutation_rate_per_locus: f64,
    pub model: FitnessModel<S>,
    pub replicates: usize,
    pub resources: Vector<S>,
    pub landscapes: [usize; 2],
    pub null_model: bool,
    pub load_landscape: bool,
    pub folder_name: String
}

impl<const S: usize> Parameters<S> {
    pub fn from_command_line() -> Self {
        let rn: Vec<String> = (0..S).map(|i| format!("res {}", i+1)).collect();
        let resource_names: Vec<&str> = rn.iter().map(|s| s.as_str()).collect();

        let matches = App::new("")
              .author("André Amado <andre.amado@pm.me>")
              .setting(AppSettings::AllowNegativeNumbers)
              // General arguments
              .arg(Arg::with_name("population_size").help("List of population sizes").short("s").long("size").takes_value(true).multiple(true).required(true))
              .arg(Arg::with_name("mutation_rate_per_locus").help("Mutation rate per locus per generation").short("m").long("mutation_rate").value_name("rate").takes_value(true).required(true))
              .arg(Arg::with_name("resources").help("Amount of each resource").short("r").long("resources").takes_value(true).value_names(&resource_names[..]).required(true))

              .arg(Arg::with_name("landscapes").long("landscapes").short("l").value_names(&["first_landscape", "last_landscape"]).help("Range of landscapes to analize").required(true))
              .arg(Arg::with_name("replicates").long("replicates").takes_value(true).help("Number of replicates per landscapes").required(true))
              .arg(Arg::with_name("load_landscape").long("load").help("Flag loading existing landscape"))

              // Models
              .arg(Arg::with_name("HoC").long("hoc").help("House of Cards model").takes_value(true).value_names(&["cb_diag", "cb_offdiag"]))
              .arg(Arg::with_name("additive").long("add").help("Additive model").takes_value(true).value_names(&["mu", "ca_diag", "ca_offdiag"]))
              .arg(Arg::with_name("RMF").long("rmf").help("Rough Mount Fuji model").takes_value(true).value_names(&["mu", "ca_diag", "ca_offdiag", "cb_diag", "cb_offdiag"]))
              .group(ArgGroup::with_name("model").args(&["HoC", "additive", "RMF"]).required(true))

              .arg(Arg::with_name("null_model").long("null").help("Flags the usage of the null model"))

              .get_matches();

        let model = if matches.is_present("HoC") {
            let model_params = values_t!(matches.values_of("HoC"), f64).unwrap();
            FitnessModel::<S>::new_hoc(model_params)
        } else if matches.is_present("additive") {
            let model_params = values_t!(matches.values_of("additive"), f64).unwrap();
            FitnessModel::<S>::new_additive(model_params)
        } else if matches.is_present("RMF") {
            let model_params = values_t!(matches.values_of("RMF"), f64).unwrap();
            FitnessModel::<S>::new_rmf(model_params)
        } else {
            panic!("No model found!")
        };

        let resources_v = values_t!(matches.values_of("resources"), f64).unwrap();
        let mut resources = Vector::<S>::new();
        for i in 0..S {
            resources[i] = resources_v[i];
        }

        let null_model = matches.is_present("null_model");
        let load_landscape = matches.is_present("load_landscape");

        let landscapes: [usize; 2] = values_t!(matches.values_of("landscapes"), usize).unwrap().try_into().unwrap();

        Self {
            pop_size: values_t!(matches.values_of("population_size"), usize).unwrap(),
            mutation_rate_per_locus: value_t!(matches.value_of("mutation_rate_per_locus"), f64).unwrap(),
            model,
            replicates: value_t!(matches.value_of("replicates"), usize).unwrap(),
            resources,
            landscapes,
            null_model,
            load_landscape,
            folder_name: "".to_string()
        }
    }

    pub fn from_command_line_landscape() -> Self {
        let rn: Vec<String> = (0..S).map(|i| format!("res {}", i+1)).collect();
        let _resource_names: Vec<&str> = rn.iter().map(|s| s.as_str()).collect();

        let matches = App::new("")
              .author("André Amado <andre.amado@pm.me>")
              .setting(AppSettings::AllowNegativeNumbers)
              // General arguments

              .arg(Arg::with_name("landscapes").long("landscapes").short("l").takes_value(true).help("Number of landscapes to analize").required(true))

              // Models
              .arg(Arg::with_name("HoC").long("hoc").help("House of Cards model").takes_value(true).value_names(&["cb_diag", "cb_offdiag"]))
              .arg(Arg::with_name("additive").long("add").help("Additive model").takes_value(true).value_names(&["mu", "ca_diag", "ca_offdiag"]))
              .arg(Arg::with_name("RMF").long("rmf").help("Rough Mount Fuji model").takes_value(true).value_names(&["mu", "ca_diag", "ca_offdiag", "cb_diag", "cb_offdiag"]))
              .group(ArgGroup::with_name("model").args(&["HoC", "additive", "RMF"]).required(true))

              .get_matches();

        let model = if matches.is_present("HoC") {
            let model_params = values_t!(matches.values_of("HoC"), f64).unwrap();
            FitnessModel::<S>::new_hoc(model_params)
        } else if matches.is_present("additive") {
            let model_params = values_t!(matches.values_of("additive"), f64).unwrap();
            FitnessModel::<S>::new_additive(model_params)
        } else if matches.is_present("RMF") {
            let model_params = values_t!(matches.values_of("RMF"), f64).unwrap();
            FitnessModel::<S>::new_rmf(model_params)
        } else {
            panic!("No model found!")
        };

        let mut resources = Vector::<S>::new();
        for i in 0..S {
            resources[i] = 1.;
        }

        Self {
            pop_size: vec![0],
            mutation_rate_per_locus: 0.,
            model,
            replicates: 0,
            resources,
            landscapes: [value_t!(matches.value_of("landscapes"), usize).unwrap(), 0],
            null_model: false,
            load_landscape: false,
            folder_name: "".to_string()
        }
    }

    pub fn from_command_line_convergence() -> Self {
        let rn: Vec<String> = (0..S).map(|i| format!("res {}", i+1)).collect();
        let resource_names: Vec<&str> = rn.iter().map(|s| s.as_str()).collect();

        let matches = App::new("")
              .author("André Amado <andre.amado@pm.me>")
              .setting(AppSettings::AllowNegativeNumbers)
              // General arguments
              .arg(Arg::with_name("population_size").help("List of population sizes").short("s").long("size").takes_value(true).multiple(true).required(true))
              .arg(Arg::with_name("mutation_rate_per_locus").help("Mutation rate per locus per generation").short("m").long("mutation_rate").value_name("rate").takes_value(true).required(true))
              .arg(Arg::with_name("resources").help("Amount of each resource").short("r").long("resources").takes_value(true).value_names(&resource_names[..]).required(true))

              .arg(Arg::with_name("landscape").long("landscape").short("l").takes_value(true).help("Index of the landscape to analize").required(true))
              .arg(Arg::with_name("folder").long("folder").short("f").takes_value(true).help("Name of the folder where to store the results").required(true))

              // Models
              .arg(Arg::with_name("HoC").long("hoc").help("House of Cards model").takes_value(true).value_names(&["cb_diag", "cb_offdiag"]))
              .arg(Arg::with_name("additive").long("add").help("Additive model").takes_value(true).value_names(&["mu", "ca_diag", "ca_offdiag"]))
              .arg(Arg::with_name("RMF").long("rmf").help("Rough Mount Fuji model").takes_value(true).value_names(&["mu", "ca_diag", "ca_offdiag", "cb_diag", "cb_offdiag"]))
              .group(ArgGroup::with_name("model").args(&["HoC", "additive", "RMF"]).required(true))

              .arg(Arg::with_name("null_model").long("null").help("Flags the usage of the null model"))

              .get_matches();

        let model = if matches.is_present("HoC") {
            let model_params = values_t!(matches.values_of("HoC"), f64).unwrap();
            FitnessModel::<S>::new_hoc(model_params)
        } else if matches.is_present("additive") {
            let model_params = values_t!(matches.values_of("additive"), f64).unwrap();
            FitnessModel::<S>::new_additive(model_params)
        } else if matches.is_present("RMF") {
            let model_params = values_t!(matches.values_of("RMF"), f64).unwrap();
            FitnessModel::<S>::new_rmf(model_params)
        } else {
            panic!("No model found!")
        };

        let resources_v = values_t!(matches.values_of("resources"), f64).unwrap();
        let mut resources = Vector::<S>::new();
        for i in 0..S {
            resources[i] = resources_v[i];
        }

        let null_model = matches.is_present("null_model");

        let mut folder_name = value_t!(matches.value_of("folder"), String).unwrap();
        if folder_name.len() > 0 && !folder_name.ends_with("/") {
            folder_name.push('/');
            println!("Warning: / appended to folder name ({})", folder_name);
        }

        Self {
            pop_size: values_t!(matches.values_of("population_size"), usize).unwrap(),
            mutation_rate_per_locus: value_t!(matches.value_of("mutation_rate_per_locus"), f64).unwrap(),
            model,
            replicates: 0,
            resources,
            landscapes: [value_t!(matches.value_of("landscape"), usize).unwrap(), 0],
            null_model,
            load_landscape: true,
            folder_name
        }
    }

}
