use configparser::ini::Ini;
use crate::my_structs::UserOptions;
use std::io::{Error, ErrorKind};

pub fn read_input(file_name: &str) -> Result<UserOptions, std::io::Error> {
    println!("Reading input from {file_name}");

    let mut config = Ini::new();
    let map = config.load(file_name);
    //println!("{:?}", map);
    //let conf = ini::Ini::load_from_file(file_name).unwrap();
    let val = config.get("Database", "Host").unwrap();
    //println!("Host: {}", val);

    //let host = map["Database"]["Host"];

}