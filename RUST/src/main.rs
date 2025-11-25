mod my_structs;
mod input_util;

use std::process;
fn main() {
    let args: Vec<String> = std::env::args().collect();
    dbg!(&args);
    println!("len = {}", args.len());
    // let query = &args[1];
    // let file_path = &args[2];
    // println!("Searching for {query} in {file_path}");
    if args.len() != 3 {
        eprintln!("🛑 Invalid number of arguments.");
        eprintln!("Usage: {} --c config.ini", args[0]);
        process::exit(1);
    }
    else{
        input_util::read_input(&args[2]);
    }
    println!("Hello, world!");
    let user_opt : my_structs::UserOptions = input_util::read_input(&args[2]).unwrap();
}
