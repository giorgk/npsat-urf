fn main() {
    let args: Vec<String> = std::env::args().collect();
    dbg!(args.len());
    if (args.len() < 2){
        panic!("Not enough arguments");
    }
    else{
        let query = &args[1];
        println!("Input file {query}");
    }
    dbg!(&args);
    // let query = &args[1];
    // let file_path = &args[2];
    // println!("Searching for {query} in {file_path}");
    println!("Hello, world!");
}
