extern crate sdl2;

use std::io;
use std::io::Write;

fn main() {
    let sdl_context = match sdl2::init() {
        Ok(v) => v,
        Err(str) => {
            writeln!(&mut io::stderr(), "Failed to initialise SDL: {}", str)
                .unwrap();
            std::process::exit(70);
        }
    };

    println!("hello sdl");
}
