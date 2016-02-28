extern crate sdl2;

use std::io;
use std::io::Write;

fn main() {
    fn die<T>(message: String) -> T {
        writeln!(&mut io::stderr(), "Failed to initialise SDL: {}", message)
            .unwrap();
        std::process::exit(70)
    }

    fn to_window_size(d: i32, dflt: u32) -> u32 {
        if d > 0 {
            d as u32
        } else {
            dflt
        }
    }

    let sdl_context = sdl2::init().unwrap_or_else(die);
    let _sdl_audio = sdl_context.audio().unwrap_or_else(die);
    let _sdl_event = sdl_context.event().unwrap_or_else(die);
    let _sdl_game_controller = sdl_context.game_controller()
        .unwrap_or_else(die);
    let sdl_video = sdl_context.video().unwrap_or_else(die);
    let _sdl_event_pump = sdl_context.event_pump().unwrap_or_else(die);

    sdl_video.gl_attr().set_context_profile(
        sdl2::video::GLProfile::Core);

    let current_mode = sdl_video.current_display_mode(0)
        .unwrap_or_else(die);
    let window = sdl_video.window(
        "Hexeline",
        to_window_size(current_mode.w, 640),
        to_window_size(current_mode.h, 480))
        .opengl()
        .build().unwrap_or_else(|error| {
            match error {
                sdl2::video::WindowBuildResult::SdlError(message) =>
                    die(message),
                any => die(format!("Unexpected error: {:?}", any)),
            }
        });

    std::thread::sleep(std::time::Duration::from_millis(5000));

    println!("hello sdl");
}
