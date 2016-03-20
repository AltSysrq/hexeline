//-
// Copyright (c) 2016, Jason Lingle
//
// Permission to  use, copy,  modify, and/or distribute  this software  for any
// purpose  with or  without fee  is hereby  granted, provided  that the  above
// copyright notice and this permission notice appear in all copies.
//
// THE SOFTWARE  IS PROVIDED "AS  IS" AND  THE AUTHOR DISCLAIMS  ALL WARRANTIES
// WITH  REGARD   TO  THIS  SOFTWARE   INCLUDING  ALL  IMPLIED   WARRANTIES  OF
// MERCHANTABILITY AND FITNESS. IN NO EVENT  SHALL THE AUTHOR BE LIABLE FOR ANY
// SPECIAL,  DIRECT,   INDIRECT,  OR  CONSEQUENTIAL  DAMAGES   OR  ANY  DAMAGES
// WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
// OF  CONTRACT, NEGLIGENCE  OR OTHER  TORTIOUS ACTION,  ARISING OUT  OF OR  IN
// CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

extern crate libc;
extern crate glium;
extern crate sdl2;

use std::io;
use std::io::Write;

mod glium_sdl2;

mod praef;
mod physics;

use glium_sdl2::DisplayBuild;

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
    let mut sdl_event_pump = sdl_context.event_pump().unwrap_or_else(die);

    sdl_video.gl_attr().set_context_profile(
        sdl2::video::GLProfile::Core);

    let current_mode = sdl_video.current_display_mode(0)
        .unwrap_or_else(die);
    let display = DisplayBuild::build_glium(
        &mut sdl_video.window(
            "Hexeline",
            to_window_size(current_mode.w, 640),
            to_window_size(current_mode.h, 480)))
        .unwrap_or_else(|error| {
            match error {
                glium::GliumCreationError::BackendCreationError(
                    sdl2::video::WindowBuildError::SdlError(message)) =>
                    die(message),
                glium::GliumCreationError::IncompatibleOpenGl(message) =>
                    die(format!("Incompatible OpenGL: {}", message)),
                any => die(format!("Unexpected error: {:?}", any)),
            }
        });

    'main_loop: loop {
        let mut target = display.draw();
        target.finish().unwrap();

        for event in sdl_event_pump.poll_iter() {
            use sdl2::event::Event;

            match event {
                Event::Quit { .. } => {
                    break 'main_loop;
                },
                _ => (),
            }
        }
    }
}
