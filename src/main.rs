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

#![allow(dead_code)]

extern crate sdl2;
extern crate gl;
extern crate cgmath as cg;
extern crate simd;

mod graphic;
mod physics;

use std::io;
use std::io::Write;
use gl::types::*;

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
    sdl_video.gl_attr().set_context_version(3, 0);

    let current_mode = sdl_video.current_display_mode(0)
        .unwrap_or_else(die);
    let screen =
        sdl_video.window(
            "Hexeline",
            to_window_size(current_mode.w, 640),
            to_window_size(current_mode.h, 480))
        .opengl()
        .build().unwrap_or_else(|error| {
            match error {
                sdl2::video::WindowBuildError::SdlError(message) =>
                    die(message),
                any => die(format!("Unexpected error: {:?}", any)),
            }
        });
    let _gl_context = screen.gl_create_context().unwrap_or_else(die);

    // Load all the GL functions, since this doesn't happen on-demand
    gl::load_with(|s| sdl_video.gl_get_proc_address(s) as
                  *const std::os::raw::c_void);

    unsafe {
        gl::Viewport(0, 0, screen.drawable_size().0 as GLsizei,
                     screen.drawable_size().1 as GLsizei);
    }

    let projection_matrix =
        // Use real pixel coordinates. Origin is at the top left, down is
        // positive Y, one unit is one pixel. We'll always be drawing at Z=0,
        // so just position the Z boundaries on either side.
        cg::ortho::<f32>(0.0, screen.drawable_size().0 as f32,
                         screen.drawable_size().1 as f32, 0.0,
                         -1.0, 1.0);

    let shader_programs = graphic::ShaderPrograms::new().unwrap_or_else(die);
    let shaders = graphic::Shaders::new(&shader_programs).unwrap_or_else(die);

    'main_loop: loop {
        draw(&projection_matrix, &shaders);
        screen.gl_swap_window();

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

fn draw(matrix: &cg::Matrix4<f32>, shaders: &graphic::Shaders) {
    use graphic::*;

    unsafe {
        gl::Clear(gl::COLOR_BUFFER_BIT);
    }

    let vertices = [
        vert::Pos2 { v: cg::Vector2 { x: 1280.0 / 2.0, y: 1024.0 / 2.0 } },
        vert::Pos2 { v: cg::Vector2 { x: 1280.0, y: 1024.0 / 2.0 } },
        vert::Pos2 { v: cg::Vector2 { x: 1280.0 / 2.0, y: 0.0 } },
    ];
    let vbo = Vbo::<vert::Pos2>::new(gl::ARRAY_BUFFER).unwrap();
    let vbo_active = vbo.activate();
    vbo_active.data(&vertices, gl::STREAM_DRAW);

    let vao = Vao::new(&shaders.flat, &vbo_active).unwrap();
    let uniform = uni::MColour {
        matrix: *matrix,
        colour: cg::Vector4 { x: 1.0, y: 0.5, z: 0.0, w: 1.0 },
    };
    vao.activate(&uniform);
    unsafe {
        gl::DrawArrays(gl::TRIANGLES, 0, 3);
    }
}
