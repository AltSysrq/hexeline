//-
// Copyright (c) 2016, 2017, Jason Lingle
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

#![cfg_attr(not(test), allow(unused_features))]
#![feature(test, platform_intrinsics, core_intrinsics, repr_simd, stdsimd)]
#![allow(dead_code)]

extern crate arc_swap;
extern crate arrayvec;
extern crate bit_set;
#[macro_use] extern crate bitflags;
extern crate cgmath as cg;
extern crate env_logger;
extern crate gl;
#[macro_use] extern crate lazy_static;
#[macro_use] extern crate log;
extern crate odds;
extern crate sdl2;
extern crate packed_simd as simd;
extern crate smallvec;

#[cfg(test)] extern crate fnv;
#[cfg(test)] #[macro_use] extern crate proptest;
#[cfg(test)] extern crate test;

use std::env;
use std::os::raw::c_char;
use std::ffi::CStr;
use std::io::{self, Write};

use gl::types::*;

#[macro_use]
mod boolext;
mod intext;
mod simdext;
#[macro_use]
mod graphic;
mod physics;
mod rpn;

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

    env_logger::init();

    match env::args().nth(1).as_ref().map(|s| &**s) {
        Some("rpn") => {
            rpn::run_rpn();
            return;
        },
        Some("circle") => {
            let mut theta = ::std::num::Wrapping(0);
            loop {
                let expected = (theta.0 as f64 / 32768.0 * 3.14159).cos();
                let actual =
                    physics::xform::Affine2d::rotate(theta).repr().extract(0)
                    as f64 / 1024.0;
                println!("{:5}: {} vs {} (rel {})",
                         theta.0, expected, actual, actual / expected);

                theta += ::std::num::Wrapping(256);
                if 0 == theta.0 { break; }
            }
            return;
        },
        _ => (),
    }

    let sdl_context = sdl2::init().unwrap_or_else(die);
    let _sdl_audio = sdl_context.audio().unwrap_or_else(die);
    let _sdl_event = sdl_context.event().unwrap_or_else(die);
    let _sdl_game_controller = sdl_context.game_controller()
        .unwrap_or_else(die);
    let sdl_video = sdl_context.video().unwrap_or_else(die);
    let mut sdl_event_pump = sdl_context.event_pump().unwrap_or_else(die);

    sdl_video.gl_attr().set_context_profile(
        sdl2::video::GLProfile::GLES);
    sdl_video.gl_attr().set_context_version(2, 0);

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

    fn get_string(val: GLenum) -> &'static CStr {
        unsafe {
            CStr::from_ptr(gl::GetString(val) as *const c_char)
        }
    }

    info!("GL Version: {} by {}",
          get_string(gl::VERSION).to_string_lossy(),
          get_string(gl::VENDOR).to_string_lossy());
    info!("GLSL Version: {}",
          get_string(gl::SHADING_LANGUAGE_VERSION).to_string_lossy());
    info!("GL Renderer: {}",
          get_string(gl::RENDERER).to_string_lossy());

    unsafe {
        gl!(Viewport,
            0, 0, screen.drawable_size().0 as GLsizei,
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

    let (sw, sh) = screen.drawable_size();
    let mut texture = graphic::Texture::new();

    let mut mouse_x = 0;
    let mut mouse_y = 0;
    sdl_context.mouse().show_cursor(false);

    'main_loop: loop {
        test_coords(&mut texture, sw, sh, mouse_x, mouse_y);
        draw(&projection_matrix, &shaders, &texture, sw, sh);
        screen.gl_swap_window();

        for event in sdl_event_pump.poll_iter() {
            use sdl2::event::Event;

            match event {
                Event::Quit { .. } => {
                    break 'main_loop;
                },
                Event::MouseMotion { x, y, .. } => {
                    mouse_x = x;
                    mouse_y = y;
                },
                _ => (),
            }
        }
    }
}

fn draw(matrix: &cg::Matrix4<f32>, shaders: &graphic::Shaders,
        texture: &graphic::Texture, sw: u32, sh: u32) {
    use graphic::*;

    unsafe {
        gl!(Clear, gl::COLOR_BUFFER_BIT);
    }

    let vertices = [
        vert::Pos2Tc2 { v: cg::Vector2 { x: 0.0, y: 0.0 },
                        tc: cg::Vector2 { x: 0.0, y: 0.0 } },
        vert::Pos2Tc2 { v: cg::Vector2 { x: sw as f32, y: 0.0 },
                        tc: cg::Vector2 { x: 1.0, y: 0.0 } },
        vert::Pos2Tc2 { v: cg::Vector2 { x: 0.0, y: sh as f32 },
                        tc: cg::Vector2 { x: 0.0, y: 1.0 } },
        vert::Pos2Tc2 { v: cg::Vector2 { x: sw as f32, y: sh as f32 },
                        tc: cg::Vector2 { x: 1.0, y: 1.0 } },
    ];
    let vbo = Vbo::<vert::Pos2Tc2>::new(gl::ARRAY_BUFFER).unwrap();
    let vbo_active = vbo.activate();
    vbo_active.data(&vertices, gl::STREAM_DRAW);

    let vao = Vao::new(&shaders.texture, &vbo_active).unwrap();
    texture.bind();
    let uniform = uni::MTex {
        matrix: *matrix,
        tex: 0,
    };
    vao.activate(&uniform);
    unsafe {
        gl!(DrawArrays, gl::TRIANGLE_STRIP, 0, 4);
    }
}

fn test_coords(tex: &mut graphic::Texture, w: u32, h: u32,
               mouse_x: i32, mouse_y: i32) {
    use std::cmp::{max, min};
    use std::collections::HashSet;
    use std::num::Wrapping;

    use physics::coords::*;
    use physics::xform::Affine2dH;

    let mouse_pos = Vos(mouse_x * 65536 / w as i32, mouse_y * 65536 / w as i32)
        .to_vhr();
    let (mouse_covering_a, mouse_covering_b) =
        mouse_pos.single().to_grid_overlap();
    let mouse_covering: HashSet<(i32,i32)> = (0..4).into_iter()
        .map(|i| (mouse_covering_a.extract(i), mouse_covering_b.extract(i)))
        .filter(|&(a, _)| -32768 != a)
        .collect();

    let mut data = vec![0u32; w as usize * h as usize];
    for y in 0..h {
        for x in 0..w {
            let cart = Vos((x * 65536 / w) as i32, (y * 65536 / w) as i32);
            let hexa = cart.to_vhr();
            let hex = hexa.single().to_index();

            let mut rg = (hex.1 * 16) as u8;
            let mut b = (hex.0 * 16) as u8;

            if (hexa.a() & CELL_HEX_MASK) <= 32 {
                rg ^= 255;
            }
            if (hexa.b() & CELL_HEX_MASK) <= 32 {
                b ^= 255;
            }

            if mouse_covering.contains(&hex) {
                rg ^= 255;
                b ^= 255;
            }

            let recart = hexa.dual().to_vod();
            if recart.x() * w as i32 / 65536 & 255 < 16 ||
                recart.y() * w as i32 / 65536 & 255 < 16
            {
                rg = 255;
                b = 255;
            }

            let mut r = rg;
            let g = rg;

            let mouse_dist = mouse_pos - hexa;
            let mouse_dist = max(max(mouse_dist.a(), mouse_dist.b()),
                                 mouse_dist.c()) -
                min(min(mouse_dist.a(), mouse_dist.b()),
                    mouse_dist.c());
            if mouse_dist.abs() <= CELL_HEX_SIZE {
                r = (Wrapping(r) + Wrapping(128)).0;
            }

            data[(y * w + x) as usize] = graphic::texture::rgba(r, g, b, 255);
        }
    }

    let radius = Vhd(100, 0);
    for theta in -32768i32..32768i32 {
        let affine = Affine2dH::rotate_hex(Wrapping(theta as i16));
        let xformed = (affine * radius).to_vod();
        let px = (xformed + Vod(256, 256)) * Vod(2,2);
        data[(px.y()*(w as i32) + px.x()) as usize] |=
            graphic::texture::rgba(255, 0, 0, 0);
        data[(px.y()*(w as i32) + px.x() + 1) as usize] |=
            graphic::texture::rgba(255, 0, 0, 0);
        data[(px.y()*(w as i32) + px.x() + w as i32) as usize] |=
            graphic::texture::rgba(255, 0, 0, 0);
        data[(px.y()*(w as i32) + px.x() + w as i32 + 1) as usize] |=
            graphic::texture::rgba(255, 0, 0, 0);
    }

    tex.blit_rgba(&data, w, w, h, false);
}
