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

use glium;

#[derive(Copy,Clone,Debug)]
pub struct Vertex {
    pub v: [f32; 2],
}

#[derive(Copy,Clone,Debug)]
pub struct Uniform {
    pub colour: [f32; 4],
}

implement_vertex!(Vertex, v);
implement_uniform_block!(Uniform, colour);

pub fn make(glf: &::Glf) -> glium::program::Program {
    glium::program::Program::from_source(
        glf, include_str!("v_flat.glsl"), include_str!("f_flat.glsl"), None)
        .unwrap()
}
