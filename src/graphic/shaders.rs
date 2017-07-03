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

use super::program::ProgramHandle;
use super::shader::Shader;

pub mod vert {
    use cg;

    vertex! {
        Pos2, Pos2Binding;
        v: cg::Vector2<f32>,
    }

    vertex! {
        Pos2Tc2, Pos2Tc2Binding;
        v: cg::Vector2<f32>,
        tc: cg::Vector2<f32>,
    }
}

pub mod uni {
    use gl::types::*;

    use cg;

    uniform! {
        MColour, MColourBinding;
        matrix: cg::Matrix4<f32>,
        colour: cg::Vector4<f32>,
    }

    uniform! {
        MTex, MTexBinding;
        matrix: cg::Matrix4<f32>,
        tex: GLuint,
    }
}

pub struct ShaderPrograms {
    flat: ProgramHandle,
    texture: ProgramHandle,
}

static F_FLAT_SRC : &'static str = include_str!("f_flat.glsl");
static F_TEX_SRC  : &'static str = include_str!("f_tex.glsl");
static V_POS2_SRC : &'static str = include_str!("v_pos2.glsl");
static V_POS2TC2_SRC : &'static str = include_str!("v_pos2tc2.glsl");

impl ShaderPrograms {
    pub fn new() -> Result<ShaderPrograms,String> {
        Ok(ShaderPrograms {
            flat: ProgramHandle::of_vf(
                "flat", "v_pos2", V_POS2_SRC, "f_flat", F_FLAT_SRC)?,
            texture: ProgramHandle::of_vf(
                "texture", "v_pos2tc2", V_POS2TC2_SRC, "f_tex", F_TEX_SRC)?,
        })
    }
}

pub struct Shaders<'a> {
    pub flat: Shader<'a, uni::MColour, vert::Pos2>,
    pub texture: Shader<'a, uni::MTex, vert::Pos2Tc2>,
}

impl<'a> Shaders<'a> {
    pub fn new(progs: &'a ShaderPrograms) -> Result<Shaders<'a>,String> {
        Ok(Shaders {
            flat: Shader::of(&progs.flat)?,
            texture: Shader::of(&progs.texture)?,
        })
    }
}
