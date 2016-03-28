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

use gl;
use gl::types::*;

use super::vbo::{Vbo,ActiveVbo};
use super::uniform::Uniform;
use super::vertex::Vertex;
use super::shader::{Shader,ActiveShader};

#[derive(Debug)]
pub struct Vao<'a, 'b: 'a, U: Uniform + 'b, V: Vertex + 'b> {
    handle: GLuint,
    shader: &'a Shader<'b,U,V>,
    vbo: &'a Vbo<V>,
}

impl<'a, 'b: 'a, U: Uniform + 'b, V: Vertex + 'b> Drop for Vao<'a,'b,U,V> {
    fn drop(&mut self) {
    unsafe {
        gl::DeleteVertexArrays(1, &self.handle);
    } }
}

impl<'a, 'b: 'a, U: Uniform + 'b, V: Vertex + 'b> Vao<'a,'b,U,V> {
    pub fn new(shader: &'a Shader<'b,U,V>,
               vbo: &ActiveVbo<'a,V>) -> Result<Self,String> {
        let this = try!(Vao::alloc(shader, vbo.0));
        this.make_current();
        this.set_up();
        Ok(this)
    }

    fn alloc(shader: &'a Shader<'b,U,V>,
             vbo: &'a Vbo<V>) -> Result<Self,String> {
    unsafe {
        let mut raw = 0 as GLuint;
        gl::GenVertexArrays(1, &mut raw);
        if 0 == raw {
            Err("Failed to allocate VAO".to_string())
        } else {
            Ok(Vao { handle: raw, vbo: vbo, shader: shader })
        }
    } }

    fn set_up(&self) {
        V::install(self.shader.vertex_binding());
    }

    pub fn make_current(&self) {
    unsafe {
        gl::BindVertexArray(self.handle);
    } }

    pub fn activate<'c>(&'c self, uniform: &U)
                        -> ActiveShader<'c,'b,U,V> {
        let ret = self.shader.activate(uniform);
        self.make_current();
        ret
    }
}
