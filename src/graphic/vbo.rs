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

use std::marker::PhantomData;
use std::mem;
use gl;
use gl::types::*;

#[derive(Debug)]
pub struct Vbo<V> {
    target: GLenum,
    handle: GLuint,
    _vertices: PhantomData<V>,
}

pub struct ActiveVbo<'a,V:'a>(pub &'a Vbo<V>);

impl<V> Vbo<V> {
    pub fn new(target: GLenum) -> Result<Vbo<V>,String> {
    unsafe {
        let mut raw = 0 as GLuint;
        gl::GenBuffers(1, &mut raw);

        if 0 == raw {
            Err("Failed to allocate VBO".to_string())
        } else {
            Ok(Vbo { target: target, handle: raw, _vertices: PhantomData })
        }
    } }

    pub fn activate<'a>(&'a self) -> ActiveVbo<'a,V> {
        unsafe {
            gl::BindBuffer(self.target, self.handle);
        }
        ActiveVbo(self)
    }
}

impl<V> Drop for Vbo<V> {
    fn drop(&mut self) {
    unsafe {
        gl::DeleteBuffers(1, &self.handle);
    } }
}

impl<'a,V:'a> ActiveVbo<'a,V> {
    pub fn data(&self, data: &[V], usage: GLenum) {
    unsafe {
        gl::BufferData(self.0.target,
                       (data.len() * mem::size_of::<V>()) as GLsizeiptr,
                       &data[0] as *const V as *const GLvoid, usage);
    } }
}
