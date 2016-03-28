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

use super::program::*;
use super::uniform::Uniform;
use super::vertex::Vertex;

#[derive(Clone,Copy,Debug)]
pub struct Shader<'a, U: Uniform, V: Vertex> {
    program: &'a ProgramHandle,
    uniform_binding: U::Binding,
    vertex_binding: V::Binding,
}

pub struct ActiveShader<'b,'a: 'b, U: Uniform + 'b, V: Vertex + 'b>(
    pub &'b Shader<'a, U, V>);

impl<'a, U: Uniform, V: Vertex> Shader<'a,U,V> {
    pub fn of(program: &'a ProgramHandle) -> Result<Self,String> {
        Ok(Shader {
            program: program,
            uniform_binding: try!(U::bind(program)),
            vertex_binding: try!(V::bind(program)),
        })
    }

    pub fn vertex_binding(&self) -> &V::Binding {
        &self.vertex_binding
    }

    pub fn activate<'b >(&'b self, uniform: &U)
                         -> ActiveShader<'b,'a,U,V> {
        self.program.activate();
        uniform.put(&self.uniform_binding);
        ActiveShader(self)
    }
}
