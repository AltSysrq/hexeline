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

use gl::types::*;
use gl;
use std;

#[derive(Debug)]
struct ShaderHandle(GLuint);

impl ShaderHandle {
    pub fn from_source(typ: GLenum, name: &str,
                       source: &str) -> Result<ShaderHandle,String> {
    unsafe {
        let raw = gl::CreateShader(typ);
        if 0 == raw {
            return Err(format!("Unable to allocate shader {}: {}",
                               name, gl::GetError()));
        }

        let wrapped = ShaderHandle(raw);

        let mut composed_source = String::from("#version 130\n");
        composed_source.push_str(source);
        composed_source.push('\x00');
        let srcptr = composed_source.as_ptr() as *const GLchar;
        gl::ShaderSource(wrapped.0, 1, &srcptr, std::ptr::null());
        gl::CompileShader(wrapped.0);

        let mut status = 0 as GLint;
        gl::GetShaderiv(wrapped.0, gl::COMPILE_STATUS, &mut status);
        if 0 == status {
            let mut log_data = [0 as u8;1024];
            let mut log_size = 0 as GLsizei;
            gl::GetShaderInfoLog(
                wrapped.0, std::mem::size_of_val(&log_data) as GLsizei,
                &mut log_size, &mut log_data[0] as *mut u8 as *mut GLchar);

            let log = std::str::from_utf8(&log_data[0..(log_size as usize)])
                .unwrap_or("(compile log is garbage)");
            return Err(format!("Unable to compile shader {}: {}\n{}",
                               name, gl::GetError(), log));
        }

        Ok(wrapped)
    } }
}

impl Drop for ShaderHandle {
    fn drop(&mut self) {
    unsafe {
        gl::DeleteShader(self.0);
    } }
}

#[derive(Debug)]
pub struct ProgramHandle {
    components: Vec<ShaderHandle>,
    handle: GLuint,
}

impl ProgramHandle {
    fn link(name: &str, components: Vec<ShaderHandle>)
            -> Result<ProgramHandle,String> {
    unsafe {
        let raw = gl::CreateProgram();
        if 0 == raw {
            return Err(format!("Unable to allocate shader program: {}",
                               gl::GetError()));
        }

        let wrapped = ProgramHandle {
            components: components,
            handle: raw,
        };

        for component in wrapped.components.iter() {
            gl::AttachShader(wrapped.handle, component.0);
        }
        gl::LinkProgram(wrapped.handle);

        let mut status = 0 as GLint;
        gl::GetProgramiv(wrapped.handle, gl::LINK_STATUS, &mut status);
        if 0 == status {
            let mut log_data = [0 as u8; 1024];
            let mut log_size = 0 as GLsizei;
            gl::GetProgramInfoLog(
                wrapped.handle, std::mem::size_of_val(&log_data) as GLsizei,
                &mut log_size, &mut log_data[0] as *mut u8 as *mut GLchar);

            let log = std::str::from_utf8(&log_data[0..(log_size as usize)])
                .unwrap_or("(compile log is garbage)");
            return Err(format!("Unable to link shader {}: {}", name, log));
        };

        Ok(wrapped)
    } }

    pub fn of_vf(name: &str, vname: &str, vs: &str,
                 fname: &str, fs: &str)
                 -> Result<ProgramHandle,String> {
        let v = try!(ShaderHandle::from_source(
            gl::VERTEX_SHADER, vname, vs));
        let f = try!(ShaderHandle::from_source(
            gl::FRAGMENT_SHADER, fname, fs));
        ProgramHandle::link(name, vec![v, f])
    }

    pub fn raw(&self) -> GLuint { self.handle }
    pub fn activate(&self) {
    unsafe {
        gl::UseProgram(self.handle);
    } }
}

impl Drop for ProgramHandle {
    fn drop(&mut self) {
    unsafe {
        gl::DeleteProgram(self.handle);
    } }
}
