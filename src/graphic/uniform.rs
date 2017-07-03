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

use std::fmt::Debug;
use std::ffi::CString;
use gl;
use gl::types::*;
use cg;

use super::program::ProgramHandle;

pub trait UniformFieldType {
    unsafe fn put(&self, ix: GLint);
}

impl UniformFieldType for f32 {
    unsafe fn put(&self, ix: GLint) {
        gl::Uniform1f(ix, *self);
    }
}

impl UniformFieldType for cg::Vector2<f32> {
    unsafe fn put(&self, ix: GLint) {
        gl::Uniform2fv(ix, 1, self as *const Self as *const f32);
    }
}

impl UniformFieldType for cg::Vector3<f32> {
    unsafe fn put(&self, ix: GLint) {
        gl::Uniform3fv(ix, 1, self as *const Self as *const f32);
    }
}

impl UniformFieldType for cg::Vector4<f32> {
    unsafe fn put(&self, ix: GLint) {
        gl::Uniform4fv(ix, 1, self as *const Self as *const f32);
    }
}

impl UniformFieldType for cg::Matrix4<f32> {
    unsafe fn put(&self, ix: GLint) {
        gl::UniformMatrix4fv(ix, 1, 0, self as *const Self as *const f32);
    }
}

impl UniformFieldType for GLuint {
    unsafe fn put(&self, ix: GLint) {
        gl::Uniform1ui(ix, *self);
    }
}

pub trait Uniform {
    type Binding : Copy + Debug;

    fn bind(program: &ProgramHandle) -> Result<Self::Binding,String>;
    fn put(&self, binding: &Self::Binding);
}

pub fn try_bind_uniform_field(program: &ProgramHandle, name: &str)
                              -> Result<GLint,String> {
unsafe {
    let ix = gl::GetUniformLocation(
        program.raw(), CString::new(name).unwrap().as_ptr());
    if -1 == ix {
        Err(format!("Uniform {} could not be bound", name))
    } else {
        Ok(ix)
    }
} }

#[macro_export]
macro_rules! uniform {
    ($name:ident, $binding:ident;
     $($field:ident: $typ:ty,)*) => {
        #[derive(Copy,Clone,Debug)]
        pub struct $name {
            $(pub $field: $typ),*
        }

        #[derive(Copy,Clone,Debug)]
        pub struct $binding {
            $($field: ::gl::types::GLint,)*
        }

        impl $crate::graphic::uniform::Uniform for $name {
            type Binding = $binding;

            fn bind(program: &$crate::graphic::program::ProgramHandle)
                    -> Result<$binding,String> {
                Ok($binding {
                    $($field: try!(
                        $crate::graphic::uniform::try_bind_uniform_field(
                            program, stringify!($field))),)*
                })
            }

            fn put(&self, binding: &$binding) {
            unsafe {
                use $crate::graphic::uniform::UniformFieldType;
                $(self.$field.put(binding.$field);)*
            } }
        }
    }
}

uniform! {
    TestUniform, TestUniformBinding;
    xform: cg::Matrix4<f32>,
    t: f32,
}
