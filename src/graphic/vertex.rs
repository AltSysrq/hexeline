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

use std::ffi::CString;
use gl;
use gl::types::*;

use super::program::ProgramHandle;

pub unsafe trait VertexFieldType {
    fn element_count() -> GLint;
    fn element_type() -> GLenum;
    fn normalized() -> GLboolean { gl::FALSE }
}

unsafe impl VertexFieldType for f32 {
    fn element_count() -> GLint { 1 }
    fn element_type() -> GLenum { gl::FLOAT }
}

unsafe impl VertexFieldType for [f32;2] {
    fn element_count() -> GLint { 2 }
    fn element_type() -> GLenum { gl::FLOAT }
}

unsafe impl VertexFieldType for [f32;3] {
    fn element_count() -> GLint { 3 }
    fn element_type() -> GLenum { gl::FLOAT }
}

unsafe impl VertexFieldType for [f32;4] {
    fn element_count() -> GLint { 4 }
    fn element_type() -> GLenum { gl::FLOAT }
}

pub trait Vertex {
    type Binding;

    fn bind(program: &ProgramHandle) -> Result<Self::Binding,String>;
    fn install(binding: &Self::Binding);
}

pub fn try_bind_vertex_attrib(program: &ProgramHandle, name: &str)
                              -> Result<GLuint,String> {
unsafe {
    let ix = gl::GetAttribLocation(
        program.raw(),  CString::new(name).unwrap().as_ptr());
    if ix < 0 {
        Err(format!("Vertex attribute {} could not be located", name))
    } else {
        Ok(ix as GLuint)
    }
} }

pub unsafe fn install_vertex_attrib<T: VertexFieldType>(
    ix: GLuint, offset: usize, stride: usize)
{
    gl::VertexAttribPointer(ix, T::element_count(), T::element_type(),
                            T::normalized(), stride as GLsizei,
                            offset as *const GLvoid);
    gl::EnableVertexAttribArray(ix);
}

#[macro_export]
macro_rules! vertex {
    ($name:ident, $binding:ident;
     $($field:ident: $typ:ty,)+) => {
        #[derive(Copy,Clone,Debug)]
        #[repr(C)]
        pub struct $name {
            $(pub $field: $typ,)+
        }

        #[derive(Copy,Clone,Debug)]
        pub struct $binding {
            $($field: ::gl::types::GLuint,)*
        }

        impl $crate::graphic::vertex::Vertex for $name {
            type Binding = $binding;

            fn bind(program: &$crate::graphic::program::ProgramHandle)
                    -> Result<$binding,String> {
                Ok($binding {
                    $($field: try!(
                        $crate::graphic::vertex::try_bind_vertex_attrib(
                            program, stringify!($field))),)*
                })
            }

            fn install(binding: &$binding) {
            unsafe {
                fn addr<T>(t: &T) -> usize { t as *const T as usize }
                let base : $name = ::std::mem::uninitialized();

                $($crate::graphic::vertex::install_vertex_attrib::<$typ>(
                    binding.$field, addr(&base.$field) - addr(&base),
                    ::std::mem::size_of::<$name>()));*;

                ::std::mem::forget(base);
            } }
        }
    }
}

vertex! {
    TestVertex, TestVertexBinding;
    v: [f32;2],
    t: f32,
}
