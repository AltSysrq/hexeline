//-
// Copyright (c) 2017 Jason Lingle
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

use std::borrow::Cow;
use std::cmp::max;
use std::ptr;

use gl;
use gl::types::*;

/// A simple RAII wrapper around a GL texture.
#[derive(Debug)]
pub struct Texture(GLuint, *const ());

impl Texture {
    pub fn new() -> Self {
        let mut unit = 0;
        unsafe {
            gl::GenTextures(1, &mut unit);
        }

        Texture(unit, ptr::null())
    }

    /// Blit `data` (an ARGB texture in machine byte order) onto the texture.
    ///
    /// `w` and `h` specify the dimensions of the texture; `pitch` is the
    /// physical row width.
    ///
    /// If the implementation does not support the data as presented, it is
    /// scaled down until it does. If `mipmap` is true, all mipmap levels will
    /// be generated automatically as well.
    pub fn blit_argb(&mut self, data: &[u32],
                     mut w: u32, mut pitch: u32, mut h: u32,
                     mipmap: bool) {
        assert!(h as usize * pitch as usize <= data.len());

        let mut data = Cow::Borrowed(data);

        // Scale down until the implementation will accept it
        loop {
            unsafe {
                gl::TexImage2D(gl::PROXY_TEXTURE_2D, 0, gl::RGBA as GLint,
                               w as i32, h as i32, 0,
                               gl::BGRA, gl::UNSIGNED_BYTE,
                               data.as_ptr() as *const GLvoid);
                let mut image_is_supported = 0;
                gl::GetTexLevelParameteriv(gl::PROXY_TEXTURE_2D, 0,
                                           gl::TEXTURE_WIDTH,
                                           &mut image_is_supported);
                if 0 != image_is_supported { break; }

                if 1 == w && 1 == h {
                    panic!("Scaled image down to 1x1 and still not accepted \
                            by OpenGL");
                }

                let data2 = data;
                data = Cow::Owned(scale_half(
                    &data2, &mut w, &mut pitch, &mut h));
            }
        }

        unsafe {
            let mut level = 0;
            gl::BindTexture(gl::TEXTURE_2D, self.id());
            loop {
                gl::PixelStorei(gl::UNPACK_ROW_LENGTH, pitch as i32);
                gl::TexImage2D(gl::TEXTURE_2D, level, gl::RGBA as GLint,
                               w as i32, h as i32, 0,
                               // This is presumably sensitive to machine byte
                               // order. TODO Address if we encounter any
                               // big-endian systems we care to support.
                               gl::BGRA, gl::UNSIGNED_BYTE,
                               data.as_ptr() as *const GLvoid);

                if !mipmap || 1 == w && 1 == h { break; }

                let data2 = data;
                data = Cow::Owned(scale_half(
                    &data2, &mut w, &mut pitch, &mut h));
                level += 1;
            }

            gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_MAG_FILTER,
                              gl::NEAREST as i32);
            gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_MIN_FILTER,
                              gl::NEAREST as i32);
            gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_WRAP_S,
                              gl::CLAMP_TO_EDGE as i32);
            gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_WRAP_T,
                              gl::CLAMP_TO_EDGE as i32);

            gl::PixelStorei(gl::UNPACK_ROW_LENGTH, 0);
        }
    }

    pub fn bind(&self) {
        unsafe {
            gl::BindTexture(gl::TEXTURE_2D, self.id());
        }
    }

    /// Return the identifier for this texture.
    pub fn id(&self) -> GLuint {
        self.0
    }
}

impl Drop for Texture {
    fn drop(&mut self) {
        unsafe {
            gl::DeleteTextures(1, &self.0)
        }
    }
}

// We may want something more sophisticated than nearest scaling. For now, use
// the simplest thing that works.
fn scale_half(px: &[u32], w: &mut u32, pitch: &mut u32, h: &mut u32)
              -> Vec<u32> {
    let nw = max(1, *w / 2);
    let nh = max(1, *h / 2);
    let mut dst = Vec::with_capacity(nw as usize * nh as usize);

    for row in 0..nh {
        for col in 0..nw {
            dst.push(px[row as usize * 2 * *pitch as usize + col as usize * 2]);
        }
    }

    *w = nw;
    *pitch = nw;
    *h = nh;
    dst
}
