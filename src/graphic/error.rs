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

use gl;

#[cfg(debug_assertions)]
#[macro_export]
macro_rules! gl {
    ($name:ident $(, $arg:expr)*) => { {
        let r = ::gl::$name($($arg),*);
        ::graphic::error::check_error(
            concat!(file!(), ":", line!(), " ",
                    stringify!($name)));
        r
    } }
}

#[cfg(not(debug_assertions))]
#[macro_export]
macro_rules! gl {
    ($name:ident $(, $arg:expr)*) => {
        ::gl::$name($($arg),*)
    }
}

#[macro_export]
macro_rules! check_gl_error {
    () => {
        ::graphic::error::check_error(
            concat!(file!(), ":", line!()))
    }
}

pub fn check_error(wo: &str) {
    loop {
        let name = match unsafe { gl::GetError() } {
            gl::NO_ERROR => break,
            gl::INVALID_ENUM => "GL_INVALID_ENUM".to_owned(),
            gl::INVALID_VALUE => "GL_INVALID_VALUE".to_owned(),
            gl::INVALID_OPERATION => "GL_INVALID_OPERATION".to_owned(),
            gl::STACK_OVERFLOW => "GL_STACK_OVERFLOW".to_owned(),
            gl::STACK_UNDERFLOW => "GL_STACK_UNDERFLOW".to_owned(),
            gl::OUT_OF_MEMORY => "GL_OUT_OF_MEMORY".to_owned(),
            code => format!("Unknown, code {}", code),
        };

        error!("{}: {}", wo, name);
    }
}

pub fn ignore_errors() {
    while gl::NO_ERROR != unsafe { gl::GetError() } { }
}
