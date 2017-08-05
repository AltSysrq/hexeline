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

//! A primitive interactive RPN calculator with specific knowledge of certain
//! Hexeline concepts.
//!
//! Each whole line of input is interpreted as either an integer or a built-in
//! function.

use std::cmp::max;
use std::fmt;
use std::io::{self, BufRead};
use std::num::Wrapping;

use simdext::*;

use physics::coords::*;
use physics::xform::*;

#[derive(Clone, Copy, Debug)]
enum Value {
    Scalar(i32),
    Vo(Vos),
    Vh(Vhr),
    Affine(Affine2dO),
}
use self::Value::*;

impl fmt::Display for Value {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Scalar(i) => write!(f, "{}", i),
            Vo(v) => write!(f, "{:?}", v),
            Vh(v) => write!(f, "{:?}", v),
            Affine(matrix) =>
                write!(f, "[{:5} {:5}\n\
                     \x20   {:5} {:5} ] / {}",
                         matrix.repr().extract(0),
                         matrix.repr().extract(1),
                         matrix.repr().extract(2),
                         matrix.repr().extract(3),
                       1 << AFFINE_POINT),
        }
    }
}

pub fn run_rpn() {
    let mut stack = Vec::new();

    for line in io::BufReader::new(io::stdin()).lines() {
        let line = line.unwrap();

        macro_rules! pop1 {
            () => {{
                let f = match stack.pop() {
                    Some(f) => f,
                    None => {
                        println!("Stack underflow");
                        continue;
                    }
                };
                f
            }}
        }
        macro_rules! pop2 {
            () => {{
                if stack.len() < 2 {
                    println!("Stack underflow");
                    continue;
                }
                let a = stack.pop().unwrap();
                let b = stack.pop().unwrap();
                (a, b)
            }}
        }
        macro_rules! domain_error {
            () => {{
                println!("Domain error");
                continue;
            }}
        }

        match &*line {
            "drop" => { pop1!(); },
            "vo2" => match pop2!() {
                (Scalar(y), Scalar(x)) => {
                    stack.push(Vo(Vos(x, y)));
                },
                _ => domain_error!(),
            },
            "vh2" => match pop2!() {
                (Scalar(b), Scalar(a)) => {
                    stack.push(Vh(Vhs(a, b).redundant()));
                },
                _ => domain_error!(),
            },
            "rot" => match pop1!() {
                Scalar(theta) => {
                    stack.push(Affine(
                        Affine2d::rotate(Wrapping(theta as i16))));
                },
                _ => domain_error!(),
            },
            "+" => match pop2!() {
                (Scalar(rhs), Scalar(lhs)) => stack.push(Scalar(lhs + rhs)),
                (Vo(rhs), Vo(lhs)) => stack.push(Vo(lhs + rhs)),
                (Vh(rhs), Vh(lhs)) => stack.push(Vh(lhs + rhs)),
                _ => domain_error!(),
            },
            "-" => match pop2!() {
                (Scalar(rhs), Scalar(lhs)) => stack.push(Scalar(lhs - rhs)),
                (Vo(rhs), Vo(lhs)) => stack.push(Vo(lhs - rhs)),
                (Vh(rhs), Vh(lhs)) => stack.push(Vh(lhs - rhs)),
                _ => domain_error!(),
            },
            "*" => match pop2!() {
                (Scalar(rhs), Scalar(lhs)) => stack.push(Scalar(lhs * rhs)),
                (Vo(rhs), Vo(lhs)) => stack.push(Vo(lhs * rhs)),
                (Vh(rhs), Vh(lhs)) => stack.push(Vh(lhs * rhs)),
                (Vo(rhs), Scalar(lhs)) => stack.push(Vo(Vos(lhs, lhs) * rhs)),
                (Scalar(rhs), Vo(lhs)) => stack.push(Vo(lhs * Vos(rhs, rhs))),
                (Vh(rhs), Scalar(lhs)) => stack.push(Vh(Vhr(lhs, lhs, lhs) * rhs)),
                (Vo(rhs), Affine(lhs)) =>
                    stack.push(Vo((lhs * rhs.dual()).single())),
                (Vh(rhs), Affine(lhs)) =>
                    stack.push(Vh((lhs.to_hexagonal() * rhs.dual()).redundant())),
                (Affine(rhs), Affine(lhs)) => stack.push(Affine(lhs * rhs)),
                _ => domain_error!(),
            },
            "/" => match pop2!() {
                (Scalar(rhs), Scalar(lhs)) => stack.push(Scalar(lhs / rhs)),
                _ => domain_error!(),
            },
            "<<" => match pop2!() {
                (Scalar(rhs), Scalar(lhs)) => stack.push(Scalar(lhs << rhs as u32)),
                (Scalar(rhs), Vo(lhs)) => stack.push(Vo(lhs << rhs as u32)),
                (Scalar(rhs), Vh(lhs)) => stack.push(Vh(lhs << rhs as u32)),
                _ => domain_error!(),
            },
            ">>" => match pop2!() {
                (Scalar(rhs), Scalar(lhs)) => stack.push(Scalar(lhs >> rhs as u32)),
                (Scalar(rhs), Vo(lhs)) => stack.push(Vo(lhs >> rhs as u32)),
                (Scalar(rhs), Vh(lhs)) => stack.push(Vh(lhs >> rhs as u32)),
                _ => domain_error!(),
            },
            "ortho" => match pop1!() {
                Vh(v) => stack.push(Vo(v.dual().to_vod().single())),
                _ => domain_error!(),
            },
            "hex" => match pop1!() {
                Vo(v) => stack.push(Vh(v.to_vhr())),
                _ => domain_error!(),
            },
            "->cell" => match pop1!() {
                Scalar(v) => stack.push(Scalar(v >> CELL_HEX_SHIFT)),
                Vh(v) => stack.push(Vh(v >> CELL_HEX_SHIFT)),
                _ => domain_error!(),
            },
            "->cont" => match pop1!() {
                Scalar(v) => stack.push(Scalar(v << CELL_HEX_SHIFT)),
                Vh(v) => stack.push(Vh(v << CELL_HEX_SHIFT)),
                _ => domain_error!(),
            },
            "x" => {
                let (a, b) = pop2!();
                stack.push(a);
                stack.push(b);
            },
            "?" => {
                for elem in &stack {
                    println!(".. {}", elem);
                }
            },
            "dup" => {
                let a = pop1!();
                stack.push(a);
                stack.push(a);
            },
            "CELL_HEX_SHIFT" => stack.push(Scalar(CELL_HEX_SHIFT as i32)),
            "CELL_HEX_MASK" => stack.push(Scalar(CELL_HEX_MASK)),
            "CELL_HEX_SIZE" => stack.push(Scalar(CELL_HEX_SIZE)),
            "CELL_HEX_LINF_EDGE" => stack.push(Scalar(CELL_HEX_LINF_EDGE)),
            "CELL_HEX_LINF_VERTEX" => stack.push(Scalar(CELL_HEX_LINF_VERTEX)),
            "CELL_L2_EDGE" => stack.push(Scalar(CELL_L2_EDGE)),
            "CELL_L2_VERTEX" => stack.push(Scalar(CELL_L2_VERTEX)),
            "l1" => match pop1!() {
                Vo(v) => stack.push(Scalar(v.repr().abs().hsum_2())),
                Vh(v) => stack.push(Scalar(v.repr().abs().hsum_3())),
                _ => domain_error!(),
            },
            "l2" => match pop1!() {
                Vo(v) => {
                    let x = v.x() as i64;
                    let y = v.y() as i64;
                    let dist = ((x*x + y*y) as f64).sqrt().ceil() as i32;
                    stack.push(Scalar(dist));
                },
                Vh(v) => {
                    let a = v.a() as i64;
                    let b = v.b() as i64;
                    let c = v.c() as i64;
                    let dist = ((a*a + b*b + c*c) as f64).sqrt().ceil() as i32;
                    stack.push(Scalar(dist));
                },
                _ => domain_error!(),
            },
            "linf" => match pop1!() {
                Vo(v) => stack.push(Scalar(max(v.x().abs(), v.y().abs()))),
                Vh(v) => stack.push(Scalar(
                    max(v.a().abs(), max(v.b().abs(), v.c().abs())))),
                _ => domain_error!(),
            },
            _ => if let Ok(val) = line.parse() {
                stack.push(Scalar(val))
            } else {
                println!("?");
            },
        }

        if let Some(last) = stack.last() {
            println!("=> {} (+ {})", last, stack.len() - 1);
        } else {
            println!("=> empty stack");
        }
    }
}
