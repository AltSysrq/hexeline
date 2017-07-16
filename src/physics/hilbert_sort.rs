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

use std::i32;
use std::mem;

use simd::*;

use physics::common_object::CommonObject;

const STRIDE_BITS_BITS: u32 = 2;
const STRIDE_BITS: u32 = 1 << STRIDE_BITS_BITS;

/// Precalculated information for a single stride for Hilbert coordinate
/// conversion.
struct HilbertEntry {
    /// Amount to add to the accumulated coordinate. `2 * STRIDE_BITS` bits are
    /// set.
    addend: u8,
    /// Control information. The low `STRIDE_BITS` bits are to be shifted into
    /// the current stride and XORed with the coordinates. All bits below the
    /// current stride are to be XORed with bit 7. If bit 6 is set, the X and Y
    /// coordinates are to be swapped.
    control: u8,
}

/// Output from `xy_to_hilbert_one_at_a_time`.
#[derive(Clone, Copy, Debug, Default)]
struct HilbertOaat {
    /// The actual Hilbert coordinate.
    coord: u32,
    /// Whether the bits below the ones in question in the input coordinates
    /// were inverted.
    not: u32,
    /// Whether the input coordinates were net swapped at the end of the
    /// function.
    swap: bool,
}

/// Convert (X,Y) coordinates to Hilbert coordinates using the simple
/// one-bit-at-a-time algorithm.
///
/// Besides the Hilbert coordinate, it also returns accumulated internal state
/// used to build the nybble-at-a-time table.
///
/// This is not intended to be particularly fast since it is only used for code
/// generation and testing, so more focus is made on being obviously correct.
///
/// `s` specifies the number of bits to evaluate.
#[allow(dead_code)]
fn xy_to_hilbert_one_at_a_time(
    mut x: u32,
    mut y: u32,
    mut s: u32
) -> HilbertOaat {
    let mut out = HilbertOaat::default();
    // Adapted from
    // https://en.wikipedia.org/wiki/Hilbert_curve#Applications_and_mapping_algorithms
    while s > 0 {
        s -= 1;

        let rx = 0 != x & (1 << s);
        let ry = 0 != y & (1 << s);
        out.coord += ((3 * rx as u32) ^ ry as u32) << s << s;

        if !ry {
            if rx {
                // The Wikipedia algorithm would have us do
                //
                // x = (1 << s) - 1 - x;
                // y = (1 << s) - 1 - y;
                //
                // This makes it less than obvious what's going on, and is
                // prone to overflow.
                //
                // First, note that we do not care about the bits at or above
                // `s` anywhere after this. Below that point, all bits on the
                // LHS of the subtraction are 1. This means that we're really
                // just XORing the coordinates with `(1 << s) - 1`.
                x ^= (1 << s) - 1;
                y ^= (1 << s) - 1;
                out.not = !out.not;
            }

            mem::swap(&mut x, &mut y);
            out.swap = !out.swap;
        }
    }

    out
}

#[allow(dead_code)]
fn gen_hilbert_table() {
    println!("static HILBERT_TABLE: [u16;{}] = [",
             1 << STRIDE_BITS << STRIDE_BITS);
    for y in 0..(1 << STRIDE_BITS) {
        print!("   ");
        for x in 0..(1 << STRIDE_BITS) {
            let out = xy_to_hilbert_one_at_a_time(x, y, STRIDE_BITS);
            let val = ((out.not as u32) << 15) | ((out.swap as u32) << 14) |
                out.coord;
            print!(" 0x{:04X},", val);
            if 7 == x % 8 && x + 1 < (1 << STRIDE_BITS) {
                print!("\n   ");
            }
        }
        println!();
    }
    println!("];");
}

/// Precomputed information for quickly converting to Hilbert coordinates.
///
/// Bits 0 through `2 * STRIDE_BITS` are the value to be shifted and ORed into
/// the accumulator.
///
/// If bit 15 is set, the X and Y coordinates should be binary inverted before
/// the next iteration.
///
/// If bit 14 is set, the X and Y coordinates should be swapped before the next
/// iteration.
static HILBERT_TABLE: [u16;256] = [
    0x0000, 0x4001, 0xC00E, 0x000F, 0x4010, 0xC013, 0x0014, 0x4015,
    0xC0EA, 0x00EB, 0x40EC, 0xC0EF, 0x00F0, 0x40F1, 0xC0FE, 0x00FF,
    0x8003, 0x4002, 0xC00D, 0x800C, 0x0011, 0x0012, 0x8017, 0x4016,
    0xC0E9, 0x80E8, 0x00ED, 0x00EE, 0x80F3, 0x40F2, 0xC0FD, 0x80FC,
    0x4004, 0xC007, 0x4008, 0xC00B, 0x801E, 0x801D, 0x0018, 0x4019,
    0xC0E6, 0x00E7, 0x80E2, 0x80E1, 0x40F4, 0xC0F7, 0x40F8, 0xC0FB,
    0x0005, 0x0006, 0x0009, 0x000A, 0x401F, 0xC01C, 0x801B, 0x401A,
    0xC0E5, 0x80E4, 0x40E3, 0xC0E0, 0x00F5, 0x00F6, 0x00F9, 0x00FA,
    0x803A, 0x8039, 0x8036, 0x8035, 0x4020, 0xC023, 0x0024, 0x4025,
    0xC0DA, 0x00DB, 0x40DC, 0xC0DF, 0x80CA, 0x80C9, 0x80C6, 0x80C5,
    0x403B, 0xC038, 0x4037, 0xC034, 0x0021, 0x0022, 0x8027, 0x4026,
    0xC0D9, 0x80D8, 0x00DD, 0x00DE, 0x40CB, 0xC0C8, 0x40C7, 0xC0C4,
    0x003C, 0x403D, 0xC032, 0x0033, 0x802E, 0x802D, 0x0028, 0x4029,
    0xC0D6, 0x00D7, 0x80D2, 0x80D1, 0x00CC, 0x40CD, 0xC0C2, 0x00C3,
    0x803F, 0x403E, 0xC031, 0x8030, 0x402F, 0xC02C, 0x802B, 0x402A,
    0xC0D5, 0x80D4, 0x40D3, 0xC0D0, 0x80CF, 0x40CE, 0xC0C1, 0x80C0,
    0x4040, 0xC043, 0x0044, 0x4045, 0xC07A, 0x007B, 0x407C, 0xC07F,
    0x4080, 0xC083, 0x0084, 0x4085, 0xC0BA, 0x00BB, 0x40BC, 0xC0BF,
    0x0041, 0x0042, 0x8047, 0x4046, 0xC079, 0x8078, 0x007D, 0x007E,
    0x0081, 0x0082, 0x8087, 0x4086, 0xC0B9, 0x80B8, 0x00BD, 0x00BE,
    0x804E, 0x804D, 0x0048, 0x4049, 0xC076, 0x0077, 0x8072, 0x8071,
    0x808E, 0x808D, 0x0088, 0x4089, 0xC0B6, 0x00B7, 0x80B2, 0x80B1,
    0x404F, 0xC04C, 0x804B, 0x404A, 0xC075, 0x8074, 0x4073, 0xC070,
    0x408F, 0xC08C, 0x808B, 0x408A, 0xC0B5, 0x80B4, 0x40B3, 0xC0B0,
    0x0050, 0x4051, 0xC05E, 0x005F, 0x0060, 0x4061, 0xC06E, 0x006F,
    0x0090, 0x4091, 0xC09E, 0x009F, 0x00A0, 0x40A1, 0xC0AE, 0x00AF,
    0x8053, 0x4052, 0xC05D, 0x805C, 0x8063, 0x4062, 0xC06D, 0x806C,
    0x8093, 0x4092, 0xC09D, 0x809C, 0x80A3, 0x40A2, 0xC0AD, 0x80AC,
    0x4054, 0xC057, 0x4058, 0xC05B, 0x4064, 0xC067, 0x4068, 0xC06B,
    0x4094, 0xC097, 0x4098, 0xC09B, 0x40A4, 0xC0A7, 0x40A8, 0xC0AB,
    0x0055, 0x0056, 0x0059, 0x005A, 0x0065, 0x0066, 0x0069, 0x006A,
    0x0095, 0x0096, 0x0099, 0x009A, 0x00A5, 0x00A6, 0x00A9, 0x00AA,
];

#[allow(unused_assignments)]
fn xy_to_hilbert(pos: i32x4) -> u64 {
    // Deal with negative values by just offsetting everything by (1<<31) and
    // treating the vector as unsigned.
    let pos = pos - i32x4::splat(i32::MIN);
    // Drop some precision to shave off a couple iterations
    let pos = pos >> 2*STRIDE_BITS;

    let mut d = 0u64;
    // The rest of the algorithm is not really amenable to SSE, so move to
    // scalar registers.
    let x = pos.extract(0) as u32;
    let y = pos.extract(1) as u32;
    let mut invert_mask = 0u32;
    let mut xs = 0;
    let mut ys = STRIDE_BITS;

    // This is a branch-free nybble-at-a-time version of the earlier
    // one-at-a-time algorithm. This works on the observation that we work
    // strictly left to right over the coordinates, and a prior operation can
    // have a combination of two net effects:
    //
    // - Swap the coordinates.
    // - Invert all the bits we care about.
    //
    // Instead of swapping the coordinates themselves, we swap the shifts used
    // to combine X and Y into a table index. Instead of inverting the bits in
    // the coordinates, we simply track a mask to XOR into the table index.
    macro_rules! stride {
        ($s:expr) => { {
            let suby = (y >> $s) & ((1 << STRIDE_BITS) - 1);
            let subx = (x >> $s) & ((1 << STRIDE_BITS) - 1);
            let ix = ((suby << ys) | (subx << xs)) ^ invert_mask;
            // The AND here is necessary since invert_mask is either 0 or
            // u32::MAX. We thus need the AND *somewhere*, and putting it here
            // ensures the optimiser will be able to remove the array bounds
            // check.
            let val = HILBERT_TABLE[ix as usize & ((1 << 2*STRIDE_BITS) - 1)];

            d += ((val  & ((1 << 2*STRIDE_BITS) - 1)) as u64) << 2*$s;
            invert_mask ^= ((val as i16 as i32) >> 31) as u32;
            let swap_mask = (val as u32) >> 14 << STRIDE_BITS_BITS;
            xs ^= swap_mask;
            ys ^= swap_mask;
        } }
    }

    stride!(24);
    stride!(20);
    stride!(16);
    stride!(12);
    stride!( 8);
    stride!( 4);
    stride!( 0);
    d
}

pub fn hilbert_sort(array: &mut [CommonObject]) {
    array.sort_by_key(|o| xy_to_hilbert(o.p));
}

#[cfg(test)]
mod test {
    use test::{Bencher, black_box};

    use physics::common_object::UnpackedCommonObject;
    use super::*;

    #[bench]
    fn bench_hilbert_sort(b: &mut Bencher) {
        let mut data: Vec<CommonObject> =
            (0..10000).into_iter().map(|i| UnpackedCommonObject {
                a: i * 256, .. UnpackedCommonObject::default()
            }.pack()).collect();
        b.iter(|| hilbert_sort(&mut data));
    }

    #[bench]
    fn bench_xy_to_hilbert(b: &mut Bencher) {
        b.iter(|| xy_to_hilbert(black_box(i32x4::splat(65536))));
    }
}
