
use crate::params::*;
/* helper functions from the paper, bits and hints */

const Q_INV: i32 = 62209; // q^(-1) mod 2^16


/* for finite field element a with -2^{31}Q <= a <= Q*2^31,
 compute a' = a*2^{-32} (mod Q) such that -Q < a' < Q.
 https://en.wikipedia.org/wiki/Montgomery_modular_multiplication#The_REDC_algorithm
 here R = 2^32 */
pub fn montgomery(a: i32) -> i16
{
    let mut m = a.wrapping_mul(Q_INV); 
    // a.wrapping_mul(b:T) returns (a * b) mod 2N, where N is the width of T in bits.
    m = (a as i32 - m * Q as i32) >> 16;
    m as i16
}


/* central reduction, returns r mod +- a */
pub fn cmod(r : i16, a: i16) -> i16
{
    let mut n = r % a;
    let mut t = a;
    if n > (t >> 1)
    {
        n -= t;
    }
    assert!(n <= (a+1)>>1);
    n 
}

pub fn _compress(x: i16, d: i16) -> i16
{
    let tmp: f32 = ((x<<d) as f32)/(Q as f32);
    let tmp: i16 = (0.000001 + tmp.trunc()) as i16;
    return cmod(tmp, 1<<d);
}

pub fn _decompress(x: i16, d: i16) -> i16
{
    let tmp: f32 = (((x*Q as i16) as f32)/((1<<d) as f32)) as f32;
    return (0.000001+tmp).ceil() as i16
}

