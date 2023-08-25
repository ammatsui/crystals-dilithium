
use crate::params::*;
/* helper functions from the paper, bits and hints */

const Q_INV: i32 = 58728449; // q^(-1) mod 2^32
const Q_i32: i32 = Q as i32;


/* for finite field element a with -2^{31}Q <= a <= Q*2^31,
 compute a' = a*2^{-32} (mod Q) such that -Q < a' < Q.
 https://en.wikipedia.org/wiki/Montgomery_modular_multiplication#The_REDC_algorithm
 here R = 2^32 */
pub fn montgomery(a: i64) -> i32
{
    let mut m = (a as i32).wrapping_mul(Q_INV) as i64; 
    // a.wrapping_mul(b:T) returns (a * b) mod 2N, where N is the width of T in bits.
    m = (a as i64 - m * Q as i64) >> 32;
    m as i32
}

/* centered reduce modulo Q */
pub fn reduce32(a: i32) -> i32
{
  let mut t = (a + (1 << 22)) >> 23;
  t = a - t * Q as i32;
  t
}

/* central reduction, returns r mod +- a */
pub fn cmod(rr : i32, a: i32) -> i32
{
    let r = ((rr % Q_i32) + Q_i32) % Q_i32;
    if a == Q as i32
    {
        reduce32(rr);
    }
    let mut n = r % a;
    let mut t = a;
    if n > (t >> 1)
    {
        n -= t;
    }
    assert!(n <= (a+1)>>1);
    n 
}


/* the straightforward bit-wise way to break up an element r = r1 · 2^d + r0 
 where r0 = r mod± 2^d and r1 = (r − r0)/2^d
 returns (r1, r0) */
 pub fn power2round_(rr: i32) -> (i32, i32)
 {
    let r = ((rr % Q_i32) + Q_i32) % Q_i32; 
    let mut r0 = cmod(r, 1<<D); 
    if r0 > (1<<(D-1))
    {
        r0 = r0 - (1<<(D-1)) ;
    }
    let r1 = (r - r0) >> D;
    assert_eq!(r, (r1 << D) + r0);
    assert!(r0 <= (1<<(D-1)));
    assert!(r1 >= 0);
    (r1, r0) 
 }

/* for finite field element r, compute high and low bits r1, r0 such
that r = r1*ALPHA + r0 with -ALPHA/2 < r0 <= ALPHA/2 except
if r1 = (Q-1)/ALPHA set r1 = 0  
ALPHA is an even divisor of Q-1 (either (q − 1)/88, (q − 1)/32, (q − 1)/3 in different modes)
returns (r1, r0) */
pub fn decompose_(rr: i32) -> (i32, i32)
{
    let r = ((rr % Q_i32) + Q_i32) % Q_i32; 
    let mut r0 = cmod(r, ALPHA as i32);
    if r0 > (ALPHA >> 1) as i32
    {
      r0 = r0 - ((ALPHA >> 1) as i32);
    }
    let mut r1 = (r - r0) / (ALPHA as i32);
    /* check it is correct */
    assert!(r1 >= 0);
    assert_eq!(r, r1*(ALPHA as i32) + r0); 
    /* fix */
    if r1 * (ALPHA as i32) == (Q - 1) as i32
    {
        r1 = 0; r0 -= 1;
    }
    assert!(r1 < (Q_i32 - 1)/(ALPHA as i32));
    assert!(cmod(r0, Q_i32).abs() <= (ALPHA >> 1) as i32);
    (r1, r0)
}


pub fn highBits_(r: i32) -> i32
{
    let (r1, _) = decompose_(r);
    r1
}


pub fn lowBits_(r: i32) -> i32
{
    let (_, r0) = decompose_(r);
    r0
}

// pub fn makeHint_(z: i32, r: i32) -> u8
// {
//   const GAMMA2_I32 : i32 = ALPHA as i32;
//   if z > GAMMA2_I32 || z < Q_i32-GAMMA2_I32 || (z == Q_i32-GAMMA2_I32 && r != 0) {
//     assert_eq!(useHint_(1, r), highBits_(z+r));
//     return 1;
//   }
//   assert_eq!(useHint_(0, r), highBits_(z+r));
//   return 0;
// }

pub fn makeHint_(z: i32, r: i32) -> u8
{
    let r1 = highBits_(r);
    let v1 = highBits_(r+z);
    if r1 == v1
    {
        assert_eq!(useHint_(0, r), highBits_(z+r));
        return 0;
    }
    assert_eq!(useHint_(1, r), highBits_(z+r));
    return 1
}


pub fn useHint_(h: u8, r: i32) -> i32
{
    let m = ((Q-1)/ ALPHA) as i32;
    let (r1, r0) = decompose_(r);
    if h == (1 as u8)
    {
        if r0 >= 0
        {
            return (r1+1 + m) % m;
        }
        return (r1-1 + m) % m;
    }
    return r1;
}


