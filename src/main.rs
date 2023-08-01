pub mod params;
pub mod ntt;
pub mod utils;
pub mod poly;
pub mod expand;
pub mod bit_packing;
pub mod dilithium;

use crate::params::*;
use crate::utils::*;

pub fn cmod(a: i32) -> i32
{
   let mut t = a >> 31;
   t = a - (t & 2 * a);
   t
}

fn main()
{
   // let z = (ALPHA>>1 ) as i32 ;
   // assert_eq!(crate::utils::cmod(z, Q as i32), cmod(z));
   // assert!(cmod(z).abs() <= (ALPHA>>1) as i32);
   // assert!(z <= (ALPHA>>1) as i32);
   //  let r = 2<<12;  //r < 2<<13 - 1
   //  let (r1,  r0) = decompose_(r);
   //  assert!(r0 <= (ALPHA>>1) as i32);
   //  assert!(cmod(r0).abs() <= (ALPHA>>1) as i32);
   //  if r0 > 0 
   //  {
   //    // println!("{}", r0);
   //    // println!("{}", z);
   //    assert!(-(ALPHA as i32)/2 < r0 + z);
   //    assert!(r0 + z <= ALPHA as i32);
   //  }
   //  else 
   //  {
   //    assert!(-(ALPHA as i32) <= r0 + z);
   //    assert!(r0 + z <= (ALPHA as i32)/2);
   //  }
   //  assert_eq!(useHint_(makeHint_(z, r), r), highBits_(z+r));


   let (pk, sk) = dilithium::keyGen();
   let sig = dilithium::sign(&sk, &pk[..10]);
   let b = dilithium::verify(&pk, &pk[..10], &sig);
   assert_eq!(b, true);
}
