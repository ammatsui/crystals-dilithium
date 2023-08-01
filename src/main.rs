pub mod params;
pub mod ntt;
pub mod utils;
pub mod poly;
pub mod expand;
pub mod bit_packing;
pub mod dilithium;

use sha3::{Shake128, Shake256, digest::{Update, ExtendableOutput, XofReader}};
use crate::params::*;



fn main()
{
   // let rho = [121, 251, 146, 198, 80, 156, 206, 181, 87, 3, 96, 206, 220, 19, 239, 221, 81, 50, 248, 130, 149, 192, 182, 73, 237, 144, 243, 178, 29, 140, 187, 220];
   // let nonce = 4;
   //  let mut seed = [0u8; SEEDBYTES + 2];
   //  seed[..SEEDBYTES].copy_from_slice(&rho[0..SEEDBYTES]);
   //  seed[SEEDBYTES] = nonce as u8;
   //  seed[SEEDBYTES+1] = (nonce>>8) as u8;
   //  let mut hasher = Shake128::default();
   //  hasher.update(&seed);
   //  let mut reader = hasher.finalize_xof();
   //  let mut b = [0u8; 32];
   //  reader.read(&mut b);


   //  let mut seed = [0u8; SEEDBYTES + 2];
   //  seed[..SEEDBYTES].copy_from_slice(&rho[0..SEEDBYTES]);
   //  seed[SEEDBYTES] = nonce as u8;
   //  seed[SEEDBYTES+1] = (nonce>>8) as u8;
   //  let mut hasher = Shake128::default();
   //  hasher.update(&seed);
   //  let mut reader = hasher.finalize_xof();
   //  let mut bb = [0u8; 32];
   //  reader.read(&mut bb);

   //  assert_eq!(b, bb);

   let (pk, sk) = dilithium::keyGen();
   let sig = dilithium::sign(&sk, &pk[..10]);
   let b = dilithium::verify(&pk, &pk[..10], &sig);
   assert_eq!(b, true);
}
