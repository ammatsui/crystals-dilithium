pub mod params;
pub mod ntt;
pub mod utils;
pub mod poly;
pub mod expand;
pub mod bit_packing;
pub mod dilithium;




fn main()
{

   let (pk, sk) = dilithium::keyGen();
   let sig = dilithium::sign(&sk, &pk[..10]);
   let b = dilithium::verify(&pk, &pk[..10], &sig);
   assert_eq!(b, true);
   let sig = dilithium::sign(&sk, &pk[..12]);
   let b = dilithium::verify(&pk, &pk[..10], &sig);
   assert_eq!(b, false);
}
