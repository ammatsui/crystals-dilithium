#[cfg(feature = "mode512")]
mod mode512;
#[cfg(not(any(feature = "mode512", feature = "mode1024")))]
mod mode768;
#[cfg(feature = "mode1024")]
mod mode1024;

#[cfg(feature = "mode512")]
pub use mode512::*;
#[cfg(not(any(feature = "mode512", feature = "mode1024")))]
pub use mode768::*;
#[cfg(feature = "mode1024")]
pub use mode1024::*;



pub const MESSAGE: usize = 32; //max size of the message


pub const Q: usize = 3329; 
pub const N: usize = 256;


pub const KYBER_SYMBYTES: usize = 32;
pub const SEEDBYTES: usize = 32;
pub const PACKED_KEYS: usize = 12; 

/// Size of the shared key
pub const KYBER_SSBYTES: usize = 32;

pub const KYBER_POLYBYTES: usize = 384;
pub const KYBER_POLYVECBYTES: usize = K * KYBER_POLYBYTES;

#[cfg(not(feature = "kyber1024"))]
pub const KYBER_POLYCOMPRESSEDBYTES: usize = 128;
#[cfg(not(feature = "kyber1024"))]
pub const KYBER_POLYVECCOMPRESSEDBYTES: usize = K * 320;

#[cfg(feature = "kyber1024")]
pub const KYBER_POLYCOMPRESSEDBYTES: usize = 160;
#[cfg(feature = "kyber1024")]
pub const KYBER_POLYVECCOMPRESSEDBYTES: usize = K * 352;

pub const KYBER_INDCPA_PUBLICKEYBYTES: usize = KYBER_POLYVECBYTES + KYBER_SYMBYTES;
pub const KYBER_INDCPA_SECRETKEYBYTES: usize = KYBER_POLYVECBYTES;
pub const KYBER_INDCPA_BYTES: usize = KYBER_POLYVECCOMPRESSEDBYTES + KYBER_POLYCOMPRESSEDBYTES;

/// Size in bytes of the Kyber public key
pub const PK_BYTES: usize = KYBER_INDCPA_PUBLICKEYBYTES;
/// Size in bytes of the Kyber secret key
pub const SK_BYTES: usize =
    KYBER_INDCPA_SECRETKEYBYTES + KYBER_INDCPA_PUBLICKEYBYTES + 2 * KYBER_SYMBYTES;
/// Size in bytes of the Kyber ciphertext
pub const KYBER_CIPHERTEXTBYTES: usize = KYBER_INDCPA_BYTES;
