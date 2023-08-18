#[cfg(feature = "mode2")]
mod mode2;
#[cfg(not(any(feature = "mode2", feature = "mode5")))]
mod mode3;
#[cfg(feature = "mode5")]
mod mode5;

#[cfg(feature = "mode2")]
pub use mode2::*;
#[cfg(not(any(feature = "mode2", feature = "mode5")))]
pub use mode3::*;
#[cfg(feature = "mode5")]
pub use mode5::*;






pub const MESSAGE: usize = 32; //max size of the message



pub const Q: usize = 8380417;
pub const N: usize = 256;
pub const D: usize = 13;
pub const ROOT: usize = 1753;

//pub const ALPHA: usize = 2*GAMMA2;

pub const SEEDBYTES: usize = 32;
pub const CRHBYTES: usize = 64; // 48???

pub const T1_BYTES: usize = 320;
pub const T0_BYTES: usize = 416;


pub const Z_BYTES: usize  = if cfg!(feature = "mode2") { 576 } else { 640 };
pub const W1_BYTES: usize = if cfg!(feature = "mode2") { 192 } else { 128 };
pub const S_BYTES: usize  = if cfg!(not(any(feature = "mode2", feature = "mode5"))) { 128 } else { 96 };

pub const SK_BYTES: usize = 3 * SEEDBYTES
+ L * S_BYTES
+ K * S_BYTES
+ K * T0_BYTES;
pub const PK_BYTES: usize = SEEDBYTES + K * T1_BYTES;

pub const SIGN_BYTES: usize = SEEDBYTES + L * Z_BYTES + OMEGA + K;