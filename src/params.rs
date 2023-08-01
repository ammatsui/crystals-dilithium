pub const MESSAGE: usize = 32; //max size of the message



pub const Q: usize = 8380417; //usize
pub const N: usize = 256;
pub const D: usize = 13;//32;//13;  //32
pub const ROOT: usize = 1753;

pub const SEEDBYTES: usize = 32;
pub const CRHBYTES: usize = 64; // 48


/* mode 3 */ 
pub const K: usize = 6;
pub const L: usize = 5;
pub const ETA: usize = 4;
pub const TAU: usize = 49;
pub const BETA: usize = 196;
pub const GAMMA1: usize = 1 << 19;
pub const GAMMA2: usize = (Q - 1) / 32;
pub const ALPHA: usize = 2*GAMMA2;
pub const OMEGA: usize = 55;

pub const g1_bits: usize = 20;


// TODO: add other modes
// add aes256 option (instead of shake-128)
pub const T1_BYTES: usize = 320;
pub const T0_BYTES: usize = 416;
//pub const POLYVECH_PACKEDBYTES: usize = OMEGA + K;

pub const Z_BYTES: usize = 640;
pub const W1_BYTES: usize = 192;
pub const S_BYTES: usize = 128;

pub const SK_BYTES: usize = 3 * SEEDBYTES
+ L * S_BYTES
+ K * S_BYTES
+ K * T0_BYTES;
pub const PK_BYTES: usize = SEEDBYTES + K * T1_BYTES;

pub const SIGN_BYTES: usize = SEEDBYTES + L * Z_BYTES + OMEGA + K;

// pub const Z_BYTES: usize =
//   if cfg!(feature = "mode2") { 576 } else { 640 };
// pub const W1_BYTES: usize =
//   if cfg!(feature = "mode2") { 192 } else { 128 };

// pub const S_BYTES: usize =
//   if cfg!(not(any(feature = "mode2", feature = "mode5"))) {
//     128
//   } else {
//     96
//   };