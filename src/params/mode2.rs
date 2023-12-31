/* mode 2 */

use super::Q;

pub const K: usize = 4;
pub const L: usize = 4;
pub const ETA: usize = 2;
pub const TAU: usize = 39;
pub const BETA: usize = 78;
pub const GAMMA1: usize = 1 << 17;
pub const GAMMA2: usize = (Q - 1) / 88;
pub const ALPHA: usize = 2*GAMMA2;
pub const OMEGA: usize = 80;

pub const g1_bits: usize = 18;