use crate::{params::*, utils::montgomery};
/* ntt library: forward, inverse,  multiplication */

/* roots of unity in order for the forward ntt */
pub const ZETAS: [i16; 128] = [
  -1044, -758, -359, -1517, 1493, 1422, 287, 202, -171, 622, 1577, 182, 962, -1202, -1474, 1468,
  573, -1325, 264, 383, -829, 1458, -1602, -130, -681, 1017, 732, 608, -1542, 411, -205, -1571,
  1223, 652, -552, 1015, -1293, 1491, -282, -1544, 516, -8, -320, -666, -1618, -1162, 126, 1469,
  -853, -90, -271, 830, 107, -1421, -247, -951, -398, 961, -1508, -725, 448, -1065, 677, -1275,
  -1103, 430, 555, 843, -1251, 871, 1550, 105, 422, 587, 177, -235, -291, -460, 1574, 1653, -246,
  778, 1159, -147, -777, 1483, -602, 1119, -1590, 644, -872, 349, 418, 329, -156, -75, 817, 1097,
  603, 610, 1322, -1285, -1465, 384, -1215, -136, 1218, -1335, -874, 220, -1187, -1659, -1185,
  -1530, -1278, 794, -1510, -854, -870, 478, -108, -308, 996, 991, 958, -1460, 1522, 1628,
];


pub fn ntt_(a: &mut [i16])
/* arg: array of coefficients 
output in bit-reversed order */
{
    let mut j;
    let mut k = 0;
    let mut len = N/2; //128
    let (mut t, mut zeta);
  
    while len > 0 {
      let mut i = 0;
      while i < N {
        k += 1;
        zeta = ZETAS[k] as i32;
        j = i;
        while j < (i + len) {
          t = montgomery(zeta * a[j + len] as i32);
          a[j + len] = a[j] - t;
          a[j] += t;
          j += 1;
        }
        i = j + len;
      }
      len >>= 1;
    }
  }


/* inverse ntt (and mult by Montgomery factor r = 2^32) */
pub fn inv_ntt_(a: &mut [i16])
{
    let mut j;
    let mut k = N;
    let mut len = 1;
    let (mut t, mut zeta);
    const F: i32 = 1441; // mont^2/256
  
    while len < N/2 {
      let mut i = 0;
      while i < 256 {
        k -= 1;
        zeta = -ZETAS[k] as i32;
        j = i;
        while j < (i + len) {
          t = a[j];
          a[j] = t + a[j + len];
          a[j + len] = t - a[j + len];
          a[j + len] = montgomery(zeta * a[j + len] as i32);
          j += 1
        }
        i = j + len;
      }
      len <<= 1;
    }
    for j in 0..N {
      a[j] = montgomery(F * a[j] as i32);
    }
  }
  