use crate::{params::*, poly::*,};
// use bitvec::prelude::*;

/* encode */
pub fn encode(r: &mut [u8], a: Poly) 
{
    let (mut t0, mut t1);

    for i in 0..(N / 2) {
        t0 = a.coeff[2 * i];
        t0 += (t0 >> 15) & Q as i16;
        t1 = a.coeff[2 * i + 1];
        t1 += (t1 >> 15) & Q as i16;
        r[3 * i + 0] = (t0 >> 0) as u8;
        r[3 * i + 1] = ((t0 >> 8) | (t1 << 4)) as u8;
        r[3 * i + 2] = (t1 >> 4) as u8;
    }
}


/* decode */
pub fn decode(r: &mut Poly, a: &[u8]) {
    for i in 0..(N / 2) {
        r.coeff[2 * i + 0] =
            ((a[3 * i + 0] >> 0) as u16 | ((a[3 * i + 1] as u16) << 8) & 0xFFF) as i16;
        r.coeff[2 * i + 1] =
            ((a[3 * i + 1] >> 4) as u16 | ((a[3 * i + 2] as u16) << 4) & 0xFFF) as i16;
    }
}


/* main packers */
/* secret key sk = (s)*/
pub fn pack_sk(sk: &mut [u8], s: &VecPoly<K>)
{
    for i in 0..K 
    {
        encode(&mut sk[ i * PACKED_KEYS..], s.poly[i]);
    }
}


pub fn unpack_sk( sk: &[u8], s: &mut VecPoly<K>)
{
    for i in 0..K 
    {
        decode( &mut s.poly[i], &sk[i * PACKED_KEYS..],);
    }
  }
  


/* public key pk = (t, rho) */
pub fn pack_pk(pk: &mut [u8], t: &VecPoly<K>, rho: &[u8])
{
    pk[..SEEDBYTES].copy_from_slice(&rho[..SEEDBYTES]);
    for i in 0..K 
    {
        encode(&mut pk[SEEDBYTES + i * PACKED_KEYS..], t.poly[i]);
    }
}

/* unpack public key pk = (t, rho) */
pub fn unpack_pk(pk: &[u8], t: &mut VecPoly<K>, rho: &mut [u8])
{
    rho[..SEEDBYTES].copy_from_slice(&pk[..SEEDBYTES]);
    for i in 0..K 
    {
        decode( &mut t.poly[i], &pk[SEEDBYTES + i * PACKED_KEYS..],);
    }
}


// /* pack signature sig = (c_hat, z, h) */
// pub fn pack_sign(sig: &mut [u8], c_hat: &[u8], z: &VecPoly<L>, h: &VecPoly<K>)
// {
//     let mut ctr = 0usize;
    
//     sig[..SEEDBYTES].copy_from_slice(&c_hat[..SEEDBYTES]);
//     ctr += SEEDBYTES;

//     for i in 0..L 
//     {
//         pack_z(&z.poly[i], &mut sig[ctr + i*Z_BYTES..]);
//     }
//     ctr += L * Z_BYTES;

//     /* pack h */
//     sig[ctr..ctr + OMEGA + K].copy_from_slice(&[0u8; OMEGA + K]);

//   let mut k = 0;
//   for i in 0..K {
//     for j in 0..N {
//       if h.poly[i].coeff[j] != 0 {
//         sig[ctr + k] = j as u8;
//         k += 1;
//       }
//     }
//     sig[ctr + OMEGA + i] = k as u8;
//     }
// }


// pub fn unpack_sign(sig: &[u8], c_hat: & mut [u8], z: &mut VecPoly<L>, h: &mut VecPoly<K>)
// {
//     let mut ctr = 0usize;
    
//     c_hat[..SEEDBYTES].copy_from_slice(&sig[..SEEDBYTES]);
//     ctr += SEEDBYTES;

//     for i in 0..L 
//     {
//         unpack_z( & sig[ctr + i*Z_BYTES..], &mut z.poly[i]);
//     }
//     ctr += L * Z_BYTES;

//     /* unpack h */
//     let mut k = 0usize;
//   for i in 0..K {
//     if sig[ctr + OMEGA + i] < k as u8 || sig[ctr + OMEGA + i] > (OMEGA as u8) {
//       return ;
//     }
//     for j in k..sig[ctr + OMEGA + i] as usize {
//       // Coefficients are ordered for strong unforgeability
//       if j > k && sig[ctr + j as usize] <= sig[ctr + j as usize - 1] {
//         return ;
//       }
//       h.poly[i].coeff[sig[ctr + j] as usize] = 1;
//     }
//     k = sig[ctr + OMEGA + i] as usize;
//   }

//   // Extra indices are zero for strong unforgeability
//   for j in k..OMEGA {
//     if sig[ctr + j as usize] > 0 {
//       return ;
//     }
//   }
// }