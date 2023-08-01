use crate::{params::*, poly::*,};
use bitvec::prelude::*;


/* general packer */
/* pack integers into bytestring with `n_bits` per integer */
pub fn pack(st: &mut [u8], n_bits: usize, coef: & [i32])
{
    for (slot, c) in st
  .view_bits_mut::<Lsb0>()
  .chunks_mut(n_bits)
  .zip(coef.iter().copied())
{
  slot.store_le(c as u32); // because the coeffs are supposed to be modified to be positive
  assert_eq!(slot.load_le::<u32>(), c as u32); 
}
}

pub fn unpack(st: &[u8], n_bits: usize, coef: &mut [i32])
{
    let bits = st.view_bits::<Lsb0>();
    for i in 0..N
    {
        coef[i] = (bits[n_bits*i..n_bits*(i+1)]).load_le::<u32>() as i32;
    }
}

/* helper packers of intermediate polys */

/* pack poly of vectors s1, s2 */
pub fn pack_s(s: &Poly, out: &mut [u8])
{
    let mut coeffs = [0i32; N];
    for i in 0..coeffs.len()
    {
        coeffs[i] = (ETA as i32) - s.coeff[i];
    }
    if ETA == 2
    {
        pack( out, 3, &coeffs);
    }
    if ETA == 4
    {
        pack( out, 4, &coeffs);
    }
}

pub fn unpack_s(st: &[u8], s: &mut Poly)
{
    if ETA == 2
    {
        unpack(&st, 3, &mut s.coeff);
    }
    if ETA == 4
    {
        unpack(&st, 4, &mut s.coeff);
    }
    for i in 0..s.coeff.len()
    {
        s.coeff[i] = (ETA as i32) - s.coeff[i];
    }
    s.ntt = false;
}

pub fn pack_t0(s: &Poly, out: &mut [u8])
{
    let mut coeffs = [0i32; N];
    for i in 0..coeffs.len()
    {
        coeffs[i] = (1<<(D-1)) - s.coeff[i];
    } 
    pack( out, 13, &coeffs);
}

pub fn unpack_t0(st: &[u8], s: &mut Poly)
{
    unpack(st, 13, &mut s.coeff);
    for i in 0..s.coeff.len()
    {
        s.coeff[i] = (1<<(D-1)) - s.coeff[i];
    }
    s.ntt = false;
}

pub fn pack_t1(s: &Poly, out: &mut [u8])
{
    pack( out, 10, &s.coeff);
}

pub fn unpack_t1(st: &[u8], s: &mut Poly)
{
    unpack(st, 10, &mut s.coeff);
    s.ntt = false;
}


pub fn pack_w(s: &Poly, out: &mut [u8])
{
    if GAMMA2 == (Q - 1) / 88
    {
        pack( out, 6, &s.coeff);
    }
    if GAMMA2 == (Q - 1) / 32
    {
        pack( out, 4, &s.coeff);
    }
}

pub fn unpack_w(st: &[u8], s: &mut Poly)
{
    if GAMMA2 == (Q - 1) / 88
    {
        unpack(st, 6, &mut s.coeff);
    }
    if GAMMA2 == (Q - 1) / 32
    {
        unpack(st, 4, &mut s.coeff);
    }
    s.ntt = false;
}


pub fn pack_z(s: &Poly, out: &mut [u8])
{
    let mut coeffs = [0i32; N];
    for i in 0..coeffs.len()
    {
        coeffs[i] = (GAMMA1 as i32) - s.coeff[i];
    }
    if GAMMA1 == (1<<17)
    {
        pack(out, 18, &coeffs);
    }
    if GAMMA1 == (1<<19)
    {
        pack(out, 20, &coeffs);
    }
}

pub fn unpack_z(st: &[u8], s: &mut Poly)
{
    if GAMMA1 == (1<<17)
    {
        unpack(st, 18, &mut s.coeff);
    }
    if GAMMA1 == (1<<19)
    {
        unpack(st, 20, &mut s.coeff);
    }
    for i in 0..s.coeff.len()
    {
        s.coeff[i] = (GAMMA1 as i32) - s.coeff[i];
    }
    s.ntt = false;
}



/* main packers */
/* secret key sk = (rho, K, tr, s1, s2, t0)*/
pub fn pack_sk(sk: &mut [u8],
               rho: &[u8],
               tr: &[u8],
               key: &[u8],
               t0: &VecPoly<K>,
               s1: &VecPoly<L>,
               s2: &VecPoly<K>)
{
    let mut ctr = 0usize;
    
    sk[ctr..SEEDBYTES].copy_from_slice(&rho[0..SEEDBYTES]);
    ctr += SEEDBYTES;

    sk[ctr..ctr + SEEDBYTES].copy_from_slice(&key[0..SEEDBYTES]);
    ctr += SEEDBYTES;

    sk[ctr..ctr + SEEDBYTES].copy_from_slice(&tr[0..SEEDBYTES]);
    ctr += SEEDBYTES;

    for i in 0..L 
    {
        pack_s(&s1.poly[i], &mut sk[ctr + i * S_BYTES..]);
    }
    ctr += L * S_BYTES;

    for i in 0..K 
    {
        pack_s(&s2.poly[i], &mut sk[ctr + i * S_BYTES..]);
    }
    ctr += K * S_BYTES;

    for i in 0..K 
    {
      pack_t0(&t0.poly[i], &mut sk[ctr + i * T0_BYTES..]);
    }
}


pub fn unpack_sk( sk: &[u8],
                  rho: &mut [u8],
                  tr: &mut [u8],
                  key: &mut [u8],
                  t0: &mut VecPoly<K>,
                  s1: &mut VecPoly<L>,
                  s2: &mut VecPoly<K>)
{
    let mut ctr = 0usize;
  
    rho[..SEEDBYTES].copy_from_slice(&sk[..SEEDBYTES]);
    ctr += SEEDBYTES;
  
    key[..SEEDBYTES].copy_from_slice(&sk[ctr..ctr + SEEDBYTES]);
    ctr += SEEDBYTES;
  
    tr[..SEEDBYTES].copy_from_slice(&sk[ctr..ctr + SEEDBYTES]);
    ctr += SEEDBYTES;
  
    for i in 0..L 
    {
       unpack_s(&sk[ctr + i * S_BYTES..], &mut s1.poly[i], );
    }
    ctr += L * S_BYTES;
  
    for i in 0..K
    {
       unpack_s(&sk[ctr + i * S_BYTES..], &mut s2.poly[i], );
    }
    ctr += K * S_BYTES;
  
    for i in 0..K 
    {
        unpack_t0(&sk[ctr + i * T0_BYTES..], &mut t0.poly[i]);
    }
  }
  


/* public key pk = (rho, t1) */
pub fn pack_pk(pk: &mut [u8], rho: &[u8], t1: &VecPoly<K>)
{
    pk[..SEEDBYTES].copy_from_slice(&rho[..SEEDBYTES]);
    for i in 0..K 
    {
        pack_t1(&t1.poly[i], &mut pk[SEEDBYTES + i * T1_BYTES..]);
    }
}

/* unpack public key pk = (rho, t1) */
pub fn unpack_pk(pk: &[u8], rho: &mut [u8], t1: &mut VecPoly<K>)
{
    rho[..SEEDBYTES].copy_from_slice(&pk[..SEEDBYTES]);
    for i in 0..K 
    {
        unpack_t1(&pk[SEEDBYTES + i * T1_BYTES..], &mut t1.poly[i]);
    }
}


/* pack signature sig = (c_hat, z, h) */
pub fn pack_sign(sig: &mut [u8], c_hat: &[u8], z: &VecPoly<L>, h: &VecPoly<K>)
{
    let mut ctr = 0usize;
    
    sig[..SEEDBYTES].copy_from_slice(&c_hat[..SEEDBYTES]);
    ctr += SEEDBYTES;

    for i in 0..L 
    {
        pack_z(&z.poly[i], &mut sig[ctr + i*Z_BYTES..]);
    }
    ctr += L * Z_BYTES;

    /* pack h */
    sig[ctr..ctr + OMEGA + K].copy_from_slice(&[0u8; OMEGA + K]);

  let mut k = 0;
  for i in 0..K {
    for j in 0..N {
      if h.poly[i].coeff[j] != 0 {
        sig[ctr + k] = j as u8;
        k += 1;
      }
    }
    sig[ctr + OMEGA + i] = k as u8;
    }
}


pub fn unpack_sign(sig: &[u8], c_hat: & mut [u8], z: &mut VecPoly<L>, h: &mut VecPoly<K>)
{
    let mut ctr = 0usize;
    
    c_hat[..SEEDBYTES].copy_from_slice(&sig[..SEEDBYTES]);
    ctr += SEEDBYTES;

    for i in 0..L 
    {
        unpack_z( & sig[ctr + i*Z_BYTES..], &mut z.poly[i]);
    }
    ctr += L * Z_BYTES;

    /* unpack h */
    let mut k = 0usize;
  for i in 0..K {
    if sig[ctr + OMEGA + i] < k as u8 || sig[ctr + OMEGA + i] > (OMEGA as u8) {
      return ;
    }
    for j in k..sig[ctr + OMEGA + i] as usize {
      // Coefficients are ordered for strong unforgeability
      if j > k && sig[ctr + j as usize] <= sig[ctr + j as usize - 1] {
        return ;
      }
      h.poly[i].coeff[sig[ctr + j] as usize] = 1;
    }
    k = sig[ctr + OMEGA + i] as usize;
  }

  // Extra indices are zero for strong unforgeability
  for j in k..OMEGA {
    if sig[ctr + j as usize] > 0 {
      return ;
    }
  }
}