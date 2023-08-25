use crate::{params::*, poly::*, utils::*, packing::*};
use sha3::{Shake128, Shake256, Sha3_256, Sha3_512, Digest, };


pub fn xof(seed: &[u8], out: &mut [u8])
{
    use sha3::digest::{Update, ExtendableOutput, XofReader};

    let mut hasher = Shake128::default();
    hasher.update(&seed);
    let mut reader = hasher.finalize_xof();
    reader.read(out);
}

pub fn H(seed: &[u8], out: &mut [u8])
{
    let mut hasher = Sha3_256::default();
    hasher.update(&seed);
    let res = hasher.finalize(); 
    out.copy_from_slice(&res[..]);
}

pub fn G(seed: &[u8], out: &mut [u8])
{
    let mut hasher = Sha3_512::default();
    hasher.update(&seed);
    let mut res = hasher.finalize();
    out.copy_from_slice(&res[..]);
}

pub fn prf(seed: &[u8], out: &mut [u8])
{
    use sha3::digest::{Update, ExtendableOutput, XofReader};

    let mut hasher = Shake256::default();
    hasher.update(&seed);
    let mut reader = hasher.finalize_xof();
    reader.read(out);
}


pub fn parse(b: &[u8]) -> Poly
{
    let mut res = Poly::default();
    res.ntt = true;
    let mut i = 0;
    let mut j: usize = 0;
    while j < N 
    {
        let mut d: i16 = (b[i] as i16) + 256*(b[i+1] as i16);
        d = cmod(d, 1<<13);
        if d < 19*Q as i16
        {
            res.coeff[j] = d;
            j += 1;
        }
        i += 2;
    }
    res
}

/* load 4 bytes into a 32-bit integer in little-endian order */
fn load32_LE(x: &[u8]) -> u32 {
    let mut r = x[0] as u32;
    r |= (x[1] as u32) << 8;
    r |= (x[2] as u32) << 16;
    r |= (x[3] as u32) << 24;
    r
}

/* load 3 bytes into a 32-bit integer in little-endian order */
fn load24_LE(x: &[u8]) -> u32 {
    let mut r = x[0] as u32;
    r |= (x[1] as u32) << 8;
    r |= (x[2] as u32) << 16;
    r
}


fn cbd2(r: &mut Poly, buf: &[u8]) 
{
    let (mut d, mut t, mut a, mut b);
    for i in 0..(N / 8) {
        t = load32_LE(&buf[4 * i..]);
        d = t & 0x55555555;
        d += (t >> 1) & 0x55555555;
        for j in 0..8 {
            a = ((d >> (4 * j)) & 0x3) as i16;
            b = ((d >> (4 * j + 2)) & 0x3) as i16;
            r.coeff[8 * i + j] = a - b;
        }
    }
}



pub fn cbd3(r: &mut Poly, buf: &[u8]) {
    let (mut d, mut t, mut a, mut b);
    for i in 0..(N / 4) {
        t = load24_LE(&buf[3 * i..]);
        d = t & 0x00249249;
        d += (t >> 1) & 0x00249249;
        d += (t >> 2) & 0x00249249;
        for j in 0..4 {
            a = ((d >> (6 * j)) & 0x7) as i16;
            b = ((d >> (6 * j + 3)) & 0x7) as i16;
            r.coeff[4 * i + j] = a - b;
        }
    }
}


pub fn CBD(b: &[u8], ETA: usize) -> Poly
{
    let mut res = Poly::default();
    res.ntt = false;
    if ETA == 2 
    {
        cbd2(&mut res, b);
    }
    cbd3(&mut res, b);
    res
}


pub fn get_matrix(rho: &[u8]) -> Mat<K, K>
{
    let mut A: Mat<K, K> = Mat::<K, K>::default();
    for i in 0..K
    {
        for j in 0..K 
        {
            let mut tmp = [0u8; 32];
            let mut seed = [0u8; 32+4];
            seed[..32].copy_from_slice(&rho[..32]);
            seed[32..].copy_from_slice(&[j as u8, (j>>8) as u8]);
            seed[34..].copy_from_slice(&[i as u8, (i>>8) as u8]);
            xof(&seed, &mut tmp);
            A.vec[i].poly[j] = parse(&tmp);
        }
    }
    A
}

pub fn get_noise(sigma: &[u8]) -> VecPoly<K>
{
    let mut s = VecPoly::<K>::default();
    for i in 0..K
    {
        let mut tmp = [0u8; 32];
        let mut seed = [0u8; 32+2];
        seed[..32].copy_from_slice(&sigma[..32]);
        seed[32..].copy_from_slice(&[i as u8, (i>>8) as u8]);
        prf(&seed, &mut tmp);
        s.poly[i] = CBD(&tmp, ETA1);
    }
    s
}


pub fn get_error(sigma: &[u8]) -> VecPoly<K>
{
    let mut s = VecPoly::<K>::default();
    for i in 0..K
    {
        let mut tmp = [0u8; 32];
        let mut seed = [0u8; 32+2];
        seed[..32].copy_from_slice(&sigma[..32]);
        seed[32..].copy_from_slice(&[(i+K) as u8, ((i+K)>>8) as u8]);
        prf(&seed, &mut tmp);
        s.poly[i] = CBD(&tmp, ETA1);
    }
    s
}