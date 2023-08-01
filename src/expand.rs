use crate::{params::*, poly::*, bit_packing::*};
use sha3::{Shake128, Shake256, digest::{Update, ExtendableOutput, XofReader}};


// todo: fix rejection sampling where needed
// it is supposed to be like in sampleInBall


pub fn crh(seed: &[u8], xof: &mut [u8])
{
    let mut hasher = Shake256::default();
    hasher.update(&seed);
    let mut reader = hasher.finalize_xof();
    reader.read(xof);
}


/* generate polynomials for the matrix */
fn p_uniform(a: &mut Poly, rho: &[u8], nonce: u16)
{  
    /* convert nonce to two u8, add to rho */
    let mut seed = [0u8; SEEDBYTES + 2];
    seed[..SEEDBYTES].copy_from_slice(&rho[0..SEEDBYTES]);
    seed[SEEDBYTES] = nonce as u8;
    seed[SEEDBYTES+1] = (nonce>>8) as u8;
    let mut hasher = Shake128::default();
    hasher.update(&seed);
    let mut reader = hasher.finalize_xof();
    /*rejection sample 256 coefficients, each is of the form b2' << 16 + b1 << 8 + b0 */
    let mut t: u32;
    let mut i = 0usize;
    while i < N
    {
        let mut b = [0u8; 3];
        reader.read(&mut b);
        t = b[0] as u32;
        t |= (b[1] as u32) << 8;
        t |= (b[2] as u32) << 16;
        t &= 0x7FFFFF;

        if t < Q as u32
        {
            a.coeff[i] = t as i32;
            i += 1;
        }
    }
    a.ntt = true;
}

/* compute an integer from bitvector */
// fn from_bits(bits: &BitVec, start: usize, end: usize) -> i32
// {
//     let mut t = 0usize;
//     for i in 0..end-start
//     {
//         if bits[start+i]
//         {
//             t |= 1 << (end-start-1-i);
//         }
//     }
//     t as i32
// }


/* generate coefficients for the mask polynomial */
fn poly_mask(a: &mut Poly, rho: &[u8], nonce: u16)
{
    let mut seed = [0u8; CRHBYTES + 2]; 
    seed[..CRHBYTES].copy_from_slice(&rho[0..CRHBYTES]);
    seed[CRHBYTES] = nonce as u8;
    seed[CRHBYTES+1] = (nonce>>8) as u8;
    let mut b = [0; N*g1_bits/8];
    crh(&seed, &mut b);
    // let mut hasher = Shake256::default();
    // hasher.update(&seed);
    // let mut reader = hasher.finalize_xof();
    /* read all the data */
   // reader.read(&mut b);
    /* no rejection sampling */
    unpack(&b, g1_bits, &mut a.coeff);
    for i in 0..N
    {
        a.coeff[i] = (GAMMA1 as i32) - a.coeff[i];
    }
    a.ntt = false;
}

/* generate coefficients for the secret error polynomials */
fn poly_s(s: &mut Poly, rho: &[u8], nonce: u16)
{
    const ETA_i32 : i32 = ETA as i32;
    let mut seed = [0u8; CRHBYTES + 2]; 
    seed[..CRHBYTES].copy_from_slice(&rho[0..CRHBYTES]);
    seed[CRHBYTES] = nonce as u8;
    seed[CRHBYTES+1] = (nonce>>8) as u8;
    let mut hasher = Shake256::default();
    hasher.update(&seed);
    let mut reader = hasher.finalize_xof();
    /* rejection sample coefficients, one byte is interpreted as 2 integers */
    let mut t: u32;
    let mut i = 0usize;
    while i < N
    {
        let mut b = [0u8; 1];
        reader.read(&mut b);
        let mut t0 = ((b[0] & 0x0F) as u32) as i32;
        let mut t1 = ((b[0] >> 4) as u32) as i32;
        /* rejection */
        if ETA == 2
        {
            if t0 < 15
            {
                s.coeff[i] = ETA_i32 - t0 % 5;
                i+= 1;
            }
            if (t1 < 15 && i < N)
            {
                s.coeff[i] = ETA_i32 - t1 % 5;
                i+= 1;
            }
        }
        if ETA == 4
        {
            if t0 < 9
            {
                s.coeff[i] = ETA_i32 - t0;
                i+= 1;
            }
            if (t1 < 9 && i < N)
            {
                s.coeff[i] = ETA_i32 - t1;
                i+= 1;
            }
        }
    }
    s.ntt = false;
}


/* computes each coefficient aˆi,j ∈ Rq of Aˆ separately. 
For the coefficient aˆi,j it absorbs the 32 bytes of ρ immediately followed by
two bytes representing 0 ≤ 256 ∗ i + j < 216 in little-endian byte order into SHAKE-128.
Then ExpandA performs rejection sampling on these 23-bit integers to sample the 256 coefficients of A^ */
pub fn expandA(rho: &[u8])  -> Mat<K, L>// returns matrix k x l with coefficients of random polys in NTT domain
{
    let mut A: Mat<K, L> = Mat::<K, L>::default();
    for i in 0..K
    {
        for j in 0..L 
        {
            p_uniform(&mut A.vec[i].poly[j], rho, ((i<<8) + j) as u16);
        }
    }
    A
}


pub fn expandMask(rho: &[u8], kappa: u16) -> VecPoly<L>
{
    let mut y = VecPoly::<L>::default();
    for i in 0..L 
    {
        poly_mask(&mut y.poly[i], rho, kappa + i as u16);
    }
    y
}



pub fn expandS(rho: &[u8]) -> (VecPoly<L>, VecPoly<K>)
{
    let mut s1 : VecPoly<L> = VecPoly::<L>::default();
    let mut s2 : VecPoly<K> = VecPoly::<K>::default();
    for i in 0..K+L 
    {
        if i < L 
        {
            poly_s(&mut s1.poly[i], rho, i as u16);
        }
        else 
        {
            poly_s(&mut s2.poly[i-L], rho, i as u16);
        }
    }
    (s1, s2)
}


/* create a random N-element array with τ ±1’s and (256 − τ) 0′s using 
the input seed ρ (and an SHAKE256) to generate the randomness needed */
pub fn SampleInBall(c_hat: &[u8]) -> Poly
{
    let mut c = Poly::default();
    let mut buf = [0u8; N*SEEDBYTES]; //why
    let mut hasher = Shake256::default();
    hasher.update(&c_hat);
    let mut reader = hasher.finalize_xof();
    reader.read(&mut buf);
    let mut _signs = 0u64;
    /* generate tau +- one's */
    for i in 0..8 
    {
        _signs |= (buf[i] as u64) << 8 * i;
    }
    let mut ctr: usize = 8;
    let mut j;
    for i in N - TAU..N 
    {
        /* rejection sampling j <= i*/
        loop // infinite loop
        {
            /* sample random bytes from hash output and
            interpret them as integers in {0, ..., 255}
            rejects values until a value j <= i is found */
            j = buf[ctr] as usize;
            ctr += 1;
            if j <= i 
            {
                break;
            }
        }
        /* save the coefficient */
        c.coeff[i] = c.coeff[j as usize];
        c.coeff[j as usize] = 1i32 - 2 * (_signs & 1) as i32;
        _signs >>= 1;
    }
    c
}