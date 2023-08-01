use crate::{params::*, utils::*, ntt::*};

/* polynomial ring (make a struct?) arithm op modulo q, ntt domain */

#[derive(Copy, Clone)]
pub struct Poly
{
    pub coeff: [i32; N],
    pub ntt: bool,
}

impl Default for Poly
{
  fn default() -> Self
  {
    Poly { coeff: [0i32; N] , ntt: false}
  }
}

pub fn neg(a: &Poly) -> Poly
{
    let mut res = Poly::default();
    for i in 0..N 
    {
        res.coeff[i] = cmod(-a.coeff[i], Q as i32);
    }
    res.ntt = a.ntt;
    res    
}

pub fn add(a: &mut Poly, b: &Poly)
{
    assert_eq!(a.ntt, b.ntt);
    if a.ntt == b.ntt
    {
        for i in 0..N 
        {
            a.coeff[i] = cmod(a.coeff[i] + b.coeff[i], Q as i32);
        }
    }
}

// pub fn sub(a: &mut Poly, b: &Poly)
// {
//     assert_eq!(a.ntt, b.ntt);
//     if a.ntt == b.ntt
//     {
//         for i in 0..N 
//         {
//             a.coeff[i] = a.coeff[i] - b.coeff[i];
//             a.coeff[i] = reduce32(a.coeff[i]);
//         }
//     }
// }


pub fn ntt(a: &Poly) -> Poly
{
    let mut res = Poly::default();
    res.coeff = a.coeff;
    if ! a.ntt
    {
        ntt_(&mut res.coeff);
        res.ntt = true;
    }
    res
}

pub fn reduce(a: &mut Poly)// -> i32
{
    for i in 0..N 
    {
        a.coeff[i] = cmod(montgomery(a.coeff[i] as i64), Q as i32);
       //a.coeff[i] = cmod(montgomery((a.coeff[i] as i64)<<16), Q as i32);
    }
}

pub fn inv_ntt(a: &Poly) -> Poly
{
    let mut res = Poly::default();
    res.coeff = a.coeff;
    if a.ntt
    {
        inv_ntt_(&mut res.coeff);
        reduce(&mut res);
        res.ntt = false;
    }
    res
}

/* scalar multiplication */
pub fn smult(s: i32, a: &Poly) -> Poly
{
    let mut res = Poly::default();
    for i in 0..N 
    {
        res.coeff[i] = cmod(s * a.coeff[i], Q as i32);
    }
    res.ntt = a.ntt;
    res
}

pub fn mult(a: &Poly, b: &Poly) -> Poly
{
    assert_eq!(a.ntt, b.ntt);
    assert_eq!(a.ntt, true);
    let mut res = Poly::default();
    for i in 0..N 
    {
        res.coeff[i] = montgomery((a.coeff[i] as i64) * b.coeff[i] as i64);
        let t = (((Q as i64 + res.coeff[i] as i64)<<32 ) % (Q as i64)) as i32; 
        res.coeff[i] = crate::utils::cmod(t, Q as i32); 
    }
    //reduce(& mut res);
    res.ntt = true;
    res
}

pub fn slow_mult(a: &Poly, b: &Poly) -> Poly
{
    assert_eq!(a.ntt, b.ntt);
    assert_eq!(a.ntt, false);
    let mut res = Poly::default();
    for i in 0..N 
    {
        for j in 0..N-i
        {
            res.coeff[i+j] += a.coeff[i] * b.coeff[j];           
        }
    } 
    for j in 1..N 
    {
        for i in N-j .. N 
        {
            res.coeff[i+j-N] -= a.coeff[i] * b.coeff[j];
        }
    }
    res.ntt = false;
    res
}

// pub fn poly_chknorm(a: &Poly, b: i32) -> u8
// {
//   // It is ok to leak which coefficient violates the bound since
//   // the probability for each coefficient is independent of secret
//   // data but we must not leak the sign of the centralized representative.
//   let mut t;

//   if b > ((Q as i32) - 1) / 8 {
//     return 1;
//   }
//   for i in 0..N {
//     // Absolute value of centralized representative
//     t = a.coeff[i] >> 31;
//     t = a.coeff[i] - (t & 2 * a.coeff[i]);

//     if t >= b {
//       return 1;
//     }
//   }
//   return 0;
// }

pub fn p_infnorm(p: &Poly) -> i32
{
    let mut norm = 0i32;
    for i in 0..N 
    {
        let t = p.coeff[i] - ((p.coeff[i] >> 31) & 2 * p.coeff[i]);
        if t >= norm 
        {
            norm = t;
           //norm = cmod(p.coeff[i], Q as i32).abs();
        }
    }
    norm
}
// pub fn p_infnorm(a: &Poly) -> usize
// {  
//     let mut t;
//     let mut norm = 0i32;

//     for i in 0..N {
//       // Absolute value of centralized representative
//       t = a.coeff[i] >> 31;
//       t = a.coeff[i] - (t & 2 * a.coeff[i]);
  
//       if t >= norm {
//         norm = t;
//       }
//     }
//     return norm as usize;
//   }


/* vector of polynomials */
#[derive(Copy, Clone)]
pub struct VecPoly<const l: usize>
{
    pub poly: [Poly; l],
}


impl<const l:usize> Default for VecPoly<{l}>
{
  fn default() -> Self
  {
    VecPoly { poly: [Poly::default(); l] }
  }
}


pub fn inf_norm<const k: usize>(p: &VecPoly<{k}>) -> i32
{
    let mut norm = 0i32;
    for i in 0..k
    {
        if p_infnorm(&p.poly[i]) >= norm 
        {
            norm = p_infnorm(&p.poly[i]);
        }
    }
    norm
}

/* matrix of polynomials */
#[derive(Copy, Clone)]
pub struct Mat<const k: usize, const l:usize>
{
    pub vec: [VecPoly<l>; k],
}


impl<const k: usize, const l:usize> Default for Mat<{k}, {l}>
{
  fn default() -> Self
  {
    Mat { vec: [VecPoly::<l>::default(); k] }
  }
}


pub fn Ntt<const k: usize>(a: &VecPoly<{k}>) -> VecPoly<{k}>
{
    let mut res = *a;//.copy();
    for i in 0..res.poly.len()
    {
        res.poly[i] = ntt(&a.poly[i]);
    }
    res
}


pub fn inv_Ntt<const k: usize>(a: &VecPoly<{k}>) -> VecPoly<{k}>
{
    let mut res = VecPoly::<{k}>::default();
    for i in 0..res.poly.len()
    {
        res.poly[i] = inv_ntt(&a.poly[i]);
    }
    res
}


pub fn Neg<const k: usize>(a: &VecPoly<{k}>) -> VecPoly<{k}>
{
    let mut res = VecPoly::<{k}>::default();
    for i in 0..k 
    {
        res.poly[i] = neg(&a.poly[i]);
    }
    res    
}

pub fn sMult<const k: usize>(s: i32, a: &VecPoly<{k}>) -> VecPoly<{k}>
{
    let mut res = VecPoly::<{k}>::default();
    for i in 0..k 
    {
        res.poly[i] = smult(s, &a.poly[i]);
    }
    res  
}


pub fn v_add<const k: usize>(a: &VecPoly<k>, b: &VecPoly<k>) -> VecPoly<k>
{
    let mut res = VecPoly::<k>::default();
    for i in 0..k 
    {
        res.poly[i] = a.poly[i];
        add(&mut res.poly[i], &b.poly[i]);
    }
    //Caddq(&res)
    res
}


pub fn v_sub<const k: usize>(a: &VecPoly<k>, b: &VecPoly<k>) -> VecPoly<k>
{
    let mut res = VecPoly::<k>::default();
    for i in 0..k 
    {
        res.poly[i] = a.poly[i];
        add(&mut res.poly[i], &smult(-1, &b.poly[i]));
    }
    res
}


pub fn p_mult_v<const k: usize>(p: &Poly, v: &VecPoly<{k}>) -> VecPoly<{k}>
{
    let mut res = VecPoly::<k>::default();
    for i in 0..k 
    {
        res.poly[i] = mult(&p, &v.poly[i]); 
    }
    Caddq(&res);

    res
}


/* scalar vector multiplication */
pub fn v_mult_v<const k:usize>(v: &VecPoly<{k}>, u: &VecPoly<{k}>) -> Poly
{
    let mut res = Poly::default();
    res.ntt = true;
    for i in 0..k 
    {
        add(&mut res, &mult(&v.poly[i], &u.poly[i])); // res += mult(v.poly[i], u.poly[i])
    }
    //c_addq(&res)
    res
}


/* matrix and vector multiplication */
pub fn m_mult_v<const k:usize, const l:usize>(A: &Mat<{k}, {l}>, s: &VecPoly<{l}>) -> VecPoly<{k}>
{
    let mut res = VecPoly::<{k}>::default();
    for i in 0..k 
    {
        res.poly[i] = v_mult_v(&A.vec[i], &s);
    }
    //Caddq(&res)
    res
}


/* utilities */
pub fn c_addq(a: &Poly) -> Poly
{
    let mut res = Poly::default();
    for i in 0..N 
    {
        res.coeff[i] = reduce32(cmod(a.coeff[i], Q as i32));
    }
    res.ntt = a.ntt;
    res
}


pub fn Caddq<const k: usize>(a: &VecPoly<k>) -> VecPoly<k>
{
    let mut res = VecPoly::<k>::default();
    for i in 0..k 
    {
        res.poly[i] = c_addq(&a.poly[i]);
    }
    res
}


pub fn power2round(r: &Poly) -> (Poly, Poly)
{
    let mut t1 = Poly::default();
    let mut t0 = Poly::default();
    for i in 0..r.coeff.len()
    {
        (t1.coeff[i], t0.coeff[i]) = power2round_(r.coeff[i]);
    }
    t1.ntt = r.ntt;
    t0.ntt = r.ntt;
    (t1, t0)
}


fn decompose(r: &Poly) -> (Poly, Poly)
{
    let mut t1 = Poly::default();
    let mut t0 = Poly::default();
    for i in 0..r.coeff.len()
    {
        (t1.coeff[i], t0.coeff[i]) = decompose_(r.coeff[i]);
    }
    t1.ntt = r.ntt;
    t0.ntt = r.ntt;
    (t1, t0)
}


pub fn highBits(r: &Poly) -> Poly
{
    let (r1, _) = decompose(r);
    r1
}


pub fn lowBits(r: &Poly) -> Poly
{
    let (_, r0) = decompose(&r);
    r0
}


pub fn makeHint(z: &Poly, r: &Poly) -> (Poly, usize)
{
  let mut cnt = 0usize;
  let mut h = Poly::default();
  for i in 0..N {
    h.coeff[i] = makeHint_(z.coeff[i], r.coeff[i]) as i32;
    cnt += h.coeff[i] as usize;
  }
  (h, cnt)
}

pub fn make_hint<const k: usize>(v0: &VecPoly<{k}>, v1: &VecPoly<{k}>)
  -> (VecPoly<{k}>, i32)
{
    let mut h = VecPoly::<{k}>::default();
  let mut s = 0i32;
  let mut cnt = 0usize;
  for i in 0..k 
  {
    (h.poly[i], cnt) = makeHint( &v0.poly[i], &v1.poly[i]);
     s += cnt as i32;
  }
  (h, s)
}

pub fn useHint(h: &Poly, r: &Poly) -> Poly
{
    let mut t = Poly::default();
    t.ntt = r.ntt;
    for i in 0..N 
    {
        t.coeff[i] = useHint_(h.coeff[i] as u8, r.coeff[i]);
    }
    t
}