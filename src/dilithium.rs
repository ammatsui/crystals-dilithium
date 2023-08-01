use crate::{params::*, expand::*, poly::*, bit_packing::*};
use rand::Rng;

/* dilithium procedures */

pub fn keyGen() -> ([u8; PK_BYTES], [u8; SK_BYTES])
{
    let ksi  = rand::thread_rng().gen::<[u8; 32]>();
    
    /* initialisation */
    let mut sk = [0u8; SK_BYTES];
    let mut pk = [0u8; PK_BYTES];

    let (mut rho, mut rhoprime, mut key) =
    ([0u8; SEEDBYTES], [0u8; CRHBYTES], [0u8; SEEDBYTES]);

    let mut tr = [0u8; SEEDBYTES];

    /* generate seeds from hash of random ksi*/
    let mut seed = [0u8; 2 * SEEDBYTES + CRHBYTES];
    crh(&ksi, &mut seed);
    rho.copy_from_slice(&seed[..SEEDBYTES]);
    rhoprime.copy_from_slice(&seed[SEEDBYTES..SEEDBYTES + CRHBYTES]);
    key.copy_from_slice(&seed[SEEDBYTES + CRHBYTES..]);
    /* generation algorithm */
    let A = expandA( &rho);
    let (s1, s2) = expandS(&rhoprime);
    let mut t: VecPoly<K> = inv_Ntt(&m_mult_v(&A, &Ntt(&s1)));
    t = v_add(&t, &s2);
    let (mut t1, mut t0) = (VecPoly::<K>::default(), VecPoly::<K>::default());
    for i in 0..K 
    {
        (t1.poly[i], t0.poly[i]) = power2round(&t.poly[i]);
    }
    /* pack public key = rho || t1 */
    pack_pk(&mut pk, &rho, &t1);
    /* hash */
    crh(&pk, &mut tr);
    /* pack secret key */
    pack_sk(&mut sk, &rho, &tr, &key, &t0, &s1, &s2);
    
    (pk, sk)
}


pub fn sign(sk: &[u8], m: &[u8]) -> [u8; SIGN_BYTES] 
{
    /* initialisation */
    let mut sig = [0u8; SIGN_BYTES];
    let mut s1 : VecPoly<L> = VecPoly::<L>::default();
    let mut s2 : VecPoly<K> = VecPoly::<K>::default();
    let mut z  : VecPoly<L> = VecPoly::<L>::default();
    let (mut rho, mut rhoprime, mut key) =
    ([0u8; SEEDBYTES], [0u8; CRHBYTES], [0u8; SEEDBYTES]);
    let mut t0 = VecPoly::<K>::default();
    let mut tr = [0u8; SEEDBYTES];

    let mut mu = [0u8; CRHBYTES];

    let mut c_hat = [0u8; SEEDBYTES];
    
    /* unpack keys */
    unpack_sk(sk, &mut rho, &mut tr, &mut key, &mut t0, &mut s1, &mut s2);

    let A = expandA(&rho);

    /* CRH(tr+message) */
    let mut seed = [0u8; MESSAGE+SEEDBYTES];
    seed[..SEEDBYTES].copy_from_slice(&tr[..]);
    seed[SEEDBYTES..SEEDBYTES+m.len()].copy_from_slice(&m[..m.len()]);
    crh(&seed[..m.len()+SEEDBYTES], &mut mu);


    /* CRH(key + CRH(tr+message)) */
    let mut seed = [0u8; SEEDBYTES+CRHBYTES];
    seed[..SEEDBYTES].copy_from_slice(&key[..]);
    seed[SEEDBYTES..].copy_from_slice(&mu[..]);
    crh(&seed, &mut rhoprime);
   
    s1 = Ntt(&s1);
    s2 = Ntt(&s2);
    t0 = Ntt(&t0);

    let mut kappa = 0usize;

    loop 
    {
        let y = expandMask(&rhoprime, kappa as u16);
        kappa += 1;

        let w =  inv_Ntt(&m_mult_v(&A, &Ntt(&y)));
        /* highbits */
        let mut w1 = VecPoly::<K>::default();
        for i in 0..K
        {
            w1.poly[i] = highBits(&w.poly[i]);
        }

        /* pack w1 */
        let mut wp = [0u8; W1_BYTES*K];
        for i in 0..K
        {
            pack_w(&w1.poly[i], &mut wp[i*W1_BYTES..]);
        }

        /* CRH(mu+w1) */
        let mut seed = [0u8; CRHBYTES + K*W1_BYTES];
        seed[..CRHBYTES].copy_from_slice(&mu[..]);
        seed[CRHBYTES..].copy_from_slice(&wp[..]);
        crh(&seed, &mut c_hat);

        let mut c = SampleInBall(&c_hat);
        c = ntt(&c);
          
        /* z = y + c s1*/
        z = v_add(&y, &inv_Ntt(&p_mult_v(&c, &s1)));
        /* tmp = w - c s2 */
        let tmp = &inv_Ntt(&p_mult_v(&c, &s2));
        let tmp = &Neg(tmp);
        let tmp = v_add(&w, &tmp);
        /* r0 = lowbits(tmp) */
        let mut r0 = VecPoly::<K>::default();
        for i in 0..K
        {
            r0.poly[i] = lowBits(&tmp.poly[i]);
        }

        /* check norm (check if it reveals secret) */ 
        if inf_norm(&z) >= (GAMMA1 as i32 - BETA as i32) || inf_norm(&r0) >= (GAMMA2 as i32 - BETA as i32) 
        { continue;}

        /* make hint */
        let ct0 = inv_Ntt(&p_mult_v(&c, &t0));
        let (h, cnt) = make_hint(&Neg(&ct0), &v_add(&tmp, &ct0));
       
        
        /* check norm */
        if inf_norm(&ct0) >= GAMMA2 as i32 || cnt > OMEGA as i32 { continue;}  

         /* assert the hint works */
        /* use(h, tmp+ct0) = high(tmp) */
        let mut r1 = VecPoly::<K>::default();
        for i in 0..K
        {
            r1.poly[i] = highBits(&tmp.poly[i]);
        }
        
        for i in 0..K 
        {
            let res = useHint(&h.poly[i], &v_add(&tmp, &ct0).poly[i]);
            assert_eq!(res.coeff, r1.poly[i].coeff);
        } 

        /* pack signature = c_hat, z, h */ 
        pack_sign(&mut sig, &c_hat, &z, &h);  
        return sig;
    }
}


pub fn verify(pk: &[u8], m: &[u8], sig: &[u8]) -> bool
{
    /* initialisation */
    let mut z  : VecPoly<L> = VecPoly::<L>::default();
    let mut ch = [0u8; SEEDBYTES];
    let mut h  : VecPoly<K> = VecPoly::<K>::default();

    let mut rho = [0u8; SEEDBYTES];
    let mut t1 = VecPoly::<K>::default();

    let mut mu = [0u8; CRHBYTES];

    /* unboxing */
    unpack_pk(&pk, &mut rho, &mut t1);
    unpack_sign(&sig, &mut ch, &mut z, &mut h); 

    let A = expandA(&rho);
    
    /* generate mu */
    crh(&pk, &mut mu);
    let mut seed = [0u8; CRHBYTES+MESSAGE];
    seed[..CRHBYTES].clone_from_slice(&mu[..]);
    seed[CRHBYTES..CRHBYTES+m.len()].copy_from_slice(&m[..m.len()]);
    crh(&seed[..CRHBYTES + m.len()], &mut mu);

    let mut c = SampleInBall(&ch);
    c = ntt(&c);

    t1 = sMult(1<<D, &t1);
    t1 = Ntt(&t1);

    let ct1 = p_mult_v(&c, &t1);
    let az = m_mult_v(&A, &Ntt(&z));

    let tmp = v_sub(&az, &ct1);

    let mut w1 = VecPoly::<K>::default();
    for i in 0..K
    {
        w1.poly[i] = useHint(&h.poly[i], &tmp.poly[i]);
    }

    /* check */
    assert!(inf_norm(&z) < (GAMMA1 as i32 - BETA as i32));
    /* pack w1 */
    let mut wp = [0u8; W1_BYTES*K];
    for i in 0..K
    {
        pack_w(&w1.poly[i], &mut wp[i*W1_BYTES..]);
    }

    /* CRH(mu+w1) */
    let mut c_hat = [0u8; SEEDBYTES];
    let mut seed = [0u8; CRHBYTES + K*W1_BYTES];
    seed[..CRHBYTES].copy_from_slice(&mu[..]);
    seed[CRHBYTES..].copy_from_slice(&wp[..]);
    crh(&seed, &mut c_hat);
    
    assert_eq!(c_hat, ch);
    for i in 0..SEEDBYTES 
    {
        if c_hat[i] != ch[i] {return false;}
    }
    /* number of ones */
    

    return true;



}