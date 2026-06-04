#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Semiprime small-factor diagnostic, S9 (same shell as the enrichment curve).
Among centres with 6N-1 prime and 6N+1 semiprime (Omega(6N+1)=2), measure the
fraction whose smaller prime factor is small (<=97), binned by omega_{>3}(N).
The mechanism for the slower semiprime enrichment predicts this fraction FALLS
with omega: a factor-rich N's small-prime avoidance removes small*large
semiprimes, pushing survivors toward large*large. Memory-light segmented sieve.
Default S9. Requires: numpy.
"""
import numpy as np, math, os, time
def primes_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(math.isqrt(n))+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0].astype(np.int64)
MAXK=int(os.environ.get("MAXK",9))
LO=10**(MAXK-1)//6+1; HI=10**MAXK//6; SEG=2_000_000
vmax=6*HI+2; PB=int(math.isqrt(vmax))+1; BP=primes_upto(PB)
SMALL=97
OMAX=7
def Omega_spf(lo,hi):
    # Omega(v) with multiplicity, and smallest prime factor, for v in [lo,hi)
    sz=hi-lo; rem=np.arange(lo,hi,dtype=np.int64)
    Om=np.zeros(sz,np.int16); spf=np.zeros(sz,np.int64)
    for p in BP:
        if p*p>hi-1: break
        f=((lo+p-1)//p)*p
        if f>=hi: continue
        idx=np.arange(f-lo,sz,p)
        if idx.size==0: continue
        first=idx[spf[idx]==0]
        if first.size: spf[first]=p
        sub=rem[idx]
        while True:
            m=(sub%p)==0
            if not m.any(): break
            sub[m]//=p; Om[idx[m]]+=1
        rem[idx]=sub
    big=rem>1
    Om[big]+=1
    needs=big&(spf==0)
    spf[needs]=rem[needs]   # prime (or large prime factor) is its own spf
    return Om,spf
def omega3(lo,hi):
    sz=hi-lo; rem=np.arange(lo,hi,dtype=np.int64); ob=np.zeros(sz,np.int16)
    for p in BP:
        if p>HI or p*p>hi-1: break
        f=((lo+p-1)//p)*p
        if f>=hi: continue
        idx=np.arange(f-lo,sz,p)
        if idx.size==0: continue
        sub=rem[idx]; m=(sub%p)==0
        while m.any(): sub[m]//=p; m=(sub%p)==0
        rem[idx]=sub
        if p>3: ob[idx]+=1
    ob[rem>1]+=1
    return ob
tot=np.zeros(OMAX+2); semic=np.zeros(OMAX+2); small=np.zeros(OMAX+2)
t0=time.time(); n=LO
while n<=HI:
    nh=min(n+SEG,HI+1); Narr=np.arange(n,nh,dtype=np.int64)
    omc=np.clip(omega3(n,nh),0,OMAX+1)
    vlo=6*n-1; vhi=6*(nh-1)+1; Om,spf=Omega_spf(vlo,vhi+1)
    Oa=Om[(6*Narr-1)-vlo]; Ob=Om[(6*Narr+1)-vlo]; spfb=spf[(6*Narr+1)-vlo]
    leftp=(Oa==1); semi=leftp&(Ob==2)
    for om in range(1,OMAX+1):
        s=leftp&(omc==om); tot[om]+=s.sum()
        ss=semi&(omc==om); semic[om]+=ss.sum()
        small[om]+=(ss&(spfb<=SMALL)).sum()
    n=nh
print(f"S{MAXK}: semiprime small-factor diagnostic ({time.time()-t0:.0f}s)")
oms=[om for om in range(1,OMAX+1) if tot[om]>=20000]
sr1=semic[oms[0]]/tot[oms[0]]
print(f"\n{'omega':>6}{'semirate':>10}{'semi/o1':>9}{'frac spf<=97':>14}")
for om in oms:
    sr=semic[om]/tot[om]; fr=small[om]/semic[om] if semic[om]>0 else 0
    print(f"{om:>6}{sr:>10.5f}{sr/sr1:>9.3f}{fr:>14.4f}")
print("\nIf frac spf<=97 falls monotonically with omega, the mechanism is confirmed")
print("on S9 (same shell as the enrichment): factor-rich N pushes semiprime 6N+1")
print("toward large*large. Compare to S8: 0.6025,0.5761,0.5423,0.4968 (omega 1-4).")

# ---- emit CSV (semiprime_spf_S{K}.csv: per-omega small-factor fraction) ----
import csv as _csv
with open(f'semiprime_spf_S{MAXK}.csv','w',newline='') as _f:
    _w=_csv.writer(_f); _w.writerow(['omega','semirate','semi_norm','frac_spf_le97'])
    for om in oms:
        sr=semic[om]/tot[om]; fr=small[om]/semic[om] if semic[om]>0 else 0
        _w.writerow([om,f'{sr:.5f}',f'{sr/sr1:.3f}',f'{fr:.4f}'])
print(f"\n[ok] wrote semiprime_spf_S{MAXK}.csv")
