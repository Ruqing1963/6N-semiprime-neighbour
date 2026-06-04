#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Semiprime-neighbour conditional density, S8/S9, memory-light segmented sieve.
Twin position (6N-1, 6N+1). Stratify by omega_{>3}(N) and compare two shapes:
  twin:      6N-1 prime AND 6N+1 prime
  semi-nbr:  6N-1 prime AND 6N+1 prime-or-semiprime (Omega(6N+1) in {1,2})
where Omega counts prime factors with multiplicity. Question (recon): is the
semiprime omega-enrichment clean and distinct from the twin (S7 showed it rises
more slowly)? This confirms the shape on larger shells and extends to higher omega.
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
vmax=6*HI+2
PB=int(math.isqrt(vmax))+1; BP=primes_upto(PB)
OMAX=8
from collections import defaultdict
tot=np.zeros(OMAX+2,dtype=np.int64)
twin=np.zeros(OMAX+2,dtype=np.int64)
semi=np.zeros(OMAX+2,dtype=np.int64)
t0=time.time()
def bigOmega_interval(lo,hi):
    # Omega(v) with multiplicity for v in [lo,hi)
    sz=hi-lo; rem=np.arange(lo,hi,dtype=np.int64); Om=np.zeros(sz,np.int16)
    for p in BP:
        if p*p>hi-1: break
        f=((lo+p-1)//p)*p
        if f>=hi: continue
        idx=np.arange(f-lo,sz,p)
        if idx.size==0: continue
        sub=rem[idx]
        while True:
            m=(sub%p)==0
            if not m.any(): break
            sub[m]//=p; Om[idx[m]]+=1
        rem[idx]=sub
    Om[rem>1]+=1  # residual large prime factor
    return Om
def omega3_interval(lo,hi):
    # omega_{>3}(N) distinct prime factors >3, for N in [lo,hi)
    sz=hi-lo; rem=np.arange(lo,hi,dtype=np.int64); ob=np.zeros(sz,np.int16)
    for p in BP:
        if p>HI: break
        if p*p>hi-1: break
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
n=LO
while n<=HI:
    nh=min(n+SEG,HI+1)
    Narr=np.arange(n,nh,dtype=np.int64)
    omN=np.clip(omega3_interval(n,nh),0,OMAX+1)
    # compute Omega for the exact values 6N-1 and 6N+1
    a_vals=6*Narr-1; b_vals=6*Narr+1
    # interval covering both a and b
    vlo=6*n-1; vhi=6*(nh-1)+1
    Omv=bigOmega_interval(vlo,vhi+1)
    Oa=Omv[a_vals-vlo]; Ob=Omv[b_vals-vlo]
    left_prime=(Oa==1)
    for om in range(1,OMAX+1):
        sel=left_prime&(omN==om); ns=sel.sum()
        if ns==0: continue
        tot[om]+=ns
        ob_sel=Ob[sel]
        twin[om]+=int((ob_sel==1).sum())
        semi[om]+=int(((ob_sel==1)|(ob_sel==2)).sum())
    n=nh
print(f"S{MAXK}: done ({time.time()-t0:.0f}s)")
oms=[om for om in range(1,OMAX+1) if tot[om]>=20000]
tw1=twin[oms[0]]/tot[oms[0]]; se1=semi[oms[0]]/tot[oms[0]]
print(f"\n{'omega':>6}{'(6N-1 pr)':>11}{'twinrate':>10}{'semirate':>10}{'tw/o1':>8}{'semi/o1':>9}")
for om in oms:
    tr=twin[om]/tot[om]; sr=semi[om]/tot[om]
    print(f"{om:>6}{int(tot[om]):>11}{tr:>10.5f}{sr:>10.5f}{tr/tw1:>8.3f}{sr/se1:>9.3f}")
print("\ntw/o1 = twin enrichment; semi/o1 = semiprime-relaxed enrichment.")
print("S7 showed semi rises more slowly than twin; confirm here on a larger shell")
print("and across more omega. Smooth, monotone, distinct shape => worth modelling.")

# ---- emit CSV (semiprime_enrich_S{K}.csv: per-omega twin/semi/semionly enrichment) ----
import csv as _csv
with open(f'semiprime_enrich_S{MAXK}.csv','w',newline='') as _f:
    _w=_csv.writer(_f); _w.writerow(['omega','tot','twinrate','semirate','semionly','twin_norm','semi_norm','semionly_norm'])
    t1=twin[oms[0]]/tot[oms[0]]; se1=semi[oms[0]]/tot[oms[0]]
    so1=(semi[oms[0]]-twin[oms[0]])/tot[oms[0]]
    for om in oms:
        tr=twin[om]/tot[om]; sr=semi[om]/tot[om]; so=(semi[om]-twin[om])/tot[om]
        _w.writerow([om,int(tot[om]),f'{tr:.5f}',f'{sr:.5f}',f'{so:.5f}',f'{tr/t1:.3f}',f'{sr/se1:.3f}',f'{so/so1:.3f}'])
print(f"\n[ok] wrote semiprime_enrich_S{MAXK}.csv")
