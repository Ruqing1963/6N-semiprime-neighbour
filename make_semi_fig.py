#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build the 2-panel semiprime figure from ../data (both panels on S9, same shell).
Left: enrichment curves — twin (fastest), semiprime-neighbour (intermediate),
exactly-semiprime part (slowest). Right: the small-factor diagnostic — fraction of
semiprime 6N+1 with a prime factor <=97 falls with omega (factor-rich N pushes
6N+1 toward large*large).
"""
import csv, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
en=list(csv.DictReader(open('../data/semiprime_enrich.csv')))
sp=list(csv.DictReader(open('../data/semiprime_spf.csv')))
om=np.array([int(r['omega']) for r in en])
tw=np.array([float(r['twin_norm']) for r in en])
se=np.array([float(r['semi_norm']) for r in en])
so=np.array([float(r['semionly_norm']) for r in en])
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(14,5.4))
ax1.plot(om,tw,'o-',color='#c0392b',lw=2.2,ms=9,label='twin (6N+1 prime)',zorder=4)
ax1.plot(om,se,'s-',color='#185FA5',lw=2,ms=8,label='semiprime-neighbour (6N+1 prime or semiprime)',zorder=3)
ax1.plot(om,so,'^--',color='#2ca25f',lw=1.8,ms=8,label='exactly-semiprime part (6N+1 semiprime)',zorder=3)
ax1.axhline(1,color='gray',ls=':',lw=1)
ax1.set_xlabel(r'$\omega_{>3}(N)$',fontsize=11)
ax1.set_ylabel(r'rate, normalised to $\omega=1$',fontsize=11)
ax1.set_title('Semiprime-neighbour enrichment is slower than the twin (S9)',fontsize=12)
ax1.legend(fontsize=8.5,loc='upper left'); ax1.grid(alpha=.25); ax1.set_xticks(om)
omS=np.array([int(r['omega']) for r in sp])
fr=np.array([float(r['frac_spf_le97']) for r in sp])
ax2.plot(omS,fr,'D-',color='#8e44ad',lw=2,ms=9)
ax2.set_xlabel(r'$\omega_{>3}(N)$',fontsize=11)
ax2.set_ylabel(r'fraction of semiprimes $6N{+}1$ with a factor $\leq97$',fontsize=11)
ax2.set_title('Mechanism: factor-rich $N$ pushes $6N{+}1$ toward\nlarge$\\times$large semiprimes (S9)',fontsize=12)
ax2.set_ylim(0.38,0.60); ax2.grid(alpha=.25); ax2.set_xticks(omS)
ax2.annotate('small-factor fraction\nfalls with $\\omega$',xy=(3.5,0.47),xytext=(1.6,0.42),
             fontsize=9,color='#8e44ad',arrowprops=dict(arrowstyle='->',color='#8e44ad'))
plt.suptitle('Semiprime-neighbour conditional density on the $6N$ skeleton: enrichment and its mechanism',fontsize=12.5,y=1.02)
plt.tight_layout()
plt.savefig('fig_paper16_semiprime.pdf',bbox_inches='tight')
plt.savefig('fig_paper16_semiprime.png',dpi=160,bbox_inches='tight')
print("figure saved")
