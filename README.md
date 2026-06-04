# 6N Semiprime Neighbour: a Slower omega-Enrichment and Its Mechanism (Part XVI)

Relax the twin's right wing from prime to **prime-or-semiprime**: count centres
with 6N-1 prime and 6N+1 either prime or a product of exactly two primes
(Omega(6N+1) in {1,2}). This neighbour count enriches with omega_{>3}(N), but
**more slowly than the twin** — and the mechanism is that a factor-rich N's
small-prime avoidance helps a prime more than a semiprime.

**A slower enrichment (S9).** Normalised to omega=1, on S9 (3.5M centres with
6N-1 prime at omega=1):

| omega | #(6N-1 pr) | twin | semi-nbr | semiprime-only |
|------:|-----------:|-----:|---------:|---------------:|
| 1 | 3,549,831 | 1.000 | 1.000 | 1.000 |
| 2 | 8,643,794 | 1.125 | 1.082 | 1.066 |
| 3 | 7,372,158 | 1.297 | 1.191 | 1.152 |
| 4 | 2,613,621 | 1.522 | 1.327 | 1.254 |
| 5 |   351,717 | 1.814 | 1.492 | 1.372 |

The twin (right wing prime) rises fastest; the exactly-semiprime part (right wing
semiprime) rises slowest; the relaxed neighbour lies between. The enrichment ratio
(semiprime-only / twin) falls monotonically 1.00 -> 0.76 over omega=1..5, so the
semiprime follows a genuinely flatter omega-law — not a scaling of the twin.

**The mechanism.** The twin enrichment (Part XII) is a small-prime avoidance
effect: a factor-rich N forces, through dead(q) = {+-6^-1 mod q}, that 6N+-1 avoid
divisibility by the small primes q|N. For the right wing to be *prime* this is
all-or-nothing and maximally helpful. For it to be *semiprime*, it is less helpful:
a semiprime may carry one small factor, so suppressing small factors removes some
of the very semiprimes being counted. The enrichment is therefore weaker.

**Diagnostic confirmation (S9, same shell).** As omega grows, surviving semiprimes
6N+1 should shift toward two *large* factors. Confirmed: among semiprime 6N+1
(with 6N-1 prime), the fraction with a prime factor <=97 falls monotonically.

| omega | semi/o1 | fraction with factor <=97 |
|------:|--------:|--------------------------:|
| 1 | 1.000 | 0.5599 |
| 2 | 1.066 | 0.5368 |
| 3 | 1.152 | 0.5056 |
| 4 | 1.254 | 0.4668 |
| 5 | 1.372 | 0.4136 |

Factor-rich N pushes its semiprime neighbour toward the large x large regime,
exactly as the mechanism requires.

> **Attribution.** The single-prime and semiprime local densities are classical
> (Hardy-Littlewood for the prime wing; Landau's semiprime density for the relaxed
> wing). This paper introduces NO new density formula. The contribution is the
> measured fact that the semiprime neighbour enriches more slowly than the twin,
> and the factor-structure mechanism for it, confirmed by the falling small-factor
> fraction.
>
> **Open: the exact double-factor closed form.** The twin enrichment reduces to
> the clean single-prime law 1/(q-2) (Part XII). The semiprime condition
> Omega(6N+1)=2 couples TWO prime factors, and the omega-dependence does not
> separate into a single per-prime factor; a faithful model would be a
> double-factor convolution that does NOT reduce to a clean product at the percent
> level. We report the enrichment and mechanism empirically and leave the exact
> closed form open, rather than force an approximate product.
>
> **Scope.** No statement about the infinitude of twins, semiprime neighbours, or
> any pattern. This is a measured, factor-resolved account of a conditional rate.

Part I: doi:10.5281/zenodo.20470367 . Part XII: doi:10.5281/zenodo.20528446

---

## Layout

```
.
├── README.md
├── LICENSE                 (MIT)
├── CITATION.cff
├── data/
│   ├── semiprime_enrich.csv    omega, tot, twin/semi/semionly rates + normalised (S9)
│   └── semiprime_spf.csv       omega, semirate, semi_norm, frac_spf_le97 (S9)
├── code/
│   ├── semiprime_enrich.py     enrichment: by omega(N), rate of {6N-1 prime,
│   │                           6N+1 prime / semiprime}; emits CSV. Default S9.
│   ├── semiprime_diagnostic.py small-factor diagnostic: among semiprime 6N+1,
│   │                           fraction with smallest prime factor <=97, by omega;
│   │                           emits CSV. Default S9, same shell as enrichment.
│   └── make_semi_fig.py        2-panel figure, reads ../data
├── figures/                fig_paper16_semiprime.{pdf,png}
└── paper/                  Chen_6N_Paper16.{tex,pdf} + figure
```

## Reproducing

Requirements: Python 3.8+, `numpy`, `matplotlib`.

```bash
pip install numpy matplotlib

# 1. Enrichment curves (twin / semi-neighbour / semiprime-only by omega). Default S9.
#    Emits semiprime_enrich_S9.csv. Memory-light segmented sieve.
python code/semiprime_enrich.py
MAXK=8 python code/semiprime_enrich.py     # S8, faster

# 2. Small-factor diagnostic (same shell). Emits semiprime_spf_S9.csv.
python code/semiprime_diagnostic.py

# 3. Figure (reads ../data).
cd code && python make_semi_fig.py
```

### Conventions (same as Parts I-XV)

- Twin position (6N-1, 6N+1); relax right wing to prime-or-semiprime.
- Omega(v) = prime-factor count with multiplicity; semiprime iff Omega=2.
- omega_{>3}(N) = number of distinct prime factors >3 of the left centre N.
- dead(q) = {+-6^-1 mod q} : residues of N making a wing divisible by q.
- Rates normalised to omega=1 for shape comparison.
- Shell S_k: 6N in [10^(k-1), 10^k). Both enrichment and diagnostic on S9.
- Engine: segmented-sieve factorisation + Omega/spf interval sieve; S10 twin count
  23,988,173 matches Part I.

## License

MIT — see `LICENSE`.
