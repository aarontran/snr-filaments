X-ray line review
=================

Aaron Tran (2014 Nov 6)

Quick review/primer on atomic lines, because this trips me up whenever I come
back to emission lines.  Mainly for my own reference.

Note: AtomDB numbering of energy states refers to each distinct energy level,
but the numbering doesn't correspond to n, l, m, or anything.  So don't get
tripped up... (I know I came across this before)

### Nomenclature

* K,L,M,... shells/levels indicate n = 1,2,3
* alpha, beta, gamma indicates n=2,3,4,... -> 1 transitions
* Ion notation: I indicates neutral

Ly(alpha/beta/...) == K(alpha/beta/...) transition in  H-like atom
He(alpha/beta/...) == K(alpha/beta/...) transition in He-like atom

Roughly, energies scale as (Z-1)^2 for H/He-like atoms (Moseley's law)
E.g., Si XIII He alpha is ~1.85 keV, S XV He alpha is ~2.45 keV
Si's He-like core has Z=13, S has 15, then (14/12)^2 * 1.85 = 2.52 is not a bad
estimate -- good to ~3%.

Looking at Hayato et al. (2010) -- it appears that naming lines by K suggests
that we can't differentiate between He/Ly alpha, so better to just call it a
single K alpha line.  Same idea for L-shell.

### SNR elements (near fully-ionized)

Relevant elements' ionization states

    Element     He,H-like states
    ----------  ----------------
    O   (Z=8)   VII, VIII
    Ne  (Z=10)  IX, X
    Mg  (Z=12)  XI, XII
    Si  (Z=14)  XIII, XIV
    S   (Z=16)  XV, XVI
    Ar  (Z=18)  XVII, XVIII
    Ca  (Z=20)  XIX, XX
    Ti  (Z=22)
    Cr  (Z=24)
    Fe  (Z=26)  XXV, XXVI
    Ni  (Z=28)

alpha-group elements: O,Ne,Mg,Si,S,Ar,Ca
Fe-group elements: Z=21-28 (Sc,Ti,V,Cr,Mn,Fe,Co,Ni)

Fe-group decay chains (Vink 2012):
    Ti-44 -> Sc-44 -> Ca-44 (time 86 yr, 6 hr)
    Ni-57 -> Co-57 -> Fe-57 (time 52 hr, 1 yr)
    Ni-56 -> Co-56 -> Fe-56 (time  8 d, 110 d)

### SNR relevant transitions

In general we consider just He and H-like ions, especially above ~1.5 keV.
Around and below 1.5 keV, we getting swamped by Fe L-shell emission and a huge
melange of other lines (O, Ne, Mg, etc).
Reference: http://www.nist.gov/pml/data/asd.cfm

For H,He-like ions, the most important emission lines are:
* He alpha: 1s2p(S=0), 1s2s(S=1)
* Ly alpha: 2p(J=1/2,3/2) (fine structure doublet)
* He beta: 1s3p(S=0) -> 1s2
The He alpha lines are smeared across ~0.03 keV (near 1-2 keV), 0.05 keV (near
4-6 keV) -- due to the splitting between 1s2s -> 1s2 and 1s2p -> 1s2.

He alpha: 1s2s(S=0,L=0,J=0) and 1s2p(S=1,L=1,J=0/1/2) are forbidden/disfavored
Ordering of energy levels looks something like:
(I'm not sure I have the J orderings correct)

    1s^2(S=0,L=0,J=0) ground

    1s2s(S=1,L=0,J=1) 2nd strongest (triplet from S=1, degenerate)
    1s2s(S=0,L=0,J=0) forbidden

    1s2p(S=1,L=1,J=0) disfavored
    1s2p(S=1,L=1,J=1) disfavored
    1s2p(S=1,L=1,J=2) disfavored
    1s2p(S=0,L=1,J=1) strongest

Ly alpha, 2s(J=1/2) is forbidden.  Ordering looks like:

    1s(S=1/2,L=0,J=1/2) ground

    2s(S=1/2,L=0,J=1/2) forbidden
    2p(S=1/2,L=1,J=1/2) strong (J=1/2,3/2 form fine structure doublet)
    2p(S=1/2,L=1,J=3/2) strong

Then He beta, Ly beta look similar of course -- same angular momentum rules.

E.g., the whole set for Si is:
* 1.85 keV, Si XIII He alpha
* 2.00 keV, Si XIV  Ly alpha
* 2.18 keV, Si XIII He beta
* 2.38 keV, Si XIV  Ly beta

For Sulfur:
* 2.45 keV, S XV He alpha
* 2.62 keV, S XVI Ly alpha
* 2.89 keV, S XV He beta
* 3.11 keV, S XVI Ly beta

### Actually useful stuff -- list of transitions

Fe L-shell transitions arises from Fe XXIV to Fe XVII (16+) (Li to Ne-like)
Strongest lines are: 2p5,3d1 -> 2p6 at 0.83 keV, and 2p5,3s1 -> 2p6 at 0.73 keV
(these are Fe XVII lines).  For all lines below ~1.6 keV there's a risk of
getting swamped out by Fe emission.

Relative strengths depend on ionization age, abundances, etc.
E.g., SNR forward shocks should not have much Fe.  So no worries here...

* 0.425 keV = N VI, He alpha
* 0.436 keV = C VI, Ly beta
* 0.500 keV = N VII, Ly alpha

* 0.57 keV = O VII, He alpha
* 0.654 keV = O VIII, Ly alpha
* 0.67 keV = O VII, He beta
* 0.775 keV = O VIII, Ly beta

* 0.91 keV = Ne IX, He alpha
* 1.02 keV = Ne X, Ly alpha
* 1.07 keV = Ne IX, He beta

* 1.34 keV = Mg XI, He alpha
* 1.47 keV = Mg XII, Ly alpha
* 1.58 keV = Mg XI, He beta
* 1.75 keV = Mg XII, Ly beta

* 1.85 keV = Si XIII, He alpha
* 2.00 keV = Si XIV, Ly alpha
* 2.183 keV = Si XIII, He beta (He gamma, delta faint at ~2.3 keV)
* 2.377 keV = Si XIV, Ly beta

* 2.45 keV = S XV, He alpha
* 2.62 keV = S XVI, Ly alpha
* 2.89 keV = S XV, He beta
* 3.11 keV = S XVI, Ly beta

* 3.12 keV = Ar XVII, He alpha
* 3.32 keV = Ar XVIII, Ly alpha
* (DARK MATTER?!?!?!?!?!?!)
* 3.69 keV = Ar XVII, He beta

* 3.88 keV = Ca XIX, He alpha
* 4.10 keV = Ca XX, Ly alpha

(not much around here)

* 6.68 keV = Fe XXV, He alpha (smeared across 6.64-6.70)
* 6.96 keV = Fe XXVI, Ly alpha

In Tycho, we care about:
* 1.34 keV (Mg He alpha, maybe)
* 1.85 keV (Si He alpha)
* 2.18 keV (Si He beta)
* 2.45 keV (S He alpha)
* 2.89 keV (S He beta)
* 3.12 keV (Ar He alpha)
* 3.88 keV (Ca He alpha)
* 6.8 keV (Fe He/Ly alpha = K alpha)



