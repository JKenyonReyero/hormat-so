DATAP_TEMPLATE = """{a}   {z}    {elab}    {qval}
0.1  {nstep}   {lmax}
60
49.0303    1.1998   0.69
7.4458    1.2363   0.69
1.9857     1.2363   0.69
0     1.25   0.65
49.0303    1.1998   0.69
7.4458    1.2363   0.69
1.9857     1.2363   0.69
0     1.25   0.65
1.26677728


This is CH89 from twofnr with no SO.
Notes: if qvalue.ne.0 then Elab is that of the
deuteron in (d,p) and proton step length is
also adjusted for use in twofnr. If qvalue=0
then Elab is proton lab energy and the input 
step length is used - e.g. for use in (p,d).

The input target mass and charge are those in
the proton-target system, whatever qvalue used.

Potential 1 is that from which iterations start.
Potential 2 must be input but is not used if 
non-local potential choice is made.
----------------------------------------
target mass   charge    lab energy   qvalue
step          points    max L-value
number of iterations - 60 usually OK

--- starting potential iterations
v1r    r1     a1      real volume  1
w1i    r1i    a1i     imag surface 1
wv1i   rv1i   av1i    imag volume  1
v1so   r1so   a1so    spin orbit   1

--- local potential if local option used
v2r    r2     a2      real volume  2
w2i    r2i    a2i     imag surface 2
wv2i   rv2i   av2i    imag volume  2
v2so   r2so   a2so    spin orbit   2

rc (Coulomb radius parameter)
----------------------------------------
"""