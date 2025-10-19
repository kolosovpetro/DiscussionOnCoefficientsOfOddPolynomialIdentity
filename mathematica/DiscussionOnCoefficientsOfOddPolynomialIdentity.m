(* ::Package:: *)

BeginPackage["DiscussionOnCoefficientsOfOddPolynomialIdentity`"]

A::usage= "A[n, k] returns the real coefficient A of non-negative integers n, k such that n <= k."
L::usage= "L[m, n, k] returns the polynomial L of integers m, n, k."
P::usage= "P[m, x, b] returns the polynomial P of m, x and b."
BivariateSum::usage="Returns bivariate sum F(n,r) = sum_{k=1}^{n} k^r (n-k)^r"
BivariateFaulhabersFormula::usage="Returns bivariate Faulhabers formula."

Begin["`Private`"]

Unprotect[Power];
Power[0|0., 0|0.] = 1;
Protect[Power];

A[n_, k_] := 0;
A[n_, k_] := (2k + 1) * Binomial[2k, k] * Sum[A[n, j] * Binomial[j, 2k + 1] * (-1)^(j - 1) / (j - k) * BernoulliB[2j - 2k], {j, 2k + 1, n}] /; 0 <= k < n;
A[n_, k_] := (2n + 1) * Binomial[2n, n] /; k == n;

L[m_, n_, k_] := Sum[A[m, r] * k^r * (n - k)^r, {r, 0, m}];
P[m_, n_, b_] := Sum[L[m, n, k], {k, 0, b - 1}];
BivariateSum[n_, r_] := Sum[k^r * (n-k)^r, {k, 1, n}];
BivariateFaulhabersFormula[n_, r_]:= 1/((2r+1) * Binomial[2r,r]) * n^(2r+1) + Sum[((-1)^r)/(r-k) * Binomial[r, 2k+1] * BernoulliB[2r-2k] * n^(2k+1), {k, 0, r-1}];

End[ ]
EndPackage[ ]



