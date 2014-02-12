# Chebfun examples inventory

Please help with this! For each of the examples, verify
that the code runs and the output looks reasonable.
This list is sorted by category.

### approx

- *AbsoluteValue.m* - good
- ***AbsoluteValueScaled.m*** - stalls at line 60:

  ```
  hold on, semilogy(abs(rs-abs(x)),'r',LW,lw)
  ```
  The first `abs` is the problem.
- *BernsteinPolys.m* - good
- *BestApprox.m* - good
- *BestL2Approximation.m* - good
- *BSplineConv.m* - good
- ***BumpFunction.m*** - this one is going seriously wrong,
  I can't even describe. It's not adding the functions
  together properly, evidently, because the final result
  is not random and it has rank 100.
- ***ChebCoeffs.m*** - strange error in `chebpolyplot`, which
  appears to fail only when given arguments for linewidth, etc
  (line 44):

  ```
  >> chebpolyplot(exp(x),'.-',LW,1,MS,16)
  Too many outputs requested.  Most likely cause
  is missing [] around left hand side that has a
  comma separated list expansion.
  Error in chebfun/chebpolyplot (line 39)
      col = col{:};
  ```
- *ChebfunFFT.m* - good
- *ChebPolysHigham.m* - good
- ***CommunicationSystem.m*** - after `feature-deltafun` is
  merged this one will probably work, but `chebtune` still
  gives an admittedly strange warning:

  ```
  Warning: File: chebtune.m Line: 68 Column: 10

  The expression on this line will generate an error when executed.
  The error will be: Error using vertcat:
  Dimensions of matrices being concatenated are not consistent.
  ```
- ***Convergence.m*** - this one works, but singfun spits out a
  lot of warnings related to `singfun.plus`.
- ***DivergentSeries.m*** - `unbndfun` not in development
- *Entire.m* - good
- *ExactChebCoeffs.m* - good
- *FastChebyshevLegendreTransform.m* - good
- ***GammaFun.m*** - errors on `minandmax` on line 49 or so:

  ```
  Undefined function 'issmooth' for input
  arguments of type 'chebtech2'.
  Error in singfun/diff (line 74)
  if ( issmooth(f) )
  Error in bndfun/diff (line 34)
      f.onefun = diff(f.onefun, k,
      dim)/rescaleFactork;
  Error in chebfun/diff>diffContinuousDim (line
  122)
          funs{k} = diff(funs{k});
  Error in chebfun/diff (line 62)
          F(k) = diffContinuousDim(F(k), n);
  Error in chebfun/minandmax>localMinAndMax (line
  138)
  df = diff(f);
  Error in chebfun/minandmax (line 37)
      [y, x] = localMinAndMax(f);
  ```
- *GreedyInterp.m* - good
- *Halphen.m* - good
- ***HermiteBasis.m*** - requires `linop`
- ***InteractiveInterp.m*** - not working, not sure why (no errors)
- *LebesgueConst.m* - good
- ***Local.m*** - requires some things with `chebop` and `chebpref`
  that I'm not sure how to work
- ***OrthPolys.m*** - complains that you can't use `subsref` to set
  a column of a quasimatrix equal to given chebfun. I don't know
  what this syntax is supposed to look like now.