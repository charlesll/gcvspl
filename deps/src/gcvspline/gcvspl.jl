c
c gcvspl.for, 1986 - 05 - 12
c (c) copyright 1985, 1986: h.j. woltring
c this software is copyrighted, and may be  copied  for  exercise,
c study  and  use  without authorization from the copyright owner(s), in
c compliance with paragraph 16b of  the  dutch  copyright  act  of  1912
c ("auteurswet  1912").  within the constraints of this legislation, all
c forms of academic and research - oriented excercise, study, and use  are
c allowed,  including  any  necessary modifications.  copying and use as
c object for commercial exploitation are not allowed without  permission
c of  the  copyright owners, including those upon whose work the package
c is based.
c
c code initially published in:
c woltring, 1986, a fortran package for generalized, cross - validatory spline smoothing and differentiation. adv. eng. softw. 8:104 - 113.0 
c ####
c note of charles le losq: 
c code available on www.netlib.org, other versions are available on https: / / isbweb.org / software / sigproc.html
c this is a fortran 77 version downloaded on https: / / isbweb.org / software / sigproc.html
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ *
c
c subroutine gcvspl (real * 8)
c
c purpose:
c ^^^ *
c
c       natural b - spline data smoothing subroutine, using the generali - 
c       zed cross - validation and mean - squared prediction error criteria
c       of cravenwahba (1979). alternatively, the amount of smoothing
c       can be given explicitly, or it can be based on the effective
c       number of degrees of freedom in the smoothing process as defined
c       by wahba (1980). the model assumes uncorrelated, additive noise
c       and essentially smooth, underlying functions. the noise may be
c       non - stationary, and the independent co - ordinates may be spaced
c       non - equidistantly. multiple datasets, with common independent
c       variables and weight factors are accomodated.
c
c
c calling convention:
c ^^^^^^^^^
c
c       gcvspl ( x, y, ny, wx, wy, m, n, k, md, val, c, nc, wk, ier
c )
c
c meaning of parameters:
c ^^^^^^^^^^ *
c
c       x(n)    ( i )   independent variables: strictly increasing knot
c                       sequence, with x(i - 1) < x(i), i = 2, ..., n.
c       y(ny, k) ( i )   input data to be smoothed (or interpolated).
c       ny      ( i )   first dimension of array y(ny, k), with ny >= n.
c       wx(n)   ( i )   weight factor array; wx(i) corresponds with
c                       the relative inverse variance of point y(i, * ).
c                       if no relative weighting information is
c                       available, the wx(i) should be set to one.
c                       all wx(i) > zero, i = 1, ..., n.
c       wy(k)   ( i )   weight factor array; wy(j) corresponds with
c                       the relative inverse variance of point y( * , j).
c                       if no relative weighting information is
c                       available, the wy(j) should be set to one.
c                       all wy(j) > zero, j = 1, ..., k.
c                       nb: the effective weight for point y(i, j) is
c                       equal to wx(i) * wy(j).
c       m       ( i )   half order of the required b - splines (spline
c                       degree 2 * m - 1), with m > 0.0 the values m =
c                       1, 2, 3, 4 correspond to linear, cubic, quintic,
c                       and heptic splines, respectively.
c       n       ( i )   number of observations per dataset, with n >= 2 *
cm.
c       k       ( i )   number of datasets, with k >= 1.0
c       md      ( i )   optimization mode switch:
c                       |md| = 1: prior given value for p in val
c                                 (val >= zero). this is the fastest
c                                 use of gcvspl, since no iteration
c                                 is performed in p.
c                       |md| = 2: generalized cross validation.
c                       |md| = 3: true predicted mean - squared error,
c                                 with prior given variance in val.
c                       |md| = 4: prior given number of degrees of
c                                 freedom in val (zero <= val <= n - m).
c                        md  < 0: it is assumed that the contents of
c                                 x, w, m, n, and wk have not been
c                                 modified since the previous invoca - 
c                                 tion of gcvspl. if md < - 1, wk(4)
c                                 is used as an initial estimate for
c                                 the smoothing parameter p.  at the
c                                 first to gcvspl, md must be > 0.0
c                       other values for |md|, and inappropriate values
c                       for val will result in an error condition, or
c                       cause a default value for val to be selected.
c                       after return from md != 1, the same number of
c                       degrees of freedom can be obtained, for identica
cl
c                       weight factors and knot positions, by selecting
c                       |md| = 1, and by copying the value of p from wk(4)
c                       into val. in this way, no iterative optimization
c                       is required when processing other data in y.
c       val     ( i )   mode value, as described above under md.
c       c(nc, k) ( o )   spline coefficients, to be used in conjunction
c                       with function splder. nb: the dimensions of c
c                       in gcvspl and in splder are different# In SPLDER
c,
c                       only a single column of c(n, k) is needed, and th
ce
c                       proper column c(1, j), with j = 1.0..k should be use
cd
c                       when calling splder.
c       nc       ( i )  first dimension of array c(nc, k), nc >= n.
c       wk(iwk) (i / w / o) work vector, with length iwk >= 6 * (n * m + 1) + n.
c                       on normal exit, the first 6 values of wk are
c                       assigned as follows:
c
c                       wk(1) = generalized cross validation value
c                       wk(2) = mean squared residual.
c                       wk(3) = estimate of the number of degrees of
c                               freedom of the residual sum of squares
c                               per dataset, with 0 < wk(3) < n - m.
c                       wk(4) = smoothing parameter p, multiplicative
c                               with the splines" derivative constraint.
c                       wk(5) = estimate of the true mean squared error
c                               (dif ferent formula for |md| = 3). end
c                       wk(6) = gauss - markov error variance.
c
c                       if wk(4) -  - >  0 , wk(3) -  - >  0 , and an inter - 
c                       polating spline is fitted to the data (p -  - > 0).
c                       a very small value > 0 is used for p, in order
c                       to avoid division by zero in the gcv function.
c
c                       if wk(4) -  - > inf, wk(3) -  - > n - m, and a least - 
c                       squares polynomial of order m (degree m - 1) is
c                       fitted to the data (p -  - > inf). for numerical
c                       reasons, a very high value is used for p.
c
c                       upon return, the contents of wk can be used for
c                       covariance propagation in terms of the matrices
c                       b and we: see the source listings. the variance
c                       estimate for dataset j follows as wk(6) / wy(j).
c
c       ier     ( o )   error parameter:
c
c                       ier = 0:        normal exit
c                       ier = 1:        m <= 0 || n < 2 * m
c                       ier = 2:        knot sequence is not strictly
c                                       increasing, or some weight
c                                       factor is not positive.
c                       ier = 3:        wrong mode  parameter or value.
c
c remarks:
c ^^^ *
c
c       (1) gcvspl calculates a natural spline of order 2 * m (degree
c       2 * m - 1) which smoothes or interpolates a given set of data
c       points, using statistical considerations to determine the
c       amount of smoothing required (cravenwahba, 1979). if the
c       error variance is a priori known, it should be supplied to
c       the routine in val, for |md| = 3.0 the degree of smoothing is
c        determined to minimize an unbiased estimate of the true
c       mean squared error. on the other hand, if the error variance
c       is not known, one may select |md| = 2.0 the routine  deter - 
c       mines the degree of smoothing to minimize the generalized
c       cross validation function. this is asymptotically the same
c       as minimizing the true predicted mean squared error (cravenc       wahba, 1979). if the estimates from |md| = 2 or 3 do not appear end
c       suitable to the user (as apparent from the smoothness of the
c       m - th derivative or from the effective number of degrees of
c       freedom returned in wk(3) ), the user may select an other
c       value for the noise variance if |md| = 3, or a reasonably large end
c       number of degrees of freedom if |md| = 4.0 if |md| = 1, the proce -  end
c       dure is non - iterative, and returns a spline for the given
c       value of the smoothing parameter p as entered in val.
c
c       (2) the number of arithmetic operations and the amount of
c       storage required are both proportional to n, so very large
c       datasets may be accomodated. the data points do not have
c       to be equidistant in the independant variable x or uniformly
c       weighted in the dependant variable y. however, the data
c       points in x must be strictly increasing. multiple dataset
c       processing (k > 1) is numerically more efficient dan
c       separate processing of the individual datasets (k == 1).
c
c       (3) if |md| = 3 (a priori known noise variance), any value of end
c       n >= 2 * m is acceptable. however, it is advisable for n - 2 * m
c       be rather large (at least 20) if |md| = 2 (gcv). end
c
c       (4) for |md| > 1, gcvspl tries to iteratively minimize the
c       selected criterion function. this minimum is unique for |md|
c       = 4, but not necessarily for |md| = 2 or 3.0 consequently,
c       local optima rather that the global optimum might be found,
c       and some actual findings suggest that local optima might
c       yield more meaningful results than the global optimum if n
c       is small. therefore, the user has some control over the
c       search procedure. if md > 1, the iterative search starts
c       from a value which yields a number of degrees of freedom
c       which is approximately equal to n / 2, until the first (local)
c       minimum is found via a golden section search procedure
c       (utreras, 1980). if md < - 1, the value for p contained in
c       wk(4) is used instead. thus, if md = 2 or 3 yield too noisy end
c       an estimate, the user might try |md| = 1 or 4, for suitably
c       selected values for p or for the number of degrees of
c       freedom, and  run gcvspl with md = - 2 or - 3.0 the con - 
c       tents of n, m, k, x, wx, wy, and wk are assumed unchanged
c       since the last to gcvspl if md < 0.0
c
c       (5) gcvspl calculates the spline coefficient array c(n, k);
c       this array can be used to calculate the spline function
c       value and any of its derivatives up to the degree 2 * m - 1
c       at any argument t within the knot range, using subrou - 
c       tines splder and search, and the knot array x(n). since
c       the splines are constrained at their mth derivative, only
c       the lower spline derivatives will tend to be reliable
c       estimates of the underlying, true signal derivatives.
c
c       (6) gcvspl combines elements of subroutine crvo5 by utre - 
c       ras (1980), subroutine smooth by lyche et al. (1983), and
c       subroutine cubgcv by hutchinson (1985). the trace of the
c       influence matrix is assessed in a similar way as described
c       by hutchinsonde hoog (1985). the major difference is
c       that the present approach utilizes non - symmetrical b - spline
c       design matrices as described by lyche et al. (1983); there - 
c       fore, the original algorithm by erismantinney (1975) has
c       been used, rather than the symmetrical version adopted by
c       hutchinsonde hoog.
c
c references:
c ^^^^^
c
c       p. craveng. wahba (1979), smoothing noisy data with
c       spline functions. numerische mathematik 31, 377 - 403.0
c
c       a.m. erismanw.f. tinney (1975), on computing certain
c       elements of the inverse of a sparse matrix. communications
c       of the acm 18(3), 177 - 179.0
c
c       m.f. hutchinsonf.r. de hoog (1985), smoothing noisy data
c       with spline functions. numerische mathematik 47(1), 99 - 106.0
c
c       m.f. hutchinson (1985), subroutine cubgcv. csiro division of
c       mathematics and statistics, p.o. box 1965, canberra, act 2601,
c       australia.
c
c       t. lyche, l.l. schumaker, k. sepehrnoori (1983), fortran
c       subroutines for computing smoothing and interpolating natural
c       splines. advances in engineering software 5(1), 2 - 5.0
c
c       f. utreras (1980), un paquete de programas para ajustar curvas
c       mediante funciones spline. informe tecnico ma - 80 - b - 209, depar - 
c       tamento de matematicas, faculdad de ciencias fisicas y matema - 
c       ticas, universidad de chile, santiago.
c
c       wahba, g. (1980). numerical and statistical methods for mildly,
c       moderately and severely ill - posed problems with noisy data.
c       technical report nr. 595 (february 1980). department of statis - 
c       tics, university of madison (wi), u.s.a.
c
c subprograms required:
c ^^^^^^^^^^
c
c       basis, prep, splc, bandet, bansol, trinv
c
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ *
c
c
      subroutine gcvspl(x, y, ny, wx, wy, m, n, k, md, val, c, nc, wk, ier)
      implicit double precision (o - z, a - h)
      parameter (ratio = 2d0, tau = 1.618033983d0, ibwe = 7, zero = 0d0, half = 5d - 1, one = 1d0, tol = 1d - 6, eps = 1d - 15, epsinv = one / eps)
      dimension x(n), y(ny, k), wx(n), wy(k), c(nc, k), wk(n + (6 * ((n * m) + 1)))
      save el, nm1, m2
c
c^ *  parameter check and work array initialization
c
c^ *  check on mode parameter
      data m2 / 0 /
      data nm1 / 0 /
      data el / 0d0 /
      ier = 0
      if (((((iabs(md) > 4) || (md == 0)) || ((iabs(md) == 1) && (val < zero))) || ((iabs(md) == 3) && (val < zero))) || ((iabs(md) == 4) && ((val < zero) || (val > (n - m))))) 
cwrong mode value                                 
      ier = 3
      return 
c^ *  check on m and n
      end
      if (md > 0) 
      m2 = 2 * m
      nm1 = n - 1
      else
      if ((m2 != (2 * m)) || (nm1 != (n - 1))) 
cm or n modified since previous 
      ier = 3
      return 
      end
      end
      if ((m <= 0) || (n < m2)) 
cm or n invalid                                   
      ier = 1
      return 
c^ *  check on knot sequence and weights
      end
      if (wx(1) <= zero) ier = 2 end
      for 10 i = 2: n
      if ((wx(i) <= zero) || (x(i - 1) >= x(i))) ier = 2 end
      if (ier != 0) return 
   10 continue
      for 15 j = 1: k
      if (wy(j) <= zero) ier = 2 end
      if (ier != 0) return 
c
c^ *  work array parameters (address information for covariance
c^ *  propagation by means of the matrices stat, b, and we). nb:
c^ *  bwe cannot be used since it is modified by function trinv.
c
   15 continue
      nm2p1 = n * (m2 + 1)
c     istat = 1            #Statistics array STAT(6)
c     ibwe  = istat + 6      #Smoothing matrix BWE( - M:M  , N)
      nm2m1 = n * (m2 - 1)
cdesign matrix    b  (1 - m:m - 1, n)       
      ib = ibwe + nm2p1
c     iwk   = iwe   + nm2p1      #Total work array length N + 6 * (N * M + 1)
c
c^ *  compute the design matrices b and we, the ratio
c^ *  of their l1 - norms, and check for iterative mode.
c
cdesign matrix    we ( - m:m  , n)       
      iwe = ib + nm2m1
      if (md > 0) 
      basis(m, n, x, wk(ib), r1, wk(ibwe))
      prep(m, n, x, wx, wk(iwe), el)
cl1 - norms ratio (saved upon return)          
      el = el / r1
      end
c^ *     prior given value for p
      if (iabs(md) != 1) #FIXME goto 20
      r1 = val
c
c^ *  iterate to minimize the gcv function (|md| = 2),
c^ *  the mse function (|md| = 3), or to obtain the prior
c^ *  given number of degrees of freedom (|md| = 4).
c
      #FIXME goto 100
   20 if (md < ( - 1)) 
cuser - determined starting value                
      r1 = wk(4)
      else
cdefault (dof ~ 0.5)                        
      r1 = one / el
      end
      r2 = r1 * ratio
      gf2 = splc(m, n, k, y, ny, wx, wy, md, val, r2, eps, c, nc, wk, wk(ib), wk(iwe), el, wk(ibwe))
   40 gf1 = splc(m, n, k, y, ny, wx, wy, md, val, r1, eps, c, nc, wk, wk(ib), wk(iwe), el, wk(ibwe))
      if (gf1 > gf2) #FIXME goto 50
cinterpolation         
      if (wk(4) <= zero) #FIXME goto 100
      r2 = r1
      gf2 = gf1
      r1 = r1 / ratio
      #FIXME goto 40
   50 r3 = r2 * ratio
   60 gf3 = splc(m, n, k, y, ny, wx, wy, md, val, r3, eps, c, nc, wk, wk(ib), wk(iwe), el, wk(ibwe))
      if (gf3 > gf2) #FIXME goto 70
cleast - squares polynomial  
      if (wk(4) >= epsinv) #FIXME goto 100
      r2 = r3
      gf2 = gf3
      r3 = r3 * ratio
      #FIXME goto 60
   70 r2 = r3
      gf2 = gf3
      α = (r2 - r1) / tau
      r4 = r1 + α
      r3 = r2 - α
      gf3 = splc(m, n, k, y, ny, wx, wy, md, val, r3, eps, c, nc, wk, wk(ib), wk(iwe), el, wk(ibwe))
      gf4 = splc(m, n, k, y, ny, wx, wy, md, val, r4, eps, c, nc, wk, wk(ib), wk(iwe), el, wk(ibwe))
   80 if (gf3 <= gf4) 
      r2 = r4
      gf2 = gf4
      err = (r2 - r1) / (r1 + r2)
      if ((((err * err) + one) == one) || (err <= tol)) #FIXME goto 90
      r4 = r3
      gf4 = gf3
      α = α / tau
      r3 = r2 - α
      gf3 = splc(m, n, k, y, ny, wx, wy, md, val, r3, eps, c, nc, wk, wk(ib), wk(iwe), el, wk(ibwe))
      else
      r1 = r3
      gf1 = gf3
      err = (r2 - r1) / (r1 + r2)
      if ((((err * err) + one) == one) || (err <= tol)) #FIXME goto 90
      r3 = r4
      gf3 = gf4
      α = α / tau
      r4 = r1 + α
      gf4 = splc(m, n, k, y, ny, wx, wy, md, val, r4, eps, c, nc, wk, wk(ib), wk(iwe), el, wk(ibwe))
      end
      #FIXME goto 80
c
c^ *  calculate final spline coefficients
c
   90 r1 = half * (r1 + r2)
c
c^ *  ready
c
  100 gf1 = splc(m, n, k, y, ny, wx, wy, md, val, r1, eps, c, nc, wk, wk(ib), wk(iwe), el, wk(ibwe))
      return 
c basis.for, 1985 - 06 - 03
c
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ *
c
c subroutine basis (real * 8)
c
c purpose:
c ^^^ *
c
c       subroutine to assess a b - spline tableau, stored in vectorized
c       form.
c
c calling convention:
c ^^^^^^^^^
c
c       basis ( m, n, x, b, bl, q )
c
c meaning of parameters:
c ^^^^^^^^^^ *
c
c       m               ( i )   half order of the spline (degree 2 * m - 1),
c                               m > 0.0
c       n               ( i )   number of knots, n > = 2 * m.
c       x(n)            ( i )   knot sequence, x(i - 1) < x(i), i = 2, n.
c       b(1 - m:m - 1, n)    ( o )   output tableau. element b(j, i) of array
c                               b corresponds with element b(i, i + j) of
c                               the tableau matrix b.
c       bl              ( o )   l1 - norm of b.
c       q(1 - m:m)        ( w )   internal work array.
c
c remark:
c ^^^
c
c       this subroutine is an adaptation of subroutine basis from the
c       paper by lyche et al. (1983). no checking is performed on the
c       validity of m and n. if the knot sequence is not strictly in - 
c       creasing, division by zero may occur.
c
c reference:
c ^^^^ *
c
c       t. lyche, l.l. schumaker, k. sepehrnoori, fortran subroutines
c       for computing smoothing and interpolating natural splines.
c       advances in engineering software 5(1983)1, pp. 2 - 5.0
c
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ *
c
      end
c
      subroutine basis(m, n, x, b, bl, q)
      implicit double precision (o - z, a - h)
      parameter (zero = 0d0, one = 1d0)
c
      dimension x(n), b(1 - m:m - 1, n), q(1 - m:m)
c^ *         linear spline
      if (m == 1) 
      for 3 i = 1: n
      b(0, i) = one
    3 continue
      bl = one
      return 
c
c^ *  general splines
c
      end
      mm1 = m - 1
      mp1 = m + 1
      m2 = 2 * m
c^ *     1st row
      for 15 l = 1: n
      for 5 j = - mm1: m
      q(j) = zero
    5 continue
      q(mm1) = one
c^ *     successive rows
      if ((l != 1) && (l != n)) q(mm1) = one / (x(l + 1) - x(l - 1)) end
      arg = x(l)
      for 13 i = 3: m2
      ir = mp1 - i
      v = q(ir)
c^ *               left - hand b - splines
      if (l < i) 
      for 6 j = l + 1: i
      u = v
      v = q(ir + 1)
      q(ir) = u + ((x(j) - arg) * v)
      ir = ir + 1
    6 continue
      end
      j1 = max0((l - i) + 1, 1)
      j2 = min0(l - 1, n - i)
c^ *               ordinary b - splines
      if (j1 <= j2) 
      if (i < m2) 
      for 8 j = j1: j2
      y = x(i + j)
      u = v
      v = q(ir + 1)
      q(ir) = u + (((v - u) * (y - arg)) / (y - x(j)))
      ir = ir + 1
    8 continue
      else
      for 10 j = j1: j2
      u = v
      v = q(ir + 1)
      q(ir) = ((arg - x(j)) * u) + ((x(i + j) - arg) * v)
      ir = ir + 1
   10 continue
      end
      end
      nmip1 = (n - i) + 1
c^ *           right - hand b - splines
      if (nmip1 < l) 
      for 12 j = nmip1: l - 1
      u = v
      v = q(ir + 1)
      q(ir) = ((arg - x(j)) * u) + v
      ir = ir + 1
   12 continue
      end
   13 continue
      for 14 j = - mm1: mm1
      b(j, l) = q(j)
   14 continue
c
c^ *  zero unused parts of b
c
   15 continue
      for 17 i = 1: mm1
      for 16 k = i: mm1
      b( - k, i) = zero
      b(k, (n + 1) - i) = zero
   16 continue
c
c^ *  assess l1 - norm of b
c
   17 continue
      bl = 0d0
      for 19 i = 1: n
      for 18 k = - mm1: mm1
      bl = bl + abs(b(k, i))
   18 continue
   19 continue
c
c^ *  ready
c
      bl = bl / n
      return 
c prep.for, 1985 - 07 - 04
c
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ *
c
c subroutine prep (real * 8)
c
c purpose:
c ^^^ *
c
c       to compute the matrix we of weighted divided difference coeffi - 
c       cients needed to set up a linear system of equations for sol - 
c       ving b - spline smoothing problems, and its l1 - norm el. the matrix
c       we is stored in vectorized form.
c
c calling convention:
c ^^^^^^^^^
c
c       prep ( m, n, x, w, we, el )
c
c meaning of parameters:
c ^^^^^^^^^^ *
c
c       m               ( i )   half order of the b - spline (degree
c                               2 * m - 1), with m > 0.0
c       n               ( i )   number of knots, with n > = 2 * m.
c       x(n)            ( i )   strictly increasing knot array, with
c                               x(i - 1) < x(i), i = 2, n.
c       w(n)            ( i )   weight matrix (diagonal), with
c                               w(i) > 0.0, i = 1, n.
c       we( - m:m, n)      ( o )   array containing the weighted divided
c                               difference terms in vectorized format.
c                               element we(j, i) of array e corresponds
c                               with element e(i, i + j) of the matrix
c                               w^ - 1 * e.
c       el              ( o )   l1 - norm of we.
c
c remark:
c ^^^
c
c       this subroutine is an adaptation of subroutine prep from the pap
cer
c       by lyche et al. (1983). no checking is performed on the validity
c       of m and n. division by zero may occur if the knot sequence is
c       not strictly increasing.
c
c reference:
c ^^^^ *
c
c       t. lyche, l.l. schumaker, k. sepehrnoori, fortran subroutines
c       for computing smoothing and interpolating natural splines.
c       advances in engineering software 5(1983)1, pp. 2 - 5.0
c
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ *
c
      end
c
      subroutine prep(m, n, x, w, we, el)
      implicit double precision (o - z, a - h)
      parameter (zero = 0d0, one = 1d0)
c
c^ *  calculate the factor f1
c
cwe( - m:m, n)              
      dimension x(n), w(n), we(((2 * m) + 1) * n)
      m2 = 2 * m
      mp1 = m + 1
      m2m1 = m2 - 1
      m2p1 = m2 + 1
      nm = n - m
      f1 = - one
      if (m != 1) 
      for 5 i = 2: m
      f1 = - (f1 * i)
    5 continue
      for 6 i = mp1: m2m1
      f1 = f1 * i
    6 continue
c
c^ *  columnwise evaluation of the unweighted design matrix e
c
      end
      i1 = 1
      i2 = m
      jm = mp1
      for 17 j = 1: n
      inc = m2p1
      if (j > nm) 
      f1 = - f1
      f = f1
      else
      if (j < mp1) 
      inc = 1
      f = f1
      else
      f = f1 * (x(j + m) - x(j - m))
      end
      end
      if (j > mp1) i1 = i1 + 1 end
      if (i2 < n) i2 = i2 + 1 end
c^ *     loop for divided difference coefficients
      jj = jm
      ff = f
      y = x(i1)
      i1p1 = i1 + 1
      for 11 i = i1p1: i2
      ff = ff / (y - x(i))
   11 continue
      we(jj) = ff
      jj = jj + m2
      i2m1 = i2 - 1
      if (i1p1 <= i2m1) 
      for 14 l = i1p1: i2m1
      ff = f
      y = x(l)
      for 12 i = i1: l - 1
      ff = ff / (y - x(i))
   12 continue
      for 13 i = l + 1: i2
      ff = ff / (y - x(i))
   13 continue
      we(jj) = ff
      jj = jj + m2
   14 continue
      end
      ff = f
      y = x(i2)
      for 16 i = i1: i2m1
      ff = ff / (y - x(i))
   16 continue
      we(jj) = ff
      jj = jj + m2
      jm = jm + inc
c
c^ *  zero the upper left and lower right corners of e
c
   17 continue
      kl = 1
      n2m = (m2p1 * n) + 1
      for 19 i = 1: m
      ku = (kl + m) - i
      for 18 k = kl: ku
      we(k) = zero
      we(n2m - k) = zero
   18 continue
      kl = kl + m2p1
c
c^ *  weighted matrix we = w^ - 1 * e and its l1 - norm
c
   19 continue
   20 jj = 0
      el = 0d0
      for 22 i = 1: n
      wi = w(i)
      for 21 j = 1: m2p1
      jj = jj + 1
      we(jj) = we(jj) / wi
      el = el + abs(we(jj))
   21 continue
   22 continue
c
c^ *  ready
c
      el = el / n
      return 
c splc.for, 1985 - 12 - 12
c
c author: h.j. woltring
c
c organizations: university of nijmegen, and
c                philips medical systems, eindhoven
c                (the netherlands)
c
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ *
c
c function splc (real * 8)
c
c purpose:
c ^^^ *
c
c       to assess the coefficients of a b - spline and various statistical
c       parameters, for a given value of the regularization parameter p.
c
c calling convention:
c ^^^^^^^^^
c
c       fv = splc ( m, n, k, y, ny, wx, wy, mode, val, p, eps, c, nc,
c       1           stat, b, we, el, bwe)
c
c meaning of parameters:
c ^^^^^^^^^^ *
c
c       splc            ( o )   gcv function value if |mode| == 2,
c                               mse value if |mode| == 3, and absolute
c                               difference with the prior given number o
cf
c                               degrees of freedom if |mode| == 4.0
c       m               ( i )   half order of the b - spline (degree 2 * m - 1
c),
c                               with m > 0.0
c       n               ( i )   number of observations, with n > = 2 * m.
c       k               ( i )   number of datasets, with k > = 1.0
c       y(ny, k)         ( i )   observed measurements.
c       ny              ( i )   first dimension of y(ny, k), with ny >= n
c.
c       wx(n)           ( i )   weight factors, corresponding to the
c                               relative inverse variance of each measur
ce - 
c                               ment, with wx(i) > 0.0.0
c       wy(k)           ( i )   weight factors, corresponding to the
c                               relative inverse variance of each datase
ct,
c                               with wy(j) > 0.0.0
c       mode            ( i )   mode switch, as described in gcvspl.
c       val             ( i )   prior variance if |mode| == 3, and
c                               prior number of degrees of freedom if
c                               |mode| == 4.0 for other values of mode,
c                               val is not used.
c       p               ( i )   smoothing parameter, with p > = 0.0.0 if
c                               p == 0.0, an interpolating spline is
c                               calculated.
c       eps             ( i )   relative rounding tolerance * 10.0.0 eps is
c                               the smallest positive number such that
c                               eps / 10.0 + 1.0 != 1.0.0
c       c(nc, k)         ( o )   calculated spline coefficient arrays. nb
c:
c                               the dimensions of in gcvspl and in splde
cr
c                               are different# In SPLDER, only a single
c                               column of c(n, k) is needed, and the prop
cer
c                               column c(1, j), with j = 1.0..k, should be u
csed
c                               when calling splder.
c       nc              ( i )   first dimension of c(nc, k), with nc >= n
c.
c       stat(6)         ( o )   statistics array. see the description in
c                               subroutine gcvspl.
c       b (1 - m:m - 1, n)   ( i )   b - spline tableau as evaluated by subrout
cine
c                               basis.
c       we( - m:m  , n)   ( i )   weighted b - spline tableau (w^ - 1 * e) as
c                               evaluated by subroutine prep.
c       el              ( i )   l1 - norm of the matrix we as evaluated by
c                               subroutine prep.
c       bwe( - m:m, n)     ( o )   central 2 * m + 1 bands of the inverted
c                               matrix ( b  +  p * w^ - 1 * e )^ - 1
c
c remarks:
c ^^^ *
c
c       this subroutine combines elements of subroutine splc0 from the
c       paper by lyche et al. (1983), and of subroutine spfit1 by
c       hutchinson (1985).
c
c references:
c ^^^^^
c
c       m.f. hutchinson (1985), subroutine cubgcv. csiro division of
c       mathematics and statistics, p.o. box 1965, canberra, act 2601,
c       australia.
c
c       t. lyche, l.l. schumaker, k. sepehrnoori, fortran subroutines
c       for computing smoothing and interpolating natural splines.
c       advances in engineering software 5(1983)1, pp. 2 - 5.0
c
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ *
c
      end
c
      function splc(m, n, k, y, ny, wx, wy, mode, val, p, eps, c, nc, stat, b, we, el, bwe)
      implicit double precision (o - z, a - h)
      parameter (zero = 0d0, one = 1d0, two = 2d0)
c
c^ *  check on p - value
c
      dimension y(ny, k), wx(n), wy(k), c(nc, k), stat(6), b(1 - m:m - 1, n), we( - m:m, n), bwe( - m:m, n)
      dp = p
      stat(4) = p
c^ *  pseudo - interpolation if p is too small
      pel = p * el
      if (pel < eps) 
      dp = eps / el
      stat(4) = zero
c^ *  pseudo least - squares polynomial if p is too large
      end
      if ((pel * eps) > one) 
      dp = one / (el * eps)
      stat(4) = dp
c
c^ *  calculate  bwe  =  b  +  p * w^ - 1 * e
c
      end
      for 40 i = 1: n
      km = - min0(m, i - 1)
      kp = min0(m, n - i)
      for 30 l = km: kp
      if (iabs(l) == m) 
      bwe(l, i) = dp * we(l, i)
      else
      bwe(l, i) = b(l, i) + (dp * we(l, i))
      end
   30 continue
c
c^ *  solve bwe * c = y, and assess trace [ b * bwe^ - 1 ]
c
   40 continue
      bandet(bwe, m, n)
      bansol(bwe, y, ny, c, nc, m, n, k)
ctrace * p = res. d.o.
      stat(3) = trinv(we, bwe, m, n) * dp
c
c^ *  compute mean - squared weighted residual
c
      trn = stat(3) / n
      esn = zero
      for 70 j = 1: k
      for 60 i = 1: n
      dt = - y(i, j)
      km = - min0(m - 1, i - 1)
      kp = min0(m - 1, n - i)
      for 50 l = km: kp
      dt = dt + (b(l, i) * c(i + l, j))
   50 continue
      esn = esn + (((dt * dt) * wx(i)) * wy(j))
   60 continue
   70 continue
c
c^ *  calculate statistics and function value
c
      esn = esn / (n * k)
cestimated variance               
      stat(6) = esn / trn
cgcv function value               
      stat(1) = stat(6) / trn
c     stat(3) = trace [p * b * bwe^ - 1] #Estimated residuals" d.o.f.
c     stat(4) = p                     #Normalized smoothing factor
cmean squared residual            
      stat(2) = esn
c^ *     unknown variance: gcv
      if (iabs(mode) != 3) 
      stat(5) = stat(6) - esn
      if (iabs(mode) == 1) splc = zero end
      if (iabs(mode) == 2) splc = stat(1) end
      if (iabs(mode) == 4) splc = dabs(stat(3) - val) end
c^ *     known variance: estimated mean squared error
      else
      stat(5) = esn - (val * ((two * trn) - one))
      splc = stat(5)
c
      end
      return 
c bandet.for, 1985 - 06 - 03
c
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ *
c
c subroutine bandet (real * 8)
c
c purpose:
c ^^^ *
c
c       this subroutine computes the lu decomposition of an n * n matrix
c       e. it is assumed that e has m bands above and m bands below the
c       diagonal. the decomposition is returned in e. it is assumed that
c       e can be decomposed without pivoting. the matrix e is stored in
c       vectorized form in the array e( - m:m, n), where element e(j, i) of
c       the array e corresponds with element e(i, i + j) of the matrix e.
c
c calling convention:
c ^^^^^^^^^
c
c       bandet ( e, m, n )
c
c meaning of parameters:
c ^^^^^^^^^^ *
c
c       e( - m:m, n)       (i / o)   matrix to be decomposed.
c       m, n            ( i )   matrix dimensioning parameters,
c                               m > = 0, n > = 2 * m.
c
c remark:
c ^^^
c
c       no checking on the validity of the input data is performed.
c       if (m <= 0), no action is taken.
c
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ *
c
      end
c
      subroutine bandet(e, m, n)
      implicit double precision (o - z, a - h)
c
      dimension e( - m:m, n)
      if (m <= 0) return 
      for 40 i = 1: n
      di = e(0, i)
      mi = min0(m, i - 1)
      if (mi >= 1) 
      for 10 k = 1: mi
      di = di - (e( - k, i) * e(k, i - k))
   10 continue
      e(0, i) = di
      end
      lm = min0(m, n - i)
      if (lm >= 1) 
      for 30 l = 1: lm
      dl = e( - l, i + l)
      km = min0(m - l, i - 1)
      if (km >= 1) 
      du = e(l, i)
      for 20 k = 1: km
      du = du - (e( - k, i) * e(l + k, i - k))
      dl = dl - (e(( - l) - k, l + i) * e(k, i - k))
   20 continue
      e(l, i) = du
      end
      e( - l, i + l) = dl / di
   30 continue
      end
c
c^ *  ready
c
   40 continue
      return 
c bansol.for, 1985 - 12 - 12
c
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ *
c
c subroutine bansol (real * 8)
c
c purpose:
c ^^^ *
c
c       this subroutine solves systems of linear equations given an lu
c       decomposition of the design matrix. such a decomposition is pro - 
c       vided by subroutine bandet, in vectorized form. it is assumed
c       that the design matrix is not singular.
c
c calling convention:
c ^^^^^^^^^
c
c       bansol ( e, y, ny, c, nc, m, n, k )
c
c meaning of parameters:
c ^^^^^^^^^^ *
c
c       e( - m:m, n)       ( i )   input design matrix, in lu - decomposed,
c                               vectorized form. element e(j, i) of the
c                               array e corresponds with element
c                               e(i, i + j) of the n * n design matrix e.
c       y(ny, k)         ( i )   right hand side vectors.
c       c(nc, k)         ( o )   solution vectors.
c       ny, nc, m, n, k ( i )   dimensioning parameters, with m > = 0,
c                               n > 2 * m, and k > = 1.0
c
c remark:
c ^^^
c
c       this subroutine is an adaptation of subroutine bansol from the
c       paper by lyche et al. (1983). no checking is performed on the
c       validity of the input parameters and data. division by zero may
c       occur if the system is singular.
c
c reference:
c ^^^^ *
c
c       t. lyche, l.l. schumaker, k. sepehrnoori, fortran subroutines
c       for computing smoothing and interpolating natural splines.
c       advances in engineering software 5(1983)1, pp. 2 - 5.0
c
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ *
c
      end
c
      subroutine bansol(e, y, ny, c, nc, m, n, k)
      implicit double precision (o - z, a - h)
c
c^ *  check on special cases: m = 0, m = 1, m>1
c
      dimension e( - m:m, n), y(ny, k), c(nc, k)
      nm1 = n - 1
c
c^ *  m = 0: diagonal system
c
      if (m - 1) 10, 40, 80
   10 for 30 i = 1: n
      for 20 j = 1: k
      c(i, j) = y(i, j) / e(0, i)
   20 continue
   30 continue
c
c^ *  m = 1: tridiagonal system
c
      return 
   40 for 70 j = 1: k
      c(1, j) = y(1, j)
cforward sweep                          
      for 50 i = 2: n
      c(i, j) = y(i, j) - (e( - 1, i) * c(i - 1, j))
   50 continue
      c(n, j) = c(n, j) / e(0, n)
cbackward sweep                          
      for 60 i = nm1, 1: - 1
      c(i, j) = (c(i, j) - (e(1, i) * c(i + 1, j))) / e(0, i)
   60 continue
   70 continue
c
c^ *  m > 1: general system
c
      return 
   80 for 130 j = 1: k
      c(1, j) = y(1, j)
cforward sweep                         
      for 100 i = 2: n
      mi = min0(m, i - 1)
      d = y(i, j)
      for 90 l = 1: mi
      d = d - (e( - l, i) * c(i - l, j))
   90 continue
      c(i, j) = d
  100 continue
      c(n, j) = c(n, j) / e(0, n)
cbackward sweep                         
      for 120 i = nm1, 1: - 1
      mi = min0(m, n - i)
      d = c(i, j)
      for 110 l = 1: mi
      d = d - (e(l, i) * c(i + l, j))
  110 continue
      c(i, j) = d / e(0, i)
  120 continue
  130 continue
c
      return 
c trinv.for, 1985 - 06 - 03
c
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ *
c
c function trinv (real * 8)
c
c purpose:
c ^^^ *
c
c       to calculate trace [ b * e^ - 1 ], where b and e are n * n
c       matrices with bandwidth 2 * m + 1, and where e is a regular matrix
c       in lu - decomposed form. b and e are stored in vectorized form,
c       compatible with subroutines bandet and bansol.
c
c calling convention:
c ^^^^^^^^^
c
c       trace = trinv ( b, e, m, n )
c
c meaning of parameters:
c ^^^^^^^^^^ *
c
c       b( - m:m, n)       ( i ) input array for matrix b. element b(j, i)
c                             corresponds with element b(i, i + j) of the
c                             matrix b.
c       e( - m:m, n)       (i / o) input array for matrix e. element e(j, i)
c                             corresponds with element e(i, i + j) of the
c                             matrix e. this matrix is stored in lu - 
c                             decomposed form, with l unit lower tri - 
c                             angular, and u upper triangular. the unit
c                             diagonal of l is not stored. upon return,
c                             the array e holds the central 2 * m + 1 bands
c                             of the inverse e^ - 1, in similar ordering.
c       m, n            ( i ) array and matrix dimensioning parameters
c                             (m > 0, n >= 2 * m + 1).
c       trinv           ( o ) output function value trace [ b * e^ - 1 ]
c
c reference:
c ^^^^ *
c
c       a.m. erismanw.f. tinney, on computing certain elements of the
c       inverse of a sparse matrix. communications of the acm 18(1975),
c       nr. 3, pp. 177 - 179.0
c
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ *
c
      end
c
      double precision function trinv(b, e, m, n)
      implicit double precision (o - z, a - h)
      parameter (zero = 0d0, one = 1d0)
c
c^ *  assess central 2 * m + 1 bands of e^ - 1 and store in array e
c
      dimension b( - m:m, n), e( - m:m, n)
cnth pivot                             
      e(0, n) = one / e(0, n)
      for 40 i = n - 1, 1: - 1
      mi = min0(m, n - i)
c^ *     save ith column of l and ith row of u, and normalize u row
cith pivot                             
      dd = one / e(0, i)
      for 10 k = 1: mi
cith row of u (normalized)    
      e(k, n) = e(k, i) * dd
cith column of l                   
      e( - k, 1) = e( - k, k + i)
   10 continue
c^ *     invert around ith pivot
      dd = dd + dd
      for 30 j = mi, 1: - 1
      du = zero
      dl = zero
      for 20 k = 1: mi
      du = du - (e(k, n) * e(j - k, i + k))
      dl = dl - (e( - k, 1) * e(k - j, i + j))
   20 continue
      e(j, i) = du
      e( - j, j + i) = dl
      dd = dd - ((e(j, n) * dl) + (e( - j, 1) * du))
   30 continue
      e(0, i) = 5d - 1 * dd
c
c^ *  assess trace [ b * e^ - 1 ] and clear working storage
c
   40 continue
      dd = zero
      for 60 i = 1: n
      mn = - min0(m, i - 1)
      mp = min0(m, n - i)
      for 50 k = mn: mp
      dd = dd + (b(k, i) * e( - k, k + i))
   50 continue
   60 continue
      trinv = dd
      for 70 k = 1: m
      e(k, n) = zero
      e( - k, 1) = zero
c
c^ *  ready
c
   70 continue
      return 
c splder.for, 1985 - 06 - 11
c
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ *
c
c function splder (real * 8)
c
c purpose:
c ^^^ *
c
c       to produce the value of the function (ider == 0) or of the
c       iderth derivative (ider > 0) of a 2m - th order b - spline at
c       the point t. the spline is described in terms of the half
c       order m, the knot sequence x(n), n >= 2 * m, and the spline
c       coefficients c(n).
c
c calling convention:
c ^^^^^^^^^
c
c       svider = splder ( ider, m, n, t, x, c, l, q )
c
c meaning of parameters:
c ^^^^^^^^^^ *
c
c       splder  ( o )   function or derivative value.
c       ider    ( i )   derivative order required, with 0 <= ider
c                       and ider <= 2 * m. if ider == 0, the function
c                       value is returned; otherwise, the ider - th
c                       derivative of the spline is returned.
c       m       ( i )   half order of the spline, with m > 0.0
c       n       ( i )   number of knots and spline coefficients,
c                       with n >= 2 * m.
c       t       ( i )   argument at which the spline or its deri - 
c                       vative is to be evaluated, with x(1) <= t
c                       and t <= x(n).
c       x(n)    ( i )   strictly increasing knot sequence array,
c                       x(i - 1) < x(i), i = 2, ..., n.
c       c(n)    ( i )   spline coefficients, as evaluated by
c                       subroutine gvcspl.
c       l       (i / o)   l contains an integer such that:
c                       x(l) <= t and t < x(l + 1) if t is within
c                       the range x(1) <= t and t < x(n). if
c                       t < x(1), l is set to 0, and if t >= x(n),
c                       l is set to n. the search for l is facili - 
c                       tated if l has approximately the right
c                       value on entry.
c       q(2 * m)  ( w )   internal work array.
c
c remark:
c ^^^
c
c       this subroutine is an adaptation of subroutine splder of
c       the paper by lyche et al. (1983). no checking is performed
c       on the validity of the input parameters.
c
c reference:
c ^^^^ *
c
c       t. lyche, l.l. schumaker, k. sepehrnoori, fortran subroutines
c       for computing smoothing and interpolating natural splines.
c       advances in engineering software 5(1983)1, pp. 2 - 5.0
c
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ *
c
      end
c
      double precision function splder(ider, m, n, t, x, c, l, q)
      implicit double precision (o - z, a - h)
      parameter (zero = 0d0, one = 1d0)
cf2py intent(out) :: splder
cf2py intent(hide) :: q
c
c^ *  derivatives of ider >= 2 * m are alway zero
c
      dimension x(n), c(n), q(2 * m)
      m2 = 2 * m
      k = m2 - ider
      if (k < 1) 
      splder = zero
      return 
c
c^ *  search for the interval value l
c
      end
c
c^ *  initialize parameters and the 1st row of the b - spline
c^ *  coefficients tableau
c
      search(n, x, t, l)
      tt = t
      mp1 = m + 1
      npm = n + m
      m2m1 = m2 - 1
      k1 = k - 1
      nk = n - k
      lk = l - k
      lk1 = lk + 1
      lm = l - m
      jl = l + 1
      ju = l + m2
      ii = n - m2
      ml = - l
      for 2 j = jl: ju
      if ((j >= mp1) && (j <= npm)) 
      q(j + ml) = c(j - m)
      else
      q(j + ml) = zero
      end
c
c^ *  the following loop computes differences of the b - spline
c^ *  coefficients. if the value of the spline is required,
c^ *  differencing is not necessary.
c
    2 continue
      if (ider > 0) 
      jl = jl - m2
      ml = ml + m2
      for 6 i = 1: ider
      jl = jl + 1
      ii = ii + 1
      j1 = max0(1, jl)
      j2 = min0(l, ii)
      mi = m2 - i
      j = j2 + 1
      if (j1 <= j2) 
      for 3 jin = j1: j2
      j = j - 1
      jm = ml + j
      q(jm) = (q(jm) - q(jm - 1)) / (x(j + mi) - x(j))
    3 continue
      end
      if (jl >= 1) #FIXME goto 6
      i1 = i + 1
      j = ml + 1
      if (i1 <= ml) 
      for 5 jin = i1: ml
      j = j - 1
      q(j) = - q(j - 1)
    5 continue
      end
    6 continue
      for 7 j = 1: k
      q(j) = q(j + ider)
    7 continue
c
c^ *  compute lower half of the evaluation tableau
c
      end
ctableau ready if ider == 2 * m - 1            
      if (k1 >= 1) 
      for 14 i = 1: k1
      nki = nk + i
      ir = k
      jj = l
      ki = k - i
c^ *        right - hand b - splines
      nki1 = nki + 1
      if (l >= nki1) 
      for 9 j = nki1: l
      q(ir) = q(ir - 1) + ((tt - x(jj)) * q(ir))
      jj = jj - 1
      ir = ir - 1
    9 continue
c^ *        middle b - splines
      end
      lk1i = lk1 + i
      j1 = max0(1, lk1i)
      j2 = min0(l, nki)
      if (j1 <= j2) 
      for 11 j = j1: j2
      xjki = x(jj + ki)
      z = q(ir)
      q(ir) = z + (((xjki - tt) * (q(ir - 1) - z)) / (xjki - x(jj)))
      ir = ir - 1
      jj = jj - 1
   11 continue
c^ *        left - hand b - splines
      end
      if (lk1i <= 0) 
      jj = ki
      lk1i1 = 1 - lk1i
      for 13 j = 1: lk1i1
      q(ir) = q(ir) + ((x(jj) - tt) * q(ir - 1))
      jj = jj - 1
      ir = ir - 1
   13 continue
      end
   14 continue
c
c^ *  compute the return value
c
      end
c^ *  multiply with factorial if ider > 0
      z = q(k)
      if (ider > 0) 
      for 16 j = k: m2m1
      z = z * j
   16 continue
      end
c
c^ *  ready
c
      splder = z
      return 
c search.for, 1985 - 06 - 03
c
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ *
c
c subroutine search (real * 8)
c
c purpose:
c ^^^ *
c
c       given a strictly increasing knot sequence x(1) < ... < x(n),
c       where n > = 1, and a real number t, this subroutine finds the
c       value l such that x(l) < = t < x(l + 1).  if t < x(1), l = 0; end
c       if x(n) < = t, l = n. end
c
c calling convention:
c ^^^^^^^^^
c
c       search ( n, x, t, l )
c
c meaning of parameters:
c ^^^^^^^^^^ *
c
c       n       ( i )   knot array dimensioning parameter.
c       x(n)    ( i )   stricly increasing knot array.
c       t       ( i )   input argument whose knot interval is to
c                       be found.
c       l       (i / o)   knot interval parameter. the search procedure
c                       is facilitated if l has approximately the
c                       right value on entry.
c
c remark:
c ^^^
c
c       this subroutine is an adaptation of subroutine search from
c       the paper by lyche et al. (1983). no checking is performed
c       on the input parameters and data; the algorithm may fail if
c       the input sequence is not strictly increasing.
c
c reference:
c ^^^^ *
c
c       t. lyche, l.l. schumaker, k. sepehrnoori, fortran subroutines
c       for computing smoothing and interpolating natural splines.
c       advances in engineering software 5(1983)1, pp. 2 - 5.0
c
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ *
c
      end
c
      subroutine search(n, x, t, l)
      implicit double precision (o - z, a - h)
c
      dimension x(n)
c^ *     out of range to the left
      if (t < x(1)) 
      l = 0
      return 
      end
c^ *     out of range to the right
      if (t >= x(n)) 
      l = n
      return 
c^ *  validate input value of l
      end
      l = max0(l, 1)
c
c^ *  often l will be in an interval adjoining the interval found
c^ *  in a previous to search
c
      if (l >= n) l = n - 1 end
      if (t >= x(l)) #FIXME goto 5
      l = l - 1
c
c^ *  perform bisection
c
      if (t >= x(l)) return 
      il = 1
    3 iu = l
    4 l = (il + iu) / 2
      if ((iu - il) <= 1) return 
      if (t < x(l)) #FIXME goto 3
      il = l
      #FIXME goto 4
    5 if (t < x(l + 1)) return 
      l = l + 1
      if (t < x(l + 1)) return 
      il = l + 1
      iu = n
c
      #FIXME goto 4
      end
