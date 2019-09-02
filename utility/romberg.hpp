#ifndef ROMBERG_HPP
#define ROMBERG_HPP

#include <vector>

#define Doub double
#define Int int
#define VecDoub std::vector<double>
#define Bool bool

#define MIN std::min
#define MAX std::max

#ifndef ITERMAX
#define ITERMAX 20
#endif

#ifndef MULT_STEP
#define MULT_STEP (1.0/9.0)
#endif

#ifndef ORDER_K
#define ORDER_K 5
#endif

typedef const VecDoub VecDoub_I;


// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// # # # INTERPOLATION # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

struct Poly_interp
{

  Int n, mm, jsav, cor, dj;
  const Doub *xx, *yy;
  Doub dy;

  Poly_interp(VecDoub_I &xv, VecDoub_I &yv, Int m):
    n(xv.size()), mm(m), jsav(0), cor(0), xx(&xv[0]), yy(&yv[0]), dy(0.)
    {
      dj = MIN(1,(int)pow((Doub)n,0.25));
    }

  Doub interp(Doub x)
  {
    Int jlo = cor ? hunt(x) : locate(x);
    return rawinterp(jlo,x);
  }

  Int locate(const Doub x)
  {
    Int ju,jm,jl;
    if (n < 2 || mm < 2 || mm > n) throw("locate size error");
    Bool ascnd=(xx[n-1] >= xx[0]);
    jl=0;
    ju=n-1;
    while (ju-jl > 1)
    {
      jm = (ju+jl) >> 1;
      if (x >= xx[jm] == ascnd)
        jl=jm;
      else
        ju=jm;
    }
    cor = abs(jl-jsav) > dj ? 0 : 1;
    jsav = jl;
    return MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
  }

  Int hunt(const Doub x)
  {
    Int jl=jsav, jm, ju, inc=1;
    if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
    Bool ascnd=(xx[n-1] >= xx[0]);
    if (jl < 0 || jl > n-1)
    {
      jl=0;
      ju=n-1;
    }
    else
    {
      if (x >= xx[jl] == ascnd)
      {
        for (;;)
        {
          ju = jl + inc;
          if (ju >= n-1) { ju = n-1; break;}
          else if (x < xx[ju] == ascnd) break;
          else
          {
            jl = ju;
            inc += inc;
          }
        }
      }
      else
      {
        ju = jl;
        for (;;)
        {
          jl = jl - inc;
          if (jl <= 0) { jl = 0; break;}
          else if (x >= xx[jl] == ascnd) break;
          else
          {
            ju = jl;
            inc += inc;
          }
        }
      }
    }
    while (ju-jl > 1)
    {
      jm = (ju+jl) >> 1;
      if (x >= xx[jm] == ascnd)
      jl=jm; else
      ju=jm;
    }
    cor = abs(jl-jsav) > dj ? 0 : 1;
    jsav = jl;
    return MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
  }

  Doub rawinterp(Int jl, Doub x)
  {
    Int i,m,ns=0;
    Doub y,den,dif,dift,ho,hp,w;
    const Doub *xa = &xx[jl], *ya = &yy[jl];
    VecDoub c(mm),d(mm);
    dif=abs(x-xa[0]);
    for (i=0;i<mm;i++)
    {
      if ((dift=abs(x-xa[i])) < dif)
      {
        ns=i;
        dif=dift;
      }
      c[i]=ya[i];
      d[i]=ya[i];
    }
    y=ya[ns--];
    for (m=1;m<mm;m++)
    {
      for (i=0;i<mm-m;i++)
      {
        ho=xa[i]-x;
        hp=xa[i+m]-x;
        w=c[i+1]-d[i];
        if ((den=ho-hp) == 0.0) throw("Poly_interp error");
        den=w/den;
        d[i]=hp*den;
        c[i]=ho*den;
      }
      y += (dy=(2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--]));
    }
    return y;
  }

};

// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// # # # QUADRATURE  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

struct Quadrature
{
  Int n;
  virtual Doub next() = 0;
};

template <class T>
struct Midpnt : Quadrature
{
  //Routine implementing the extended midpoint rule.
  Doub a,b,s;   // Limits of integration and current value of inte-
  T &funk;      // gral.
  Midpnt(T &funcc, const Doub aa, const Doub bb) :
    funk(funcc), a(aa), b(bb) { n=0; }
  // The constructor takes as inputs func, the function or functor to be integrated between
  // limits a and b, also input.
  Doub next()
  {
  Int it,j;
    Doub x,tnm,sum,del,ddel;
    n++;
    if (n == 1)
    {
        return (s=(b-a)*func(0.5*(a+b)));
    }
    else
    {
      for(it=1,j=1;j<n-1;j++) it *= 3;
      tnm=it;
      del=(b-a)/(3.0*tnm);
      ddel=del+del;
      x=a+0.5*del;
      sum=0.0;
      for (j=0;j<it;j++)
      {
        sum += func(x);
        x += ddel;
        sum += func(x);
        x += del;
      }
      s=(s+(b-a)*sum/tnm)/3.0;
      return s;
    }
  }
  // The added points alternate in spacing be- tween del and ddel.
  // The new sum is combined with the old inte- gral to give a refined integral.
  virtual Doub func(const Doub x)
  {
    return funk(x);
  }
};

template <class T>
struct Midinf : Midpnt<T>
{
  // This routine is an exact replacement for midpnt, i.e., returns the nth stage
  // of refinement of the integral of funcc from aa to bb, except that the
  // function is evaluated at evenly spaced points in 1=x rather than in x. This
  // allows the upper limit bb to be as large and positive as the computer allows,
  // or the lower limit aa to be as large and negative, but not both. aa and bb
  // must have the same sign.
  Doub func(const Doub x)
  {
    return Midpnt<T>::funk(1.0/x)/(x*x);
  }
  Midinf(T &funcc, const Doub aa, const Doub bb) :
    Midpnt<T>(funcc, aa, bb)
    {
      Midpnt<T>::a = 1.0/bb;
      Midpnt<T>::b = 1.0/aa;
    }
};

template<class T>
struct Trapzd : Quadrature
{
  T &func;
  Doub a,b,s;
  Trapzd() {};
  Trapzd(T &funcc, const Doub aa, const Doub bb) :
    func(funcc), a(aa), b(bb) {n=0;}
  Doub next()
  {
    Doub x,tnm,sum,del;
    Int it,j;
    n++;
    if (n == 1)
    {
      return (s=0.5*(b-a)*(func(a)+func(b)));
    }
    else
    {
      for (it=1,j=1;j<n-1;j++) it <<= 1;
      tnm=it;
      del=(b-a)/tnm;
      x=a+0.5*del;
      for (sum=0.0,j=0;j<it;j++,x+=del) sum += func(x);
      s=0.5*(s+(b-a)*sum/tnm);
      return s;
    }
  }
};

// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// # # # ROMBERG TEMPLATE CLASS  # # # # # # # # # # # # # # # # # # # # # # # #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

template <class T>
Doub qromb(T &func, Doub a, Doub b, const Doub eps=1.0e-10)
{
  const Int JMAX=ITERMAX, JMAXP=JMAX+1, K=ORDER_K;
  VecDoub s(JMAX), h(JMAXP);
  Poly_interp polint(h,s,K);
  h[0]=1.0;
  Trapzd<T> t(func,a,b);
  Doub ss = 0.0;
  for (Int j=1; j<=JMAX; j++)
  {
    s[j-1]=t.next();
    if (j >= K)
    {
      ss=polint.rawinterp(j-K,0.0);
      if (abs(polint.dy) <= eps*abs(ss)) return ss;
    }
    h[j]=MULT_STEP*h[j-1];
  }
  // std::cout << "WARNING: Too many steps in routine qromb" << std::endl;
  return ss;
}

template<class T>
Doub qromo(Midpnt<T> &q, const Doub eps=3.0e-9)
{
  const Int JMAX=ITERMAX, JMAXP=JMAX+1, K=ORDER_K;
  VecDoub h(JMAXP), s(JMAX);
  Poly_interp polint(h,s,K);
  h[0]=1.0;
  Doub ss = 0.0;
  for (Int j=1;j<=JMAX;j++)
  {
    s[j-1]=q.next();
    if (j >= K)
    {
      ss=polint.rawinterp(j-K,0.0);
      if (abs(polint.dy) <= eps*abs(ss)) return ss;
    }
    h[j]=MULT_STEP*h[j-1];
  }
  std::cout << "WARNING: Too many steps in routine qromb" << std::endl;
  return ss;
}

// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#endif /* ROMBERG_HPP */
