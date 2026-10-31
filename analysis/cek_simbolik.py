import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_HERE, _os.path.join(_HERE, "..", "experiments")]
DATA = _os.path.join(_HERE, "..", "data") + _os.sep
MANUSCRIPT = _os.path.join(_HERE, "..", "manuscript")
_os.chdir(MANUSCRIPT)  # figures and JSON statistics are written to manuscript/gambar and manuscript/gambar_en
import sympy as sp
th, al, psi, C, R, c, a, rho, p, q = sp.symbols('theta alpha psi C R c a rho p q', real=True)
E = lambda f: sp.simplify(sp.integrate(f, (th, 0, 2*sp.pi)) / (2*sp.pi))
print("E cos, cos^2, cos^3, cos^4:", E(sp.cos(th)), E(sp.cos(th)**2), E(sp.cos(th)**3), E(sp.cos(th)**4))
y = C + R*sp.cos(th)
mu2 = sp.expand(E(y**2)); mu4 = sp.expand(E(y**4))
print("mu2 =", mu2, "  mu4 =", mu4)
b = sp.simplify((mu2**2/mu4).subs(R, sp.sqrt(rho)*C))
print("b(rho) =", sp.factor(b), " target equal:", sp.simplify(b - 2*(2+rho)**2/(3*rho**2+24*rho+8)) == 0)
db = sp.factor(sp.diff(2*(2+rho)**2/(3*rho**2+24*rho+8), rho))
print("b'(rho) =", db, " roots:", sp.solve(sp.numer(sp.together(db)), rho))
print("b(8/3) =", sp.nsimplify((2*(2+rho)**2/(3*rho**2+24*rho+8)).subs(rho, sp.Rational(8,3))),
      " lim0 =", sp.limit(2*(2+rho)**2/(3*rho**2+24*rho+8), rho, 0), " limInf =", sp.limit(2*(2+rho)**2/(3*rho**2+24*rho+8), rho, sp.oo))
# Prop 1 general form, m=2 check: R^2 = a^2 |u1 + u2 e^{-i psi}|^2
u1, u2 = sp.cos(al), sp.sin(al)
R2 = sp.simplify(sp.expand((u1 + u2*sp.cos(psi))**2 + (u2*sp.sin(psi))**2))
print("R^2/a^2 sesayap =", sp.simplify(sp.trigsimp(R2)), " vs 1+sin2a cos psi:", sp.simplify(R2 - (1+sp.sin(2*al)*sp.cos(psi))) == 0)
R2x = sp.simplify(sp.expand((u1 - u2*sp.cos(psi))**2 + (u2*sp.sin(psi))**2))
print("R'^2/a^2 beda sayap:", sp.simplify(R2x - (1-sp.sin(2*al)*sp.cos(psi))) == 0)
# direct check of projection identity: u1 cos th + u2 cos(th-psi) == |Z| cos(th - delta)
# mixture
K2, K4 = sp.symbols('K2 K4', positive=True)
m2 = p*K2*(1+q) + (1-p)*K2*(1-q); m4 = p*K4*(1+q)**2 + (1-p)*K4*(1-q)**2
print("mu2mix:", sp.factor(sp.expand(m2)), " mu4mix:", sp.expand(m4/K4))
bmix = m2**2/m4
print("b at q=-1:", sp.simplify(bmix.subs(q,-1)), " -> (1-p) b0 with b0=K2^2/K4:", sp.simplify(bmix.subs(q,-1) - (1-p)*K2**2/K4) == 0)
print("b at q=+1:", sp.simplify(bmix.subs(q,1)))
print("K2^2/K4 == b((a/c)^2):", sp.simplify(((c**2+a**2/2)**2/(c**4+3*c**2*a**2+3*a**4/8)) - (2*(2+rho)**2/(3*rho**2+24*rho+8)).subs(rho, a**2/c**2)) == 0)
# Sarle finite sample: limit of correction 3(n-1)^2/((n-2)(n-3)) -> 3
n = sp.symbols('n', positive=True)
print("koreksi Sarle ->", sp.limit(3*(n-1)**2/((n-2)*(n-3)), n, sp.oo))
