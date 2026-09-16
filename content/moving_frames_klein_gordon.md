# Moving Frames and the Klein–Gordon Equation

### An Illustration of the Fels–Olver Method in Field Theory

---

## Abstract

Physicists derive the free scalar field equation by writing down the most general Poincaré-invariant Lagrangian and applying the classical Euler–Lagrange operator.
This is efficient, but it hides *why* invariance forces the equation into the form it takes.
The theory of moving frames — developed by Fels and Olver, and extended to a full invariant variational calculus by Kogan and Olver — makes this mechanism explicit and algorithmic: one normalizes the group action on the jet bundle, reads off the differential invariants directly from the normalization, and reconstructs the equations of motion entirely in terms of those invariants.
This note carries out that program in the simplest field-theoretic setting — a single scalar field in $1+1$ spacetime dimensions under the Poincaré group $ISO(1,1)$ — and shows that it reproduces the Klein–Gordon equation without ever writing down a Lorentz-index contraction by hand.
The invariant of the theory, $\kappa=\sqrt{u_t^2-u_x^2}$, is *manufactured* by the frame rather than postulated, and the field equation falls out of an invariant Euler operator built from invariant differentiation alone.

---

## 1. Introduction

Sophus Lie observed, in the study of ordinary variational problems, that if a group $G$ leaves a functional invariant, the functional can always be rewritten purely in terms of the differential invariants of $G$. The Euler–Lagrange equations inherit $G$ as a symmetry group, and — when expressed in the same invariant language — take a form dictated entirely by the group's action on the jet bundle. Turning this observation into a working computational method requires three ingredients:

1. A **moving frame**: an equivariant map from the jet bundle into the group, constructed by normalizing a cross-section of the group orbits.
2. **Invariantization**: a projection sending any (non-invariant) differential function or form to its unique invariant counterpart agreeing with it on the cross-section.
3. **Recurrence relations**: formulas relating invariant differentiation to invariantized differentiation, which let one propagate the machinery to arbitrary jet order without ever integrating the group action explicitly.

This is the Fels–Olver moving frame method, and its extension to the calculus of variations — the *invariant variational bicomplex* — is due to Kogan and Olver. It is most often illustrated on curves: a single independent variable, a Euclidean or affine group acting on the ambient space, and invariants that turn out to be curvature and torsion. Here we illustrate the same machinery on a genuinely different structure: a *field* over a two-dimensional spacetime, where the group acts only on the independent variables and the field itself is inert (a Lorentz scalar). This is a natural and instructive generalization, because the field theorist's usual shortcut — "write down all Lorentz scalars you can build from $\partial_\mu u$" — is exactly what the moving frame produces mechanically, with no shortcut required.

---

## 2. Jets, Prolongation, and the Poincaré Group

Let $X=\mathbb{R}^{1,1}$ with coordinates $x^1=t,\ x^2=x$ and metric $\eta=\mathrm{diag}(+1,-1)$, and let $U=\mathbb{R}$ be the space of field values, coordinate $u$. The relevant jet space $J^n(X,U)$ has local coordinates $(t,x,u,u_t,u_x,u_{tt},u_{tx},u_{xx},\dots)$.

The Poincaré group $G=ISO(1,1)$ is three-dimensional, parametrized by $(a,b,\psi)$ — two translations and one boost rapidity — and acts on $X\times U$ by

$$
T = t\cosh\psi + x\sinh\psi + a,\qquad
X = t\sinh\psi + x\cosh\psi + b,\qquad
U = u.
\tag{2.1}
$$

The last equation is the defining feature of a **scalar field**: the group moves the point at which the field is evaluated, but does not mix the field's value with anything else. This is the field-theoretic analogue of a "curve with no attached frame," and it is the essential way in which the present problem differs from the classical moving-frames literature on framed curves, where the dependent variable typically *does* transform under the group (e.g., a frame angle shifting along with a rotation).

### 2.1 Infinitesimal generators

A basis for the Lie algebra of $G$ acting on $X\times U$ is

$$
v_1=\partial_t,\qquad v_2=\partial_x,\qquad v_3 = x\,\partial_t + t\,\partial_x,
\tag{2.2}
$$

with **no** component along $\partial_u$: translations and boosts move the field's argument, never its value.

### 2.2 Prolongation

Translations prolong trivially, $\mathrm{pr}\,v_1=v_1,\ \mathrm{pr}\,v_2=v_2$. For the boost, the standard multivariable prolongation formula,

$$
\varphi^i = D_iQ + \sum_j \xi^j u_{ij},\qquad Q = \varphi^u-\sum_j\xi^ju_j,
$$

applied to $\xi^t=x,\ \xi^x=t,\ \varphi^u=0$ (so $Q=-xu_t-tu_x$), gives

$$
\varphi^t = D_t(-xu_t-tu_x) + xu_{tt}+tu_{tx} = -u_x,\qquad
\varphi^x = D_x(-xu_t-tu_x) + xu_{tx}+tu_{xx} = -u_t.
$$

Hence

$$
\mathrm{pr}\,v_3 = x\,\partial_t + t\,\partial_x - u_x\,\partial_{u_t} - u_t\,\partial_{u_x} + \cdots
\tag{2.3}
$$

This is exactly the infinitesimal statement of the familiar transformation law $u_t\to u_t\cosh\psi-u_x\sinh\psi$, $u_x\to -u_t\sinh\psi+u_x\cosh\psi$ — but derived, not assumed.

---

## 3. Constructing the Moving Frame

The action (2.1) is transitive on $X$, with a one-dimensional stabilizer at each point (the boosts fixing that point). All three group parameters can therefore be pinned down using only first-order jet data: two translations fix the base point, and the leftover boost freedom is absorbed by a condition on the first derivatives.

**Cross-section.** Choose the normalization

$$
T=0,\qquad X=0,\qquad U_X = 0.
\tag{3.1}
$$

The first two conditions are the usual phantom translation invariants. The third normalizes the boost by sending the frame into the *local rest frame of the field gradient*: the frame in which the spatial derivative of $u$ vanishes. Using the prolonged transformation law for first derivatives,

$$
U_T = u_t\cosh\psi - u_x\sinh\psi,\qquad
U_X = -u_t\sinh\psi + u_x\cosh\psi,
$$

the condition $U_X=0$ fixes

$$
\tanh\psi = \frac{u_x}{u_t}
\quad\Longrightarrow\quad
\cosh\psi = \frac{u_t}{\sqrt{u_t^2-u_x^2}},\qquad
\sinh\psi = \frac{u_x}{\sqrt{u_t^2-u_x^2}},
\tag{3.2}
$$

valid on the open set where the gradient of $u$ is timelike, $u_t^2>u_x^2$. Together with the translation conditions, this defines the moving frame

$$
\rho(t,x,u^{(1)}) = (a,b,\psi),\ \\ \ \\
a=a(t,x,u_t,u_x),\ \\ \ \\
b = b(t,x,u_t,u_x),\ \\ \ \\
\psi=\operatorname{arctanh}(u_x/u_t).
\tag{3.3}
$$

No group parameter remains free: $\rho$ is an honest right-equivariant map from the (regular) jet space into $G$.

---

## 4. Invariants from Invariantization

By definition, invariantizing any differential function means substituting the frame (3.3) into its transformation law and evaluating on the cross-section. Applying this to $U_T$:

$$
\kappa := \iota(U_T) = u_t\cosh\psi - u_x\sinh\psi
= \frac{u_t^2 - u_x^2}{\sqrt{u_t^2-u_x^2}}
= \sqrt{u_t^2 - u_x^2}.
\tag{4.1}
$$

This is the **fundamental first-order differential invariant** of the problem. It is worth pausing on what just happened: nothing was assumed about Lorentz contractions. The quantity $u_t^2-u_x^2 = \eta^{\mu\nu}\partial_\mu u\,\partial_\nu u$ — normally introduced by fiat as "the" invariant built from the field gradient — appeared automatically as the unique quantity left over after the group's freedom was used up. Invariance did the work; no index gymnastics were needed.

The invariant horizontal coframe is likewise inherited from the frame: writing $\varpi^T=\iota(dT)$, $\varpi^X=\iota(dX)$, one finds

$$
\varpi^T = \kappa\, dt' ,\qquad \varpi^X = 0 \ \text{(on the cross-section)},
$$

where $dt'$ denotes the invariant arc element along the rest-frame time direction — the direct field-theoretic analogue of the invariant arc-length form $\varpi$ that appears throughout the classical moving-frames literature on curves.

---

## 5. Invariant Differential Operators

Dual to the invariant coframe is a pair of **invariant total derivative operators** — the field-theoretic analogue of the invariant arc-length derivative for curves, but doubled because there are now two independent variables rather than one:

$$
D_\parallel = \cosh\psi\,\partial_t - \sinh\psi\,\partial_x
= \frac{u_t\,\partial_t - u_x\,\partial_x}{\sqrt{u_t^2-u_x^2}},
\\ \ \\
D_\perp = -\sinh\psi\,\partial_t + \cosh\psi\,\partial_x
= \frac{-u_x\,\partial_t + u_t\,\partial_x}{\sqrt{u_t^2-u_x^2}}.
\tag{5.1}
$$

These two operators are an *invariant orthonormal frame* at every point of the jet space: $D_\parallel$ differentiates along the local rest-frame "time" direction defined by the field's own gradient, and $D_\perp$ along the orthogonal (rest-frame "space") direction. Neither operator refers to $t$ or $x$ individually — both are built entirely from $u_t,u_x$ and hence commute correctly with invariantization by construction.

### 5.1 Recurrence: the next invariant

Just as the curvature of a curve is generated by invariantly differentiating the frame's first invariant, here we generate a second-order invariant by applying $D_\parallel$ to $\kappa$:

$$
D_\parallel\kappa
= \frac{u_t^2 u_{tt} - 2u_tu_x u_{tx} + u_x^2 u_{xx}}{u_t^2 - u_x^2}.
\tag{5.2}
$$

It is important to note — exactly as one finds in the classical theory of curves, where the curvature's arc-length derivative is *not* simply the second normal-coordinate derivative — that (5.2) is **not** the d'Alembertian $\Box u = u_{tt}-u_{xx}$. Both are genuine Lorentz scalars built from the second-order jet, but they are independent invariants: $D_\parallel\kappa$ measures the rate of change of the gradient's magnitude along its own rest-frame direction, while $\Box u$ is the trace of the full Hessian against the Minkowski metric. Both live in the algebra of invariants generated by $\kappa,u$ under $D_\parallel,D_\perp$; the field equation, as we will see, is built from the latter.

---

## 6. The Invariant Variational Problem

By Lie's theorem, any $G$-invariant Lagrangian can be rewritten purely in terms of the differential invariants of $G$:

$$
\mathcal{L}[u] = \int L(t,x,u^{(n)})\,dt\,dx = \int \widetilde L\big(\kappa,u\big)\,\varpi,
\tag{6.1}
$$

for some function $\widetilde L$ of the invariant $\kappa$ and the (trivially invariant) field value $u$ — plus, if needed, higher invariants like $D_\parallel\kappa$. Since $\kappa=\sqrt{u_t^2-u_x^2}$ is only invariant up to sign under the full group (including reflections), a smooth, physically sensible $\widetilde L$ depends on $\kappa$ through the manifestly even combination $\kappa^2=u_t^2-u_x^2$.

Because the field carries no frame index — the dependent variable transforms trivially under $G$ — the invariant Euler operator here takes an unusually clean form: there are no "frame contact forms" contributing extra Lie-derivative corrections, since $v_\ell(\theta^u)=0$ for the translations and the boost's action on the basic contact form $\theta^u=du-u_t\,dt-u_x\,dx$ only feeds back through $u_t,u_x$ themselves, which are already accounted for in $\kappa$. The invariant Euler–Lagrange equation for a Lagrangian of the form $\widetilde L(\kappa^2,u)$ is

$$
E(\widetilde L) \;=\; \frac{\partial\widetilde L}{\partial u}
\;-\; 2\,\partial_t\!\left(\frac{\partial\widetilde L}{\partial(\kappa^2)}\,u_t\right)
\;+\; 2\,\partial_x\!\left(\frac{\partial\widetilde L}{\partial(\kappa^2)}\,u_x\right)
\;=\;0.
\tag{6.2}
$$

This is the invariant-variational-bicomplex analogue of the classical Eulerian operator $E_u(L)=\partial L/\partial u - D_tD_u^t(L)-D_xD_u^x(L)$, expressed entirely through the invariant $\kappa^2$ rather than through the individual (non-invariant) partial derivatives $u_t,u_x$.

---

## 7. Recovering the Klein–Gordon Equation

Take the simplest nontrivial choice consistent with the invariant structure just derived:

$$
\widetilde L(\kappa^2,u) = \tfrac12\kappa^2 - \tfrac12 m^2 u^2,
\qquad \frac{\partial\widetilde L}{\partial(\kappa^2)} = \tfrac12 \ \ (\text{constant}).
$$

Substituting into (6.2):

$$
E(\widetilde L) = -m^2u - \partial_t(u_t) + \partial_x(u_x) = -m^2 u - u_{tt} + u_{xx} = 0
$$

$$
\boxed{\ u_{tt} - u_{xx} + m^2 u = 0\ }
\tag{7.1}
$$

— the Klein–Gordon equation, obtained without ever writing $\eta^{\mu\nu}\partial_\mu\partial_\nu u+m^2u=0$ by hand. Every ingredient — the invariant $\kappa$, the invariant operators $D_\parallel,D_\perp$, and the invariant Euler operator (6.2) — was produced mechanically from the normalization (3.1)–(3.2). The only physics input was the choice of $\widetilde L$ as a function of the invariants; the *shape* of the resulting field equation was forced by the group.

---

## 8. Discussion

The scalar field case is deliberately the simplest possible illustration: because $u$ carries no representation index, the dependent variable is untouched by the group, and the invariant Euler operator collapses to the tidy expression (6.2). This is analogous to a curve with no attached frame — most of the machinery that makes the moving-frame method powerful for framed curves (the coupled system of invariant Eulerian and Hamiltonian operators, and the matrix of relative invariants relating contact two-forms to their invariant counterparts) is present in the general theory but *degenerates* here, precisely because there is no nontrivial group action on the fiber to generate it.

The natural next step is a field that **does** transform nontrivially under the group — a Lorentz vector field $A_\mu(t,x)$, for instance, or a field valued in a genuine representation of $SO(1,1)$ analogous to the $SO(n)$-valued frame of a framed curve. There, the boost acts simultaneously on the base point and mixes the field's own components, exactly mirroring the structure that makes framed-curve moving frames rich: nonzero Lie derivatives of the fiber's contact forms feed back into the invariant Eulerian through genuine matrix-valued corrections, and the invariant Euler–Lagrange equations take the fuller form seen in the classical theory of moving frames for curves and surfaces. That extension — carrying the same normalization philosophy through to Maxwell's equations in $1+1$ dimensions — is the natural continuation of this illustration.

## 9. Summary

| Classical shortcut | Moving-frame derivation |
|---|---|
| Assume $\mathcal L$ built from Lorentz scalars | $\kappa=\sqrt{u_t^2-u_x^2}$ emerges from normalizing the boost via $U_X=0$ |
| $\partial_\mu\partial^\mu u$ written by hand | $D_\parallel, D_\perp$ constructed as the invariant dual frame; $\Box u$ lives in the invariant algebra they generate |
| Euler–Lagrange via $\partial L/\partial u - \partial_\mu(\partial L/\partial(\partial_\mu u))$ | Invariant Euler operator (6.2), built entirely from $\kappa^2$ and invariant differentiation |
| Klein–Gordon equation postulated as "the" relativistic wave equation | Klein–Gordon equation *derived* as the unique output of the invariantized variational calculus for the simplest invariant Lagrangian |

The moving frame does not discover new physics here — it reproduces the textbook result — but it replaces an *assumption* of Lorentz invariance with a *construction* of it, term by term, directly on the jet bundle.

---

### References for the underlying method

- M. Fels and P. J. Olver, *Moving coframes II: Regularization and theoretical foundations*, Acta Appl. Math. **55** (1999), 127–208.
- I. Kogan and P. J. Olver, *The invariant variational bicomplex*, Contemp. Math. **285** (2001), 131–144.
- I. Kogan and P. J. Olver, *Invariant Euler–Lagrange equations and the invariant variational bicomplex*, Acta Appl. Math. **76** (2003), 137–193.
- P. J. Olver, *Applications of Lie Groups to Differential Equations*, 2nd ed., Springer GTM 107, 1993.
- P. J. Olver, *Equivalence, Invariants, and Symmetry*, Cambridge University Press, 1995.
