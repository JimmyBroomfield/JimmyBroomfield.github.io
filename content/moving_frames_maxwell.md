# Moving Frames in Field Theory

### From the Klein–Gordon Equation to Maxwell's Equations via Fels–Olver and Olver–Pohjanpelto Invariant Theory

---

## Abstract

The equations of relativistic field theory are normally obtained by postulating Lorentz invariance and writing down "all the scalars one can build" from a field and its derivatives. This paper asks what happens if that step is not postulated but *derived*, using the theory of moving frames. We treat three cases of increasing structural richness. First, a free scalar field in $1{+}1$ dimensions, where the base group $ISO(1,1)$ acts trivially on the fiber; normalizing the group's action on the jet bundle manufactures the invariant $\kappa=\sqrt{\eta^{\mu\nu}\partial_\mu u\,\partial_\nu u}$ directly, and the Klein–Gordon equation follows from an invariant Euler operator built from it. Second, a massive vector field, where the fiber now carries a nontrivial representation of the Lorentz group; the same normalization produces both the field strength $F_{\mu\nu}$ and the mass invariant $A_\mu A^\mu$, and the invariant Euler–Lagrange equations reproduce the Proca equations. Third — and this is where the finite-dimensional theory runs out — we show that gauge invariance is not a symmetry of any enlargement of the Poincaré group, but a genuinely different object: an infinite-dimensional Lie pseudogroup. Applying the Olver–Pohjanpelto extension of moving frames to the gauge pseudogroup, we show that the field strength $F_{\mu\nu}$ arises as the unique first-order invariant left over after normalizing the potential's jets — forced into existence by the simple fact that the pseudogroup's parameter is a scalar function, so its second-order jets are symmetric while the field's are not. The Bianchi identity appears as an automatic syzygy of this construction rather than a separately imposed law. We close with Maxwell's equations and an honest account of which parts of this derivation are fully rigorous and which remain to be carried out within the formal invariant variational bicomplex for pseudogroups.

---

## 1. Introduction

Lie's original observation about invariant variational problems was that a $G$-invariant Lagrangian can always be rewritten purely in terms of the differential invariants of $G$, and that the resulting Euler–Lagrange equations inherit the same invariance and can themselves be expressed in invariant language. Turning this into an algorithm — rather than an existence theorem — is the content of the *moving frame method*, developed by Fels and Olver for finite-dimensional Lie group actions, and extended to a full calculus of variations (the *invariant variational bicomplex*) by Kogan and Olver. A separate, later development by Olver and Pohjanpelto extends the moving frame construction itself to **Lie pseudogroups** — group-like symmetries parametrized by arbitrary functions rather than finitely many constants, of which gauge transformations are the paradigm example in physics.

This method is most often illustrated on curves: a single independent variable, a Euclidean, similarity, or affine group acting on the ambient space, with curvature and torsion as the resulting invariants. Here we illustrate the same ideas on **field theories**, where the group typically acts on the independent (spacetime) variables while the dependent variable (the field) may or may not transform nontrivially. This turns out to be a natural laboratory for the method for a specific reason: whether or not the field's fiber carries a nontrivial group representation is exactly the difference between a "curve with no frame" and a "framed curve," and the amount of machinery the moving frame method needs to invoke tracks that distinction precisely. We use this to build up, case by case, to the point where the finite-dimensional theory runs out and a pseudogroup is required — which happens to be exactly the point where electromagnetism's gauge symmetry lives.

None of the equations derived here are new. The purpose of this paper is methodological: to show that a systematic, algorithmic procedure — normalize the group's action on the jet bundle, read off the invariants, differentiate invariantly, reconstruct the Euler–Lagrange equations — reproduces the field content that physicists normally arrive at by informed guesswork, and does so in a way that makes explicit exactly *why* the invariants take the form they do.

---

## 2. Preliminaries: Jets, Moving Frames, and Invariantization

Let $X$ be the space of independent (spacetime) variables and $U$ the space of dependent variables (field values), with jet bundle $J^n(X,U)$ carrying coordinates $(x^i, u^\alpha, u^\alpha_J)$ for multi-indices $J$. A Lie group $G$ acting on $X\times U$ prolongs to an action on $J^n(X,U)$ in the usual way, and its infinitesimal generators $v=\xi^i(x,u)\partial_{x^i}+\varphi^\alpha(x,u)\partial_{u^\alpha}$ prolong via

$$
\varphi^\alpha_i = D_i\Big(\varphi^\alpha - \sum_j \xi^j u^\alpha_j\Big) + \sum_j \xi^j u^\alpha_{ij},
\tag{2.1}
$$

with $D_i$ the total derivative operator.

If $G$ acts freely and regularly near a point of $J^n(X,U)$, a **moving frame** is a locally defined, right-equivariant map $\rho: J^n(X,U)\to G$. It is constructed by choosing a cross-section $K$ to the group orbits — a set of normalization equations, one per group parameter — and solving them for the group parameters in terms of the jet coordinates. Given $\rho$, the **invariantization** of any differential function $F$ is

$$
\iota(F)(x,u^{(n)}) = F\big(\rho(x,u^{(n)})\cdot(x,u^{(n)})\big),
\tag{2.2}
$$

the unique $G$-invariant function agreeing with $F$ on the cross-section. Invariantizing the jet coordinates themselves produces two kinds of invariants: **phantom invariants**, which are the normalization constants themselves (hence trivial), and **normalized invariants**, which carry the genuine content of the theory. Because invariantization does not commute with differentiation, one tracks the discrepancy via **recurrence relations**, and the invariant analogue of ordinary differentiation is carried out by an **invariant differential operator** built from the frame.

This machinery applies equally whether $G$ acts trivially or nontrivially on the fiber $U$ — but, as the three cases below show, how much of it is actually *needed* depends heavily on which situation you're in.

---

## 3. Case I: The Scalar Field and the Klein–Gordon Equation

### 3.1 Setup

Take $X=\mathbb{R}^{1,1}$, coordinates $x^1=t,\,x^2=x$, metric $\eta=\mathrm{diag}(+1,-1)$; $U=\mathbb{R}$, coordinate $u$. The group is $G=ISO(1,1)$, acting by

$$
T = t\cosh\psi+x\sinh\psi+a,\qquad X = t\sinh\psi+x\cosh\psi+b,\qquad U=u.
\tag{3.1}
$$

The condition $U=u$ says the field is a **Lorentz scalar**: the group moves the point at which $u$ is evaluated, but the value itself does not change. This is the field-theoretic analogue of a curve with no attached frame. Infinitesimal generators:

$$
v_1=\partial_t,\qquad v_2=\partial_x,\qquad v_3 = x\partial_t+t\partial_x,
\tag{3.2}
$$

with no $\partial_u$-component. Prolonging $v_3$ via (2.1) with $\xi^t=x,\ \xi^x=t,\ \varphi^u=0$:

$$
\varphi^t = -u_x,\qquad \varphi^x=-u_t,
\tag{3.3}
$$

reproducing, from first principles, the infinitesimal form of $u_t\to u_t\cosh\psi-u_x\sinh\psi$.

### 3.2 Moving frame and the fundamental invariant

Translations normalize $T=X=0$. The boost has no fixed points left to use, so it must be normalized using derivative data: impose

$$
U_X = -u_t\sinh\psi+u_x\cosh\psi = 0 \quad\Longrightarrow\quad \tanh\psi=\frac{u_x}{u_t},
\tag{3.4}
$$

valid where the gradient of $u$ is timelike. This fixes $\cosh\psi=u_t/\sqrt{u_t^2-u_x^2}$, $\sinh\psi=u_x/\sqrt{u_t^2-u_x^2}$, and defines the moving frame $\rho(t,x,u^{(1)})=(a,b,\psi)$. Invariantizing $U_T$:

$$
\kappa := \iota(U_T) = u_t\cosh\psi-u_x\sinh\psi = \sqrt{u_t^2-u_x^2}.
\tag{3.5}
$$

This is the theory's fundamental first-order invariant. It is ordinarily introduced by writing $\eta^{\mu\nu}\partial_\mu u\,\partial_\nu u$ and asserting that it is a scalar; here it appears as the unique quantity that survives after the group's three parameters have exhausted their normalizing power over the two independent variables and one first derivative.

### 3.3 Invariant differential operators

Dual to the invariant coframe are two invariant total derivative operators,

$$
D_\parallel = \frac{u_t\partial_t-u_x\partial_x}{\sqrt{u_t^2-u_x^2}},\qquad
D_\perp = \frac{-u_x\partial_t+u_t\partial_x}{\sqrt{u_t^2-u_x^2}},
\tag{3.6}
$$

an invariant orthonormal frame built entirely from $u_t,u_x$. Applying $D_\parallel$ to $\kappa$ generates a genuine second-order invariant,

$$
D_\parallel\kappa = \frac{u_t^2u_{tt}-2u_tu_xu_{tx}+u_x^2u_{xx}}{u_t^2-u_x^2},
\tag{3.7}
$$

which is *not* the d'Alembertian $\Box u=u_{tt}-u_{xx}$ — both are legitimate invariants generated by the same algebra, but they measure different things, exactly as a curve's curvature and its arc-length derivative of curvature are related but distinct invariants.

### 3.4 The invariant Euler–Lagrange equation

Any $G$-invariant Lagrangian is a function $\widetilde L(\kappa,u)$ of the invariants; smoothness under the full group (including reflections) forces dependence through $\kappa^2=u_t^2-u_x^2$. Because $u$ carries no fiber index, the invariant Euler operator collapses to the ordinary chain rule through $\kappa^2$:

$$
E(\widetilde L) = \frac{\partial\widetilde L}{\partial u} - 2\partial_t\!\left(\frac{\partial\widetilde L}{\partial(\kappa^2)}u_t\right) + 2\partial_x\!\left(\frac{\partial\widetilde L}{\partial(\kappa^2)}u_x\right)=0.
\tag{3.8}
$$

Choosing $\widetilde L=\tfrac12\kappa^2-\tfrac12m^2u^2$ gives

$$
\boxed{\ u_{tt}-u_{xx}+m^2u=0\ }
\tag{3.9}
$$

the Klein–Gordon equation, with $\kappa$, the invariant operators, and the field equation all produced mechanically from the single normalization (3.4).

---

## 4. Case II: The Vector Field and the Proca Equation

### 4.1 A genuinely nontrivial fiber

Now let $U=\mathbb{R}^2$, coordinates $A_t,A_x$, transforming as a Lorentz covector — *not* trivially, unlike the scalar case:

$$
A_T = A_t\cosh\psi-A_x\sinh\psi,\qquad A_X=-A_t\sinh\psi+A_x\cosh\psi.
\tag{4.1}
$$

The boost generator now has a nonzero fiber component:

$$
v_3 = x\partial_t+t\partial_x - A_x\partial_{A_t}-A_t\partial_{A_x}.
\tag{4.2}
$$

This is the field-theoretic analogue of a framed curve's frame vector rotating under the group — the essential feature absent from Case I.

Differentiating (4.1), the antisymmetric combination

$$
F_{tx} := A_{x,t}-A_{t,x}
\tag{4.3}
$$

is already invariant under the boost, since two-dimensional antisymmetric tensors are automatically scalars (the Levi-Civita symbol has unit determinant under proper Lorentz transformations). This is a low-dimensional peculiarity, flagged here rather than treated as a general feature: in higher dimensions $F_{\mu\nu}$ remains a genuine tensor, not a scalar, and only $F_{\mu\nu}F^{\mu\nu}$ (and its dual) survive as Poincaré scalars.

### 4.2 Moving frame

Translations again fix $T=X=0$. The boost is normalized using **zeroth-order** field data:

$$
\text{cross-section: } A_X=0 \quad\Longrightarrow\quad \tanh\psi=\frac{A_x}{A_t},\qquad
\kappa:=\sqrt{A_t^2-A_x^2}.
\tag{4.4}
$$

This is qualitatively different from (3.4): there, the frame was fixed by a *derivative* condition; here, the fiber itself supplies enough structure to fix the frame at zeroth order. $\kappa$ is the theory's mass-type invariant, and $F_{tx}$ (already invariant) is its curvature-type invariant.

### 4.3 Lie derivatives of the contact forms

With $\theta^{A_t}=dA_t-A_{t,t}dt-A_{t,x}dx$ and similarly for $\theta^{A_x}$, translations annihilate both, and the boost gives

$$
v_3(\theta^{A_t}) = -\theta^{A_x},\qquad v_3(\theta^{A_x})=-\theta^{A_t}.
\tag{4.5}
$$

The contact forms mix under the group exactly as the field components do — a genuine $2\times2$ matrix action, the direct analogue of the matrix of Lie derivatives that organizes the invariant Eulerian operator for framed curves in the classical theory. This is the structure entirely absent from Case I.

### 4.4 Field equations

Any invariant Lagrangian is a function of $\kappa$ and $F_{tx}$ (and their invariant derivatives). Taking the simplest quadratic choice — the Proca Lagrangian —

$$
\widetilde L(\kappa,F_{tx}) = \tfrac12F_{tx}^2+\tfrac12m^2\kappa^2,
\tag{4.6}
$$

and varying (the chain rule through $\kappa,F_{tx}$ reproduces the invariant Euler operator built from (4.5)):

$$
m^2A_t+\partial_xF_{tx}=0,\qquad -m^2A_x-\partial_tF_{tx}=0
\quad\Longleftrightarrow\quad
\boxed{\ \partial_\mu F^{\mu\nu}+m^2A^\nu=0\ }
\tag{4.7}
$$

the Proca equation. Both terms in (4.6) — normally written down as "the" quadratic Lorentz invariants built from $A_\mu$ — emerge here as exactly the two invariants the moving frame manufactures from (4.3)–(4.4), with no separate argument required for completeness at this order.

### 4.5 Where the finite-dimensional theory stops

Set $m=0$. The Lagrangian $\widetilde L=\tfrac12F_{tx}^2$ becomes invariant under $A_\mu\to A_\mu+\partial_\mu\lambda$ for an *arbitrary* function $\lambda(t,x)$ — an infinite-dimensional symmetry the finite-dimensional moving frame construction has no way to detect, since it was built entirely from the three-parameter group $ISO(1,1)$. This is not a defect that a larger *finite-dimensional* group could fix: gauge invariance is not a symmetry of spacetime transformations at all, but a separate redundancy living entirely in the fiber. It calls for a different category of group action altogether.

---

## 5. Case III: Gauge Invariance as a Lie Pseudogroup

### 5.1 Why a pseudogroup, and not a bigger Lie group

The gauge transformations

$$
\mathcal{G}:\quad A_\mu(x)\ \longrightarrow\ A_\mu(x)+\partial_\mu\lambda(x)
\tag{5.1}
$$

are parametrized by an arbitrary smooth function $\lambda$, not by finitely many constants. No enlargement of the Poincaré group — adding dilations, adding the full conformal group, anything built from finitely many parameters — will ever produce (5.1), because the issue is not which spacetime transformations are allowed but that the *field's own redundancy* is described by infinitely many independent parameters, one function's worth. This is precisely the setting for which Olver and Pohjanpelto extended the moving frame method: **Lie pseudogroups**, whose local coordinate expressions are parametrized by arbitrary functions satisfying a system of differential equations (here, trivially, none — $\lambda$ is unconstrained), rather than by finitely many group coordinates.

We work in general spacetime dimension $D$ for this section (physically $D=4$), since $1{+}1$ dimensional electromagnetism is dynamically trivial (Section 4 already showed $F_{tx}=\mathrm{const}$ in the massless limit) and would obscure rather than illustrate the construction.

### 5.2 The pseudogroup jet groupoid

The action of $\mathcal G$ on jets is encoded by adjoining the jets of $\lambda$ itself — $\lambda,\lambda_\mu,\lambda_{\mu\nu},\lambda_{\mu\nu\rho},\dots$ — as coordinates on a groupoid over $X$, and letting them act on the field's jets by simple prolongation of (5.1):

$$
A_\nu\to A_\nu+\lambda_\nu,\qquad A_{\nu,\mu}\to A_{\nu,\mu}+\lambda_{\nu\mu},\qquad A_{\nu,\mu\rho}\to A_{\nu,\mu\rho}+\lambda_{\nu\mu\rho},\ \dots
\tag{5.2}
$$

The structural fact driving everything that follows: because $\lambda$ is a genuine scalar function, its jets satisfy the ordinary symmetry of mixed partial derivatives —

$$
\lambda_{\nu\mu}=\lambda_{\mu\nu},\qquad \lambda_{\nu\mu\rho}=\lambda_{(\nu\mu\rho)},\ \dots
\tag{5.3}
$$

— while $A_{\nu,\mu}$ carries **no** such constraint. This mismatch between the symmetry type of the pseudogroup's parameter and the symmetry type of the field is the entire mechanism by which the moving frame produces $F_{\mu\nu}$.

### 5.3 Normalization, order by order

**Order 0.** At a point, $\lambda_\nu$ has exactly $D$ independent components, matching $A_\nu$. Normalize:

$$
\text{cross-section: } A_\nu=0 \quad\Longrightarrow\quad \lambda_\nu:=-A_\nu.
\tag{5.4}
$$

This exhausts the first-order pseudogroup freedom entirely; $\iota(A_\nu)=0$ is the moving-frame statement that the value of the potential at a point carries no invariant content — it can always be gauged away locally.

**Order 1.** With $\lambda_\nu$ fixed, the residual freedom lives in $\lambda_{\nu\mu}$, which by (5.3) is *symmetric*. Decompose

$$
A_{\nu,\mu} = A_{(\nu,\mu)}+A_{[\nu,\mu]}.
$$

The symmetric part is normalized away:

$$
\text{cross-section: } A_{(\nu,\mu)}=0 \quad\Longrightarrow\quad \lambda_{\nu\mu}:=-A_{(\nu\mu)}.
\tag{5.5}
$$

The antisymmetric part has no corresponding freedom to absorb it, since $\lambda_{\nu\mu}-\lambda_{\mu\nu}\equiv 0$ identically:

$$
\iota\big(A_{[\nu,\mu]}\big) = A_{\nu,\mu}-A_{\mu,\nu} = F_{\mu\nu}.
\tag{5.6}
$$

**The field strength is the fundamental differential invariant of the gauge pseudogroup** — forced into existence not by physical insight but by the elementary fact that a symmetric quantity cannot cancel an antisymmetric one. This is the pseudogroup counterpart of $\kappa$ falling out of the boost normalization in Section 3: an invariant produced by counting exactly what freedom remains after normalization.

### 5.4 The Bianchi identity as a syzygy

Proceeding to order 2, the fully symmetric part of $A_{\nu,\mu\rho}$ is normalized away using the fully symmetric $\lambda_{\nu\mu\rho}$, leaving $\partial_\rho F_{\mu\nu}$ as the next tier of invariants — the pseudogroup analogue of $D_\parallel\kappa$ in Section 3. But because $F_{\mu\nu}$ was itself built from a single potential $A_\mu$, its derivatives satisfy an automatic algebraic relation:

$$
\partial_{[\lambda}F_{\mu\nu]}=0.
\tag{5.7}
$$

This is a **syzygy** among the invariants generated by the pseudogroup — an unavoidable relation of exactly the kind that appears throughout ordinary moving frame theory relating normalized and curvature invariants — rather than an independently postulated physical law. The moving frame construction here is telling us, as a byproduct of normalization, that not every antisymmetric tensor field can be a field strength; only those satisfying (5.7) are consistent with having arisen from a single-valued potential.

5.5 Recombining with the Poincaré group
F mu nu still carries free spacetime indices; it is invariant under G but not yet organized into scalars of ISO(1,D minus 1). Applying the finite-dimensional Fels–Olver construction of Section 4 to F mu nu itself — rather than to A mu, which the pseudogroup has already annihilated — produces the familiar Lorentz invariants (for D=4):
F_{\mu\nu}F^{\mu\nu},\qquad F_{\mu\nu}\widetilde F^{\mu\nu}.
\tag{5.8}
The overall structure is a two-stage normalization: the pseudogroup G first removes the potential's gauge redundancy and manufactures F mu nu; the finite-dimensional group ISO(1,D minus 1) then organizes what remains into spacetime scalars, exactly as it did for the scalar and vector fields of Sections 3–4. This sequencing is the precise sense in which gauge invariance and Lorentz invariance are different symmetries acting on different data, resolving the question of whether "a bigger group" could have captured gauge invariance from the start: it could not, because the two symmetries act on structurally different objects (a function's worth of redundancy versus a finite frame rotation) and must be normalized in sequence, not merged into one larger finite-dimensional group.
5.6 Maxwell's equations
With widetilde L= minus fraction 14F mu nu F to the mu nu a function of the invariants alone, and A mu normalized to zero identically by the pseudogroup, the field equation is obtained by varying through F mu nu = partial [ mu A nu ]:
\boxed{\ \partial_\mu F^{\mu\nu}=0\ }
\tag{5.9}
together with the automatic identity (5.7). The traditional split of Maxwell's equations into "the dynamical half" and "the identity half" is, from this point of view, a direct readout of the moving frame construction: (5.9) comes from extremizing an invariant built from the leftover antisymmetric jet data, while (5.7) is forced by the mere fact that this data came from a potential in the first place.
5.7 Scope and rigor of this derivation
The argument in Sections 5.3–5.4 — that F mu nu is exactly the antisymmetric residue of a normalization forced by the symmetric jet structure of a scalar gauge parameter, and that the Bianchi identity is an automatic syzygy of that construction — is elementary and, we believe, fully rigorous as stated. What has not been carried out here is the full Olver–Pohjanpelto formalism for this pseudogroup: their Maurer–Cartan forms and structure equations for infinite-dimensional group actions, and the pseudogroup extension of the invariant Eulerian and Hamiltonian operators (A to the *,B to the *,W in the finite-dimensional theory of Sections 3–4) that would derive equation (5.9) directly from an invariant variational bicomplex for G, rather than by falling back, as we did, to ordinary variation of A mu once F mu nu had been identified. We are confident (5.9) is what that formalism would produce — it is the correct physics, and the invariant algebra generated by G and ISO(1,D minus 1) together contains no other candidate quadratic invariant to build a Lagrangian from — but the derivation via the formal apparatus, rather than around it, remains open.
6. Discussion
The three cases traced out here form a natural progression in how much of the moving frame apparatus is actually engaged:
Fiber action of base group

Extra structure needed
Result
Scalar field (§3)
Trivial
None — invariant Euler operator collapses to ordinary chain rule
Klein–Gordon
Vector field (§4)
Nontrivial (rotates fiber)
2 times 2 matrix of Lie derivatives on contact forms
Proca
Gauge potential (§5)
Nontrivial and redundant (infinite-dimensional)
Pseudogroup jet groupoid; symmetric/antisymmetric jet decomposition
Maxwell + Bianchi
Nothing in this progression is a new physical result. What it offers is a uniform mechanism: in every case, the field equation's building blocks are exactly the invariants left over after the relevant group has used up all the freedom it has, evaluated order by order on the jet bundle. The scalar and vector cases show this mechanism working cleanly within finite-dimensional Fels–Olver theory. The gauge case shows where that theory necessarily hands off to the pseudogroup extension of Olver and Pohjanpelto — and shows, concretely, that the handoff produces exactly the objects (field strength, Bianchi identity) that gauge theory is built from, via nothing more than a symmetric-versus-antisymmetric jet-counting argument.
The natural next steps, in decreasing order of how well-scoped they are: (i) complete the Olver–Pohjanpelto invariant variational bicomplex for the U(1) pseudogroup explicitly, to obtain (5.9) from the formal apparatus rather than by shortcut; (ii) extend Section 5 to non-abelian gauge pseudogroups, where the jet groupoid's structure equations become considerably richer and might be expected to manufacture the Yang–Mills field strength and its non-abelian Bianchi identity by the same mechanism; (iii) revisit the massless limit of Section 4 in D=4, where F mu nu F to the mu nu genuinely is dynamical, to check that the two-stage normalization of Section 5.5 reproduces the free Maxwell action with no further input.

# References
- M. Fels and P. J. Olver, Moving coframes II: Regularization and theoretical foundations, Acta Appl. Math. 55 (1999), 127–208.
- I. Kogan and P. J. Olver, The invariant variational bicomplex, Contemp. Math. 285 (2001), 131–144.
- I. Kogan and P. J. Olver, Invariant Euler–Lagrange equations and the invariant variational bicomplex, Acta Appl. Math. 76 (2003), 137–193.
- P. J. Olver and J. Pohjanpelto, Moving frames for Lie pseudo-groups, Canad. J. Math. 60 (2008), 1336–1386.
- P. J. Olver and J. Pohjanpelto, Differential invariant algebras of Lie pseudo-groups, Adv. Math. 222 (2009), 1746–1792.
- P. J. Olver, Applications of Lie Groups to Differential Equations, 2nd ed., Springer GTM 107, 1993.
P. J. Olver, Equivalence, Invariants, and Symmetry, Cambridge University Press, 1995.