### New Parametrization

Consider general diffusion equation, $J=A_\alpha P - \partial_\alpha (DP)$ where we assume $D$ is just homoegeous (a scalar constant). Then the equilibrium solution satisfies, $A_\alpha P = \partial_\alpha D P + D\partial_\alpha P$. If we assume a potential form $P=e^{-\phi}$ then, we have 
$$
\begin{align*}
A_\alpha P &= \partial_\alpha D P - D \partial_\alpha \phi P \\
\partial_\alpha \phi &= \frac{\partial_\alpha D - A_\alpha}{D}
\end{align*}
$$
So the goal is to find a parametrization whereby this produces a potential where in the wells it's easily approximated by second order taylor expansion (meaning $\partial_\alpha \phi$ should be near linear around those points). One proxy may be to get a simple expression for $\partial_\alpha \phi$. 

Note that $D = \frac{1}{2}B^2$ where $A, B$ are determined by the Ito SDE expression (assume $\alpha=\beta$),
$$
\begin{align*}
\frac{d\alpha}{d\tau} &= b -\alpha + \lambda \alpha - g^2 \alpha^3 + \sqrt{\lambda - g^2 \alpha^2} W
\end{align*}
$$
Now assume a general parametrization $u(\alpha), \alpha(u)$, we have,
$$
\begin{align*}
\frac{du}{d\tau} &= u' \left(b -\alpha + \lambda \alpha - g^2 \alpha^3 \right) + u'\sqrt{\lambda - g^2 \alpha^2} W
\end{align*}
$$
basically multiplying $A$, $B$ by $u'$. And $D= \frac{1}{2}u'^2 (\lambda -g^2 \alpha^2)$. Now we find $\partial_u \phi$, 
$$
\begin{align*}
\partial_u \phi &= \frac{\partial_u D - A_u}{D} \\
&= \frac{u' \partial_u u'(\lambda - g^2 \alpha^2)  - u'^2g^2\alpha \alpha' - u'(b-\alpha + \lambda\alpha - g^2\alpha^3)}{\frac{1}{2}u'^2 (\lambda -g^2\alpha^2)} \\
&= 2\frac{\partial_u u'(\lambda - g^2 \alpha^2)  - g^2\alpha - (b-\alpha + \lambda\alpha - g^2\alpha^3)}{u' (\lambda -g^2\alpha^2)} \\
\end{align*}
$$
now suppose for convenience that $u' = (\lambda - g^2 \alpha^2)^{-1}$. Then $u = \frac{1}{g\sqrt{\lambda}}\tanh^{-1}(\frac{g\alpha}{\sqrt{\lambda}})$ IT DOESN'T PRODUCE ANYTHING NICE.


### Dot Product

$$
\begin{align*}
V_s(u, v) &= \lambda (\cos u - \cos v) - 2\ln\left(\cos u + \cos v\right) \\
V_a(u, v) &= \frac{2gb}{\sqrt{\lambda}} \ln\left(\cos u + \cos v\right) - \frac{2gb}{\sqrt{\lambda}} \ln\left( 2 + 2\sin\left(\frac{u+v}{2}\right) + 2\sin\left(\frac{u-v}{2}\right)  + \cos v - \cos u \right)
\end{align*}
$$

Now let us attempt approximation via dot product. Consider a one dimensional approximation $c_0 + c_1 u + \frac{1}{2}c_2 u^2$. 

Note $\int \ln(\cos \frac{u}{2}) = \int_0^{u_A} \ln(\sin(\frac{\pi}{2} - \frac{u}{2})) = -\int_{\pi-u_A}^{\pi} \ln(2\sin \frac{x}{2})+ u_A\ln 2 = \text{Cl}(\pi) -\text{Cl}(\pi-u_A)+ u_A \ln 2 = -\text{Cl}(\pi-u_A)+u_A \ln 2$. For now let's denote this $K$. We have the following
$$
\begin{align*}
&\int_0^{u_A} \lambda(\cos u - 1) c_0 + \int_0^{u_A} \lambda(\cos u - 1) c_1 u + \int_0^{u_A} \lambda(\cos u - 1) \frac{c_2}{2} u^2 \\
-&\int_0^{u_A} 2\ln\left(\cos u + 1\right) c_0 - \int_0^{u_A} 2\ln\left(\cos u + 1\right) c_1 u - \int_0^{u_A} 2\ln\left(\cos u + 1\right)\frac{c_2}{2} u^2 \\
=& \lambda c_0 (\sin u_A - u_A) + \lambda c_1 u_A \sin u_A + \lambda c_1 (\cos u_A - 1) - \lambda c_1 \frac{u_A^2}{2} + \lambda \frac{c_2}{2} \sin u_A u_A^2 \\
&- \lambda c_2 u_A (\cos u_A - 1) - \lambda \frac{c_2}{6} u_A^3  - 2 c_0 u_A \ln 2  - 4c_0 K -\cdots
\end{align*}
$$
This is too messy, let's start over. Let $\int \lambda(\cos u - 1) = M$. Let $\int 2\ln(\cos u+1) = N$. Then the integration by part gives 
$$
\begin{align*}
M c_0 + M c_1 u_A - \int  M c_1 + M \frac{c_2}{2}u_A^2 - \int M c_2 u_A + \int M c_2
\end{align*}
$$ 
and same for $N$. $M= \lambda (\sin u_A - u_A)$, and $\int M = -\lambda (\cos u_A + \frac{u_A^2}{2})$. So in total we got for the $M$ terms,
$$
\begin{align*}
\lambda \sin u_A ( c_0 + c_1 u_A + \frac{c_2}{2}u_A^2) + \lambda \cos u_A (c_1 +c_2 u_A - c_2)
\end{align*}
$$
