import jax
import jax.numpy as jnp
from jax import grad, hessian, jit
from scipy.optimize import minimize

jax.config.update("jax_enable_x64", True)  # use float64, important for numerics

def f(x):
  return jnp.sum(x**2) + (x[0] - 1) * x[1]

x0 = jnp.array([1.0, 2.0, 3.0])

# Value, gradient, Hessian
f(x0)
grad(f)(x0)
hessian(f)(x0)

# Optimize — scipy.optimize.minimize works well with JAX via .numpy() conversion
grad_f    = jit(grad(f))
hessian_f = jit(hessian(f))

res = minimize(
  fun=lambda x: float(f(jnp.array(x))),
  x0=[1.0, 2.0, 0.0],
  jac=lambda x: jax.device_get(grad_f(jnp.array(x))),
  method="BFGS",
  options={"gtol": 1e-8}
)

res.x
hessian_f(jnp.array(res.x))