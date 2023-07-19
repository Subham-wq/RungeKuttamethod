def runge_kutta(f, y0, t0, t_end, h):
    """
    Solve the first-order ordinary differential equation (ODE) using the Runge-Kutta method.

    Parameters:
        f: Function representing the ODE dy/dt = f(t, y), where t is the independent variable and y is the dependent variable.
        y0: Initial value of the dependent variable (y) at t = t0.
        t0: Initial value of the independent variable (t).
        t_end: End value of the independent variable (t) for the solution.
        h: Step size for the Runge-Kutta method.

    Returns:
        Two lists: t_values and y_values, representing the time steps and corresponding values of the dependent variable y.
    """
    t_values = [t0]
    y_values = [y0]

    while t0 < t_end:
        k1 = h * f(t0, y0)
        k2 = h * f(t0 + h/2, y0 + k1/2)
        k3 = h * f(t0 + h/2, y0 + k2/2)
        k4 = h * f(t0 + h, y0 + k3)

        y0 = y0 + (k1 + 2*k2 + 2*k3 + k4) / 6
        t0 = t0 + h

        t_values.append(t0)
        y_values.append(y0)

    return t_values, y_values

# Example usage:
def example_ode(t, y):
    # Define the ODE dy/dt = f(t, y) here
    return t * y  # Example: t * y

# Initial conditions
initial_t = 0
initial_y = 1

# End value for t and step size
end_t = 5
step_size = 0.1

# Solving the ODE using the Runge-Kutta method
time_steps, y_values = runge_kutta(example_ode, initial_y, initial_t, end_t, step_size)

# Printing the results
for t, y in zip(time_steps, y_values):
    print(f"t = {t:.2f}, y = {y:.5f}")
