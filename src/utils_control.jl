#= Adaptive-DDTO -- Utility functions for generalized control capabilities.
Author: Samuel Buckner (UW-ACL)
=#

function optimal_controller(t::CReal, x::CVector, sol::Solution)::CVector
    """
    Compute the optimal control input (quadcopter thrust vector) at time t.

    Args:
        t (CReal): the time at which to compute the optimal control.
        x (CVector): the current state (not currently used, no feedback control).
        sol (Solution): the optimized solution to track.
    
    Returns:
        U (CVector): the optimal input.
    """

    # Get current optimal thrust (ZOH interpolation)
    i = findlast(τ->τ<=t,sol.t)
    if typeof(i)==Nothing || i>=size(sol.T,2)
        T = sol.T[:,end]
    else
        T = sol.T[:,i]
    end
   
    # Create the input vector for the state-space dynamics
    U = CVector(vcat(T,norm(T,2)))

    return U
end

function c2d_zoh(lander::Lander, Δt::CReal)::Tuple{CMatrix, CMatrix, CVector}
    """
    Discretize the continuous-time dynamics of the lander using zeroth-order hold (ZOH)
    for a discrete-time state-space system of the form:
        x_{k+1} = A*x_k + B*u_k + p

    Args:
        lander (Lander): the lander object.
        Δt (CReal): the discretization time step.

    Returns:
        A (CMatrix): the discrete-time state transition matrix.
        B (CMatrix): the discrete-time input matrix.
        p (CVector): the discrete-time additive affine vector.
    """

    A_c,B_c,p_c,n,m = lander.A_c,lander.B_c,lander.p_c,lander.n,lander.m

    _M = exp(CMatrix([
        A_c B_c p_c;
        zeros(m+1,n+m+1)
    ])*Δt)

    A = _M[1:n,1:n]
    B = _M[1:n,n+1:n+m]
    p = _M[1:n,n+m+1]
    return (A,B,p)
end

function rk4_step(x_cur::CVector, f::Function, t_cur::CReal, Δt::CReal)::CVector
    """
    Integrate a system of ordinary differential equations (ODE)
    one time-step forward using RK4 (updates x_cur in place).

    Args:
        x_cur (CVector): the current state.
        f (Function): the function defining the ODE, dx/dt=f(t,x).
        t_cur (CReal): the current time (in DDTO solution).
        Δt (CReal): the integration time step.

    Returns:
        x_new (CVector): the new state.
    """

    # ..:: Integrate one time-step forward ::..
    y = x_cur
    h = Δt
    t_ = t_cur
    k1 = f(t_,y)
    k2 = f(t_+h/2,y+h*k1/2)
    k3 = f(t_+h/2,y+h*k2/2)
    k4 = f(t_+h,y+h*k3)
    x_cur = y+h/6*(k1+2*k2+2*k3+k4)

    return x_cur
end

function simulate_cont(lander::Lander, sol::Solution, control::Function)::Solution
    """
    Simulate the quadcopter dynamics using a predefined control input
    trajectory in continuous time with RK4 integration.

    Args:
        lander (Lander): the vehicle object.
        sol (Solution): Optimization solution object.
        control (Function): the control trajectory, expressed as a function
                            u=control(t,x,r) where u is the control input, t
                            is the time, x is the state, and r is the lander
                            object.

    Returns:
        sim (Solution): the simulation output as a Solution object.
    """

    # ..:: Simulate ::..

    dynamics = (t,x) -> lander.A_c*x + lander.B_c*control(t,x,sol) + lander.p_c
    x0 = CVector(vcat(sol.r[:,1], sol.v[:,1]))
    Δt = 1e-2
    tf = sol.t[end]
    t = 0:Δt:tf
    N = length(t)
    X = zeros(lander.n,N)
    X[:,1] = x0
    for i = 1:N-1
        X[:,i+1] = rk4_step(X[:,i],dynamics,t[i],Δt)
    end
    U = CMatrix(hcat([control(t[n],X[:,n],sol) for n = 1:N]...))

    # ..:: Save solution ::..

    r = X[1:3,:]
    v = X[4:6,:]
    T = U[1:3,:]
    Γ = U[4,:]

    T_nrm = CVector([norm(T[:,i],2) for i=1:N])
    γ = CVector([acos(dot(T[:,k],e_z)/norm(T[:,k],2)) for k=1:N])

    sim = Solution(t,r,v,T,Γ,0.0,T_nrm,γ)

    return sim
end