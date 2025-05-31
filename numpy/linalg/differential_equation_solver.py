from sympy import symbols, integrate, sin, cos, exp,simplify,solve,Eq,sqrt
import numpy as np

x = symbols('x')
y = symbols('y')
c = symbols('c')
f_y = symbols('f(y)')
g_x = symbols('g(x)')
q_x = symbols('q(x)')
v = symbols('v')
f_v = symbols('f(v)')
mu = symbols('mu')
n = symbols('n')
a = symbols('a')
b = symbols('b')
c = symbols('c')
c1 = symbols('c1')
c2 = symbols('c2')
i = symbols('i')
def solve_separable_eqn(f_y ,g_x , x0=None , y0=None ,solution_type='general'):
    if solution_type=='general':
        f1 = integrate((1/f_y) , y)
        f2 = integrate(g_x , x)
        print(f"{simplify(f1)} = {simplify(f2)} +{c}") ##arbitary constant
        return f1,f2
    elif solution_type=='particular':
        f1 = integrate((1/f_y) , (y,y0))
        f2 = integrate(g_x , (x,x0))
        print(f"{f1} = {f2}")
        return f1 , f2
    

def solve_homogenous_separable_eqn(f_v , x0=None , y0=None , solution_type='general'):
    if solution_type=='general':
        f1 = integrate((1/x) , x)
        f2 = integrate((1/(f_v -v)) ,v)
        print(f"{simplify(f1)} = {simplify(f2)} + {c}")
        return f2,f1
    
    elif solution_type=='particular':
        f1 = integrate((1/x) ,(x,x0))
        f2 = integrate((1/(f_v - v)) , (y , (y0/x0)))
        print(f"{simplify(f2)} = {simplify(f1)}")
        return f2,f1
    

def solve_linear_ode(g_x ,y,q_x ,x0=None , y0 = None , solution_type = 'general'):
    if solution_type=='general':
        mu = exp(integrate(g_x,x))
        sol = (integrate((mu * q_x)) + c)/mu
        sol = simplify(sol)
        print("y =" , sol)
        return sol
    elif solution_type=='particular':
        mu = exp(integrate(g_x, x))  # integrating factor
        general_sol = (integrate(mu * q_x, x) + c) / mu
        particular_eq = general_sol.subs(x, x0) - y0
        c_val = solve(particular_eq, c)[0]
        final_sol = simplify(general_sol.subs(c, c_val))
        return final_sol
    

def solve_bernoulli(g_x,y,q_x,n,x0=None,y0=None,solution_type='general'):
    if solution_type=='general':
        mu = exp(integrate(((1-n) * g_x) , x))
        v = ((integrate((mu* (1-n) * q_x) ,x)) + c)/mu
        v = v**(1/(1-n))
        sol = simplify(v)
        print("y = " , sol)
        return sol
    
    elif solution_type=='particular':
        mu = exp(integrate((g_x * (1-n)), x))  # integrating factor
        general_sol = (integrate((mu * q_x *(1-n)), x) + c) / mu
        particular_eq = general_sol.subs(x, x0) - y0
        c_val = solve(particular_eq, c)[0]
        final_sol = simplify(general_sol.subs(c, c_val))
        return final_sol
    


def second_order_homogenous(a,b,c,x0=None , y0=None ,x1=None ,y1=None ,solution_type='general'):
    D = b**2 - 4*a*c
    if D > 0:
        m1 = (-b + sqrt(D)) / (2*a)
        m2 = (-b - sqrt(D)) / (2*a)
        sol = c1 * exp(m1 * x) + c2 * exp(m2 * x)
    
    elif D == 0:
        m = -b / (2*a)
        sol = (c1 + c2 * x) * exp(m * x)

    
    else:
        alpha = -b / (2*a)
        beta = sqrt(-D) / (2*a)
        sol = exp(alpha * x) * (c1 * cos(beta * x) + c2 * sin(beta * x))
        
    if solution_type == 'general':
        return simplify(sol)

    
    elif solution_type == 'particular':
        if None in (x0, y0, x1, y1):
            raise ValueError("Initial conditions (x0, y0, x1, y1) must be provided for particular solution.")
        eq1 = Eq(sol.subs(x, x0), y0)
        eq2 = Eq(sol.subs(x, x1), y1)
        constants = solve((eq1, eq2), (c1, c2))
        sol = sol.subs(constants)
        return simplify(sol)
    



             
                
            
            
            
        
       




            





    













        





