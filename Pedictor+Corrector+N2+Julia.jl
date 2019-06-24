using LinearAlgebra
using PolynomialRoots
import Distances
using Plots


setprecision(1000)

function LW_ineq(r,t)

    A = zeros(BigFloat, 3r+1,2r)
    b = zeros(BigFloat, 3r+1,1)
    c = zeros(BigFloat, 2r,1)
# defining A
    A[1,1], A[2,2] = 1 , 1
    A[3r,2r-1], A[3r+1,2r] = -1 , -1

    for j = 1:r-1
        A[3j,2j-1] ,   A[3j,2j+1] =  -t , 1
        A[3j+1 ,2j],   A[3j+1 ,2j+1] =  -t , 1
        A[3j+2 ,2j-1], A[3j+2,2j] , A[3j+2 ,2j+2] =  -t^(1 - 2.0^(-j)) , -t^(1 - 2.0 ^(-j)) , 1
    end
#defining b

    b[1],b[2] = t^2.0, t
#defining c
    c[1] = 1

    return (A,b,c)
end




function LW_eq(r,t)

    A,b,c = LW_ineq(r,t)

    big_A = zeros(BigFloat, 3r+1,5r+1)
    big_c = zeros(BigFloat, 5r+1,1)

    big_A[1:3r+1,1:2r] = A
    big_A[1:3r+1,2r+1:5r+1] = Matrix{BigFloat}(I, 3r+1, 3r+1)

    big_c[1:2r] = c

    return (big_A,b,big_c)

end


function tropical_central_path_x_y(r::Int64,lambda::BigFloat)

    x = zeros(BigFloat,2r)
    y = zeros(BigFloat,3r + 1)

    x[1], x[2]  = min(2,lambda) , 1
    y[1], y[2]  = lambda - 2, lambda - 1


    for j = 1:r-1

        x[2j+1] = 1 + min(x[2j-1] , x[2j])
        x[2j+2] = (1 - 2.0^(-j)) + max(x[2j-1],x[2j])

        y[3j] = lambda - 1 - x[2j-1]
        y[3j+1] = lambda - 1 - x[2j]
        y[3j+2] = lambda - x[2j + 2]
    end

    y[3r], y[3r + 1] = lambda - x[2r-1], lambda - x[2r]

    return x , y

end

function mu_bar(z,n::Int64,m::Int64)

    x,s = z[1:n], z[n+m+1:2n + m]

    return dot(x,s)/n

end

function Newton_matrix(n::Int64,m::Int64,A,z)

    X , S = Diagonal(z[1:n]),  Diagonal(z[n+m+1:2n + m] );

    M = zeros(BigFloat,2n+m, 2n+m);

    M[1:n, n+1:n+m] = transpose(A)
    M[1:n, n+m+1: 2n + m] = Matrix{BigFloat}(I, n, n)

    M[n+1:n+m, 1:n] = A

    M[m+n+1:2n+m, 1:n] = S
    M[m+n+1:2n+m, m+n+1: 2n+m] = X

    return M

end

function get_Newton_direction(n,m,A,z,mu)

    x,s = z[1:n], z[n+m+1:2n + m]
    M = Newton_matrix(n,m,A,z)
    X = Diagonal(x)

    v = zeros(BigFloat,2n+m)
    v[n+m+1:2n+m] = mu*ones(BigFloat,n,1) - X*s

    d = M\v

    return d
end

function get_max_step(n,m,z,dz,theta_p)

    dx, ds = dz[1:n], dz[n+m+1:2n + m]
    dX, dS = Diagonal(dx), Diagonal(ds)

    x, s = z[1:n], z[n+m+1:2n + m]
    X, S = Diagonal(x), Diagonal(s)

    v_1 = X*s - mu_bar(z,n,m)*ones(BigFloat,n,1)
    v_2 = X*ds + dX*s  - ( (dot(x,ds) + dot(s,dx))/n )*ones(BigFloat,n,1)
    v_3 = dX*ds - mu_bar(dz,n,m)*ones(BigFloat,n,1)

    a_0 = dot(v_1,v_1) - (theta_p*mu_bar(z,n,m))^2
    a_1 = 2*( dot(v_1,v_2) - (theta_p^2) * mu_bar(z,n,m) * ( (dot(x,ds) + dot(s,dx))/n ) )
    a_2 = dot(v_2,v_2)+ 2*dot(v_1,v_3) - (theta_p^2)* (2*mu_bar(z,n,m)*mu_bar(dz,n,m) +  ( (dot(x,ds) + dot(s,dx))/n  )^2 )
    a_3 = 2*dot(v_2,v_3)   - 2 * theta_p * mu_bar(dz,n,m) * ((dot(x,ds) + dot(s,dx))/n)
    a_4 = dot(v_3,v_3) - (theta_p*mu_bar(dz,n,m))^2

    P = [a_0,a_1,a_2,a_3,a_4]

    alpha = - Inf
    bool = false

    for a in PolynomialRoots.roots(P)
        if abs(Complex(a).im) <= 10.0^(-10)
            bool = true
            alpha = max(alpha, Complex(a).re)
        end
    end

    if bool == false
        print("Roots are:", roots(P))
        throw("No real soution found")
    else
        return alpha
    end

end


function lifting_coefs(r,l_alpha,l_beta)

    alpha = zeros(BigFloat,2r)
    beta = zeros(BigFloat,3r+1)

    alpha[2r-1:2r] = [1/2 - l_alpha^(-7), 1/2 - l_alpha^(-7)]
    beta[1:2] = [l_beta^r,l_beta^r]


    for j = 1:r-1

        alpha[2j-1:2j] = [1/2 - l_alpha^(-r+j-7) , 1/2 - l_alpha^(-r+j-7)]
        beta[3j:3j+2] =[l_beta^(r-j),l_beta^(r-j),l_beta^(r-j)]

    end

    beta[3r:3r+1] = [1,1]

    beta = beta/(l_beta^(r+1))

    return alpha, beta

end


function lifter(r,t,lambda,l_alpha, l_beta)

    alpha, beta = lifting_coefs(r,l_alpha,l_beta)
    A,b,c = LW_ineq(r,t)

    x_trop, y_trop = tropical_central_path_x_y(r,lambda)

    f(trop_element) = t ^ trop_element

    x = broadcast(*,alpha,  map(f, x_trop) )
    y = broadcast(*,beta ,  map(f, y_trop) )
    w = b - A*x
    s = transpose(A)*y + c

    z = Array{BigFloat}(undef,13r+3,1)

    z[1:2r] = x
    z[2r+1:5r+1] = w
    z[5r+2:8r+2] = -y
    z[8r+3:10r+2] = s
    z[10r+3:13r+3] = y

    return z
end



function initializer(n,m,z_0,A,theta)

    z = copy(z_0)
    x, s = z[1:n], z[n+m+1:2n + m]
    X, S = Diagonal(x), Diagonal(s)

    mu = mu_bar(z,n,m)
    cr = Distances.euclidean(X*s , mu*ones(BigFloat,n,1) ) / mu


    while cr >= theta

        dz = get_Newton_direction(n,m,A,z,mu)
        z = z + dz

        x, s = z[1:n], z[n+m+1:2n + m]
        X, S = Diagonal(x), Diagonal(s)
        mu = mu_bar(z,n,m)
        cr = Distances.euclidean(X*s , mu*ones(BigFloat,n,1) ) / mu


    end

    return z

end


function solve_until_mu_1(n,m,A,z,theta,theta_p)

    predictions = []
    corrections = []

    iteration = 0

    z_p = copy(z)

    mu = mu_bar(z_p,n,m)

    while mu > 0.9
        iteration += 1

        #### Prediction ####
        mu = 0.0
        d = get_Newton_direction(n,m,A,z_p,mu)
        alpha = get_max_step(n,m,z_p,d,theta_p)
        z_p = z_p + alpha * d
        append!(predictions, [z_p[n-1:n]])


        #### Correction ####
        mu = mu_bar(z_p,n,m)
        d = get_Newton_direction(n,m,A,z_p,mu)
        z_p = z_p + d
        append!(corrections,[z_p[n-1:n]])

        #################
        mu = mu_bar(z_p,n,m)

        println("iteration number: ", iteration)

    end

    return predictions,corrections

end




r = 6
t = big(10.0^30)
lambda = big(2.0)

n = 5r+1
m = 3r+1

theta = big(0.25)
theta_p = big(0.5)

A,b,c = LW_eq(r,t)

z_0 = lifter(r,t,lambda,big(2.0),big(2.0001))
z = initializer(n,m,z_0,A,theta)


predictions, corrections = solve_until_mu_1(n,m,A,z,theta,theta_p)

q = size(predictions)[1]

Prediction_Data = Array{BigFloat}(undef,2,q)
Correction_Data = Array{BigFloat}(undef,2,q)

for i = 1:q
    Prediction_Data[1:2,i] = predictions[i]
    Correction_Data[1:2,i] = corrections[i]
end


g(elem) = log(elem)/log(t)

Prediction_Data = map(g,Prediction_Data)
Correction_Data = map(g,Correction_Data)


plotly()

plot(Correction_Data[1,1:q], Correction_Data[1,1:q])