# NumericalAnalysisMethod


INTRODUCTION

This Respository tries to familiarize the beginner with numerical methods. I am working a lot with numerical analysis and methods, and I want to share with you some of my experiences and the results that I encountered. This is intended to be the first article in a series of Numerical Analysis Methods and Their Implementation in C++.

In this Respository, we are going to focus on two famous equations in Mathematics, which are of huge importance in almost every domain of science: F(x) = 0 and F(x) = x. Many problems in Mathematics, Chemistry, Physics are reduced to solving one of these two equations, and that's why scientists over time have tried hard to offer good solutions to resolve these equations. We will go through some numerical methods that try to solve them, and we will try to make a short analysis of those methods, showing the advantages/disadvantages of every one of them.
The methods and their implementation

The first issue of our article is the problem of solving the equation F(x) = 0, where F(x) can be any kind of function. For instance, we can have F(x) = x^2 + 6x + 6, or F(x) = cos(x). This is a well known problem in Mathematics, and it is known also that there is no way of finding out exactly the solution or solutions of every equation of this form. That's why mathematicians have been trying to offer approximate solutions for this equation. Please also note that the solutions of the equation F(x) = 0 are usually called the roots of the function F(x).

#Newton's method

Newton's approach is the following: we start with an initial value for the solution (also called initial approximation), then we replace the function by its tangent, and we compute the root of this tangent which will be a better approximation for the function's root. We repeat this process until we find a suitable solution (one that is close enough to the actual solution and fits very well the equation F(x) = 0). It is obvious that this process is, in fact, an iterative process. Note also that the function F must be a real valued, differentiable function in order to apply Newton's algorithm.

In detail, if we have a current approximation xCrt, the next approximation nNxt will be computed using the following formula:

xNxt = xCrt - (F(xCrt) / F`(xCrt))

where F denotes the derivative of the function F. The iteration process stops when we have gone through a maximum permitted number of iterations and we still can't find the solution, or when we have found an approximation which is close enough to the actual solution of the equation. Here is the code that implements 

Newton's method:

---------------------------------------------------------------------------------------------------------------------------------
/* 
 * Newton's method for solving equation F(x) = 0
 * Output: 
 * x - the resulted approximation of the solution
 * Return:
 * The number of iterations passed
 */
int NewtonMethodForEquation(double& x)
{
    int n = 1;

    while( ( fabs(F(x)) > error ) && ( n <= MAXITER ) )
    {
        x = x - ( F(x) / Fd(x) );

        n++;
    }

    return n;
}

In the above code snippet, Fd denotes the derivative of the function F.

#Secant method

The secant method is another approach for solving the equation F(x) = 0. The method is almost identical with Newton's method, except the fact that we choose two initial approximations instead of one before we start the iteration process. Suppose we have the current approximations xCrt0 and xCrt1. The next approximation xNxt will be computed this time using the following formula:

xNxt = xCrt1 - (F(xCrt1)(xCrt1 - xCrt0)) / (F(xCrt1) - F(xCrt0))

Note that this method doesn't require the derivative of the function F, like Newton's method did. Here is the code that implements the Secant method.

/* 
 * Secant method for solving equation F(x) = 0
 * Input:
 * x0 - the first initial approximation of the solution
 * x1 - the second initial approximation of the solution
 * Output:
 * x - the resulted approximation of the solution
 * Return:
 * The number of iterations passed
 */
int SecantMethodForEquation(double& x, double x0, double x1)
{
    int n = 2;

    while( ( fabs(F(x1)) > error ) && ( n <= MAXITER ) )
    {
        x = x1 - (F(x1) * (x1 - x0)) / (F(x1) - F(x0));
        x0 = x1;
        x1 = x;

        n++;
    }

    return n;
}

Another problem that comes into attention some times is solving the equation F(x) = x. If we write the equation like this: F(x) - x = 0 and we note G(x) = F(x) - x, then the equation becomes G(x) = 0. But the equation in the form F(x) = x presents a particular interest for mathematicians. It is said that if x0 is a solution of the equation F(x) = x, then x0 is called a fixed point of the function F(x). Of course, we can apply the methods learned before for the equation G(x) = 0, but our interest is to present methods for solving the equation F(x) = x.

#Successive approximations method

This method, as simple as it may be, is of huge importance in Mathematics, being widely used in many fixed point theories. Let's see how the method works. First, like before, we choose an initial approximation x0, and we start the iterative process. If xCrt denotes the current approximation, then we compute xNxt like this:

xNxt = F(xCrt)

This is a pretty simple formula, and against all odds, it has been proved that the method usually converges after a number of iterations, leading to a good approximation of the equation solution. The source code is pretty straightforward:

/* 
 * Successive approximations method for solving equation F(x) = x
 * Output:
 * x - the resulted approximation of the solution
 * Return:
 * The number of iterations passed
 */
int SuccessiveApproxForEquation(double& x)
{
    int n = 1;

    while( ( fabs(x - F1(x)) > error ) && ( n <= MAXITER ))
    {
        x = F1(x);

        n++;
    }

    return n;
}

This method, as simple and useful as it may be, has one big disadvantage: its convergence rate is very slow, meaning that the number of iterations passed until we get a solution will be pretty high. If the initial approximation is chosen close to the actual solution of the equation, the iteration process will be fast enough and the algorithm will find a suitable solution. But, if the initial approximation is chosen at random, the process may not find a solution at all (depending on the number of maximum iterations permitted, and how close to the actual solution the initial approximation is).

To overcome this problem, some mathematicians tried to speed up the process of finding the right solution for the iterative methods. So some algorithms were developed to do just that (those algorithms are usually called convergence acceleration algorithms). The most important algorithms of convergence acceleration are Aitken's algorithm and Overholt's algorithm.

The idea behind a convergence acceleration algorithm is the following: if we look closer at the iteration process, we see that if we choose every approximation value at each step, we can form a real valued number sequence: x0, x1, x2, ... ,xn which, after a number of iterations, converges to the solution (x) of the equation F(x) = x. A convergence acceleration algorithm transforms the number sequence x0, x1, x2, ... ,xn into another real valued number sequence y0, y1, y2, ... yn, which has the very important property that it converges faster to the solution x.

We present in the following section, three such algorithms for convergence acceleration, and we stick to our problem of solving the equation F(x) = x.

#Aitken's method

Aitken's method is an iterative process similar to the ones presented. I will not go into the mathematical details again. Instead, I prefer to present the code for Aitken's method:
Hide   Shrink Image 3 for Some simple numerical methods in C++   Copy Code

/* 
 * Aitken's method for solving equation F(x) = x
 * Input:
 * x0 - the initial approximation of the solution
 * Output:
 * y - the resulted approximation of the solution
 * Return:
 * The number of iterations passed
 */
int AitkenMethodForEquation(double& y, double x0)
{
    int n = 1;
    double x;

    do
    {
        x = F1(x0);

        y = x + 1 / ((1 / (F1(x) - x)) - (1 / (x - x0)));

        n++;

        x0 = x;
    }
    while((fabs(y - F1(y)) > error) && (n <= MAXITER));

    return n;
}

---------------------------------------------------------------------------------------------------------------------------------
#Steffenson's method

This method is a simplified version of Aitken's method, observing that if we apply Aitken's formula for the values xCrt, F(xCrt), F(F(xCrt)), we obtain:

xNxt = (xCrt F(F(xCrt)) - (F(xCrt)) ^ 2) / (F(F(xCrt)) - 2F(xCrt) + xCrt)

xCrt = F(xCrt) + 1 / ((1 / (F(F(xCrt)) - F(xCrt))) - (1 / (F(xCrt) - xCrt)))

/* 
 * Steffensen's method for solving equation F(x) = x
 * Output:
 * x - the resulted approximation of the solution
 * Return:
 * The number of iterations passed
 */
int SteffensenMethodForEquation(double& x)
{
    int n = 0;

    do
    {
        x = F1(x) + 1 / ( (1 / (F1(F1(x)) - F1(x)) ) - (1 / (F1(x) - x) ) );

        n++;
    }
    while((fabs(x - F1(x)) > error) && (n <= MAXITER));

    return n;
}

#Overholt's Method

Overholt's method is the fastest method for solving the equation F(x) = x. I'm going to present the code for the algorithm, the method itself being pretty straightforward from the source code.

/* 
 * Overholt's method for solving equation F(x) = x
 * Input:
 * x0 - the initial approximation of the solution
 * s - the convergence order (s must be  >= 2)
 * Output:
 * x - the resulted approximation of the solution
 * Return:
 * The number of iterations passed
 */
int OverholtMethodForEquation(double &x, double x0, int s)
{
    int m = 0;
    
    x = x0;

    do
    {
        V[0][0] = x;

        for(int n = 1; n <= s; n++)
            V[0][n] = F1(V[0][n - 1]);

        for(int k = 0; k <= s - 2; k++)
            for(int n = 0; n <= s - k - 2; n++)
                V[k + 1][n] = (pow(V[0][n + k + 2] - V[0][n + k + 1], k+1) * V[k][n] - 
                          pow(V[0][n + k + 1] - V[0][n + k], k + 1) * V[k][n + 1]) / 
                         (pow(V[0][n + k + 2] - V[0][n + k + 1], k + 1) - 
                          pow(V[0][n + k + 1] - V[0][n + k], k + 1));

        m++;

        x = V[s - 1][0];
    }
    while((fabs(x - F1(x)) > error) && (m <= MAXITER));

    return m;
}
