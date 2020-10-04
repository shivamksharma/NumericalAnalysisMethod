# **Numerical Analysis in Six Methods**

In This Repository, I'm gonna show you'all on two famous equations in Mathematics which are of huge importance in almost every domain of Science F(x) = 0 & F(x) = x. Some of Equations Problems in Mathematics, Chemistry, Physics are reduced to solving one of these two quations & thats why Mathematicians/Scientists over time have tried hard to offer good solutions to resolve these equations. We'll go thorugh some of Numerical Methods that try to solve them & We'll try to make a short analysis of those Methods & also showing the advantages or disadvantages of every one of them.



Methods Whhich We'll Discuss On :

- Newton's Method
- Secant Method
- Successive Approximations Method
- Aitken's Method
- Steffenson's Method
- Overholt's Method



# **Newton's Method**

Starts with an initial value for the solution then we replace the function by its tangent & we compute the root of this tangent which will be a better approximation for the function's root. I repeat this process until we find a suitable solution (one that is close enough to the actual solution and fits very well the equation F(x) = 0). It is obvious that this process is, in fact, an iterative process. Note also that the function F must be a real valued, differentiable function in order to apply Newton's algorithm.

If We've a current approximation xCrt, The next appoximation nNxt will be computed using the following formula down below.

```
 xNxt = xCrt - (F(xCrt) / F`(xCrt))
```
F` Denotes the Derivatice of the function F & The Iteration process stops when we've gone through a maximum permitted number of iteration & we still can't find the solution or when we've found an appoximation which is close enough to the actual solution of the equation.

Code Down Below which Implements The Newton's Method:

```C++
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
```
In the above code snippet, Fd denotes the Derivative of the Function F.



# **Secant Method**

Same approach for solving the quation F(x) = 0. This method is almost identical with Newton's Method except the fact that we choose two initial approximations instead ofone before we start the Iteration Process. Suppose We've the current approximations xCrt0 and xCrt1. The next approximation xNxt will be computed this time using the following formula down below.

```
xNxt = xCrt1 - (F(xCrt1)(xCrt1 - xCrt0)) / (F(xCrt1) - F(xCrt0))
```

This Methos doesn't require the Derivatice of the Funtion F like Newton's Method. Code Down Below which Implements The Secant Method.

```
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
```

Another problem that comes into attention some times is solving the equation F(x) = x. If we write the equation like this: F(x) - x = 0 and we note G(x) = F(x) - x, then the equation becomes G(x) = 0. But the equation in the form F(x) = x presents a particular interest for mathematicians. It is said that if x0 is a solution of the equation F(x) = x, then x0 is called a fixed point of the function F(x). Of course, we can apply the methods learned before for the equation G(x) = 0, but our interest is to present methods for solving the equation F(x) = x.



# **Successive Approximations Method**

This method is as simple as it may be, is of huge importance in Mathematics, being widely used in many fixed point theories. Let's see how the method works. First, like before, we choose an initial approximation x0, and we start the iterative process. If xCrt denotes the current approximation. We Compute xNxt Like as This:

```
xNxt = F(xCrt)
```
This is a pretty simple formula & against all odds, it has been proved that the method usually converges after a number of iterations, leading to a good approximation of the equation solution & the source code is pretty straightforward.

```C++
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
```

This method, as simple and useful as it may be, has one big disadvantage: its convergence rate is very slow, meaning that the number of iterations passed until we get a solution will be pretty high. If the initial approximation is chosen close to the actual solution of the equation, the iteration process will be fast enough and the algorithm will find a suitable solution. But, if the initial approximation is chosen at random, the process may not find a solution at all (depending on the number of maximum iterations permitted, and how close to the actual solution the initial approximation is).

To overcome this problem, some mathematicians tried to speed up the process of finding the right solution for the iterative methods. So some algorithms were developed to do just that (those algorithms are usually called convergence acceleration algorithms). The most important algorithms of convergence acceleration are Aitken's algorithm and Overholt's algorithm.

The idea behind a convergence acceleration algorithm is the following: if we look closer at the iteration process, we see that if we choose every approximation value at each step, we can form a real valued number sequence: x0, x1, x2, ... ,xn which, after a number of iterations, converges to the solution (x) of the equation F(x) = x. A convergence acceleration algorithm transforms the number sequence x0, x1, x2, ... ,xn into another real valued number sequence y0, y1, y2, ... yn, which has the very important property that it converges faster to the solution x.

We present in the following section, three such algorithms for convergence acceleration, and we stick to our problem of solving the equation F(x) = x.



# **Aitken's Method**

Aitken's method is an iterative process similar to the ones presented. I won't go into the mathematical details again. Instead, I prefer to present the code for Aitken's method.

```C++
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
```



# **Steffenson's Method**

This method is a simplified version of Aitken's method, observing that if we apply Aitken's formula for the values xCrt, F(xCrt), F(F(xCrt)), we obtain:

```
xNxt = (xCrt F(F(xCrt)) - (F(xCrt)) ^ 2) / (F(F(xCrt)) - 2F(xCrt) + xCrt)
```

In a Simplified form, This is written as:

```
xCrt = F(xCrt) + 1 / ((1 / (F(F(xCrt)) - F(xCrt))) - (1 / (F(xCrt) - xCrt)))
```

The Code for The Algorithm is:

```C++
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
```



# **Overholt's Method**

Overholt's method is the fastest method for solving the equation F(x) = x. I'm going to present the code for the algorithm, the method itself being pretty straightforward from the source code:

```C++
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
```

In the source code above, V is a global variable declared as:

```
// matrix used in Overholt algorithm
// this is declared as a global variable to avoid stack overflow
    double V[100][100];
```




## **USING THE CODE**

Now, let's see how to effectively use the methods described above, and on the way, we will make a short analysis of these methods.

In order to use the methods, we must first define some constants:

```C++
// the maximum number of iterations
const int MAXITER    = 100;

// the accepted error
const double error    = 0.0001;
```

MAXITER represents the maximum number of iterations an algorithm is permitted to pass, meaning that if one of our algorithms has passed MAXITER iterations and it still didn't find a suitable solution, the algorithm will stop nevertheless. error represents the minimum accepted error for the solution, meaning that the approximate solution must be close enough to the real solution of the equation with respect to this error. The smaller the value of this error, the closer to the real solution the approximate solution will be.

Let's take an example equation of the first kind, for instance, x*e^x - 1 = 0. So, we have F(x) = x*e^x - 1. We define the function in our code, and also its derivative, because we will need it in order to apply Newton's method:

```C++
// the Euler constant
const double e        = 2.718281828459;

// the function F(x)
#define F(x) (  x * pow(e, x) - 1  )

// the derivative of the function F(x), meaning F`(x)
#define Fd(x) (  (x + 1) * pow(e, x)  )
```

You can change the function F & its derivative Fd & you'll solve any kind of equation you like.

```C++
double x;
int n;

/* Example usage for Newton Method */
cout << "Newton's method: " << endl << endl;

cout << "Give the initial approximation: ";
cin >> x;


// now we apply Newton's method
n = NewtonMethodForEquation(x);

if(n > MAXITER)
    cout << "In " << MAXITER << " iterations no solution was found!" << endl;
else
    cout << "The solution is: " << x << " and it was found in " 
         << n << " iterations" << endl;

double x0, x1;

/* Example usage for Secant Method */
cout << "Secant method: " << endl << endl;

cout << "Give the first initial approximation: ";
cin >> x0;

cout << "Give the second initial approximation: ";
cin >> x1;

// now we apply the Secant method
n = SecantMethodForEquation(x, x0, x1);

if(n > MAXITER)
    cout << "In " << MAXITER << " iterations no solution was found!" << endl;
else
    cout << "The solution is: " << x << " and it was found in " 
         << n << " iterations" << endl;
```

Now you can play with the algorithms, giving various initial approximations and decreasing or increasing the error value to see how they are behaving. We can see that, overall, Newton's method is faster than Secant method.

For the second type of equation, let's take, for example, F(x) = e^(-x), and the equation becomes e^(-x) = x.

```C++
#define F1(x) ( pow(e, -x) )
```

And Here is how we apply our algorithms:

```C++
/* Example usage for Successive Approximations Method */
cout << "Successive approximations method: " << endl << endl;

cout << "Give the initial approximation: ";
cin >> x;

n = SuccessiveApproxForEquation(x);

if(n > MAXITER)
    cout << "In " << MAXITER << " iterations no solution was found!" << endl;
else
    cout << "The solution is: " << x << " and it was found in " 
         << n << " iterations" << endl;

/* Example usage for Aitken Method */
cout << "Aitken's method: " << endl << endl;

cout << "Give the initial approximation: ";
cin >> x;

double y;

n = AitkenMethodForEquation(y, x);

if(n > MAXITER)
    cout << "In " << MAXITER 
         << " iterations no solution was found!" << endl;
else
    cout << "The solution is: " << y << " and it was found in " 
         << n << " iterations" << endl;


/* Example usage for Steffensen Method */
cout << "Steffensen's method: " << endl << endl;

cout << "Give the initial approximation: ";
cin >> x;

n = SteffensenMethodForEquation(x);

if(n > MAXITER)
    cout << "In " << MAXITER 
         << " iterations no solution was found!" << endl;
else
    cout << "The solution is: " << x << " and it was found in " 
         << n << " iterations" << endl;

/* Example usage for Overholt Method */
cout << "Overholt's method: " << endl << endl;

cout << "Give the initial approximation: ";
cin >> x;

double rez;

n = OverholtMethodForEquation(rez, x, 2);

if(n > MAXITER)
    cout << "In " << MAXITER 
         << " iterations no solution was found!" << endl;
else
    cout << "The solution is: " << rez << " and it was found in " 
         << n << " iterations" << endl;
```

Try and play with these algorithms too, and you will see that, indeed, Overholt's method is a stable and fast method for solving equations of type F(x) = x.





