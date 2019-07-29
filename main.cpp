#include <iostream>
#include <math.h>


using namespace std;


// the Euler constant
const double e		= 2.718281828459;

// the maximum number of iterations
const int MAXITER	= 100;

// the accepted error
const double error	= 0.0001;

// matrix used in Overholt algorithm
// this is declared as a global variable to avoid stack overflow
double V[100][100];


// the function F(x)
#define F(x) (  x * pow(e, x) - 1  )

// the derivative of the function F(x), meaning F`(x)
#define Fd(x) (  (x + 1) * pow(e, x)  )

// the function used in successive approximations method
#define F1(x) ( pow(e, -x) )


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


int main()
{
	double x;
	int n;


	/* Example usage for Newton Method */
	cout << "Newton's method: " << endl << endl;

	cout << "Give the initial approximation: ";
	cin >> x;


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

	n = SecantMethodForEquation(x, x0, x1);

	if(n > MAXITER)
		cout << "In " << MAXITER << " iterations no solution was found!" << endl;
	else
		cout << "The solution is: " << x << " and it was found in " 
			 << n << " iterations" << endl;


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
		cout << "In " << MAXITER << " iterations no solution was found!" << endl;
	else
		cout << "The solution is: " << y << " and it was found in " 
			 << n << " iterations" << endl;


	/* Example usage for Steffensen Method */
	cout << "Steffensen's method: " << endl << endl;

	cout << "Give the initial approximation: ";
	cin >> x;

	n = SteffensenMethodForEquation(x);

	if(n > MAXITER)
		cout << "In " << MAXITER << " iterations no solution was found!" << endl;
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
		cout << "In " << MAXITER << " iterations no solution was found!" << endl;
	else
		cout << "The solution is: " << rez << " and it was found in " 
			 << n << " iterations" << endl;


	return 0;
}
