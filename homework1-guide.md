# Homework 1 Grading Explanations
 - APC 524
 - Fall 2019
 - Princeton University
 - by Dev Dabke & Emily Walters

## Introduction
This document provides a brief overview to some common errors we saw while grading assignment 1.
Additionally, we've made some notes to explain parts of the rubric in more detail.
As we came across them, we also made note of snippets of code that are paragon examples of good code and of these common errors.

Since the rubric has two main components (design and style), we've broken down this guide into those two components.

## Design

### Mistake D1: Code Repetition
The first element of design that we considered was the level of code repetition.
One particular place we saw code repetition was in the ab2 file, where we tended to see something like this:
```c
  /* For the first time step */
  if (t == 0) {
    /* Run the RHS function */
    int rhserr = integrator->rhs(n, t, x, fx);
    if (rhserr != 0)
    {
      return rhserr;
    }

    /* Now run forward Euler */
    for (int i = 0; i < n; ++i)
    {
      x[i] += integrator->dt * fx[i];
    }

  } else { /* For the rest of the time steps */
    /* Run the RHS function */
    int rhserr = integrator->rhs(n, t, x, fx);
    if (rhserr != 0)
    {
      return rhserr;
    }
    
    /* Adams-Bashforth method */
    for (int i = 0; i < n; ++i) {

      x[i] += (3/2)*(integrator->dt)*fx[i] - (1/2)*(integrator->dt)*fxp[i];

    }

  }
```
However, this code exemplifies "structural repetition" in that we see a lot of redundancy in code patterns:
 - the conditional check on the return code of `rhserr`
 - the loop
This code can be condensed to
```c
  /* Bubble up any errors from RHS function */
  int rhserr = integrator->rhs(n, t, x, fx);
  if (rhserr != 0)
  {
    return rhserr;
  }

  /* Check if we are at the first step.
   * If so, we calculate this x with Euler */
  if (integrator->initial_step == 0)
  {
    /* Forward Euler algorithm for dx */
    for (int i = 0; i < n; ++i)
    {
      x[i] += dt * fx[i];
      integrator->fx0[i] = fx[i];
    }

    integrator->initial_step = 1;

    /* Successful exit */
    return 0;
  }

  /* Adams-Bashworth algorithm for dx */
  for (int i = 0; i < n; ++i)
  {
    x[i] += 1.5*dt * fx[i] - 0.5*dt * integrator->fx0[i];
    integrator->fx0[i] = fx[i];
  }
```
to avoid the repetition of checking `rhserr`.

Or, to avoid the repetition of the loop:
```c
  /* Bubble up any errors from RHS function */
  int rhserr = integrator->rhs(n, t, x, fx);
  if (rhserr != 0)
  {
    return rhserr;
  }

  /* ab2 algorithm for dx */
  for (int i = 0; i < n; ++i)
  {
    x[i] += (3/2) * integrator->dt * fx[i];
  
    if(fx[i] != 0)
    {
      x[i] += x[i] - (1/2) * integrator->dt *fx[i];
    }
  }
```

### Mistake D2: Lack of Modularity
Sometimes, repetition comes from a lack of "modularity" or in other words, a function is doing too many things at once.
There are (at least) two benefits to using subfunctions:
 - you can reuse these subfunctions
 - they provide a nice label to a section of your code

Here is an example of an `rk4` implemention that involves a lot of repetition, but its root cause is a lack of modularity
```c
int integrator_step(Integrator *integrator, double t, double *x)
{
  assert(integrator);
  const int n = integrator->n;	/* Shorthand `n` for use
				   in-function */
  double fx[n];			/* NB: This is a VLA!! */

  /* Defining variables to be used with RK */
  double fxk2[n];
  double fxk3[n];
  double fxk4[n];

  double xk2[n];
  double xk3[n];
  double xk4[n];

  double k1[n];
  double k2[n];
  double k3[n];
  double k4[n];


  /* Bubble up any errors from RHS function */
  int rhserr = integrator->rhs(n, t, x, fx);
  if (rhserr != 0)
  {
    return rhserr;
  }

  /* Runge-Kutta algorithm for dx */
  for (int i = 0; i < n; ++i)
  {
    /*k1*/
    k1[i] = integrator->dt * fx[i];

    /*k2*/
    xk2[i] = x[i] + k1[i]/2;
  }

/* calculate fx for k2 value */
  double tk2 = t + (integrator->dt / 2);
  int rhserr2 = integrator->rhs(n, tk2, xk2, fxk2);

/* calculate k2 */
  for (int i = 0; i < n; ++i)
  {
    k2[i] = integrator->dt * fxk2[i];
  }

/* calculate xk3 used to find k3 */
  for (int i = 0; i < n; ++i)
  {
    xk3[i] = x[i] + k2[i]/2;
  }

/* calculate fx used to calculate k3 */
  double tk3 = t + (integrator->dt / 2);
  int rhserr3 = integrator->rhs(n, tk3, xk3, fxk3);

/* calculate k3 */
  for (int i = 0; i < n; ++i)
  {
    k3[i] = integrator->dt * fxk3[i];
  }

/* calculate x used to calculate k4 */
  for (int i = 0; i < n; ++i)
  {
    xk4[i] = x[i] + k3[i];
  }

/* calculate fx used to find k4 */
  double tk4 = t + (integrator->dt);
  int rhserr4 = integrator->rhs(n, tk4, xk4, fxk4);

/* calculate k4 */
  for (int i = 0; i < n; ++i)
  {
    k4[i] = integrator->dt * fxk4[i];
  }

  /*final (wooo!)*/
  for (int i = 0; i < n; ++i)
  {
    x[i] += (k1[i]/6 + (k2[i]+k3[i])/3 + k4[i]/6);
  }

  /* Successful exit */
  return 0;


}
```

Contrast this with
```c
  integrator->rhs(n, t, x, k1);
  k_vector_product(n, dt, k1);
  x_vector_sum(n, x1, k1, mul_h);
```
Although this is just a snippet, we see that by delegating to helper functions `k_vector_product` and `x_vector_sum`, the code is both more readable and more terse.

### Mistake D3: Lack of Error Checking
Code may not always function as expected.
It is important to include basic error checks as code is written, saving yourself time and frustration in the future.
By adding a few assertions to check for malformed input values, `NULL` pointers, or other unexpected program states, you can catch bugs much more easily and quickly.
For example, a simple check that the user passed an appropriate set of command line arguments is a good idea.
```c
  int main( int argc, char *argv[] )
  {
    if (argc != 3) /* user did not pass correct number of command line arguments */
    {
      /* explain correct command line syntax */
      printf( "to use, use the bash syntax %s omega N_steps\n",argv[0] );
      printf( "omega is the desired frequency for the damped driven oscillator\n");
      printf( "Nsteps is the number of iteration steps desired\n");
    }
    else
    ...
```
Additionally writing code in a way (say by keeping your interfaces clean) helps prevent possible sources of error.


## Style

### Mistake S1: Magic Numbers, etc.
Code is written for humans.
For this reason, we want to make code as understandable as possible.
Therefore, we want to avoid things like "magic" numbers
```c
  /* Runge_Kutta methods for dx */
  for (int i = 0; i < n; ++i)
  {
    x[i] += k1[i] / 6 + k2[i] / 3 + k3[i] / 3 + k4[i] / 6;
  }
```
Immediately, looking at this code: what is `6` and what is `3`?
Where do these come from?
Although, this particular example comes from rk4, which perhaps the reader is more familiar with.

A more egregious example:
```c
  dt = 18.849555921538 / Sn;

  for (i=0; i < Sn; ++i)
  {
    double xold = x;
    double dxold = dx;
    
    x += dt * dxold;
    dx += dt * (cos(Fq*t) - xold - 0.5*dxold);
    t += dt;
    printf("%15.8f %15.8f %15.8f\n", t, x, dx);
  }
  printf("\n");
```
What is `18.849555921538` supposed to mean?
Compare this code with the following:
```c
#define pi 3.14159265358979323846
...
const double    T = 6*pi;							// Total time 
```
Now, we can more semantically understand where this number is coming from.

More generally, pick variable and function names that are easy to understand.
A reader should not have to guess as to what a variable represents.
More pertinently, avoid "magic" numbers and constants.
One way to think of a variable is that it acts like a label on your code.

### Mistake S2: Lacking Well-Documented Intent (and Poor Readbility)
Although you may have good code design, it is important for your code to be readable.
Moreover, we want to clearly understand the intent and purpose of every single line of code without thinking too hard about it.
```c
  frequency= atof(argv[1]);
  Nsteps= atoi(argv[2]);

  double *x0 = malloc(Nsteps*sizeof(double));
  double *x1 = malloc(Nsteps*sizeof(double));

  /*dt is time total divided by number of steps*/
  dt = totaltime / Nsteps;
  printf("\n dt is %lf", dt);
  printf("\nForcing frequency is %lf and Nsteps is %d \n", frequency, Nsteps); 
```
In the above example, the code is fine, but does a lot of thinkings one after another.
It's not obvious why we jump from something to do with parsing arguments, then initializing some other variables, and then printing out some computed value.
Organized code should read like an essay: we should almost be able to guess what you are going to say next.

The way to document intent and organize code is both through the structure of the code and the strategic use of comments.
The following code block is organized and uses minimal comments to help clarify what is happening.
```c
  	//Variables
 	double t[Nsteps+1];         						// Arrays to store results 
  	double x[Nsteps+1];       	
	double dx[Nsteps+1];      

	// Initial Conditions  
	t[0] = 0.0;											
	x[0] = 0.0;
	dx[0] = 0.0;

	// Solve & Step
  	stepper(t, dt, x, dx, omg, Nsteps);

  	// Print
  	printer(Nsteps, t, x, dx);
```
**Beware**: comments are not always good and they can promote clutter.
Good variable names and good structure can help obviate the need for comments:
```c
/* Function that verify the number of arguments. */
void verify_arguments(int argc);
```
In the above snippet, we have a fairly well-named function which makes the comment above it both unhelpful and cluttering.

## References
All of these code snippets are from real assignments.
We have not greatly modified them (only perhaps cut out subsections).
